#include "HPXsugarscape.h"

#include <hpx/hpx_init.hpp>
#include <hpx/program_options.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/async_combinators/wait_all.hpp>
#include <hpx/runtime_distributed/find_all_localities.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>

namespace po = hpx::program_options;

static std::vector<GridCell> localGrid;
static std::vector<Agent> localAgents;
static std::vector<std::pair<int, int>> targetCells;
static std::vector<char> moveWon;
static std::vector<int> growingSugarCells;
static std::vector<unsigned char> sugarIsGrowing;

static int localityId = 0;
static int localityCount = 1;
static int localStartRow = 0;
static int localEndRow = 0;
static int gridHeight = 0;
static int gridWidth = 0;

static std::vector<int> bestLocalAgentIndex;
static std::vector<int> bestAgentId;
static std::vector<int> cellStamp;
static int currentStamp = 1;

// Expand one seed into full RNG state
uint64_t Xoshiro256::splitmix64(
    uint64_t& seedValue
) {
    seedValue += 0x9e3779b97f4a7c15ULL;
    uint64_t mixedValue = seedValue;

    mixedValue = (mixedValue ^ (mixedValue >> 30)) *
        0xbf58476d1ce4e5b9ULL;
    mixedValue = (mixedValue ^ (mixedValue >> 27)) *
        0x94d049bb133111ebULL;

    return mixedValue ^ (mixedValue >> 31);
}

// Rotate bits for xoshiro state update
uint64_t Xoshiro256::rotateLeft(
    uint64_t value,
    int shift
) {
    return (value << shift) | (value >> (64 - shift));
}

// Initialize RNG state
Xoshiro256::Xoshiro256(
    uint64_t seed
) {
    for (int i = 0; i < 4; i++) {
        state[i] = splitmix64(
            seed
        );
    }
}

// Generate next RNG value
uint64_t Xoshiro256::next() {
    uint64_t result = rotateLeft(
        state[1] * 5,
        7
    ) * 9;
    uint64_t shiftedValue = state[1] << 17;

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];

    state[2] ^= shiftedValue;
    state[3] = rotateLeft(
        state[3],
        45
    );

    return result;
}

// Generate unbiased value in inclusive range
uint64_t Xoshiro256::nextInRange(
    uint64_t minValue,
    uint64_t maxValue
) {
    if (maxValue < minValue) {
        return next();
    }

    uint64_t range = maxValue - minValue + 1;
    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);
    uint64_t randomValue;

    do {
        randomValue = next();
    } while (randomValue >= limit);

    return minValue + (randomValue % range);
}

// Mark a cell for future sugar regrowth
static void markSugarGrowing(
    int localCellIndex
) {
    if (
        localCellIndex < 0 ||
        localCellIndex >= (int)localGrid.size()
    ) {
        return;
    }

    if (!sugarIsGrowing[(size_t)localCellIndex]) {
        sugarIsGrowing[(size_t)localCellIndex] = 1;
        growingSugarCells.push_back(
            localCellIndex
        );
    }
}

// Map a global row to owning locality
static int partitionForRow(
    int globalRow,
    int height,
    int localities
) {
    if (localities <= 1) {
        return 0;
    }

    int rowsPerLocality = (height + localities - 1) / localities;
    int owner = globalRow / rowsPerLocality;

    if (owner >= localities) {
        return localities - 1;
    }

    return owner;
}

// Restore stable per-locality agent order
static void sortAgentsById() {
    std::sort(
        localAgents.begin(),
        localAgents.end(),
        [](
            const Agent& leftAgent,
            const Agent& rightAgent
        ) {
            return leftAgent.id < rightAgent.id;
        }
    );
}

// Find agent index after sorting by id
static int findAgentIndex(
    int agentId
) {
    auto iterator = std::lower_bound(
        localAgents.begin(),
        localAgents.end(),
        agentId,
        [](
            const Agent& agent,
            int searchId
        ) {
            return agent.id < searchId;
        }
    );

    if (
        iterator == localAgents.end() ||
        iterator->id != agentId
    ) {
        return -1;
    }

    return (int)(iterator - localAgents.begin());
}

// Rank directions for tie breaking
static int directionIndex(
    int targetRow,
    int targetCol,
    int agentRow,
    int agentCol
) {
    int rowDelta = targetRow - agentRow;
    int colDelta = targetCol - agentCol;

    if (rowDelta < 0 && colDelta == 0) {
        return 0;
    }

    if (rowDelta > 0 && colDelta == 0) {
        return 1;
    }

    if (rowDelta == 0 && colDelta > 0) {
        return 2;
    }

    if (rowDelta == 0 && colDelta < 0) {
        return 3;
    }

    return 4;
}

// Reset conflict cell state for current timestep
static void touchConflictCell(
    int localCellIndex
) {
    if (cellStamp[(size_t)localCellIndex] != currentStamp) {
        cellStamp[(size_t)localCellIndex] = currentStamp;
        bestLocalAgentIndex[(size_t)localCellIndex] = -1;
        bestAgentId[(size_t)localCellIndex] =
            std::numeric_limits<int>::max();
    }
}

// Build sugar capacity from two peak landscape
static int sugarCapacity(
    int row,
    int col,
    int height,
    int width,
    int minCapacity,
    int maxCapacity
) {
    double firstPeakRow = height * 0.25;
    double firstPeakCol = width * 0.25;
    double secondPeakRow = height * 0.75;
    double secondPeakCol = width * 0.75;
    double maxRadius = std::min(
        (double)height,
        (double)width
    ) * 0.4;

    double firstDistance = std::sqrt(
        (row - firstPeakRow) * (row - firstPeakRow) +
        (col - firstPeakCol) * (col - firstPeakCol)
    );
    double secondDistance = std::sqrt(
        (row - secondPeakRow) * (row - secondPeakRow) +
        (col - secondPeakCol) * (col - secondPeakCol)
    );
    double distance = std::min(
        firstDistance,
        secondDistance
    );

    if (distance <= maxRadius * 0.25) {
        return maxCapacity;
    }

    if (distance <= maxRadius * 0.50) {
        return std::min(
            maxCapacity,
            std::max(
                minCapacity,
                maxCapacity - 1
            )
        );
    }

    if (distance <= maxRadius * 0.75) {
        return std::min(
            maxCapacity,
            std::max(
                minCapacity,
                maxCapacity - 2
            )
        );
    }

    return minCapacity;
}

// Initialize worker locality state
void initWorkerImpl(
    int locid,
    int nlocs,
    int rowStart,
    int rowEnd,
    std::vector<GridCell> gridChunk,
    std::vector<Agent> agents,
    int height,
    int width
) {
    localityId = locid;
    localityCount = nlocs;
    localStartRow = rowStart;
    localEndRow = rowEnd;
    gridHeight = height;
    gridWidth = width;

    localGrid = std::move(
        gridChunk
    );
    localAgents = std::move(
        agents
    );

    sugarIsGrowing.assign(
        localGrid.size(),
        0
    );
    growingSugarCells.clear();

    size_t localCellCount = localGrid.size();
    bestLocalAgentIndex.assign(
        localCellCount,
        -1
    );
    bestAgentId.assign(
        localCellCount,
        std::numeric_limits<int>::max()
    );
    cellStamp.assign(
        localCellCount,
        0
    );
    currentStamp = 1;

    sortAgentsById();
}
HPX_PLAIN_ACTION(
    initWorkerImpl,
    initWorkerAction
);

// Regrow changed sugar cells
void regrowLocalImpl(
    int growthRate
) {
    if (growingSugarCells.empty()) {
        return;
    }

    size_t writeIndex = 0;

    for (
        size_t readIndex = 0;
        readIndex < growingSugarCells.size();
        readIndex++
    ) {
        int cellIndex = growingSugarCells[readIndex];
        int& sugar = localGrid[(size_t)cellIndex].currentSugar;
        int capacity = localGrid[(size_t)cellIndex].sugarCapacity;
        int grownSugar = sugar + growthRate;

        if (grownSugar >= capacity) {
            sugar = capacity;
            sugarIsGrowing[(size_t)cellIndex] = 0;
        }
        else {
            sugar = grownSugar;
            growingSugarCells[writeIndex] = cellIndex;
            writeIndex++;
        }
    }

    growingSugarCells.resize(
        writeIndex
    );
}
HPX_PLAIN_ACTION(
    regrowLocalImpl,
    regrowLocalAction
);

// Return top or bottom sugar rows for halo reads
std::vector<int> getSugarHaloImpl(
    int fromTop,
    int rowCount
) {
    int localRows = localEndRow - localStartRow;

    if (rowCount > localRows) {
        rowCount = localRows;
    }

    std::vector<int> halo(
        (size_t)rowCount * (size_t)gridWidth,
        0
    );

    if (fromTop) {
        for (int row = 0; row < rowCount; row++) {
            for (int col = 0; col < gridWidth; col++) {
                halo[(size_t)row * gridWidth + col] =
                    localGrid[(size_t)row * gridWidth + col].currentSugar;
            }
        }
    }
    else {
        int startRow = localRows - rowCount;

        for (int row = 0; row < rowCount; row++) {
            for (int col = 0; col < gridWidth; col++) {
                halo[(size_t)row * gridWidth + col] =
                    localGrid[
                        (size_t)(startRow + row) * gridWidth + col
                    ].currentSugar;
            }
        }
    }

    return halo;
}
HPX_PLAIN_ACTION(
    getSugarHaloImpl,
    getSugarHaloAction
);

// Compute preferred movement target for each local agent
void computeTargetsLocalImpl(
    std::vector<int> topHalo,
    std::vector<int> bottomHalo
) {
    int topHaloRows = gridWidth > 0
        ? (int)(topHalo.size() / (size_t)gridWidth)
        : 0;
    int bottomHaloRows = gridWidth > 0
        ? (int)(bottomHalo.size() / (size_t)gridWidth)
        : 0;
    int topStartRow = localStartRow - topHaloRows;

    auto sugarAt = [
        &
    ](
        int row,
        int col
    ) -> int {
        if (row < localStartRow) {
            return topHalo[(size_t)(row - topStartRow) * gridWidth + col];
        }

        if (row >= localEndRow) {
            return bottomHalo[(size_t)(row - localEndRow) * gridWidth + col];
        }

        return localGrid[
            (size_t)(row - localStartRow) * gridWidth + col
        ].currentSugar;
    };

    (void)bottomHaloRows;

    int agentCount = (int)localAgents.size();
    targetCells.assign(
        (size_t)agentCount,
        {0, 0}
    );
    moveWon.assign(
        (size_t)agentCount,
        0
    );

    std::vector<std::pair<int, int>> candidates;
    std::vector<std::pair<int, int>> closestCandidates;
    candidates.reserve(
        64
    );
    closestCandidates.reserve(
        64
    );

    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        const Agent& agent = localAgents[(size_t)agentIndex];
        candidates.clear();

        int maxSugar = sugarAt(
            agent.row,
            agent.col
        );
        candidates.push_back(
            {agent.row, agent.col}
        );

        for (const auto& direction : directions) {
            for (int vision = 1; vision <= agent.vision; vision++) {
                int targetRow = agent.row + direction.first * vision;
                int targetCol = agent.col + direction.second * vision;

                if (
                    targetRow < 0 ||
                    targetRow >= gridHeight ||
                    targetCol < 0 ||
                    targetCol >= gridWidth
                ) {
                    break;
                }

                int sugar = sugarAt(
                    targetRow,
                    targetCol
                );

                if (sugar > maxSugar) {
                    maxSugar = sugar;
                    candidates.clear();
                    candidates.push_back(
                        {targetRow, targetCol}
                    );
                }
                else if (sugar == maxSugar) {
                    candidates.push_back(
                        {targetRow, targetCol}
                    );
                }
            }
        }

        int chosenRow;
        int chosenCol;

        if (candidates.size() > 1) {
            int minDistance = std::numeric_limits<int>::max();
            closestCandidates.clear();

            for (const auto& candidate : candidates) {
                int distance = std::abs(candidate.first - agent.row) +
                    std::abs(candidate.second - agent.col);

                if (distance < minDistance) {
                    minDistance = distance;
                    closestCandidates.clear();
                    closestCandidates.push_back(
                        candidate
                    );
                }
                else if (distance == minDistance) {
                    closestCandidates.push_back(
                        candidate
                    );
                }
            }

            std::sort(
                closestCandidates.begin(),
                closestCandidates.end(),
                [&agent](
                    const std::pair<int, int>& leftCell,
                    const std::pair<int, int>& rightCell
                ) {
                    return directionIndex(
                        leftCell.first,
                        leftCell.second,
                        agent.row,
                        agent.col
                    ) < directionIndex(
                        rightCell.first,
                        rightCell.second,
                        agent.row,
                        agent.col
                    );
                }
            );

            chosenRow = closestCandidates[0].first;
            chosenCol = closestCandidates[0].second;
        }
        else {
            chosenRow = candidates[0].first;
            chosenCol = candidates[0].second;
        }

        targetCells[(size_t)agentIndex] = {
            chosenRow,
            chosenCol
        };
    }
}
HPX_PLAIN_ACTION(
    computeTargetsLocalImpl,
    computeTargetsLocalAction
);

// Pack candidates that target neighboring localities
BoundaryCandidatePack packBoundaryCandidatesImpl() {
    BoundaryCandidatePack pack;
    int agentCount = (int)localAgents.size();

    pack.toAbove.reserve(
        (size_t)(agentCount / 16 + 1)
    );
    pack.toBelow.reserve(
        (size_t)(agentCount / 16 + 1)
    );

    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        const Agent& agent = localAgents[(size_t)agentIndex];
        int targetRow = targetCells[(size_t)agentIndex].first;
        int targetCol = targetCells[(size_t)agentIndex].second;

        if (targetRow >= localStartRow && targetRow < localEndRow) {
            continue;
        }

        MoveCandidate candidate;
        candidate.agentId = agent.id;
        candidate.agentRow = agent.row;
        candidate.agentCol = agent.col;
        candidate.targetRow = targetRow;
        candidate.targetCol = targetCol;

        if (targetRow < localStartRow) {
            pack.toAbove.push_back(
                candidate
            );
        }
        else {
            pack.toBelow.push_back(
                candidate
            );
        }
    }

    return pack;
}
HPX_PLAIN_ACTION(
    packBoundaryCandidatesImpl,
    packBoundaryCandidatesAction
);

// Resolve conflicts for owned target cells
BoundaryWinPack resolveConflictsLocalImpl(
    std::vector<MoveCandidate> fromAbove,
    std::vector<MoveCandidate> fromBelow
) {
    int agentCount = (int)localAgents.size();

    currentStamp++;

    if (currentStamp == std::numeric_limits<int>::max()) {
        std::fill(
            cellStamp.begin(),
            cellStamp.end(),
            0
        );
        currentStamp = 1;
    }

    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        moveWon[(size_t)agentIndex] = 0;
        int targetRow = targetCells[(size_t)agentIndex].first;
        int targetCol = targetCells[(size_t)agentIndex].second;

        if (targetRow < localStartRow || targetRow >= localEndRow) {
            continue;
        }

        int localCellIndex = (targetRow - localStartRow) * gridWidth +
            targetCol;
        touchConflictCell(
            localCellIndex
        );

        if (
            localAgents[(size_t)agentIndex].id <
            bestAgentId[(size_t)localCellIndex]
        ) {
            bestLocalAgentIndex[(size_t)localCellIndex] = agentIndex;
            bestAgentId[(size_t)localCellIndex] =
                localAgents[(size_t)agentIndex].id;
        }
    }

    for (const auto& candidate : fromAbove) {
        if (
            candidate.targetRow < localStartRow ||
            candidate.targetRow >= localEndRow
        ) {
            continue;
        }

        int localCellIndex = (candidate.targetRow - localStartRow) *
            gridWidth + candidate.targetCol;
        touchConflictCell(
            localCellIndex
        );

        if (candidate.agentId < bestAgentId[(size_t)localCellIndex]) {
            bestLocalAgentIndex[(size_t)localCellIndex] = -1;
            bestAgentId[(size_t)localCellIndex] = candidate.agentId;
        }
    }

    for (const auto& candidate : fromBelow) {
        if (
            candidate.targetRow < localStartRow ||
            candidate.targetRow >= localEndRow
        ) {
            continue;
        }

        int localCellIndex = (candidate.targetRow - localStartRow) *
            gridWidth + candidate.targetCol;
        touchConflictCell(
            localCellIndex
        );

        if (candidate.agentId < bestAgentId[(size_t)localCellIndex]) {
            bestLocalAgentIndex[(size_t)localCellIndex] = -1;
            bestAgentId[(size_t)localCellIndex] = candidate.agentId;
        }
    }

    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        int targetRow = targetCells[(size_t)agentIndex].first;
        int targetCol = targetCells[(size_t)agentIndex].second;

        if (targetRow < localStartRow || targetRow >= localEndRow) {
            continue;
        }

        int localCellIndex = (targetRow - localStartRow) * gridWidth +
            targetCol;

        if (
            cellStamp[(size_t)localCellIndex] == currentStamp &&
            bestLocalAgentIndex[(size_t)localCellIndex] == agentIndex
        ) {
            moveWon[(size_t)agentIndex] = 1;
        }
    }

    BoundaryWinPack replies;
    replies.winsForAbove.reserve(
        fromAbove.size()
    );
    replies.winsForBelow.reserve(
        fromBelow.size()
    );

    for (const auto& candidate : fromAbove) {
        if (
            candidate.targetRow >= localStartRow &&
            candidate.targetRow < localEndRow
        ) {
            int localCellIndex = (candidate.targetRow - localStartRow) *
                gridWidth + candidate.targetCol;

            if (
                cellStamp[(size_t)localCellIndex] == currentStamp &&
                bestAgentId[(size_t)localCellIndex] == candidate.agentId
            ) {
                replies.winsForAbove.push_back(
                    candidate.agentId
                );
            }
            else {
                replies.winsForAbove.push_back(
                    -1
                );
            }
        }
        else {
            replies.winsForAbove.push_back(
                -1
            );
        }
    }

    for (const auto& candidate : fromBelow) {
        if (
            candidate.targetRow >= localStartRow &&
            candidate.targetRow < localEndRow
        ) {
            int localCellIndex = (candidate.targetRow - localStartRow) *
                gridWidth + candidate.targetCol;

            if (
                cellStamp[(size_t)localCellIndex] == currentStamp &&
                bestAgentId[(size_t)localCellIndex] == candidate.agentId
            ) {
                replies.winsForBelow.push_back(
                    candidate.agentId
                );
            }
            else {
                replies.winsForBelow.push_back(
                    -1
                );
            }
        }
        else {
            replies.winsForBelow.push_back(
                -1
            );
        }
    }

    return replies;
}
HPX_PLAIN_ACTION(
    resolveConflictsLocalImpl,
    resolveConflictsLocalAction
);

// Apply boundary conflict replies to source agents
void applyBoundaryWinsImpl(
    std::vector<int> returnedWinsFromAbove,
    std::vector<int> returnedWinsFromBelow,
    std::vector<MoveCandidate> sentAbove,
    std::vector<MoveCandidate> sentBelow
) {
    for (size_t i = 0; i < sentAbove.size(); i++) {
        if (
            i < returnedWinsFromAbove.size() &&
            returnedWinsFromAbove[i] == sentAbove[i].agentId
        ) {
            int agentIndex = findAgentIndex(
                sentAbove[i].agentId
            );

            if (agentIndex >= 0) {
                moveWon[(size_t)agentIndex] = 1;
            }
        }
    }

    for (size_t i = 0; i < sentBelow.size(); i++) {
        if (
            i < returnedWinsFromBelow.size() &&
            returnedWinsFromBelow[i] == sentBelow[i].agentId
        ) {
            int agentIndex = findAgentIndex(
                sentBelow[i].agentId
            );

            if (agentIndex >= 0) {
                moveWon[(size_t)agentIndex] = 1;
            }
        }
    }
}
HPX_PLAIN_ACTION(
    applyBoundaryWinsImpl,
    applyBoundaryWinsAction
);

// Move winners and harvest owned target sugar
void applyMovesAndTakeSugarImpl() {
    for (size_t agentIndex = 0; agentIndex < localAgents.size(); agentIndex++) {
        if (!moveWon[agentIndex]) {
            continue;
        }

        Agent& agent = localAgents[agentIndex];
        int targetRow = targetCells[agentIndex].first;
        int targetCol = targetCells[agentIndex].second;

        if (agent.row >= localStartRow && agent.row < localEndRow) {
            int oldLocalCell = (agent.row - localStartRow) * gridWidth +
                agent.col;

            if (localGrid[(size_t)oldLocalCell].agentID == agent.id) {
                localGrid[(size_t)oldLocalCell].agentID = -1;
            }
        }

        agent.row = targetRow;
        agent.col = targetCol;

        if (targetRow >= localStartRow && targetRow < localEndRow) {
            int newLocalCell = (targetRow - localStartRow) * gridWidth +
                targetCol;
            int harvestedSugar = localGrid[(size_t)newLocalCell].currentSugar;
            agent.wealth += harvestedSugar;

            if (harvestedSugar > 0) {
                localGrid[(size_t)newLocalCell].currentSugar = 0;
                markSugarGrowing(
                    newLocalCell
                );
            }

            localGrid[(size_t)newLocalCell].agentID = agent.id;
        }
    }
}
HPX_PLAIN_ACTION(
    applyMovesAndTakeSugarImpl,
    applyMovesAndTakeSugarAction
);

// Pack agents whose new row belongs to a neighbor
BoundaryAgentPack packMigrationsImpl() {
    BoundaryAgentPack pack;

    for (const Agent& agent : localAgents) {
        int owner = partitionForRow(
            agent.row,
            gridHeight,
            localityCount
        );

        if (owner == localityId) {
            continue;
        }

        if (owner < localityId) {
            pack.toUp.push_back(
                agent
            );
        }
        else {
            pack.toDown.push_back(
                agent
            );
        }
    }

    return pack;
}
HPX_PLAIN_ACTION(
    packMigrationsImpl,
    packMigrationsAction
);

// Integrate received migrants and rebuild occupancy.
//
// PM2/MASS parity: when a foreign migrant Z wins a boundary cell C currently
// occupied by a local loser-stayer S (S targeted C but lost the conflict to
// Z), Z is the cell's rightful owner. MASS::integrateMigrations harvests on
// arrival and immediately overwrites occupancy with the migrant id BEFORE
// the second-pass take_sugar runs. The prior HPX implementation cleared and
// rebuilt the entire occupancy map by iterating sort-by-id agents, so when
// Z.id < S.id the rebuild overwrote grid[C] with S.id, take_sugar then gave
// the cell's sugar to S instead of Z. Order here mirrors MASS:
//   1. Filter out agents that left this strip.
//   2. Rebuild occupancy from remaining local agents.
//   3. Absorb each migrant: harvest sugar at arrival cell, claim occupancy.
//   4. Sort by id
//      This keeps binary search valid for next apply boundary wins step
void integrateMigrationsImpl(
    std::vector<Agent> fromUp,
    std::vector<Agent> fromDown
) {
    std::vector<Agent> keptAgents;
    keptAgents.reserve(
        localAgents.size() + fromUp.size() + fromDown.size()
    );

    for (const auto& agent : localAgents) {
        if (
            partitionForRow(
                agent.row,
                gridHeight,
                localityCount
            ) == localityId
        ) {
            keptAgents.push_back(
                agent
            );
        }
    }

    localAgents = std::move(
        keptAgents
    );

    for (auto& cell : localGrid) {
        cell.agentID = -1;
    }

    for (const auto& agent : localAgents) {
        if (agent.row >= localStartRow && agent.row < localEndRow) {
            int localCell = (agent.row - localStartRow) * gridWidth +
                agent.col;
            localGrid[(size_t)localCell].agentID = agent.id;
        }
    }

    auto absorbMigrant = [&](Agent& agent) {
        if (agent.row < localStartRow || agent.row >= localEndRow) {
            return;
        }
        int localCell = (agent.row - localStartRow) * gridWidth + agent.col;
        int harvestedSugar = localGrid[(size_t)localCell].currentSugar;
        agent.wealth += harvestedSugar;
        if (harvestedSugar > 0) {
            localGrid[(size_t)localCell].currentSugar = 0;
            markSugarGrowing(
                localCell
            );
        }
        localGrid[(size_t)localCell].agentID = agent.id;
        localAgents.push_back(
            agent
        );
    };

    for (auto& agent : fromUp) {
        absorbMigrant(
            agent
        );
    }
    for (auto& agent : fromDown) {
        absorbMigrant(
            agent
        );
    }

    sortAgentsById();
}
HPX_PLAIN_ACTION(
    integrateMigrationsImpl,
    integrateMigrationsAction
);

// Harvest sugar under agents after migration
void takeSugarAtCurrentLocalImpl() {
    for (auto& agent : localAgents) {
        int localCell = (agent.row - localStartRow) * gridWidth + agent.col;

        if (localGrid[(size_t)localCell].agentID != agent.id) {
            continue;
        }

        int harvestedSugar = localGrid[(size_t)localCell].currentSugar;

        if (harvestedSugar > 0) {
            agent.wealth += harvestedSugar;
            localGrid[(size_t)localCell].currentSugar = 0;
            markSugarGrowing(
                localCell
            );
        }
    }
}
HPX_PLAIN_ACTION(
    takeSugarAtCurrentLocalImpl,
    takeSugarAtCurrentLocalAction
);

// Apply metabolism and remove starved agents
void metabolismAndDeathLocalImpl() {
    for (auto& agent : localAgents) {
        agent.wealth -= agent.metabolism;
    }

    std::vector<Agent> aliveAgents;
    aliveAgents.reserve(
        localAgents.size()
    );

    for (const auto& agent : localAgents) {
        if (agent.wealth <= 0) {
            if (agent.row >= localStartRow && agent.row < localEndRow) {
                int localCell = (agent.row - localStartRow) * gridWidth +
                    agent.col;
                localGrid[(size_t)localCell].agentID = -1;
            }
        }
        else {
            aliveAgents.push_back(
                agent
            );
        }
    }

    localAgents = std::move(
        aliveAgents
    );
}
HPX_PLAIN_ACTION(
    metabolismAndDeathLocalImpl,
    metabolismAndDeathLocalAction
);

// Return local agents for final summary and optional CSV
std::vector<Agent> finalAgentsSnapshotImpl() {
    return localAgents;
}
HPX_PLAIN_ACTION(
    finalAgentsSnapshotImpl,
    finalAgentsSnapshotAction
);

// Return local sugar total for final summary
int finalLocalSugarImpl() {
    int totalSugar = 0;

    for (const auto& cell : localGrid) {
        totalSugar += cell.currentSugar;
    }

    return totalSugar;
}
HPX_PLAIN_ACTION(
    finalLocalSugarImpl,
    finalLocalSugarAction
);
// Initialize global grid and seed agents on owner localities
static void initializeLocalities(
    int height,
    int width,
    int agentCount,
    const std::vector<hpx::id_type>& localities,
    int localityTotal
) {
    const int sugarCapMin = 1;
    const int sugarCapMax = 4;
    const int metabolismMin = 1;
    const int metabolismMax = 4;
    const int visionMin = 1;
    const int visionMax = HPX_VISION_MAX;
    const unsigned int initialWealth = (unsigned int)-1;

    Xoshiro256 sugarRandom(
        1
    );
    Xoshiro256 metabolismRandom(
        2
    );
    Xoshiro256 visionRandom(
        3
    );
    Xoshiro256 coordRandom(
        4
    );
    Xoshiro256 wealthRandom(
        100
    );

    int totalCells = height * width;
    std::vector<GridCell> grid(
        (size_t)totalCells
    );

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            sugarRandom.next();
            int cellIndex = row * width + col;
            int capacity = sugarCapacity(
                row,
                col,
                height,
                width,
                sugarCapMin,
                sugarCapMax
            );

            grid[(size_t)cellIndex].currentSugar = capacity;
            grid[(size_t)cellIndex].sugarCapacity = capacity;
            grid[(size_t)cellIndex].agentID = -1;
        }
    }

    std::vector<Agent> agents(
        (size_t)agentCount
    );

    for (int i = 0; i < agentCount; i++) {
        Agent& agent = agents[(size_t)i];
        agent.id = i;
        agent.metabolism = (int)metabolismRandom.nextInRange(
            metabolismMin,
            metabolismMax
        );
        agent.vision = (int)visionRandom.nextInRange(
            visionMin,
            visionMax
        );

        if (initialWealth == (unsigned int)-1) {
            agent.wealth = (int)wealthRandom.nextInRange(
                (uint64_t)(agent.metabolism * 5),
                (uint64_t)(agent.metabolism * 10)
            );
        }
        else {
            agent.wealth = (int)initialWealth;
        }

        agent.age = 0;
    }

    std::vector<int> shuffledRows(
        (size_t)totalCells
    );
    std::vector<int> shuffledCols(
        (size_t)totalCells
    );

    int writeIndex = 0;

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            shuffledRows[(size_t)writeIndex] = row;
            shuffledCols[(size_t)writeIndex] = col;
            writeIndex++;
        }
    }

    for (int i = totalCells - 1; i > 0; i--) {
        int j = (int)coordRandom.nextInRange(
            0,
            (uint64_t)i
        );

        if (j != i) {
            std::swap(
                shuffledRows[(size_t)i],
                shuffledRows[(size_t)j]
            );
            std::swap(
                shuffledCols[(size_t)i],
                shuffledCols[(size_t)j]
            );
        }
    }

    for (int i = 0; i < agentCount; i++) {
        Agent& agent = agents[(size_t)i];
        agent.row = shuffledRows[(size_t)i];
        agent.col = shuffledCols[(size_t)i];
        grid[(size_t)(agent.row * width + agent.col)].agentID = agent.id;
    }

    int rowsPerLocality = (height + localityTotal - 1) / localityTotal;

    auto startRowFor = [rowsPerLocality](
        int localityIndex
    ) {
        return localityIndex * rowsPerLocality;
    };

    auto endRowFor = [
        height,
        rowsPerLocality
    ](
        int localityIndex
    ) {
        return std::min(
            localityIndex * rowsPerLocality + rowsPerLocality,
            height
        );
    };

    std::vector<hpx::future<void>> futures;
    futures.reserve(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        int rowStart = startRowFor(
            localityIndex
        );
        int rowEnd = endRowFor(
            localityIndex
        );

        std::vector<GridCell> chunk(
            (size_t)(rowEnd - rowStart) * (size_t)width
        );

        for (int row = rowStart; row < rowEnd; row++) {
            for (int col = 0; col < width; col++) {
                chunk[(size_t)(row - rowStart) * width + col] =
                    grid[(size_t)row * width + col];
            }
        }

        std::vector<Agent> chunkAgents;
        chunkAgents.reserve(
            agents.size() / (size_t)std::max(
                localityTotal,
                1
            ) + 64
        );

        for (const auto& agent : agents) {
            if (
                partitionForRow(
                    agent.row,
                    height,
                    localityTotal
                ) == localityIndex
            ) {
                chunkAgents.push_back(
                    agent
                );
            }
        }

        std::sort(
            chunkAgents.begin(),
            chunkAgents.end(),
            [](
                const Agent& leftAgent,
                const Agent& rightAgent
            ) {
                return leftAgent.id < rightAgent.id;
            }
        );

        futures.push_back(
            hpx::async<initWorkerAction>(
                localities[(size_t)localityIndex],
                localityIndex,
                localityTotal,
                rowStart,
                rowEnd,
                std::move(
                    chunk
                ),
                std::move(
                    chunkAgents
                ),
                height,
                width
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

// Exchange sugar halo rows between neighboring localities
static void exchangeHaloRows(
    const std::vector<hpx::id_type>& localities,
    int localityTotal,
    std::vector<std::vector<int>>& topHalo,
    std::vector<std::vector<int>>& bottomHalo
) {
    topHalo.assign(
        (size_t)localityTotal,
        std::vector<int>()
    );
    bottomHalo.assign(
        (size_t)localityTotal,
        std::vector<int>()
    );

    if (localityTotal <= 1) {
        return;
    }

    std::vector<hpx::future<std::vector<int>>> futures;
    futures.reserve(
        (size_t)localityTotal * 2
    );

    std::vector<int> topRequest(
        (size_t)localityTotal,
        -1
    );
    std::vector<int> bottomRequest(
        (size_t)localityTotal,
        -1
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        if (localityIndex > 0) {
            topRequest[(size_t)localityIndex] = (int)futures.size();
            futures.push_back(
                hpx::async<getSugarHaloAction>(
                    localities[(size_t)(localityIndex - 1)],
                    0,
                    HPX_VISION_MAX
                )
            );
        }

        if (localityIndex < localityTotal - 1) {
            bottomRequest[(size_t)localityIndex] = (int)futures.size();
            futures.push_back(
                hpx::async<getSugarHaloAction>(
                    localities[(size_t)(localityIndex + 1)],
                    1,
                    HPX_VISION_MAX
                )
            );
        }
    }

    hpx::wait_all(
        futures
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        if (topRequest[(size_t)localityIndex] >= 0) {
            topHalo[(size_t)localityIndex] = futures[
                (size_t)topRequest[(size_t)localityIndex]
            ].get();
        }

        if (bottomRequest[(size_t)localityIndex] >= 0) {
            bottomHalo[(size_t)localityIndex] = futures[
                (size_t)bottomRequest[(size_t)localityIndex]
            ].get();
        }
    }
}

// Run target computation on all localities
static void computeTargetsForAllLocalities(
    const std::vector<hpx::id_type>& localities,
    int localityTotal,
    std::vector<std::vector<int>>& topHalo,
    std::vector<std::vector<int>>& bottomHalo
) {
    std::vector<hpx::future<void>> futures(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        futures[(size_t)localityIndex] = hpx::async<computeTargetsLocalAction>(
            localities[(size_t)localityIndex],
            std::move(
                topHalo[(size_t)localityIndex]
            ),
            std::move(
                bottomHalo[(size_t)localityIndex]
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

// Resolve boundary moves and apply returned win flags
static void resolveBoundaryMoves(
    const std::vector<hpx::id_type>& localities,
    int localityTotal
) {
    std::vector<BoundaryCandidatePack> candidatePacks(
        (size_t)localityTotal
    );

    {
        std::vector<hpx::future<BoundaryCandidatePack>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<packBoundaryCandidatesAction>(
                    localities[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            candidatePacks[(size_t)localityIndex] =
                futures[(size_t)localityIndex].get();
        }
    }

    std::vector<std::vector<MoveCandidate>> fromAbove(
        (size_t)localityTotal
    );
    std::vector<std::vector<MoveCandidate>> fromBelow(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        if (localityIndex > 0) {
            fromAbove[(size_t)localityIndex] = candidatePacks[
                (size_t)(localityIndex - 1)
            ].toBelow;
        }

        if (localityIndex < localityTotal - 1) {
            fromBelow[(size_t)localityIndex] = candidatePacks[
                (size_t)(localityIndex + 1)
            ].toAbove;
        }
    }

    std::vector<BoundaryWinPack> winPacks(
        (size_t)localityTotal
    );

    {
        std::vector<hpx::future<BoundaryWinPack>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<resolveConflictsLocalAction>(
                    localities[(size_t)localityIndex],
                    fromAbove[(size_t)localityIndex],
                    fromBelow[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            winPacks[(size_t)localityIndex] =
                futures[(size_t)localityIndex].get();
        }
    }

    std::vector<std::vector<int>> returnedWinsAbove(
        (size_t)localityTotal
    );
    std::vector<std::vector<int>> returnedWinsBelow(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        if (localityIndex > 0) {
            returnedWinsAbove[(size_t)localityIndex] = winPacks[
                (size_t)(localityIndex - 1)
            ].winsForBelow;
        }

        if (localityIndex < localityTotal - 1) {
            returnedWinsBelow[(size_t)localityIndex] = winPacks[
                (size_t)(localityIndex + 1)
            ].winsForAbove;
        }
    }

    std::vector<hpx::future<void>> futures(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        futures[(size_t)localityIndex] = hpx::async<applyBoundaryWinsAction>(
            localities[(size_t)localityIndex],
            std::move(
                returnedWinsAbove[(size_t)localityIndex]
            ),
            std::move(
                returnedWinsBelow[(size_t)localityIndex]
            ),
            std::move(
                candidatePacks[(size_t)localityIndex].toAbove
            ),
            std::move(
                candidatePacks[(size_t)localityIndex].toBelow
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

// Apply local winning moves
static void applyMovesOnAllLocalities(
    const std::vector<hpx::id_type>& localities,
    int localityTotal
) {
    std::vector<hpx::future<void>> futures(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        futures[(size_t)localityIndex] =
            hpx::async<applyMovesAndTakeSugarAction>(
                localities[(size_t)localityIndex]
            );
    }

    hpx::wait_all(
        futures
    );
}

// Move agents to owning localities after crossing boundaries
static void migrateAgentsBetweenLocalities(
    const std::vector<hpx::id_type>& localities,
    int localityTotal
) {
    std::vector<BoundaryAgentPack> migrationPacks(
        (size_t)localityTotal
    );

    {
        std::vector<hpx::future<BoundaryAgentPack>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<packMigrationsAction>(
                    localities[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            migrationPacks[(size_t)localityIndex] =
                futures[(size_t)localityIndex].get();
        }
    }

    std::vector<std::vector<Agent>> fromUp(
        (size_t)localityTotal
    );
    std::vector<std::vector<Agent>> fromDown(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        if (localityIndex > 0) {
            fromUp[(size_t)localityIndex] = migrationPacks[
                (size_t)(localityIndex - 1)
            ].toDown;
        }

        if (localityIndex < localityTotal - 1) {
            fromDown[(size_t)localityIndex] = migrationPacks[
                (size_t)(localityIndex + 1)
            ].toUp;
        }
    }

    std::vector<hpx::future<void>> futures(
        (size_t)localityTotal
    );

    for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
        futures[(size_t)localityIndex] = hpx::async<integrateMigrationsAction>(
            localities[(size_t)localityIndex],
            std::move(
                fromUp[(size_t)localityIndex]
            ),
            std::move(
                fromDown[(size_t)localityIndex]
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

// Harvest after migration and apply metabolism
static void finishTimestepOnAllLocalities(
    const std::vector<hpx::id_type>& localities,
    int localityTotal
) {
    {
        std::vector<hpx::future<void>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<takeSugarAtCurrentLocalAction>(
                    localities[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );
    }

    {
        std::vector<hpx::future<void>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<metabolismAndDeathLocalAction>(
                    localities[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );
    }
}

// Write optional final agent CSV
static void writeFinalAgentCsv(
    const std::string& outputPrefix,
    int timesteps,
    std::vector<Agent>& finalAgents
) {
    if (outputPrefix.empty()) {
        return;
    }

    std::sort(
        finalAgents.begin(),
        finalAgents.end(),
        [](
            const Agent& leftAgent,
            const Agent& rightAgent
        ) {
            return leftAgent.id < rightAgent.id;
        }
    );

    std::ostringstream timestepText;
    timestepText
        << std::setfill(
            '0'
        )
        << std::setw(
            6
        )
        << timesteps;

    std::ofstream agentFile(
        outputPrefix + "_" + timestepText.str() + "_AGENT.csv"
    );

    agentFile << "agentID,row,col,wealth,metabolism,vision,age\n";

    for (const auto& agent : finalAgents) {
        agentFile
            << agent.id << ","
            << agent.row << ","
            << agent.col << ","
            << agent.wealth << ","
            << agent.metabolism << ","
            << agent.vision << ","
            << agent.age << "\n";
    }

    agentFile.close();
}

// Collect final stats after timed loop
static void collectAndPrintFinalStats(
    const std::vector<hpx::id_type>& localities,
    int localityTotal,
    int timesteps,
    const std::string& outputPrefix,
    std::chrono::steady_clock::time_point totalStartTime,
    std::chrono::steady_clock::time_point simulationStartTime,
    std::chrono::steady_clock::time_point totalEndTime
) {
    std::vector<Agent> finalAgents;

    {
        std::vector<hpx::future<std::vector<Agent>>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<finalAgentsSnapshotAction>(
                    localities[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );

        std::vector<std::vector<Agent>> perLocality(
            (size_t)localityTotal
        );
        size_t totalAgentCount = 0;

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            perLocality[(size_t)localityIndex] =
                futures[(size_t)localityIndex].get();
            totalAgentCount += perLocality[(size_t)localityIndex].size();
        }

        finalAgents.reserve(
            totalAgentCount
        );

        for (auto& agentChunk : perLocality) {
            finalAgents.insert(
                finalAgents.end(),
                std::make_move_iterator(
                    agentChunk.begin()
                ),
                std::make_move_iterator(
                    agentChunk.end()
                )
            );
        }
    }

    long long totalWealth = 0;

    for (const auto& agent : finalAgents) {
        totalWealth += agent.wealth;
    }

    int totalSugar = 0;

    {
        std::vector<hpx::future<int>> futures(
            (size_t)localityTotal
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            futures[(size_t)localityIndex] =
                hpx::async<finalLocalSugarAction>(
                    localities[(size_t)localityIndex]
                );
        }

        hpx::wait_all(
            futures
        );

        for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
            totalSugar += futures[(size_t)localityIndex].get();
        }
    }

    std::cout
        << "BENCH_STAT agents=" << finalAgents.size()
        << " wealth=" << totalWealth
        << " sugar=" << totalSugar
        << std::endl;

    long initializationMs = std::chrono::duration_cast<
        std::chrono::milliseconds
    >(
        simulationStartTime - totalStartTime
    ).count();

    long simulationMs = std::chrono::duration_cast<
        std::chrono::milliseconds
    >(
        totalEndTime - simulationStartTime
    ).count();

    long totalMs = std::chrono::duration_cast<
        std::chrono::milliseconds
    >(
        totalEndTime - totalStartTime
    ).count();

    std::cout << "INIT_MS=" << initializationMs << std::endl;
    std::cout << "WALL_TIME_MS=" << simulationMs << std::endl;
    std::cout << "TOTAL_MS=" << totalMs << std::endl;

    writeFinalAgentCsv(
        outputPrefix,
        timesteps,
        finalAgents
    );
}

// Run distributed SugarScape simulation
static void runSimulation(
    int height,
    int width,
    int agentCount,
    int timesteps,
    const std::string& outputPrefix,
    const std::vector<hpx::id_type>& localities,
    int localityTotal
) {
    auto totalStartTime = std::chrono::steady_clock::now();

    initializeLocalities(
        height,
        width,
        agentCount,
        localities,
        localityTotal
    );

    auto simulationStartTime = std::chrono::steady_clock::now();

    for (int timestep = 0; timestep < timesteps; timestep++) {
        if (timestep > 0) {
            std::vector<hpx::future<void>> futures(
                (size_t)localityTotal
            );

            for (
        int localityIndex = 0;
        localityIndex < localityTotal;
        localityIndex++
    ) {
                futures[(size_t)localityIndex] =
                    hpx::async<regrowLocalAction>(
                        localities[(size_t)localityIndex],
                        1
                    );
            }

            hpx::wait_all(
                futures
            );
        }

        std::vector<std::vector<int>> topHalo;
        std::vector<std::vector<int>> bottomHalo;

        exchangeHaloRows(
            localities,
            localityTotal,
            topHalo,
            bottomHalo
        );

        computeTargetsForAllLocalities(
            localities,
            localityTotal,
            topHalo,
            bottomHalo
        );

        resolveBoundaryMoves(
            localities,
            localityTotal
        );

        applyMovesOnAllLocalities(
            localities,
            localityTotal
        );

        migrateAgentsBetweenLocalities(
            localities,
            localityTotal
        );

        finishTimestepOnAllLocalities(
            localities,
            localityTotal
        );
    }

    auto totalEndTime = std::chrono::steady_clock::now();

    collectAndPrintFinalStats(
        localities,
        localityTotal,
        timesteps,
        outputPrefix,
        totalStartTime,
        simulationStartTime,
        totalEndTime
    );
}

// HPX entry point
int hpx_main(
    po::variables_map& variables
) {
    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    int localityTotal = (int)localities.size();

    int height = variables["height"].as<int>();
    int width = variables["width"].as<int>();
    int agents = variables["agents"].as<int>();
    int timesteps = variables["timesteps"].as<int>();
    std::string outputPrefix = variables.count(
        "output"
    ) ? variables["output"].as<std::string>() : "";

    if (
        height <= 0 ||
        width <= 0 ||
        agents > height * width
    ) {
        std::cerr << "Invalid height, width, or agentCount\n";
        return hpx::finalize();
    }

    if (
        localityTotal > 1 &&
        height / localityTotal < 2 * HPX_VISION_MAX
    ) {
        std::cerr
            << "Error: partition height " << (height / localityTotal)
            << " < 2*visionMax=" << (2 * HPX_VISION_MAX)
            << ". Reduce nlocs or increase gridHeight.\n";
        return hpx::finalize();
    }

    runSimulation(
        height,
        width,
        agents,
        timesteps,
        outputPrefix,
        localities,
        localityTotal
    );

    return hpx::finalize();
}

// Create HPX command line options and start runtime
int runSugarScapeBenchmark(
    int argc,
    char* argv[]
) {
    po::options_description description(
        "HPX SugarScape"
    );

    description.add_options()
        (
            "height",
            po::value<int>()->default_value(
                32
            ),
            "Grid height"
        )
        (
            "width",
            po::value<int>()->default_value(
                32
            ),
            "Grid width"
        )
        (
            "agents",
            po::value<int>()->default_value(
                50
            ),
            "Initial agent count"
        )
        (
            "timesteps",
            po::value<int>()->default_value(
                100
            ),
            "Simulation steps"
        )
        (
            "output,o",
            po::value<std::string>(),
            "Output prefix for CSV"
        );

    hpx::init_params parameters;
    parameters.desc_cmdline = description;

    return hpx::init(
        argc,
        argv,
        parameters
    );
}
