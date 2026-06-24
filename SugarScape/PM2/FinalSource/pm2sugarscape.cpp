#include "pm2sugarscape.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// SugarScape movement only benchmark
// Reproduction and age-based death disabled

MPI_Datatype createAgentType() {
    MPI_Datatype agentType;
    int blockLengths[7] = { 1, 1, 1, 1, 1, 1, 1 };
    MPI_Aint offsets[7];
    MPI_Datatype fieldTypes[7];
    Agent dummyAgent;
    MPI_Aint baseAddress;

    // Build MPI layout from struct field addresses
    MPI_Get_address(
        &dummyAgent,
        &baseAddress
    );
    MPI_Get_address(
        &dummyAgent.id,
        &offsets[0]
    );
    MPI_Get_address(
        &dummyAgent.row,
        &offsets[1]
    );
    MPI_Get_address(
        &dummyAgent.col,
        &offsets[2]
    );
    MPI_Get_address(
        &dummyAgent.wealth,
        &offsets[3]
    );
    MPI_Get_address(
        &dummyAgent.metabolism,
        &offsets[4]
    );
    MPI_Get_address(
        &dummyAgent.vision,
        &offsets[5]
    );
    MPI_Get_address(
        &dummyAgent.age,
        &offsets[6]
    );

    // Convert absolute addresses into relative offsets
    for (int i = 0; i < 7; i++) {
        offsets[i] -= baseAddress;
        fieldTypes[i] = MPI_INT;
    }

    // Commit packed Agent layout used by sends and gathers
    MPI_Type_create_struct(
        7,
        blockLengths,
        offsets,
        fieldTypes,
        &agentType
    );
    MPI_Type_commit(
        &agentType
    );

    return agentType;
}

// Build MPI datatype for move candidates
MPI_Datatype createMoveCandidateType() {
    MPI_Datatype moveCandidateType;
    int blockLengths[5] = { 1, 1, 1, 1, 1 };
    MPI_Aint offsets[5];
    MPI_Datatype fieldTypes[5];
    MoveCandidate dummyCandidate;
    MPI_Aint baseAddress;

    MPI_Get_address(
        &dummyCandidate,
        &baseAddress
    );
    // Match MoveCandidate field order exactly
    MPI_Get_address(
        &dummyCandidate.agentId,
        &offsets[0]
    );
    MPI_Get_address(
        &dummyCandidate.agentRow,
        &offsets[1]
    );
    MPI_Get_address(
        &dummyCandidate.agentCol,
        &offsets[2]
    );
    MPI_Get_address(
        &dummyCandidate.targetRow,
        &offsets[3]
    );
    MPI_Get_address(
        &dummyCandidate.targetCol,
        &offsets[4]
    );

    for (int i = 0; i < 5; i++) {
        offsets[i] -= baseAddress;
        fieldTypes[i] = MPI_INT;
    }

    // Commit candidate layout
    MPI_Type_create_struct(
        5,
        blockLengths,
        offsets,
        fieldTypes,
        &moveCandidateType
    );
    MPI_Type_commit(
        &moveCandidateType
    );

    return moveCandidateType;
}

// Expand one seed into full RNG state
uint64_t Xoshiro256::splitmix64(
    uint64_t& seedValue
) {
    seedValue += 0x9e3779b97f4a7c15ULL;

    uint64_t mixedValue = seedValue;
    mixedValue = (mixedValue ^ (mixedValue >> 30)) * 0xbf58476d1ce4e5b9ULL;
    mixedValue = (mixedValue ^ (mixedValue >> 27)) * 0x94d049bb133111ebULL;

    return mixedValue ^ (mixedValue >> 31);
}

// Rotate bits for xoshiro state update
uint64_t Xoshiro256::rotateLeft(
    uint64_t value,
    int shift
) {
    return (value << shift) | (value >> (64 - shift));
}

// Initialize deterministic RNG state
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

// Generate unbiased value in an inclusive range
uint64_t Xoshiro256::nextInRange(
    uint64_t minValue,
    uint64_t maxValue
) {
    if (maxValue < minValue) {
        return next();
    }

    // Rejection bound to avoid bias from modulo
    uint64_t range = maxValue - minValue + 1;
    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    uint64_t randomValue;
    do {
        randomValue = next();
    } while (randomValue >= limit);

    return minValue + (randomValue % range);
}

// Map global row to slab owner
int rankForRow(
    int globalRow,
    int gridHeight,
    int processCount
) {
    if (processCount <= 1) {
        return 0;
    }

    // Use ceiling division 
    int rowsPerRank = (gridHeight + processCount - 1) / processCount;
    int ownerRank = globalRow / rowsPerRank;

    if (ownerRank >= processCount) {
        return processCount - 1;
    }

    return ownerRank;
}

// Get first owned global row
int firstRow(
    int rank,
    int gridHeight,
    int processCount
) {
    // Use ceiling division 
    int rowsPerRank = (gridHeight + processCount - 1) / processCount;
    return rank * rowsPerRank;
}

// Get one-past-last owned global row
int endRow(
    int rank,
    int gridHeight,
    int processCount
) {
    // Use ceiling division 
    int rowsPerRank = (gridHeight + processCount - 1) / processCount;
    return min(
        rank * rowsPerRank + rowsPerRank,
        gridHeight
    );
}


// Find sugar capacity in two peak landscape
int sugarCapacity(
    int row,
    int col,
    int gridHeight,
    int gridWidth,
    int minCapacity,
    int maxCapacity
) {
    // Fixed sugar peaks
    double firstPeakRow = gridHeight * 0.25;
    double firstPeakCol = gridWidth * 0.25;
    double secondPeakRow = gridHeight * 0.75;
    double secondPeakCol = gridWidth * 0.75;
    double maxRadius = min(
        (double)gridHeight,
        (double)gridWidth
    ) * 0.4;

    double firstDistance = sqrt(
        (row - firstPeakRow) * (row - firstPeakRow) +
        (col - firstPeakCol) * (col - firstPeakCol)
    );

    double secondDistance = sqrt(
        (row - secondPeakRow) * (row - secondPeakRow) +
        (col - secondPeakCol) * (col - secondPeakCol)
    );

    // Capacity depends on closest peak
    double nearestDistance = min(
        firstDistance,
        secondDistance
    );

    if (nearestDistance <= maxRadius * 0.25) {
        return maxCapacity;
    }

    if (nearestDistance <= maxRadius * 0.50) {
        return min(
            maxCapacity,
            max(
                minCapacity,
                maxCapacity - 1
            )
        );
    }

    if (nearestDistance <= maxRadius * 0.75) {
        return min(
            maxCapacity,
            max(
                minCapacity,
                maxCapacity - 2
            )
        );
    }

    return minCapacity;
}


// Exchange halo rows with neighbors
int postHaloExchange(
    const vector<int>& localSugar,
    vector<int>& haloRowsAbove,
    vector<int>& haloRowsBelow,
    int localStartRow,
    int localRowCount,
    int gridWidth,
    int maxVision,
    int rank,
    int processCount,
    int gridHeight,
    vector<MPI_Request>& requests
) {
    // Reuse request vector 
    requests.clear();

    // Only exchange rows that actually exist
    int haloRowCountAbove = min(
        maxVision,
        localStartRow
    );
    int haloRowCountBelow = min(
        maxVision,
        gridHeight - (localStartRow + localRowCount)
    );

    haloRowsAbove.resize(
        haloRowCountAbove * gridWidth
    );
    haloRowsBelow.resize(
        haloRowCountBelow * gridWidth
    );

    int previousRank = rank - 1;
    int nextRank = rank + 1;

    // Send top owned rows upward and receive upper halo
    if (previousRank >= 0 && haloRowCountAbove > 0) {
        int sendRowCount = min(
            maxVision,
            localRowCount
        );
        MPI_Request sendRequest;
        MPI_Request recvRequest;

        MPI_Isend(
            localSugar.data(),
            sendRowCount * gridWidth,
            MPI_INT,
            previousRank,
            0,
            MPI_COMM_WORLD,
            &sendRequest
        );

        MPI_Irecv(
            haloRowsAbove.data(),
            haloRowCountAbove * gridWidth,
            MPI_INT,
            previousRank,
            1,
            MPI_COMM_WORLD,
            &recvRequest
        );

        requests.push_back(
            sendRequest
        );
        requests.push_back(
            recvRequest
        );
    }

    // Send bottom owned rows downward and receive lower halo
    if (nextRank < processCount && haloRowCountBelow > 0) {
        int sendRowCount = min(
            maxVision,
            localRowCount
        );
        int sendOffset = (localRowCount - sendRowCount) * gridWidth;
        MPI_Request sendRequest;
        MPI_Request recvRequest;

        MPI_Isend(
            localSugar.data() + sendOffset,
            sendRowCount * gridWidth,
            MPI_INT,
            nextRank,
            1,
            MPI_COMM_WORLD,
            &sendRequest
        );

        MPI_Irecv(
            haloRowsBelow.data(),
            haloRowCountBelow * gridWidth,
            MPI_INT,
            nextRank,
            0,
            MPI_COMM_WORLD,
            &recvRequest
        );

        requests.push_back(
            sendRequest
        );
        requests.push_back(
            recvRequest
        );
    }

    return (int)requests.size();
}

// Resolve movement conflicts
void resolveConflictsDistributed(
    vector<Agent>& agents,
    vector<pair<int, int>>& targetCells,
    vector<int>& localSugar,
    vector<int>& occupantIds,
    int localStartRow,
    int localRowCount,
    int gridWidth,
    int gridHeight,
    int rank,
    int processCount,
    int maxVision,
    vector<MoveResult>& moveResults
) {
    (void)localSugar;
    (void)occupantIds;
    (void)gridHeight;
    (void)maxVision;

    int localEndRow = localStartRow + localRowCount;
    int agentCount = (int)agents.size();

    // One result per local agent
    moveResults.resize(
        agentCount
    );

    vector<MoveCandidate> candidatesToAbove;
    vector<MoveCandidate> candidatesToBelow;
    vector<int> candidateIndexesToAbove;
    vector<int> candidateIndexesToBelow;

    // Most moves are still local, so reserve a small boundary space
    candidatesToAbove.reserve(
        agentCount / 16 + 1
    );
    candidatesToBelow.reserve(
        agentCount / 16 + 1
    );
    candidateIndexesToAbove.reserve(
        agentCount / 16 + 1
    );
    candidateIndexesToBelow.reserve(
        agentCount / 16 + 1
    );

    // Timestamp arrays avoid clearing every step
    static vector<int> bestLocalAgentIndex;
    static vector<int> bestAgentId;
    static vector<int> cellStamp;
    static int currentStamp = 1;

    size_t localCellCount = (size_t)localRowCount * (size_t)gridWidth;

    // Resize stamps if local slab size changed
    if (bestLocalAgentIndex.size() != localCellCount) {
        bestLocalAgentIndex.assign(
            localCellCount,
            -1
        );
        bestAgentId.assign(
            localCellCount,
            numeric_limits<int>::max()
        );
        cellStamp.assign(
            localCellCount,
            0
        );
        currentStamp = 1;
    }
    else {
        currentStamp++;

        if (currentStamp == numeric_limits<int>::max()) {
            fill(
                cellStamp.begin(),
                cellStamp.end(),
                0
            );
            currentStamp = 1;
        }
    }

    // Reset conflict state for touched cells
    auto touchCell = [&](int localCellIndex) {
        if (cellStamp[localCellIndex] != currentStamp) {
            cellStamp[localCellIndex] = currentStamp;
            bestLocalAgentIndex[localCellIndex] = -1;
            bestAgentId[localCellIndex] = numeric_limits<int>::max();
        }
    };

    // Build local winners and outbound boundary candidates
    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        moveResults[agentIndex].agentId = agents[agentIndex].id;
        moveResults[agentIndex].targetRow = targetCells[agentIndex].first;
        moveResults[agentIndex].targetCol = targetCells[agentIndex].second;
        moveResults[agentIndex].won = 0;

        int targetRow = targetCells[agentIndex].first;
        int targetCol = targetCells[agentIndex].second;

        // Boundary targets are resolved by owner rank
        if (targetRow < localStartRow || targetRow >= localEndRow) {
            MoveCandidate candidate;
            candidate.agentId = agents[agentIndex].id;
            candidate.agentRow = agents[agentIndex].row;
            candidate.agentCol = agents[agentIndex].col;
            candidate.targetRow = targetRow;
            candidate.targetCol = targetCol;

            if (targetRow < localStartRow) {
                candidatesToAbove.push_back(
                    candidate
                );
                candidateIndexesToAbove.push_back(
                    agentIndex
                );
            }
            else {
                candidatesToBelow.push_back(
                    candidate
                );
                candidateIndexesToBelow.push_back(
                    agentIndex
                );
            }
        }
        else {
            int localCellIndex =
                (targetRow - localStartRow) * gridWidth + targetCol;
            touchCell(
                localCellIndex
            );

            // Lowest ID wins (Sorry to high IDs!)
            if (agents[agentIndex].id < bestAgentId[localCellIndex]) {
                bestLocalAgentIndex[localCellIndex] = agentIndex;
                bestAgentId[localCellIndex] = agents[agentIndex].id;
            }
        }
    }

    static MPI_Datatype moveCandidateType = MPI_DATATYPE_NULL;

    if (moveCandidateType == MPI_DATATYPE_NULL) {
        moveCandidateType = createMoveCandidateType();
    }

    int previousRank = rank - 1;
    int nextRank = rank + 1;

    vector<MoveCandidate> candidatesFromAbove;
    vector<MoveCandidate> candidatesFromBelow;

    // Send remote candidates to neighboring owner ranks
    exchangeNeighborVectors(
        candidatesToAbove,
        candidatesToBelow,
        candidatesFromAbove,
        candidatesFromBelow,
        moveCandidateType,
        previousRank,
        nextRank,
        processCount,
        10,
        11,
        12,
        13
    );

    // Upper neighbor candidates -> Owned cell winners
    for (const auto& moveCandidate : candidatesFromAbove) {
        int targetRow = moveCandidate.targetRow;
        int targetCol = moveCandidate.targetCol;

        if (targetRow >= localStartRow && targetRow < localEndRow) {
            int localCellIndex =
                (targetRow - localStartRow) * gridWidth + targetCol;
            touchCell(
                localCellIndex
            );

            if (moveCandidate.agentId < bestAgentId[localCellIndex]) {
                bestLocalAgentIndex[localCellIndex] = -1;
                bestAgentId[localCellIndex] = moveCandidate.agentId;
            }
        }
    }

    // Lower neighbor candidates -> Owned cell winners
    for (const auto& moveCandidate : candidatesFromBelow) {
        int targetRow = moveCandidate.targetRow;
        int targetCol = moveCandidate.targetCol;

        if (targetRow >= localStartRow && targetRow < localEndRow) {
            int localCellIndex =
                (targetRow - localStartRow) * gridWidth + targetCol;
            touchCell(
                localCellIndex
            );

            if (moveCandidate.agentId < bestAgentId[localCellIndex]) {
                bestLocalAgentIndex[localCellIndex] = -1;
                bestAgentId[localCellIndex] = moveCandidate.agentId;
            }
        }
    }

    // Mark local agents that won conflicts
    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        int targetRow = targetCells[agentIndex].first;
        int targetCol = targetCells[agentIndex].second;

        if (targetRow < localStartRow || targetRow >= localEndRow) {
            continue;
        }

        int localCellIndex =
            (targetRow - localStartRow) * gridWidth + targetCol;

        if (
            cellStamp[localCellIndex] == currentStamp &&
            bestLocalAgentIndex[localCellIndex] == agentIndex
        ) {
            moveResults[agentIndex].won = 1;
        }
    }

    // Return boundary winner IDs to requesting ranks
    vector<int> winnerIdsForAbove;
    vector<int> winnerIdsForBelow;

    winnerIdsForAbove.reserve(
        candidatesFromAbove.size()
    );
    winnerIdsForBelow.reserve(
        candidatesFromBelow.size()
    );

    // Upper neighbor candidates -> Owned cell winners
    for (const auto& moveCandidate : candidatesFromAbove) {
        int targetRow = moveCandidate.targetRow;
        int targetCol = moveCandidate.targetCol;

        if (targetRow >= localStartRow && targetRow < localEndRow) {
            int localCellIndex =
                (targetRow - localStartRow) * gridWidth + targetCol;

            if (
                cellStamp[localCellIndex] == currentStamp &&
                bestAgentId[localCellIndex] == moveCandidate.agentId
            ) {
                winnerIdsForAbove.push_back(
                    moveCandidate.agentId
                );
            }
            else {
                winnerIdsForAbove.push_back(
                    -1
                );
            }
        }
        else {
            winnerIdsForAbove.push_back(
                -1
            );
        }
    }

    // Lower neighbor candidates -> Owned cell winners
    for (const auto& moveCandidate : candidatesFromBelow) {
        int targetRow = moveCandidate.targetRow;
        int targetCol = moveCandidate.targetCol;

        if (targetRow >= localStartRow && targetRow < localEndRow) {
            int localCellIndex =
                (targetRow - localStartRow) * gridWidth + targetCol;

            if (
                cellStamp[localCellIndex] == currentStamp &&
                bestAgentId[localCellIndex] == moveCandidate.agentId
            ) {
                winnerIdsForBelow.push_back(
                    moveCandidate.agentId
                );
            }
            else {
                winnerIdsForBelow.push_back(
                    -1
                );
            }
        }
        else {
            winnerIdsForBelow.push_back(
                -1
            );
        }
    }

    vector<int> returnedWinsFromAbove(
        candidatesToAbove.size(),
        -1
    );
    vector<int> returnedWinsFromBelow(
        candidatesToBelow.size(),
        -1
    );

    // Send winner decisions to candidate source ranks
    exchangeNeighborDataKnownCounts(
        winnerIdsForAbove,
        winnerIdsForBelow,
        returnedWinsFromAbove,
        returnedWinsFromBelow,
        MPI_INT,
        previousRank,
        nextRank,
        processCount,
        14,
        15
    );

    // Apply returned decisions for upward moves
    for (size_t i = 0; i < candidatesToAbove.size(); i++) {
        if (returnedWinsFromAbove[i] == candidatesToAbove[i].agentId) {
            moveResults[candidateIndexesToAbove[i]].won = 1;
        }
    }

    // Apply returned decisions for downward moves
    for (size_t i = 0; i < candidatesToBelow.size(); i++) {
        if (returnedWinsFromBelow[i] == candidatesToBelow[i].agentId) {
            moveResults[candidateIndexesToBelow[i]].won = 1;
        }
    }
}

// Move agents to ranks that own their new rows
void migrateAgents(
    vector<Agent>& agents,
    vector<int>& occupantIds,
    int localStartRow,
    int localRowCount,
    int gridWidth,
    int gridHeight,
    int rank,
    int processCount
) {
    static MPI_Datatype agentType = MPI_DATATYPE_NULL;

    if (agentType == MPI_DATATYPE_NULL) {
        agentType = createAgentType();
    }

    int localEndRow = localStartRow + localRowCount;
    (void)localEndRow;

    int previousRank = rank - 1;
    int nextRank = rank + 1;

    // Split agents by row ownership
    vector<Agent> agentsToSendUp;
    vector<Agent> agentsToSendDown;
    vector<Agent> localAgentsToKeep;

    // Movement is bounded by vision, but owner check is general still
    for (const auto& agent : agents) {
        int ownerRank = rankForRow(
            agent.row,
            gridHeight,
            processCount
        );

        if (ownerRank == rank) {
            localAgentsToKeep.push_back(
                agent
            );
        }
        else if (ownerRank < rank) {
            agentsToSendUp.push_back(
                agent
            );
        }
        else {
            agentsToSendDown.push_back(
                agent
            );
        }
    }

    // Exchange migration counts first
    int sendUpCount = (int)agentsToSendUp.size();
    int sendDownCount = (int)agentsToSendDown.size();
    int recvUpCount = 0;
    int recvDownCount = 0;

    MPI_Request sizeRequests[4];
    int sizeRequestCount = 0;

    if (previousRank >= 0) {
        MPI_Isend(
            &sendUpCount,
            1,
            MPI_INT,
            previousRank,
            20,
            MPI_COMM_WORLD,
            &sizeRequests[sizeRequestCount++]
        );
        MPI_Irecv(
            &recvUpCount,
            1,
            MPI_INT,
            previousRank,
            21,
            MPI_COMM_WORLD,
            &sizeRequests[sizeRequestCount++]
        );
    }

    if (nextRank < processCount) {
        MPI_Isend(
            &sendDownCount,
            1,
            MPI_INT,
            nextRank,
            21,
            MPI_COMM_WORLD,
            &sizeRequests[sizeRequestCount++]
        );
        MPI_Irecv(
            &recvDownCount,
            1,
            MPI_INT,
            nextRank,
            20,
            MPI_COMM_WORLD,
            &sizeRequests[sizeRequestCount++]
        );
    }

    if (sizeRequestCount > 0) {
        MPI_Waitall(
            sizeRequestCount,
            sizeRequests,
            MPI_STATUSES_IGNORE
        );
    }

    // Receive migrants from near slabs
    vector<Agent> receivedAgentsFromUp(
        recvUpCount
    );
    vector<Agent> receivedAgentsFromDown(
        recvDownCount
    );

    MPI_Request dataRequests[4];
    int dataRequestCount = 0;

    if (previousRank >= 0) {
        if (sendUpCount > 0) {
            MPI_Isend(
                agentsToSendUp.data(),
                sendUpCount,
                agentType,
                previousRank,
                22,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }

        if (recvUpCount > 0) {
            MPI_Irecv(
                receivedAgentsFromUp.data(),
                recvUpCount,
                agentType,
                previousRank,
                23,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }
    }

    if (nextRank < processCount) {
        if (sendDownCount > 0) {
            MPI_Isend(
                agentsToSendDown.data(),
                sendDownCount,
                agentType,
                nextRank,
                23,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }

        if (recvDownCount > 0) {
            MPI_Irecv(
                receivedAgentsFromDown.data(),
                recvDownCount,
                agentType,
                nextRank,
                22,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }
    }

    if (dataRequestCount > 0) {
        MPI_Waitall(
            dataRequestCount,
            dataRequests,
            MPI_STATUSES_IGNORE
        );
    }

    // Append migrants, then rebuild cell state for occupied cells
    agents = move(
        localAgentsToKeep
    );
    agents.insert(
        agents.end(),
        receivedAgentsFromUp.begin(),
        receivedAgentsFromUp.end()
    );
    agents.insert(
        agents.end(),
        receivedAgentsFromDown.begin(),
        receivedAgentsFromDown.end()
    );

    rebuildLocalOccupantIds(
        agents,
        occupantIds,
        localStartRow,
        localRowCount,
        gridWidth
    );
}

// Gather full state and write optional CSV output on rank 0
void gatherAndWriteCsv(
    const vector<Agent>& agents,
    const vector<int>& localSugar,
    const vector<int>& localCapacity,
    int localStartRow,
    int localRowCount,
    int gridWidth,
    int gridHeight,
    int rank,
    int processCount,
    const SimConfig& config,
    const char* outputPrefix,
    int stepNumber
) {
    (void)localCapacity;
    (void)localStartRow;

    // Skip expensive full gather unless CSV output is needed / was asked for
    if (!outputPrefix) {
        return;
    }

    static MPI_Datatype agentType = MPI_DATATYPE_NULL;

    if (agentType == MPI_DATATYPE_NULL) {
        agentType = createAgentType();
    }

    // Gather sugar grid by row slab
    vector<int> sugarReceiveCounts(
        processCount
    );
    vector<int> sugarDisplacements(
        processCount
    );

    int localSugarCount = localRowCount * gridWidth;

    MPI_Gather(
        &localSugarCount,
        1,
        MPI_INT,
        sugarReceiveCounts.data(),
        1,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );

    vector<int> fullSugarGrid;

    // Print benchmark summary once
    if (rank == 0) {
        sugarDisplacements[0] = 0;

        for (int i = 1; i < processCount; i++) {
            sugarDisplacements[i] =
                sugarDisplacements[i - 1] + sugarReceiveCounts[i - 1];
        }

        fullSugarGrid.resize(
            gridHeight * gridWidth
        );
    }

    MPI_Gatherv(
        localSugar.data(),
        localSugarCount,
        MPI_INT,
        fullSugarGrid.data(),
        sugarReceiveCounts.data(),
        sugarDisplacements.data(),
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );

    // Gather agents after sugar grid
    int localAgentCount = (int)agents.size();

    vector<int> agentReceiveCounts(
        processCount
    );
    vector<int> agentDisplacements(
        processCount
    );

    MPI_Gather(
        &localAgentCount,
        1,
        MPI_INT,
        agentReceiveCounts.data(),
        1,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );

    vector<Agent> allAgents;

    // Print benchmark summary once
    if (rank == 0) {
        agentDisplacements[0] = 0;

        for (int i = 1; i < processCount; i++) {
            agentDisplacements[i] =
                agentDisplacements[i - 1] + agentReceiveCounts[i - 1];
        }

        int totalAgentCount =
            agentDisplacements[processCount - 1] +
            agentReceiveCounts[processCount - 1];

        allAgents.resize(
            totalAgentCount
        );
    }

    MPI_Gatherv(
        agents.data(),
        localAgentCount,
        agentType,
        allAgents.data(),
        agentReceiveCounts.data(),
        agentDisplacements.data(),
        agentType,
        0,
        MPI_COMM_WORLD
    );

    if (rank != 0) {
        return;
    }

    // Write full sugar grid snapshot
    ostringstream gridFilePath;
    gridFilePath
        << outputPrefix
        << "_"
        << setfill(
            '0'
        )
        << setw(
            6
        )
        << stepNumber
        << "_GRID.csv";

    ofstream gridFile(
        gridFilePath.str()
    );
    gridFile << "row,col,sugarCount\n";

    for (int row = 0; row < gridHeight; row++) {
        for (int col = 0; col < gridWidth; col++) {
            gridFile
                << row
                << ","
                << col
                << ","
                << fullSugarGrid[row * gridWidth + col]
                << "\n";
        }
    }

    gridFile.close();

    // Write gathered agent snapshot
    ostringstream agentFilePath;
    agentFilePath
        << outputPrefix
        << "_"
        << setfill(
            '0'
        )
        << setw(
            6
        )
        << stepNumber
        << "_AGENT.csv";

    ofstream agentFile(
        agentFilePath.str()
    );
    agentFile << "agentID,row,col,wealth,metabolism,vision,age\n";

    for (const auto& agent : allAgents) {
        agentFile
            << agent.id << ","
            << agent.row << ","
            << agent.col << ","
            << agent.wealth << ","
            << agent.metabolism << ","
            << agent.vision << ","
            << agent.age
            << "\n";
    }

    agentFile.close();
}


int runSugarScapeBenchmark(
    int argc,
    char* argv[]
) {
    // Start MPI runtime
    MPI_Init(
        &argc,
        &argv
    );

    int rank;
    int processCount;

    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &rank
    );
    MPI_Comm_size(
        MPI_COMM_WORLD,
        &processCount
    );

    // Need height, width, and agent count
    if (argc < 4) {
        if (rank == 0) {
            cerr
                << "usage: mpirun -np P ./mpisugarscape height width agentCount"
                << " [timesteps [outputPrefix [csvInterval"
                << " [growthRate [metabolismMax]]]]]\n";
        }

        // Shut down MPI runtime
        MPI_Finalize();
        return -1;
    }

    // Parse simulation config
    SimConfig config;
    config.height = atoi(
        argv[1]
    );
    config.width = atoi(
        argv[2]
    );
    config.agentCount = atoi(
        argv[3]
    );
    config.timesteps = (
        argc > 4
            ? atoi(
                argv[4]
            )
            : 100
    );

    const char* outputPrefix = (
        argc > 5
            ? argv[5]
            : nullptr
    );

    config.csvInterval = (
        argc > 6
            ? atoi(
                argv[6]
            )
            : 1
    );

    if (config.csvInterval <= 0) {
        config.csvInterval = 1;
    }

    config.growthRate = (
        argc > 7
            ? atoi(
                argv[7]
            )
            : 1
    );
    config.sugarCapMin = 1;
    config.sugarCapMax = 4;
    config.metabolismMin = 1;
    config.metabolismMax = (
        argc > 8
            ? atoi(
                argv[8]
            )
            : 4
    );
    config.visionMin = 1;
    config.visionMax = 6;
    config.initialWealth = (unsigned int)-1;

    int gridHeight = config.height;
    int gridWidth = config.width;

    // Validate grid and population size
    if (
        gridHeight <= 0 ||
        gridWidth <= 0 ||
        config.agentCount > gridHeight * gridWidth
    ) {
        if (rank == 0) {
            cerr << "Invalid height, width, or agentCount\n";
        }

        // Shut down MPI runtime
        MPI_Finalize();
        return -1;
    }

    // Require enough rows per rank for halo movement
    if (processCount > 1 && gridHeight / processCount < 2 * config.visionMax) {
        if (rank == 0) {
            cerr
                << "Error: partition height "
                << gridHeight / processCount
                << " < 2*visionMax="
                << 2 * config.visionMax
                << ". Reduce process count or increase grid height.\n";
        }

        // Shut down MPI runtime
        MPI_Finalize();
        return -1;
    }

    auto totalStartTime = chrono::steady_clock::now();

    // Compute this rank's row slab
    int localStartRow = firstRow(
        rank,
        gridHeight,
        processCount
    );
    int localEndRow = endRow(
        rank,
        gridHeight,
        processCount
    );
    int localRowCount = localEndRow - localStartRow;

    // Allocate owned grid state
    vector<int> localSugar(
        localRowCount * gridWidth
    );
    vector<int> localCapacity(
        localRowCount * gridWidth
    );
    vector<int> occupantIds(
        localRowCount * gridWidth,
        -1
    );
    vector<Agent> localAgents;

    // Initialize sugar capacity and starting sugar
    for (int row = localStartRow; row < localEndRow; row++) {
        for (int col = 0; col < gridWidth; col++) {
            int localIndex = (row - localStartRow) * gridWidth + col;
            int capacity = sugarCapacity(
                row,
                col,
                gridHeight,
                gridWidth,
                config.sugarCapMin,
                config.sugarCapMax
            );

            localSugar[localIndex] = capacity;
            localCapacity[localIndex] = capacity;
        }
    }

    {
        // Rank 0 creates initial population
        vector<Agent> allAgents;

        if (rank == 0) {
            // Separate RNG streams preserve initialization order
            Xoshiro256 metabolismRng(
                2
            );
            Xoshiro256 visionRng(
                3
            );
            Xoshiro256 coordinateRng(
                4
            );
            Xoshiro256 wealthRng(
                100
            );

            allAgents.resize(
                config.agentCount
            );

            for (int i = 0; i < config.agentCount; i++) {
                Agent& agent = allAgents[i];
                agent.id = i;
                agent.metabolism = (int)metabolismRng.nextInRange(
                    config.metabolismMin,
                    config.metabolismMax
                );
                agent.vision = (int)visionRng.nextInRange(
                    config.visionMin,
                    config.visionMax
                );
                if (config.initialWealth == (unsigned int)-1) {
                    agent.wealth = (int)wealthRng.nextInRange(
                        agent.metabolism * 5,
                        (uint64_t)(agent.metabolism * 10)
                    );
                }
                else {
                    agent.wealth = (int)config.initialWealth;
                }

                agent.age = 0;
            }

            // Shuffle cell IDs to place agents without duplicates
            int totalCellCount = gridHeight * gridWidth;
            vector<int> shuffledCells(
                totalCellCount
            );

            for (int cellIndex = 0; cellIndex < totalCellCount; cellIndex++) {
                shuffledCells[cellIndex] = cellIndex;
            }

            for (int i = totalCellCount - 1; i > 0; i--) {
                int j = (int)coordinateRng.nextInRange(
                    0,
                    (uint64_t)i
                );

                if (j != i) {
                    swap(
                        shuffledCells[i],
                        shuffledCells[j]
                    );
                }
            }

            for (int i = 0; i < config.agentCount; i++) {
                int cellIndex = shuffledCells[i];
                allAgents[i].row = cellIndex / gridWidth;
                allAgents[i].col = cellIndex % gridWidth;
            }
        }

        // Distribute initialized agents by owning row slab
        static MPI_Datatype agentType = createAgentType();

        if (rank == 0) {
            vector<vector<Agent>> agentsByRank(
                processCount
            );

            for (const auto& agent : allAgents) {
                agentsByRank[rankForRow(
                    agent.row,
                    gridHeight,
                    processCount
                )].push_back(
                    agent
                );
            }

            for (
                int destinationRank = 1;
                destinationRank < processCount;
                destinationRank++
            ) {
                int sendCount = (int)agentsByRank[destinationRank].size();

                MPI_Send(
                    &sendCount,
                    1,
                    MPI_INT,
                    destinationRank,
                    40,
                    MPI_COMM_WORLD
                );

                if (sendCount > 0) {
                    MPI_Send(
                        agentsByRank[destinationRank].data(),
                        sendCount,
                        agentType,
                        destinationRank,
                        41,
                        MPI_COMM_WORLD
                    );
                }
            }

            localAgents = move(
                agentsByRank[0]
            );
        }
        else {
            int recvCount;

            MPI_Recv(
                &recvCount,
                1,
                MPI_INT,
                0,
                40,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE
            );

            localAgents.resize(
                recvCount
            );

            if (recvCount > 0) {
                MPI_Recv(
                    localAgents.data(),
                    recvCount,
                    agentType,
                    0,
                    41,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
            }
        }
    }

    // Initial deterministic processing order
    sortAgentsById(
        localAgents
    );

    rebuildLocalOccupantIds(
        localAgents,
        occupantIds,
        localStartRow,
        localRowCount,
        gridWidth
    );

    auto simulationStartTime = chrono::steady_clock::now();

    int maxVision = config.visionMax;

    vector<int> haloRowsAbove;
    vector<int> haloRowsBelow;
    vector<MPI_Request> haloRequests;
    vector<pair<int, int>> targetCells;
    vector<MoveResult> moveResults;

    // Track only harvested cells that still need regrowth
    vector<int> growingSugarCells;
    vector<unsigned char> sugarIsGrowing(
        (size_t)localRowCount * (size_t)gridWidth,
        0
    );

    // Add cell to sparse regrowth set once
    auto markSugarGrowing = [&](int localCellIndex) {
        if (!sugarIsGrowing[localCellIndex]) {
            sugarIsGrowing[localCellIndex] = 1;
            growingSugarCells.push_back(
                localCellIndex
            );
        }
    };

    // Main timestep loop
    for (int timestep = 0; timestep < config.timesteps; timestep++) {
        // Regrow only cells harvested in prior steps
        if (timestep > 0 && !growingSugarCells.empty()) {
            size_t writeIndex = 0;

            for (
                size_t readIndex = 0;
                readIndex < growingSugarCells.size();
                readIndex++
            ) {
                int localCellIndex = growingSugarCells[readIndex];
                int grownSugar = localSugar[localCellIndex] + config.growthRate;

                if (grownSugar >= localCapacity[localCellIndex]) {
                    localSugar[localCellIndex] = localCapacity[localCellIndex];
                    sugarIsGrowing[localCellIndex] = 0;
                }
                else {
                    localSugar[localCellIndex] = grownSugar;
                    growingSugarCells[writeIndex++] = localCellIndex;
                }
            }

            growingSugarCells.resize(
                writeIndex
            );
        }

        // Exchange sugar halos before deciding on movement 
        postHaloExchange(
            localSugar,
            haloRowsAbove,
            haloRowsBelow,
            localStartRow,
            localRowCount,
            gridWidth,
            maxVision,
            rank,
            processCount,
            gridHeight,
            haloRequests
        );

        int haloRowCountAbove = min(
            maxVision,
            localStartRow
        );
        int haloRowCountBelow = min(
            maxVision,
            gridHeight - localEndRow
        );

        if (!haloRequests.empty()) {
            MPI_Waitall(
                (int)haloRequests.size(),
                haloRequests.data(),
                MPI_STATUSES_IGNORE
            );
        }

        // Find one target per local agent
        targetCells.resize(
            localAgents.size()
        );

        for (
            size_t agentIndex = 0;
            agentIndex < localAgents.size();
            agentIndex++
        ) {
            const Agent& agent = localAgents[agentIndex];

            int targetRow;
            int targetCol;

            computeTargetWithHaloRows(
                agent.row,
                agent.col,
                agent.vision,
                gridHeight,
                gridWidth,
                localSugar,
                haloRowsAbove,
                haloRowsBelow,
                localStartRow,
                localRowCount,
                haloRowCountAbove,
                haloRowCountBelow,
                &targetRow,
                &targetCol
            );

            targetCells[agentIndex] = { targetRow, targetCol };
        }

        // Resolve target conflicts across local and boundary cells
        resolveConflictsDistributed(
            localAgents,
            targetCells,
            localSugar,
            occupantIds,
            localStartRow,
            localRowCount,
            gridWidth,
            gridHeight,
            rank,
            processCount,
            maxVision,
            moveResults
        );

        // Apply winning moves and harvest target sugar
        for (
            size_t agentIndex = 0;
            agentIndex < localAgents.size();
            agentIndex++
        ) {
            if (!moveResults[agentIndex].won) {
                continue;
            }

            Agent& agent = localAgents[agentIndex];
            int targetRow = moveResults[agentIndex].targetRow;
            int targetCol = moveResults[agentIndex].targetCol;

            int oldLocalRow = agent.row - localStartRow;

            if (
                oldLocalRow >= 0 &&
                oldLocalRow < localRowCount &&
                occupantIds[oldLocalRow * gridWidth + agent.col] == agent.id
            ) {
                occupantIds[oldLocalRow * gridWidth + agent.col] = -1;
            }

            agent.row = targetRow;
            agent.col = targetCol;

            int newLocalRow = targetRow - localStartRow;

            if (newLocalRow >= 0 && newLocalRow < localRowCount) {
                int localCellIndex = newLocalRow * gridWidth + targetCol;
                int harvestedSugar = localSugar[localCellIndex];

                agent.wealth += harvestedSugar;

                if (harvestedSugar > 0) {
                    localSugar[localCellIndex] = 0;
                    markSugarGrowing(
                        localCellIndex
                    );
                }

                occupantIds[localCellIndex] = agent.id;
            }
        }

        // Move agents to owner ranks
        migrateAgents(
            localAgents,
            occupantIds,
            localStartRow,
            localRowCount,
            gridWidth,
            gridHeight,
            rank,
            processCount
        );

        // Restore order before second harvest pass
        sortAgentsById(
            localAgents
        );

        // Harvest sugar under agents after migration
        for (auto& agent : localAgents) {
            int localRow = agent.row - localStartRow;

            if (
                localRow >= 0 &&
                localRow < localRowCount &&
                occupantIds[localRow * gridWidth + agent.col] == agent.id
            ) {
                int localCellIndex = localRow * gridWidth + agent.col;
                int harvestedSugar = localSugar[localCellIndex];

                if (harvestedSugar > 0) {
                    agent.wealth += harvestedSugar;
                    localSugar[localCellIndex] = 0;
                    markSugarGrowing(
                        localCellIndex
                    );
                }
            }
        }

        // Apply metabolism and compact surviving agents
        size_t survivorCount = 0;

        for (
            size_t agentIndex = 0;
            agentIndex < localAgents.size();
            agentIndex++
        ) {
            Agent& agent = localAgents[agentIndex];
            agent.wealth -= agent.metabolism;

            if (agent.wealth > 0) {
                if (survivorCount != agentIndex) {
                    localAgents[survivorCount] = agent;
                }

                survivorCount++;
            }
            else {
                int localRow = agent.row - localStartRow;

                if (localRow >= 0 && localRow < localRowCount) {
                    occupantIds[localRow * gridWidth + agent.col] = -1;
                }
            }
        }

        // Drop starved agents from local vector
        localAgents.resize(
            survivorCount
        );

        sortAgentsById(
            localAgents
        );

        // Optional debug/output path
        if (outputPrefix) {
            int outputStepNumber = timestep + 1;

            if (
                outputStepNumber == config.timesteps ||
                outputStepNumber % config.csvInterval == 0
                ) {
                gatherAndWriteCsv(
                    localAgents,
                    localSugar,
                    localCapacity,
                    localStartRow,
                    localRowCount,
                    gridWidth,
                    gridHeight,
                    rank,
                    processCount,
                    config,
                    outputPrefix,
                    outputStepNumber
                );
            }
        }
    }

    auto totalEndTime = chrono::steady_clock::now();

    // Compute local summary statistics
    int localLiveAgentCount = (int)localAgents.size();

    long long localTotalWealth = 0;

    for (const auto& agent : localAgents) {
        localTotalWealth += agent.wealth;
    }

    int localTotalSugar = 0;

    for (int i = 0; i < localRowCount * gridWidth; i++) {
        localTotalSugar += localSugar[i];
    }

    int globalLiveAgentCount = 0;
    long long globalTotalWealth = 0;
    int globalTotalSugar = 0;

    // Reduce final stats to rank 0
    MPI_Reduce(
        &localLiveAgentCount,
        &globalLiveAgentCount,
        1,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    MPI_Reduce(
        &localTotalWealth,
        &globalTotalWealth,
        1,
        MPI_LONG_LONG,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    MPI_Reduce(
        &localTotalSugar,
        &globalTotalSugar,
        1,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    // Print benchmark summary only once
    if (rank == 0) {
        cout
            << "BENCH_STAT agents=" << globalLiveAgentCount
            << " wealth=" << globalTotalWealth
            << " sugar=" << globalTotalSugar
            << endl;

        // Split initialization from timestep loop timing
        long initializationMs =
            chrono::duration_cast<chrono::milliseconds> (
                simulationStartTime - totalStartTime
            ).count();

        long simulationMs =
            chrono::duration_cast<chrono::milliseconds> (
                totalEndTime - simulationStartTime
            ).count();

        long totalMs =
            chrono::duration_cast<chrono::milliseconds> (
                totalEndTime - totalStartTime
            ).count();

        cout << "INIT_MS=" << initializationMs << endl;
        cout << "WALL_TIME_MS=" << simulationMs << endl;
        cout << "TOTAL_MS=" << totalMs << endl;

        // Optional debug/output path
        if (outputPrefix) {
            ostringstream summaryFilePath;
            summaryFilePath
                << outputPrefix
                << "_"
                << setfill(
                    '0'
                )
                << setw(
                    6
                )
                << config.timesteps
                << "_TIMESTEPCONSTANTS.csv";

            ofstream summaryFile(
                summaryFilePath.str()
            );

            summaryFile
                << "timestep,liveAgentCount,totalWealth,"
                << "averageWealth,totalSugarOnGrid\n";

            summaryFile
                << config.timesteps << ","
                << globalLiveAgentCount << ","
                << globalTotalWealth << ","
                << (
                    globalLiveAgentCount
                        ? (double)globalTotalWealth / globalLiveAgentCount
                        : 0
                )
                << ","
                << globalTotalSugar
                << "\n";

            summaryFile.close();
        }
    }

    // Shut down MPI runtime
    MPI_Finalize();
    return 0;
}
