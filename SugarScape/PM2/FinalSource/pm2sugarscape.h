/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                              *
 *                                                                            *
 * Benchmark: SugarScape                                                      *
 * Library: PM2/MPI                                                           *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Final Formatting and Cleanup: Kyryll Kotyk                                 *
 *                                                                            *
 *****************************************************************************/

#ifndef PM2_SUGARSCAPE_H_
#define PM2_SUGARSCAPE_H_

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

using namespace std;

// Movement directions: N, S, E, W
static const int directions[4][2] = {
    {-1, 0},
    {1, 0},
    {0, 1},
    {0, -1}
};

// One agent in the simulation
struct Agent {
    int id;
    int row;
    int col;

    int wealth;
    int metabolism;
    int vision;

    int age;
};


// Candidate move sent to other ranks
struct MoveCandidate {
    int agentId;
    int agentRow;
    int agentCol;

    int targetRow;
    int targetCol;
};

// Result of a move after conflict resolution
struct MoveResult {
    int agentId;
    int targetRow;
    int targetCol;

    // 1 = won, 0 = lost
    int won;
};

// Simulation config (same for all ranks)
struct SimConfig {
    int height;
    int width;
    int agentCount;
    int timesteps;

    int growthRate;

    int sugarCapMin;
    int sugarCapMax;

    int metabolismMin;
    int metabolismMax;

    int visionMin;
    int visionMax;

    // -1 = random initialization
    unsigned int initialWealth;

    int csvInterval;
};

// RNG for determinism regardless of rank and node configuration
struct Xoshiro256 {
    uint64_t state[4];

    // Seed expansion
    static uint64_t splitmix64(
        uint64_t& seedValue
    );

    // Bit rotation
    static uint64_t rotateLeft(
        uint64_t value,
        int shift
    );

    explicit Xoshiro256(
        uint64_t seed
    );

    // Next random value
    uint64_t next();

    // Uniform random in [minValue, maxValue]
    uint64_t nextInRange(
        uint64_t minValue,
        uint64_t maxValue
    );
};

// Run benchmark from command line arguments
int runSugarScapeBenchmark(
    int argc,
    char* argv[]
);

// Map a global row to owning rank
int rankForRow(
    int globalRow,
    int gridHeight,
    int processCount
);

// First row owned by the rank
int firstRow(
    int rank,
    int gridHeight,
    int processCount
);

// One-past-last row owned by the rank
int endRow(
    int rank,
    int gridHeight,
    int processCount
);

// Sugar capacity based on position
int sugarCapacity(
    int row,
    int col,
    int gridHeight,
    int gridWidth,
    int minCapacity,
    int maxCapacity
);

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
);

// Resolve move conflicts
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
);

// Move agents to correct rank
void migrateAgents(
    vector<Agent>& agents,
    vector<int>& occupantIds,
    int localStartRow,
    int localRowCount,
    int gridWidth,
    int gridHeight,
    int rank,
    int processCount
);

// Gather full state and write CSV on rank 0
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
);

// MPI datatypes for structs
MPI_Datatype createAgentType();
MPI_Datatype createMoveCandidateType();

// Read sugar from owned slab or halo rows
inline int sugarAtWithHaloRows(
    int globalRow,
    int globalCol,
    int localStartRow,
    int localRowCount,
    int gridWidth,
    const vector<int>& localSugar,
    const vector<int>& haloRowsAbove,
    const vector<int>& haloRowsBelow,
    int haloRowCountAbove,
    int haloRowCountBelow
) {
    // Reject invalid column before row lookup
    if (globalCol < 0 || globalCol >= gridWidth) {
        return -1;
    }

    int localEndRow = localStartRow + localRowCount;

    // Fast path for owned rows
    if (globalRow >= localStartRow && globalRow < localEndRow) {
        return localSugar[(globalRow - localStartRow) * gridWidth + globalCol];
    }

    // Check upper halo when row is above slab
    if (globalRow < localStartRow) {
        int haloStartRow = localStartRow - haloRowCountAbove;
        int haloRow = globalRow - haloStartRow;

        if (haloRow >= 0 && haloRow < haloRowCountAbove) {
            return haloRowsAbove[haloRow * gridWidth + globalCol];
        }

        return -1;
    }

    // Otherwise check lower halo
    int haloRow = globalRow - localEndRow;

    if (haloRow >= 0 && haloRow < haloRowCountBelow) {
        return haloRowsBelow[haloRow * gridWidth + globalCol];
    }

    return -1;
}

// Find best visible target without copying full halo buffer
inline void computeTargetWithHaloRows(
    int agentRow,
    int agentCol,
    int vision,
    int gridHeight,
    int gridWidth,
    const vector<int>& localSugar,
    const vector<int>& haloRowsAbove,
    const vector<int>& haloRowsBelow,
    int localStartRow,
    int localRowCount,
    int haloRowCountAbove,
    int haloRowCountBelow,
    int* outTargetRow,
    int* outTargetCol
) {
    // Staying in place is valid 
    // Staying in place starts as best target
    int bestRow = agentRow;
    int bestCol = agentCol;
    int bestSugar = sugarAtWithHaloRows(
        agentRow,
        agentCol,
        localStartRow,
        localRowCount,
        gridWidth,
        localSugar,
        haloRowsAbove,
        haloRowsBelow,
        haloRowCountAbove,
        haloRowCountBelow
    );
    int bestDistance = 0;
    int bestDirection = 4;

    // Check each movement direction in priority order
    for (int direction = 0; direction < 4; direction++) {
        int rowStep = directions[direction][0];
        int colStep = directions[direction][1];

        // Scan until vision limit or grid edge
        for (int distance = 1; distance <= vision; distance++) {
            int targetRow = agentRow + rowStep * distance;
            int targetCol = agentCol + colStep * distance;

            if (targetRow < 0 || targetRow >= gridHeight || targetCol < 0 ||
                targetCol >= gridWidth) {
                break;
            }

            int sugarValue = sugarAtWithHaloRows(
                targetRow,
                targetCol,
                localStartRow,
                localRowCount,
                gridWidth,
                localSugar,
                haloRowsAbove,
                haloRowsBelow,
                haloRowCountAbove,
                haloRowCountBelow
            );

            if (sugarValue < 0) {
                break;
            }

            // Pick highest sugar, then nearest cell, then earlier direction
            if (
                sugarValue > bestSugar ||
                (sugarValue == bestSugar && distance < bestDistance) ||
                (sugarValue == bestSugar && distance == bestDistance &&
                    direction < bestDirection)
                ) {
                bestSugar = sugarValue;
                bestDistance = distance;
                bestDirection = direction;
                bestRow = targetRow;
                bestCol = targetCol;
            }
        }
    }

    *outTargetRow = bestRow;
    *outTargetCol = bestCol;
}

// Restore local processing order
inline void sortAgentsById(
    vector<Agent>& agents
) {
    sort(
        agents.begin(),
        agents.end(),
        [](const Agent& leftAgent, const Agent& rightAgent) {
            return leftAgent.id < rightAgent.id;
        }
    );
}

// Rebuild occupied cell lookup for owned 
inline void rebuildLocalOccupantIds(
    const vector<Agent>& localAgents,
    vector<int>& occupantIds,
    int localStartRow,
    int localRowCount,
    int gridWidth
) {
    // Full clear for correctness 
    fill(
        occupantIds.begin(),
        occupantIds.end(),
        -1
    );

    for (const auto& agent : localAgents) {
        int localRow = agent.row - localStartRow;

        // Register only agents owned by this rank
        if (localRow >= 0 && localRow < localRowCount) {
            occupantIds[localRow * gridWidth + agent.col] = agent.id;
        }
    }
}

// Exchange vectors with adjacent ranks
template<typename T>
inline void exchangeNeighborVectors(
    const vector<T>& sendUp,
    const vector<T>& sendDown,
    vector<T>& recvUp,
    vector<T>& recvDown,
    MPI_Datatype dataType,
    int previousRank,
    int nextRank,
    int processCount,
    int countUpTag,
    int countDownTag,
    int dataUpTag,
    int dataDownTag
) {
    int sendUpCount = (int)sendUp.size();
    int sendDownCount = (int)sendDown.size();
    int recvUpCount = 0;
    int recvDownCount = 0;

    // Exchange sizes before payloads so they know what to expect
    MPI_Request countRequests[4];
    int countRequestCount = 0;

    if (previousRank >= 0) {
        MPI_Isend(
            &sendUpCount,
            1,
            MPI_INT,
            previousRank,
            countUpTag,
            MPI_COMM_WORLD,
            &countRequests[countRequestCount++]
        );

        MPI_Irecv(
            &recvUpCount,
            1,
            MPI_INT,
            previousRank,
            countDownTag,
            MPI_COMM_WORLD,
            &countRequests[countRequestCount++]
        );
    }

    if (nextRank < processCount) {
        MPI_Isend(
            &sendDownCount,
            1,
            MPI_INT,
            nextRank,
            countDownTag,
            MPI_COMM_WORLD,
            &countRequests[countRequestCount++]
        );
        MPI_Irecv(
            &recvDownCount,
            1,
            MPI_INT,
            nextRank,
            countUpTag,
            MPI_COMM_WORLD,
            &countRequests[countRequestCount++]
        );
    }

    if (countRequestCount > 0) {
        MPI_Waitall(
            countRequestCount,
            countRequests,
            MPI_STATUSES_IGNORE
        );
    }

    // Allocate receive buffers exactly as needed
    recvUp.resize(
        recvUpCount
    );
    recvDown.resize(
        recvDownCount
    );

    // Post data sends and receives now that the sizes are known
    MPI_Request dataRequests[4];
    int dataRequestCount = 0;

    if (previousRank >= 0) {
        if (sendUpCount > 0) {
            MPI_Isend(
                sendUp.data(),
                sendUpCount,
                dataType,
                previousRank,
                dataUpTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }

        if (recvUpCount > 0) {
            MPI_Irecv(
                recvUp.data(),
                recvUpCount,
                dataType,
                previousRank,
                dataDownTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }
    }

    if (nextRank < processCount) {
        if (sendDownCount > 0) {
            MPI_Isend(
                sendDown.data(),
                sendDownCount,
                dataType,
                nextRank,
                dataDownTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }

        if (recvDownCount > 0) {
            MPI_Irecv(
                recvDown.data(),
                recvDownCount,
                dataType,
                nextRank,
                dataUpTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }
    }

    // Wait...
    if (dataRequestCount > 0) {
        MPI_Waitall(
            dataRequestCount,
            dataRequests,
            MPI_STATUSES_IGNORE
        );
    }
}

// Exchange data given known receive buffer sizes
template<typename T>
inline void exchangeNeighborDataKnownCounts(
    const vector<T>& sendUp,
    const vector<T>& sendDown,
    vector<T>& recvUp,
    vector<T>& recvDown,
    MPI_Datatype dataType,
    int previousRank,
    int nextRank,
    int processCount,
    int dataUpTag,
    int dataDownTag
) {
    // Reuse caller receive buffers
    MPI_Request dataRequests[4];
    int dataRequestCount = 0;

    if (previousRank >= 0) {
        if (!sendUp.empty()) {
            MPI_Isend(
                sendUp.data(),
                (int)sendUp.size(),
                dataType,
                previousRank,
                dataUpTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }

        if (!recvUp.empty()) {
            MPI_Irecv(
                recvUp.data(),
                (int)recvUp.size(),
                dataType,
                previousRank,
                dataDownTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }
    }

    if (nextRank < processCount) {
        if (!sendDown.empty()) {
            MPI_Isend(
                sendDown.data(),
                (int)sendDown.size(),
                dataType,
                nextRank,
                dataDownTag,
                MPI_COMM_WORLD,
                &dataRequests[dataRequestCount++]
            );
        }

        if (!recvDown.empty()) {
            MPI_Irecv(
                recvDown.data(),
                (int)recvDown.size(),
                dataType,
                nextRank,
                dataUpTag,
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
}

#endif // PM2_SUGARSCAPE_H_
