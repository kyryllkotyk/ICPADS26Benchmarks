/******************************************************************************
 *                                                                            *
 * ICPADS2 Benchmarks Collection                                              *
 *                                                                            *
 * Benchmark: SugarScape                                                      *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef SUGAR_CHUNK_H_
#define SUGAR_CHUNK_H_

#include "MethodRegistry.h"
#include "Place.h"

#include <cstdint>
#include <utility>
#include <vector>

#define SUGARCHUNK_VISION_MAX 6

static const int SUGARCHUNK_DIRS[4][2] = {
    {-1, 0},
    {1, 0},
    {0, 1},
    {0, -1}
};

struct SugarXoshiro256 {
    uint64_t state[4];

    static inline uint64_t splitmix64(
        uint64_t& seedValue
    ) {
        seedValue += 0x9e3779b97f4a7c15ULL;
        uint64_t z = seedValue;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    static inline uint64_t rotateLeft(
        const uint64_t value,
        const int shift
    ) {
        return (value << shift) | (value >> (64 - shift));
    }

    inline explicit SugarXoshiro256(
        uint64_t seed
    ) {
        for (int i = 0; i < 4; i++) {
            state[i] = splitmix64(
                seed
            );
        }
    }

    inline uint64_t next() {
        const uint64_t result = rotateLeft(
            state[1] * 5,
            7
        ) * 9;
        const uint64_t nextState = state[1] << 17;

        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];
        state[2] ^= nextState;
        state[3] = rotateLeft(
            state[3],
            45
        );

        return result;
    }

    inline uint64_t nextInRange(
        const uint64_t minValue,
        const uint64_t maxValue
    ) {
        if (maxValue < minValue) {
            return next();
        }

        const uint64_t range = maxValue - minValue + 1;
        const uint64_t limit = UINT64_MAX - (UINT64_MAX % range);
        uint64_t value;

        do {
            value = next();
        } while (value >= limit);

        return minValue + value % range;
    }
};

struct AgentRec {
    int id;
    int row;
    int col;
    int wealth;
    int metabolism;
    int vision;
    int age;
};

struct MoveCandidateRec {
    int agentId;
    int agentRow;
    int agentCol;
    int targetRow;
    int targetCol;
};

struct SugarChunkConfig {
    int height;
    int width;
    int numChunks;
    int agentCount;

    int growthRate;
    int sugarCapMin;
    int sugarCapMax;
    int metabolismMin;
    int metabolismMax;
    int visionMin;
    int visionMax;

    // -1 means random initial wealth
    unsigned int initialWealth;

    int outMessageBytes;
    int reportAgentsBytes;
    int reportSugarInts;
};

struct SugarChunkSeedHeader {
    int numAgents;
};

struct SugarChunkStats {
    int localAgentCount;
    long long localWealth;
    int localSugar;
    int localMaxAgentId;
};

struct SugarChunkStepArg {
    int step;
};

class SugarChunk : public Place {
public:
    explicit SugarChunk(
        void* argument
    );

    MASS_DISPATCH_TABLE(
        SugarChunk,
        MASS_METHOD(SugarChunk, init),
        MASS_METHOD(SugarChunk, growSugar),
        MASS_METHOD(SugarChunk, packHaloSugar),
        MASS_METHOD(SugarChunk, computeTargets),
        MASS_METHOD(SugarChunk, packMoveCandidates),
        MASS_METHOD(SugarChunk, resolveConflictsLocal),
        MASS_METHOD(SugarChunk, packWinnerReplies),
        MASS_METHOD(SugarChunk, applyMovesAndTakeSugar),
        MASS_METHOD(SugarChunk, packMigrations),
        MASS_METHOD(SugarChunk, integrateMigrations),
        MASS_METHOD(SugarChunk, takeSugarAtCurrent),
        MASS_METHOD(SugarChunk, metabolismAndDeath),
        MASS_METHOD(SugarChunk, sortLocalAgents),
        MASS_METHOD(SugarChunk, setTimestep),
        MASS_METHOD(SugarChunk, advanceTimestep),
        MASS_METHOD(SugarChunk, reportStats),
        MASS_METHOD(SugarChunk, reportAgentsCount),
        MASS_METHOD(SugarChunk, reportAgents),
        MASS_METHOD(SugarChunk, reportSugar),
        MASS_METHOD(SugarChunk, seedInitialAgents),
        MASS_METHOD(SugarChunk, redistributeSeedAgents)
    )

    void* init(
        void* argument
    );

    void* growSugar(
        void* argument
    );

    void* packHaloSugar(
        void* argument
    );

    void* computeTargets(
        void* argument
    );

    void* packMoveCandidates(
        void* argument
    );

    void* resolveConflictsLocal(
        void* argument
    );

    void* packWinnerReplies(
        void* argument
    );

    void* applyMovesAndTakeSugar(
        void* argument
    );

    void* packMigrations(
        void* argument
    );

    void* integrateMigrations(
        void* argument
    );

    void* takeSugarAtCurrent(
        void* argument
    );

    void* metabolismAndDeath(
        void* argument
    );

    void* sortLocalAgents(
        void* argument
    );

    void* setTimestep(
        void* argument
    );

    void* advanceTimestep(
        void* argument
    );

    void* reportStats(
        void* argument
    );

    void* reportAgentsCount(
        void* argument
    );

    void* reportAgents(
        void* argument
    );

    void* reportSugar(
        void* argument
    );

    void* seedInitialAgents(
        void* argument
    );

    void* redistributeSeedAgents(
        void* argument
    );

private:
    void ensureOutMessage();

    void resetOutMessageHeader(
        int phase
    );

    void writeOutU8(
        int offset,
        int value
    );

    int globalRowOf(
        int localRow
    ) const {
        return localStartRow_ + localRow;
    }

    int localRowOf(
        int globalRow
    ) const {
        return globalRow - localStartRow_;
    }

    bool ownsRow(
        int globalRow
    ) const;

    int ownerChunkOfRow(
        int globalRow
    ) const;

    void rebuildOccupantIndex();

    void clearOccupantAt(
        int globalRow,
        int col
    );

    void setOccupantAt(
        int globalRow,
        int col,
        int agentId
    );

    int getOccupantAt(
        int globalRow,
        int col
    ) const;

    void computeTargetFor(
        const AgentRec& agent,
        const std::vector<int>& sugarBuffer,
        int sugarBufferStartRow,
        int sugarBufferRowCount,
        int& outRow,
        int& outCol
    ) const;

    void markSugarGrowing(
        int localCellIndex
    );

    void touchConflictCell(
        int localCellIndex
    );

    int chunkIndex_;
    int numChunks_;
    int gridHeight_;
    int gridWidth_;
    int localStartRow_;
    int localRowCount_;
    int previousChunk_;
    int nextChunk_;

    int growthRate_;
    int sugarCapMin_;
    int sugarCapMax_;
    int metabolismMin_;
    int metabolismMax_;
    int visionMin_;
    int visionMax_;
    unsigned int initialWealth_;
    int outMessageBytes_;

    std::vector<int> localSugar_;
    std::vector<uint8_t> localCapacity_;
    std::vector<int> occupantIds_;
    std::vector<int> growingSugarCells_;
    std::vector<unsigned char> sugarIsGrowing_;

    std::vector<AgentRec> localAgents_;

    int currentStep_;
    bool compoundStepPrimed_;

    std::vector<int> bestLocalAgentIndex_;
    std::vector<int> bestAgentId_;
    std::vector<int> cellStamp_;
    int currentStamp_;

    std::vector<std::pair<int, int>> targetCells_;
    std::vector<int> moveWon_;
    std::vector<MoveCandidateRec> sentCandidatesAbove_;
    std::vector<MoveCandidateRec> sentCandidatesBelow_;
    std::vector<int> sentCandidateIdxAbove_;
    std::vector<int> sentCandidateIdxBelow_;
    std::vector<int> winsForFromAbove_;
    std::vector<int> winsForFromBelow_;
    std::vector<MoveCandidateRec> recvCandidatesFromAbove_;
    std::vector<MoveCandidateRec> recvCandidatesFromBelow_;
    std::vector<AgentRec> pendingMigrationsAbove_;
    std::vector<AgentRec> pendingMigrationsBelow_;

    SugarChunkStats reportStatsBuf_{};
    int reportCountBuf_{0};
    std::vector<char> reportAgentsBuf_;
    std::vector<int> reportSugarBuf_;
    int reportAgentsCapacityBytes_{0};
    int reportSugarCapacityInts_{0};
};

#endif
