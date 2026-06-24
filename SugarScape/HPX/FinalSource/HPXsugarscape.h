/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: SugarScape                                                      *
 * Library: HPX                                                               *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                    *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                           *
 *                                                                            *
 *****************************************************************************/

#ifndef HPX_SUGARSCAPE_H_
#define HPX_SUGARSCAPE_H_

#include <cstdint>
#include <utility>
#include <vector>

#define HPX_VISION_MAX 6

// Movement directions: N, S, E, W
static const std::pair<int, int> directions[4] = {
    {-1, 0},
    {1, 0},
    {0, 1},
    {0, -1}
};

// One local grid cell
struct GridCell {
    int currentSugar = 0;
    int sugarCapacity = 0;
    int agentID = -1;

    template <class Archive>
    void serialize(
        Archive& archive,
        unsigned int version
    ) {
        (void)version;
        archive & currentSugar & sugarCapacity & agentID;
    }
};

// One SugarScape agent
struct Agent {
    int id = 0;
    int row = 0;
    int col = 0;

    int wealth = 0;
    int metabolism = 0;
    int vision = 0;

    int age = 0;

    template <class Archive>
    void serialize(
        Archive& archive,
        unsigned int version
    ) {
        (void)version;
        archive & id & row & col & wealth & metabolism & vision & age;
    }
};

// Candidate movement target
struct MoveCandidate {
    int agentId = 0;
    int agentRow = 0;
    int agentCol = 0;

    int targetRow = 0;
    int targetCol = 0;

    template <class Archive>
    void serialize(
        Archive& archive,
        unsigned int version
    ) {
        (void)version;
        archive & agentId & agentRow & agentCol & targetRow & targetCol;
    }
};

// Boundary candidates sent between neighboring localities
struct BoundaryCandidatePack {
    std::vector<MoveCandidate> toAbove;
    std::vector<MoveCandidate> toBelow;

    template <class Archive>
    void serialize(
        Archive& archive,
        unsigned int version
    ) {
        (void)version;
        archive & toAbove & toBelow;
    }
};

// Boundary conflict replies sent back to source localities
struct BoundaryWinPack {
    std::vector<int> winsForAbove;
    std::vector<int> winsForBelow;

    template <class Archive>
    void serialize(
        Archive& archive,
        unsigned int version
    ) {
        (void)version;
        archive & winsForAbove & winsForBelow;
    }
};

// Migrating agents sent between neighboring localities
struct BoundaryAgentPack {
    std::vector<Agent> toUp;
    std::vector<Agent> toDown;

    template <class Archive>
    void serialize(
        Archive& archive,
        unsigned int version
    ) {
        (void)version;
        archive & toUp & toDown;
    }
};

// RNG for stable initialization
struct Xoshiro256 {
    uint64_t state[4];

    static uint64_t splitmix64(
        uint64_t& seedValue
    );

    static uint64_t rotateLeft(
        uint64_t value,
        int shift
    );

    explicit Xoshiro256(
        uint64_t seed
    );

    uint64_t next();

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

#endif
