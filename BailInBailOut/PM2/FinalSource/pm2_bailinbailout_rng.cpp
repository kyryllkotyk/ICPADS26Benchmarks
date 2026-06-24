/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Bail-In Bail-Out                                                *
 * Library: PM2/MPI                                                           *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "pm2_bailinbailout.h"

// Xoshiro random number generator helpers

// Expand one seed value into high quality generator state
uint64_t BailInBailOut::Xoshiro256::splitmix64(
    uint64_t& x
) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// Initialize xoshiro state from splitmix64 stream
BailInBailOut::Xoshiro256::Xoshiro256(
    uint64_t seed
) {
    for (int i = 0; i < 4; i++) {
        s[i] = splitmix64(
            seed
        );
    }
}

// Rotate helper used by xoshiro transition
uint64_t BailInBailOut::Xoshiro256::rotl(
    uint64_t x,
    int k
) {
    return (
        x << k
    ) | (x >> (
        64 - k
    ));
}

// Advance generator state and return next 64 bit value
uint64_t BailInBailOut::Xoshiro256::next() {
    uint64_t result = rotl(
        s[1] * 5,
        7
    ) * 9;
    uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rotl(
        s[3],
        45
    );

    return result;
}

// Draw unbiased integer value in inclusive range
uint64_t BailInBailOut::Xoshiro256::nextInRange(
    uint64_t min,
    uint64_t max
) {
    uint64_t range = max - min + 1;
    uint64_t x;
    uint64_t limit = numeric_limits<uint64_t>::max() - (numeric_limits
        <uint64_t>::max() % range);

    // Rejection sampling avoids modulo bias for non power of 2 ranges
    do {
        x = next();
    } while (x >= limit);

    return min + (x % range);
}

// Convert high bits into double in [0, 1)
double BailInBailOut::Xoshiro256::nextDouble01() {
    return (
        next() >> 11
    ) * (1.0 / 9007199254740992.0); // 2^53
}

// Draw deterministic floating point value in range
double BailInBailOut::Xoshiro256::nextInRangeDouble(
    double min,
    double max
) {
    return min + (max - min) * nextDouble01();
}

// Deterministic seeding and selection helpers

// One shot splitmix64 hash for deterministic entity mapping
uint64_t BailInBailOut::splitmix64Hash(
    uint64_t x
) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// Assign worker to bank or firm without central coordination
unsigned int BailInBailOut::selectForWorker(
    uint64_t baseSeed,
    unsigned int workerGlobalId,
    unsigned int countTotal,
    uint64_t streamTag
) {

    // Worker assignment uses streamTag to separate bank and firm mapping
    uint64_t seed = makeSeed(
        baseSeed,
        0,
        0,
        workerGlobalId,
        streamTag
    );
    BailInBailOut::Xoshiro256 rng(
        seed
    );

    // Bank IDs are global
    return (
        unsigned int
    )rng.nextInRange(
        0,
        (uint64_t)countTotal - 1
    );
}


// Combine run, timestep, entity id, and stream into one seed
uint64_t BailInBailOut::makeSeed(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int timestep,
    unsigned int globalId,
    uint64_t streamTag
) {

    uint64_t key = baseSeed;

    // Mix in run number
    key ^= 0xD6E8FEB86659FD93ULL * (uint64_t)run;

    // Mix in timestep (0 for static parameters like graph/structure)
    key ^= 0xBF58476D1CE4E5B9ULL * (uint64_t)timestep;

    // Mix in entity global ID (bank/worker/firm)
    key ^= 0x9E3779B97F4A7C15ULL * (uint64_t)globalId;

    // Mix in stream tag to separate RNG
    key ^= 0x94D049BB133111EBULL * streamTag;

    return splitmix64Hash(
        key
    );
}

unsigned short BailInBailOut::percent0to100FromMixedSeed(
    uint64_t seed
) {
    // Assumes mixed seed

    uint64_t x = seed;

    const uint64_t range = 101;
    const uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    // Retry by adding an odd constant
    while (x >= limit) {
        x += 0x9e3779b97f4a7c15ULL;
    }

    return (
        unsigned short
    )(x % range);
}

// Deterministically draw integer value from seed without storing RNG
uint64_t BailInBailOut::randomInRangeFromSeed(
    uint64_t seed,
    uint64_t minVal,
    uint64_t maxVal
) {
    // Assumes mixed seed

    if (minVal > maxVal) {
        swap(
            minVal,
            maxVal
        );
    }

    uint64_t range = maxVal - minVal + 1;
    uint64_t x = seed;

    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    // Keep going until retry
    while (x >= limit) {
        x += 0x9e3779b97f4a7c15ULL;
    }

    return minVal + (x % range);
}
