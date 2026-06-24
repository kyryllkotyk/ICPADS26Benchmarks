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

// Compute floor(x * pct / 100) using integer arithmetic
uint64_t BailInBailOut::percentFloorU64(
    uint64_t x,
    unsigned short pct
) {
    return (
        x * (uint64_t)pct
    ) / 100ULL;
}

// Validate user parameters and normalize percent based ranges
bool BailInBailOut::errorDetectionAndClamping(
    const unsigned int runs,
    const unsigned int timesteps,
    const uint64_t baseSeed,

    const unsigned int bankCountTotal,
    const unsigned int firmCountTotal,
    const unsigned int workerCountTotal,

    const unsigned int initialBankLiquidity,
    const unsigned int initialFirmLiquidity,
    const unsigned int initialProductionCost,
    const unsigned int wage,

    unsigned short& wageConsumptionPercent,
    unsigned short& profitMultiplierMin,
    unsigned short& profitMultiplierMax,
    unsigned short& shockMultiplierMin,
    unsigned short& shockMultiplierMax,
    unsigned short& minInterestRate,
    unsigned short& maxInterestRate,
    unsigned short& firmShockMinPercent,
    unsigned short& firmShockMaxPercent,
    unsigned short& firmProfitMinPercent,
    unsigned short& firmProfitMaxPercent,
    unsigned short& firmOperatingCostMinPercent,
    unsigned short& firmOperatingCostMaxPercent,
    unsigned short& firmShockProfitMultiplierMinPercent,
    unsigned short& firmShockProfitMultiplierMaxPercent,
    unsigned short& firmInitialLiquidityMultiplierMinPercent,
    unsigned short& firmInitialLiquidityMultiplierMaxPercent,
    unsigned short& firmRareEventProbabilityPercent,
    unsigned short& firmRareEventImpactMinPercent,
    unsigned short& firmRareEventImpactMaxPercent,

    unsigned short& interbankDensity,
    const unsigned short maxInterbankLenderSamplingK,
    unsigned short& maxInterbankLoanPercent,
    unsigned short& maxFirmLoanPercent,
    unsigned short& firmLenderDegree,
    unsigned short& firmRepayPercent,
    unsigned short& bankRepayPercent,
    unsigned int interventionDelay,

    unsigned short& bailInCoveragePercent,
    unsigned short& bailOutCoveragePercent
) {
    (void)baseSeed;

    // Clamp all percentage settings into valid 0 to 100 range
    auto clampPercent = [](
        unsigned short value
    ) -> unsigned short {
        return (
            value > 100
        ) ? 100 : value;
        };

    wageConsumptionPercent = clampPercent(
        wageConsumptionPercent
    );
    profitMultiplierMin = clampPercent(
        profitMultiplierMin
    );
    profitMultiplierMax = clampPercent(
        profitMultiplierMax
    );
    shockMultiplierMin = clampPercent(
        shockMultiplierMin
    );
    shockMultiplierMax = clampPercent(
        shockMultiplierMax
    );

    interbankDensity = clampPercent(
        interbankDensity
    );
    maxInterbankLoanPercent = clampPercent(
        maxInterbankLoanPercent
    );
    maxFirmLoanPercent = clampPercent(
        maxFirmLoanPercent
    );

    firmRepayPercent = clampPercent(
        firmRepayPercent
    );
    bankRepayPercent = clampPercent(
        bankRepayPercent
    );

    bailInCoveragePercent = clampPercent(
        bailInCoveragePercent
    );
    bailOutCoveragePercent = clampPercent(
        bailOutCoveragePercent
    );

    firmShockMinPercent = clampPercent(
        firmShockMinPercent
    );
    firmShockMaxPercent = clampPercent(
        firmShockMaxPercent
    );
    firmProfitMinPercent = clampPercent(
        firmProfitMinPercent
    );
    firmProfitMaxPercent = clampPercent(
        firmProfitMaxPercent
    );
    firmOperatingCostMinPercent = clampPercent(
        firmOperatingCostMinPercent
    );
    firmOperatingCostMaxPercent = clampPercent(
        firmOperatingCostMaxPercent
    );
    firmShockProfitMultiplierMinPercent =
        clampPercent(
            firmShockProfitMultiplierMinPercent
        );
    firmShockProfitMultiplierMaxPercent =
        clampPercent(
            firmShockProfitMultiplierMaxPercent
        );
    firmInitialLiquidityMultiplierMinPercent =
        clampPercent(
            firmInitialLiquidityMultiplierMinPercent
        );
    firmInitialLiquidityMultiplierMaxPercent =
        clampPercent(
            firmInitialLiquidityMultiplierMaxPercent
        );
    firmRareEventProbabilityPercent =
        clampPercent(
            firmRareEventProbabilityPercent
        );
    firmRareEventImpactMinPercent = clampPercent(
        firmRareEventImpactMinPercent
    );
    firmRareEventImpactMaxPercent = clampPercent(
        firmRareEventImpactMaxPercent
    );

    // Normalize min/max inputs that may have been passed in reversed order
    auto swapIfMinGreaterThanMax = [](
        unsigned short& minVal,
        unsigned short& maxVal
    ) {
        if (minVal > maxVal) {
            unsigned short tmp = minVal;
            minVal = maxVal;
            maxVal = tmp;
        }
    };


    swapIfMinGreaterThanMax(
        profitMultiplierMin,
        profitMultiplierMax
    );
    swapIfMinGreaterThanMax(
        shockMultiplierMin,
        shockMultiplierMax
    );
    swapIfMinGreaterThanMax(
        minInterestRate,
        maxInterestRate
    );

    swapIfMinGreaterThanMax(
        firmShockMinPercent,
        firmShockMaxPercent
    );
    swapIfMinGreaterThanMax(
        firmProfitMinPercent,
        firmProfitMaxPercent
    );
    swapIfMinGreaterThanMax(
        firmOperatingCostMinPercent,
        firmOperatingCostMaxPercent
    );
    swapIfMinGreaterThanMax(
        firmShockProfitMultiplierMinPercent,
        firmShockProfitMultiplierMaxPercent
    );
    swapIfMinGreaterThanMax(
        firmInitialLiquidityMultiplierMinPercent,
        firmInitialLiquidityMultiplierMaxPercent
    );
    swapIfMinGreaterThanMax(
        firmRareEventImpactMinPercent,
        firmRareEventImpactMaxPercent
    );

    // Reject invalid sizes after clamping non fatal percent ranges
    if (runs < 1) {
        cerr << "runs variable out of bounds (Must be at least 1)";
        return false;
    }

    if (timesteps < 1) {
        cerr << "timesteps variable out of bounds (Must be at least 1)";
        return false;
    }

    if (bankCountTotal < 1) {
        cerr << "bankCountTotal variable out of bounds (Must be at least 1)";
        return false;
    }

    if (firmCountTotal < 1) {
        cerr << "firmCountTotal variable out of bounds (Must be at least 1)";
        return false;
    }

    if (workerCountTotal < 1) {
        cerr << "workerCountTotal variable out of bounds (Must be at least 1)";
        return false;
    }

    if (initialBankLiquidity < 1) {
        cerr << "initialBankLiquidity variable out of bounds (Must be at least 1)";
        return false;
    }

    if (initialFirmLiquidity < 1) {
        cerr << "initialFirmLiquidity variable out of bounds (Must be at least 1)";
        return false;
    }

    if (initialProductionCost < 1) {
        cerr << "initialProductionCost variable out of bounds (Must be at least 1)";
        return false;
    }

    if (wage < 1) {
        cerr << "wage variable out of bounds (Must be at least 1)";
        return false;
    }

    if (maxInterbankLenderSamplingK < 1) {
        cerr << "maxInterbankLenderSamplingK variable out of bounds (Must be at least 1)";
        return false;
    }

    if (firmLenderDegree < 1) {
        cerr << "firmLenderDegree variable out of bounds (Must be at least 1)";
        return false;
    }

    // Avoid obviously unreasonable allocation sizes
    if (bankCountTotal > 10000000) {
        cerr << "bankCountTotal too large (Unreasonable allocation size)";
        return false;
    }

    if (timesteps > 1000000) {
        cerr << "timesteps too large (Not in scope of the simulation)";
        return false;
    }

    // Reject configurations where no money can move through the model
    if (wageConsumptionPercent == 0 &&
        interbankDensity == 0 &&
        maxInterbankLoanPercent == 0 &&
        maxFirmLoanPercent == 0 &&
        firmRepayPercent == 0 &&
        bankRepayPercent == 0 &&
        bailInCoveragePercent == 0 &&
        bailOutCoveragePercent == 0) {
        cerr << "configuration out of bounds (No transfers possible. Simulation will be degenerate)";
        return false;
    }

    return true;
}

// Return owned count and global start for one rank
pair<unsigned int, unsigned int> BailInBailOut::computeRange(
    unsigned int totalCount,
    unsigned int rank,
    unsigned int totalRanks
) {

    unsigned int base = totalCount / totalRanks;
    unsigned int remainder = totalCount % totalRanks;

    pair<unsigned int, unsigned int> range;

    if (rank < remainder) {
        // Count
        range.first = base + 1;
        // Start
        range.second = rank * (base + 1);
    }
    else {
        range.first = base;
        range.second = remainder * (base + 1) + (rank - remainder) * base;
    }

    return range;
}


// Invert computeRange ownership mapping from global id to rank
unsigned int BailInBailOut::ownerRankFromGlobalId(
    unsigned int countTotal,
    unsigned int totalRanks,
    unsigned int globalId
) {

    // Inverse mapping (global ID to owning rank, matches computeRange splitting)

    unsigned int base = countTotal / totalRanks;
    unsigned int remainder = countTotal % totalRanks;

    // Give +1 remainder to the first ranks until no remainder is left

    unsigned int cutoff = (base + 1) * remainder;

    if (globalId < cutoff) {
        return globalId / (base + 1);
    }
    else {
        return remainder + (globalId - cutoff) / base;
    }
}

// Reconstruct deterministic lender interest rate from id and seed
unsigned int BailInBailOut::getInterestRate(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int globalId,
    uint64_t minVal,
    uint64_t maxVal
) {
    uint64_t seed = makeSeed(
        baseSeed,
        run,
        0,
        globalId,
        STREAM_INTEREST_RATE_SELECTION
    );

    // Rebuild interest rates from seed to avoid MPI calls
    unsigned short rate = (unsigned short)randomInRangeFromSeed
    (
        seed,
        minVal,
        maxVal
    );

    return rate;
}

// Generate deterministic percent value for one entity and stream
unsigned short BailInBailOut::randomPercentFromStream(
    const uint64_t baseSeed,
    const unsigned int run,
    const unsigned int step,
    const unsigned int globalId,
    const unsigned int streamTag,
    const unsigned short minPercent,
    const unsigned short maxPercent
) {
    uint64_t seed = makeSeed(
        baseSeed,
        run,
        step,
        globalId,
        streamTag
    );

    return (
        unsigned short
    )randomInRangeFromSeed(
        seed,
        minPercent,
        maxPercent
    );
}

// Debt aggregation helper

// Sum all outstanding debt entries for one borrower
uint64_t BailInBailOut::sumDebtU64(
    const vector<BailInBailOut::DebtEntry>& debts
) {
    uint64_t s = 0;
    for (const auto& e : debts) {
        s += e.amount;
    }
    return s;
}
