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

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <mpi.h>

// Check whether a boolean command line flag is present
static bool hasFlag(
    int argc,
    char** argv,
    const char* flag
) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(
            argv[i],
            flag
        ) == 0) {
            return true;
        }
    }
    return false;
}

// Return value following a command line option
static const char* getArgValue(
    int argc,
    char** argv,
    const char* flag
) {
    for (int i = 1; i < argc - 1; i++) {
        if (strcmp(
            argv[i],
            flag
        ) == 0) {
            return argv[i + 1];
        }
    }
    return nullptr;
}

// Parse unsigned integer option or keep default
static unsigned int getUIntArg(
    int argc,
    char** argv,
    const char* flag,
    unsigned int defaultValue
) {
    const char* value = getArgValue(
        argc,
        argv,
        flag
    );
    if (!value) {
        return defaultValue;
    }
    return static_cast<unsigned int>(
        strtoul(
            value,
            nullptr,
            10
        )
    );
}

// Parse 64 bit unsigned option or keep default
static uint64_t getU64Arg(
    int argc,
    char** argv,
    const char* flag,
    uint64_t defaultValue
) {
    const char* value = getArgValue(
        argc,
        argv,
        flag
    );
    if (!value) {
        return defaultValue;
    }
    return static_cast<uint64_t>(
        strtoull(
            value,
            nullptr,
            10
        )
    );
}

// Parse unsigned short option with range checking
static unsigned short getUShortArg(
    int argc,
    char** argv,
    const char* flag,
    unsigned short defaultValue
) {
    const char* value = getArgValue(
        argc,
        argv,
        flag
    );
    if (!value) {
        return defaultValue;
    }

    errno = 0;
    char* endPtr = nullptr;
    unsigned long parsed = strtoul(
        value,
        &endPtr,
        10
    );

    if (errno != 0 || endPtr == value || *endPtr != '\0' ||
        parsed > numeric_limits<unsigned short>::max()) {
        cerr << "Invalid value for " << flag << ": " << value << "\n";
        MPI_Abort(
            MPI_COMM_WORLD,
            1
        );
    }

    return static_cast<unsigned short>(
        parsed
    );
}

// Parse bool option from common string forms
static bool getBoolArg(
    int argc,
    char** argv,
    const char* flag,
    bool defaultValue
) {
    const char* value = getArgValue(
        argc,
        argv,
        flag
    );

    if (!value) {
        return hasFlag(
            argc,
            argv,
            flag
        ) ? true : defaultValue;
    }

    if (
        strcmp(
            value,
            "1"
        ) == 0 ||
        strcmp(
            value,
            "true"
        ) == 0 ||
        strcmp(
            value,
            "yes"
        ) == 0 ||
        strcmp(
            value,
            "on"
        ) == 0
        ) {
        return true;
    }

    if (
        strcmp(
            value,
            "0"
        ) == 0 ||
        strcmp(
            value,
            "false"
        ) == 0 ||
        strcmp(
            value,
            "no"
        ) == 0 ||
        strcmp(
            value,
            "off"
        ) == 0
        ) {
        return false;
    }

    return defaultValue;
}

// Print launcher usage on rank 0 only
static void printUsage(
    int mpiRank
) {
    if (mpiRank != 0) {
        return;
    }

    printf(
        "Usage:\n"
        "  mpirun ... ./pm2bailinbailout [options]\n\n"
        "Core options:\n"
        "  --runs <uint>\n"
        "  --timesteps <uint>\n"
        "  --seed <uint64>\n\n"
        "Structure:\n"
        "  --banks <uint>\n"
        "  --firms <uint>\n"
        "  --workers <uint>\n"
        "  --bankWorkers <uint>\n\n"
        "System:\n"
        "  --initialBankLiquidity <uint>\n"
        "  --initialFirmLiquidity <uint>\n"
        "  --initialProductionCost <uint>\n"
        "  --wage <uint>\n"
        "  --bankEmployeeWage <uint>\n\n"
        "General percentages:\n"
        "  --wageConsumptionPercent <ushort>\n"
        "  --profitMultiplierMin <ushort>\n"
        "  --profitMultiplierMax <ushort>\n"
        "  --shockMultiplierMin <ushort>\n"
        "  --shockMultiplierMax <ushort>\n"
        "  --minInterestRate <ushort>\n"
        "  --maxInterestRate <ushort>\n\n"
        "Firm heterogeneity:\n"
        "  --firmShockMinPercent <ushort>\n"
        "  --firmShockMaxPercent <ushort>\n"
        "  --firmProfitMinPercent <ushort>\n"
        "  --firmProfitMaxPercent <ushort>\n"
        "  --firmOperatingCostMinPercent <ushort>\n"
        "  --firmOperatingCostMaxPercent <ushort>\n"
        "  --firmShockProfitMultiplierMinPercent <ushort>\n"
        "  --firmShockProfitMultiplierMaxPercent <ushort>\n"
        "  --firmInitialLiquidityMultiplierMinPercent <ushort>\n"
        "  --firmInitialLiquidityMultiplierMaxPercent <ushort>\n"
        "  --firmRareEventProbabilityPercent <ushort>\n"
        "  --firmRareEventImpactMinPercent <ushort>\n"
        "  --firmRareEventImpactMaxPercent <ushort>\n\n"
        "Banking/policy:\n"
        "  --interbankDensity <ushort>\n"
        "  --maxInterbankLenderSamplingK <ushort>\n"
        "  --maxInterbankLoanPercent <ushort>\n"
        "  --maxFirmLoanPercent <ushort>\n"
        "  --firmLenderDegree <ushort>\n"
        "  --firmRepayPercent <ushort>\n"
        "  --bankRepayPercent <ushort>\n"
        "  --interventionDelay <uint>\n"
        "  --bailInCoveragePercent <ushort>\n"
        "  --bailOutCoveragePercent <ushort>\n\n"
        "Other:\n"
        "  --help\n"
    );
}

// Print compact start summary before measured run
static void printRunSummary(
    int mpiRank,
    int mpiSize,
    unsigned int runs,
    unsigned int timesteps,
    uint64_t baseSeed,
    unsigned int bankCountTotal,
    unsigned int firmCountTotal,
    unsigned int workerCountTotal
) {
    if (mpiRank != 0) {
        return;
    }

    printf(
        "BailInBailOut start: mpiSize=%d runs=%u timesteps=%u seed=%llu banks=%u firms=%u workers=%u\n",
        mpiSize,
        runs,
        timesteps,
        (unsigned long long)baseSeed,
        bankCountTotal,
        firmCountTotal,
        workerCountTotal
    );
    fflush(
        stdout
    );
}

int main(
    int argc,
    char** argv
) {
    MPI_Init(
        &argc,
        &argv
    );

    int mpiRank = 0;
    int mpiSize = 1;
    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &mpiRank
    );
    MPI_Comm_size(
        MPI_COMM_WORLD,
        &mpiSize
    );

    if (hasFlag(
        argc,
        argv,
        "--help"
    )) {
        printUsage(
            mpiRank
        );
        MPI_Finalize();
        return 0;
    }

    /* Simulation Specifiers */
    unsigned int runs = getUIntArg(
        argc,
        argv,
        "--runs",
        1
    );
    unsigned int timesteps = getUIntArg(
        argc,
        argv,
        "--timesteps",
        8
    );
    uint64_t baseSeed = getU64Arg(
        argc,
        argv,
        "--seed",
        123456789ULL
    );

    /* Structural Parameters */
    unsigned int bankCountTotal = getUIntArg(
        argc,
        argv,
        "--banks",
        4
    );
    unsigned int firmCountTotal = getUIntArg(
        argc,
        argv,
        "--firms",
        8
    );
    unsigned int workerCountTotal = getUIntArg(
        argc,
        argv,
        "--workers",
        16
    );
    unsigned int bankWorkerCount = getUIntArg(
        argc,
        argv,
        "--bankWorkers",
        2
    );

    /* System Parameters */
    unsigned int initialBankLiquidity = getUIntArg(
        argc,
        argv,
        "--initialBankLiquidity",
        5000
    );
    unsigned int initialFirmLiquidity = getUIntArg(
        argc,
        argv,
        "--initialFirmLiquidity",
        100
    );
    unsigned int initialProductionCost = getUIntArg(
        argc,
        argv,
        "--initialProductionCost",
        200
    );
    unsigned int wage = getUIntArg(
        argc,
        argv,
        "--wage",
        500
    );
    unsigned int bankEmployeeWage = getUIntArg(
        argc,
        argv,
        "--bankEmployeeWage",
        100
    );

    /* Multipliers / rates */
    unsigned short wageConsumptionPercent = getUShortArg(
        argc,
        argv,
        "--wageConsumptionPercent",
        50
    );
    unsigned short profitMultiplierMin = getUShortArg(
        argc,
        argv,
        "--profitMultiplierMin",
        100
    );
    unsigned short profitMultiplierMax = getUShortArg(
        argc,
        argv,
        "--profitMultiplierMax",
        100
    );
    unsigned short shockMultiplierMin = getUShortArg(
        argc,
        argv,
        "--shockMultiplierMin",
        0
    );
    unsigned short shockMultiplierMax = getUShortArg(
        argc,
        argv,
        "--shockMultiplierMax",
        0
    );
    unsigned short minInterestRate = getUShortArg(
        argc,
        argv,
        "--minInterestRate",
        7
    );
    unsigned short maxInterestRate = getUShortArg(
        argc,
        argv,
        "--maxInterestRate",
        7
    );

    /* Firm heterogeneity */
    unsigned short firmShockMinPercent = getUShortArg(
        argc,
        argv,
        "--firmShockMinPercent",
        0
    );
    unsigned short firmShockMaxPercent = getUShortArg(
        argc,
        argv,
        "--firmShockMaxPercent",
        0
    );
    unsigned short firmProfitMinPercent = getUShortArg(
        argc,
        argv,
        "--firmProfitMinPercent",
        100
    );
    unsigned short firmProfitMaxPercent = getUShortArg(
        argc,
        argv,
        "--firmProfitMaxPercent",
        100
    );
    unsigned short firmOperatingCostMinPercent = getUShortArg(
        argc,
        argv,
        "--firmOperatingCostMinPercent",
        100
    );
    unsigned short firmOperatingCostMaxPercent = getUShortArg(
        argc,
        argv,
        "--firmOperatingCostMaxPercent",
        100
    );
    unsigned short firmShockProfitMultiplierMinPercent = getUShortArg(
        argc,
        argv,
        "--firmShockProfitMultiplierMinPercent",
        100
    );
    unsigned short firmShockProfitMultiplierMaxPercent = getUShortArg(
        argc,
        argv,
        "--firmShockProfitMultiplierMaxPercent",
        100
    );
    unsigned short firmInitialLiquidityMultiplierMinPercent = getUShortArg(
        argc,
        argv,
        "--firmInitialLiquidityMultiplierMinPercent",
        100
    );
    unsigned short firmInitialLiquidityMultiplierMaxPercent = getUShortArg(
        argc,
        argv,
        "--firmInitialLiquidityMultiplierMaxPercent",
        100
    );
    unsigned short firmRareEventProbabilityPercent = getUShortArg(
        argc,
        argv,
        "--firmRareEventProbabilityPercent",
        0
    );
    unsigned short firmRareEventImpactMinPercent = getUShortArg(
        argc,
        argv,
        "--firmRareEventImpactMinPercent",
        0
    );
    unsigned short firmRareEventImpactMaxPercent = getUShortArg(
        argc,
        argv,
        "--firmRareEventImpactMaxPercent",
        0
    );

    /* Banking Interaction Parameters */
    unsigned short interbankDensity = getUShortArg(
        argc,
        argv,
        "--interbankDensity",
        0
    );
    unsigned short maxInterbankLenderSamplingK = getUShortArg(
        argc,
        argv,
        "--maxInterbankLenderSamplingK",
        1
    );
    unsigned short maxInterbankLoanPercent = getUShortArg(
        argc,
        argv,
        "--maxInterbankLoanPercent",
        0
    );
    unsigned short maxFirmLoanPercent = getUShortArg(
        argc,
        argv,
        "--maxFirmLoanPercent",
        20
    );
    unsigned short firmLenderDegree = getUShortArg(
        argc,
        argv,
        "--firmLenderDegree",
        1
    );
    unsigned short firmRepayPercent = getUShortArg(
        argc,
        argv,
        "--firmRepayPercent",
        25
    );
    unsigned short bankRepayPercent = getUShortArg(
        argc,
        argv,
        "--bankRepayPercent",
        0
    );

    unsigned int interventionDelay = getUIntArg(
        argc,
        argv,
        "--interventionDelay",
        1
    );

    /* Policy */
    unsigned short bailInCoveragePercent = getUShortArg(
        argc,
        argv,
        "--bailInCoveragePercent",
        0
    );
    unsigned short bailOutCoveragePercent = getUShortArg(
        argc,
        argv,
        "--bailOutCoveragePercent",
        0
    );


    // Print run summary before timing begins
    printRunSummary(
        mpiRank,
        mpiSize,
        runs,
        timesteps,
        baseSeed,
        bankCountTotal,
        firmCountTotal,
        workerCountTotal
    );

    MPI_Barrier(
        MPI_COMM_WORLD
    );
    auto t0 = chrono::high_resolution_clock::now();

    BailInBailOut sim;

    // Run benchmark using parsed command line parameters
    sim.runBenchmarkInAndOut(
        /* Simulation Specifiers */
        runs,
        timesteps,
        baseSeed,

        /* Structural Parameters */
        bankCountTotal,
        firmCountTotal,
        workerCountTotal,
        bankWorkerCount,

        /* System Parameters */
        initialBankLiquidity,
        initialFirmLiquidity,
        initialProductionCost,
        wage,
        bankEmployeeWage,

        /* Multipliers (general) */
        wageConsumptionPercent,
        profitMultiplierMin,
        profitMultiplierMax,
        shockMultiplierMin,
        shockMultiplierMax,
        minInterestRate,
        maxInterestRate,
        firmShockMinPercent,
        firmShockMaxPercent,
        firmProfitMinPercent,
        firmProfitMaxPercent,
        firmOperatingCostMinPercent,
        firmOperatingCostMaxPercent,
        firmShockProfitMultiplierMinPercent,
        firmShockProfitMultiplierMaxPercent,
        firmInitialLiquidityMultiplierMinPercent,
        firmInitialLiquidityMultiplierMaxPercent,
        firmRareEventProbabilityPercent,
        firmRareEventImpactMinPercent,
        firmRareEventImpactMaxPercent,

        /* Banking Interaction Parameters */
        interbankDensity,
        maxInterbankLenderSamplingK,
        maxInterbankLoanPercent,
        maxFirmLoanPercent,
        firmLenderDegree,
        firmRepayPercent,
        bankRepayPercent,

        interventionDelay,

        /* Bail-In Parameters */
        bailInCoveragePercent,

        /* Bail-Out Parameters */
        bailOutCoveragePercent
    );

    MPI_Barrier(
        MPI_COMM_WORLD
    );
    auto t1 = chrono::high_resolution_clock::now();
    const long long elapsedMs = chrono::duration_cast<chrono::milliseconds>(
        t1 - t0
    ).count();

    long long maxElapsedMs = 0;
    MPI_Reduce(
        &elapsedMs,
        &maxElapsedMs,
        1,
        MPI_LONG_LONG,
        MPI_MAX,
        0,
        MPI_COMM_WORLD
    );

    if (mpiRank == 0) {
        printf(
            "BailInBailOut done: elapsedMs(END TO END)=%lld\n",
            maxElapsedMs
        );
        fflush(
            stdout
        );
    }

    MPI_Finalize();
    return 0;
}
