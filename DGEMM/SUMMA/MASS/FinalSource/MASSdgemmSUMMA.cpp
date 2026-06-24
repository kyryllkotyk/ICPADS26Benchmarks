/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: DGEMM (SUMMA)                                                   *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "MASSdgemmSUMMA.h"

#ifndef LOGGING
static const bool printOutput = false;
#else
static const bool printOutput = true;
#endif

// Compute greatest common divisor for SUMMA step count
int gcd(
    const int leftValue,
    const int rightValue
)
{
    int currentLeft = leftValue;
    int currentRight = rightValue;

    while (currentRight != 0) {
        const int remainder = currentLeft % currentRight;

        currentLeft = currentRight;
        currentRight = remainder;
    }

    return currentLeft;
}

// Create Places, run SUMMA pipeline, and report timing/checksum results
void runMassDgemmSummaBenchmark(
    char* massArgs[4],
    const int processes,
    const int threads,
    const MassDgemmSummaConfig& config
)
{
    const int stepCount = (
        config.processGridRows / gcd(
            config.processGridRows,
            config.processGridCols
        )
    ) * config.processGridCols;

    const int totalPlaces =
        config.processGridRows * config.processGridCols;

    if (printOutput) {
        std::cerr
            << "Initializing MASS ..."
            << std::endl;
    }

    MASS::init(
        massArgs,
        processes,
        threads
    );

    DGEMMConfig placeConfig{};
    placeConfig.matrixRows = config.matrixRows;
    placeConfig.matrixCols = config.matrixCols;
    placeConfig.sharedDimension = config.sharedDimension;
    placeConfig.runs = config.runs;
    placeConfig.baseSeed = config.baseSeed;

    std::strncpy(
        placeConfig.initType,
        config.initType.c_str(),
        sizeof(
            placeConfig.initType
        ) - 1
    );

    placeConfig.initType[
        sizeof(
            placeConfig.initType
        ) - 1
    ] = '\0';

    DGEMMConfig* const heapConfig = new DGEMMConfig(
        placeConfig
    );

    if (printOutput) {
        std::cerr
            << "Creating "
            << config.processGridRows
            << "x"
            << config.processGridCols
            << " Places"
            << std::endl;
    }

    Places* const dgemm = new Places(
        DGEMM_PLACES_HANDLE,
        "DGEMMPlace",
        0,
        heapConfig,
        sizeof(
            DGEMMConfig
        ),
        2,
        config.processGridRows,
        config.processGridCols
    );

    // Initialize each Place with its tile size and generation settings
    dgemm->callAll(
        "DGEMMPlace::init",
        heapConfig,
        sizeof(
            DGEMMConfig
        )
    );

    mass::IterationConfig pipeline;

    // Each iteration publishes panels, exchanges them, then computes
    pipeline
        .iterations(
            stepCount
        )
        .placeCompute(
            "DGEMMPlace::publishPanels"
        )
        .placeExchangeAll(
            DGEMM_PLACES_HANDLE,
            "DGEMMPlace::recvPanels"
        )
        .placeCompute(
            "DGEMMPlace::accumulate"
        );

    std::cout
        << "--------------DGEMM MASS SUMMA--------------\n"
        << "Matrix: "
        << config.matrixRows
        << "x"
        << config.matrixCols
        << " = "
        << config.matrixRows
        << "x"
        << config.sharedDimension
        << " * "
        << config.sharedDimension
        << "x"
        << config.matrixCols
        << "\n"
        << "Process grid: "
        << config.processGridRows
        << "x"
        << config.processGridCols
        << "  Steps: "
        << stepCount
        << "  Init: "
        << config.initType
        << "  Runs: "
        << config.runs
        << "\n";

    for (int run = 1; run <= config.runs; run++) {
        // Reset local C tiles before each measured run
        dgemm->callAll(
            "DGEMMPlace::resetC"
        );

        const auto timedStart =
            std::chrono::high_resolution_clock::now();

        // Run all SUMMA steps as one MASS iteration pipeline
        MASS::runIterations(
            DGEMM_PLACES_HANDLE,
            pipeline
        );

        const auto timedEnd =
            std::chrono::high_resolution_clock::now();

        double checksum = 0.0;

        try {
            // Collect 1 checksum value from each Place
            double* const localSums = reinterpret_cast<double*>(
                dgemm->callAll(
                    "DGEMMPlace::getChecksum",
                    static_cast<void**>(
                        nullptr
                    ),
                    0,
                    static_cast<int>(
                        sizeof(
                            double
                        )
                    )
                )
            );

            if (localSums != nullptr) {
                for (
                    int placeIndex = 0;
                    placeIndex < totalPlaces;
                    placeIndex++
                ) {
                    checksum += localSums[placeIndex];
                }
            }
        }
        catch (...) {
            std::cerr
                << "Warning: checksum collection failed"
                << std::endl;
        }

        const double elapsed = std::chrono::duration<double>(
            timedEnd - timedStart
        ).count();

        dgemm_report::printRunResults(
            run,
            checksum,
            elapsed,
            config.matrixRows,
            config.matrixCols,
            config.sharedDimension
        );
    }

    MASS::finish();
}
