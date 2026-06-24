/****************************************************************************
 *                                                                          *
 * ICPADS26 Benchmarks Collection                                           *
 *                                                                          *
 * Benchmark: DGEMM (Cannon's Algorithm)                                    *
 * Library: MASS                                                            *
 *                                                                          *
 * Author: Ahmed Bera Pay                                                   *
 * Faculty Advisor: Munehiro Fukuda                                         *
 * Code Finalization: Kyryll Kotyk                                          *
 *                                                                          *
 ****************************************************************************/

#include "DGEMMCannonPlace.h"
#include "IterationConfig.h"
#include "MASS.h"

#include "../../include/dgemm_run_report.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#ifndef LOGGING
#define PRINT_OUTPUT false
#else
#define PRINT_OUTPUT true
#endif

#define PLACES_HANDLE 1
#define MASS_DGEMM_DEFAULT_DIMENSION 2048
#define MASS_DGEMM_DEFAULT_GRID_SIZE 1
#define MASS_DGEMM_DEFAULT_RUNS 1
#define MASS_DGEMM_DEFAULT_SEED 1
#define MASS_ARGS_COUNT 7

namespace
{

struct CliConfig
{
    int M = MASS_DGEMM_DEFAULT_DIMENSION;
    int N = MASS_DGEMM_DEFAULT_DIMENSION;
    int K = MASS_DGEMM_DEFAULT_DIMENSION;
    int P = MASS_DGEMM_DEFAULT_GRID_SIZE;
    int runs = MASS_DGEMM_DEFAULT_RUNS;
    uint64_t baseSeed = MASS_DGEMM_DEFAULT_SEED;
    std::string init_type = "hash";
};

// Parse benchmark options after required MASS launch arguments
// Unknown options are ignored so shared runners can pass extra flags
void parseArgs(
    const int argc,
    char *argv[],
    CliConfig &config
) {
    for (
        int argumentIndex = MASS_ARGS_COUNT;
        argumentIndex < argc;
        argumentIndex++
    ) {
        const char *key = argv[argumentIndex];

        // Every supported option expects a value after it
        if (argumentIndex + 1 >= argc) {
            break;
        }

        const char *value = argv[argumentIndex + 1];

        if (!std::strcmp(
            key,
            "--M"
        )) {
            config.M = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--N"
        )) {
            config.N = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--K"
        )) {
            config.K = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--P"
        )) {
            config.P = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--Px"
        )) {
            config.P = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--Py"
        )) {
            config.P = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--runs"
        )) {
            config.runs = std::atoi(
                value
            );
        }
        else if (!std::strcmp(
            key,
            "--seed"
        )) {
            config.baseSeed = std::strtoull(
                value,
                nullptr,
                10
            );
        }
        else if (!std::strcmp(
            key,
            "--init-type"
        )) {
            config.init_type = value;
        }
        else if (!std::strcmp(
            key,
            "--block-size"
        )) {
            // Block size is accepted for runner compatibility
            // Cannon tile sizes are determined by matrix and grid size
        }
        else {
            // Unknown options are ignored to preserve runner compatibility
        }

        argumentIndex++;
    }
}

// Print MASS launcher and benchmark option format
void usage(
    const char *executableName
) {
    std::cerr
        << "Usage: "
        << executableName
        << " user pass machinefile port nProc nThr [options]\n\n"
        << "Cannon requires square Place grid P×P (--Px/--Py).\n"
        << "nProc = SSH processes (nodes); nThr = threads per process (ppn).\n\n"
        << "DGEMM options:\n"
        << "  --M N            rows of A/C (default 2048)\n"
        << "  --N N            cols of B/C (default 2048)\n"
        << "  --K N            cols of A / rows of B (default 2048)\n"
        << "  --P N            grid dimension\n"
        << "  --Px N --Py N    grid dimension aliases\n"
        << "  --runs N         benchmark repetitions (default 1)\n"
        << "  --seed N         RNG base seed (default 1)\n"
        << "  --init-type T    hash | deterministic | identity | ones\n";
}

}

int main(
    int argc,
    char *argv[]
) {
    if (argc < MASS_ARGS_COUNT) {
        usage(
            argv[0]
        );

        return 1;
    }

    // MASS expects only login, password, machinefile, and port here
    char *massArgs[4] = {
        argv[1],
        argv[2],
        argv[3],
        argv[4]
    };

    // Remaining required launch arguments choose process and thread counts
    const int processCount = std::atoi(
        argv[5]
    );

    const int threadCount = std::atoi(
        argv[6]
    );

    // Benchmark options start after MASS launch arguments
    CliConfig config;

    parseArgs(
        argc,
        argv,
        config
    );

    // Invalid grid size falls back to 1 Place
    if (config.P <= 0) {
        config.P = MASS_DGEMM_DEFAULT_GRID_SIZE;
    }

    const int gridSize = config.P;
    const int totalPlaces = gridSize * gridSize;

    // Cannon Place grid is P×P; Hermes runners pass P for nodes×ppn localities.
    if (totalPlaces < 1) {
        std::cerr
            << "Error: invalid Cannon grid P="
            << config.P
            << std::endl;

        return 1;
    }

    if (PRINT_OUTPUT) {
        std::cerr
            << "Initialising MASS ..."
            << std::endl;
    }

    MASS::init(
        massArgs,
        processCount,
        threadCount
    );

    // Pack configuration so each Place can initialize local ownership
    DGEMMCannonConfig placeConfig{};

    placeConfig.M = config.M;
    placeConfig.N = config.N;
    placeConfig.K = config.K;
    placeConfig.runs = config.runs;
    placeConfig.baseSeed = config.baseSeed;

    std::strncpy(
        placeConfig.init_type,
        config.init_type.c_str(),
        sizeof(
            placeConfig.init_type
        ) - 1
    );

    placeConfig.init_type[
        sizeof(
            placeConfig.init_type
        ) - 1
    ] = '\0';

    // MASS keeps pointer data during Place construction and init call
    DGEMMCannonConfig *heapConfig = new DGEMMCannonConfig(
        placeConfig
    );

    if (PRINT_OUTPUT) {
        std::cerr
            << "Creating "
            << gridSize
            << "x"
            << gridSize
            << " Places (Cannon)"
            << std::endl;
    }

    // Create 1 MASS Place for each Cannon tile
    Places *dgemm = new Places(
        PLACES_HANDLE,
        "DGEMMCannonPlace",
        0,
        heapConfig,
        sizeof(
            DGEMMCannonConfig
        ),
        2,
        gridSize,
        gridSize
    );

    // Initialize local tile ranges, panels, and neighbor metadata
    dgemm->callAll(
        "DGEMMCannonPlace::init",
        heapConfig,
        sizeof(
            DGEMMCannonConfig
        )
    );

    // Pipeline models repeated Cannon shift and multiply steps
    mass::IterationConfig pipeline;

    if (gridSize > 1) {
        // Steps 1 through P - 1 shift panels then multiply
        pipeline.iterations(
            gridSize - 1
        ).placeCompute(
            "DGEMMCannonPlace::publishShift"
        ).placeExchangeAll(
            PLACES_HANDLE,
            "DGEMMCannonPlace::recvShift"
        ).placeCompute(
            "DGEMMCannonPlace::shiftFinish"
        ).placeCompute(
            "DGEMMCannonPlace::multiplyOnly"
        );
    }

    // Print config once before timed benchmark runs
    std::cout
        << "--------------DGEMM MASS "
        << "(Cannon, placeExchangeAll ring shifts)----\n"
        << "Matrix: "
        << config.M
        << "x"
        << config.N
        << " = "
        << config.M
        << "x"
        << config.K
        << " * "
        << config.K
        << "x"
        << config.N
        << "\n"
        << "Grid: "
        << gridSize
        << "x"
        << gridSize
        << " (Cannon, square)"
        << "  Steps: "
        << gridSize
        << "  Init: "
        << config.init_type
        << "  Runs: "
        << config.runs
        << "\n";

    for (int run = 1; run <= config.runs; run++) {
        // Reset C outside timed region before each run
        dgemm->callAll(
            "DGEMMCannonPlace::resetC"
        );

        // Timed phase includes initial fill, multiply, shifts, and updates
        const auto startTime =
            std::chrono::high_resolution_clock::now();

        // Fill and step 0 multiply happen before shift pipeline starts
        dgemm->callAll(
            "DGEMMCannonPlace::initialFill"
        );

        dgemm->callAll(
            "DGEMMCannonPlace::multiplyStep0"
        );

        if (gridSize > 1) {
            // Remaining steps shift panels then multiply each new pair
            MASS::runIterations(
                PLACES_HANDLE,
                pipeline
            );
        }

        const auto endTime =
            std::chrono::high_resolution_clock::now();

        // Checksum is collected after timed region
        double checksum = 0.0;

        try {
            double *sums = reinterpret_cast<double *>(
                dgemm->callAll(
                    "DGEMMCannonPlace::getChecksum",
                    static_cast<void **>(
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

            if (sums) {
                // Sum 1 checksum contribution from each Place
                for (
                    int placeIndex = 0;
                    placeIndex < totalPlaces;
                    placeIndex++
                ) {
                    checksum += sums[placeIndex];
                }
            }
        }
        catch (...) {
            std::cerr
                << "Warning: checksum collection failed"
                << std::endl;
        }

        // Convert measured wall time for report output
        const double elapsedSeconds =
            std::chrono::duration<double>(
                endTime - startTime
            ).count();

        dgemm_report::printRunResults(
            run,
            checksum,
            elapsedSeconds,
            config.M,
            config.N,
            config.K
        );
    }

    // Shut down MASS runtime after all benchmark runs finish
    MASS::finish();

    return 0;
}
