/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Heat3D                                                          *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "MASS.h"
#include "IterationConfig.h"
#include "Heat3DPlace.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>

#define PLACES_HANDLE 1

namespace
{
    // Store MASS launcher options and Heat3D benchmark parameters
    struct CliConfig
    {
        int gridX = 128;
        int gridY = 128;
        int gridZ = 128;
        int decompX = 1;
        int decompY = 1;
        int decompZ = 1;
        int timesteps = 100;
        int totalRuns = 1;
        uint64_t seed = 1;
        double initMin = 0.0;
        double initMax = 1.0;
        double alpha = 0.5;
        double beta = 0.1;
        bool debug = false;
    };

    // Parse Heat3D options after the required MASS runtime arguments
    void parseArgs(
        int argc,
        char* argv[],
        CliConfig& config
    ) {
        // First 6 arguments belong to MASS startup
        for (int argumentIndex = 7; argumentIndex < argc; argumentIndex++) {
            const char* key = argv[argumentIndex];

            if (argumentIndex + 1 >= argc) {
                break;
            }

            // Options are parsed as key value pairs
            const char* value = argv[argumentIndex + 1];

            if (!std::strcmp(
                key,
                "--gridX"
            )) {
                config.gridX = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--gridY"
            )) {
                config.gridY = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--gridZ"
            )) {
                config.gridZ = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--decompX"
            )) {
                config.decompX = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--decompY"
            )) {
                config.decompY = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--decompZ"
            )) {
                config.decompZ = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--timesteps"
            )) {
                config.timesteps = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--runs"
            )) {
                config.totalRuns = std::atoi(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--seed"
            )) {
                config.seed = std::strtoull(
                    value,
                    nullptr,
                    10
                );
            }
            else if (!std::strcmp(
                key,
                "--initMin"
            )) {
                config.initMin = std::atof(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--initMax"
            )) {
                config.initMax = std::atof(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--alpha"
            )) {
                config.alpha = std::atof(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--beta"
            )) {
                config.beta = std::atof(
                    value
                );
            }
            else if (!std::strcmp(
                key,
                "--debug"
            )) {
                config.debug = std::atoi(
                    value
                ) != 0;
            }
            else {
                continue;
            }

            argumentIndex++;
        }
    }

    // Print MASS launcher syntax and benchmark specific options
    void usage(
        const char* executable
    ) {
        std::cerr
            << "Usage: "
            << executable
            << " user pass machinefile port nProc nThr [options]\n"
            << "Heat3D options:\n"
            << "  --gridX/Y/Z N     global grid\n"
            << "  --decompX/Y/Z N   decomposition\n"
            << "  --timesteps N     iteration count\n"
            << "  --runs N          benchmark repeats\n"
            << "  --seed N          init seed\n"
            << "  --initMin/Max F   init range\n"
            << "  --alpha/beta F    stencil weights\n"
            << "  --debug 0|1       per-run output table\n";
    }
}

int main(
    int argc,
    char* argv[]
) {
    if (argc < 7) {
        usage(
            argv[0]
        );

        return 1;
    }

    // MASS expects only login, password, machine file, and port here
    char* massArgs[4] = {
        argv[1],
        argv[2],
        argv[3],
        argv[4]
    };

    const int processes = std::atoi(
        argv[5]
    );

    const int threads = std::atoi(
        argv[6]
    );

    CliConfig config;

    parseArgs(
        argc,
        argv,
        config
    );

    // One MASS Place is created for each 3D decomposition chunk
    const int totalChunks = config.decompX * config.decompY * config.decompZ;

    // End-to-end timing includes MASS setup, all runs, and cleanup
    const auto endToEndStart = std::chrono::high_resolution_clock::now();

    MASS::init(
        massArgs,
        processes,
        threads
    );

    // Package benchmark settings passed to every Place
    Heat3DConfig placeConfig{};
    placeConfig.gridX = config.gridX;
    placeConfig.gridY = config.gridY;
    placeConfig.gridZ = config.gridZ;
    placeConfig.seed = config.seed;
    placeConfig.initMin = config.initMin;
    placeConfig.initMax = config.initMax;
    placeConfig.alpha = config.alpha;
    placeConfig.beta = config.beta;

    // Pass dimensions as z,y,x so x varies fastest in Place index order
    Places* heat = new Places(
        PLACES_HANDLE,
        "Heat3DPlace",
        &placeConfig,
        sizeof(
            Heat3DConfig
            ),
        3,
        config.decompZ,
        config.decompY,
        config.decompX
    );

    // Initialize all Places before building the iteration pipeline
    heat->callAll(
        "Heat3DPlace::init",
        &placeConfig,
        sizeof(
            Heat3DConfig
            )
    );

    // Build the timestep pipeline once and reuse it for every run
    mass::IterationConfig pipeline;

    pipeline.iterations(
        config.timesteps
    ).placeCompute(
        "Heat3DPlace::packFaces"
    ).placeExchangeAll(
        PLACES_HANDLE,
        "Heat3DPlace::recvHalo"
    ).placeCompute(
        "Heat3DPlace::computeStep"
    );

    if (config.debug) {
        std::cout
            << "--------------CONFIG-----------------\n"
            << "  [MASS chunk-per-Place + IterationConfig]\n"
            << "Chunks: "
            << totalChunks
            << " ("
            << config.decompX
            << "x"
            << config.decompY
            << "x"
            << config.decompZ
            << ")\n"
            << "Global Grid: "
            << config.gridX
            << " x "
            << config.gridY
            << " x "
            << config.gridZ
            << "\n"
            << "Timesteps: "
            << config.timesteps
            << "\n"
            << "Runs: "
            << config.totalRuns
            << "\n"
            << "Alpha: "
            << config.alpha
            << " | Beta: "
            << config.beta
            << "\n"
            << "Seed: "
            << config.seed
            << "\n"
            << "Init range: ["
            << config.initMin
            << ", "
            << config.initMax
            << "]\n"
            << "-------------------------------------\n"
            << "Run | Simulation(ms) | Checksum\n"
            << "-----------------------------------------------\n";
    }

    // Run the same timestep pipeline for each measured repeat
    for (int run = 0; run < config.totalRuns; run++) {
        if (run > 0) {
            // Reset Place grids while preserving allocated communication state
            heat->callAll(
                "Heat3DPlace::reInit"
            );
        }

        const auto wallStart = std::chrono::high_resolution_clock::now();

        MASS::runIterations(
            PLACES_HANDLE,
            pipeline
        );

        const auto wallEnd = std::chrono::high_resolution_clock::now();

        // Gather one local checksum from each Place
        double* partials = reinterpret_cast<double*>(
            heat->callAll(
                "Heat3DPlace::collectChecksum",
                static_cast<void**>(
                    nullptr
                    ),
                0,
                sizeof(
                    double
                    )
            )
            );

        double globalChecksum = 0.0;

        // Sum Place checksums into one global validation value
        for (int chunk = 0; chunk < totalChunks; chunk++) {
            globalChecksum += partials[chunk];
        }

        const long simulationMs = std::chrono::duration_cast<
            std::chrono::milliseconds
        >(
            wallEnd - wallStart
        ).count();

        if (config.debug) {
            std::cout
                << std::setw(
                    3
                )
                << run
                << " | "
                << std::setw(
                    14
                )
                << simulationMs
                << " | "
                << std::setprecision(
                    15
                )
                << globalChecksum
                << "\n";
        }
    }

    // Shut down MASS after all Place work is complete
    MASS::finish();

    const auto endToEndEnd = std::chrono::high_resolution_clock::now();

    const long endToEndMs = std::chrono::duration_cast<
        std::chrono::milliseconds
    >(
        endToEndEnd - endToEndStart
    ).count();

    std::cout
        << "elapsedMs(END TO END)="
        << endToEndMs
        << "\n";

    return 0;
}
