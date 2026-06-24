/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Heat3D                                                          *
 * Library: PM2/MPI                                                           *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "pm2heat3d.h"

 // Parse --option value pairs into a simple lookup table
static unordered_map<string, string> parseArgs(
    const int argc,
    char** argv
) {
    unordered_map<string, string> args;

    for (int argumentIndex = 1; argumentIndex < argc; argumentIndex++) {
        const string key = argv[argumentIndex];

        if (key.rfind(
            "--",
            0
        ) == 0) {
            if (argumentIndex + 1 >= argc) {
                cerr
                    << "Missing value for "
                    << key
                    << "\n";

                exit(
                    1
                );
            }

            args[key.substr(
                2
            )] = argv[++argumentIndex];
        }
        else {
            cerr
                << "Invalid argument: "
                << key
                << "\n";

            exit(
                1
            );
        }
    }

    return args;
}

template <typename T>
static T getArg(
    const unordered_map<string, string>& args,
    const string& key,
    const T& defaultValue
) {
    const auto iterator = args.find(
        key
    );

    if (iterator == args.end()) {
        return defaultValue;
    }

    // Convert stored string to requested argument type
    if constexpr (is_same<T, int>::value) {
        return stoi(
            iterator->second
        );
    }
    else if constexpr (is_same<T, double>::value) {
        return stod(
            iterator->second
        );
    }
    else if constexpr (is_same<T, float>::value) {
        return stof(
            iterator->second
        );
    }
    else if constexpr (is_same<T, uint64_t>::value) {
        return stoull(
            iterator->second
        );
    }
    else {
        return iterator->second;
    }
}

int main(
    int argc,
    char** argv
) {
    // Initialize MPI before benchmark setup or communication
    MPI_Init(
        &argc,
        &argv
    );

    // Read command line options once before extracting defaults
    const unordered_map<string, string> args = parseArgs(
        argc,
        argv
    );

    // Parse benchmark parameters with same defaults as original driver
    const bool debug = getArg<int>(
        args,
        "debug",
        0
    );

    const short timesteps = static_cast<short>(
        getArg<int>(
            args,
            "timesteps",
            100
        )
        );

    const short totalRuns = static_cast<short>(
        getArg<int>(
            args,
            "runs",
            1
        )
        );

    const short screenshotEvery = static_cast<short>(
        getArg<int>(
            args,
            "screenshotEvery",
            0
        )
        );

    const uint64_t seed = getArg<uint64_t>(
        args,
        "seed",
        1
    );

    const short globalGridX = static_cast<short>(
        getArg<int>(
            args,
            "gridX",
            128
        )
        );

    const short globalGridY = static_cast<short>(
        getArg<int>(
            args,
            "gridY",
            128
        )
        );

    const short globalGridZ = static_cast<short>(
        getArg<int>(
            args,
            "gridZ",
            128
        )
        );

    const double initMin = getArg<double>(
        args,
        "initMin",
        0.0
    );

    const double initMax = getArg<double>(
        args,
        "initMax",
        1.0
    );

    const double alpha = getArg<double>(
        args,
        "alpha",
        0.5
    );

    const double beta = getArg<double>(
        args,
        "beta",
        0.1
    );

    const short decompX = static_cast<short>(
        getArg<int>(
            args,
            "decompX",
            1
        )
        );

    const short decompY = static_cast<short>(
        getArg<int>(
            args,
            "decompY",
            1
        )
        );

    const short decompZ = static_cast<short>(
        getArg<int>(
            args,
            "decompZ",
            1
        )
        );

    PM2_Heat3D solver;

    // Hand parsed benchmark configuration to implementation file
    solver.runBenchmark(
        debug,
        timesteps,
        totalRuns,
        screenshotEvery,
        seed,
        globalGridX,
        globalGridY,
        globalGridZ,
        initMin,
        initMax,
        alpha,
        beta,
        decompX,
        decompY,
        decompZ
    );

    // Finalize MPI after benchmark output is complete
    MPI_Finalize();
    return 0;
}