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

// Print required MASS arguments and optional benchmark arguments
static void printUsage(
    const char* executableName
)
{
    std::cerr
        << "Usage: "
        << executableName
        << " user pass machinefile port nProc nThr [options]\n\n"
        << "MASS arguments:\n"
        << "  user         SSH username\n"
        << "  pass         SSH password, or '' for key based auth\n"
        << "  machinefile  file with cluster hostnames\n"
        << "  port         MASS communication port\n"
        << "  nProc        SSH processes (typically cluster nodes)\n"
        << "  nThr         threads per process (Hermes: ppn)\n"
        << "  Px * Py      process-grid Places (may exceed nProc; see Heat3D model)\n\n"
        << "DGEMM options:\n"
        << "  --M N            rows of A / C            (default 2048)\n"
        << "  --N N            cols of B / C            (default 2048)\n"
        << "  --K N            cols of A / rows of B    (default 2048)\n"
        << "  --Px N           process grid rows        (default 1)\n"
        << "  --Py N           process grid cols        (default 1)\n"
        << "  --runs N         benchmark repetitions    (default 1)\n"
        << "  --seed N         RNG base seed            (default 1)\n"
        << "  --init-type T    hash | deterministic | identity | ones\n";
}

// Parse optional DGEMM arguments after required MASS launch arguments
static void parseArgs(
    const int argc,
    char* argv[],
    MassDgemmSummaConfig& config
)
{
    for (int argumentIndex = 7; argumentIndex < argc; argumentIndex++) {
        const char* const optionName = argv[argumentIndex];

        if (argumentIndex + 1 >= argc) {
            break;
        }

        const char* const optionValue = argv[argumentIndex + 1];

        if (std::strcmp(
            optionName,
            "--M"
        ) == 0) {
            config.matrixRows = std::atoi(
                optionValue
            );
        }
        else if (std::strcmp(
            optionName,
            "--N"
        ) == 0) {
            config.matrixCols = std::atoi(
                optionValue
            );
        }
        else if (std::strcmp(
            optionName,
            "--K"
        ) == 0) {
            config.sharedDimension = std::atoi(
                optionValue
            );
        }
        else if (std::strcmp(
            optionName,
            "--Px"
        ) == 0) {
            config.processGridRows = std::atoi(
                optionValue
            );
        }
        else if (std::strcmp(
            optionName,
            "--Py"
        ) == 0) {
            config.processGridCols = std::atoi(
                optionValue
            );
        }
        else if (std::strcmp(
            optionName,
            "--runs"
        ) == 0) {
            config.runs = std::atoi(
                optionValue
            );
        }
        else if (std::strcmp(
            optionName,
            "--seed"
        ) == 0) {
            config.baseSeed = static_cast<uint64_t>(
                std::strtoull(
                    optionValue,
                    nullptr,
                    10
                )
            );
        }
        else if (std::strcmp(
            optionName,
            "--init-type"
        ) == 0) {
            config.initType = optionValue;
        }
        else if (std::strcmp(
            optionName,
            "--block-size"
        ) == 0) {
            // Accepted for runner compatibility
            // SUMMA Place implementation does not use this value
        }
        else {
            // Unknown options are ignored to match existing runner behavior
        }

        argumentIndex++;
    }
}

int main(
    int argc,
    char* argv[]
)
{
    if (argc < 7) {
        printUsage(
            argv[0]
        );

        return 1;
    }

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

    MassDgemmSummaConfig config;

    parseArgs(
        argc,
        argv,
        config
    );

    const int gridPlaces =
        config.processGridRows * config.processGridCols;

    if (gridPlaces < 1) {
        std::cerr
            << "Error: invalid process grid Px="
            << config.processGridRows
            << " Py="
            << config.processGridCols
            << std::endl;

        return 1;
    }

    runMassDgemmSummaBenchmark(
        massArgs,
        processes,
        threads,
        config
    );

    return 0;
}
