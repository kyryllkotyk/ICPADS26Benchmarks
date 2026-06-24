#include "pm2dgemmSUMMA.h"

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <string>

using namespace std;

static void printUsage(
    const char* const programName
) {
    cout
        << "Usage:\n"
        << "  "
        << programName
        << " <algorithm>"
        << " [--M rowsOfAAndC]"
        << " [--N colsOfBAndC]"
        << " [--K colsOfA_rowsOfB]"
        << " [--Px processGridRows]"
        << " [--Py processGridCols]"
        << " [--runs numberOfRuns]"
        << " [--seed baseSeed]\n\n"
        << "Algorithms:\n"
        << "  summa | regular | s   -> runBenchmark\n";
}

static bool isSummaName(
    const string& algorithmName
) {
    return algorithmName == "summa"
        || algorithmName == "regular"
        || algorithmName == "s";
}

int main(
    int argc,
    char* argv[]
) {
    MPI_Init(
        &argc,
        &argv
    );

    int rank = 0;

    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &rank
    );

    if (argc < 2) {
        if (rank == 0) {
            cerr << "Missing algorithm name.\n";

            printUsage(
                argv[0]
            );
        }

        MPI_Finalize();
        return 1;
    }

    const string algorithmName = argv[1];

    if (algorithmName == "--help" || algorithmName == "-h") {
        if (rank == 0) {
            printUsage(
                argv[0]
            );
        }

        MPI_Finalize();
        return 0;
    }

    if (!isSummaName(
        algorithmName
    )) {
        if (rank == 0) {
            cerr
                << "Unknown algorithm: "
                << algorithmName
                << '\n';

            printUsage(
                argv[0]
            );
        }

        MPI_Finalize();
        return 1;
    }

    int rowsOfAAndC = 2048;
    int colsOfBAndC = 2048;
    int sharedDimension = 2048;
    int processGridRows = 1;
    int processGridCols = 1;
    int runs = 1;
    uint64_t seed = 1;

    bool parseError = false;

    for (int argumentIndex = 2; argumentIndex < argc; argumentIndex++) {
        const string argument = argv[argumentIndex];

        // Fetch value that belongs to current command line option
        auto requireValue = [
            &
        ](
            const string& optionName
            ) -> const char* {
                if (argumentIndex + 1 >= argc) {
                    if (rank == 0) {
                        cerr
                            << "Missing value for "
                            << optionName
                            << ".\n";
                    }

                    parseError = true;
                    return nullptr;
                }

                argumentIndex++;
                return argv[argumentIndex];
            };

            if (argument == "--M") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    rowsOfAAndC = atoi(
                        value
                    );
                }
            }
            else if (argument == "--N") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    colsOfBAndC = atoi(
                        value
                    );
                }
            }
            else if (argument == "--K") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    sharedDimension = atoi(
                        value
                    );
                }
            }
            else if (argument == "--Px") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    processGridRows = atoi(
                        value
                    );
                }
            }
            else if (argument == "--Py") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    processGridCols = atoi(
                        value
                    );
                }
            }
            else if (argument == "--runs") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    runs = atoi(
                        value
                    );
                }
            }
            else if (argument == "--seed") {
                const char* value = requireValue(
                    argument
                );

                if (value != nullptr) {
                    seed = static_cast<uint64_t>(
                        strtoull(
                            value,
                            nullptr,
                            10
                        )
                        );
                }
            }
            else if (argument == "--help" || argument == "-h") {
                if (rank == 0) {
                    printUsage(
                        argv[0]
                    );
                }

                MPI_Finalize();
                return 0;
            }
            else {
                if (rank == 0) {
                    cerr
                        << "Unknown argument: "
                        << argument
                        << '\n';

                    printUsage(
                        argv[0]
                    );
                }

                parseError = true;
            }
    }

    if (parseError) {
        MPI_Finalize();
        return 1;
    }

    PM2_DGEMM benchmark;

    benchmark.runBenchmark(
        static_cast<short>(
            runs
            ),
        rowsOfAAndC,
        sharedDimension,
        colsOfBAndC,
        seed,
        processGridRows,
        processGridCols
    );

    MPI_Finalize();
    return 0;
}