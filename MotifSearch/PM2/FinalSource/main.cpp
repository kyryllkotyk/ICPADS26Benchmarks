#include "pm2motifsearch.h"

#include <chrono>
#include <cctype>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

// Store 1 benchmark row from CLI defaults or batch CSV
struct BatchRun
{
    int vertices = 0;
    double probability = 0.0;
    unsigned seed = 1;
    string motifName = "triangle";
    string inputFile = "";
    bool debug = false;
    int runs = 3;
};

// Remove surrounding whitespace from CSV fields
static string trim(
    const string& text
) {
    const size_t start = text.find_first_not_of(
        " \t\r\n"
    );

    if (start == string::npos) {
        return "";
    }

    const size_t end = text.find_last_not_of(
        " \t\r\n"
    );

    return text.substr(
        start,
        end - start + 1
    );
}

// Accept simple true values from batch CSV
static bool parseBoolValue(
    const string& text
) {
    string lowered = text;

    for (char& character : lowered) {
        character = static_cast<char>(
            tolower(
                static_cast<unsigned char>(
                    character
                )
            )
        );
    }

    return lowered == "1" || lowered == "true" || lowered == "yes";
}

// Split 1 CSV row while preserving quoted commas
static vector<string> splitCsvLine(
    const string& line
) {
    vector<string> fields;
    string current;
    bool inQuotes = false;

    for (const char character : line) {
        if (character == '"') {
            inQuotes = !inQuotes;
        }
        else if (character == ',' && !inQuotes) {
            fields.push_back(
                trim(
                    current
                )
            );
            current.clear();
        }
        else {
            current += character;
        }
    }

    fields.push_back(
        trim(
            current
        )
    );

    return fields;
}

// Make input path safe for generated log file names
static string sanitizeFileStem(
    const string& path
) {
    const size_t slashPosition = path.find_last_of(
        "/\\"
    );
    string name = (
        slashPosition == string::npos
            ? path
            : path.substr(
                slashPosition + 1
            )
    );

    const size_t dotPosition = name.find_last_of(
        '.'
    );

    if (dotPosition != string::npos) {
        name = name.substr(
            0,
            dotPosition
        );
    }

    for (char& character : name) {
        const bool isSafeCharacter = isalnum(
            static_cast<unsigned char>(
                character
            )
        ) || character == '_' || character == '-';

        if (!isSafeCharacter) {
            character = '_';
        }
    }

    return name;
}

// Build result log name from benchmark parameters
static string deriveLogfileName(
    const string& implementation,
    const string& motif,
    const int vertices,
    const double edgeProbability,
    const unsigned seed,
    const int nodes,
    const string& inputFile
) {
    ostringstream output;
    output
        << implementation
        << "_"
        << motif
        << "_";

    if (!inputFile.empty()) {
        output
            << "F"
            << sanitizeFileStem(
                inputFile
            );
    }
    else {
        output
            << "V"
            << vertices
            << "_p"
            << edgeProbability
            << "_s"
            << seed;
    }

    output
        << "_N"
        << nodes;

    return output.str();
}

// Create timestamp suffix for default CSV name
static string makeTimestampString() {
    using namespace chrono;

    const auto now = system_clock::now();
    const auto millisecondsPart = duration_cast<milliseconds>(
        now.time_since_epoch()
    ) % 1000;

    const time_t nowTime = system_clock::to_time_t(
        now
    );
    tm localTime{};

#if defined(_WIN32)
    localtime_s(
        &localTime,
        &nowTime
    );
#else
    localtime_r(
        &nowTime,
        &localTime
    );
#endif

    ostringstream output;
    output
        << put_time(
            &localTime,
            "%Y%m%d%H%M%S"
        )
        << setw(
            3
        )
        << setfill(
            '0'
        )
        << millisecondsPart.count();

    return output.str();
}

// Build default output CSV file name
static string makeOutputCsvName() {
    return "pm2_motifsearch_" + makeTimestampString() + ".csv";
}

// Write common result CSV columns
static void writeCsvHeader(
    ofstream& output
) {
    output
        << "implementation,motif,vertices,edge_prob,seed,total_processes,"
        << "nodes_used,processes_per_node,edges,raw_embeddings,"
        << "automorphisms,instance_count,search_time_ms,"
        << "total_wall_time_ms,status,logfile\n";
}

// Write 1 averaged benchmark result row
static void writeCsvRow(
    ofstream& output,
    const string& motifName,
    const int requestedVertices,
    const double probability,
    const unsigned seed,
    const int totalProcesses,
    const int nodesUsed,
    const int processesPerNode,
    const string& inputFile,
    const vector<PM2_MotifSearch::MotifSearchRunResult>& results
) {
    double totalSearchTimeMs = 0.0;

    for (const auto& result : results) {
        totalSearchTimeMs += result.searchTimeMs;
    }

    const double averageSearchTimeMs = totalSearchTimeMs / results.size();
    const auto& lastResult = results.back();
    const int outputVertices = (
        lastResult.vertices > 0
            ? lastResult.vertices
            : requestedVertices
    );

    output
        << "pm2"
        << ","
        << motifName
        << ","
        << outputVertices
        << ","
        << probability
        << ","
        << seed
        << ","
        << totalProcesses
        << ","
        << nodesUsed
        << ","
        << processesPerNode
        << ","
        << lastResult.edges
        << ","
        << lastResult.rawEmbeddings
        << ","
        << lastResult.automorphisms
        << ","
        << lastResult.instanceCount
        << ","
        << fixed
        << setprecision(
            3
        )
        << averageSearchTimeMs
        << ","
        << averageSearchTimeMs
        << ","
        << "OK"
        << ","
        << deriveLogfileName(
            "pm2",
            motifName,
            outputVertices,
            probability,
            seed,
            totalProcesses,
            inputFile
        )
        << "\n";
}

// Convert parsed CSV fields into 1 batch run config
static BatchRun parseBatchRun(
    const vector<string>& fields
) {
    BatchRun run;
    run.motifName = fields[0];
    run.vertices = (
        fields[1].empty()
            ? 0
            : stoi(
                fields[1]
            )
    );
    run.probability = (
        fields[2].empty()
            ? 0.0
            : stod(
                fields[2]
            )
    );
    run.seed = (
        fields[3].empty()
            ? 1u
            : static_cast<unsigned>(
                stoul(
                    fields[3]
                )
            )
    );
    run.inputFile = fields[4];
    run.debug = parseBoolValue(
        fields[5]
    );

    if (fields.size() >= 7 && !fields[6].empty()) {
        run.runs = stoi(
            fields[6]
        );
    }

    return run;
}

// Load all valid benchmark rows from batch CSV
static vector<BatchRun> loadBatchFile(
    const string& batchFilePath
) {
    ifstream input(
        batchFilePath
    );

    if (!input.is_open()) {
        throw runtime_error(
            "Could not open batch file: " + batchFilePath
        );
    }

    vector<BatchRun> runs;
    string line;
    bool firstNonEmptyLine = true;

    // Parse optional CSV header and each benchmark row
    while (getline(
        input,
        line
    )) {
        line = trim(
            line
        );

        if (line.empty() || line[0] == '#') {
            continue;
        }

        vector<string> fields = splitCsvLine(
            line
        );

        if (firstNonEmptyLine) {
            firstNonEmptyLine = false;

            if (!fields.empty()) {
                string firstField = fields[0];

                for (char& character : firstField) {
                    character = static_cast<char>(
                        tolower(
                            static_cast<unsigned char>(
                                character
                            )
                        )
                    );
                }

                if (firstField == "motif" || firstField == "vertices") {
                    continue;
                }
            }
        }

        if (fields.size() < 6) {
            throw runtime_error(
                "Batch row must have at least 6 columns: "
                "motif,vertices,edge_prob,seed,input,debug"
            );
        }

        runs.push_back(
            parseBatchRun(
                fields
            )
        );
    }

    return runs;
}

// Launch PM2/MPI motif search benchmark
int main(
    int argc,
    char** argv
) {
    // Initialize MPI before reading rank specific behavior
    MPI_Init(
        &argc,
        &argv
    );

    int mpiRank;
    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &mpiRank
    );

    int mpiSize;
    MPI_Comm_size(
        MPI_COMM_WORLD,
        &mpiSize
    );

    int vertices = 0;
    double probability = 0.0;
    unsigned seed = 1;
    string motifName = "triangle";
    string inputFile = "";
    bool debug = false;
    int runs = 3;
    string batchFile = "";
    string outputFile = makeOutputCsvName();
    int nodesUsed = 0;
    int processesPerNode = 0;

    // Parse command line options used by benchmark scripts
    // Parse CLI options on every rank so benchmark inputs match
    for (int index = 1; index < argc; ++index) {
        const string argument = argv[index];

        if (argument == "--vertices" && index + 1 < argc) {
            vertices = stoi(
                argv[++index]
            );
        }
        else if (argument == "--edge-prob" && index + 1 < argc) {
            probability = stod(
                argv[++index]
            );
        }
        else if (argument == "--seed" && index + 1 < argc) {
            seed = static_cast<unsigned>(
                stoul(
                    argv[++index]
                )
            );
        }
        else if (argument == "--motif" && index + 1 < argc) {
            motifName = argv[++index];
        }
        else if (argument == "--input" && index + 1 < argc) {
            inputFile = argv[++index];
        }
        else if (argument == "--debug" && index + 1 < argc) {
            debug = stoi(
                argv[++index]
            ) != 0;
        }
        else if (argument == "--runs" && index + 1 < argc) {
            runs = stoi(
                argv[++index]
            );
        }
        else if (argument == "--batch-file" && index + 1 < argc) {
            batchFile = argv[++index];
        }
        else if (argument == "--output" && index + 1 < argc) {
            outputFile = argv[++index];
        }
        else if (argument == "--nodes-used" && index + 1 < argc) {
            nodesUsed = stoi(
                argv[++index]
            );
        }
        else if (argument == "--processes-per-node" && index + 1 < argc) {
            processesPerNode = stoi(
                argv[++index]
            );
        }
    }

    if (nodesUsed < 1) {
        nodesUsed = mpiSize;
    }

    if (processesPerNode < 1) {
        processesPerNode = mpiSize / nodesUsed;
    }

    if (runs < 1) {
        if (mpiRank == 0) {
            cerr
                << "Error: --runs must be at least 1\n";
        }

        MPI_Finalize();
        return 1;
    }

    PM2_MotifSearch runner;

    try {
        if (mpiRank == 0) {
            ofstream output(
                outputFile
            );

            if (!output.is_open()) {
                throw runtime_error(
                    "Could not open output file: " + outputFile
                );
            }

            writeCsvHeader(
                output
            );
        }

        // Batch mode runs each CSV row and appends 1 result row
        if (!batchFile.empty()) {
            vector<BatchRun> batchRuns = loadBatchFile(
                batchFile
            );

            for (const auto& batchRun : batchRuns) {
                if (batchRun.runs < 1) {
                    if (mpiRank == 0) {
                        cerr
                            << "Error: batch file run count must be at "
                            << "least 1\n";
                    }

                    MPI_Finalize();
                    return 1;
                }

                vector<PM2_MotifSearch::MotifSearchRunResult> results;
                results.reserve(
                    batchRun.runs
                );

                for (int run = 0; run < batchRun.runs; ++run) {
                    auto result = runner.runMotifSearchBenchmark(
                        batchRun.vertices,
                        batchRun.probability,
                        batchRun.seed,
                        batchRun.motifName,
                        batchRun.inputFile,
                        batchRun.debug
                    );

                    if (!result.success) {
                        MPI_Finalize();
                        return 1;
                    }

                    results.push_back(
                        result
                    );
                }

                if (mpiRank == 0) {
                    ofstream output(
                        outputFile,
                        ios::app
                    );

                    if (!output.is_open()) {
                        throw runtime_error(
                            "Could not append to output file: " + outputFile
                        );
                    }

                    writeCsvRow(
                        output,
                        batchRun.motifName,
                        batchRun.vertices,
                        batchRun.probability,
                        batchRun.seed,
                        mpiSize,
                        nodesUsed,
                        processesPerNode,
                        batchRun.inputFile,
                        results
                    );
                }
            }
        }
        else {
            // Single config mode repeats same run count
            vector<PM2_MotifSearch::MotifSearchRunResult> results;
            results.reserve(
                runs
            );

            for (int run = 0; run < runs; ++run) {
                auto result = runner.runMotifSearchBenchmark(
                    vertices,
                    probability,
                    seed,
                    motifName,
                    inputFile,
                    debug
                );

                if (!result.success) {
                    MPI_Finalize();
                    return 1;
                }

                results.push_back(
                    result
                );
            }

            if (mpiRank == 0) {
                double totalSearchTimeMs = 0.0;

                for (const auto& result : results) {
                    totalSearchTimeMs += result.searchTimeMs;
                }

                const double averageSearchTimeMs = totalSearchTimeMs
                    / results.size();
                const auto& lastResult = results.back();
                const int outputVertices = (
                    lastResult.vertices > 0
                        ? lastResult.vertices
                        : vertices
                );

                cout
                    << "\n-- Results "
                    << "------------------------------------------\n"
                    << "  Raw embeddings:   "
                    << lastResult.rawEmbeddings
                    << "\n"
                    << "  Automorphisms:    "
                    << lastResult.automorphisms
                    << "\n"
                    << "  Instance count:   "
                    << lastResult.instanceCount
                    << "\n"
                    << "  Edges:            "
                    << lastResult.edges
                    << "\n"
                    << "  Runs:             "
                    << runs
                    << "\n"
                    << "  Avg search time:  "
                    << fixed
                    << setprecision(
                        3
                    )
                    << averageSearchTimeMs
                    << " ms\n";

                for (int index = 0;
                     index < static_cast<int>(
                         results.size()
                     );
                     ++index) {
                    cout
                        << "  Run "
                        << index + 1
                        << ":           "
                        << fixed
                        << setprecision(
                            3
                        )
                        << results[index].searchTimeMs
                        << " ms\n";
                }

                cout
                    << "\n-- Configuration "
                    << "------------------------------------\n"
                    << "  Motif:            "
                    << motifName
                    << "\n"
                    << "  Vertices:         "
                    << outputVertices
                    << "\n"
                    << "  Edge prob:        "
                    << probability
                    << "\n"
                    << "  Seed:             "
                    << seed
                    << "\n"
                    << "  Ranks:            "
                    << mpiSize
                    << "\n"
                    << "  Logfile:          "
                    << deriveLogfileName(
                        "pm2",
                        motifName,
                        outputVertices,
                        probability,
                        seed,
                        mpiSize,
                        inputFile
                    )
                    << "\n"
                    << "  Output CSV:       "
                    << outputFile
                    << "\n"
                    << "========================================\n";

                ofstream output(
                    outputFile,
                    ios::app
                );

                if (!output.is_open()) {
                    throw runtime_error(
                        "Could not append to output file: " + outputFile
                    );
                }

                writeCsvRow(
                    output,
                    motifName,
                    vertices,
                    probability,
                    seed,
                    mpiSize,
                    nodesUsed,
                    processesPerNode,
                    inputFile,
                    results
                );
            }
        }
    }
    catch (const exception& exception) {
        if (mpiRank == 0) {
            cerr
                << "Error: "
                << exception.what()
                << "\n";
        }

        MPI_Finalize();
        return 1;
    }

    MPI_Finalize();
    return 0;
}
