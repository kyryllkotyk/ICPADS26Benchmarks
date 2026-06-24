#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using namespace std;

static void printUsage(
    const char* programName
) {
    cerr
        << "Usage: "
        << programName
        << " --vertices N --edge-prob P --seed S --output FILE\n";
}

static string requireValue(
    int& index,
    int argc,
    char** argv,
    const string& argument
) {
    if (index + 1 >= argc) {
        throw runtime_error(
            "Missing value for " + argument
        );
    }

    ++index;
    return argv[index];
}

int main(
    int argc,
    char** argv
) {
    int vertices = 0;
    double edgeProbability = 0.0;
    unsigned seed = 1;
    string outputFile;

    for (int index = 1; index < argc; ++index) {
        const string argument = argv[index];

        if (argument == "--vertices") {
            vertices = stoi(
                requireValue(
                    index,
                    argc,
                    argv,
                    argument
                )
            );
        }
        else if (argument == "--edge-prob") {
            edgeProbability = stod(
                requireValue(
                    index,
                    argc,
                    argv,
                    argument
                )
            );
        }
        else if (argument == "--seed") {
            seed = static_cast<unsigned>(
                stoul(
                    requireValue(
                        index,
                        argc,
                        argv,
                        argument
                    )
                )
            );
        }
        else if (argument == "--output") {
            outputFile = requireValue(
                index,
                argc,
                argv,
                argument
            );
        }
        else {
            throw runtime_error(
                "Unknown argument: " + argument
            );
        }
    }

    if (vertices <= 0 || edgeProbability < 0.0 || edgeProbability > 1.0
            || outputFile.empty()) {
        printUsage(
            argv[0]
        );
        return 1;
    }

    ofstream output(
        outputFile
    );

    if (!output.is_open()) {
        cerr
            << "Error: could not open output file: "
            << outputFile
            << "\n";
        return 1;
    }

    mt19937 randomGenerator(
        seed
    );
    uniform_real_distribution<double> distribution(
        0.0,
        1.0
    );

    int64_t edgeCount = 0;

    // Write each undirected edge once
    for (int sourceVertex = 0;
            sourceVertex < vertices;
            ++sourceVertex) {
        for (int targetVertex = sourceVertex + 1;
                targetVertex < vertices;
                ++targetVertex) {
            if (distribution(
                    randomGenerator
                ) < edgeProbability) {
                output
                    << sourceVertex
                    << ' '
                    << targetVertex
                    << '\n';

                ++edgeCount;
            }
        }
    }

    cout
        << "vertices="
        << vertices
        << " edge_prob="
        << edgeProbability
        << " seed="
        << seed
        << " edges="
        << edgeCount
        << " output="
        << outputFile
        << "\n";

    return 0;
}
