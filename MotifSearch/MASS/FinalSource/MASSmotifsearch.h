/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Motif Search                                                    *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef MASS_MOTIFSEARCH_
#define MASS_MOTIFSEARCH_

#include "MASS.h"
#include "MASSmotifsearchMotifPlace.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#define MASS_MOTIF_PLACES_HANDLE 1

class MASS_MotifSearch
{
public:
    // Run full MASS motif search benchmark from command line arguments
    int run(
        int argc,
        char* argv[]
    );

private:
    // Store motif pattern metadata and adjacency matrix
    struct MotifPattern
    {
        std::string name;
        int motifK = 0;

        // Flattened motif adjacency matrix
        bool adjacency[MAX_MOTIF_K * MAX_MOTIF_K] = {};

        // Number of equivalent labelings for this motif
        int automorphisms = 0;

        // Check whether 2 motif vertices are connected
        bool hasEdge(
            const int leftMotifVertex,
            const int rightMotifVertex
        ) const;
    };

    // Store graph in compressed sparse row format
    struct CSRGraph
    {
        int vertices = 0;
        int64_t directedEdges = 0;

        // offsets[v] gives start of vertex v adjacency list
        std::vector<int64_t> offsets;

        // Flattened directed adjacency lists
        std::vector<int> adjacency;
    };

    // Store benchmark options parsed from command line
    struct RunConfig
    {
        int vertices = 0;
        double edgeProbability = 0.0;
        unsigned seed = 42;
        std::string inputFile;
        std::string motifName = "triangle";
        int chunks = 0;
    };

    // Print expected MASS and benchmark argument format
    static void printUsage(
        const char* programName
    );

    // Parse graph, motif, seed, and chunk options
    static void parseArguments(
        int argc,
        char* argv[],
        RunConfig& config
    );

    // Generate an undirected Erdos Renyi graph in CSR form
    static CSRGraph generateErdosRenyi(
        const int vertices,
        const double edgeProbability,
        const unsigned seed
    );

    // Load graph edge list and convert it to CSR form
    static CSRGraph loadGraphFromFile(
        const std::string& filename
    );

    // Count label permutations that preserve motif structure
    static int computeAutomorphisms(
        const MotifPattern& motif
    );

    // Build a motif pattern from named edge pairs
    static MotifPattern makeMotif(
        const std::string& name,
        const int motifVertices,
        const std::vector<std::pair<int, int>>& edges
    );

    // Create supported motif patterns used by benchmark runs
    static std::unordered_map<std::string, MotifPattern> builtinMotifs();

    // Split vertices into chunks assigned to MASS Places
    static std::vector<int> buildChunkStarts(
        const CSRGraph& graph,
        int chunkCount
    );

    // Pack graph, motif, and chunk metadata for Place initialization
    static std::vector<char> serializeBenchmarkData(
        const CSRGraph& graph,
        const MotifPattern& motif,
        const std::vector<int>& chunkStarts
    );
};

#endif
