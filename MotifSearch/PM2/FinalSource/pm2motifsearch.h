/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Motif Search                                                    *
 * Library: PM2/MPI                                                           *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef PM2_MOTIF_SEARCH_
#define PM2_MOTIF_SEARCH_

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <mpi.h>

#define PM2_MOTIF_SEARCH_MAX_K 6

using namespace std;

class PM2_MotifSearch
{
public:
    struct MotifSearchRunResult
    {
        int64_t edges = 0;
        int64_t rawEmbeddings = 0;
        int automorphisms = 0;
        int64_t instanceCount = 0;
        double searchTimeMs = 0.0;
        bool success = true;
        int vertices = 0;
    };

    MotifSearchRunResult runMotifSearchBenchmark(
        const int vertices,
        const double probability,
        const unsigned seed,
        const string& motifName,
        const string& inputFile,
        const bool debug
    );

private:
    struct CSRGraph
    {
        int vertices = 0;
        uint64_t edges = 0;

        // Store start offsets for each vertex neighbor list
        vector<int64_t> offsets;

        // Store all neighbor lists in 1 contiguous array
        vector<int> adjacency;

        // Return neighbor count for 1 vertex
        int degree(
            const int vertex
        ) const;

        // Return pointer to first neighbor for 1 vertex
        const int* firstNeighborOfVertex(
            const int vertex
        ) const;

        // Return pointer after last neighbor for 1 vertex
        const int* afterLastNeighborOfVertex(
            const int vertex
        ) const;

        // Check whether 2 vertices share an edge
        bool hasEdge(
            const int vertex1,
            const int vertex2
        ) const;
    };

    struct MotifPattern
    {
        string name;
        int motifVertices = 0;
        bool adjacency[PM2_MOTIF_SEARCH_MAX_K * PM2_MOTIF_SEARCH_MAX_K] = {};
        int automorphisms = 0;

        // Check whether 2 motif vertices share an edge
        bool hasEdge(
            const int motifVertex1,
            const int motifVertex2
        ) const;
    };

    // Build all supported motif patterns
    static unordered_map<string, MotifPattern> builtinMotifs();

    // Generate random graph from seeded Erdos Renyi model
    static CSRGraph generateErdosRenyi(
        const int verticesCount,
        const double edgeProbability,
        const unsigned seed
    );

    // Load graph from edge pair file
    static CSRGraph loadGraphFromFile(
        const string& filename
    );

    // Count edge preserving permutations of pattern
    static int computeAutomorphisms(
        const MotifPattern& motifPattern
    );

    // Build 1 motif pattern from vertex count and edge list
    static MotifPattern makeMotif(
        const string& name,
        const int motifVertices,
        const vector<pair<int, int>>& edges
    );

    // Search all embeddings rooted at 1 graph vertex
    static int64_t searchFromRoot(
        const CSRGraph& graph,
        const MotifPattern& pattern,
        const int root
    );
};

#endif
