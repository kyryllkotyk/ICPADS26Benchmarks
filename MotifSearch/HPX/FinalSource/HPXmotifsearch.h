/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Motif Search                                                    *
 * Library: HPX                                                               *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef HPX_MOTIFSEARCH_
#define HPX_MOTIFSEARCH_

#ifndef HPX_UTIL_FROM_STRING_INCLUDED
#define HPX_UTIL_FROM_STRING_INCLUDED
#include <hpx/util/bad_lexical_cast.hpp>
#include <hpx/util/from_string.hpp>
#endif

#include <hpx/algorithm.hpp>
#include <hpx/execution.hpp>
#include <hpx/executors/parallel_executor.hpp>
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#define HPX_MOTIF_MAX_K 6

namespace hpx
{
    extern char const hpx_check_boost_version_108900[];
    extern char const hpx_check_boost_version_107500[];
}

// Store graph in compressed sparse row format
struct CSRGraph
{
    int vertices = 0;
    int64_t directedEdges = 0;

    // offsets[v] gives start of vertex v adjacency list
    std::vector<int64_t> offsets;

    // Flattened directed adjacency lists
    std::vector<int> adjacency;

    int degree(
        const int vertex
    ) const noexcept;

    const int* firstNeighborOfVertex(
        const int vertex
    ) const noexcept;

    const int* afterLastNeighborOfVertex(
        const int vertex
    ) const noexcept;

    // Check graph edge existence using sorted adjacency rows
    bool hasEdge(
        const int leftVertex,
        const int rightVertex
    ) const noexcept;

    template <typename Archive>
    void serialize(
        Archive& archive,
        unsigned version
    );
};

// Store motif pattern metadata and adjacency matrix
struct MotifPattern
{
    std::string name;
    int motifK = 0;
    int automorphisms = 0;

    // Flattened motif adjacency matrix
    bool adjacency[HPX_MOTIF_MAX_K * HPX_MOTIF_MAX_K] = {};

    // Check motif edge existence using flattened adjacency matrix
    bool hasEdge(
        const int leftMotifVertex,
        const int rightMotifVertex
    ) const noexcept;

    template <typename Archive>
    void serialize(
        Archive& archive,
        unsigned version
    );
};

// Precomputed motif dependencies for each search depth
struct RequiredTable
{
    // Prior motif vertices that must connect to current depth
    std::array<int8_t, HPX_MOTIF_MAX_K* HPX_MOTIF_MAX_K> data{};

    // Number of required prior vertices for each depth
    std::array<int8_t, HPX_MOTIF_MAX_K> length{};

    int motifK = 0;
    bool connectedAtDepth1 = false;

    const int8_t* row(
        const int depth
    ) const noexcept;

    int count(
        const int depth
    ) const noexcept;
};

// Store one worker shard's search result
struct WorkerRunResult
{
    int64_t count = 0;
    double searchMs = 0.0;

    template <typename Archive>
    void serialize(
        Archive& archive,
        unsigned version
    );
};

using SearchClock = std::chrono::high_resolution_clock;

// Generate an undirected Erdos Renyi graph in CSR form
CSRGraph generateErdosRenyi(
    const int vertices,
    const double edgeProbability,
    const unsigned seed
);

// Load graph edge list and convert it to CSR form
CSRGraph loadGraphFromFile(
    const std::string& filename
);

// Count label permutations that preserve motif structure
int computeAutomorphisms(
    const MotifPattern& motif
);

// Build a motif pattern from named edge pairs
MotifPattern makeMotif(
    const std::string& name,
    const int motifVertices,
    const std::vector<std::pair<int, int>>& edges
);

// Create supported motif patterns used by benchmark runs
std::unordered_map<std::string, MotifPattern> builtinMotifs();

// Build compact table of prior motif vertices required by depth
RequiredTable buildRequiredTable(
    const MotifPattern& motif
);

// Split root vertices into worker chunks
std::vector<int> buildChunkStarts(
    const int vertices,
    int chunkCount
);

// Convert clock interval to milliseconds
double elapsedMs(
    const SearchClock::time_point start,
    const SearchClock::time_point end
);

// Search all root vertices in one local vertex range
int64_t searchRootsSequential(
    const CSRGraph& graph,
    const RequiredTable& requiredTable,
    const int startVertex,
    const int endVertex
);

// Search worker shard range assigned to one HPX locality
WorkerRunResult searchLocalBlockShards(
    const CSRGraph& graph,
    const RequiredTable& requiredTable,
    const std::vector<int>& chunkStarts,
    const int firstWorker,
    const int lastWorker
);

// Run full HPX motif search benchmark
int hpx_main(
    hpx::program_options::variables_map& variablesMap
);

#endif