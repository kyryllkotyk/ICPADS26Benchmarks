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

#ifndef MASS_MOTIFSEARCH_MOTIFPLACE_
#define MASS_MOTIFSEARCH_MOTIFPLACE_

#include "MASS_base.h"
#include "MethodRegistry.h"
#include "Place.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#define MAX_MOTIF_K 6

 // Serialized benchmark metadata received by each Place
struct MotifSearchConfig
{
    int n;
    int64_t m;
    int motifK;
    int motifAutomorphisms;
    int chunkCount;
    int pad;

    // Flattened motif adjacency matrix
    bool motifAdjacency[MAX_MOTIF_K * MAX_MOTIF_K];

    // Motif name is used only for debug output
    char motifName[32];
};

// Precomputed motif dependencies for each search depth
struct RequiredTable
{
    // Prior motif vertices that must connect to current depth
    std::array<int8_t, MAX_MOTIF_K* MAX_MOTIF_K> data{};

    // Number of required prior vertices for each depth
    std::array<int8_t, MAX_MOTIF_K> length{};

    int motifK = 0;
    bool connectedAtDepth1 = false;

    const int8_t* row(
        const int depth
    ) const noexcept;

    int count(
        const int depth
    ) const noexcept;
};

class MotifPlace : public Place
{
public:
    MotifPlace(
        void* argument
    );

    // Unpack graph data and assign this Place's vertex range
    void* init(
        void* argument
    );

    // Search all root vertices owned by this Place
    void* search(
        void* argument
    );

    // Return local raw embedding count
    void* getCount(
        void* argument
    );

    // Return local search time in milliseconds
    void* getSearchMs(
        void* argument
    );

    MASS_DISPATCH_TABLE(
        MotifPlace,
        MASS_METHOD(
            MotifPlace,
            init
        ),
        MASS_METHOD(
            MotifPlace,
            search
        ),
        MASS_METHOD(
            MotifPlace,
            getCount
        ),
        MASS_METHOD(
            MotifPlace,
            getSearchMs
        )
    )

private:
    // Per Place chunk ownership and result state
    int chunkId;
    int vertexStart;
    int vertexEnd;
    int64_t localCount;
    double localSearchMs;

    // Shared graph state is loaded once per process
    static std::mutex initMutex;
    static bool graphReady;
    static int vertexCount;
    static int64_t directedEdgeCount;
    static std::vector<int64_t> offsets;
    static std::vector<int> adjacency;

    // Shared motif state is reused by all local Places
    static int motifK;
    static int motifAutomorphisms;
    static bool motifAdjacency[MAX_MOTIF_K * MAX_MOTIF_K];
    static RequiredTable requiredTable;

    // Chunk table maps Place index to root vertex range
    static std::vector<int> chunkStarts;

    static int degree(
        const int vertex
    );

    static const int* firstNeighborOfVertex(
        const int vertex
    );

    static const int* afterLastNeighborOfVertex(
        const int vertex
    );

    // Check graph edge existence using sorted adjacency rows
    static bool hasEdge(
        const int leftVertex,
        const int rightVertex
    );

    // Check motif edge existence using flattened adjacency matrix
    static bool motifEdge(
        const int leftMotifVertex,
        const int rightMotifVertex
    );
};

#endif