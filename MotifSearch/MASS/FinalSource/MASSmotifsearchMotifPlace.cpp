#include "MASSmotifsearchMotifPlace.h"

#ifndef LOGGING
static const bool printOutput = false;
#else
static const bool printOutput = true;
#endif

// MASS factory hook used when Places creates MotifPlace objects
extern "C" Place* instantiate(
    void* argument
) {
    return new MotifPlace(
        argument
    );
}

// MASS cleanup hook for Place objects
extern "C" void destroy(
    Place* object
) {
    delete object;
}

// Shared graph state is loaded once per process
std::mutex MotifPlace::initMutex;
bool MotifPlace::graphReady = false;
int MotifPlace::vertexCount = 0;
int64_t MotifPlace::directedEdgeCount = 0;
std::vector<int64_t> MotifPlace::offsets;
std::vector<int> MotifPlace::adjacency;

// Shared motif state is reused by all local Places
int MotifPlace::motifK = 0;
int MotifPlace::motifAutomorphisms = 0;
bool MotifPlace::motifAdjacency[MAX_MOTIF_K * MAX_MOTIF_K] = {};
RequiredTable MotifPlace::requiredTable;

// Chunk table maps each Place index to a root vertex range
std::vector<int> MotifPlace::chunkStarts;

namespace
{
    // Reuse search buffers per thread to avoid repeated allocation
    thread_local std::vector<int> threadMatch;
    thread_local std::vector<uint64_t> threadUsedEpoch;
    thread_local uint64_t threadCurrentEpoch = 0;

    // Store mutable state used by recursive backtracking
    struct BacktrackSearch
    {
        const std::vector<int64_t>& offsets;
        const std::vector<int>& adjacency;
        const RequiredTable& requiredTable;
        int* match;
        uint64_t* usedEpoch;
        uint64_t currentEpoch;
        int vertexCount;
        int motifK;
        int64_t count = 0;

        int degree(
            const int vertex
        ) const {
            return static_cast<int>(
                offsets[vertex + 1] - offsets[vertex]
                );
        }

        const int* firstNeighborOfVertex(
            const int vertex
        ) const {
            return adjacency.data() + offsets[vertex];
        }

        const int* afterLastNeighborOfVertex(
            const int vertex
        ) const {
            return adjacency.data() + offsets[vertex + 1];
        }

        // Use sorted adjacency row for fast edge existence check
        bool hasEdge(
            const int leftVertex,
            const int rightVertex
        ) const {
            return std::binary_search(
                firstNeighborOfVertex(
                    leftVertex
                ),
                afterLastNeighborOfVertex(
                    leftVertex
                ),
                rightVertex
            );
        }

        // Pick required matched vertex with smallest degree as pivot
        void pickPivot(
            const int depth,
            const int8_t* requiredVertices,
            const int requiredCount,
            int& pivotIndex,
            int& pivotVertex
        ) const {
            pivotIndex = requiredVertices[0];
            pivotVertex = match[pivotIndex];

            int pivotDegree = degree(
                pivotVertex
            );

            // Lower degree pivot reduces candidate fanout
            for (int index = 1; index < requiredCount; index++) {
                const int requiredVertex = match[requiredVertices[index]];
                const int requiredDegree = degree(
                    requiredVertex
                );

                if (requiredDegree < pivotDegree) {
                    pivotIndex = requiredVertices[index];
                    pivotVertex = requiredVertex;
                    pivotDegree = requiredDegree;
                }
            }
        }

        // Search connected motifs by expanding from matched neighbors
        void runConnected(
            const int depth
        ) {
            if (depth == motifK) {
                count++;
                return;
            }

            // Required table lists prior motif vertices connected to depth
            const int8_t* requiredVertices = requiredTable.row(
                depth
            );

            const int requiredCount = requiredTable.count(
                depth
            );

            int pivotIndex = 0;
            int pivotVertex = 0;

            pickPivot(
                depth,
                requiredVertices,
                requiredCount,
                pivotIndex,
                pivotVertex
            );

            int requiredGraphVertices[MAX_MOTIF_K];

            // Store already matched graph vertices needed at this depth
            for (int index = 0; index < requiredCount; index++) {
                requiredGraphVertices[index] = match[requiredVertices[index]];
            }

            // Candidate vertices come from pivot adjacency list
            for (
                const int* neighbor = firstNeighborOfVertex(
                    pivotVertex
                );
                neighbor != afterLastNeighborOfVertex(
                    pivotVertex
                );
                neighbor++
                ) {
                const int candidate = *neighbor;

                if (usedEpoch[candidate] == currentEpoch) {
                    continue;
                }

                bool valid = true;

                // Candidate must connect to every required matched vertex
                for (
                    int index = 0;
                    index < requiredCount && valid;
                    index++
                    ) {
                    if (requiredVertices[index] == pivotIndex) {
                        continue;
                    }

                    if (!hasEdge(
                        requiredGraphVertices[index],
                        candidate
                    )) {
                        valid = false;
                    }
                }

                if (!valid) {
                    continue;
                }

                usedEpoch[candidate] = currentEpoch;
                match[depth] = candidate;

                runConnected(
                    depth + 1
                );

                // Unmark by changing epoch value instead of clearing array
                usedEpoch[candidate] = currentEpoch - 1;
            }
        }

        // Search motifs that may have disconnected expansion points
        void runGeneral(
            const int depth
        ) {
            if (depth == motifK) {
                count++;
                return;
            }

            const int8_t* requiredVertices = requiredTable.row(
                depth
            );

            const int requiredCount = requiredTable.count(
                depth
            );

            if (requiredCount == 0) {
                // Disconnected motifs may choose any graph vertex not used
                for (int candidate = 0; candidate < vertexCount; candidate++) {
                    if (usedEpoch[candidate] == currentEpoch) {
                        continue;
                    }

                    usedEpoch[candidate] = currentEpoch;
                    match[depth] = candidate;

                    runGeneral(
                        depth + 1
                    );

                    usedEpoch[candidate] = currentEpoch - 1;
                }

                return;
            }

            int pivotIndex = 0;
            int pivotVertex = 0;

            pickPivot(
                depth,
                requiredVertices,
                requiredCount,
                pivotIndex,
                pivotVertex
            );

            int requiredGraphVertices[MAX_MOTIF_K];

            // Cache matched vertices that candidate must connect to
            for (int index = 0; index < requiredCount; index++) {
                requiredGraphVertices[index] = match[requiredVertices[index]];
            }

            // Connected part can still expand from pivot adjacency list
            for (
                const int* neighbor = firstNeighborOfVertex(
                    pivotVertex
                );
                neighbor != afterLastNeighborOfVertex(
                    pivotVertex
                );
                neighbor++
                ) {
                const int candidate = *neighbor;

                if (usedEpoch[candidate] == currentEpoch) {
                    continue;
                }

                bool valid = true;

                // Enforce all motif edges to already matched vertices
                for (
                    int index = 0;
                    index < requiredCount && valid;
                    index++
                    ) {
                    if (requiredVertices[index] == pivotIndex) {
                        continue;
                    }

                    if (!hasEdge(
                        requiredGraphVertices[index],
                        candidate
                    )) {
                        valid = false;
                    }
                }

                if (!valid) {
                    continue;
                }

                usedEpoch[candidate] = currentEpoch;
                match[depth] = candidate;

                runGeneral(
                    depth + 1
                );

                usedEpoch[candidate] = currentEpoch - 1;
            }
        }
    };

    // Count raw embeddings that start from 1 assigned root vertex
    int64_t searchFromRoot(
        const int root,
        const int vertexCount,
        const std::vector<int64_t>& offsets,
        const std::vector<int>& adjacency,
        const RequiredTable& requiredTable
    ) {
        const int motifK = requiredTable.motifK;

        if (static_cast<int>(threadMatch.size()) < motifK) {
            threadMatch.assign(
                motifK,
                -1
            );
        }

        if (static_cast<int>(threadUsedEpoch.size()) < vertexCount) {
            threadUsedEpoch.assign(
                vertexCount,
                0
            );
        }

        // New epoch makes prior marks irrelevant without clearing vector
        threadCurrentEpoch++;

        const uint64_t currentEpoch = threadCurrentEpoch;

        // Root is fixed as first motif vertex
        threadMatch[0] = root;
        threadUsedEpoch[root] = currentEpoch;

        BacktrackSearch searchState{
            offsets,
            adjacency,
            requiredTable,
            threadMatch.data(),
            threadUsedEpoch.data(),
            currentEpoch,
            vertexCount,
            motifK
        };

        // Connected motifs can immediately search from neighbor lists
        if (requiredTable.connectedAtDepth1) {
            searchState.runConnected(
                1
            );
        }
        else {
            searchState.runGeneral(
                1
            );
        }

        return searchState.count;
    }
}

const int8_t* RequiredTable::row(
    const int depth
) const noexcept {
    return data.data() + depth * MAX_MOTIF_K;
}

int RequiredTable::count(
    const int depth
) const noexcept {
    return length[depth];
}

MotifPlace::MotifPlace(
    void* argument
) :
    Place(
        argument
    ),
    chunkId(
        -1
    ),
    vertexStart(
        0
    ),
    vertexEnd(
        0
    ),
    localCount(
        0
    ),
    localSearchMs(
        0.0
    ) {
    outMessage_size = 0;
    inMessage_size = 0;
}

int MotifPlace::degree(
    const int vertex
) {
    return static_cast<int>(
        offsets[vertex + 1] - offsets[vertex]
        );
}

const int* MotifPlace::firstNeighborOfVertex(
    const int vertex
) {
    return adjacency.data() + offsets[vertex];
}

const int* MotifPlace::afterLastNeighborOfVertex(
    const int vertex
) {
    return adjacency.data() + offsets[vertex + 1];
}

bool MotifPlace::hasEdge(
    const int leftVertex,
    const int rightVertex
) {
    return std::binary_search(
        firstNeighborOfVertex(
            leftVertex
        ),
        afterLastNeighborOfVertex(
            leftVertex
        ),
        rightVertex
    );
}

bool MotifPlace::motifEdge(
    const int leftMotifVertex,
    const int rightMotifVertex
) {
    return motifAdjacency[
        leftMotifVertex * motifK + rightMotifVertex
    ];
}

// Initialize Place chunk state and shared graph state
void* MotifPlace::init(
    void* argument
) {
    chunkId = index[0];
    localCount = 0;

    if (argument == nullptr) {
        return nullptr;
    }

    // Only 1 Place per process should unpack shared graph arrays
    std::lock_guard<std::mutex> lock(
        initMutex
    );

    if (!graphReady) {
        const char* buffer = static_cast<const char*>(
            argument
            );

        // Config header is stored at start of serialized buffer
        const MotifSearchConfig* config =
            reinterpret_cast<const MotifSearchConfig*>(
                buffer
                );

        vertexCount = config->n;
        directedEdgeCount = config->m;
        motifK = config->motifK;
        motifAutomorphisms = config->motifAutomorphisms;

        std::memcpy(
            motifAdjacency,
            config->motifAdjacency,
            sizeof(
                motifAdjacency
                )
        );

        // Pointer advances through serialized payload sections
        const char* pointer = buffer + sizeof(
            MotifSearchConfig
            );

        // Read chunk boundaries before graph arrays
        chunkStarts.resize(
            config->chunkCount + 1
        );

        std::memcpy(
            chunkStarts.data(),
            pointer,
            (config->chunkCount + 1) * sizeof(
                int
                )
        );

        pointer += (config->chunkCount + 1) * sizeof(
            int
            );

        // Read CSR offsets shared by all Places on this process
        offsets.resize(
            vertexCount + 1
        );

        std::memcpy(
            offsets.data(),
            pointer,
            (vertexCount + 1) * sizeof(
                int64_t
                )
        );

        pointer += (vertexCount + 1) * sizeof(
            int64_t
            );

        // Read flattened adjacency list shared by all Places
        adjacency.resize(
            directedEdgeCount
        );

        std::memcpy(
            adjacency.data(),
            pointer,
            directedEdgeCount * sizeof(
                int
                )
        );

        // Build compact table of prior motif vertices required by depth
        requiredTable = {};
        requiredTable.motifK = motifK;

        for (int depth = 1; depth < motifK; depth++) {
            int requiredCount = 0;

            for (int prior = 0; prior < depth; prior++) {
                if (!motifEdge(
                    depth,
                    prior
                )) {
                    continue;
                }

                requiredTable.data[
                    depth * MAX_MOTIF_K + requiredCount
                ] = static_cast<int8_t>(
                    prior
                    );

                    requiredCount++;
            }

            requiredTable.length[depth] = static_cast<int8_t>(
                requiredCount
                );
        }

        // Depth 1 connection determines which search path is cheaper
        requiredTable.connectedAtDepth1 = requiredTable.length[1] > 0;
        graphReady = true;

        if (printOutput) {
            std::cerr
                << "[MotifPlace] graph initialized n="
                << vertexCount
                << " m="
                << directedEdgeCount
                << " motif="
                << config->motifName
                << " k="
                << motifK
                << " automorphisms="
                << motifAutomorphisms
                << " chunks="
                << config->chunkCount
                << "\n";
        }
    }

    // Each Place maps its linear index to 1 root range
    if (
        chunkId >= 0
        && chunkId + 1 < static_cast<int>(
            chunkStarts.size()
            )
        ) {
        vertexStart = chunkStarts[chunkId];
        vertexEnd = chunkStarts[chunkId + 1];
    }
    else {
        // Invalid chunk id gets empty range instead of failing search
        vertexStart = 0;
        vertexEnd = 0;
    }

    return nullptr;
}

// Search all root vertices assigned to this Place
void* MotifPlace::search(
    void* argument
) {
    (void)argument;

    localCount = 0;
    localSearchMs = 0.0;

    const auto searchStart = std::chrono::high_resolution_clock::now();

    for (int vertex = vertexStart; vertex < vertexEnd; vertex++) {
        localCount += searchFromRoot(
            vertex,
            vertexCount,
            offsets,
            adjacency,
            requiredTable
        );
    }

    const auto searchEnd = std::chrono::high_resolution_clock::now();

    localSearchMs = std::chrono::duration<double, std::milli>(
        searchEnd - searchStart
    ).count();

    return nullptr;
}

void* MotifPlace::getCount(
    void* argument
) {
    (void)argument;

    int64_t* output = new int64_t;
    *output = localCount;

    return output;
}

void* MotifPlace::getSearchMs(
    void* argument
) {
    (void)argument;

    double* output = new double;
    *output = localSearchMs;

    return output;
}