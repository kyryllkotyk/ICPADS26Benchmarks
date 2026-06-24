#include "pm2motifsearch.h"

PM2_MotifSearch::MotifSearchRunResult
PM2_MotifSearch::runMotifSearchBenchmark(
    const int vertices,
    const double probability,
    const unsigned seed,
    const string& motifName,
    const string& inputFile,
    const bool debug
) {
    MotifSearchRunResult result;

    // Get rank information before splitting root vertices
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

    // Build supported motifs and select requested pattern
    const auto motifLibrary = builtinMotifs();
    const auto motifIterator = motifLibrary.find(
        motifName
    );

    if (motifIterator == motifLibrary.end()) {
        if (mpiRank == 0) {
            cerr
                << "Error: unknown motif '"
                << motifName
                << "'\n";

            for (const auto& motifEntry : motifLibrary) {
                cerr
                    << "  "
                    << motifEntry.first
                    << "\n";
            }
        }

        result.success = false;
        return result;
    }

    const MotifPattern& motif = motifIterator->second;
    CSRGraph graph;

    // Use file input when caller gives an edge list
    if (!inputFile.empty()) {
        if (debug && mpiRank == 0) {
            cout
                << "Loading graph from: "
                << inputFile
                << "\n";
        }

        graph = loadGraphFromFile(
            inputFile
        );
    }

    // Otherwise generate same random graph on every rank
    else if (vertices > 0) {
        if (debug && mpiRank == 0) {
            cout
                << "Generating Erdos-Renyi G("
                << vertices
                << ", "
                << probability
                << ") seed="
                << seed
                << "\n";
        }

        graph = generateErdosRenyi(
            vertices,
            probability,
            seed
        );
    }
    else {
        if (mpiRank == 0) {
            cerr
                << "Error: specify --vertices N --edge-prob P "
                << "or --input FILE\n";
        }

        result.success = false;
        return result;
    }

    // Assign a contiguous root vertex range to this rank
    const int baseChunk = graph.vertices / mpiSize;
    const int remainder = graph.vertices % mpiSize;
    const int vertexStart = mpiRank * baseChunk + min(
        mpiRank,
        remainder
    );
    const int vertexEnd = vertexStart + baseChunk + (
        mpiRank < remainder
            ? 1
            : 0
    );

    // Start timing after every rank has built its local graph copy
    MPI_Barrier(
        MPI_COMM_WORLD
    );
    const double startTime = MPI_Wtime();

    int64_t localCount = 0;

    // Search only root vertices owned by this rank
    for (int vertex = vertexStart; vertex < vertexEnd; ++vertex) {
        localCount += searchFromRoot(
            graph,
            motif,
            vertex
        );
    }

    // Sum raw embeddings across all ranks
    int64_t globalCount = 0;
    MPI_Reduce(
        &localCount,
        &globalCount,
        1,
        MPI_LONG_LONG,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    const double endTime = MPI_Wtime();
    const double localSearchMs = (endTime - startTime) * 1000.0;

    // Use slowest rank time as benchmark search time
    double maxSearchMs = 0.0;
    MPI_Reduce(
        &localSearchMs,
        &maxSearchMs,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        0,
        MPI_COMM_WORLD
    );

    // Only rank 0 owns final printable benchmark result fields
    if (mpiRank == 0) {
        result.vertices = graph.vertices;
        result.edges = graph.edges / 2;
        result.rawEmbeddings = globalCount;
        result.automorphisms = motif.automorphisms;
        result.instanceCount = globalCount / motif.automorphisms;
        result.searchTimeMs = maxSearchMs;
    }

    return result;
}

int PM2_MotifSearch::CSRGraph::degree(
    const int vertex
) const {
    return static_cast<int>(
        offsets[vertex + 1] - offsets[vertex]
    );
}

const int* PM2_MotifSearch::CSRGraph::firstNeighborOfVertex(
    const int vertex
) const {
    return adjacency.data() + offsets[vertex];
}

const int* PM2_MotifSearch::CSRGraph::afterLastNeighborOfVertex(
    const int vertex
) const {
    return adjacency.data() + offsets[vertex + 1];
}

bool PM2_MotifSearch::CSRGraph::hasEdge(
    const int vertex1,
    const int vertex2
) const {
    // Neighbor lists are sorted during graph construction
    return binary_search(
        firstNeighborOfVertex(
            vertex1
        ),
        afterLastNeighborOfVertex(
            vertex1
        ),
        vertex2
    );
}

bool PM2_MotifSearch::MotifPattern::hasEdge(
    const int motifVertex1,
    const int motifVertex2
) const {
    const int adjacencyIndex = motifVertex1 * motifVertices + motifVertex2;

    return adjacency[adjacencyIndex];
}

unordered_map<string, PM2_MotifSearch::MotifPattern>
PM2_MotifSearch::builtinMotifs() {
    unordered_map<string, MotifPattern> motifLibrary;

    // 3 vertex complete graph
    motifLibrary["triangle"] = makeMotif(
        "triangle",
        3,
        {
            { 0, 1 },
            { 0, 2 },
            { 1, 2 }
        }
    );

    // 3 vertex path
    motifLibrary["path3"] = makeMotif(
        "path3",
        3,
        {
            { 0, 1 },
            { 1, 2 }
        }
    );

    // 4 vertex cycle
    motifLibrary["square"] = makeMotif(
        "square",
        4,
        {
            { 0, 1 },
            { 1, 2 },
            { 2, 3 },
            { 3, 0 }
        }
    );

    // 4 vertex complete graph
    motifLibrary["4clique"] = makeMotif(
        "4clique",
        4,
        {
            { 0, 1 },
            { 0, 2 },
            { 0, 3 },
            { 1, 2 },
            { 1, 3 },
            { 2, 3 }
        }
    );

    // 4 vertex graph with 1 missing clique edge
    motifLibrary["diamond"] = makeMotif(
        "diamond",
        4,
        {
            { 0, 1 },
            { 0, 2 },
            { 0, 3 },
            { 1, 2 },
            { 2, 3 }
        }
    );

    // Triangle with 1 extra tail edge
    motifLibrary["tailed_triangle"] = makeMotif(
        "tailed_triangle",
        4,
        {
            { 0, 1 },
            { 0, 2 },
            { 1, 2 },
            { 2, 3 }
        }
    );

    // 3 leaf star around vertex 0
    motifLibrary["star3"] = makeMotif(
        "star3",
        4,
        {
            { 0, 1 },
            { 0, 2 },
            { 0, 3 }
        }
    );

    // 4 vertex path
    motifLibrary["path4"] = makeMotif(
        "path4",
        4,
        {
            { 0, 1 },
            { 1, 2 },
            { 2, 3 }
        }
    );

    // 4 leaf star around vertex 0
    motifLibrary["star4"] = makeMotif(
        "star4",
        5,
        {
            { 0, 1 },
            { 0, 2 },
            { 0, 3 },
            { 0, 4 }
        }
    );

    // 5 vertex cycle
    motifLibrary["5cycle"] = makeMotif(
        "5cycle",
        5,
        {
            { 0, 1 },
            { 1, 2 },
            { 2, 3 },
            { 3, 4 },
            { 4, 0 }
        }
    );

    // 5 vertex complete graph
    motifLibrary["5clique"] = makeMotif(
        "5clique",
        5,
        {
            { 0, 1 },
            { 0, 2 },
            { 0, 3 },
            { 0, 4 },
            { 1, 2 },
            { 1, 3 },
            { 1, 4 },
            { 2, 3 },
            { 2, 4 },
            { 3, 4 }
        }
    );

    return motifLibrary;
}

PM2_MotifSearch::CSRGraph PM2_MotifSearch::generateErdosRenyi(
    const int verticesCount,
    const double edgeProbability,
    const unsigned seed
) {
    CSRGraph graph;
    graph.vertices = verticesCount;

    // Store temporary adjacency lists before CSR flattening
    vector<vector<int>> temporaryAdjacency(
        verticesCount
    );

    // Generate repeatable edge decisions from given seed
    mt19937 randomGenerator(
        seed
    );
    uniform_real_distribution<double> distribution(
        0.0,
        1.0
    );

    // Visit each unique undirected vertex pair once
    for (int sourceVertex = 0;
         sourceVertex < verticesCount;
         ++sourceVertex) {
        for (int targetVertex = sourceVertex + 1;
             targetVertex < verticesCount;
             ++targetVertex) {
            if (distribution(
                    randomGenerator
                ) < edgeProbability) {
                temporaryAdjacency[sourceVertex].push_back(
                    targetVertex
                );
                temporaryAdjacency[targetVertex].push_back(
                    sourceVertex
                );
            }
        }
    }

    graph.offsets.resize(
        verticesCount + 1
    );
    graph.offsets[0] = 0;

    // Sort each neighbor list before binary edge checks
    for (int vertex = 0; vertex < verticesCount; ++vertex) {
        sort(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end()
        );

        graph.offsets[vertex + 1] = graph.offsets[vertex]
            + static_cast<int64_t>(
                temporaryAdjacency[vertex].size()
            );
    }

    graph.edges = graph.offsets[verticesCount];
    graph.adjacency.resize(
        graph.edges
    );

    // Flatten sorted neighbor lists into CSR adjacency storage
    for (int vertex = 0; vertex < verticesCount; ++vertex) {
        copy(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            graph.adjacency.begin() + graph.offsets[vertex]
        );
    }

    return graph;
}

PM2_MotifSearch::CSRGraph PM2_MotifSearch::loadGraphFromFile(
    const string& filename
) {
    ifstream inputFile(
        filename
    );

    if (!inputFile.is_open()) {
        cerr
            << "Error: cannot open "
            << filename
            << "\n";
        exit(
            1
        );
    }

    vector<pair<int, int>> fileEdges;
    int maxId = -1;
    string line;

    // Read noncomment edge pairs from file
    while (getline(
        inputFile,
        line
    )) {
        if (line.empty() || line[0] == '#' || line[0] == '%') {
            continue;
        }

        istringstream lineStream(
            line
        );

        int sourceVertex;
        int targetVertex;

        if (!(lineStream >> sourceVertex >> targetVertex)) {
            continue;
        }

        // Keep original edge list until vertex count is known
        fileEdges.push_back(
            {
                sourceVertex,
                targetVertex
            }
        );

        // Track largest vertex id to size graph storage
        maxId = max(
            maxId,
            max(
                sourceVertex,
                targetVertex
            )
        );
    }

    const int verticesCount = maxId + 1;
    vector<vector<int>> temporaryAdjacency(
        verticesCount
    );

    // Store every file edge as an undirected pair
    for (const auto& edge : fileEdges) {
        const int sourceVertex = edge.first;
        const int targetVertex = edge.second;

        temporaryAdjacency[sourceVertex].push_back(
            targetVertex
        );
        temporaryAdjacency[targetVertex].push_back(
            sourceVertex
        );
    }

    CSRGraph graph;
    graph.vertices = verticesCount;
    graph.offsets.resize(
        verticesCount + 1
    );
    graph.offsets[0] = 0;

    // Normalize each neighbor list for binary edge checks
    for (int vertex = 0; vertex < verticesCount; ++vertex) {
        sort(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end()
        );

        const auto uniqueEnd = unique(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end()
        );
        temporaryAdjacency[vertex].erase(
            uniqueEnd,
            temporaryAdjacency[vertex].end()
        );

        const auto selfLoopEnd = remove(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            vertex
        );
        temporaryAdjacency[vertex].erase(
            selfLoopEnd,
            temporaryAdjacency[vertex].end()
        );

        graph.offsets[vertex + 1] = graph.offsets[vertex]
            + static_cast<int64_t>(
                temporaryAdjacency[vertex].size()
            );
    }

    graph.edges = graph.offsets[verticesCount];
    graph.adjacency.resize(
        graph.edges
    );

    // Flatten normalized neighbor lists into CSR adjacency storage
    for (int vertex = 0; vertex < verticesCount; ++vertex) {
        copy(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            graph.adjacency.begin() + graph.offsets[vertex]
        );
    }

    return graph;
}

int PM2_MotifSearch::computeAutomorphisms(
    const MotifPattern& motifPattern
) {
    vector<int> permutation(
        motifPattern.motifVertices
    );

    // Start with identity permutation
    for (int vertex = 0;
         vertex < motifPattern.motifVertices;
         ++vertex) {
        permutation[vertex] = vertex;
    }

    int count = 0;

    // Check every vertex permutation because motif size is small
    do {
        bool isAutomorphism = true;

        // Compare each original motif edge with permuted edge
        for (int sourceVertex = 0;
             sourceVertex < motifPattern.motifVertices && isAutomorphism;
             ++sourceVertex) {
            for (int targetVertex = sourceVertex + 1;
                 targetVertex < motifPattern.motifVertices && isAutomorphism;
                 ++targetVertex) {
                const bool originalEdge = motifPattern.hasEdge(
                    sourceVertex,
                    targetVertex
                );
                const bool permutedEdge = motifPattern.hasEdge(
                    permutation[sourceVertex],
                    permutation[targetVertex]
                );

                if (originalEdge != permutedEdge) {
                    isAutomorphism = false;
                }
            }
        }

        if (isAutomorphism) {
            count++;
        }
    } while (next_permutation(
        permutation.begin(),
        permutation.end()
    ));

    return count;
}

PM2_MotifSearch::MotifPattern PM2_MotifSearch::makeMotif(
    const string& name,
    const int motifVertices,
    const vector<pair<int, int>>& edges
) {
    MotifPattern pattern;
    pattern.name = name;
    pattern.motifVertices = motifVertices;

    // Clear fixed adjacency matrix before setting motif edges
    memset(
        pattern.adjacency,
        0,
        sizeof(
            pattern.adjacency
        )
    );

    // Store each motif edge in both directions
    for (const auto& edge : edges) {
        const int sourceVertex = edge.first;
        const int targetVertex = edge.second;
        const int forwardIndex = sourceVertex * motifVertices + targetVertex;
        const int reverseIndex = targetVertex * motifVertices + sourceVertex;

        pattern.adjacency[forwardIndex] = true;
        pattern.adjacency[reverseIndex] = true;
    }

    // Save automorphism count for final instance normalization
    pattern.automorphisms = computeAutomorphisms(
        pattern
    );

    return pattern;
}

int64_t PM2_MotifSearch::searchFromRoot(
    const CSRGraph& graph,
    const MotifPattern& pattern,
    const int root
) {
    const int motifVertices = pattern.motifVertices;
    int64_t count = 0;

    vector<int> match(
        motifVertices,
        -1
    );
    vector<bool> used(
        graph.vertices,
        false
    );
    vector<vector<int>> required(
        motifVertices
    );

    // Precompute earlier motif vertices required by each depth
    for (int depth = 1; depth < motifVertices; ++depth) {
        for (int previousDepth = 0;
             previousDepth < depth;
             ++previousDepth) {
            if (pattern.hasEdge(
                    depth,
                    previousDepth
                )) {
                required[depth].push_back(
                    previousDepth
                );
            }
        }
    }

    function<void(int)> backtrack = [&](
        const int depth
    ) {
        if (depth == motifVertices) {
            count++;
            return;
        }

        const auto& requiredNeighbors = required[depth];

        // No required neighbor means any not used graph vertex can fit
        if (requiredNeighbors.empty()) {
            for (int candidate = 0;
                 candidate < graph.vertices;
                 ++candidate) {
                if (!used[candidate]) {
                    match[depth] = candidate;
                    used[candidate] = true;
                    backtrack(
                        depth + 1
                    );
                    used[candidate] = false;
                }
            }
        }
        else {
            int pivot = requiredNeighbors[0];

            // Choose required neighbor with smallest graph degree
            for (const int requiredIndex : requiredNeighbors) {
                if (graph.degree(
                        match[requiredIndex]
                    ) < graph.degree(
                        match[pivot]
                    )) {
                    pivot = requiredIndex;
                }
            }

            const int* neighborBegin = graph.firstNeighborOfVertex(
                match[pivot]
            );
            const int* neighborEnd = graph.afterLastNeighborOfVertex(
                match[pivot]
            );

            // Enumerate only neighbors of selected pivot match
            for (const int* iterator = neighborBegin;
                 iterator != neighborEnd;
                 ++iterator) {
                const int candidate = *iterator;

                if (used[candidate]) {
                    continue;
                }

                bool valid = true;

                // Candidate must connect to every required match
                for (const int requiredIndex : requiredNeighbors) {
                    if (requiredIndex == pivot) {
                        continue;
                    }

                    if (!graph.hasEdge(
                            match[requiredIndex],
                            candidate
                        )) {
                        valid = false;
                        break;
                    }
                }

                if (!valid) {
                    continue;
                }

                match[depth] = candidate;
                used[candidate] = true;
                backtrack(
                    depth + 1
                );
                used[candidate] = false;
            }
        }
    };

    // Root vertex is fixed by rank ownership loop
    match[0] = root;
    used[root] = true;
    backtrack(
        1
    );

    return count;
}
