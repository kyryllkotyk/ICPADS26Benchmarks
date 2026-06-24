#include "HPXmotifsearch.h"

namespace hpx
{
    char const hpx_check_boost_version_108900[] = "";
    char const hpx_check_boost_version_107500[] = "";
}

namespace
{
    // Reuse search buffers per HPX worker thread
    thread_local std::vector<int> threadMatch;
    thread_local std::vector<uint64_t> threadEpoch;
    thread_local uint64_t threadCurrentEpoch = 0;

    // Replicated benchmark state stored once per locality
    std::unique_ptr<CSRGraph> globalGraph;
    MotifPattern globalMotif;
    RequiredTable globalRequiredTable;
    std::vector<int> globalChunkStarts;

    // Mutable state for one recursive motif search
    struct BacktrackSearch
    {
        const CSRGraph& graph;
        const RequiredTable& requiredTable;
        int* match;
        uint64_t* epoch;
        uint64_t currentEpoch;
        int motifK;
        int64_t count = 0;

        // Pick lowest degree matched neighbor to reduce candidate fanout
        void pickPivot(
            const int depth,
            const int8_t* requiredVertices,
            const int requiredCount,
            int& pivotIndex,
            int& pivotVertex
        ) const noexcept {
            pivotIndex = requiredVertices[0];
            pivotVertex = match[pivotIndex];

            int pivotDegree = graph.degree(
                pivotVertex
            );

            for (int index = 1; index < requiredCount; index++) {
                const int requiredVertex = match[requiredVertices[index]];
                const int requiredDegree = graph.degree(
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

            // Required vertices are prior motif vertices adjacent to depth
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

            int requiredGraphVertices[HPX_MOTIF_MAX_K];

            // Cache already matched graph vertices needed at this depth
            for (int index = 0; index < requiredCount; index++) {
                requiredGraphVertices[index] = match[requiredVertices[index]];
            }

            // Candidate vertices come from pivot adjacency list
            for (
                const int* neighbor = graph.firstNeighborOfVertex(
                    pivotVertex
                );
                neighbor != graph.afterLastNeighborOfVertex(
                    pivotVertex
                );
                neighbor++
                ) {
                const int candidate = *neighbor;

                if (epoch[candidate] == currentEpoch) {
                    continue;
                }

                bool valid = true;

                // Candidate must connect to all required matched vertices
                for (
                    int index = 0;
                    index < requiredCount && valid;
                    index++
                    ) {
                    if (requiredVertices[index] == pivotIndex) {
                        continue;
                    }

                    if (!graph.hasEdge(
                        requiredGraphVertices[index],
                        candidate
                    )) {
                        valid = false;
                    }
                }

                if (!valid) {
                    continue;
                }

                epoch[candidate] = currentEpoch;
                match[depth] = candidate;

                runConnected(
                    depth + 1
                );

                // Unmark without clearing full epoch array
                epoch[candidate] = currentEpoch - 1;
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
                // Disconnected depth can choose any graph vertex not used
                for (
                    int candidate = 0;
                    candidate < graph.vertices;
                    candidate++
                    ) {
                    if (epoch[candidate] == currentEpoch) {
                        continue;
                    }

                    epoch[candidate] = currentEpoch;
                    match[depth] = candidate;

                    runGeneral(
                        depth + 1
                    );

                    epoch[candidate] = currentEpoch - 1;
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

            int requiredGraphVertices[HPX_MOTIF_MAX_K];

            // Cache matched vertices that candidate must connect to
            for (int index = 0; index < requiredCount; index++) {
                requiredGraphVertices[index] = match[requiredVertices[index]];
            }

            for (
                const int* neighbor = graph.firstNeighborOfVertex(
                    pivotVertex
                );
                neighbor != graph.afterLastNeighborOfVertex(
                    pivotVertex
                );
                neighbor++
                ) {
                const int candidate = *neighbor;

                if (epoch[candidate] == currentEpoch) {
                    continue;
                }

                bool valid = true;

                // Enforce motif edges to all required prior vertices
                for (
                    int index = 0;
                    index < requiredCount && valid;
                    index++
                    ) {
                    if (requiredVertices[index] == pivotIndex) {
                        continue;
                    }

                    if (!graph.hasEdge(
                        requiredGraphVertices[index],
                        candidate
                    )) {
                        valid = false;
                    }
                }

                if (!valid) {
                    continue;
                }

                epoch[candidate] = currentEpoch;
                match[depth] = candidate;

                runGeneral(
                    depth + 1
                );

                epoch[candidate] = currentEpoch - 1;
            }
        }
    };

    // Count raw embeddings that start from one root vertex
    int64_t searchFromRoot(
        const CSRGraph& graph,
        const RequiredTable& requiredTable,
        const int root
    ) {
        const int motifK = requiredTable.motifK;

        if (static_cast<int>(threadMatch.size()) < motifK) {
            threadMatch.assign(
                motifK,
                -1
            );
        }

        if (static_cast<int>(threadEpoch.size()) < graph.vertices) {
            threadEpoch.assign(
                graph.vertices,
                0
            );
        }

        // New epoch avoids clearing the used marker vector
        threadCurrentEpoch++;

        threadMatch[0] = root;
        threadEpoch[root] = threadCurrentEpoch;

        BacktrackSearch searchState{
            graph,
            requiredTable,
            threadMatch.data(),
            threadEpoch.data(),
            threadCurrentEpoch,
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

// Store graph and motif state on one locality
void initLocalityImpl(
    CSRGraph graph,
    MotifPattern motif,
    std::vector<int> chunkStarts
) {
    globalGraph = std::make_unique<CSRGraph>(
        std::move(
            graph
        )
    );

    globalMotif = motif;

    // Build locality copy of required table from motif metadata
    globalRequiredTable = buildRequiredTable(
        motif
    );

    globalChunkStarts = std::move(
        chunkStarts
    );
}

HPX_PLAIN_ACTION(
    initLocalityImpl,
    InitLocalityAction
)

// Empty action used to synchronize localities before timed search
void searchBarrierImpl() {
}

HPX_PLAIN_ACTION(
    searchBarrierImpl,
    SearchBarrierAction
)

// Run worker chunks assigned to one locality
WorkerRunResult localityWorkerImpl(
    const int localityId,
    const int workersPerLocality
) {
    if (
        !globalGraph
        || globalChunkStarts.size() < 2
        ) {
        return {};
    }

    const int firstWorker = localityId * workersPerLocality;

    const int lastWorker = std::min(
        firstWorker + workersPerLocality,
        static_cast<int>(
            globalChunkStarts.size()
            ) - 1
    );

    // Search only worker chunks assigned to this locality
    return searchLocalBlockShards(
        *globalGraph,
        globalRequiredTable,
        globalChunkStarts,
        firstWorker,
        lastWorker
    );
}

HPX_PLAIN_ACTION(
    localityWorkerImpl,
    LocalityWorkerAction
)

int CSRGraph::degree(
    const int vertex
) const noexcept {
    return static_cast<int>(
        offsets[vertex + 1] - offsets[vertex]
        );
}

const int* CSRGraph::firstNeighborOfVertex(
    const int vertex
) const noexcept {
    return adjacency.data() + offsets[vertex];
}

const int* CSRGraph::afterLastNeighborOfVertex(
    const int vertex
) const noexcept {
    return adjacency.data() + offsets[vertex + 1];
}

bool CSRGraph::hasEdge(
    const int leftVertex,
    const int rightVertex
) const noexcept {
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

template <typename Archive>
void CSRGraph::serialize(
    Archive& archive,
    unsigned version
) {
    (void)version;

    archive
        & vertices
        & directedEdges
        & offsets
        & adjacency;
}

bool MotifPattern::hasEdge(
    const int leftMotifVertex,
    const int rightMotifVertex
) const noexcept {
    return adjacency[
        leftMotifVertex * motifK + rightMotifVertex
    ];
}

template <typename Archive>
void MotifPattern::serialize(
    Archive& archive,
    unsigned version
) {
    (void)version;

    archive
        & name
        & motifK
        & automorphisms;

    for (int index = 0; index < HPX_MOTIF_MAX_K * HPX_MOTIF_MAX_K; index++) {
        archive& adjacency[index];
    }
}

const int8_t* RequiredTable::row(
    const int depth
) const noexcept {
    return data.data() + depth * HPX_MOTIF_MAX_K;
}

int RequiredTable::count(
    const int depth
) const noexcept {
    return length[depth];
}

template <typename Archive>
void WorkerRunResult::serialize(
    Archive& archive,
    unsigned version
) {
    (void)version;

    archive
        & count
        & searchMs;
}

CSRGraph generateErdosRenyi(
    const int vertices,
    const double edgeProbability,
    const unsigned seed
) {
    CSRGraph graph;
    graph.vertices = vertices;

    std::vector<std::vector<int>> temporaryAdjacency(
        vertices
    );

    std::mt19937 generator(
        seed
    );

    std::uniform_real_distribution<double> distribution(
        0.0,
        1.0
    );

    // Build undirected graph by testing each upper triangle edge
    for (int source = 0; source < vertices; source++) {
        for (
            int destination = source + 1;
            destination < vertices;
            destination++
            ) {
            if (distribution(
                generator
            ) < edgeProbability) {
                temporaryAdjacency[source].push_back(
                    destination
                );

                temporaryAdjacency[destination].push_back(
                    source
                );
            }
        }
    }

    graph.offsets.resize(
        vertices + 1
    );

    graph.offsets[0] = 0;

    // Sort adjacency rows and build CSR offsets
    for (int vertex = 0; vertex < vertices; vertex++) {
        std::sort(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end()
        );

        graph.offsets[vertex + 1] = graph.offsets[vertex]
            + static_cast<int64_t>(
                temporaryAdjacency[vertex].size()
                );
    }

    graph.directedEdges = graph.offsets[vertices];

    graph.adjacency.resize(
        graph.directedEdges
    );

    // Flatten temporary adjacency lists into CSR adjacency storage
    for (int vertex = 0; vertex < vertices; vertex++) {
        std::copy(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            graph.adjacency.begin() + graph.offsets[vertex]
        );
    }

    return graph;
}

CSRGraph loadGraphFromFile(
    const std::string& filename
) {
    std::ifstream input(
        filename
    );

    if (!input.is_open()) {
        std::cerr
            << "Error: cannot open "
            << filename
            << "\n";

        std::exit(
            1
        );
    }

    std::vector<std::pair<int, int>> edges;
    int maxVertexId = -1;
    std::string line;

    // Read whitespace separated undirected edge list
    while (std::getline(
        input,
        line
    )) {
        if (
            line.empty()
            || line[0] == '#'
            || line[0] == '%'
            ) {
            continue;
        }

        std::istringstream stream(
            line
        );

        int source = 0;
        int destination = 0;

        if (!(stream >> source >> destination)) {
            continue;
        }

        edges.push_back({
            source,
            destination
            });

        maxVertexId = std::max(
            maxVertexId,
            std::max(
                source,
                destination
            )
        );
    }

    const int vertices = maxVertexId + 1;

    std::vector<std::vector<int>> temporaryAdjacency(
        vertices
    );

    // Store both directions before CSR conversion
    for (const auto& edge : edges) {
        temporaryAdjacency[edge.first].push_back(
            edge.second
        );

        temporaryAdjacency[edge.second].push_back(
            edge.first
        );
    }

    CSRGraph graph;
    graph.vertices = vertices;

    graph.offsets.resize(
        vertices + 1
    );

    graph.offsets[0] = 0;

    // Remove duplicates and self loops before CSR conversion
    for (int vertex = 0; vertex < vertices; vertex++) {
        std::sort(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end()
        );

        temporaryAdjacency[vertex].erase(
            std::unique(
                temporaryAdjacency[vertex].begin(),
                temporaryAdjacency[vertex].end()
            ),
            temporaryAdjacency[vertex].end()
        );

        temporaryAdjacency[vertex].erase(
            std::remove(
                temporaryAdjacency[vertex].begin(),
                temporaryAdjacency[vertex].end(),
                vertex
            ),
            temporaryAdjacency[vertex].end()
        );

        graph.offsets[vertex + 1] = graph.offsets[vertex]
            + static_cast<int64_t>(
                temporaryAdjacency[vertex].size()
                );
    }

    graph.directedEdges = graph.offsets[vertices];

    graph.adjacency.resize(
        graph.directedEdges
    );

    // Flatten cleaned adjacency rows into CSR storage
    for (int vertex = 0; vertex < vertices; vertex++) {
        std::copy(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            graph.adjacency.begin() + graph.offsets[vertex]
        );
    }

    return graph;
}

int computeAutomorphisms(
    const MotifPattern& motif
) {
    std::vector<int> permutation(
        motif.motifK
    );

    std::iota(
        permutation.begin(),
        permutation.end(),
        0
    );

    int count = 0;

    do {
        bool valid = true;

        // Count permutations that preserve all motif edges
        for (
            int left = 0;
            left < motif.motifK && valid;
            left++
            ) {
            for (
                int right = left + 1;
                right < motif.motifK && valid;
                right++
                ) {
                if (motif.hasEdge(
                    left,
                    right
                ) != motif.hasEdge(
                    permutation[left],
                    permutation[right]
                )) {
                    valid = false;
                }
            }
        }

        if (valid) {
            count++;
        }
    } while (std::next_permutation(
        permutation.begin(),
        permutation.end()
    ));

    return count;
}

MotifPattern makeMotif(
    const std::string& name,
    const int motifVertices,
    const std::vector<std::pair<int, int>>& edges
) {
    MotifPattern motif;
    motif.name = name;
    motif.motifK = motifVertices;

    std::memset(
        motif.adjacency,
        0,
        sizeof(
            motif.adjacency
            )
    );

    // Add undirected motif edges to flattened adjacency matrix
    for (const auto& edge : edges) {
        motif.adjacency[edge.first * motifVertices + edge.second] = true;
        motif.adjacency[edge.second * motifVertices + edge.first] = true;
    }

    motif.automorphisms = computeAutomorphisms(
        motif
    );

    return motif;
}

std::unordered_map<std::string, MotifPattern> builtinMotifs() {
    std::unordered_map<std::string, MotifPattern> library;

    // Small built in motif library used by benchmark configurations
    library["triangle"] = makeMotif(
        "triangle",
        3,
        { {0, 1}, {0, 2}, {1, 2} }
    );

    library["path3"] = makeMotif(
        "path3",
        3,
        { {0, 1}, {1, 2} }
    );

    library["square"] = makeMotif(
        "square",
        4,
        { {0, 1}, {1, 2}, {2, 3}, {3, 0} }
    );

    library["4clique"] = makeMotif(
        "4clique",
        4,
        { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3} }
    );

    library["diamond"] = makeMotif(
        "diamond",
        4,
        { {0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3} }
    );

    library["tailed_triangle"] = makeMotif(
        "tailed_triangle",
        4,
        { {0, 1}, {0, 2}, {1, 2}, {2, 3} }
    );

    library["star3"] = makeMotif(
        "star3",
        4,
        { {0, 1}, {0, 2}, {0, 3} }
    );

    library["path4"] = makeMotif(
        "path4",
        4,
        { {0, 1}, {1, 2}, {2, 3} }
    );

    library["star4"] = makeMotif(
        "star4",
        5,
        { {0, 1}, {0, 2}, {0, 3}, {0, 4} }
    );

    library["5cycle"] = makeMotif(
        "5cycle",
        5,
        { {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0} }
    );

    library["5clique"] = makeMotif(
        "5clique",
        5,
        {
            {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2},
            {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}
        }
    );

    return library;
}

RequiredTable buildRequiredTable(
    const MotifPattern& motif
) {
    RequiredTable table;
    table.motifK = motif.motifK;

    // Each row stores prior motif vertices that must be adjacent
    for (int depth = 1; depth < motif.motifK; depth++) {
        int requiredCount = 0;

        for (int prior = 0; prior < depth; prior++) {
            if (motif.hasEdge(
                depth,
                prior
            )) {
                table.data[depth * HPX_MOTIF_MAX_K + requiredCount]
                    = static_cast<int8_t>(
                        prior
                        );

                requiredCount++;
            }
        }

        table.length[depth] = static_cast<int8_t>(
            requiredCount
            );
    }

    // Depth 1 connection selects faster connected search path
    table.connectedAtDepth1 = table.length[1] > 0;

    return table;
}

std::vector<int> buildChunkStarts(
    const int vertices,
    int chunkCount
) {
    if (chunkCount <= 0) {
        chunkCount = 1;
    }

    std::vector<int> starts(
        static_cast<std::size_t>(
            chunkCount
            ) + 1,
        0
    );

    const int baseChunkSize = vertices / chunkCount;
    const int extraVertices = vertices % chunkCount;
    int start = 0;

    // Split root vertices into contiguous equal count chunks
    for (int chunk = 0; chunk < chunkCount; chunk++) {
        starts[static_cast<std::size_t>(
            chunk
            )] = start;

        start += baseChunkSize + (
            chunk < extraVertices
            ? 1
            : 0
            );
    }

    starts[static_cast<std::size_t>(
        chunkCount
        )] = vertices;

    return starts;
}

double elapsedMs(
    const SearchClock::time_point start,
    const SearchClock::time_point end
) {
    return std::chrono::duration<double, std::milli>(
        end - start
    ).count();
}

int64_t searchRootsSequential(
    const CSRGraph& graph,
    const RequiredTable& requiredTable,
    const int startVertex,
    const int endVertex
) {
    int64_t sum = 0;

    // Search each root in this worker shard sequentially
    for (int vertex = startVertex; vertex < endVertex; vertex++) {
        sum += searchFromRoot(
            graph,
            requiredTable,
            vertex
        );
    }

    return sum;
}

WorkerRunResult searchLocalBlockShards(
    const CSRGraph& graph,
    const RequiredTable& requiredTable,
    const std::vector<int>& chunkStarts,
    const int firstWorker,
    const int lastWorker
) {
    WorkerRunResult result;

    const int shardCount = lastWorker - firstWorker;

    if (shardCount <= 0) {
        return result;
    }

    const int threadCount = static_cast<int>(
        std::max<std::size_t>(
            1,
            hpx::get_num_worker_threads()
        )
        );

    std::vector<int64_t> partialCounts(
        static_cast<std::size_t>(
            shardCount
            ),
        0
    );

    std::vector<double> partialTimes(
        static_cast<std::size_t>(
            shardCount
            ),
        0.0
    );

    // Search one chunk and record its local count and elapsed time
    const auto searchShard = [
        &
    ](
        const int shardIndex
        ) {
            const int worker = firstWorker + shardIndex;

            const int startVertex = chunkStarts[static_cast<std::size_t>(
                worker
                )];

            const int endVertex = chunkStarts[static_cast<std::size_t>(
                worker
                ) + 1];

            const auto searchStart = SearchClock::now();

            partialCounts[static_cast<std::size_t>(
                shardIndex
                )] = searchRootsSequential(
                    graph,
                    requiredTable,
                    startVertex,
                    endVertex
                );

            partialTimes[static_cast<std::size_t>(
                shardIndex
                )] = elapsedMs(
                    searchStart,
                    SearchClock::now()
                );
        };

    if (
        threadCount == 1
        || shardCount == 1
        ) {
        for (int shardIndex = 0; shardIndex < shardCount; shardIndex++) {
            searchShard(
                shardIndex
            );
        }
    }
    else {
        // Parallelize local shards when multiple HPX worker threads exist
        hpx::experimental::for_loop(
            hpx::execution::par,
            0,
            shardCount,
            searchShard
        );
    }

    // Combine shard counts and keep slowest shard time
    for (int shardIndex = 0; shardIndex < shardCount; shardIndex++) {
        result.count += partialCounts[static_cast<std::size_t>(
            shardIndex
            )];

        result.searchMs = std::max(
            result.searchMs,
            partialTimes[static_cast<std::size_t>(
                shardIndex
                )]
        );
    }

    return result;
}

int hpx_main(
    hpx::program_options::variables_map& variablesMap
) {
    // Only root locality performs setup and result reporting
    if (hpx::get_locality_id() != 0) {
        return hpx::finalize();
    }

    const int vertices = variablesMap["vertices"].as<int>();
    const double edgeProbability = variablesMap["edge_prob"].as<double>();
    const unsigned seed = variablesMap["seed"].as<unsigned>();

    const std::string motifName = variablesMap["motif"].as<std::string>();
    const std::string inputFile = variablesMap["input"].as<std::string>();

    std::unordered_map<std::string, MotifPattern> library = builtinMotifs();
    auto motifIterator = library.find(
        motifName
    );

    if (motifIterator == library.end()) {
        std::cerr
            << "Error: unknown motif '"
            << motifName
            << "'\n";

        for (const auto& entry : library) {
            std::cerr
                << "  "
                << entry.first
                << "\n";
        }

        return hpx::finalize();
    }

    MotifPattern& motif = motifIterator->second;

    std::cout
        << "========================================\n"
        << "  Graph Motif Search - HPX\n"
        << "========================================\n\n";

    CSRGraph graph;

    if (!inputFile.empty()) {
        std::cout
            << "Loading graph from: "
            << inputFile
            << "\n";

        graph = loadGraphFromFile(
            inputFile
        );
    }
    else if (vertices > 0) {
        std::cout
            << "Generating Erdos-Renyi G("
            << vertices
            << ", "
            << edgeProbability
            << ") seed="
            << seed
            << "\n";

        // Graph generation isn't included in search timing
        const auto generationStart = SearchClock::now();

        graph = generateErdosRenyi(
            vertices,
            edgeProbability,
            seed
        );

        const double generationMs = elapsedMs(
            generationStart,
            SearchClock::now()
        );

        std::cout
            << "Graph generation: "
            << std::fixed
            << std::setprecision(
                2
            )
            << generationMs
            << " ms\n";
    }
    else {
        std::cerr
            << "Error: specify --vertices N or --input FILE\n";

        return hpx::finalize();
    }

    // Discover localities that will receive graph and worker ranges
    std::vector<hpx::id_type> localities = hpx::find_all_localities();

    const int localityCount = static_cast<int>(
        localities.size()
        );

    const int threadOption = variablesMap[
        "threads_per_locality"
    ].as<int>();

    const std::size_t polledWorkers = hpx::get_num_worker_threads();

    int threadsPerLocality = 1;

    if (threadOption > 0) {
        threadsPerLocality = threadOption;
    }
    else if (localityCount > 0) {
        threadsPerLocality = std::max(
            1,
            static_cast<int>(
                polledWorkers / static_cast<std::size_t>(
                    localityCount
                    )
                )
        );

        if (
            localityCount > 1
            && polledWorkers % static_cast<std::size_t>(
                localityCount
                ) != 0
            ) {
            std::cerr
                << "Warning: get_num_worker_threads()="
                << polledWorkers
                << " not divisible by localities="
                << localityCount
                << "; set --threads-per-locality explicitly.\n";
        }
    }

    const int workers = localityCount * threadsPerLocality;

    // Build one root vertex chunk per logical worker
    const std::vector<int> starts = buildChunkStarts(
        graph.vertices,
        workers
    );

    std::cout
        << "  Vertices:     "
        << graph.vertices
        << "\n"
        << "  Edges:        "
        << graph.directedEdges / 2
        << "\n"
        << "  Motif:        "
        << motif.name
        << " (k="
        << motif.motifK
        << ", automorphisms="
        << motif.automorphisms
        << ")\n"
        << "  HPX threads:  "
        << threadsPerLocality
        << "\n"
        << "  Localities:   "
        << localityCount
        << "\n"
        << "  Workers:      "
        << workers
        << "\n";

    if (localityCount > 1) {
        std::cout
            << "  Polled workers: "
            << polledWorkers
            << "\n";
    }

    std::cout
        << "\n";

    int64_t totalRawEmbeddings = 0;
    double searchMs = 0.0;

    // Run empty action on all localities as a timing barrier
    const auto runBarrier = [
        &
    ](
        const std::vector<hpx::id_type>& targets
        ) {
            std::vector<hpx::future<void>> futures;

            futures.reserve(
                targets.size()
            );

            for (const hpx::id_type& locality : targets) {
                futures.push_back(
                    hpx::async<SearchBarrierAction>(
                        locality
                    )
                );
            }

            hpx::wait_all(
                futures
            );

            for (hpx::future<void>& future : futures) {
                future.get();
            }
        };

    // Collect all locality results into global count and max search time
    const auto reduceWorkerResults = [
        &
    ](
        std::vector<hpx::future<WorkerRunResult>>& futures
        ) {
            for (hpx::future<WorkerRunResult>& future : futures) {
                const WorkerRunResult run = future.get();

                totalRawEmbeddings += run.count;

                searchMs = std::max(
                    searchMs,
                    run.searchMs
                );
            }
        };

    std::cout
        << "Running motif search...\n"
        << std::flush;

    if (localityCount == 1) {
        RequiredTable requiredTable = buildRequiredTable(
            motif
        );

        runBarrier(
            localities
        );

        // Single locality can search directly without remote actions
        const WorkerRunResult run = searchLocalBlockShards(
            graph,
            requiredTable,
            starts,
            0,
            workers
        );

        totalRawEmbeddings = run.count;
        searchMs = run.searchMs;
    }
    else {
        std::vector<hpx::future<void>> initFutures;

        initFutures.reserve(
            localityCount
        );

        // Send replicated graph and chunk table to every locality
        for (
            int localityIndex = 0;
            localityIndex < localityCount;
            localityIndex++
            ) {
            initFutures.push_back(
                hpx::async<InitLocalityAction>(
                    localities[localityIndex],
                    graph,
                    motif,
                    starts
                )
            );
        }

        hpx::wait_all(
            initFutures
        );

        for (hpx::future<void>& future : initFutures) {
            future.get();
        }

        runBarrier(
            localities
        );

        std::vector<hpx::future<WorkerRunResult>> workerFutures;

        workerFutures.reserve(
            localityCount
        );

        // Launch one worker action per locality
        for (
            int localityIndex = 0;
            localityIndex < localityCount;
            localityIndex++
            ) {
            workerFutures.push_back(
                hpx::async<LocalityWorkerAction>(
                    localities[localityIndex],
                    localityIndex,
                    threadsPerLocality
                )
            );
        }

        hpx::wait_all(
            workerFutures
        );

        reduceWorkerResults(
            workerFutures
        );
    }

    // Unique instance count divides raw embeddings by automorphisms
    const double wallMs = searchMs;
    const int64_t instanceCount = totalRawEmbeddings / motif.automorphisms;

    std::cout
        << "\n-- Results ------------------------------------------\n"
        << "  Raw embeddings:   "
        << totalRawEmbeddings
        << "\n"
        << "  Automorphisms:    "
        << motif.automorphisms
        << "\n"
        << "  Instance count:   "
        << instanceCount
        << "\n"
        << "  Search time:      "
        << std::fixed
        << std::setprecision(
            3
        )
        << searchMs
        << " ms\n"
        << "  Total wall time:  "
        << wallMs
        << " ms\n"
        << "\n-- Configuration ------------------------------------\n"
        << "  Vertices:         "
        << graph.vertices
        << "\n"
        << "  Edges:            "
        << graph.directedEdges / 2
        << "\n"
        << "  Motif:            "
        << motif.name
        << " (k="
        << motif.motifK
        << ")\n"
        << "  Localities:       "
        << localityCount
        << "\n"
        << "  Threads/locality: "
        << threadsPerLocality
        << "\n"
        << "  Workers:          "
        << workers
        << "\n"
        << "========================================\n";

    return hpx::finalize();
}