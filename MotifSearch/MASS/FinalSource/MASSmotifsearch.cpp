#include "MASSmotifsearch.h"

#ifndef LOGGING
static const bool printOutput = false;
#else
static const bool printOutput = true;
#endif

void MASS_MotifSearch::printUsage(
    const char* programName
) {
    std::cerr
        << "Usage: "
        << programName
        << " user pass machinefile port nProc nThr [options]\n\n"
        << "MASS arguments:\n"
        << "  user          SSH username\n"
        << "  pass          SSH password\n"
        << "  machinefile   file with cluster hostnames\n"
        << "  port          MASS communication port\n"
        << "  nProc         number of processes\n"
        << "  nThr          threads per process\n\n"
        << "Options:\n"
        << "  --vertices N      Number of vertices\n"
        << "  --edge-prob P     Edge probability\n"
        << "  --seed S          Random seed\n"
        << "  --input FILE      Load graph from edge list file\n"
        << "  --motif NAME      Motif pattern\n"
        << "  --chunks N        Chunk count\n";
}

void MASS_MotifSearch::parseArguments(
    int argc,
    char* argv[],
    RunConfig& config
) {
    // First 6 arguments belong to MASS startup
    for (int argumentIndex = 6; argumentIndex < argc; argumentIndex++) {
        if (
            !std::strcmp(
                argv[argumentIndex],
                "--vertices"
            )
            && argumentIndex + 1 < argc
        ) {
            config.vertices = std::atoi(
                argv[++argumentIndex]
            );
        }
        else if (
            !std::strcmp(
                argv[argumentIndex],
                "--edge-prob"
            )
            && argumentIndex + 1 < argc
        ) {
            config.edgeProbability = std::atof(
                argv[++argumentIndex]
            );
        }
        else if (
            !std::strcmp(
                argv[argumentIndex],
                "--seed"
            )
            && argumentIndex + 1 < argc
        ) {
            config.seed = static_cast<unsigned>(
                std::atoi(
                    argv[++argumentIndex]
                )
            );
        }
        else if (
            !std::strcmp(
                argv[argumentIndex],
                "--input"
            )
            && argumentIndex + 1 < argc
        ) {
            config.inputFile = argv[++argumentIndex];
        }
        else if (
            !std::strcmp(
                argv[argumentIndex],
                "--motif"
            )
            && argumentIndex + 1 < argc
        ) {
            config.motifName = argv[++argumentIndex];
        }
        else if (
            !std::strcmp(
                argv[argumentIndex],
                "--chunks"
            )
            && argumentIndex + 1 < argc
        ) {
            config.chunks = std::atoi(
                argv[++argumentIndex]
            );
        }
    }
}

MASS_MotifSearch::CSRGraph MASS_MotifSearch::generateErdosRenyi(
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

    // Build an undirected graph by testing each upper triangle edge
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

    // Sort adjacency lists before flattening into CSR
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

    for (int vertex = 0; vertex < vertices; vertex++) {
        std::copy(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            graph.adjacency.begin() + graph.offsets[vertex]
        );
    }

    return graph;
}

MASS_MotifSearch::CSRGraph MASS_MotifSearch::loadGraphFromFile(
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

    for (int vertex = 0; vertex < vertices; vertex++) {
        std::copy(
            temporaryAdjacency[vertex].begin(),
            temporaryAdjacency[vertex].end(),
            graph.adjacency.begin() + graph.offsets[vertex]
        );
    }

    return graph;
}

int MASS_MotifSearch::computeAutomorphisms(
    const MotifPattern& motif
) {
    std::vector<int> permutation(
        motif.motifK
    );

    for (int index = 0; index < motif.motifK; index++) {
        permutation[index] = index;
    }

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

MASS_MotifSearch::MotifPattern MASS_MotifSearch::makeMotif(
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

    for (const auto& edge : edges) {
        motif.adjacency[edge.first * motifVertices + edge.second] = true;
        motif.adjacency[edge.second * motifVertices + edge.first] = true;
    }

    motif.automorphisms = computeAutomorphisms(
        motif
    );

    return motif;
}

std::unordered_map<std::string, MASS_MotifSearch::MotifPattern>
MASS_MotifSearch::builtinMotifs() {
    std::unordered_map<std::string, MotifPattern> library;

    library["triangle"] = makeMotif(
        "triangle",
        3,
        {{0, 1}, {0, 2}, {1, 2}}
    );

    library["path3"] = makeMotif(
        "path3",
        3,
        {{0, 1}, {1, 2}}
    );

    library["square"] = makeMotif(
        "square",
        4,
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}}
    );

    library["4clique"] = makeMotif(
        "4clique",
        4,
        {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}
    );

    library["diamond"] = makeMotif(
        "diamond",
        4,
        {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}}
    );

    library["tailed_triangle"] = makeMotif(
        "tailed_triangle",
        4,
        {{0, 1}, {0, 2}, {1, 2}, {2, 3}}
    );

    library["star3"] = makeMotif(
        "star3",
        4,
        {{0, 1}, {0, 2}, {0, 3}}
    );

    library["path4"] = makeMotif(
        "path4",
        4,
        {{0, 1}, {1, 2}, {2, 3}}
    );

    library["star4"] = makeMotif(
        "star4",
        5,
        {{0, 1}, {0, 2}, {0, 3}, {0, 4}}
    );

    library["5cycle"] = makeMotif(
        "5cycle",
        5,
        {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}
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

bool MASS_MotifSearch::MotifPattern::hasEdge(
    const int leftMotifVertex,
    const int rightMotifVertex
) const {
    return adjacency[
        leftMotifVertex * motifK + rightMotifVertex
    ];
}

std::vector<int> MASS_MotifSearch::buildChunkStarts(
    const CSRGraph& graph,
    int chunkCount
) {
    if (chunkCount <= 0) {
        chunkCount = 1;
    }

    std::vector<int> starts(
        chunkCount + 1,
        0
    );

    const int baseChunkSize = graph.vertices / chunkCount;
    const int extraVertices = graph.vertices % chunkCount;
    int start = 0;

    // Split root vertices into contiguous equal count chunks
    for (int chunk = 0; chunk < chunkCount; chunk++) {
        starts[chunk] = start;
        start += baseChunkSize + (
            chunk < extraVertices
                ? 1
                : 0
        );
    }

    starts[chunkCount] = graph.vertices;

    return starts;
}
std::vector<char> MASS_MotifSearch::serializeBenchmarkData(
    const CSRGraph& graph,
    const MotifPattern& motif,
    const std::vector<int>& chunkStarts
) {
    // Last chunk start is sentinel end index
    const int chunkCount = static_cast<int>(
        chunkStarts.size()
        ) - 1;

    // Header stores fixed config and motif metadata
    const size_t headerBytes = sizeof(
        MotifSearchConfig
        );

    // Chunk table maps each Place to its vertex range
    const size_t chunkBytes = chunkStarts.size() * sizeof(
        int
        );

    // CSR offsets need 1 extra entry for end offset
    const size_t offsetBytes = (graph.vertices + 1) * sizeof(
        int64_t
        );

    // Directed adjacency list stores both directions
    const size_t adjacencyBytes = graph.directedEdges * sizeof(
        int
        );

    // Single buffer keeps MASS initialization argument simple
    const size_t totalBytes = headerBytes
        + chunkBytes
        + offsetBytes
        + adjacencyBytes;

    std::vector<char> buffer(
        totalBytes
    );

    // Pointer walks through serialized sections in order
    char* pointer = buffer.data();

    MotifSearchConfig config;

    // Clear padding and string fields before copying values
    std::memset(
        &config,
        0,
        sizeof(
            config
            )
    );

    // Store graph and motif sizes for each Place
    config.n = graph.vertices;
    config.m = graph.directedEdges;
    config.motifK = motif.motifK;
    config.motifAutomorphisms = motif.automorphisms;
    config.chunkCount = chunkCount;

    // Copy fixed size motif adjacency matrix into config
    std::memcpy(
        config.motifAdjacency,
        motif.adjacency,
        sizeof(
            config.motifAdjacency
            )
    );

    // Copy motif name and leave final byte available for null end
    std::strncpy(
        config.motifName,
        motif.name.c_str(),
        sizeof(
            config.motifName
            ) - 1
    );

    // Write config header first so Places can parse buffer layout
    std::memcpy(
        pointer,
        &config,
        headerBytes
    );

    pointer += headerBytes;

    // Write chunk start table after fixed header
    std::memcpy(
        pointer,
        chunkStarts.data(),
        chunkBytes
    );

    pointer += chunkBytes;

    // Write CSR offset table after chunk starts
    std::memcpy(
        pointer,
        graph.offsets.data(),
        offsetBytes
    );

    pointer += offsetBytes;

    // Write CSR adjacency data last
    std::memcpy(
        pointer,
        graph.adjacency.data(),
        adjacencyBytes
    );

    return buffer;
}

int MASS_MotifSearch::run(
    int argc,
    char* argv[]
) {
    // MASS needs 4 runtime arguments before benchmark arguments
    if (argc < 7) {
        printUsage(
            argv[0]
        );

        return 1;
    }

    // Pass only MASS runtime arguments into MASS initialization
    char* massArguments[4] = {
        argv[1],
        argv[2],
        argv[3],
        argv[4]
    };

    // Benchmark process count is used for chunk planning and printing
    const int processes = std::atoi(
        argv[5]
    );

    // Thread count helps choose default chunk count
    const int threads = std::atoi(
        argv[6]
    );

    RunConfig config;

    // Parse graph, motif, seed, and chunk options
    parseArguments(
        argc,
        argv,
        config
    );

    // Require either generated graph parameters or input graph file
    if (
        config.inputFile.empty()
        && config.vertices <= 0
        ) {
        std::cerr
            << "Error: specify --vertices N --edge-prob P or --input FILE\n";

        return 1;
    }

    // Load built in motif patterns by name
    std::unordered_map<std::string, MotifPattern> library = builtinMotifs();
    auto motifIterator = library.find(
        config.motifName
    );

    // Stop early when motif name isn't known
    if (motifIterator == library.end()) {
        std::cerr
            << "Error: unknown motif '"
            << config.motifName
            << "'\n";

        // Print supported motif names to make command fix obvious
        for (const auto& entry : library) {
            std::cerr
                << "  "
                << entry.first
                << "\n";
        }

        return 1;
    }

    // Use selected motif for search and normalization
    MotifPattern& motif = motifIterator->second;

    std::cout
        << "========================================\n"
        << "  Graph Motif Search - MASS\n"
        << "========================================\n\n";

    CSRGraph graph;

    if (!config.inputFile.empty()) {
        std::cout
            << "Loading graph from: "
            << config.inputFile
            << "\n";

        // Load graph once before broadcasting to Places
        graph = loadGraphFromFile(
            config.inputFile
        );
    }
    else {
        std::cout
            << "Generating Erdos-Renyi G("
            << config.vertices
            << ", "
            << config.edgeProbability
            << ") seed="
            << config.seed
            << "\n";

        // Graph generation isn't part of search timing
        const auto generationStart = std::chrono::high_resolution_clock::now();

        graph = generateErdosRenyi(
            config.vertices,
            config.edgeProbability,
            config.seed
        );

        const auto generationEnd = std::chrono::high_resolution_clock::now();

        const double generationMs = std::chrono::duration<double, std::milli>(
            generationEnd - generationStart
        ).count();

        std::cout
            << "Graph generation: "
            << std::fixed
            << std::setprecision(
                2
            )
            << generationMs
            << " ms\n";
    }

    // Default to 1 chunk per process thread when chunks aren't specified
    int chunkCount = config.chunks > 0
        ? config.chunks
        : processes * std::max(
            1,
            threads
        );

    // Do not create more chunks than vertices
    if (chunkCount > graph.vertices) {
        chunkCount = graph.vertices;
    }

    // Keep at least 1 chunk for valid Place creation
    if (chunkCount < 1) {
        chunkCount = 1;
    }

    // Partition vertices into graph chunks with similar work
    const std::vector<int> chunkStarts = buildChunkStarts(
        graph,
        chunkCount
    );

    // Graph stores undirected edges as 2 directed adjacency entries
    const int64_t undirectedEdges = graph.directedEdges / 2;

    std::cout
        << "  Vertices:     "
        << graph.vertices
        << "\n"
        << "  Edges:        "
        << undirectedEdges
        << "\n"
        << "  Motif:        "
        << motif.name
        << " (k="
        << motif.motifK
        << ", automorphisms="
        << motif.automorphisms
        << ")\n"
        << "  Chunks:       "
        << chunkCount
        << " (processes="
        << processes
        << " threads="
        << threads
        << ")\n\n";

    // Serialize graph and motif once for MASS Place initialization
    std::vector<char> buffer = serializeBenchmarkData(
        graph,
        motif,
        chunkStarts
    );

    const int bufferSize = static_cast<int>(
        buffer.size()
        );

    std::cout
        << "Serialized buffer: "
        << bufferSize
        << " bytes ("
        << std::fixed
        << std::setprecision(
            2
        )
        << bufferSize / 1048576.0
        << " MB)\n\n";

    if (printOutput) {
        std::cerr
            << "Initializing MASS...\n";
    }

    // Start MASS runtime after graph preparation is complete
    MASS::init(
        massArguments,
        processes,
        threads
    );

    std::cout
        << "Creating "
        << chunkCount
        << " Places ("
        << processes
        << " processes x "
        << threads
        << " threads)...\n";

    // Create 1 Place per graph chunk
    Places* graphPlaces = new Places(
        MASS_MOTIF_PLACES_HANDLE,
        "MotifPlace",
        0,
        static_cast<void*>(
            nullptr
            ),
        0,
        1,
        chunkCount
    );

    std::cout
        << "Broadcasting graph and chunk table...\n";

    // Initialize every Place with same serialized graph buffer
    graphPlaces->callAll(
        "MotifPlace::init",
        buffer.data(),
        bufferSize
    );

    std::cout
        << "Running motif search...\n"
        << std::flush;

    // Each Place searches only vertices assigned to its chunk
    graphPlaces->callAll(
        "MotifPlace::search"
    );

    int64_t totalRawEmbeddings = 0;
    double searchMs = 0.0;

    try {
        // Collect raw embedding count from every Place
        int64_t* counts = reinterpret_cast<int64_t*>(
            graphPlaces->callAll(
                "MotifPlace::getCount",
                static_cast<void**>(
                    nullptr
                    ),
                0,
                sizeof(
                    int64_t
                    )
            )
            );

        if (counts != nullptr) {
            // Sum all Place counts into global raw embedding count
            for (int chunk = 0; chunk < chunkCount; chunk++) {
                totalRawEmbeddings += counts[chunk];
            }
        }
    }
    catch (...) {
        std::cerr
            << "Warning: count collection failed\n";
    }

    try {
        // Collect measured search time from every Place
        double* searchTimes = reinterpret_cast<double*>(
            graphPlaces->callAll(
                "MotifPlace::getSearchMs",
                static_cast<void**>(
                    nullptr
                    ),
                0,
                sizeof(
                    double
                    )
            )
            );

        if (searchTimes != nullptr) {
            // Search time is slowest Place time
            for (int chunk = 0; chunk < chunkCount; chunk++) {
                searchMs = std::max(
                    searchMs,
                    searchTimes[chunk]
                );
            }
        }
    }
    catch (...) {
        std::cerr
            << "Warning: search time collection failed\n";
    }

    // MASS search wall time is represented by slowest Place search
    const double wallMs = searchMs;

    // Divide by automorphisms to report unique motif instances
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
        << std::fixed
        << std::setprecision(
            3
        )
        << wallMs
        << " ms\n";

    std::cout
        << "\n-- Configuration ------------------------------------\n"
        << "  Vertices:         "
        << graph.vertices
        << "\n"
        << "  Edges:            "
        << undirectedEdges
        << "\n"
        << "  Motif:            "
        << motif.name
        << " (k="
        << motif.motifK
        << ")\n"
        << "  Processes:        "
        << processes
        << "\n"
        << "  Threads/process:  "
        << threads
        << "\n"
        << "  Chunks:           "
        << chunkCount
        << "\n"
        << "========================================\n";

    // Shut down MASS runtime after Place results are collected
    MASS::finish();

    return 0;
}