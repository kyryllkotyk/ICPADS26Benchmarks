/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: DGEMM (Cannon's Algorithm)                                      *
 * Library: HPX                                                               *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "HPXdgemmCannon.h"

// Split one global dimension into balanced chunks
std::pair<int, int> split1D(
    const int totalDimensionSize,
    const int processesInDimension,
    const int processIndex
) {
    const int baseSize = totalDimensionSize / processesInDimension;
    const int extraCells = totalDimensionSize % processesInDimension;

    const int localCount = baseSize + (
        processIndex < extraCells
            ? 1
            : 0
    );

    const int localStart = processIndex * baseSize + (
        processIndex < extraCells
            ? processIndex
            : extraCells
    );

    return {
        localStart,
        localCount
    };
}

// Generate matrix values without storing full matrices
double valueAt(
    const int rowIndex,
    const int colIndex,
    const std::uint64_t seed
) {
    std::uint64_t hashValue = seed;

    hashValue ^= static_cast<std::uint64_t>(
        rowIndex
    ) * ROW_HASH_MULTIPLIER;

    hashValue ^= static_cast<std::uint64_t>(
        colIndex
    ) * COL_HASH_MULTIPLIER;

    hashValue ^= (
        hashValue >> 30
    );

    hashValue *= FIRST_MIX_MULTIPLIER;

    hashValue ^= (
        hashValue >> 27
    );

    hashValue *= SECOND_MIX_MULTIPLIER;

    hashValue ^= (
        hashValue >> 31
    );

    return (
        hashValue >> 11
    ) * DOUBLE_SCALE;
}

// Convert command line initialization option to enum
InitType parseInitType(
    const std::string& raw
) {
    if (raw == "ones") {
        return InitType::ONES;
    }

    if (raw == "identity") {
        return InitType::IDENTITY;
    }

    if (raw == "deterministic") {
        return InitType::DETERMINISTIC;
    }

    return InitType::HASH;
}

HPXDgemmCannonTile::HPXDgemmCannonTile(
    const int rowsOfAAndCInput,
    const int colsOfBAndCInput,
    const int sharedDimensionInput,
    const int processGridSizeInput,
    const int processRowInput,
    const int processColInput,
    const std::uint64_t seed,
    const std::string& init
) :
    rowsOfAAndC(
        rowsOfAAndCInput
    ),
    colsOfBAndC(
        colsOfBAndCInput
    ),
    sharedDimension(
        sharedDimensionInput
    ),
    processGridSize(
        processGridSizeInput
    ),
    processRow(
        processRowInput
    ),
    processCol(
        processColInput
    ),
    baseSeed(
        seed
    )
{
    // Use separate seed stream for B values
    bMatrixSeed = baseSeed ^ B_MATRIX_SEED_XOR;

    // Cache selected generation mode once per tile
    initType = parseInitType(
        init
    );

    // Determine local output row range from tile row
    const std::pair<int, int> rowRange = split1D(
        rowsOfAAndC,
        processGridSize,
        processRow
    );

    localRowStart = rowRange.first;
    localRowCount = rowRange.second;

    // Determine local output column range from tile column
    const std::pair<int, int> colRange = split1D(
        colsOfBAndC,
        processGridSize,
        processCol
    );

    localColStart = colRange.first;
    localColCount = colRange.second;

    int maxSharedBlockWidth = 0;

    // Find largest shared block so panel storage fits every step
    for (int step = 0; step < processGridSize; step++) {
        const int kWidth = split1D(
            sharedDimension,
            processGridSize,
            step
        ).second;

        maxSharedBlockWidth = std::max(
            maxSharedBlockWidth,
            kWidth
        );
    }

    maxKWidth = maxSharedBlockWidth;

    // A panel holds local rows for 1 shared dimension block
    aPanel.assign(
        static_cast<std::size_t>(
            localRowCount
        ) * maxKWidth,
        0.0
    );

    // B panel holds 1 shared dimension block for local columns
    bPanel.assign(
        static_cast<std::size_t>(
            maxKWidth
        ) * localColCount,
        0.0
    );

    // Receive buffers stage incoming panels before swapping
    aReceive.assign(
        aPanel.size(),
        0.0
    );

    bReceive.assign(
        bPanel.size(),
        0.0
    );

    // C tile stores final local output block
    cTile.assign(
        static_cast<std::size_t>(
            localRowCount
        ) * localColCount,
        0.0
    );
}

// Generate one A value for selected initialization mode
double HPXDgemmCannonTile::genA(
    const int globalRow,
    const int globalCol
) const {
    switch (initType) {
    case InitType::ONES:
        return 1.0;

    case InitType::IDENTITY:
        return globalRow == globalCol
            ? 1.0
            : 0.0;

    case InitType::DETERMINISTIC:
        return static_cast<double>(
            globalRow * sharedDimension + globalCol
        ) / static_cast<double>(
            rowsOfAAndC * sharedDimension
        );

    default:
        return valueAt(
            globalRow,
            globalCol,
            baseSeed
        );
    }
}

// Generate one B value for selected initialization mode
double HPXDgemmCannonTile::genB(
    const int globalRow,
    const int globalCol
) const {
    switch (initType) {
    case InitType::ONES:
        return 1.0;

    case InitType::IDENTITY:
        return globalRow == globalCol
            ? 1.0
            : 0.0;

    case InitType::DETERMINISTIC:
        return static_cast<double>(
            globalRow * colsOfBAndC + globalCol
        ) / static_cast<double>(
            sharedDimension * colsOfBAndC
        );

    default:
        return valueAt(
            globalRow,
            globalCol,
            bMatrixSeed
        );
    }
}

// Generate initial skewed A and B panels for Cannon step 0
void HPXDgemmCannonTile::initialFill()
{
    const int initialBlock = (
        processRow + processCol
    ) % processGridSize;

    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        processGridSize,
        initialBlock
    );

    const int kStart = kRange.first;
    const int kWidth = kRange.second;

    // Clear panels before generating initial skewed blocks
    std::fill(
        aPanel.begin(),
        aPanel.end(),
        0.0
    );

    std::fill(
        bPanel.begin(),
        bPanel.end(),
        0.0
    );

    // Generate local A rows for initial Cannon skew
    for (int localRow = 0; localRow < localRowCount; localRow++) {
        const int globalRow = localRowStart + localRow;

        double* const row = &aPanel[
            static_cast<std::size_t>(
                localRow
            ) * maxKWidth
        ];

        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            row[kOffset] = genA(
                globalRow,
                kStart + kOffset
            );
        }
    }

    // Generate local B columns for initial Cannon skew
    for (int kOffset = 0; kOffset < kWidth; kOffset++) {
        const int globalK = kStart + kOffset;

        double* const row = &bPanel[
            static_cast<std::size_t>(
                kOffset
            ) * localColCount
        ];

        for (int localCol = 0; localCol < localColCount; localCol++) {
            row[localCol] = genB(
                globalK,
                localColStart + localCol
            );
        }
    }

    currentStep = 0;
}

// Multiply current local A and B panels into C tile
void HPXDgemmCannonTile::multiplyCurrent(
    const int step
) {
    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        processGridSize,
        (
            processRow + processCol + step
        ) % processGridSize
    );

    const int currentKWidth = kRange.second;
    const int innerCol = localColCount;

    double* __restrict__ cBase = cTile.data();
    const double* __restrict__ aBase = aPanel.data();
    const double* __restrict__ bBase = bPanel.data();

    int localRow = 0;

    // Reuse each B row across 4 C rows
    for (
        ;
        localRow + MICRO_KERNEL_ROWS - 1 < localRowCount;
        localRow += MICRO_KERNEL_ROWS
    ) {
        // Get 4 output rows from local C tile
        double* __restrict__ cRow0 = cBase +
            static_cast<std::size_t>(
                localRow + 0
            ) * innerCol;

        double* __restrict__ cRow1 = cBase +
            static_cast<std::size_t>(
                localRow + 1
            ) * innerCol;

        double* __restrict__ cRow2 = cBase +
            static_cast<std::size_t>(
                localRow + 2
            ) * innerCol;

        double* __restrict__ cRow3 = cBase +
            static_cast<std::size_t>(
                localRow + 3
            ) * innerCol;

        // Get 4 input rows from local A panel
        const double* __restrict__ aRow0 = aBase +
            static_cast<std::size_t>(
                localRow + 0
            ) * maxKWidth;

        const double* __restrict__ aRow1 = aBase +
            static_cast<std::size_t>(
                localRow + 1
            ) * maxKWidth;

        const double* __restrict__ aRow2 = aBase +
            static_cast<std::size_t>(
                localRow + 2
            ) * maxKWidth;

        const double* __restrict__ aRow3 = aBase +
            static_cast<std::size_t>(
                localRow + 3
            ) * maxKWidth;

        // Walk through current shared dimension block
        for (int kOffset = 0; kOffset < currentKWidth; kOffset++) {
            // Load A scalars once for reuse across local columns
            const double aValue0 = aRow0[kOffset];
            const double aValue1 = aRow1[kOffset];
            const double aValue2 = aRow2[kOffset];
            const double aValue3 = aRow3[kOffset];

            // Current B row is reused by all 4 C rows
            const double* __restrict__ bRow = bBase +
                static_cast<std::size_t>(
                    kOffset
                ) * innerCol;

            // Update 1 column position in 4 C rows
            for (int localCol = 0; localCol < innerCol; localCol++) {
                const double bValue = bRow[localCol];

                cRow0[localCol] += aValue0 * bValue;
                cRow1[localCol] += aValue1 * bValue;
                cRow2[localCol] += aValue2 * bValue;
                cRow3[localCol] += aValue3 * bValue;
            }
        }
    }

    // Handle remaining rows that do not fill a 4 row block
    for (; localRow < localRowCount; localRow++) {
        // Get output row for tail computation
        double* __restrict__ cRow = cBase +
            static_cast<std::size_t>(
                localRow
            ) * innerCol;

        // Get matching A row for scalar tail computation
        const double* __restrict__ aRow = aBase +
            static_cast<std::size_t>(
                localRow
            ) * maxKWidth;

        // Accumulate current k block into 1 C row
        for (int kOffset = 0; kOffset < currentKWidth; kOffset++) {
            const double aValue = aRow[kOffset];

            // Current B row contributes to all local C columns
            const double* __restrict__ bRow = bBase +
                static_cast<std::size_t>(
                    kOffset
                ) * innerCol;

            for (int localCol = 0; localCol < innerCol; localCol++) {
                cRow[localCol] += aValue * bRow[localCol];
            }
        }
    }
}

void HPXDgemmCannonTile::resetC()
{
    // Clear output tile before each measured run
    std::fill(
        cTile.begin(),
        cTile.end(),
        0.0
    );

    currentStep = 0;
}

double HPXDgemmCannonTile::checksum() const
{
    double sum = 0.0;

    // Sum local C tile for final validation value
    for (const double value : cTile) {
        sum += value;
    }

    return sum;
}

// Locality state used by HPX plain actions
static std::unique_ptr<HPXDgemmCannonTile> globalTile;
static std::vector<hpx::id_type> globalLocalities;

void initTileImpl(
    const int rowsOfAAndC,
    const int colsOfBAndC,
    const int sharedDimension,
    const int processGridSize,
    const int processRow,
    const int processCol,
    const std::uint64_t seed,
    const std::string init
) {
    // Create 1 Cannon tile on each HPX locality
    globalTile = std::make_unique<HPXDgemmCannonTile>(
        rowsOfAAndC,
        colsOfBAndC,
        sharedDimension,
        processGridSize,
        processRow,
        processCol,
        seed,
        init
    );
}
HPX_PLAIN_ACTION(
    initTileImpl,
    initTileAction
);

void installLocalityTableImpl(
    std::vector<hpx::id_type> localities
) {
    // Cache locality table for later neighbor shift actions
    globalLocalities = std::move(
        localities
    );
}
HPX_PLAIN_ACTION(
    installLocalityTableImpl,
    installLocalityTableAction
);

void receiveAImpl(
    std::vector<double> data
) {
    if (globalTile) {
        // Stage received A panel until all shifts finish
        globalTile->aReceive = std::move(
            data
        );
    }
}
HPX_PLAIN_ACTION(
    receiveAImpl,
    receiveAAction
);

void receiveBImpl(
    std::vector<double> data
) {
    if (globalTile) {
        // Stage received B panel until all shifts finish
        globalTile->bReceive = std::move(
            data
        );
    }
}
HPX_PLAIN_ACTION(
    receiveBImpl,
    receiveBAction
);

void initialFillImpl()
{
    if (globalTile) {
        // Build initial skewed A and B panels inside timed phase
        globalTile->initialFill();
    }
}
HPX_PLAIN_ACTION(
    initialFillImpl,
    initialFillAction
);

void multiplyStep0Impl()
{
    if (!globalTile) {
        return;
    }

    // Step 0 uses panels from initial skew without shifting first
    globalTile->multiplyCurrent(
        0
    );

    globalTile->currentStep = 1;
}
HPX_PLAIN_ACTION(
    multiplyStep0Impl,
    multiplyStep0Action
);

void shiftPanelsImpl()
{
    if (!globalTile) {
        return;
    }

    HPXDgemmCannonTile& tile = *globalTile;

    const int leftId = tile.processRow * tile.processGridSize + (
        (
            tile.processCol - 1 + tile.processGridSize
        ) % tile.processGridSize
    );

    const int upId = (
        (
            tile.processRow - 1 + tile.processGridSize
        ) % tile.processGridSize
    ) * tile.processGridSize + tile.processCol;

    std::vector<hpx::future<void>> futures;

    futures.reserve(
        2
    );

    // Send A panel to left neighbor for Cannon shift
    futures.push_back(
        hpx::async<receiveAAction>(
            globalLocalities[static_cast<std::size_t>(
                leftId
            )],
            tile.aPanel
        )
    );

    // Send B panel to upper neighbor for Cannon shift
    futures.push_back(
        hpx::async<receiveBAction>(
            globalLocalities[static_cast<std::size_t>(
                upId
            )],
            tile.bPanel
        )
    );

    hpx::wait_all(
        futures
    );
}
HPX_PLAIN_ACTION(
    shiftPanelsImpl,
    shiftPanelsAction
);

void shiftFinishImpl()
{
    if (!globalTile) {
        return;
    }

    // Swap staged A panel into active panel storage
    std::swap(
        globalTile->aPanel,
        globalTile->aReceive
    );

    // Swap staged B panel into active panel storage
    std::swap(
        globalTile->bPanel,
        globalTile->bReceive
    );
}
HPX_PLAIN_ACTION(
    shiftFinishImpl,
    shiftFinishAction
);

void multiplyOnlyImpl()
{
    if (!globalTile) {
        return;
    }

    // Multiply current shifted panels for this locality
    globalTile->multiplyCurrent(
        globalTile->currentStep
    );

    globalTile->currentStep = (
        globalTile->currentStep + 1
    ) % globalTile->processGridSize;
}
HPX_PLAIN_ACTION(
    multiplyOnlyImpl,
    multiplyOnlyAction
);

void resetCImpl()
{
    if (globalTile) {
        // Reset output tile outside measured timed phase
        globalTile->resetC();
    }
}
HPX_PLAIN_ACTION(
    resetCImpl,
    resetCAction
);

// Empty action used to synchronize localities before timing
void timedPhaseBarrierImpl()
{
}
HPX_PLAIN_ACTION(
    timedPhaseBarrierImpl,
    timedPhaseBarrierAction
);

// Return checksum for this locality tile
double getChecksumImpl()
{
    return globalTile
        ? globalTile->checksum()
        : 0.0;
}
HPX_PLAIN_ACTION(
    getChecksumImpl,
    getChecksumAction
);

static void waitForExpectedLocalities(
    std::vector<hpx::id_type>& foundLocalities,
    const int expectedLocalities
) {
    // Give HPX runtime time to discover all requested localities
    for (int attempt = 0; attempt < 300; attempt++) {
        foundLocalities = hpx::find_all_localities();

        if (static_cast<int>(
            foundLocalities.size()
        ) >= expectedLocalities) {
            break;
        }

        std::this_thread::sleep_for(
            std::chrono::milliseconds(
                100
            )
        );
    }
}

static bool buildLocalityTable(
    const std::vector<hpx::id_type>& foundLocalities,
    std::vector<hpx::id_type>& localities,
    const int expectedLocalities
) {
    std::vector<char> slotOk(
        static_cast<std::size_t>(
            expectedLocalities
        ),
        0
    );

    // Place each locality into slot matching its HPX locality id
    for (const hpx::id_type& locality : foundLocalities) {
        const std::uint32_t rawId =
            hpx::naming::get_locality_id_from_id(
                locality
            );

        if (rawId >= static_cast<std::uint32_t>(
            expectedLocalities
        )) {
            continue;
        }

        const int localityId = static_cast<int>(
            rawId
        );

        if (slotOk[static_cast<std::size_t>(
            localityId
        )]) {
            std::cerr
                << "Error: duplicate HPX locality id "
                << localityId
                << "\n";

            return false;
        }

        slotOk[static_cast<std::size_t>(
            localityId
        )] = 1;

        localities[static_cast<std::size_t>(
            localityId
        )] = locality;
    }

    // Require contiguous locality ids so Cannon grid mapping is valid
    for (int localityId = 0; localityId < expectedLocalities; localityId++) {
        if (!slotOk[static_cast<std::size_t>(
            localityId
        )]) {
            std::cerr
                << "Error: missing HPX locality id "
                << localityId
                << " (need contiguous ids 0.."
                << expectedLocalities - 1
                << " among "
                << foundLocalities.size()
                << " reported)\n";

            return false;
        }
    }

    return true;
}

static void initializeTiles(
    const std::vector<hpx::id_type>& localities,
    const int rowsOfAAndC,
    const int colsOfBAndC,
    const int sharedDimension,
    const int processGridSize,
    const std::uint64_t seed,
    const std::string& init
) {
    std::vector<hpx::future<void>> futures;

    futures.reserve(
        localities.size()
    );

    // Create matching Cannon tile on every locality
    for (int localityId = 0;
         localityId < static_cast<int>(
             localities.size()
         );
         localityId++) {
        futures.push_back(
            hpx::async<initTileAction>(
                localities[static_cast<std::size_t>(
                    localityId
                )],
                rowsOfAAndC,
                colsOfBAndC,
                sharedDimension,
                processGridSize,
                localityId / processGridSize,
                localityId % processGridSize,
                seed,
                init
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

static void installLocalityTables(
    const std::vector<hpx::id_type>& localities
) {
    std::vector<hpx::future<void>> futures;

    futures.reserve(
        localities.size()
    );

    // Share complete locality table with every tile
    for (const hpx::id_type& locality : localities) {
        futures.push_back(
            hpx::async<installLocalityTableAction>(
                locality,
                localities
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

// Launch same action on every locality and wait for completion
template <typename ActionType>
static void runActionOnAllLocalities(
    const std::vector<hpx::id_type>& localities
) {
    std::vector<hpx::future<void>> futures;

    futures.reserve(
        localities.size()
    );

    for (const hpx::id_type& locality : localities) {
        futures.push_back(
            hpx::async<ActionType>(
                locality
            )
        );
    }

    hpx::wait_all(
        futures
    );
}

static double collectChecksum(
    const std::vector<hpx::id_type>& localities
) {
    std::vector<hpx::future<double>> futures;

    futures.reserve(
        localities.size()
    );

    // Request each local C tile checksum
    for (const hpx::id_type& locality : localities) {
        futures.push_back(
            hpx::async<getChecksumAction>(
                locality
            )
        );
    }

    double checksum = 0.0;

    // Combine local checksums on current locality
    for (hpx::future<double>& future : futures) {
        checksum += future.get();
    }

    return checksum;
}

int hpx_main(
    hpx::program_options::variables_map& variables
) {
    // Read matrix, process grid, run, and initialization options
    const int rowsOfAAndC = variables["M"].as<int>();
    const int colsOfBAndC = variables["N"].as<int>();
    const int sharedDimension = variables["K"].as<int>();

    int processGridSize = variables["P"].as<int>();

    // Allow Px/Py spelling as an alias for square Cannon grid size
    if (processGridSize <= 0) {
        const int processGridRows = variables["Px"].as<int>();
        const int processGridCols = variables["Py"].as<int>();

        processGridSize = processGridRows > 0
            && processGridRows == processGridCols
            ? processGridRows
            : 1;
    }

    const int runs = variables["runs"].as<int>();
    const std::uint64_t seed = variables["seed"].as<std::uint64_t>();
    const std::string init = variables["init_type"].as<std::string>();

    // Cannon grid maps one tile to each HPX locality
    const int expectedLocalities = processGridSize * processGridSize;

    std::vector<hpx::id_type> foundLocalities;

    // Wait for HPX runtime to report all expected localities
    waitForExpectedLocalities(
        foundLocalities,
        expectedLocalities
    );

    if (static_cast<int>(
        foundLocalities.size()
    ) != expectedLocalities) {
        std::cerr
            << "Error: HPX locality count must exactly match P*P. Found "
            << foundLocalities.size()
            << ", expected "
            << expectedLocalities
            << "\n";

        return hpx::finalize();
    }

    // Reorder locality ids into Cannon grid order
    std::vector<hpx::id_type> localities(
        static_cast<std::size_t>(
            expectedLocalities
        )
    );

    if (!buildLocalityTable(
        foundLocalities,
        localities,
        expectedLocalities
    )) {
        return hpx::finalize();
    }

    const int localityCount = expectedLocalities;

    std::cout
        << "\n--------------DGEMM HPX "
        << "(Cannon, neighbor shifts + installed localities)----\n"
        << "Matrix: "
        << rowsOfAAndC
        << "x"
        << colsOfBAndC
        << " = "
        << rowsOfAAndC
        << "x"
        << sharedDimension
        << " * "
        << sharedDimension
        << "x"
        << colsOfBAndC
        << "\n"
        << "Grid: "
        << processGridSize
        << "x"
        << processGridSize
        << "  Localities: "
        << localityCount
        << "  Steps: "
        << processGridSize
        << "  Init: "
        << init
        << "  Runs: "
        << runs
        << "\n";

    // Create one local Cannon tile on each locality
    initializeTiles(
        localities,
        rowsOfAAndC,
        colsOfBAndC,
        sharedDimension,
        processGridSize,
        seed,
        init
    );

    // Install full locality table so each tile can contact neighbors
    installLocalityTables(
        localities
    );

    for (int run = 1; run <= runs; run++) {
        // Each run repeats same Cannon pipeline from reset state
        // Clear output before starting timed phase
        runActionOnAllLocalities<resetCAction>(
            localities
        );

        // Make sure every locality is ready before timing starts
        runActionOnAllLocalities<timedPhaseBarrierAction>(
            localities
        );

        // Timed phase includes initial fill, shifts, and multiply work
        const auto startTime = std::chrono::high_resolution_clock::now();

        // Initial skew generation is part of measured Cannon work
        runActionOnAllLocalities<initialFillAction>(
            localities
        );

        // First multiply uses initially skewed panels
        runActionOnAllLocalities<multiplyStep0Action>(
            localities
        );

        for (int step = 1; step < processGridSize; step++) {
            // Send panels to Cannon neighbors for next step
            runActionOnAllLocalities<shiftPanelsAction>(
                localities
            );

            // Activate received panels after all sends finish
            runActionOnAllLocalities<shiftFinishAction>(
                localities
            );

            // Multiply current shifted panels
            runActionOnAllLocalities<multiplyOnlyAction>(
                localities
            );
        }

        const auto endTime = std::chrono::high_resolution_clock::now();

        // Collect checksum after final Cannon step
        const double checksum = collectChecksum(
            localities
        );

        const double elapsedSeconds = std::chrono::duration<double>(
            endTime - startTime
        ).count();

        // Print shared DGEMM run format
        dgemm_report::printRunResults(
            run,
            checksum,
            elapsedSeconds,
            rowsOfAAndC,
            colsOfBAndC,
            sharedDimension
        );
    }

    return hpx::finalize();
}
