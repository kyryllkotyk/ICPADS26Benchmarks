#include "HPXdgemmSUMMA.h"

// Compute greatest common divisor for SUMMA step schedule
int gcdInt(
    int leftValue,
    int rightValue
) {
    while (rightValue != 0) {
        const int remainder = leftValue % rightValue;

        leftValue = rightValue;
        rightValue = remainder;
    }

    return leftValue;
}

// Split 1 global dimension into balanced chunks
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

// Convert CLI init type into compact local enum
InitType parseInitType(
    const std::string& rawInput
) {
    if (rawInput == "ones") {
        return InitType::ONES;
    }

    if (rawInput == "identity") {
        return InitType::IDENTITY;
    }

    if (rawInput == "deterministic") {
        return InitType::DETERMINISTIC;
    }

    return InitType::HASH;
}

LocalTile::LocalTile(
    const int rowsOfAAndCValue,
    const int colsOfBAndCValue,
    const int sharedDimensionValue,
    const int processGridRowsValue,
    const int processGridColsValue,
    const int tileRowValue,
    const int tileColValue,
    const std::uint64_t baseSeedValue,
    const std::string& initTypeName
) : rowsOfAAndC(
    rowsOfAAndCValue
),
colsOfBAndC(
    colsOfBAndCValue
),
sharedDimension(
    sharedDimensionValue
),
processGridRows(
    processGridRowsValue
),
processGridCols(
    processGridColsValue
),
tileRow(
    tileRowValue
),
tileCol(
    tileColValue
),
localRowStart(
    0
),
localRowCount(
    0
),
localColStart(
    0
),
localColCount(
    0
),
stepCount(
    0
),
maxKWidth(
    0
),
baseSeed(
    baseSeedValue
),
bMatrixSeed(
    baseSeedValue^ B_MATRIX_SEED_XOR
),
initType(
    parseInitType(
        initTypeName
    )
) {
    // Each tile owns 1 row chunk and 1 column chunk of C
    const std::pair<int, int> rowRange = split1D(
        rowsOfAAndC,
        processGridRows,
        tileRow
    );

    localRowStart = rowRange.first;
    localRowCount = rowRange.second;

    const std::pair<int, int> colRange = split1D(
        colsOfBAndC,
        processGridCols,
        tileCol
    );

    localColStart = colRange.first;
    localColCount = colRange.second;

    // Step count covers all row and column broadcast roots
    stepCount = (
        processGridRows / gcdInt(
            processGridRows,
            processGridCols
        )
        ) * processGridCols;

    // Largest k block sets fixed panel buffer capacity
    maxKWidth = (
        sharedDimension + stepCount - 1
        ) / stepCount;

    // Store only active panels and local C tile
    aPanel.assign(
        static_cast<std::size_t>(
            localRowCount
            ) * maxKWidth,
        0.0
    );

    bPanel.assign(
        static_cast<std::size_t>(
            maxKWidth
            ) * localColCount,
        0.0
    );

    cTile.assign(
        static_cast<std::size_t>(
            localRowCount
            ) * localColCount,
        0.0
    );
}

double LocalTile::genA(
    const int globalRow,
    const int globalK
) const {
    // Init mode controls whether A is generated analytically or by hash
    switch (initType) {
    case InitType::ONES:
        return 1.0;

    case InitType::IDENTITY:
        return globalRow == globalK
            ? 1.0
            : 0.0;

    case InitType::DETERMINISTIC:
        return static_cast<double>(
            globalRow * sharedDimension + globalK
            ) / static_cast<double>(
                rowsOfAAndC * sharedDimension
                );

    default:
        return valueAt(
            globalRow,
            globalK,
            baseSeed
        );
    }
}

double LocalTile::genB(
    const int globalK,
    const int globalCol
) const {
    // B uses matching init mode but different hash seed when needed
    switch (initType) {
    case InitType::ONES:
        return 1.0;

    case InitType::IDENTITY:
        return globalK == globalCol
            ? 1.0
            : 0.0;

    case InitType::DETERMINISTIC:
        return static_cast<double>(
            globalK * colsOfBAndC + globalCol
            ) / static_cast<double>(
                sharedDimension * colsOfBAndC
                );

    default:
        return valueAt(
            globalK,
            globalCol,
            bMatrixSeed
        );
    }
}

void LocalTile::fillAPanelForStep(
    const int step,
    std::vector<double>& outputPanel
) const {
    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        stepCount,
        step
    );

    const int kStart = kRange.first;
    const int kWidth = kRange.second;

    // Owner builds compact A panel for this shared dimension block
    outputPanel.assign(
        static_cast<std::size_t>(
            localRowCount
            ) * kWidth,
        0.0
    );

    for (int localRow = 0; localRow < localRowCount; localRow++) {
        const int globalRow = localRowStart + localRow;

        double* const outputRow = &outputPanel[
            static_cast<std::size_t>(
                localRow
                ) * kWidth
        ];

        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            outputRow[kOffset] = genA(
                globalRow,
                kStart + kOffset
            );
        }
    }
}

void LocalTile::fillBPanelForStep(
    const int step,
    std::vector<double>& outputPanel
) const {
    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        stepCount,
        step
    );

    const int kStart = kRange.first;
    const int kWidth = kRange.second;

    // Owner builds compact B panel for this shared dimension block
    outputPanel.assign(
        static_cast<std::size_t>(
            kWidth
            ) * localColCount,
        0.0
    );

    for (int kOffset = 0; kOffset < kWidth; kOffset++) {
        const int globalK = kStart + kOffset;

        double* const outputRow = &outputPanel[
            static_cast<std::size_t>(
                kOffset
                ) * localColCount
        ];

        for (int localCol = 0; localCol < localColCount; localCol++) {
            outputRow[localCol] = genB(
                globalK,
                localColStart + localCol
            );
        }
    }
}

void LocalTile::storeAFromPayload(
    const int step,
    const std::vector<double>& payload
) {
    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        stepCount,
        step
    );

    const int kWidth = kRange.second;

    // Clear padded storage before installing compact received panel
    std::fill(
        aPanel.begin(),
        aPanel.end(),
        0.0
    );

    const std::size_t rowBytes = static_cast<std::size_t>(
        kWidth
        ) * sizeof(
            double
            );

    // Copy compact row data into fixed width local panel storage
    for (int localRow = 0; localRow < localRowCount; localRow++) {
        std::memcpy(
            &aPanel[
                static_cast<std::size_t>(
                    localRow
                    ) * maxKWidth
            ],
            &payload[
                static_cast<std::size_t>(
                    localRow
                    ) * kWidth
            ],
            rowBytes
        );
    }
}

void LocalTile::storeBFromPayload(
    const int step,
    const std::vector<double>& payload
) {
    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        stepCount,
        step
    );

    const int kWidth = kRange.second;

    // Clear padded storage before installing compact received panel
    std::fill(
        bPanel.begin(),
        bPanel.end(),
        0.0
    );

    // B panel width is localColCount, so rows copy directly
    for (int kOffset = 0; kOffset < kWidth; kOffset++) {
        std::memcpy(
            &bPanel[
                static_cast<std::size_t>(
                    kOffset
                    ) * localColCount
            ],
            &payload[
                static_cast<std::size_t>(
                    kOffset
                    ) * localColCount
            ],
            static_cast<std::size_t>(
                localColCount
                ) * sizeof(
                    double
                    )
        );
    }
}

void LocalTile::accumulatePanel(
    const int step
) {
    const std::pair<int, int> kRange = split1D(
        sharedDimension,
        stepCount,
        step
    );

    const int kWidth = kRange.second;
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
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
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

        // Get matching A row for tail computation
        const double* __restrict__ aRow = aBase +
            static_cast<std::size_t>(
                localRow
                ) * maxKWidth;

        // Accumulate current k block into 1 C row
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
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

void LocalTile::resetC() {
    std::fill(
        cTile.begin(),
        cTile.end(),
        0.0
    );
}

double LocalTile::checksum() const {
    double localSum = 0.0;

    // Sum local C tile for global validation
    for (const double value : cTile) {
        localSum += value;
    }

    return localSum;
}

// One tile and locality lookup table are stored on each locality
static std::unique_ptr<LocalTile> localTile;
static std::vector<hpx::id_type> localityTable;

void initTileImpl(
    const int rowsOfAAndC,
    const int colsOfBAndC,
    const int sharedDimension,
    const int processGridRows,
    const int processGridCols,
    const int tileRow,
    const int tileCol,
    const std::uint64_t seed,
    const std::string initTypeName
) {
    // Create 1 tile object on each locality
    localTile = std::make_unique<LocalTile>(
        rowsOfAAndC,
        colsOfBAndC,
        sharedDimension,
        processGridRows,
        processGridCols,
        tileRow,
        tileCol,
        seed,
        initTypeName
    );
}

HPX_PLAIN_ACTION(
    initTileImpl,
    initTileAction
);

void installLocalityTableImpl(
    std::vector<hpx::id_type> localities
) {
    // Install shared locality lookup table on each locality
    localityTable = std::move(
        localities
    );
}

HPX_PLAIN_ACTION(
    installLocalityTableImpl,
    installLocalityTableAction
);

void receiveAImpl(
    const int step,
    const std::vector<double> payload
) {
    if (localTile) {
        // Receiver installs row broadcast panel for this step
        localTile->storeAFromPayload(
            step,
            payload
        );
    }
}

HPX_PLAIN_ACTION(
    receiveAImpl,
    receiveAAction
);

void receiveBImpl(
    const int step,
    const std::vector<double> payload
) {
    if (localTile) {
        // Receiver installs column broadcast panel for this step
        localTile->storeBFromPayload(
            step,
            payload
        );
    }
}

HPX_PLAIN_ACTION(
    receiveBImpl,
    receiveBAction
);

void broadcastPanelsImpl(
    const int step
) {
    if (!localTile) {
        return;
    }

    LocalTile& tile = *localTile;

    const int aOwnerCol = step % tile.processGridCols;
    const int bOwnerRow = step % tile.processGridRows;

    std::vector<hpx::future<void>> sendFutures;

    // A owner pushes panel across its process row
    if (tile.tileCol == aOwnerCol) {
        std::vector<double> aPayload;

        tile.fillAPanelForStep(
            step,
            aPayload
        );

        // Owner also installs its own local A panel
        tile.storeAFromPayload(
            step,
            aPayload
        );

        for (
            int otherCol = 0;
            otherCol < tile.processGridCols;
            otherCol++
            ) {
            if (otherCol == tile.tileCol) {
                continue;
            }

            const int peerId = tile.tileRow * tile.processGridCols + otherCol;

            sendFutures.push_back(
                hpx::async<receiveAAction>(
                    localityTable[peerId],
                    step,
                    aPayload
                )
            );
        }
    }

    // B owner pushes panel across its process column
    if (tile.tileRow == bOwnerRow) {
        std::vector<double> bPayload;

        tile.fillBPanelForStep(
            step,
            bPayload
        );

        // Owner also installs its own local B panel
        tile.storeBFromPayload(
            step,
            bPayload
        );

        for (
            int otherRow = 0;
            otherRow < tile.processGridRows;
            otherRow++
            ) {
            if (otherRow == tile.tileRow) {
                continue;
            }

            const int peerId = otherRow * tile.processGridCols + tile.tileCol;

            sendFutures.push_back(
                hpx::async<receiveBAction>(
                    localityTable[peerId],
                    step,
                    bPayload
                )
            );
        }
    }

    if (!sendFutures.empty()) {
        // Owners wait until remote receives have installed panels
        hpx::wait_all(
            sendFutures
        );
    }
}

HPX_PLAIN_ACTION(
    broadcastPanelsImpl,
    broadcastPanelsAction
);

void accumulatePanelImpl(
    const int step
) {
    if (localTile) {
        localTile->accumulatePanel(
            step
        );
    }
}

HPX_PLAIN_ACTION(
    accumulatePanelImpl,
    accumulatePanelAction
);

void resetCImpl() {
    if (localTile) {
        localTile->resetC();
    }
}

HPX_PLAIN_ACTION(
    resetCImpl,
    resetCAction
);

void timedPhaseBarrierImpl() {
    // Empty remote action acts as a pre timing synchronization point
}

HPX_PLAIN_ACTION(
    timedPhaseBarrierImpl,
    timedPhaseBarrierAction
);

double getChecksumImpl() {
    if (localTile) {
        return localTile->checksum();
    }

    return 0.0;
}

HPX_PLAIN_ACTION(
    getChecksumImpl,
    getChecksumAction
);

int hpx_main(
    hpx::program_options::variables_map& variablesMap
) {
    // Read benchmark dimensions and runtime options
    const int rowsOfAAndC = variablesMap["M"].as<int>();
    const int colsOfBAndC = variablesMap["N"].as<int>();
    const int sharedDimension = variablesMap["K"].as<int>();
    const int processGridRows = variablesMap["Px"].as<int>();
    const int processGridCols = variablesMap["Py"].as<int>();
    const int runs = variablesMap["runs"].as<int>();
    const std::uint64_t seed = variablesMap[
        "seed"
    ].as<std::uint64_t>();
    const std::string initTypeName = variablesMap[
        "init_type"
    ].as<std::string>();

    const int expectedLocalities = processGridRows * processGridCols;

    std::vector<hpx::id_type> foundLocalities;

    // Wait for expected HPX localities to become visible
    for (
        int attempt = 0;
        attempt < LOCALITY_DISCOVERY_ATTEMPTS;
        attempt++
        ) {
        foundLocalities = hpx::find_all_localities();

        if (static_cast<int>(
            foundLocalities.size()
            ) == expectedLocalities) {
            break;
        }

        std::this_thread::sleep_for(
            std::chrono::milliseconds(
                LOCALITY_DISCOVERY_SLEEP_MS
            )
        );
    }

    if (static_cast<int>(
        foundLocalities.size()
        ) != expectedLocalities) {
        std::cerr
            << "Error: HPX locality count must exactly match Px*Py. Found "
            << foundLocalities.size()
            << ", expected "
            << expectedLocalities
            << std::endl;

        return hpx::finalize();
    }

    std::vector<hpx::id_type> localities(
        static_cast<std::size_t>(
            expectedLocalities
            )
    );

    std::vector<char> slotFilled(
        static_cast<std::size_t>(
            expectedLocalities
            ),
        0
    );

    // Map HPX locality id to geometric tile id
    for (const hpx::id_type& locality : foundLocalities) {
        const std::uint32_t localityId =
            hpx::naming::get_locality_id_from_id(
                locality
            );

        if (localityId >= static_cast<std::uint32_t>(
            expectedLocalities
            )) {
            continue;
        }

        const int geometricId = static_cast<int>(
            localityId
            );

        if (slotFilled[
            static_cast<std::size_t>(
                geometricId
                )
        ]) {
            std::cerr
                << "Error: duplicate HPX locality id "
                << geometricId
                << std::endl;

            return hpx::finalize();
        }

        slotFilled[
            static_cast<std::size_t>(
                geometricId
                )
        ] = 1;

        localities[
            static_cast<std::size_t>(
                geometricId
                )
        ] = locality;
    }

    // Require contiguous locality ids so grid coordinates are unambiguous
    for (int geometricId = 0; geometricId < expectedLocalities; geometricId++) {
        if (!slotFilled[
            static_cast<std::size_t>(
                geometricId
                )
        ]) {
            std::cerr
                << "Error: missing HPX locality id "
                << geometricId
                << " among "
                << foundLocalities.size()
                << " reported localities\n";

            return hpx::finalize();
        }
    }

    const int stepCount = (
        processGridRows / gcdInt(
            processGridRows,
            processGridCols
        )
        ) * processGridCols;

    std::cout
        << "\n--------------DGEMM HPX SUMMA---------------------\n"
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
        << "\nGrid: "
        << processGridRows
        << "x"
        << processGridCols
        << "  Localities: "
        << expectedLocalities
        << "  Steps: "
        << stepCount
        << "  Init: "
        << initTypeName
        << "  Runs: "
        << runs
        << "\n";

    // Create tile state on each HPX locality
    {
        std::vector<hpx::future<void>> futures;

        futures.reserve(
            static_cast<std::size_t>(
                expectedLocalities
                )
        );

        for (int id = 0; id < expectedLocalities; id++) {
            futures.push_back(
                hpx::async<initTileAction>(
                    localities[id],
                    rowsOfAAndC,
                    colsOfBAndC,
                    sharedDimension,
                    processGridRows,
                    processGridCols,
                    id / processGridCols,
                    id % processGridCols,
                    seed,
                    initTypeName
                )
            );
        }

        hpx::wait_all(
            futures
        );
    }

    // Install locality table so tile actions can push panels to peers
    {
        std::vector<hpx::future<void>> futures;

        futures.reserve(
            static_cast<std::size_t>(
                expectedLocalities
                )
        );

        for (int id = 0; id < expectedLocalities; id++) {
            futures.push_back(
                hpx::async<installLocalityTableAction>(
                    localities[id],
                    localities
                )
            );
        }

        hpx::wait_all(
            futures
        );
    }

    for (int run = 1; run <= runs; run++) {
        // Reset local C before timed phase
        {
            std::vector<hpx::future<void>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            for (int id = 0; id < expectedLocalities; id++) {
                futures.push_back(
                    hpx::async<resetCAction>(
                        localities[id]
                    )
                );
            }

            hpx::wait_all(
                futures
            );
        }

        // Synchronize localities before total timer starts
        {
            std::vector<hpx::future<void>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            for (int id = 0; id < expectedLocalities; id++) {
                futures.push_back(
                    hpx::async<timedPhaseBarrierAction>(
                        localities[id]
                    )
                );
            }

            hpx::wait_all(
                futures
            );
        }

        const auto startTime = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < stepCount; step++) {
            // Step phase 1 broadcasts A and B panels
            {
                std::vector<hpx::future<void>> futures;

                futures.reserve(
                    static_cast<std::size_t>(
                        expectedLocalities
                        )
                );

                for (int id = 0; id < expectedLocalities; id++) {
                    futures.push_back(
                        hpx::async<broadcastPanelsAction>(
                            localities[id],
                            step
                        )
                    );
                }

                hpx::wait_all(
                    futures
                );
            }

            // Step phase 2 multiplies received panels into local C
            {
                std::vector<hpx::future<void>> futures;

                futures.reserve(
                    static_cast<std::size_t>(
                        expectedLocalities
                        )
                );

                for (int id = 0; id < expectedLocalities; id++) {
                    futures.push_back(
                        hpx::async<accumulatePanelAction>(
                            localities[id],
                            step
                        )
                    );
                }

                hpx::wait_all(
                    futures
                );
            }
        }

        const auto endTime = std::chrono::high_resolution_clock::now();

        double globalChecksum = 0.0;

        // Pull local checksums and combine them on the root locality
        {
            std::vector<hpx::future<double>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            for (int id = 0; id < expectedLocalities; id++) {
                futures.push_back(
                    hpx::async<getChecksumAction>(
                        localities[id]
                    )
                );
            }

            for (hpx::future<double>& future : futures) {
                globalChecksum += future.get();
            }
        }

        const double elapsedSeconds = std::chrono::duration<double>(
            endTime - startTime
        ).count();

        dgemm_report::printRunResults(
            run,
            globalChecksum,
            elapsedSeconds,
            rowsOfAAndC,
            colsOfBAndC,
            sharedDimension
        );
    }

    return hpx::finalize();
}