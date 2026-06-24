#include "pm2dgemmCANNONS.h"

// Split 1 global dimension into balanced chunks
pair<int, int> PM2_DGEMM::split1D(
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
double PM2_DGEMM::valueAt(
    const int rowIndex,
    const int colIndex,
    const uint64_t seed
) {
    uint64_t hashValue = seed;

    hashValue ^= uint64_t(
        rowIndex
    ) * ROW_HASH_MULTIPLIER;

    hashValue ^= uint64_t(
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

// Validate inputs before starting benchmark work
bool PM2_DGEMM::inputValidation(
    const int rank,
    const short runs,
    const int rowsOfAAndC,
    const int sharedDimension,
    const int colsOfBAndC,
    const int processGridRows,
    const int processGridCols
) {
    const bool invalidInput =
        runs <= 0
        || rowsOfAAndC <= 0
        || sharedDimension <= 0
        || colsOfBAndC <= 0
        || processGridRows <= 0
        || processGridCols <= 0;

    if (rank == 0 && invalidInput) {
        cerr
            << "ERROR: Inputs out of bounds, ensure all inputs "
            << "are greater than 0\n";

        return false;
    }

    return true;
}

// Run Cannon's DGEMM benchmark and report timing/checksum data
void PM2_DGEMM::runBenchmarkCannons(
    const short runs,
    const int rowsOfAAndC,
    const int sharedDimension,
    const int colsOfBAndC,
    const uint64_t baseSeed,
    const int processGridRows,
    const int processGridCols
) {
    int rank = 0;

    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &rank
    );

    // Stop early if the config isn't valid
    if (!inputValidation(
        rank,
        runs,
        rowsOfAAndC,
        sharedDimension,
        colsOfBAndC,
        processGridRows,
        processGridCols
    )) {
        exit(
            1
        );
    }

    int processes = 0;

    MPI_Comm_size(
        MPI_COMM_WORLD,
        &processes
    );

    // Validate that process grid dimensions match total process count
    const bool processGridMismatch =
        processGridRows * processGridCols != processes;

    if (processGridMismatch) {
        if (rank == 0) {
            cerr
                << "ERROR: Process grid dimensions must match "
                << "MPI rank count\n";
        }

        exit(
            1
        );
    }

    // Cannon's algorithm needs a square process grid
    if (processGridRows != processGridCols) {
        if (rank == 0) {
            cerr
                << "ERROR: Cannon's algorithm requires a square "
                << "process grid (Px == Py), got "
                << processGridRows
                << " x "
                << processGridCols
                << "\n";
        }

        exit(
            1
        );
    }

    // Store square grid size once because Cannon's requires Px == Py
    const int processGridSize = processGridRows;
    const int myProcessRow = rank / processGridSize;
    const int myProcessCol = rank % processGridSize;

    // Use different seed from A to reduce A/B value correlation
    const uint64_t bMatrixSeed = baseSeed ^ B_MATRIX_SEED_XOR;

    // Get local A row range for this process row
    const pair<int, int> aRowRange = split1D(
        rowsOfAAndC,
        processGridSize,
        myProcessRow
    );

    // Get local B column range for this process column
    const pair<int, int> bColRange = split1D(
        colsOfBAndC,
        processGridSize,
        myProcessCol
    );

    // Store local block sizes for loop bounds and buffer sizes
    const int myARows = aRowRange.second;
    const int myBCols = bColRange.second;

    // Store each shared dimension block range for reuse
    vector<int> kStarts(
        processGridSize
    );

    vector<int> kWidths(
        processGridSize
    );

    int maxKWidth = 0;

    for (
        int blockIndex = 0;
        blockIndex < processGridSize;
        blockIndex++
    ) {
        // Get shared dimension block assigned to this Cannon's step
        const pair<int, int> kBlock = split1D(
            sharedDimension,
            processGridSize,
            blockIndex
        );

        // Save start and width so each run can reuse block info
        kStarts[blockIndex] = kBlock.first;
        kWidths[blockIndex] = kBlock.second;

        // Track largest k block width for shift buffers with constant size
        if (kBlock.second > maxKWidth) {
            maxKWidth = kBlock.second;
        }
    }

    // Use padded A/B panels so shifted buffers have fixed message sizes
    vector<double> aPanel(
        myARows * maxKWidth
    );

    vector<double> bPanel(
        maxKWidth * myBCols
    );

    vector<double> aRecv(
        myARows * maxKWidth
    );

    vector<double> bRecv(
        maxKWidth * myBCols
    );

    vector<double> cLocal(
        myARows * myBCols,
        0.0
    );

    MPI_Comm rowComm;
    MPI_Comm colComm;

    // Row communicator shifts A panels horizontally
    MPI_Comm_split(
        MPI_COMM_WORLD,
        myProcessRow,
        myProcessCol,
        &rowComm
    );

    // Column communicator shifts B panels vertically
    MPI_Comm_split(
        MPI_COMM_WORLD,
        myProcessCol,
        myProcessRow,
        &colComm
    );

    // A panels shift left across each process row
    const int aSendTo = (
        myProcessCol - 1 + processGridSize
    ) % processGridSize;

    // A panels are received from right neighbor
    const int aRecvFrom = (
        myProcessCol + 1
    ) % processGridSize;

    // B panels shift upward across each process column
    const int bSendTo = (
        myProcessRow - 1 + processGridSize
    ) % processGridSize;

    // B panels are received from lower neighbor
    const int bRecvFrom = (
        myProcessRow + 1
    ) % processGridSize;

    // Repeat full benchmark for given run count
    for (int run = 1; run <= runs; run++) {
        // Reset local C block before each run
        fill(
            cLocal.begin(),
            cLocal.end(),
            0.0
        );

        // Track local timing before reducing to max rank time
        double localFillSeconds = 0.0;
        double localComputeSeconds = 0.0;
        double localShiftWaitSeconds = 0.0;

        // Synchronize ranks for fair timing
        MPI_Barrier(
            MPI_COMM_WORLD
        );

        // Start relevant clocks
        const auto totalStart =
            chrono::high_resolution_clock::now();
        const auto fillStart =
            chrono::high_resolution_clock::now();

        // Initial skew chooses which k block this process starts with
        const int initialBlock = (
            myProcessRow + myProcessCol
        ) % processGridSize;

        // Initial shared dimension block range for A and B panels
        const int initialKStart = kStarts[initialBlock];
        const int initialKWidth = kWidths[initialBlock];

        // Build initial skewed A panel
        for (int localRow = 0; localRow < myARows; localRow++) {
            // Convert local A row index to global matrix row
            const int globalRow = aRowRange.first + localRow;

            // Fill only active k width
            // Remaining padded entries aren't used
            for (
                int kOffset = 0;
                kOffset < initialKWidth;
                kOffset++
            ) {
                aPanel[localRow * maxKWidth + kOffset] = valueAt(
                    globalRow,
                    initialKStart + kOffset,
                    baseSeed
                );
            }
        }

        // Build initial skewed B panel
        for (
            int kOffset = 0;
            kOffset < initialKWidth;
            kOffset++
        ) {
            // Convert local k offset to global shared dimension index
            const int globalK = initialKStart + kOffset;

            // Fill B values for this process column's owned columns
            for (int localCol = 0; localCol < myBCols; localCol++) {
                bPanel[kOffset * myBCols + localCol] = valueAt(
                    globalK,
                    bColRange.first + localCol,
                    bMatrixSeed
                );
            }
        }

        // Stop fill timer
        const auto fillEnd =
            chrono::high_resolution_clock::now();

        localFillSeconds = chrono::duration<double>(
            fillEnd - fillStart
        ).count();

        // Each step computes 1 local k block contribution
        for (
            int cannonStep = 0;
            cannonStep < processGridSize;
            cannonStep++
        ) {
            const int currentBlock = (
                myProcessRow + myProcessCol + cannonStep
            ) % processGridSize;

            // Only active width participates in multiplication
            const int currentKWidth = kWidths[currentBlock];

            MPI_Request requests[4];

            // Last step doesn't shift because it isn't needed
            const bool shifting = cannonStep < processGridSize - 1;

            if (shifting) {
                // Start A send before compute to overlap communication
                MPI_Isend(
                    aPanel.data(),
                    myARows * maxKWidth,
                    MPI_DOUBLE,
                    aSendTo,
                    0,
                    rowComm,
                    &requests[0]
                );

                // Receive next A panel from row neighbor
                MPI_Irecv(
                    aRecv.data(),
                    myARows * maxKWidth,
                    MPI_DOUBLE,
                    aRecvFrom,
                    0,
                    rowComm,
                    &requests[1]
                );

                // Start B send before compute to overlap communication
                MPI_Isend(
                    bPanel.data(),
                    maxKWidth * myBCols,
                    MPI_DOUBLE,
                    bSendTo,
                    0,
                    colComm,
                    &requests[2]
                );

                // Receive next B panel from column neighbor
                MPI_Irecv(
                    bRecv.data(),
                    maxKWidth * myBCols,
                    MPI_DOUBLE,
                    bRecvFrom,
                    0,
                    colComm,
                    &requests[3]
                );
            }

            // Start compute timer after shift requests are posted
            const auto computeStart =
                chrono::high_resolution_clock::now();

            int localRow = 0;

            // Reuse each B row across 4 C rows
            for (
                ;
                localRow + MICRO_KERNEL_ROWS - 1 < myARows;
                localRow += MICRO_KERNEL_ROWS
            ) {
                // Get 4 output rows from local C block
                double* cRow0 =
                    &cLocal[(localRow + 0) * myBCols];

                double* cRow1 =
                    &cLocal[(localRow + 1) * myBCols];

                double* cRow2 =
                    &cLocal[(localRow + 2) * myBCols];

                double* cRow3 =
                    &cLocal[(localRow + 3) * myBCols];

                // Get 4 input rows from local A panel
                const double* aRow0 =
                    &aPanel[(localRow + 0) * maxKWidth];

                const double* aRow1 =
                    &aPanel[(localRow + 1) * maxKWidth];

                const double* aRow2 =
                    &aPanel[(localRow + 2) * maxKWidth];

                const double* aRow3 =
                    &aPanel[(localRow + 3) * maxKWidth];

                // Walk through current shared dimension block
                for (
                    int kOffset = 0;
                    kOffset < currentKWidth;
                    kOffset++
                ) {
                    // Load A scalars once for reuse across local B columns
                    const double aValue0 = aRow0[kOffset];
                    const double aValue1 = aRow1[kOffset];
                    const double aValue2 = aRow2[kOffset];
                    const double aValue3 = aRow3[kOffset];

                    // Current B row is reused by all 4 C rows
                    const double* bRow =
                        &bPanel[kOffset * myBCols];

                    // Update 1 column position in 4 C rows
                    for (
                        int localCol = 0;
                        localCol < myBCols;
                        localCol++
                    ) {
                        const double bValue = bRow[localCol];

                        cRow0[localCol] += aValue0 * bValue;
                        cRow1[localCol] += aValue1 * bValue;
                        cRow2[localCol] += aValue2 * bValue;
                        cRow3[localCol] += aValue3 * bValue;
                    }
                }
            }

            // Handle remaining rows that don't fill a 4 row block
            for (; localRow < myARows; localRow++) {
                // Get output row for tail computation
                double* cRow =
                    &cLocal[localRow * myBCols];

                // Get matching A row for scalar tail computation
                const double* aRow =
                    &aPanel[localRow * maxKWidth];

                // Accumulate current k block into 1 C row
                for (
                    int kOffset = 0;
                    kOffset < currentKWidth;
                    kOffset++
                ) {
                    const double aValue = aRow[kOffset];

                    // Current B row contributes to all local C columns
                    const double* bRow =
                        &bPanel[kOffset * myBCols];

                    for (
                        int localCol = 0;
                        localCol < myBCols;
                        localCol++
                    ) {
                        cRow[localCol] +=
                            aValue * bRow[localCol];
                    }
                }
            }

            // Stop compute timer before waiting for unfinished shifts
            const auto computeEnd =
                chrono::high_resolution_clock::now();

            localComputeSeconds += chrono::duration<double>(
                computeEnd - computeStart
            ).count();

            if (shifting) {
                // Wait timer captures only communication left after compute
                const auto waitStart =
                    chrono::high_resolution_clock::now();

                MPI_Waitall(
                    4,
                    requests,
                    MPI_STATUSES_IGNORE
                );

                const auto waitEnd =
                    chrono::high_resolution_clock::now();

                localShiftWaitSeconds += chrono::duration<double>(
                    waitEnd - waitStart
                ).count();

                // Newly received panels become active for next step
                swap(
                    aPanel,
                    aRecv
                );

                swap(
                    bPanel,
                    bRecv
                );
            }
        }

        // End total timer after everything is done
        const auto totalEnd =
            chrono::high_resolution_clock::now();

        double localSum = 0.0;

        // Reduce local C block into 1 checksum value
        for (const double value : cLocal) {
            localSum += value;
        }

        double globalSum = 0.0;

        // Combine local checksums on rank 0
        MPI_Reduce(
            &localSum,
            &globalSum,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0,
            MPI_COMM_WORLD
        );

        double maxFillSeconds = 0.0;
        double maxComputeSeconds = 0.0;
        double maxShiftWaitSeconds = 0.0;
        double maxTotalSeconds = 0.0;

        // Local total is reduced with max to represent parallel wall time
        const double localTotalSeconds = chrono::duration<double>(
            totalEnd - totalStart
        ).count();

        // Report max rank time so output reflects parallel bottleneck time
        MPI_Allreduce(
            &localFillSeconds,
            &maxFillSeconds,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );

        // Compute time uses slowest rank for fairness
        MPI_Allreduce(
            &localComputeSeconds,
            &maxComputeSeconds,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );

        // Shift wait time shows uncovered communication after overlap
        MPI_Allreduce(
            &localShiftWaitSeconds,
            &maxShiftWaitSeconds,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );

        // Total time uses slowest rank's timed phase
        MPI_Allreduce(
            &localTotalSeconds,
            &maxTotalSeconds,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );

        if (rank == 0) {
            if (run == 1) {
                // Print configuration once per benchmark config
                cout
                    << "--------------NEW CONFIG-----------------\n";

                cout
                    << "RUNNING CANNON'S BENCHMARK "
                    << "FOR CONFIGURATION:\n"
                    << "RANKS: "
                    << processes
                    << "\nA (row) size: "
                    << rowsOfAAndC
                    << "\nShared dim size: "
                    << sharedDimension
                    << "\nB (col) size: "
                    << colsOfBAndC
                    << "\nRuns: "
                    << runs
                    << '\n'
                    << "MPI Process Grid: "
                    << processGridSize
                    << " x "
                    << processGridSize
                    << '\n'
                    << "Cannon's Steps: "
                    << processGridSize
                    << '\n'
                    << "Reported comparison metric: "
                    << "Compute Time only\n"
                    << "MPI overheads reported separately: "
                    << "Fill Time, Shift Wait Time\n"
                    << "Shift Wait Time = time blocked in "
                    << "MPI_Waitall after compute "
                    << "(0 = fully hidden)\n";
            }

            // Print 1 timing block per completed run
            cout
                << "RUN: "
                << run
                << '\n';

            cout
                << "Checksum: "
                << fixed
                << setprecision(
                    15
                )
                << globalSum
                << "\n";

            cout
                << "Compute Time: "
                << maxComputeSeconds
                << "s\n";

            cout
                << "Shift Wait Time: "
                << maxShiftWaitSeconds
                << "s\n";

            cout
                << "Fill Time: "
                << maxFillSeconds
                << "s\n";

            cout
                << "Total MPI Timed Phase: "
                << maxTotalSeconds
                << "s\n";

            cout
                << "<=========================================================>"
                << '\n';
        }
    }

    MPI_Comm_free(
        &rowComm
    );

    MPI_Comm_free(
        &colComm
    );
}
