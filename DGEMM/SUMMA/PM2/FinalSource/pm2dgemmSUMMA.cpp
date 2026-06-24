#include "pm2dgemmSUMMA.h"

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

// Run SUMMA benchmark
void PM2_DGEMM::runBenchmark(
    const short runs,
    const int rowsOfAAndC,
    const int sharedDimension,
    const int colsOfBAndC,
    const uint64_t baseSeed,
    const int processGridRows,
    const int processGridCols
) {
    // MPI rank for this process
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

    // Total MPI process count
    int processes = 0;

    MPI_Comm_size(
        MPI_COMM_WORLD,
        &processes
    );

    // Validate that process grid matches total processes
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

    // Compute greatest common divisor for SUMMA step count
    const auto gcdInt = [](
        int leftValue,
        int rightValue
        ) {
            while (rightValue != 0) {
                const int remainder = leftValue % rightValue;

                leftValue = rightValue;
                rightValue = remainder;
            }

            return leftValue;
        };

    // Step count covers all row and column broadcast roots
    const int stepCount = (
        processGridRows / gcdInt(
            processGridRows,
            processGridCols
        )
        ) * processGridCols;

    const int myProcessRow = rank / processGridCols;
    const int myProcessCol = rank % processGridCols;

    // Use a different B seed from A
    const uint64_t bMatrixSeed = baseSeed ^ B_MATRIX_SEED_XOR;

    // Get local A row range for this process row
    const pair<int, int> aRowRange = split1D(
        rowsOfAAndC,
        processGridRows,
        myProcessRow
    );

    // Get local B column range for this process column
    const pair<int, int> bColRange = split1D(
        colsOfBAndC,
        processGridCols,
        myProcessCol
    );

    // Store local block sizes for loop bounds and buffer sizes
    const int myARows = aRowRange.second;
    const int myBCols = bColRange.second;

    // Store each shared dimension block range for reuse
    vector<int> kStarts(
        stepCount
    );

    vector<int> kWidths(
        stepCount
    );

    int maxKWidth = 0;

    for (
        int step = 0;
        step < stepCount;
        step++
        ) {
        // Get shared dimension block for this SUMMA step
        const pair<int, int> kBlock = split1D(
            sharedDimension,
            stepCount,
            step
        );

        // Save start and width so each run can reuse block info
        kStarts[step] = kBlock.first;
        kWidths[step] = kBlock.second;

        // Track largest k block width for panel buffer sizing
        if (kBlock.second > maxKWidth) {
            maxKWidth = kBlock.second;
        }
    }

    // Use 2 A panel buffers to overlap next broadcast with compute
    vector<double> aPanels[2] = {
        vector<double>(
            myARows * maxKWidth
        ),
        vector<double>(
            myARows * maxKWidth
        )
    };

    // Use 2 B panel buffers to overlap next broadcast with compute
    vector<double> bPanels[2] = {
        vector<double>(
            maxKWidth * myBCols
        ),
        vector<double>(
            maxKWidth * myBCols
        )
    };

    vector<double> cLocal(
        myARows * myBCols,
        0.0
    );

    MPI_Comm rowComm;
    MPI_Comm colComm;

    // Row communicator broadcasts A panels horizontally
    MPI_Comm_split(
        MPI_COMM_WORLD,
        myProcessRow,
        myProcessCol,
        &rowComm
    );

    // Column communicator broadcasts B panels vertically
    MPI_Comm_split(
        MPI_COMM_WORLD,
        myProcessCol,
        myProcessRow,
        &colComm
    );

    // Fill local A and B panel data for a SUMMA step root
    const auto fillPanelForStep = [
        &
    ](
        const int step,
        const int bufferIndex
        ) {
            const int aRoot = step % processGridCols;
            const int bRoot = step % processGridRows;
            const int kStart = kStarts[step];
            const int kWidth = kWidths[step];

            if (myProcessCol == aRoot) {
                double* const aPanel = aPanels[bufferIndex].data();

                // A root fills its local A rows for this k block
                for (int localRow = 0; localRow < myARows; localRow++) {
                    const int globalRow = aRowRange.first + localRow;

                    double* const aRow = &aPanel[localRow * kWidth];

                    for (
                        int kOffset = 0;
                        kOffset < kWidth;
                        kOffset++
                        ) {
                        aRow[kOffset] = valueAt(
                            globalRow,
                            kStart + kOffset,
                            baseSeed
                        );
                    }
                }
            }

            if (myProcessRow == bRoot) {
                double* const bPanel = bPanels[bufferIndex].data();

                // B root fills its local B columns for this k block
                for (
                    int kOffset = 0;
                    kOffset < kWidth;
                    kOffset++
                    ) {
                    const int globalK = kStart + kOffset;

                    double* const bRow = &bPanel[kOffset * myBCols];

                    for (int localCol = 0; localCol < myBCols; localCol++) {
                        bRow[localCol] = valueAt(
                            globalK,
                            bColRange.first + localCol,
                            bMatrixSeed
                        );
                    }
                }
            }
        };

    // Start nonblocking row and column broadcasts for a panel pair
    const auto startBroadcastForStep = [
        &
    ](
        const int step,
        const int bufferIndex,
        MPI_Request requests[2]
        ) {
            const int aRoot = step % processGridCols;
            const int bRoot = step % processGridRows;
            const int kWidth = kWidths[step];

            MPI_Ibcast(
                aPanels[bufferIndex].data(),
                myARows * kWidth,
                MPI_DOUBLE,
                aRoot,
                rowComm,
                &requests[0]
            );

            MPI_Ibcast(
                bPanels[bufferIndex].data(),
                kWidth * myBCols,
                MPI_DOUBLE,
                bRoot,
                colComm,
                &requests[1]
            );
        };

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
        double localBroadcastWaitSeconds = 0.0;

        // Synchronize ranks for fair timing
        MPI_Barrier(
            MPI_COMM_WORLD
        );

        // Start relevant clocks
        const auto totalStart =
            chrono::high_resolution_clock::now();

        MPI_Request activeRequests[2];

        const auto firstFillStart =
            chrono::high_resolution_clock::now();

        // Fill first panel before entering compute loop
        fillPanelForStep(
            0,
            0
        );

        const auto firstFillEnd =
            chrono::high_resolution_clock::now();

        localFillSeconds += chrono::duration<double>(
            firstFillEnd - firstFillStart
        ).count();

        // Start first broadcasts before entering compute loop
        startBroadcastForStep(
            0,
            0,
            activeRequests
        );

        const auto firstWaitStart =
            chrono::high_resolution_clock::now();

        // Wait for first panels before first compute step
        MPI_Waitall(
            2,
            activeRequests,
            MPI_STATUSES_IGNORE
        );

        const auto firstWaitEnd =
            chrono::high_resolution_clock::now();

        localBroadcastWaitSeconds += chrono::duration<double>(
            firstWaitEnd - firstWaitStart
        ).count();

        // Each step computes 1 local k block contribution
        for (
            int step = 0;
            step < stepCount;
            step++
            ) {
            const int currentBuffer = step % 2;
            const int currentKWidth = kWidths[step];

            MPI_Request nextRequests[2];

            const bool hasNextBroadcast =
                step + 1 < stepCount;

            if (hasNextBroadcast) {
                const int nextBuffer = (
                    step + 1
                    ) % 2;

                const auto fillStart =
                    chrono::high_resolution_clock::now();

                // Fill next panel before computing current panel
                fillPanelForStep(
                    step + 1,
                    nextBuffer
                );

                const auto fillEnd =
                    chrono::high_resolution_clock::now();

                localFillSeconds += chrono::duration<double>(
                    fillEnd - fillStart
                ).count();

                // Launch next broadcasts so they can overlap compute
                startBroadcastForStep(
                    step + 1,
                    nextBuffer,
                    nextRequests
                );
            }

            // Start compute timer after next broadcasts are launched
            const auto computeStart =
                chrono::high_resolution_clock::now();

            const double* const aPanel =
                aPanels[currentBuffer].data();

            const double* const bPanel =
                bPanels[currentBuffer].data();

            int localRow = 0;

            // Reuse each B row across 4 C rows
            for (
                ;
                localRow + MICRO_KERNEL_ROWS - 1 < myARows;
                localRow += MICRO_KERNEL_ROWS
                ) {
                // Get 4 output rows from local C block
                double* const cRow0 =
                    &cLocal[(localRow + 0) * myBCols];

                double* const cRow1 =
                    &cLocal[(localRow + 1) * myBCols];

                double* const cRow2 =
                    &cLocal[(localRow + 2) * myBCols];

                double* const cRow3 =
                    &cLocal[(localRow + 3) * myBCols];

                // Get 4 input rows from local A panel
                const double* const aRow0 =
                    &aPanel[(localRow + 0) * currentKWidth];

                const double* const aRow1 =
                    &aPanel[(localRow + 1) * currentKWidth];

                const double* const aRow2 =
                    &aPanel[(localRow + 2) * currentKWidth];

                const double* const aRow3 =
                    &aPanel[(localRow + 3) * currentKWidth];

                // Walk through current shared dimension block
                for (
                    int kOffset = 0;
                    kOffset < currentKWidth;
                    kOffset++
                    ) {
                    // Load A scalars once for reuse across local columns
                    const double aValue0 = aRow0[kOffset];
                    const double aValue1 = aRow1[kOffset];
                    const double aValue2 = aRow2[kOffset];
                    const double aValue3 = aRow3[kOffset];

                    // Current B row is reused by all 4 C rows
                    const double* const bRow =
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

            // Handle remaining rows that do not fill a 4 row block
            for (; localRow < myARows; localRow++) {
                // Get output row for tail computation
                double* const cRow =
                    &cLocal[localRow * myBCols];

                // Get matching A row for tail computation
                const double* const aRow =
                    &aPanel[localRow * currentKWidth];

                // Accumulate current k block into 1 C row
                for (
                    int kOffset = 0;
                    kOffset < currentKWidth;
                    kOffset++
                    ) {
                    const double aValue = aRow[kOffset];

                    // Current B row contributes to all local C columns
                    const double* const bRow =
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

            // Stop compute timer before waiting for next broadcast
            const auto computeEnd =
                chrono::high_resolution_clock::now();

            localComputeSeconds += chrono::duration<double>(
                computeEnd - computeStart
            ).count();

            if (hasNextBroadcast) {
                const auto waitStart =
                    chrono::high_resolution_clock::now();

                // Wait timer captures only broadcast left after compute
                MPI_Waitall(
                    2,
                    nextRequests,
                    MPI_STATUSES_IGNORE
                );

                const auto waitEnd =
                    chrono::high_resolution_clock::now();

                localBroadcastWaitSeconds += chrono::duration<double>(
                    waitEnd - waitStart
                ).count();
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
        double maxBroadcastWaitSeconds = 0.0;
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

        // Broadcast wait time shows uncovered communication after overlap
        MPI_Allreduce(
            &localBroadcastWaitSeconds,
            &maxBroadcastWaitSeconds,
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
                    << "RUNNING SUMMA BENCHMARK "
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
                    << "MPI Process Split: "
                    << processGridRows
                    << " x "
                    << processGridCols
                    << '\n'
                    << "SUMMA Steps: "
                    << stepCount
                    << '\n'
                    << "Reported comparison metric: "
                    << "Compute Time only\n"
                    << "MPI overheads reported separately: "
                    << "Fill Time, Broadcast Wait Time\n"
                    << "Broadcast Wait Time = time blocked in "
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
                << "Broadcast Wait Time: "
                << maxBroadcastWaitSeconds
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
                << "<================================================"
                << "=========>\n";
        }
    }

    MPI_Comm_free(
        &rowComm
    );

    MPI_Comm_free(
        &colComm
    );
}