/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: DGEMM (SUMMA)                                                   *
 * Library: PM2/MPI (MadMPI)                                                  *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef PM2_DGEMM_SUMMA_
#define PM2_DGEMM_SUMMA_

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include <mpi.h>

#define ROW_HASH_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define COL_HASH_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL
#define B_MATRIX_SEED_XOR 0x9e3779b9ULL
#define DOUBLE_SCALE (1.0 / (1ULL << 53))
#define MICRO_KERNEL_ROWS 4

using namespace std;

class PM2_DGEMM
{
public:
    /**
     * @brief SUMMA matrix multiplication C = A x B
     * @param runs Number of timed runs
     * @param rowsOfAAndC Row count of A and C
     * @param sharedDimension Column count of A and row count of B
     * @param colsOfBAndC Column count of B and C
     * @param baseSeed RNG seed for matrix generation
     * @param processGridRows Process grid rows
     * @param processGridCols Process grid columns
     */
    void runBenchmark(
        const short runs,
        const int rowsOfAAndC,
        const int sharedDimension,
        const int colsOfBAndC,
        const uint64_t baseSeed,
        const int processGridRows,
        const int processGridCols
    );

private:
    static pair<int, int> split1D(
        const int totalDimensionSize,
        const int processesInDimension,
        const int processIndex
    );

    static double valueAt(
        const int rowIndex,
        const int colIndex,
        const uint64_t seed
    );

    static bool inputValidation(
        const int rank,
        const short runs,
        const int rowsOfAAndC,
        const int sharedDimension,
        const int colsOfBAndC,
        const int processGridRows,
        const int processGridCols
    );
};

#endif