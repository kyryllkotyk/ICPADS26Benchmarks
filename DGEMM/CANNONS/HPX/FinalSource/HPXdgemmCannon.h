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

#ifndef HPX_DGEMM_CANNON_
#define HPX_DGEMM_CANNON_

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/async.hpp>
#include <hpx/naming_base/id_type.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "../../include/dgemm_run_report.hpp"

// Hash constants and micro-kernel size used by tile generation and multiply
#define ROW_HASH_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define COL_HASH_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL
#define B_MATRIX_SEED_XOR 0x9e3779b9ULL
#define DOUBLE_SCALE (1.0 / (1ULL << 53))
#define MICRO_KERNEL_ROWS 4

// Supported local matrix initialization modes
enum class InitType : int
{
    HASH = 0,
    DETERMINISTIC = 1,
    IDENTITY = 2,
    ONES = 3
};

// Split one matrix dimension across one grid dimension
std::pair<int, int> split1D(
    const int totalDimensionSize,
    const int processesInDimension,
    const int processIndex
);

// Generate deterministic hash based matrix value
double valueAt(
    const int rowIndex,
    const int colIndex,
    const std::uint64_t seed
);

// Parse initialization mode from command line text
InitType parseInitType(
    const std::string& raw
);

// Local Cannon tile owned by one HPX locality
class HPXDgemmCannonTile
{
public:
    HPXDgemmCannonTile(
        const int rowsOfAAndC,
        const int colsOfBAndC,
        const int sharedDimension,
        const int processGridSize,
        const int processRow,
        const int processCol,
        const std::uint64_t seed,
        const std::string& init
    );

    // Generate A value for local panel fill
    double genA(
        const int globalRow,
        const int globalCol
    ) const;

    // Generate B value for local panel fill
    double genB(
        const int globalRow,
        const int globalCol
    ) const;

    // Fill initial skewed panels for Cannon step 0
    void initialFill();

    // Multiply current panels into local C tile
    void multiplyCurrent(
        const int step
    );

    // Clear output tile before a measured run
    void resetC();

    // Sum local C tile for validation
    double checksum() const;

    // Global matrix dimensions and Cannon grid location
    int rowsOfAAndC = 0;
    int colsOfBAndC = 0;
    int sharedDimension = 0;
    int processGridSize = 0;
    int processRow = 0;
    int processCol = 0;
    int localRowStart = 0;
    int localRowCount = 0;
    int localColStart = 0;
    int localColCount = 0;
    int maxKWidth = 0;

    // Matrix generation state
    std::uint64_t baseSeed = 0;
    std::uint64_t bMatrixSeed = 0;
    InitType initType = InitType::HASH;

    // Active panels, receive buffers, and local output tile
    std::vector<double> aPanel;
    std::vector<double> bPanel;
    std::vector<double> aReceive;
    std::vector<double> bReceive;
    std::vector<double> cTile;

    int currentStep = 0;
};

// HPX runtime entry point for benchmark execution
int hpx_main(
    hpx::program_options::variables_map& variables
);

#endif
