/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: DGEMM (SUMMA)                                                   *
 * Library: HPX                                                               *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef HPX_DGEMM_SUMMA_
#define HPX_DGEMM_SUMMA_

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/util.hpp>
#include <hpx/naming_base/id_type.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "../../include/dgemm_run_report.hpp"

 // Hash constants used for coordinate based matrix generation
#define ROW_HASH_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define COL_HASH_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL
#define B_MATRIX_SEED_XOR 0x9e3779b9ULL
#define DOUBLE_SCALE (1.0 / (1ULL << 53))

// Local multiply processes 4 C rows per inner kernel block
#define MICRO_KERNEL_ROWS 4

// Locality discovery retry settings for multi locality startup
#define LOCALITY_DISCOVERY_ATTEMPTS 300
#define LOCALITY_DISCOVERY_SLEEP_MS 100

// Matrix initialization modes supported by benchmark runs
enum class InitType : int
{
    HASH = 0,
    DETERMINISTIC = 1,
    IDENTITY = 2,
    ONES = 3
};

// Compute greatest common divisor for SUMMA step count
int gcdInt(
    int leftValue,
    int rightValue
);

// Split 1 global dimension into balanced local ranges
std::pair<int, int> split1D(
    const int totalDimensionSize,
    const int processesInDimension,
    const int processIndex
);

// Generate deterministic matrix value from global coordinates
double valueAt(
    const int rowIndex,
    const int colIndex,
    const std::uint64_t seed
);

// Convert command line init mode string into enum value
InitType parseInitType(
    const std::string& rawInput
);

// Store one logical SUMMA tile and its local panel buffers
struct LocalTile
{
    // Global matrix dimensions
    int rowsOfAAndC;
    int colsOfBAndC;
    int sharedDimension;

    // Logical 2D process grid dimensions
    int processGridRows;
    int processGridCols;

    // This tile's logical grid coordinate
    int tileRow;
    int tileCol;

    // Owned output row and column ranges
    int localRowStart;
    int localRowCount;
    int localColStart;
    int localColCount;

    // SUMMA panel schedule information
    int stepCount;
    int maxKWidth;

    // Matrix generation settings
    std::uint64_t baseSeed;
    std::uint64_t bMatrixSeed;
    InitType initType;

    // Reused panel and output buffers for this tile
    std::vector<double> aPanel;
    std::vector<double> bPanel;
    std::vector<double> cTile;

    // Create local tile metadata and allocate panel storage
    LocalTile(
        const int rowsOfAAndC,
        const int colsOfBAndC,
        const int sharedDimension,
        const int processGridRows,
        const int processGridCols,
        const int tileRow,
        const int tileCol,
        const std::uint64_t baseSeed,
        const std::string& initTypeName
    );

    // Generate A value for selected initialization mode
    double genA(
        const int globalRow,
        const int globalK
    ) const;

    // Generate B value for selected initialization mode
    double genB(
        const int globalK,
        const int globalCol
    ) const;

    // Fill local A panel when this tile is row broadcast root
    void fillAPanelForStep(
        const int step,
        std::vector<double>& outputPanel
    ) const;

    // Fill local B panel when this tile is column broadcast root
    void fillBPanelForStep(
        const int step,
        std::vector<double>& outputPanel
    ) const;

    // Store received A panel payload for current SUMMA step
    void storeAFromPayload(
        const int step,
        const std::vector<double>& payload
    );

    // Store received B panel payload for current SUMMA step
    void storeBFromPayload(
        const int step,
        const std::vector<double>& payload
    );

    // Accumulate current A and B panels into local C tile
    void accumulatePanel(
        const int step
    );

    // Clear local C tile before a measured run
    void resetC();

    // Sum local C tile values for final validation
    double checksum() const;
};

// Run full HPX SUMMA benchmark after runtime initialization
int hpx_main(
    hpx::program_options::variables_map& variablesMap
);

#endif