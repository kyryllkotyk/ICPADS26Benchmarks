/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Heat3D                                                          *
 * Library: PM2/MPI                                                           *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef PM2_HEAT_3D_
#define PM2_HEAT_3D_

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <mpi.h>

 // Message tags for each halo direction
#define TAG_NEG_X 100
#define TAG_POS_X 101
#define TAG_NEG_Y 102
#define TAG_POS_Y 103
#define TAG_NEG_Z 104
#define TAG_POS_Z 105

// Hash constants used for coordinate based grid initialization
#define HEAT_HASH_X_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define HEAT_HASH_Y_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define HEAT_HASH_Z_MULTIPLIER 0x94d049bb133111ebULL
#define HEAT_HASH_FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define HEAT_HASH_SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL

using namespace std;

class PM2_Heat3D
{
public:
    PM2_Heat3D() = default;

    // Run full Heat3D benchmark for one grid and decomposition config
    void runBenchmark(
        const bool debug,
        const short timesteps,
        const short totalRuns,
        const short screenshotEvery,
        const uint64_t seed,
        const short globalGridX,
        const short globalGridY,
        const short globalGridZ,
        const double initTemperatureMin,
        const double initTemperatureMax,
        const double alpha,
        const double beta,
        const short decompX,
        const short decompY,
        const short decompZ
    );

private:
    // Split 1 global dimension across 1 decomposition axis
    void computeLocal(
        int& localSize,
        int& globalStart,
        const int globalGridDim,
        const int rankDim,
        const int decompDim
    );

    // Fill owned cells using global coordinate based values
    void generateLocalGrid(
        const uint64_t seed,
        const double initTemperatureMin,
        const double initTemperatureMax,
        vector<double>& localGrid,
        const int totalLocalCells,
        const int globalStartX,
        const int globalStartY,
        const int globalStartZ,
        const int localSizeX,
        const int localSizeY,
        const int localSizeZ
    );

    // Convert 3D block coordinate into flat MPI rank
    int blockToRank(
        const int blockX,
        const int blockY,
        const int blockZ,
        const int decompX,
        const int decompY,
        const int decompZ
    );

    // Check whether block coordinate is inside decomposition grid
    bool inBoundsBlock(
        const int blockX,
        const int blockY,
        const int blockZ,
        const int decompX,
        const int decompY,
        const int decompZ
    );

    // Return neighboring rank for 1 direction or no neighbor sentinel
    int getNeighborRank(
        const int rankX,
        const int rankY,
        const int rankZ,
        const int deltaX,
        const int deltaY,
        const int deltaZ,
        const int decompX,
        const int decompY,
        const int decompZ
    );
};

#endif