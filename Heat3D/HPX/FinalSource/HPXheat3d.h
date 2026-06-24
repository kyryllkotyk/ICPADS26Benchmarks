/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Heat3D                                                          *
 * Library: HPX                                                               *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef HPX_HEAT3D_
#define HPX_HEAT3D_

#ifndef HPX_UTIL_FROM_STRING_INCLUDED
#define HPX_UTIL_FROM_STRING_INCLUDED
#include <hpx/util/bad_lexical_cast.hpp>
#include <hpx/util/from_string.hpp>
#endif

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/util.hpp>
#include <hpx/naming_base/id_type.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

#define HEAT_DIRECTION_COUNT 6
#define NO_HEAT_NEIGHBOR -1

#define HEAT_NEG_X 0
#define HEAT_POS_X 1
#define HEAT_NEG_Y 2
#define HEAT_POS_Y 3
#define HEAT_NEG_Z 4
#define HEAT_POS_Z 5

#define HEAT_HASH_X_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define HEAT_HASH_Y_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define HEAT_HASH_Z_MULTIPLIER 0x94d049bb133111ebULL
#define HEAT_HASH_FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define HEAT_HASH_SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL
#define LOCALITY_DISCOVERY_ATTEMPTS 50
#define LOCALITY_DISCOVERY_SLEEP_MS 200

namespace hpx
{
    extern char const hpx_check_boost_version_108900[];
    extern char const hpx_check_boost_version_107500[];
}

// Halo face direction used by remote face transfer actions
enum class HeatDirection : int
{
    NEG_X = HEAT_NEG_X,
    POS_X = HEAT_POS_X,
    NEG_Y = HEAT_NEG_Y,
    POS_Y = HEAT_POS_Y,
    NEG_Z = HEAT_NEG_Z,
    POS_Z = HEAT_POS_Z
};

// One HPX locality's owned block, ghost layers, and halo buffers
struct HeatSubdomain
{
    std::vector<double> currentGrid;
    std::vector<double> nextGrid;

    int localSizeX = 0;
    int localSizeY = 0;
    int localSizeZ = 0;

    int localAllocX = 0;
    int localAllocY = 0;
    int localAllocZ = 0;
    int sliceSize = 0;

    int globalStartX = 0;
    int globalStartY = 0;
    int globalStartZ = 0;

    int globalGridX = 0;
    int globalGridY = 0;
    int globalGridZ = 0;

    int decompX = 0;
    int decompY = 0;
    int decompZ = 0;

    int localityId = 0;
    int rankX = 0;
    int rankY = 0;
    int rankZ = 0;

    double alpha = 0.5;
    double beta = 0.1;

    int neighborIds[HEAT_DIRECTION_COUNT] = {
        NO_HEAT_NEIGHBOR,
        NO_HEAT_NEIGHBOR,
        NO_HEAT_NEIGHBOR,
        NO_HEAT_NEIGHBOR,
        NO_HEAT_NEIGHBOR,
        NO_HEAT_NEIGHBOR
    };

    // Send, receive, and inbound buffers are reused every timestep
    std::vector<double> sendBuffer[HEAT_DIRECTION_COUNT];
    std::vector<double> receiveBuffer[HEAT_DIRECTION_COUNT];
    std::vector<double> receiveInbound[HEAT_DIRECTION_COUNT];
};

int hpx_main(
    hpx::program_options::variables_map& variablesMap
);

#endif
