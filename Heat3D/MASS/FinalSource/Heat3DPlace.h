/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Heat3D                                                          *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef HEAT3D_PLACE_
#define HEAT3D_PLACE_

#include "MASS_base.h"
#include "MethodRegistry.h"
#include "Place.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#define HEAT3D_HEADER_INTS 4
#define HEAT3D_HEADER_BYTES \
    (HEAT3D_HEADER_INTS * static_cast<int>(sizeof(int)))
#define HEAT3D_AXIAL_DIRECTIONS 6

 // Benchmark configuration copied to every Place during initialization
struct Heat3DConfig
{
    int gridX;
    int gridY;
    int gridZ;
    uint64_t seed;
    double initMin;
    double initMax;
    double alpha;
    double beta;
};

// MASS Place that owns 1 Heat3D subdomain
class Heat3DPlace : public Place
{
public:
    explicit Heat3DPlace(
        void* argument
    );

    ~Heat3DPlace() override;

    // Initialize local grid, chunk metadata, and neighbor face layout
    void* init(
        void* argument
    );

    // Pack owned boundary faces into MASS outMessage buffer
    void* packFaces(
        void* argument
    );

    // Receive neighbor face data and unpack it into ghost layers
    void* recvHalo(
        void* argument
    );

    // Apply 7 point stencil and swap local grid buffers
    void* computeStep(
        void* argument
    );

    // Return checksum over owned cells only
    void* collectChecksum(
        void* argument
    );

    // Reset local grid before next measured run
    void* reInit(
        void* argument
    );

    MASS_DISPATCH_TABLE(
        Heat3DPlace,
        MASS_METHOD(
            Heat3DPlace,
            init
        ),
        MASS_METHOD(
            Heat3DPlace,
            packFaces
        ),
        MASS_METHOD(
            Heat3DPlace,
            recvHalo
        ),
        MASS_METHOD(
            Heat3DPlace,
            computeStep
        ),
        MASS_METHOD(
            Heat3DPlace,
            collectChecksum
        ),
        MASS_METHOD(
            Heat3DPlace,
            reInit
        )
    )

private:
    // Metadata needed to find one neighbor face inside received payload
    struct NbrInfo
    {
        bool exists;
        int faceToExtract;
        int offsetInNeighbor;
        int faceSize;
    };

    // Convert local x,y,z coordinate into flat grid index
    int idx(
        const int x,
        const int y,
        const int z
    ) const;

    // Fill owned cells from deterministic global coordinate values
    void seedLocalGrid();

    // Allocate MASS outMessage storage for all packed faces
    void allocateOutMessage();

    // Pack all owned faces that neighbors may request
    void packAllFaces();

    // Unpack 1 received neighbor face into matching ghost layer
    void unpackOneFace(
        const int direction,
        const double* neighborPayload
    );

    // Mirror owned boundary cells into ghost layers at physical edges
    void applyNeumann();

    // Convert sender chunk delta into local receive direction
    int directionFromSenderDelta(
        const int deltaX,
        const int deltaY,
        const int deltaZ
    ) const;

    // Double buffered local grid with ghost layers included
    std::vector<double> currentGrid_;
    std::vector<double> nextGrid_;

    // Local owned sizes and allocation sizes including ghost layers
    int localSizeX_;
    int localSizeY_;
    int localSizeZ_;
    int localAllocX_;
    int localAllocY_;
    int localAllocZ_;
    int sliceSize_;

    // Global start coordinate for this Place's owned subdomain
    int globalStartX_;
    int globalStartY_;
    int globalStartZ_;

    // Global grid dimensions
    int globalGridX_;
    int globalGridY_;
    int globalGridZ_;

    // This Place's chunk coordinate in decomposition grid
    int chunkX_;
    int chunkY_;
    int chunkZ_;

    // Decomposition dimensions
    int decompX_;
    int decompY_;
    int decompZ_;

    // Stencil weights and deterministic initialization settings
    double alpha_;
    double beta_;
    uint64_t seed_;
    double initMin_;
    double initMax_;

    // Offsets and sizes for this Place's packed face payload
    int ownFaceOffset_[HEAT3D_AXIAL_DIRECTIONS];
    int ownFaceSize_[HEAT3D_AXIAL_DIRECTIONS];
    int totalFaceDoubles_;

    // Neighbor payload metadata for all 6 axial directions
    NbrInfo nbrInfo_[HEAT3D_AXIAL_DIRECTIONS];
};

#endif