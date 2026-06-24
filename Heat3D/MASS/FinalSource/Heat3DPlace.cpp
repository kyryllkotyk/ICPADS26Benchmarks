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

#include "Heat3DPlace.h"

 // MASS factory hook used when Places creates Heat3DPlace objects
extern "C" Place* instantiate(
    void* argument
) {
    return new Heat3DPlace(
        argument
    );
}

// MASS cleanup hook for Place objects
extern "C" void destroy(
    Place* object
) {
    delete object;
}

namespace
{
    // Generate deterministic initial value from global cell coordinates
    double valueAt(
        const uint64_t seed,
        const int globalX,
        const int globalY,
        const int globalZ,
        const double initMin,
        const double initMax
    ) {
        uint64_t value = seed;

        value ^= static_cast<uint64_t>(
            globalX
            ) * 0x9e3779b97f4a7c15ULL;

        value ^= static_cast<uint64_t>(
            globalY
            ) * 0xbf58476d1ce4e5b9ULL;

        value ^= static_cast<uint64_t>(
            globalZ
            ) * 0x94d049bb133111ebULL;

        value ^= value >> 30;
        value *= 0xbf58476d1ce4e5b9ULL;
        value ^= value >> 27;
        value *= 0x94d049bb133111ebULL;
        value ^= value >> 31;

        const double normalized = static_cast<double>(
            value
            ) / static_cast<double>(
                UINT64_MAX
                );

        return initMin + normalized * (initMax - initMin);
    }

    // Split 1 global dimension across 1 decomposition axis
    void computeLocal(
        int& localSize,
        int& globalStart,
        const int globalDimension,
        const int rankDimension,
        const int decompositionDimension
    ) {
        const int base = globalDimension / decompositionDimension;
        const int remainder = globalDimension % decompositionDimension;

        localSize = base + (
            rankDimension < remainder
            ? 1
            : 0
            );

        globalStart = rankDimension * base + std::min(
            rankDimension,
            remainder
        );
    }

    // Compute packed face offsets and sizes for 6 axial directions
    void faceLayout(
        const int localSizeX,
        const int localSizeY,
        const int localSizeZ,
        int offsets[HEAT3D_AXIAL_DIRECTIONS],
        int sizes[HEAT3D_AXIAL_DIRECTIONS],
        int& totalFaceDoubles
    ) {
        // X faces store one value for each y,z coordinate
        sizes[0] = localSizeY * localSizeZ;
        sizes[1] = localSizeY * localSizeZ;
        // Y faces store one value for each x,z coordinate
        sizes[2] = localSizeX * localSizeZ;
        sizes[3] = localSizeX * localSizeZ;
        // Z faces store one value for each x,y coordinate
        sizes[4] = localSizeX * localSizeY;
        sizes[5] = localSizeX * localSizeY;

        offsets[0] = 0;

        for (int direction = 1; direction < HEAT3D_AXIAL_DIRECTIONS;
            direction++) {
            offsets[direction] = offsets[direction - 1]
                + sizes[direction - 1];
        }

        totalFaceDoubles = offsets[5] + sizes[5];
    }
}

Heat3DPlace::Heat3DPlace(
    void* argument
) :
    Place(
        argument
    ),
    localSizeX_(
        0
    ),
    localSizeY_(
        0
    ),
    localSizeZ_(
        0
    ),
    localAllocX_(
        0
    ),
    localAllocY_(
        0
    ),
    localAllocZ_(
        0
    ),
    sliceSize_(
        0
    ),
    globalStartX_(
        0
    ),
    globalStartY_(
        0
    ),
    globalStartZ_(
        0
    ),
    globalGridX_(
        0
    ),
    globalGridY_(
        0
    ),
    globalGridZ_(
        0
    ),
    chunkX_(
        0
    ),
    chunkY_(
        0
    ),
    chunkZ_(
        0
    ),
    decompX_(
        0
    ),
    decompY_(
        0
    ),
    decompZ_(
        0
    ),
    alpha_(
        0.5
    ),
    beta_(
        0.1
    ),
    seed_(
        0
    ),
    initMin_(
        0.0
    ),
    initMax_(
        0.0
    ),
    totalFaceDoubles_(
        0
    ) {
    // Clear face layout and neighbor metadata before init() fills it
    std::memset(
        ownFaceOffset_,
        0,
        sizeof(
            ownFaceOffset_
            )
    );

    std::memset(
        ownFaceSize_,
        0,
        sizeof(
            ownFaceSize_
            )
    );

    std::memset(
        nbrInfo_,
        0,
        sizeof(
            nbrInfo_
            )
    );
}

Heat3DPlace::~Heat3DPlace() {
    // MASS owns Place object, but this Place owns its outMessage buffer
    if (outMessage) {
        std::free(
            outMessage
        );

        outMessage = nullptr;
    }
}

// Convert local x,y,z coordinate into flat grid index
int Heat3DPlace::idx(
    const int x,
    const int y,
    const int z
) const {
    return z * sliceSize_ + y * localAllocX_ + x;
}

// Convert sender chunk delta into local receive direction
int Heat3DPlace::directionFromSenderDelta(
    const int deltaX,
    const int deltaY,
    const int deltaZ
) const {
    if (
        deltaY == 0
        && deltaZ == 0
        ) {
        if (deltaX == -1) {
            return 0;
        }

        if (deltaX == 1) {
            return 1;
        }
    }

    if (
        deltaX == 0
        && deltaZ == 0
        ) {
        if (deltaY == -1) {
            return 2;
        }

        if (deltaY == 1) {
            return 3;
        }
    }

    if (
        deltaX == 0
        && deltaY == 0
        ) {
        if (deltaZ == -1) {
            return 4;
        }

        if (deltaZ == 1) {
            return 5;
        }
    }

    return -1;
}

// Fill owned cells and clear ghost layers for a new measured run
void Heat3DPlace::seedLocalGrid() {
    const int totalCells = localAllocX_ * localAllocY_ * localAllocZ_;

    currentGrid_.assign(
        totalCells,
        0.0
    );

    nextGrid_.assign(
        totalCells,
        0.0
    );

    // Fill owned cells from global coordinates
    // Compute only owned cells, leaving ghost layers as communication state
    for (int z = 1; z <= localSizeZ_; z++) {
        const int globalZ = globalStartZ_ + z - 1;

        for (int y = 1; y <= localSizeY_; y++) {
            const int globalY = globalStartY_ + y - 1;

            for (int x = 1; x <= localSizeX_; x++) {
                const int globalX = globalStartX_ + x - 1;

                currentGrid_[idx(
                    x,
                    y,
                    z
                )] = valueAt(
                    seed_,
                    globalX,
                    globalY,
                    globalZ,
                    initMin_,
                    initMax_
                );
            }
        }
    }
}

// Allocate reusable MASS outMessage payload for packed faces
void Heat3DPlace::allocateOutMessage() {
    if (outMessage) {
        std::free(
            outMessage
        );

        outMessage = nullptr;
    }

    const std::size_t bytes = static_cast<std::size_t>(
        HEAT3D_HEADER_BYTES
        ) + static_cast<std::size_t>(
            totalFaceDoubles_
            ) * sizeof(
                double
                );

    outMessage = std::malloc(
        bytes
    );

    if (!outMessage) {
        std::cerr
            << "Heat3DPlace::init failed to allocate outMessage ("
            << bytes
            << " bytes)\n";

        std::exit(
            1
        );
    }

    outMessage_size = static_cast<int>(
        bytes
        );

    // recvHalo returns 1 double to satisfy framework bookkeeping
    inMessage_size = static_cast<int>(
        sizeof(
            double
            )
        );

    int* header = static_cast<int*>(
        outMessage
        );

    // Header identifies sender chunk during placeExchangeAll
    header[0] = chunkX_;
    header[1] = chunkY_;
    header[2] = chunkZ_;
    header[3] = 0;
}

void* Heat3DPlace::init(
    void* argument
) {
    // Last Place dimension varies fastest, so index[2] is x chunk
    chunkX_ = index[2];
    chunkY_ = index[1];
    chunkZ_ = index[0];

    decompX_ = size[2];
    decompY_ = size[1];
    decompZ_ = size[0];

    seed_ = 1;
    initMin_ = 0.0;
    initMax_ = 1.0;
    alpha_ = 0.5;
    beta_ = 0.1;
    globalGridX_ = decompX_;
    globalGridY_ = decompY_;
    globalGridZ_ = decompZ_;

    if (argument) {
        // Runtime config overrides safe constructor defaults
        const Heat3DConfig* config = static_cast<const Heat3DConfig*>(
            argument
            );

        globalGridX_ = config->gridX;
        globalGridY_ = config->gridY;
        globalGridZ_ = config->gridZ;
        seed_ = config->seed;
        initMin_ = config->initMin;
        initMax_ = config->initMax;
        alpha_ = config->alpha;
        beta_ = config->beta;
    }

    // Compute owned size and global start for each dimension
    computeLocal(
        localSizeX_,
        globalStartX_,
        globalGridX_,
        chunkX_,
        decompX_
    );

    computeLocal(
        localSizeY_,
        globalStartY_,
        globalGridY_,
        chunkY_,
        decompY_
    );

    computeLocal(
        localSizeZ_,
        globalStartZ_,
        globalGridZ_,
        chunkZ_,
        decompZ_
    );

    // Allocate 1 ghost layer on every side of local owned block
    localAllocX_ = localSizeX_ + 2;
    localAllocY_ = localSizeY_ + 2;
    localAllocZ_ = localSizeZ_ + 2;
    sliceSize_ = localAllocX_ * localAllocY_;

    seedLocalGrid();

    // Build this Place's packed face layout
    faceLayout(
        localSizeX_,
        localSizeY_,
        localSizeZ_,
        ownFaceOffset_,
        ownFaceSize_,
        totalFaceDoubles_
    );

    // Direction vectors follow x-, x+, y-, y+, z-, z+ order
    static const int deltaX[HEAT3D_AXIAL_DIRECTIONS] = {
        -1,
        1,
        0,
        0,
        0,
        0
    };

    static const int deltaY[HEAT3D_AXIAL_DIRECTIONS] = {
        0,
        0,
        -1,
        1,
        0,
        0
    };

    static const int deltaZ[HEAT3D_AXIAL_DIRECTIONS] = {
        0,
        0,
        0,
        0,
        -1,
        1
    };

    static const int oppositeFace[HEAT3D_AXIAL_DIRECTIONS] = {
        1,
        0,
        3,
        2,
        5,
        4
    };

    // Precompute neighbor payload offsets for every existing face
    for (int direction = 0; direction < HEAT3D_AXIAL_DIRECTIONS;
        direction++) {
        const int neighborX = chunkX_ + deltaX[direction];
        const int neighborY = chunkY_ + deltaY[direction];
        const int neighborZ = chunkZ_ + deltaZ[direction];

        if (
            neighborX < 0
            || neighborX >= decompX_
            || neighborY < 0
            || neighborY >= decompY_
            || neighborZ < 0
            || neighborZ >= decompZ_
            ) {
            nbrInfo_[direction].exists = false;
            continue;
        }

        nbrInfo_[direction].exists = true;
        nbrInfo_[direction].faceToExtract = oppositeFace[direction];

        int neighborSizeX = 0;
        int neighborSizeY = 0;
        int neighborSizeZ = 0;
        int ignoredStart = 0;

        computeLocal(
            neighborSizeX,
            ignoredStart,
            globalGridX_,
            neighborX,
            decompX_
        );

        computeLocal(
            neighborSizeY,
            ignoredStart,
            globalGridY_,
            neighborY,
            decompY_
        );

        computeLocal(
            neighborSizeZ,
            ignoredStart,
            globalGridZ_,
            neighborZ,
            decompZ_
        );

        int neighborOffsets[HEAT3D_AXIAL_DIRECTIONS];
        int neighborSizes[HEAT3D_AXIAL_DIRECTIONS];
        int ignoredTotal = 0;

        faceLayout(
            neighborSizeX,
            neighborSizeY,
            neighborSizeZ,
            neighborOffsets,
            neighborSizes,
            ignoredTotal
        );

        nbrInfo_[direction].offsetInNeighbor =
            neighborOffsets[oppositeFace[direction]];
        nbrInfo_[direction].faceSize =
            neighborSizes[oppositeFace[direction]];
    }

    allocateOutMessage();

    // Register 6 axial neighbors for placeExchangeAll
    addNeighbors(
        neighborPattern::VON_NEUMANN3D
    );

    packAllFaces();

    return nullptr;
}

// Pack all 6 owned boundary faces into outMessage
void Heat3DPlace::packAllFaces() {
    double* output = reinterpret_cast<double*>(
        static_cast<char*>(
            outMessage
            ) + HEAT3D_HEADER_BYTES
        );

    double* pointer = nullptr;

    // Pack x negative face
    pointer = output + ownFaceOffset_[0];

    for (int z = 1; z <= localSizeZ_; z++) {
        for (int y = 1; y <= localSizeY_; y++) {
            *pointer = currentGrid_[idx(
                1,
                y,
                z
            )];
            pointer++;
        }
    }

    // Pack x positive face
    pointer = output + ownFaceOffset_[1];

    for (int z = 1; z <= localSizeZ_; z++) {
        for (int y = 1; y <= localSizeY_; y++) {
            *pointer = currentGrid_[idx(
                localSizeX_,
                y,
                z
            )];
            pointer++;
        }
    }

    // Pack y negative face
    pointer = output + ownFaceOffset_[2];

    for (int z = 1; z <= localSizeZ_; z++) {
        for (int x = 1; x <= localSizeX_; x++) {
            *pointer = currentGrid_[idx(
                x,
                1,
                z
            )];
            pointer++;
        }
    }

    // Pack y positive face
    pointer = output + ownFaceOffset_[3];

    for (int z = 1; z <= localSizeZ_; z++) {
        for (int x = 1; x <= localSizeX_; x++) {
            *pointer = currentGrid_[idx(
                x,
                localSizeY_,
                z
            )];
            pointer++;
        }
    }

    // Pack z negative face
    pointer = output + ownFaceOffset_[4];

    for (int y = 1; y <= localSizeY_; y++) {
        for (int x = 1; x <= localSizeX_; x++) {
            *pointer = currentGrid_[idx(
                x,
                y,
                1
            )];
            pointer++;
        }
    }

    // Pack z positive face
    pointer = output + ownFaceOffset_[5];

    for (int y = 1; y <= localSizeY_; y++) {
        for (int x = 1; x <= localSizeX_; x++) {
            *pointer = currentGrid_[idx(
                x,
                y,
                localSizeZ_
            )];
            pointer++;
        }
    }
}

// Unpack 1 neighbor payload into the matching local ghost layer
void Heat3DPlace::unpackOneFace(
    const int direction,
    const double* neighborPayload
) {
    if (
        !nbrInfo_[direction].exists
        || !neighborPayload
        ) {
        return;
    }

    // Neighbor payload contains all 6 faces, so jump to requested face
    const double* face = neighborPayload
        + nbrInfo_[direction].offsetInNeighbor;

    switch (direction) {
    case 0:
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int y = 1; y <= localSizeY_; y++) {
                currentGrid_[idx(
                    0,
                    y,
                    z
                )] = *face;
                face++;
            }
        }
        break;

    case 1:
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int y = 1; y <= localSizeY_; y++) {
                currentGrid_[idx(
                    localSizeX_ + 1,
                    y,
                    z
                )] = *face;
                face++;
            }
        }
        break;

    case 2:
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    0,
                    z
                )] = *face;
                face++;
            }
        }
        break;

    case 3:
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    localSizeY_ + 1,
                    z
                )] = *face;
                face++;
            }
        }
        break;

    case 4:
        for (int y = 1; y <= localSizeY_; y++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    y,
                    0
                )] = *face;
                face++;
            }
        }
        break;

    case 5:
        for (int y = 1; y <= localSizeY_; y++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    y,
                    localSizeZ_ + 1
                )] = *face;
                face++;
            }
        }
        break;

    default:
        break;
    }
}

// Mirror owned boundary cells into ghost layers at physical boundaries
void Heat3DPlace::applyNeumann() {
    if (!nbrInfo_[0].exists) {
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int y = 1; y <= localSizeY_; y++) {
                currentGrid_[idx(
                    0,
                    y,
                    z
                )] = currentGrid_[idx(
                    1,
                    y,
                    z
                )];
            }
        }
    }

    if (!nbrInfo_[1].exists) {
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int y = 1; y <= localSizeY_; y++) {
                currentGrid_[idx(
                    localSizeX_ + 1,
                    y,
                    z
                )] = currentGrid_[idx(
                    localSizeX_,
                    y,
                    z
                )];
            }
        }
    }

    if (!nbrInfo_[2].exists) {
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    0,
                    z
                )] = currentGrid_[idx(
                    x,
                    1,
                    z
                )];
            }
        }
    }

    if (!nbrInfo_[3].exists) {
        for (int z = 1; z <= localSizeZ_; z++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    localSizeY_ + 1,
                    z
                )] = currentGrid_[idx(
                    x,
                    localSizeY_,
                    z
                )];
            }
        }
    }

    if (!nbrInfo_[4].exists) {
        for (int y = 1; y <= localSizeY_; y++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    y,
                    0
                )] = currentGrid_[idx(
                    x,
                    y,
                    1
                )];
            }
        }
    }

    if (!nbrInfo_[5].exists) {
        for (int y = 1; y <= localSizeY_; y++) {
            for (int x = 1; x <= localSizeX_; x++) {
                currentGrid_[idx(
                    x,
                    y,
                    localSizeZ_ + 1
                )] = currentGrid_[idx(
                    x,
                    y,
                    localSizeZ_
                )];
            }
        }
    }
}

// MASS method called before placeExchangeAll
void* Heat3DPlace::packFaces(
    void* argument
) {
    (void)argument;

    packAllFaces();

    return nullptr;
}

// MASS method called once per received neighbor message
void* Heat3DPlace::recvHalo(
    void* argument
) {
    if (!argument) {
        return new double(
            0.0
            );
    }

    // Header identifies which neighbor chunk sent this payload
    const int* header = static_cast<const int*>(
        argument
        );

    const int senderX = header[0];
    const int senderY = header[1];
    const int senderZ = header[2];

    const int direction = directionFromSenderDelta(
        senderX - chunkX_,
        senderY - chunkY_,
        senderZ - chunkZ_
    );

    if (
        direction >= 0
        && direction < HEAT3D_AXIAL_DIRECTIONS
        && nbrInfo_[direction].exists
        ) {
        const double* payload = reinterpret_cast<const double*>(
            static_cast<const char*>(
                argument
                ) + HEAT3D_HEADER_BYTES
            );

        unpackOneFace(
            direction,
            payload
        );
    }

    return new double(
        0.0
        );
}

// Apply one Jacobi stencil step after halo data is ready
void* Heat3DPlace::computeStep(
    void* argument
) {
    (void)argument;

    // Fill physical boundary ghost layers before stencil compute
    applyNeumann();

    const double alpha = alpha_;
    const double beta = beta_;

    double* __restrict__ next = nextGrid_.data();
    const double* __restrict__ current = currentGrid_.data();

    for (int z = 1; z <= localSizeZ_; z++) {
        const double* currentZ = current + z * sliceSize_;
        double* nextZ = next + z * sliceSize_;

        for (int y = 1; y <= localSizeY_; y++) {
            const double* currentRow = currentZ + y * localAllocX_;
            const double* currentYNeg = currentZ + (y - 1) * localAllocX_;
            const double* currentYPos = currentZ + (y + 1) * localAllocX_;
            const double* currentZNeg = currentZ
                - sliceSize_
                + y * localAllocX_;

            const double* currentZPos = currentZ
                + sliceSize_
                + y * localAllocX_;
            double* nextRow = nextZ + y * localAllocX_;

            for (int x = 1; x <= localSizeX_; x++) {
                nextRow[x] = alpha * currentRow[x]
                    + beta * (
                        currentRow[x - 1]
                        + currentRow[x + 1]
                        + currentYNeg[x]
                        + currentYPos[x]
                        + currentZNeg[x]
                        + currentZPos[x]
                        );
            }
        }
    }

    std::swap(
        currentGrid_,
        nextGrid_
    );

    return nullptr;
}

// Return checksum over owned cells only
void* Heat3DPlace::collectChecksum(
    void* argument
) {
    (void)argument;

    double* output = new double;
    *output = 0.0;

    for (int z = 1; z <= localSizeZ_; z++) {
        for (int y = 1; y <= localSizeY_; y++) {
            for (int x = 1; x <= localSizeX_; x++) {
                *output += currentGrid_[idx(
                    x,
                    y,
                    z
                )];
            }
        }
    }

    return output;
}

// Reset local grid state before another measured run
void* Heat3DPlace::reInit(
    void* argument
) {
    (void)argument;

    // Reuse Place state and reseed only the local grids between runs
    seedLocalGrid();
    packAllFaces();

    return nullptr;
}
