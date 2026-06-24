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

#include "pm2heat3d.h"

 // Sentinel used for faces on global domain boundary
#define NO_NEIGHBOR -1

// MPI rank of each direct neighbor
struct NeighborRanks
{
    int negX = NO_NEIGHBOR;
    int posX = NO_NEIGHBOR;
    int negY = NO_NEIGHBOR;
    int posY = NO_NEIGHBOR;
    int negZ = NO_NEIGHBOR;
    int posZ = NO_NEIGHBOR;
};

// Packed send and receive buffers for each halo face
struct HaloBuffers
{
    vector<double> sendNegX;
    vector<double> sendPosX;
    vector<double> sendNegY;
    vector<double> sendPosY;
    vector<double> sendNegZ;
    vector<double> sendPosZ;

    vector<double> recvNegX;
    vector<double> recvPosX;
    vector<double> recvNegY;
    vector<double> recvPosY;
    vector<double> recvNegZ;
    vector<double> recvPosZ;
};

// Allocate matching send and receive buffers for 1 face
// Empty for faces with no neighbor
static void resizeHaloBufferPair(
    vector<double>& sendBuffer,
    vector<double>& recvBuffer,
    const int neighborRank,
    const int elementCount
) {
    if (neighborRank == NO_NEIGHBOR) {
        sendBuffer.resize(
            0
        );

        recvBuffer.resize(
            0
        );

        return;
    }

    sendBuffer.resize(
        elementCount
    );

    recvBuffer.resize(
        elementCount
    );
}

// Prepare all halo buffer pairs once before persistent requests are created
// After this call, send and receive buffer pointers do not change
static void initializeHaloBuffers(
    HaloBuffers& buffers,
    const NeighborRanks& neighbors,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ
) {
    resizeHaloBufferPair(
        buffers.sendNegX,
        buffers.recvNegX,
        neighbors.negX,
        localSizeY * localSizeZ
    );

    resizeHaloBufferPair(
        buffers.sendPosX,
        buffers.recvPosX,
        neighbors.posX,
        localSizeY * localSizeZ
    );

    resizeHaloBufferPair(
        buffers.sendNegY,
        buffers.recvNegY,
        neighbors.negY,
        localSizeX * localSizeZ
    );

    resizeHaloBufferPair(
        buffers.sendPosY,
        buffers.recvPosY,
        neighbors.posY,
        localSizeX * localSizeZ
    );

    resizeHaloBufferPair(
        buffers.sendNegZ,
        buffers.recvNegZ,
        neighbors.negZ,
        localSizeX * localSizeY
    );

    resizeHaloBufferPair(
        buffers.sendPosZ,
        buffers.recvPosZ,
        neighbors.posZ,
        localSizeX * localSizeY
    );
}

// Create 1 persistent send and 1 persistent receive for 1 neighbor face
// Both bind to long lived halo buffers so they survive grid swaps
static void initSendRecvPair(
    vector<double>& sendBuffer,
    vector<double>& recvBuffer,
    const int neighborRank,
    const int sendTag,
    const int recvTag,
    MPI_Request requests[],
    int& requestCount
) {
    if (neighborRank == NO_NEIGHBOR) {
        return;
    }

    MPI_Send_init(
        sendBuffer.data(),
        static_cast<int>(
            sendBuffer.size()
            ),
        MPI_DOUBLE,
        neighborRank,
        sendTag,
        MPI_COMM_WORLD,
        &requests[requestCount]
    );

    requestCount++;

    MPI_Recv_init(
        recvBuffer.data(),
        static_cast<int>(
            recvBuffer.size()
            ),
        MPI_DOUBLE,
        neighborRank,
        recvTag,
        MPI_COMM_WORLD,
        &requests[requestCount]
    );

    requestCount++;
}

// Create up to 12 persistent halo requests for this rank
// Returns active count used by MPI_Startall and MPI_Waitall
static int createPersistentRequests(
    HaloBuffers& buffers,
    const NeighborRanks& neighbors,
    MPI_Request requests[12]
) {
    int requestCount = 0;

    initSendRecvPair(
        buffers.sendNegX,
        buffers.recvNegX,
        neighbors.negX,
        TAG_NEG_X,
        TAG_POS_X,
        requests,
        requestCount
    );

    initSendRecvPair(
        buffers.sendPosX,
        buffers.recvPosX,
        neighbors.posX,
        TAG_POS_X,
        TAG_NEG_X,
        requests,
        requestCount
    );

    initSendRecvPair(
        buffers.sendNegY,
        buffers.recvNegY,
        neighbors.negY,
        TAG_NEG_Y,
        TAG_POS_Y,
        requests,
        requestCount
    );

    initSendRecvPair(
        buffers.sendPosY,
        buffers.recvPosY,
        neighbors.posY,
        TAG_POS_Y,
        TAG_NEG_Y,
        requests,
        requestCount
    );

    initSendRecvPair(
        buffers.sendNegZ,
        buffers.recvNegZ,
        neighbors.negZ,
        TAG_NEG_Z,
        TAG_POS_Z,
        requests,
        requestCount
    );

    initSendRecvPair(
        buffers.sendPosZ,
        buffers.recvPosZ,
        neighbors.posZ,
        TAG_POS_Z,
        TAG_NEG_Z,
        requests,
        requestCount
    );

    return requestCount;
}

// Free all active persistent requests
// Inactive entries were never created
static void freePersistentRequests(
    MPI_Request requests[],
    const int requestCount
) {
    for (int index = 0; index < requestCount; index++) {
        MPI_Request_free(
            &requests[index]
        );
    }
}

// Pack owned boundary faces into send buffers
// Copy owned edge into ghost layer when no neighbor exists
static void packFacesOrApplyBoundaries(
    double* grid,
    HaloBuffers& buffers,
    const NeighborRanks& neighbors,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize
) {
    const size_t rowBytes = static_cast<size_t>(
        localSizeX
        ) * sizeof(
            double
            );

    if (neighbors.negX == NO_NEIGHBOR) {
        // Domain boundary mirrors first owned x plane into ghost layer
        for (int z = 1; z <= localSizeZ; z++) {
            double* row = grid + z * sliceSize + localAllocX;

            for (int y = 1; y <= localSizeY; y++) {
                row[0] = row[1];
                row += localAllocX;
            }
        }
    }
    else {
        double* out = buffers.sendNegX.data();

        // Pack x negative face as 1 value per y,z coordinate
        for (int z = 1; z <= localSizeZ; z++) {
            const double* row = grid + z * sliceSize + localAllocX;

            for (int y = 1; y <= localSizeY; y++) {
                *out = row[1];
                out++;
                row += localAllocX;
            }
        }
    }

    if (neighbors.posX == NO_NEIGHBOR) {
        const int ghostX = localSizeX + 1;

        // Domain boundary mirrors last owned x plane into ghost layer
        for (int z = 1; z <= localSizeZ; z++) {
            double* row = grid + z * sliceSize + localAllocX;

            for (int y = 1; y <= localSizeY; y++) {
                row[ghostX] = row[localSizeX];
                row += localAllocX;
            }
        }
    }
    else {
        double* out = buffers.sendPosX.data();

        // Pack x positive face as 1 value per y,z coordinate
        for (int z = 1; z <= localSizeZ; z++) {
            const double* row = grid + z * sliceSize + localAllocX;

            for (int y = 1; y <= localSizeY; y++) {
                *out = row[localSizeX];
                out++;
                row += localAllocX;
            }
        }
    }

    if (neighbors.negY == NO_NEIGHBOR) {
        // Domain boundary mirrors first owned y row into ghost row
        for (int z = 1; z <= localSizeZ; z++) {
            double* ghostRow = grid + z * sliceSize;
            const double* firstRow = ghostRow + localAllocX;

            memcpy(
                ghostRow + 1,
                firstRow + 1,
                rowBytes
            );
        }
    }
    else {
        double* out = buffers.sendNegY.data();

        // Pack y negative face as contiguous x rows over z
        for (int z = 1; z <= localSizeZ; z++) {
            const double* row = grid + z * sliceSize + localAllocX;

            memcpy(
                out,
                row + 1,
                rowBytes
            );

            out += localSizeX;
        }
    }

    if (neighbors.posY == NO_NEIGHBOR) {
        const int ghostY = localSizeY + 1;

        // Domain boundary mirrors last owned y row into ghost row
        for (int z = 1; z <= localSizeZ; z++) {
            double* ghostRow = grid
                + z * sliceSize
                + ghostY * localAllocX;
            const double* lastRow = grid
                + z * sliceSize
                + localSizeY * localAllocX;

            memcpy(
                ghostRow + 1,
                lastRow + 1,
                rowBytes
            );
        }
    }
    else {
        double* out = buffers.sendPosY.data();

        // Pack y positive face as contiguous x rows over z
        for (int z = 1; z <= localSizeZ; z++) {
            const double* row = grid
                + z * sliceSize
                + localSizeY * localAllocX;

            memcpy(
                out,
                row + 1,
                rowBytes
            );

            out += localSizeX;
        }
    }

    if (neighbors.negZ == NO_NEIGHBOR) {
        double* ghostSlice = grid;
        const double* firstSlice = grid + sliceSize;

        // Domain boundary mirrors first owned z slice into ghost slice
        for (int y = 1; y <= localSizeY; y++) {
            double* ghostRow = ghostSlice + y * localAllocX;
            const double* firstRow = firstSlice + y * localAllocX;

            memcpy(
                ghostRow + 1,
                firstRow + 1,
                rowBytes
            );
        }
    }
    else {
        double* out = buffers.sendNegZ.data();
        const double* slice = grid + sliceSize;

        // Pack z negative face as contiguous x rows over y
        for (int y = 1; y <= localSizeY; y++) {
            const double* row = slice + y * localAllocX;

            memcpy(
                out,
                row + 1,
                rowBytes
            );

            out += localSizeX;
        }
    }

    if (neighbors.posZ == NO_NEIGHBOR) {
        const int ghostZ = localSizeZ + 1;
        double* ghostSlice = grid + ghostZ * sliceSize;
        const double* lastSlice = grid + localSizeZ * sliceSize;

        // Domain boundary mirrors last owned z slice into ghost slice
        for (int y = 1; y <= localSizeY; y++) {
            double* ghostRow = ghostSlice + y * localAllocX;
            const double* lastRow = lastSlice + y * localAllocX;

            memcpy(
                ghostRow + 1,
                lastRow + 1,
                rowBytes
            );
        }
    }
    else {
        double* out = buffers.sendPosZ.data();
        const double* slice = grid + localSizeZ * sliceSize;

        // Pack z positive face as contiguous x rows over y
        for (int y = 1; y <= localSizeY; y++) {
            const double* row = slice + y * localAllocX;

            memcpy(
                out,
                row + 1,
                rowBytes
            );

            out += localSizeX;
        }
    }
}

// Copy received halo values into local ghost layers
// Each unpack path mirrors the packing order for the matching face
static void unpackFaces(
    double* grid,
    const HaloBuffers& buffers,
    const NeighborRanks& neighbors,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize
) {
    const size_t rowBytes = static_cast<size_t>(
        localSizeX
        ) * sizeof(
            double
            );

    if (neighbors.negX != NO_NEIGHBOR) {
        // Unpack x negative face into ghost column x = 0
        const double* in = buffers.recvNegX.data();

        for (int z = 1; z <= localSizeZ; z++) {
            double* row = grid + z * sliceSize + localAllocX;

            for (int y = 1; y <= localSizeY; y++) {
                row[0] = *in;
                in++;
                row += localAllocX;
            }
        }
    }

    if (neighbors.posX != NO_NEIGHBOR) {
        // Unpack x positive face into ghost column x = localSizeX + 1
        const int ghostX = localSizeX + 1;
        const double* in = buffers.recvPosX.data();

        for (int z = 1; z <= localSizeZ; z++) {
            double* row = grid + z * sliceSize + localAllocX;

            for (int y = 1; y <= localSizeY; y++) {
                row[ghostX] = *in;
                in++;
                row += localAllocX;
            }
        }
    }

    if (neighbors.negY != NO_NEIGHBOR) {
        // Unpack y negative face as contiguous x rows
        const double* in = buffers.recvNegY.data();

        for (int z = 1; z <= localSizeZ; z++) {
            double* row = grid + z * sliceSize;

            memcpy(
                row + 1,
                in,
                rowBytes
            );

            in += localSizeX;
        }
    }

    if (neighbors.posY != NO_NEIGHBOR) {
        // Unpack y positive face as contiguous x rows
        const int ghostY = localSizeY + 1;
        const double* in = buffers.recvPosY.data();

        for (int z = 1; z <= localSizeZ; z++) {
            double* row = grid
                + z * sliceSize
                + ghostY * localAllocX;

            memcpy(
                row + 1,
                in,
                rowBytes
            );

            in += localSizeX;
        }
    }

    if (neighbors.negZ != NO_NEIGHBOR) {
        // Unpack z negative face into ghost slice z = 0
        const double* in = buffers.recvNegZ.data();
        double* slice = grid;

        for (int y = 1; y <= localSizeY; y++) {
            double* row = slice + y * localAllocX;

            memcpy(
                row + 1,
                in,
                rowBytes
            );

            in += localSizeX;
        }
    }

    if (neighbors.posZ != NO_NEIGHBOR) {
        // Unpack z positive face into final ghost slice
        const int ghostZ = localSizeZ + 1;
        const double* in = buffers.recvPosZ.data();
        double* slice = grid + ghostZ * sliceSize;

        for (int y = 1; y <= localSizeY; y++) {
            double* row = slice + y * localAllocX;

            memcpy(
                row + 1,
                in,
                rowBytes
            );

            in += localSizeX;
        }
    }
}

// Apply heat stencil to 1 rectangular owned cell box
// Inner x loop is branch free so the compiler can vectorize it
static void computeStencilBox(
    const double* __restrict__ currentGrid,
    double* __restrict__ nextGrid,
    const int xBegin,
    const int xEnd,
    const int yBegin,
    const int yEnd,
    const int zBegin,
    const int zEnd,
    const int localAllocX,
    const int sliceSize,
    const double alpha,
    const double beta
) {
    if (xBegin > xEnd
        || yBegin > yEnd
        || zBegin > zEnd) {
        return;
    }

    for (int z = zBegin; z <= zEnd; z++) {
        const double* currentZ = currentGrid + z * sliceSize;
        const double* currentZNeg = currentZ - sliceSize;
        const double* currentZPos = currentZ + sliceSize;
        double* nextZ = nextGrid + z * sliceSize;

        for (int y = yBegin; y <= yEnd; y++) {
            const double* currentRow = currentZ + y * localAllocX;
            const double* currentYNeg = currentRow - localAllocX;
            const double* currentYPos = currentRow + localAllocX;
            const double* currentLowerZ = currentZNeg + y * localAllocX;
            const double* currentUpperZ = currentZPos + y * localAllocX;
            double* nextRow = nextZ + y * localAllocX;

            for (int x = xBegin; x <= xEnd; x++) {
                nextRow[x] = alpha * currentRow[x]
                    + beta * (
                        currentRow[x - 1]
                        + currentRow[x + 1]
                        + currentYNeg[x]
                        + currentYPos[x]
                        + currentLowerZ[x]
                        + currentUpperZ[x]
                        );
            }
        }
    }
}

// Compute cells that do not depend on received ghost layers
// Used while halo exchange is in flight
static void computeInteriorStencil(
    const double* __restrict__ currentGrid,
    double* __restrict__ nextGrid,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize,
    const double alpha,
    const double beta
) {
    computeStencilBox(
        currentGrid,
        nextGrid,
        2,
        localSizeX - 1,
        2,
        localSizeY - 1,
        2,
        localSizeZ - 1,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );
}

// Compute 6 boundary slabs that depend on received ghost layers
// Slabs are disjoint and cover the 1 cell shell around the interior box
static void computeBoundaryStencil6Slabs(
    const double* __restrict__ currentGrid,
    double* __restrict__ nextGrid,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize,
    const double alpha,
    const double beta
) {
    // Z negative face: full XY area at z = 1
    computeStencilBox(
        currentGrid,
        nextGrid,
        1,
        localSizeX,
        1,
        localSizeY,
        1,
        1,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    // Z positive face: full XY area at z = localSizeZ
    computeStencilBox(
        currentGrid,
        nextGrid,
        1,
        localSizeX,
        1,
        localSizeY,
        localSizeZ,
        localSizeZ,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    // Y negative band: full X row at y = 1, interior z only
    computeStencilBox(
        currentGrid,
        nextGrid,
        1,
        localSizeX,
        1,
        1,
        2,
        localSizeZ - 1,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    // Y positive band: full X row at y = localSizeY, interior z only
    computeStencilBox(
        currentGrid,
        nextGrid,
        1,
        localSizeX,
        localSizeY,
        localSizeY,
        2,
        localSizeZ - 1,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    // X negative edge: x = 1, interior y and z only
    computeStencilBox(
        currentGrid,
        nextGrid,
        1,
        1,
        2,
        localSizeY - 1,
        2,
        localSizeZ - 1,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    // X positive edge: x = localSizeX, interior y and z only
    computeStencilBox(
        currentGrid,
        nextGrid,
        localSizeX,
        localSizeX,
        2,
        localSizeY - 1,
        2,
        localSizeZ - 1,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );
}

// Compute every owned cell in 1 pass
// Used when the local block is too small for overlap
static void computeFullStencil(
    const double* __restrict__ currentGrid,
    double* __restrict__ nextGrid,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize,
    const double alpha,
    const double beta
) {
    computeStencilBox(
        currentGrid,
        nextGrid,
        1,
        localSizeX,
        1,
        localSizeY,
        1,
        localSizeZ,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );
}

// 1 timestep with compute and halo overlap
// Requires non empty interior box in every dimension
static void doTimestepOverlap(
    double* currentGrid,
    double* nextGrid,
    HaloBuffers& buffers,
    const NeighborRanks& neighbors,
    MPI_Request requests[12],
    const int requestCount,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize,
    const double alpha,
    const double beta,
    double& haloTimeMs,
    double& computeTimeMs
) {
    const auto haloStart = chrono::steady_clock::now();

    // Pack boundary faces before posting persistent transfers
    packFacesOrApplyBoundaries(
        currentGrid,
        buffers,
        neighbors,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize
    );

    if (requestCount > 0) {
        // Start all persistent halo sends and receives together
        MPI_Startall(
            requestCount,
            requests
        );
    }

    const auto haloPostEnd = chrono::steady_clock::now();
    const auto interiorStart = chrono::steady_clock::now();

    // Interior work overlaps with in flight halo exchange
    computeInteriorStencil(
        currentGrid,
        nextGrid,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    const auto interiorEnd = chrono::steady_clock::now();
    const auto haloWaitStart = chrono::steady_clock::now();

    if (requestCount > 0) {
        // Wait only after interior work has had chance to overlap
        MPI_Waitall(
            requestCount,
            requests,
            MPI_STATUSES_IGNORE
        );
    }

    unpackFaces(
        currentGrid,
        buffers,
        neighbors,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize
    );

    const auto haloWaitEnd = chrono::steady_clock::now();
    const auto boundaryStart = chrono::steady_clock::now();

    // Boundary shell is safe only after ghost layers are unpacked
    computeBoundaryStencil6Slabs(
        currentGrid,
        nextGrid,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    const auto boundaryEnd = chrono::steady_clock::now();

    haloTimeMs += chrono::duration<double, milli>(
        haloPostEnd - haloStart
    ).count()
        + chrono::duration<double, milli>(
            haloWaitEnd - haloWaitStart
        ).count();

    computeTimeMs += chrono::duration<double, milli>(
        interiorEnd - interiorStart
    ).count()
        + chrono::duration<double, milli>(
            boundaryEnd - boundaryStart
        ).count();
}

// 1 timestep without overlap
// Used when any local dimension is too small for non empty interior
static void doTimestepSerial(
    double* currentGrid,
    double* nextGrid,
    HaloBuffers& buffers,
    const NeighborRanks& neighbors,
    MPI_Request requests[12],
    const int requestCount,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize,
    const double alpha,
    const double beta,
    double& haloTimeMs,
    double& computeTimeMs
) {
    const auto haloStart = chrono::steady_clock::now();

    // Serial path completes halo exchange before any stencil work
    packFacesOrApplyBoundaries(
        currentGrid,
        buffers,
        neighbors,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize
    );

    if (requestCount > 0) {
        MPI_Startall(
            requestCount,
            requests
        );

        MPI_Waitall(
            requestCount,
            requests,
            MPI_STATUSES_IGNORE
        );
    }

    unpackFaces(
        currentGrid,
        buffers,
        neighbors,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize
    );

    const auto haloEnd = chrono::steady_clock::now();
    const auto computeStart = chrono::steady_clock::now();

    computeFullStencil(
        currentGrid,
        nextGrid,
        localSizeX,
        localSizeY,
        localSizeZ,
        localAllocX,
        sliceSize,
        alpha,
        beta
    );

    const auto computeEnd = chrono::steady_clock::now();

    haloTimeMs += chrono::duration<double, milli>(
        haloEnd - haloStart
    ).count();

    computeTimeMs += chrono::duration<double, milli>(
        computeEnd - computeStart
    ).count();
}

// Get checksum without including ghost layers
static double computeLocalChecksum(
    const vector<double>& grid,
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int localAllocX,
    const int sliceSize
) {
    double checksum = 0.0;

    for (int z = 1; z <= localSizeZ; z++) {
        for (int y = 1; y <= localSizeY; y++) {
            const double* row = grid.data()
                + z * sliceSize
                + y * localAllocX;

            for (int x = 1; x <= localSizeX; x++) {
                checksum += row[x];
            }
        }
    }

    return checksum;
}

// Gather timing and checksum values to rank 0 and print 1 run line
static void reduceAndPrintRun(
    const int mpiRank,
    const int mpiSize,
    const int run,
    const int runCount,
    const int blocksX,
    const int blocksY,
    const int blocksZ,
    const int gridX,
    const int gridY,
    const int gridZ,
    const int timestepCount,
    const double alpha,
    const double beta,
    const bool debug,
    const double localHaloTimeMs,
    const double localComputeTimeMs,
    const double localWallTimeMs,
    const double localChecksum
) {
    // Gather wall, halo, compute, and checksum in 1 compact record
    double localMetrics[4] = {
        localWallTimeMs,
        localHaloTimeMs,
        localComputeTimeMs,
        localChecksum
    };

    vector<double> gatheredMetrics;

    if (mpiRank == 0) {
        gatheredMetrics.resize(
            static_cast<size_t>(
                mpiSize
                ) * 4
        );
    }

    // Rank 0 receives one metric record from each rank
    MPI_Gather(
        localMetrics,
        4,
        MPI_DOUBLE,
        mpiRank == 0
        ? gatheredMetrics.data()
        : nullptr,
        4,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    if (mpiRank != 0) {
        return;
    }

    int wallRank = 0;
    double wallMaxMs = gatheredMetrics[0];
    double globalChecksum = 0.0;
    double haloMaxMs = gatheredMetrics[1];
    double computeMaxMs = gatheredMetrics[2];

    // Track slowest wall rank and max phase time across ranks
    for (int rank = 0; rank < mpiSize; rank++) {
        const int offset = rank * 4;
        const double wallMs = gatheredMetrics[offset];
        const double haloMs = gatheredMetrics[offset + 1];
        const double computeMs = gatheredMetrics[offset + 2];
        const double checksum = gatheredMetrics[offset + 3];

        if (wallMs > wallMaxMs) {
            wallMaxMs = wallMs;
            wallRank = rank;
        }

        if (haloMs > haloMaxMs) {
            haloMaxMs = haloMs;
        }

        if (computeMs > computeMaxMs) {
            computeMaxMs = computeMs;
        }

        globalChecksum += checksum;
    }

    // Report phase times from the rank that determines wall time
    const int wallOffset = wallRank * 4;
    const double haloOnWallRankMs = gatheredMetrics[wallOffset + 1];
    const double computeOnWallRankMs = gatheredMetrics[wallOffset + 2];

    // Print benchmark configuration once before first run result
    if (run == 0) {
        cout
            << "--------------CONFIG-----------------\n";

        cout
            << "MPI Ranks: "
            << mpiSize
            << "\n";

        cout
            << "Decomposition: "
            << blocksX
            << " x "
            << blocksY
            << " x "
            << blocksZ
            << "\n";

        cout
            << "Global Grid: "
            << gridX
            << " x "
            << gridY
            << " x "
            << gridZ
            << "\n";

        cout
            << "Timesteps: "
            << timestepCount
            << "\n";

        cout
            << "Runs: "
            << runCount
            << "\n";

        cout
            << "Alpha: "
            << alpha
            << " | Beta: "
            << beta
            << "\n";

        cout
            << "-------------------------------------\n";
    }

    cout
        << "Run "
        << run
        << " HaloWallRankMs="
        << fixed
        << setprecision(
            3
        )
        << haloOnWallRankMs
        << " ComputeWallRankMs="
        << computeOnWallRankMs
        << " WallMaxMs="
        << wallMaxMs
        << " HaloMaxMs="
        << haloMaxMs
        << " ComputeMaxMs="
        << computeMaxMs
        << " WallRank="
        << wallRank
        << " Checksum="
        << setprecision(
            15
        )
        << globalChecksum
        << "\n";

    if (debug) {
        cout
            << "-------------------------------------\n";
    }
}

// Guard against invalid decompositions where a rank owns no cells
static bool hasValidLocalBlock(
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ
) {
    return localSizeX > 0
        && localSizeY > 0
        && localSizeZ > 0;
}

// True if interior box has at least 1 cell in every dimension
// Interior is [2..Sx-1] x [2..Sy-1] x [2..Sz-1]
static bool overlapIsBeneficial(
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ
) {
    return localSizeX > 2
        && localSizeY > 2
        && localSizeZ > 2;
}

// Benchmark runner with persistent halo requests and compute/halo overlap
void PM2_Heat3D::runBenchmark(
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
) {
    // Screenshots aren't used in PM2 timing runs
    (void)screenshotEvery;

    int mpiRank = 0;

    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &mpiRank
    );

    int mpiSize = 0;

    MPI_Comm_size(
        MPI_COMM_WORLD,
        &mpiSize
    );

    const int timestepCount = static_cast<int>(
        timesteps
        );

    const int runCount = static_cast<int>(
        totalRuns
        );

    const int gridX = static_cast<int>(
        globalGridX
        );

    const int gridY = static_cast<int>(
        globalGridY
        );

    const int gridZ = static_cast<int>(
        globalGridZ
        );

    const int blocksX = static_cast<int>(
        decompX
        );

    const int blocksY = static_cast<int>(
        decompY
        );

    const int blocksZ = static_cast<int>(
        decompZ
        );

    const double alphaValue = static_cast<double>(
        alpha
        );

    const double betaValue = static_cast<double>(
        beta
        );

    if (blocksX * blocksY * blocksZ != mpiSize) {
        if (mpiRank == 0) {
            cerr
                << "decompX * decompY * decompZ must equal process count\n";
        }

        MPI_Abort(
            MPI_COMM_WORLD,
            1
        );

        return;
    }

    // Convert flat rank into 3D block coordinates
    const int rankX = mpiRank % blocksX;
    const int rankY = (
        mpiRank / blocksX
        ) % blocksY;
    const int rankZ = mpiRank / (
        blocksX * blocksY
        );

    int localSizeX = 0;
    int localSizeY = 0;
    int localSizeZ = 0;

    int globalStartX = 0;
    int globalStartY = 0;
    int globalStartZ = 0;

    // Split each global dimension across matching decomposition axis
    computeLocal(
        localSizeX,
        globalStartX,
        gridX,
        rankX,
        blocksX
    );

    computeLocal(
        localSizeY,
        globalStartY,
        gridY,
        rankY,
        blocksY
    );

    computeLocal(
        localSizeZ,
        globalStartZ,
        gridZ,
        rankZ,
        blocksZ
    );

    if (!hasValidLocalBlock(
        localSizeX,
        localSizeY,
        localSizeZ
    )) {
        if (mpiRank == 0) {
            cerr
                << "Each rank must own at least 1 cell in every dimension\n";
        }

        MPI_Abort(
            MPI_COMM_WORLD,
            1
        );

        return;
    }

    // Local allocation includes 1 ghost layer on each side
    const int localAllocX = localSizeX + 2;
    const int localAllocY = localSizeY + 2;
    const int localAllocZ = localSizeZ + 2;
    const int sliceSize = localAllocX * localAllocY;

    const size_t totalLocalCells = static_cast<size_t>(
        localAllocX
        ) * static_cast<size_t>(
            localAllocY
            ) * static_cast<size_t>(
                localAllocZ
                );

    NeighborRanks neighbors;

    // Resolve 6 direct neighbor ranks from 3D block coordinates
    neighbors.negX = getNeighborRank(
        rankX,
        rankY,
        rankZ,
        -1,
        0,
        0,
        blocksX,
        blocksY,
        blocksZ
    );

    neighbors.posX = getNeighborRank(
        rankX,
        rankY,
        rankZ,
        1,
        0,
        0,
        blocksX,
        blocksY,
        blocksZ
    );

    neighbors.negY = getNeighborRank(
        rankX,
        rankY,
        rankZ,
        0,
        -1,
        0,
        blocksX,
        blocksY,
        blocksZ
    );

    neighbors.posY = getNeighborRank(
        rankX,
        rankY,
        rankZ,
        0,
        1,
        0,
        blocksX,
        blocksY,
        blocksZ
    );

    neighbors.negZ = getNeighborRank(
        rankX,
        rankY,
        rankZ,
        0,
        0,
        -1,
        blocksX,
        blocksY,
        blocksZ
    );

    neighbors.posZ = getNeighborRank(
        rankX,
        rankY,
        rankZ,
        0,
        0,
        1,
        blocksX,
        blocksY,
        blocksZ
    );

    HaloBuffers buffers;

    // Halo buffers reach final size here
    // They are never resized after persistent requests bind to them
    initializeHaloBuffers(
        buffers,
        neighbors,
        localSizeX,
        localSizeY,
        localSizeZ
    );

    MPI_Request requests[12];

    const int requestCount = createPersistentRequests(
        buffers,
        neighbors,
        requests
    );

    // Use overlap only when a non empty interior box exists
    const bool useOverlap = overlapIsBeneficial(
        localSizeX,
        localSizeY,
        localSizeZ
    );

    // Double buffer owned cells plus ghost layers
    vector<double> currentGrid(
        totalLocalCells,
        0.0
    );

    vector<double> nextGrid(
        totalLocalCells,
        0.0
    );

    for (int run = 0; run < runCount; run++) {
        // Reinitialize grid for each measured run
        generateLocalGrid(
            seed,
            initTemperatureMin,
            initTemperatureMax,
            currentGrid,
            static_cast<int>(
                totalLocalCells
                ),
            globalStartX,
            globalStartY,
            globalStartZ,
            localSizeX,
            localSizeY,
            localSizeZ
        );

        double totalHaloTimeMs = 0.0;
        double totalComputeTimeMs = 0.0;

        const auto runWallStart = chrono::steady_clock::now();

        // Timed loop includes halo exchange and stencil compute
        for (int step = 0; step < timestepCount; step++) {
            if (useOverlap) {
                doTimestepOverlap(
                    currentGrid.data(),
                    nextGrid.data(),
                    buffers,
                    neighbors,
                    requests,
                    requestCount,
                    localSizeX,
                    localSizeY,
                    localSizeZ,
                    localAllocX,
                    sliceSize,
                    alphaValue,
                    betaValue,
                    totalHaloTimeMs,
                    totalComputeTimeMs
                );
            }
            else {
                doTimestepSerial(
                    currentGrid.data(),
                    nextGrid.data(),
                    buffers,
                    neighbors,
                    requests,
                    requestCount,
                    localSizeX,
                    localSizeY,
                    localSizeZ,
                    localAllocX,
                    sliceSize,
                    alphaValue,
                    betaValue,
                    totalHaloTimeMs,
                    totalComputeTimeMs
                );
            }

            // Next grid becomes current state for following timestep
            swap(
                currentGrid,
                nextGrid
            );
        }

        const auto runWallEnd = chrono::steady_clock::now();

        const double localWallTimeMs = chrono::duration<double, milli>(
            runWallEnd - runWallStart
        ).count();

        // Checksum validates owned cell state after final timestep
        const double localChecksum = computeLocalChecksum(
            currentGrid,
            localSizeX,
            localSizeY,
            localSizeZ,
            localAllocX,
            sliceSize
        );

        reduceAndPrintRun(
            mpiRank,
            mpiSize,
            run,
            runCount,
            blocksX,
            blocksY,
            blocksZ,
            gridX,
            gridY,
            gridZ,
            timestepCount,
            alphaValue,
            betaValue,
            debug,
            totalHaloTimeMs,
            totalComputeTimeMs,
            localWallTimeMs,
            localChecksum
        );
    }

    // Persistent requests must be freed after all runs finish
    freePersistentRequests(
        requests,
        requestCount
    );
}

// Split 1 global dimension across 1 decomposition dimension
void PM2_Heat3D::computeLocal(
    int& localSize,
    int& globalStart,
    const int globalGridDim,
    const int rankDim,
    const int decompDim
) {
    const int base = globalGridDim / decompDim;
    const int remainder = globalGridDim % decompDim;

    if (rankDim < remainder) {
        localSize = base + 1;
    }
    else {
        localSize = base;
    }

    globalStart = rankDim * base + min(
        rankDim,
        remainder
    );
}

// Fill owned cells using global coordinate based values
// Same global cell gets same value in every decomposition
void PM2_Heat3D::generateLocalGrid(
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
) {
    (void)totalLocalCells;

    const auto valueAt = [](
        const uint64_t seed,
        const int globalX,
        const int globalY,
        const int globalZ,
        const double initTemperatureMin,
        const double initTemperatureMax
        ) -> double {
            uint64_t value = seed;

            value ^= static_cast<uint64_t>(
                globalX
                ) * HEAT_HASH_X_MULTIPLIER;

            value ^= static_cast<uint64_t>(
                globalY
                ) * HEAT_HASH_Y_MULTIPLIER;

            value ^= static_cast<uint64_t>(
                globalZ
                ) * HEAT_HASH_Z_MULTIPLIER;

            value ^= value >> 30;
            value *= HEAT_HASH_FIRST_MIX_MULTIPLIER;
            value ^= value >> 27;
            value *= HEAT_HASH_SECOND_MIX_MULTIPLIER;
            value ^= value >> 31;

            const double normalized = static_cast<double>(
                value
                ) / static_cast<double>(
                    numeric_limits<uint64_t>::max()
                    );

            return initTemperatureMin
                + normalized * (initTemperatureMax - initTemperatureMin);
        };

    const int localAllocX = localSizeX + 2;
    const int localAllocY = localSizeY + 2;
    const int sliceSize = localAllocX * localAllocY;

    // Fill only owned cells and leave ghost layers untouched until packing
    for (int z = 1; z <= localSizeZ; z++) {
        const int globalZ = globalStartZ + z - 1;

        for (int y = 1; y <= localSizeY; y++) {
            const int globalY = globalStartY + y - 1;
            double* row = localGrid.data()
                + z * sliceSize
                + y * localAllocX;

            for (int x = 1; x <= localSizeX; x++) {
                const int globalX = globalStartX + x - 1;

                row[x] = valueAt(
                    seed,
                    globalX,
                    globalY,
                    globalZ,
                    initTemperatureMin,
                    initTemperatureMax
                );
            }
        }
    }
}

// Convert a 3D block coordinate into flat MPI rank
int PM2_Heat3D::blockToRank(
    const int blockX,
    const int blockY,
    const int blockZ,
    const int decompX,
    const int decompY,
    const int decompZ
) {
    (void)decompZ;

    return blockZ * (
        decompX * decompY
        ) + blockY * decompX + blockX;
}

// Check whether block coordinate is inside decomposition grid
bool PM2_Heat3D::inBoundsBlock(
    const int blockX,
    const int blockY,
    const int blockZ,
    const int decompX,
    const int decompY,
    const int decompZ
) {
    return blockX >= 0
        && blockX < decompX
        && blockY >= 0
        && blockY < decompY
        && blockZ >= 0
        && blockZ < decompZ;
}

// Return MPI rank of neighboring block or NO_NEIGHBOR if outside
int PM2_Heat3D::getNeighborRank(
    const int rankX,
    const int rankY,
    const int rankZ,
    const int deltaX,
    const int deltaY,
    const int deltaZ,
    const int decompX,
    const int decompY,
    const int decompZ
) {
    const int neighborX = rankX + deltaX;
    const int neighborY = rankY + deltaY;
    const int neighborZ = rankZ + deltaZ;

    if (!inBoundsBlock(
        neighborX,
        neighborY,
        neighborZ,
        decompX,
        decompY,
        decompZ
    )) {
        return NO_NEIGHBOR;
    }

    return blockToRank(
        neighborX,
        neighborY,
        neighborZ,
        decompX,
        decompY,
        decompZ
    );
}
