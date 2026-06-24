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

#include "HPXheat3d.h"

namespace hpx
{
    char const hpx_check_boost_version_108900[] = "";
    char const hpx_check_boost_version_107500[] = "";
}

namespace
{
    // State owned by the current HPX locality
    std::unique_ptr<HeatSubdomain> localSubdomain;

    // Geometric rank to HPX locality id lookup table
    std::vector<hpx::id_type> localityTable;

    // Generate deterministic value from global cell coordinates
    double valueAt(
        const std::uint64_t seed,
        const int globalX,
        const int globalY,
        const int globalZ,
        const double initMin,
        const double initMax
    ) {
        std::uint64_t value = seed;

        value ^= static_cast<std::uint64_t>(
            globalX
            ) * HEAT_HASH_X_MULTIPLIER;

        value ^= static_cast<std::uint64_t>(
            globalY
            ) * HEAT_HASH_Y_MULTIPLIER;

        value ^= static_cast<std::uint64_t>(
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
                std::numeric_limits<std::uint64_t>::max()
                );

        return initMin + normalized * (initMax - initMin);
    }

    // Split one global dimension across one decomposition axis
    void computeLocalRange(
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

    // Convert local coordinates, including ghost cells, to flat index
    int gridIndex(
        const int x,
        const int y,
        const int z,
        const int localAllocX,
        const int localAllocY
    ) {
        return z * localAllocX * localAllocY
            + y * localAllocX
            + x;
    }

    HeatDirection oppositeDirection(
        const int direction
    ) {
        // Direction ids are paired as 0/1, 2/3, and 4/5
        return static_cast<HeatDirection>(
            direction ^ 1
            );
    }

    // Return geometric neighbor id or sentinel for physical boundary
    int getNeighborId(
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

        if (
            neighborX < 0
            || neighborX >= decompX
            || neighborY < 0
            || neighborY >= decompY
            || neighborZ < 0
            || neighborZ >= decompZ
            ) {
            return NO_HEAT_NEIGHBOR;
        }

        return neighborZ * (
            decompX * decompY
            ) + neighborY * decompX + neighborX;
    }

    // Resolve the 6 direct neighbors for this subdomain
    void computeNeighborIds(
        HeatSubdomain& subdomain
    ) {
        subdomain.neighborIds[HEAT_NEG_X] = getNeighborId(
            subdomain.rankX,
            subdomain.rankY,
            subdomain.rankZ,
            -1,
            0,
            0,
            subdomain.decompX,
            subdomain.decompY,
            subdomain.decompZ
        );

        subdomain.neighborIds[HEAT_POS_X] = getNeighborId(
            subdomain.rankX,
            subdomain.rankY,
            subdomain.rankZ,
            1,
            0,
            0,
            subdomain.decompX,
            subdomain.decompY,
            subdomain.decompZ
        );

        subdomain.neighborIds[HEAT_NEG_Y] = getNeighborId(
            subdomain.rankX,
            subdomain.rankY,
            subdomain.rankZ,
            0,
            -1,
            0,
            subdomain.decompX,
            subdomain.decompY,
            subdomain.decompZ
        );

        subdomain.neighborIds[HEAT_POS_Y] = getNeighborId(
            subdomain.rankX,
            subdomain.rankY,
            subdomain.rankZ,
            0,
            1,
            0,
            subdomain.decompX,
            subdomain.decompY,
            subdomain.decompZ
        );

        subdomain.neighborIds[HEAT_NEG_Z] = getNeighborId(
            subdomain.rankX,
            subdomain.rankY,
            subdomain.rankZ,
            0,
            0,
            -1,
            subdomain.decompX,
            subdomain.decompY,
            subdomain.decompZ
        );

        subdomain.neighborIds[HEAT_POS_Z] = getNeighborId(
            subdomain.rankX,
            subdomain.rankY,
            subdomain.rankZ,
            0,
            0,
            1,
            subdomain.decompX,
            subdomain.decompY,
            subdomain.decompZ
        );
    }

    // Pack one owned boundary face into reusable send buffer
    void packFace(
        const std::vector<double>& grid,
        const HeatDirection direction,
        const HeatSubdomain& subdomain,
        std::vector<double>& buffer
    ) {
        const int localSizeX = subdomain.localSizeX;
        const int localSizeY = subdomain.localSizeY;
        const int localSizeZ = subdomain.localSizeZ;
        const int localAllocX = subdomain.localAllocX;
        const int localAllocY = subdomain.localAllocY;

        switch (direction) {
        case HeatDirection::NEG_X:
            // x negative face stores one value for each y,z coordinate
            for (int z = 1; z <= localSizeZ; z++) {
                for (int y = 1; y <= localSizeY; y++) {
                    buffer[(z - 1) * localSizeY + (y - 1)] = grid[
                        gridIndex(
                            1,
                            y,
                            z,
                            localAllocX,
                            localAllocY
                        )
                    ];
                }
            }
            break;

        case HeatDirection::POS_X:
            // x positive face mirrors same y,z layout
            for (int z = 1; z <= localSizeZ; z++) {
                for (int y = 1; y <= localSizeY; y++) {
                    buffer[(z - 1) * localSizeY + (y - 1)] = grid[
                        gridIndex(
                            localSizeX,
                            y,
                            z,
                            localAllocX,
                            localAllocY
                        )
                    ];
                }
            }
            break;

        case HeatDirection::NEG_Y:
            // y faces store contiguous x rows across z
            for (int z = 1; z <= localSizeZ; z++) {
                for (int x = 1; x <= localSizeX; x++) {
                    buffer[(z - 1) * localSizeX + (x - 1)] = grid[
                        gridIndex(
                            x,
                            1,
                            z,
                            localAllocX,
                            localAllocY
                        )
                    ];
                }
            }
            break;

        case HeatDirection::POS_Y:
            for (int z = 1; z <= localSizeZ; z++) {
                for (int x = 1; x <= localSizeX; x++) {
                    buffer[(z - 1) * localSizeX + (x - 1)] = grid[
                        gridIndex(
                            x,
                            localSizeY,
                            z,
                            localAllocX,
                            localAllocY
                        )
                    ];
                }
            }
            break;

        case HeatDirection::NEG_Z:
            // z faces store contiguous x rows across y
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    buffer[(y - 1) * localSizeX + (x - 1)] = grid[
                        gridIndex(
                            x,
                            y,
                            1,
                            localAllocX,
                            localAllocY
                        )
                    ];
                }
            }
            break;

        case HeatDirection::POS_Z:
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    buffer[(y - 1) * localSizeX + (x - 1)] = grid[
                        gridIndex(
                            x,
                            y,
                            localSizeZ,
                            localAllocX,
                            localAllocY
                        )
                    ];
                }
            }
            break;
        }
    }

    // Unpack one received face into the matching ghost layer
    void unpackFace(
        std::vector<double>& grid,
        const std::vector<double>& buffer,
        const HeatDirection direction,
        const HeatSubdomain& subdomain
    ) {
        const int localSizeX = subdomain.localSizeX;
        const int localSizeY = subdomain.localSizeY;
        const int localSizeZ = subdomain.localSizeZ;
        const int localAllocX = subdomain.localAllocX;
        const int localAllocY = subdomain.localAllocY;

        switch (direction) {
        case HeatDirection::NEG_X:
            for (int z = 1; z <= localSizeZ; z++) {
                for (int y = 1; y <= localSizeY; y++) {
                    grid[gridIndex(
                        0,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )] = buffer[(z - 1) * localSizeY + (y - 1)];
                }
            }
            break;

        case HeatDirection::POS_X:
            for (int z = 1; z <= localSizeZ; z++) {
                for (int y = 1; y <= localSizeY; y++) {
                    grid[gridIndex(
                        localSizeX + 1,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )] = buffer[(z - 1) * localSizeY + (y - 1)];
                }
            }
            break;

        case HeatDirection::NEG_Y:
            for (int z = 1; z <= localSizeZ; z++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        0,
                        z,
                        localAllocX,
                        localAllocY
                    )] = buffer[(z - 1) * localSizeX + (x - 1)];
                }
            }
            break;

        case HeatDirection::POS_Y:
            for (int z = 1; z <= localSizeZ; z++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        localSizeY + 1,
                        z,
                        localAllocX,
                        localAllocY
                    )] = buffer[(z - 1) * localSizeX + (x - 1)];
                }
            }
            break;

        case HeatDirection::NEG_Z:
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        y,
                        0,
                        localAllocX,
                        localAllocY
                    )] = buffer[(y - 1) * localSizeX + (x - 1)];
                }
            }
            break;

        case HeatDirection::POS_Z:
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        y,
                        localSizeZ + 1,
                        localAllocX,
                        localAllocY
                    )] = buffer[(y - 1) * localSizeX + (x - 1)];
                }
            }
            break;
        }
    }

    // Mirror owned boundary values into physical ghost layers
    void applyNeumannBoundaries(
        std::vector<double>& grid,
        const HeatSubdomain& subdomain
    ) {
        const int localSizeX = subdomain.localSizeX;
        const int localSizeY = subdomain.localSizeY;
        const int localSizeZ = subdomain.localSizeZ;
        const int localAllocX = subdomain.localAllocX;
        const int localAllocY = subdomain.localAllocY;
        const int* neighborIds = subdomain.neighborIds;

        if (neighborIds[HEAT_NEG_X] == NO_HEAT_NEIGHBOR) {
            for (int z = 1; z <= localSizeZ; z++) {
                for (int y = 1; y <= localSizeY; y++) {
                    grid[gridIndex(
                        0,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )] = grid[gridIndex(
                        1,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }

        if (neighborIds[HEAT_POS_X] == NO_HEAT_NEIGHBOR) {
            for (int z = 1; z <= localSizeZ; z++) {
                for (int y = 1; y <= localSizeY; y++) {
                    grid[gridIndex(
                        localSizeX + 1,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )] = grid[gridIndex(
                        localSizeX,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }

        if (neighborIds[HEAT_NEG_Y] == NO_HEAT_NEIGHBOR) {
            for (int z = 1; z <= localSizeZ; z++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        0,
                        z,
                        localAllocX,
                        localAllocY
                    )] = grid[gridIndex(
                        x,
                        1,
                        z,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }

        if (neighborIds[HEAT_POS_Y] == NO_HEAT_NEIGHBOR) {
            for (int z = 1; z <= localSizeZ; z++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        localSizeY + 1,
                        z,
                        localAllocX,
                        localAllocY
                    )] = grid[gridIndex(
                        x,
                        localSizeY,
                        z,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }

        if (neighborIds[HEAT_NEG_Z] == NO_HEAT_NEIGHBOR) {
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        y,
                        0,
                        localAllocX,
                        localAllocY
                    )] = grid[gridIndex(
                        x,
                        y,
                        1,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }

        if (neighborIds[HEAT_POS_Z] == NO_HEAT_NEIGHBOR) {
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    grid[gridIndex(
                        x,
                        y,
                        localSizeZ + 1,
                        localAllocX,
                        localAllocY
                    )] = grid[gridIndex(
                        x,
                        y,
                        localSizeZ,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }
    }

    // Initialize owned cells from global coordinates
    void seedGrid(
        HeatSubdomain& subdomain,
        const std::uint64_t seed,
        const double initMin,
        const double initMax
    ) {
        const int localSizeX = subdomain.localSizeX;
        const int localSizeY = subdomain.localSizeY;
        const int localSizeZ = subdomain.localSizeZ;
        const int localAllocX = subdomain.localAllocX;
        const int localAllocY = subdomain.localAllocY;

        const int totalCells = subdomain.localAllocX
            * subdomain.localAllocY
            * subdomain.localAllocZ;

        subdomain.currentGrid.assign(
            totalCells,
            0.0
        );

        subdomain.nextGrid.assign(
            totalCells,
            0.0
        );

        for (int z = 1; z <= localSizeZ; z++) {
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    const int globalX = subdomain.globalStartX + x - 1;
                    const int globalY = subdomain.globalStartY + y - 1;
                    const int globalZ = subdomain.globalStartZ + z - 1;

                    subdomain.currentGrid[gridIndex(
                        x,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )] = valueAt(
                        seed,
                        globalX,
                        globalY,
                        globalZ,
                        initMin,
                        initMax
                    );
                }
            }
        }
    }

    // Apply 7-point stencil over one owned rectangular box
    void stencilSlab(
        double* __restrict__ nextGrid,
        const double* __restrict__ currentGrid,
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

    // Sum owned cells only so ghost layers do not affect validation
    double computeLocalChecksum(
        const HeatSubdomain& subdomain
    ) {
        const int localSizeX = subdomain.localSizeX;
        const int localSizeY = subdomain.localSizeY;
        const int localSizeZ = subdomain.localSizeZ;
        const int localAllocX = subdomain.localAllocX;
        const int localAllocY = subdomain.localAllocY;

        double sum = 0.0;

        for (int z = 1; z <= localSizeZ; z++) {
            for (int y = 1; y <= localSizeY; y++) {
                for (int x = 1; x <= localSizeX; x++) {
                    sum += subdomain.currentGrid[gridIndex(
                        x,
                        y,
                        z,
                        localAllocX,
                        localAllocY
                    )];
                }
            }
        }

        return sum;
    }
}

// Install locality table used by later halo send actions
void installLocalityTableImpl(
    std::vector<hpx::id_type> table
) {
    localityTable = std::move(
        table
    );
}

HPX_PLAIN_ACTION(
    installLocalityTableImpl,
    InstallLocalityTableAction
);

// Remote action target for one incoming halo face
void putFaceRemoteImpl(
    const HeatDirection direction,
    std::vector<double> data,
    const int step
) {
    (void)step;

    if (!localSubdomain) {
        return;
    }

    const int directionIndex = static_cast<int>(
        direction
        );

    if (
        directionIndex >= 0
        && directionIndex < HEAT_DIRECTION_COUNT
        ) {
        localSubdomain->receiveInbound[directionIndex] = std::move(
            data
        );
    }
}

HPX_PLAIN_ACTION(
    putFaceRemoteImpl,
    PutFaceAction
);

// Initialize one locality subdomain and allocate all reusable buffers
void initWorkerImpl(
    const int localSizeX,
    const int localSizeY,
    const int localSizeZ,
    const int globalStartX,
    const int globalStartY,
    const int globalStartZ,
    const int globalGridX,
    const int globalGridY,
    const int globalGridZ,
    const int decompX,
    const int decompY,
    const int decompZ,
    const int localityId,
    const std::uint64_t seed,
    const double initMin,
    const double initMax,
    const double alpha,
    const double beta
) {
    localSubdomain = std::make_unique<HeatSubdomain>();
    HeatSubdomain& subdomain = *localSubdomain;

    // Store owned size, ghost allocation, and global coordinate metadata

    subdomain.localSizeX = localSizeX;
    subdomain.localSizeY = localSizeY;
    subdomain.localSizeZ = localSizeZ;
    subdomain.localAllocX = localSizeX + 2;
    subdomain.localAllocY = localSizeY + 2;
    subdomain.localAllocZ = localSizeZ + 2;
    subdomain.sliceSize = subdomain.localAllocX * subdomain.localAllocY;

    subdomain.globalStartX = globalStartX;
    subdomain.globalStartY = globalStartY;
    subdomain.globalStartZ = globalStartZ;
    subdomain.globalGridX = globalGridX;
    subdomain.globalGridY = globalGridY;
    subdomain.globalGridZ = globalGridZ;

    subdomain.decompX = decompX;
    subdomain.decompY = decompY;
    subdomain.decompZ = decompZ;
    subdomain.localityId = localityId;

    // Convert geometric locality id into 3D decomposition coordinate
    subdomain.rankX = localityId % decompX;
    subdomain.rankY = (
        localityId / decompX
        ) % decompY;
    subdomain.rankZ = localityId / (
        decompX * decompY
        );

    subdomain.alpha = alpha;
    subdomain.beta = beta;

    computeNeighborIds(
        subdomain
    );

    seedGrid(
        subdomain,
        seed,
        initMin,
        initMax
    );

    // Allocate one send, receive, and inbound buffer per face direction
    const int yzFaceSize = localSizeY * localSizeZ;
    const int xzFaceSize = localSizeX * localSizeZ;
    const int xyFaceSize = localSizeX * localSizeY;

    subdomain.sendBuffer[HEAT_NEG_X].resize(
        yzFaceSize
    );
    subdomain.sendBuffer[HEAT_POS_X].resize(
        yzFaceSize
    );
    subdomain.sendBuffer[HEAT_NEG_Y].resize(
        xzFaceSize
    );
    subdomain.sendBuffer[HEAT_POS_Y].resize(
        xzFaceSize
    );
    subdomain.sendBuffer[HEAT_NEG_Z].resize(
        xyFaceSize
    );
    subdomain.sendBuffer[HEAT_POS_Z].resize(
        xyFaceSize
    );

    subdomain.receiveBuffer[HEAT_NEG_X].resize(
        yzFaceSize
    );
    subdomain.receiveBuffer[HEAT_POS_X].resize(
        yzFaceSize
    );
    subdomain.receiveBuffer[HEAT_NEG_Y].resize(
        xzFaceSize
    );
    subdomain.receiveBuffer[HEAT_POS_Y].resize(
        xzFaceSize
    );
    subdomain.receiveBuffer[HEAT_NEG_Z].resize(
        xyFaceSize
    );
    subdomain.receiveBuffer[HEAT_POS_Z].resize(
        xyFaceSize
    );

    subdomain.receiveInbound[HEAT_NEG_X].resize(
        yzFaceSize
    );
    subdomain.receiveInbound[HEAT_POS_X].resize(
        yzFaceSize
    );
    subdomain.receiveInbound[HEAT_NEG_Y].resize(
        xzFaceSize
    );
    subdomain.receiveInbound[HEAT_POS_Y].resize(
        xzFaceSize
    );
    subdomain.receiveInbound[HEAT_NEG_Z].resize(
        xyFaceSize
    );
    subdomain.receiveInbound[HEAT_POS_Z].resize(
        xyFaceSize
    );
}

HPX_PLAIN_ACTION(
    initWorkerImpl,
    InitWorkerAction
);

// Reseed local grid before a repeated measured run
void reInitWorkerImpl(
    const std::uint64_t seed,
    const double initMin,
    const double initMax
) {
    if (!localSubdomain) {
        return;
    }

    seedGrid(
        *localSubdomain,
        seed,
        initMin,
        initMax
    );
}

HPX_PLAIN_ACTION(
    reInitWorkerImpl,
    ReInitWorkerAction
);

// Move newly arrived halos into buffers consumed by next timestep
void swapHaloReceiveBuffersImpl() {
    if (!localSubdomain) {
        return;
    }

    HeatSubdomain& subdomain = *localSubdomain;

    for (int direction = 0; direction < HEAT_DIRECTION_COUNT; direction++) {
        subdomain.receiveBuffer[direction].swap(
            subdomain.receiveInbound[direction]
        );
    }
}

HPX_PLAIN_ACTION(
    swapHaloReceiveBuffersImpl,
    SwapHaloReceiveBuffersAction
);

// Send current boundary faces before first timestep starts
void sendInitialFacesImpl() {
    if (
        !localSubdomain
        || localityTable.empty()
        ) {
        return;
    }

    HeatSubdomain& subdomain = *localSubdomain;
    const int localityCount = static_cast<int>(
        localityTable.size()
        );

    // Fill initial receive buffers for timestep 0
    for (int direction = 0; direction < HEAT_DIRECTION_COUNT; direction++) {
        const int target = subdomain.neighborIds[direction];

        if (
            target < 0
            || target >= localityCount
            ) {
            continue;
        }

        packFace(
            subdomain.currentGrid,
            static_cast<HeatDirection>(
                direction
                ),
            subdomain,
            subdomain.sendBuffer[direction]
        );

        const HeatDirection receiveDirection = oppositeDirection(
            direction
        );

        hpx::async<PutFaceAction>(
            localityTable[static_cast<std::size_t>(
                target
                )],
            receiveDirection,
            subdomain.sendBuffer[direction],
            -1
        ).get();
    }
}

HPX_PLAIN_ACTION(
    sendInitialFacesImpl,
    SendInitialFacesAction
);

// Execute one bulk synchronous Jacobi timestep on this locality
void runTimestepImpl(
    const int step
) {
    if (
        !localSubdomain
        || localityTable.empty()
        ) {
        return;
    }

    HeatSubdomain& subdomain = *localSubdomain;
    const int localityCount = static_cast<int>(
        localityTable.size()
        );

    // Existing receive buffer holds halo values from prior step exchange
    if (localityCount > 1) {
        for (int direction = 0; direction < HEAT_DIRECTION_COUNT; direction++) {
            if (
                subdomain.neighborIds[direction] >= 0
                && !subdomain.receiveBuffer[direction].empty()
                ) {
                unpackFace(
                    subdomain.currentGrid,
                    subdomain.receiveBuffer[direction],
                    static_cast<HeatDirection>(
                        direction
                        ),
                    subdomain
                );
            }
        }
    }

    // Physical boundaries are handled after remote halos are unpacked
    applyNeumannBoundaries(
        subdomain.currentGrid,
        subdomain
    );

    // Compute full owned block using current ghost layer values
    stencilSlab(
        subdomain.nextGrid.data(),
        subdomain.currentGrid.data(),
        1,
        subdomain.localSizeX,
        1,
        subdomain.localSizeY,
        1,
        subdomain.localSizeZ,
        subdomain.localAllocX,
        subdomain.sliceSize,
        subdomain.alpha,
        subdomain.beta
    );

    if (localityCount > 1) {
        // Send next-grid faces for the following timestep
        std::vector<hpx::future<void>> sends;

        sends.reserve(
            HEAT_DIRECTION_COUNT
        );

        for (int direction = 0; direction < HEAT_DIRECTION_COUNT; direction++) {
            const int target = subdomain.neighborIds[direction];

            if (
                target < 0
                || target >= localityCount
                ) {
                continue;
            }

            packFace(
                subdomain.nextGrid,
                static_cast<HeatDirection>(
                    direction
                    ),
                subdomain,
                subdomain.sendBuffer[direction]
            );

            const HeatDirection receiveDirection = oppositeDirection(
                direction
            );

            sends.push_back(
                hpx::async<PutFaceAction>(
                    localityTable[static_cast<std::size_t>(
                        target
                        )],
                    receiveDirection,
                    subdomain.sendBuffer[direction],
                    step
                )
            );
        }

        for (hpx::future<void>& future : sends) {
            future.get();
        }
    }

    std::swap(
        subdomain.currentGrid,
        subdomain.nextGrid
    );
}

HPX_PLAIN_ACTION(
    runTimestepImpl,
    RunTimestepAction
);

// Return local checksum to root locality
double getRemoteChecksumImpl() {
    if (!localSubdomain) {
        return 0.0;
    }

    return computeLocalChecksum(
        *localSubdomain
    );
}

HPX_PLAIN_ACTION(
    getRemoteChecksumImpl,
    GetChecksumAction
);

// HPX benchmark driver after runtime initialization
int hpx_main(
    hpx::program_options::variables_map& variablesMap
) {
    // Read benchmark parameters from HPX program options
    const int globalGridX = variablesMap["gridX"].as<int>();
    const int globalGridY = variablesMap["gridY"].as<int>();
    const int globalGridZ = variablesMap["gridZ"].as<int>();
    const int decompX = variablesMap["decompX"].as<int>();
    const int decompY = variablesMap["decompY"].as<int>();
    const int decompZ = variablesMap["decompZ"].as<int>();
    const int timesteps = variablesMap["timesteps"].as<int>();
    const int totalRuns = variablesMap["runs"].as<int>();
    const bool debug = variablesMap["debug"].as<int>() != 0;
    const std::uint64_t seed = variablesMap["seed"].as<std::uint64_t>();
    const double initMin = variablesMap["initMin"].as<double>();
    const double initMax = variablesMap["initMax"].as<double>();
    const double alpha = variablesMap["alpha"].as<double>();
    const double beta = variablesMap["beta"].as<double>();

    // End to end timer includes setup, all runs, and checksum collection
    const auto endToEndStart = std::chrono::high_resolution_clock::now();
    const int expectedLocalities = decompX * decompY * decompZ;

    std::vector<hpx::id_type> foundLocalities;

    // HPX localities can appear slightly after runtime startup
    for (
        int attempt = 0;
        attempt < LOCALITY_DISCOVERY_ATTEMPTS;
        attempt++
        ) {
        foundLocalities = hpx::find_all_localities();

        if (static_cast<int>(
            foundLocalities.size()
            ) >= expectedLocalities) {
            break;
        }

        std::this_thread::sleep_for(
            std::chrono::milliseconds(
                LOCALITY_DISCOVERY_SLEEP_MS
            )
        );
    }

    if (static_cast<int>(
        foundLocalities.size()
        ) < expectedLocalities) {
        std::cerr
            << "Error: need "
            << expectedLocalities
            << " localities, got "
            << foundLocalities.size()
            << "\n";

        return hpx::finalize();
    }

    std::vector<hpx::id_type> localities(
        static_cast<std::size_t>(
            expectedLocalities
            )
    );

    std::vector<char> slotFilled(
        static_cast<std::size_t>(
            expectedLocalities
            ),
        0
    );

    // Map HPX locality ids into geometric decomposition order
    for (const hpx::id_type& locality : foundLocalities) {
        const std::uint32_t geometricId =
            hpx::naming::get_locality_id_from_id(
                locality
            );

        if (geometricId >= static_cast<std::uint32_t>(
            expectedLocalities
            )) {
            continue;
        }

        const int localityIndex = static_cast<int>(
            geometricId
            );

        if (slotFilled[static_cast<std::size_t>(
            localityIndex
            )]) {
            std::cerr
                << "Error: duplicate HPX locality id "
                << localityIndex
                << "\n";

            return hpx::finalize();
        }

        slotFilled[static_cast<std::size_t>(
            localityIndex
            )] = 1;

        localities[static_cast<std::size_t>(
            localityIndex
            )] = locality;
    }

    // Require a complete contiguous geometric locality table
    for (int index = 0; index < expectedLocalities; index++) {
        if (!slotFilled[static_cast<std::size_t>(
            index
            )]) {
            std::cerr
                << "Error: missing HPX locality id "
                << index
                << " among "
                << foundLocalities.size()
                << " reported localities\n";

            return hpx::finalize();
        }
    }

    // Print either detailed config table or compact progress line
    if (debug) {
        std::cout
            << "--------------CONFIG-----------------\n"
            << "Localities: "
            << expectedLocalities
            << "\n"
            << "Decomposition: "
            << decompX
            << " x "
            << decompY
            << " x "
            << decompZ
            << "\n"
            << "Global Grid: "
            << globalGridX
            << " x "
            << globalGridY
            << " x "
            << globalGridZ
            << "\n"
            << "Timesteps: "
            << timesteps
            << "\n"
            << "Runs: "
            << totalRuns
            << "\n"
            << "Alpha: "
            << alpha
            << " | Beta: "
            << beta
            << "\n"
            << "Seed: "
            << seed
            << "\n"
            << "Init range: ["
            << initMin
            << ", "
            << initMax
            << "]\n"
            << "-------------------------------------\n"
            << "Run | Simulation(ms) | Checksum\n"
            << "-----------------------------------------------\n"
            << std::flush;
    }
    else {
        std::cout
            << "HPX Heat3D: "
            << expectedLocalities
            << " localities; "
            << totalRuns
            << " run(s) x "
            << timesteps
            << " timesteps; output appears once each run finishes.\n"
            << std::flush;
    }

    // Create subdomain state on every locality
    {
        std::vector<hpx::future<void>> futures;

        futures.reserve(
            static_cast<std::size_t>(
                expectedLocalities
                )
        );

        for (int index = 0; index < expectedLocalities; index++) {
            const int rankX = index % decompX;
            const int rankY = (
                index / decompX
                ) % decompY;
            const int rankZ = index / (
                decompX * decompY
                );

            int localSizeX = 0;
            int localSizeY = 0;
            int localSizeZ = 0;
            int globalStartX = 0;
            int globalStartY = 0;
            int globalStartZ = 0;

            computeLocalRange(
                localSizeX,
                globalStartX,
                globalGridX,
                rankX,
                decompX
            );

            computeLocalRange(
                localSizeY,
                globalStartY,
                globalGridY,
                rankY,
                decompY
            );

            computeLocalRange(
                localSizeZ,
                globalStartZ,
                globalGridZ,
                rankZ,
                decompZ
            );

            futures.push_back(
                hpx::async<InitWorkerAction>(
                    localities[static_cast<std::size_t>(
                        index
                        )],
                    localSizeX,
                    localSizeY,
                    localSizeZ,
                    globalStartX,
                    globalStartY,
                    globalStartZ,
                    globalGridX,
                    globalGridY,
                    globalGridZ,
                    decompX,
                    decompY,
                    decompZ,
                    index,
                    seed,
                    initMin,
                    initMax,
                    alpha,
                    beta
                )
            );
        }

        hpx::wait_all(
            futures
        );
    }

    // Install locality table once so halo actions can target neighbors
    {
        std::vector<hpx::future<void>> futures;

        futures.reserve(
            static_cast<std::size_t>(
                expectedLocalities
                )
        );

        for (int index = 0; index < expectedLocalities; index++) {
            futures.push_back(
                hpx::async<InstallLocalityTableAction>(
                    localities[static_cast<std::size_t>(
                        index
                        )],
                    localities
                )
            );
        }

        hpx::wait_all(
            futures
        );
    }

    if (debug) {
        std::cout
            << "Grid init complete; starting timestep loop.\n"
            << std::flush;
    }

    for (int run = 0; run < totalRuns; run++) {
        // Runs after the first reuse locality state and only reseed grids
        if (run > 0) {
            std::vector<hpx::future<void>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            for (int index = 0; index < expectedLocalities; index++) {
                futures.push_back(
                    hpx::async<ReInitWorkerAction>(
                        localities[static_cast<std::size_t>(
                            index
                            )],
                        seed,
                        initMin,
                        initMax
                    )
                );
            }

            hpx::wait_all(
                futures
            );
        }

        if (expectedLocalities > 1) {
            std::vector<hpx::future<void>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            // Initial halo fill prepares recvBuffer for first unpack
            for (int index = 0; index < expectedLocalities; index++) {
                futures.push_back(
                    hpx::async<SendInitialFacesAction>(
                        localities[static_cast<std::size_t>(
                            index
                            )]
                    )
                );
            }

            hpx::wait_all(
                futures
            );

            futures.clear();

            for (int index = 0; index < expectedLocalities; index++) {
                futures.push_back(
                    hpx::async<SwapHaloReceiveBuffersAction>(
                        localities[static_cast<std::size_t>(
                            index
                            )]
                    )
                );
            }

            hpx::wait_all(
                futures
            );
        }

        // Timed simulation section excludes setup and checksum collection
        const auto simulationStart = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < timesteps; step++) {
            std::vector<hpx::future<void>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            for (int index = 0; index < expectedLocalities; index++) {
                futures.push_back(
                    hpx::async<RunTimestepAction>(
                        localities[static_cast<std::size_t>(
                            index
                            )],
                        step
                    )
                );
            }

            hpx::wait_all(
                futures
            );

            if (expectedLocalities > 1) {
                // Swap after all timestep actions finish their sends
                futures.clear();

                // Make halos produced by this step visible to next step
                for (int index = 0; index < expectedLocalities; index++) {
                    futures.push_back(
                        hpx::async<SwapHaloReceiveBuffersAction>(
                            localities[static_cast<std::size_t>(
                                index
                                )]
                        )
                    );
                }

                hpx::wait_all(
                    futures
                );
            }
        }

        const auto simulationEnd = std::chrono::high_resolution_clock::now();
        double globalChecksum = 0.0;

        {
            // Pull one checksum from each locality and reduce on root
            std::vector<hpx::future<double>> futures;

            futures.reserve(
                static_cast<std::size_t>(
                    expectedLocalities
                    )
            );

            for (int index = 0; index < expectedLocalities; index++) {
                futures.push_back(
                    hpx::async<GetChecksumAction>(
                        localities[static_cast<std::size_t>(
                            index
                            )]
                    )
                );
            }

            for (hpx::future<double>& future : futures) {
                globalChecksum += future.get();
            }
        }

        const long simulationMs = std::chrono::duration_cast<
            std::chrono::milliseconds
        >(
            simulationEnd - simulationStart
        ).count();

        if (debug) {
            std::cout
                << std::setw(
                    3
                )
                << run
                << " | "
                << std::setw(
                    14
                )
                << simulationMs
                << " | "
                << std::setprecision(
                    15
                )
                << globalChecksum
                << "\n"
                << std::flush;
        }
    }

    const auto endToEndEnd = std::chrono::high_resolution_clock::now();

    const long endToEndMs = std::chrono::duration_cast<
        std::chrono::milliseconds
    >(
        endToEndEnd - endToEndStart
    ).count();

    std::cout
        << "elapsedMs(END TO END)="
        << endToEndMs
        << "\n";

    return hpx::finalize();
}
