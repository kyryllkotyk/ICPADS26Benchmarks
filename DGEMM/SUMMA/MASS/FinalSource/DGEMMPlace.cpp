/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: DGEMM (SUMMA)                                                   *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "DGEMMPlace.h"

extern "C" Place* instantiate(
    void* argument
)
{
    return new DGEMMPlace(
        argument
    );
}

extern "C" void destroy(
    Place* object
)
{
    delete object;
}

DGEMMPlace::DGEMMPlace(
    void* argument
) :
    Place(
        argument
    ),
    matrixRows_(
        0
    ),
    matrixCols_(
        0
    ),
    sharedDimension_(
        0
    ),
    processGridRows_(
        0
    ),
    processGridCols_(
        0
    ),
    tileRow_(
        0
    ),
    tileCol_(
        0
    ),
    localRowStart_(
        0
    ),
    localRowCount_(
        0
    ),
    localColStart_(
        0
    ),
    localColCount_(
        0
    ),
    stepCount_(
        0
    ),
    maxKWidth_(
        0
    ),
    maxRowCount_(
        0
    ),
    maxColCount_(
        0
    ),
    aPanelMaxDoubles_(
        0
    ),
    bPanelMaxDoubles_(
        0
    ),
    baseSeed_(
        0
    ),
    bMatrixSeed_(
        0
    ),
    initType_(
        InitType::HASH
    ),
    aPanelOffsetDoubles_(
        0
    ),
    bPanelOffsetDoubles_(
        0
    ),
    currentStep_(
        0
    )
{
    inMessage_size = 0;
    outMessage_size = 0;
}

DGEMMPlace::~DGEMMPlace()
{
    if (outMessage != nullptr) {
        std::free(
            outMessage
        );

        outMessage = nullptr;
    }
}

// Compute greatest common divisor for SUMMA step count
int DGEMMPlace::gcd(
    const int leftValue,
    const int rightValue
)
{
    int currentLeft = leftValue;
    int currentRight = rightValue;

    while (currentRight != 0) {
        const int remainder = currentLeft % currentRight;

        currentLeft = currentRight;
        currentRight = remainder;
    }

    return currentLeft;
}

// Split 1 global dimension into balanced chunks
std::pair<int, int> DGEMMPlace::split1D(
    const int totalDimensionSize,
    const int parts,
    const int partIndex
)
{
    const int baseSize = totalDimensionSize / parts;
    const int extraCells = totalDimensionSize % parts;

    const int localCount = baseSize + (
        partIndex < extraCells
            ? 1
            : 0
    );

    const int localStart = partIndex * baseSize + (
        partIndex < extraCells
            ? partIndex
            : extraCells
    );

    return {
        localStart,
        localCount
    };
}

// Generate matrix values without storing full matrices
// This keeps initialization tied to global matrix coordinates
// The same coordinate always produces the same value
// Process count and tile ownership do not affect generated data
double DGEMMPlace::valueAt(
    const int rowIndex,
    const int colIndex,
    const uint64_t seed
)
{
    uint64_t hashValue = seed;

    hashValue ^= static_cast<uint64_t>(
        rowIndex
    ) * ROW_HASH_MULTIPLIER;

    hashValue ^= static_cast<uint64_t>(
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

// Convert string config into a cheap enum checked by panel generation
void DGEMMPlace::parseInitType(
    const char* rawInitType
)
{
    if (rawInitType == nullptr) {
        initType_ = InitType::HASH;
        return;
    }

    if (std::strncmp(
        rawInitType,
        "ones",
        INIT_TYPE_NAME_BYTES
    ) == 0) {
        initType_ = InitType::ONES;
    }
    else if (std::strncmp(
        rawInitType,
        "identity",
        INIT_TYPE_NAME_BYTES
    ) == 0) {
        initType_ = InitType::IDENTITY;
    }
    else if (std::strncmp(
        rawInitType,
        "deterministic",
        INIT_TYPE_NAME_BYTES
    ) == 0) {
        initType_ = InitType::DETERMINISTIC;
    }
    else {
        initType_ = InitType::HASH;
    }
}

// Generate A panel values for the selected initialization mode
double DGEMMPlace::generateA(
    const int globalRow,
    const int globalCol
) const
{
    switch (initType_) {
        case InitType::ONES:
            return 1.0;

        case InitType::IDENTITY:
            return globalRow == globalCol
                ? 1.0
                : 0.0;

        case InitType::DETERMINISTIC:
            return static_cast<double>(
                globalRow * sharedDimension_ + globalCol
            ) / static_cast<double>(
                matrixRows_ * sharedDimension_
            );

        default:
            return valueAt(
                globalRow,
                globalCol,
                baseSeed_
            );
    }
}

// Generate B panel values for the selected initialization mode
double DGEMMPlace::generateB(
    const int globalRow,
    const int globalCol
) const
{
    switch (initType_) {
        case InitType::ONES:
            return 1.0;

        case InitType::IDENTITY:
            return globalRow == globalCol
                ? 1.0
                : 0.0;

        case InitType::DETERMINISTIC:
            return static_cast<double>(
                globalRow * matrixCols_ + globalCol
            ) / static_cast<double>(
                sharedDimension_ * matrixCols_
            );

        default:
            return valueAt(
                globalRow,
                globalCol,
                bMatrixSeed_
            );
    }
}

// Get A panel owner column for a SUMMA step
int DGEMMPlace::ownerColForStep(
    const int step
) const
{
    return step % processGridCols_;
}

// Get B panel owner row for a SUMMA step
int DGEMMPlace::ownerRowForStep(
    const int step
) const
{
    return step % processGridRows_;
}

// Allocate reusable outgoing message buffer for panel exchange
void DGEMMPlace::allocateOutMessage()
{
    if (outMessage != nullptr) {
        std::free(
            outMessage
        );

        outMessage = nullptr;
    }

    aPanelOffsetDoubles_ = 0;
    bPanelOffsetDoubles_ = aPanelMaxDoubles_;

    const std::size_t doublesTotal =
        static_cast<std::size_t>(
            aPanelMaxDoubles_
        )
        + static_cast<std::size_t>(
            bPanelMaxDoubles_
        );

    const std::size_t bytes =
        static_cast<std::size_t>(
            DGEMM_MESSAGE_HEADER_BYTES
        )
        + doublesTotal * sizeof(
            double
        );

    outMessage = std::malloc(
        bytes
    );

    if (outMessage == nullptr) {
        std::cerr
            << "DGEMMPlace::init: malloc failed for outMessage ("
            << bytes
            << " bytes)"
            << std::endl;

        std::exit(
            1
        );
    }

    std::memset(
        outMessage,
        0,
        bytes
    );

    outMessage_size = static_cast<int>(
        bytes
    );

    // recvPanels returns 1 double for framework bookkeeping
    inMessage_size = static_cast<int>(
        sizeof(
            double
        )
    );

    int* const header = static_cast<int*>(
        outMessage
    );

    header[0] = tileRow_;
    header[1] = tileCol_;
    header[2] = 0;
    header[3] = 0;
}

// Initialize local tile ranges, panel buffers, and neighbor list
void* DGEMMPlace::init(
    void* argument
)
{
    // MASS creates the Place grid as row dimension then column dimension
    tileRow_ = index[0];
    tileCol_ = index[1];
    processGridRows_ = size[0];
    processGridCols_ = size[1];

    if (argument != nullptr) {
        const DGEMMConfig* const config = static_cast<const DGEMMConfig*>(
            argument
        );

        matrixRows_ = config->matrixRows;
        matrixCols_ = config->matrixCols;
        sharedDimension_ = config->sharedDimension;
        baseSeed_ = config->baseSeed;
        bMatrixSeed_ = baseSeed_ ^ B_MATRIX_SEED_XOR;

        parseInitType(
            config->initType
        );
    }
    else {
        matrixRows_ = 2048;
        matrixCols_ = 2048;
        sharedDimension_ = 2048;
        baseSeed_ = 1;
        bMatrixSeed_ = baseSeed_ ^ B_MATRIX_SEED_XOR;
        initType_ = InitType::HASH;
    }

    const std::pair<int, int> rowRange = split1D(
        matrixRows_,
        processGridRows_,
        tileRow_
    );

    localRowStart_ = rowRange.first;
    localRowCount_ = rowRange.second;

    const std::pair<int, int> colRange = split1D(
        matrixCols_,
        processGridCols_,
        tileCol_
    );

    localColStart_ = colRange.first;
    localColCount_ = colRange.second;

    // Step count covers all row and column broadcast roots
    stepCount_ = (
        processGridRows_ / gcd(
            processGridRows_,
            processGridCols_
        )
    ) * processGridCols_;

    // Largest k block bounds panel storage for all steps
    maxKWidth_ = (
        sharedDimension_ + stepCount_ - 1
    ) / stepCount_;

    // Split index 0 receives any extra cells in split1D
    maxRowCount_ = split1D(
        matrixRows_,
        processGridRows_,
        0
    ).second;

    maxColCount_ = split1D(
        matrixCols_,
        processGridCols_,
        0
    ).second;

    aPanelMaxDoubles_ = maxRowCount_ * maxKWidth_;
    bPanelMaxDoubles_ = maxKWidth_ * maxColCount_;

    // Local A panel holds owned or received A data for current step
    aPanel_.assign(
        static_cast<std::size_t>(
            localRowCount_
        ) * maxKWidth_,
        0.0
    );

    // Local B panel holds owned or received B data for current step
    bPanel_.assign(
        static_cast<std::size_t>(
            maxKWidth_
        ) * localColCount_,
        0.0
    );

    // Local C tile is accumulated across all SUMMA steps
    cTile_.assign(
        static_cast<std::size_t>(
            localRowCount_
        ) * localColCount_,
        0.0
    );

    allocateOutMessage();
    setupRowColNeighbors();

    currentStep_ = 0;
    return nullptr;
}

// Register row and column peers used by placeExchangeAll
void DGEMMPlace::setupRowColNeighbors()
{
    std::vector<int*> offsets;

    offsets.reserve(
        static_cast<std::size_t>(
            processGridCols_ - 1 + processGridRows_ - 1
        )
    );

    // Row peers receive A panels from this tile when it owns A
    for (
        int otherCol = 0;
        otherCol < processGridCols_;
        otherCol++
    ) {
        if (otherCol == tileCol_) {
            continue;
        }

        int* const offset = new int[2];
        offset[0] = 0;
        offset[1] = otherCol - tileCol_;

        offsets.push_back(
            offset
        );
    }

    // Column peers receive B panels from this tile when it owns B
    for (
        int otherRow = 0;
        otherRow < processGridRows_;
        otherRow++
    ) {
        if (otherRow == tileRow_) {
            continue;
        }

        int* const offset = new int[2];
        offset[0] = otherRow - tileRow_;
        offset[1] = 0;

        offsets.push_back(
            offset
        );
    }

    addNeighbors(
        offsets
    );
}

// Pack owned A and B panels into outMessage for current SUMMA step
void* DGEMMPlace::publishPanels(
    void* argument
)
{
    (void)argument;

    const int step = currentStep_;

    const std::pair<int, int> kRange = split1D(
        sharedDimension_,
        stepCount_,
        step
    );

    const int kStart = kRange.first;
    const int kWidth = kRange.second;

    int* const header = static_cast<int*>(
        outMessage
    );

    header[0] = tileRow_;
    header[1] = tileCol_;

    double* const outDoubles = reinterpret_cast<double*>(
        static_cast<char*>(
            outMessage
        ) + DGEMM_MESSAGE_HEADER_BYTES
    );

    const bool ownsA = tileCol_ == ownerColForStep(
        step
    );

    const bool ownsB = tileRow_ == ownerRowForStep(
        step
    );

    if (ownsA) {
        double* const aOut = outDoubles + aPanelOffsetDoubles_;

        // A owner fills its local A rows for this k block
        for (int localRow = 0; localRow < localRowCount_; localRow++) {
            const int globalRow = localRowStart_ + localRow;

            double* const outRow = aOut + static_cast<std::size_t>(
                localRow
            ) * maxKWidth_;

            double* const privateRow = aPanel_.data()
                + static_cast<std::size_t>(
                    localRow
                ) * maxKWidth_;

            for (int kOffset = 0; kOffset < kWidth; kOffset++) {
                const double value = generateA(
                    globalRow,
                    kStart + kOffset
                );

                outRow[kOffset] = value;
                privateRow[kOffset] = value;
            }

            if (kWidth < maxKWidth_) {
                const std::size_t tailBytes =
                    static_cast<std::size_t>(
                        maxKWidth_ - kWidth
                    ) * sizeof(
                        double
                    );

                std::memset(
                    outRow + kWidth,
                    0,
                    tailBytes
                );

                std::memset(
                    privateRow + kWidth,
                    0,
                    tailBytes
                );
            }
        }

        header[2] = 1;
    }
    else {
        header[2] = 0;
    }

    if (ownsB) {
        double* const bOut = outDoubles + bPanelOffsetDoubles_;

        // B owner fills its local B columns for this k block
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            const int globalK = kStart + kOffset;

            double* const outRow = bOut + static_cast<std::size_t>(
                kOffset
            ) * localColCount_;

            double* const privateRow = bPanel_.data()
                + static_cast<std::size_t>(
                    kOffset
                ) * localColCount_;

            for (int localCol = 0; localCol < localColCount_; localCol++) {
                const double value = generateB(
                    globalK,
                    localColStart_ + localCol
                );

                outRow[localCol] = value;
                privateRow[localCol] = value;
            }
        }

        if (kWidth < maxKWidth_) {
            const std::size_t tailBytes =
                static_cast<std::size_t>(
                    maxKWidth_ - kWidth
                ) * localColCount_ * sizeof(
                    double
                );

            std::memset(
                bOut + static_cast<std::size_t>(
                    kWidth
                ) * localColCount_,
                0,
                tailBytes
            );

            std::memset(
                bPanel_.data() + static_cast<std::size_t>(
                    kWidth
                ) * localColCount_,
                0,
                tailBytes
            );
        }

        header[3] = 1;
    }
    else {
        header[3] = 0;
    }

    return nullptr;
}

// Unpack received row and column panels into local panel buffers
void* DGEMMPlace::recvPanels(
    void* argument
)
{
    if (argument == nullptr) {
        return new double(
            0.0
        );
    }

    const int* const header = static_cast<const int*>(
        argument
    );

    const int senderTileRow = header[0];
    const int senderTileCol = header[1];
    const int hasAPanel = header[2];
    const int hasBPanel = header[3];

    const double* const payload = reinterpret_cast<const double*>(
        static_cast<const char*>(
            argument
        ) + DGEMM_MESSAGE_HEADER_BYTES
    );

    // Accept A only from a sender in the same tile row
    if (hasAPanel != 0 && senderTileRow == tileRow_) {
        const double* const aIn = payload + aPanelOffsetDoubles_;

        std::memcpy(
            aPanel_.data(),
            aIn,
            static_cast<std::size_t>(
                localRowCount_ * maxKWidth_
            ) * sizeof(
                double
            )
        );
    }

    // Accept B only from a sender in the same tile column
    if (hasBPanel != 0 && senderTileCol == tileCol_) {
        const double* const bIn = payload + bPanelOffsetDoubles_;

        std::memcpy(
            bPanel_.data(),
            bIn,
            static_cast<std::size_t>(
                maxKWidth_ * localColCount_
            ) * sizeof(
                double
            )
        );
    }

    return new double(
        0.0
    );
}

// Accumulate current A and B panels into the local C tile
void* DGEMMPlace::accumulate(
    void* argument
)
{
    (void)argument;

    const int step = currentStep_;

    const std::pair<int, int> kRange = split1D(
        sharedDimension_,
        stepCount_,
        step
    );

    const int kWidth = kRange.second;
    const int localCols = localColCount_;

    double* __restrict__ const cBase = cTile_.data();
    const double* __restrict__ const aBase = aPanel_.data();
    const double* __restrict__ const bBase = bPanel_.data();

    int localRow = 0;

    // Reuse each B row across 4 C rows
    for (
        ;
        localRow + MICRO_KERNEL_ROWS - 1 < localRowCount_;
        localRow += MICRO_KERNEL_ROWS
    ) {
        // Get 4 output rows from local C block
        double* __restrict__ const cRow0 = cBase
            + static_cast<std::size_t>(
                localRow + 0
            ) * localCols;

        double* __restrict__ const cRow1 = cBase
            + static_cast<std::size_t>(
                localRow + 1
            ) * localCols;

        double* __restrict__ const cRow2 = cBase
            + static_cast<std::size_t>(
                localRow + 2
            ) * localCols;

        double* __restrict__ const cRow3 = cBase
            + static_cast<std::size_t>(
                localRow + 3
            ) * localCols;

        // Get 4 input rows from local A panel
        const double* __restrict__ const aRow0 = aBase
            + static_cast<std::size_t>(
                localRow + 0
            ) * maxKWidth_;

        const double* __restrict__ const aRow1 = aBase
            + static_cast<std::size_t>(
                localRow + 1
            ) * maxKWidth_;

        const double* __restrict__ const aRow2 = aBase
            + static_cast<std::size_t>(
                localRow + 2
            ) * maxKWidth_;

        const double* __restrict__ const aRow3 = aBase
            + static_cast<std::size_t>(
                localRow + 3
            ) * maxKWidth_;

        // Walk through current shared dimension block
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            // Load A scalars once for reuse across local columns
            const double aValue0 = aRow0[kOffset];
            const double aValue1 = aRow1[kOffset];
            const double aValue2 = aRow2[kOffset];
            const double aValue3 = aRow3[kOffset];

            // Current B row is reused by all 4 C rows
            const double* __restrict__ const bRow = bBase
                + static_cast<std::size_t>(
                    kOffset
                ) * localCols;

            // Update 1 column position in 4 C rows
            for (int localCol = 0; localCol < localCols; localCol++) {
                const double bValue = bRow[localCol];

                cRow0[localCol] += aValue0 * bValue;
                cRow1[localCol] += aValue1 * bValue;
                cRow2[localCol] += aValue2 * bValue;
                cRow3[localCol] += aValue3 * bValue;
            }
        }
    }

    // Handle remaining rows that do not fill a 4 row block
    for (; localRow < localRowCount_; localRow++) {
        // Get output row for tail computation
        double* __restrict__ const cRow = cBase
            + static_cast<std::size_t>(
                localRow
            ) * localCols;

        // Get matching A row for tail computation
        const double* __restrict__ const aRow = aBase
            + static_cast<std::size_t>(
                localRow
            ) * maxKWidth_;

        // Accumulate current k block into 1 C row
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            const double aValue = aRow[kOffset];

            // Current B row contributes to all local C columns
            const double* __restrict__ const bRow = bBase
                + static_cast<std::size_t>(
                    kOffset
                ) * localCols;

            for (int localCol = 0; localCol < localCols; localCol++) {
                cRow[localCol] += aValue * bRow[localCol];
            }
        }
    }

    // Advance Place local step counter for next pipeline iteration
    currentStep_ = (
        currentStep_ + 1
    ) % stepCount_;

    return nullptr;
}

// Clear C tile before each measured benchmark run
void* DGEMMPlace::resetC(
    void* argument
)
{
    (void)argument;

    std::fill(
        cTile_.begin(),
        cTile_.end(),
        0.0
    );

    currentStep_ = 0;
    return nullptr;
}

// Return local C tile checksum for run validation
void* DGEMMPlace::getChecksum(
    void* argument
)
{
    (void)argument;

    double* const out = new double;
    *out = 0.0;

    for (const double value : cTile_) {
        *out += value;
    }

    return out;
}
