#include "DGEMMCannonPlace.h"

#define MICRO_KERNEL_ROWS 4

// MASS creates 1 Place object for each Cannon grid cell
extern "C" Place *instantiate(
    void *argument
) {
    return new DGEMMCannonPlace(
        argument
    );
}

// MASS calls this when Place objects are released
extern "C" void destroy(
    Place *object
) {
    delete object;
}

// Start with empty matrix state until init receives config data
DGEMMCannonPlace::DGEMMCannonPlace(
    void *argument
) :
    Place(
        argument
    ),
    M_(0),
    N_(0),
    K_(0),
    P_(0),
    pi_(0),
    pj_(0),
    myRowStart_(0),
    myRowCount_(0),
    myColStart_(0),
    myColCount_(0),
    maxKWidth_(0),
    maxRowCount_(0),
    maxColCount_(0),
    aPanelMaxDoubles_(0),
    bPanelMaxDoubles_(0),
    baseSeed_(0),
    bSeed_(0),
    initType_(InitType::HASH),
    aPanelOffsetDoubles_(0),
    bPanelOffsetDoubles_(0),
    currentStep_(0)
{
    inMessage_size = 0;
    outMessage_size = 0;
}

// Release shift payload owned by this Place
DGEMMCannonPlace::~DGEMMCannonPlace()
{
    // Free payload if it was allocated
    if (outMessage) {
        std::free(
            outMessage
        );

        outMessage = nullptr;
    }
}

// Split 1 global dimension into balanced chunks
std::pair<int, int> DGEMMCannonPlace::split1D(
    const int totalSize,
    const int numParts,
    const int index
) {
    const int baseSize = totalSize / numParts;
    const int extraCells = totalSize % numParts;

    const int localCount = baseSize + (
        index < extraCells
            ? 1
            : 0
    );

    const int localStart = index * baseSize + (
        index < extraCells
            ? index
            : extraCells
    );

    return {
        localStart,
        localCount
    };
}

// Generate matrix values without storing full matrices
double DGEMMCannonPlace::valueAt(
    const int rowIndex,
    const int colIndex,
    const uint64_t seed
) {
    uint64_t hashValue = seed;

    hashValue ^= static_cast<uint64_t>(
        rowIndex
    ) * MASS_DGEMM_ROW_HASH_MULTIPLIER;

    hashValue ^= static_cast<uint64_t>(
        colIndex
    ) * MASS_DGEMM_COL_HASH_MULTIPLIER;

    hashValue ^= (
        hashValue >> 30
    );

    hashValue *= MASS_DGEMM_FIRST_MIX_MULTIPLIER;

    hashValue ^= (
        hashValue >> 27
    );

    hashValue *= MASS_DGEMM_SECOND_MIX_MULTIPLIER;

    hashValue ^= (
        hashValue >> 31
    );

    return (
        hashValue >> 11
    ) * MASS_DGEMM_DOUBLE_SCALE;
}

// Select matrix initialization mode from command line config
void DGEMMCannonPlace::parseInitType(
    const char *raw
) {
    // Missing or invalid init type falls back to hash values
    if (raw == nullptr) {
        initType_ = InitType::HASH;
        return;
    }

    if (std::strncmp(
        raw,
        "ones",
        32
    ) == 0) {
        initType_ = InitType::ONES;
    }
    else if (std::strncmp(
        raw,
        "identity",
        32
    ) == 0) {
        initType_ = InitType::IDENTITY;
    }
    else if (std::strncmp(
        raw,
        "deterministic",
        32
    ) == 0) {
        initType_ = InitType::DETERMINISTIC;
    }
    else {
        initType_ = InitType::HASH;
    }
}

// Generate A value for requested global coordinate
double DGEMMCannonPlace::generateA(
    const int globalRow,
    const int globalCol
) const {
    // Keep alternate modes for validation and quick sanity runs
    switch (initType_) {
    case InitType::ONES:
        return 1.0;

    case InitType::IDENTITY:
        return globalRow == globalCol
            ? 1.0
            : 0.0;

    case InitType::DETERMINISTIC:
        return static_cast<double>(
            globalRow * K_ + globalCol
        ) / static_cast<double>(
            M_ * K_
        );

    default:
        return valueAt(
            globalRow,
            globalCol,
            baseSeed_
        );
    }
}

// Generate B value for requested global coordinate
double DGEMMCannonPlace::generateB(
    const int globalRow,
    const int globalCol
) const {
    // Keep B generation separate so rectangular dimensions are handled
    switch (initType_) {
    case InitType::ONES:
        return 1.0;

    case InitType::IDENTITY:
        return globalRow == globalCol
            ? 1.0
            : 0.0;

    case InitType::DETERMINISTIC:
        return static_cast<double>(
            globalRow * N_ + globalCol
        ) / static_cast<double>(
            K_ * N_
        );

    default:
        return valueAt(
            globalRow,
            globalCol,
            bSeed_
        );
    }
}

// Allocate fixed size MASS shift message for A and B panels
void DGEMMCannonPlace::allocateOutMessage()
{
    if (outMessage) {
        std::free(
            outMessage
        );

        outMessage = nullptr;
    }

    // Store A first and B immediately after it in payload doubles
    aPanelOffsetDoubles_ = 0;
    bPanelOffsetDoubles_ = aPanelMaxDoubles_;

    const std::size_t doublesTotal =
        static_cast<std::size_t>(
            aPanelMaxDoubles_
        ) + static_cast<std::size_t>(
            bPanelMaxDoubles_
        );

    const std::size_t bytes =
        static_cast<std::size_t>(
            MASS_DGEMM_CANNON_HEADER_BYTES
        ) + doublesTotal * sizeof(
            double
        );

    // MASS exchange uses outMessage as raw byte storage
    outMessage = std::malloc(
        bytes
    );

    if (!outMessage) {
        std::cerr
            << "DGEMMCannonPlace::init malloc failed for outMessage ("
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

    inMessage_size = static_cast<int>(
        sizeof(
            double
        )
    );

    // Header stores sender tile coordinates and spare metadata slots
    int *header = static_cast<int *>(
        outMessage
    );

    header[0] = pi_;
    header[1] = pj_;
    header[2] = 0;
    header[3] = 0;
}

// Register left and up neighbors for Cannon panel shifts
void DGEMMCannonPlace::setupShiftNeighbors()
{
    if (P_ == 1) {
        return;
    }

    // MASS neighbor offsets are relative to current Place index
    std::vector<int *> offsets;

    offsets.reserve(
        2
    );

    // A shifts left across each Cannon row
    int *leftOffset = new int[2];

    leftOffset[0] = 0;
    leftOffset[1] = pj_ == 0
        ? P_ - 1
        : -1;

    offsets.push_back(
        leftOffset
    );

    // B shifts up across each Cannon column
    int *upOffset = new int[2];

    upOffset[0] = pi_ == 0
        ? P_ - 1
        : -1;

    upOffset[1] = 0;

    offsets.push_back(
        upOffset
    );

    addNeighbors(
        offsets
    );
}

// Initialize local tile ownership, buffers, and MASS shift layout
void *DGEMMCannonPlace::init(
    void *argument
) {
    // MASS Place index gives this tile coordinate in square grid
    pi_ = index[0];
    pj_ = index[1];
    P_ = size[0];

    // Driver passes benchmark dimensions and value generation settings
    if (argument) {
        const DGEMMCannonConfig *config =
            static_cast<const DGEMMCannonConfig *>(
                argument
            );

        M_ = config->M;
        N_ = config->N;
        K_ = config->K;
        baseSeed_ = config->baseSeed;
        bSeed_ = baseSeed_ ^ MASS_DGEMM_B_MATRIX_SEED_XOR;

        parseInitType(
            config->init_type
        );
    }
    else {
        M_ = 2048;
        N_ = 2048;
        K_ = 2048;
        baseSeed_ = 1;
        bSeed_ = baseSeed_ ^ MASS_DGEMM_B_MATRIX_SEED_XOR;
        initType_ = InitType::HASH;
    }

    // Split global C rows across Cannon grid rows
    const std::pair<int, int> rowRange = split1D(
        M_,
        P_,
        pi_
    );

    myRowStart_ = rowRange.first;
    myRowCount_ = rowRange.second;

    // Split global C columns across Cannon grid columns
    const std::pair<int, int> colRange = split1D(
        N_,
        P_,
        pj_
    );

    myColStart_ = colRange.first;
    myColCount_ = colRange.second;

    // Find largest shared dimension block for padded panel storage
    int maxK = 0;

    for (int step = 0; step < P_; step++) {
        const int kWidth = split1D(
            K_,
            P_,
            step
        ).second;

        if (kWidth > maxK) {
            maxK = kWidth;
        }
    }

    maxKWidth_ = maxK;

    maxRowCount_ = split1D(
        M_,
        P_,
        0
    ).second;

    maxColCount_ = split1D(
        N_,
        P_,
        0
    ).second;

    // Padded message sizes use largest possible local panel shapes
    aPanelMaxDoubles_ = maxRowCount_ * maxKWidth_;
    bPanelMaxDoubles_ = maxKWidth_ * maxColCount_;

    // Panels are padded to max size so each MASS message has fixed size
    A_panel_.assign(
        static_cast<std::size_t>(
            aPanelMaxDoubles_
        ),
        0.0
    );

    B_panel_.assign(
        static_cast<std::size_t>(
            bPanelMaxDoubles_
        ),
        0.0
    );

    A_recv_.assign(
        A_panel_.size(),
        0.0
    );

    B_recv_.assign(
        B_panel_.size(),
        0.0
    );

    C_tile_.assign(
        static_cast<std::size_t>(
            myRowCount_
        ) * myColCount_,
        0.0
    );

    // Initial skew fill happens inside timed region
    allocateOutMessage();
    setupShiftNeighbors();

    currentStep_ = 0;

    return nullptr;
}

// Fill initial skewed A and B panels inside timed region
void *DGEMMCannonPlace::initialFill(
    void *argument
) {
    (void)argument;

    // Initial skew chooses the k block already aligned to this tile
    const int initialBlock = (
        pi_ + pj_
    ) % P_;

    const std::pair<int, int> kRange = split1D(
        K_,
        P_,
        initialBlock
    );

    const int kStart = kRange.first;
    const int kWidth = kRange.second;

    // kStart maps local offsets back to global shared dimension

    // Clear padded panel tails before filling active values
    std::fill(
        A_panel_.begin(),
        A_panel_.end(),
        0.0
    );

    std::fill(
        B_panel_.begin(),
        B_panel_.end(),
        0.0
    );

    // A panel stores local C rows and current shared dimension block
    for (int localRow = 0; localRow < myRowCount_; localRow++) {
        const int globalRow = myRowStart_ + localRow;

        double *row = &A_panel_[
            static_cast<std::size_t>(
                localRow
            ) * maxKWidth_
        ];

        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            row[kOffset] = generateA(
                globalRow,
                kStart + kOffset
            );
        }
    }

    // B panel stores current shared dimension block and local C columns
    for (int kOffset = 0; kOffset < kWidth; kOffset++) {
        const int globalK = kStart + kOffset;

        double *row = &B_panel_[
            static_cast<std::size_t>(
                kOffset
            ) * myColCount_
        ];

        for (int localCol = 0; localCol < myColCount_; localCol++) {
            row[localCol] = generateB(
                globalK,
                myColStart_ + localCol
            );
        }
    }

    // Step 0 is consumed by multiplyStep0
    currentStep_ = 0;

    return nullptr;
}

// Multiply initial skewed panels before shift pipeline starts
void *DGEMMCannonPlace::multiplyStep0(
    void *argument
) {
    (void)argument;

    const std::pair<int, int> kRange = split1D(
        K_,
        P_,
        (
            pi_ + pj_
        ) % P_
    );

    const int kWidth = kRange.second;
    const int innerCol = myColCount_;

    // innerCol is local C column count and B panel row stride
    // Raw pointers keep inner multiply loops simple and fast
    double *__restrict__ cBase = C_tile_.data();
    const double *__restrict__ aBase = A_panel_.data();
    const double *__restrict__ bBase = B_panel_.data();

    int localRow = 0;

    // Reuse each B row across 4 C rows
    for (
        ;
        localRow + MICRO_KERNEL_ROWS - 1 < myRowCount_;
        localRow += MICRO_KERNEL_ROWS
    ) {
        // Get 4 output rows from local C block
        double *__restrict__ cRow0 = cBase +
            static_cast<std::size_t>(
                localRow + 0
            ) * innerCol;

        double *__restrict__ cRow1 = cBase +
            static_cast<std::size_t>(
                localRow + 1
            ) * innerCol;

        double *__restrict__ cRow2 = cBase +
            static_cast<std::size_t>(
                localRow + 2
            ) * innerCol;

        double *__restrict__ cRow3 = cBase +
            static_cast<std::size_t>(
                localRow + 3
            ) * innerCol;

        // Get 4 input rows from local A panel
        const double *__restrict__ aRow0 = aBase +
            static_cast<std::size_t>(
                localRow + 0
            ) * maxKWidth_;

        const double *__restrict__ aRow1 = aBase +
            static_cast<std::size_t>(
                localRow + 1
            ) * maxKWidth_;

        const double *__restrict__ aRow2 = aBase +
            static_cast<std::size_t>(
                localRow + 2
            ) * maxKWidth_;

        const double *__restrict__ aRow3 = aBase +
            static_cast<std::size_t>(
                localRow + 3
            ) * maxKWidth_;

        // Walk through current shared dimension block
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            // Load A scalars once for reuse for all local B columns
            const double aValue0 = aRow0[kOffset];
            const double aValue1 = aRow1[kOffset];
            const double aValue2 = aRow2[kOffset];
            const double aValue3 = aRow3[kOffset];

            // Current B row is reused by all 4 C rows
            const double *__restrict__ bRow = bBase +
                static_cast<std::size_t>(
                    kOffset
                ) * innerCol;

            // Update one column position in 4 C rows
            for (int localCol = 0; localCol < innerCol; localCol++) {
                const double bValue = bRow[localCol];

                cRow0[localCol] += aValue0 * bValue;
                cRow1[localCol] += aValue1 * bValue;
                cRow2[localCol] += aValue2 * bValue;
                cRow3[localCol] += aValue3 * bValue;
            }
        }
    }

    // Handle remaining rows that do not fill a 4 row block
    for (; localRow < myRowCount_; localRow++) {
        // Get output row for tail computation
        double *__restrict__ cRow = cBase +
            static_cast<std::size_t>(
                localRow
            ) * innerCol;

        // Get matching A row for scalar tail computation
        const double *__restrict__ aRow = aBase +
            static_cast<std::size_t>(
                localRow
            ) * maxKWidth_;

        // Accumulate current k block into 1 C row
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            const double aValue = aRow[kOffset];

            // Current B row contributes to all local C columns
            const double *__restrict__ bRow = bBase +
                static_cast<std::size_t>(
                    kOffset
                ) * innerCol;

            for (int localCol = 0; localCol < innerCol; localCol++) {
                cRow[localCol] += aValue * bRow[localCol];
            }
        }
    }

    currentStep_ = 1;

    return nullptr;
}

// Pack current panels so MASS can send them to shift neighbors
void *DGEMMCannonPlace::publishShift(
    void *argument
) {
    (void)argument;

    // Header stores sender tile coordinates and spare metadata slots
    int *header = static_cast<int *>(
        outMessage
    );

    header[0] = pi_;
    header[1] = pj_;

    // Double payload starts after integer header bytes
    double *outDoubles = reinterpret_cast<double *>(
        static_cast<char *>(
            outMessage
        ) + MASS_DGEMM_CANNON_HEADER_BYTES
    );

    std::memcpy(
        outDoubles + aPanelOffsetDoubles_,
        A_panel_.data(),
        static_cast<std::size_t>(
            aPanelMaxDoubles_
        ) * sizeof(
            double
        )
    );

    std::memcpy(
        outDoubles + bPanelOffsetDoubles_,
        B_panel_.data(),
        static_cast<std::size_t>(
            bPanelMaxDoubles_
        ) * sizeof(
            double
        )
    );

    return nullptr;
}

// Receive shifted A or B panel based on sender tile position
void *DGEMMCannonPlace::recvShift(
    void *argument
) {
    if (!argument) {
        return new double(
            0.0
        );
    }

    // Decode sender tile to decide whether payload is A or B
    const int *header = static_cast<const int *>(
        argument
    );

    const int senderTileI = header[0];
    const int senderTileJ = header[1];

    // Payload layout matches publishShift packing order
    const double *payload = reinterpret_cast<const double *>(
        static_cast<const char *>(
            argument
        ) + MASS_DGEMM_CANNON_HEADER_BYTES
    );

    // Sender on right provides next A panel after left shift
    if (senderTileI == pi_ && senderTileJ == (
        (
            pj_ + 1
        ) % P_
    )) {
        std::memcpy(
            A_recv_.data(),
            payload + aPanelOffsetDoubles_,
            static_cast<std::size_t>(
                aPanelMaxDoubles_
            ) * sizeof(
                double
            )
        );
    }
    // Sender below provides next B panel after upward shift
    else if (senderTileI == (
        (
            pi_ + 1
        ) % P_
    ) && senderTileJ == pj_) {
        std::memcpy(
            B_recv_.data(),
            payload + bPanelOffsetDoubles_,
            static_cast<std::size_t>(
                bPanelMaxDoubles_
            ) * sizeof(
                double
            )
        );
    }

    return new double(
        0.0
    );
}

// Make received panels active for next multiply step
void *DGEMMCannonPlace::shiftFinish(
    void *argument
) {
    (void)argument;

    std::swap(
        A_panel_,
        A_recv_
    );

    std::swap(
        B_panel_,
        B_recv_
    );

    return nullptr;
}

// Multiply shifted panels for current Cannon step
void *DGEMMCannonPlace::multiplyOnly(
    void *argument
) {
    (void)argument;

    const std::pair<int, int> kRange = split1D(
        K_,
        P_,
        (
            pi_ + pj_ + currentStep_
        ) % P_
    );

    const int kWidth = kRange.second;
    const int innerCol = myColCount_;

    double *__restrict__ cBase = C_tile_.data();
    const double *__restrict__ aBase = A_panel_.data();
    const double *__restrict__ bBase = B_panel_.data();

    int localRow = 0;

    // Reuse each B row across 4 C rows
    for (
        ;
        localRow + MICRO_KERNEL_ROWS - 1 < myRowCount_;
        localRow += MICRO_KERNEL_ROWS
    ) {
        // Get 4 output rows from local C block
        double *__restrict__ cRow0 = cBase +
            static_cast<std::size_t>(
                localRow + 0
            ) * innerCol;

        double *__restrict__ cRow1 = cBase +
            static_cast<std::size_t>(
                localRow + 1
            ) * innerCol;

        double *__restrict__ cRow2 = cBase +
            static_cast<std::size_t>(
                localRow + 2
            ) * innerCol;

        double *__restrict__ cRow3 = cBase +
            static_cast<std::size_t>(
                localRow + 3
            ) * innerCol;

        // Get 4 input rows from local A panel
        const double *__restrict__ aRow0 = aBase +
            static_cast<std::size_t>(
                localRow + 0
            ) * maxKWidth_;

        const double *__restrict__ aRow1 = aBase +
            static_cast<std::size_t>(
                localRow + 1
            ) * maxKWidth_;

        const double *__restrict__ aRow2 = aBase +
            static_cast<std::size_t>(
                localRow + 2
            ) * maxKWidth_;

        const double *__restrict__ aRow3 = aBase +
            static_cast<std::size_t>(
                localRow + 3
            ) * maxKWidth_;

        // Walk through current shared dimension block
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            // Load A scalars once for reuse for all local B columns
            const double aValue0 = aRow0[kOffset];
            const double aValue1 = aRow1[kOffset];
            const double aValue2 = aRow2[kOffset];
            const double aValue3 = aRow3[kOffset];

            // Current B row is reused by all 4 C rows
            const double *__restrict__ bRow = bBase +
                static_cast<std::size_t>(
                    kOffset
                ) * innerCol;

            // Update one column position in 4 C rows
            for (int localCol = 0; localCol < innerCol; localCol++) {
                const double bValue = bRow[localCol];

                cRow0[localCol] += aValue0 * bValue;
                cRow1[localCol] += aValue1 * bValue;
                cRow2[localCol] += aValue2 * bValue;
                cRow3[localCol] += aValue3 * bValue;
            }
        }
    }

    // Handle remaining rows that do not fill a 4 row block
    for (; localRow < myRowCount_; localRow++) {
        // Get output row for tail computation
        double *__restrict__ cRow = cBase +
            static_cast<std::size_t>(
                localRow
            ) * innerCol;

        // Get matching A row for scalar tail computation
        const double *__restrict__ aRow = aBase +
            static_cast<std::size_t>(
                localRow
            ) * maxKWidth_;

        // Accumulate current k block into 1 C row
        for (int kOffset = 0; kOffset < kWidth; kOffset++) {
            const double aValue = aRow[kOffset];

            // Current B row contributes to all local C columns
            const double *__restrict__ bRow = bBase +
                static_cast<std::size_t>(
                    kOffset
                ) * innerCol;

            for (int localCol = 0; localCol < innerCol; localCol++) {
                cRow[localCol] += aValue * bRow[localCol];
            }
        }
    }

    // Advance logical Cannon step after contribution is accumulated
    currentStep_ = (
        currentStep_ + 1
    ) % P_;

    return nullptr;
}

// Reset C tile outside timed region before each benchmark run
void *DGEMMCannonPlace::resetC(
    void *argument
) {
    (void)argument;

    std::fill(
        C_tile_.begin(),
        C_tile_.end(),
        0.0
    );

    currentStep_ = 0;

    return nullptr;
}

// Return local checksum so driver can sum all Place results
void *DGEMMCannonPlace::getChecksum(
    void *argument
) {
    (void)argument;

    // MASS callAll returns 1 allocated value per Place
    double *out = new double;

    *out = 0.0;

    for (const double value : C_tile_) {
        *out += value;
    }

    return out;
}
