/****************************************************************************
 *                                                                          *
 * ICPADS26 Benchmarks Collection                                           *
 *                                                                          *
 * Benchmark: DGEMM (Cannon's Algorithm)                                    *
 * Library: MASS                                                            *
 *                                                                          *
 * Author: Ahmed Bera Pay                                                   *
 * Faculty Advisor: Munehiro Fukuda                                         *
 * Code Finalization: Kyryll Kotyk                                          *
 *                                                                          *
 ****************************************************************************/

#ifndef MASS_DGEMM_CANNON_PLACE_
#define MASS_DGEMM_CANNON_PLACE_

#include "MASS_base.h"
#include "MethodRegistry.h"
#include "Place.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>
#include <vector>

#define MASS_DGEMM_ROW_HASH_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define MASS_DGEMM_COL_HASH_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define MASS_DGEMM_FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define MASS_DGEMM_SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL
#define MASS_DGEMM_B_MATRIX_SEED_XOR 0x9e3779b9ULL
#define MASS_DGEMM_DOUBLE_SCALE (1.0 / (1ULL << 53))
#define MASS_DGEMM_CANNON_HEADER_INTS 4
#define MASS_DGEMM_CANNON_HEADER_BYTES \
    (MASS_DGEMM_CANNON_HEADER_INTS * static_cast<int>(sizeof(int)))

struct DGEMMCannonConfig
{
    int M;
    int N;
    int K;
    int runs;
    uint64_t baseSeed;
    char init_type[32];
};

class DGEMMCannonPlace : public Place
{
public:
    explicit DGEMMCannonPlace(
        void *argument
    );

    ~DGEMMCannonPlace() override;

    MASS_DISPATCH_TABLE(
        DGEMMCannonPlace,
        MASS_METHOD(
            DGEMMCannonPlace,
            init
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            initialFill
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            multiplyStep0
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            publishShift
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            recvShift
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            shiftFinish
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            multiplyOnly
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            resetC
        ),
        MASS_METHOD(
            DGEMMCannonPlace,
            getChecksum
        )
    )

    void *init(
        void *argument
    );

    void *initialFill(
        void *argument
    );

    void *multiplyStep0(
        void *argument
    );

    void *publishShift(
        void *argument
    );

    void *recvShift(
        void *argument
    );

    void *shiftFinish(
        void *argument
    );

    void *multiplyOnly(
        void *argument
    );

    void *resetC(
        void *argument
    );

    void *getChecksum(
        void *argument
    );

private:
    enum class InitType : int
    {
        HASH = 0,
        DETERMINISTIC = 1,
        IDENTITY = 2,
        ONES = 3
    };

    static std::pair<int, int> split1D(
        const int totalSize,
        const int numParts,
        const int index
    );

    static double valueAt(
        const int rowIndex,
        const int colIndex,
        const uint64_t seed
    );

    void parseInitType(
        const char *raw
    );

    double generateA(
        const int globalRow,
        const int globalCol
    ) const;

    double generateB(
        const int globalRow,
        const int globalCol
    ) const;

    void allocateOutMessage();

    void setupShiftNeighbors();

    // Global matrix dimensions
    int M_;
    int N_;
    int K_;

    // Square Cannon grid position
    int P_;
    int pi_;
    int pj_;

    // Local C tile ownership
    int myRowStart_;
    int myRowCount_;
    int myColStart_;
    int myColCount_;

    // Padded panel dimensions used for fixed shift message sizes
    int maxKWidth_;
    int maxRowCount_;
    int maxColCount_;
    int aPanelMaxDoubles_;
    int bPanelMaxDoubles_;

    // Matrix value generation state
    uint64_t baseSeed_;
    uint64_t bSeed_;
    InitType initType_;

    // Active panels and shift staging buffers
    std::vector<double> A_panel_;
    std::vector<double> B_panel_;
    std::vector<double> A_recv_;
    std::vector<double> B_recv_;

    // Running local C tile accumulator
    std::vector<double> C_tile_;

    // Payload offsets inside outMessage after integer header
    int aPanelOffsetDoubles_;
    int bPanelOffsetDoubles_;

    // Current Cannon step after initial skewed multiply
    int currentStep_;
};

#endif
