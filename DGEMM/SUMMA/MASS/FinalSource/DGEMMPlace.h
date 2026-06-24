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

#ifndef MASS_DGEMM_PLACE_
#define MASS_DGEMM_PLACE_

#include "Place.h"
#include "MethodRegistry.h"
#include "MASS_base.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>
#include <vector>

#define DGEMM_MESSAGE_HEADER_INTS 4
#define DGEMM_MESSAGE_HEADER_BYTES \
    (DGEMM_MESSAGE_HEADER_INTS * static_cast<int>(sizeof(int)))
#define ROW_HASH_MULTIPLIER 0x9e3779b97f4a7c15ULL
#define COL_HASH_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define FIRST_MIX_MULTIPLIER 0xbf58476d1ce4e5b9ULL
#define SECOND_MIX_MULTIPLIER 0x94d049bb133111ebULL
#define B_MATRIX_SEED_XOR 0x9e3779b9ULL
#define DOUBLE_SCALE (1.0 / (1ULL << 53))
#define MICRO_KERNEL_ROWS 4
#define INIT_TYPE_NAME_BYTES 32

struct DGEMMConfig
{
    int matrixRows;
    int matrixCols;
    int sharedDimension;
    int runs;
    uint64_t baseSeed;
    char initType[INIT_TYPE_NAME_BYTES];
};

class DGEMMPlace : public Place
{
public:
    explicit DGEMMPlace(
        void* argument
    );

    ~DGEMMPlace() override;

    // Register Place methods for MASS string dispatch
    MASS_DISPATCH_TABLE(
        DGEMMPlace,
        MASS_METHOD(
            DGEMMPlace,
            init
        ),
        MASS_METHOD(
            DGEMMPlace,
            publishPanels
        ),
        MASS_METHOD(
            DGEMMPlace,
            recvPanels
        ),
        MASS_METHOD(
            DGEMMPlace,
            accumulate
        ),
        MASS_METHOD(
            DGEMMPlace,
            resetC
        ),
        MASS_METHOD(
            DGEMMPlace,
            getChecksum
        )
    )

    void* init(
        void* argument
    );

    void* publishPanels(
        void* argument
    );

    void* recvPanels(
        void* argument
    );

    void* accumulate(
        void* argument
    );

    void* resetC(
        void* argument
    );

    void* getChecksum(
        void* argument
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
        const int totalDimensionSize,
        const int parts,
        const int partIndex
    );

    static double valueAt(
        const int rowIndex,
        const int colIndex,
        const uint64_t seed
    );

    static int gcd(
        const int leftValue,
        const int rightValue
    );

    void parseInitType(
        const char* rawInitType
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

    void setupRowColNeighbors();

    int ownerColForStep(
        const int step
    ) const;

    int ownerRowForStep(
        const int step
    ) const;

    int matrixRows_;
    int matrixCols_;
    int sharedDimension_;
    int processGridRows_;
    int processGridCols_;
    int tileRow_;
    int tileCol_;
    int localRowStart_;
    int localRowCount_;
    int localColStart_;
    int localColCount_;
    int stepCount_;
    int maxKWidth_;
    int maxRowCount_;
    int maxColCount_;
    int aPanelMaxDoubles_;
    int bPanelMaxDoubles_;

    uint64_t baseSeed_;
    uint64_t bMatrixSeed_;
    InitType initType_;

    std::vector<double> aPanel_;
    std::vector<double> bPanel_;
    std::vector<double> cTile_;

    int aPanelOffsetDoubles_;
    int bPanelOffsetDoubles_;
    int currentStep_;
};

#endif
