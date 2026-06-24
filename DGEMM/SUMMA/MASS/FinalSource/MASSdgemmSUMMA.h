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

#ifndef MASS_DGEMM_SUMMA_
#define MASS_DGEMM_SUMMA_

#include "DGEMMPlace.h"
#include "IterationConfig.h"
#include "MASS.h"

#include "../../include/dgemm_run_report.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#define DGEMM_PLACES_HANDLE 1

struct MassDgemmSummaConfig
{
    int matrixRows = 2048;
    int matrixCols = 2048;
    int sharedDimension = 2048;
    int processGridRows = 1;
    int processGridCols = 1;
    int runs = 1;
    uint64_t baseSeed = 1;
    std::string initType = "hash";
};

int gcd(
    const int leftValue,
    const int rightValue
);

void runMassDgemmSummaBenchmark(
    char* massArgs[4],
    const int processes,
    const int threads,
    const MassDgemmSummaConfig& config
);

#endif
