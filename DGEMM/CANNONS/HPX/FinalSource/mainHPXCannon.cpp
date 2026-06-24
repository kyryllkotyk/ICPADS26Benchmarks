/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: DGEMM (Cannon's Algorithm)                                      *
 * Library: HPX                                                               *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#include "HPXdgemmCannon.h"

int main(
    int argc,
    char* argv[]
) {
    namespace programOptions = hpx::program_options;

    // Define command line options consumed by hpx_main
    programOptions::options_description description(
        "DGEMM HPX Cannon options"
    );

    // Register matrix, grid, run, and initialization options
    description.add_options()
        (
            "M",
            programOptions::value<int>()->default_value(
                2048
            ),
            "rows of A/C"
        )
        (
            "N",
            programOptions::value<int>()->default_value(
                2048
            ),
            "cols of B/C"
        )
        (
            "K",
            programOptions::value<int>()->default_value(
                2048
            ),
            "cols of A / rows of B"
        )
        (
            "P",
            programOptions::value<int>()->default_value(
                0
            ),
            "grid dimension (P*P localities)"
        )
        (
            "Px",
            programOptions::value<int>()->default_value(
                0
            ),
            "alias of P"
        )
        (
            "Py",
            programOptions::value<int>()->default_value(
                0
            ),
            "must equal Px"
        )
        (
            "runs",
            programOptions::value<int>()->default_value(
                1
            ),
            "benchmark runs"
        )
        (
            "seed",
            programOptions::value<std::uint64_t>()->default_value(
                1
            ),
            "base seed"
        )
        (
            "init_type",
            programOptions::value<std::string>()->default_value(
                "hash"
            ),
            "hash | deterministic | identity | ones"
        )
        (
            "block_size",
            programOptions::value<int>()->default_value(
                128
            ),
            "not used"
        );

    hpx::init_params initParams;

    // Pass command line description into HPX runtime startup
    initParams.desc_cmdline = description;

    // Start HPX runtime and transfer control to hpx_main
    return hpx::init(
        argc,
        argv,
        initParams
    );
}
