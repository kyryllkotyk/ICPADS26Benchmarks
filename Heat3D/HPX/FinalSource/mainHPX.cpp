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

int main(
    int argc,
    char* argv[]
) {
    namespace programOptions = hpx::program_options;

    // Register command line options consumed by hpx_main
    programOptions::options_description description(
        "Heat3D HPX options"
    );

    description.add_options()
        (
            "gridX",
            programOptions::value<int>()->default_value(
                128
            ),
            "global grid X"
        )
        (
            "gridY",
            programOptions::value<int>()->default_value(
                128
            ),
            "global grid Y"
        )
        (
            "gridZ",
            programOptions::value<int>()->default_value(
                128
            ),
            "global grid Z"
        )
        (
            "decompX",
            programOptions::value<int>()->default_value(
                1
            ),
            "decomposition X"
        )
        (
            "decompY",
            programOptions::value<int>()->default_value(
                1
            ),
            "decomposition Y"
        )
        (
            "decompZ",
            programOptions::value<int>()->default_value(
                1
            ),
            "decomposition Z"
        )
        (
            "timesteps",
            programOptions::value<int>()->default_value(
                100
            ),
            "timesteps"
        )
        (
            "runs",
            programOptions::value<int>()->default_value(
                1
            ),
            "benchmark runs"
        )
        (
            "debug",
            programOptions::value<int>()->default_value(
                0
            ),
            "debug output"
        )
        (
            "seed",
            programOptions::value<std::uint64_t>()->default_value(
                1
            ),
            "init seed"
        )
        (
            "initMin",
            programOptions::value<double>()->default_value(
                0.0
            ),
            "minimum initial temperature"
        )
        (
            "initMax",
            programOptions::value<double>()->default_value(
                1.0
            ),
            "maximum initial temperature"
        )
        (
            "alpha",
            programOptions::value<double>()->default_value(
                0.5
            ),
            "center weight"
        )
        (
            "beta",
            programOptions::value<double>()->default_value(
                0.1
            ),
            "neighbor weight"
        );

    hpx::init_params initParameters;

    // Keep launcher setup here and benchmark logic in hpx_main
    initParameters.desc_cmdline = description;

    return hpx::init(
        argc,
        argv,
        initParameters
    );
}
