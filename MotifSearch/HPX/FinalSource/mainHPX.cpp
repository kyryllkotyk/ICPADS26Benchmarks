#include "HPXmotifsearch.h"

int main(
    int argc,
    char* argv[]
) {
    namespace programOptions = hpx::program_options;

    // Register command line options consumed by hpx_main
    programOptions::options_description description(
        "Motif Search HPX options"
    );

    description.add_options()
        (
            "vertices",
            programOptions::value<int>()->default_value(
                0
            ),
            "number of vertices"
            )
        (
            "edge_prob,edge-prob",
            programOptions::value<double>()->default_value(
                0.0
            ),
            "edge probability for Erdos-Renyi generation"
            )
        (
            "seed",
            programOptions::value<unsigned>()->default_value(
                42
            ),
            "random seed"
            )
        (
            "motif",
            programOptions::value<std::string>()->default_value(
                "triangle"
            ),
            "motif name"
            )
        (
            "input",
            programOptions::value<std::string>()->default_value(
                ""
            ),
            "graph edge list file"
            )
        (
            "threads_per_locality,threads-per-locality",
            programOptions::value<int>()->default_value(
                0
            ),
            "HPX worker threads per locality"
            );

    hpx::init_params initParameters;

    // Keep launcher setup here and benchmark logic in hpx_main
    initParameters.desc_cmdline = description;

    // Start HPX runtime and hand execution to hpx_main
    return hpx::init(
        argc,
        argv,
        initParameters
    );
}