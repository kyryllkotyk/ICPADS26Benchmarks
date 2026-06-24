#include "HPXdgemmSUMMA.h"

int main(
    int argc,
    char* argv[]
) {
    namespace options = hpx::program_options;

    options::options_description description(
        "DGEMM HPX SUMMA options"
    );

    // Register matrix, process grid, and run count options
    description.add_options()
        (
            "M",
            options::value<int>()->default_value(
                2048
            ),
            "rows of A/C"
        )
        (
            "N",
            options::value<int>()->default_value(
                2048
            ),
            "cols of B/C"
        )
        (
            "K",
            options::value<int>()->default_value(
                2048
            ),
            "cols of A / rows of B"
        )
        (
            "Px",
            options::value<int>()->default_value(
                1
            ),
            "process grid rows"
        )
        (
            "Py",
            options::value<int>()->default_value(
                1
            ),
            "process grid cols"
        )
        (
            "runs",
            options::value<int>()->default_value(
                1
            ),
            "number of benchmark runs"
        )
        (
            "seed",
            options::value<std::uint64_t>()->default_value(
                1
            ),
            "RNG base seed"
        )
        (
            "init_type",
            options::value<std::string>()->default_value(
                "hash"
            ),
            "matrix init: hash | deterministic | identity | ones"
        )
        (
            "block_size",
            options::value<int>()->default_value(
                128
            ),
            "not used, kept for old run scripts"
        );

    hpx::init_params initParams;

    initParams.desc_cmdline = description;

    return hpx::init(
        argc,
        argv,
        initParams
    );
}
