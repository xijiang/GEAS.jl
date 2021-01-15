"""
    function workflow()
---
This records the framework of the package.
"""
function workflow(;debug = true)
    @info "run function `binUpdate()` to update binaries if necessary"
    # binUpdate()
    if debug
        @info "Running in debug mode"
        @info join(["Preparing data",
                    "  1. some real data with ID masked",
                    "  2. simulation data"], "\n")
        @info "Phase the genotypes of 600 ID with beagle"
        gt600()
    else
        @info "Running in release mode"
    end
end
