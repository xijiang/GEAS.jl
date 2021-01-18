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
        @info join(["STEP I: Preparing data",
                    "  1. some real data with ID masked",
                    "  2. simulation data",
                    "  ToDo: simulation of base"], "\n")
        @info "Phase the genotypes of 600 ID with beagle"
        # @time gt600()           # -> c.vcf.gz
        # base = vcf2dic(joinpath(dat_dir, "run/c.vcf.gz"))
        # serialize(joinpath(dat_dir, "run/ns.ser"), base)
        base = deserize(joinpath(dat_dir, "run/ns.ser")) # to save time
        nqtl = 500
        QTL = sim_QTL(base, nqtl)
    else
        @info "Running in release mode"
    end
end
