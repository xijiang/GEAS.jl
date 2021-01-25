"""
    function workflow()
---
This records the framework of the package.
It can also serve as a template pipeline.
"""
function workflow(; debug = true)
    # Par = Parameter(500, 600, 100, .5, .5, .5)
    prd = Trait(500, .5)        # production trait
    dsr = Trait(500, .5)        # disease resisitance trait
    nbs, nsb = 600, 100         # base population size, full sibship size
    @info "run function `Update()` to update binaries if necessary"
    if debug
        @info "Running in debug mode"
        @info join(["",
                    "STEP I: Preparing data",
                    "  1. some real data with ID masked",
                    "  2. simulation data",
                    "  ToDo: simulation of base"], "\n")
        @info "Phase the genotypes of 600 ID with beagle"
        # read back real data
        base = deserialize(joinpath(dat_dir, "run/ns.ser")) # to save time
        qtls = sim_QTL(base, prd.nqtl, dsr.nqtl)
        ped  = random_mate(nbs, nsb)
        r    = haldane(base[:pos])
        goff = gdrop(base[:hap], ped, r)

        @info "Start to breeding"
    else
        @info "Running in release mode"
        @time gt600()           # -> c.vcf.gz
        @time base = vcf2dic(joinpath(dat_dir, "run/c.vcf.gz"))
        @time serialize(joinpath(dat_dir, "run/ns.ser"), base)
    end
end
