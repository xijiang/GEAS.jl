"""
    function workflow()
---
This records the framework of the package.
"""
function workflow(; debug = true)
    # Par = Parameter(500, 600, 100, .5, .5, .5)
    prd = Trait(500, .5)
    dsr = Trait(500, .5)
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
        ############################################################
        # Below are a few steps to organize some real data from Nofima
        # @time gt600()           # -> c.vcf.gz
        # base = vcf2dic(joinpath(dat_dir, "run/c.vcf.gz"))
        # serialize(joinpath(dat_dir, "run/ns.ser"), base)
        ############################################################
        # read back real data
        base = deserialize(joinpath(dat_dir, "run/ns.ser")) # to save time
        Qinf = sim_QTL(base, prd.nqtl, dsr.nqtl)
        ped = random_mate(nbs, nsb)
        r = haldane(base[:pos])
        goff = gdrop(base[:hap], ped, r)
        TBV(goff, Qinf[1])
    else
        @info "Running in release mode"
    end
end
