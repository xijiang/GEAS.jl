"""
    function workflow()
---
This records the framework of the package.
It can also serve as a template pipeline.
"""
function workflow(; debug = true)
    copyright()
    
    if debug
        @info join(["Debug mode",
                    "      ==========",
                    "",
                    "-- run `GEAS.Update()` to update binaries if necessary"],
                   "\n")
        # read back real data
        @load "dat/run/base.jld2"                           # =>base
        #base = deserialize(joinpath(dat_dir, "run/ns.ser")) # to save time
        nQTL = [1000, 1000]
        qtl = sim_QTL(base, nQTL...)

        println("\n")
        for weight in [3, 6, 9, 12]
            for nsib in [10, 20, 40]
                test_2_breeding(base, qtl, nsib, weight)
            end
        end
    else
        @info join(["Running in release mode",
                    "STEP I: Preparing base",
                    "  1. base population",
                    "     - some real Nofima data",
                    "     - simulate with `macs`",
                    "  2. If real data",
                    "     - remove duplicates",
                    "     - phase with `beagle.jar`",
                    "  3. serialization to ease later load"],
                   "\n")
        @time gt600()           # -> c.vcf.gz
        @time base = vcf2dic(joinpath(dat_dir, "run/c.vcf.gz"))
        @time serialize(joinpath(dat_dir, "run/ns.ser"), base)
        @save joinpath(dat_dir, "run/base.jld2") base
        @info join(["",
                    "STEP II: The breeding program",
                    "-----------------------------"],
                   "\n")
        @time breeding_program(base, par)
    end
end
