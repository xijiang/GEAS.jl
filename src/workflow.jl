"""
    function workflow()
---
This records the framework of the package.
It can also serve as a template pipeline.
"""
function workflow(par; debug = true)
    #par = Breeding(375, 2, 10, 100, [500, 500], [.5, .5], 40, 0.)
    if debug
        @info join(["Debug mode",
                    "      ==========",
                    "",
                    "-- run `GEAS.Update()` to update binaries if necessary"],
                   "\n")
        # read back real data
        base = deserialize(joinpath(dat_dir, "run/ns.ser")) # to save time

        println("\n")
        @info join(["",
                    "STEP II: The breeding program",
                    "-----------------------------"],
                   "\n")
        breeding_program(base, par)
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
    end
end
