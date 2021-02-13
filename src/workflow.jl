"""
    function workflow()
---
This records the framework of the package.
It can also serve as a template pipeline.
"""
function workflow(; debug = true)
    copyright()
    # It is better to setup simulation scenarios with a dictionary
    # Then pass it as a named tuple, so that we can conveniantly use `x.y`
    # using a dictionary also makes it easy to modify parameter structure.
    # the keys in this dictionry must exist to setup a simulation
    Parameters = Dict(:nSire=> 125,
                      :nDam => 250, # male:female = 1:2
                      :nSib => 100,
                      :nC7e => 40, # number of sib to be challenged
                      :nG8n => 10, # number of generations
                      :nQTL => [1000, 1000],
                      :h²   => [.5, .5],
                      :p8e  => .5, # percentage of affected in the binary trait
                      )
    par = (; Parameters...)     # named tuple.  contents as above
    
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
