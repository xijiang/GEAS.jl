module GEAS

using LinearAlgebra, Dates, HTTP, DataFrames, Random, Distributions,
    Serialization, Statistics, StatsPlots, JLD2
import Octavian: matmul, matmul!
export workflow

prj_dir = begin
    t = splitdir(pathof(GEAS))[1]
    l = findlast('/', t) - 1
    t[1:l]
end

dat_dir = joinpath(prj_dir, "dat")
bin_dir = joinpath(prj_dir, "bin")
plink   = joinpath(bin_dir, "plink")
beagle  = joinpath(bin_dir, "beagle.jar")
macs    = joinpath(bin_dir, "macs")

include("workflow.jl")
include("genotypes.jl")
include("b_Update.jl")
include("simQTL.jl")
include("dropping.jl")
include("structs.jl")
include("phenotypes.jl")
include("breeding.jl")
include("evaluation.jl")
include("sim-base.jl")
include("gedit.jl")

# Test files
# all functions in files below have a prefix `test_`.
# These can be excluded in the release version
include("tst/breeding.jl")
include("tst/evaluation.jl")
include("tst/dropping.jl")
include("tst/editing.jl")

copyright()

end # module
