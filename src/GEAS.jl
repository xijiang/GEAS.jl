module GEAS

using LinearAlgebra, Dates, HTTP, DataFrames, Random, Distributions,
    Serialization, Statistics, StatsPlots, JLD2
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

include("workflow.jl")
include("genotypes.jl")
include("b_Update.jl")
include("simQTL.jl")
include("dropping.jl")
include("structs.jl")
include("phenotypes.jl")
include("breeding.jl")
include("evaluation.jl")

# Test files
# all functions in files below have a prefix `test_`.
include("test-breeding.jl")
include("test-evaluation.jl")

copyright()

end # module
