module GEAS

using LinearAlgebra, Dates, HTTP, DataFrames, Random, Distributions,
    Serialization, Statistics
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
include("estimation.jl")

copyright()

end # module
