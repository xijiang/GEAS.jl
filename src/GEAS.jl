module GEAS

using LinearAlgebra, Dates, HTTP, DataFrames, Random, Distributions,
    Serialization, Statistics, StatsPlots, JLD2, LaTeXStrings
import Octavian: matmul, matmul!

const prj_dir = begin
    t = splitdir(pathof(GEAS))[1]
    l = findlast('/', t) - 1
    t[1:l]
end

dat_dir = joinpath(prj_dir, "dat")
bin_dir = joinpath(prj_dir, "bin")
plink   = joinpath(bin_dir, "plink")
beagle  = joinpath(bin_dir, "beagle.jar")
macs    = joinpath(bin_dir, "macs")
qmsim   = joinpath(bin_dir, "qmsim")

include("breeding.jl")
include("b-Update.jl")
include("dropping.jl")
include("evaluation.jl")
include("gedit.jl")
include("genotypes.jl")
include("inbreeding.jl")
include("phenotypes.jl")
include("sim-base.jl")
include("sim-QTL.jl")
include("structs.jl")
include("summarize.jl")
include("workflow.jl")
include("check-par.jl")
include("utils.jl")
include("gwas.jl")
include("bayes-alphabet.jl")

copyright()

end # module
