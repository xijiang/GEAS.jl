using JLD2, DataFrames #, Distributions, Serialization, Dates, Statistics, DataFrames
import GEAS:create_storage, sim_QTL
function t_2021_08_06()
    @load "dat/run/base.jld2" base
    nqtl = 100
    #qtl = sim_QTL(base, nqtl, nqtl)
    par = begin
        Parameters = Dict(
            :nSire=> 100,
            :nDam => 200,   # male:female = 1:2
            :nSib => 40,
            :nC7e => 30,    # number of sib to be challenged
            :nG8n => 10,    # number of generations
            :nQTL => [nqtl, nqtl],
            :h²   => [.5, .5],
            :p8e  => .5,    # death rate after challenge
            :e19e => 1.,    # edit_successsful_rate
            :w4t  => 1,     # weight on the binary trait EBV
            :nk3n => 10,    # number of known QTL on binary trait
            :dd   => 0.01,  # A small value to be added to diagonals
            #:log  => slog
        )
        (; Parameters...)       # named tuple.  contents as above
    end
    create_storage(base, par)
end
