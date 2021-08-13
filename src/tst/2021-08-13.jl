using JLD2, DataFrames
#, Distributions, Serialization, Dates, Statistics, DataFrames
import GEAS:breeding, sim_QTL, summarize
function t_2021_08_06()
    @load "dat/run/base.jld2" base
    nqtl = 100
    par = begin
        Parameters = Dict(
            :nSire=> 100,
            :nDam => 200,   # male:female = 1:2
            :nSib => 40,
            :nC7e => 30,    # number of sib to be challenged
            :nG8n => 10,    # number of generations
            :nQTL => [nqtl, nqtl],
            :h²   => [.5, .5],
            :p8e  => .5,    # percentage of death after challenge
            :e19e => 1.,    # edit_successsful_rate
            :w4t  => 1,     # weight on the binary trait EBV
            :nK3n => 10,    # number of known QTL on binary trait
            :dd   => 0.01,  # A small value to be added to diagonals
            #:log  => slog
        )
        (; Parameters...)       # named tuple.  contents as above
    end
    #=
    Below is to test selection on phenotypes. To achieve this, insert
    prd[(g8n=ig,)][:, :val] = prd[(g8n=ig,)][:, :p7e]
    below calc_idx in function breeding
    =#
    nrpt = 20
    rst = DataFrame()
    for i in 1:nrpt
        @info "====== Repeat $i of $nrpt ======"
        qtl = sim_QTL(base, nqtl, nqtl)
        ped, snp = breeding(base, par, qtl, fixed=true)
        df = summarize(ped)
        df.rpt = ones(Int32, nrow(df)) .* i
        append!(rst, df)
    end
    rst
end
