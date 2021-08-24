using JLD2, DataFrames
#, Distributions, Serialization, Dates, Statistics, DataFrames
import GEAS:breeding, sim_QTL, summarize

function vnqtl(nqtl)
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
end

function t_2021_08_14()
    @load "dat/run/base.jld2" base
    nqtl = 500
    par = vnqtl(nqtl)
    nrpt = 20

    rst_fixed = DataFrame()
    for i in 1:nrpt
        @info "====== Repeat with fixed QTL, $i of $nrpt ======"
        qtl = sim_QTL(base, nqtl, nqtl)
        ped, snp = breeding(base, par, qtl, fixed=true)
        df = summarize(ped)
        df.rpt = ones(Int32, nrow(df)) .* i
        append!(rst_fixed, df)
    end
    
    rst_edit = DataFrame()
    for i in 1:nrpt
        @info "====== Repeat with editing, $i of $nrpt ======"
        qtl = sim_QTL(base, nqtl, nqtl)
        ped, snp = breeding(base, par, qtl, edit=true)
        df = summarize(ped)
        df.rpt = ones(Int32, nrow(df)) .* i
        append!(rst_edit, df)
    end
    rst_fixed, rst_edit
end
