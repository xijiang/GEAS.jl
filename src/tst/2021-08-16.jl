using JLD2, DataFrames, Serialization
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

function append_rst(rst, ped, nqtl, method, i)
    df = summarize(ped)
    ng = nrow(df)
    df.nqtl = ones(Int, ng) .* nqtl
    df.method = ones(Int, ng) .* method
    df.repeat = ones(Int, ng) .* i
    append!(rst, df)
end

function t_2021_08_16(file)
    @load "dat/run/base.jld2" base
    nrpt = 30
    rst = DataFrame()
    for nqtl in [100, 500]
        par = vnqtl(nqtl)
        for i in 1:nrpt
            @info "==== Repeat $i of $nrpt ===="
            qtl = sim_QTL(base, nqtl, nqtl)
            
            ped, snp = breeding(base, par, qtl)
            append_rst(rst, ped, nqtl, 1, i)
            
            ped, snp = breeding(base, par, qtl, fixed=true)
            append_rst(rst, ped, nqtl, 2, i)

            ped, snp = breeding(base, par, qtl, edit=true)
            append_rst(rst, ped, nqtl, 3, i)
        end
    end
    serialize(file, rst)
end
