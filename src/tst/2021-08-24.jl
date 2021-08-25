using JLD2, DataFrames, Serialization
#, Distributions, Serialization, Dates, Statistics, DataFrames
import GEAS:breeding, sim_QTL, summarize, create_storage, simple_breeding

function gpar(nqtl, edit, fixed)
    par = begin
        Parameters = Dict(
            :nSire=> 100,
            :nDam => 200,   # male:female = 1:2
            :nSib => 40,
            :nC7e => 30,    # number of sib to be challenged
            :nG8n => 11,    # number of generations
            :nQTL => [nqtl, nqtl],
            :h²   => [.5, .5],
            :p8e  => .5,    # percentage of death after challenge
            :e19e => 1.,    # edit_successsful_rate
            :w4t  => 1,     # weight on the binary trait EBV
            :nK3n => 10,    # number of known QTL on binary trait
            :dd   => 0.01,  # A small value to be added to diagonals
            :edit => edit,
            :fix  => fixed,
            #:log  => slog
        )
        (; Parameters...)       # named tuple.  contents as above
    end
end

function t_2021_08_24()
    isdir("dat/tmp") || mkpath("dat/tmp")
    @load "dat/run/base.jld2" base
    par = gpar(100, false, false)
    qtl = sim_QTL(base, 100, 100)
    #@info "Creating storage and 1st generation"
    #snp, ped = create_storage(base, par, qtl)
    #serialize("dat/tmp/snp.ser", snp)
    #serialize("dat/tmp/ped.ser", ped)
    @info "===== Breeding generation 2-$(par.nG8n) ====="
    #breeding(ped, snp, base, qtl, par)
    snp = deserialize("dat/tmp/snp.ser")
    ped = deserialize("dat/tmp/ped.ser")
    simple_breeding(ped, snp, base, qtl, par, 3)
    snp, ped
end

"""
    function t_one_trait()
---
Test selection on `TBV`, `phenotype`, and `SNP-BLUP` on production trait only.
"""
function t_one_trait()
    nrpt = 1
    isdir("dat/tmp") || mkpath("dat/tmp")
    @load "dat/run/base.jld2" base
    nqtl = 500
    par = gpar(nqtl, false, false)
    mtd = ["TBV", "phenotype", "EBV"]
    
    rst = DataFrame()
    for ir in 1:nrpt
        qtl = sim_QTL(base, par.nQTL...)
        @info "Creating storage and 1st generation, repeat $ir of $nrpt"
        snp, ped = create_storage(base, par, qtl)
        serialize("dat/tmp/snp.ser", snp)
        serialize("dat/tmp/ped.ser", ped)
        @info "===== Breeding generation 2-$(par.nG8n), on TBV ====="
        for op in 1:3
            println("----- Selection on $(mtd[op]) -----")
            snp = deserialize("dat/tmp/snp.ser")
            ped = deserialize("dat/tmp/ped.ser")
            simple_breeding(ped, snp, base, qtl, par, op)
        end
        df = summarize(ped.prd, snp, qtl[2])
        df.method = ones(Int, nrow(df)) .* op
        df.repeat = ones(Int, Nrow(df)) .* ir
        append!(rst, df)
    end
    serialize("dat/rst/one-trait.ser", rst)
    rst
end


#=
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
=#
