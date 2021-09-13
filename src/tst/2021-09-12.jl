using JLD2, DataFrames, Serialization, LinearAlgebra, Distributions

import GEAS:breeding, sim_QTL, summarize, create_storage

function gpar(nqtl)
    Dict(
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
        :edit => false,
        :fix  => false
    )
end

function t_2021_09_12(nrpt = 10)
    BLAS.set_num_threads(12)
    nqtl = [100, 500]
    isdir("dat/tmp") || mkpath("dat/tmp")
    @load "dat/run/base.jld2" base
    qdist = [Normal(), Laplace()]
    rst = DataFrame()
    ng8n = 11

    for nq in nqtl
        @info "With n_QTL = $nq"
        dpr = gpar(nq)  # parameter dictionary
        for qd in qdist
            @info "QTL distribution $(string(qd))"
            for rpt in 1:nrpt
                @info "Repeat $rpt of $nrpt"
                @info "Preparing common base + generation 1"
                dpr[:fix], dpr[:edit] = false, false
                par = (; dpr...)
                qtl = sim_QTL(base, par.nQTL..., d=qd)
                snp, ped = create_storage(base, par, qtl)
                serialize("dat/tmp/snp.ser", snp)
                serialize("dat/tmp/ped.ser", ped)
                
                @info "Selection with SNP-BLUP"
                par = (; dpr...)
                breeding(ped, snp, base, qtl, par)
                df = summarize(ped.prd, snp.prd, qtl[2])
                df.repeat = repeat([rpt], ng8n)
                df.dist = repeat([string(qd)[1:3]], ng8n)
                df.nq = repeat([nq], ng8n)
                df.method = repeat(["SBLP"], ng8n)
                append!(rst, df)
                
                @info "Selection with emphasis on top $(par.nK3n) QTL"
                dpr[:fix], dpr[:edit] = true, false
                par = (; dpr...)
                snp = deserialize("dat/tmp/snp.ser")
                ped = deserialize("dat/tmp/ped.ser")
                breeding(ped, snp, base, qtl, par)
                df = summarize(ped.prd, snp.prd, qtl[2])
                df.repeat = repeat([rpt], ng8n)
                df.dist = repeat([string(qd)[1:3]], ng8n)
                df.nq = repeat([nq], ng8n)
                df.method = repeat(["Fix"], ng8n)
                append!(rst, df)
                
                dpr[:fix], dpr[:edit] = false, true
                par = (; dpr...)
                snp = deserialize("dat/tmp/snp.ser")
                ped = deserialize("dat/tmp/ped.ser")
                breeding(ped, snp, base, qtl, par)
                df = summarize(ped.prd, snp.prd, qtl[2])
                df.repeat = repeat([rpt], ng8n)
                df.dist = repeat([string(qd)[1:3]], ng8n)
                df.nq = repeat([nq], ng8n)
                df.method = repeat(["Edit"], ng8n)
                append!(rst, df)
                serialize("dat/rst/3w-cmp.ser", rst)
            end
        end
    end
end

function shortcut_2021_09_12()
    BLAS.set_num_threads(12)
    dpr = gpar(100)
    par = (; dpr...)
    @load "dat/run/base.jld2" base
    qtl = sim_QTL(base, par.nQTL...)
    snp, ped = create_storage(base, par, qtl)
    breeding(ped, snp, base, qtl, par)
    serialize("dat/tmp/sample.ser", Dict(:snp=>snp, :ped=>ped, :qtl=>qtl))
end
