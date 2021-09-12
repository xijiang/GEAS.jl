using JLD2, DataFrames, Serialization, LinearAlgebra, Distributions, LaTeXStrings
#, Distributions, Serialization, Dates, Statistics, DataFrames
import GEAS:breeding, sim_QTL, summarize, create_storage, simple_breeding,
    sim_gamma_QTL, rank_QTL, sum_qtl_frq

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
    BLAS.set_num_threads(12)
    nrpt = 30
    
    # to save the common base + generation 1
    isdir("dat/tmp") || mkpath("dat/tmp")
    @load "dat/run/base.jld2" base
    nqtl = 500
    par = gpar(nqtl, false, false)
    mtd = ["TBV", "phenotype", "EBV"]
    d = Normal()
    
    rst = DataFrame()
    for ir in 1:nrpt
        qtl = sim_QTL(base, par.nQTL..., d=d)
        rkq = rank_QTL(base.hap, qtl[1]) # QTL ranking in ↓ order
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
            df = summarize(ped.prd, snp.prd, qtl[2])
            df.method = ones(Int, nrow(df)) .* op
            df.repeat = ones(Int, nrow(df)) .* ir
            append!(rst, df)
        end
    end
    serialize("dat/rst/one-trait.ser", rst)
    rst
end

"""
    function freq_change()
---
This to test frequency changes of the top QTL with different selection methods
and different QTL distribution.

Before with Laplacian distribution, GS was already very effective at changing
the QTL frequencies, which is suspicious.
"""
function freq_change()
    BLAS.set_num_threads(12)
    dq = [Laplace(), Normal(), Gamma(0.4)]
    isdir("dat/tmp") || mkpath("dat/tmp")
    mtd = ["TBV", "phenotype", "EBV"]
    @load "dat/run/base.jld2" base
    
    for nq in [100, 500]
        for d in dq
            qtl = sim_QTL(base, nq, nq, d=d)
            tpq = rank_QTL(base.hap, qtl[1])
            par = gpar(nq, false, false)
            snp, ped = create_storage(base, par, qtl)
            serialize("dat/tmp/snp.ser", snp)
            serialize("dat/tmp/ped.ser", ped)
            for op in 1:3
                println("----- Selection on $(mtd[op]) -----")
                snp = deserialize("dat/tmp/snp.ser")
                ped = deserialize("dat/tmp/ped.ser")
                simple_breeding(ped, snp, base, qtl, par, op)
                sp = snp.prd[abs.(tpq), :]
                sc = snp.prd[qtl[2].pos, :]
                prd = ped.prd
                @save "dat/tmp/q-$nq-$(string(d)[1:3])-$(mtd[op][1:3])" prd sp sc
            end
        end
    end
end

function org_one_trait_rst()
    df = deserialize("dat/rst/one-trait.ser")
    tbv, phe, ebv = groupby(df, :method)
    a = combine(groupby(tbv, :g8n),
            :mpbv => mean => :mpbv,
            :pvg => mean => :mpvg,
            :mcbv => mean => :mcbv,
            :cvg => mean => :mcvg,
                   :mica => mean => :mica)
    b = combine(groupby(phe, :g8n),
            :mpbv => mean => :mpbv,
            :pvg => mean => :mpvg,
            :mcbv => mean => :mcbv,
            :cvg => mean => :mcvg,
            :mica => mean => :mica)
    c = combine(groupby(ebv, :g8n),
            :mpbv => mean => :mpbv,
            :pvg => mean => :mpvg,
            :mcbv => mean => :mcbv,
            :cvg => mean => :mcvg,
                :mica => mean => :mica)
    
    p1 = plot(1:11, a.mpbv, label = "on TBV",
              ylabel = L"\textrm{Production\ }\overline{\mathrm{TBV}}",
              dpi = 300,
              legendfontsize=6,
              legend=:bottomright);
    plot!(p1, 1:11, b.mpbv, label = "on Phenotype");
    plot!(p1, 1:11, c.mpbv, label = "on EBV");
    
    p2 = plot(1:11, a.mpvg, label = "on TBV",
              ylabel = L"\textrm{Production\ }\overline{V_g}",
              dpi = 300,
              legendfontsize=6,
              legend=:topright);
    plot!(p2, 1:11, b.mpvg, label = "on Phenotype");
    plot!(p2, 1:11, c.mpvg, label = "on EBV");

    p3 = plot(1:11, a.mcbv, label = "on TBV",
              ylabel = L"\textrm{Challenge\ }\overline{\mathrm{TBV}}",
              dpi = 300,
              legendfontsize=6,
              legend=:bottomright);
    plot!(p3, 1:11, b.mcbv, label = "on Phenotype");
    plot!(p3, 1:11, c.mcbv, label = "on EBV");

    p4 = plot(1:11, a.mcvg, label = "on TBV",
              xlabel = "Generation",
              ylabel = L"\textrm{Challenge\ }\overline{V_g}",
              dpi = 300,
              legendfontsize=6,
              legend=:bottomleft);
    plot!(p4, 1:11, b.mcvg, label = "on Phenotype");
    plot!(p4, 1:11, c.mcvg, label = "on EBV");

    p5 = plot(1:11, a.mica, label = "on TBV",
              ylabel = L"\overline{\mathrm{Inbreeding}}",
              dpi = 300,
              legendfontsize=6,
              legend=:topleft);
    plot!(p5, 1:11, b.mica, label = "on Phenotype");
    plot!(p5, 1:11, c.mica, label = "on EBV");
    p6 = plot(p1, p2, p5, p3, p4);
    savefig(p6, "dat/rst/one-trait.pdf")
end
