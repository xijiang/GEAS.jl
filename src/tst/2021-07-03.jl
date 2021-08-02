using JLD2, Distributions, Serialization, Dates, Statistics, DataFrames
import GEAS:sim_QTL, breeding_program, top_QTL, sim_gamma_QTL, sim_pt_QTL

function current_qtl(base, nqtl, op)
    if op == 1                  # gamma
        return sim_gamma_QTL(base, nqtl, nqtl)
    elseif op == 2              # mv normal
        d = MvNormal(zeros(2), [1 0; 0 1])
        return sim_pt_QTL(base, nqtl, d)
    elseif op == 3              # Laplacian
        return sim_QTL(base, nqtl, nqtl)
    end
end
    
"""
    function the_simulation(jld, nqtl, nrpt; rst = "")
---
## Description
- 3 comparisons
  - Normal genomic selection, no known QTL
  - Known 5 QTL as fixed effects
  - Edit the biggest of known 5 QTL in generation 2-6

## parameters
1. jld: the base
2. nqtl: e.g., 500
4. nrpt: number of repeats
5. rst: result path, in `2021-06-30/`
"""
function the_simulation(jld, nqtl, nrpt, desc, op; rst = "")
    @load jld base
    if length(rst) > 1
        isdir(rst) || mkpath(rst)
    end
    slog = joinpath(rst, "simu.log") # simulation log
    write(slog, desc)


    w = ndigits(nrpt)
    for r in 1:nrpt
        #qtl = sim_QTL(base, nqtl, nqtl)
        qtl = current_qtl(base, nqtl, op)

        suffix = lpad(string(r), w, '0')
        serialize(joinpath(rst, "qtl-$suffix.ser"), qtl)
        par = begin
            Parameters = Dict(
                :nSire=> 100,
                :nDam => 200,   # male:female = 1:2
                :nSib => 40,
                :nC7e => 30,    # number of sib to be challenged
                :nG8n => 11,    # number of generations
                :nQTL => [nqtl, nqtl],
                :h²   => [.5, .5],
                :p8e  => .5,    # death rate after challenge
                :e19e => 1.,    # edit_successsful_rate
                :w4t  => 1,     # weight on the binary trait EBV
                :nk3n => 10,    # number of known QTL on binary trait
                :dd   => 0.01, # A small value to be added to diagonals
                :log  => slog
            )
            (; Parameters...)       # named tuple.  contents as above
        end

        # genomic selection with no known QTL
        prd, snp = breeding_program(base, par, qtl)
        serialize(joinpath(rst, "prd-a-$suffix.ser"), prd)
        serialize(joinpath(rst, "snp-a-$suffix.ser"), snp)

        # genomic selection with known QTL as fixed effects
        prd, snp = breeding_program(base, par, qtl, fixed=true)
        serialize(joinpath(rst, "prd-b-$suffix.ser"), prd)
        serialize(joinpath(rst, "snp-b-$suffix.ser"), snp)

        # genomic selection, edit biggest known QTL each generation
        prd, snp = breeding_program(base, par, qtl, edit=true)
        serialize(joinpath(rst, "prd-c-$suffix.ser"), prd)
        serialize(joinpath(rst, "snp-c-$suffix.ser"), snp)
    end
end

function simulation_01()
    rst = joinpath("dat", string(DateTime(now())))
    @info "Results are in $rst"
    
    desc = "Uses a Γ(0.4) distritubion.
The Gamma distribution was modified to have half positive and half negative.
Total 100 QTL.
Number of known QTL 10.
Rank QTL fist of MAF(0.1, 0.9), then on abs(effects).

The top QTL of each round of simulation are as below:
"
    the_simulation("dat/run/base.jld2", 100, 10, desc, 1, rst=rst)
end

function simulation_02()
    rst = joinpath("dat", string(DateTime(now())))
    @info "Results are in $rst"
    
    desc = "Uses a pleiotropic MvNormal distritubion.
Total 100 QTL.
Number of known QTL 10.
Rank QTL fist of MAF(0.1, 0.9), then on abs(effects).

The top QTL of each round of simulation are as below:
"
    the_simulation("dat/run/base.jld2", 100, 10, desc, 2, rst=rst)
end

function simulation_03()
    rst = joinpath("dat", string(DateTime(now())))
    @info "Results are in $rst"
    
    desc = "Uses a independent Laplacian distribution. 
Total 100 QTL.
Number of known QTL 10.
Rank QTL fist of MAF(0.1, 0.9), then on abs(effects).

The top QTL of each round of simulation are as below:
"
    the_simulation("dat/run/base.jld2", 100, 10, desc, 3, rst=rst)
end

function simulation_04()
    rst = joinpath("dat", string(DateTime(now())))
    @info "Results are in $rst"
    
    desc = "Uses a independent Laplacian distribution. 
Total 500 QTL.
Number of known QTL 10.
Rank QTL fist of MAF(0.1, 0.9), then on abs(effects).

The top QTL of each round of simulation are as below:
"
    the_simulation("dat/run/base.jld2", 500, 10, desc, 3, rst=rst)
end

function qfrq(snp)
    fq = zeros(11, 10)
    start = 1
    i = 1
    for stop = 4000:4000:44000
        f = mean(snp[:, start:stop], dims=2)
        fq[i, :] = f'
        i += 1
        start = stop + 1
    end
    fq
end

"""
Generation mean and variance of `prd-abc` on `c`, e.g., `:tbv`.
"""
function mv(c, prd...)
    m = zeros(11, 3)
    v = zeros(11, 3)
    for i in 1:3
        gp = groupby(prd[i], :g8n)
        t = combine(gp, c => mean => :m, c => var => :v)
        m[:, i] = t.m
        v[:, i] = t.v
    end
    return m, v
end

function summarize_202107(dir, rpt)
    ir = lpad(string(rpt), 2, '0')
    @info "Read top QTL of the first repeat"
    log = readlines(joinpath(dir, "simu.log"))
    t = 1
    for i in 1:length(log)
        t += 1
        length(log[i]) > 0 && log[i][end] == ':' && break
    end
    t += 3(rpt - 1)
    tq = abs.(parse.(Int, split(log[t])))

    @info "Read SNP of the first repeat, ordinary GS"
    snp = deserialize(joinpath(dir, "snp-a-$ir.ser"))
    fq = qfrq(snp[tq, :])
    plot(fq, leg=false, title = "ogs")
    savefig(joinpath(dir, "afq.pdf"))

    @info "Read SNP of the first repeat, top QTL as fixed effects"
    snp = deserialize(joinpath(dir, "snp-b-$ir.ser"))
    fq = qfrq(snp[tq, :])
    plot(fq, leg=false, title = "fxq")
    savefig(joinpath(dir, "bfq.pdf"))
    
    @info "Read SNP of the first repeat, top QTL edited in turn"
    snp = deserialize(joinpath(dir, "snp-c-$ir.ser"))
    fq = qfrq(snp[tq, :])
    plot(fq, leg=false, title = "edt")
    savefig(joinpath(dir, "cfq.pdf"))

    @info "Read the production trait"
    lb = ["ogs" "fxq" "edt"]
    pa = deserialize(joinpath(dir, "prd-a-$ir.ser"))
    pb = deserialize(joinpath(dir, "prd-b-$ir.ser"))
    pc = deserialize(joinpath(dir, "prd-c-$ir.ser"))
    m, v = mv(:t2c, pa, pb, pc)
    plot(m, label = lb, xlabel="generation", ylabel="Bin mean")
    savefig(joinpath(dir, "bin-mean.pdf"))
    plot(v, label = lb, xlabel="generation", ylabel="Bin var")
    savefig(joinpath(dir, "bin-var.pdf"))
    m, v = mv(:tbv, pa, pb, pc)
    plot(m, label = lb, xlabel="generation", ylabel="Prd mean")
    savefig(joinpath(dir, "prd-mean.pdf"))
    plot(v, label = lb, xlabel="generation", ylabel="Prd var")
    savefig(joinpath(dir, "prd-var.pdf"))
end

function sum_mmv(dir, nr)
    mm = zeros(11, 3)
    vm = zeros(11, 3)
    for i in 1:nr
        r = lpad(string(i), 2, '0')
        pa = deserialize(joinpath(dir, "prd-a-$r.ser"))
        pb = deserialize(joinpath(dir, "prd-b-$r.ser"))
        pc = deserialize(joinpath(dir, "prd-c-$r.ser"))
        m, v = mv(:t2c, pa, pb, pc)
        mm += m
        vm += v
    end
    mm ./= nr
    vm ./= nr
    lb = ["ogs" "fxq" "edt"]
    plot(mm, label = lb, xlabel="generation", ylabel="Average bin mean")
    savefig(joinpath(dir, "ave-bin-mean.pdf"))
    plot(vm, label = lb, xlabel="generation", ylabel="Average bin var")
    savefig(joinpath(dir, "ave-bin-var.pdf"))
end

#=
!!!!! Obsolated !!!!!

"""
This is to use Laplacian distribution to have a few large QTL.
Only 100 QTL were simulated to check the editing effects.
Before, 500 QTL scenario show little advantage of the editing.
And, emphasize known QTL effects give negative selection effects.
"""
function simulation_2021_07_03()
    the_simulation("dat/run/base.jld2", 100, 10, "test", rst="dat/2021-07-03")
end

function simulation_2021_07_04()
    rst = joinpath("dat", string(DateTime(now())))
    @info "Results are in $rst"
    
    desc = "The previous simulation just recorded the top QTL positions in the last repeat.
This simulation changed `w` to `a` for QTL log.
This also uses a Gamma distritubion.
The Gamma distribution was modified to have half positive and half negative.
Number of known QTL changed from 5 to 10.

The top QTL of each round of simulation are as below:
"
    the_simulation("dat/run/base.jld2", 100, 10, desc, rst=rst)
end

function simulation_2021_07_05()
    rst = joinpath("dat", string(DateTime(now())))
    @info "Results are in $rst"
    
    desc = "The previous simulation just recorded the top QTL positions in the last repeat.
This simulation changed `desc` writing position out of the loop.
This also uses a Gamma distritubion.
The Gamma distribution was modified to have half positive and half negative.
Number of known QTL changed from 5 to 10.
QTL ranking was changed from var to effect only.
Only MAF of [0.1, 0.9] are to be considered.

The top QTL of each round of simulation are as below:
"
    the_simulation("dat/run/base.jld2", 100, 10, desc, rst=rst)
end
=#
