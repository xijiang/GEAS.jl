using JLD2, Distributions, Serialization, Dates
import GEAS:sim_QTL, breeding_program, top_QTL, sim_gamma_QTL

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
        qtl = sim_gamma_QTL(base, nqtl, op)
        
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

function summarize_05()
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
