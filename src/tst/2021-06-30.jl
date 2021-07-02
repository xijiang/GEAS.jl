using JLD2, Distributions, Serialization
import GEAS:sim_pt_QTL, breeding_program, top_QTL

"""
    function the_simulation(jld, nqtl, rt, nrpt; rst = "")
---
## Description
- 3 comparisons
  - Normal genomic selection, no known QTL
  - Known 5 QTL as fixed effects
  - Edit the biggest of known 5 QTL in generation 2-6

## parameters
1. jld: the base
2. nqtl: e.g., 500
3. rt: correlation of QTL effects between 2 traits
4. nrpt: number of repeats
5. rst: result path, in `2021-06-30/`
"""
function the_simulation(jld, nqtl, rt, nrpt; rst = "")
    @load jld base
    if length(rst) > 1
        isdir(rst) || mkpath(rst)
    end

    w = ndigits(nrpt)
    for r in 1:nrpt
        d = MvNormal(zeros(2), [1 rt; rt 1])
        qtl = sim_pt_QTL(base, nqtl, d)
        
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
                :nk3n => 5,     # number of known QTL on binary trait
                :dd   => 0.01, # A small value to be added to diagonals
                :log  => joinpath(rst, "qtl.log")
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

function simulation_2021_06_30()
    the_simulation("dat/run/base.jld2", 500, 0, 10, rst="dat/2021-06-30")
end

function simulation_2021_07_01()
    the_simulation("dat/run/base.jld2", 500, -.2, 10, rst="dat/2021-07-01")
end

function simulation_2021_07_02()
    the_simulation("dat/run/base.jld2", 500, .2, 10, rst="dat/2021-07-02")
end
