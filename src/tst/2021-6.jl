using JLD2, Distributions, Serialization
import GEAS:sim_pt_QTL, breeding_program, top_QTL, QTL

"""
    function run_2021_06
---
This test is the first formal test on editint.
I will use the real genotype first.
This is because that the genotypes from `QMSimu` and `ms` still look wierd to me.

## Scenario to be considered
- Number of QTL, 100:450:1000
- Number of full sibs per family, 60:20:100.  Half of them are to be challenged
- number of QTL to be edited, 3, or 5 biggest in the base
  - I send a subset of the simulated QTL for editing
  - The subset is determined in the breeding program, not here.
- Correlation `r` between the QTL of the 2 traits.
- Number of sibs per full sib family
  - Always challenge half of them
- Editing vs. non editing
- Number of repeats will be as number of instances to be submitted on `Saga`.
"""
function run_2021_06()
    @load "dat/run/base.jld2" base
    rst = "dat/editing"
    isdir(rst) || mkpath(rst)
    
    for nqtl in 100:450:1000
        for r in -.3:.15:.3
            d = MvNormal(zeros(2), [1 r; r 1])
            qtl = sim_pt_QTL(base, nqtl, d)
            RR = Int(r*100)
            serialize("$rst/qtl-$nqtl-$RR.ser", qtl)
            for nsib in 30:10:50
                par = begin
                    Parameters = Dict(
                        :nSire=> 100,
                        :nDam => 200, # male:female = 1:2
                        :nSib => 2nsib,
                        :nC7e => nsib, # number of sib to be challenged
                        :nG8n => 10,   # number of generations
                        :nQTL => [nqtl, nqtl],
                        :h²   => [.5, .5],
                        :p8e  => .5, # percentage of affected in the binary trait
                        :e19e => .8, # edit_successsful_rate
                        :w4t  => 1, # weight on the binary trait EBV
                        :nk3n => 3  # number of known QTL
                    )
                    (; Parameters...) # named tuple.  contents as above
                end
                prd1, snp1 = breeding_program(base, par, qtl)
                serialize("$rst/prd-$nqtl-$RR-$nsib-n.ser", prd1)
                serialize("$rst/snp-$nqtl-$RR-$nsib-n.ser", snp1)
                prd2, snp2 = breeding_program(base, par, qtl, edit=true)
                serialize("$rst/prd-$nqtl-$RR-$nsib-y.ser", prd2)
                serialize("$rst/snp-$nqtl-$RR-$nsib-y.ser", snp2)
            end
        end
    end
end


function my_test()
    @load "dat/run/base.jld2" base
    r = .5
    d = MvNormal(zeros(2), [1 r; r 1])
    qtl = sim_pt_QTL(base, 100, d)
    q = QTL(qtl[2].pos[1:3], qtl[2].effect[1:3], 0, 0)
end

