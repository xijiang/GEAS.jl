import GEAS:qmsim_base, sim_pt_QTL, breeding_program, bp_summary
using Distributions

"""
    function my_sim_pt_qtl()
---
1. Simulate a base population with `QMSim`.
2. go through a breeding program
3. summarize
"""
function my_sim_pt_qtl()
    dicpar = Dict(
        :hg_size => "1200 [0]\n1200 [2000]",
        :male => 600,
        :female => 600,         # 600 + 600 => 1200
        :nchr => 29,
        :lchr => 100,
        :nsnp => 1800,          # to simulate  1725*29 => 50025 SNP
    )
    bpar = (; dicpar...)
    base = qmsim_base(bpar, dir="dat/test/QMSim_Linux")
    #d = MvNormal(zeros(2), [1 -.1; -.1 1])
    #qtl = sim_pt_QTL(base, 1000, d)
    #
    #nqtl = [length(qtl[1].pos), length(qtl[2].pos)]
    #par = begin
    #    Parameters = Dict(:nSire=> 100,
    #                      :nDam => 200, # male:female = 1:2
    #                      :nSib => 60,
    #                      :nC7e => 30, # number of sib to be challenged
    #                      :nG8n => 10, # number of generations
    #                      :nQTL => nqtl, # usually [1000, 1000]
    #                      :h²   => [.5, .5],
    #                      :p8e  => .5, # percentage of affected in the binary trait
    #                      :e19e => .8, # edit_successsful_rate
    #                      :w4t  => 1 # weight on the binary trait EBV
    #                      )
    #    (; Parameters...)       # named tuple.  contents as above
    #end
    #prd, snp = breeding_program(base, par, qtl)
    #bp_summary(prd, snp, qtl, dir="dat/test/QMSim_Linux")
end
