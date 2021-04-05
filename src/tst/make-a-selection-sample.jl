"""
    function make_a_result()
---
Since one run with 50k snps and 10 generations needs one hour with this package,
this function runs once and store the results to `{prd,snp,qtl}.ser`
to write codes for summary.
"""
function make_a_result()
    par = begin
        nsib = 50
        pp = Dict(:nSire=> 100,
                  :nDam => 200, # male:female = 1:2
                  :nSib => 2nsib,
                  :nC7e => nsib, # number of sib to be challenged
                  :nG8n => 10,   # number of generations
                  :nQTL => [1000, 1000],
                  :h²   => [.5, .5],
                  :p8e  => .5, # percentage of affected in the binary trait
                  :e19e => .8, # edit_successsful_rate
                  :w4t  => 2   # weight on the binary trait EBV
                  )
        (; pp...)               # named tuple.  contents as above
    end
    @load "dat/run/base.jld2" base
    qtl = sim_QTL(base, par.nQTL...)
    prd, snp = breeding_program(base, par, qtl; edit=false)
    serialize("dat/run/prd.ser", prd)
    serialize("dat/run/snp.ser", snp)
    serialize("dat/run/qtl.ser", qtl)
end
