"""
    function test_2021_03_26()
---
To test:
- if the right loci is found
- if the editing is right
"""
function test_2021_03_26()
    @load "dat/run/base.jld2" base
    qtl = sim_QTL(base, 10)[1]
    println(qtl.pos)
    snp = copy(base[:hap])
    ep = gedit(base[:hap], qtl, .8)
    println(sum(snp[ep, :]), ' ', sum(base[:hap][ep, :]))
end
