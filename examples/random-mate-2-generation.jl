# to use this, `using Revise; using GEAS; includet("this file")`
# then `ped, snp = two_random_generations(dir)`
# use `ped, snp` for your test.

using JLD2, Serialization, DataFrames
import GEAS:sim_QTL, create_storage, breeding, alleles2gt, snp_blup
"""
    function two_random_generations()
---
This is to generate a random mating two-generation pedigree to test about SNP-BLUP.
One-trait Results genotypes, pedigree and QTL info are written to `dir`.
"""
function two_random_generations()
    par = begin
        nsib = 10
        pp = Dict(:nSire=> 100,
                  :nDam => 200, # male:female = 1:2
                  :nSib => 2nsib,
                  :nC7e => nsib, # number of sib to be challenged
                  :nG8n => 1,    # number of generations
                  :nQTL => [100, 100],
                  :h²   => (a = .5, b = .5),
                  :p8h  => .5,  # percentage-of-death after challenge
                  :e19e => 1.,  # edit_successsful_rate
                  :b4y  => false,
                  :fix  => false,
                  :w4t  => 2    # weight on the binary trait EBV
                  )
        (; pp...)               # named tuple.  contents as above
    end
    @load "dat/run/base.jld2" base
    qtl = sim_QTL(base, par.nQTL...)
    ped, snp = create_storage(base, par, qtl)
    base, ped, snp
end


function test_trg()
    base, ped, snp = two_random_generations()
    g = alleles2gt(snp.prd)
    p = ped.prd.p7e
    snp_blup(g, p, .5, base.twop)
end
