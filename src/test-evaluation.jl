"""
    function test_evaluation(base)
---
Create F₁ to test genomic seleciton methods.

## Usage
```julia
using Serialization
base = deserialize("dat/run/ns.ser")
GEAS.test_evaluation(base)
```
"""
function test_evaluation()
    base = deserialize("dat/run/ns.ser")
    
    nSire, nDam, nClg = 200, 400, 10
    h² = [.5, .5]
    par = begin
        nqtl = repeat([100, 500, 1000], inner=3)
        nsib = repeat([5, 10, 30], 3) .+ nClg
        [nqtl nsib]
    end
    # below can actually just simulate a large sibship size,
    # then sample QTL and sub-sibships
    # but since the speed is OK.  I ran 9x time, to avoid re-programming
    for (nqtl, nSib) in eachrow(par)
        # simulation
        qtl = sim_QTL(base, [nqtl, nqtl]...)
        ped = random_mate(nSire, nDam, nSib)
        nxt = gdrop(base[:hap], ped ,base[:r])
        ip, ic, hp, hc = test_sets(nDam, nSib, nClg)
        bv₁, p₁ = phenotype(nxt, qtl[1], h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], h²[2], 0.5) # 50% death

        # genotypes and phenotypes
        gp = alleles2gt(nxt[:, hp])                  # genotypes
        gc = alleles2gt(nxt[:, hc])
        pp, pc = p₁[ip], p₂[ic] # phenotypes

        results = zeros(6)
        results[1:2] = [cor(bv₁, p₁), cor(bv₂, p₂)]
        # evaluation on phenotypes
        _, sp = snp_blup(gp, pp)
        gebv = gp'sp
        bv = bv₁[ip]
        results[3:4] = [cor(gebv, bv), cor(gebv, pp)]

        # evaluation on TBV
        _, sp = snp_blup(gp, bv)
        gebv = gp'sp
        results[5:6] = [cor(gebv, bv), cor(gebv, pp)]
        print(round.(results; digits = 3), "\n")
    end
end
