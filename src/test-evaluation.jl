"""
    function test_eval()
---
Create F₁ to test genomic seleciton methods.
"""
function test_eval(; just_one = false)
    #base = deserialize("dat/run/ns.ser")
    @load "dat/run/base.jld2"
    
    nSire, nDam, nSib = 200, 400, 42
    h² = [.5, .5]
    par = begin
        nqtl = repeat([100, 500, 1000], inner=5)
        nclg = repeat([37, 32, 22, 12, 2], 3)
        [nqtl nclg]
    end
    ped = random_mate(nSire, nDam, nSib)
    nxt = gdrop(base[:hap], ped, base[:r])

    rst = open(joinpath(dat_dir, "run/test-evaluation.txt"), "w")
    lhs = nothing
    for (nqtl, nclg) in eachrow(par)
        qtl = sim_QTL(base, [nqtl, nqtl]...)
        ip, ic, hp, hc = test_sets(nDam, nSib, nclg)
        bv₁, p₁ = phenotype(nxt, qtl[1], h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], h²[2], 0.5) # 50% death

        # genotypes and phenotypes
        gp = alleles2gt(nxt[:, hp]) # genotypes
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
        nprd = nSib - nclg
        print("$nqtl, $nprd: ", round.(results; digits = 3), "\n")
        print(rst, "$nqtl, $nprd: ", round.(results; digits = 3), "\n")

        just_one && break
    end
    close(rst)
end
