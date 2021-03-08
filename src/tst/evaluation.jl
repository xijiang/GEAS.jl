"""
    function test_2_eval()
---
This is to test the `evaluation` function on two traits.
Only `base` genotypes were used.
The first 400 are for production.
The rest 200 are for challenge.
"""
function test_2_eval()
    @load "dat/run/base.jld2"
    qtl = sim_QTL(base, 500, 500)
    snp1 = base[:hap][:, 1:600] # used for production
    snp2 = base[:hap][:, 601:end] # used for challenge

    # bv and phenotypes for the production pop
    b1, p1 = phenotype(snp1, qtl[1], .5)

    # and also simulate bv and phenotypes as if challenged
    b3, p3 = phenotype(snp1, qtl[2], .5, .5)

    # bv and phenotypes for the challenge pop.
    b2, p2 = phenotype(base[:hap][:, 601:end], qtl[2], .5, .5)

    # SNP effects for these two pops.
    g1 = alleles2gt(snp1)       # genotypes of 0, 1, and 2
    g2 = alleles2gt(snp2)
    m1, s1 = snp_blup(g1, p1)
    m2, s2 = snp_blup(g2, p2)

    ebv1 = g1's1 .+ m1
    ebv2 = g2's2 .+ m2
    ebv3 = g1's2 .+ m2
    println(cor(ebv1, b1), ' ', cor(ebv1, p1))
    println(cor(ebv2, b2), ' ', cor(b2, p2))
    println(cor(ebv3, b1), ' ', cor(ebv3, p1))
end

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
        nqtl = repeat([100, 500, 1000], inner=3)
        nclg = repeat([32, 22, 2], 3)
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
