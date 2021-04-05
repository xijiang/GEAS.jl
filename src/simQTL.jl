"""
    function sim_QTL(base, nqtl...; d = Laplace())
---
Given the population, a dictionary of SNP loci and SNP genotypes,
this function tries to:
1. sample QTL from the loci
2. sampe QTL effects
   - by default from `Laplace()`
   - an alternative is `Γ(0.4, 1)`,
   - or from `N(0, 1)`, if `shape < 0`
   - effects are scaled such that σₐ² of base = 1.
3. the potential, or the upper limit of the population, or the ideal population is also calculated.

Note:
- A problem with `Γ(0.4, 1)` is that two traits will have an expectation correlation ~0.2.
"""
function sim_QTL(base, nqtl...; d = Laplace())
    @info join(["",
                "Sample QTL locations and effects",
                "QTL effects were simulated such that, var(BV) ≈ 1"],
                "\n")
    pos = base[:pos]
    GT  = base[:hap]            # short hands
    tsnp = size(pos)[1]
    ntrt = length(nqtl)

    qinfo = QTL[]
    for n in nqtl
        loci = sort(randperm(tsnp)[1:n])
        Q  = qtl_gt(GT, loci)
        #a = (shape < 0) ? randn(n) : rand(Gamma(0.4, 1), n)
        a = rand(d, n)
        vₐ = var(Q'a)           # variance of TBV
        a ./= sqrt(vₐ)          # scale allele effect, such that vₐ = 1
        m = mean(Q'a)           # base population expectation
        t = 2sum(a[a.>0]) - m   # expectation of an ideal ID
        push!(qinfo, QTL(loci, a, m, t))
    end
    qinfo
end

"""
    function qtl_gt(snp, qtl)
---
Given population genotypes `snp`, this function returns the QTL genotypes `Q`,
according to the loci specified in `qtl`.
"""
function qtl_gt(snp, qtl)
    nid = size(snp)[2] ÷ 2
    nlc = length(qtl)
    Q = zeros(Int8, nlc, nid)   # QTL genotypes 0/1 => 0, 1, 2
    for id in 1:nid
        Q[:, id] = snp[qtl, id*2 - 1] + snp[qtl, id*2]
    end
    Q
end
