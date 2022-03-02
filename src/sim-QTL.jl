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
- Randomly flip the sign of the effects may be a solution.
"""
function sim_QTL(base, nqtl...; d = Laplace(), ϵ = 1e-5)
    @debug "Sample QTL locations and effects" "With the QTL effects, mean ≈ 0, var(BV) ≈ 1"
    pos = base.pos
    GT  = base.hap              # short hands
    tsnp = size(pos)[1]
    ntrt = length(nqtl)

    qinfo = QTL[]
    for n in nqtl
        loci = sort(randperm(tsnp)[1:n])
        Q  = qtl_gt(GT, loci)
        a = rand(d, n) .* rand([-1, 1], n) # random effect directions
        y = Q'a
        m, s = mean(y), std(y)
        while abs(m) > ϵ || abs(s-1) > ϵ
            a .-= 2m/n
            a ./= s
            y = Q'a
            m, s = mean(y), std(y)
        end
        push!(qinfo, QTL(loci, a))
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

"""
    function sim_pt_QTL(base, nqtl, d)
---
Simulation of pleiotropic QTL, where `pt` is a shorthand of pleiotropic.
That is, sample QTL for two traits that are all of the same positions.
The effects of this set of QTL for the two traits,
however, are correlated.
For example, they can be reasonably negative correlated.

An example of `d` is `MvNormal(zeros(2), [1 0; 0 1.])`.
"""
function sim_pt_QTL(base, nqtl, d)
    @debug "Sample pleiotropic QTL for two trait" "QTL effects are define by `d`"
    pos = base.pos
    GT  = base.hap              # short hands
    tsnp = size(pos)[1]
    qinfo = QTL[]
    # two traits share below QTL
    lqtl = sort(randperm(tsnp)[1:nqtl])
    eqtl = rand(d, nqtl)'       # effects for each trait in column
    Q = qtl_gt(GT, lqtl)        # QTL genotypes
    for i in 1:2
        a = eqtl[:, i]
        y = Q'a
        m = mean(y)             # base population mean
        v = var(y)
        a ./= sqrt(v)
        push!(qinfo, QTL(lqtl, a))
    end
    qinfo
end
