"""
    function alleles2gt(snp)
---
Merge alleles to genotypes.

- SNP haplotypes are SNP × (Nid × 2)
- SNP alleles are 0 & 1.
- 00 → 0; 01/10 → 1; 11 → 2
"""
function alleles2gt(snp)
    nlc, nid = size(snp)
    nid ÷= 2
    g = Array{Int8, 2}(undef, nlc, nid)
    @inbounds for i in 1:nid
        g[:, i] = snp[:, 2i-1] + snp[:, 2i]
    end
    return g
end

"""
    function hap2gt(hap; GBLUP = false)
Note `hap` here are ID majored.  The two haplotypes of one ID are in two columns.
This is for ease of gene-dropping.
"""
function hap2gt(hap; GBLUP = false)
    nlc, nid = size(hap)
    nid ÷= 2
    if GBLUP
        gt = Matrix{Int8}(undef, nlc, nid)
        for i in 1:nid
            gt[:, i] = hap[:, 2i-1] + hap[:, 2i]
        end
    else
        gt = Matrix{Int8}(undef, nid, nlc)
        for i in 1:nid
            gt[i, :] = hap[:, 2i-1] + hap[:, 2i]
        end
    end
    gt
end

"""
    function breeding_value(snp, qtl)
Given SNP genotypes `snp`, and QTL loci and effects in `qtl`, this function return
true breeding values of each individual.
"""
function breeding_value(snp, qtl)
    l = qtl.pos
    e = qtl.effect             # shorthands
    n = size(snp)[2] ÷ 2       # number of individuals
    g = alleles2gt(snp[l, :])
    return g'e
end

"""
    function phenotype(bv, h²)
---
Given breeding values `bv` and `h²`, and assuming `σₐ²` = 1,
this function returns a phenotype array.
"""
function phenotype(bv, h²)
    σₑ = √(1/h² - 1)            # σₐ² = 1
    ni = length(bv)             # number of individuals
    er = randn(ni) .* σₑ
    return bv, bv + er
end

"""
    function phenotype(snp, qtl, h²)
---
This function returns phenotypes of a trait.
`σₑ²` is determined assuming `σₐ²` = 1 in the base population.
"""
function phenotype(snp, qtl, h²)
    σₑ = √(1/h² - 1)            # σₐ² = 1
    bv = breeding_value(snp, qtl)
    ni = size(bv)
    er = randn(ni) .* σₑ
    return bv, bv + er
end

"""
    function phenotype(snp, qtl, h², threshold)
---
This function returns phenotypes of a binary trait.
`σₑ²` is determined assuming `σₐ² = 1` in the base population.
This result can also serve as a challenge result.
Assuming the liability is for resisitance.
ID with a `phenotype < qtl.mean + threshold` are incident/dead.
That is an ID is alive and phenotype is 1.
"""
function phenotype(snp, qtl, h², percentage)
    bv, ph = phenotype(snp, qtl, h²)
    th = quantile(ph, percentage)
    ni = length(ph)
    bn = zeros(Int8, ni)
    @inbounds for i in 1:ni
        (ph[i] > th) && (bn[i] = 1)
    end
    return bv, bn
end
