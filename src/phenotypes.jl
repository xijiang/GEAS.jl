"""
Given SNP genotypes `snp`, and QTL loci and effects in `qtl`, this function return
true breeding values of each individual.
"""
function breeding_value(snp, qtl)
    l = qtl.pos
    e = qtl.effect             # shorthands
    n = size(snp)[2] ÷ 2       # number of individuals
    g = zeros(Int8, length(l), n)      # QTL genotypes (012)
    @inbounds for i in 1:n
        g[:, i] = snp[l, 2i-1] + snp[l, 2i]
    end
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
σₑ² is determined assuming σₐ²=1 in the base population.
This result can also serve as a challenge result.
Assuming the liability is for resisitance.
ID with a `phenotype < qtl.mean + threshold` are incident.
That is an ID is affected and phenotype is 1.
"""
function phenotype(snp, qtl, h², percentage)
    bv, ph = phenotype(snp, qtl, h²)
    th = quantile(ph, percentage)
    ni = length(ph)
    bn = zeros(Int8, ni)
    @inbounds for i in 1:ni
        (ph[i] < th) && (bn[i] = 1)
    end
    return bv, bn
end
