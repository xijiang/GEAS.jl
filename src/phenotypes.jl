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
    return bv + er
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
    ni = size(snp)[2] ÷ 2       # number of individuals
    er = randn(ni) .* σₑ
    return bv + er
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
function phenotype(snp, qtl, h², threshold)
    ph = phenotype(snp, qtl, h²)
    ni = length(ph)
    bn = zeros(Int8, ni)
    th = qtl.mean + threshold
    @inbounds for i in 1:ni
        bn[i] = (ph[i] < threshold) ? 1 : 0
    end
    return bn
end

#=
function TBV(snp, qtl)
    nid = size(snp)[2] ÷ 2
    pos = qtl.pos
    nql = length(pos)
    
    g = zeros(Int8, nql, nid)
    e = qtl.effect
    for id in 1:nid
        g[:, id] = snp[pos, id*2 - 1] + snp[pos, id*2]
    end

    bv = g'e
    tmp = DataFrame(bv = bv, id = 1:nid)
    sort!(tmp, :bv, rev=true)
    rank = tmp.id               # rank offspring in BV order, high->low
    @inbounds for s in 50:50:600
        x = g[:, rank[1:s]]     # QTL genotypes of the selected ones
        b = bv[rank[1:s]]       # BV of the selected ones
        q = sum(x, dims=2)      # frequency: 0 <--> 2s
        m = mean(b)
        v = var(b)
        t = qtl.max
        f = 0
        @inbounds for i in 1:length(e)
            if q[i] == 0
                if e[i] > 0
                    t -= e[i]
                end
                f += 1
            elseif q[i] == 2s
                if e[i] < 0
                    t += e[i]
                end
                f += 1
            end
        end
        println("$m, $v, $t, $f")
    end
end
=#
