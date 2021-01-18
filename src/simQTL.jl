"""
    function sim_QTL(base, nqtl...)
---
Given the population, a dictionary of SNP loci and SNP genotypes,
this function tries to:
1. sample QTL from the loci
2. sampe QTL effects such that μ ≈ 0.
3. the potential, or the upper limit of the population, or the ideal population is also calculated.
"""
function sim_QTL(base, nqtl...)
    @info "Sample QTL locations and effects"
    tsnp = length(base[:pos])
    ntrt = length(nqtl)

    loci = Any[]
    for n in nqtl               # sample locations
        push!(loci, sort(randperm(tsnp)[1:n]))
        push!(loci, sort(randperm(tsnp)[1:n]))
    end

    effect = Any[]
    for n in nqtl
        eb = randn(n)
    end
    Dict(:loci => loci)
end
