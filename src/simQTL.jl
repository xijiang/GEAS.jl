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
    @info join(["",
                "Sample QTL locations and effects",
                "QTL effects were simulated such that, var(BV) ≈ 1"],
                "\n")
    pos = base[:pos]
    GT  = base[:hap]            # short hands
    tsnp = size(pos)[1]
    ntrt = length(nqtl)
    nid  = size(GT)[2] ÷ 2

    qinfo = QTL[]
    for n in nqtl
        loci = sort(randperm(tsnp)[1:n])
        
        q = zeros(Int8, n, nid) # QTL genotypes
        for id in 1:nid
            q[:, id] = base[:hap][loci, id*2 - 1] + base[:hap][loci, id*2]
        end
        e = randn(n)
        v = var(q'e)
        e ./= sqrt(v)
        m = mean(q'e)
        t = 2sum(e[e.>0]) - m
        push!(qinfo, QTL(loci, e, m, t))
    end
    qinfo
end
