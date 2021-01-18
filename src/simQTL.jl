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

    loci  = Any[]
    efft  = Any[]
    ave   = Float64[]
    ideal = Float64[]
    for n in nqtl
        lc = sort(randperm(tsnp)[1:n])
        push!(loci, lc)         # sampled locations
        
        q = zeros(Int8, n, nid) # QTL genotypes
        for id in 1:nid
            q[:, id] = base[:hap][lc, id*2 - 1] + base[:hap][lc, id*2]
        end
        e = randn(n)
        v = var(q'e)
        e ./= sqrt(v)
        push!(efft, e)
        m = mean(q'e)
        push!(ave, m)
        t = sum(e[e.>0]) - m
        push!(ideal, t)
    end
    Dict(:loci => loci,
         :effects => efft,
         :mean => ave,
         :ideal => ideal)
end
