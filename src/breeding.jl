"""
    function test_sets(nfam, nsib, nclg)
---
Given number of full sib families `nfam`, full sibship size `nsib`,
and number of ID to be challenged `nclg`,
this function return two vectors.
They indicate ID index for challenge and for production.
There are no need to random sample again as they were randomly dropped
from the previous generation.
"""
function test_sets(nfam, nsib, nclg)
    nprd = nsib - nclg
    clg = repeat(1:nclg, nfam)
    prd = repeat(nclg+1:nsib, nfam)
    a, b, c, d, inc = 1, nclg, 1, nprd, 0
    for _ in 1:nfam
        clg[a:b] .+= inc
        prd[c:d] .+= inc
        a += nclg
        b += nclg
        c += nprd
        d += nprd
        inc += nsib
    end
    prd, clg
end

"""
    function evaluate_n_select(obs, nSire, nDam)
---
Given observation `obs`, return genotypes of `nSire + nDam`.
"""
function evaluate_n_select(obs, nCur, nSire, nDam)
    # Remember to write/return the numbers of ID selected 
    @warn "evaluation: under construction"
    nh = (nSire + nDam) * 2
    obs.g₁[:, 1:nh]            # only select from the production group
end

"""
    function generation_one(base, par, qtl)
---
The base population size varies.
The rest of the generations have a constant size.
This function is to scale generation one, such that,
the production population and the challenge population are
approximately equal to the later generations.
This returns SNP genotypes of the nuclear population of generation 2.
"""
function generation_one(base, par, qtl)
    R   = par.nDam ÷ par.nSire     # female:male ratio
    if R < 1
        @warn join(["",
                    "- Sires are more than dams.",
                    "- Set female:male to 1:1 in base."],
                   "\n")
        R = 1
    end
    n₀  = size(base[:hap])[2] ÷ 2  # base population size
    nd  = n₀ ÷ (R + 1) * R         # number of dams
    ns  = par.nDam * par.nSib ÷ nd # sibship size
    nc  = par.nDam * par.nC7e ÷ nd # number challenged

    # breed and measure
    obs, nCur = begin
        x₁, x₂ = test_sets(nd, ns, nc)
        ped = random_mate(n₀ - nd, nd, ns)
        nxt = gdrop(base[:hap], ped, base[:r])
        p₁  = phenotype(nxt, qtl[1], par.h²[1])
        p₂  = phenotype(nxt, qtl[2], par.h²[2], par.t7d)
        tmp = Dict(:g₁ => nxt[:, x₁],
                   :g₂ => nxt[:, x₂],
                   :p₁ => p₁[x₁],
                   :p₂ => p₂[x₂],
                   :p₃ => ped[x₁, :])
        tmp, length(x₁)
    end
    serialize(joinpath(dat_dir, "run/pre.ser"), obs)
    
    # select nuclear of next generation
    evaluate_n_select((; obs...), nCur, par.nSire, par.nDam)
end

"""
    function breeding_program(base, par)
---
Do the breeding with the given parameter `par`.
"""
function breeding_program(base, par)
    # Parameters for this simulation
    qtl = sim_QTL(base, par.nQTL...)
    println()                   # breeding cycles begin below

    start = 1
    snp = begin                 # genotypes for nuclear population
        if size(base[:hap])[2] ÷ 2 ≠ par.nSire + par.nDam
            start = 2
            @info "F-1: Scaled from base"
            generation_one(base, par, qtl)
        else
            base[:hap]
        end
    end

    x₁, x₂ = test_sets(par.nDam, par.nSib, par.nC7e)
    for ig in start:par.nG8n    # breeding and seletion cycles
        @info "F-$ig"
        cur = begin
            ped = random_mate(par.nSire, par.nDam, par.nSib)
            nxt = gdrop(snp, ped, base[:r])
            p₁ = phenotype(nxt, qtl[1], par.h²[1])
            p₂ = phenotype(nxt, qtl[2], par.h²[2], par.t7d)
            Dict(:g₁ => nxt[:, x₁],
                 :g₂ => nxt[:, x₂],
                 :p₁ => p₁[x₁],
                 :p₂ => p₂[x₂],
                 :p₃ => ped[x₁, :])
        end
        obs = begin
            if ig > 1
                tmp = deserialize(joinpath(dat_dir, "run/pre.ser"))
                Dict(:g₁ => [cur[:g₁]  tmp[:g₁]],
                     :g₂ => [cur[:g₂]  tmp[:g₂]],
                     :p₁ => [cur[:p₁]; tmp[:p₁]],
                     :p₂ => [cur[:p₂]; tmp[:p₂]],
                     :p₃ => [cur[:p₃]; tmp[:p₃]])
            else
                cur
            end
        end
        serialize(joinpath(dat_dir, "run/pre.ser"), obs)
        snp = evaluate_n_select((; obs...), length(x₁), par.nSire, par.nDam)
    end
end
