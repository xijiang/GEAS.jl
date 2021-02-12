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
    hp = 2repeat(prd, inner=2)
    hp[1:2:end] .-= 1
    hc = 2repeat(clg, inner=2)
    hc[1:2:end] .-= 1
    prd, clg, hp, hc
end

"""
    function select_nuclear(df, nSire, nDam)
---
This is to select `nSire`, and `nDam` from `df` on df.val.
DataFrame `df` has 3 columns:
1. ID number
2. sex, 0 for female, 1 for male.
3. val

This is constructed from a evaluation procedure, 
where a method can be chosen to calculate `df.val`.
Be careful that no range check is performed.
"""
function select_nuclear(df, nSire, nDam)
    gp = groupby(df, :sex)
    sires = last(sort(gp[1], :val), nSire).ID
    dams  = last(sort(gp[2], :val), nDam ).ID
    return [sires; dams]
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
    R   = par.nDam ÷ par.nSire  # female:male ratio
    if R < 1
        @warn join(["",
                    "- Sires are more than dams.",
                    "- female:male set to 1:1 in base."],
                   "\n")
        R = 1
    end
    n₀ = size(base[:hap])[2] ÷ 2  # base population size
    nd = n₀ ÷ (R + 1) * R         # number of dams, or families
    ns = par.nDam * par.nSib ÷ nd # sibship size
    nc = par.nDam * par.nC7e ÷ nd # number challenged
    nn = par.nSire + par.nDam     # nuclear population size

    # breed and measure
    x₁, x₂, i₁, i₂ = test_sets(nd, ns, nc)
    np₁ = length(x₁)            # nID for production in generation 1
    ped = random_mate(n₀ - nd, nd, ns)
    nxt = gdrop(base[:hap], ped, base[:r])
    bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
    bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
    sex = rand(0:1, np₁)
    obs = Dict(:g₁  => nxt[:, i₁],
               :g₂  => nxt[:, i₂],
               :bv₁ => bv₁[x₁],
               :bv₂ => bv₂[x₂],
               :p₁  => p₁[x₁],
               :p₂  => p₂[x₂],
               :n₂  => [length(x₂)],
               :ped => [ones(Int16, length(x₁)) sex ped[x₁, :]])
    
    # select nuclear of next generation
    # !!!!!!!!!! evaluation line should be inserted here.
    df = DataFrame(ID = 1:np₁, sex = rand(0:1, np₁), val = obs[:p₁])
    nuclear = select_nuclear(df, par.nSire, par.nDam)
    obs, nuclear
end

"""
    function determine_storage(base, par)
---
Determine array sizes for one-go genotype storage.
This is to avoid too much memory allocation, 
and garbage collection.
Above also cause program crash.
"""
function determine_storage(base, par)
    n₀ = size(base[:hap])[2] ÷ 2  # base population size
    nSNP = size(base[:hap])[1]
    r  = par.nDam / par.nSire
    nSir = Int(floor(n₀/(1+r)))
    nDam = n₀ - nSir
    nSib = par.nDam * par.nSib ÷ nDam # sibship size
    nClg = par.nDam * par.nC7e ÷ nDam # number challenged
    nPrd = nSib - nClg
    tPrd = nDam * nPrd + (par.nG8n - 1) * (par.nSib - par.nC7e) * par.nDam
    tClg = nDam * nClg + (par.nG8n - 1) *  par.nC7e             * par.nDam
    gOne = Dict(:nSire => nSir,
                :nDam  => nDam,
                :nSib  => nSib,
                :nP8n  => nPrd, # n-production
                :nC7e  => nClg) # n-challenged
    nSNP, tPrd, tClg, (; gOne...)
end

"""
    function breeding_program(base, par)
---
Do the breeding with the given parameter `par`.
I always breed the first generation separately.
One is to scale, if necessary, the base to the size of nuclear population.
The other is for ease of merging generation data.
Disk I/O is avoided as much as possible.

## Memo
A full pedigree should include size(base) rows of `0 0`.
"""
function breeding_program(base, par)
    # Storage in memory
    nSnp, nPrd, nClg, pOne = determine_storage(base, par)
    gPrd = zeros(Int8, nSnp, nPrd)
    gClg = zeros(Int8, nSnp, nClg)
    Prd = DataFrame(g8n  = ones(Int8, nPrd), # generation
                    sex  = rand(0:1, nPrd),  # determine sex in very good advance
                    Sire = zeros(Int, nPrd),
                    Dam  = zeros(Int, nPrd),
                    tbv  = zeros(nPrd),
                    p7e  = zeros(nPrd)) # phenotype
    Clg = DataFrame(tbv = zeros(nClg),
                    p7e = zeros(nClg))

    ## QTL setup for this simulation
    qtl = sim_QTL(base, par.nQTL...)
    println()                   # breeding cycles begin below
    @info "F-1: Scaled from base"
    begin                       # generation one
        ped = random_mate(pOne.nSire, pOne.nDam, pOne.nSib)
        nxt = gdrop(base[:hap], ped, base[:r])
        #fill_storage(Prd, Clg, cur, qtl, par
        bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
    end

    # The rest generations
    pPls = begin
        tmp = Dict(:nSire => par.nSire,
                   :nDam  => par.nDam,
                   :nSib  => par.nSib,
                   :nP8n  => par.nSib - par.nC7e,
                   :nC7e  => par.nC7e)
        (; tmp...)
    end
    #for ig in 2:par.nG8n
    #    @info "F-$ig"
    #end
    #obs, nuclear = generation_one(base, par, qtl)
    #
    ## other generations
    #x₁, x₂, i₁, i₂ = test_sets(par.nDam, par.nSib, par.nC7e)
    #ns = par.nSire + par.nDam   # starts from base population size
    #np = length(x₁)             # No. of sibs for production in g8n 2+.
    #nn = par.nSire + par.nDam   # nuclear population size
    #for ig in 2:par.nG8n        # breeding and selection cycles
    #    @info "F-$ig"
    #    cur = length(obs[:p₁]) # ID of new g8n starts after this.
    #    obs = begin
    #        ped = random_mate(nuclear, par.nSire, par.nDam, par.nSib)
    #        nxt = gdrop(obs[:g₁], ped, base[:r])
    #        ped .+= ns          # to make it start from generation 0
    #        sex = rand(0:1, np)
    #        bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
    #        bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
    #        Dict(:g₁  => [obs[:g₁]   nxt[:, i₁]],
    #             :g₂  => [obs[:g₂]   nxt[:, i₂]],
    #             :bv₁ => [obs[:bv₁]; bv₁[x₁]],
    #             :bv₂ => [obs[:bv₂]; bv₂[x₂]],
    #             :p₁  => [obs[:p₁];  p₁[x₁]],
    #             :p₂  => [obs[:p₂];  p₂[x₂]],
    #             :n₂  => [obs[:n₂];  length(x₂)],
    #             :ped => [obs[:ped]; [ones(Int16, np).*ig sex ped[x₁, :]]]
    #             )
    #    end
    #    df = DataFrame(ID = cur+1:cur+np,
    #                   sex = rand(0:1, np),
    #                   val = obs[:p₁][cur+1:cur+np])
    #    nuclear = select_nuclear(df, par.nSire, par.nDam)
    #end
    #obs
end
