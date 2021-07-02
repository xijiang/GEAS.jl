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
    # Index for ID in, e.g., phenotypes
    ic = begin
        a = repeat(0:nsib:(nfam-1)*nsib, inner=nclg)
        a .+= repeat(1:nclg, nfam)
    end
    ip = begin
        n = nsib - nclg
        a = repeat(0:nsib:(nfam-1)*nsib, inner=n)
        a .+= repeat(nclg+1:nsib, nfam)
    end
    
    # Index for haplotypes
    hp = 2repeat(ip, inner=2)
    hp[1:2:end] .-= 1
    hc = 2repeat(ic, inner=2)
    hc[1:2:end] .-= 1
    ip, ic, hp, hc
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
    sires = last(sort(gp[1], :val), nSire).id
    dams  = last(sort(gp[2], :val), nDam ).id
    return [sires; dams]
end

"""
    function create_storage(base, par)
---
Determine array sizes for one-go genotype storage.
This is to avoid too much memory allocation, 
and garbage collection.
GC may also cause program crash.
"""
function create_storage(base, par, qtl)
    # Generation 0:
    n₀   = size(base[:hap])[2] ÷ 2 # base population size
    nSnp = size(base[:hap])[1]
    nSire = begin
        r = par.nDam / par.nSire
        Int(floor(n₀/(1+r)))
    end
    nDam = n₀ - nSire
    
    # Generation 1:
    nSib = par.nDam * par.nSib ÷ nDam # sibship size
    nClg = par.nDam * par.nC7e ÷ nDam # number challenged
    nPrd = nSib - nClg
    ped = random_mate(nSire, nDam, nSib)
    nxt = gdrop(base[:hap], ped, base[:r])
    ip, ic, hp, hc = test_sets(nDam, nSib, nClg)
    bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
    bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)

    # Storage, with simulated generation one included.
    tPrd = nDam * nPrd + (par.nG8n - 1) * (par.nSib - par.nC7e) * par.nDam
    tClg = nDam * nClg + (par.nG8n - 1) *  par.nC7e             * par.nDam

    #- SNP genotypes
    snp₁ = BitArray(undef, nSnp, tPrd*2)
    snp₂ = BitArray(undef, nSnp, tClg*2)
    snp₁[:, 1:length(hp)] = nxt[:, hp]
    snp₂[:, 1:length(hc)] = nxt[:, hc]

    #- Produciton population
    prd = begin
        t = (par.nSib - par.nC7e) * par.nDam
        g8n = Int8.([ones(nPrd * nDam); repeat(2:par.nG8n, inner=t)])
        id  = 1:tPrd
        sex = rand(Int8.(0:1), tPrd) # determine sex in very good advance
        Sir = Array{Int32, 1}(undef, tPrd)
        Dam = Array{Int32, 1}(undef, tPrd)
        TBV = Array{Float32, 1}(undef, tPrd)
        T2C = Array{Float32, 1}(undef, tPrd) # TBV for disease trait
        p7e = Array{Float32, 1}(undef, tPrd)
        n₁  = length(ip)        # fill in generation one info
        # !!!! OBS !!! Pedigree below
        #     if to calculate A matrix, add n_base to 2+ generations.
        Sir[1:n₁] = ped[ip, 1]
        Dam[1:n₁] = ped[ip, 2]
        TBV[1:n₁] = bv₁[ip]
        T2C[1:n₁] = bv₂[ip]
        p7e[1:n₁] = p₁[ip]
        DataFrame(g8n = g8n,
                  id  = id,
                  sex = sex,
                  sir = Sir,
                  dam = Dam,
                  tbv = TBV,    # TBV for the production trait
                  t2c = T2C,    # TBV for the binary trait
                  p7e = p7e)    # phenotype for the production trait
    end

    #- The challenged ones
    clg = begin
        t = par.nC7e * par.nDam
        g8n = Int8.([ones(nClg * nDam); repeat(2:par.nG8n, inner = t)])
        TBV = Array{Float32, 1}(undef, tClg)
        p7e = Array{Float32, 1}(undef, tClg)
        n₂ = length(ic)         # fill in info of generation one
        TBV[1:n₂] = bv₂[ic]
        p7e[1:n₂] = p₂[ic]
        DataFrame(g8n = g8n, tbv = TBV, p7e = p7e)
    end
    snp₁, snp₂, prd, clg, length(ip), length(ic), length(hp), length(hc)
end

"""
    function breeding_program(base, par, nsib; edit=false)
---
Do the breeding with the given parameter `par`.
I always breed the first generation separately.
One is to scale, if necessary, the base to the size of nuclear population.
The other is for ease of merging generation data.
Disk I/O is avoided as much as possible.

## Memo
A full pedigree should include size(base) rows of `0 0`.
"""
function breeding_program(base, par, qtl; edit=false, fixed=false)
    flog = open(par.log, "w")
    # determine which QTL are known in the beginning.
    loci_2_edit = top_QTL(base[:hap], qtl[2], n=par.nk3n)
    println(flog, join(loci_2_edit, ' '))
    ed_qtl = QTL(qtl[2].pos[loci_2_edit], qtl[2].effect[loci_2_edit], 0, 0)
    Q = fixed ? loci_2_edit : Int[] # whether to deem the known QTL as fixed effects
    
    # general parameters
    @info "Create storage and generation one"
    snp₁, snp₂, prd, clg, f1, f2, f3, f4 = # f: from
        create_storage(base, par, qtl)
    # shorthands
    nSire, nDam, nSib, nChallenge, percentage, generation =
        par.nSire, par.nDam, par.nSib, par.nC7e, par.p8e, par.nG8n
    ip, ic, hp, hc = test_sets(nDam, nSib, nChallenge)
    l1, l2, l3, l4 = length.((ip, ic, hp, hc)) # lengths of generation 2+
    h2 = begin                                 # for the binary trait
        p = par.p8e             # percentage to be killed in challenge
        z = pdf(Normal(), p)
        par.h²[2]*z*z/(p*(1-p)) # ~0.25 here
    end
    for ig in 2:generation
        @info "generation $ig"
        # prepare for evaluation
        
        # OBS!!!
        # select ID for nuclear from ig-1 on data of all previous generations
        val = begin             # scores for selective candidates
            # edit or emphasize on known QTL
            edit && gedit(snp₁, ed_qtl, par.e19e)
            # the production trait
            g1 = alleles2gt(view(snp₁, :, 1:f3))
            p1 = select(prd[prd.g8n .< ig, :], :p7e).p7e
            m1, s1 = snp_blup(g1, p1, par.h²[1], dd = par.dd)
            # the challenge trait
            g2 = alleles2gt(view(snp₂, :, 1:f4))
            p2 = select(clg[clg.g8n .< ig, :], :p7e).p7e
            m2, s2 = begin
                # generation as fixed effects
                f = select(clg[clg.g8n .< ig, :], :g8n).g8n
                snp_blup(g2, p2, h2, Q=Q, F=f, dd = par.dd)
            end
            g1's1 + g1's2 .* par.w4t
        end
        
        df = select(prd[(prd.g8n .== ig-1), :], :id, :sex)
        tl = length(df.id) - 1
        df.val = vec(val)[end-tl:end]
        # or, select(prd[(prd.g8n .== ig-1), :], :id, :sex, :tbv)
        nuclear = select_nuclear(df, nSire, nDam)
        
        # simulate current generation
        ped = random_mate(nuclear, nSire, nDam, nSib)
        nxt = gdrop(snp₁, ped, base[:r])
        bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
        
        # update storage.  codes below are ugly.
        prd.tbv[f1+1:f1+l1] = bv₁[ip]
        prd.t2c[f1+1:f1+l1] = bv₂[ip]
        prd.p7e[f1+1:f1+l1] =  p₁[ip]
        prd.sir[f1+1:f1+l1] = ped[ip, 1]
        prd.dam[f1+1:f1+l1] = ped[ip, 2]
        clg.tbv[f2+1:f2+l2] = bv₂[ic]
        clg.p7e[f2+1:f2+l2] =  p₂[ic]
        snp₁[:, f3+1:f3+l3] = nxt[:, hp]
        snp₂[:, f4+1:f4+l4] = nxt[:, hc]
        f1 += l1
        f2 += l2
        f3 += l3
        f4 += l4
    end
    close(flog)
    return prd, snp₁            # bitset snp is only of a few GiB.
end
