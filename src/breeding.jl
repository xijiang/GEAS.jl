"""
    function rand_mate(Sires, Dams, nSib)
---
Ramdomly mate `Sires` and `Dams`, which are name vectors.
Each family will have `nSib` full sibs.
If the bigger to smaller ratio of length of `Sires` and `Dams` is
not integer, some ID may be used more often than others.
"""
function rand_mate(Sires, Dams, nSib)
    np, nm = length(Sires), length(Dams) # number of pa and ma
    nf = max(np, nm)                     # number of full sib families
    nid = nSib * nf
    pa = repeat(shuffle(Sires), inner=nSib, outer=Int(ceil(nf/np)))
    ma = repeat(shuffle(Dams), inner=nSib, outer=Int(ceil(nf/nm)))
    pa[1:nid], ma[1:nid]
end

"""
    function create_storage(base, par)
---
This is a refactorized function for creating storage.
My function before was too long, making it difficult for debugging.
This function is still very long, but well blockwised.

The storage is created in one-go.
This is to avoid too many times of memory allocations,
which can crash the program.

Finally knew how to use `view`.
Use `view`s for gene-dropping later.

This function:
1. create storage for SNP of production and challenge group
   - base SNP are copied to the production group as generation 0.
2. create pedigree generation 0 → 1 for the production group
   - which include the `base` population.
3. create pedigree for the challenge group.
4. generations are marked
5. sex sampling determined in good advance
6. determine parents of the 1st generation.
"""
function create_storage(base, par)
    # the pedigree are of 3 parts: the base, 1st generation, 2:end generations
    nSnp = size(base.hap)[1]

    # Generation 0/base: (for simplicity, I copy base SNP into `snp_prd` later)
    n₀    = size(base.hap)[2] ÷ 2   # base population size
    fmr   = par.nDam / par.nSire      # female:male ratio
    n₀Sir = Int(ceil(n₀ / (1 + fmr))) # valid even when nSires > nDams
    n₀Dam = n₀ - n₀Sir
    g₀Sex = shuffle([ones(Int8, n₀Sir); zeros(Int8, n₀Dam)])
    g₀Sir = (1:n₀)[g₀Sex .== 1]
    g₀Dam = (1:n₀)[g₀Sex .== 0]
    
    # Generation 2:end
    n₂Fam = max(par.nDam, par.nSire) # in generation 2:end
    # - nID in each of generation 2:end
    n₂Prd, n₂Clg = n₂Fam .* [par.nSib - par.nC7e, par.nC7e]
    
    # Generation 1: generate similar number of offspring in 2+ generation
    n₁Fam = max(n₀Sir, n₀Dam)   # full sib groups from base
    # - sibship sizes in 1st generation
    n₁P8n, n₁C7e = Int.(ceil.([n₂Prd, n₂Clg]./n₁Fam))
    # - total offspring in 1st generation
    n₁Prd, n₁Clg = [n₁P8n, n₁C7e] .* n₁Fam
    
    # Overall number of offspring
    nPrd, nClg = [n₁Prd + n₀, n₁Clg] + [n₂Prd, n₂Clg] .* (par.nG8n - 1)

    # storage For SNP, after the sizes are determined above
    snp = begin
        prd = BitArray(undef, nSnp, 2nPrd)
        snp₀ = view(prd, :, 1:2n₀)
        copy!(snp₀, base.hap)     # the SNP copy
        clg = BitArray(undef, nSnp, 2nClg)
        dic = Dict(:prd => prd, :clg => clg)
        (; dic...)
    end

    # For pedigree
    ped = begin
        prd = begin
            pa, ma = rand_mate(g₀Sir, g₀Dam, n₁P8n)
            DataFrame(
                g8n = Int16.([zeros(n₀);
                              ones(n₁Prd);
                              repeat(2:par.nG8n, inner = n₂Prd)]),
                id  = Int32.(1:nPrd),
                sex = [g₀Sex; rand(Int8.(0:1), nPrd - n₀)],
                sir = Int32.([zeros(Int32, n₀); pa; zeros(Int32, nPrd-n₀-n₁Prd)]),
                dam = Int32.([zeros(Int32, n₀); ma; zeros(Int32, nPrd-n₀-n₁Prd)]),
                tbv = zeros(Float32, nPrd), # TBV for the production trait
                p7e = zeros(Float32, nPrd), # production phenotypes
                val = zeros(Float32, nPrd)  # selection index
            )
        end
        clg = begin
            pa, ma = rand_mate(g₀Sir, g₀Dam, n₁C7e)
            DataFrame(
                g8n = Int16.([ones(n₁Clg); repeat(2:par.nG8n, inner = n₂Clg)]),
                id  = Int32.(1:nClg),
                sir = Int32.([pa; zeros(Int32, nClg - n₁Clg)]),
                dam = Int32.([ma; zeros(Int32, nClg - n₁Clg)]),
                tbv = zeros(Float32, nClg),
                p7e = zeros(Float32, nClg)
            )
        end
        dic = Dict(:prd => prd, :clg => clg)
        (; dic...)
    end
    snp, ped
end

"""
    function sample_alleles(r)
---
This function returns a vector of 1:2 of the same size of `r`,
which is the vector recombination rates.
`r` for the first locus of each chromosome is 0.5,
such that they are independent.
"""
function sample_alleles(r)
    n = length(r)
    t = rand(n)
    v = ones(Int, n)
    k = 1
    @inbounds for i in 1:n
        (t[i] < r[i]) && (k = 3-k)
        v[i] = k
    end
    v
end

"""
    function gene_drop(df, r, psnp, osnp, qtl, h²)
---
Given a dataframe `df` of a generation, 
this function drop genes in `prt_snp` specified by `df.sir` and `df.dam`
int `off_snp` specified by `df.id`.
When the genotypes of this genration are created,
the breeding values, and phenotypes are also calculated.
"""
function gene_drop(df, r, psnp, osnp, qtl, h²)
    Threads.@threads for x in eachrow(df)
        gp = view(psnp, :, 2x.sir-1:2x.sir) # pa's genotypes
        gm = view(psnp, :, 2x.dam-1:2x.dam) # ma's genotypes
        op = view(osnp, :, 2x.id-1)         # id's haplotype from pa
        om = view(osnp, :, 2x.id)           # id's haplotype from ma
        ix = sample_alleles(r)
        copy!(op, [gp[i, ix[i]] for i in 1:length(op)])
        ix = sample_alleles(r)
        copy!(om, [gm[i, ix[i]] for i in 1:length(om)])
    end
    fra = 2first(df.id) - 1      # from column in SNP genotypes
    til = 2last(df.id)           # end column
    snp = view(osnp, :, fra:til)
    df[:, :tbv], df[:, :p7e] = phenotype(snp, qtl, h²)
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
    sires = last(sort(gp[(sex=1,)], :val), nSire).id
    dams  = last(sort(gp[(sex=0,)], :val), nDam ).id
    return sort(sires), sort(dams)
end

"""
    function calc_idx(ped, snp, cur, lqtl, par)
---
Calculate selection index SNP and phenotype information up to current generation.
Store the results in `cur[:, :val]`.

- Generations are treated as a fixed effect in the challenge group.
- If `lqtl` is not empty, then its specified loci are fixed as fixed effects
"""
function calc_idx(ped, snp, cur, Q, par)
    fra = findlast(x -> x == 0, ped.prd.g8n)
    til = last(cur).id
    gp = alleles2gt(view(snp.prd, :, 2fra+1:2til))
    mp, sp = snp_blup(gp, ped.prd[:, :p7e][fra+1:til], par.h²[1], dd = par.dd)

    til = findlast(x -> x == first(cur).g8n, ped.clg.g8n)
    gc = alleles2gt(view(snp.clg, :, 1:2til))
    F = ped.clg.g8n[1:til]
    mc, sc = snp_blup(gc, ped.clg[:,:p7e][1:til], par.h²[2], Q=Q, F=F, dd=par.dd)

    tid = size(gp)[2]
    cgp = view(gp, :, tid-nrow(cur)+1:tid) # genotypes of current generation
    cur[:, :val] = cgp'sp + cgp'sc .* par.w4t
end

"""
    function assign_pama(ped, sires, dams, nsib)
---
Randomly mate `sires` and `dams` according to the sibship sizes `nsib`,
and assign them into the parent columns of `ped`.
"""
function assign_pama(ped, sires, dams, nsib)
    pa, ma = rand_mate(sires, dams, nsib)
    ped[:, :sir] = pa
    ped[:, :dam] = ma
end

"""
    function breeding(base, par, qtl; edit = false, fixed = false)
---
The breeding program.
- `edit` is to instruct whether to edit top QTL or not.
- `fixed` is to determine whether to trait the top QTL as fixed effect or not.
"""
function breeding(base, par, qtl; edit = false, fixed = false)
    rkq = rank_QTL(base.hap, qtl[2]) # QTL rank in ↓ order
    Q = fixed ? abs.(rkq[1:par.nK3n]) : Int[] # whether to emphaize known QTL or not
    @debug rkq[1:10]                    # print top QTL at debug level
    # prepare pedigree and storage
    snp, ped = create_storage(base, par)
    prd, clg = groupby(ped.prd, :g8n), groupby(ped.clg, :g8n)
    for ig in 1:par.nG8n-1
        @info "Breeding generation $ig of $(par.nG8n)"
        # dropping and create geno-, pheno- types.
        gene_drop(prd[(g8n=ig,)], base.r, snp.prd, snp.prd, qtl[1], par.h²[1])
        gene_drop(clg[(g8n=ig,)], base.r, snp.prd, snp.clg, qtl[2], par.h²[2])

        edit && gedit(snp.prd, prd[(g8n=ig,)], popfirst!(rkq), par.e19e)

        # - evaluation, saving the results in prd.val
        calc_idx(ped, snp, prd[(g8n=ig,)], Q, par)
        # some shortcuts below
        #prd[(g8n=ig,)][:, :val] = prd[(g8n=ig,)][:, :tbv]

        # - select breeders for next generation, and drop genes
        sires, dams = select_nuclear(prd[(g8n=ig,)], par.nSire, par.nDam)

        # - assign parents of the next generation
        assign_pama(prd[(g8n=ig+1,)], sires, dams, par.nSib - par.nC7e)
        assign_pama(clg[(g8n=ig+1,)], sires, dams, par.nC7e)
    end
    @info "Generation: $(par.nG8n), dropping only in the last generation"
    gene_drop(prd[(g8n=par.nG8n,)], base.r, snp.prd, snp.prd, qtl[1], par.h²[1])
    return ped.prd, snp.prd         # only needed for summary
end
