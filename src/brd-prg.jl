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
    function scale_generation_one(par, n₀)
---
The base population size varies.
The rest of the generations have a constant size.
This function is to scale generation one, such that,
the production population and the binary population are
approximately equal to the later generations.
"""
function scale_generation_one(par, n₀)
    nfam = par.nn ÷ (1 + par.fm) * par.fm # full sib families in g 2+
    nclg = nfam * par.nc                  # challenged ID
    noff = nfam * par.ss                  # selection and challenge canddt.
    nf   = n₀ ÷ (1 + par.fm) * par.fm     # full sib families in g₁
    ss   = noff ÷ nf                      # full sibship size in g₁
    nc   = nclg ÷ nf            # sibs to challege in a sibship
    nf, ss, nc
end

"""
    function breed_n_select(snp, par, qtl, r)
---
Given the genotypes of a nuclear population `snp`, breeding parameter `par`, 
qtl loci and effects `qtl`, and recombination rates `r`,
this function return genotypes of its next generation,
and two vectors of phenotypes.
"""
function breed_n_measure(snp, par, qtl, r)
    ped  = random_mate(par.nn, par.ss) # sample parents for offspring
    goff = gdrop(snp, ped, r)          # dropping -> offspring SNP
    p₁ = phenotype(goff, qtl[1], par.h²[1])
    p₂ = phenotype(goff, qtl[2], par.h²[2], par.th)
    goff, p₁, p₂
end

"""
    function evaluate_n_select()
---
"""
function evaluate_n_select(snp, ix₁, ix₂, prd, bin, nn, fm)
    @warn "under construction"
    snp[:, 1:2nn]
end

"""
    function breeding_program(base, par)
---
Do the breeding with the given parameter `par`.
"""
function breeding_program(base, par)
    # Parameters for this simulation
    qtl = sim_QTL(base, par.nq...)
    r   = haldane(base[:pos])
    println()
    
    @info join(["",
                "generation 1",
                "  - scale nuclear population size to par.nn",
                "  - Assuming sizes for challenging, genotyping are similar"],
               "\n")
    snp = begin                      # nuclear pop genotypes of g-2
        n₀ = size(base[:hap])[2] ÷ 2 # base population size
        nf, ss, nc = scale_generation_one(par, n₀)
        tp = Breeding(n₀, par.fm, par.ng, ss, par.nq, par.h², nc, par.th)
        goff, prd, bin = breed_n_measure(base[:hap], tp, qtl, r)
        i₁, i₂ = test_sets(nf, ss, nc)
        evaluate_n_select(goff, i₁, i₂, prd, bin, par.nn, par.fm)
    end

    nf = par.nn ÷ (1 + par.fm) * par.fm
    ix₁, ix₂ = test_sets(nf, par.ss, par.nc)
    for ig in 2:par.ng
        println()
        @info join(["",
                    "generation $ig"],
                   "\n")
        goff, prd, bin = breed_n_measure(snp, par, qtl, r)
        snp = evaluate_n_select(goff, ix₁, ix₂, prd, bin, par.nn, par.fm)
        # gene editing goes here?
        @warn "record something"
    end
end
