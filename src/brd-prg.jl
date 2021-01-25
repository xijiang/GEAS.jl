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
    snps = begin                              # nuclear pop genotypes of g-2
        # determine pedigree struct
        nfam = par.nn ÷ (1 + par.fm) * par.fm # full sib families
        nclg = nfam * par.nc                  # challenged ID
        noff = nfam * par.ss                  # selection and challenge canddt.
        n₀   = size(base[:hap])[2] ÷ 2        # base population size
        f1   = n₀ ÷ (1 + par.fm) * par.fm # full sib families in g₁
        ss   = noff ÷ f1                  # sibship size in g₁
        nc   = nclg ÷ f1           # sibs to challege in a sibship
        # simulation of pedigree
        ped  = random_mate(n₀, ss) # sample parents for offspring
        goff = gdrop(base[:hap], ped, r) # dropping -> offspring SNP
        tbv₁ = breeding_value(goff, qtl[1]) # breeding value of the prd trait
        tbv₂ = breeding_value(goff, qtl[2]) # breeding value of the bin trait
        pht₁ = phenotype(tbv₁, par.h²[1])
        pht₂ = phenotype(goff, qtl[2], par.h²[2], par.th)
        # estimation procedure
        # -- sample ID to challenge
        chg = zeros(Int, nc*f1)
        prd = zeros(Int, noff-chg)
    end

    #for ig in 2:par.ng
    #    println()
    #    @info join(["",
    #                "generation $ig"],
    #               "\n")
    #end
end
