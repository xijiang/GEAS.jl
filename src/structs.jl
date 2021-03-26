"""
QTL parameters.
- **pos**: QTL location indices of the SNPs
- **effect**: effects of allele `1`, of each QTL
- **mean::Float64**: mean of the base population
- **max::Float64**: TBV of an ideal ID from the base
"""
struct QTL
    pos
    effect
    mean::Float64               # mean and max are about the base population
    max::Float64
end

#=
"""
    struct Breeding
---
    nn::Int                     # size of nuclear population, a multiple of (1+fm)
    fm::Int                     # female:male in nuclear, default 2
    ng::Int                     # number of generations
    ss::Int                     # full sibship size
    nq                          # number of QTL vector
    h²                          # heritability vector
    # for the binary trait
    nc::Int                     # number of sibs for challenge test
    th::Float64                 # threshold for the binary trait. e.g., 0 -> half incdn.

## Example
- `par = GEAS.Breeding(375, 2, 10, 100, [500, 500], [.5, .5], 40, 0.)`
- `par.nn` to retrieve the nuclear population size
"""
struct Breeding
    nn::Int                     # size of nuclear population, a multiple of (1+fm)
    fm::Int                     # female:male in nuclear, default 2
    ng::Int                     # number of generations
    ss::Int                     # full sibship size
    nq                          # number of QTL vector
    h²                          # heritability vector
    # for the binary trait
    nc::Int                     # number of sibs for challenge test
    th::Float64                 # threshold for the binary trait. e.g., 0 -> half incdn.
end

"""
Data for genomic selection, SNP BLUP method.
- **g**: genotypes of value 0, 1, or 2, of nSNP × nID
- **p**: phenotypes, nID × 1
- **a**: overall genetic variance
- **e**: overall environmental variance
- **s**: SNP loci to be fitted as fixed effects
- **f**: fixed effect vector with integer levels
- **w**: weight on this trait
"""
struct GS_dat
    g                           # genotypes of value 0, 1, or 2, of nSNP × nID
    p                           # phenotypes, nID × 1
    a                           # overall genetic variance
    e                           # overall environmental variance
    s                           # SNP loci to be fitted as fixed effects
    f                           # fixed effect vector with integer levels
    w                           # weight on this trait
end
=#
