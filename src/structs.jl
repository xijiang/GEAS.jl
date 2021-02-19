struct QTL
    pos
    effect
    rank                        # rand high -> low on 2pqa^2
    mean::Float64               # mean and max are about the base population
    max::Float64
end

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
