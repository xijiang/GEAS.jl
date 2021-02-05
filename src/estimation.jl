"""
    function grm(gt)
---
Genotypes are of value 0, 1 or 2, and of ID column majored.
They should be stored per ID per column.
This function uses vanRaden method I.
I do no genotype and frequency check here,
As it should be done somewhere else.
"""
function grm(g)
    twop = mean(g, dims = 2)
    s2pq = (1 .- .5 .* twop)'twop # one element array
    r2pq = 1. / s2pq[1]
    Z = g .- twop
    Z'Z .* r2pq
end

"""
    function grm_yang(gt)
---
GRM used by Yang et al. (2010), Nat Genet 42:565-571.
`gt` should be stored ID column majored.
"""
function grm_yang(g)
    nlc, nid = size(g)
    p = mean(g, dims = 2) ./ 2
    q = 1 .- p
    r = 1 ./ (2 .* p .* q)      # '×' is faster than '÷'

    d = zeros(nid)              # diagonals
    t₁ = 1 .+ 2 .* p
    t₂ = 2 .* p .* p
    for i in 1:nid
        w = view(g, :, i)
        d[i] = 1 + sum((w .* w - t₁ .* w + t₂) .* r) / nlc
    end

    t = 2 .* p
    r = sqrt.(r)
    z = (g .- t) .* r
    grm = z'z ./ nlc
    for i in 1:nid
        grm[i, i] = d[i]
    end
    grm
end

# below all return GEBV corresponding to `g`, or `g₁`.
function g_blup(g, p)
    @warn "Under construction"
end

function g_blup(g₁, g₂, p)
    @warn "Under construction"
end

function snp_blup(g, p; σₐ² = 1, σₑ² = 1)
    @warn "Under construction"
end

function snp_blup(g₁, g₂, p)
    @warn "Under construction"
end
# Other methods, e.g., bayes-α, β, γ, π, may be added.
