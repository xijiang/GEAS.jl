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

"""
    function evaluation(g₁, p₁, g₂, p₂, qlc, h²)
---
- ₁: for production trait
- ₂: for binary trait

Vector `qlc` specify known QTL to be fitted as fixed effects.

The function returns a vector of score for ID in `₁`.
"""
function evaluation(g₁, p₁, q₁, h₁², g₂, p₂, q₂, h₂²)
end

# below all return GEBV corresponding to `g`, or `g₁`.
function g_blup(g, p)
    @warn "Under construction"
end

function g_blup(g₁, g₂, p)
    @warn "Under construction"
end

"""
    function snp_blup(g, p, q; σₐ² = 1, σₑ² = 1)
---
With SNP genotypes `g`, phenotypes `p`, and QTL SNP `q`,
this function return a vector of SNP effects with SNP-BLUP method.

SNP specified in `q` is fitted as fixed effects.
"""
function snp_blup(g, p, q; σₐ² = 1, σₑ² = 1)
    
end

"""
    function snp_blup(g, p, q, f; σₐ² = 1, σₑ² = 1)
---
With SNP genotypes `g`, phenotypes `p`, QTL SNP `q`, and filial `f`,
this function return a vector of SNP effects with SNP-BLUP method.

SNP specified in `q` is fitted as fixed effects.
Filial/generation `f` is fitted as fixed effects.
"""
function snp_blup(g, p, q, f; σₐ² = 1, σₑ² = 1)
    @warn "Under construction"
end

"""
    function snp_blup(g, p; σₐ² = 1, σₑ² = 1)
---
The plain SNP-BLUP procedure.
Returns E(population) and  a vector of SNP effects.
"""
function snp_blup(g, p; σₐ² = 1, σₑ² = 1)
    nlc, nid = size(g)
    @info join(["",
                "SNP-BLUP",
                "  $nlc SNP",
                "  $nid ID"], "\n")

    λ = (σₑ²/σₐ²) * nlc

    lhs = begin
        edge = sum(g, dims=2)
        tmp = Array{Float32, 2}(undef, nlc+1, nlc+1)
        tmp[1, 1] = nid
        tmp[1, 2:end] = edge'
        tmp[2:end, 1] = edge
        tmp[2:end, 2:end] = g * g' + I .* λ
        tmp
    end

    rhs = begin
        tmp = Array{Float32, 1}(undef, nlc + 1)
        tmp[1] = sum(p)
        tmp[2:end] = g * p
        tmp
    end
    sol = lhs \ rhs
    sol[1], sol[2:end]
end
# Other methods, e.g., bayes-α, β, γ, π, may be added.
