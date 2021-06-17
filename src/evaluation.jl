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
    function fxmat(F)
---
If the length of categorical vector `F` is `n`, which has `v` leveles,
this function returns a `n × v` fixed effect matrix.
`μ` is thus not fitted.

In this simulation, the fixed effects are just the generation effects.
As selection on the binary trait goes on, the genetic values are elavated.
The challenge, however, is still killing half of the challenge population.

I may later replace this part using `MixedModels`.
"""
function fxmat(F)
    levels = unique(F)
    dic = begin
        tmp = Dict()
        col = 1
        for l in levels
            tmp[l] = col
            col += 1
        end
        tmp
    end
    nlvl = length(levels)
    nobs = length(F)
    X = zeros(Int16, nobs, nlvl) # I bet levels are no more than 32,767
    for i in 1:length(F)
        X[i, dic[F[i]]] = 1
    end
    X
end

"""
    function snp_blup(g, p; σₐ² = 1, σₑ² = 1, Q = [], F = [])
---
The plain SNP-BLUP procedure.
Returns E(population) and  a vector of SNP effects.

**Warning**: I use Float32 to save memory.
This will be changed if to run on a bigger machine.

`function snp_blup(g, p; σₐ² = 1, σₑ² = 1, Q = [], F = [])`
changed to
`function snp_blup(g, p, h²; Q = [], F = [])`, on 2021-Apr.-1
"""
function snp_blup(g, p, h²; Q = [], F = [])
    nlc, nid = size(g)
    @info join(["",
                "SNP-BLUP",
                "  $nlc SNP",
                "  $nid ID"], "\n")

    # λ = (σₑ²/σₐ²) * nlc  # changed to below
    λ = (1. - h²)/h² * nlc
    X = begin           # incident matrix for fixed effects, include μ
        if length(F) > 0
            fxmat(F)
        else
            ones(Int16, nid, 1)
        end
    end
    nf = size(X)[2]
    sl = nlc + nf               # size of lhs
    #α, β = Float32[1, 0]        # scales for GEMM
    lhs = begin
        tmp = Array{Float32, 2}(undef, sl, sl)

        # upper left block
        ul = view(tmp, 1:nf, 1:nf)
        matmul!(ul, X', X)

        # lower left block
        ll = view(tmp, nf+1:sl, 1:nf)
        matmul!(ll, g, X)

        # upper right block
        for i in 1:nf
            for j in nf+1:sl
                tmp[i, j] = tmp[j, i]
            end
        end

        # lower right block
        lr = view(tmp, nf+1:sl, nf+1:sl)
        matmul!(lr, g, g')
        for i in nf+1:sl
            tmp[i, i] += λ
        end

        # if some QTL as fixed effects
        if length(Q) > 0
            l = Q .+ nf
            for x in l
                tmp[x, x] -= λ
            end
        end
        tmp
    end
    rhs = begin
        tmp = zeros(Float32, sl)
        y = reshape(p, length(p), :)
        tmp[1:nf]     = matmul(X', y)
        tmp[nf+1:end] = matmul(g, y)
        vec(tmp)
    end

    # Solving
    # ToDo: may try `preconditioned conjugate gradient` later
    LAPACK.posv!('L', lhs, rhs)

    # return fixed effects and SNP effects separately
    rhs[1:nf], rhs[nf+1:end]
end

# Other methods, e.g., bayes-α, β, γ, π, may be added.
#=
"""
    function evaluate(gsdat)
---
Vector `gsdat` has data of `1+` traits for evaluation with SNP BLUP method.
Each trait has `g`enotypes, `p`henotypes, trait `σₐ²`, `σₑ²`, 
`w`eight for a select index, `f`ixed effect vector, and `s`NP to be
fitted as fixed effects.

Returns `EBV(1)⋅w₁ + EBV(2)⋅w₂ + ...` of ID in `g₁` as a vector
"""
function evaluate(gsdat)
    t = gsdat[1]                # shorthand for production data
    f1, s1 = snp_blup(t.g, t.p, σₐ² = t.a, σₑ² = t.e, Q = t.s, F = t.f)
    w1 = t.w
    t = gsdat[2]                # shorthand for challenge data
    f2, s2 = snp_blup(t.g, t.p, σₐ² = t.a, σₑ² = t.e, Q = t.s, F = t.f)
    t = gsdat[2]
    w2 = t.w
    
end
=#
