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
this function returns a `n × (v+1)` fixed effect matrix.
`+1` is for mean, `μ`.

In this simulation, the fixed effects are just the generation effects.
As selection on the binary trait goes on, the genetic values are elavated.
The challenge, however, is still killing half of the challenge population.

I may later replace this part using `MixedModels`.
"""
function fxmat(F)
    levels = unique(F)
    dic = begin
        tmp = Dict()
        col = 2
        for l in levels
            tmp[l] = col
            col += 1
        end
        tmp
    end
    nlvl = length(levels) + 1
    nobs = length(F)
    X = zeros(Float32, nobs, nlvl)
    X[:, 1] .= 1
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
"""
function snp_blup(g, p; σₐ² = 1, σₑ² = 1, Q = [], F = [])
    nlc, nid = size(g)
    @info join(["",
                "SNP-BLUP",
                "  $nlc SNP",
                "  $nid ID"], "\n")

    λ = (σₑ²/σₐ²) * nlc
    X = begin           # incident matrix for fixed effects, include μ
        if length(F) > 0
            fxmat(F)
        else
            ones(Float32, nid, 1)
        end
    end
    nf = size(X)[2]
    sl = nlc + nf               # size of lhs
    α, β = Float32[1, 0]        # scales for GEMM
    lhs = begin
        tmp = Array{Float32, 2}(undef, sl, sl)

        # upper left block
        ul = view(tmp, 1:nf, 1:nf)
        BLAS.gemm!('T', 'N', α, X, X, β, ul)

        # lower left block
        ll = view(tmp, nf+1:sl, 1:nf)
        BLAS.gemm!('N', 'N', α, g, X, β, ll)

        # upper right block
        for i in 1:nf
            for j in nf+1:sl
                tmp[i, j] = tmp[j, i]
            end
        end

        # lower right block
        lr = view(tmp, nf+1:sl, nf+1:sl)
        BLAS.gemm!('N', 'T', α, g, g, β, lr)
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
        tmp = Array{Float32, 1}(undef, sl)
        y = convert.(Float32, p)
        # Fixed effects
        uv = view(tmp, 1:nf)
        BLAS.gemv!('T', α, X, y, β, uv)

        # Random effecs (may be some fixed QTL also)
        lv = view(tmp, nf+1:sl)
        BLAS.gemv!('N', α, g, y, β, lv)
        tmp
    end

    # Solving
    LAPACK.posv!('L', lhs, rhs)

    # return fixed effects and SNP effects separately
    rhs[1:nf], rhs[nf+1:end]
end

# Other methods, e.g., bayes-α, β, γ, π, may be added.
