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
    function snp_blup(g, p, h²; Q = [], F = [], dd=0)
---
The plain SNP-BLUP procedure.
Returns E(population) and  a vector of SNP effects.

**Warning**: I use Float32 to save memory.
This will be changed if to run on a bigger machine.
"""
function snp_blup(g, p, h²; Q = [], F = [], dd=0)
    nlc, nid = size(g)
    λ = (1. - h²)/h² * nlc + dd
    X = begin           # incident matrix for fixed effects, include μ
        if length(F) > 0
            fxmat(F)
        else
            ones(Int16, nid, 1)
        end
    end
    nf = size(X)[2]
    sl = nlc + nf               # size of lhs

    lhs = zeros(Float32, sl, sl)
    # -- upper left block
    ul = view(lhs, 1:nf, 1:nf)
    matmul!(ul, X', X)
    # -- lower left block
    ll = view(lhs, nf+1:sl, 1:nf)
    matmul!(ll, g, X)
    # -- upper right block
    ur = view(lhs, 1:nf, nf+1:sl)
    copyto!(ur, ll')
    # -- lower right block
    lr = view(lhs, nf+1:sl, nf+1:sl)
    matmul!(lr, g, g')
    for i in nf+1:sl
        lhs[i, i] += λ
    end
    # -- if some QTL as fixed effects
    if length(Q) > 0
        l = Q .+ nf
        for x in l
            lhs[x, x] -= λ
        end
    end
    
    rhs = zeros(Float32, sl)
    y = reshape(p, length(p), :)
    rhs[1:nf]     = matmul(X', y)
    rhs[nf+1:end] = matmul(g, y)

    # Solving
    LAPACK.posv!('L', lhs, rhs)

    # return fixed effects and SNP effects separately
    rhs[1:nf], rhs[nf+1:end]
end
