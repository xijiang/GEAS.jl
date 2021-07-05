"""
    function gedit(snp, pos, sr)
---
With the SNP genotypes `snp` on position `pos`.
If `pos`>0, make them alleles at `pos` 1, or 0, with a successful rate, `sr`.
"""
function gedit(snp, pos, sr)
    nid = size(snp)[2] ÷ 2
    l = abs(pos)
    br = rand(Bernoulli(sr), 2nid) # succeeded editing
    if pos > 0
        snp[l, br] .= 1
    else
        snp[l, br] .= 0
    end
end

"""
    function top_QTL(snp, qtl, n, onVar=false)
---
Return the row number of top `n` QTL `qtl` in `snp`, on variance (`2pqa^2`)
or on `effect`.

## Amendment
- I realized that it's better to decide the QTL and their order to edit in the beginning.
- will return positions of the QTL
- if a position is positive, then allele `1` has positive effect.  Otherwise, negagive.
"""
function top_QTL(snp, qtl, n; onVar=false)
    p = vec(mean(snp, dims=2))
    x = qtl.pos                 # short hands
    e = qtl.effect

    if onVar
        q = 1 .- p
        v = 2 .* p[x] .* q[x] .* (e .^ 2) # 2pqa^2
        o = sortperm(v, rev=true)
        return Int.(x[o[1:n]] .* sign.(e[o[1:n]]))
    else
        o = sortperm(abs.(e), rev=true)
        r = Int[]               # results
        for i in 1:length(o)
            j = o[i]
            k = x[j]            # real SNP position
            p[k] <= 0.1 && continue
            p[k] >= 0.9 && continue
            push!(r, Int(sign(e[j]) * k))
            (length(r) == n) && (return r)
        end
    end
end
