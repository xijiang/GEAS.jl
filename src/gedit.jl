"""
    function gedit(svc, sr, tov)
---
Change alleles in SNP vector `svc` to `tov` (0/1),
with a successful rate, `sr`.
"""
function gedit(snp, prd, l, sr)
    @debug "editing" l
    tov = l > 0 ? 1 : 0
    fra, til = 2first(prd.id)-1, 2last(prd.id)
    svc = view(snp, abs(l), fra:til)
    br = rand(Bernoulli(sr), length(svc)) # succeeded editing
    svc[br] .= tov
end

"""
    function rank_QTL(snp, qtl; onVar=false)
---
Rank QTL on either `2pqa^2` or their absolute effects, in descending order.
If ranking on effects, allele frequecies of ≤0.1, or ≥0.9 are put
at the bottom.
If a position is positive, then allele `1` has positive effect.
Otherwise, negagive.
"""
function rank_QTL(snp, qtl; onVar=false)
    t = 0.1                     # threshold when ranking on abs effect
    p = vec(mean(snp, dims=2))
    e = abs.(qtl.effect)
    if onVar
        q = 1 .- p
        v = 2 .* p[qtl.pos] .* q[qtl.pos] .* (e .^ 2)
        o = sortperm(v, rev=true) # order
        return Int.(sign.(qtl.effect[o])) .* qtl.pos[o]
    else
        for i in length(e)
            xp = p[qtl.pos[i]]
            xp <= t || xp >= 1-t && (e[i] = 0)
        end
        o = sortperm(e, rev=true)
        return Int.(sign.(qtl.effect[o])) .* qtl.pos[o]
    end
end
