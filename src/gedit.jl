"""
    function gedit(snp, qtl, sr)
---
With the SNP genotypes `snp`, and `qtl` information,
this function find the the biggest QTL of this genotype array,
according to `2pqa^2`.
The editing will have a successful rate, `sr` ∈ (0, 1).

For simplicity, the whole generation genotypes can be sent here for editing.
"""
function gedit(snp, qtl, sr)
    nid = size(snp)[2] ÷ 2
    nq  = length(qtl.pos)
    pos = qtl.pos
    ep = top_QTL(snp, qtl)[1]   # position for editing
    br = rand(Bernoulli(sr), 2nid) # succeeded editing
    if qtl.effect[ep] < 0       # negative, 1 -> 0
        snp[ep, br] .= 0
    else                        # positive, 0 -> 1
        snp[ep, br] .= 1
    end
end

"""
    function top_QTL(base, qtl; n = 1)
---
Return the subset of top QTL `qtl` according to their variance, `2pqa^2`, in `base`.
Only one is returned by default.
`n` can be changed to other value, as number of known QTL in the base population.
"""
function top_QTL(snp, qtl; n = 1)
    gt = qtl_gt(snp, qtl.pos)
    p = mean(gt, dims=2)
    q = 1 .- p
    v = 2 .* p .* q .* (qtl.effect .^ 2) # 2pqa^2
    r = sortperm(vec(v), rev=true)[1:n]
end
