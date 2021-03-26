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
    ep = begin                  # position for editing
        Q = zeros(Int8, nq, nid) # QTL genotypes
        for id in 1:nid
            Q[:, id] = snp[pos, id*2 - 1] + snp[pos, id*2]
        end
        p = mean(Q, dims=2) ./2
        q = 1 .- p
        a = qtl.effect          # shorthand
        σ² = 2 .* p .* q .* a .* a
        qtl.pos[sortperm(vec(σ²))[end]]
    end
    for i in 1:2nid
        rand(Bernoulli(sr)) && (snp[ep, i] = 1)
    end
    ep                          # return the edited locus for testing
end
