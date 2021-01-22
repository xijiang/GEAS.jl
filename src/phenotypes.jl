"""
    function phenotype(gt, QTL, h²)
---
This function returns phenotypes a trait.
"""
function phenotype(gt, qtl, h²)
    @warn "under construction"
end

function TBV(snp, qtl)
    nid = size(snp)[2] ÷ 2
    pos = qtl.pos
    nql = length(pos)
    
    g = zeros(Int8, nql, nid)
    e = qtl.effect
    for id in 1:nid
        g[:, id] = snp[pos, id*2 - 1] + snp[pos, id*2]
    end

    bv = g'e
    tmp = DataFrame(bv = bv, id = 1:nid)
    sort!(tmp, :bv, rev=true)
    rank = tmp.id               # rank offspring in BV order, high->low
    @inbounds for s in 50:50:600
        x = g[:, rank[1:s]]     # QTL genotypes of the selected ones
        b = bv[rank[1:s]]       # BV of the selected ones
        q = sum(x, dims=2)      # frequency: 0 <--> 2s
        m = mean(b)
        v = var(b)
        t = qtl.max
        f = 0
        @inbounds for i in 1:length(e)
            if q[i] == 0
                if e[i] > 0
                    t -= e[i]
                end
                f += 1
            elseif q[i] == 2s
                if e[i] < 0
                    t += e[i]
                end
                f += 1
            end
        end
        println("$m, $v, $t, $f")
    end
end
