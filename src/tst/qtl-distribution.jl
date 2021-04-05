"""
The codes in this file test the correlation of the simulated 2 traits.
It seems that there is a constant positive correlation of the 2 traits.
This made it unfair for the selection of the 2 traits.
As this is like to select one, the other trait is automatically selected.
"""
cc = zeros(1000)

function qc(n)
    nqtl = [n, n]
    qtl = sim_QTL(base, nqtl...)
    qa = qtl_gt(base[:hap], qtl[1].pos)
    qb = qtl_gt(base[:hap], qtl[2].pos)
    aa = qa'qtl[1].effect
    ab = qb'qtl[2].effect
    cc = cor(aa, ab)
    println(cc)
    cc
end

for n in 100:100:1000
    for i in 1:1000
        c = qc(n)
        cc[i] = c
    end
    histogram(cc, bins=19, label="$n QTL")
    savefig("$n.pdf")
end

#=
for these 10 pdf figures:

```bash
sudo dnf install texlive-pdfjam # if necessary
pdfjam {100..1000..100}.pdf --nup 5x2 --landscape -o p.pdf
evince p.pdf
```
=#