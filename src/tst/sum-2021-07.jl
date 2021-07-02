using Serialization, DataFrames, JLD2, Statistics, Plots
using GEAS:QTL, top_QTL

function sum_2021_07(rpt)
    @load "dat/run/base.jld2" base
    qtl = deserialize("dat/2021-06-30/qtl-$rpt.ser")
    tq = top_QTL(base[:hap], qtl[2], n=10) # top QTL loci
    prd = deserialize("dat/2021-06-30/prd-c-$rpt.ser")
    snp = deserialize("dat/2021-06-30/snp-c-$rpt.ser")
    st = snp[tq, :]             # only the top QTL
    frq = zeros(11, 10)
    i = start = 1
    for stop in 4000:4000:44000
        frq[i, :] = mean(st[:, start:stop], dims=2)'
        start = stop + 1
        i += 1
    end
    plot(frq, leg=false)
end

function improvement()
    m2 = zeros(11, 3)
    tt = ['a', 'b', 'c']
    
    for i in 1:3
        for j in 1:10
            w = lpad(string(j), 2, '0')
            prd = deserialize("dat/2021-06-30/prd-$(tt[i])-$w.ser")
            bt = combine(groupby(prd, :g8n), :t2c => mean => :m)
            m2[:, i] += bt.m
        end
    end
    m2 ./= 10
    plot(m2, dpi=300)
end
