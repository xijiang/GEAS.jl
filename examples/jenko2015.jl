using Serialization, Statistics, DataFrames
import GEAS:MaCS, macs_2_base, macs_args, qtl_breeding, sim_QTL

function jenko_fork(fork, qtl, dic, nedit, esires, df, class)
    for t in nedit
        cc = "$class, $t edit"
        @info cc
        dic[:nEdit] = t         # how many (top) QTN to edit
        arg = (; dic...)
        base = deserialize(fork)
        ped, _ = qtl_breeding(base, arg, qtl, eSires = esires)
        tdf = combine(ped, :tbv => mean => :mtbv)
        tdf.sel .= cc
        append!(df, tdf)
    end
end

"""
    function jenko2015(; nid = 1000, nchr = 10, nsnp = 1000, maf = .05, dir = "dat/jenko2015")
---
## Description
This is to repeat Jenko 2015.
"""
function jenko2015(; nid = 1000, nchr = 10, nsnp = 1000, maf = .05, dir = "dat/jenko2015")
    @info "Simulation of a base with `MaCS`"
    par = macs_args(nid)
    tmp = joinpath(dir, "macs")
    isdir(tmp) || mkpath(tmp)
    
    # MaCS(par.cattle, nchr, dir = tmp)
    # base = macs_2_base(10, nsnp, dir = tmp)
    # serialize("$dir/base.ser", base)
    base = deserialize("$dir/base.ser")
    qtl = sim_QTL(base, nchr * nsnp)[1]

    dic = Dict(
        :nSire=> 25,
        :nDam => 500,
        :nSib => 2,             # such that selection is only on sires
        :nG8n => 21,
        :h²   => .5)            # dummy, only to utilize GEAS
    
    @info "Genomic selection on generation -20 → 0"
    
    arg = (; dic...)
    # Note, orginal work was just using TBV. No evaluation method envolved.
    ped, snp = qtl_breeding(base, arg, qtl)
    df = combine(ped, :tbv => mean => :mtbv)
    df.g8n[:] = arg.nG8n-1 : -1 : 0
    df.sel .= "pre-fork"        # 0 for generation -20 → 0
    
    # create new base, and insert editing option
    base.hap[:, :] = snp    # the new base, only to replace genotypes
    serialize("$dir/fork.ser", base) # this base will go through many different scenarios

    @info "Generation 1 → 20"
    dic[:nG8n] = 20
    
    @info "Continue selection with no editing"
    arg = (; dic...)
    ped, _ = qtl_breeding(base, arg, qtl)
    tdf = combine(ped, :tbv => mean => :mtbv)
    tdf.sel .= "GS only"     # 1 for continued normal selection from 0
    append!(df, tdf)
    
    # Editing the top 5 bulls, 25, 50, 100 QTN edits
    nedit = [25, 50, 100]
    esires = [ones(Int, 5); zeros(Int, 20)]
    jenko_fork("$dir/fork.ser", qtl, dic, nedit, esires, df, "top 5 sires")
    
    # Editing the top 10 bulls, 1, 5, 10, 20 QTN edits
    nedit = [1, 5, 10, 20]
    esires = [ones(Int, 10); zeros(Int, 15)]
    jenko_fork("$dir/fork.ser", qtl, dic, nedit, esires, df, "top 10 sires")
    
    # Editing the bottom 10 bulls, 1, 5, 10, 20 QTN edits
    esires = [zeros(Int, 15); ones(Int, 10)]
    jenko_fork("$dir/fork.ser", qtl, dic, nedit, esires, df, "bottom 10 sires")
    
    @info "Editing the all 25 bulls, 1, 5, 10 20 edits"
    esires = ones(Int, 25)
    jenko_fork("$dir/fork.ser", qtl, dic, nedit, esires, df, "all sires")

    df
end

function mytest(dic)
    dic[:tmp] => [1, 2, 3]
end
