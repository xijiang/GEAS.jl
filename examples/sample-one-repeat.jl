# Below packages can be used for simulation, data storage and summary
using JLD2, DataFrames, LinearAlgebra, Distributions, Serialization, StatsPlots

# These functions from `GEAS` are to be used in this simulation
import GEAS:breeding, sim_QTL, summarize, create_storage, sim_pt_QTL

"""
    function gpar(; edit = false, emph = false, cv = 0., g8n = 3)
---
Generate parameters for a simulation scenario.
Note, once the named parameter tuple is set, its value can't be changed.
You can modify the arguments of this function to change elements of `par`.
The return is a named tuple.
That is you can quote a value like `par.nSire`.
To avoid too long tuple names, I used numbers as shortcuts, e.g.,
`generation` → `g8n`.

`cv` is for pleiotropic QTL distribution, which follow a distribution of
`MvNormal(zeros(2), [1 cc; cc 1])`.  `cv` is the covariance of the QTL
effects for traits.
"""
function gpar(; binary = false, edit = false, emph = false, cv = 0., g8n = 3)
    par = Dict(
        :nSire=> 100,
        :nDam => 200,   # male:female = 1:2
        :nSib => 40,
        :nC7e => 30,    # number of sib of `nSib` to be challenged
        :nG8n => g8n,   # number of generations to select
        
        # Note: Using tuples are more efficient than to use vectors
        :h²   => (a = .5, b =.5),
        :w4t  => 1,     # weight on the binary trait EBV
        :cv   => cv,    # covariance of pleiotripic QTL for two traits
        
        :p17h => .5,    # percentage-of-death after challenge
        :e19e => 1.,    # edit_successsful_rate (e19e), here, 1.0 = 100%
        :nK3n => 10,    # number of known (k3n) QTL on binary trait
        
        :dd   => 0.01,  # A small value to be added to diagonals
        :edit => edit,  # whether to edit (the biggest) known QTL
        :fix  => emph, # if to fit (the largest) known QTL as fixed effects
        :b4y  => binary # for debug only, remove the threshold of binary trait
    )
    (; par...)
end


"""
    function example_simple_simulation(base, qtl, g8n; dir = ".")
---
# Objectives
- Do a simple breeding simulation with a production and a challenge trait.
- Write intermediate and final results in `dir`, which, by default, is the `pwd`.

# Introduction
This module is designed to run on a large mainframe machine, 
where many repeats can be ran simultaneously.
The run is only one repeat of some certain scenario.

This function uses all the set threads for matrix calculation.
Run:
```bash
export JULIA_NUM_THREADS=12
```
before starting Julia.
Or, put above line into your `~/.bashrc`.
You'd better set this number ≤ your available physical CPU cores.
"""
function example_simple_simulation(base, qtl, g8n; dir = ".")
    rst = DataFrame()
    
    par = gpar(g8n = g8n)
    @info "prepare common base + generation 1 of $(par.nG8n)"
    # by default Lapace below.
    ped, snp = create_storage(base, par, qtl)
    serialize("$dir/cmbase.ser", (ped=ped, snp=snp, qtl=qtl))

    @info "Breeding with SNP-BLUP"
    breeding(ped, snp, base, qtl, par)
    df = summarize(ped.prd, snp.prd, qtl[2])
    df.method = repeat(["SNP-BLUP"], nrow(df))
    append!(rst, df)

    @info "Breeding with emphasis on known QTL"
    ped, snp, qtl = deserialize("$dir/cmbase.ser")
    par = gpar(emph = true, g8n = g8n)
    breeding(ped, snp, base, qtl, par)
    df = summarize(ped.prd, snp.prd, qtl[2])
    df.method = repeat(["Fixed"], nrow(df))
    append!(rst, df)
    
    @info "Breeding with edit on known QTL"
    ped, snp, qtl = deserialize("$dir/cmbase.ser")
    par = gpar(edit = true, g8n = g8n)
    breeding(ped, snp, base, qtl, par)
    df = summarize(ped.prd, snp.prd, qtl[2])
    df.method = repeat(["Edit"], nrow(df))
    append!(rst, df)

    oo = tempname(dir)
    @save "$oo-r.jld2" rst
end

"""
    function one_repeat(nqtl, g8n, jld; dir = ".")
---
Load a base population in `jld`, either real or simulated.
This function simulate two sets of QTL of Laplace distribution for two traits.
The `example_simple_simulation` will do breeding with `SNP-BLUP`, `Fixed`
and `Edit` methods.
Results (`*.jld2`) and intermediate data (`*ser`) are stored in `dir`.
"""
function one_repeat(nqtl, g8n, jld; dir = ".")
    isdir(dir) || mkpath(dir)
    ntd = Threads.nthreads()
    BLAS.set_num_threads(ntd)   # use all specified threads (before REPL)
    @load jld base  # load population, which can be real, or simulated
    qtl = sim_QTL(base, nqtl, nqtl) # can change this to sim_pt_QTL()
    example_simple_simulation(base, qtl, g8n, dir = dir)
end

"""
    function pt_repeat(nqtl, g8n, jld; dir = ".")
---
Load a base population in `jld`, either real or simulated.
This function simulate one set of QTL of `MvNormal([0, 0], [1 0; 0 1])`
distribution for two traits. QTL are pleiotropic.
The `example_simple_simulation` will do breeding with `SNP-BLUP`, `Fixed`
and `Edit` methods.
Results (`*.jld2`) and intermediate data (`*ser`) are stored in `dir`.
"""
function pt_repeat(nqtl, g8n, jld; dir = ".")
    isdir(dir) || mkpath(dir)
    ntd = Threads.nthreads()
    BLAS.set_num_threads(ntd)
    @load jld base
    qtl = sim_pt_QTL(base, 100, MvNormal([0., 0.], [1 0; 0 1.]))
    example_simple_simulation(base, qtl, g8n, dir = dir)
end


function collect_repeats(; dir = ".")
    df = DataFrame()
    for file in readdir(dir)
        if file[1:3] == "jl_"
            @load "$dir/$file" rst
            append!(df, rst)
        end
    end
    rst = DataFrame()
    for gp in groupby(df, :method)
        for gn in groupby(gp, :g8n)
            tmp = combine(gn,
                          :mpbv => mean => :mpbv,
                          :pvg  => mean => :pvg,
                          :mcbv => mean => :mcbv,
                          :cvg  => mean => :cvg,
                          :mica => mean => :mica)
            tmp.method .= gp.method[1]
            tmp.g8n .= gn.g8n[1]
            append!(rst, tmp)
        end
    end
    bymethod = groupby(rst, :method)
    p1 = plot(bymethod[1].g8n, bymethod[1].mpbv,
         xlabel = "Generation",
         ylabel = "Mean TBV, production",
         legend = :topleft,
         label = bymethod[1].method[1],
         xticks = 0:11);
    plot!(p1, bymethod[2].g8n, bymethod[2].mpbv, label = bymethod[2].method[1]);
    plot!(p1, bymethod[3].g8n, bymethod[3].mpbv, label = bymethod[3].method[1]);

    p2 = plot(bymethod[1].g8n, bymethod[1].mcbv,
         xlabel = "Generation",
         ylabel = "Mean TBV, challenge (cont.)",
         legend = :topleft,
         label = bymethod[1].method[1],
         xticks = 0:11);
    plot!(p2, bymethod[2].g8n, bymethod[2].mcbv, label = bymethod[2].method[1]);
    plot!(p2, bymethod[3].g8n, bymethod[3].mcbv, label = bymethod[3].method[1]);

    p3 = plot(bymethod[1].g8n, bymethod[1].pvg,
         xlabel = "Generation",
         ylabel = "Var(TBV), production",
         legend = :topright,
         label = bymethod[1].method[1],
         xticks = 0:11);
    plot!(p3, bymethod[2].g8n, bymethod[2].pvg, label = bymethod[2].method[1]);
    plot!(p3, bymethod[3].g8n, bymethod[3].pvg, label = bymethod[3].method[1]);

    p4 = plot(bymethod[1].g8n, bymethod[1].cvg,
         xlabel = "Generation",
         ylabel = "Var(TBV), challenge (cont.)",
         legend = :topright,
         label = bymethod[1].method[1],
         xticks = 0:11);
    plot!(p4, bymethod[2].g8n, bymethod[2].cvg, label = bymethod[2].method[1]);
    plot!(p4, bymethod[3].g8n, bymethod[3].cvg, label = bymethod[3].method[1]);

    p5 = plot(bymethod[1].g8n, bymethod[1].mica,
              xlabel = "Generation",
              ylabel = "Mean inbreeding, production",
              legend = :topleft,
              label = bymethod[1].method[1],
              xticks = 0:11);
    plot!(p5, bymethod[2].g8n, bymethod[2].mica, label = bymethod[2].method[1]);
    plot!(p5, bymethod[3].g8n, bymethod[3].mica, label = bymethod[3].method[1]);

    plot(p1, p3, p5, p2, p4);
    savefig("$dir/rst.pdf")
end
