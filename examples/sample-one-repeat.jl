# Below packages can be used for simulation, data storage and summary
using JLD2, DataFrames, LinearAlgebra, Distributions, Serialization, StatsPlots

# These functions from `GEAS` are to be used in this simulation
import GEAS:breeding, sim_QTL, summarize, create_storage, sim_pt_QTL
import GEAS:rank_QTL, calc_idx

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
function gpar(;
              binary = false,   # whether to convert challenge trait binary
              edit   = false,
              emph   = false,   # if to emphasize known QTL
              cv     = 0.,      # for pleiotropic trait
              g8n    = 3,       # number of generation to simulate
              u38y   = false)
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
        :b4y  => binary, # for debug only, remove the threshold of binary trait
        :u38y => u38y # "use genotypes of current generation only" :40char
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
    
    # par = gpar(g8n = g8n, u38y = true)
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

You may find this `example` file by:
```julia
dir = pathof(GEAS)
til = findlast("/src", dir)[1]
example = joinpath(tmp[1:til], "examples/sample-one-repeat.jl")
```
"""
function one_repeat(nqtl, g8n, jld; dir = ".")
    isdir(dir) || mkpath(dir)
    ntd = Threads.nthreads()
    BLAS.set_num_threads(ntd)   # use all specified threads (before REPL)
    # Note `twop` was added to base.jld2, Jan, 2022
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

function test_codes(u38y = false)
    ntd = Threads.nthreads()
    BLAS.set_num_threads(ntd)   # use all specified threads (before REPL)
    @load "dat/run/base.jld2" base
    qtl = sim_QTL(base, 100, 100) # can change this to sim_pt_QTL()
    par = gpar(g8n = 3, u38y = u38y)
    @info "prepare common base + generation 1 of $(par.nG8n)"
    # by default Lapace below.
    ped, snp = create_storage(base, par, qtl)

    rkq = rank_QTL(base.hap, qtl[2]) # QTL for binary trait ranking in ↓ order
    Q = par.fix ? abs.(rkq[1:par.nK3n]) : Int[]
    prd, clg = groupby(ped.prd, :g8n), groupby(ped.clg, :g8n) # shorthands

    # note, there are `nG8n+1` generations, including the base.
    @info "Evaluating generation 1"
    calc_idx(ped, snp, prd[2], Q, par) # evaluation of generation one
end
