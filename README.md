# GEAS.jl

## Installation
```julia
]
add https://github.com/xijiang/GEAS.jl
```

This package currently runs best on a Linux system.
Running on Windows were not tested.
It is most probably that the package runs on a Windows Subsystem for Linux(WSL), though.

## System requirement
This package paralleled calculations wherever possible.
It is hence better to specify as many threads as possible to speedup the simulation.
For example with an AMD 3900x CPU, you can add a line

```bash
export NUM_OF_THREADS=12
```

into your `~/.bashrc`.
Numbers that are greater than available physical cores won't help the speed.

## Base population preparation
The package now only deals with SNP data.  Multiallelic genotypes may be considered later.
The base population can either from some real data-sets, or from simulation.
A `base` is a named tuple with two components.
The `hap` part is a `BitArray` of phased SNP genotypes.
Its dimension is `n_loci` by `n_ID x 2`.
The `pos` part indicates the SNP locations.
It is a two-column matrix.
The first column contains chromosome numbers.
The second column holds the base pair positions.

`GEAS` has a function that can pull, compile and run `macs` and `QMSim` to make a base.

## Simulation

Note, below codes can be written into a driver file.
You may want to add these two lines in the file:
```julia
#!/usr/bin/env julia
using JLD2, DataFrames, LinearAlgebra, Distributions, Serialization, StatsPlots
```

### Setup scenarios
A simulation is defined by a named parameter tuple, which can be derived from a dictionary.
Below is an example:
```julia
par = Dict(
    :nSire=> 100,
    :nDam => 200,   # male:female = 1:2
    :nSib => 40,
    :nC7e => 30,    # number of sib of `nSib` to be challenged
    :nG8n => 10,   # number of generations to select
    
    # Note: Using tuples are more efficient than to use vectors
    :h²   => (a = .5, b =.5),
    :w4t  => 1,     # weight on the binary trait EBV
    :cv   => .2,    # covariance of pleiotripic QTL for two traits
    
    :p17h => .5,    # percentage-of-death after challenge
    :e19e => 1.,    # edit_successsful_rate (e19e), here, 1.0 = 100%
    :nK3n => 10,    # number of known (k3n) QTL on binary trait
    
    :dd   => 0.05,  # A small value to be added to diagonals
    :edit => false,  # whether to edit (the biggest) known QTL
    :fix  => false,  # if to fit (the largest) known QTL as fixed effects
    :b4y  => false, # for debug only, remove the threshold of binary trait
    :u38y => false  # use-genetypes-of-current-generation-only(u38y)
)
(; par...)
```

### An example simulation
Three scenarios can be ran with the following function. 
They are:
1. Breeding with ordinary SNP-BLUP
2. Breeding with emphasis on known QTL
3. Breeding with editing on known QTL

All three scenarios are using the same base (randomly mated) + first generation.
The results are saved to a `DataFrame`.
The `DataFrame` is then written to a random file as a repeat.

```julia
function example_simple_simulation(base, qtl, g8n)
    rst = DataFrame()
    
    par = gpar(g8n = g8n)
    @info "prepare common base + generation 1 of $(par.nG8n)"
    # by default Lapace below.
    ped, snp = create_storage(base, par, qtl)
    serialize("cmbase.ser", (ped=ped, snp=snp, qtl=qtl))

    @info "Breeding with SNP-BLUP"
    breeding(ped, snp, base, qtl, par)
    df = summarize(ped.prd, snp.prd, qtl[2])
    df.method = repeat(["SNP-BLUP"], nrow(df))
    append!(rst, df)

    @info "Breeding with emphasis on known QTL"
    ped, snp, qtl = deserialize("cmbase.ser")
    par = gpar(emph = true, g8n = g8n)
    breeding(ped, snp, base, qtl, par)
    df = summarize(ped.prd, snp.prd, qtl[2])
    df.method = repeat(["Fixed"], nrow(df))
    append!(rst, df)
    
    @info "Breeding with editing on known QTL"
    ped, snp, qtl = deserialize("cmbase.ser")
    par = gpar(edit = true, g8n = g8n)
    breeding(ped, snp, base, qtl, par)
    df = summarize(ped.prd, snp.prd, qtl[2])
    df.method = repeat(["Edit"], nrow(df))
    append!(rst, df)

    oo = basename(tempname())
    @save "$oo.jld2" rst
end
```

The arguments `base`, `qtl` and `g8n` can be generated as below:
```julia
function one_repeat(nqtl, g8n, jld)
    ntd = Threads.nthreads()
    BLAS.set_num_threads(ntd)   # use all specified threads (before REPL)
    @load jld base  # load population, which can be real, or simulated
    qtl = sim_QTL(base, nqtl, nqtl) # can change this to sim_pt_QTL()
    example_simple_simulation(base, qtl, g8n)
end
```
In above, `base` was prepared either from real data or simulation as described above.
It is saved in a `JLD2` file.
QTL `qtl` here are from two independent Laplace distributions.
QTL effects were scaled, such that the mean and additive genetic variance are approximate 0 and 1,
respectively.

One can replace `qtl = sim_QTL(base, nqtl, nqtl)` with
`qtl = sim_QTL(base, nqtl, nqtl, d = some_distribution)` to create scenarios of 
QTL effect of different distributions.

Or, with the following function:
```julia
function pt_repeat(nqtl, g8n, jld)
    ntd = Threads.nthreads()
    BLAS.set_num_threads(ntd)
    @load jld base
    qtl = sim_pt_QTL(base, 100, MvNormal([0., 0.], [1 0; 0 1.]))
    example_simple_simulation(base, qtl, g8n)
end
```
for MvNormally distributed pleiotropic QTL.

Then one can simply call, e.g.,  `one_repeat(500, 10, "base.jld2")` to have one repeat of simulation.

On a mainframe computer, you can simultaneously run many instances of such `one_repeats`. 
For one repeat of above scenarios, which contains 3 different breeding schemes,
it takes about 150 miutes with 12 threads of an AMD 3900x to finish.
With enough computation resources, you may obtain hundreds of repeats in similar time.

## Citing
- Xijiang Yu, Binyam S. Dagnachew, and Theo HE Meuwissen, A simulator to evaluate gene editing assested slelection program. WCGALP 2022, Rotterdam.
