"""
    function read_macs(file)
---
Read genotypes and physical positions from simulation results of `macs`.
Returns
- genotypes of Array{Int8, 2}
- physical positions and
- allele frequencies
"""
function read_macs(file)
    gt = Int8[]
    ps = Float64[]
    open(file, "r") do io
        for line in eachline(io)
            (line[1:4] ≠ "SITE") && continue
            _, _, p, _, g = split(line)
            push!(ps, parse(Float64, p))
            append!(gt, parse.(Int8, collect(g)))
        end
    end
    nlc = length(ps)
    nhp = length(gt) ÷ nlc
    gt = reshape(gt, nhp, :)'
    fq = mean(gt, dims=2)
    return gt, ps, fq
end

"""
    function sample_macs(frq, maf, n)
---
Sample `n` loci of `maf` according `frq`. 
A warning message is given if not enough loci.
"""
function sample_macs(frq, maf, n)
    nlc = length(frq)
    df = DataFrame(ix = 1:nlc, freq = vec(frq))
    ix = filter(row -> maf < row.freq < 1 - maf, df).ix
    (length(ix) < n) && return ix
    return sort(shuffle(ix)[1:n])
end

"""
    function sim_base(h, c, l, s; t = 1e-5, r = 1e-5, maf=0.05)
---
Simulate a base population with [`macs`](https://github.com/gchen98/macs).
This function will simulate a population of 
- `h` haplotypes, or `h/2` individuals
- `c` chromosomes, each individual
- each chromosome of `l` base pairs, 
- each chromosome of `s` number of SNP
- loci of maf<0.05, by default, are excluded

`s` = min(s, number of simulated SNP fitered with `maf`,
which will be sampled from the simulated results.

## Example
```julia
GEAS.sim_base(1200, 29, 1e8, 2000)
```
will simulate 600 ID.  Each ID has 29 chromosome.
Each chromsome is of 1e8 base pairs.
Returns 
"""
function sim_base(h, c, l, s; t = 1e-5, r = 1e-5, maf=0.05)
    @info "Simulate a base population with `macs`"
    cmd = `$macs $h $l -t $t -r $r`
    sim = joinpath(dat_dir, "run/sim")
    hps = joinpath(sim, "haps.txt")
    isdir(sim) || mkpath(sim)
    nlc = c * s
    # base is not big, BitArray can be slow.
    # the real base is also of Int8. so
    gt  = Array{Int8, 2}(undef, nlc, h)
    # other simulation results
    chr = Array{Int, 1}(undef, nlc)
    bp  = Array{Int, 1}(undef, nlc)
    r   = Array{Float64, 1}(undef, nlc)
    from = 1
    totl = 0                    # in case not enough SNP simulated
    print("Simulating chromosome: ")
    for ic in 1:c
        print(" $ic")
        run(pipeline(cmd,
                     stderr = joinpath(sim, "debug.log"),
                     stdout = hps))
        # chromsome genotype, physical positions, frequences
        cg, cp, cf = read_macs(hps) # c: current chromosome
        ix = sample_macs(cf, maf, s)
        cs = length(ix)         # check if enough loci sampled
        cs < s && @warn "Less than expected SNP sampled on chromosome $ic"
        totl += cs
        gt[from:totl, :] = cg[ix, :] # add GT block
        chr[from:totl] = repeat([ic], cs)
        bp[from:totl] = Int.(floor.(cp[ix] .* l))
        from += cs
    end
    println()
    begin                       # create the dictionary and return
        pos = [chr bp]
        r   = haldane(pos)
        Dict(:pos => pos,
             :r => r,
             :hap => Bool.(gt))
    end
end
