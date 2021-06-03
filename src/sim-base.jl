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
    isdir(sim) || mkpath(sim)
    print("Simulate chromosome:")
    Threads.@threads for ic in 1:c
        print(" $ic")
        run(pipeline(cmd,
                     stderr = devnull,
                     stdout = "$sim/$ic.txt"))
    end
    println()
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
    print("Collect chromosome:")
    for ic in 1:c
        print(" $ic")
        # chromsome genotype (cg), physical positions (cp), frequences (cf)
        cg, cp, cf = read_macs("$sim/$ic.txt") # c: current chromosome
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
        pos = [chr[1:totl] bp[1:totl]]
        r   = haldane(pos)
        Dict(:pos => pos,
             :r => r,
             :hap => Bool.(gt[1:totl, :]))
    end
end


"""
    function qmsim()
--
Simulate a base with `QMSim`.  Note, this program has many functions.
We only simulate SNP of a base population here.
Only the Linux version is used here.
Cross-platform running will be considered later.
"""
function qmsim_base(par; dir = ".")
    @info "Simulate a base with QMSim"
    isdir(dir) || mkpath(dir)
    tmp = mktempdir(dir)
    prm = [
        # Global parameters
        "title = \"A base with QMSim\";",
        "nrep = 1;",
        "h2 = 0;",
        "qtlh2 = 0;",
        "phvar = 1.0;",
        # Historical population
        # - I use a string here, two rows, the 2nd row specify 2Ne generations
        "begin_hp;",
        "    hg_size = $(par.hg_size);",
        "end_hp;",
        # Population
        "begin_pop = \"base\";",
        "    begin_founder;",
        "        male   [n = $(par.male), pop = \"hp\"];",
        "        female [n = $(par.female), pop = \"hp\"];",
        "    end_founder;",
        "    ls = 1;",          # litter size, to be discussed
        "    ng  = 1;",
        "    begin_popoutput;",
        "        genotype /gen 1;",
        "    end_popoutput;",
        "end_pop;",
        # Genome
        "begin_genome;",
        "    begin_chr = $(par.nchr);",
        "        chrlen = $(par.lchr);",
        "        nmloci = $(par.nsnp);",
        "        mpos   = rnd;",
        "        nma    = all 2;",
        "        maf    = eql;",
        "        nqloci = 0;",
        "        qpos   = rnd;",
        "        nqa    = all 2;",
        "        qaf    = eql;",
        "        qae    = rndg 0.4;",
        "    end_chr;",
        "end_genome;",
        # Output options
        "begin_output;",
        "    linkage_map;",
        "    output_folder = \"$tmp\";",
        "end_output;"
    ]
    open("$dir/qmsim.prm", "w") do io
        println(io, join(prm, "\n"))
    end
    run(pipeline(`$qmsim $dir/qmsim.prm -o`, devnull))
    read_qmsim_gt(tmp)
end

"""
    function read_qmsim_gt(dir)
---
Read `QMSim` genotype results into a `base`.
"""
function read_qmsim_gt(dir)
    chr = Int[]
    mgn = Float64[]             # Morgan positions
    open("$dir/lm_mrk_001.txt", "r") do io
        _ = readline(io)        # skip head line
        for line in eachline(io)
            c, m = split(line)[2:3]
            push!(chr, parse(Int, c))
            push!(mgn, parse(Float64, m))
        end
    end
    bp = Int.(floor.(mgn .* 1e6)) # suppose 1cM = 1e6 bp
    pos = [chr bp]
    r = haldane(pos)

    gt = begin
        t = Int[]
        open("$dir/base_mrk_001.txt", "r") do io
            _ = readline(io)
            for line in eachline(io)
                g = split(line)[2:end]
                append!(t, parse.(Int, g))
            end
        end
        t .-= 1                 # collected vector
        nlc = length(chr)
        nhp = length(t) ÷ nlc
        u = reshape(t, 2nlc, :)
        v = BitArray(undef, nlc, nhp)
        for i in 1:nlc
            h = 2i
            for j in 2:2:nhp
                id = j ÷ 2
                v[i, j-1:j] = u[h-1:h, id]
            end
        end
        v
    end
    
    Dict(:pos => pos,
         :r => r,
         :hap => gt
         )
end
