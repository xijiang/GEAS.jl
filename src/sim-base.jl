"""
    function qsdgp(nid, nlc, μ, σ, ϵ; maf = .2)
`Q`uick (and `d`irty) `s`imulation of (012) `g`enotypes and `p`henotypes.
The results are to evaluate SNP effects, so genotypes are ID majored,
as there will be a `z'z` for `lhs`.
This function also simulate an intercept for individuals.
Allele frequencies were firstly sampled from `U(maf, 1 - maf)`.
Genotypes `T` were then sampled from `Binomial(2, pᵢ)`.
In the end, it returns a `(n_Loci + 1) × n_ID X` matrix, 
an observation vector `y = X'b + e` (`e ~ MVN(0, Iϵ)`), 
and true values `b` from `MVN(0, Iσ)`.
For simplicity, they are all of `Float64`.
The number of QTL are also affected by the QTL probability `pr` of simulated SNP.
"""
function qsdgp(nid, nlc, μ, σ, ϵ; maf = .2, pr = 1.)
    # uniform distributed allele freq, can be U-shaped
    X = ones(Int8, nlc+1, nid)
    for i in 2:nlc+1
        p = rand(Uniform(maf, 1 - maf)) # p~Uniform(maf, 1-maf)
        rand!(Binomial(2, p), view(X, i, :))
    end

    b = rand(Normal(μ, σ), nlc + 1)
    q = rand(Binomial(1, pr), nlc) # which ones are QTL
    @info "$(sum(q)) out of $nlc simulated SNP are QTL"
    b .*= [1; q]
    tbv = zeros(nid)
    matmul!(tbv, X', b)         # true breeding values
    vg = var(tbv)
    y = tbv + rand(Normal(0., ϵ), nid)
    vp = var(y)
    X, y, b, vg, vp
end

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
    function MaCS(args, nchr; dir = ".")
---
## Description
This function simulate `nchr` chromosomes with `MaCS`.

## Example
```julia
par = GEAS.macs_args(100)
GEAS.MaCS(par.cattle, 29, dir = "where-you-want-to-store-the-results"
```
Above will simulate `29` chromosomes. The results are `dir/chr.{1..nchr}`, and
`dir/info.{1..chr}`.

## Note
- Each chromosome might take about 1G disk space.
- The function runs the chromosome simulation in parallel.
- The number of threads is defined by `JULIA_NUM_THREADS` of your bash environment.
"""
function MaCS(args, nchr; dir = ".")
    isdir(dir) || mkpath(dir)
    seed = rand(Int, nchr)

    @info "Simulate a base population with `macs`" 
    Threads.@threads for chr in 1:nchr
        cmd = `$macs $(split(args)) -s $(seed[chr])`
        run(pipeline(cmd,
                     stderr = joinpath(dir, "info.$chr"),
                     stdout = joinpath(dir, "chr.$chr"))
            )
    end
end

"""
    function macs_2_base(nchr, nsnp; maf = .05, dir = ".", bppm = 1e8)
---
## Description
Sample SNP of `maf` from simulated `nchr` MaCS results in `dir`.
Each chromosome has `nsnp` SNP.
The extracted results are written in `base` format required by `GEAS`.

## Arguments
- `nchr`: Number of chromosome to be read.
  - Chromsomes were created from function `MaCS`
  - Error if `dir/chr.{1..nchr}` not exists.
- `nsnp`: Number of SNP to be sampled on each chromosome
  - An error is thrown When there aren't enough SNP
- `bppm`: base pairs per morgan.
"""
function macs_2_base(nchr, nsnp; maf = .05, dir = ".", bppm = 1e8)
    nhp = begin                 # Determine number of ID simulated
        file = joinpath(dir, "chr.1")
        isfile(file) || error("File $dir/chr.1 not exists")
        tmp = 0
        for line in eachline(file)
            if line[1:4] == "SITE"
                tmp = length(split(line)[5])
                break
            end
        end
        tmp
    end
    nlc = nchr * nsnp
    
    # Create storage for base
    base = (; Dict(
        :hap => BitArray(undef, nlc, nhp),
        :pos => zeros(Int, nlc, 2),
        :r => zeros(nlc))...)

    print("Dealing with chromsome:")
    for ic in 1:nchr            # Input the genotypes
        print(" $ic")
        file = joinpath(dir, "chr.$ic")
        isfile(file) || error("File $file not exist")
        gt, ps, fq = read_macs(file)
        ix = sample_macs(fq, maf, nsnp)
        length(ix) < nsnp && error("Not enough SNP sampled")
        fra, til = (ic - 1) * nsnp + 1, ic * nsnp
        base.hap[fra:til, :] = gt[ix, shuffle(1:nhp)]
        base.pos[fra:til, 1] .= ic
        base.pos[fra:til, 2] = Int.(floor.(ps[ix] .* bppm))
    end
    base.r[:] = haldane(base.pos, M = bppm)
    println()
    base
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
