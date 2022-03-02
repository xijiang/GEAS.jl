"""
    function workflow()
---
This records the framework of the package.
It can also serve as a template pipeline.

## Parameter example
- The keys must be of these names

```julia
    Parameters = Dict(:nSire=> 100,
                      :nDam => 200, # male:female = 1:2
                      :nSib => 2nsib,
                      :nC7e => nsib, # number of sib to be challenged
                      :nG8n => 10, # number of generations
                      :nQTL => [1000, 1000],
                      :h²   => [.5, .5],
                      :p8e  => .5, # percentage of affected in the binary trait
                      :e19e => .8, # edit_successsful_rate
                      :w4t  => weight # weight on the binary trait EBV
                      )
    par = (; Parameters...)     # named tuple.  contents as above
```
"""
function workflow(; debug = true)
    #copyright()
    if_normalization()
end

"""
    function data_collection()
This function collects some real data for downstream simulation.  It
1. remove duplicates
2. phase with `beagle`.
3. serialize the data to disk
4. also save a `jld2` file
"""
function data_collection()
    @time gt600()           # -> c.vcf.gz
    @time base = vcf2dic(joinpath(dat_dir, "run/c.vcf.gz"))
    @time serialize(joinpath(dat_dir, "run/ns.ser"), base)
    @save joinpath(dat_dir, "run/base.jld2") base
end

"""
    function a_simple_snp_blup_test()
This function using a real dataset from Nofima, which has 600 ID and 50,793 loci.
It consists 4 scenarios:
1. SNP-BLUP with 012 genotypes.
2. SNP-BLUP with normalized genotypes.
"""
function if_normalization()
    @load "dat/run/base.jld2" base
    nb = begin                  # a smaller new base
        loci = sort(randperm(size(base.hap)[1])[1:400])
        (hap = base.hap[loci, :], pos = base.pos[loci])
    end
    p = mean(nb.hap, dims=2)
    qtl = sim_QTL(nb, 350)[1]       # one QTL
    gt = hap2gt(nb.hap)             # → nlc × nid
    Q = qtl_gt(nb.hap, qtl.pos)     # QTL genotypes
    tbv = Q'qtl.effect
    z = zeros(size(gt))
    copyto!(z, gt)

    # snp-blup with no normalization
    lhs, rhs = begin
        x = [ones(600) z]
        x'x, x'tbv
    end
    esnp = lhs \ rhs
    ebv = gt * esnp[2:end]
    cor(ebv, tbv), norm(ebv - tbv)
    # returns usually close to (1, 0), meaning not necessary to normalize genotypes.
end
