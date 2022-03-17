#!/usr/bin/env julia
#using Plots, StatsPlots, Distributions, LaTeXStrings, Turing

"""
    function qsdgwas(nid, nlc, μ, σ, ϵ; maf = .1)
`Q`uick (and dirty) `s`imulation `d`ata for GWAS.
This function uses a uniform distribution of allele frequencies.
It is closer to reality than using just one unique frequency.

Here I used some smart, may be stupid later, one line code to generate genotypes.
- Ca = `rand.(Binomial.(2, p)` will generate genotypes of one individual.
- Cb = `[Ca for i in 1:nid]` then generate genotypes for all ID.
- `reduce(hcat, Cb)` will convert the results to a matrix of dimension `nid × nlc`.
"""
function qsdgwas(nid, nlc, μ, σ, ϵ; maf = .2)
    p = rand(Uniform(maf, 1-maf), nlc)              # allele frequencies,
    training_snp   = reduce(hcat, [rand.(Binomial.(2, p)) for i in 1:nid])'
    validation_snp = reduce(hcat, [rand.(Binomial.(2, p)) for i in 1:nid])'
    T = [ones(nid) training_snp]
    V = [ones(nid) validation_snp]
    b = rand(Normal(μ, σ), nlc + 1)                 # true SNP effects
    y = T * b + rand(Normal(0., ϵ), nid)            # phenotypes
    T, y, V, b, p
end

function naive_mcmc(T, y, σₑ; nit = 10_000)
    np = size(T, 2)             # number of parameters to estimate
    b = zeros(np)
    mcmcSample = zeros(nit, np)
    sp, tt = zeros(np), zeros(np)
    for i in 1:np
        tt[i] = T[:, i]'T[:, i]
        sp[i] = σₑ / sqrt(tt[i])
    end
    for i in 1:nit
        for t in 1:np
            b[t] = .0           # nullify this parameter
            yadj = y - T * b    # adjust y with every one of other parameters
            tcol = view(T, :, t)
            bhat = tcol'yadj / tt[t] # ̂b conditional on other parameters
            b[t] = rand(Normal(bhat, sp[t]))
        end
        mcmcSample[i, :] = b
    end
    mcmcSample
end

#=
function using_turing(T, y; nit = 2_000, nthreads = 8)
    # must using Turing to be valid
    @model linreg(X, y; predictors=size(X, 2)) = begin
        α ~ Normal(mean(y), 2.5std(y))
        β ~ filldist(TDist(3), predictors)
        σ ~ Exponential(1)
        y ~ MvNormal(α .+ X * β, σ)
    end
    model = linreg(T, y)
    chain = sample(model, NUTS(), MCMCThreads(), nit, nthreads)
    rst = summarystats(chain)
    rst[:, :mean][1:end-1]
end
=#

function _tst_gwas()
    nid, nlc, μ, σ, ϵ = 1000, 5, 0., .5, 1.0
    T, y, V, b, p = qsdgwas(nid, nlc, μ, σ, ϵ)
    sample = naive_mcmc(T, y, σ^2)

    # show the estimates side by side
    [b mean(sample, dims=1)' rst[:, :mean][1:1+nlc]]
end

"""
    function qsgp(nlc, nid, μ, σ, ϵ; maf = .2)
`Q`uick (and dirty) `s`imulation of (012) `g`enotypes and `p`henotypes.
The results are to evaluate SNP effects, so genotypes are locus majored.
That is, it is of `n_Loci × n_ID`.
"""
function qsgp(nlc, nid, μ, σ, ϵ; maf = .2)
    ps = rand(Uniform(maf, 1-maf), nlc) # allele frequencies
    T = reduce(hcat, [rand(Binomial(2, p), nid) for p in ps])'
    b = rand(Normal(μ, σ), nlc)
    y = T'b + rand(Normal(0., ϵ), nid)
    Int8.(T), y, b
end
