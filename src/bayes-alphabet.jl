"""
    function rr_blup(X, gt, y, par)
## Description
Estimate breeding values and locus effects with method RR-BLUP.
Where,
- `X`: is `n × p` incidence matrix for the non-genetic fixed effects.
- `gt`: are 0, 1, 2 genotypes.
- `y`: is a vector of phenotypes.
Argument `par` is a named tuple of parameters for the methods.
- `nit`: number of iterations, or chain length.
- `dflv`: hyperparameter (degree of freedom) for locus effect variance.
- `dfrv`: hyperparameter (degree of freedom) for residual variance.
- `vg`: used to derive hyper parameter (scale) for locus effect variance.
- `vr`: used to derive hyper parameter (scale) for residual variance.
- `window`: number of consecutive markers in a genomic window.

## Notes
This function is based on chapter 10, Bayesian methods applied to GWAS 
by Fernando and Garrick of Gondro, van der Werf, and Heyes, Genome-wide 
assocation studies and genomic prediction, 2013, Humana Press.
"""
function rr_blup(X, gt, y, par)
    num = (loci    = size(gt)[1],
           ID      = size(gt)[2],
           windows = Int(ceil(size(gt)[1] / par.window)),
           samples = Int(floor(par.n/par.ofrq))
           )
    p = mean(gt, dims = 2) ./ 2 # realized allele frequency
    Z = zeros(Float64, size(gt))
    copyto!(Z, gt)
    Z .-= 2p
    
    m2pq = mean(2 .* p .* (1 .- p))
    vq = par.var.genotype / num.loci / m2pq  # QTL variance
    # Scale parameters for the scaled inverse Χ² priors
    S² = (genotype = vq * (par.df.genotype - 2) / par.df.genotype,
          residue  = par.var.residue * (par.df.residue - 2) / par.df.residue)

    α, mα, μ, mμ = zeros(num.loci), zeros(num.loci), mean(y), 0
    mdfrq = zeros(num.loci)
    vg = zeros(num.samples)
    wvp = zeros(num.samples, num.windows) # window var prop.
    sc = 0                      # sample counts
    yc = y .- μ
    dχ² = Chisq(num.ID + par.df.residue)

    # the MCMC part
    for i in 1:par.n
        # sample residue variance
        ve = (yc'yc + par.df.residue * S².residue) / rand(dχ²)
        # sample intercept
        yc .+= μ
        rhs = sum(yc)
        
    end
end

function _test_ba(dat)
    T, y, b = begin
        if isfile(dat)
            @debug "Load formerly generated data"
            deserialize(dat)
        else
            @debug "Create data and share them with R"
            nlc, nid, μ, σ, ϵ = 5, 10, 0., .5, 1.0
            T, y, b = qsgp(nlc, nid, μ, σ, ϵ)
            dir = dirname(dat)
            open("$dir/T.csv", "w") do io
                for i in 1:nid
                    println(io, join(T[:, i], ", "))
                end
            end
            open("$dir/y.csv", "w") do io
                for x in y
                    println(io, x)
                end
            end
            open("$dir/b.csv", "w") do io
                for x in b
                    println(io, b)
                end
            end
            serialize(dat, (T = T, y = y, b = b))
            T, y, b
        end
    end

    par = (; Dict(:n => 10_000,
                  :df => (genotype = 4., residue = 4.),
                  :var => (genotype = 1., residue = 1.),
                  :window => 40, # windows size, n markers
                  :ofrq => 10     # output frequency
                  )...)
    nlc, nid = size(T)
    rr_blup(ones(nid), T, y, par)
end

function bayes_a()
end

function bayes_b(; nit = 10_000)
end

function bayes_c(; nit = 10_000)
end

function bayes_cπ(; nit = 10_000)
end
