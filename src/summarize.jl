"""
    function compact_ped(ped)
---
This function is to simplify inbreeding calculation.
With *A* matrix, full-sibs have the same inbreeding value.
This function returns a dataframe, taking the first ID of a full sibship.
There is also a column counting the sibship size.
"""
function compact_ped(ped)
    t = DataFrame(g8n=Int16[], id=Int32[], ic=Float64[])
    pa = ma = count = 0
    for x in eachrow(ped)
        if x.sir != pa && x.dam != ma
            pa, ma = x.sir, x.dam
            push!(t, (x.g8n, x.id, 0.))
        end
    end
    v = combine(groupby(ped, [:g8n, :sir, :dam]), nrow => :count)
    t.count = v.count[2:end]
    t
end


"""
    function summarize(ped, snp, qtl)
---
# Objective
Summarize the results from function `breeding`.
Return results in a DataFrame

The following are calculated
1. mean TBV of each trait
2. genetic variance of each trait
3. mean inbreeding coefficient of each generation based on A matrix

# Arguments
Here,
- `ped` is the pedigree of the production/1st trait, where selection is
also taking place.
- `snp` is the SNP genotype matrix of the production trait.
- `qtl` is the `QTL` vector for the binary (2nd) traits.

# Inbreeding using A matrix
- Recursive method is used, as here we have large full sibship, 
    only one calculation is needed.
- Inbreeding calculation is also multiple threaded

# Returns
- mean inbreeding coefficients by generation
- mean TBV by generation
"""
function summarize(ped, snp, qtl)
    bgt = alleles2gt(view(snp, qtl.pos, :))
    ped.clg = bgt'qtl.effect
    rst = combine(groupby(ped, :g8n),
                  :tbv => mean => :mpbv,
                  :tbv => var => :pvg,
                  :clg => mean => :mcbv,
                  :clg => var => :cvg)
    mpd = Matrix(select(ped, :sir, :dam)) # my pedigree, for func `kinship`
    
    # Inbreeding using A matrix, create a compact pedigree below
    cpd = compact_ped(ped)
    @Threads.threads for i in 1:nrow(cpd)
        cpd[i, :ic] = kinship(mpd, cpd[i, :id], cpd[i,:id]) - 1
    end
    # mean inbreeding coefficients with A matrix
    mica = combine(groupby(cpd, :g8n),
                   [:ic, :count] => ((x, y) -> x'y/sum(y)) => :mica).mica
    rst.mica = [0; mica]
    rst
end


"""
    function sumsum(df)
---
This `sumsum` is to summarize some summarized results.
The results, in a `DataFrame`, are grouped by repeats.
Then, the `combine` function is used.
"""
function sumsum(df)
    gp = groupby(df, :rpt)
    m = copy(gp[1])
    for j in 2:ncol(m)-1
        for i in 2:length(gp)
            m[:, j] += gp[i][:, j]
        end
        m[:, j] ./= last(df).rpt
    end
    m
end

"""
    function sum_qtl_frq(prd, snp, qtl)
---
QTL frequency changes in the production 
"""
function sum_qtl_frq(prd, snp, qtl)
    gg = groupby(prd, :g8n)     # group on generation
    fq = Float64[]
    for grp in gg
        fra = 2first(grp).id - 1
        til = 2last(grp).id
        qa = view(snp, qtl.pos, fra:til)
        append!(fq, mean(qa, dims=2))
    end
    reshape(fq, length(gg), :)
end

"""
    function allele_frq_flux(prd, snp, loci)
---
Return the frequency change of `loci` over generations.
"""
function allele_frq_flux(prd, snp, loci)
    gp = groupby(prd, :g8n)
    ng = length(gp)
    df = DataFrame()
    ft = []                     # [[from, to]]
    for df in gp
        fra, til = 2first(df.id) - 1, 2last(df.id)
        push!(ft, [fra, til])
    end
    for l in loci
        frq = Float64[]
        for (fra, til) in ft
            push!(frq, mean(view(snp, l, fra:til)))
        end
        df[!, string(l)] = frq
    end
    df
end

"""
    function bp_summary(prd, snp, qtl, dir)
---
Summarize simulation results of the `b`reeding `p`rogram.
The figures and tables will be saved in `dir`.
"""
function bp_summary(prd, snp, qtl; dir=".")
    isdir(dir) || mkpath(dir)
    # construct pedigree
    n₀ = length(
        unique(
            Matrix(
                select(
                    filter(row ->row.g8n ==1, prd), :sir, :dam))))
    ped = begin                 # this was the way I recorded the pedigree
        n = n₀ + nrow(prd)      # total ID number, including base
        ped = zeros(Int32, n, 2) # the first `n₀` are automatically 0s
        # parents of generation 1 refer to base
        pm₁ = Matrix(select(filter(row -> row.g8n == 1, prd), :sir, :dam))
        n₁ = size(pm₁)[1]
        ped[n₀+1:n₀+n₁, :] = pm₁
        # parents of generation 2+ are about descendants
        pm₂ = Matrix(select(filter(row -> row.g8n > 1, prd), :sir, :dam))
        pm₂ .+= n₀
        ped[n₀+n₁+1:end, :] = pm₂
        ped
    end

    @info "Inbreeding using A matrix"
    begin                 # generation average inbreeding coefficients
        # only need to record generation, and one ID number of its many full sibs.
        n = size(ped)[1]
        g8n = zeros(Int, n)
        g8n[n₀+1:end] = prd.g8n
        df = DataFrame(g8n=Int[], id=Int[], cnt=Int[], ic=Float64[])
        a = b = 0
        for (g, id, pa, ma) in eachrow([g8n 1:n ped])
            g == 0 && continue
            if pa == a && ma == b
                df.cnt[end] += 1
                continue
            end
            push!(df, [g, id, 1, 0])
            a, b = pa, ma
        end
        Threads.@threads for i in 1:nrow(df)
            df.ic[i] = kinship(ped, df.id[i], df.id[i]) - 1
        end
        boxplot(df.g8n, df.ic, label="Inbreeding measured by A")
        savefig("$dir/inbreeding-changes-with-A.pdf")
    end

    @info "Inbreeding using homozygous loci"
    begin
        m = size(snp)[1]
        df = DataFrame(g8n = prd.g8n, ic = zeros(nrow(prd)))
        Threads.@threads for k in 1:nrow(prd)
            df.ic[k] = m - sum(snp[:, 2k-1] .⊻ snp[:, 2k])
        end
        df.ic ./= m
        t = mean(filter(row -> row.g8n == 1, df).ic)
        df.ic ./= t
        df.ic .-= 1
        boxplot(df.g8n, df.ic, label="Inbreeding measured with homozygotes ratio")
        savefig("$dir/inbreeding-changes-with-homozygous-loci.pdf")
    end

    # Something not right below.  Will visit this later.
    #@info "Inbreeding using G matrix"
    #begin
    #    p = mean(snp, dims = 2)
    #    tp = 2p
    #    s2pq = sum((1 .- p) .* tp)
    #    r2pq = 1 / s2pq
    #    df = DataFrame(g8n = prd.g8n, ic = zeros(nrow(prd)))
    #    for i in 1:nrow(prd)
    #        z = (snp[:, 2i-1] + snp[:, 2i]) - tp
    #        df.ic[i] = (z'z)[1]*r2pq - 1
    #    end
    #    boxplot(df.g8n, df.ic, label="Inbreeding measured with G")
    #    savefig("doc/fig/inbreeding-changes-with-G.pdf")
    #end
    
    gp = groupby(prd, :g8n)
    @info "Genetic variance and inbreeding over generations"
    begin                                 # plot vₐ² trends of 2 traits
        p = combine(gp, :tbv => var => :vprd)
        b = combine(gp, :t2c => var => :vbin)
        plot(p.vprd, label="Production", dpi=300)
        plot!(b.vbin, label="Binary")
        savefig("$dir/Va-changes.pdf")
    end
    @info join(["",
                "QTL frequency changes",
                "  - if negtive a, switch p to 1-p",
                "  - bubble plot, sizes are with respect to initial frq"], '\n')
    begin
        qg = GEAS.alleles2gt(snp[qtl[2].pos, :])
        n1 = nrow(filter(row -> row.g8n == 1, prd))
        p1 = mean(qg[:, 1:n1], dims=2)./2
        n2 = nrow(filter(row -> row.g8n == prd.g8n[end], prd))
        p2 = mean(qg[:, end-n2:end], dims=2)./2
        scatter(p1, p2 .- p1, ms=p1.*10,
                label="Frequency changes",
                xlabel="MAF",
                ylabel="frequency change",
                dpi=300)
        savefig("$dir/frequency-changes.pdf")
    end
end

function sum_improve(prd)
    gp = groupby(prd, :g8n)
    @info "Genetic variance and inbreeding over generations"
    p = combine(gp, :tbv => mean => :mprd)
    b = combine(gp, :t2c => mean => :mbin)
    p, b
end
