"""
    function bp_summary(prd, snp, qtl)
---
Summarize simulation results of the `b`reeding `p`rogram.
The figures and tables will be saved with an initial 'bps-'.
"""
function bp_summary(prd, snp, qtl)
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
        savefig("doc/fig/inbreeding-changes-with-A.pdf")
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
        savefig("doc/fig/inbreeding-changes-with-homozygous-loci.pdf")
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
        savefig("doc/fig/Va-changes.pdf")
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
        savefig("doc/fig/frequency-changes.pdf")
    end
end
