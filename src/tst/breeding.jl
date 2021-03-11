"""
    function test_2_breeding(base, qtl, nsib)
---
- 100 sires, 200 dams
- sibship size: 10, 20, 40 for both produciton and challenge group
- edit 20,000 fertilized eggs with a rate
  * also chech Hicky's editing method
- compare
  * no editing
  * editing
- order QTL on `2pqa²`.
- equal weights on both traits.
- 2021-March-06
"""
function test_2_breeding(base, qtl, nsib, weight)
    # It is better to setup simulation scenarios with a dictionary
    # Then pass it as a named tuple, so that we can conveniantly use `x.y`
    # using a dictionary also makes it easy to modify parameter structure.
    # the keys in this dictionry must exist to setup a simulation
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
    prd = breeding_program(base, par, qtl)
    serialize(joinpath(dat_dir, "run/breeding/w/$weight-$nsib.ser"), prd)
end

"""
    function sum_w_breeding()
---
Summarize genetic improvement using different weight on GEBV of binary trait.
The conclusion is `weight` should be some where between 2 and 3.
- 2021-Mar-11
"""
function sum_w_breeding()
    sib    = [10, 20, 40]
    weight = [3, 6, 9, 12]
    ave = Array{Float64, 4}(undef, 3, 4, 2, 10)
    for i in 1:3
        for j in 1:4
            rst = deserialize("dat/run/breeding/w/$(weight[j])-$(sib[i]).ser")
            k = 1
            for d in groupby(rst, :g8n)
                ave[i, j, 1, k] = mean(d.tbv)
                ave[i, j, 2, k] = mean(d.t2c)
                k += 1
            end
        end
    end
    for i in 1:3
        ave[i, 1, 1, :] .-= ave[i, 1, 1, 1]
        plot(ave[i, 1, 1, :], label="prd-$(sib[i])s-3w", dpi=300, legend=:topleft)
        for j in 2:4
            ave[i, j, 1, :] .-= ave[i, j, 1, 1]
            plot!(ave[i, j, 1, :], label="prd-$(sib[i])s-$(weight[j])w")
        end
        for j in 1:4
            ave[i, j, 2, :] .-= ave[i, j, 2, 1]
            plot!(ave[i, j, 2, :], label="tst-$(sib[i])s-$(weight[j])w", linestyle=:dash)
        end
        savefig(joinpath(dat_dir, "run/breeding/w/sib-$(sib[i]).pdf"))
    end
end

"""
    function sum_2_breeding()
---
Function `test_2_breeding` of v0.5.2 resulted in `breeding/{1,2,3}'.
This function is to summarize the results in 9 files.
Below was done Mar. 08, 2021
"""
function sum_2_breeding()
    for i in 1:3
        m = Float64[]
        for sib in [10, 20, 40]
            dat = deserialize(joinpath(dat_dir, "run/breeding/$i/$sib.ser"))
            for d in groupby(dat, :g8n)
                push!(m, mean(d.tbv))
                push!(m, mean(d.t2c))
            end
        end
        m = reshape(m, 6, :)
        t = reshape(m[1:20], 2, :)
        t .-= t[:, 1]           # set the first generation as 0.
        plot(t[1, :], label="Production-10-sibs", dpi=300, legend=:topleft)
        plot!(t[2, :], label="Challenge-10-sibs")
        t = reshape(m[21:40], 2, :)
        t .-= t[:, 1]           # set the first generation as 0.
        plot!(t[1, :], label="Production-20-sibs")
        plot!(t[2, :], label="Challenge-20-sibs")
        t = reshape(m[41:60], 2, :)
        t .-= t[:, 1]           # set the first generation as 0.
        plot!(t[1, :], label="Production-40-sibs")
        plot!(t[2, :], label="Challenge-40-sibs")
        savefig(joinpath(dat_dir, "run/breeding/fig/$i.pdf"))
    end
end

"""
    function test_breeding()
---
This is a independent function.  It tests selection on
TBV,
phenotypes,
and random variable (no selection).
The result is a pdf figure.
So far so good.
"""
function test_breeding()
    Parameters = Dict(:nSire=> 125,
                      :nDam => 250, # male:female = 1:2
                      :nSib => 100,
                      :nC7e => 40, # number of sib to be challenged
                      :nG8n => 40, # number of generations
                      :nQTL => [1000, 1000],
                      :h²   => [.2, .2],
                      :p8e  => .5, # percentage of affected in the binary trait
                      )
    par = (; Parameters...)     # named tuple.  contents as above
    
    # read back real data
    #base = deserialize(joinpath(dat_dir, "run/ns.ser")) # to save time
    @load "dat/run/base.jld2"
    
    # general parameters
    qtl = sim_QTL(base, par.nQTL...)
    @info "Create storage and generation one"
    snp₁, snp₂, prd, clg, f1, f2, f3, f4 = # f: from
        create_storage(base, par, qtl)
    # shorthands
    nSire, nDam, nSib, nChallenge, percentage, generation =
        par.nSire, par.nDam, par.nSib, par.nC7e, par.p8e, par.nG8n
    ip, ic, hp, hc = test_sets(nDam, nSib, nChallenge)
    l1, l2, l3, l4 = length.((ip, ic, hp, hc)) # lengths of generation 2+
    mp = zeros(generation)      # mean phenotype on generation
    average(ig) = begin
        tmp = select(prd, :g8n, :p7e)
        mean(tmp.p7e[tmp.g8n .== ig])
    end
    mp[1] = average(1)
    
    @info "!!!!! Select on TBV"
    a, b, c, d = f1, f2, f3, f4
    for ig in 2:generation
        @info "generation $ig"
        df = select(prd[(prd.g8n .== ig-1), :], :id, :sex, :tbv => :val)
        nuclear = select_nuclear(df, nSire, nDam)
        
        # simulate current generation
        ped = random_mate(nuclear, nSire, nDam, nSib)
        nxt = gdrop(snp₁, ped, base[:r])
        bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
        
        # update storage.  codes below are ugly.
        prd.tbv[f1+1:f1+l1] = bv₁[ip]
        prd.p7e[f1+1:f1+l1] =  p₁[ip]
        prd.sir[f1+1:f1+l1] = ped[ip, 1]
        prd.dam[f1+1:f1+l1] = ped[ip, 2]
        clg.tbv[f2+1:f2+l2] = bv₂[ic]
        clg.p7e[f2+1:f2+l2] =  p₂[ic]
        snp₁[:, f3+1:f3+l3] = nxt[:, hp]
        snp₂[:, f4+1:f4+l4] = nxt[:, hc]
        f1 += l1
        f2 += l2
        f3 += l3
        f4 += l4
        mp[ig] = average(ig)
    end
    plot(mp, label="TBV", legend=:topleft)
    
    @info "!!!!! Select on phenotypes"
    f1, f2, f3, f4 = a, b, c, d
    for ig in 2:generation
        @info "generation $ig"
        df = select(prd[(prd.g8n .== ig-1), :], :id, :sex, :p7e => :val)
        nuclear = select_nuclear(df, nSire, nDam)
        
        # simulate current generation
        ped = random_mate(nuclear, nSire, nDam, nSib)
        nxt = gdrop(snp₁, ped, base[:r])
        bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
        
        # update storage.  codes below are ugly.
        prd.tbv[f1+1:f1+l1] = bv₁[ip]
        prd.p7e[f1+1:f1+l1] =  p₁[ip]
        prd.sir[f1+1:f1+l1] = ped[ip, 1]
        prd.dam[f1+1:f1+l1] = ped[ip, 2]
        clg.tbv[f2+1:f2+l2] = bv₂[ic]
        clg.p7e[f2+1:f2+l2] =  p₂[ic]
        snp₁[:, f3+1:f3+l3] = nxt[:, hp]
        snp₂[:, f4+1:f4+l4] = nxt[:, hc]
        f1 += l1
        f2 += l2
        f3 += l3
        f4 += l4
        mp[ig] = average(ig)
    end
    plot!(mp, label="Phenotype")

    @info "!!!!! Random mating"
    f1, f2, f3, f4 = a, b, c, d
    for ig in 2:generation
        @info "generation $ig"
        # simulate current generation
        df = select(prd[(prd.g8n .== ig-1), :], :id, :sex)
        df.val = rand(length(df.id))
        nuclear = select_nuclear(df, nSire, nDam)
        ped = random_mate(nuclear, nSire, nDam, nSib)
        nxt = gdrop(snp₁, ped, base[:r])
        bv₁, p₁ = phenotype(nxt, qtl[1], par.h²[1])
        bv₂, p₂ = phenotype(nxt, qtl[2], par.h²[2], par.p8e)
        
        # update storage.  codes below are ugly.
        prd.tbv[f1+1:f1+l1] = bv₁[ip]
        prd.p7e[f1+1:f1+l1] =  p₁[ip]
        prd.sir[f1+1:f1+l1] = ped[ip, 1]
        prd.dam[f1+1:f1+l1] = ped[ip, 2]
        clg.tbv[f2+1:f2+l2] = bv₂[ic]
        clg.p7e[f2+1:f2+l2] =  p₂[ic]
        snp₁[:, f3+1:f3+l3] = nxt[:, hp]
        snp₂[:, f4+1:f4+l4] = nxt[:, hc]
        f1 += l1
        f2 += l2
        f3 += l3
        f4 += l4
        mp[ig] = average(ig)
    end
    plot!(mp, label="Random mating")
    savefig("doc/fig/test-breeding-v3.pdf")
end
