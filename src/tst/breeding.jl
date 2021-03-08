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
function test_2_breeding(base, qtl, nsib)
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
                      :e19e => .8  # edit_successsful_rate
                      )
    par = (; Parameters...)     # named tuple.  contents as above
    prd = breeding_program(base, par, qtl)
    serialize(joinpath(dat_dir, "run/$nsib.ser"))
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
