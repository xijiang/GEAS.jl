"""
    function test_2021_03_26()
---
To test:
- if the right loci is found
- if the editing is right
"""
function test_2021_03_26()
    @load "dat/run/base.jld2" base
    qtl = sim_QTL(base, 10)[1]
    println(qtl.pos)
    snp = copy(base[:hap])
    ep = gedit(base[:hap], qtl, .8)
    println(sum(snp[ep, :]), ' ', sum(base[:hap][ep, :]))
end

"""
    function edit_2021_03_26()
---
1. simulate a base with `macs`, of 29 chromosomes
2. simulate a QTL setup
3. compare
   - random mating
   - selection with GS
   - selection with GEAS
"""
function edit_2021_03_26()
    Parameters = Dict(:nSire=> 125,
                      :nDam => 250, # male:female = 1:2
                      :nSib => 100,
                      :nC7e => 40, # number of sib to be challenged
                      :nG8n => 5,  # number of generations
                      :nQTL => [1000, 1000],
                      :h²   => [.5, .5],
                      :p8e  => .5, # percentage of affected in the binary trait
                      :e19e => .8, # edit_successsful_rate
                      :w4t  => 2   # weight on the binary trait EBV
                      )
    par = (; Parameters...)     # named tuple.  contents as above

    base = sim_base(1200, 29, 1e8, 1725)
    qtl = sim_QTL(base, par.nQTL...)
    
    p1, _ = breeding_program(base, par, qtl, edit=false)
    p2, _ = breeding_program(base, par, qtl, edit=true)
end
