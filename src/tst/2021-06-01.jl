#!/usr/bin/env julia

using JLD2, Distributions, Serialization
import GEAS:sim_pt_QTL, breeding_program, top_QTL

"""
    function run_2021_06_01(jld; rst = "")
---
## Description
- 20 repeats
- 3 comparisons
  - Normal genomic selection, no known QTL
  - Known 5 QTL as fixed effects
  - Edit the biggest of known 5 QTL in generation 2-6

## parameters
1. jld: the base
2. nqtl: e.g., 500
3. rt: correlation of QTL effects between 2 traits
4. nsib: number of sibs challenged and in production. total 2sib
5. nrpt: number of repeats
6. rst: result path
"""
function run_2021_06_01(jld, nqtl, rt, nsib, nrpt; rst = "")
    @load jld base
    if length(rst) > 1
        isdir(rst) || mkpath(rst)
    end

    w = ndigits(nrpt)
    for r in 1:nrpt
        d = MvNormal(zeros(2), [1 rt; rt 1])
        qtl = sim_pt_QTL(base, nqtl, d)
        
        suffix = lpad(string(r), w, '0')
        serialize(joinpath(rst, "qtl-$suffix.ser"), qtl)
        par = begin
            Parameters = Dict(
                :nSire=> 100,
                :nDam => 200,       # male:female = 1:2
                :nSib => 2nsib,
                :nC7e => nsib,      # number of sib to be challenged
                :nG8n => 6,         # number of generations
                :nQTL => [nqtl, nqtl],
                :h²   => [.5, .5],
                :p8e  => .5,        # death rate after challenge
                :e19e => 1.,        # edit_successsful_rate
                :w4t  => 1,         # weight on the binary trait EBV
                :nk3n => 5          # number of known QTL on binary trait
            )
            (; Parameters...)       # named tuple.  contents as above
        end

        # genomic selection with no known QTL
        prd, snp = breeding_program(base, par, qtl)
        serialize(joinpath(rst, "prd-a-$suffix.ser"), prd)
        serialize(joinpath(rst, "snp-a-$suffix.ser"), snp)

        # genomic selection with known QTL as fixed effects
        prd, snp = breeding_program(base, par, qtl, fixed=true)
        serialize(joinpath(rst, "prd-b-$suffix.ser"), prd)
        serialize(joinpath(rst, "snp-b-$suffix.ser"), snp)

        # genomic selection, edit biggest known QTL each generation
        prd, snp = breeding_program(base, par, qtl, edit=true)
        serialize(joinpath(rst, "prd-c-$suffix.ser"), prd)
        serialize(joinpath(rst, "snp-c-$suffix.ser"), snp)
    end
end

# run_2021_06_01("dat/run/base.jld2", 500, -0.2, 40, 20, rst="dat")

#=
Note:
- this is a task discription to be submitted on slurm system
- below is slurm task bash script.
- copy below to `example.sh`
- uncomment the function call
- `sbatch example.sh`

#!/usr/bin/env bash
############################################################
# Environment setup
############################################################
#SBATCH --job-name=GEAS
#SBATCH --account=nn9891k
#SBATCH --job-name=GEAS
#SBATCH --time=0-10:00:00       # d-hh:mm:ss
#SBATCH --mem-per-cpu=2G        # 10x2 => 20G memory
#SBATCH --cpus-per-task=10      # 80 threads available
##sbatch --ntasks=6             # 190/20 => 9 independent tasks
##sbatch --ntasks-per-node=2    # not sure for the moment

# other job-controls may be needed:
#  - nodes=;
#  - ntasks-per-node=;
#  - cpus-per-task=;

log=results.geas

############################################################
# copy to workdirectory on the computation node
# i always write the working scripts in Julia
############################################################
cp $SUBMITDIR/geas.{sh,jl} $SCRATCH/
cp $SUBMITDIR/../dat/base.jld2 $SCRATCH

############################################################
# run the task
############################################################
cd $SCRATCH
echo "Time start: $(date)" > $log
./geas.jl
echo "Time stop: $(date)" >> $log

############################################################
# collect results
############################################################
cp *ser $log $SUBMITDIR/
=#
