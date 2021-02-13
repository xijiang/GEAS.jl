import GZip
"""
    function gt600()
---
1. Phase the genotypes with beagle
2. Convert results to array
3.
"""
function gt600()
    rst = joinpath(dat_dir, "run")
    fra = joinpath(dat_dir, "real/nofima/for_CMSEdit")
    Nₑ  = 100
    
    @info join(["",
                "Remove SNP of duplicates and unknown postions",
                "  - See data/real/nofima/plink.sh"], "\n")

    @info "Phasing with beagle.jar"
    
    # make alleles to C/T, and replace positions, so that beagle can phase them
    GZip.open("$rst/b.vcf.gz", "w") do io
        i = 0
        for line in eachline(GZip.open("$rst/a.vcf.gz", "r"))
            if line[1] == '#'
                write(io, line, "\n")
            else
                i += 1
                t  = split(line)
                t[4] = "T"
                t[5] = "C"
                write(io, join(t, "\t"), "\n")
            end
        end
    end

    # phasing
    run(`java -jar $beagle gt=$rst/b.vcf.gz ne=$Nₑ out=$rst/c`)
end

"""
    function vcf2dic(vcf)
---
Read chromosome and haplotypes from a `vcf` file of **gzip** format.
Returns a Dictionary: pos=>[chr, bp], hap=>hap of chr8[][].
Each column of hap is a haplotype across genome.
"""
function vcf2dic(vcf)
    GZip.open(vcf, "r") do io
        # skip vcf header and determine N_ID
        nid = 0
        for line in eachline(io)
            if line[2] != '#'
                nid = length(split(line)) - 9
                break
            end
        end
        pos = Int[]
        hap = Int8[]
        for line in eachline(io)
            f = split(line)
            append!(pos, parse.(Int, f[1:2]))
            append!(hap, parse.(Int8, collect(join(f[10:end], ' ')[1:2:end])))
        end
        pos = reshape(pos, (2, :))'
        r   = haldane(pos)
        hap = reshape(hap, (nid*2, :))'
        Dict(:pos => pos,
             :r => r,
             :hap => Bool.(hap))
    end
end

"""
    function sim_base(nChr)
---
This procedure were copied from _Jenko J et. al_ 2015.
The command included will simulate a population of N_e = 100.
The history is to mimic cattle history.
I may think some other history structure later.

- mt: mutation rate per site per 4N generations
- rc: recombination rate per site per 4N generations
"""
function sim_base(nID, nChr, nSNP; nBP = 100_000_000, mt=0.00001, rc=0.000004)
    pcmd = pipeline(`$bin_dir/macs $nID $nBP -t $mt -r $rc
				-eN    0.03   1.75
				-eN    0.06   2.00
				-eN    0.13   3.50
				-eN    0.25   5.00
				-eN    0.50   7.00
				-eN    0.75   8.20
				-eN    1.00   8.50
				-eN    1.25   9.00
				-eN    1.50  10.00
				-eN    1.75  11.00
				-eN    2.00  12.75
				-eN    2.25  13.00
				-eN    2.50  12.00
				-eN    5.00  20.00
				-eN    7.50  25.00
				-eN   10.00  30.00
				-eN   12.50  32.00
				-eN   15.00  35.00
				-eN   17.50  38.00
				-eN   20.00  40.00
				-eN   22.50  42.00
				-eN   25.00  45.00
				-eN   50.00  54.56
				-eN  100.00  73.67
				-eN  150.00  92.78
				-eN  200.00 111.90
				-eN  250.00 131.01
				-eN  500.00 226.58
				-eN 1000.00 417.72
				-eN 1500.00 608.86
				-eN 2000.00 800.00`,
                    `$bin_dir/msformatter`,
                    "$dat_dir/run/hap.txt")
    for i in 1:nChr
    #    println("Simulate chromosome $i")
        run(pipeline(pcmd, stderr="/dev/null")) # supress stderr message.
        # sample SNP & write.
    end
end
