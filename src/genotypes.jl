import GZip
"""
    function macs()
---
Simulation of genotypes with macs.
"""
function macs(;nchr = 26, nid = 300, ns=600)
    
    @info join(["Simulation with ms",
                "1",
                ], "\n")
end

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
    
    @info  "Phasing with beagle.jar"
    
    # convert to VCF for beagle to impute and phase
    run(`$plink --bfile $fra --chr-set 30 --recode vcf-iid bgz --out $rst/a`)
    run(pipeline(`zcat $rst/a.vcf.gz`,
                 `grep -v \#`,
                 `gawk '{print $1, $2}'`,
                 "$rst/bp"))
    # make positions unique
    # - if two SNP have same chr and bp, make the later bp += 1000
    positions = begin
        chr = Int[]
        bp = Int[]
        for line in eachline("$rst/bp")
            c, b = parse.(Int, split(line))
            push!(chr, c)
            push!(bp, b)
        end
        DataFrame(chr=chr, bp=bp)
    end
    for df in groupby(positions, :chr)
        for i in 2:length(df.bp)
            if df.bp[i] == df.bp[i-1]
                df.bp[i:end] .+= 1000
            end
        end
    end

    # make alleles to C/T, and replace positions, so that beagle can phase them
    GZip.open("$rst/b.vcf.gz", "w") do io
        i = 0
        for line in eachline(GZip.open("$rst/a.vcf.gz", "r"))
            if line[1] == '#'
                write(io, line, "\n")
            else
                i += 1
                t  = split(line)
                t[2] = string(positions.bp[i])
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
    function simBase()
---
Simulate a base population with `macs`.
"""
function simBase()
end
