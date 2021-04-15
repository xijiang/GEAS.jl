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
