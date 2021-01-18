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
        Dict(:pos => reshape(pos, (2, :))', :hap => reshape(hap, (nid*2, :))')
    end
end

"""
    function macs()
---
Simulate a base population with `macs`.
"""
function macs()
    ## Let the time and the population size at the time are
    #t = [50., 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000,
    #     4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 40000, 60000,
    #     80000, 100000, 200000, 400000, 600000, 800000, 1000000]
    ## this is suppose to give an Nₑ = 100
    #mut = 1e-8
    #nid = 600
    #chr = 1e8
    #length(t)/sum(1. ./ t)
    @warn "To be written"
end
