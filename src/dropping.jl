"""
    function random_mate(nprt, noff; ratio::Int = 1)
---
Assume the first half of the ID are males, and the second half are females,
this function randomize parent pairs.
Each pair will have `noff` offspring.
The few odds of `nprt` might be ommitted.
A warning message will be given if this happens.
The `ratio` here is female:male, which should always be an integer.
"""
function random_mate(nprt, noff; ratio::Int = 1)
    np = Int(nprt ÷ (1 + ratio)) # number of pa
    nm = Int(floor(np * ratio))
    omt = nprt - np - nm
    omt > 0 && @warn "The last $omt ID will ommitted"
    pa = repeat(shuffle(1:np),
                inner = ratio * noff)
    ma = repeat(shuffle(np+1:np+nm),
                inner = noff)
    [pa ma]
end
    
"""
    function haldane(pos; M = 1e10)
---
Return recommbination rate column given bp distances.

Haldane function:
- Recombination rate: r = .5(1 - exp(-2d))
- Distance in Morgan: d = -.5log(1 - 2r)
"""
function haldane(pos; M = 1e8)
    chr = bp = 0
    nl = size(pos)[1]
    r = zeros(Float64, nl)
    @inbounds for i in 1:nl
        if pos[i, 1] == chr
            d = (pos[i, 2] - bp) / M # distance in cM
            r[i] = .5(1 - exp(-2d))
        else
            chr = pos[i, 1]
            r[i] = 0.5
        end
        bp = pos[i, 2]
    end
    r
end

"""
    function gdrop(base, ped)
---
Drop genotypes of `base`, through `ped`.

ToDo:
Needs optimization
"""
function gdrop(snp, ped, r)
    noff = size(ped)[1]
    nsnp = size(snp)[1]
    osnp = Array{Int8, 2}(undef, nsnp, 2noff)
    i = 1
    sample = zeros(nsnp)
    @inbounds for prt in eachrow(ped)
        @inbounds for p in prt
            hap = view(snp, :, 2p-1:2p)
            off = view(osnp, :, i)
            k   = 1
            rand!(sample)
            @inbounds for j in 1:nsnp
                sample[j] < r[j] && (k = 3 - k)
                off[j] = hap[j, k]
            end
            i += 1
        end
    end
    osnp     
end

"""
Seems OK
"""
function test_drop()
    snp = rand(0:1, 10, 8)
    ped = [1 2
           3 4
           ]
    r = ones(10).* 0.1
    r[1] = .5
    off = gdrop(snp, ped, r)
    [snp off]
end
