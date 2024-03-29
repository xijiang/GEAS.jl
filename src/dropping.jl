"""
    function random_mate(nsire, ndam, nsib)
---
Random mate `nsire` with `ndam`, each dam will have `nsib`.
Assuming the first `nsire` are sires.
The rest are dams.
This is for generation one.
And can be replaced by call `random_mate(shuffle(1:nid), nsire, ndam, nsib)`.
"""
function random_mate(nsire, ndam, nsib)
    r  = Int(ceil(ndam/nsire))
    pa = repeat(repeat(shuffle(1:nsire), r), inner=nsib)
    ma = repeat(1:ndam, inner = nsib) .+ nsire
    [pa[1:length(ma)] ma]
end

"""
    function random_mate(ID, nsire, ndam, nsib)
---
Random mate `ID`.
Its first `nsire` are sires, last `ndam` are dams.

**Note**: `nsire + ndam` must be of `length(ID)`.
"""
function random_mate(ID, nsire, ndam, nsib)
    r  = Int(ceil(ndam/nsire))
    pa = repeat(repeat(shuffle(ID[1:nsire]), r), inner = nsib)
    ma = repeat(ID[(nsire+1):(ndam+nsire)], inner = nsib)
    [pa[1:length(ma)] ma]
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
Drop genotypes of a subset or all of `snp`, through `ped` into `osnp`.
The last will be returned.

## Note:
If to use `osnp = BitArray(undef, nsnp, 2noff)` to save memory,
the speed is greatly reduced.

## ToDo:
Needs optimization.
May be parallel into 10 tasks.
Least urgent.
"""
function gdrop(snp, ped, r)
    noff = size(ped)[1]
    nsnp = size(snp)[1]
    osnp = zeros(Bool, nsnp, 2noff)
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
