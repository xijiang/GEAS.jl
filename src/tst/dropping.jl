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
