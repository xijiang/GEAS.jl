"""
    function A_matrix(ped)
---
Give a matrix of 2 columns:
- Sire
- Dam

With row number as ID number, this function return a **A** matrix.
It is assumed that row number is the `ID` number.
`Sire`s and `Dam`s are of integers less than their offspring `ID`.
This function cannot deal very large pedigree, i.e.,
with 32G momory, this function can only calculate a pedigree with
~90k ID.
With half precision, this can be ~130k.
With double precision, this function can only handle ~65k ID.
"""
function A_matrix(ped)
    n = size(ped, 1)
    A = zeros(Float32, n, n)
    for i in 1:n
        ip, im = ped[i, :]      # ID i's pa and ma
        A[i, i] = 1.            # below avoid swap
        ip >0 && im >0 && (A[i, i] += .5(A[ip, im] + A[im, ip]))
        for j in i+1:n          # upper triangular
            jp, jm = ped[j, :]
            jp > 0 && (A[i, j] += .5A[i, jp])
            jm > 0 && (A[i, j] += .5A[i, jm])
            A[j, i] = A[i, j]
        end
        i % 100 == 0 && (print("\t$i"))
    end
    A
end

"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + .5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return .5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return .5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end
