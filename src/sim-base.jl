"""
    function sim_base()
---
Simulate a base population with [`macs`](https://github.com/gchen98/macs).
"""
function sim_base()
    ne1d = `macs 4000 100000000 -t 0.00001 -r 0.000004
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
                -eN 2000.00 800.00`
    ne1k = `macs 4000 100000000 -t 0.0001 -r 0.00004
                -eN   0.50  2.00
                -eN   0.75  2.50
                -eN   1.00  3.00
                -eN   1.25  3.20
                -eN   1.50  3.50
                -eN   1.75  3.80
                -eN   2.00  4.00
                -eN   2.25  4.20
                -eN   2.50  4.50
                -eN   5.00  5.46
                -eN  10.00  7.37
                -eN  15.00  9.28
                -eN  20.00 11.19
                -eN  25.00 13.10
                -eN  50.00 22.66
                -eN 100.00 41.77
                -eN 150.00 60.89
                -eN 200.00 80.00`
    run(pipeline(ne1d,
                 stdout = joinpath(dat_dir, "run/sim/t.txt"),
                 stderr = joinpath(dat_dir, "run/sim/d.txt")))
end

    #=
    2>debug.txt | msformatter > haplotypes.txt
tail --lines +6 haplotypes.txt >Intermediate.txt
tail --lines +2 Intermediate.txt > MacsHaplotypes.txt
head --lines 1 Intermediate.txt > TempMap.txt
tail --bytes=+11 TempMap.txt > PhysicalMapInput.txt
grep segsites haplotypes.txt > TempSegSites.txt
tail --bytes=+10 TempSegSites.txt > SegSites.txt
rm Intermediate.txt
rm TempMap.txt
rm TempSegSites.txt
sed -i 's/0/ 0/g' MacsHaplotypes.txt
sed -i 's/1/ 1/g' MacsHaplotypes.txt
cp SegSites.txt FinishedMacs.txt


#### For $N_e=1000$
./ \
       2>debug.txt |
    ./msformatter > haplotypes.txt
tail --lines +6 haplotypes.txt >Intermediate.txt
tail --lines +2 Intermediate.txt > MacsHaplotypes.txt
head --lines 1 Intermediate.txt > TempMap.txt
tail --bytes=+11 TempMap.txt > PhysicalMapInput.txt
grep segsites haplotypes.txt > TempSegSites.txt
tail --bytes=+10 TempSegSites.txt > SegSites.txt
rm Intermediate.txt
rm TempMap.txt
rm TempSegSites.txt
sed -i 's/0/ 0/g' MacsHaplotypes.txt
sed -i 's/1/ 1/g' MacsHaplotypes.txt
cp SegSites.txt FinishedMacs.txt
=#
