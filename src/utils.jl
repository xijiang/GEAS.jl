function macs_args(nid)
    par = Dict(
        :generic => "$(2nid) 1e8 -t 1e-5 -r 4e-6 -eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0",
        :cattle => "$(2nid) 1e8 -t 9E-6 -r 3.6E-6 -eN 0.011 1.33 -eN 0.019 2.78 -eN 0.036 3.89 -eN 0.053 11.11 -eN 0.069 16.67 -eN 0.431 22.22 -eN 1.264 27.78 -eN 1.819 38.89 -eN 4.875 77.78 -eN 6.542 111.11 -eN 9.319 188.89 -eN 92.097 688.89 -eN 2592.097 688.89",
        :wheat => "$(2nid) 8E8 -t 4E-7 -r 3.6E-7 -eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640",
        :maize => "$(2nid) 2E8 -t 5E-6 -r 4E-6-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100",
        :european => "$(2nid) 1.3E8 -t 0.0483328 -r 0.02054849 -G 1.0195 -eG 0.0001000977 1.0031 -eN 0.0004492188 0.002015625 -eN 0.000449707 0.003634766")
    
    (; par...)
end
