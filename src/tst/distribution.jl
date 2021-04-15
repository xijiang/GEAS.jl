"""
    function distribution()
---
This function is to test how to simulate correlated traits.
"""
function distribution(;nqtl = 1000)
    @info "Simulate a larege population"
    # to simulate
    #base = sim_base(20000, 24, 1e8, 4000, maf=0.01)
    #serialize("dat/run/large-base.ser", base)
    # to use simulate results above
    base = deserialize("dat/run/large-base.ser")
    ds = [Gamma(0.4),
          Laplace(),
          Normal(),
          Uniform(),
          Weibull(),
          Chi(1)]
    label = [L"\Gamma(0.4)", "Laplace()", L"N(0, 1)", L"U(0, 1)", "Weibull()", L"\chi^2(1)"]
    for i in 1:length(ds)
        qtl = sim_QTL(base, nqtl, d=ds[i])[1]
        qgt = qtl_gt(base[:hap], qtl.pos)
        tbv = qgt'qtl.effect
        histogram(tbv, normalize=true, label=label[i], dpi=300)
        savefig("dat/run/$i.pdf")
    end
end
