using Serialization, Statistics, DataFrames, Plots
function qfrq(snp)
    fq = zeros(11, 10)
    start = 1
    i = 1
    for stop = 4000:4000:44000
        f = mean(snp[:, start:stop], dims=2)
        fq[i, :] = f'
        i += 1
        start = stop + 1
    end
    fq
end

"""
Generation mean and variance of `prd-abc` on `c`, e.g., `:tbv`.
"""
function mv(c, prd...)
    m = zeros(11, 3)
    v = zeros(11, 3)
    for i in 1:3
        gp = groupby(prd[i], :g8n)
        t = combine(gp, c => mean => :m, c => var => :v)
        m[:, i] = t.m
        v[:, i] = t.v
    end
    return m, v
end

function sum_rpt(dir, rpt)
    ir = lpad(string(rpt), 2, '0')
    @info "Read top QTL of the first repeat"
    log = readlines(joinpath(dir, "simu.log"))
    t = 1
    for i in 1:length(log)
        t += 1
        length(log[i]) > 0 && log[i][end] == ':' && break
    end
    t += 3(rpt - 1)
    tq = abs.(parse.(Int, split(log[t])))

    @info "Read SNP of the first repeat, ordinary GS"
    snp = deserialize(joinpath(dir, "snp-a-$ir.ser"))
    fq = qfrq(snp[tq, :])
    plot(fq, leg=false, title = "ogs")
    savefig(joinpath(dir, "afq.pdf"))

    @info "Read SNP of the first repeat, top QTL as fixed effects"
    snp = deserialize(joinpath(dir, "snp-b-$ir.ser"))
    fq = qfrq(snp[tq, :])
    plot(fq, leg=false, title = "fxq")
    savefig(joinpath(dir, "bfq.pdf"))
    
    @info "Read SNP of the first repeat, top QTL edited in turn"
    snp = deserialize(joinpath(dir, "snp-c-$ir.ser"))
    fq = qfrq(snp[tq, :])
    plot(fq, leg=false, title = "edt")
    savefig(joinpath(dir, "cfq.pdf"))

    @info "Read the production trait"
    lb = ["ogs" "fxq" "edt"]
    pa = deserialize(joinpath(dir, "prd-a-$ir.ser"))
    pb = deserialize(joinpath(dir, "prd-b-$ir.ser"))
    pc = deserialize(joinpath(dir, "prd-c-$ir.ser"))
    m, v = mv(:t2c, pa, pb, pc)
    plot(m, label = lb, xlabel="generation", ylabel="Bin mean")
    savefig(joinpath(dir, "bin-mean.pdf"))
    plot(v, label = lb, xlabel="generation", ylabel="Bin var")
    savefig(joinpath(dir, "bin-var.pdf"))
    m, v = mv(:tbv, pa, pb, pc)
    plot(m, label = lb, xlabel="generation", ylabel="Prd mean")
    savefig(joinpath(dir, "prd-mean.pdf"))
    plot(v, label = lb, xlabel="generation", ylabel="Prd var")
    savefig(joinpath(dir, "prd-var.pdf"))
end

function sum_mmv(dir, nr)
    mm = zeros(11, 3)
    vm = zeros(11, 3)
    for i in 1:nr
        r = lpad(string(i), 2, '0')
        pa = deserialize(joinpath(dir, "prd-a-$r.ser"))
        pb = deserialize(joinpath(dir, "prd-b-$r.ser"))
        pc = deserialize(joinpath(dir, "prd-c-$r.ser"))
        m, v = mv(:t2c, pa, pb, pc)
        mm += m
        vm += v
    end
    mm ./= nr
    vm ./= nr
    lb = ["ogs" "fxq" "edt"]
    plot(mm, label = lb, xlabel="generation", ylabel="Average bin mean")
    savefig(joinpath(dir, "ave-bin-mean.pdf"))
    plot(vm, label = lb, xlabel="generation", ylabel="Average bin var")
    savefig(joinpath(dir, "ave-bin-var.pdf"))
end
