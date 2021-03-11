"""
    function copyright()
---
Show copyright information.
"""
function copyright()
    w = 60
    year = Dates.year(now())
    msg = ["",
           "GEAS: Genome Editing Assisted Selection",
           "===========",
           "Developed by Xijiang Yu @NMBU",
           "Copyright © $year",
           "MIT license"]
    @info join(lpad.(msg, w), "\n")
end

function isLaterBeagle(a, b)
    mon = Dict("Jan" => 1, "Feb" => 2, "Mar" => 3,
               "Apr" => 4, "May" => 5, "Jun" => 6,
               "Jul" => 7, "Aug" => 8, "Sep" => 9,
               "Oct" => 10, "Nov" => 11, "Dec" => 12, "mmm" => 0)
    da, db = a[8:9], b[8:9]
    ma, mb = mon[a[10:12]], mon[b[10:12]]
    ya, yb = a[13:14], b[13:14]
    if ya == yb
        if ma == mb
            return da > db ? a : b
        else
            return ma > mb ? a : b
        end
    else
        return ya > yb ? a : b
    end
end

"""
    function update_beagle()
---
Return the latest beagle URL
"""
function update_beagle()
    BeagleURL = "https://faculty.washington.edu/browning/beagle/"
    # Note Beagle is not version controlled by git.
    # It seems that names for new versions are following the pattern like
    # beagle.ddMmmyy.xxx.jar
    # so I just grab all the file names in $BeagleURL
    # There also must be better way to get the URL, like XML?
    beagle = "beagle.00mmm00.xxx.jar"
    web = HTTP.request("GET", BeagleURL)
    txt = split(String(web.body), "\n")
    for line in txt
        r = match(r"beagle\.\d{2}\w{3}\d{2}.{5}jar", line)
        r == nothing || (beagle = isLaterBeagle(beagle, r.match))
    end
    @info "Updating Beagle to $beagle."
    local_beagle = joinpath(bin_dir, "beagle.jar")
    isfile(local_beagle) && rm(local_beagle, force=true)
    download(joinpath(BeagleURL, beagle), local_beagle)

    beagle2vcfURL = "https://faculty.washington.edu/browning/beagle_utilities/beagle2vcf.jar"
    beagle2vcf = joinpath(bin_dir, "beagle2vcf.jar")

    if !isfile(joinpath(bin_dir, beagle2vcf))
        @info "Downloading beagle2vcf.jar"
        download(beagle2vcfURL, beagle2vcf)
    end
end

"""
    function update_plink()
---
Return URL of the latest plink version 1 for Linux x86_64.
"""
function update_plink()
    @info "Updating plink to the latest"
    plinkURL = "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip"
    cd(joinpath(dat_dir, "test"))
    download(plinkURL, "plink.zip")
    run(`unzip -f plink.zip`)
    cp("plink", joinpath(bin_dir, "plink"), force = true)
    cd(prj_dir)
end


"""
    function update_macs()
---
Clone and compile marcs into `bin`.
"""
function update_macs()
    @info "Updating macs"
    marcs = "https://github.com/gchen98/macs"
    cd(dat_dir)
    if isdir("macs")
        cd("macs")
        run(`git pull --no-rebase`)
        cd("..")
    else
        run(`git clone $marcs`)
    end
    cd("macs")
    run(`g++ -o macs -O3 -Wall simulator.cpp algorithm.cpp datastructures.cpp`)
    run(`g++ -o msformatter -O3 -Wall msformat.cpp`)
    cp("macs", joinpath(bin_dir, "macs"), force = true)
    cp("msformatter", joinpath(bin_dir, "msformatter"), force = true)
    cd(prj_dir)
end


"""
    Update()
---
Every time one starts this package, the package will automatically check if `plink` and `beagle.jar` exist.
If not, the package will download the latest version of these two files.
The package will update, or download them on 17th every month anyway.
"""
function Update()
    @info "Updating binaries"
    isdir(bin_dir) || mkdir(bin_dir)
    update_beagle()
    update_plink()
    update_macs()
end
