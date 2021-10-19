"""
    function check_par(par)
---
## Description
Check if enough and right.
"""
function check_par(par)
    @info "Check if the simulation parameters are valid"
    must_keys = Dict(
        :nSire => "Number of sires",
        :nDam  => "Number of dams",
        :nSib  => "Number of sibs",
        :nC7e  => "Number of sibs for challenge(c7e)",
        :nG8n  => "Number of generation(g8n) to breed",
        :h²    => "Heritability for traits",
        :w4t   => "Weight(w4t) for the binary trait",
        :cv    => "Covariance of pleiotropic QTL effects for 2 traits",
        :p17h  => "Percentage-of-death(p17h) after challenge",
        :e19e  => "edit-successful-rate(e19e), [0, 1]",
        :nK3n  => "Number of known(k3n) QTL of the binary trait",
        :dd    => "An ϵ to be added to diag(LHS matrix) to avoid singular",
        :edit  => "Whether to edit known QTL",
        :fix   => "Whether to emphasize known QTL",
        :b4y   => "Whether the binary(b4y) trait to be recorded as binary",
        :u38y  => "use-genetypes-of-current-generation-only(u38y)"
    )
    missingkey = []
    for key in keys(must_keys)
        haskey(par, key) || push!(missingkey, key)
    end
    isempty(missingkey) || begin
        for key in missingkey
            println("Key $key \"$(must_keys[key])\" is missing")
        end
        error("Some must-have key(s) is/are missing")
    end
    haskey(par, :b4y) && !par.b4y && 
        @warn "Continuous disease phenotype was for debugging" par.b4y
    # More to be added
end
