struct Trait
    nqtl::Int
    h²  ::Float64
end

struct QTL
    pos
    effect
    mean::Float64               # mean and max are about the base population
    max::Float64
end
