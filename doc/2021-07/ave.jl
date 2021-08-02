#!/usr/bin/env julia
include("sum.jl")
sum_mmv(".", parse(Int, ARGS[1]))
