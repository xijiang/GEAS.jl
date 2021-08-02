#!/usr/bin/env julia
include("sum.jl")
sum_rpt(".", parse(Int, ARGS[1]))
