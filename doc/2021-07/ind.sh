#!/usr/bin/env bash
julia ind.jl $1
pandoc -o ind.pdf -t beamer ind.md
