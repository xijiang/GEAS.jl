#!/usr/bin/env bash
julia ave.jl $1
pandoc -o ave.pdf -t beamer ave.md
