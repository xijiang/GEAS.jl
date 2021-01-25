# A very brief manual for GEAS
- Xijiang Yu

## Make a local copy of GEAS
As this project is curretly set as private,
there is only one way to have a local copy.
In your working directory, say, `~/pkg`, do
```bash
# as a co-developer
# and assume that you have setup public/private keys on github
git pull git@github.com:xijiang/GEAS.jl GEAS
```

## Prepare data

In, say `~/pkg/GEAS`:
```bash
mkdir -p dat/real/nofima dat/run
```
Then put Binyam's plink files in `~/pkg/GEAS/dat/real/nofima`.

## Run the functions
In, say, `~/pkg/GEAS`:
```bash
julia
```
```julia
]
activate .
instantiate # for first the time to install packages if not installed.
<backspace> # to REPL
using GEAS
workflow()
```

## To develop the package
I usually open 3 terminals:
- julia REPL
- emacs -nw
- bash

Above is much easier with `tmux` or `screen`.

In the `REPL`:
```bash
julia
```
```julia
]
activate .
<backspace> # to REPL
using Revise
using GEAS
```

Note the difference between above and its previous one.
Package `Revise` is a must-have for `julia` development.
You edit the `julia` files.  The changes immediately reflected in the `REPL`.

It is better to `branch` my codes if you want to make some modifications.
I will then `merge` them to the `main`.
