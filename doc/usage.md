# A very brief manual for GEAS
- Xijiang Yu

## Make a key pair for GitHub
### In a terminal:
```bash
ssh-keygen
# it is better to set a password for this key pair.
```
### On github.com
1. Click your icon on the top-right corner
2. Select `settings`
3. Click the `SSH and GPG keys`
4. Click `New SSH key`.
5. Paste the content of `id_rsa.pub` here.

There is also a guide on this page.
In the local terminal, you can use `keychain` (http://www.funtoo.org) to manage your ssh-agent.  Please refer its manual.

If you are running a Linux machine, this can usually be automatically managed.

Once you have setup this key pair, you can use command lines to develop your packages in GitHub repo.
For example the manipulations in later sections.

## Make a local copy of GEAS
As this project is curretly set as private,
there is only one way to have a local copy.
In a terminal and in your working directory,
say, `~/workplace`, do

```bash
# as a co-developer
# and assume that you have setup public/private keys on github
git clone git@github.com:xijiang/GEAS.jl GEAS
```

## Prepare data

In, say `~/workspace/GEAS`:
```bash
mkdir -p dat/run
wget https://nmbu.org/tmp/ns.ser
```

These data were imputed and `serialized` for ease of later simulations and tweaks.

## Run the functions
In, say, `~/workspace/GEAS`:
```bash
julia
```
```julia
]
activate .
instantiate # for first the time to install packages if not installed.
<backspace> # to REPL
using Revise			# a must have package for julia devel.
using GEAS
?GEAS.breeding.program		# docs attached to this function
@time workflow()		# the only exported function of GEAS
```

The `workflow` function shows an example on how to use the package. 

## To develop the package
I usually open 3 terminals:
- julia REPL
- emacs -nw
- bash

Above is much easier with `tmux` or `screen`.

Package `Revise` is a must-have for `julia` development.
You edit the `julia` files.  The changes immediately reflected in the `REPL`.

It is better to `branch` my codes if to make some modifications.
I will then `merge` them to the `main`.

## Branching/flow

Below is a workflow, which roughly shows how to branch a git repo.
Please refer https://guides.github.com/introduction/flow/.

## Issues
We can discuss the package with `issues` on the `GEAS` webpage.