#!/bin/bash

PS4='\e[35mRUNNING: \e[m'
set -x

salloc -t 1 --mem 16G --ntasks-per-node 32 --nodes 1 --partition m12 -- julia -t 32 --project=LifeGame.jl time-thread.jl
salloc -t 1 --mem  8G --ntasks-per-node  4 --nodes 8 --partition m12 -- mpirun julia --project=LifeGameMPI.jl time-mpi.jl
