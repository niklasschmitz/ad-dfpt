#!/bin/sh
#=
BN=$(basename "$0" .jl)
WORKDIR=$1
LAM=$2
mkdir $WORKDIR
julia --project -t4 $BN.jl $WORKDIR $LAM 2>&1 | tee ${WORKDIR}/${BN}.log
exit $?
=#
include("lib.jl")
DFTK.setup_threading()
workdir = ARGS[1]
λ = parse(Float64, ARGS[2])
train(workdir; λ)
