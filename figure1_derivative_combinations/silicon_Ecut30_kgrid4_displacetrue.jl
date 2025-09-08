#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t4 $BN.jl 2>&1 | tee $BN.log
exit $?
=#

include("lib.jl")
DFTK.setup_threading()
prefix, _ = splitext(@__FILE__)
run(prefix; Ecut=30, nkpt=4, displace=true)
