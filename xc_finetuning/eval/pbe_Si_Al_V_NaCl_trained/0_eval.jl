#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t4 $BN.jl 2>&1 | tee $BN.log
exit $?
=#
using AtomsIO
using Dates
using DFTK
using DftFunctionals
using JSON3
using JLD2
DFTK.setup_threading()
include(joinpath(@__DIR__, "../../sol58lc/load.jl"))
include(joinpath(@__DIR__, "../../../lattice_relaxation.jl"))

file = ARGS[1]

name, ext = splitext(basename(file))
system = load_system(file)
case = convert_system_to_case(system)
(; Ecut, kgrid) = case()

make_pbe(θ) = [PbeExchange(θ), DftFunctional(:gga_c_pbe)]
θ = ComponentVector(κ=0.5274266054802829, μ=0.2225734371774011)

# Use smearing with a small temperature for all systems
temperature = 1e-3
smearing = Smearing.Gaussian()

# Lattice relaxation
@info now() name Ecut kgrid
p_opt = @time "$name" optimal_lattice(case, make_pbe, θ; Ecut, kgrid, smearing, temperature)
a_opt = p_opt.a

# For convenience: one final SCF
(; lattice, atoms, positions) = case(p_opt)
model = model_DFT(lattice, atoms, positions; functionals=make_pbe(θ),
                  temperature, smearing)
basis = PlaneWaveBasis(model; Ecut, kgrid, kshift=(0, 0, 0))
scfres = self_consistent_field(basis; tol=1e-8)
extra_data = Dict("a_opt"=>a_opt, "name"=>name)
save_scfres("$(name).jld2", scfres; extra_data, save_ψ=false)
save_scfres("$(name).json", scfres; extra_data)

