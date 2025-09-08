#!/bin/sh
## To be called like `bash 0_run.jl Si_diamond.extxyz`.
#=
BN=$(basename "$0" .jl)
INPUTFILE=$1
julia --project -t4 $BN.jl $INPUTFILE >runs/$INPUTFILE.log 2>&1
exit $?
=#
using AtomsIO
using Dates
using DFTK
using JLD2
using JSON3
using LinearAlgebra
using Unitful, UnitfulAtomic
using PseudoPotentialData
DFTK.setup_threading()
DFTK.versioninfo()
DFTK.reset_timer!(DFTK.timer)

cd("runs")

data_path = joinpath(@__DIR__, "../xc_finetuning/sol58lc/structures")
pseudopotentials = PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf")
tol_scf = 1e-8
Ecut = 20

file = ARGS[1]
name, _ = splitext(file)
system = load_system(joinpath(data_path, file))
model = model_DFT(system; functionals=PBE(), pseudopotentials)

kspacing = 0.15u"Å^-1"
kgrid = kgrid_from_maximal_spacing(model.lattice, kspacing)

# Cartesian displacement of first atom along (0.1, 0.2, 0.3) direction
# to obtain a force of about 0.5 eV/Å.
system_rescalings = JSON3.read(read("../system_rescalings.json", String))
system_rescalings = Dict(d["name"] => d["rescaling"] for d in system_rescalings)
cartesian_displacement = normalize([0.1, 0.2, 0.3]) * austrip(0.1u"Å") * system_rescalings[split(name, "_")[1]]

updated_positions = copy(model.positions)
updated_positions[1] += DFTK.vector_cart_to_red(model, cartesian_displacement)
model = Model(
    model;
    positions=updated_positions,
    symmetries=false,
)

"""
Run the SCF algorithm.
If `start_with_temperature` is true, run a first SCF with temperature,
then use the resulting density and eigenfunctions to start a second SCF without temperature.
"""
function run_scf(basis; tol)
    # Run a first SCF with temperature
    model_init  = Model(model;
                        temperature=1e-3, smearing=Smearing.Gaussian())
    basis_init  = PlaneWaveBasis(model_init; basis.Ecut, basis.kgrid, basis.fft_size)
    scfres_init = self_consistent_field(basis_init; tol)

    # Compute suitable starting diagtol
    diagtol_first = DFTK.determine_diagtol(AdaptiveDiagtol(), scfres_init)
    diagtolalg = AdaptiveDiagtol(; diagtol_first)

    self_consistent_field(basis; tol, scfres_init.ρ, scfres_init.ψ, diagtolalg)
end

basis = PlaneWaveBasis(model; Ecut, kgrid)
scfres = run_scf(basis; tol=tol_scf)
forces = compute_forces(scfres)
save_scfres("$(name)_kgrid$(tuple(kgrid.kgrid_size...))_Ecut$(Ecut)_scfres.jld2", 
            scfres; extra_data=Dict("forces"=>forces))

println("SCF", DFTK.timer)
DFTK.reset_timer!(DFTK.timer)


# Run refinement
# Determine Ecut from pseudopotential recommendations
Ecuts_atoms = [
    pseudometa(pseudopotentials, Symbol(at.species))["cutoffs_high"]["Ecut"]
    for at in unique(model.atoms)
]
Ecut_max = maximum(Ecuts_atoms)
Ecut_large = 1.5 * Ecut_max
basis_large = PlaneWaveBasis(model; Ecut=Ecut_large, kgrid)
refinement = refine_scfres(scfres, basis_large)
forces_refined = refine_forces(refinement)
jldopen("$(name)_kgrid$(tuple(kgrid.kgrid_size...))_Ecut$(Ecut)_refined_Ecutlarge$(Ecut_large).jld2", "w") do file
    file["ρ"] = refinement.ρ
    file["δρ"] = refinement.δρ
    file["ΩpK_res"] = refinement.ΩpK_res

    (; F, dF) = forces_refined
    file["F_base"] = F
    file["F_refined"] = F + dF
end
println("REFINEMENT", DFTK.timer)
DFTK.reset_timer!(DFTK.timer)


# Verification: run at large cutoff
scfres_large = run_scf(basis_large; tol=tol_scf)
forces_large = compute_forces(scfres_large)
save_scfres("$(name)_kgrid$(tuple(kgrid.kgrid_size...))_Ecut$(Ecut_large)_scfres.jld2", 
            scfres_large; extra_data=Dict("forces"=>forces_large))
println("SCFLARGE", DFTK.timer)

