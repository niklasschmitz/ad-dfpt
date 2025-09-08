#!/bin/sh
## To be called like `bash 1_relax.jl runs/Si.extxyz 8`.
#=
BN=$(basename "$0" .jl)
INPUTFILE=$1
KGRID=$2
julia --project -t4 $BN.jl $INPUTFILE $KGRID 2>&1 | tee ${INPUTFILE}_kgrid${KGRID}_relax.log
exit $?
=#
include("lib.jl")
DFTK.setup_threading()
DFTK.versioninfo()
file = ARGS[1]
nkpt = parse(Int, ARGS[2])
Ecut = :auto
kgrid = (nkpt, nkpt, nkpt)
prefix = "$(file)_kgrid$(nkpt)_relax"

system = load_system(file)
pseudopotentials = PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf")
kinetic_blowup = BlowupIdentity()  # Use no kinetic energy smearing
smearing = Smearing.Gaussian()
temperature = 1e-3
θ = zeros(1)  # Dummy variable to make optimal_lattice happy... TODO avoid this
functional(θ) = PBE()
model_init = model_DFT(system; pseudopotentials, functionals=functional(θ),
                       kinetic_blowup, smearing, temperature)

if Ecut == :auto
    Ecut = recommended_cutoff(model_init).Ecut
end

# Phase 1: Lattice relaxation
function case(p=ComponentVector(; scale=0.))
    lattice = (1 + p.scale) * model_init.lattice
    (; name=prefix, p, lattice, model_init.atoms, model_init.positions, Ecut, kgrid)
end
p_opt = optimal_lattice(case, functional, θ; Ecut, kgrid, temperature, smearing, kinetic_blowup, g_abstol=1e-8)  # TODO How reasonable?

# Construct model at relaxed geometry
case_opt = case(p_opt)
model0 = Model(model_init; lattice=case_opt.lattice, positions=case_opt.positions)

# Phase 2: check stresses here as sanity check, to see that xy yz xz stresses are small too
@info now()
basis0 = @time "Basis" PlaneWaveBasis(model0; Ecut, kgrid)
display(basis0)
scfres0 = @time "SCF" self_consistent_field(basis0; tol=1e-10)
stresses0 = @time "Stress" DFTK.full_stress_to_voigt(compute_stresses_cart(scfres0))
@info now() stresses0
save_scfres(prefix * "_scfres0.jld2", scfres0; extra_data=Dict("stresses0"=>stresses0), save_ψ=false)
save_scfres(prefix * "_scfres0.json", scfres0; extra_data=Dict("stresses0"=>stresses0))
