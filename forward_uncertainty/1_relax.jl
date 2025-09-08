#!/bin/sh
## To be called like `bash 1_relax.jl mp-134.cif 4`.
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
# kgrid = :auto
kgrid = (nkpt, nkpt, nkpt)
prefix = "$(basename(file))_kgrid$(nkpt)_relax"

system = load_system(file)
case = convert_system_to_case(system)
(; Ecut) = case()  # Ecut recommended by pseudopotentials

temperature = 1e-3
smearing = Smearing.Gaussian()
kinetic_blowup = BlowupCHV()

# Load functional and mean parameters
functional = make_beef
θ = params_beef_2005().θ_bf

# functional = make_beefvdw
# θ = params_beefvdw_2012().θ_bf

# Lattice relaxation
@info now() functional θ
p_opt = optimal_lattice(case, functional, θ; Ecut, kgrid, temperature, smearing,
                        kinetic_blowup, g_abstol=1e-8)

# For convenience: one final SCF
(; lattice, atoms, positions) = case(p_opt)
model0 = model_DFT(lattice, atoms, positions; functionals=functional(θ),
                   temperature, smearing, kinetic_blowup)
@info now()
basis0 = @time "Basis" PlaneWaveBasis(model0; Ecut, kgrid)
display(basis0)
scfres0 = @time "SCF" self_consistent_field(basis0; tol=1e-10)
stresses0 = @time "Stress" DFTK.full_stress_to_voigt(compute_stresses_cart(scfres0))
@info now() stresses0 p_opt.a
extra_data=Dict("stresses0"=>stresses0, "a_opt"=>p_opt.a)
save_scfres(prefix * "_scfres0.jld2", scfres0; extra_data, save_ψ=false)
save_scfres(prefix * "_scfres0.json", scfres0; extra_data)

