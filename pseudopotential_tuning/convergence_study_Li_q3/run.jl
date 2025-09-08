#!/bin/sh
#=
BN=$(basename "$0" .jl)
FILE=$1
ECUT=$2
KGRID=$3
LOGFILE=$(basename "$FILE" .extxyz)_ecut${ECUT}_kgrid${KGRID}.log
julia --project -t4 $BN.jl $FILE $ECUT $KGRID 2>&1 | tee $LOGFILE
exit $?
=#
using AtomsIO
using Dates
using DFTK
using PseudoPotentialData
using JSON3
using JLD2
DFTK.setup_threading()


file = ARGS[1]
Ecut = parse(Int, ARGS[2])
kgrid = parse(Int, ARGS[3])
workdir = "."

# Load system from file
system = load_system(file)
name = system.system_data.name
resname = "$(name)_Ecut$(Ecut)_kgrid$(kgrid)"

# Load GTH PBE pseudopotential table
pseudopotentials = PseudoFamily("cp2k.nc.sr.pbe.v0_1.semicore.gth")

# DFT settings
functionals = PBE()
smearing    = Smearing.Gaussian()
temperature = 0.00225
tol         = 1e-4
kgrid       = (kgrid, kgrid, kgrid)

model0 = model_DFT(system; pseudopotentials, functionals, smearing, temperature)
basis0 = PlaneWaveBasis(model0; Ecut, kgrid)
show(stdout, "text/plain", basis0)

start_time = now()

volume_scaling_range = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]
energies = map(volume_scaling_range) do scale
    model = Model(model0; lattice=cbrt(scale) * model0.lattice)
    basis = PlaneWaveBasis(model; Ecut, kgrid)
    scfres = self_consistent_field(basis; tol)

    scf_path = joinpath(workdir, "$(resname)_v$(scale)_scfres.jld2")
    save_scfres(scf_path, scfres; save_ψ=false)

    scfres.energies.total / length(model.atoms)
end

end_time = now()
results = (;
    file,
    functionals=string(functionals),
    pseudopotentials=pseudopotentials.identifier,
    Ecut,
    kgrid,
    temperature,
    smearing,
    volume_scaling_range,
    energies,
    date=end_time,
    time=string(canonicalize(end_time - start_time))
)
save_path = joinpath(workdir, "$(resname).json")
open(save_path, "w") do io
    JSON3.pretty(io, results)
end


