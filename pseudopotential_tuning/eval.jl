#!/bin/sh
#=
BN=$(basename "$0" .jl)
TRAINWORKDIR=$1
WORKDIR=${TRAINWORKDIR}/eval
mkdir $WORKDIR
julia --project -t4 $BN.jl $TRAINWORKDIR $WORKDIR 2>&1 | tee ${WORKDIR}/${BN}.log
exit $?
=#
include("lib.jl")
DFTK.setup_threading()
datadir = @__DIR__
trainworkdir = ARGS[1]
workdir = ARGS[2]
trajectory_file = joinpath(trainworkdir, "trajectory.json")
# trajectory_file = "train0.0/trajectory.json"

# Load GTH PBE pseudopotential table
pseudopotentials = PseudoFamily("cp2k.nc.sr.pbe.v0_1.largecore.gth")

# Load initial parameters for Li-q1 (only LDA available for q1)
Li_q1_init = ElementPsp(:Li, PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"))

# Expose optimisable parameters
make_li_q1(θ) = unflatten(ComponentVector(; θ..., rloc=Li_q1_init.psp.rloc))
params = flatten(Li_q1_init)
tunable_param_keys = [:cloc1, :cloc2, :rp1, :rp2, :h1, :h2]
θ_init = params[tunable_param_keys]

# Load trained parameters
trajectory = JSON3.read(trajectory_file)
θ = deepcopy(θ_init)
θ .= trajectory[end]["theta"]
Li_q1 = make_li_q1(θ)
@info "Loaded params" trajectory_file θ_init θ

# DFT settings
functionals = PBE()
smearing    = Smearing.Gaussian()
temperature = 0.00225
tol         = 1e-7

# Load systems
systems = [
    (; file="Li-BCC.extxyz", Ecut= 20, kgrid=(8,8,8)),
    (; file="Li-XO.extxyz",  Ecut=120, kgrid=(8,8,8)),
]
volume_scaling_range = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]

for (; file, Ecut, kgrid) in systems
    start_time = now()
    system = load_system(joinpath(datadir, file))
    name = system.system_data.name
    @info start_time name Ecut kgrid

    # Update pseudopotential on each Li atom
    (; lattice, atoms, positions) = DFTK.parse_system(system, pseudopotentials)
    atoms = map(atoms) do atom
        (element_symbol(atom) == :Li) ? Li_q1 : atom
    end
    model0 = model_DFT(lattice, atoms, positions; functionals, temperature, smearing)

    energies = map(volume_scaling_range) do scale
        model = Model(model0; lattice=cbrt(scale) * model0.lattice)
        basis = PlaneWaveBasis(model; Ecut, kgrid)
        scfres = self_consistent_field(basis; tol)

        scf_path = joinpath(workdir, "$(name)_v$(scale)_Ecut$(Ecut)_kgrid$(kgrid)_scfres.jld2")
        save_scfres(scf_path, scfres; save_ψ=false)

        energy_per_atom = scfres.energies.total / length(model.atoms)
    end

    end_time = now()
    results = (;
        theta=θ,
        trajectory_file,
        file,
        functionals=string(functionals),
        Ecut,
        kgrid,
        temperature,
        smearing,
        volume_scaling_range,
        energies,
        date=end_time,
        time=string(canonicalize(end_time - start_time))
    )
    save_path = joinpath(workdir, "$name.json")
    open(save_path, "w") do io
        JSON3.pretty(io, results)
    end
end
