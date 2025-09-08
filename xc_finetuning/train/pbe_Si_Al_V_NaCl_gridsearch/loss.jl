#!/bin/sh
#=
BN=$(basename "$0" .jl)
TS=$(date +"%Y%m%d_%H%M%S")
LOGFILE="$1/${BN}_${TS}_kappa$2_mu$3.log"
julia --project -t4 "$BN.jl" "$1 $2 $3 2>&1 | tee -a "$LOGFILE"
exit $?
=#
using AtomsIO
using Dates
using DFTK
using JSON3
using JLD2
using Statistics
using DftFunctionals
using ForwardDiff
using Optim
DFTK.setup_threading()
include(joinpath(@__DIR__, "../../sol58lc/load.jl"))
include(joinpath(@__DIR__, "../../../lattice_relaxation.jl"))

data_path = joinpath(@__DIR__, "../../sol58lc/structures")

trainset = [
    "Si_diamond",
    # "GaAs_b3",
    "NaCl_b1",
    "Al_fcc",
    "V_bcc",
]
files_train = ["$x.extxyz" for x in trainset]
@assert issubset(files_train, readdir(data_path))

use_smearing = Dict(
    "Si"=>true,
    "NaCl"=>false,  # Large band gap -> unstable Fermi-level derivative under temperature
    "Al"=>true,
    "V"=>true,
)

function main(workdir, κ, μ)
    mkpath(workdir)
    save_path = joinpath(workdir, "kappa$(κ)_mu$(μ).json")

    # Make PBE Ansatz
    make_pbe(θ) = [PbeExchange(θ), DftFunctional(:gga_c_pbe)]

    # Load training set
    x_train = map(files_train) do file
        path = joinpath(data_path, file)
        system = load_system(path)
        convert_system_to_case(system)
    end

    # Storage for extracting predictions from loss as a side effect
    trajectory = []

    function loss_fn(θ::AbstractVector{T}, cases) where T
        # Side effect: Collect state
        state = Dict()
        state["theta"] = ForwardDiff.value.(θ)
        state["preds"] = Dict()
        push!(trajectory, state)
        open(save_path, "w") do io
            JSON3.pretty(io, trajectory)
        end

        loss = zero(T)
        batch_size = length(cases)
        for case in cases  # TODO: might parallelize if memory allows
            (; name, Ecut, kgrid, a0_exp) = case()
            
            temperature = 0
            smearing = Smearing.None()
            if use_smearing[name]
                temperature = 1e-3
                smearing = Smearing.Gaussian()
            end

            # Initialize lattice, either from PBE or from previous iteration
            p0 = case().p
            if length(trajectory) > 1
                p0.a = trajectory[end-1]["preds"][name]
            end

            @info now() name Ecut kgrid p0 temperature smearing
            a_opt = @time "$name" optimal_lattice(case, make_pbe, θ; Ecut, kgrid, p0, temperature, smearing).a

            # Side effect: save predictions
            state["preds"][name] = ForwardDiff.value(a_opt)
            open(save_path, "w") do io
                JSON3.pretty(io, trajectory)
            end

            loss += ((a_opt - a0_exp) /  a0_exp)^2
        end
        loss /= batch_size

        # Side effect: Save loss
        state["loss"] = ForwardDiff.value(loss)
        state["grad"] = ForwardDiff.partials(loss)
        state["time"] = now()
        open(save_path, "w") do io
            JSON3.pretty(io, trajectory)
        end
        
        loss
    end

    θ = ComponentVector(; κ, μ)
    loss_fn(θ, x_train)
end

if !isinteractive()
    if length(ARGS) != 3
        error("This script requires three command-line arguments: workdir and kappa and mu.")
    end
    workdir = ARGS[1]
    κ = parse(Float64, ARGS[2])
    μ = parse(Float64, ARGS[3])
    main(workdir, κ, μ)
end
