#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t4 $BN.jl 2>&1 | tee $BN.log
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
using LineSearches
DFTK.setup_threading()
include(joinpath(@__DIR__, "../../sol58lc/load.jl"))
include(joinpath(@__DIR__, "../../../lattice_relaxation.jl"))

data_path = joinpath(@__DIR__, "../../sol58lc/structures")

trainset = [
    "Si_diamond",
    "NaCl_b1",
    "Al_fcc",
    "V_bcc",
]
files_train = ["$x.extxyz" for x in trainset]
@assert issubset(files_train, readdir(data_path))

function train(workdir=".", num_train_steps=20)
    mkpath(workdir)
    save_path = joinpath(workdir, "trajectory.json")

    # Load initial parameters
    make_pbe(θ) = [PbeExchange(θ), DftFunctional(:gga_c_pbe)]
    θ_init = parameters(DftFunctional(:gga_x_pbe))

    # Load previous trajectory if it exists
    if isfile(save_path)
        trajectory = JSON3.read(save_path)
        if !isempty(trajectory) && haskey(trajectory[end], "theta")
            θ_loaded = trajectory[end]["theta"]
            θ_init = ComponentVector(; θ_init..., κ=θ_loaded[1], μ=θ_loaded[2])
            println("Found previous trajectory, loading θ_init from it")
            @show θ_init
        end
    end

    # Load training set
    x_train = map(files_train) do file
        path = joinpath(data_path, file)
        system = load_system(path)
        convert_system_to_case(system)
    end

    # Use smearing with a small temperature
    temperature = 1e-3
    smearing = Smearing.Gaussian()

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
        @info now() θ

        loss = zero(T)
        batch_size = length(cases)
        for case in cases  # TODO: might parallelize if memory allows
            (; name, Ecut, kgrid, a0_exp) = case()

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

    res = Optim.optimize(
        θ -> loss_fn(θ, x_train),
        θ_init,
        BFGS(linesearch=LineSearches.BackTracking(order=3)),
        autodiff = :forward,
        Optim.Options(;
            allow_f_increases=true,
            show_trace=true,
            extended_trace=true,
            iterations=num_train_steps,
        )
    )

    println(res)
end

train()
