using AtomsIO
using Dates
using DFTK
using PseudoPotentialData
using Unitful
using UnitfulAtomic
using JSON3
using JLD2
using Statistics
using DifferentiationInterface
using ForwardDiff
using ComponentArrays
using Optim
using LineSearches


"""Helper function for ForwardDiff for 1st-order Hellman-Feynman derivatives."""
function total_energy(basis::PlaneWaveBasis; kwargs...)
    scfres = self_consistent_field(basis; kwargs...)
    scfres.energies.total
end

"""ForwardDiff rule for 1st-order Hellman-Feynman total energy derivatives."""
function total_energy(basis::PlaneWaveBasis{T}; kwargs...) where {T <: ForwardDiff.Dual}
    basis_primal = DFTK.construct_value(basis)
    scfres = self_consistent_field(basis_primal; kwargs...)
    ρ = compute_density(basis, scfres.ψ, scfres.occupation)
    (; energies) = DFTK.energy(basis, scfres.ψ, scfres.occupation; ρ, scfres.eigenvalues, scfres.εF)
    energies.total
end

function flatten(Li::ElementPsp)
    params = ComponentArray(
        rloc=Li.psp.rloc,
        cloc1=Li.psp.cloc[1],
        cloc2=Li.psp.cloc[2],
        rp1=Li.psp.rp[1],
        rp2=Li.psp.rp[2],
        h1=Li.psp.h[1][1],
        h2=Li.psp.h[2][1],
    )
    params
end

function unflatten(params)
    (; rloc, cloc1, cloc2, rp1, rp2, h1, h2) = params
    T = eltype(params)
    psp = PspHgh{T}(
        1, rloc, [cloc1, cloc2, 0., 0.], 1, [rp1, rp2],
        [[h1;;], [h2;;]], "", "Li GTH-q1",
    )
    ElementPsp(:Li, psp)
end

function train(workdir="."; datadir=@__DIR__, num_train_steps=50, λ=0.0)
    mkpath(workdir)
    save_path = joinpath(workdir, "trajectory.json")

    # Load GTH PBE pseudopotential table
    pseudopotentials = PseudoFamily("cp2k.nc.sr.pbe.v0_1.semicore.gth")

    # Load initial parameters for Li-q1 (only LDA available for q1)
    Li_q1_init = ElementPsp(:Li, PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"))

    # Expose optimisable parameters
    make_li_q1(θ) = unflatten(ComponentVector(; θ..., rloc=Li_q1_init.psp.rloc))
    params = flatten(Li_q1_init)
    tunable_param_keys = [:cloc1, :cloc2, :rp1, :rp2, :h1, :h2]
    θ_init = params[tunable_param_keys]

    temperature = 0.00225  # 0.0045 Ry
    # smearing = Smearing.FermiDirac()  # TODO Hunt down NaN in Entropy for ForwardDiff
    smearing = Smearing.Gaussian()    
    scf_tol = 1e-7

    # Load training set and attach relative weights
    systems = [
        (; weight=λ,     file="Li-BCC.extxyz", Ecut= 20, kgrid=(8,8,8)),
        (; weight=(1-λ), file="Li-XO.extxyz",  Ecut=120, kgrid=(8,8,8)),
    ]

    volume_scaling_range = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]
    cases_train = map(systems) do (; weight, file, Ecut, kgrid)
        system = load_system(joinpath(datadir, file))
        name = system.system_data.name

        # # Load reference total energies from all-electron data
        # data = JSON3.read(joinpath(datadir, "all_electron_data/wien2k_$(name).json"))
        # @assert isapprox(
        #     data[:volumes_Ang3] ./ data[:volumes_Ang3][4],
        #     volume_scaling_range,
        #     atol=1e-6
        # ) "Inconsistent volumes found for $name."
        # # Convert total energies to atomic units and per-atom
        # num_atoms = length(system)
        # reference_energies = austrip.(data[:energies_eV] * u"eV") / num_atoms

        # Load reference total energies from Li-q3 results on same kgrid
        # (Note that the Li-q3 results are already energies-per-atom)
        data = JSON3.read(joinpath(datadir, "convergence_study_Li_q3/$(name)_Ecut120_kgrid8.json"))
        reference_energies = Vector(data[:energies])

        (; weight, system, Ecut, kgrid, reference_energies)
    end

    # Storage for extracting predictions from loss as a side effect
    trajectory = []

    # Storage for passing previous SCF results as initial guesses
    scf_results = Dict(system.system_data.name=>Dict() for (; system) in cases_train)

    function loss_and_gradient!(L, G, θ::AbstractVector)
        start_time = now()
        @info start_time θ

        # Side effect: Collect state
        state = Dict()
        state["theta"] = deepcopy(θ)
        state["cases"] = Dict()
        push!(trajectory, state)
        open(save_path, "w") do io
            JSON3.pretty(io, trajectory)
        end

        function loss_fn(θ, case)
            T = eltype(θ)

            (; system, Ecut, kgrid, reference_energies) = case
            name = system.system_data.name
            (; lattice, atoms, positions) = DFTK.parse_system(system, pseudopotentials)

            @info now() name Ecut kgrid

            # Update pseudopotential on each Li atom
            Li_q1 = make_li_q1(θ)
            atoms = map(atoms) do atom
                (element_symbol(atom) == :Li) ? Li_q1 : atom
            end
            model0 = model_DFT(T.(lattice), atoms, positions; functionals=PBE(), temperature, smearing)

            energies = map(volume_scaling_range) do scale
                model = Model(model0; lattice=cbrt(scale) * model0.lattice)
                basis = PlaneWaveBasis(model; Ecut, kgrid)
                ρ = guess_density(DFTK.construct_value(basis))
                ψ = nothing
                
                # Load SCF results from previous iteration as initial guesses
                if haskey(scf_results[name], scale)
                    scfres = scf_results[name][scale]
                    ρ = scfres.ρ
                    # ψ = scfres.ψ  # TODO re-enable psi initial guess
                end

                # Save new scfres (in the primal) by a callback
                function save_callback(scfres)
                    if scfres.stage == :finalize
                        scf_path = joinpath(workdir, "scfres_$(name)_v$(scale).jld2")
                        @time "I/O" save_scfres(scf_path, scfres; save_ψ=false)
                        scf_results[name][scale] = scfres
                    end
                    scfres
                end
                callback = ScfDefaultCallback(; show_damping=false) ∘ save_callback

                # Run SCF (with good initial guesses) and compute total energy
                energy = total_energy(basis; tol=scf_tol, ρ, ψ, callback)
                energy_per_atom = energy / length(model.atoms)
                energy_per_atom
            end

            # TODO can try other loss functions here (e.g. Δ or ν metric)
            # Compute per-system loss
            @assert length(energies) == 7
            @assert length(reference_energies) == 7
            Ea = energies .- energies[4]
            Eb = reference_energies .- reference_energies[4]
            loss = sqrt(
                sum(abs2, Ea - Eb) / (
                    sqrt(sum(abs2, Ea .- mean(Ea)) * sum(abs2, Eb .- mean(Eb)))
                )
            )

            @info now() name ForwardDiff.value.(Ea) Eb ForwardDiff.value(loss)
            state["cases"][name] = Dict()
            state["cases"][name]["energies"] = ForwardDiff.value.(Ea)
            state["cases"][name]["energies_grads"] = ForwardDiff.partials.(Ea)
            state["cases"][name]["loss"] = ForwardDiff.value(loss)
            state["cases"][name]["grad"] = ForwardDiff.partials(loss)

            loss
        end

        L, grad = value_and_gradient(AutoForwardDiff(), θ) do θ
            weighted_losses = [
                case.weight * loss_fn(θ, case) for case in cases_train
                if case.weight > 0
            ]
            sum(weighted_losses)
        end

        end_time = now()
        state["loss"] = L
        state["grad"] = deepcopy(grad)
        state["date"] = end_time
        state["time"] = string(canonicalize(end_time - start_time))
        open(save_path, "w") do io
            JSON3.pretty(io, trajectory)
        end

        if !isnothing(G)
            G .= grad
        end

        L
    end

    res = Optim.optimize(
        Optim.only_fg!(loss_and_gradient!),
        θ_init,
        BFGS(linesearch=LineSearches.BackTracking(order=3, maxstep=1e-1)),  # TODO maxstep?
        Optim.Options(;
            allow_f_increases=true,
            show_trace=true,
            extended_trace=true,
            iterations=num_train_steps,
        )
    )

    println(res)
end
