using ForwardDiff
using DFTK
using DifferentiationInterface
include("optimise.jl")


function optimal_lattice(case, functional, θ;
                         Ecut, kgrid, temperature=1e-3, smearing=Smearing.Gaussian(),
                         kinetic_blowup=BlowupCHV(), magnetic_moments=[],
                         g_abstol=1e-6, f_abstol=1e-8, p0=case().p)    
    function dftk_fg!(E, G, p, θ)
        T = promote_type(eltype(p), eltype(θ))
        (; lattice, atoms, positions) = case(p)
        model = model_DFT(Matrix{T}(lattice), atoms, positions;
                          functionals=functional(θ), temperature, smearing,
                          kinetic_blowup, magnetic_moments)
        basis = PlaneWaveBasis(model; Ecut, kgrid)
        ρ0 = guess_density(basis, magnetic_moments)
        ρ0 = ForwardDiff.value.(ρ0)  # TODO Fix this upstream in DFTK to not have to do this here
        scfres = self_consistent_field(basis; tol=g_abstol/10, ρ=ρ0)
        if !isnothing(G)
            G .= ForwardDiff.gradient(p) do p
                # Use Hellman-Feynman explicitly to compute dE / dp
                (; lattice, atoms, positions) = case(p)
                model_p = Model(model; lattice, atoms, positions)
                basis_p = PlaneWaveBasis(model_p; Ecut, kgrid)
                ρ = DFTK.compute_density(basis_p, scfres.ψ, scfres.occupation)
                (; energies) = DFTK.energy(basis_p, scfres.ψ, scfres.occupation;
                                           ρ, scfres.eigenvalues, scfres.εF)
                energies.total
            end
        end
        E = scfres.energies.total
    end
    optimise(dftk_fg!, p0, θ; f_abstol, x_abstol=-1, g_abstol, show_trace=true, extended_trace=true)
end


function lattice_sensitivity(case, functional, θ; kwargs...)
    function compute_lattice_constant(θ)
        optimal_lattice(case, functional, θ; kwargs...).a
    end
    value_and_gradient(compute_lattice_constant, AutoForwardDiff(), θ)
end
