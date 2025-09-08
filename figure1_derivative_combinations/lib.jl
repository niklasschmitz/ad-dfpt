using AtomsBuilder
using Brillouin
using ComponentArrays
using Dates
using DftFunctionals
using DFTK
using ForwardDiff
using Functors
using JLD2
using LinearAlgebra
using PseudoPotentialData
using StaticArrays
using Unitful
using UnitfulAtomic
include(joinpath(@__DIR__, "compute_bands_forwarddiff.jl"))


"""Helper for injecting differentiable params into PspHgh for silicon."""
function construct_psp(cloc1::T) where {T}
    PspHgh{T}(4, 0.44, SVector(cloc1, 0.0, 0.0, 0.0), 1, [0.43563383, 0.49794218], [[8.9517415 -2.70627082; -2.70627082 3.4937806], [2.43127673;;]], "hgh/pbe/si-q4", "Si GTH-PBE-q4")
end
function destruct_psp(psp::PspHgh)
    psp.cloc[1]
end


function make_model(θ; symmetries)
    T = eltype(θ)

    # Construct XC from parameters
    params_x = parameters(DftFunctional(:gga_x_pbe))
    functionals = [
        PbeExchange(ComponentVector(; params_x..., θ.κ)),
        DftFunctional(:gga_c_pbe)
    ]

    # Construct pseudopotential from parameters
    Si_psp = construct_psp(θ.cloc1)
    pseudopotentials = fill(Si_psp, 8)

    # Bulk silicon in the cubic cell
    silicon = model_DFT(bulk(:Si; cubic=true); functionals, pseudopotentials)

    # Scale the lattice and displace positions
    lattice = (1 + θ.strain) * T.(silicon.lattice)
    positions = silicon.positions + θ.displacement
    Model(silicon; lattice, positions, symmetries)
end

function compute_quantities(θ; Ecut, kgrid, kgrid_bs, tol=1e-10, symmetries=true, return_scfres=false)
    model = make_model(θ; symmetries)
    basis = PlaneWaveBasis(model; Ecut, kgrid)
    scfres = self_consistent_field(basis; tol)

    # Compute quantities of interest
    stress = compute_stresses_cart(scfres)
    forces = compute_forces_cart(scfres)
    band_data = compute_bands_forwarddiff(scfres.basis, kgrid_bs; scfres.ρ, tol=1e-6)
    y = ComponentArray(;
        energies=collect(values(scfres.energies)),
        forces,
        stress,
        scfres.ρ,
        scfres.eigenvalues,
        scfres.εF,
        occupation=scfres.occupation,
        eigenvalues_bs=band_data.eigenvalues,
    )

    if return_scfres
        return (; y, scfres, band_data)
    end

    y
end

function run(prefix; Ecut=30, nkpt=10, displace=false)
    kgrid = (nkpt, nkpt, nkpt)

    # Uniform shift to better visualize atoms in x-y slice later
    displacement = fill([0.25, 0.25, 0.0], 8)

    kinter = let
        model = model_atomic(bulk(:Si; cubic=true); pseudopotentials=fill(nothing, 8))
        model = Model(model; positions=model.positions + displacement)
        kpath = irrfbz_path(model)
        kpath = KPath(kpath.points, [[:K, :Γ, :L]], kpath.basis, kpath.setting)  # Pick smaller KPath
        Brillouin.interpolate(kpath, density=austrip(40u"bohr"))
    end
    kgrid_bs = ExplicitKpoints(kinter)

    if displace
        displacement[7] += [0.1, 0.02, 0.0]
        symmetries = false
    else
        symmetries = true
    end

    psp = load_psp(PseudoFamily("cp2k.nc.sr.pbe.v0_1.largecore.gth")[:Si])

    θ = ComponentVector(
        κ=parameters(DftFunctional(:gga_x_pbe)).κ,
        cloc1=destruct_psp(psp),
        strain=0.0,
        displacement=displacement,
    )

    @show now() Ecut nkpt θ

    @show now() "Run basic SCF without derivatives"
    (; y, scfres, band_data) = compute_quantities(θ; Ecut, kgrid, kgrid_bs, return_scfres=true)
    save_scfres(
        prefix * "_scfres.jld2",
        scfres;
        save_ψ=false,
        extra_data=Dict(
            "y" => y,
            "θ" => θ,
            "bands" => DFTK.data_for_plotting(merge(band_data, (; kinter))),
        )
    )

    @show now() "Run Xc derivative (no symmetry-breaking)"
    jac_xc = ForwardDiff.derivative(θ.κ) do κ
        compute_quantities(ComponentVector(; θ..., κ); Ecut, kgrid, kgrid_bs, symmetries)
    end
    jldsave(prefix * "_derivatives_xc_kappa.jld2"; jac=jac_xc, θ)

    @show now() "Run Psp derivative (no symmetry-breaking)"
    jac_psp = ForwardDiff.derivative(θ.cloc1) do cloc1
        compute_quantities(ComponentVector(; θ..., cloc1); Ecut, kgrid, kgrid_bs, symmetries)
    end
    jldsave(prefix * "_derivatives_psp_cloc1.jld2"; jac=jac_psp, θ)

    @show now() "Run isotropic strain derivative (no symmetry-breaking)"
    jac_strain = ForwardDiff.derivative(θ.strain) do strain
        compute_quantities(ComponentVector(; θ..., strain); Ecut, kgrid, kgrid_bs, symmetries)
    end
    jldsave(prefix * "_derivatives_strain.jld2"; jac=jac_strain, θ)

    @show now() "Finished."
end


