using Test
using FiniteDifferences
using LinearAlgebra
include("lib.jl")


@testset "total_energy and ForwardDiff"  begin
    # Load GTH PBE pseudopotential table
    pseudopotentials = PseudoFamily("cp2k.nc.sr.pbe.v0_1.largecore.gth")

    # Load initial parameters for Li-q1 (only LDA available for q1)
    Li_q1_init = ElementPsp(:Li, PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"))

    # Expose optimisable parameters
    make_li_q1(θ) = unflatten(ComponentVector(; θ..., rloc=Li_q1_init.psp.rloc))
    params = flatten(Li_q1_init)
    tunable_param_keys = [:cloc1, :cloc2, :rp1, :rp2, :h1, :h2]
    θ_init = params[tunable_param_keys]

    system      = load_system("Li-BCC.extxyz")
    smearing    = Smearing.Gaussian()
    temperature = 0.00225
    kgrid       = (2, 2, 2)
    Ecut        = 20
    tol         = 1e-8
    
    (; lattice, atoms, positions) = DFTK.parse_system(system, pseudopotentials)

    # Update pseudopotential on each Li atom
    atoms = map(atoms) do atom
        (element_symbol(atom) == :Li) ? Li_q1_init : atom
    end
    model0 = model_DFT(lattice, atoms, positions; functionals=PBE(), temperature, smearing)

    basis0 = PlaneWaveBasis(model0; Ecut, kgrid)
    scfres = self_consistent_field(basis0; tol)

    E1 = total_energy(basis0; tol)
    @test isapprox(E1, scfres.energies.total; atol=1e-12)

    function lithium_psp_total_energy(θ)
        T = eltype(θ)
        Li_q1 = make_li_q1(θ)
        atoms = map(model0.atoms) do atom
            (element_symbol(atom) == :Li) ? Li_q1 : atom
        end
        model = Model(model0; lattice=T.(model0.lattice), atoms)
        basis = PlaneWaveBasis(model; Ecut, kgrid)
        total_energy(basis; tol)
    end

    E2 = lithium_psp_total_energy(θ_init)
    @test isapprox(E2, scfres.energies.total; atol=1e-12)

    E3, grad3 = value_and_gradient(lithium_psp_total_energy, AutoForwardDiff(), θ_init)

    E4, grad4 = value_and_gradient(lithium_psp_total_energy, 
        AutoFiniteDifferences(; fdm=central_fdm(5, 1, max_range=1e-3)), 
        θ_init)
    
    @show grad3
    @show grad4
    @show grad3 - grad4
    @test isapprox(grad3, grad4, atol=1e-7, norm=x -> norm(x, Inf))
end

