using ComponentArrays
using DFTK
using PseudoPotentialData
using Unitful, UnitfulAtomic


function convert_system_to_case(
    system; 
    pseudopotentials=PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf"),
    kspacing=0.15u"Å^-1",
)
    name = system.system_data.name

    # Extract lattice constant and convert from Angstrom to atomic units
    a0_pbe = austrip(system.system_data.a0_pbe * u"Å")
    a0_exp = austrip(system.system_data.a0_exp * u"Å")
    
    parsed = DFTK.parse_system(system, pseudopotentials)
    atoms = parsed.atoms
    positions = parsed.positions

    # Generate k-point grid from PBE relaxed lattice and keep fixed
    kgrid = kgrid_from_maximal_spacing(parsed.lattice, kspacing)

    # Determine Ecut from pseudopotential
    Ecuts_atoms = [
        recommended_cutoff(pseudopotentials, Symbol(at.species)).Ecut
        for at in unique(atoms)
    ]
    Ecut = maximum(Ecuts_atoms)

    # Normalize lattice by dividing out the lattice constant. We do this so we
    # can later below multiply the lattice by a chosen trial lattice constant
    # during optimization.
    # TODO: This approach is a bit hacky, could be nicer if AtomsBuilder or AtomsBase
    #       had helper functions to rescale / strain a system.
    normalized_lattice = parsed.lattice / a0_pbe

    function case(p=ComponentArray(; a=a0_pbe))
        lattice = p.a * normalized_lattice
        (; name, lattice, atoms, positions, p, Ecut, kgrid, a0_pbe, a0_exp)
    end
end
