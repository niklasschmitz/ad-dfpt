# julia --project 0_generate_structures.jl  # Already included in runs/

structure_list=(
    "Si.extxyz"
    "C.extxyz"
    "CsCl.extxyz"
)

for structure in "${structure_list[@]}"; do
    bash 1_relax.jl runs/${structure} 8  # Run volume relaxation with 8x8x8 k-grid
done

# SCF tolerances to test
tol_list=(
    1e-2
    1e-4
    1e-6
    1e-8
    1e-10
    1e-12
)

for structure in "${structure_list[@]}"; do
    for tol in "${tol_list[@]}"; do
        bash 2_elastic_dfpt.jl runs/${structure}_kgrid8_relax_scfres0.jld2 $tol
        bash 2_elastic_finitediff.jl runs/${structure}_kgrid8_relax_scfres0.jld2 $tol
    done
done

for structure in "${structure_list[@]}"; do
    julia --project 3_collect_summary.jl "runs/${structure}_kgrid8_relax_scfres0.jld2"
done
