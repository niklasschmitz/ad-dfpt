# KGRID=4  # For debugging
KGRID=14  # As in XC finetuning experiment (ensures kspacing \leq 0.15 Å⁻¹)
bash 1_relax.jl Si_diamond.extxyz $KGRID
bash 2_linear_pushforward.jl Si_diamond.extxyz_kgrid${KGRID}_relax_scfres0.jld2
bash 2_ensemble.jl           Si_diamond.extxyz_kgrid${KGRID}_relax_scfres0.jld2
