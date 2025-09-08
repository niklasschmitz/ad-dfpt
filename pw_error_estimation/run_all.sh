#!/bin/bash

# Define the structures to run
structures="LiF_b1 LiCl_b1 NaF_b1 NaCl_b1 BN_b3 MgO_b1 C_diamond CaO_b1 AlN_b3 LiH_b1 MgS_b1 AlP_b3 GaN_b3 GaP_b3 AlAs_b3 SiC_b3 BP_b3 BAs_b3 InAs_b3 Si_diamond InP_b3"

# Maximum number of parallel jobs
MAX_JOBS=5

echo "Structures: $structures"
structure_count=$(echo $structures | wc -w)
echo "Total jobs: $structure_count"
echo "Max concurrent jobs: $MAX_JOBS"
echo "-------------------------------------"

mkdir -p runs

# Loop over structures and launch jobs in the background
for structure in $structures; do
    # Enforce max concurrent jobs
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 1
    done
    
    # Launch job
    (
        echo "Launching $structure"
        bash 0_run.jl $structure.extxyz
    ) &
done

# Wait for all remaining background jobs to finish
wait

echo "All plane-wave error estimation jobs completed."
