#!/bin/bash

# Define the ranges for kappa and mu
kappa_range=$(seq 0.5 0.1 1.5)
mu_range=$(seq 0.15 0.01 0.25)

# Directory to store the results
TS=$(date +"%Y%m%d_%H%M%S")
workdir="grid_search_results_${TS}"
mkdir -p "$workdir"

# Maximum number of parallel jobs
MAX_JOBS=25

echo "Work directory: $workdir"
echo "Kappa range: $kappa_range"
echo "Mu range: $mu_range"
# Calculate product of lengths of kappa and mu ranges
kappa_count=$(echo $kappa_range | wc -w)
mu_count=$(echo $mu_range | wc -w)
total_jobs=$((kappa_count * mu_count))
echo "Total jobs: $total_jobs"
echo "Max concurrent jobs: $MAX_JOBS"
echo "-------------------------------------"

# Loop over kappa and mu values
# and launch jobs in the background
for kappa in $kappa_range; do
    for mu in $mu_range; do
        # Enforce max concurrent jobs
        while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
            sleep 1
        done

        LOGFILE="$workdir/kappa${kappa}_mu${mu}.log"
        
        # Launch job
        (
            echo "Launching $LOGFILE"
            julia --project -t1 loss.jl $workdir $kappa $mu > "$LOGFILE" 2>&1
        ) &
    done
done

# Wait for all remaining background jobs to finish
wait

echo "All grid search jobs completed. Results are in the '$workdir' directory."