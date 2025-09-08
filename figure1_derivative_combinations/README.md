## Basic: Plots
To reproduce plots from JSON files, see `analysis.ipynb`.

## Advanced: Re-run DFT scripts

1. To re-run all DFT calculations (<1 hour in total on a decent workstation):
```
bash silicon_Ecut30_kgrid4_displacefalse.jl  # used in the band structure plot
bash silicon_Ecut30_kgrid4_displacetrue.jl   # used in the density and forces plot
```
2. To extract JSON3 summaries out of the results:
```
bash extract_json.jl
```
3. Plot results with `analysis.ipynb`

