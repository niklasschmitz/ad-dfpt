
# Algorithmic differentiation for plane-wave DFT: materials design, error control and learning model parameters

Supporting information containing structures, raw data and computational scripts for the paper:

> Niklas Frederik Schmitz, Bruno Ploumhans, Michael F. Herbst.  
> *Algorithmic differentiation for plane-wave DFT:
materials design, error control and learning model parameters*  
> Preprint on arxiv: **TODO add link**

The code in this repository has been used to run all calculations and produce all plots of the above paper. It relies on [DFTK](dftk.org) version 0.7.16.

## Running the code and reproducing the plots

Running the code requires an installation of [Julia 1.11.5](https://www.julialang.org/downloads/). To install all further necessary dependencies, run the following command in this directory:
```sh
julia --project=@. -e "import Pkg; Pkg.instantiate()"  # Install dependencies
```

For each example, we provide the code to reproduce all results and plots. Raw computational data is also included, so you can regenerate all plots without re-running the DFT calculations. For full reproducibility, all computations can be run from scratch, with running times ranging from seconds to several hours on a modern workstation. For reference, we also include log files from our original runs. Note that performance could be further improved by tweaking threading settings.


The examples sections from the paper map as follows to this repository:
1. **Figure 1: Systematic DFT derivatives**: See [figure1_derivative_combinations/](figure1_derivative_combinations/)
2. **Elasticity**: See [elastic_constants/](elastic_constants/)
   - For the plot, see [analysis.ipynb](elastic_constants/analysis.ipynb)
   - To regenerate the convergence study data, run `run.sh`
3. **Inverse materials design**: See [inverse_design/](inverse_design/)
   - For the plot, see [analysis.ipynb](inverse_design/analysis.ipynb).
   - To rerun the optimization, run `bash run.jl`
4. **Learning the exchange-correlation functional**: See [xc_finetuning/](xc_finetuning/).
   - For the final plot, see [eval/pbe_Si_Al_V_NaCl_trained/analysis.ipynb](xc_finetuning/eval/pbe_Si_Al_V_NaCl_trained/analysis.ipynb)
   - For the data preprocessing, see [sol58lc](xc_finetuning/sol58lc/)
   - For the training, see [train/pbe_Si_Al_V_NaCl/train.jl](xc_finetuning/train/pbe_Si_Al_V_NaCl/train.jl)
   - For the evaluation, see [eval/](xc_finetuning/eval/)
5. **Property-driven pseudopotential optimization**: See [pseudopotential_tuning/](pseudopotential_tuning/).
   - For the final plot, see [analysis.ipynb](pseudopotential_tuning/analysis.ipynb)
   - For the convergence study for GTH Li-q3 target results, see [convergence_study_Li_q3](pseudopotential_tuning/convergence_study_Li_q3/)
   - For the training code of Li-q1, see [lib.jl](pseudopotential_tuning/lib.jl)
6. **Propagating XC functional uncertainty**: See [forward_uncertainty/](forward_uncertainty/).
   - For the final plot, see [analysis.ipynb](forward_uncertainty/analysis.ipynb).
   - For rerunning all computations, see [run.sh](forward_uncertainty/run.sh)  
7. **Estimation of the plane-wave basis error**: See [pw_error_estimation](pw_error_estimation/).
   - For the final plot, see [2_analysis.ipynb](pw_error_estimation/2_analysis.ipynb)
   - For rerunning all computations, see [run_all.sh](pw_error_estimation/run_all.sh).

