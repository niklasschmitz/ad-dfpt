
# Algorithmic differentiation for plane-wave DFT: materials design, error control and learning model parameters

Supporting information containing structures, raw data and computational scripts for the paper:

> Niklas Frederik Schmitz, Bruno Ploumhans, Michael F. Herbst.  
> *Algorithmic differentiation for plane-wave DFT:
materials design, error control and learning model parameters*  
> Preprint on arxiv: https://arxiv.org/abs/2509.07785

The code in this repository has been used to run all calculations and produce all plots of the above paper. It relies on [DFTK version 0.7.16](https://github.com/JuliaMolSim/DFTK.jl/releases/tag/v0.7.16).

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

## License
All code in this repository is licensed under the MIT license.

> MIT License
>
> Copyright (c) 2025 Niklas F. Schmitz and collaborators
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

All data in this repository, except where noted below, is licensed under [the **Creative Commons Attribution 4.0 (CC-BY 4.0)** license](https://creativecommons.org/licenses/by/4.0/legalcode.en).

**Exceptions**:
- The Sol58LC dataset in `xc_finetuning/sol58lc`. See [the corresponding README](xc_finetuning/sol58lc/README.md).
- The all-electron data in `pseudopotential_tuning/all_electron_data`. See [the corresponding README](pseudopotential_tuning/all_electron_data/README.md).

## Funding

This research was supported by
the Swiss National Science Foundation (SNSF, Grant No. 221186)
as well as the NCCR MARVEL, a National Centre of Competence in Research,
funded by the SNSF (Grant No. 205602).
