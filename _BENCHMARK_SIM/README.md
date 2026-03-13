
Welcome to the folder containing the benchmarks for simulated data of the `fixest` paper.

# Running the benchmarks

If you have `make` installed, you can run `make` and it should work (provided you have R/python/Julia installed, see below). Given that the estimation tools use multithreading, `make` should not be run in multithreaded mode or else the performance of each run would be artificially degraded.

If you do not have `make`, you would need to run these scripts in that order (only the order of the first and last matter) using the appropriate software for each language:
- `code/data.R`: generates the data in `./data/`
- `code/ols.R`: generates `results/results_ols_R.csv`
- `code/ols.jl`: generates `results/results_ols_jl.csv`
- `code/ols.py`: generates `results/results_ols_py.csv`
- `code/poisson.R`: generates `results/results_poisson_R.csv`
- `code/poisson.jl`: generates `results/results_poisson_jl.csv`
- `code/poisson.py`: generates `results/results_poisson_py.csv`
- `code/combine-results.R`: generates `results/all_results.csv`

The results of the benchmarks will be stored in `results/all_results.csv`.

# Installation of the software

### R

- install R from https://cran.r-project.org/

- you will need the following dependencies:
  - data.table
  - fixest
  - lfe
  - alpaca

- the code to install these dependencies is within the scripts and should run automatically

- you can also install them manually with:
  - `R> install.packages(c("data.table", "fixest", "lfe", "alpaca"))`

- IMPORTANT NOTE for MacOS users. The software `fixest` uses multithreading via OpenMP. OpenMP may not be activated on MacOS, please see this discussion to resolve the problem https://github.com/lrberge/fixest/issues/63

### Python

We strongly advise to use the `uv` package manager.

- install `uv` from https://docs.astral.sh/uv/getting-started/installation/

- to sync the project and download the required dependencies (including python), run this command within the current folder:
  - `> uv sync`

- if it does not work, you can install the depdendencies manually:
  - `> uv add pandas pyfixest`

- IMPORTANT NOTE: as of 2026-01-15 you need `python<=3.13` for `pyfixest` to work, this is due to a depdency issue with `scipy` (which may have issues to be installed on `python>=3.14`). This should be automatically taken care of with `uv`.

### Julia

- install Julia from https://github.com/JuliaLang/juliaup

- to download the dependencies, in this folder run Julia, then:
  - type `]` to access the package manager within Julia
  - Run the following code:
    - `pkg> activate .`
    - `pkg> instantiate`

- to install the dependencies manually, type in Julia's package manager:
  - `pkg> add CSV DataFrames FixedEffectModels GLFixedEffectModels GLM Distributions`
