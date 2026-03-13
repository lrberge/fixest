
# This makefile generates the data, runs all the estimations and gathers the results

ifeq ($(OS),Windows_NT)
  R:=Rscript.exe
else
  R:=Rscript
endif

# pick your preferred way to run python
# PYTHON:=python
PYTHON:=uv run

.PHONY: all clean all_benchmarks
all: all_benchmarks

all_benchmarks: results/all_results.csv

# Data generations
main_data: data/data_n-1e7_rep-3.gz
data/data_n-1e7_rep-3.gz: code/data.R
	$(R) code/data.R

# OLS benchmarks
results/results_ols_R.csv: main_data code/ols.R
	$(R) code/ols.R

results/results_ols_py.csv: main_data code/ols.py
	$(PYTHON) code/ols.py

results/results_ols_jl.csv: main_data code/ols.jl
	julia --project=. code/ols.jl

# Poisson benchmarks
results/results_poisson_R.csv: main_data code/poisson.R
	$(R) code/poisson.R

results/results_poisson_py.csv: main_data code/poisson.py
	$(PYTHON) code/poisson.py

results/results_poisson_jl.csv: main_data code/poisson.jl
	julia --project=. code/poisson.jl

# gathering the results
results/all_results.csv: results/results_ols_R.csv results/results_ols_py.csv results/results_ols_jl.csv results/results_poisson_R.csv results/results_poisson_py.csv results/results_poisson_jl.csv
	$(R) code/combine-results.R


#
# clean
#

clean:
	rm results/*
