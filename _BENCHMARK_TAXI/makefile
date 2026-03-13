# This makefile runs all the estimations and gathers the results

ifeq ($(OS),Windows_NT)
  R:=Rscript.exe
else
  R:=Rscript
endif

# pick your preferred way to run python
# PYTHON:=python
PYTHON:=uv run

.PHONY: all clean benchmark_taxi
all: benchmark_taxi

benchmark_taxi: results/all_results_taxi.csv

# Data download and formatting
data/nyc_taxi.parquet: code/data.R
	$(R) code/data.R

# Estimations
results/results_taxi_R.csv: code/taxi.R
	$(R) code/taxi.R

results/results_taxi_py.csv: code/taxi.py
	$(PYTHON) code/taxi.py

results/results_taxi_jl.csv: code/taxi.jl
	julia --project=. code/taxi.jl

# gathering the results
results/all_results_taxi.csv: results/results_taxi_R.csv results/results_taxi_py.csv results/results_taxi_jl.csv 
	$(R) code/combine-taxi-results.R


#
# clean
#

clean:
	rm results/*


