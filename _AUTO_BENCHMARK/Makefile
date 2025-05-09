.PHONY: all

DATA_DIR=DATA
JULIA=julia
R=Rscript --vanilla
STATA=/Applications/Stata/StataSE.app/Contents/MacOS/StataSE -b do

# Data generation

$(DATA_DIR)/base_all_simulations.Rdata: src/10_data_generation.R
	$(R) src/10_data_generation.R

$(DATA_DIR)/base_10M.csv: $(DATA_DIR)/base_all_simulations.Rdata
$(DATA_DIR)/base_all_diff.Rdata: $(DATA_DIR)/base_all_simulations.Rdata

# R Benchmarks

$(DATA_DIR)/results_bench_R.txt: $(DATA_DIR)/base_all_simulations.Rdata $(DATA_DIR)/base_10M.csv src/20_benchmark_r.R
	$(R) src/20_benchmark_r.R

$(DATA_DIR)/results_diff_bench_R.txt: $(DATA_DIR)/results_bench_R.txt

# Julia Benchmarks

$(DATA_DIR)/julia_bench_1FE.txt: $(DATA_DIR)/base_10M.csv base_all_diff.Rdata src/30_benchmark.jl
	$(JULIA) src/30_benchmark.jl

$(DATA_DIR)/julia_bench_2FE.txt: $(DATA_DIR)/julia_bench_1FE.txt
$(DATA_DIR)/julia_bench_3FE.txt: $(DATA_DIR)/julia_bench_2FE.txt
$(DATA_DIR)/julia_bench_diff.txt: $(DATA_DIR)/julia_bench_3FE.txt

# Stata Results (n.b. need to enumerate still)

_STATA/reghdfe_diff.txt: $(DATA_DIR)/base_10M.csv src/40_benchmark_stata.do
	$(STATA) src/40_benchmark_stata.do


rstats: $(DATA_DIR)/results_bench_R.txt $(DATA_DIR)/results_diff_bench_R.txt
julia: $(DATA_DIR)/julia_bench_diff.txt $(DATA_DIR)/julia_bench_3FE.txt $(DATA_DIR)/julia_bench_2FE.txt $(DATA_DIR)/julia_bench_1FE.txt
stata: _STATA/reghdfe_diff.txt
all: rstats julia stata

clean:
	rm -f \
		$(DATA_DIR)/base_10M.csv \
		$(DATA_DIR)/base_all_diff.Rdata \
		$(DATA_DIR)/base_all_simulations.Rdata \
		$(DATA_DIR)/julia_bench_1FE.txt \
		$(DATA_DIR)/julia_bench_2FE.txt \
		$(DATA_DIR)/julia_bench_3FE.txt \
		$(DATA_DIR)/julia_bench_diff.txt \
		$(DATA_DIR)/results_bench_R.txt \
		$(DATA_DIR)/results_diff_bench_R.txt \
		_STATA/reghdfe_diff.txt