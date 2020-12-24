Replication code and data for the benchmarking of the package:
    ____________
    || fixest ||
    ************

***********
* SCRIPTS *
***********

    - src/10_data_generation.R: R code for data simulation
    - src/20_benchmakr_r.R: benchmarking code in R
    - src/30_benchmark.jl: Benchmarking code in Julia
    - src/40_benchmark_stata.do: benchmarking code in Stata
    - src/50_collecting_results.R: compilation of R and Stata results

***************
* REPLICATING *
***************

We ran these benchmarks on a 2019 15-inch MacBook Pro with a 2.6GHz 6-Core Intel Core i7 and 32GB of RAM.

R (version 4.0.3) and Julia (1.5.3) were installed using Homebrew (commit 2c77a540b). Stata (version 16.0 Revision 24 Jul 2019) was installed from the DMG.

See the `setup` folder for scripts to setup your environment.

*************
* DATA SETS *
*************

Main data sets used in replication:
    - [DATA/] base_10M.csv: simulated data set for 10M obs.
    - [DATA/] base_all_simulations.Rdata: array of all simulated data sets
    - [_STATA/] base_s1_r1.dta to base_square_s4_g2_r10.dta: same data sets as in base_all_simulations.Rdata but in Stata format

Timing results of the replication:
    - [DATA/] results_bench_R.txt: Results of the benchmarking of R functions
    - [DATA/] results_diff_bench_R.txt: Results of the benchmarking of fixest and lfe for the difficult data set
    - [DATA/] julia_bench_[nbFixef]FE.txt & julia_bench_diff.txt: results of the Julia benchmarking for nbFixef FEs, and for the difficult data set
    - [_STATA/] a2reg_G2.txt to xtpoisson_G3.txt ([funname]_G[nbFixef].txt): benchmarking results for each Stata function, for each number of fixed-effects (from 1 to 3)
    - [DATA/] results_all.txt & results_diff_all.txt: compilation of R, Stata and Julia results -- same for the difficult data set

**************************
* ADDITIONAL INFORMATION *
**************************

    - All packages for R, Stata and Julia were updated in February 2020.
    - You need to provide the right path to the data sets for the code to run


******************
* LOG OF CHANGES *
******************

    - [7 Fev 2020]:
        * The benchmark folder containing all replication files has been moved to GitHub
        * New "difficult" benchmark has been added for OLS

	- [18 Nov 2019]:
		* name harmonizaiton (dropping _rect suffix [coming from a former battery of benchmarks])
		* variable ln_y (=log(y+1)) now directly exists in the data sets (would otherwise lead to unfair advantage to Stata where the variable was computed before the estimation). Code updated accordingly
		* functions plm, a2reg, xtpoisson and poi2hdfe removed from the benchmark.
		* reordering of the models in Benchmark_Stata.do





















































