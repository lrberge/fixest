#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Oct 18 16:36:41 2019
# ~: Benchmark: collecting results
#----------------------------------------------#

library(data.table)
library(here)

####
#### Main results ####
####

RESULTS_DIR <- "results"

results_R <- fread(here(RESULTS_DIR, "results_bench_R.txt"))

# Poisson
results_poisson <- results_R[family == "poisson", -1]

# We add the stata results
ppmlhdfe_G1_raw <- fread(here("_STATA", "ppmlhdfe_G1.txt"))
ppmlhdfe_G2_raw <- fread(here("_STATA", "ppmlhdfe_G2.txt"))
ppmlhdfe_G3_raw <- fread(here("_STATA", "ppmlhdfe_G3.txt"))

ppmlhdfe_G1 <- ppmlhdfe_G1_raw[-1, .(method = "ppmlhdfe", n_obs = (10**(3:6))[c1], G = 1, rep = 1:10, time = c2)]
ppmlhdfe_G2 <- ppmlhdfe_G2_raw[-1, .(method = "ppmlhdfe", n_obs = (10**(3:6))[c1], G = 2, rep = 1:10, time = c2)]
ppmlhdfe_G3 <- ppmlhdfe_G3_raw[-1, .(method = "ppmlhdfe", n_obs = (10**(3:6))[c1], G = 3, rep = 1:10, time = c2)]

results_poisson_all <- rbindlist(list(results_poisson, ppmlhdfe_G1, ppmlhdfe_G2, ppmlhdfe_G3))
results_poisson_all$model <- "Poisson"

#
# Gaussian
#

results_gaussian <- results_R[family == "gaussian", -1]

# We add the stata results
reghdfe_G1_raw <- fread(here("_STATA", "reghdfe_G1.txt"))
reghdfe_G2_raw <- fread(here("_STATA", "reghdfe_G2.txt"))
reghdfe_G3_raw <- fread(here("_STATA", "reghdfe_G3.txt"))

reghdfe_G1 <- reghdfe_G1_raw[-1, .(method = "reghdfe", n_obs = (10**(3:7))[c1], G = 1, rep = 1:10, time = c2)]
reghdfe_G2 <- reghdfe_G2_raw[-1, .(method = "reghdfe", n_obs = (10**(3:7))[c1], G = 2, rep = 1:10, time = c2)]
reghdfe_G3 <- reghdfe_G3_raw[-1, .(method = "reghdfe", n_obs = (10**(3:7))[c1], G = 3, rep = 1:10, time = c2)]

# The Julia results
julia_G1_raw <- fread(here(RESULTS_DIR, "julia_bench_1FE.txt"))
julia_G2_raw <- fread(here(RESULTS_DIR, "julia_bench_2FE.txt"))
julia_G3_raw <- fread(here(RESULTS_DIR, "julia_bench_3FE.txt"))

n_obs <- c(rep(10**(3:6), 10), rep(1e7, 10))
julia_G1 <- julia_G1_raw[, .(method = "FixedEffectModels", n_obs = n_obs, G = 1, rep = 1:10, time = V1)]
julia_G2 <- julia_G2_raw[, .(method = "FixedEffectModels", n_obs = n_obs, G = 2, rep = 1:10, time = V1)]
julia_G3 <- julia_G3_raw[, .(method = "FixedEffectModels", n_obs = n_obs, G = 3, rep = 1:10, time = V1)]

# merging
results_gaussian_all <- rbindlist(list(results_gaussian, reghdfe_G1, reghdfe_G2, reghdfe_G3, julia_G1, julia_G2, julia_G3))
results_gaussian_all$model <- "Gaussian"

#
# NEGBIN
#

results_negbin <- results_R[family == "negbin", -1]

# We add the stata results
nbreg_G1_raw <- fread(here("_STATA", "nbreg_G1.txt"))
nbreg_G2_raw <- fread(here("_STATA", "nbreg_G2.txt"))
nbreg_G3_raw <- fread(here("_STATA", "nbreg_G3.txt"))

nbreg_G1 <- nbreg_G1_raw[-1, .(method = "nbreg", n_obs = (10**(3:6))[c1], G = 1, rep = 1:10, time = c2)]
nbreg_G2 <- nbreg_G2_raw[-1, .(method = "nbreg", n_obs = (10**(3:6))[c1], G = 2, rep = 1:10, time = c2)]
nbreg_G3 <- nbreg_G3_raw[-1, .(method = "nbreg", n_obs = (10**(3:6))[c1], G = 3, rep = 1:10, time = c2)]

results_negbin_all <- rbindlist(list(results_negbin, nbreg_G1, nbreg_G2, nbreg_G3))
results_negbin_all$model <- "Negative Binomial"

#
# LOGIT
#

results_logit <- results_R[family == "logit", -1]

# We add the stata results
logit_G1_raw <- fread(here("_STATA", "logit_G1.txt"))
logit_G2_raw <- fread(here("_STATA", "logit_G2.txt"))
logit_G3_raw <- fread(here("_STATA", "logit_G3.txt"))

logit_G1 <- logit_G1_raw[-1, .(method = "logit", n_obs = (10**(3:6))[c1], G = 1, rep = 1:10, time = c2)]
logit_G2 <- logit_G2_raw[-1, .(method = "logit", n_obs = (10**(3:6))[c1], G = 2, rep = 1:10, time = c2)]
logit_G3 <- logit_G3_raw[-1, .(method = "logit", n_obs = (10**(3:6))[c1], G = 3, rep = 1:10, time = c2)]

results_logit_all <- rbindlist(list(results_logit, logit_G1, logit_G2, logit_G3))
results_logit_all$model <- "Logit"

#
# Merge all
#

results_all <- rbindlist(list(results_poisson_all, results_negbin_all, results_logit_all, results_gaussian_all))

fwrite(results_all, here(RESULTS_DIR, "results_all.txt"))


####
#### "Difficult" Benchmark ####
####

# Only estimations

results_diff_R <- fread(here(RESULTS_DIR, "results_diff_bench_R.txt"))[, -1]

results_diff_reghdfe <- fread(here("_STATA", "reghdfe_diff.txt"))
reghdfe_diff <- results_diff_reghdfe[-1, .(method = "reghdfe", n_obs = 10**c1, G = rep(1:3, each = 10), rep = 1:10, time = c2)]
reghdfe_diff <- reghdfe_diff[n_obs >= 1e4] # only 10K+ obs

results_diff_julia <- fread(here(RESULTS_DIR, "julia_bench_diff.txt"))
julia_diff <- results_diff_julia[, .(method = "FixedEffectModels", n_obs = rep(10**(4:7), each = 30), G = rep(1:3, each = 10), rep = 1:10, time = V1)]

results_diff_all <- rbindlist(list(results_diff_R, reghdfe_diff, julia_diff))
results_diff_all <- results_diff_all[order(n_obs, G, rep, method)]

fwrite(results_diff_all, here(RESULTS_DIR, "results_diff_all.txt"))
