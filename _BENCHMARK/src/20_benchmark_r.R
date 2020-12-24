#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Wed Jun 26 16:59:14 2019
# Purpose: Benchmarking -- full
#----------------------------------------------#

library(data.table)

library(MASS)
library(lfe)
library(glmmML)
library(alpaca)
library(plm)
library(fixest)
library(here)

####
#### Estimations ####
####

# Methods:
# - OLS
# - poisson
# - negbin
# - logit

# The simulated data sets (see script "Data generation.R" for details on the generation process)
load(here("DATA", "base_all_simulations.Rdata"))

# Loading the packages
# install.packages(c("fixest", "lfe", "glmmML", "alpaca", "plm"))

# custom functions for time monitoring
getime <- function(x) (proc.time() - x)[[3]]
# same + error handling (e.g. out of memory problems)
getime_errcheck <- function(x, y) ifelse("try-error" %in% class(y), NA, (proc.time() - x)[[3]])

# observations + replication number
all_n <- 1000 * 10**(0:3)
all_rep <- 1:10

results_all <- data.frame(family = NA, method = NA, n_obs = NA, G = NA, rep = NA, time = NA)

all_fixef <- c("dum_1", "dum_2", "dum_3")
all_families <- c("poisson", "gaussian", "negbin", "logit")

setFixest_notes(FALSE)

for (fam in all_families) {
  cat("-------\nFAMILY: ", fam, "\n-------\n")
  for (i in 1:length(all_n)) {
    cat("---- i =", i, "\n")
    n_obs <- all_n[i]

    for (g in 1:3) {
      cat("g =", g)
      for (r in all_rep) {
        cat(".")

        # we retrieve the right index
        base <- base_all[i, r][[1]]

        #
        # Negative binomial
        #

        if (fam == "negbin") {

          # fixest
          pt <- proc.time()
          res_fixest <- fenegbin(y ~ X1, base, fixef = all_fixef[1:g])
          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "fenegbin", n_obs = n_obs, G = g, rep = r,
              time = getime(pt)
            )
          )

          # glm.nb
          if (i <= 2) {
            pt <- proc.time()
            fml_glmnb <- switch(g, "1" = y ~ X1 + factor(dum_1), "2" = y ~ X1 + factor(dum_1) + factor(dum_2), "3" = y ~ X1 + factor(dum_1) + factor(dum_2) + factor(dum_3))
            res_glmnb <- try(glm.nb(fml_glmnb, base), silent = TRUE)

            results_all <- rbind(
              results_all,
              data.frame(
                family = fam, method = "glmnb", n_obs = n_obs, G = g, rep = r,
                time = getime_errcheck(pt, res_glmnb)
              )
            )
          }
        }

        #
        # Gaussian
        #

        if (fam == "gaussian") {

          # fixest
          pt <- proc.time()
          res_fixest <- feols(ln_y ~ X1, base, fixef = all_fixef[1:g])
          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "feols", n_obs = n_obs, G = g, rep = r,
              time = getime(pt)
            )
          )

          # lfe
          pt <- proc.time()
          fml_lfe <- switch(g, "1" = ln_y ~ X1 | dum_1, "2" = ln_y ~ X1 | dum_1 + dum_2, "3" = ln_y ~ X1 | dum_1 + dum_2 + dum_3)
          res_lfe <- try(felm(fml_lfe, base), silent = TRUE)

          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "lfe", n_obs = n_obs, G = g, rep = r,
              time = getime_errcheck(pt, res_lfe)
            )
          )

          # plm => Too slow, dropped
          if (FALSE && g == 1) {
            pt <- proc.time()
            res_plm <- try(plm(ln_y ~ X1, base, index = "dum_1", model = "within"), silent = TRUE)
            results_all <- rbind(
              results_all,
              data.frame(
                family = fam, method = "plm", n_obs = n_obs, G = g, rep = r,
                time = getime_errcheck(pt, res_plm)
              )
            )
          }
        }

        #
        # Logit
        #

        if (fam == "logit") {

          # fixest
          pt <- proc.time()
          res_fixest <- fixest::feglm(sign(y) ~ X1, base, fixef = all_fixef[1:g], family = binomial())
          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "feglm (fixest)", n_obs = n_obs, G = g, rep = r,
              time = getime(pt)
            )
          )

          # glmmboot
          if (g == 1) {
            pt <- proc.time()
            res_glmmML <- try(glmmboot(sign(y) ~ X1, base, family = binomial, cluster = base$dum_1), silent = FALSE)

            results_all <- rbind(
              results_all,
              data.frame(
                family = fam, method = "glmmboot", n_obs = n_obs, G = g, rep = r,
                time = getime_errcheck(pt, res_glmmML)
              )
            )
          }

          # alpaca
          base$sy <- sign(base$y) # alpaca does not support LHS formulas!!!
          pt <- proc.time()
          fml_alpaca <- switch(g, "1" = sy ~ X1 | dum_1, "2" = sy ~ X1 | dum_1 + dum_2, "3" = sy ~ X1 | dum_1 + dum_2 + dum_3)
          res_alpaca <- try(alpaca::feglm(fml_alpaca, base), silent = TRUE)
          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "feglm (alpaca)", n_obs = n_obs, G = g, rep = r,
              time = getime_errcheck(pt, res_alpaca)
            )
          )
        }

        #
        # Poisson
        #

        if (fam == "poisson") {

          # fixest
          pt <- proc.time()
          res_fixest <- fepois(y ~ X1, base, fixef = all_fixef[1:g])
          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "fepois", n_obs = n_obs, G = g, rep = r,
              time = getime(pt)
            )
          )

          # glmmboot
          if (g == 1) {
            pt <- proc.time()
            res_glmmML <- try(glmmboot(y ~ X1, base, family = poisson, cluster = base$dum_1), silent = FALSE)

            results_all <- rbind(
              results_all,
              data.frame(
                family = fam, method = "glmmboot", n_obs = n_obs, G = g, rep = r,
                time = getime_errcheck(pt, res_glmmML)
              )
            )
          }

          # alpaca
          pt <- proc.time()
          fml_alpaca <- switch(g, "1" = y ~ X1 | dum_1, "2" = y ~ X1 | dum_1 + dum_2, "3" = y ~ X1 | dum_1 + dum_2 + dum_3)
          res_alpaca <- try(alpaca::feglm(fml_alpaca, base, family = poisson()), silent = TRUE)
          results_all <- rbind(
            results_all,
            data.frame(
              family = fam, method = "feglm (alpaca)", n_obs = n_obs, G = g, rep = r,
              time = getime_errcheck(pt, res_alpaca)
            )
          )
        }
      }
      cat("\n")
    }
  }
}

# removing first NA line
results_all <- results_all[-1, ]

#
# Adding the 10M obs OLS
#

base <- fread(here("DATA", "base_10M.csv"))

for (g in 1:3) {
  cat("g =", g)
  for (r in 1:10) {
    cat("|")
    # fixest
    pt <- proc.time()
    res_fixest <- feols(ln_y ~ X1, base, fixef = all_fixef[1:g])
    results_all <- rbind(
      results_all,
      data.frame(
        family = "gaussian", method = "feols", n_obs = 1e7, G = g, rep = r,
        time = getime(pt)
      )
    )
    cat(".")
    # lfe
    pt <- proc.time()
    fml_lfe <- switch(g, "1" = ln_y ~ X1 | dum_1, "2" = ln_y ~ X1 | dum_1 + dum_2, "3" = ln_y ~ X1 | dum_1 + dum_2 + dum_3)
    res_lfe <- felm(fml_lfe, base)

    results_all <- rbind(
      results_all,
      data.frame(
        family = "gaussian", method = "lfe", n_obs = 1e7, G = g, rep = r,
        time = getime_errcheck(pt, res_lfe)
      )
    )
  }
  cat("\n")
}


data.table::fwrite(results_all, file = here("DATA", "results_bench_R.txt"))


####
#### "Difficult" Benchmark ####
####

# We benchmark only in OLS, otherwise, too long

library(fixest)
library(lfe)

load(here("DATA", "base_all_diff.Rdata"))

results_all_diff <- data.frame(family = NA, method = NA, n_obs = NA, G = NA, rep = NA, time = NA)

for (i in 1:4) {
  cat("------\nOBS = ", 10**(3 + i), "\n------\n")

  base <- base_all_diff[[i]]

  for (g in 1:3) {
    cat("g =", g)

    for (r in 1:10) {
      fml <- switch(g, "1" = y ~ x1 + x2 | id_indiv,
        "2" = y ~ x1 + x2 | id_indiv + id_firm,
        "3" = y ~ x1 + x2 | id_indiv + id_firm + id_year
      )

      cat("|")
      # fixest
      pt <- proc.time()
      res_fixest <- feols(fml, base)
      results_all_diff <- rbind(
        results_all_diff,
        data.frame(
          family = "gaussian", method = "feols", n_obs = 10**(3 + i), G = g, rep = r,
          time = getime(pt)
        )
      )
      cat(".")

      # lfe
      if (i < 4) {
        pt <- proc.time()

        res_lfe <- felm(fml, base)

        results_all_diff <- rbind(
          results_all_diff,
          data.frame(
            family = "gaussian", method = "lfe", n_obs = 10**(3 + i), G = g, rep = r,
            time = getime_errcheck(pt, res_lfe)
          )
        )
      }
    }
    cat("\n")
  }
}

results_all_diff <- results_all_diff[-1, ]

data.table::fwrite(results_all_diff, file = here("DATA", "results_diff_bench_R.txt"))
