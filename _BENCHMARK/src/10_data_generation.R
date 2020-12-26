#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Oct 18 17:05:15 2019
# ~: Benchmarking: data generation
#----------------------------------------------#

####
#### SIMULATION ####
####

library(MASS)
library(here)
library(data.table)
library(haven)

# Some constants

DATA_DIR <- "data"
STATA_DIR <- "_STATA"
RESULTS_DIR <- "results"

dir.create(DATA_DIR, showWarnings = FALSE)
dir.create(STATA_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)

# We simulate databases

all_n <- 1000 * 10**(0:3)
all_rep <- 1:10

# META parameters
a <- 1
b <- 0.05

# Array of lists to store the results
base_all <- array(data.frame(), dim = c(length(all_n), length(all_rep)))

for (i in 1:length(all_n)) {
  cat("i =", i)
  n <- all_n[i]

  dum_all <- list()

  nb_dum <- c(n / 20, floor(sqrt(n)), floor(n**.33))
  N <- nb_dum**3
  dum_all[[1]] <- sample(nb_dum[1], n, TRUE)
  dum_all[[2]] <- sample(nb_dum[2], n, TRUE)
  dum_all[[3]] <- sample(nb_dum[3], n, TRUE)

  for (r in all_rep) {
    cat(".")

    X1 <- rnorm(n)
    X2 <- X1**2

    mu <- a * X1 + b * X2

    for (m in 1:3) {
      coef_dum <- rnorm(nb_dum[m])
      mu <- mu + coef_dum[dum_all[[m]]]
    }

    mu <- exp(mu)
    y <- rnegbin(mu, theta = 0.5)

    base <- data.frame(y, X1, ln_y = log(y + 1))

    for (m in 1:3) {
      base[[paste0("dum_", m)]] <- dum_all[[m]]
    }

    base_all[i, r][[1]] <- base
  }
  cat("\n")
}

save(base_all, file = here(DATA_DIR, "base_all_simulations.Rdata"))

#
# Data with 10M observation for OLS (just one to save size)
#

n <- 1e7
dum_all <- list()

nb_dum <- c(n / 20, floor(sqrt(n)), floor(n**.33))
N <- nb_dum**3
dum_all[[1]] <- sample(nb_dum[1], n, TRUE)
dum_all[[2]] <- sample(nb_dum[2], n, TRUE)
dum_all[[3]] <- sample(nb_dum[3], n, TRUE)

X1 <- rnorm(n)
X2 <- X1**2

mu <- a * X1 + b * X2

for (m in 1:3) {
  coef_dum <- rnorm(nb_dum[m])
  mu <- mu + coef_dum[dum_all[[m]]]
}

mu <- exp(mu)
y <- rnegbin(mu, theta = 0.5)

base <- data.frame(ln_y = log(y + 1), X1)

for (m in 1:3) {
  base[[paste0("dum_", m)]] <- dum_all[[m]]
}

fwrite(base, here(DATA_DIR, "base_10M.csv"))



#
# Exportation to stata
#

load(here(DATA_DIR, "base_all_simulations.Rdata"))

# base_all: 4 sizes, 2 groups, 10 replications

for (size in 1:4) {
  for (replication in 1:10) {
    cat(".")
    stata_name <- paste0("_STATA/base_s", size, "_r", replication, ".dta")
    write_dta(as.data.frame(base_all[size, replication]), stata_name)
  }
}


####
#### Difficult Data ####
####

# This benchmark data set is an adaptation of a benchmark of
# the authors of the Julia FixedEffectModels.jl software

set.seed(1) # for replication
base_all_diff <- list()

for (pow in 4:7) {
  cat(".")
  n <- 10**pow
  nb_indiv <- n / 20
  nb_firm <- round(n / 160)
  nb_year <- round(n**.3)

  id_indiv <- sample(1:nb_indiv, n, TRUE)
  id_firm <- pmin(sample(0:20, n, TRUE) + pmax(1, id_indiv %/% 8 - 10), nb_firm)
  id_year <- sample(nb_year, n, TRUE)

  x1 <- 5 * cos(id_indiv) + 5 * sin(id_firm) + 5 * sin(id_year) + runif(n)
  x2 <- cos(id_indiv) + sin(id_firm) + sin(id_year) + rnorm(n)
  y <- 3 * x1 + 5 * x2 + cos(id_indiv) + cos(id_firm)^2 + sin(id_year) + rnorm(n)
  df <- data.frame(id_indiv = id_indiv, id_firm = id_firm, id_year = id_year, x1 = x1, x2 = x2, y = y)

  base_all_diff[[length(base_all_diff) + 1]] <- df

  if (pow < 7) {
    stata_name <- paste0("_STATA/base_diff_e", pow, ".dta")
    haven::write_dta(df, stata_name)
  }
}

save(base_all_diff, file = here(DATA_DIR, "base_all_diff.Rdata"))
