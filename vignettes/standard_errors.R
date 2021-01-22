## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      eval = TRUE,
                      comment = "#>")
Sys.setenv(lang = "en")

library(fixest)

if(requireNamespace("plm", quietly = TRUE)) library(plm)

if(requireNamespace("sandwich", quietly = TRUE)) library(sandwich)

setFixest_nthreads(1)

## -----------------------------------------------------------------------------
library(fixest)
data(trade)
# OLS estimation
gravity = feols(log(Euros) ~ log(dist_km) | Destination + Origin + Product + Year, trade)
# Two-way clustered SEs
summary(gravity, se = "twoway")
# Two-way clustered SEs, without DOF correction
summary(gravity, se = "twoway", dof = dof(adj = FALSE, cluster.adj = FALSE))

## -----------------------------------------------------------------------------
# Data generation
set.seed(0)
N = 20 ; n_id = N/5; n_time = N/n_id
base = data.frame(y = rnorm(N), x = rnorm(N), id = rep(1:n_id, n_time), 
                  time = rep(1:n_time, each = n_id))


## ---- echo = FALSE------------------------------------------------------------
if(!requireNamespace("sandwich", quietly = TRUE)){
    knitr::opts_chunk$set(eval = FALSE)
    cat("Evaluation of the next chunks requires 'sandwich', which is not present.")
} else {
    knitr::opts_chunk$set(eval = TRUE)
}

## -----------------------------------------------------------------------------
library(sandwich)

# Estimations
res_lm    = lm(y ~ x, base)
res_feols = feols(y ~ x, base)

# Same standard-errors
rbind(se(res_lm), se(res_feols))

# Heteroskedasticity-robust covariance
se_lm_hc    = sqrt(diag(vcovHC(res_lm, type = "HC1")))
se_feols_hc = se(res_feols, se = "hetero")
rbind(se_lm_hc, se_feols_hc)

## ---- echo = FALSE------------------------------------------------------------
if(!requireNamespace("plm", quietly = TRUE)){
    knitr::opts_chunk$set(eval = FALSE)
    cat("Evaluation of the next chunks requires 'plm', which is not present.")
} else {
    knitr::opts_chunk$set(eval = TRUE)
}

## -----------------------------------------------------------------------------
#  library(plm)
#  
#  # Estimations
#  est_lm    = lm(y ~ x + as.factor(id) + as.factor(time), base)
#  est_plm   = plm(y ~ x + as.factor(time), base, index = c("id", "time"), model = "within")
#  est_feols = feols(y ~ x | id + time, base)
#  
#  #
#  # "Standard" standard-errors
#  #
#  
#  # By default fixest clusters the SEs when FEs are present,
#  #  so we need to ask for standard SEs explicitly.
#  rbind(se(est_lm)["x"], se(est_plm)["x"], se(est_feols, se = "standard"))
#  
#  # p-values:
#  rbind(pvalue(est_lm)["x"], pvalue(est_plm)["x"], pvalue(est_feols, se = "standard"))
#  

## -----------------------------------------------------------------------------
#  # Clustered by id
#  se_lm_id    = sqrt(vcovCL(est_lm, cluster = base$id, type = "HC1")["x", "x"])
#  se_plm_id   = sqrt(vcovHC(est_plm, cluster = "group")["x", "x"])
#  se_stata_id = 0.165385      # vce(cluster id)
#  se_feols_id = se(est_feols) # By default: clustered according to id
#  
#  rbind(se_lm_id, se_plm_id, se_stata_id, se_feols_id)

## -----------------------------------------------------------------------------
#  # How to get the lm version
#  se_feols_id_lm = se(est_feols, dof = dof(fixef.K = "full"))
#  rbind(se_lm_id, se_feols_id_lm)
#  
#  # How to get the plm version
#  se_feols_id_plm = se(est_feols, dof = dof(fixef.K = "none", cluster.adj = FALSE))
#  rbind(se_plm_id, se_feols_id_plm)

## ---- eval = FALSE------------------------------------------------------------
#  library(lfe)
#  
#  # lfe: clustered by id
#  est_lfe = felm(y ~ x | id + time | 0 | id, base)
#  se_lfe_id = se(est_lfe)
#  
#  # The two are different, and it cannot be directly replicated by feols
#  rbind(se_lfe_id, se_feols_id)
#  #>                     x
#  #> se_lfe_id   0.1458559
#  #> se_feols_id 0.1653850
#  
#  # You have to provide a custom VCOV to replicate lfe's VCOV
#  my_vcov = vcov(est_feols, dof = dof(adj = FALSE))
#  se(est_feols, .vcov = my_vcov * 19/18) # Note that there are 20 observations
#  #>         x
#  #> 0.1458559
#  
#  # Differently from feols, the SEs in lfe are different if time is not a FE:
#  # => now SEs are identical.
#  rbind(se(felm(y ~ x + factor(time) | id | 0 | id, base))["x"],
#        se(feols(y ~ x + factor(time) | id, base))["x"])
#  #>             x
#  #> [1,] 0.165385
#  #> [2,] 0.165385
#  
#  # Now with two-way clustered standard-errors
#  est_lfe_2way = felm(y ~ x | id + time | 0 | id + time, base)
#  se_lfe_2way  = se(est_lfe_2way)
#  se_feols_2way = se(est_feols, se = "twoway")
#  rbind(se_lfe_2way, se_feols_2way)
#  #>                       x
#  #> se_lfe_2way   0.3268584
#  #> se_feols_2way 0.3080378
#  
#  # To obtain the same SEs, use cluster.df = "conventional"
#  sum_feols_2way_conv = summary(est_feols, se = "twoway", dof = dof(cluster.df = "conv"))
#  rbind(se_lfe_2way, se(sum_feols_2way_conv))
#  #>                     x
#  #> se_lfe_2way 0.3268584
#  #>             0.3268584
#  
#  # We also obtain the same p-values
#  rbind(pvalue(est_lfe_2way), pvalue(sum_feols_2way_conv))
#  #>              x
#  #> [1,] 0.3347851
#  #> [2,] 0.3347851

## -----------------------------------------------------------------------------
#  setFixest_dof(dof(adj = FALSE))

## -----------------------------------------------------------------------------
#  setFixest_se(no_FE = "standard", one_FE = "standard", two_FE = "standard")

