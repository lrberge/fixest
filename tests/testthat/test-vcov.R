## Ctrl + Shift + S to run all the script

library(fixest)
vcovClust <- fixest:::vcovClust

### https://lrberge.github.io/fixest/articles/standard_errors.html

### Standard SEs
base <- datab3()
est <- feols(y ~ x | fe1 + fe2, base)
n <- est$nobs
VCOV_raw <- vcov(est, se = "standard") / ((n - 1) / (n - est$nparams))

patrick::with_parameters_test_that("Standard VarCov estimations are correct",
  {
    vcov_est <- vcov(est, se = "standard", dof = dof(adj = adj, fixef.K = k_val))
    vcov_raw <- VCOV_raw * my_adj
    testthat::expect_equal(vcov_est, vcov_raw)
  },
  .cases = vcov_cases1()
)

### Clustered SEs
VCOV_raw <- vcov(est, se = "standard") / est$sigma2
H <- vcovClust(est$fixef_id$fe1, VCOV_raw, scores = est$scores, dof = FALSE)

patrick::with_parameters_test_that("Clustered VarCov estimations are correct",
  {
    vcov_est <- vcov(est, se = "cluster", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj))
    vcov_raw <- H * cluster_factor * my_adj
    testthat::expect_equal(as.numeric(vcov_est), as.numeric(vcov_raw))
  },
  .cases = vcov_cases2()
)

### Two-way clustered se
VCOV_raw <- vcov(est, se = "standard") / est$sigma2
M_i <- vcovClust(est$fixef_id$fe1, VCOV_raw, scores = est$scores, dof = FALSE)
M_t <- vcovClust(est$fixef_id$fe2, VCOV_raw, scores = est$scores, dof = FALSE)
M_it <- vcovClust(paste(base$fe1, base$fe2), VCOV_raw, scores = est$scores, dof = FALSE, do.unclass = TRUE)


patrick::with_parameters_test_that("Two-Way Clustered VarCov estimations are correct",
  {
    vcov_est <- vcov(est, se = "two", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj, cluster.df = cdf))
    vcov_raw <- V_matrix(Mi, M_t, M_it, c_adj, cdf) * my_adj
    testthat::expect_equal(as.numeric(vcov_est), as.numeric(vcov_raw))
  },
  .cases = vcov_cases3()
)

