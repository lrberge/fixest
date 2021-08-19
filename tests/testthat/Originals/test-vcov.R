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
    vcov_raw <- V_matrix(M_i, M_t, M_it, c_adj, cdf) * my_adj
    testthat::expect_equal(as.numeric(vcov_est), as.numeric(vcov_raw))
  },
  .cases = vcov_cases3()
)



base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$clu <- sample(6, 150, TRUE)
base$clu[1:5] <- NA

test_that("vcov estimation from different sources are equal", {
  est <- feols(y ~ x1 | species, base, cluster = ~clu, dof = dof(adj = FALSE))
  v1 <- est$cov.scaled
  v1b <- vcov(est)
  v1c <- summary(est)$cov.scaled

  expect_equal(as.numeric(v1), as.numeric(v1b))
  expect_equal(as.numeric(v1), as.numeric(v1c))

  # Only dof change
  v2 <- summary(est, dof = dof())$cov.scaled
  v2b <- vcov(est, cluster = ~clu, dof = dof())

  expect_equal(as.numeric(v2), as.numeric(v2b))
  expect_true(max(abs(v1 - v2)) != 0)

  # SE change only
  v3 <- summary(est, se = "hetero")$cov.scaled
  v3b <- vcov(est, se = "hetero", dof = dof(adj = FALSE))

  expect_equal(as.numeric(v3), as.numeric(v3b))
  expect_true(max(abs(v1 - v3)) != 0)
  expect_true(max(abs(v2 - v3)) != 0)
})
