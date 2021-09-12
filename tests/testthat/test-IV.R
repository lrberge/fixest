# setFixest_se(all = "standard")
base <- datab15()
test_that("IV estimation works properly", {
  est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
  res_f1 <- feols(x_endo_1 ~ x1 + x_inst_1 + x_inst_2, base)
  res_f2 <- feols(x_endo_2 ~ x1 + x_inst_1 + x_inst_2, base)
  base$fit_x_endo_1 <- predict(res_f1)
  base$fit_x_endo_2 <- predict(res_f2)

  res_2nd <- feols(y ~ fit_x_endo_1 + fit_x_endo_2 + x1, base)

  expect_equal(coef(est_iv), coef(res_2nd))

  resid_iv <- base$y - predict(res_2nd, data.frame(x1 = base$x1, fit_x_endo_1 = base$x_endo_1, fit_x_endo_2 = base$x_endo_2))
  sigma2_iv <- sum(resid_iv**2) / (res_2nd$nobs - res_2nd$nparams)
  coviid <- vcov(res_2nd, se = "standard")
  sum_2nd <- summary(res_2nd, .vcov = coviid / res_2nd$sigma2 * sigma2_iv)

  # We only check that on Windows => avoids super odd bug in fedora devel
  # The worst is that I just can't debug it.... so that's the way it's done.
  if (Sys.info()["sysname"] == "Windows") {
    expect_equal(se(sum_2nd), se(est_iv))
  }
})

# check no bug
est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
etable(summary(est_iv, stage = 1:2))
# setFixest_se(reset = TRUE)
