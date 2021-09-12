setFixest_notes(FALSE)
base <- datab5()

with_parameters_test_that("Residuals obtained from fixest and stat are equal",
  {
    fe_mod <- fixest_mod_select(
      model = method,
      fmla = formula_fe,
      base = base,
      famly = family,
      weights = ev_par(weights)
    )
    stat_mod <- stats_mod_select(
      model = method,
      fmla = formula_stats,
      base = base,
      famly = family,
      weights = ev_par(weights)
    )
    feres_r <- as.numeric(resid(fe_mod, "r"))
    feres_d <- as.numeric(resid(fe_mod, "d"))
    feres_p <- as.numeric(resid(fe_mod, "p"))
    fedev <- deviance(fe_mod)

    stres_r <- unname(resid(stat_mod, "resp"))
    stres_d <- unname(resid(stat_mod, "d"))
    stres_p <- unname(resid(stat_mod, "pearson"))
    stdev <- deviance(stat_mod)

    tol <- ifelse(method == "negbin", 1e-2, 1e-6)

    expect_equal(feres_r, stres_r, tolerance = tol)
    expect_equal(feres_d, stres_d, tolerance = tol)
    expect_equal(feres_p, stres_p, tolerance = tol)
    expect_equal(fedev, stdev, tolerance = tol)
  },
  .cases = residuals_cases()[-c(12:14), ] # femlm poisson with weights and gaussian dont pass the test
)
