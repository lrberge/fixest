
setFixest_notes(FALSE)
base <- datab5()

patrick::with_parameters_test_that("Fitting without intercept works properly",
  {
    fe_mod <- fixest_mod_select(
      model = method,
      fmla = formula_fe,
      base = base,
      famly = family,
      weights = ev_par(wghts)
    )
    stat_mod <- stats_mod_select(
      model = method,
      fmla = formula_stats,
      base = base,
      famly = family,
      weights = ev_par(wghts)
    )
    feres_r <- resid(fe_mod, "r")
    feres_d <- resid(fe_mod, "d")
    feres_p <- resid(fe_mod, "p")
    fedev <- deviance(fe_mod)

    stres_r <- unname(resid(stat_mod, "resp"))
    stres_d <- unname(resid(stat_mod, "d"))
    stres_p <- unname(resid(stat_mod, "pearson"))
    stdev <- deviance(stat_mod)

    expect_equal(feres_r, stres_r, tolerance = tol)
    expect_equal(feres_d, stres_d, tolerance = tol)
    expect_equal(feres_p, stres_p, tolerance = tol)
    expect_equal(fedev, stdev, tolerance = tol)
  },
  .cases = residuals_cases()
)
