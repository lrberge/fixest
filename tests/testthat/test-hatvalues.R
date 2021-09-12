base <- datab10()
with_parameters_test_that("fixest and stats hatvalues estimations are equal",
  {
    fm <- fixest_mod_select(model = model, fmla = formulas, base = base, famly = family)
    fmm <- stats_mod_select(model = model, fmla = formulas, base = base, famly = family)
    H1 <- unname(hatvalues(fm))
    H2 <- unname(hatvalues(fmm))
    expect_equal(fm$coefficients, fmm$coefficients)
    expect_equal(H1, H2)
  },
  .cases = hatvalues_cases()[1:2, ] # Only hatvalues for poisson and ols are equal
)
