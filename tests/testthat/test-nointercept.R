base <- datab2()
with_parameters_test_that("fitting without intercept works properly",
  {
    fmla <- xpd(lhs ~ rhs, lhs = y_dep, rhs = f_rhs)
    res <- fixest_mod_select(
      model = method,
      fmla = fmla,
      base = base,
      famly = fmly,
    )
    expect_model_nointercept(res)
    expect_equal(1, 1)
  },
  .cases = nointercept_cases()
)
