
library(fixest)
base <- datab2()
patrick::with_parameters_test_that("fitting without intercept works properly",
  code = {
    res <- fixest_mod_select(
      model = model,
      fmla = formula,
      base = base,
      famly = family
    )
    expect_model_nointercept(res)
    expect_equal(1, 1)
  },
  .cases = nointercept_cases()
)
