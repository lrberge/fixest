
library(fixest)
setFixest_notes(FALSE)
base <- datab2()
patrick::with_parameters_test_that("fitting without intercept works properly",
                                   code = {
                                       res = fixest_mod_select(model = model,
                                                        fmla = formula,
                                                        data = base,
                                                        famly = family)
                                       expect_model_nointercept(res)
                                       expect_equal(1,1)
                                   },
                                   .cases = nointercept_cases()
)
warnings() ## The `code` argument to `test_that()` must be a braced expression to get accurate file-line information for failures.
setFixest_notes(TRUE)

