
library(fixest)
<<<<<<< HEAD
=======
setFixest_notes(FALSE)
>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110
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
<<<<<<< HEAD

=======
warnings() ## The `code` argument to `test_that()` must be a braced expression to get accurate file-line information for failures.
setFixest_notes(TRUE)
>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110

