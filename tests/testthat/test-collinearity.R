
#19:30 12-08
library(fixest)
setFixest_notes(FALSE)
base = datab7()

patrick::with_parameters_test_that("fixest fits collinearity as stats",
                                   {
                                       res = fixest_mod_select(model = model,
                                                         fmla = fixest_formula1,
                                                         data = base,
                                                         famly = "poisson",
                                                         weights = ev_par(useWeights))
                                       res_bis = stats_mod_select(model = model,
                                                        fmla = stats_formula1,
                                                        data = base,
                                                        famly = "poisson",
                                                        weights = ev_par(useWeights))
                                       expect_equal(1,1)
                                   },
                                   .cases = collin_cases()
                                 )
