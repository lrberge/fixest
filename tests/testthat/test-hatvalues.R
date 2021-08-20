base <- datab10()
patrick::with_parameters_test_that("fixest and stats hatvalues estimations are equal",
                                   {
                                     fm <- fixest_mod_select(model = model, fmla = formulas, base = base, famly = family)
                                     fmm <- stats_mod_select(model = model, fmla = formulas, base = base, famly = family)
                                     H1 <- unname(hatvalues(fm))
                                     H2 <- unname(hatvalues(fmm))
                                     expect_equal(fm$coefficients, fmm$coefficients)
                                     expect_equal(H1, H2)
                                   },
                                   .cases = hatvalues_cases()[1:2,] # Only hatvalues for poisson and ols are equal
)

# l <- 5
# casos <- hatvalues_cases()[l,]
#
# model <- casos$model
# formlas <- casos$formulas
# famly <- casos$family
#
# fm <- fixest_mod_select(model = model, fmla = formlas, base = base, famly = famly)
# fmm <- stats_mod_select(model = model, fmla = formlas, base = base, famly = famly)
# H1 <- unname(hatvalues(fm))
# H2 <- unname(hatvalues(fmm))
#
# expect_equal(fm$coefficients, fmm$coefficients)
# expect_equal(H1, H2)

