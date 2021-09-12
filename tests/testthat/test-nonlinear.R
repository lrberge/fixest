# base <- datab17()
# test_that("estimation of a linear model between feNmlm and feols are equal", {
#   est_nl <- feNmlm(y ~ x1, data = base, NL.fml = ~ fun_nl(a, b, tab2), NL.start = 1, family = "gaussian")
#   est_lin <- feols(y ~ x1 + var_spec + I(var_spec^2), base)
#
#   coef_nl <- coef(est_nl)
#   coef_lin <- coef(est_lin)[c(3, 4, 1, 2)]
#   expect_equal2(unname(coef_nl), unname(coef_lin))
# })
