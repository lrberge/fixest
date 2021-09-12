# base <- datab18()
# test_that("there is no bug in corner cases for feols", {
#   res <- feols(y ~ 1 | csw(fe1, fe1^fe2), base)
#   res <- feols(y ~ 1 + csw(x1, i(fe1)) | fe2, base)
#   res <- feols(y ~ csw(f(x1, 1:2), x2) | sw0(fe2, fe2^fe3), base, panel.id = ~ fe1 + period)
#   res <- feols(d(y) ~ -1 + d(x2), base, panel.id = ~ fe1 + period)
#   expect_equal(length(coef(res)), 1)
#   res <- feols(c(y, x1) ~ 1 | fe1 | x2 ~ x3, base)
#   res <- feols(y ~ x1 | fe1[x2] + fe2[x2], base)
# })
#
# test_that("Error when warn = TRUE",
#           {
#             res = feols(y ~ x_cst | fe1, base, warn = FALSE)
#             expect_silent(res)         # => no error
#             expect_silent(etable(res)) # => no error
#             expect_error(feols(y ~ x_cst | fe1, base))
#
#             res = feols(c(y, x1) ~ x_cst | fe1, base, warn = FALSE) # warn doesnt work with multiple estimation
#             expect_silent(res)         # => no error
#             expect_silent(etable(res)) # => no error
#           })
