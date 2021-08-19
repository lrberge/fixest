
library(fixest)
setFixest_notes(FALSE)
base <- datab7()

patrick::with_parameters_test_that("fixest fits collinearity as stats",
  {
    res <- fixest_mod_select(
      model = model,
      fmla = fixest_formula1,
      base = base,
      famly = "poisson",
      weights = ev_par(useWeights)
    )
    res_bis <- stats_mod_select(
      model = model,
      fmla = stats_formula1,
      base = base,
      famly = "poisson",
      weights = ev_par(useWeights)
    )

    # test(coef(res)["x1"], coef(res_bis)["x1"], "~")
    # test(se(res, se = "st", dof = dof(adj = adj))["x1"], se(res_bis)["x1"], "~")
    #
    expect_model_equal(res, res_bis, method = model)
    expect_equal(1, 1)
  },
  .cases = collin_cases()
)
#
# K = 1
# casos = collin_cases()[K,]
#
# model = casos$model
# fixest_formula1 = casos$fixest_formula1
# stats_formula1 = casos$stats_formula1
# useWeights = casos$useWeights
#
# res <- fixest_mod_select(
#   model = model,
#   fmla = fixest_formula1,
#   base = base,
#   famly = "poisson",
#   weights = ev_par(useWeights)
# )
# res_bis <- stats_mod_select(
#   model = model,
#   fmla = stats_formula1,
#   base = base,
#   famly = "poisson",
#   weights = ev_par(useWeights)
# )
# expect_model_equal(res,res_bis, method = model)
#
#
#
#
# object = res
# reference = res_bis
# tol = 1e-5
#
# testthat::test_that("fixest and stats have equal x1 coefficient", {
#   expect_equal2(coef(object)["x1"],
#                 coef(reference)["x1"],
#                 tolerance = tol,
#                 scale = 1 # Absolute difference
#   )
# })
#
# testthat::test_that("fixest and stats have equal standard errors", {
#   expect_equal2(se(object, se = "st", dof = dof(adj = adj))["x1"],
#                 se(reference)["x1"],
#                 tolerance = tol,
#                 scale = 1 # Absolute difference
#   )
# })
#
# testthat::test_that("fixest and stats have equal p-values", {
#   expect_equal2(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
#                 pvalue(reference)["x1"],
#                 tolerance = tol,
#                 scale = 1 # Absolute difference
#   )
# })
