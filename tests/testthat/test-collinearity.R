# setFixest_notes(FALSE)
# base <- datab7()
# with_parameters_test_that("fixest fits collinearity as stats", {
#     fix_fmla <- paste(y_dep, " ~ ", fixest_formula1)
#     stat_fmla <- paste(y_dep, " ~ ", stats_formula2)
#
#     res <- fixest_mod_select(
#         model = model,
#         fmla = fix_fmla,
#         base = base,
#         famly = fmly,
#         weights = ev_par(useWeights)
#     )
#     res_bis <- stats_mod_select(
#         model = model,
#         fmla = stat_fmla,
#         base = base,
#         famly = fmly,
#         weights = ev_par(useWeights)
#     )
#     expect_model_equal(res, res_bis, method = model)
#     expect_equal(1, 1)
# },
# .cases = collin_cases2()[1:12,] # femlm does not omit constant on collinearity
# )
