# base2 <- datab4()
#
# with_parameters_test_that("Santandard, Clustered, HC and tw-way Clustered SE estimations are correct between fixest and stats",
#                           {
#                               est <- fixest_mod_select(model = model,
#                                                        fmla = xpd(lhs ~ x | grp + tm, lhs = y_dep),
#                                                        base = base2,
#                                                        famly = fmly)
#
#                               est_bis <- stats_mod_select(model = model,
#                                                           fmla = xpd(lhs ~ x + as.factor(grp) + as.factor(tm), lhs = y_dep),
#                                                           base = base2,
#                                                           famly = fmly)
#
#                               test_that("Clustered Standard Errors are equal between lm and sandwich", {
#                                   # Standard Errors
#                                   # vcov_fxt <- se(est, se = "standard")["x"]
#                                   vcov_fxt <- est$coeftable[,2]
#                                   vcov_sw <- unname(se(est_bis)["x"])
#                                   expect_equal2(vcov_fxt, vcov_sw, tol = 1e-6)
#
#                                   # Clustered SE type HC0
#                                   vcov_fxt <- unname(se(est, dof = dof(adj = FALSE, fixef.K = "full")))
#                                   vcov_sw <- sqrt(vcovCL(est_bis, cluster = base2$grp, type = "HC0")["x", "x"])
#                                   expect_equal2(vcov_fxt, vcov_sw, tol = 1e-6) # Almost equal
#
#                                   # Clusteres SE type HC1
#                                   vcov_fxt <- unname(se(est, dof = dof(fixef.K = "full")))
#                                   vcov_sw <- sqrt(vcovCL(est_bis, cluster = base2$grp, type = "HC1")["x", "x"])
#                                   expect_equal2(vcov_fxt, vcov_sw, tol = 1e-6) # Almost equal
#                               })
#
#                               ### Heteroskedasticity-robust
#                               test_that("Heterscoedastcity-robus Standard Errors are equal between feols and sandwich", {
#                                   vcov_fxt <- unname(se(est, se = "hetero", dof = dof(adj = FALSE, cluster.adj = FALSE)))
#                                   vcov_sw <- sqrt(vcovHC(est_bis, type = "HC0")["x", "x"])
#                                   expect_equal2(vcov_fxt, vcov_sw, tol = 1e-5)
#
#                                   vcov_fxt <- unname(se(est, se = "hetero"))
#                                   vcov_sw <- sqrt(vcovHC(est_bis, type = "HC1")["x", "x"])
#                                   expect_equal2(vcov_fxt, vcov_sw, tol = 1e-5)
#                               })
#
#                               ### Two-way clustered standard errors
#                               test_that("Two-way Clustered Standard Errors are equal between feols and sandwich", {
#                                   vcov_fxt <- unname(se(est, se = "twoway", dof = dof(fixef.K = "full", cluster.df = "conv")))
#                                   vcov_sw <- sqrt(vcovCL(est_bis, cluster = ~ grp + tm, type = "HC1")["x", "x"])
#                                   expect_equal2(vcov_fxt, vcov_sw, tol = 1e-5)
#                               })
#                               expect_equal(1,1)
#                           },
#                           .cases = SE_cases()[1:6,]
# )
#
#
#
# ## cluster = data vs cluster = formula
# test_that("cluster = data vs cluster = formula works properly", {
#     for (k in 1:3) {
#         est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, vcov_db(k))
#         expect_equal_vcov(est_pois)
#     }
# })
#
#
# ## Testing Errors
# test_that("SE is reporting errors correctly", {
#     est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, vcov_db(3))
#
#     expect_error(se(est_pois, cluster = "Origin_na"))
#     expect_error(se(est_pois, cluster = base$Origin_na))
#     expect_error(se(est_pois, cluster = list(datab4()$Origin_na)))
#     expect_error(se(est_pois, cluster = ~ Origin_na^Destination))
#     expect_error(se(est_pois, se = "cluster", cluster = ~ Origin_na^not_there))
#     expect_error(se(est_pois, se = "cluster", cluster = ~ Origin_na^not_there))
#     expect_error(se(est_pois, se = "twoway", cluster = c("Origin^Destination", "Product", "error")))
#     expect_error(se(est_pois, se = "twoway", cluster = datab4()[, 1:4]))
#     expect_error(se(est_pois, se = "twoway", cluster = ~ Origin + Destination + Product))
#     expect_error(se(est_pois, se = "fourway", cluster = ~ Origin + Destination + Product))
# })
#
# ## Testing aliases
#
# test_that("Aliases works properly", {
#     est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, vcov_db(3))
#
#     se_hetero <- se(est_pois, se = "hetero")
#     se_white <- se(est_pois, se = "white")
#     # se_hc1    = se(est_pois, se = "hc1") # Hc1 gives error -> no matches found for hc1
#
#     expect_equal(se_hetero, se_white)
#     # expect_equal(se_hetero, se_hc1)
# })
#
#
#
# ## Testing compatibility between fixest and sandwich
#
# ## Only vcovCL is tested. Extend to se = hetero, twoway, threeway
# ## Parametrize for future extension
#
# data(base_did)
# base <- base_did
# base$y_int <- as.integer(base$y) + 20
#
# with_parameters_test_that("fixest is compatible with sandwich's vcov",
#                           {
#                               est <- fixest_mod_select(
#                                   model = model,
#                                   fmla = fmlas,
#                                   base = base,
#                                   famly = family
#                               )
#
#                               if (isFALSE(FE)) {
#                                   expect_equal(vcov(est, cluster = ~id), vcovCL(est, cluster = ~id, type = "HC1"))
#                               } else {
#                                   expect_equal(vcov(est, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est, cluster = ~id))
#                               }
#                           },
#                           .cases = sandwcomp_cases()
# )
#
# # K = 1
# # casos = sandwcomp_cases()[K,]
# # model = casos$model
# # fmlas = casos$fmlas
# # fmla = as.formula(fmlas)
# # # data = base
# # famly = casos$family
# # Data = base
# # est <- fixest_mod_select(
# #   model = model,
# #   fmla = fmlas,
# #   base = base,
# #   famly = famly,
# #   weights = NULL
# # )
# #
# # expect_equal(vcov(est, cluster = ~id), vcovCL(est, cluster = ~id, type = "HC1"))
# # expect_equal(vcov(est, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est, cluster = ~id))
# #
# # ## I must call my database with the same name of fixest_mod_select Data argument
# # ## For some reason vcovCL doesnt work without the previous detail
#
# #
# # est_lm <- lm(y ~ x + as.factor(grp) + as.factor(tm), data = base)
# # est_feols <- feols(y ~ x | grp + tm, data = base)
# #
# # est_lm <- glm(y_int ~ x + as.factor(grp) + as.factor(tm), data = base, family = "poisson")
# # est_feols <- feglm(y_int ~ x | grp + tm, data = base, family = "poisson")
# #
# # # est_lm <- glm(y_int ~ x + as.factor(grp) + as.factor(tm), data = base, family = "poisson")
# # # est_feols <- femlm(y_int ~ x | grp + tm, data = base, family = "poisson")
# #
# # summary(est_lm)$coefficients[2, ]
# # est_feols$coeftable
# # est_feols # Notice that the standard error shown here is clustered
# # se(est_feols, se = "standard")
# # coeftable(est_feols, se = "standard")
# # coeftable(est_feols)
