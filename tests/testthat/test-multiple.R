# setFixest_notes(FALSE)
# base <- datab14()
# species_indx <- list(
#   "setosa" = which(base$species == "setosa"),
#   "versicolor" = which(base$species == "versicolor"),
#   "virginica" = which(base$species == "virginica")
# )
#
#
# ### Feols multiple estimation
# Model <- "ols"
# est_multi <- fixest_mod_select(
#   model = Model,
#   fmla = c(y1, y2) ~ x1 + sw(x2, x3),
#   base = base,
#   split = ~species
# )
# with_parameters_test_that("multiple estimations works properly for feols",
#                                    {
#                                      est <- est_multi[[num_fmla]]
#
#                                      base_aux <- base[species_indx[[s]], ]
#                                      fmlas <- xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs)
#                                      res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#
#                                      expect_equal(coef(est), coef(res))
#                                      expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
#                                    },
#                                    .cases = multiple_cases(Model)
# )
#
# #### Feglm multiple estimation
# Model <- "glm"
# famlys <- c("poisson", "quasipoisson", "gaussian", "quasi") # cant add Gamma and inverse.gaussian beucase min(base$y2) == 0
# for(l in 1:length(famlys)){ ## Testing different families
#   est_multi <- fixest_mod_select(
#     model = Model,
#     fmla = c(y1, y2) ~ x1 + sw(x2, x3),
#     base = base,
#     split = ~species,
#     famly = famlys[l]
#   )
#   with_parameters_test_that("multiple estimations works properly for feglm",
#                                      {
#                                        est <- est_multi[[num_fmla]]
#
#                                        base_aux <- base[species_indx[[s]], ]
#                                        fmlas <- xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs)
#                                        res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#
#                                        expect_equal(coef(est), coef(res))
#                                        expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
#                                      },
#                                      .cases = multiple_cases(Model, famly = famlys[l])
#   )
# }
#
# #### Femlm multiple estimation
# Model <- "femlm"
# famlys <- c("poisson", "gaussian", "negbin")
# for(l in 1:length(famlys)){
#   est_multi <- fixest_mod_select(
#     model = Model,
#     fmla = c(y1, y2) ~ x1 + sw(x2, x3),
#     base = base,
#     split = ~species,
#     famly = famlys[l]
#   )
#   with_parameters_test_that("multiple estimations works properly for feglm",
#                                      {
#                                        est <- est_multi[[num_fmla]]
#
#                                        base_aux <- base[species_indx[[s]], ]
#                                        fmlas <- xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs)
#                                        res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#
#                                        expect_equal(coef(est), coef(res))
#                                        expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
#                                      },
#                                      .cases = multiple_cases(Model, famly = famlys[l])
#   )
# }
#
#
# #### FeNmlm multiple estimation
# ## Is it necessary to test the following knowing that femlm has been already tested?
# Model <- "feNmlm"
# famlys <- c("poisson", "gaussian", "negbin")
# for(l in 1:length(famlys)){
#   est_multi <- fixest_mod_select(
#     model = Model,
#     fmla = c(y1, y2) ~ x1 + sw(x2, x3),
#     base = base,
#     split = ~species,
#     famly = famlys[l]
#   )
#   with_parameters_test_that("multiple estimations works properly for feglm",
#                                      {
#                                        est <- est_multi[[num_fmla]]
#
#                                        base_aux <- base[species_indx[[s]], ]
#                                        fmlas <- xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs)
#                                        res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#
#                                        expect_equal(coef(est), coef(res))
#                                        expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
#                                      },
#                                      .cases = multiple_cases(Model, famly = famlys[l])
#   )
# }
#
#
# ############ Test that multiple estimation with fixed effects
# ############ works properly for feols, feglm, femlm and feNmlm
#
# ##### feols
# Model <- "ols"
# all_rhs <- c("", "x2", "x3")
# est_multi <- fixest_mod_select(
#   model = Model,
#   fmla = c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2,
#   base = base,
#   fsplit = ~species
# )
# with_parameters_test_that("multiple estimations with fixed effects works properly for feols",
#                                    {
#                                      fmlas <- xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs])
#                                      if (s == "all") {
#                                        res <- fixest_mod_select(model = Model, fmla = fmlas, base = base, famly = fmly, notes = FALSE)
#                                      } else {
#                                        base_aux <- base[species_indx[[s]], ]
#                                        res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#                                      }
#
#                                      est <- est_multi[[num_fmla]]
#
#                                      vname <- names(coef(res))
#                                      expect_equal2(coef(est_multi[[num_fmla]])[vname], coef(res), tolerance = 1e-6)
#                                      expect_equal2(se(est_multi[[num_fmla]], cluster = "fe3")[vname], se(res, cluster = "fe3"), tolerance = 1e-6)
#                                    },
#                                    .cases = multiple_cases2(Model)
# )
#
#
# ##### feglm
# Model <- "glm"
# all_rhs <- c("", "x2", "x3")
# famlys <- c("poisson", "quasipoisson", "gaussian", "quasi") # cant add Gamma and inverse.gaussian beucase min(base$y2) == 0
# for(l in 1:length(famlys)){
#   est_multi <- fixest_mod_select(
#     model = Model,
#     fmla = c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2,
#     base = base,
#     fsplit = ~species,
#     famly = famlys[l]
#   )
#   with_parameters_test_that("multiple estimations with fixed effects works properly for feglm",
#                                      {
#                                        fmlas <- xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs])
#                                        if (s == "all") {
#                                          res <- fixest_mod_select(model = Model, fmla = fmlas, base = base, famly = fmly, notes = FALSE)
#                                        } else {
#                                          base_aux <- base[species_indx[[s]], ]
#                                          res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#                                        }
#
#                                        est <- est_multi[[num_fmla]]
#
#                                        vname <- names(coef(res))
#                                        expect_equal2(coef(est_multi[[num_fmla]])[vname], coef(res), tolerance = 1e-6)
#                                        expect_equal2(se(est_multi[[num_fmla]], cluster = "fe3")[vname], se(res, cluster = "fe3"), tolerance = 1e-6)
#                                      },
#                                      .cases = multiple_cases2(Model, famly = famlys[l])
#   )
# }
#
# ##### femlm
# Model <- "femlm"
# famlys <- c("poisson", "gaussian", "negbin")
# for(l in 1:length(famlys)){
#   est_multi <- fixest_mod_select(
#     model = Model,
#     fmla = c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2,
#     base = base,
#     fsplit = ~species,
#     famly = famlys[l]
#   )
#   with_parameters_test_that("multiple estimations with fixed effects works properly for femlm",
#                                      {
#                                        fmlas <- xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs])
#                                        if (s == "all") {
#                                          res <- fixest_mod_select(model = Model, fmla = fmlas, base = base, famly = fmly, notes = FALSE)
#                                        } else {
#                                          base_aux <- base[species_indx[[s]], ]
#                                          res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#                                        }
#
#                                        est <- est_multi[[num_fmla]]
#
#                                        vname <- names(coef(res))
#                                        expect_equal2(coef(est_multi[[num_fmla]])[vname], coef(res), tolerance = 1e-6)
#                                        expect_equal2(se(est_multi[[num_fmla]], cluster = "fe3")[vname], se(res, cluster = "fe3"), tolerance = 1e-6)
#                                      },
#                                      .cases = multiple_cases2(Model, famly = famlys[l])
#   )
# }
#
# ##### feNmlm
# Model <- "femlm"
# famlys <- c("poisson", "gaussian", "negbin")
# for(l in 1:length(famlys)){
#   est_multi <- fixest_mod_select(
#     model = Model,
#     fmla = c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2,
#     base = base,
#     fsplit = ~species,
#     famly = famlys[l]
#   )
#   with_parameters_test_that("multiple estimations with fixed effects works properly for feNmlm",
#                                      {
#                                        fmlas <- xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs])
#                                        if (s == "all") {
#                                          res <- fixest_mod_select(model = Model, fmla = fmlas, base = base, famly = fmly, notes = FALSE)
#                                        } else {
#                                          base_aux <- base[species_indx[[s]], ]
#                                          res <- fixest_mod_select(model = Model, fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#                                        }
#
#                                        est <- est_multi[[num_fmla]]
#
#                                        vname <- names(coef(res))
#                                        expect_equal2(coef(est_multi[[num_fmla]])[vname], coef(res), tolerance = 1e-6)
#                                        expect_equal2(se(est_multi[[num_fmla]], cluster = "fe3")[vname], se(res, cluster = "fe3"), tolerance = 1e-6)
#                                      },
#                                      .cases = multiple_cases2(Model, famly = famlys[l])
#   )
# }
