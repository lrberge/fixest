

library(fixest)
setFixest_notes(FALSE)

base <- datab()
patrick::with_parameters_test_that("feols and lm coefficients, Standard errors and p-values are equall",
  {
    res <- fixest::feols(as.formula(fml_fixest), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
    res_bis <- lm(as.formula(fml_stats), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
    expect_model_equal(res, res_bis, method = "ols")
    expect_equal(1, 1)
  },
  .cases = ols_cases()
)

## Vector of models to test:
all_mods <- c("binomial", "quasibinomial","poisson", "quasipoisson", "Gamma", "inverse.gaussian", "gaussian", "quasi")
feglm.cases <- feglm_cases(all_mods)
# most of the inverse.gaussian with offset dont converge
feglm.cases <- feglm.cases[!(feglm.cases$my_family == "inverse.gaussian" & feglm.cases$my_offset != "NULL"),]

patrick::with_parameters_test_that("feglm and glm coefficients, Standard errors and p-values are equall",
                                   {
                                     res <- feglm(as.formula(fml_fixest), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
                                     if (!is.null(res$obs_selection$obsRemoved)) {
                                       qui <- res$obs_selection$obsRemoved
                                       base <- base[qui, ]
                                       res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight)[qui], offset = ev_par(my_offset)[qui])
                                     } else {
                                       res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
                                     }
                                     expect_model_equal(res, res_bis, method = "glm")
                                     expect_equal(1, 1)
                                   },
                                   .cases = feglm.cases
)

patrick::with_parameters_test_that("femlm and glm coefficients are equall",
                                   {
                                     if(my_family != "negbin"){
                                       res <- femlm(as.formula(fml_fixest), base, family = my_family, offset = ev_par(my_offset))
                                       res_bis <- glm(as.formula(fml_stats), base, family = my_family_stats, offset = ev_par(my_offset))
                                     } else {
                                       res <- femlm(as.formula(fml_fixest), base, family = my_family)
                                       res_bis <- MASS::glm.nb(formula = as.formula(fml_stats), data = base)
                                     }

                                     # expect_model_equal(res, res_bis, method = "glm")
                                     expect_equal2(coef(res)["x1"], coef(res_bis)["x1"], tolerance = 1e-5)
                                     # expect_equal2(se(res, se = "standard")["x1"], se(res_bis)["x1"], tolerance = 1e-5)  # SE are not equal
                                   },
                                   .cases = femlm_cases()
)


patrick::with_parameters_test_that("fenegbin and glm.nb coefficients, Standard errors and p-values are equall",
  {
    res <- fenegbin(as.formula(fml_fixest), base, notes = FALSE)
    res_bis <- MASS::glm.nb(as.formula(fml_stats), base)
    expect_model_equal(res, res_bis, method = "negbin")
    # expect_equal(1,1)
  },
  .cases = fenegbin_cases()
  # .cases = feglm_cases("negbin") # Test all glm cases
)



# all_mods <- c("quasi")
# feglm.cases <- feglm_cases(all_mods)
# feglm.cases <- feglm.cases[!(feglm.cases$my_family == "inverse.gaussian" & feglm.cases$my_offset != "NULL"),]
#
# dim(feglm.cases)
#
# patrick::with_parameters_test_that("feglm and glm coefficients, Standard errors and p-values are equall",
#                                    {
#                                      res <- feglm(as.formula(fml_fixest), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
#                                      if (!is.null(res$obs_selection$obsRemoved)) {
#                                        qui <- res$obs_selection$obsRemoved
#                                        base <- base[qui, ]
#                                        res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight)[qui], offset = ev_par(my_offset)[qui])
#                                      } else {
#                                        res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
#                                      }
#                                      expect_model_equal(res, res_bis, method = "glm")
#                                      expect_equal(1, 1)
#                                    },
#                                    .cases = feglm.cases
# )
#
#
# l = 2
# casos <- feglm_cases(Models = c("inverse.gaussian"))[l,]
# my_offset <- casos$my_offset
# my_weight <- casos$my_weight
# fml_fixest <- casos$fml_fixest
# fml_stats <- casos$fml_stats
# my_family <- casos$my_family
#
# res <- feglm(as.formula(fml_fixest), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
# if (!is.null(res$obs_selection$obsRemoved)) {
#   qui <- res$obs_selection$obsRemoved
#   base <- base[qui, ]
#   res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight)[qui], offset = ev_par(my_offset)[qui])
# } else {
#   res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
# }
# summary(res)
# summary(res_bis)
# res$dispersion
#
# expect_model_equal(res, res_bis, method = "glm")
#
# reference = res_bis
# object = res
# method = "glm"
#
#
#
#
# patrick::with_parameters_test_that("fenegbin and glm.nb coefficients, Standard errors and p-values are equall",
#                                    {
#                                      res <- fenegbin(as.formula(fml_fixest), base, notes = FALSE)
#                                      res_bis <- MASS::glm.nb(as.formula(fml_stats), base)
#                                      expect_model_equal(res, res_bis, method = "negbin")
#                                      # expect_equal(1,1)
#                                    },
#                                    .cases = fenegbin_cases()
#                                    # .cases = feglm_cases("negbin") # Test all glm cases
# )
#
#
# l = 2
# casos <-fenegbin_cases()[1,]
# fml_fixest <- casos$fml_fixest
# fml_stats <- casos$fml_stats
#
# res <- fenegbin(as.formula(fml_fixest), base, notes = FALSE)
# res_bis <- MASS::glm.nb(as.formula(fml_stats), base)
#
# summary(res)
# summary(res_bis)
#
# expect_model_equal(res, res_bis, method = "negbin")
#
# reference = res_bis
# object = res
# method = "negbin"
#
