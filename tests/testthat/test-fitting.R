

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
  .cases = feglm_cases()
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
