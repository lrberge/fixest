
library(fixest)
<<<<<<< HEAD
setFixest_notes(FALSE)
=======
base <- datab()
>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110

base <- datab()
patrick::with_parameters_test_that("feols and lm coefficients, Standard errors and p-values are equall",
  {
<<<<<<< HEAD
    res = fixest::feols(as.formula(fml_fixest), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
    res_bis = lm(as.formula(fml_stats), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
=======
    res <- fixest::feols(as.formula(fml_fixest), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
    res_bis <- lm(as.formula(fml_stats), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110
    expect_ols_equal(res, res_bis)
    expect_equal(1, 1)
  },
  .cases = ols_cases()
)

patrick::with_parameters_test_that("feglm and glm coefficients, Standard errors and p-values are equall",
  {
    res <- feglm(as.formula(fml_fixest), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
<<<<<<< HEAD
    if (!is.null(res$obs_selection$obsRemoved)) {
      qui <- res$obs_selection$obsRemoved
      base <- base[qui, ]
      res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight)[qui], offset = ev_par(my_offset)[qui])
    } else {
      res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
    }
=======
    res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110
    expect_glm_equal(res, res_bis)
    expect_equal(1, 1)
  },
  .cases = feglm_cases()
)

<<<<<<< HEAD
=======
#### fixest_test.R line 126 occurs a filtering of the database, which is included
#### in the following test, but there are errors in certain cases (e.g. feglm_cases()[19:20,]), which
#### doesnt occur in the former test

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
    expect_glm_equal(res, res_bis)
    expect_equal(1, 1)
  },
  .cases = feglm_cases()
)


>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110
patrick::with_parameters_test_that("fenegbin and glm.nb coefficients, Standard errors and p-values are equall",
  {
    res <- fenegbin(as.formula(fml_fixest), base, notes = FALSE)
    res_bis <- MASS::glm.nb(as.formula(fml_stats), base)
    expect_negbin_equal(res, res_bis)
    # expect_equal(1,1)
  },
  .cases = fenegbin_cases()
  # .cases = feglm_cases("negbin") # Test all glm cases
)
