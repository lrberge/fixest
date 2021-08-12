## Ctrl + Alt + R to run all the script


expect_ols_equal <- function(object, reference) {
  testthat::test_that("Equal x1 coefficient between feols and lm", {
    local_edition(2)
    testthat::expect_equal(
      coef(object)["x1"],
      coef(reference)["x1"],
      scale = 1
    )
  })

  testthat::test_that("Equal standard errors between feols and lm", {
    adj <- 1
    local_edition(2)
    testthat::expect_equal(
      se(object, se = "st", dof = dof(adj = adj))["x1"],
      se(reference)["x1"],
      scale = 1
    )
  })

  testthat::test_that("Equal p-values between feols and lm", {
    adj <- 1
    local_edition(2)
    testthat::expect_equal(
      pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
      pvalue(reference)["x1"],
      scale = 1
    )
  })
}


### expect function for glms
expect_glm_equal <- function(object, reference) {
  model <- reference$family$family
  tol <- switch(model,
    "binomial" = 3e-5,
    1e-5
  )
  adj <- 0

  if (model == "binomial" & (reference$formula == (y_01 ~ x1 + species + i(species, x2) + factor(fe_2) + i(fe_2, x3) + factor(fe_3)) |
    reference$formula == (y_01 ~ x1 + species + factor(fe_2) + i(fe_2, x2) + i(fe_2, x3) + factor(fe_3)))) {
    tol <- 0.5
  }

  testthat::test_that("feglm and glm have equal x1 coefficient", {
    local_edition(2)
    expect_equal(coef(object)["x1"],
      coef(reference)["x1"],
      tolerance = tol,
      scale = 1 # Absolute difference
    )
  })

  testthat::test_that("feglm and glm have equal standard errors", {
    local_edition(2)
    expect_equal(se(object, se = "st", dof = dof(adj = adj))["x1"],
      se(reference)["x1"],
      tolerance = tol,
      scale = 1 # Absolute difference
    )
  })

  testthat::test_that("Equal p-values between feglm and glm", {
    local_edition(2)
    expect_equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
      pvalue(reference)["x1"],
      tolerance = tol,
      scale = 1 # Absolute difference
    )
  })
}


### Expectation function for Negbin models

expect_negbin_equal <- function(object, reference) {
  tol <- 1e-2
  adj <- 1

  testthat::test_that("Equal x1 coefficient between fenegbin and MASS::glm.nb", {
    local_edition(2)
    expect_equal(coef(object)["x1"],
      coef(reference)["x1"],
      tolerance = tol,
      scale = 1
    )
  })

  testthat::test_that("Equal standard errors between fenegbin and MASS::glm.nb", {
    local_edition(2)
    expect_equal(se(object, se = "st", dof = dof(adj = adj))["x1"],
      se(reference)["x1"],
      tolerance = tol,
      scale = 1
    )
  })

  testthat::test_that("Equal p-values between between fenegbin and MASS::glm.nb", {
    local_edition(2)
    expect_equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
      pvalue(reference)["x1"],
      tolerance = tol * 10,
      scale = 1
    )
  })
}


expect_model_nointercept <- function(res) {
  testthat::test_that("feols omits correctly the intercept", {
    testthat::expect_false("(Intercept)" %in% names(res$coefficients))
  })
}


## Helper function for test-sandwich.R

expect_equal_vcov <- function(est_pois) {
  testthat::test_that("cluster = data vs cluster = formula works properly", {
    ## Clustered SE
    vcov1 <- se(est_pois, se = "cluster", cluster = "Product")
    vcov2 <- se(est_pois, cluster = trade$Product)
    vcov3 <- se(est_pois, cluster = ~Product)
    testthat::expect_equal(vcov1, vcov2)
    testthat::expect_equal(vcov1, vcov3)

    ## Two-way SE
    vcov1 <- se(est_pois, se = "twoway", cluster = trade[, c("Product", "Destination")])
    vcov2 <- se(est_pois, cluster = c("Product", "Destination"))
    vcov3 <- se(est_pois, cluster = ~ Product + Destination)
    testthat::expect_equal(vcov1, vcov2)
    testthat::expect_equal(vcov1, vcov3)

    ## Combined clusters
    vcov1 <- se(est_pois, cluster = "Product^Destination")
    vcov2 <- se(est_pois, cluster = paste(trade$Product, trade$Destination))
    vcov3 <- se(est_pois, cluster = ~ Product^Destination)
    testthat::expect_equal(vcov1, vcov2)
    testthat::expect_equal(vcov1, vcov3)

    ## Two-way Combined clusters
    vcov1 <- se(est_pois, cluster = c("Origin^Destination", "Product"))
    vcov2 <- se(est_pois, cluster = list(paste(trade$Origin, trade$Destination), trade$Product))
    vcov3 <- se(est_pois, cluster = ~ Origin^Destination + Product)
    testthat::expect_equal(vcov1, vcov2)
    testthat::expect_equal(vcov1, vcov3)
  })
}
