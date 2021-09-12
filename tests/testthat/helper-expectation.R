## Ctrl + Alt + R to run all the script

# 2nd edition of expect equal allows to use scale argument to switch between absolute and relative differences
expect_equal2 <- function(object, expected, tolerance = if (edition_get() >= 3) testthat_tolerance(), scale = NULL) {
  local_edition(2)
  expect_equal(object, expected, tolerance = tolerance, scale = scale)
}

expect_model_equal <- function(object, reference, method) {
  model <- ifelse(method == "ols" | method == "feNmlm", method, reference$family$family)
  tol <- ifelse(model == "binomial", 3e-5, 1e-5)
  tol <- ifelse(method == "negbin", 1e-2, tol)
  tol <- ifelse(model == "quasibinomial", 1e-2, tol)

  if (model == "binomial") {
    tol <- ifelse((reference$formula == (y_01 ~ x1 + species + i(species, x2) + factor(fe_2) + i(fe_2, x3) + factor(fe_3)) |
      reference$formula == (y_01 ~ x1 + species + factor(fe_2) + i(fe_2, x2) + i(fe_2, x3) + factor(fe_3))), 0.5, tol)
  }

  do_adj <<- ifelse(method == "glm", FALSE, TRUE)

  # fixest and stats have equal x1 coefficient
  expect_equal2(
    coef(object)["x1"],
    coef(reference)["x1"],
    tolerance = tol,
    scale = 1 # Absolute difference
  )

  # fixest and stats have equal standard errors
  expect_equal2(
    se(object, se = "st", ssc = ssc(adj = do_adj))["x1"],
    se(reference)["x1"],
    tolerance = tol,
    scale = 1 # Absolute difference
  )

  # fixest and stats have equal p-values
  expect_equal2(
    pvalue(object, se = "st", ssc = ssc(adj = do_adj))["x1"],
    pvalue(reference)["x1"],
    tolerance = tol,
    scale = 1 # Absolute difference
  )

  if (!is.null(object$dispersion)) {
    # fixest and stats have the same dispersion parameter
    if (!is.null(reference$family$family)) {
      tol <- ifelse(grepl("quasi", reference$family$family), 1e-3, tol)
    }
    expect_equal2(object$dispersion, summary(reference)$dispersion, tolerance = tol, scale = 1)
  }
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


## Helper function for test-fixef.R
expect_fixef <- function(all_coef, m_fe, k, str1, str2) {
  M <- switch(as.character(k),
    "1" = 2,
    "2" = 3,
    "3" = 3,
    "4" = 4,
    "5" = 4
  )
  for (i in 1:M) {
    Cs <- get_coef(all_coef, str2[i])
    test(var(Cs - m_fe[[str1[i]]][names(Cs)]), 0, "~")
  }
}

expect_fixef <- function(all_coef, m_fe, k, str1, str2) {
  M <- switch(as.character(k),
    "1" = 2,
    "2" = 3,
    "3" = 3,
    "4" = 4,
    "5" = 4
  )
  for (i in 1:M) {
    Cs <- get_coef(all_coef, str2[i])
    expect_equal2(var(Cs - m_fe[[str1[i]]][names(Cs)]), 0, scale = 1)
  }
}

### Function for test-nonlinear.R
# a and b are te parameters to be estimated
fun_nl <- function(a, b, tab2) {
  res <- as.numeric(tab2)
  return(a * res + b * res^2)
}
