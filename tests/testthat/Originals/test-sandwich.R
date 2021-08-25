
## http://thestatsgeek.com/2013/10/12/the-robust-sandwich-variance-estimator-for-linear-regression/

library(sandwich)

base <- datab4()

est_lm <- lm(y ~ x + as.factor(grp) + as.factor(tm), data = datab4())
est_feols <- feols(y ~ x | grp + tm, data = datab4())

######### Testing sandwich covariance estimations

### Clustered Standard Errors
testthat::test_that("Clustered Standard Errors are equal between lm and sandwich", {
  # Standard Errors
  vcov_fxt <- se(est_feols, se = "st")["x"]
  vcov_sw <- se(est_lm)["x"]
  testthat::expect_equal(vcov_fxt, vcov_sw)

  # Clustered SE type HC0
  vcov_fxt <- unname(se(est_feols, dof = dof(adj = FALSE, fixef.K = "full")))
  vcov_sw <- sqrt(vcovCL(est_lm, cluster = datab4()$grp, type = "HC0")["x", "x"])
  testthat::expect_equal(vcov_fxt, vcov_sw)

  # Clusteres SE type HC1
  vcov_fxt <- unname(se(est_feols, dof = dof(fixef.K = "full")))
  vcov_sw <- sqrt(vcovCL(est_lm, cluster = datab4()$grp, type = "HC1")["x", "x"])
  testthat::expect_equal(vcov_fxt, vcov_sw)
})

### Heteroskedasticity-robust
testthat::test_that("Heterscoedastcity-robus Standard Errors are equal between feols and sandwich", {
  vcov_fxt <- unname(se(est_feols, se = "hetero", dof = dof(adj = FALSE, cluster.adj = FALSE)))
  vcov_sw <- sqrt(vcovHC(est_lm, type = "HC0")["x", "x"])
  testthat::expect_equal(vcov_fxt, vcov_sw)

  vcov_fxt <- unname(se(est_feols, se = "hetero"))
  vcov_sw <- sqrt(vcovHC(est_lm, type = "HC1")["x", "x"])
  testthat::expect_equal(vcov_fxt, vcov_sw)
})

### Two-way clustered standard errors
testthat::test_that("Two-way Clustered Standard Errors are equal between feols and sandwich", {
  vcov_fxt <- unname(se(est_feols, se = "twoway", dof = dof(fixef.K = "full", cluster.df = "conv")))
  vcov_sw <- sqrt(vcovCL(est_lm, cluster = ~ grp + tm, type = "HC1")["x", "x"])
  testthat::expect_equal(vcov_fxt, vcov_sw)
})


## cluster = data vs cluster = formula
testthat::test_that("cluster = data vs cluster = formula works properly", {
  for (k in 1:3) {
    datab4() <- vcov_db(k)
    est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, datab4())
    expect_equal_vcov(est_pois)
  }
})


## Testing Errors
testthat::test_that("SE is reporting errors correctly", {
  base <- vcov_db(3)
  est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, vcov_db(3))

  testthat::expect_error(se(est_pois, cluster = "Origin_na"))
  testthat::expect_error(se(est_pois, cluster = vcov_db(3)$Origin_na))
  testthat::expect_error(se(est_pois, cluster = list(vcov_db(3)$Origin_na)))
  testthat::expect_error(se(est_pois, cluster = ~ Origin_na^Destination))
  testthat::expect_error(se(est_pois, se = "cluster", cluster = ~ Origin_na^not_there))
  testthat::expect_error(se(est_pois, se = "cluster", cluster = ~ Origin_na^not_there))
  testthat::expect_error(se(est_pois, se = "twoway", cluster = c("Origin^Destination", "Product", "error")))
  testthat::expect_error(se(est_pois, se = "twoway", cluster = vcov_db(3)[, 1:4]))
  testthat::expect_error(se(est_pois, se = "twoway", cluster = ~ Origin + Destination + Product))
  testthat::expect_error(se(est_pois, se = "fourway", cluster = ~ Origin + Destination + Product))
})

## Testing aliases

testthat::test_that("Aliases works properly", {
  est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, vcov_db(3))

  se_hetero <- se(est_pois, se = "hetero")
  se_white <- se(est_pois, se = "white")
  # se_hc1    = se(est_pois, se = "hc1") # Hc1 gives error -> no matches found for hc1

  testthat::expect_equal(se_hetero, se_white)
  # testthat::expect_equal(se_hetero, se_hc1)
})



## Testing compatibility between fixest and sandwich

## Only vcovCL is tested. Extend to se = hetero, twoway, threeway
## Parametrize for future extension

data(base_did)
base_did <- base_did
base_did$y_int <- as.integer(base_did$y) + 20

patrick::with_parameters_test_that("fixest is compatible with sandwich's vcov",
  {
    est <- fixest_mod_select(
      model = model,
      fmla = fmlas,
      base = base,
      famly = family,
      weights = NULL
    )

    if (isFALSE(FE)) {
      # test(vcov(est, cluster = ~id), vcovCL(est, cluster = ~id, type = "HC1"))
      expect_equal(vcov(est, cluster = ~id), vcovCL(est, cluster = ~id, type = "HC1"))
    } else {
      # test(vcov(est, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est, cluster = ~id))
      expect_equal(vcov(est, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est, cluster = ~id))
    }
  },
  .cases = sandwcomp_cases()
)

# K = 1
# casos = sandwcomp_cases()[K,]
# model = casos$model
# fmlas = casos$fmlas
# fmla = as.formula(fmlas)
# # data = base
# famly = casos$family
# Data = base
# est <- fixest_mod_select(
#   model = model,
#   fmla = fmlas,
#   base = base,
#   famly = famly,
#   weights = NULL
# )
#
# expect_equal(vcov(est, cluster = ~id), vcovCL(est, cluster = ~id, type = "HC1"))
# expect_equal(vcov(est, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est, cluster = ~id))
#
# ## I must call my database with the same name of fixest_mod_select Data argument
# ## For some reason vcovCL doesnt work without the previous detail
