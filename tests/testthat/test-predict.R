base <- datab13()
with_parameters_test_that("feols, feglm and femlm predict correctly with different families",
  {
    fmla <- xpd(lhs ~ rhs, lhs = y_dep, rhs = fmlas)
    res <- fixest_mod_select(model = method, fmla = fmla, base = base, famly = fmly)
    pred_res <- as.numeric(predict(res))
    pred_resb <- predict(res, base)
    expect_equal2(pred_res, pred_resb, tolerance = tol)
  },
  .cases = predict_cases()
)

test_that("Predict with factors works properly", {
  res <- feols(y ~ x1 + i(species) + i(fe_bis), base)
  expect_equal(predict(res), predict(res, base))

  quoi <- head(base[, c("y", "x1", "species", "fe_bis")])
  # test(head(predict(res)), predict(res, quoi))
  expect_equal(head(predict(res)), predict(res, quoi))
  quoi$species <- as.character(quoi$species)
  quoi$species[1:3] <- "zz"
  # test(head(predict(res)), predict(res, quoi))
  expect_equal(head(predict(res)), predict(res, quoi))
})

## HR: predict(sample = "original") doesnt work anymore. sample is not an argument of predict
## The following test may not be evaluating if prediction works properly
## due to the treatment needed with na.omits to obtain the TRUE's
test_that("prediction with lags works properly", {
  data(base_did)
  res <- feols(y ~ x1 + l(x1), base_did, panel.id = ~ id + period)
  pred_1 <- predict(res)
  pred_2 <- as.numeric(na.omit(predict(res, base_did))) # <<---
  expect_equal(pred_1, pred_2)

  qui <- sample(which(base_did$id %in% 1:5))
  base_bis <- base_did[qui, ]
  pred_1 <- predict(res)
  pred_2 <- as.numeric(na.omit(predict(res, base_bis)))

  expect_equal(pred_2 %in% pred_1, rep(TRUE, length(pred_2)))
})


test_that("prediction with poly() works properly", {
  res_poly <- feols(y ~ poly(x1, 2), base)
  pred_all <- predict(res_poly)
  pred_head <- predict(res_poly, head(base, 20))
  pred_tail <- predict(res_poly, tail(base, 20))
  expect_equal(head(pred_all, 20), pred_head)
  expect_equal(tail(pred_all, 20), pred_tail)
})

#
# "Predicting" fixed-effects
#

### HR: predict(fixef = TRUE) doesnt work anymore
#
# test_that(" 'predicting' fixed effects works properly", {
#   res <- feols(y ~ x1 | species^fe_bis[x2], base, combine.quick = FALSE)
#
#   obs_fe <- predict(res, fixef = TRUE)
#   fe_coef_all <- fixef(res, sorted = FALSE)
#
#   coef_fe <- fe_coef_all[[1]]
#   coef_vs <- fe_coef_all[[2]]
#
#   fe_names <- paste0(base$species, "_", base$fe_bis)
#
#   expect_equal(coef_fe[fe_names], obs_fe[, 1])
#   expect_equal(coef_vs[fe_names], obs_fe[, 2])
# })
