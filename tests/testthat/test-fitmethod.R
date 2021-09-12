fitmethod.cases <- fitmethod_cases()[-c(4:6), ] # Eliminating ols with fmly (makes no sense)

with_parameters_test_that("feols.fit works properly",
  {
    fmla <- paste(y_dep, "-1 + x1 + x2 + x3", sep = " ~ ")
    res <- feols.fit(y = ev_par(paste0("base$", y_dep)), X = base[, 2:4])
    res_bis <- feols(fml = as.formula(fmla), data = base)

    expect_equal(coef(res), coef(res_bis))
  },
  .cases = fitmethod.cases[1:3, ]
)

with_parameters_test_that("feglm.fit works properly",
  {
    fmla <- paste(y_dep, "-1 + x1 + x2 + x3", sep = " ~ ")
    if (isTRUE(with_fmly)) {
      res <- feglm.fit(y = ev_par(paste0("base$", y_dep)), X = base[, 2:4], family = fmly)
      res_bis <- feglm(fml = as.formula(fmla), data = base, family = fmly)
    } else {
      res <- feglm.fit(y = ev_par(paste0("base$", y_dep)), X = base[, 2:4])
      res_bis <- feglm(fml = as.formula(fmla), data = base)
    }

    expect_equal(coef(res), coef(res_bis))
  },
  .cases = fitmethod.cases[4:6, ]
)
