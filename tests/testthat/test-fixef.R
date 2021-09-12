base <- datab6()
w <- 1
AuxL1 <- fixef.strings()[[1]]
AuxL2 <- fixef.strings()[[2]]
with_parameters_test_that("manually extracting fixed effects coefficients are equal to nlme::fixef",
  {
    if (K == 5) w <- 3 * (as.integer(base$species) - 0.95)
    m1 <- feols(ev_par(formula1), base, weights = w)
    m2 <- feols(ev_par(formula2), base, weights = w)
    m_fe <- fixef(m1)
    all_coef <- coef(m2)

    test_that("extracted fixed effects are equal to nlme's fixef", {
      expect_fixef(all_coef, m_fe, K, str1 = AuxL1[[K]], str2 = AuxL2[[K]])
    })
    expect_equal(1, 1)
  },
  .cases = fixef_cases()
)
