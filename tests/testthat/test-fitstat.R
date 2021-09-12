test_that("fitstat works properly", {
  base <- iris
  names(base) <- c("y", "x1", "x_endo_1", "x_inst_1", "fe")
  set.seed(2)
  base$x_inst_2 <- 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
  base$x_endo_2 <- 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)

  # Checking a basic estimation
  est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
  fitstat(est_iv, ~ f + ivf + ivf2 + wald + ivwald + ivwald2 + wh + sargan + rmse + g + n + ll + sq.cor + r2)

  est_fe <- feols(y ~ x1 | fe, base)
  fitstat(est_fe, ~wf)

  expect_equal(1, 1) # Omit "Empty test" message
})
