base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$fe_2 <- round(rnorm(150))

test_that("the following works without error", {
  m <- feols(y ~ x1 + i(fe_2), base)
  coefplot(m)
  etable(m, dict = c("0" = "zero"))

  m <- feols(y ~ x1 + i(fe_2) + i(fe_2, x2), base)
  coefplot(m)
  etable(m, dict = c("0" = "zero"))

  a <- i(base$fe_2)
  b <- i(base$fe_2, ref = 0:1)
  d <- i(base$fe_2, keep = 0:1)

  expect_equal(ncol(a), ncol(b) + 2)
  expect_equal(ncol(d), 2)
})
dev.off()
