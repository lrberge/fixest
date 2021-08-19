
test_that("only.env works properly (without error)", {
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")

  env <- feols(y ~ x1 + x2 | species, base, only.env = TRUE)
  feols(env = env)

  env <- feglm(y ~ x1 + x2 | species, base, only.env = TRUE)
  feglm(env = env)

  env <- fepois(y ~ x1 + x2 | species, base, only.env = TRUE)
  fepois(env = env)

  env <- fenegbin(y ~ x1 + x2 | species, base, only.env = TRUE)
  fenegbin(env = env)

  env <- femlm(y ~ x1 + x2 | species, base, only.env = TRUE)
  femlm(env = env)

  env <- feNmlm(y ~ x1 + x2 | species, base, only.env = TRUE)
  feNmlm(env = env)
})
