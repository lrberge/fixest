Base <- datab9()
X <- Base[, c("ln_euros", "ln_dist")]
fe <- Base[, c("Origin", "Destination")]
base_new <- demean(X, fe)

a <- feols(ln_euros ~ ln_dist, base_new)
b <- feols(ln_euros ~ ln_dist | Origin + Destination, Base, demeaned = TRUE)

test_that("feols estimates correctly demean'd data", {
  expect_equal2(coef(a)[-1], coef(b), tolerance = 1e-12, scale = 1)
})

test_that("feols and demean produce same demean'd data", {
  expect_equal(base_new$ln_euros, as.numeric(b$y_demeaned))
  expect_equal(base_new$ln_dist, as.numeric(b$X_demeaned))
})

# Now we just check there's no error

# NAs
X_NA <- X
fe_NA <- fe
X_NA[1:5, 1] <- NA
fe_NA[6:10, 1] <- NA

test_that("demean removing NAs works properly", {
  X_demean <- demean(X_NA, fe_NA, na.rm = FALSE)
  expect_equal(nrow(X_demean), nrow(X))
})

# integer
X_int <- X
X_int[[1]] <- as.integer(X_int[[1]])
X_demean <- demean(X_int, fe)

test_that("demean's as.matrix", {
  X_demean1 <- demean(X_int, fe, as.matrix = TRUE)
  X_demean2 <- demean(as.matrix(X_int), fe, as.matrix = FALSE)
  expect_true(is.matrix(X_demean1))
  expect_true(!is.matrix(X_demean2))
})

test_that("demean with formula and slopes", {
  X_dm_slopes <- demean(ln_dist ~ Origin + Destination[ln_euros], data = Base)
  X_dm_slopes_bis <- demean(Base$ln_dist, fe, slope.vars = Base$ln_euros, slope.flag = c(0, 1))
  expect_equal(X_dm_slopes[[1]], as.numeric(X_dm_slopes_bis))
})
