base <- iris
names(base) <- c("y", "x1", "x2", "x3", "fe1")
base$fe2 <- rep(1:5, 30)
base$y[1:5] <- NA
base$x1[4:8] <- NA
base$x2[4:21] <- NA
base$x3[110:111] <- NA
base$fe1[110:118] <- NA
base$fe2[base$fe2 == 1] <- 0
base$fe3 <- sample(letters[1:5], 150, TRUE)
base$period <- rep(1:50, 3)
base$x_cst <- 1


test_that("there is no bug in corner cases for feols",
          {
              res <- feols(y ~ 1 | csw(fe1, fe1^fe2), base)
              res <- feols(y ~ 1 + csw(x1, i(fe1)) | fe2, base)
              res <- feols(y ~ csw(f(x1, 1:2), x2) | sw0(fe2, fe2^fe3), base, panel.id = ~ fe1 + period)
              res <- feols(d(y) ~ -1 + d(x2), base, panel.id = ~ fe1 + period)
              expect_equal(length(coef(res)), 1)
              res <- feols(c(y, x1) ~ 1 | fe1 | x2 ~ x3, base)
              res <- feols(y ~ x1 | fe1[x2] + fe2[x2], base)
          })
