Bases = datab16()
base = Bases[[1]]
base_bis = Bases[[2]]

test_that("model.matrix removes NA correctly",
          {
              res <- feols(y1 ~ x1 + x2 + x3, base)
              m1 <- model.matrix(res, type = "lhs")
              expect_equal(length(m1), res$nobs)

              m1_na <- model.matrix(res, type = "lhs", na.rm = FALSE)
              expect_equal(length(m1_na), res$nobs_origin)
              expect_equal(max(abs(m1_na - base$y1), na.rm = TRUE), 0)

              y <- model.matrix(res, type = "lhs", data = base, na.rm = FALSE)
              X <- model.matrix(res, type = "rhs", data = base, na.rm = FALSE)
              obs_rm <- res$obs_selection$obsRemoved
              res_bis <- lm.fit(X[obs_rm, ], y[obs_rm])
              expect_equal(res_bis$coefficients, res$coefficients)
          })
test_that("model.matrix works properly with lagged variables",
          {
              res_lag <- feols(y1 ~ l(x1, 1:2) + x2 + x3, base, panel = ~ id + time)
              m_lag <- model.matrix(res_lag)
              expect_equal(nrow(m_lag), nobs(res_lag))

              # Laged subset
              m_lag_x1 <- model.matrix(res_lag, subset = "x1")
              expect_equal(ncol(m_lag_x1), 2)

              # lag with subset, new data
              mbis_lag_x1 <- model.matrix(res_lag, base_bis[, c("x1", "x2", "id", "time")], subset = TRUE)
              # l(x1, 1) + l(x1, 2) + x2
              expect_equal(ncol(mbis_lag_x1), 3)
              # 13 NAs: 2 per ID for the lags, 3 for x2
              expect_equal(nrow(mbis_lag_x1), 37)
          }
          )
test_that("model.matrix works properly with poly()",
          {
              res_poly <- feols(y1 ~ poly(x1, 2), base)
              m_poly_old <- model.matrix(res_poly)
              m_poly_new <- model.matrix(res_poly, base_bis)
              expect_equal(m_poly_old[1:50, 3], m_poly_new[, 3])
          }
          )
test_that("model.matrix works proprly with fixef",
          {
              res <- feols(y1 ~ x1 + x2 + x3 | species + fe2, base)
              m_fe <- model.matrix(res, type = "fixef") # model has fixed effects but model.matrix is not recognizing it
              expect_equal(ncol(m_fe), 2)
          })
test_that("type lhs argument from model.matrix works properly",
          {
              m_lhs <- model.matrix(res, type = "lhs", na.rm = FALSE)
              expect_equal(m_lhs, base$y1)
          })
test_that("model.matrix works properly with a IV model",
          {
              res_iv <- feols(y1 ~ x1 | x2 ~ x3, base)

              m_rhs1 <- model.matrix(res_iv, type = "iv.rhs1")
              expect_equal(colnames(m_rhs1)[-1], c("x3", "x1"))

              m_rhs2 <- model.matrix(res_iv, type = "iv.rhs2")
              expect_equal(colnames(m_rhs2)[-1], c("fit_x2", "x1"))

              m_endo <- model.matrix(res_iv, type = "iv.endo")
              expect_equal(colnames(m_endo), "x2")

              m_exo <- model.matrix(res_iv, type = "iv.exo")
              expect_equal(colnames(m_exo)[-1], "x1")

              m_inst <- model.matrix(res_iv, type = "iv.inst")
              expect_equal(colnames(m_inst), "x3")

              # Several
              res_mult <- feols(y1 ~ x1 | species | x2 ~ x3, base)
              m_lhs_rhs_fixef <- model.matrix(res_mult, type = c("lhs", "iv.rhs2", "fixef"), na.rm = FALSE) #Again model.matrix is not working
              expect_equal(names(m_lhs_rhs_fixef), c("y1", "fit_x2", "x1", "species"))
          })
