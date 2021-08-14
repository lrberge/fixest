####
#### Lagging ####
####

# Different types of lag
# 1) check no error in wide variety of situations
# 2) check consistency

base = datab12()

#
# Checks consistency
#

test_that("fixtest::lag function works properly",{
    expect_equal(lag(x1 ~ id + period, data = base), base$x1_lag)
    expect_equal(lag(x1 ~ id + period, -1, data = base), base$x1_lead)

    expect_equal(lag(x1 ~ id + period_bis, data = base), base$x1_lag_hole)
    expect_equal(lag(x1 ~ id + period_bis, -1, data = base), base$x1_lead_hole)

    expect_equal(lag(x1 ~ id + period_txt, data = base), base$x1_lag)
    expect_equal(lag(x1 ~ id + period_txt, -1, data = base), base$x1_lead)

    expect_equal(lag(x1 ~ id + period_date, data = base), base$x1_lag)
    expect_equal(lag(x1 ~ id + period_date, -1, data = base), base$x1_lead)

})

for (depvar in c("y", "y_na", "y_0")) {
    for (p in c("period", "period_txt", "period_date")) {
        base$per <- base[[p]]

        cat(".")

        base$y_dep <- base[[depvar]]
        pdat <- panel(base, ~ id + period)

        if (depvar == "y_0") {
            estfun <- fepois
        } else {
            estfun <- feols
        }

        est_raw <- estfun(y_dep ~ x1 + x1_lag + x1_lead, base)
        est <- estfun(y_dep ~ x1 + l(x1) + f(x1), base, panel.id = "id,per")
        est_pdat <- estfun(y_dep ~ x1 + l(x1, 1) + f(x1, 1), pdat)
        test(coef(est_raw), coef(est))
        test(coef(est_raw), coef(est_pdat))


        # Now diff
        est_raw <- estfun(y_dep ~ x1 + x1_diff, base)
        est <- estfun(y_dep ~ x1 + d(x1), base, panel.id = "id,per")
        est_pdat <- estfun(y_dep ~ x1 + d(x1, 1), pdat)
        test(coef(est_raw), coef(est))
        test(coef(est_raw), coef(est_pdat))

        # Now we just check that calls to l/f works without checking coefs

        est <- estfun(y_dep ~ x1 + l(x1) + f(x1), base, panel.id = "id,per")
        est <- estfun(y_dep ~ l(x1, -1:1) + f(x1, 2), base, panel.id = c("id", "per"))
        est <- estfun(y_dep ~ l(x1, -1:1, fill = 1), base, panel.id = ~ id + per)
        if (depvar == "y") test(est$nobs, n)
        est <- estfun(f(y_dep) ~ f(x1, -1:1), base, panel.id = ~ id + per)
    }
}

lagging_cases()


patrick::with_parameters_test_that("lagging fits works",{
    base$per <- base[[p]]
    base$y_dep <- base[[depvar]]
    pdat <- panel(base, ~ id + period)

    est_raw <- estfun(y_dep ~ x1 + x1_lag + x1_lead, base)
    est <- estfun(y_dep ~ x1 + l(x1) + f(x1), base, panel.id = "id,per")
    est_pdat <- estfun(y_dep ~ x1 + l(x1, 1) + f(x1, 1), pdat)

    est_raw <- fixest_mod_select(model = method,fmla = y_dep ~ x1 + x1_lag + x1_lead, base = base)
    est <- fixest_mod_select(model = method,fmla = y_dep ~ x1 + l(x1) + f(x1), base = base, panel.id = "id,per")
    est_pdat <- estfun(y_dep ~ x1 + l(x1, 1) + f(x1, 1), pdat)



    expect_equal(coef(est_raw), coef(est))
    expect_equal(coef(est_raw), coef(est_pdat))

    # Now diff
    est_raw <- estfun(y_dep ~ x1 + x1_diff, base)
    est <- estfun(y_dep ~ x1 + d(x1), base, panel.id = "id,per")
    est_pdat <- estfun(y_dep ~ x1 + d(x1, 1), pdat)
    test(coef(est_raw), coef(est))
    test(coef(est_raw), coef(est_pdat))

    # Now we just check that calls to l/f works without checking coefs

    est <- estfun(y_dep ~ x1 + l(x1) + f(x1), base, panel.id = "id,per")
    est <- estfun(y_dep ~ l(x1, -1:1) + f(x1, 2), base, panel.id = c("id", "per"))
    est <- estfun(y_dep ~ l(x1, -1:1, fill = 1), base, panel.id = ~ id + per)
    if (depvar == "y") test(est$nobs, n)
    est <- estfun(f(y_dep) ~ f(x1, -1:1), base, panel.id = ~ id + per)
},
.cases = lagging_cases()
)



#
# Data table
#

cat("data.table...")
# We just check there is no bug (consistency should be OK)

library(data.table)

base_dt <- data.table(
    id = c("A", "A", "B", "B"),
    time = c(1, 2, 1, 3),
    x = c(5, 6, 7, 8)
)

base_dt <- panel(base_dt, ~ id + time)

base_dt[, x_l := l(x)]
test(base_dt$x_l, c(NA, 5, NA, NA))

lag_creator <- function(dt) {
    dt2 <- panel(dt, ~ id + time)
    dt2[, x_l := l(x)]
    return(dt2)
}

base_bis <- lag_creator(base_dt)

base_bis[, x_d := d(x)]

cat("done.\n\n")


####
#### predict ####
####

chunk("PREDICT")

base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$fe_bis <- sample(letters, 150, TRUE)

#
# Same generative data
#

# Predict with fixed-effects
res <- feols(y ~ x1 | species + fe_bis, base)
test(predict(res), predict(res, base))

res <- fepois(y ~ x1 | species + fe_bis, base)
test(predict(res), predict(res, base))

res <- femlm(y ~ x1 | species + fe_bis, base)
test(predict(res), predict(res, base))


# Predict with varying slopes -- That's normal that tolerance is high (because FEs are computed with low precision)
res <- feols(y ~ x1 | species + fe_bis[x3], base)
test(predict(res), predict(res, base), "~", tol = 1e-4)

res <- fepois(y ~ x1 | species + fe_bis[x3], base)
test(predict(res), predict(res, base), "~", tol = 1e-3)


# Prediction with factors
res <- feols(y ~ x1 + i(species), base)
test(predict(res), predict(res, base))

res <- feols(y ~ x1 + i(species) + i(fe_bis), base)
test(predict(res), predict(res, base))

quoi <- head(base[, c("y", "x1", "species", "fe_bis")])
test(head(predict(res)), predict(res, quoi))

quoi$species <- as.character(quoi$species)
quoi$species[1:3] <- "zz"
test(head(predict(res)), predict(res, quoi))

# prediction with lags
data(base_did)
res <- feols(y ~ x1 + l(x1), base_did, panel.id = ~ id + period)
test(predict(res, sample = "original"), predict(res, base_did))

qui <- sample(which(base_did$id %in% 1:5))
base_bis <- base_did[qui, ]
test(predict(res, sample = "original")[qui], predict(res, base_bis))

# prediction with poly
res_poly <- feols(y ~ poly(x1, 2), base)
pred_all <- predict(res_poly)
pred_head <- predict(res_poly, head(base, 20))
pred_tail <- predict(res_poly, tail(base, 20))
test(head(pred_all, 20), pred_head)
test(tail(pred_all, 20), pred_tail)

#
# "Predicting" fixed-effects
#


res <- feols(y ~ x1 | species^fe_bis[x2], base, combine.quick = FALSE)

obs_fe <- predict(res, fixef = TRUE)
fe_coef_all <- fixef(res, sorted = FALSE)

coef_fe <- fe_coef_all[[1]]
coef_vs <- fe_coef_all[[2]]

fe_names <- paste0(base$species, "_", base$fe_bis)

test(coef_fe[fe_names], obs_fe[, 1])
test(coef_vs[fe_names], obs_fe[, 2])
