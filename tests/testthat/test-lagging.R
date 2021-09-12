####
#### Lagging ####
####

# Different types of lag
# 1) check no error in wide variety of situations
# 2) check consistency

setFixest_notes(FALSE)
base <- datab12()
n <- nrow(base)
#
# Checks consistency
#

test_that("fixtest::lag function works properly", {
  expect_equal(lag(x1 ~ id + period, data = base), base$x1_lag)
  expect_equal(lag(x1 ~ id + period, -1, data = base), base$x1_lead)

  expect_equal(lag(x1 ~ id + period_bis, data = base), base$x1_lag_hole)
  expect_equal(lag(x1 ~ id + period_bis, -1, data = base), base$x1_lead_hole)

  expect_equal(lag(x1 ~ id + period_txt, data = base), base$x1_lag)
  expect_equal(lag(x1 ~ id + period_txt, -1, data = base), base$x1_lead)

  expect_equal(lag(x1 ~ id + period_date, data = base), base$x1_lag)
  expect_equal(lag(x1 ~ id + period_date, -1, data = base), base$x1_lead)
})

with_parameters_test_that("fixest model fitting with panel data works properly",
  {
    base$per <- base[[p]]
    base$y_dep <- base[[depvar]]
    pdat <- panel(base, ~ id + period)
    est_raw <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + x1_lag + x1_lead, base = base, famly = fmly)
    est <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + l(x1) + f(x1), base = base, panel.id = "id,per", famly = fmly)
    est_pdat <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + l(x1, 1) + f(x1, 1), base = pdat, famly = fmly)

    expect_equal(unname(coef(est_raw)), unname(coef(est)))
    expect_equal(unname(coef(est_raw)), unname(coef(est_pdat)))
  },
  .cases = lagging_cases()
)


with_parameters_test_that("fixest model fitting with panel data works properly with differentiations",
  {
    base$per <- base[[p]]
    base$y_dep <- base[[depvar]]
    pdat <- panel(base, ~ id + period)

    est_raw <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + x1_diff, base = base, famly = fmly)
    est <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + d(x1), base = base, panel.id = "id,per", famly = fmly)
    est_pdat <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + d(x1, 1), base = pdat, famly = fmly)

    expect_equal(unname(coef(est_raw)), unname(coef(est)))
    expect_equal(unname(coef(est_raw)), unname(coef(est_pdat)))
  },
  .cases = lagging_cases()
)


with_parameters_test_that("fitting works with forward and lagging functions",
  {
    base$per <- base[[p]]
    base$y_dep <- base[[depvar]]
    pdat <- panel(base, ~ id + period)

    est <- fixest_mod_select(model = method, fmla = y_dep ~ x1 + l(x1) + f(x1), base = base, panel.id = "id,per", famly = fmly)
    est <- fixest_mod_select(model = method, fmla = y_dep ~ l(x1, -1:1) + f(x1, 2), base = base, panel.id = c("id", "per"), famly = fmly)
    est <- fixest_mod_select(model = method, fmla = y_dep ~ l(x1, -1:1, fill = 1), base = base, panel.id = ~ id + per, famly = fmly)
    if (depvar == "y") expect_equal(est$nobs, n)
    est <- fixest_mod_select(model = method, fmla = f(y_dep) ~ f(x1, -1:1), base = base, panel.id = ~ id + per, famly = fmly)
    expect_equal(1, 1)
  },
  .cases = lagging_cases()
)

# We just check there is no bug (consistency should be OK)
test_that("there is no bug using data.table for panel data", {
  base_dt <- data.table::data.table(
    id = c("A", "A", "B", "B"),
    time = c(1, 2, 1, 3),
    x = c(5, 6, 7, 8)
  )

  base_dt <- panel(base_dt, ~ id + time)

  base_dt[, x_l := l(x)]
  expect_equal(base_dt$x_l, c(NA, 5, NA, NA))

  lag_creator <- function(dt) {
    dt2 <- panel(dt, ~ id + time)
    dt2[, x_l := l(x)]
    return(dt2)
  }

  base_bis <- lag_creator(base_dt)
  base_bis[, x_d := d(x)]
})
