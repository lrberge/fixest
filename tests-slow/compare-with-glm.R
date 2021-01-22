# NEVER PUT THIS ON TESTS/ FOLDER
# THE DATA IS TOO LARGE AND TAKES +5 MINUTES TO FIT THE MODEL

library(dplyr)
library(testthat)
library(fixest)
library(yotover)
library(microbenchmark)

# data ----

ch1_application1_2 <-  yotov_data("ch1_application1") %>%
    filter(year %in% seq(1986, 2006, 4))

ch1_application1_2 <- ch1_application1_2 %>%
    mutate(
        log_trade = log(trade),
        log_dist = log(dist)
    )

ch1_application1_2 <- ch1_application1_2 %>%
    # Create Yit
    group_by(exporter, year) %>%
    mutate(
        y = sum(trade),
        log_y = log(y)
    ) %>%

    # Create Eit
    group_by(importer, year) %>%
    mutate(
        e = sum(trade),
        log_e = log(e)
    )

ch1_application1_2 <- ch1_application1_2 %>%
    # Replicate total_e
    group_by(exporter, year) %>%
    mutate(total_e = sum(e)) %>%
    group_by(year) %>%
    mutate(total_e = max(total_e)) %>%

    # Replicate rem_exp
    group_by(exporter, year) %>%
    mutate(
        remoteness_exp = sum(dist *  total_e / e),
        log_remoteness_exp = log(remoteness_exp)
    ) %>%

    # Replicate total_y
    group_by(importer, year) %>%
    mutate(total_y = sum(y)) %>%
    group_by(year) %>%
    mutate(total_y = max(total_y)) %>%

    # Replicate rem_imp
    group_by(importer, year) %>%
    mutate(
        remoteness_imp = sum(dist / (y / total_y)),
        log_remoteness_imp = log(remoteness_imp)
    )

ch1_application1_2 <- ch1_application1_2 %>%
    # This merges the columns exporter/importer with year
    mutate(
        exp_year = paste0(exporter, year),
        imp_year = paste0(importer, year)
    )

ch1_application1_2 <- ch1_application1_2 %>%
    filter(exporter != importer)

# models ----

m1 <- glm(trade ~ log_dist + cntg + lang + clny + exp_year + imp_year,
          family = quasipoisson(link = "log"),
          data = ch1_application1_2,
          y = FALSE,
          model = FALSE
)

m2 <- feglm(trade ~ log_dist + cntg + lang + clny + exp_year + imp_year,
            family = quasipoisson(link = "log"),
            data = ch1_application1_2
)

# tests working ----

# altered current to avoid
# Error: m1$coefficients not equal to m2$coefficients.
# Lengths differ: 831 is not 826
expect_equal(m1$coefficients[!is.na(m1$coefficients)], m2$coefficients)

# altered current to avoid
# Error: m1$residuals not equal to m2$residuals.
# names for target but not for current
expect_equal(as.numeric(m1$fitted.values), m2$fitted.values)
expect_equal(as.numeric(m1$linear.predictors), m2$linear.predictors)

expect_equal(m1$family$family, m2$family$family)
expect_equal(m1$family$link, m2$family$link)

expect_equal(m1$deviance, m2$deviance)

expect_equal(m1$y, m2$y)

expect_equal(m1$call$family, m2$call$family)

expect_equal(m1$call$data, m2$call$data)

# TESTS FAILING ----

expect_equal(as.numeric(m1$residuals), m2$residuals)
expect_equal(m1$aic, m2$aic)
expect_equal(m1$null.deviance, m2$null.deviance)
expect_equal(m1$df.residual, m2$df.residual)
expect_equal(m1$df.null, m2$df.null)
expect_equal(m1$call$formula, m2$call$formula)

# benchmark ----

microbenchmark(
    glm(trade ~ log_dist + cntg + lang + clny + exp_year + imp_year,
        family = quasipoisson(link = "log"),
        data = ch1_application1_2,
        y = FALSE,
        model = FALSE
    ),
    feglm(trade ~ log_dist + cntg + lang + clny + exp_year + imp_year,
           family = quasipoisson(link = "log"),
           data = ch1_application1_2
    ),
    times = 3L
)

# Unit: seconds
# expr
# glm(trade ~ log_dist + cntg + lang + clny + exp_year + imp_year,      family = quasipoisson(link = "log"), data = ch1_application1_2,      y = FALSE, model = FALSE)
# feglm(trade ~ log_dist + cntg + lang + clny + exp_year + imp_year,      family = quasipoisson(link = "log"), data = ch1_application1_2)
# min        lq      mean   median         uq        max neval
# 89.247475 94.063727 97.098759 98.87998 101.024401 103.168824     3
# 4.199365  4.300148  4.420545  4.40093   4.531135   4.661339     3
