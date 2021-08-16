## Ctrl + Alt + R to run all the script

library(fixest)

## Database used for test-fitting.R
datab <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$fe_2 <- rep(1:5, 30)
  base$fe_3 <- sample(15, 150, TRUE)
  base$constant <- 5
  base$y_int <- as.integer(base$y)
  base$w <- as.vector(unclass(base$species) - 0.95)
  base$offset_value <- unclass(base$species) - 0.95
  base$y_01 <- 1 * ((scale(base$x1) + rnorm(150)) > 0)
  # what follows to avoid removal of fixed-effects (logit is pain in the neck)
  base$y_01[1:5 + rep(c(0, 50, 100), each = 5)] <- 1
  base$y_01[6:10 + rep(c(0, 50, 100), each = 5)] <- 0
  # We enforce the removal of observations
  base$y_int_null <- base$y_int
  base$y_int_null[base$fe_3 %in% 1:5] <- 0
  return(base)
}

# database used for test-nointercept.R
datab2 <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "fe1")
  base$y_r <- round(base$y)
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
  return(base)
}

# Database used for test-vcov.R
datab3 <- function() {
  # We create "irregular" FEs
  set.seed(0)
  base <- data.frame(x = rnorm(20))
  base$y <- base$x + rnorm(20)
  base$fe1 <- rep(rep(1:3, c(4, 3, 3)), 2)
  base$fe2 <- rep(rep(1:5, each = 2), 2)
  return(base)
}

# Database used for test-sandwich.R
datab4 <- function() {
  set.seed(0)
  N <- 20
  G <- N / 5
  T <- N / G
  d <- data.frame(y = rnorm(N), x = rnorm(N), grp = rep(1:G, T), tm = rep(1:T, each = G))
  return(d)
}

datab5 <- function() {
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$y_int <- as.integer(base$y) + 1
  base$w <- as.numeric(base$species)
  return(base)
}

# database for test-fixef.R
datab6 <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$x4 <- rnorm(150) + 0.25 * base$y
  base$fe_bis <- sample(10, 150, TRUE)
  base$fe_ter <- sample(15, 150, TRUE)
  return(base)
}

# database for test-collinearity.R
datab7 <- function() {
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$constant <- 5
  base$y_int <- as.integer(base$y)
  base$w <- as.vector(unclass(base$species) - 0.95)
  return(base)
}

datab8 <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$z <- sample(5, 150, TRUE)
  base$species_na <- base$species
  base$species_na[base$species == "setosa"] <- NA
  return(base)
}

# Database for test-demean.R
datab9 <- function() {
  data(trade)
  base <- trade
  base$ln_euros <- log(base$Euros)
  base$ln_dist <- log(base$dist_km)
  return(base)
}

# Database for test-hatvalues.R
datab10 <- function() {
  set.seed(0)
  x <- sin(1:10)
  y <- rnorm(10)
  y_int <- rpois(10, 2)

  return(data.frame(y = y, y_int = y_int, x))
}

# Database for test-subset.R
datab11 <- function() {
  set.seed(5)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$fe_bis <- sample(letters, 150, TRUE)
  base$x4 <- rnorm(150)
  base$x1[sample(150, 5)] <- NA
  return(base)
}

# Database for test-lagging.R
datab12 <- function() {
  data(base_did)
  base <- base_did

  n <- nrow(base)

  set.seed(0)
  base$y_na <- base$y
  base$y_na[sample(n, 50)] <- NA
  base$period_txt <- letters[base$period]
  ten_dates <- c("1960-01-15", "1960-01-16", "1960-03-31", "1960-04-05", "1960-05-12", "1960-05-25", "1960-06-20", "1960-07-30", "1965-01-02", "2002-12-05")
  base$period_date <- as.Date(ten_dates, "%Y-%m-%d")[base$period]
  base$y_0 <- base$y**2
  base$y_0[base$id == 1] <- 0

  # We compute the lags "by hand"
  base <- base[order(base$id, base$period), ]
  base$x1_lag <- c(NA, base$x1[-n])
  base$x1_lag[base$period == 1] <- NA
  base$x1_lead <- c(base$x1[-1], NA)
  base$x1_lead[base$period == 10] <- NA
  base$x1_diff <- base$x1 - base$x1_lag

  # we create holes
  base$period_bis <- base$period
  base$period_bis[base$period_bis == 5] <- 50
  base$x1_lag_hole <- base$x1_lag
  base$x1_lag_hole[base$period %in% c(5, 6)] <- NA
  base$x1_lead_hole <- base$x1_lead
  base$x1_lead_hole[base$period %in% c(4, 5)] <- NA

  # we reshuffle the base
  base <- base[sample(n), ]

  return(base)
}

# database function for test-predict.R
datab13 <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$fe_bis <- sample(letters, 150, TRUE)
  return(base)
}

datab14 <- function() {
  set.seed(2)
  base <- iris
  names(base) <- c("y1", "x1", "x2", "x3", "species")
  base$y2 <- 10 + rnorm(150) + 0.5 * base$x1
  base$x4 <- rnorm(150) + 0.5 * base$y1
  base$fe2 <- rep(letters[1:15], 10)
  base$fe2[50:51] <- NA
  base$y2[base$fe2 == "a" & !is.na(base$fe2)] <- 0
  base$x2[1:5] <- NA
  base$x3[6] <- NA
  base$fe3 <- rep(letters[1:10], 15)
  return(base)
}

# Database function for test-IV
datab15 <- function() {
  base <- iris
  names(base) <- c("y", "x1", "x_endo_1", "x_inst_1", "fe")
  set.seed(2)
  base$x_inst_2 <- 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
  base$x_endo_2 <- 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)
  return(base)
}

# Database for test-model_matrix.R

datab16 <- function() {
  base <- iris
  names(base) <- c("y1", "x1", "x2", "x3", "species")
  base$y2 <- 10 + rnorm(150) + 0.5 * base$x1
  base$x4 <- rnorm(150) + 0.5 * base$y1
  base$fe2 <- rep(letters[1:15], 10)
  base$fe2[50:51] <- NA
  base$y2[base$fe2 == "a" & !is.na(base$fe2)] <- 0
  base$x2[1:5] <- NA
  base$x3[6] <- NA
  base$fe3 <- rep(letters[1:10], 15)
  base$id <- rep(1:15, each = 10)
  base$time <- rep(1:10, 15)

  base_bis <- base[1:50, ]
  base_bis$id <- rep(1:5, each = 10)
  base_bis$time <- rep(1:10, 5)
  return(list(base, base_bis))
}

ev_par <- function(string) {
  eval(parse(text = string))
}

fomla <- function(model, fml_fixest1, fml_stats1) {
  if (model == "binomial") {
    fml_fixest <- paste("y_01", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y_01", fml_stats1, sep = " ~ ")
  } else if (model == "poisson") {
    fml_fixest <- paste("y_int_null", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y_int_null", fml_stats1, sep = " ~ ")
  } else if (model %in% c("negbin", "Gamma")) {
    fml_fixest <- paste("y_int", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y_int", fml_stats1, sep = " ~ ")
  }
  return(list(fml_fixest, fml_stats))
}

ols_cases <- function() {
  fml_fixest <- c(
    "y ~ x1",
    "y ~ x1 | species",
    "y ~ x1 | species + fe_2",
    "y ~ x1 | species[[x2]]",
    "y ~ x1 | species[[x2]] + fe_2",
    "y ~ x1 | species[x2]",
    "y ~ x1 | species[[x2]] + fe_2[[x3]]",
    "y ~ x1 + x2 | species^fe_2",
    "y ~ x1 | species[x2] + fe_2[x3] + fe_3",
    "y ~ x1 | species + fe_2[x2,x3] + fe_3"
  )
  fml_stats <- c(
    "y ~ x1",
    "y ~ x1 + factor(species)",
    "y ~ x1 + factor(species) + factor(fe_2)",
    "y ~ x1 + x2:species",
    "y ~ x1 + x2:species + factor(fe_2)",
    "y ~ x1 + x2:species + species",
    "y ~ x1 + x2:species + x3:factor(fe_2)",
    "y ~ x1 + x2 + paste(species, fe_2)",
    "y ~ x1 + species + i(species, x2) + factor(fe_2) + i(fe_2, x3) + factor(fe_3)",
    "y ~ x1 + species + factor(fe_2) + i(fe_2, x2) + i(fe_2, x3) + factor(fe_3)"
  )

  DF_l <- list()
  cont <- 1
  for (weights in c("NULL", "base$w")) {
    for (offset in c("NULL", "base$offset_value")) {
      DF_l[[cont]] <- data.frame(
        fml_fixest = fml_fixest,
        fml_stats = fml_stats,
        my_offset = offset,
        my_weight = weights
      )
      cont <- cont + 1
    }
  }
  DF <- rbind(
    DF_l[[1]],
    DF_l[[2]],
    DF_l[[3]],
    DF_l[[4]]
  )
  DF$test_name <- 1:dim(DF)[1]
  return(DF)
}



feglm_cases <- function(Models = c("binomial", "poisson", "Gamma")) {
  fml_fixest1 <- c(
    "x1",
    "x1 | species",
    "x1 | species + fe_2",
    "x1 | species[[x2]]",
    "x1 | species[[x2]] + fe_2",
    "x1 | species[x2]",
    "x1 | species[[x2]] + fe_2[[x3]]",
    "x1 + x2 | species^fe_2",
    "x1 | species[x2] + fe_2[x3] + fe_3",
    "x1 | species + fe_2[x2,x3] + fe_3"
  )
  fml_stats1 <- c(
    "x1",
    "x1 + factor(species)",
    "x1 + factor(species) + factor(fe_2)",
    "x1 + x2:species",
    "x1 + x2:species + factor(fe_2)",
    "x1 + x2:species + species",
    "x1 + x2:species + x3:factor(fe_2)",
    "x1 + x2 + paste(species, fe_2)",
    "x1 + species + i(species, x2) + factor(fe_2) + i(fe_2, x3) + factor(fe_3)",
    "x1 + species + factor(fe_2) + i(fe_2, x2) + i(fe_2, x3) + factor(fe_3)"
  )

  DF_l <- list()
  cont <- 1
  for (model in Models) {

    ## Adding link function
    my_fam <- switch(model,
      binomial = "binomial",
      poisson = "poisson",
      negbin = "NULL",
      Gamma = "Gamma"
    )

    for (weights in c("NULL", "base$w")) {
      for (offset in c("NULL", "base$offset_value")) {
        ## Change the resposne variable accordingly to the model

        # fixest_test.R line 120. For some reason we are omitting gamma model with offset.
        if (model == "Gamma" && offset != "NULL") next # Deleting this line produce different results between feglm and glm

        aux <- fomla(model, fml_fixest1, fml_stats1)
        fml_fixest <- aux[[1]]
        fml_stats <- aux[[2]]

        DF_l[[cont]] <- data.frame(
          fml_fixest = fml_fixest,
          fml_stats = fml_stats,
          model = model,
          my_family = my_fam,
          my_offset = offset,
          my_weight = weights,
          test_name = paste(model, weights, offset, 1:length(fml_fixest), sep = "_")
        )
        cont <- cont + 1
      }
    }
  }
  DF <- DF_l[[1]]
  for (k in 2:length(DF_l)) DF <- rbind(DF, DF_l[[k]])
  return(DF)
}


# fixest_test.R doesnt include weights or more than one fixed effect
fenegbin_cases <- function() {
  fml_fixest1 <- c(
    "y_int ~ x1",
    "y_int ~ x1 | species"
  )
  fml_stats1 <- c(
    "y_int ~ x1",
    "y_int ~x1 + factor(species)"
  )
  DF <- data.frame(
    fml_fixest = fml_fixest1,
    fml_stats = fml_stats1,
    model = "negbin",
    test_name = paste("negbin", 1:length(fml_fixest1), sep = "_")
  )
  return(DF)
}


vcov_cases1 <- function() {
  DF <- data.frame(expand.grid(
    c(FALSE, TRUE),
    c("none", "nested", "full")
  ))
  names(DF) <- c("adj", "k_val")
  DF$K <- ifelse(DF$k_val == "none", 1, 8)
  DF$my_adj <- ifelse(DF$adj, (20 - 1) / (20 - DF$K), 1)

  DF$test_name <- paste(DF$adj, DF$k_val, sep = " - ")
  return(DF)
}

vcov_cases2 <- function() {
  DF <- data.frame(expand.grid(
    c("none", "nested", "full"),
    c("conventional", "min"),
    c(FALSE, TRUE),
    c(FALSE, TRUE)
  ))
  names(DF) <- c("k_val", "tdf", "c_adj", "adj")

  DF$K[DF$k_val == "none"] <- 1
  DF$K[DF$k_val == "nested"] <- 6
  DF$K[DF$k_val == "full"] <- 8
  DF$my_adj <- ifelse(DF$adj, (20 - 1) / (20 - DF$K), 1)
  DF$cluster_factor <- ifelse(DF$c_adj, 3 / 2, 1)
  DF$df <- ifelse(DF$tdf == "min", 2, 20 - 8)

  DF$test_name <- paste(DF$k_val, DF$tdf, DF$c_adj, DF$adj, sep = " - ")
  return(DF)
}

vcov_cases3 <- function() {
  DF <- data.frame(expand.grid(
    c("conventional", "min"),
    c("conventional", "min"),
    c("none", "nested", "full"),
    c(FALSE, TRUE),
    c(FALSE, TRUE)
  ))
  names(DF) <- c("cdf", "tdf", "k_val", "c_adj", "adj")

  DF$K[DF$k_val == "none"] <- 1
  DF$K[DF$k_val == "nested"] <- 2
  DF$K[DF$k_val == "full"] <- 8
  DF$my_adj <- ifelse(DF$adj, (20 - 1) / (20 - DF$K), 1)
  DF$cluster_factor <- ifelse(DF$c_adj, 3 / 2, 1)
  DF$df <- ifelse(DF$tdf == "min", 2, 20 - 8)

  DF$test_name <- paste(DF$cdf, DF$tdf, DF$k_val, DF$c_adj, DF$adj, sep = " - ")
  return(DF)
}

V_matrix <- function(Mi, M_t, M_it, c_adj, cdf) {
  if (c_adj) {
    if (cdf == "min") {
      V <- (M_i + M_t - M_it) * 3 / 2
    } else {
      V <- M_i * 3 / 2 + M_t * 5 / 4 - M_it * 6 / 5
    }
  } else {
    V <- M_i + M_t - M_it
  }
  return(V)
}

## Database function for test-sandwich.R
vcov_db <- function(k) {
  data(trade)
  if (k == 1) {
    base <- trade
    return(base)
  } else if (k == 2) {
    base <- trade
    base$Euros[base$Origin == "FR"] <- 0
    return(base)
  } else if (k == 3) {
    base <- trade
    base$Euros[base$Origin == "FR"] <- 0
    base$Euros_na <- base$Euros
    base$Euros_na[sample(nrow(base), 50)] <- NA
    base$Destination_na <- base$Destination
    base$Destination_na[sample(nrow(base), 50)] <- NA
    base$Origin_na <- base$Origin
    base$Origin_na[sample(nrow(base), 50)] <- NA
    base$Product_na <- base$Product
    base$Product_na[sample(nrow(base), 50)] <- NA
    return(base)
  }
}


## Cases function for test-nointercept.R

nointercept_cases <- function() {
  f_rhs <- c(
    "-1 + x1 + i(fe1)",
    "-1 + x1 + factor(fe1)",
    "-1 + x1 + i(fe1) + i(fe2)"
  )
  fmly <- c("gaussian", "poisson", "Gamma")
  DF <- data.frame(
    model = "ols", family = "NULL", formula = paste("y", f_rhs, sep = " ~ "),
    test_name = paste("lm", paste("formula", c(1, 2, 3)), sep = " - ")
  )
  for (k in 1:3) {
    f_lhs <- ifelse(fmly[k] == "gaussian", "y", "y_r")
    DF <- rbind(DF, data.frame(
      model = "glm",
      family = fmly[k],
      formula = paste(f_lhs, f_rhs, sep = " ~ "),
      test_name = paste("glm", fmly[k], paste("formula", c(1, 2, 3)), sep = " - ")
    ))
  }
  DF <- rbind(DF, data.frame(
    model = "negbin",
    family = "NULL",
    formula = paste("y", f_rhs, sep = " ~ "),
    test_name = paste("negbin", paste("formula", c(1, 2, 3)), sep = " - ")
  ))
  return(DF)
}

# fitting function selection for test-nointercept.R
# fixest_mod_select <- function(model, fmla, base, famly = NULL, weights = NULL, subset = NULL) {
#   fmla <- as.formula(fmla)
#   if (model == "ols") {
#     res <- feols(fml = fmla, data = base, weights = weights, subset = subset)
#     return(res)
#   } else if (model == "glm") {
#     res <- feglm(fml = fmla, data = base, family = famly, weights = weights, subset = subset)
#     return(res)
#   } else if (model == "negbin") {
#     res <- fenegbin(fml = fmla, data = base, subset = subset)
#     return(res)
#   } else if (model == "femlm") {
#     res <- femlm(fml = fmla, data = base, family = famly, subset = subset)
#     return(res)
#   } else if (model == "feNmlm"){
#     res <- feNmlm(fml = fmla, data = base, family = famly, subset = subset)
#   }
# }
#
# ## test-sandwich.R doesnt work with the subset argument (even if it is just ommited)
# ## Its a strange error! The following version of fixest_mod_slect omits subset argument
# fixest_mod_select2 <- function(model, fmla, base, famly = NULL, weights = NULL) {
#   fmla <- as.formula(fmla)
#   if (model == "ols") {
#     res <- feols(fml = fmla, data = base, weights = weights)
#     return(res)
#   } else if (model == "glm") {
#     res <- feglm(fml = fmla, data = base, family = famly, weights = weights)
#     return(res)
#   } else if (model == "negbin") {
#     res <- fenegbin(fml = fmla, data = base)
#     return(res)
#   } else if (model == "femlm") {
#     res <- femlm(fml = fmla, data = base, family = famly)
#     return(res)
#   } else if (model == "feNmlm"){
#     res <- feNmlm(fml = fmla, data = base, family = famly)
#   }
# }



stats_mod_select <- function(model, fmla, base, famly = NULL, weights = NULL) {
  fmla <- as.formula(fmla)
  if (model == "ols") {
    res <- lm(formula = fmla, data = base, weights = weights)
    return(res)
  } else if (model == "glm" | model == "femlm") {
    res <- glm(formula = fmla, data = base, family = famly, weights = weights)
    return(res)
  } else if (model == "negbin") {
    res <- MASS::glm.nb(formula = fmla, data = base)
    return(res)
  }
}

fixest_mod_select <- function(model, fmla, base, famly = NULL, ...) {
  # options(warn=-1)
  fmla <- as.formula(fmla)
  if (model == "ols") {
    res <- feols(fml = fmla, data = base, ...)
    # options(warn=0)
    return(res)
  } else if (model == "glm") {
    res <- feglm(fml = fmla, data = base, family = famly, ...)
    # options(warn=0)
    return(res)
  } else if (model == "negbin") {
    res <- fenegbin(fml = fmla, data = base, ...)
    # options(warn=0)
    return(res)
  } else if (model == "femlm") {
    res <- femlm(fml = fmla, data = base, family = famly, ...)
    # options(warn=0)
    return(res)
  } else if (model == "feNmlm") {
    res <- feNmlm(fml = fmla, data = base, family = famly, ...)
    # options(warn=0)
    return(res)
  }
}

### Cases for test-residuals.R


residuals_cases <- function() {
  method <- c(rep("ols", 2), rep("glm", 2), "femlm", "negbin")
  formula_fe <- "y_int ~ x1 | species"
  formula_stats <- "y_int ~ x1 + species"
  famly <- c(rep("NULL", 2), rep("poisson", 3), rep("NULL", 1))
  do_weight <- c(rep(c(FALSE, TRUE), 2), FALSE, FALSE)
  wghts <- c(rep(c("NULL", "base$w"), 2), rep("NULL", 2))
  tol <- c(rep(1e-6, 5), 1e-2)

  DF <- data.frame(
    method = method,
    formula_fe = formula_fe,
    formula_stats = formula_stats,
    family = famly,
    do_weight = do_weight,
    wghts = wghts,
    tol = tol,
    test_name = paste(method, paste("weight", do_weight, sep = "="), sep = " - ")
  )
  return(DF)
}

### test-fixef.R helper functinos

get_coef <- function(all_coef, x) {
  res <- all_coef[grepl(x, names(all_coef), perl = TRUE)]
  names(res) <- gsub(x, "", names(res), perl = TRUE)
  res
}

## fixef.strings

fixef.strings <- function() {
  AuxL1 <- list(
    c("species", "fe_bis"),
    c("species", "fe_bis", "fe_bis[[x3]]"),
    c("species", "fe_bis", "fe_bis[[x3]]", "species[[x2]]"),
    c("species", "fe_bis", "fe_bis[[x2]]", "fe_bis[[x3]]"),
    c("species", "fe_bis", "fe_bis[[x2]]", "fe_bis[[x3]]")
  )

  AuxL2 <- list(
    c("species", "factor\\(fe_bis\\)"),
    c("species", "factor\\(fe_bis\\)", "fe_bis::|:x3"),
    c("^species(?=[^:])", "^factor\\(fe_bis\\)", "fe_bis::|:x3", "species::|:x2"),
    c("^species", "^factor\\(fe_bis\\)", "fe_bis::(?=.+x2)|:x2", "fe_bis::(?=.+x3)|:x3"),
    c("^species", "^factor\\(fe_bis\\)", "fe_bis::(?=.+x2)|:x2", "fe_bis::(?=.+x3)|:x3")
  )

  return(list(AuxL1, AuxL2))
}
#### fixef cases for test-fixef.R
fixef_cases <- function() {
  Fmlas1 <- c(
    "y ~ x1 + x2 | species + fe_bis",
    "y ~ x1 + x2 | species + fe_bis[x3]",
    "y ~ x1 | species[x2] + fe_bis[x3] + fe_ter",
    "y ~ x1 | species + fe_bis[x2, x3] + fe_ter",
    "y ~ x1 | species + fe_bis[x2, x3] + fe_ter"
  )
  Fmlas2 <- c(
    "y ~ -1 + x1 + x2 + species + factor(fe_bis)",
    "y ~ -1 + x1 + x2 + species + factor(fe_bis) + i(fe_bis, x3)",
    "y ~ -1 + x1 + species + i(species, x2) + factor(fe_bis) + i(fe_bis, x3) + factor(fe_ter)",
    "y ~ x1 + species + factor(fe_bis) + i(fe_bis, x2) + i(fe_bis, x3) + factor(fe_ter)",
    "y ~ x1 + species + factor(fe_bis) + i(fe_bis, x2) + i(fe_bis, x3) + factor(fe_ter)"
  )
  K <- seq(1:5)
  DF <- data.frame(
    formula1 = Fmlas1,
    formula2 = Fmlas2,
    K,
    test_name = c(
      "With 2 x 1 FE",
      "With 1 FE + 1 FE 1 VS",
      "With 2 x (1 FE + 1 VS) + 1 FE",
      "With 2 x (1 FE) + 1 FE 2 VS",
      "Previous With weights"
    )
  )
  return(DF)
}



### Cases function for test-collinearity.R

collin_cases <- function() {
  useWeights <- c("NULL", "base$w") # c(FALSE, TRUE)
  model <- c("ols", "glm")
  use_fe <- c(FALSE, TRUE)
  fixest_formula1 <- c(" ~ x1 + constant", " ~ x1 + constant | species")
  stats_formula2 <- c(" ~ x1 + constant", " ~ x1 + constant + species")

  DF_l <- list()
  for (k in 1:2) {
    if (model[k] == "ols") {
      fmlas1 <- paste("y", fixest_formula1)
      fmlas2 <- paste("y", stats_formula2)
    } else if (model[k] == "glm") {
      fmlas1 <- paste("y_int", fixest_formula1)
      fmlas2 <- paste("y_int", stats_formula2)
    }
    df <- expand.grid(
      use_fe = use_fe,
      useWeights = useWeights,
      stringsAsFactors = FALSE
    )
    df_aux <- data.frame(use_fe, fixest_formula1 = fmlas1, stats_formula1 = fmlas2)
    df <- dplyr::left_join(df, df_aux, by = "use_fe")
    df$model <- model[k]
    df$test_name <- paste(model[k],
      paste("weights=", df$useWeights),
      paste("use_fe=", df$use_fe),
      sep = " - "
    )
    DF_l[[k]] <- df
  }
  DF <- rbind(DF_l[[1]], DF_l[[2]])
  return(DF)
}


## Cases function for test-hatvalues.R
hatvalues_cases <- function() {
  model <- c("ols", "glm")
  formulas <- c("y ~ x", "y_int ~ x")
  family <- c("NULL", "poisson")

  return(data.frame(
    model = model,
    formulas = formulas,
    family = family,
    test_name = paste(model, family, sep = " - ")
  ))
}

## Cases function for test-sandwhich.R
sandwcomp_cases <- function() {
  model <- c("ols", "ols", "glm", "glm")
  family <- c("NULL", "NULL", "poisson", "poisson")
  fmlas <- c(
    "y ~ x1 + I(x1**2) + factor(id)",
    "y ~ x1 + I(x1**2) | id",
    "y_int ~ x1 + I(x1**2) + factor(id)",
    "y_int ~ x1 + I(x1**2) | id"
  )
  FE <- rep(c(FALSE, TRUE), 2)
  DF <- data.frame(
    test_name = paste(model, paste("FE = ", FE)),
    model = model,
    family = family,
    fmlas = fmlas,
    FE = FE
  )
  return(DF)
}


subset_cases <- function() {
  method <- c("ols", "glm", "femlm", "negbin", "feNmlm")
  fmlas <- c(
    "y ~ x1 + x2",
    "y ~ x1 + x2 | species",
    "y ~ x1 + x2 | fe_bis",
    "y ~ x1 + x2 + i(fe_bis)",
    "y ~ x1 + x2 | fe_bis[x3]"
  )

  DF_l <- list()
  for (k in 1:5) {
    # for femlm, negbin and feNmlm omit lasa formula
    if (method[k] %in% c("femlm", "negbin", "feNmlm")) {
      Fmlas <- fmlas[-5]
    } else {
      Fmlas <- fmlas
    }
    DF_l[[k]] <- data.frame(
      test_name = paste(method[k], paste("formula", c(1:length(Fmlas))), sep = " - "),
      method = method[k],
      fmlas = Fmlas
    )
  }
  DF <- DF_l[[1]]
  for (k in 2:length(DF_l)) DF <- rbind(DF, DF_l[[k]])

  aux <- DF$method == "ols" | DF$method == "negbin"
  DF$famly[aux] <- "NULL"
  DF$famly[!aux] <- "poisson"
  return(DF)
}

# case function for test-lagging.R
lagging_cases <- function() {
  depvar <- c("y", "y_na", "y_0")
  p <- c("period", "period_txt", "period_date")
  method <- c("ols", "glm")

  DF <- expand.grid(p = p, depvar = depvar, stringsAsFactors = FALSE)
  aux <- DF$depvar != "y_0"
  DF$method[aux] <- "ols"
  DF$method[!aux] <- "glm"
  DF$fmly[aux] <- "NULL"
  DF$fmly[!aux] <- "poisson"

  DF$test_name <- paste(DF$method, paste("depvar =", DF$depvar), paste("p =", DF$p), sep = " - ")
  return(DF)
}


### cases function for test-predict.R
predict_cases <- function() {
  fmlas <- c("y ~ x1 | species + fe_bis", "y ~ x1 | species + fe_bis[x3 ]", "y ~ x1 + i(species)", "y ~ x1 + i(species) + i(fe_bis)")
  method <- c("ols", "glm", "femlm")
  fmly <- c("NULL", "poisson", "poisson")
  tol <- c(1e-12, 1e-3, 1e-12, 1e-12)
  DF_l <- list()
  for (k in 1:3) {
    Fmlas <- switch(method[k],
      "femlm" = fmlas[-2],
      fmlas
    )
    tol <- switch(method[k],
      "femlm" = tol[-2],
      tol
    )

    DF_l[[k]] <- data.frame(
      test_name = paste(method[k], fmly[k], paste("formula", seq(1, length(Fmlas))), sep = " - "),
      method = method[k],
      fmly = fmly[k],
      fmlas = Fmlas,
      tol = tol
    )
  }
  DF <- DF_l[[1]]
  for (k in 2:3) DF <- rbind(DF, DF_l[[k]])
  DF$test_name <- paste(seq(1, dim(DF)[1]), DF$test_name, sep = ") ")
  return(DF)
}


# case function for test-multiple.R
multiple_cases <- function(met = NULL) {
  method <- c("ols", "glm", "femlm", "feNmlm")
  s <- c("setosa", "versicolor", "virginica")
  fmly <- c("NULL", "poisson", "poisson", "poisson")
  df_aux <- expand.grid(rhs2 = c("x2", "x3"), rhs1 = "x1", lhs = c("y1", "y2"), stringsAsFactors = FALSE)
  # df_aux = df_aux[c(1,4,2,3),] # Ordered to match estimation order of fixtest
  fmlas <- paste(df_aux$lhs, paste(df_aux$rhs1, df_aux$rhs2, sep = " + "), sep = " ~ ")

  DF_l <- list()
  cont <- 0
  for (k in 1:4) {
    for (j in 1:3) {
      cont <- cont + 1
      test_name <- paste(method[k], paste("filter =", s[j]), paste("formula ", 1:length(fmlas)), sep = " - ")
      DF_l[[cont]] <- data.frame(
        test_name = test_name,
        method = method[k],
        fmly = fmly[k],
        s = s[j],
        fmlas = fmlas,
        num_fmla = 1:length(fmlas)
      )
    }
  }

  DF <- DF_l[[1]]
  for (k in 2:cont) DF <- rbind(DF, DF_l[[k]])
  DF <- DF[sort(DF$method, index.return = TRUE, decreasing = TRUE)$ix, ]
  DF$num_fmla <- rep(c(1:12), 4)
  DF$test_name <- paste(seq(1, dim(DF)[1]), DF$test_name, sep = ") ")

  if (!is.null(met)) {
    return(DF[DF$method == met, ])
  } else {
    return(DF)
  }
}
