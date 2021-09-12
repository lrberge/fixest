set.seed(0)
base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$y_int <- as.integer(base$y)
base$y_log <- sample(c(TRUE, FALSE), 150, TRUE)

# Data used for test-fitting.R ----

datab <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$y_int <- as.integer(base$y)
  base$y_log <- sample(c(TRUE, FALSE), 150, TRUE)
  base$fe_2 <- rep(1:5, 30)
  base$fe_3 <- sample(15, 150, TRUE)
  base$constant <- 5
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

# Data used for test-nointercept.R ----

datab2 <- function() {
  set.seed(0)
  base$y_r <- round(base$y)
  base$fe1 <- base$species
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

# Data used for test-vcov.R ----

datab3 <- function() {
  # We create "irregular" FEs
  set.seed(0)
  base <- data.frame(x = rnorm(20))
  base$y <- base$x + rnorm(20)
  base$y_int <- round(abs(base$y))
  base$fe1 <- rep(rep(1:3, c(4, 3, 3)), 2)
  base$fe2 <- rep(rep(1:5, each = 2), 2)
  return(base)
}

# Data used for test-sandwich.R ----
datab4 <- function() {
  set.seed(0)
  N <- 20
  G <- N / 5
  T <- N / G
  d <- data.frame(y = rnorm(N), x = rnorm(N), grp = rep(1:G, T), tm = rep(1:T, each = G))
  d$y_int <- round(abs(10 * d$y)) + 1
  return(d)
}

datab5 <- function() {
  base$y_int <- as.integer(base$y) + 1
  base$w <- as.numeric(base$species)
  return(base)
}

# Data used for test-fixef.R ----

datab6 <- function() {
  set.seed(0)
  base$x4 <- rnorm(150) + 0.25 * base$y
  base$fe_bis <- sample(10, 150, TRUE)
  base$fe_ter <- sample(15, 150, TRUE)
  return(base)
}

# database for test-collinearity.R ----
datab7 <- function() {
  base$constant <- 5
  base$y_int <- as.integer(base$y)
  base$w <- as.vector(unclass(base$species) - 0.95)
  return(base)
}

datab8 <- function() {
  set.seed(0)
  base$z <- sample(5, 150, TRUE)
  base$species_na <- base$species
  base$species_na[base$species == "setosa"] <- NA
  return(base)
}

# Data used for test-demean.R ----

datab9 <- function() {
  data(trade)
  base <- trade
  base$ln_euros <- log(base$Euros)
  base$ln_dist <- log(base$dist_km)
  return(base)
}

# Data used for test-hatvalues.R ----

datab10 <- function() {
  set.seed(0)
  x <- sin(1:10)
  y <- rnorm(10)
  y_01 <- c(y > 0) + 0
  y_int <- rpois(10, 2) + 2
  return(data.frame(y = y, y_int = y_int, y_01 = y_01, x))
}

# Data used for test-subset.R ----

datab11 <- function() {
  set.seed(5)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$y_int <- as.integer(base$y)
  base$y_log <- sample(c(TRUE, FALSE), 150, TRUE)
  base$fe_bis <- sample(letters, 150, TRUE)
  base$x4 <- rnorm(150)
  base$x1[sample(150, 5)] <- NA
  return(base)
}

# Data used for test-lagging.R ----

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
  base$y_00 <- base$y_0 + 1

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

# Data used for test-predict.R ----

datab13 <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$y_int <- as.integer(base$y)
  base$y_log <- sample(c(TRUE, FALSE), 150, TRUE)
  base$fe_bis <- sample(letters, 150, TRUE)
  return(base)
}

datab14 <- function() {
  set.seed(2)
  base$y1 <- base$y
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

# Data used for test-model_matrix.R ----

datab16 <- function() {
  set.seed(0)
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")
  base$y_int <- as.integer(base$y)
  base$y_log <- sample(c(TRUE, FALSE), 150, TRUE)
  base$y1 <- base$y
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

  # for glm models:
  base$y_int <- round(base$y1)
  return(list(base, base_bis))
}

# Data used for test-nonlinear.R ----

datab17 <- function() {
  base$tab <- c("versicolor" = 5, "setosa" = 0, "virginica" = -5)
  base$tab2 <- base$tab[base$species]
  base$var_spec <- as.numeric(base$tab[base$species])
  return(base)
}

datab18 <- function() {
  base$fe1 <- base$species
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

ev_par <- function(string) {
  eval(parse(text = string))
}

fomla <- function(model, fml_fixest1, fml_stats1) {
  if (model %in% c("binomial", "logit", "quasibinomial")) {
    fml_fixest <- paste("y_01", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y_01", fml_stats1, sep = " ~ ")
  } else if (model %in% c("poisson", "quasipoisson")) {
    fml_fixest <- paste("y_int_null", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y_int_null", fml_stats1, sep = " ~ ")
  } else if (model %in% c("negbin", "Gamma", "gaussian")) {
    fml_fixest <- paste("y_int", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y_int", fml_stats1, sep = " ~ ")
  } else if (model %in% c("inverse.gaussian", "quasi", "gaussian")) {
    fml_fixest <- paste("y", fml_fixest1, sep = " ~ ")
    fml_stats <- paste("y", fml_stats1, sep = " ~ ")
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

  df_aux <- expand.grid(
    ID = 1:length(fml_stats),
    my_weight = c("NULL", "base$w"),
    my_offset = c("NULL", "base$offset_value"),
    stringsAsFactors = FALSE
  )
  df_fmla <- data.frame(
    fml_fixest = fml_fixest,
    fml_stats = fml_stats,
    ID = 1:length(fml_stats)
  )
  DF <- merge(df_aux, df_fmla, by = "ID", all.x = TRUE)[-1] # deleting ID

  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", paste(paste("weight:", DF$my_weight),
    paste("offset:", DF$my_offset),
    paste("formula =", 1:length(fml_stats)),
    sep = " - "
  ))
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
  for (k in 1:length(Models)) {
    df_aux <- expand.grid(
      my_offset = c("NULL", "base$offset_value"),
      my_weight = c("NULL", "base$w"),
      ID = 1:length(fml_stats1),
      stringsAsFactors = FALSE
    )
    fmlas <- fomla(Models[k], fml_fixest1, fml_stats1)
    df_fmla <- data.frame(
      ID = 1:length(fml_stats1),
      fml_fixest = fmlas[[1]],
      fml_stats = fmlas[[2]]
    )
    DF_l[[k]] <- merge(df_aux, df_fmla, by = "ID", all.x = TRUE)[-1] # deleting ID
    DF_l[[k]]$my_family <- Models[k]
  }

  DF <- do.call("rbind", DF_l)
  DF <- DF[!(DF$my_family == "Gamma" & DF$my_offset != "NULL"), ]

  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", paste(DF$my_family,
    paste("offset:", DF$my_offset),
    paste("weights:", DF$my_weight),
    paste("formula =", 1:length(fml_stats1)),
    sep = " - "
  ))
  return(DF)
}

femlm_cases <- function(fmly = c("poisson", "negbin", "logit", "gaussian")) {
  fml_fixest1 <- c(
    "x1",
    "x1 | species",
    "x1 | species + fe_2",
    "x1 + x2 | species^fe_2"
  )
  fml_stats1 <- c(
    "x1",
    "x1 + factor(species)",
    "x1 + factor(species) + factor(fe_2)",
    "x1 + x2 + paste(species, fe_2)"
  )

  DF_l <- list()
  for (k in 1:length(fmly)) {
    aux <- fomla(fmly[k], fml_fixest1, fml_stats1)
    fml_fixest <- aux[[1]]
    fml_stats <- aux[[2]]
    df_fmla <- data.frame(
      fml_fixest = fml_fixest,
      fml_stats = fml_stats,
      ID = 1:length(fml_fixest)
    )
    df_aux <- expand.grid(
      my_offset = c("NULL", "base$offset_value"),
      my_family = fmly[k],
      ID = 1:length(fml_fixest),
      stringsAsFactors = FALSE
    )
    DF_l[[k]] <- merge(df_aux, df_fmla, by = "ID", all.x = TRUE)[-1] # deleting ID
  }

  DF <- do.call("rbind", DF_l)
  DF <- DF[sort(DF$my_family, index.return = TRUE)$ix, ]

  DF$my_family_stats <- DF$my_family
  DF$my_family_stats[DF$my_family == "logit"] <- "binomial"
  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", paste(DF$my_family,
    paste("offset:", DF$my_offset),
    paste("formula =", 1:length(fml_stats1)),
    sep = " - "
  ))
  return(DF)
}

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

# Case function for test-fitmethod.R ----

fitmethod_cases <- function() {
  y_dep <- c("base$y", "base$y_int", "base$y_log")
  y_dep <- c("y", "y_int", "y_log")
  with_fmly <- c(FALSE, TRUE)
  method <- c("ols", "glm")
  DF <- expand.grid(
    y_dep = y_dep,
    with_fmly = with_fmly,
    method = method,
    stringsAsFactors = FALSE
  )
  DF$fmly <- c(rep("NULL", 6), rep("poisson", 6))
  DF$test_name <- paste(DF$method, paste("with family =", DF$with_fmly), paste("y_dep =", DF$y_dep), sep = " - ")
  DF$test_name <- paste(seq(1, dim(DF)[1]), DF$test_name, sep = ") ")
  return(DF)
}


## Case function for test-vcov
vcov_cases1 <- function() {
  method <- c("ols", "glm", "femlm", "feNmlm")
  DF <- expand.grid(
    adj = c(FALSE, TRUE),
    k_val = c("none", "nested", "full"),
    method = method,
    stringsAsFactors = FALSE
  )
  DF$K <- ifelse(DF$k_val == "none", 1, 8)
  DF$fmly <- ifelse(DF$method == "ols", "NULL", "poisson")
  DF$my_adj <- ifelse(DF$adj, (20 - 1) / (20 - DF$K), 1)

  DF$test_name <- paste(DF$method, DF$fmly, paste("adj = ", DF$adj), paste("K =", DF$K), sep = " - ")
  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", DF$test_name)
  return(DF)
}

vcov_cases2 <- function() {
  DF <- expand.grid(
    k_val = c("none", "nested", "full"),
    tdf = c("conventional", "min"),
    c_adj = c(FALSE, TRUE),
    adj = c(FALSE, TRUE),
    method = c("ols", "glm", "femlm", "feNmlm"),
    stringsAsFactors = FALSE
  )

  DF$K[DF$k_val == "none"] <- 1
  DF$K[DF$k_val == "nested"] <- 6
  DF$K[DF$k_val == "full"] <- 8
  DF$fmly <- ifelse(DF$method == "ols", "NULL", "poisson")
  DF$my_adj <- ifelse(DF$adj, (20 - 1) / (20 - DF$K), 1)
  DF$cluster_factor <- ifelse(DF$c_adj, 3 / 2, 1)
  DF$df <- ifelse(DF$tdf == "min", 2, 20 - 8)

  DF$test_name <- paste(DF$k_val, DF$tdf, DF$c_adj, DF$adj, sep = " - ")
  DF$test_name <- paste(DF$method, DF$fmly, paste("K =", DF$K), paste("c_adj =", DF$adj), paste("tdf =", DF$tdf), sep = " - ")
  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", DF$test_name)
  return(DF)
}

vcov_cases3 <- function() {
  DF <- expand.grid(
    cdf = c("conventional", "min"),
    tdf = c("conventional", "min"),
    k_val = c("none", "nested", "full"),
    c_adj = c(FALSE, TRUE),
    adj = c(FALSE, TRUE),
    method = c("ols", "glm", "femlm", "feNmlm"),
    stringsAsFactors = FALSE
  )

  DF$K[DF$k_val == "none"] <- 1
  DF$K[DF$k_val == "nested"] <- 2
  DF$K[DF$k_val == "full"] <- 8
  DF$fmly <- ifelse(DF$method == "ols", "NULL", "poisson")
  DF$my_adj <- ifelse(DF$adj, (20 - 1) / (20 - DF$K), 1)
  DF$cluster_factor <- ifelse(DF$c_adj, 3 / 2, 1)
  DF$df <- ifelse(DF$tdf == "min", 2, 20 - 8)

  DF$test_name <- paste(DF$method, DF$fmly,
    paste("cdf =", DF$cdf),
    paste("tdf =", DF$tdf),
    paste("c_adj =", DF$c_adj),
    paste("adj =", DF$adj),
    sep = " - "
  )
  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", DF$test_name)
  return(DF)
}

vcov_cases4 <- function() {
  method <- c("ols", rep("glm", 4))
  fmly <- c("NULL", "poisson", "quasipoisson", "Gamma", "binomial")
  y_dep <- c("y", "y_int", "y_int", "y", "y_01")

  DF <- data.frame(
    method = method,
    fmly = fmly,
    y_dep = y_dep
  )
  DF$test_name <- paste(method, fmly, y_dep, sep = " - ")
  DF$test_name <- paste0(1:dim(DF)[1], "1)", DF$test_name)

  return(DF)
}

# Function for test-vcov.R ----

V_matrix <- function(M_i, M_t, M_it, c_adj, cdf) {
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

# Data function for test-sandwich.R ----

vcov_db <- function(k) {
  data(trade)
  if (k == 1) {
    BASE <- trade
    return(BASE)
  } else if (k == 2) {
    BASE <- trade
    BASE$Euros[BASE$Origin == "FR"] <- 0
    return(BASE)
  } else if (k == 3) {
    BASE <- trade
    BASE$Euros[BASE$Origin == "FR"] <- 0
    BASE$Euros_na <- BASE$Euros
    BASE$Euros_na[sample(nrow(BASE), 50)] <- NA
    BASE$Destination_na <- BASE$Destination
    BASE$Destination_na[sample(nrow(BASE), 50)] <- NA
    BASE$Origin_na <- BASE$Origin
    BASE$Origin_na[sample(nrow(BASE), 50)] <- NA
    BASE$Product_na <- BASE$Product
    BASE$Product_na[sample(nrow(BASE), 50)] <- NA
    return(BASE)
  }
}

# Case function for test-sandwich.R ----

SE_cases <- function() {
  model <- c("ols", rep("glm", 4), rep("femlm", 2))
  y_dep <- c(
    "y",
    rep("y_int", 3), "y",
    "y_int", "y"
  )

  fmly <- c(
    "NULL",
    "Gamma", "poisson", "quasipoisson", "gaussian",
    "poisson", "gaussian"
  )

  DF <- data.frame(
    y_dep = y_dep,
    model = model,
    fmly = fmly
  )

  DF$test_name <- paste(model, fmly, y_dep, sep = " - ")
  return(DF)
}

# Case function for test-nointercept.R ----

nointercept_cases <- function() {
  f_rhs <- c(
    "-1 + x1 + i(fe1)",
    "-1 + x1 + factor(fe1)",
    "-1 + x1 + i(fe1) + i(fe2)"
  )

  DF_l <- list()
  # OLS
  DF_l[[1]] <- data.frame(
    fmly = "NULL",
    f_rhs = f_rhs,
    method = "ols",
    y_dep = "y"
  )

  # GLM - FEMLM
  fmly <- c("poisson", "gaussian", "Gamma")
  df_aux <- expand.grid(
    f_rhs = f_rhs,
    fmly = fmly,
    method = c("glm", "femlm"),
    stringsAsFactors = FALSE
  )
  df_ydep <- data.frame(
    fmly = fmly,
    y_dep = c("y", "y_r", "y_r")
  )
  DF <- merge(df_aux, df_ydep, by = "fmly", all.x = TRUE)
  DF <- DF[sort(DF$method, index.return = TRUE)$ix, ]
  DF_l[[2]] <- DF[!(DF$fmly == "Gamma" & DF$method == "femlm"), ]

  ## negbin
  DF_l[[3]] <- data.frame(
    fmly = "NULL",
    f_rhs = f_rhs,
    method = "negbin",
    y_dep = "y"
  )

  DF <- do.call("rbind", DF_l)
  DF$test_name <- paste(DF$method, DF$fmly, paste("formula = ", DF$y_dep, DF$f_rhs), sep = " - ")
  return(DF)
}

stats_mod_select <- function(model, fmla, base, famly = NULL, ...) {
  fmla <- as.formula(fmla)
  if (model == "ols") {
    res <- lm(formula = fmla, data = base, ...)
    return(res)
  } else if (model == "glm" | model == "femlm") {
    res <- glm(formula = fmla, data = base, family = famly, ...)
    return(res)
  } else if (model == "negbin") {
    res <- MASS::glm.nb(formula = fmla, data = base, ...)
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

# Cases for test-residuals.R ----

residuals_cases <- function() {
  method <- c("ols", "glm", "femlm", "negbin")
  formula_fe <- "y_int ~ x1 | species"
  formula_stats <- "y_int ~ x1 + species"
  fmly <- c("poisson", "quasipoisson", "Gamma", "gaussian")
  fmly_mlm <- c("poisson", "gaussian") # not yet included logit
  wghts <- c("NULL", "base$w")

  DF_l <- list()
  for (k in 1:4) {
    Famly <- fmly
    Wghts <- wghts
    if (method[k] %in% c("negbin", "ols")) Famly <- "NULL"
    if (method[k] == "negbin") Wghts <- "NULL"
    if (method[k] == "femlm") Famly <- fmly_mlm

    DF_l[[k]] <- expand.grid(
      method = method[k],
      formula_fe = formula_fe,
      formula_stats = formula_stats,
      weights = Wghts,
      family = Famly,
      stringsAsFactors = FALSE
    )
  }

  DF <- do.call("rbind", DF_l)
  DF$test_name <- paste(DF$method,
    paste("weight:", DF$weights),
    sep = " - "
  )
  DF$test_name <- paste0(seq(1, dim(DF)[1]), ") ", DF$test_name)
  return(DF)
}

# test-fixef.R helper functions ----

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

# fixef cases for test-fixef.R ----

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

# Case function for test-collinearity.R ----

collin_cases2 <- function() {
  useWeights <- c("NULL", "base$w") # c(FALSE, TRUE)
  model <- c("ols", "glm", "femlm")
  use_fe <- c(FALSE, TRUE)
  fmly <- c("poisson", "Gamma") # , "quasipoisson", "Gamma", "binomial")
  DF_l <- list()
  for (k in 1:length(model)) {
    Fmly <- fmly
    UseWeights <- useWeights
    if (model[k] == "ols") Fmly <- "NULL"
    if (model[k] == "femlm") {
      UseWeights <- "NULL"
      Fmly <- c("poisson", "negbin")
    }
    df_aux <- expand.grid(
      useWeights = UseWeights,
      ID = c(1, 2),
      model = model[k],
      fmly = Fmly,
      stringsAsFactors = FALSE
    )

    df_fmla <- data.frame(
      ID = c(1, 2),
      use_fe = use_fe,
      fixest_formula1 = c("x1 + constant", "x1 + constant | species"),
      stats_formula2 = c("x1 + constant", "x1 + constant + species")
    )

    DF_l[[k]] <- merge(df_aux, df_fmla, by = "ID", all.x = TRUE)[-1] # deleting ID
  }
  DF <- do.call("rbind", DF_l)

  DF$y_dep <- "y_int"
  DF$y_dep[DF$model == "ols" | DF$fmly == "gaussian" | DF$fmly == "Gamma"] <- "y"
  DF$test_name <- paste(DF$model, fmly, paste("FE:", DF$use_fe), paste("weights:", DF$useWeights), paste("y_dep:", DF$y_dep), sep = " - ")
  DF$test_name <- paste0(1:dim(DF)[1], ") ", DF$test_name)
  return(DF)
}

# Case function for test-hatvalues.R ----

hatvalues_cases <- function() {
  model <- c("ols", rep("glm", 6))
  formulas <- c(
    "y ~ x",
    "y_int ~ x",
    "y_int ~ x",
    "y ~ x",
    "y_int ~ x",
    "y_01 ~ x",
    "y_01 ~ x"
  )
  family <- c(
    "NULL",
    "poisson",
    "quasipoisson",
    "gaussian",
    "Gamma",
    "binomial",
    "quasibinomial"
  )

  return(data.frame(
    model = model,
    formulas = formulas,
    family = family,
    test_name = paste(model, family, sep = " - ")
  ))
}

# Case function for test-sandwich.R ----

sandwcomp_cases <- function() {
  model <- c("ols", "ols", "glm", "glm", "glm", "glm")
  family <- c("NULL", "NULL", "poisson", "poisson", "Gamma", "Gamma")
  fmlas <- c(
    "y ~ x1 + I(x1**2) + factor(id)",
    "y ~ x1 + I(x1**2) | id",
    "y_int ~ x1 + I(x1**2) + factor(id)",
    "y_int ~ x1 + I(x1**2) | id",
    "y_int ~ x1 + I(x1**2) + factor(id)",
    "y_int ~ x1 + I(x1**2) | id"
  )
  FE <- rep(c(FALSE, TRUE), 3)

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
  # DF <- DF_l[[1]]
  # for (k in 2:length(DF_l)) DF <- rbind(DF, DF_l[[k]])
  DF <- do.call("rbind", DF_l)

  aux <- DF$method == "ols" | DF$method == "negbin"
  DF$famly[aux] <- "NULL"
  DF$famly[!aux] <- "poisson"
  return(DF)
}

# Case function for test-lagging.R ----

lagging_cases <- function() {
  p <- c("period", "period_txt", "period_date")
  DF_l <- list()
  # OLS
  DF_l[[1]] <- expand.grid(
    p = p,
    depvar = c("y", "y_na"),
    fmly = "NULL",
    method = "ols",
    stringsAsFactors = FALSE
  )

  # GLM
  fmly <- c("poisson", "quasipoisson", "Gamma", "gaussian")
  y_dep_glm <- c("y_0", "y_0", "y_00", "y")
  df_aux <- expand.grid(
    p = p,
    fmly = fmly,
    method = "glm",
    stringsAsFactors = FALSE
  )
  df_dep <- data.frame(
    fmly = fmly,
    depvar = y_dep_glm
  )
  DF_l[[2]] <- merge(df_aux, df_dep, by = "fmly", all.x = TRUE)

  # FEMLM
  fmly <- c("poisson", "gaussian")
  y_dep <- c("y_0", "y")
  df_aux <- expand.grid(
    p = p,
    fmly = fmly,
    method = "femlm",
    stringsAsFactors = FALSE
  )
  df_dep <- data.frame(
    fmly = fmly,
    depvar = y_dep
  )
  DF_l[[3]] <- merge(df_aux, df_dep, by = "fmly", all.x = TRUE)

  DF <- do.call("rbind", DF_l)

  DF$test_name <- paste(DF$method, paste("depvar =", DF$depvar), paste("p =", DF$p), sep = " - ")
  return(DF)
}

# Case function for test-predict.R ----

predict_cases <- function() {
  y_dep <- c("y")
  fmlas <- c("x1 | species + fe_bis", "x1 | species + fe_bis[x3]", "x1 + i(species)", "x1 + i(species) + i(fe_bis)")
  method <- c("ols", "glm", "femlm")
  fmly_glm <- c("poisson", "quasipoisson", "Gamma", "inverse.gaussian")
  fmly_mlm <- c("poisson", "gaussian", "negbin")
  DF_l <- list()
  for (k in 1:length(method)) {
    Fmlas <- switch(method[k],
      "femlm" = fmlas[-2], # without varying slope
      fmlas
    )
    Fmly <- switch(method[k],
      "ols" = "NULL",
      "glm" = fmly_glm,
      "femlm" = fmly_mlm
    )

    DF_l[[k]] <- expand.grid(
      fmlas = Fmlas,
      y_dep = y_dep,
      fmly = Fmly,
      method = method[k],
      stringsAsFactors = FALSE
    )
  }

  DF <- do.call("rbind", DF_l)

  DF$tol <- 1e-12
  DF$tol[DF$fmlas == "x1 | species + fe_bis[x3]"] <- 1e-3
  DF$test_name <- paste(DF$method, DF$fmly, paste("formula", DF$fmlas), sep = " - ")
  DF$test_name <- paste(seq(1, dim(DF)[1]), DF$test_name, sep = ") ")
  return(DF)
}

# Case function for test-multiple.R ----

multiple_cases <- function(met = c("ols", "glm", "femlm", "feNmlm"), famly = "poisson") {
  method <- met
  fmly <- famly
  rhs <- c("x2", "x3")
  lhs <- c("y1", "y2")
  s <- c("setosa", "versicolor", "virginica")
  df_aux <- expand.grid(
    rhs = rhs,
    lhs = lhs,
    s = s,
    stringsAsFactors = FALSE
  )

  DF_l <- list()
  for (k in 1:length(method)) {
    Fmly <- switch(method[k],
      "ols" = "NULL",
      fmly
    )

    DF_l[[k]] <- data.frame(
      method = method[k],
      fmly = Fmly,
      rhs = df_aux$rhs,
      lhs = df_aux$lhs,
      s = df_aux$s,
      num_fmla = 1:dim(df_aux)[1]
    )
  }

  DF <- do.call("rbind", DF_l)
  DF$test_name <- paste(DF$method, paste("filter =", DF$s), paste("formula ", DF$num_fmla), sep = " - ")
  DF$test_name <- paste(seq(1, dim(DF)[1]), DF$test_name, sep = ") ")
  return(DF)
}



multiple_cases2 <- function(met = c("ols", "glm", "femlm", "feNmlm"), famly = "poisson") {
  method <- met
  fmly <- famly
  all_rhs <- c("", "x2", "x3")
  s <- c("all", "setosa", "versicolor", "virginica")
  lhs <- c("y1", "y2")
  n_rhs <- 1:3

  df_aux <- expand.grid(
    n_rhs = n_rhs,
    lhs = lhs,
    s = s,
    stringsAsFactors = FALSE
  )

  DF_l <- list()
  for (k in 1:length(method)) {
    Fmly <- switch(method[k],
      "ols" = "NULL",
      fmly
    )

    DF_l[[k]] <- data.frame(
      method = method[k],
      fmly = Fmly,
      s = df_aux$s,
      lhs = df_aux$lhs,
      n_rhs = df_aux$n_rhs,
      num_fmla = 1:dim(df_aux)[1]
    )
  }
  DF <- do.call("rbind", DF_l)
  DF$test_name <- paste(DF$method, paste("filter =", DF$s), paste("formula ", 1:dim(df_aux)[1]), sep = " - ")
  DF$test_name <- paste(seq(1, dim(DF)[1]), DF$test_name, sep = ") ")
  return(DF)
}

# Case function fo test-model_matrix.R ----

model_matrix_cases <- function() {
  method <- c("ols", "glm", "femlm")
  fmly <- c("NULL", "poisson", "poisson")
  y_dep <- c("y1", "y_int", "y_int")
  DF <- data.frame(
    method = method,
    fmly = fmly,
    y_dep = y_dep
  )
  DF$test_name <- paste(method, fmly, y_dep, sep = " - ")
  return(DF)
}
