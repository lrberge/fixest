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
  base$w = as.numeric(base$species)
  return(base)
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


## Cases function for test-nointercept.R

nointercept_cases = function(){
  f_rhs = c("-1 + x1 + i(fe1)",
            "-1 + x1 + factor(fe1)",
            "-1 + x1 + i(fe1) + i(fe2)")
  fmly = c("gaussian", "poisson", "Gamma")
  DF = data.frame(model = "ols", family = "NULL", formula = paste("y",f_rhs, sep = " ~ "),
                  test_name = paste("lm", paste("formula",c(1,2,3)), sep = " - ") )
  for(k in 1:3){
    f_lhs = ifelse(fmly[k] == "gaussian","y", "y_r")
    DF = rbind(DF,data.frame(model = "glm",
                             family = fmly[k],
                             formula = paste(f_lhs,f_rhs, sep = " ~ "),
                             test_name = paste("glm", fmly[k],paste("formula",c(1,2,3)), sep = " - ")))
  }
  DF = rbind(DF, data.frame(model = "negbin",
                            family = "NULL",
                            formula = paste("y",f_rhs, sep = " ~ "),
                            test_name = paste("negbin", paste("formula",c(1,2,3)), sep = " - ")))
  return(DF)
}

# fitting function selection for test-nointercept.R
fixest_mod_select = function(model, fmla, data, famly = NULL, weights = NULL){
  fmla = as.formula(fmla)
  if(model == "ols"){
    res = feols(fml = fmla, data = data, weights = weights)
    return(res)
  }else if(model == "glm"){
    res = feglm(fml = fmla,data = base, family = famly, weights = weights)
    return(res)
  }else if(model == "negbin"){
    res = fenegbin(fml = fmla, data = data, weights = weights)
    return(res)
  }else if(model == "mlm"){
    res = femlm(fml = fmla, data = data, family = famly)
    return(res)
  }
}


stats_mod_select = function(model, fmla, data, famly = NULL, weights = NULL){
  fmla = as.formula(fmla)
  if(model == "ols"){
    res = lm(formula = fmla, data = data, weights = weights)
    return(res)
  }else if(model == "glm" | model == "mlm"){
    res = glm(formula = fmla,data = base, family = famly, weights = weights)
    return(res)
  }else if(model == "negbin"){
    res = MASS::glm.nb(formula = fmla, data = data, weights = weights)
    return(res)
  }
}

residuals_cases = function(){
  method <- c("ols", "glm", "mlm", "negbin")
  fmla_fe = "y_int ~ x1 | species"
  fmla_stats = "y_int ~ x1 + species"
  do_weight <- c(FALSE, TRUE)
  wghts <- c("NULL", "base$w")

  DF = data.frame(method = "ols",
                  formula_fe = fmla_fe,
                  formula_stats = fmla_stats,
                  family = "NULL",
                  do_weight = do_weight,
                  wghts = wghts,
                  tol = 1e-6,
                  test_name = paste("ols", paste("weight", do_weight, sep = "="), sep =" - " ))
  for(k in 2:length(method)){
    tol = ifelse(method[k] == "negbin", 1e-2, 1e-6)
    fmly = ifelse((method[k] == "glm" | method[k] == "mlm"), "poisson", "NULL")
    DF = rbind(DF,data.frame(method = method[k],
                             formula_fe = fmla_fe,
                             formula_stats = fmla_stats,
                             family = fmly,
                             do_weight = do_weight,
                             wghts = wghts,
                             tol = tol,
                             test_name = paste(method[k], paste("weight", do_weight, sep = "="), sep =" - " )))
  }

  return(DF)
}

