#----------------------------------------------#
# Author: Dianyi Yang (+ AI)
# Date creation: Wed Nov 26 23:11:00 2025
# ~: Anderson-Rubin test for IV models
#----------------------------------------------#



#' Anderson-Rubin test for IV models
#'
#' @description
#' Performs the Anderson-Rubin test for instrumental variable (IV) models. 
#' This test is robust to weak instruments and can be used to test hypotheses 
#' about the coefficients of endogenous variables.
#'
#' @param object A `fixest` object obtained from [`feols`] with an IV specification 
#'   (i.e., formula includes `| endo ~ inst`).
#' @param beta0 Numeric vector. The hypothesized values for the endogenous variable 
#'   coefficients under the null hypothesis. If `NULL` (default), tests 
#'   H0: beta = 0 for all endogenous variables. The length must match the number 
#'   of endogenous variables in the model.
#' @param vcov Versatile argument to specify the VCOV. If `NULL` (default), uses 
#'   the same VCOV specification as the original model. See [`vcov.fixest`] for 
#'   details on available options.
#' @param ... Additional arguments passed to [`summary.fixest`].
#'
#' @details
#' The Anderson-Rubin (AR) test is a weak-instrument robust test for the null 
#' hypothesis H0: beta = beta0, where beta is the vector of coefficients on the 
#' endogenous variables.
#'
#' The test proceeds as follows:
#' \enumerate{
#'   \item Construct the transformed outcome: y* = y - X * beta0, where X contains 
#'     the endogenous variables.
#'   \item Regress y* on the exogenous controls (W), instruments (Z), and any 
#'     fixed effects, using the same VCOV specification as the original model.
#'   \item Test the joint significance of the instrument coefficients (H0: pi = 0).
#' }
#'
#' The AR test statistic follows an F-distribution under the null, regardless of 
#' instrument strength. This makes it particularly useful when instruments may be 
#' weak.
#'
#' @return
#' An object of class `fixest_ar` containing:
#' \item{beta0}{The hypothesized coefficient values tested.}
#' \item{stat}{The AR test statistic (F-statistic).}
#' \item{df1}{Numerator degrees of freedom (number of instruments).}
#' \item{df2}{Denominator degrees of freedom.}
#' \item{pvalue}{The p-value of the test.}
#' \item{method}{Character string: "Anderson-Rubin test".}
#' \item{vcov_type}{The type of VCOV used for the test.}
#' \item{endo_names}{Names of the endogenous variables.}
#' \item{inst_names}{Names of the instruments.}
#' \item{model}{The original `feols` IV object.}
#'
#' @seealso 
#' [`ar_confint`] for computing AR confidence intervals, [`wald`] for standard 
#' Wald tests, [`feols`] for IV estimation.
#'
#' @references
#' Anderson, T. W. and Rubin, H. (1949). "Estimation of the Parameters of a Single 
#' Equation in a Complete System of Stochastic Equations." Annals of Mathematical 
#' Statistics, 20(1), 46-63.
#'
#' @examples
#' # Simple IV example
#' set.seed(123)
#' n <- 1000
#' z <- rnorm(n)
#' x <- 0.5 * z + rnorm(n)
#' y <- 1 + 0.8 * x + rnorm(n)
#' data <- data.frame(y = y, x = x, z = z)
#'
#' # IV estimation
#' est <- feols(y ~ 1 | x ~ z, data = data)
#'
#' # Test H0: beta = 0
#' ar_test(est)
#'
#' # Test H0: beta = 0.8 (close to true value, should not reject)
#' ar_test(est, beta0 = 0.8)
#'
#' # Test H0: beta = 0 (far from true value, should reject)
#' ar_test(est, beta0 = 0)
#'
#' # With clustered standard errors
#' data$cl <- rep(1:100, each = 10)
#' est_cl <- feols(y ~ 1 | x ~ z, data = data, vcov = ~cl)
#' ar_test(est_cl)
#'
#' @export
ar_test = function(object, beta0 = NULL, vcov = NULL, ...) {
  
  # Input validation
  check_arg(object, "class(fixest)")
  
  # Check that object is an IV model
  if (!isTRUE(object$iv)) {
    stop("The `ar_test` function requires an IV model. ",
         "The input object is not an IV estimation (no endogenous variables/instruments found). ",
         "Please provide a `feols` object estimated with an IV formula ",
         "(e.g., `feols(y ~ controls | endo ~ inst, data)`).")
  }
  
  # Check that it's a second stage estimation
  if (isTRUE(object$iv_stage == 1)) {
    stop("The `ar_test` function requires a second-stage IV estimation. ",
         "The input object appears to be a first-stage estimation. ",
         "Please use the main IV estimation result, not `summary(est, stage = 1)`.")
  }
  
  # Get endogenous variable names and instrument names
  endo_names <- object$iv_endo_names
  inst_names <- object$iv_inst_names_xpd  # expanded instrument names
  n_endo <- length(endo_names)
  n_inst <- length(inst_names)
  
  # Handle beta0
  if (is.null(beta0)) {
    beta0 <- rep(0, n_endo)
  } else {
    check_arg(beta0, "numeric vector no na")
    if (length(beta0) != n_endo) {
      stop("Argument `beta0` has length ", length(beta0), " but there are ", 
           n_endo, " endogenous variables (", paste(endo_names, collapse = ", "), "). ",
           "Please provide a vector of length ", n_endo, ".")
    }
  }
  
  # Name beta0 for clarity
  names(beta0) <- endo_names
  
  # Fetch the data
  data <- fetch_data(object, "To apply `ar_test`, ")
  
  # Get the dependent variable
  fml_linear <- formula(object, type = "linear")
  y_vec <- eval(fml_linear[[2]], data)
  
  # Get the endogenous variable matrix
  endo_mat <- model.matrix(object, data = data, type = "iv.endo", collin.rm = FALSE)
  
  # Compute transformed outcome: y* = y - X * beta0
  y_tilde <- as.numeric(y_vec - endo_mat %*% beta0)
  
  # Store y_tilde in a temporary column with unique name
  temp_y_name <- "..ar_y_tilde"
  # Ensure uniqueness if column already exists (with a maximum of 100 attempts)
  attempt <- 0
  while (temp_y_name %in% names(data) && attempt < 100) {
    temp_y_name <- paste0("..ar_y_tilde_", sample.int(1e6, 1))
    attempt <- attempt + 1
  }
  if (temp_y_name %in% names(data)) {
    stop("Could not create a unique temporary column name in the data. ",
         "This is unexpected - please report this issue.")
  }
  data[[temp_y_name]] <- y_tilde
  
  # Build the AR regression formula
  # We need: y_tilde ~ controls + instruments | fixed_effects
  
  # Get the exogenous controls part (the linear formula RHS)
  fml_linear_rhs <- object$fml_all$linear
  if (length(fml_linear_rhs) == 3) {
    controls_str <- deparse_long(fml_linear_rhs[[3]])
  } else {
    controls_str <- "1"
  }
  
  # Get the instrument formula RHS
  fml_iv <- object$fml_all$iv
  inst_str <- deparse_long(fml_iv[[3]])
  
  # Get the fixed effects part if any
  fml_fe <- object$fml_all$fixef
  if (!is.null(fml_fe)) {
    fe_str <- paste0(" | ", deparse_long(fml_fe[[2]]))
  } else {
    fe_str <- ""
  }
  
  # Build the new formula: y_tilde ~ controls + instruments | FE
  new_fml_str <- paste0(temp_y_name, " ~ ", controls_str, " + ", inst_str, fe_str)
  new_fml <- as.formula(new_fml_str)
  
  # Get weights if present
  weights_val <- NULL
  if (!is.null(object$weights)) {
    weights_val <- object$weights
  }
  
  # Get vcov specification - use original if not specified
  if (is.null(vcov)) {
    # Try to get the vcov from the summarized object
    if (!is.null(object$summary_flags$vcov)) {
      vcov <- object$summary_flags$vcov
      if (is.null(vcov)) vcov <- "iid"
    } else {
      vcov <- "iid"
    }
  }
  
  # Run the AR regression
  ar_fit <- feols(new_fml, data = data, weights = weights_val, vcov = vcov, ...)
  
  # Summarize to ensure we have the VCOV
  ar_fit_sum <- summary(ar_fit, vcov = vcov, ...)
  
  # Get the coefficient names in the AR regression that correspond to instruments
  ar_coef_names <- names(ar_fit_sum$coefficients)
  
  # Find which coefficients correspond to instruments
  # The instrument names in the AR fit should match the expanded instrument names
  inst_in_fit <- intersect(ar_coef_names, inst_names)
  
  if (length(inst_in_fit) == 0) {
    # Try a more flexible approach using the instrument formula variables
    # This handles cases where instruments are expanded (e.g., factor variables)
    inst_formula_vars <- all.vars(fml_iv[[3]])
    for (v in inst_formula_vars) {
      # Escape special regex characters in variable name
      v_escaped <- gsub("([\\^$*+?{}\\[\\]()|.])", "\\\\\\1", v)
      pattern <- paste0("^", v_escaped, "($|::)")
      matches <- grep(pattern, ar_coef_names, value = TRUE)
      inst_in_fit <- c(inst_in_fit, matches)
    }
    inst_in_fit <- unique(inst_in_fit)
  }
  
  if (length(inst_in_fit) == 0) {
    stop("Could not identify instrument coefficients in the AR regression. ",
         "This may be a bug - please report it.")
  }
  
  # Perform the Wald test on the instrument coefficients
  wald_result <- wald(ar_fit_sum, keep = inst_in_fit, print = FALSE)
  
  # Build the result object
  res <- list(
    beta0 = beta0,
    stat = wald_result$stat,
    df1 = wald_result$df1,
    df2 = wald_result$df2,
    pvalue = wald_result$p,
    method = "Anderson-Rubin test",
    vcov_type = wald_result$vcov,
    endo_names = endo_names,
    inst_names = inst_names,
    inst_tested = inst_in_fit,
    model = object
  )
  
  class(res) <- "fixest_ar"
  
  return(res)
}


#' @rdname ar_test
#' @export
print.fixest_ar = function(x, ...) {
  
  cat("Anderson-Rubin Test\n")
  cat("-------------------\n")
  
  # Show the null hypothesis
  endo_str <- paste(x$endo_names, collapse = ", ")
  beta0_str <- paste(format(x$beta0, digits = 4), collapse = ", ")
  
  if (length(x$endo_names) == 1) {
    cat("H0: ", x$endo_names, " = ", format(x$beta0, digits = 4), "\n", sep = "")
  } else {
    cat("H0: (", endo_str, ") = (", beta0_str, ")\n", sep = "")
  }
  
  cat("\n")
  
  # Show test statistics
  cat("Stat (F): ", format(x$stat, digits = 4), "\n", sep = "")
  cat("DoF:      ", x$df1, " and ", x$df2, "\n", sep = "")
  
  # Format p-value
  if (x$pvalue < 2.2e-16) {
    cat("P-value:  < 2.2e-16 ***\n")
  } else if (x$pvalue < 0.001) {
    cat("P-value:  ", format(x$pvalue, digits = 4), " ***\n", sep = "")
  } else if (x$pvalue < 0.01) {
    cat("P-value:  ", format(x$pvalue, digits = 4), " **\n", sep = "")
  } else if (x$pvalue < 0.05) {
    cat("P-value:  ", format(x$pvalue, digits = 4), " *\n", sep = "")
  } else if (x$pvalue < 0.1) {
    cat("P-value:  ", format(x$pvalue, digits = 4), " .\n", sep = "")
  } else {
    cat("P-value:  ", format(x$pvalue, digits = 4), "\n", sep = "")
  }
  
  cat("VCOV:     ", x$vcov_type, "\n", sep = "")
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  invisible(x)
}


#' Anderson-Rubin confidence interval for IV models
#'
#' @description
#' Computes Anderson-Rubin confidence intervals by inverting the AR test.
#' These confidence intervals are robust to weak instruments.
#'
#' @param object A `fixest` object obtained from [`feols`] with an IV specification.
#' @param level Numeric scalar between 0 and 1. The confidence level. Default is 0.95.
#' @param grid Numeric vector. The grid of values over which to compute the AR test. 
#'   If `NULL` (default), a grid is automatically constructed around the point estimate.
#' @param grid_range Numeric vector of length 2. If `grid` is `NULL`, the grid is 
#'   constructed from `grid_range[1]` to `grid_range[2]`. If `NULL` (default), 
#'   the range is set to the point estimate plus/minus 5 standard errors.
#' @param n_grid Integer. The number of grid points to use if `grid` is `NULL`. 
#'   Default is 100.
#' @param vcov Versatile argument to specify the VCOV. See [`ar_test`] for details.
#' @param ... Additional arguments passed to [`ar_test`].
#'
#' @details
#' This function only supports models with a single endogenous variable. For models 
#' with multiple endogenous variables, the confidence region is multi-dimensional 
#' and cannot be easily represented.
#'
#' The AR confidence interval is constructed by inverting the AR test: it includes 
#' all values of beta0 for which the AR test does not reject at the specified level.
#'
#' Note that AR confidence intervals may be:
#' \itemize{
#'   \item Empty (if the model is severely misspecified)
#'   \item Unbounded (if instruments are very weak)
#'   \item Disconnected (consisting of multiple disjoint intervals)
#' }
#'
#' @return
#' An object of class `fixest_ar_confint` containing:
#' \item{level}{The confidence level used.}
#' \item{grid}{The grid of values tested.}
#' \item{pvalues}{The p-values at each grid point.}
#' \item{accept}{The grid values that were not rejected.}
#' \item{intervals}{A matrix with columns "lower" and "upper" giving the 
#'   confidence interval(s). Each row is a separate interval (in case of 
#'   disconnected regions).}
#' \item{endo_name}{The name of the endogenous variable.}
#' \item{point_est}{The 2SLS point estimate.}
#' \item{model}{The original `feols` IV object.}
#'
#' @seealso 
#' [`ar_test`] for the Anderson-Rubin test, [`confint.fixest`] for standard 
#' confidence intervals.
#'
#' @examples
#' # Simple IV example
#' set.seed(123)
#' n <- 500
#' z <- rnorm(n)
#' x <- 0.5 * z + rnorm(n)  # moderate instrument strength
#' y <- 1 + 0.8 * x + rnorm(n)
#' data <- data.frame(y = y, x = x, z = z)
#'
#' # IV estimation
#' est <- feols(y ~ 1 | x ~ z, data = data)
#'
#' # AR confidence interval
#' ar_ci <- ar_confint(est)
#' print(ar_ci)
#'
#' # Compare with standard confidence interval
#' confint(est)
#'
#' @export
ar_confint = function(object, level = 0.95, grid = NULL, grid_range = NULL, 
                      n_grid = 100, vcov = NULL, ...) {
  
  # Input validation
  check_arg(object, "class(fixest)")
  check_arg(level, "numeric scalar GT{0} LT{1}")
  check_arg(n_grid, "integer scalar GE{10}")
  
  # Check that object is an IV model
  if (!isTRUE(object$iv)) {
    stop("The `ar_confint` function requires an IV model. ",
         "Please provide a `feols` object estimated with an IV formula.")
  }
  
  # Currently only support single endogenous variable
  endo_names <- object$iv_endo_names
  if (length(endo_names) > 1) {
    stop("The `ar_confint` function currently only supports models with a single ",
         "endogenous variable. This model has ", length(endo_names), 
         " endogenous variables: ", paste(endo_names, collapse = ", "), ".")
  }
  
  endo_name <- endo_names[1]
  
  # Get the point estimate and SE
  # Need to handle the case where the coefficient name might have "fit_" prefix
  coef_names <- names(object$coefficients)
  endo_coef_name <- paste0("fit_", endo_name)
  if (!endo_coef_name %in% coef_names) {
    endo_coef_name <- endo_name
  }
  
  if (!endo_coef_name %in% coef_names) {
    stop("Could not find the coefficient for the endogenous variable '", 
         endo_name, "' in the model.")
  }
  
  point_est <- object$coefficients[endo_coef_name]
  
  # Get SE for grid construction
  obj_sum <- summary(object, vcov = vcov, ...)
  se_est <- se(obj_sum)[endo_coef_name]
  
  # Construct grid if not provided
  if (is.null(grid)) {
    if (is.null(grid_range)) {
      # Default: point estimate +/- 5 SEs
      grid_range <- c(point_est - 5 * se_est, point_est + 5 * se_est)
    }
    check_arg(grid_range, "numeric vector len(2) no na")
    grid <- seq(grid_range[1], grid_range[2], length.out = n_grid)
  } else {
    check_arg(grid, "numeric vector no na")
  }
  
  # Compute AR test for each grid point
  pvalues <- numeric(length(grid))
  for (i in seq_along(grid)) {
    ar_res <- ar_test(object, beta0 = grid[i], vcov = vcov, ...)
    pvalues[i] <- ar_res$pvalue
  }
  
  # Identify accepted values (where p-value > 1 - level, i.e., not rejected)
  alpha <- 1 - level
  accept <- grid[pvalues > alpha]
  
  # Find contiguous intervals
  if (length(accept) == 0) {
    # Empty confidence set
    intervals <- matrix(NA_real_, nrow = 0, ncol = 2)
    colnames(intervals) <- c("lower", "upper")
  } else {
    # Find breaks in the accepted values
    # Two consecutive grid values are in different intervals if there's a 
    # rejected value between them
    accept_idx <- which(pvalues > alpha)
    breaks <- which(diff(accept_idx) > 1)
    
    # Build intervals
    n_intervals <- length(breaks) + 1
    intervals <- matrix(NA_real_, nrow = n_intervals, ncol = 2)
    colnames(intervals) <- c("lower", "upper")
    
    start_idx <- 1
    for (i in seq_len(n_intervals)) {
      if (i <= length(breaks)) {
        end_idx <- breaks[i]
      } else {
        end_idx <- length(accept_idx)
      }
      
      interval_accept_idx <- accept_idx[start_idx:end_idx]
      intervals[i, "lower"] <- grid[min(interval_accept_idx)]
      intervals[i, "upper"] <- grid[max(interval_accept_idx)]
      
      start_idx <- end_idx + 1
    }
    
    # Check if intervals extend to grid boundaries (possibly unbounded)
    # Only check the boundaries of the extreme intervals
    # First interval's lower bound: check if it touches the grid minimum
    if (intervals[1, "lower"] == min(grid)) {
      intervals[1, "lower"] <- -Inf
    }
    # Last interval's upper bound: check if it touches the grid maximum
    if (intervals[n_intervals, "upper"] == max(grid)) {
      intervals[n_intervals, "upper"] <- Inf
    }
  }
  
  res <- list(
    level = level,
    grid = grid,
    pvalues = pvalues,
    accept = accept,
    intervals = intervals,
    endo_name = endo_name,
    point_est = unname(point_est),
    model = object
  )
  
  class(res) <- "fixest_ar_confint"
  
  return(res)
}


#' @rdname ar_confint
#' @export
print.fixest_ar_confint = function(x, ...) {
  
  cat("Anderson-Rubin Confidence Interval\n")
  cat("-----------------------------------\n")
  cat("Variable: ", x$endo_name, "\n", sep = "")
  cat("Level:    ", format(x$level * 100, digits = 2), "%\n", sep = "")
  cat("\n")
  
  cat("2SLS point estimate: ", format(x$point_est, digits = 4), "\n", sep = "")
  cat("\n")
  
  if (nrow(x$intervals) == 0) {
    cat("AR Confidence Set: EMPTY\n")
    cat("  (This may indicate model misspecification)\n")
  } else if (nrow(x$intervals) == 1) {
    lower <- x$intervals[1, "lower"]
    upper <- x$intervals[1, "upper"]
    
    if (is.infinite(lower) && is.infinite(upper)) {
      cat("AR Confidence Set: (-Inf, Inf)\n")
      cat("  (Instruments may be very weak)\n")
    } else if (is.infinite(lower)) {
      cat("AR Confidence Set: (-Inf, ", format(upper, digits = 4), "]\n", sep = "")
      cat("  (Lower bound may be unbounded - consider expanding the grid)\n")
    } else if (is.infinite(upper)) {
      cat("AR Confidence Set: [", format(lower, digits = 4), ", Inf)\n", sep = "")
      cat("  (Upper bound may be unbounded - consider expanding the grid)\n")
    } else {
      cat("AR Confidence Set: [", format(lower, digits = 4), ", ", 
          format(upper, digits = 4), "]\n", sep = "")
    }
  } else {
    cat("AR Confidence Set: (disconnected)\n")
    for (i in seq_len(nrow(x$intervals))) {
      lower <- x$intervals[i, "lower"]
      upper <- x$intervals[i, "upper"]
      
      lower_bracket <- if (is.infinite(lower)) "(" else "["
      upper_bracket <- if (is.infinite(upper)) ")" else "]"
      lower_str <- if (is.infinite(lower)) "-Inf" else format(lower, digits = 4)
      upper_str <- if (is.infinite(upper)) "Inf" else format(upper, digits = 4)
      
      cat("  Interval ", i, ": ", lower_bracket, lower_str, ", ", 
          upper_str, upper_bracket, "\n", sep = "")
    }
    cat("  (Multiple intervals may indicate misspecification or non-convex confidence region)\n")
  }
  
  invisible(x)
}
