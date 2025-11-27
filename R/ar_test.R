#----------------------------------------------#
# Author: Dianyi Yang (+ AI)
# Date creation: Wed Nov 26 23:11:00 2025
# ~: Anderson-Rubin test for IV models
#----------------------------------------------#


# =============================================================================
# INTERNAL CORE FUNCTION
# =============================================================================

#' Internal core function for Anderson-Rubin test
#'
#' @description
#' This is an internal function that computes the AR test statistic and p-value.
#' It is called by both ar_test() and ar_confint().
#'
#' @param object A `fixest` object obtained from [`feols`] with an IV specification.
#' @param beta0 Numeric vector. The hypothesized values for the endogenous variable 
#'   coefficients under the null hypothesis.
#' @param vcov Versatile argument to specify the VCOV.
#' @param ... Additional arguments passed to [`feols`].
#'
#' @return A list containing:
#' \item{stat}{The test statistic (F or Chi-squared depending on vcov).}
#' \item{df1}{Numerator degrees of freedom.}
#' \item{df2}{Denominator degrees of freedom (NA for Chi-squared).}
#' \item{pvalue}{The p-value.}
#' \item{dist}{Distribution used: "F" or "Chi-squared".}
#' \item{vcov_type}{The VCOV type used.}
#' \item{endo_names}{Names of endogenous variables.}
#' \item{inst_names}{Names of instruments.}
#' \item{inst_tested}{Instrument coefficients actually tested.}
#'
#' @keywords internal
.ar_test_core = function(object, beta0, vcov = NULL, ...) {
  
  # Get endogenous variable names and instrument names
  endo_names <- object$iv_endo_names
  inst_names <- object$iv_inst_names_xpd  # expanded instrument names
  n_endo <- length(endo_names)
  
  # Name beta0 for clarity
  names(beta0) <- endo_names
  
  # Fetch the data
  data <- fetch_data(object, "To apply AR test, ")
  
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
  
  # Determine vcov specification
  # If not specified, try to get from original object
  if (is.null(vcov)) {
    if (!is.null(object$cov.scaled)) {
      vcov <- attr(object$cov.scaled, "type")
    }
    if (is.null(vcov)) {
      vcov <- "iid"
    }
  }
  
  # Run the AR regression - inherit settings from original object where possible
  ar_fit <- feols(new_fml, data = data, weights = weights_val, vcov = vcov, 
                  notes = FALSE, ...)
  
  # Summarize to ensure we have the VCOV
  ar_fit_sum <- summary(ar_fit, vcov = vcov, ...)
  
  # Get the coefficient names in the AR regression that correspond to instruments
  ar_coef_names <- names(ar_fit_sum$coefficients)
  
  # Find which coefficients correspond to instruments
  inst_in_fit <- intersect(ar_coef_names, inst_names)
  
  if (length(inst_in_fit) == 0) {
    # Try a more flexible approach using the instrument formula variables
    inst_formula_vars <- all.vars(fml_iv[[3]])
    for (v in inst_formula_vars) {
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
  
  # Determine if vcov is iid or not for choosing distribution
  vcov_type <- wald_result$vcov
  is_iid <- identical(vcov_type, "IID") || identical(vcov_type, "iid") || 
            identical(tolower(as.character(vcov_type)), "iid")
  
  # Compute statistic and p-value based on vcov type
  df1 <- wald_result$df1
  df2 <- wald_result$df2
  
  if (is_iid) {
    # For iid: use F-distribution (exact under homoskedastic errors)
    stat <- wald_result$stat
    pvalue <- pf(stat, df1, df2, lower.tail = FALSE)
    dist <- "F"
  } else {
    # For non-iid: use Chi-squared (asymptotic)
    # Convert F-statistic to Chi-squared: chi2 = F * df1
    stat <- wald_result$stat * df1
    pvalue <- pchisq(stat, df = df1, lower.tail = FALSE)
    dist <- "Chi-squared"
    df2 <- NA_real_
  }
  
  # Return results
  list(
    stat = stat,
    df1 = df1,
    df2 = df2,
    pvalue = pvalue,
    dist = dist,
    vcov_type = vcov_type,
    endo_names = endo_names,
    inst_names = inst_names,
    inst_tested = inst_in_fit,
    is_iid = is_iid
  )
}


# =============================================================================
# EXACT AR CI FOR IID (SINGLE ENDOGENOUS VARIABLE)
# =============================================================================

#' Compute exact AR confidence interval for iid errors (single endogenous variable)
#'
#' @description
#' Uses the closed-form quadratic formula for AR confidence sets when:
#' - There is exactly one endogenous variable
#' - The vcov is iid (homoskedastic errors)
#'
#' @param object A `fixest` object with IV and single endogenous variable.
#' @param level Confidence level (e.g., 0.95).
#'
#' @return A list with:
#' \item{method}{"AR"}
#' \item{exact}{TRUE}
#' \item{level}{The confidence level}
#' \item{intervals}{Matrix with columns "lower" and "upper"}
#'
#' @keywords internal
.ar_ci_exact_iid = function(object, level = 0.95) {
  
  # Fetch data
  data <- fetch_data(object, "To compute AR CI, ")
  
  # Get y, endogenous variable, and instruments
  fml_linear <- formula(object, type = "linear")
  y <- eval(fml_linear[[2]], data)
  
  endo_mat <- model.matrix(object, data = data, type = "iv.endo", collin.rm = FALSE)
  inst_mat <- model.matrix(object, data = data, type = "iv.inst", collin.rm = FALSE)
  
  # Get exogenous controls (if any)
  exo_mat <- model.matrix(object, data = data, type = "iv.exo", collin.rm = FALSE)
  
  n <- length(y)
  L <- ncol(inst_mat)  # number of instruments
  
  # Check for fixed effects - if present, we need to demean
  has_fe <- !is.null(object$fixef_vars)
  
  if (has_fe || !is.null(exo_mat)) {
    # We need to partial out controls and/or fixed effects
    # Use demean if FE present, otherwise just residualize
    
    if (has_fe) {
      # Demean all variables with respect to fixed effects
      dm_data <- demean(object)
      
      # Extract demeaned variables - the order is: y, then X (controls), then IV parts
      # We need to reconstruct the demeaned versions
      # For simplicity, re-run the demeaning manually
      
      # Create a combined matrix of all variables to demean
      all_vars <- cbind(y, endo_mat, inst_mat)
      if (!is.null(exo_mat) && ncol(exo_mat) > 0) {
        all_vars <- cbind(all_vars, exo_mat)
      }
      
      # Use fixest's demeaning
      fixef_id <- object$fixef_id
      slope_flag <- object$slope_flag
      slope_vars <- NULL
      if (!is.null(object$slope_variables_reordered)) {
        slope_vars <- object$slope_variables_reordered
      }
      
      # Demean using the internal algorithm
      all_vars_dm <- cpp_demean(all_vars, 
                                fixef_id, 
                                weights = object$weights,
                                iterMax = object$fixef.iter,
                                diffMax = object$fixef.tol,
                                nthreads = 1L,
                                fixef.algo = demeaning_algo(internal = TRUE),
                                slope.flag = if(!is.null(slope_flag)) slope_flag else integer(0),
                                slope.vars = if(!is.null(slope_vars)) slope_vars else matrix(0, 0, 0),
                                notes = FALSE)
      
      y_star <- all_vars_dm[, 1]
      d_star <- all_vars_dm[, 2]
      z_star <- all_vars_dm[, 2 + seq_len(L)]
      
      # If there are exogenous controls, also partial them out from the residuals
      if (!is.null(exo_mat) && ncol(exo_mat) > 0) {
        exo_dm <- all_vars_dm[, (2 + L + 1):ncol(all_vars_dm), drop = FALSE]
        # Further residualize y_star, d_star, z_star with respect to exo_dm
        # Using QR decomposition
        qr_exo <- qr(exo_dm)
        y_star <- qr.resid(qr_exo, y_star)
        d_star <- qr.resid(qr_exo, d_star)
        z_star <- qr.resid(qr_exo, z_star)
      }
      
      # Adjust degrees of freedom for FE
      k <- sum(object$fixef_sizes - 1) + 1  # number of FE parameters absorbed
      if (!is.null(exo_mat)) k <- k + ncol(exo_mat)
      
    } else {
      # No FE, just partial out exogenous controls
      qr_exo <- qr(cbind(1, exo_mat))  # include intercept
      y_star <- qr.resid(qr_exo, y)
      d_star <- qr.resid(qr_exo, endo_mat[, 1])
      z_star <- qr.resid(qr_exo, inst_mat)
      k <- ncol(exo_mat) + 1  # controls + intercept
    }
    
  } else {
    # No FE, no controls - just demean for intercept
    y_star <- y - mean(y)
    d_star <- endo_mat[, 1] - mean(endo_mat[, 1])
    z_star <- sweep(inst_mat, 2, colMeans(inst_mat))
    k <- 1  # just intercept
  }
  
  # Compute the quadratic coefficients
  # S_DD = sum(D*^2), S_YY = sum(Y*^2), S_DY = sum(D* * Y*)
  S_DD <- sum(d_star^2)
  S_YY <- sum(y_star^2)
  S_DY <- sum(d_star * y_star)
  
  # P_Y = projection of Y* onto Z*, P_D = projection of D* onto Z*
  qr_z <- qr(z_star)
  P_Y <- qr.fitted(qr_z, y_star)
  P_D <- qr.fitted(qr_z, d_star)
  
  S_PDD <- sum(P_D^2)
  S_PYY <- sum(P_Y^2)
  S_PDY <- sum(P_D * P_Y)
  
  # Critical value
  alpha <- 1 - level
  Fcrit <- qf(1 - alpha, df1 = L, df2 = n - k - L)
  c_val <- Fcrit * L / (n - k - L)
  
  # Quadratic coefficients for: A * beta0^2 + B * beta0 + C >= 0
  A <- c_val * S_DD - (c_val + 1) * S_PDD
  B <- -2 * c_val * S_DY + 2 * (c_val + 1) * S_PDY
  C <- c_val * S_YY - (c_val + 1) * S_PYY
  
  # Solve the quadratic inequality
  discriminant <- B^2 - 4 * A * C
  
  if (abs(A) < 1e-14) {
    # Linear case: B * beta0 + C >= 0
    if (abs(B) < 1e-14) {
      # Constant case
      if (C >= 0) {
        # Whole real line
        intervals <- matrix(c(-Inf, Inf), nrow = 1, ncol = 2)
      } else {
        # Empty set
        intervals <- matrix(NA_real_, nrow = 0, ncol = 2)
      }
    } else {
      root <- -C / B
      if (B > 0) {
        intervals <- matrix(c(root, Inf), nrow = 1, ncol = 2)
      } else {
        intervals <- matrix(c(-Inf, root), nrow = 1, ncol = 2)
      }
    }
  } else if (discriminant < 0) {
    # No real roots
    if (A > 0) {
      # Parabola opens upward, always positive
      intervals <- matrix(c(-Inf, Inf), nrow = 1, ncol = 2)
    } else {
      # Parabola opens downward, always negative
      intervals <- matrix(NA_real_, nrow = 0, ncol = 2)
    }
  } else {
    # Two real roots
    sqrt_disc <- sqrt(discriminant)
    root1 <- (-B - sqrt_disc) / (2 * A)
    root2 <- (-B + sqrt_disc) / (2 * A)
    
    # Ensure root1 < root2
    if (root1 > root2) {
      tmp <- root1
      root1 <- root2
      root2 <- tmp
    }
    
    if (A > 0) {
      # Parabola opens upward: f(beta0) >= 0 for beta0 <= root1 or beta0 >= root2
      intervals <- matrix(c(-Inf, root1, root2, Inf), nrow = 2, ncol = 2, byrow = TRUE)
    } else {
      # Parabola opens downward: f(beta0) >= 0 for root1 <= beta0 <= root2
      intervals <- matrix(c(root1, root2), nrow = 1, ncol = 2)
    }
  }
  
  colnames(intervals) <- c("lower", "upper")
  
  list(
    method = "AR",
    exact = TRUE,
    level = level,
    intervals = intervals
  )
}


# =============================================================================
# NUMERIC AR CI FOR NON-IID VCOV
# =============================================================================

#' Compute numeric AR confidence interval by test inversion
#'
#' @description
#' Inverts the AR test numerically to find the confidence set.
#' Used when vcov is not iid or when exact computation is not possible.
#'
#' @param object A `fixest` object with IV.
#' @param level Confidence level (e.g., 0.95).
#' @param vcov VCOV specification.
#' @param n_grid Number of coarse grid points. Default 200.
#' @param tol Tolerance for root finding. Default 1e-8.
#' @param ... Additional arguments.
#'
#' @return A list with:
#' \item{method}{"AR"}
#' \item{exact}{FALSE}
#' \item{level}{The confidence level}
#' \item{intervals}{Matrix with columns "lower" and "upper"}
#'
#' @keywords internal
.ar_ci_numeric = function(object, level = 0.95, vcov = NULL, n_grid = 200, 
                          tol = 1e-8, ...) {
  
  # Get point estimate and SE for grid construction
  endo_names <- object$iv_endo_names
  endo_name <- endo_names[1]
  
  coef_names <- names(object$coefficients)
  endo_coef_name <- paste0("fit_", endo_name)
  if (!endo_coef_name %in% coef_names) {
    endo_coef_name <- endo_name
  }
  
  point_est <- object$coefficients[endo_coef_name]
  obj_sum <- summary(object, vcov = vcov, ...)
  se_est <- se(obj_sum)[endo_coef_name]
  
  # Default grid range: point estimate +/- 5 SEs
  grid_range <- c(point_est - 5 * se_est, point_est + 5 * se_est)
  grid <- seq(grid_range[1], grid_range[2], length.out = n_grid)
  
  # Get the critical value
  # First, run one AR test to determine distribution type
  test_result <- .ar_test_core(object, beta0 = point_est, vcov = vcov, ...)
  
  if (test_result$dist == "F") {
    crit <- qf(level, df1 = test_result$df1, df2 = test_result$df2)
  } else {
    crit <- qchisq(level, df = test_result$df1)
  }
  
  # Evaluate g(beta0) = stat(beta0) - crit on the coarse grid
  g_values <- numeric(n_grid)
  for (i in seq_along(grid)) {
    ar_res <- .ar_test_core(object, beta0 = grid[i], vcov = vcov, ...)
    g_values[i] <- ar_res$stat - crit
  }
  
  # Find sign changes (brackets for roots)
  sign_changes <- which(diff(sign(g_values)) != 0)
  
  if (length(sign_changes) == 0) {
    # No sign changes - either all accepted or all rejected
    if (g_values[1] <= 0) {
      # All accepted (potentially unbounded)
      intervals <- matrix(c(-Inf, Inf), nrow = 1, ncol = 2)
    } else {
      # All rejected (empty set)
      intervals <- matrix(NA_real_, nrow = 0, ncol = 2)
    }
  } else {
    # Find roots using uniroot
    roots <- numeric(length(sign_changes))
    
    g_func <- function(beta0) {
      ar_res <- .ar_test_core(object, beta0 = beta0, vcov = vcov, ...)
      ar_res$stat - crit
    }
    
    for (i in seq_along(sign_changes)) {
      idx <- sign_changes[i]
      root_result <- uniroot(g_func, 
                            interval = c(grid[idx], grid[idx + 1]),
                            tol = tol)
      roots[i] <- root_result$root
    }
    
    roots <- sort(roots)
    
    # Build intervals based on which regions have g <= 0 (accepted)
    # Check the sign at the leftmost point to determine pattern
    n_roots <- length(roots)
    
    if (g_values[1] <= 0) {
      # Left edge is accepted - first interval starts at -Inf
      interval_list <- list()
      
      # First interval: from -Inf to first root
      if (n_roots >= 1) {
        interval_list[[1]] <- c(-Inf, roots[1])
      }
      
      # Subsequent intervals: pairs of roots (entry, exit)
      if (n_roots >= 2) {
        for (j in seq(2, n_roots, by = 2)) {
          if (j + 1 <= n_roots) {
            interval_list[[length(interval_list) + 1]] <- c(roots[j], roots[j + 1])
          } else {
            # Odd root at end - extends to infinity
            interval_list[[length(interval_list) + 1]] <- c(roots[j], Inf)
          }
        }
      }
      
    } else {
      # Left edge is rejected - intervals start at roots
      interval_list <- list()
      for (j in seq(1, n_roots, by = 2)) {
        if (j + 1 <= n_roots) {
          interval_list[[length(interval_list) + 1]] <- c(roots[j], roots[j + 1])
        } else {
          # Odd number of roots - last interval extends to infinity
          interval_list[[length(interval_list) + 1]] <- c(roots[j], Inf)
        }
      }
    }
    
    if (length(interval_list) == 0) {
      intervals <- matrix(NA_real_, nrow = 0, ncol = 2)
    } else {
      intervals <- do.call(rbind, interval_list)
    }
  }
  
  colnames(intervals) <- c("lower", "upper")
  
  list(
    method = "AR",
    exact = FALSE,
    level = level,
    intervals = intervals
  )
}


# =============================================================================
# MAIN EXPORTED FUNCTION: ar_test
# =============================================================================

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
#' @param level Numeric scalar between 0 and 1. The confidence level for the 
#'   optional confidence interval. Default is 0.95.
#' @param ci Logical or NULL. Whether to compute an AR confidence interval.
#'   If `NULL` (default):
#'   - Returns CI if vcov is "iid" AND there is exactly one endogenous variable
#'   - Does not return CI otherwise
#'   If `TRUE`: Attempts to compute CI (warns if not possible with multiple endo vars)
#'   If `FALSE`: Does not compute CI
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
#' For iid (homoskedastic) errors, the AR test statistic follows an exact 
#' F-distribution. For non-iid errors (heteroskedastic, clustered, etc.), 
#' a Chi-squared asymptotic distribution is used.
#'
#' When `ci = TRUE` (or by default for iid with single endo var):
#' - For iid errors with a single endogenous variable, an exact closed-form 
#'   AR confidence interval is computed using the quadratic formula.
#' - For non-iid errors, the CI is computed by numerically inverting the AR test.
#'
#' @return
#' An object of class `fixest_ar` containing:
#' \item{beta0}{The hypothesized coefficient values tested.}
#' \item{stat}{The AR test statistic.}
#' \item{dist}{Distribution used: "F" or "Chi-squared".}
#' \item{df1}{Numerator degrees of freedom.}
#' \item{df2}{Denominator degrees of freedom (NA for Chi-squared).}
#' \item{pvalue}{The p-value of the test.}
#' \item{method}{Character string: "Anderson-Rubin test".}
#' \item{vcov_type}{The type of VCOV used for the test.}
#' \item{endo_names}{Names of the endogenous variables.}
#' \item{inst_names}{Names of the instruments.}
#' \item{model}{The original `feols` IV object.}
#' \item{ci}{(Optional) A list containing the AR confidence interval, with elements:
#'   `method`, `exact`, `level`, `intervals`.}
#'
#' @seealso 
#' [`wald`] for standard Wald tests, [`feols`] for IV estimation.
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
#' # Test H0: beta = 0 (default) - includes CI since vcov is iid
#' ar_test(est)
#'
#' # Test H0: beta = 0.8 (close to true value, should not reject)
#' ar_test(est, beta0 = 0.8)
#'
#' # Test H0: beta = 0 (far from true value, should reject)
#' ar_test(est, beta0 = 0)
#'
#' # With clustered standard errors (Chi-squared distribution)
#' data$cl <- rep(1:100, each = 10)
#' est_cl <- feols(y ~ 1 | x ~ z, data = data, vcov = ~cl)
#' ar_test(est_cl)
#'
#' # Force CI computation even for clustered SEs
#' ar_test(est_cl, ci = TRUE)
#'
#' @export
ar_test = function(object, beta0 = NULL, vcov = NULL, level = 0.95, ci = NULL, ...) {
  
  # Input validation
  check_arg(object, "class(fixest)")
  check_arg(level, "numeric scalar GT{0} LT{1}")
  
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
  
  # Get endogenous variable names
  endo_names <- object$iv_endo_names
  n_endo <- length(endo_names)
  
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
  names(beta0) <- endo_names
  
  # Run the core test
  core_result <- .ar_test_core(object, beta0 = beta0, vcov = vcov, ...)
  
  # Build the result object
  res <- list(
    beta0 = beta0,
    stat = core_result$stat,
    dist = core_result$dist,
    df1 = core_result$df1,
    df2 = core_result$df2,
    pvalue = core_result$pvalue,
    method = "Anderson-Rubin test",
    vcov_type = core_result$vcov_type,
    endo_names = endo_names,
    inst_names = core_result$inst_names,
    inst_tested = core_result$inst_tested,
    model = object
  )
  
  # Decide whether to compute CI
  compute_ci <- FALSE
  
  if (is.null(ci)) {
    # Default: compute CI only if iid and single endo var
    if (core_result$is_iid && n_endo == 1) {
      compute_ci <- TRUE
    }
  } else if (isTRUE(ci)) {
    if (n_endo > 1) {
      warning("Cannot compute AR confidence interval for models with multiple ",
              "endogenous variables. Returning test result only.")
    } else {
      compute_ci <- TRUE
    }
  }
  # If ci = FALSE, compute_ci stays FALSE
  
  # Compute CI if requested
  if (compute_ci) {
    if (core_result$is_iid) {
      # Use exact closed-form CI
      ci_result <- .ar_ci_exact_iid(object, level = level)
    } else {
      # Use numeric inversion
      ci_result <- .ar_ci_numeric(object, level = level, vcov = vcov, ...)
    }
    res$ci <- ci_result
  }
  
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
  
  # Show test statistics with distribution info
  if (x$dist == "F") {
    cat("Stat (F):     ", format(x$stat, digits = 4), "\n", sep = "")
    cat("DoF:          ", x$df1, " and ", x$df2, "\n", sep = "")
  } else {
    cat("Stat (Chi2):  ", format(x$stat, digits = 4), "\n", sep = "")
    cat("DoF:          ", x$df1, "\n", sep = "")
  }
  
  # Format p-value
  if (x$pvalue < 2.2e-16) {
    cat("P-value:      < 2.2e-16 ***\n")
  } else if (x$pvalue < 0.001) {
    cat("P-value:      ", format(x$pvalue, digits = 4), " ***\n", sep = "")
  } else if (x$pvalue < 0.01) {
    cat("P-value:      ", format(x$pvalue, digits = 4), " **\n", sep = "")
  } else if (x$pvalue < 0.05) {
    cat("P-value:      ", format(x$pvalue, digits = 4), " *\n", sep = "")
  } else if (x$pvalue < 0.1) {
    cat("P-value:      ", format(x$pvalue, digits = 4), " .\n", sep = "")
  } else {
    cat("P-value:      ", format(x$pvalue, digits = 4), "\n", sep = "")
  }
  
  cat("Distribution: ", x$dist, "\n", sep = "")
  cat("VCOV:         ", x$vcov_type, "\n", sep = "")
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  # Print CI if present
  if (!is.null(x$ci)) {
    cat("\n")
    .print_ar_ci(x$ci, x$endo_names[1])
  }
  
  invisible(x)
}


#' Internal function to print AR confidence interval
#' @keywords internal
.print_ar_ci = function(ci, endo_name) {
  
  cat("Anderson-Rubin Confidence Interval\n")
  cat("-----------------------------------\n")
  cat("Variable: ", endo_name, "\n", sep = "")
  cat("Level:    ", format(ci$level * 100, digits = 2), "%\n", sep = "")
  cat("Method:   ", if (ci$exact) "Exact (closed-form)" else "Numeric inversion", "\n", sep = "")
  cat("\n")
  
  intervals <- ci$intervals
  
  if (nrow(intervals) == 0) {
    cat("AR Confidence Set: EMPTY\n")
    cat("  (This may indicate model misspecification)\n")
  } else if (nrow(intervals) == 1) {
    lower <- intervals[1, "lower"]
    upper <- intervals[1, "upper"]
    
    if (is.infinite(lower) && is.infinite(upper)) {
      cat("AR Confidence Set: (-Inf, Inf)\n")
      cat("  (Instruments may be very weak)\n")
    } else if (is.infinite(lower)) {
      cat("AR Confidence Set: (-Inf, ", format(upper, digits = 4), "]\n", sep = "")
    } else if (is.infinite(upper)) {
      cat("AR Confidence Set: [", format(lower, digits = 4), ", Inf)\n", sep = "")
    } else {
      cat("AR Confidence Set: [", format(lower, digits = 4), ", ", 
          format(upper, digits = 4), "]\n", sep = "")
    }
  } else {
    cat("AR Confidence Set: (disconnected)\n")
    for (i in seq_len(nrow(intervals))) {
      lower <- intervals[i, "lower"]
      upper <- intervals[i, "upper"]
      
      lower_bracket <- if (is.infinite(lower)) "(" else "["
      upper_bracket <- if (is.infinite(upper)) ")" else "]"
      lower_str <- if (is.infinite(lower)) "-Inf" else format(lower, digits = 4)
      upper_str <- if (is.infinite(upper)) "Inf" else format(upper, digits = 4)
      
      cat("  Interval ", i, ": ", lower_bracket, lower_str, ", ", 
          upper_str, upper_bracket, "\n", sep = "")
    }
    cat("  (Multiple intervals may indicate non-convex confidence region)\n")
  }
}


# =============================================================================
# BACKWARD COMPATIBLE ar_confint (now a wrapper)
# =============================================================================

#' Anderson-Rubin confidence interval for IV models
#'
#' @description
#' Computes Anderson-Rubin confidence intervals by inverting the AR test.
#' These confidence intervals are robust to weak instruments.
#'
#' Note: This function is provided for backward compatibility. The recommended
#' approach is to use `ar_test(object, ci = TRUE)` which integrates the test
#' and confidence interval computation.
#'
#' @param object A `fixest` object obtained from [`feols`] with an IV specification.
#' @param level Numeric scalar between 0 and 1. The confidence level. Default is 0.95.
#' @param vcov Versatile argument to specify the VCOV. See [`ar_test`] for details.
#' @param ... Additional arguments passed to [`ar_test`].
#'
#' @details
#' This function only supports models with a single endogenous variable. For models 
#' with multiple endogenous variables, the confidence region is multi-dimensional 
#' and cannot be easily represented.
#'
#' For iid errors, an exact closed-form confidence interval is computed.
#' For non-iid errors (clustered, heteroskedastic, etc.), the confidence interval
#' is computed by numerically inverting the AR test.
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
#' \item{method}{"AR".}
#' \item{exact}{Logical, TRUE if exact closed-form was used.}
#' \item{intervals}{A matrix with columns "lower" and "upper" giving the 
#'   confidence interval(s). Each row is a separate interval (in case of 
#'   disconnected regions).}
#' \item{endo_name}{The name of the endogenous variable.}
#' \item{point_est}{The 2SLS point estimate.}
#' \item{model}{The original `feols` IV object.}
#'
#' @seealso 
#' [`ar_test`] for the Anderson-Rubin test with integrated CI computation.
#'
#' @examples
#' # Simple IV example
#' set.seed(123)
#' n <- 500
#' z <- rnorm(n)
#' x <- 0.5 * z + rnorm(n)
#' y <- 1 + 0.8 * x + rnorm(n)
#' data <- data.frame(y = y, x = x, z = z)
#'
#' # IV estimation
#' est <- feols(y ~ 1 | x ~ z, data = data)
#'
#' # AR confidence interval (exact for iid)
#' ar_confint(est)
#'
#' # With clustered SEs (numeric inversion)
#' data$cl <- rep(1:50, each = 10)
#' est_cl <- feols(y ~ 1 | x ~ z, data = data, vcov = ~cl)
#' ar_confint(est_cl)
#'
#' @export
ar_confint = function(object, level = 0.95, vcov = NULL, ...) {
  
  # Input validation
  check_arg(object, "class(fixest)")
  check_arg(level, "numeric scalar GT{0} LT{1}")
  
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
  
  # Get point estimate
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
  
  # Determine vcov type
  if (is.null(vcov)) {
    if (!is.null(object$cov.scaled)) {
      vcov_type <- attr(object$cov.scaled, "type")
    } else {
      vcov_type <- "iid"
    }
  } else {
    # Run a quick test to determine the effective vcov type
    test_result <- .ar_test_core(object, beta0 = 0, vcov = vcov, ...)
    vcov_type <- test_result$vcov_type
  }
  
  is_iid <- identical(vcov_type, "IID") || identical(vcov_type, "iid") || 
            identical(tolower(as.character(vcov_type)), "iid")
  
  # Compute CI
  if (is_iid) {
    ci_result <- .ar_ci_exact_iid(object, level = level)
  } else {
    ci_result <- .ar_ci_numeric(object, level = level, vcov = vcov, ...)
  }
  
  res <- list(
    level = level,
    method = ci_result$method,
    exact = ci_result$exact,
    intervals = ci_result$intervals,
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
  
  cat("2SLS point estimate: ", format(x$point_est, digits = 4), "\n\n", sep = "")
  .print_ar_ci(list(level = x$level, exact = x$exact, intervals = x$intervals), 
               x$endo_name)
  
  invisible(x)
}
