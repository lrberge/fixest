#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Jun 01 08:21:26 2020
# ~: Most methods are here
#----------------------------------------------#


####
#### print/summary ####
####



#' A print facility for `fixest` objects.
#'
#' This function is very similar to usual `summary` functions as it
#' provides the table of coefficients along with other information on the fit of
#' the estimation. The type of output can be customized by the user (using
#' function `setFixest_print`).
#'
#' @method print fixest
#'
#' @param x A `fixest` object. Obtained using the methods
#'   [`femlm`], [`feols`] or
#'   [`feglm`].
#' @param n Integer, number of coefficients to display. By default, only the
#'   first 8 coefficients are displayed if `x` does not come from
#'   [`summary.fixest`].
#' @param type Either `"table"` (default) to display the coefficients table
#'   or `"coef"` to display only the coefficients.
#' @param fitstat A formula or a character vector representing which fit
#'   statistic to display. The types must be valid types of the function
#'   [`fitstat`]. The default fit statistics depend on the
#'   type of estimation (OLS, GLM, IV, with/without fixed-effect). Providing the
#'   argument `fitstat` overrides the default fit statistics, you can
#'   however use the point "." to summon them back. Ex 1: `fitstat = ~ . + ll` adds the log-likelihood
#'   to the default values. Ex 2: `fitstat = ~ ll + pr2` only displays the log-likelihood and the pseudo-R2.
#' @param ... Other arguments to be passed to [`vcov.fixest`].
#'
#' @details It is possible to set the default values for the arguments
#' `type` and `fitstat` by using the function `setFixest_print`.
#'
#' @seealso See also the main estimation functions [`femlm`],
#' [`feols`] or [`feglm`]. Use
#' [`summary.fixest`] to see the results with the appropriate
#' standard-errors, [`fixef.fixest`] to extract the
#' fixed-effects coefficients, and the function [`etable`] to
#' visualize the results of multiple estimations.
#'
#' @author Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' #   => we account for 3 fixed-effects (FEs)
#' est_pois = fepois(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # displaying the results
#' #  (by default SEs are clustered if FEs are used)
#' print(est_pois)
#'
#' # By default the coefficient table is displayed.
#' #  If the user wished to display only the coefficents, use option type:
#' print(est_pois, type = "coef")
#'
#' # To permanently display coef. only, use setFixest_print:
#' setFixest_print(type = "coef")
#' est_pois
#' # back to default:
#' setFixest_print(type = "table")
#'
#' #
#' # fitstat
#' #
#'
#' # We modify which fit statistic to display
#' print(est_pois, fitstat = ~ . + lr)
#'
#' # We add the LR test to the default (represented by the ".")
#'
#' # to show only the LR stat:
#' print(est_pois, fitstat = ~ . + lr.stat)
#'
#' # To modify the defaults:
#' setFixest_print(fitstat = ~ . + lr.stat + rmse)
#' est_pois
#'
#' # Back to default (NULL == default)
#' setFixest_print(fitstat = NULL)
#'
#'
print.fixest = function(x, n, type = "table", fitstat = NULL, ...){

  # checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = c("n", "type", "vcov"),
                  valid_args = c("vcov", "se", "cluster", "ssc", "forceCovariance", "keepBounded"))
  }

  # The objects from the estimation and the summary are identical, except regarding the vcov
  from_summary = isTRUE(x$summary)

  if(!missnull(fitstat)){
    fitstat = fitstat_validate(fitstat, TRUE)
  }

  # User options
  set_defaults("fixest_print")

  # if NOT from summary, we consider the argument 'type'
  if(!from_summary){
    # checking argument type
    check_set_arg(type, "match(coef, table)")

    if(type == "coef"){
      print(coef(x))
      return(invisible())
    }
  }

  isNegbin = x$method == "fenegbin" || (x$method_type == "feNmlm" && x$family == "negbin")

  x = summary(x, fromPrint = TRUE, ...)

  check_arg(n, "integer scalar GE{1}")

  msgRemaining = ""
  nb_coef = length(coef(x)) - isNegbin
  if(missing(n) && is.null(x$n_print)){
    if(from_summary && !isTRUE(x$summary_from_fit)){
      n = Inf
    } else {
      if(nb_coef <= 15){
        n = 15
      } else {
        n = 8
        msgRemaining = sma("... {nb_coef - n} coefficients remaining (display them with summary() or use argument n)\n")
      }
    }

  } else {
    if(!is.null(x$n_print)) n = x$n_print

    if(n < nb_coef){
      msgRemaining = sma("... {nb_coef - n} coefficients remaining\n")
    }
  }

  # We also add the collinearity message
  collinearity_msg = ""
  if(!is.null(x$collin.var)){
    n_collin = length(x$collin.var)
    collinearity_msg = sma("... {n_collin} variable{#s, were} removed because of collinearity ({enum.3 ? x$collin.var}{&n_collin>3; [full set in $collin.var]})\n")
    if(isTRUE(x$iv) && any(grepl("^fit_", x$collin.var))){
      if(!any(grepl("^fit_", names(x$coefficients)))){
        iv_msg = "NOTE: all endogenous regressors were removed.\n"
      } else {
        n_rm = sum(grepl("^fit_", x$collin.var))
        iv_msg = paste0("Important note: ", n_letter(n_rm), " endogenous regressor", plural(n_rm, "s.was"), " removed => IV estimation not valid.\n")
      }

      collinearity_msg = paste0(collinearity_msg, iv_msg)
    }
  }

  if(isFALSE(x$convStatus)){
    last_warn = getOption("fixest_last_warning")
    if(is.null(last_warn) || (proc.time() - last_warn)[3] > 1){
      if(x$method_type == "feNmlm"){
        warning("The optimization algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
      } else if(x$method_type == "feols"){
        warning("The demeaning algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
      } else {
        warning("The GLM algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
      }
    }

  }

  coeftable = x$coeftable

  # The type of SE
  se.type = attr(coeftable, "type")
  if(is.null(se.type)) se.type = "Custom"

  if(x$method_type == "feNmlm"){
    family_format = c(poisson="Poisson", negbin="Negative Binomial", logit="Logit", gaussian="Gaussian")
    msg = ifelse(is.null(x$call$NL.fml), "", "Non-linear ")
    half_line = paste0(msg, "ML estimation, family = ", family_format[x$family])
  } else if(x$method %in% c("feglm", "feglm.fit")) {
    fam_call = x$call$family
    if(is.null(names(fam_call))){
      half_line = paste0("GLM estimation, family = ", x$family$family)
    } else {
      half_line = paste0("GLM estimation, family = ", deparse_long(fam_call))
    }
  } else if(x$method == "fepois") {
    half_line = "Poisson estimation"
  } else if(x$method == "fenegbin") {
    half_line = "Negative Binomial ML estimation"
  } else {
    half_line = "OLS estimation"
  }

  if(isTRUE(x$iv)){
    glue = function(...) paste(..., collapse = ", ")
    first_line = paste0("TSLS estimation, Dep. Var.: ", as.character(x$fml)[[2]], ", Endo.: ", glue(get_vars(x$iv_endo_fml)), ", Instr.: ", glue(x$iv_inst_names), "\n")
    second_line = paste0(ifunit(x$iv_stage, "First", "Second"), " stage: Dep. Var.: ", as.character(x$fml)[[2]], "\n")
    cat(first_line, second_line, sep = "")
  } else {
    cat(half_line, ", Dep. Var.: ", as.character(x$fml)[[2]], "\n", sep="")
  }


  cat("Observations:", addCommas(x$nobs), "\n")

  extra_info = c("subset", "sample", "offset", "weights")
  for(i in seq_along(extra_info)){
    nm = extra_info[i]
    if(nm %in% names(x$model_info)){
      xtra = x$model_info[[extra_info[i]]]

      if(length(xtra) == 1){
        cat(dsb(".[u?nm]:"), xtra, "\n")
      } else {
        if(xtra$value != "Full sample"){
          cat(dsb(".[u?nm] (.[xtra$var]): .[xtra$value]\n"))
        }
      }
    }
  }

  if(!is.null(x$fixef_terms)){
    terms_full = extract_fe_slope(x$fixef_terms)
    fixef_vars = terms_full$fixef_vars

    if(length(fixef_vars) > 0){
      cat("Fixed-effects: ", paste0(fixef_vars, ": ", addCommas(x$fixef_sizes[fixef_vars]), collapse=",  "), "\n", sep = "")
    }

    cat("Varying slopes: ", paste0(terms_full$slope_vars, " (", terms_full$slope_fe, ": ", addCommas(x$fixef_sizes[terms_full$slope_fe]), ")", collapse = ",  "), "\n", sep = "")

  } else {
    if(!is.null(x$fixef_sizes)) cat("Fixed-effects: ", paste0(x$fixef_vars, ": ", addCommas(x$fixef_sizes), collapse = ",  "), "\n", sep = "")
  }


  if(is.null(x$onlyFixef)){

    cat("Standard-errors:", se.type, "\n")

    last_line = paste0(msgRemaining, collinearity_msg)

    # The matrix of coefficients
    if(isNegbin){
      if(nrow(coeftable) == 2){
        new_table = coeftable[1, , drop = FALSE]
      } else {
        new_table = coeftable[-nrow(coeftable), ]
      }

      print_coeftable(head(new_table, n), lastLine = last_line)

      theta = coeftable[".theta", 1]
      noDispInfo = ifelse(theta > 1000, "(theta >> 0, no sign of overdispersion, you may consider a Poisson model)", "")
      cat("Over-dispersion parameter: theta =", theta, noDispInfo, "\n")
    } else {
      print_coeftable(head(coeftable, n), lastLine = last_line)
    }
  }

  if(isTRUE(x$NA_model)){
    return(invisible())
  }

  if(!is.null(fitstat) && identical(fitstat, NA)){
    # No fitstat

  } else {

    if(is.null(fitstat) || "." %in% fitstat){
      if(x$method_type == "feols"){
        default_fit = c("rmse", "ar2")

        if(!is.null(x$fixef_sizes) && is.null(x$onlyFixef)){
          default_fit = c(default_fit, "wr2")
        }

        if(isTRUE(x$iv)){
          default_fit = c(default_fit, "ivf1", "wh", "sargan")
        }

      } else {
        default_fit = c("ll", "apr2", "bic", "cor2")
      }

      if("." %in% fitstat){
        fitstat = setdiff(c(default_fit, fitstat), ".")
      } else {
        fitstat = default_fit
      }
    }

    print(fitstat(x, fitstat), na.rm = TRUE, group.solo = TRUE)
  }

  if(isFALSE(x$convStatus)){
    iter_format = x$iterations
    if(length(iter_format)== 1){
      iter_format = paste0("lhs: ", iter_format)
    } else {
      n_iter = length(iter_format)
      iter_format = paste0("lhs: ", iter_format[n_iter], ", rhs: ", paste0(head(iter_format, min(n_iter - 1, n)), collapse = ", "))
    }
    cat("# Evaluations:", iter_format, "--", x$message, "\n")
  }

}

##

#' Summary of a `fixest` object. Computes different types of standard errors.
#'
#' This function is similar to `print.fixest`. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
#'
#' @inheritParams feNmlm
#' @inheritParams aggregate.fixest
#'
#' @method summary fixest
#' @param vcov Versatile argument to specify the VCOV. In general, it is either a character scalar equal to a VCOV type, either a formula of the form: `vcov_type ~ variables`. The VCOV types implemented are: "iid", "hetero" (or "HC1"), "cluster", "twoway", "NW" (or "newey_west"), "DK" (or "driscoll_kraay"), and "conley". It also accepts object from [`vcov_cluster`], [`vcov_NW`][fixest::vcov_hac], [`NW`][fixest::vcov_hac], [`vcov_DK`][fixest::vcov_hac], [`DK`][fixest::vcov_hac], [`vcov_conley`] and [`conley`][fixest::vcov_conley]. It also accepts covariance matrices computed externally. Finally it accepts functions to compute the covariances. See the `vcov` documentation in the [vignette](https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1).
#' @param se Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: `se = "cluster"`, otherwise `se = "iid"`. Note that this argument is deprecated, you should use `vcov` instead.
#' @param cluster Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over `var1` and `var2` contained in the data.frame `base` used for the estimation. All the following `cluster` arguments are valid and do the same thing: `cluster = base[, c("var1", "var2")]`, `cluster = c("var1", "var2")`, `cluster = ~var1+var2`. If the two variables were used as fixed-effects in the estimation, you can leave it blank with `vcov = "twoway"` (assuming `var1` \[resp. `var2`\] was the 1st \[resp. 2nd\] fixed-effect). You can interact two variables using `^` with the following syntax: `cluster = ~var1^var2` or `cluster = "var1^var2"`.
#' @param stage Can be equal to `2` (default), `1`, `1:2` or `2:1`. Only used if the object is an IV estimation: defines the stage to which `summary` should be applied. If `stage = 1` and there are multiple endogenous regressors or if `stage` is of length 2, then an object of class `fixest_multi` is returned.
#' @param object A `fixest` object. Obtained using the functions [`femlm`], [`feols`] or [`feglm`].
#' @param ssc An object of class `ssc.type` obtained with the function [`ssc`]. Represents how the degree of freedom correction should be done.You must use the function [`ssc`] for this argument. The arguments and defaults of the function [`ssc`] are: `adj = TRUE`, `fixef.K="nested"`, `cluster.adj = TRUE`, `cluster.df = "min"`, `t.df = "min"`, `fixef.force_exact=FALSE)`. See the help of the function [`ssc`] for details.
#' @param .vcov A user provided covariance matrix or a function computing this matrix. If a matrix, it must be a square matrix of the same number of rows as the number of variables estimated. If a function, it must return the previously mentioned matrix.
#' @param lean Logical, default is `FALSE`. Used to reduce the (memory) size of the summary object. If `TRUE`, then all objects of length N (the number of observations) are removed from the result. Note that some `fixest` methods may consequently not work when applied to the summary.
#' @param forceCovariance (Advanced users.) Logical, default is `FALSE`. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option to force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded (Advanced users -- `feNmlm` with non-linear part and bounded coefficients only.) Logical, default is `FALSE`. If `TRUE`, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param vcov_fix Logical scalar, default is `TRUE`. If the VCOV ends up not being positive definite, whether to "fix" it using an eigenvalue decomposition (a la Cameron, Gelbach & Miller 2011).
#' @param n Integer, default is 1000. Number of coefficients to display when the print method is used.
#' @param ... Only used if the argument `.vocv` is provided and is a function: extra arguments to be passed to that function.
#'
#' @section Compatibility with \pkg{sandwich} package:
#' The VCOVs from `sandwich` can be used with `feols`, `feglm` and `fepois` estimations. If you want to have a `sandwich` VCOV when using `summary.fixest`, you can use the argument `vcov` to specify the VCOV function to use (see examples).
#' Note that if you do so and you use a formula in the `cluster` argument, an innocuous warning can pop up if you used several non-numeric fixed-effects in the estimation (this is due to the function [`expand.model.frame`] used in `sandwich`).
#'
#' @return
#' It returns a `fixest` object with:
#' \item{cov.scaled}{The new variance-covariance matrix (computed according to the argument `se`).}
#' \item{se}{The new standard-errors (computed according to the argument `se`).}
#' \item{coeftable}{The table of coefficients with the new standard errors.}
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. Use [`fixef.fixest`] to extract the fixed-effects coefficients, and the function [`etable`] to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 fixed-effects)
#' est_pois = fepois(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Comparing different types of standard errors
#' sum_standard = summary(est_pois, vcov = "iid")
#' sum_hetero   = summary(est_pois, vcov = "hetero")
#' sum_oneway   = summary(est_pois, vcov = "cluster")
#' sum_twoway   = summary(est_pois, vcov = "twoway")
#'
#' etable(sum_standard, sum_hetero, sum_oneway, sum_twoway)
#'
#' # Alternative ways to cluster the SE:
#' summary(est_pois, vcov = cluster ~ Product + Origin)
#' summary(est_pois, vcov = ~Product + Origin)
#' summary(est_pois, cluster = ~Product + Origin)
#'
#' # You can interact the clustering variables "live" using the var1 ^ var2 syntax.#'
#' summary(est_pois, vcov = ~Destination^Product)
#'
#' #
#' # Newey-West and Driscoll-Kraay SEs
#' #
#'
#' data(base_did)
#' # Simple estimation on a panel
#' est = feols(y ~ x1, base_did)
#'
#' # --
#' # Newey-West
#' # Use the syntax NW ~ unit + time
#' summary(est, NW ~ id + period)
#'
#' # Now take a lag of 3:
#' summary(est, NW(3) ~ id + period)
#'
#' # --
#' # Driscoll-Kraay
#' # Use the syntax DK ~ time
#' summary(est, DK ~ period)
#'
#' # Now take a lag of 3:
#' summary(est, DK(3) ~ period)
#'
#' #--
#' # Implicit deductions
#' # When the estimation is done with a panel.id, you don't need to
#' # specify these values.
#'
#' est_panel = feols(y ~ x1, base_did, panel.id = ~id + period)
#'
#' # Both methods, NM and DK, now work automatically
#' summary(est_panel, "NW")
#' summary(est_panel, "DK")
#'
#' #
#' # VCOVs robust to spatial correlation
#' #
#'
#' data(quakes)
#' est_geo = feols(depth ~ mag, quakes)
#'
#' # --
#' # Conley
#' # Use the syntax: conley(cutoff) ~ lat + lon
#' # with lat/lon the latitude/longitude variable names in the data set
#' summary(est_geo, conley(100) ~ lat + long)
#'
#' # Change the cutoff, and how the distance is computed
#' summary(est_geo, conley(200, distance = "spherical") ~ lat + long)
#'
#' # --
#' # Implicit deduction
#' # By default the latitude and longitude are directly fetched in the data based
#' # on pattern matching. So you don't have to specify them.
#' # Further an automatic cutoff is computed by default.
#'
#' # The following works
#' summary(est_geo, "conley")
#'
#'
#'
#' #
#' # Compatibility with sandwich
#' #
#'
#' # You can use the VOCVs from sandwich by using the argument .vcov:
#' library(sandwich)
#' summary(est_pois, .vcov = vcovCL, cluster = trade[, c("Destination", "Product")])
#'
#'
summary.fixest = function(object, vcov = NULL, cluster = NULL, ssc = NULL, .vcov = NULL,
              stage = NULL, lean = FALSE, agg = NULL, forceCovariance = FALSE,
              se = NULL, keepBounded = FALSE, n = 1000, vcov_fix = TRUE,
              nthreads = getFixest_nthreads(), ...){

  # computes the clustered SEs and returns the modified vcov and coeftable
  # NOTA: if the object is already a summary

  if(isTRUE(object$onlyFixef) || isTRUE(object$NA_model)){
    # means that the estimation is done without variables
    return(object)
  }

  mc = match.call()

  dots = list(...)

  check_arg(n, "integer scalar GE{1}")
  if(!missing(n)){
    object$n_print = n
  }

  # we need this to save the summary flags
  # All three arguments se+cluster+.vcov are formatted into a valid vcov arg.
  vcov_in = vcov = oldargs_to_vcov(se, cluster, vcov, .vcov)

  check_arg(lean, "logical scalar")
  check_arg(stage, "NULL integer vector no na len(,2) GE{1} LE{2}")

  skip_vcov = FALSE
  if(isTRUE(object$summary)){

    if("fromPrint" %in% names(dots)){
      # From print
      return(object)

    } else if(is.null(vcov) && is.null(ssc)){
      # We return directly the object ONLY if not any other argument has been passed

      skip_vcov = TRUE
      if(is.null(agg) && is.null(stage)){
        if(!lean || (lean && isTRUE(object$lean))){
          # No modification required
          object$summary_from_fit = FALSE
          return(object)
        }
      }
    }

    assign_flags(object$summary_flags, vcov = vcov, ssc = ssc, agg = agg)
  }

  # Checking arguments in ...
  if(is_user_level_call()){
    if(!is_function_in_it(vcov)){
      validate_dots(suggest_args = c("se", "cluster", "ssc"), valid_args = "dof")
    }
  }

  if(is.null(stage)) stage = 2


  # IV
  if(isTRUE(object$iv) && !isTRUE(dots$iv)){
    stage = unique(stage)
    res = list()

    # if lean, we still compute the summary for the first stage,
    #  then we will inject it in the iv_first_stage object of the 2nd stage
    # => this is critical to get the right Wald stat (of the 1st stage),
    #  otherwise it won't be possible to get it.
    remove_stage_1 = FALSE
    if(lean && !1 %in% stage){
      remove_stage_1 = TRUE
      stage = 1:2
    }

    stage_names = c()

    for(s in seq_along(stage)){
      if(stage[s] == 1){
        for(i in seq_along(object$iv_first_stage)){
          res[[length(res) + 1]] = summary(object$iv_first_stage[[i]],
                           vcov = vcov, ssc = ssc, lean = lean,
                           forceCovariance = forceCovariance,
                           vcov_fix = vcov_fix,
                           n = n, nthreads = nthreads, iv = TRUE)

          stage_names[length(stage_names) + 1] = paste0("First stage: ", names(object$iv_first_stage)[i])
        }

      } else {
        # We keep the information on clustering => matters for wald tests of 1st stage
        my_res = summary(object, vcov = vcov, ssc = ssc, lean = lean,
                 forceCovariance = forceCovariance, vcov_fix = vcov_fix,
                 n = n, nthreads = nthreads, iv = TRUE)

        res[[length(res) + 1]] = my_res
        stage_names[length(stage_names) + 1] = "Second stage"
      }
    }

    if(lean && 2 %in% stage){
      # we inject the summary of the first stage into the iv_first_stage
      qui_1st = which(grepl("^First", stage_names))
      qui_2nd = which(stage_names == "Second stage")

      tmp_1st = res[qui_1st]
      names(tmp_1st) = names(object$iv_first_stage)

      res[[qui_2nd]][["iv_first_stage"]] = tmp_1st
    }

    if(remove_stage_1){
      qui_2nd = which(stage_names == "Second stage")
      return(res[[qui_2nd]])
    }

    if(length(res) == 1){
      return(res[[1]])
    }

    values = list("iv" = stage_names)
    res_multi = setup_multi(res, values)
    attr(res_multi, "print_request") = "long"

    return(res_multi)
  }


  # The new VCOV
  if(skip_vcov){
    vcov = object$cov.scaled

  } else {
    vcov = vcov.fixest(object, vcov = vcov, ssc = ssc, forceCovariance = forceCovariance,
               vcov_fix = vcov_fix,
               keepBounded = keepBounded, nthreads = nthreads,
               attr = TRUE, se = se, cluster = cluster, ...)
  }

  # NOTA:
  # I need to add se and cluster even if they're not needed only to ensure it
  # works fine when vcov is a function and cluster/se are arguments

  sd2 = diag(vcov)
  sd2[sd2 < 0] = NA
  se = sqrt(sd2)

  # used to handle the case of bounded parameters
  params = names(object$coefficients)
  if(length(se) != length(params)){
    se = se[params]
  }
  names(se) = params

  # The coeftable is modified accordingly
  coeftable = object$coeftable

  # th z & p values
  zvalue = object$coefficients/se
  pvalue = fixest_pvalue(object, zvalue, vcov)

  # update of se if bounded
  se_format = se
  isBounded = object$isBounded
  if(!is.null(isBounded) && any(isBounded)){
    if(!keepBounded){
      se_format[!isBounded] = decimalFormat(se_format[!isBounded])
      se_format[isBounded] = attr(isBounded, "type")
    }
  }

  # modifs of the table
  coeftable = cbind("Estimate" = object$coefficients, "Std. Error" = se_format,
            "t value" = zvalue, "Pr(>|t|)" = pvalue)
  if(object$method != "feols"){
    colnames(coeftable) = colnames(object$coeftable)
  }

  attr(coeftable, "type") = attr(se, "type") = attr(vcov, "type")

  object$cov.scaled = vcov
  object$coeftable = coeftable
  object$se = se

  if(lean){
    var2clean = c("fixef_id", "residuals", "fitted.values", "scores", "sumFE",
            "slope_variables_reordered", "y", "weights", "irls_weights",
            "obs_selection", "iv_residuals", "fitted.values_demean",
            "working_residuals", "linear.predictors")

    object[var2clean] = NULL

    object$lean = TRUE

    # We also clean the environments from the formulas
    for(i in seq_along(object$fml_all)){
      attr(object$fml_all[[i]], ".Environment") = .GlobalEnv
    }
    attr(object$fml, ".Environment") = .GlobalEnv

  }

  object$summary = TRUE

  # We save the arguments used to construct the summary
  if("summary_flags" %in% names(dots)){
    # If here => this is a call from an estimation (=fit)
    object$summary_flags = dots$summary_flags
    object$summary_from_fit = TRUE
  } else {
    # build_flags does not accept missing arguments
    if(missing(ssc)) ssc = NULL

    if(lean){

      size_KB = as.numeric(object.size(vcov)) / 8 / 1000

      if(size_KB > 100){
        # Here => means the user has manually provided a cluster => will be of size N at least
        # To respect lean = TRUE we keep no memory of this choice
        vcov_in = NULL
      }

    }

    object$summary_flags = build_flags(mc, vcov = vcov_in, ssc = ssc)
    object$summary_from_fit = NULL
  }

  # agg
  if(!missnull(agg)){
    agg_result = aggregate(object, agg, full = TRUE, from_summary = TRUE)
    object$coeftable = agg_result$coeftable
    object$model_matrix_info = agg_result$model_matrix_info
    object$is_agg = TRUE
  }

  return(object)
}


#' @rdname summary.fixest
summary.fixest_list = function(object, se, cluster, ssc = getFixest_ssc(), .vcov, stage = 2, lean = FALSE, n, ...){

  dots = list(...)

  res = list()
  for(i in seq_along(object)){
    my_res = summary(object[[i]], se = se, cluster = cluster, ssc = ssc, .vcov = .vcov, stage = stage, lean = lean, n = n)

    # we unroll in case of IV
    if("fixest_multi" %in% class(my_res)){
      for(j in seq_along(my_res)){
        res[[length(res) + 1]] = my_res[[j]]
      }
    } else {
      res[[length(res) + 1]] = my_res
    }
  }

  # We return a simple list
  class(res) = NULL

  res
}

#' Summary method for fixed-effects coefficients
#'
#' This function summarizes the main characteristics of the fixed-effects coefficients. It shows the number of fixed-effects that have been set as references and the first elements of the fixed-effects.
#'
#' @method summary fixest.fixef
#'
#' @param object An object returned by the function [`fixef.fixest`].
#' @param n Positive integer, defaults to 5. The `n` first fixed-effects for each fixed-effect dimension are reported.
#' @param ... Not currently used.
#'
#' @return
#' It prints the number of fixed-effect coefficients per fixed-effect dimension, as well as the number of fixed-effects used as references for each dimension, and the mean and variance of the fixed-effect coefficients. Finally, it reports the first 5 (arg. `n`) elements of each fixed-effect.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' [`femlm`], [`fixef.fixest`], [`plot.fixest.fixef`].
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' # => we account for 3 fixed-effects effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # obtaining the fixed-effects coefficients
#' fe_trade = fixef(est_pois)
#'
#' # printing some summary information on the fixed-effects coefficients:
#' summary(fe_trade)
#'
#'
summary.fixest.fixef = function(object, n = 5, ...){
  # This function shows some generic information on the fixed-effect coefficients

  # checking arguments in dots
  if(is_user_level_call()){
    validate_dots(suggest_args = "n")
  }

  Q = length(object)
  fixef_names = names(object)
  slope_flag = grepl("\\[", fixef_names)
  fe = gsub("\\[.+", "", fixef_names)
  slope = gsub(".+\\[|\\].+", "", fixef_names)

  isSlope = any(slope_flag)
  isFE = any(!slope_flag)
  info = as.character(10*isFE + isSlope)

  # we rework the names
  fixef_names[slope_flag] = paste0(slope[slope_flag], " (slopes: ", fe[slope_flag], ")")

  isRegular = TRUE
  if(Q > 1){
    nb_ref = attr(object, "references")
    nb_per_cluster = sapply(object, length)
    mean_per_cluster = sd_per_cluster = c()
    for(i in 1:Q){
      mean_per_cluster[i] = as.character(signif(mean(object[[i]]), 3))
      sd_per_cluster[i] = as.character(signif(sd(object[[i]]), 3))
    }
    res = as.data.frame(rbind(nb_per_cluster, nb_ref, mean_per_cluster, sd_per_cluster))

    row_1 = paste0("Number of ", switch(info, "11" = "fixed-effects/slopes", "10"="fixed-effects", "1"="slopes"))

    rownames(res) = c(row_1, "Number of references", "Mean", "Standard-deviation")

    colnames(res) = fixef_names

    if(sum(nb_ref) > Q-1){
      isRegular = FALSE
    }
  }

  # The message

  my_title = paste0(switch(info, "11" = "Fixed-effects/Slope", "10"="Fixed_effects", "1"="Slope"), " coefficients\n")
  cat(my_title)
  if(Q == 1){
    x1 = object[[1]]
    if(slope_flag){
      cat("Number of slope coefficients for variable ", slope, " (slope: ", fe, ") is ", length(x1), ".\n", sep = "")
    } else {
      cat("Number of fixed-effects for variable ", fixef_names, " is ", length(x1), ".\n", sep = "")
    }

    cat("\tMean = ", signif(mean(x1), 3), "\tVariance = ", signif(var(x1), 3), "\n", sep = "")
  } else {
    print(res)
  }

  # We print the first 5 elements of each fixed-effect
  cat("\nCOEFFICIENTS:\n")
  for(i in 1:Q){
    m = head(object[[i]], n)

    m_char = as.data.frame(t(as.data.frame(c("", as.character(signif(m, 4))))))
    names(m_char) = c(paste0(fixef_names[i], ":"), names(m))
    rownames(m_char) = " "

    n_cluster = length(object[[i]])
    if(n_cluster > n){
      m_char[["   "]] = paste0("... ", addCommas(n_cluster - n), " remaining")
    }

    print(m_char)
    if(i != Q) cat("-----\n")
  }

}

####
#### fixef ####
####


#' Extract the Fixed-Effects from a `fixest` estimation.
#'
#' This function retrieves the fixed effects from a `fixest` estimation. It is useful only when there are one or more fixed-effect dimensions.
#'
#' @inheritParams feNmlm
#'
#' @param object A `fixest` estimation (e.g. obtained using [`feols`] or [`feglm`]).
#' @param notes Logical. Whether to display a note when the fixed-effects coefficients are not regular.
#' @param sorted Logical, default is `TRUE`. Whether to order the fixed-effects by their names. If `FALSE`, then the order used in the demeaning algorithm is used.
#'
#' @details
#' If the fixed-effect coefficients are not regular, then several reference points need to be set: this means that the fixed-effects coefficients cannot be directly interpreted. If this is the case, then a warning is raised.
#'
#' @return
#' A list containing the vectors of the fixed effects.
#'
#' If there is more than 1 fixed-effect, then the attribute \dQuote{references} is created. This is a vector of length the number of fixed-effects, each element contains the number of coefficients set as references. By construction, the elements of the first fixed-effect dimension are never set as references. In the presence of regular fixed-effects, there should be Q-1 references (with Q the number of fixed-effects).
#'
#' @seealso
#' [`plot.fixest.fixef`]. See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. Use [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effect coefficients, and the function [`etable`] to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 fixed-effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Obtaining the fixed-effects coefficients:
#' fe_trade = fixef(est_pois)
#'
#' # The fixed-effects of the first fixed-effect dimension:
#' head(fe_trade$Origin)
#'
#' # Summary information:
#' summary(fe_trade)
#'
#' # Plotting them:
#' plot(fe_trade)
#'
fixef.fixest = function(object, notes = getFixest_notes(), sorted = TRUE, nthreads = getFixest_nthreads(),
            fixef.tol = 1e-5, fixef.iter = 10000, ...){

  # object is a fixest object
  # This function retrieves the dummies

  check_arg(notes, sorted, "logical scalar")

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots()
  }

  check_value(fixef.tol, "numeric scalar GT{0} LT{1}")
  check_value(fixef.iter, "strict integer scalar GT{0}")

  if(isTRUE(object$lean)){
    # LATER: recompute the FEs by extracting them from the data
    stop("Fixed-effects from 'lean' fixest objects cannot be extracted. Please re-estimate with 'lean = FALSE'.")
  }

  # Preliminary stuff
  S = object$sumFE

  if(is.null(S)){
    stop("The estimation was done without fixed-effects (FE). The FE coefficients cannot be retrieved.")
  }

  family = object$family
  fixef_names = object$fixef_vars

  fixef_id = object$fixef_id

  Q = length(fixef_id)
  N = length(S)

  # either (we need to clean its attributes for unlist to be efficient)
  id_dummies_vect = list()
  for(i in 1:Q) id_dummies_vect[[i]] = as.vector(fixef_id[[i]])

  is_ref_approx = FALSE
  isSlope = FALSE
  if(!is.null(object$fixef_terms)){
    isSlope = TRUE
    # This is an estimation with slopes
    # we apply another method => we use the demeaning function

    slope_variables = object$slope_variables_reordered
    slope_flag = object$slope_flag_reordered

    new_order = object$fe.reorder
    fixef_vars = object$fixef_vars[new_order]
    fixef_sizes = as.integer(object$fixef_sizes[new_order])

    # We reconstruct the terms
    fixef_terms = c()
    start = c(0, cumsum(abs(slope_flag)))
    for(i in seq_along(slope_flag)){
      sf = slope_flag[i]
      if(sf >= 0){
        fixef_terms = c(fixef_terms, fixef_vars[i])
      }

      if(abs(sf) > 0){
        fixef_terms = c(fixef_terms, paste0(fixef_vars[i], "[[", names(slope_variables)[start[i] + 1:abs(sf)], "]]"))
      }
    }

    fe_id_list = object$fixef_id[new_order]

    #
    # STEP 2: demeaning
    #

    nthreads = check_set_nthreads(nthreads)

    table_id_I = as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

    S_demean = cpp_demean(y = S, X_raw = 0, r_weights = 0, iterMax = as.integer(fixef.iter),
                diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
                fe_id_list = fe_id_list, table_id_I = table_id_I,
                slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
                r_init = 0, nthreads = nthreads, save_fixef = TRUE)

    fixef_coef = S_demean$fixef_coef

    names(fixef_sizes) = fixef_vars

    fe_all = c()
    for(i in seq_along(slope_flag)){
      fe_all = c(fe_all, rep(fixef_vars[i], 1 + abs(slope_flag[i]) - (slope_flag[i] < 0)))
    }

    start = 1
    i = 1
    fixef_values = list()
    for(q in seq_along(slope_flag)){
      sf = slope_flag[q]
      if(sf == 0){
        fixef_values[[i]] = fixef_coef[seq(start, length.out = fixef_sizes[q])]
        i = i + 1
        start = start + fixef_sizes[q]
      } else {
        nb = abs(sf) + (sf > 0)

        adj = 0
        if(sf > 0){
          # The fixed-effects is in the last position
          j_fe = nb - 1
          fixef_values[[i]] = fixef_coef[seq(start + j_fe, by = nb, length.out = fixef_sizes[q])]
          adj = 1
        }

        for(j in 0:(nb - 1 - adj)){
          fixef_values[[i + j + adj]] = fixef_coef[seq(start + j, by = nb, length.out = fixef_sizes[q])]
        }
        i = i + nb
        start = start + fixef_sizes[q] * nb
      }

    }

    #
    # Now the referenes
    #

    nb_ref = integer(length(fixef_terms))

    # FE references
    who_fe = slope_flag >= 0
    Q_fe = sum(who_fe)
    if(Q_fe >= 2){

      my_dum = fe_id_list[who_fe]

      dumMat = matrix(unlist(my_dum, use.names = FALSE), N, Q_fe) - 1
      orderCluster = matrix(unlist(lapply(my_dum, order), use.names = FALSE), N, Q_fe) - 1

      nbCluster = sapply(my_dum, max)

      fixef_values_tmp = cpp_get_fe_gnl(Q_fe, N, rep(1, N), dumMat, nbCluster, orderCluster)

      # the information on the references
      nb_ref_fe = fixef_values_tmp[[Q_fe+1]]
    } else {
      nb_ref_fe = integer(Q_fe)
    }

    # Slope references (if associated FE + constant)

    names(slope_flag) = fixef_vars

    Q_slope = sum(abs(slope_flag))
    nb_ref_slope = integer(Q_slope)
    i_noVS = 1
    for(i in seq_along(fixef_terms)){

      ft = fixef_terms[i]

      if(!grepl("[[", ft, fixed = TRUE)){
        # No slope => already computed
        nb_ref[i] = nb_ref_fe[i_noVS]
        i_noVS = i_noVS + 1

      } else {
        # Slope
        fe_name = gsub("\\[.+", "", ft)
        my_dum = fe_id_list[[fe_name]]

        my_order = order(my_dum)
        var_sorted = slope_variables[[gsub(".+\\[|\\]+", "", ft)]][my_order]

        # if no associated FE => we check only 0 values
        if(slope_flag[fe_name] < 0){
          nb_ref[i] = cpp_constant_dum(fixef_sizes[fe_name], var_sorted, my_dum[my_order], only_0 = TRUE)
        } else {
          nb_ref[i] = cpp_constant_dum(fixef_sizes[fe_name], var_sorted, my_dum[my_order])
        }
      }
    }

    # we recreate that to avoid conditioning on isSlope later
    fixef_id = fixef_id[fe_all]
    fixef_names = fixef_terms

  } else if(Q == 1){
    # This is the simplest case
    id = id_dummies_vect[[1]]

    myOrder = order(id)
    myDiff = c(1, diff(id[myOrder]))

    select = myOrder[myDiff == 1]

    fixef_values = list(S[select])

    # There are no references => no need to set nb_ref
  } else {
    # We apply a Rcpp script to handle complicated cases (and we don't know beforehand if the input is one)

    dumMat = matrix(unlist(id_dummies_vect), N, Q) - 1
    orderCluster = matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

    nbCluster = sapply(fixef_id, max)

    fixef_values = cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

    # The algorithm is fast but may fail on some instance. We need to check
    if(any(fixef_values[[Q + 1]] > 1) && Q >= 3){
      # we re-compute the "sumFE"
      sum_FE = fixef_values[[1]][1 + dumMat[, 1]]
      for(i in 2:Q){
        sum_FE = sum_FE + fixef_values[[i]][1 + dumMat[, i]]
      }

      if(max(abs(sum_FE - S)) > 1e-1){
        # divergence => we need to correct
        # we recompute the FEs

        is_ref_approx = TRUE

        fixef_sizes = as.integer(object$fixef_sizes)
        new_order = order(object$fixef_sizes, decreasing = TRUE)
        fixef_sizes = fixef_sizes[new_order]

        fe_id_list = object$fixef_id[new_order]
        table_id_I = as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

        slope_flag = rep(0L, Q)
        slope_variables = list()

        S_demean = cpp_demean(y = S, X_raw = 0, r_weights = 0, iterMax = as.integer(fixef.iter),
                    diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
                    fe_id_list = fe_id_list, table_id_I = table_id_I,
                    slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
                    r_init = 0, nthreads = nthreads, save_fixef = TRUE)

        fixef_coef = S_demean$fixef_coef

        end = cumsum(fixef_sizes)
        start = c(1, end + 1)
        for(i in 1:Q){
          fixef_values[[new_order[i]]] = fixef_coef[start[i]:end[i]]
        }
      }
    }

    # the information on the references
    nb_ref = fixef_values[[Q + 1]]
    fixef_values[[Q + 1]] = NULL
  }

  # now saving & adding information
  all_clust = list()
  Q_all = ifelse(isSlope, length(fixef_terms), Q)
  for(i in 1:Q_all){
    # We put it in the right order, if requested
    fn = attr(fixef_id[[i]], "fixef_names")

    if(sorted){
      if(all(!grepl("[^[:digit:]]", fn))) fn = as.numeric(fn)
      my_order = order(fn)

      cv = fixef_values[[i]][my_order]
      names(cv) = fn[my_order]
      all_clust[[fixef_names[i]]] = cv
    } else {
      cv = fixef_values[[i]]
      names(cv) = fn
      all_clust[[fixef_names[i]]] = cv
    }

  }

  class(all_clust) = c("fixest.fixef", "list")

  # Dealing with the references
  if(Q_all > 1){
    names(nb_ref) = fixef_names
    attr(all_clust, "references") = nb_ref

    if(!isSlope) slope_flag = rep(FALSE, Q)

    # warning if unbalanced
    if(notes && sum(nb_ref[!slope_flag]) > Q-1){
      if(is_ref_approx){
        message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted. The number of references is only approximate.")
      } else {
        message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted.")
      }

    }
  }

  # Family information
  attr(all_clust, "exponential") = FALSE
  if(object$method_type == "feNmlm" && object$family %in% c("poisson", "negbin")){
    attr(all_clust, "exponential") = TRUE
  } else if(object$method_type == "feglm" && object$family$link == "log"){
    attr(all_clust, "exponential") = TRUE
  }

  return(all_clust)
}

#' Functions exported from \pkg{nlme} to implement \pkg{fixest} methods
#'
#' The package \pkg{fixest} uses the `fixef` method from \pkg{nlme}. Unfortunately, re-exporting this method is required in order not to attach package \pkg{nlme}.
#'
#' * Here is the help from package \pkg{nlme}: [`fixef`][nlme::fixed.effects]. The help from package \pkg{fixest} is here: [`fixef.fixest`].
#'
#' @note
#' I could find this workaround thanks to the package \pkg{plm}.
#'
#' @name fixef_reexported
#' @keywords internal
NULL

#' @rdname fixef_reexported
#' @name fixef
NULL




#' Displaying the most notable fixed-effects
#'
#' This function plots the 5 fixed-effects with the highest and lowest values, for each of the fixed-effect dimension. It takes as an argument the fixed-effects obtained from the function [`fixef.fixest`] after an estimation using [`femlm`], [`feols`] or [`feglm`].
#'
#' @method plot fixest.fixef
#'
#' @param x An object obtained from the function [`fixef.fixest`].
#' @param n The number of fixed-effects to be drawn. Defaults to 5.
#' @param ... Not currently used.
#'
#' Note that the fixed-effect coefficients might NOT be interpretable. This function is useful only for fully regular panels.
#'
#' If the data are not regular in the fixed-effect coefficients, this means that several \sQuote{reference points} are set to obtain the fixed-effects, thereby impeding their interpretation. In this case a warning is raised.
#'
#' @seealso
#' [`fixef.fixest`] to extract clouster coefficients. See also the main estimation function [`femlm`], [`feols`] or [`feglm`]. Use [`summary.fixest`] to see the results with the appropriate standard-errors, the function [`etable`] to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' # => we account for 3 fixed-effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # obtaining the fixed-effects coefficients
#' fe_trade = fixef(est_pois)
#'
#' # plotting them
#' plot(fe_trade)
#'
#'
plot.fixest.fixef = function(x, n = 5, ...){

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = "n")
  }

  Q = length(x)

  mfrow = as.character(c(11, 12, 22, 22, 32, 32, 33, 33))

  fixef_names = names(x)
  slope_flag = grepl("\\[", fixef_names)

  if(Q > 1 && sum(attr(x, "references")[!slope_flag]) > sum(!slope_flag)-1){
    warning("The fixed-effects are not regular, they cannot be straightforwardly interpreted.", call. = FALSE)
  }

  # modification par:
  opar = par(no.readonly =TRUE)
  on.exit(par(opar))

  par(mfrow = as.numeric(strsplit(mfrow[Q], "")[[1]]), mar = c(3, 3, 2.5, 3))

  addExp = attr(x, "exponential")
  for(i in 1:Q){
    plot_single_cluster(x[[i]], n = n, addExp = addExp, fe_name = fixef_names[i])
  }

}



####
#### fixest own methods ####
####



#' Extracts the coefficients table from an estimation
#'
#' Methods to extracts the coefficients table and its sub-components from an estimation.
#'
#' @param object An estimation (fitted model object), e.g. a `fixest` object.
#' @param ... Other arguments to the methods.
#'
#' @return
#' Returns a matrix (`coeftable`) or vectors.
#'
#' @seealso
#' Please look at the [`coeftable.fixest`] page for more detailed information.
#'
#' @examples
#'
#' est = lm(mpg ~ cyl, mtcars)
#' coeftable(est)
#'
coeftable = function(object, ...){
  UseMethod("coeftable")
}

#' @rdname coeftable
se = function(object, ...){
  UseMethod("se")
}

#' @rdname coeftable
pvalue = function(object, ...){
  UseMethod("pvalue")
}

#' @rdname coeftable
tstat = function(object, ...){
  UseMethod("tstat")
}


#' Extracts the coefficients table from an estimation
#'
#' Default method to extracts the coefficients table and its sub-components from an estimation.
#'
#' @inheritParams etable
#'
#' @method coeftable default
#'
#' @param object The result of an estimation (a fitted model object). Note that this function is made to work with `fixest` objects so it may not work for the specific model you provide.
#' @param ... Other arguments that will be passed to `summary`.
#'
#' First the method summary is applied if needed, then the coefficients table is extracted from its output.
#'
#' The default method is very naive and hopes that the resulting coefficients table contained in the summary of the fitted model is well formed: this assumption is very often wrong. Anyway, there is no development intended since the coeftable/se/pvalue/tstat series of methods is only intended to work well with `fixest` objects. To extract the coefficients table from fitted models in a general way, it's better to use [tidy from broom](https://broom.tidymodels.org/).
#'
#' @return
#' Returns a matrix (`coeftable`) or vectors.
#'
#' @examples
#'
#' # NOTA: This function is really made to handle fixest objects
#' # The default methods works for simple structures, but you'd be
#' # likely better off with broom::tidy for other models
#'
#' est = lm(mpg ~ cyl, mtcars)
#' coeftable(est)
#'
#' se(est)
#'
#'
#'
#'
#'
coeftable.default = function(object, keep, drop, order, ...){
  # This function is EXTREMELY naive and I don't intend to improve it
  # => there is tidy for that which is much better
  # I just created that method to handle fixest/fixest_multi more easily

  check_arg(keep, drop, order, "NULL character vector no na")

  if(!any(grepl("summary", class(object)))){
    object_sum = try(summary(object, ...))
    if(!inherits(object_sum, "try-error")){
      object = object_sum
    }
  }

  list_mat = object[sapply(object, is.matrix)]

  ok = FALSE
  for(i in seq_along(list_mat)){
    mat = list_mat[[i]]
    if(!is.null(colnames(mat)) && any(grepl("(?i)(coef|estimate|value|Pr\\()", colnames(mat)))){
      ok = TRUE
      res = mat
    }
  }

  if(ok == FALSE){
    stop("Sorry, the coeffficients table could not be extracted.")
  }

  if(!missnull(keep) || !missnull(drop) || !missnull(order)){
    r_names = rownames(res)
    r_names = keep_apply(r_names, keep)
    r_names = drop_apply(r_names, drop)
    r_names = order_apply(r_names, order)

    if(length(r_names) == 0){
      return(NULL)
    }

    res = res[r_names, , drop = FALSE]
  }

  res
}


#' @describeIn coeftable.default Extracts the standard-errors from an estimation
se.default = function(object, keep, drop, order, ...){
  # There is NO GARANTEE that it works

  check_arg(keep, drop, order, "NULL character vector no na")

  mat = coeftable(object, keep = keep, drop = drop, order = order, ...)


  if(is.null(mat)){
    return(NULL)
  }

  res = mat[, 2]
  names(res) = row.names(mat)

  res
}

#' @describeIn coeftable.default Extracts the standard-errors from an estimation
tstat.default = function(object, keep, drop, order, ...){
  # There is NO GARANTEE that it works

  check_arg(keep, drop, order, "NULL character vector no na")

  mat = coeftable(object, keep = keep, drop = drop, order = order, ...)

  if(is.null(mat)){
    return(NULL)
  }

  if(ncol(mat) != 4){
    stop("The p-values could not be obtained. Use broom::tidy?")
  }

  res = mat[, 3]
  names(res) = row.names(mat)

  res
}

#' @describeIn coeftable.default Extracts the p-values from an estimation
pvalue.default = function(object, keep, drop, order, ...){
  # There is NO GARANTEE that it works

  check_arg(keep, drop, order, "NULL character vector no na")

  mat = coeftable(object, keep = keep, drop = drop, order = order, ...)

  if(is.null(mat)){
    return(NULL)
  }

  if(ncol(mat) != 4){
    stop("The p-values could not be obtained. Use broom::tidy?")
  }

  res = mat[, 4]
  names(res) = row.names(mat)

  res
}

#' @describeIn coeftable.default Extracts the standard-errors from a VCOV matrix
se.matrix = function(object, keep, drop, order, ...){
  # There is NO GARANTEE that it works

  check_arg(object, "square numeric matrix")
  check_arg(keep, drop, order, "NULL character vector no na")

  vcov_diag = diag(object)
  vcov_diag[vcov_diag < 0] = NA

  all_names = rownames(object)
  se = sqrt(vcov_diag)
  names(se) = all_names

  if(!missnull(keep) || !missnull(drop) || !missnull(order)){

    if(is.null(all_names)){
      stop("The VCOV matrix has no names attributes: the arguments keep, drop or order cannot be used.")
    }

    all_names = keep_apply(all_names, keep)
    all_names = drop_apply(all_names, drop)
    all_names = order_apply(all_names, order)

    if(length(all_names) == 0){
      return(NULL)
    }

    se = se[all_names]
  }

  se
}



#' Obtain various statistics from an estimation
#'
#' Set of functions to directly extract some commonly used statistics, like the p-value or the table of coefficients, from estimations. This was first implemented for `fixest` estimations, but has some support for other models.
#'
#' @inheritParams etable
#'
#' @method coeftable fixest
#'
#' @param object A `fixest` object. For example an estimation obtained from [`feols`].
#' @param cluster Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over `var1` and `var2` contained in the data.frame `base` used for the estimation. All the following `cluster` arguments are valid and do the same thing: `cluster = base[, c("var1, "var2")]`, `cluster = c("var1, "var2")`, `cluster = ~var1+var2`. If the two variables were used as clusters in the estimation, you could further use `cluster = 1:2` or leave it blank with `se = "twoway"` (assuming `var1` \[resp. `var2`\] was the 1st \[resp. 2nd\] cluster).
#' @param list Logical, default is `FALSE`. If `TRUE`, then a nested list is returned, the first layer is accessed with the coefficients names; the second layer with the following values: `coef`, `se`, `tstat`, `pvalue`. Note that the variable `"(Intercept)"` is renamed into `"constant"`.
#' @param ... Other arguments to be passed to [`summary.fixest`].
#'
#' @details
#' This set of tiny functions is primarily constructed for `fixest` estimations.
#'
#' @return
#' Returns a table of coefficients, with in rows the variables and four columns: the estimate, the standard-error, the t-statistic and the p-value.
#'
#' If `list = TRUE` then a nested list is returned, the first layer is accessed with the coefficients names; the second layer with the following values: `coef`, `se`, `tstat`, `pvalue`. For example, with `res = coeftable(est, list = TRUE)` you can access the SE of the coefficient `x1` with `res$x1$se`; and its coefficient with `res$x1$coef`, etc.
#'
#' @examples
#'
#' # Some data and estimation
#' data(trade)
#' est = fepois(Euros ~ log(dist_km) | Origin^Product + Year, trade)
#'
#' #
#' # Coeftable/se/tstat/pvalue
#' #
#'
#' # Default is clustering along Origin^Product
#' coeftable(est)
#' se(est)
#' tstat(est)
#' pvalue(est)
#'
#' # Now with two-way clustered standard-errors
#' #  and using coeftable()
#'
#' coeftable(est, cluster = ~Origin + Product)
#' se(est, cluster = ~Origin + Product)
#' pvalue(est, cluster = ~Origin + Product)
#' tstat(est, cluster = ~Origin + Product)
#'
#' # Or you can cluster only once:
#' est_sum = summary(est, cluster = ~Origin + Product)
#' coeftable(est_sum)
#' se(est_sum)
#' tstat(est_sum)
#' pvalue(est_sum)
#'
#' # You can use the arguments keep, drop, order
#' # to rearrange the results
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' est_iv = feols(y ~ x1 | x2 ~ x3, base)
#'
#' tstat(est_iv, keep = "x1")
#' coeftable(est_iv, keep = "x1|Int")
#'
#' coeftable(est_iv, order = "!Int")
#'
#' #
#' # Using lists
#' #
#'
#' # Returning the coefficients table as a list can be useful for quick
#' # reference in markdown documents.
#' # Note that the "(Intercept)" is renamed into "constant"
#'
#' res = coeftable(est_iv, list = TRUE)
#'
#' # coefficient of the constant:
#' res$constant$coef
#'
#' # pvalue of x1
#' res$x1$pvalue
#'
#'
#'
coeftable.fixest = function(object, vcov = NULL, ssc = NULL, cluster = NULL,
              keep = NULL, drop = NULL, order = NULL, list = FALSE, ...){
  # We don't explicitly refer to the other arguments

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(list, "logical scalar")


  if(!isTRUE(object$summary) || !(all_missing(vcov, ssc, cluster) && ...length() > 0)){
    object = summary(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)
  }

  # Let's find out the coefficients table
  res = object$coeftable

  if(!missnull(keep) || !missnull(drop) || !missnull(order)){
    r_names = rownames(res)
    r_names = keep_apply(r_names, keep)
    r_names = drop_apply(r_names, drop)
    r_names = order_apply(r_names, order)

    if(length(r_names) == 0){
      return(NULL)
    }

    res = res[r_names, , drop = FALSE]
  }

  if(list){
    res_list = list()
    coef_names = rownames(res)
    coef_names[coef_names == "(Intercept)"] = "constant"
    for(i in 1:nrow(res)){
      row = res[i, ]
      item = list(coef = row[1], se = row[2], tstat = row[3], pvalue = row[4])
      res_list[[coef_names[i]]] = item
    }
    res = res_list
  }

  res
}

#' @describeIn coeftable.fixest Extracts the standard-error of an estimation
se.fixest = function(object, vcov = NULL, ssc = NULL, cluster = NULL,
           keep = NULL, drop = NULL, order = NULL, ...){

  check_arg(keep, drop, order, "NULL character vector no na")

  mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster,
          keep = keep, drop = drop, order = order, ...)

  if(is.null(mat)){
    return(NULL)
  }

  res = mat[, 2]
  names(res) = row.names(mat)

  res
}

#' @describeIn coeftable.fixest Extracts the t-statistics of an estimation
tstat.fixest = function(object, vcov = NULL, ssc = NULL, cluster = NULL,
            keep = NULL, drop = NULL, order = NULL, ...){

  check_arg(keep, drop, order, "NULL character vector no na")

  mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster,
          keep = keep, drop = drop, order = order, ...)

  if(is.null(mat)){
    return(NULL)
  }

  res = mat[, 3]
  names(res) = row.names(mat)

  res
}

#' @describeIn coeftable.fixest Extracts the p-value of an estimation
pvalue.fixest = function(object, vcov = NULL, ssc = NULL, cluster = NULL,
             keep = NULL, drop = NULL, order = NULL, ...){

  check_arg(keep, drop, order, "NULL character vector no na")

  mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster,
          keep = keep, drop = drop, order = order, ...)

  if(is.null(mat)){
    return(NULL)
  }

  res = mat[, 4]
  names(res) = row.names(mat)

  res
}


####
#### stats ####
####



#' Extracts the number of observations form a `fixest` object
#'
#' This function simply extracts the number of observations form a `fixest` object, obtained using the functions [`femlm`], [`feols`] or [`feglm`].
#'
#' @inheritParams summary.fixest
#'
#' @param ... Not currently used.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. Use [`summary.fixest`] to see the results with the appropriate standard-errors, [`fixef.fixest`] to extract the fixed-effects coefficients, and the function [`etable`] to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @return
#' It returns an interger.
#'
#' @examples
#'
#' # simple estimation on iris data with "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' nobs(res)
#' logLik(res)
#'
#'
nobs.fixest = function(object, ...){
  object$nobs
}

#' Aikake's an information criterion
#'
#' This function computes the AIC (Aikake's, an information criterion) from a `fixest` estimation.
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Optionally, more fitted objects.
#' @param k A numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC (i.e. `AIC=-2*LL+k*nparams`).
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2\times LogLikelihood + k\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this criterion on [`AIC`][stats::AIC].
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. Other statictics methods: [`BIC.fixest`], [`logLik.fixest`], [`nobs.fixest`].
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'              Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' AIC(res1, res2)
#' BIC(res1, res2)
#'
#'
AIC.fixest = function(object, ..., k = 2){

  dots = list(...)
  if(length(dots) > 0){
    # we check consistency with observations
    nobs_all = c(nobs(object), sapply(dots, nobs))

    if(any(diff(nobs_all) != 0)){
      warning("Models are not all fitted to the same number of observations.")
    }

    otherAIC = sapply(dots, AIC)
  } else {
    otherAIC = c()
  }

  all_AIC = c(-2*logLik(object) + k*object$nparams, otherAIC)

  all_AIC
}

#' Bayesian information criterion
#'
#' This function computes the BIC (Bayesian information criterion) from a `fixest` estimation.
#'
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Optionally, more fitted objects.
#'
#' @details
#' The BIC is computed as follows:
#' \deqn{BIC = -2\times LogLikelihood + \log(nobs)\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this criterion on [`AIC`][stats::AIC].
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. Other statistics functions: [`AIC.fixest`], [`logLik.fixest`].
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' AIC(res1, res2)
#' BIC(res1, res2)
#'
BIC.fixest = function(object, ...){

  dots = list(...)
  if(length(dots) > 0){
    # we check consistency with observations
    nobs_all = c(nobs(object), sapply(dots, nobs))

    if(any(diff(nobs_all) != 0)){
      warning("Models are not all fitted to the same number of observations.")
    }

    otherBIC = sapply(dots, BIC)
  } else {
    otherBIC = c()
  }

  all_BIC = c(-2*logLik(object) + object$nparams*log(nobs(object)), otherBIC)

  all_BIC
}

#' Extracts the log-likelihood
#'
#' This function extracts the log-likelihood from a `fixest` estimation.
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Not currently used.
#'
#' @details
#' This function extracts the log-likelihood based on the model fit. You can have more information on the likelihoods in the details of the function [`femlm`].
#'
#' @return
#' It returns a numeric scalar.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. Other statistics functions: [`AIC.fixest`], [`BIC.fixest`].
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data with "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' nobs(res)
#' logLik(res)
#'
#'
logLik.fixest = function(object, ...){

  if(object$method_type == "feols"){
    # if the summary is 'lean', then no way we can compute that
    resid = object$residuals
    if(is.null(resid)) resid = NA

    sigma = sqrt(mean(resid^2))
    n = length(resid)
    ll = -1/2/sigma^2 * sum(resid^2) - n * log(sigma) - n * log(2*pi)/2
  } else {
    ll = object$loglik
  }

  ll
}

#' Extracts the coefficients from a `fixest` estimation
#'
#' This function extracts the coefficients obtained from a model estimated with [`femlm`], [`feols`] or [`feglm`].
#'
#' @inheritParams nobs.fixest
#' @inheritParams etable
#'
#' @param agg Logical scalar, default is `TRUE`. If the coefficients of the estimation have been aggregated, whether to report the aggregated coefficients. If `FALSE`, the raw coefficients will be returned.
#' @param collin Logical, default is `FALSE`. Whether the coefficients removed because of collinearity should be also returned as `NA`. It cannot be used when coefficients aggregation is also used.
#' @param ... Not currently used.
#'
#' @details
#' The coefficients are the ones that have been found to maximize the log-likelihood of the specified model. More information can be found on the models from the estimations help pages: [`femlm`], [`feols`] or [`feglm`].
#'
#' Note that if the model has been estimated with fixed-effects, to obtain the fixed-effect coefficients, you need to use the function [`fixef.fixest`].
#'
#' @return
#' This function returns a named numeric vector.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`summary.fixest`], [`confint.fixest`], [`vcov.fixest`], [`etable`], [`fixef.fixest`].
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # the coefficients of the variables:
#' coef(res)
#'
#' # the fixed-effects coefficients:
#' fixef(res)
#'
#'
coef.fixest = coefficients.fixest = function(object, keep, drop, order,
                       collin = FALSE, agg = TRUE, ...){

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(collin, agg, "logical scalar")

  if(isTRUE(object$is_agg) && agg){
    if(collin){
      warning("The argument 'collin = TRUE' cannot be used when there is coefficient aggregation.")
    }

    res = object$coeftable[, 1]
    names(res) = rownames(object$coeftable)

  } else if(collin && !is.null(object$collin.coef)){
    res = object$collin.coef
  } else {
    res = object$coefficients
  }

  if(!missnull(keep) || !missnull(drop) || !missnull(order)){
    cnames = names(res)
    cnames = keep_apply(cnames, keep)
    cnames = drop_apply(cnames, drop)
    cnames = order_apply(cnames, order)

    if(length(cnames) == 0){
      return(numeric(0))
    }

    res = res[cnames]
  }

  if(identical(object$family, "negbin")){
    res = res[-length(res)]
  }


  # deltaMethod tweak
  if(is_calling_fun("deltaMethod")){
    sysOrigin = sys.parent()
    mc_DM = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin))

    if("parameterNames" %in% names(mc_DM)){
      PN = eval(mc_DM$parameterNames, parent.frame())

      check_value(PN, "character vector no na len(data)", .data = res,
            .arg_name = "parameterNames", .up = 1)

      names(res) = PN
    }
  }

  res
}

#' @rdname coef.fixest
coefficients.fixest <- coef.fixest


#' Extracts fitted values from a `fixest` fit
#'
#' This function extracts the fitted values from a model estimated with [`femlm`], [`feols`] or [`feglm`]. The fitted values that are returned are the *expected predictor*.
#'
#' @inheritParams nobs.fixest
#'
#' @param type Character either equal to `"response"` (default) or `"link"`. If `type="response"`, then the output is at the level of the response variable, i.e. it is the expected predictor \eqn{E(Y|X)}. If `"link"`, then the output is at the level of the explanatory variables, i.e. the linear predictor \eqn{X\cdot \beta}.
#' @param na.rm Logical, default is `TRUE`. If `FALSE` the number of observation returned will be the number of observations in the original data set, otherwise it will be the number of observations used in the estimation.
#' @param ... Not currently used.
#'
#' @details
#' This function returns the *expected predictor* of a `fixest` fit. The likelihood functions are detailed in [`femlm`] help page.
#'
#' @return
#' It returns a numeric vector of length the number of observations used to estimate the model.
#'
#' If `type = "response"`, the value returned is the expected predictor, i.e. the expected value of the dependent variable for the fitted model: \eqn{E(Y|X)}.
#' If `type = "link"`, the value returned is the linear predictor of the fitted model, that is \eqn{X\cdot \beta} (remind that \eqn{E(Y|X) = f(X\cdot \beta)}).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`resid.fixest`], [`predict.fixest`], [`summary.fixest`], [`vcov.fixest`], [`fixef.fixest`].
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we extract the fitted values
#' y_fitted_poisson = fitted(res_poisson)
#'
#' # Same estimation but in OLS (Gaussian family)
#' res_gaussian = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris, family = "gaussian")
#'
#' y_fitted_gaussian = fitted(res_gaussian)
#'
#' # comparison of the fit for the two families
#' plot(iris$Sepal.Length, y_fitted_poisson)
#' points(iris$Sepal.Length, y_fitted_gaussian, col = 2, pch = 2)
#'
#'
fitted.fixest = fitted.values.fixest = function(object, type = c("response", "link"), na.rm = TRUE, ...){

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = "type")
  }

  type = match.arg(type)

  fit = predict(object)

  if(type == "response" || object$method_type == "feols"){
    res = fit
  } else if(!is.null(object$mu)){
    res = object$mu
  } else if(object$method_type == "feNmlm"){
    family = object$family
    famFuns = switch(family,
             poisson = ml_poisson(),
             negbin = ml_negbin(),
             logit = ml_logit(),
             gaussian = ml_gaussian())

    res = famFuns$linearFromExpected(fit)
  } else {
    res = object$family$linkfun(fit)
  }

  # Nota: obs can be removed: either because of NA, either because perfect fit
  # Shall I put perfect fit as NA since they're out of the estimation???
  # Still pondering...
  # Actually adding them means a lot of work to ensure consistency (also in predict...)
  if(!na.rm) res = fill_with_na(res, object)

  res
}

#' @rdname fitted.fixest
#' @method fitted.values fixest
fitted.values.fixest <- fitted.fixest

#' Extracts residuals from a `fixest` object
#'
#' This function extracts residuals from a fitted model estimated with [`femlm`], [`feols`] or [`feglm`].
#'
#' @inheritParams nobs.fixest
#'
#' @param type A character scalar, either `"response"` (default), `"deviance"`, `"pearson"`, or `"working"`. Note that the `"working"` corresponds to the residuals from the weighted least square and only applies to [`feglm`] models.
#' @param na.rm Logical, default is `TRUE`. Whether to remove the observations with NAs from the original data set. If `FALSE`, then the vector returned is always of the same length as the original data set.
#' @param ... Not currently used.
#'
#'
#' @return
#' It returns a numeric vector of the length the number of observations used for the estimation (if `na.rm = TRUE`) or of the length of the original data set (if `na.rm = FALSE`).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`fitted.fixest`], [`predict.fixest`], [`summary.fixest`], [`vcov.fixest`], [`fixef.fixest`].
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we plot the residuals
#' plot(resid(res_poisson))
#'
resid.fixest = residuals.fixest = function(object, type = c("response", "deviance", "pearson", "working"),
                       na.rm = TRUE, ...){

  check_set_arg(type, "match")
  check_set_arg(na.rm, "logical scalar")

  method = object$method
  family = object$family

  r = object$residuals
  w = object[["weights"]]

  if(isTRUE(object$lean)){
    stop("The method 'resid.fixest' cannot be applied to a 'lean' fixest object. Please apply reestimate with 'lean = FALSE'.")
  }

  if(method %in% c("feols", "feols.fit") || (method %in% c("feNmlm", "femlm") && family == "gaussian")){

    if(type == "working") stop("Type 'working' only applies to models fitted via feglm (thus is not valid for feols).")

    if(type %in% c("deviance", "pearson") && !is.null(w)){
      res = r * sqrt(w)
    } else {
      res = r
    }

  } else if(method %in% c("fepois", "feglm")){

    if(type == "response"){
      res = r

    } else if(type == "working"){
      res = object$working_residuals

    } else {
      mu = object$fitted.values
      if(is.null(w)) w = rep(1, length(r))

      if(type == "deviance"){
        y = r + mu

        res = sqrt(pmax((object$family$dev.resids)(y, mu, w), 0))
        qui = y < mu
        res[qui] = -res[qui]

      } else if(type == "pearson"){
        res = r * sqrt(w)/sqrt(object$family$variance(object$fitted.values))

      }
    }


  } else {

    if(type == "working") stop("Type 'working' only applies to models fitted via feglm (thus is not valid for ", method, ").")

    if(type == "response"){
      res = r

    } else {
      # deviance or pearson
      mu = object$fitted.values
      if(is.null(w)) w = rep(1, length(r))

      theta = ifelse(family == "negbin", object$theta, 1)

      if(type == "deviance"){

        # dev.resids function
        if(family == "poisson"){
          dev.resids = poisson()$dev.resids

        } else if(family == "logit"){
          dev.resids = binomial()$dev.resids

        } else if(family == "negbin"){
          dev.resids = function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))

        }

        y = object$residuals + mu

        res = sqrt(pmax(dev.resids(y, mu, w), 0))
        qui = y < mu
        res[qui] = -res[qui]

      } else if(type == "pearson"){

        # variance function
        if(family == "poisson"){
          variance = poisson()$variance

        } else if(family == "logit"){
          variance = binomial()$variance

        } else if(family == "negbin"){
          variance = function(mu) mu + mu^2/theta

        }

        res = r * sqrt(w)/sqrt(variance(mu))

      }
    }


  }

  if(!na.rm){
    res = fill_with_na(res, object)
  }

  res
}

#' @rdname resid.fixest
residuals.fixest <- resid.fixest

#' Predict method for `fixest` fits
#'
#' This function obtains prediction from a fitted model estimated with [`femlm`], [`feols`] or [`feglm`].
#'
#' @inheritParams nobs.fixest
#' @inheritParams fitted.fixest
#' @inheritParams summary.fixest
#'
#' @param newdata A data.frame containing the variables used to make the prediction. If not provided, the fitted expected (or linear if `type = "link"`) predictors are returned.
#' @param sample Either "estimation" (default) or "original". This argument is only used when arg. 'newdata' is missing, and is ignored otherwise. If equal to "estimation", the vector returned matches the sample used for the estimation. If equal to "original", it matches the original data set (the observations not used for the estimation being filled with NAs).
#' @param se.fit Logical, default is `FALSE`. If `TRUE`, the standard-error of the predicted value is computed and returned in a column named `se.fit`. This feature is only available for OLS models not containing fixed-effects.
#' @param interval Either "none" (default), "confidence" or "prediction". What type of confidence interval to compute. Note that this feature is only available for OLS models not containing fixed-effects (GLM/ML models are not covered).
#' @param level A numeric scalar in between 0.5 and 1, defaults to 0.95. Only used when the argument 'interval' is requested, it corresponds to the width of the confidence interval.
#' @param fixef Logical scalar, default is `FALSE`. If `TRUE`, a data.frame is returned, with each column representing the fixed-effects coefficients for each observation in `newdata` -- with as many columns as fixed-effects. Note that when there are variables with varying slopes, the slope coefficients are returned (i.e. they are not multiplied by the variable).
#' @param vs.coef Logical scalar, default is `FALSE`. Only used when `fixef = TRUE` and when variables with varying slopes are present. If `TRUE`, the coefficients of the variables with varying slopes are returned instead of the coefficient multiplied by the value of the variables (default).
#' @param ... Not currently used.
#'
#'
#' @return
#' It returns a numeric vector of length equal to the number of observations in argument `newdata`.
#' If `newdata` is missing, it returns a vector of the same length as the estimation sample, except if `sample = "original"`, in which case the length of the vector will match the one of the original data set (which can, but also cannot, be the estimation sample).
#' If `fixef = TRUE`, a `data.frame` is returned.
#' If `se.fit = TRUE` or `interval != "none"`, the object returned is a data.frame with the following columns: `fit`, `se.fit`, and, if CIs are requested, `ci_low` and `ci_high`.
#'
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`update.fixest`], [`summary.fixest`], [`vcov.fixest`], [`fixef.fixest`].
#'
#' @examples
#'
#' # Estimation on iris data
#' res = fepois(Sepal.Length ~ Petal.Length | Species, iris)
#'
#' # what would be the prediction if the data was all setosa?
#' newdata = data.frame(Petal.Length = iris$Petal.Length, Species = "setosa")
#' pred_setosa = predict(res, newdata = newdata)
#'
#' # Let's look at it graphically
#' plot(c(1, 7), c(3, 11), type = "n", xlab = "Petal.Length",
#'      ylab = "Sepal.Length")
#'
#' newdata = iris[order(iris$Petal.Length), ]
#' newdata$Species = "setosa"
#' lines(newdata$Petal.Length, predict(res, newdata))
#'
#' # versicolor
#' newdata$Species = "versicolor"
#' lines(newdata$Petal.Length, predict(res, newdata), col=2)
#'
#' # virginica
#' newdata$Species = "virginica"
#' lines(newdata$Petal.Length, predict(res, newdata), col=3)
#'
#' # The original data
#' points(iris$Petal.Length, iris$Sepal.Length, col = iris$Species, pch = 18)
#' legend("topleft", lty = 1, col = 1:3, legend = levels(iris$Species))
#'
#'
#' #
#' # Getting the fixed-effect coefficients for each obs.
#' #
#'
#' data(trade)
#' est_trade = fepois(Euros ~ log(dist_km) | Destination^Product +
#'                                            Origin^Product + Year, trade)
#' obs_fe = predict(est_trade, fixef = TRUE)
#' head(obs_fe)
#'
#' # can we check we get the right sum of fixed-effects
#' head(cbind(rowSums(obs_fe), est_trade$sumFE))
#'
#'
#' #
#' # Standard-error of the prediction
#' #
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#'
#' est = feols(y ~ x1 + species, base)
#'
#' head(predict(est, se.fit = TRUE))
#'
#' # regular confidence interval
#' head(predict(est, interval = "conf"))
#'
#' # adding the residual to the CI
#' head(predict(est, interval = "predi"))
#'
#' # You can change the type of SE on the fly
#' head(predict(est, interval = "conf", vcov = ~species))
#'
#'
#'
predict.fixest = function(object, newdata, type = c("response", "link"), se.fit = FALSE,
              interval = "none", level = 0.95, fixef = FALSE,
              vs.coef = FALSE, sample = c("estimation", "original"),
              vcov = NULL, ssc = NULL, ...){

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = c("newdata", "type"))
  }

  # Controls
  check_set_arg(type, sample, "match")
  check_arg(fixef, vs.coef, "logical scalar")
  check_arg(se.fit, "logical scalar")
  check_arg(level, "numeric scalar GT{.50} LT{1}")
  check_set_arg(interval, "match(none, confidence, prediction)")
  if(!se.fit && interval != "none"){
    se.fit = TRUE
  }

  if(se.fit && object$method_type != "feols"){
    stop("The standard-error of the prediction is currently only available for OLS models, sorry.")
  }

  # renaming to clarify
  fixef.return = fixef
  do_anyway = fixef.return || se.fit

  # if newdata is missing
  is_original_data = FALSE
  if(missing(newdata)){

    if(do_anyway || isTRUE(object$lean)){
      newdata = fetch_data(object, "In 'predict', ")
      is_original_data = TRUE
    } else {
      if(type == "response" || object$method_type == "feols"){
        res = object$fitted.values
      } else if(object$method_type == "feNmlm") {
        if("mu" %in% names(object)){
          res = object$mu
        } else {
          family = object$family
          famFuns = switch(family,
                   poisson = ml_poisson(),
                   negbin = ml_negbin(),
                   logit = ml_logit(),
                   gaussian = ml_gaussian())

          res = famFuns$linearFromExpected(object$fitted.values)
        }
      } else {
        res = object$family$linkfun(object$fitted.values)
      }

      if(sample == "original") res = fill_with_na(res, object)

      return(res)
    }

  }

  if(!is.matrix(newdata) && !"data.frame" %in% class(newdata)){
    stop("Argument 'newdata' must be a data.frame.")
  }

  # we ensure it really is a clean data.frame
  newdata = as.data.frame(newdata)

  mc = match.call()
  if(fixef.return){

    if(is.null(object$fixef_vars)){
      stop("The argument 'fixef=TRUE' cannot work since the estimation did not contain fixed-effects.")
    }

    if("type" %in% names(mc)){
      warning("Argument 'type' is ignored when fixef = TRUE.")
    }
  }

  # We deconstruct it in four steps:
  # 1) cluster
  # 2) linear
  # 3) non-linear
  # 4) offset

  # + step 0: panel setup

  n = nrow(newdata)

  # NOTA 2019-11-26: I'm pondering whether to include NA-related messages
  # (would it be useful???)


  # STEP 0: panel setup

  fml = object$fml
  panel__meta__info = set_panel_meta_info(object, newdata)

  #
  # 1) Fixed-effects
  #

  # init fixed-effect values
  value_fixef = 0

  fixef_vars = object$fixef_vars
  if(!is.null(fixef_vars)){

    n_fe = length(fixef_vars)

    # Extraction of the FEs
    id_fixef = list()
    for(i in 1:n_fe){
      # checking if the variable is in the newdata
      fe_var = fixef_vars[i]
      variable = all.vars(str2lang(fe_var))
      isNotHere = !variable %in% names(newdata)
      if(any(isNotHere)){
        stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a fixed-effect variable).")
      }

      # The values taken by the FE variable
      fixef_values_possible = attr(object$fixef_id[[i]], "fixef_names")

      # Checking if ^ is present
      if(grepl("\\^", fe_var)){
        # If fastCombine was used => we're screwed, impossible to recover
        if(!all(grepl("_", fixef_values_possible, fixed = TRUE))){
          stop("You cannot use predict() based on the initial regression since the fixed-effect '", fe_var, "' was combined using an algorithm dropping the FE values (but fast). Please re-run the regression using the argument 'combine.quick=FALSE'.")
        }

        fe_var = fml_combine(fe_var, fastCombine = FALSE, vars = TRUE)
      }

      # Obtaining the vector of fixed-effect
      fixef_current = eval(str2lang(fe_var), newdata)

      fixef_current_num = unclass(factor(fixef_current, levels = fixef_values_possible))
      id_fixef[[i]] = fixef_current_num
    }

    names(id_fixef) = fixef_vars

    # Value of the fixef coefficients // we don't show the notes, it's inelegant
    fixef_coef = fixef(object, sorted = FALSE, notes = FALSE)

    # We create the DF to be returned
    if(fixef.return){
      fixef_df = list()
    }

    # Adding the FEs and Slopes
    if(!is.null(object$fixef_terms)){

      terms_full = extract_fe_slope(object$fixef_terms)
      fixef_vars = terms_full$fixef_vars
      slope_fe = terms_full$slope_fe
      slope_vars = terms_full$slope_vars
      slope_terms = terms_full$slope_terms

      # We extract the slope variables
      slope_vars_unik = unique(slope_vars)

      slope_var_list = list()
      for(i in 1:length(slope_vars_unik)){
        variable = all.vars(str2lang(slope_vars_unik[i]))
        isNotHere = !variable %in% names(newdata)
        if(any(isNotHere)){
          stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a variable with varying slope).")
        }

        slope_var_list[[slope_vars_unik[i]]] = eval(str2lang(slope_vars_unik[i]), newdata)
      }

      # Adding the FE values
      for(var in fixef_vars){
        fixef_current_num = id_fixef[[var]]
        fixef_coef_current = fixef_coef[[var]]

        if(fixef.return){
          fixef_df[[var]] = fixef_coef_current[fixef_current_num]

        } else {
          value_fixef = value_fixef + fixef_coef_current[fixef_current_num]
        }

      }

      # Adding the slopes
      for(i in seq_along(slope_vars)){

        fixef_current_num = id_fixef[[slope_fe[i]]]
        fixef_coef_current = fixef_coef[[slope_terms[i]]]

        if(fixef.return){
          vname = slope_terms[i]

          # We return only the coefs OR the coef * the variable
          if(vs.coef){
            fixef_df[[vname]] = fixef_coef_current[fixef_current_num]

          } else {
            fixef_df[[vname]] = fixef_coef_current[fixef_current_num] * slope_var_list[[slope_vars[i]]]
          }


        } else {
          value_fixef = value_fixef + fixef_coef_current[fixef_current_num] * slope_var_list[[slope_vars[i]]]
        }
      }


    } else {
      # Adding only FEs
      for(i in 1:n_fe){
        fixef_current_num = id_fixef[[i]]
        fixef_coef_current = fixef_coef[[i]]

        if(fixef.return){
          fixef_df[[fixef_vars[i]]] = fixef_coef_current[fixef_current_num]

        } else {
          value_fixef = value_fixef + fixef_coef_current[fixef_current_num]
        }

      }
    }

    if(fixef.return){

      # putting the results into a DF
      res = fixef_df
      attr(res, "row.names") = .set_row_names(length(res[[1L]]))
      oldClass(res) = "data.frame"

      if(is_original_data && sample == "estimation"){
        # here we want the same nber of obs
        # as in the estimation sample
        for(i in seq_along(object$obs_selection)){
          res = res[object$obs_selection[[i]], , drop = FALSE]
        }
      }

      return(res)
    }

    # dropping names
    value_fixef = as.vector(value_fixef)
  }

  #
  # 2) Linear values
  #

  coef = object$coefficients

  value_linear = 0
  var_keep = NULL
  rhs_fml = fml_split(fml, 1)
  linear.varnames = all_vars_with_i_prefix(rhs_fml[[3]])

  if(length(linear.varnames) > 0){
    # Checking all variables are there

    if(isTRUE(object$iv) && object$iv_stage == 2){
      names(coef) = gsub("^fit_", "", names(coef))
      linear.varnames = c(linear.varnames, all_vars_with_i_prefix(object$fml_all$iv[[2]]))
      iv_fml = object$fml_all$iv
      rhs_fml = .xpd(..lhs ~ ..endo + ..rhs, ..lhs = rhs_fml[[2]], ..endo = iv_fml[[2]], ..rhs = rhs_fml[[3]])
    }

    varNotHere = setdiff(linear.varnames, names(newdata))
    if(length(varNotHere) > 0){
      stop("The variable", enumerate_items(varNotHere, "s.quote"),
         " used to estimate the model (in fml) ", ifsingle(varNotHere, "is", "are"),
         " missing in the data.frame given by the argument 'newdata'.")
    }

    # we create the matrix
    matrix_linear = error_sender(fixest_model_matrix_extra(object = object, newdata = newdata,
                                 original_data = FALSE, fml = rhs_fml,
                                 i_noref = TRUE),
                   "Error when creating the linear matrix: ")

    # Checking the levels created with i()
    mm_info_new = attr(matrix_linear, "model_matrix_info")
    if(!is.null(mm_info_new)){
      mm_info = object$model_matrix_info
      # The order of creation is exactly the same (same fun used),
      # so the two mm_info are identical in structure
      for(i in seq_along(mm_info)){
        mm_new_i = mm_info_new[[i]]
        mm_i = mm_info[[i]]
        if("coef_names_full" %in% names(mm_i)){
          pblm = setdiff(mm_new_i$coef_names_full, mm_i$coef_names_full)
          if(length(pblm) > 0){
            stop(dsb("In i(), predictions cannot be done for values that were not present at estimation time.",
                 " It concerns the value.[*s_, 3KO, C?pblm]."))
          }
        }
      }
    }

    var_keep = intersect(names(coef), colnames(matrix_linear))
    value_linear = value_linear + as.vector(matrix_linear[, var_keep, drop = FALSE] %*% coef[var_keep])
  }

  #
  # 3) Non linear terms
  #

  value_NL = 0
  NL_fml = object$NL.fml
  if(!is.null(NL_fml)){
    # controlling that we can evaluate that
    NL_vars = all.vars(NL_fml)
    varNotHere = setdiff(NL_vars, c(names(coef), names(newdata)))
    if(length(varNotHere) > 0){
      stop("Some variables used to estimate the model (in the non-linear formula) are missing from argument 'newdata': ", enumerate_items(varNotHere), ".")
    }

    var2send = intersect(NL_vars, names(newdata))
    env = new.env()
    for(var in var2send){
      assign(var, newdata[[var]], env)
    }

    coef2send = setdiff(NL_vars, names(newdata))
    for(iter_coef in coef2send){
      assign(iter_coef, coef[iter_coef], env)
    }

    # Evaluation of the NL part
    value_NL = eval(NL_fml[[2]], env)
  }

  #
  # 4) offset value
  #

  value_offset = 0
  offset = object$call$offset
  if(!is.null(offset)){
    # evaluation of the offset

    if(is.numeric(offset)){
      # simple numeric offset
      value_offset = offset

    } else {
      # offset valid only if formula
      offset_char = as.character(offset)

      if(length(offset_char) == 2 && offset_char[1] == "~"){
        offset_fml = eval(offset)
        varNotHere = setdiff(all.vars(offset_fml), names(newdata))
        if(length(varNotHere) > 0){
          stop("In the offset, the variable", enumerate_items(varNotHere, "s.is"), " not present in 'newdata'.")
        }

        value_offset = eval(offset_fml[[length(offset_fml)]], newdata)
      } else {
        stop("Predict can't be applied to this estimation because the offset (", deparse_long(offset), ") cannot be evaluated for the new data. Use a formula for the offset in the first estimation to avoid this.")
      }

    }

  }

  value_predicted = value_fixef + value_linear + value_NL + value_offset

  if(type == "link" || object$method_type == "feols"){
    res = value_predicted
  } else if(object$method_type == "feNmlm") {
    # Now the expected predictor
    family = object$family
    famFuns = switch(family,
             poisson = ml_poisson(),
             negbin = ml_negbin(),
             logit = ml_logit(),
             gaussian = ml_gaussian())

    if(family == "gaussian"){
      exp_value = 0
    } else {
      exp_value = exp(value_predicted)
    }

    res = famFuns$expected.predictor(value_predicted, exp_value)
  } else {
    res = object$family$linkinv(value_predicted)
  }


  #
  # se.fit
  #

  if(se.fit){

    if(!is.null(object$fixef_vars)){
      stop("The standard-errors (SEs) of the prediction cannot be computed in the presence of fixed-effects. To obtain the SEs, you would need to include the FEs as standard factors in the model.")
    }

    if(!is.null(NL_fml)){
      stop("The standard-errors (SEs) of the prediction cannot be computed in models containing non-linear in parameter elements.")
    }

    # The matrix has been created already

    V_raw = vcov(object, attr = TRUE, vcov = vcov, ssc = ssc)
    V = V_raw[var_keep, var_keep, drop = FALSE]
    X = matrix_linear[, var_keep, drop = FALSE]

    var.fit = rowSums((X %*% V) * X)
    se.fit = sqrt(var.fit)

    res = data.frame(fit = res, se.fit = se.fit)

    if(interval != "none"){
      fact = fixest_CI_factor(object, level, V_raw)

      if(interval == "prediction"){
        w = object$weights

        if(!is.null(w)){
          var.u = cpp_ssq(resid(object), w) / w
        } else {
          var.u = cpp_ssq(resid(object))
        }
        var.u = var.u / degrees_freedom(object, "resid")

        se_obs = sqrt(var.fit + var.u)

      } else {
        se_obs = se.fit
      }

      res$ci_low  = res$fit + fact[1] * se_obs
      res$ci_high = res$fit + fact[2] * se_obs
    }

  }

  res
}


#' Confidence interval for parameters estimated with `fixest`
#'
#' This function computes the confidence interval of parameter estimates obtained from a model estimated with [`femlm`], [`feols`] or [`feglm`].
#'
#' @inheritParams nobs.fixest
#' @inheritParams vcov.fixest
#'
#' @param parm The parameters for which to compute the confidence interval (either an integer vector OR a character vector with the parameter name). If missing, all parameters are used.
#' @param level The confidence level. Default is 0.95.
#' @param coef.col Logical, default is `FALSE`. If `TRUE` the column `coefficient` is inserted in the first position containing the coefficient names.
#'
#' @return
#' Returns a data.frame with two columns giving respectively the lower and upper bound of the confidence interval. There is as many rows as parameters.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 fixed-effects)
#' est_pois = femlm(Euros ~ log(dist_km) + log(Year) | Origin + Destination +
#'                  Product, trade)
#'
#' # confidence interval with "normal" VCOV
#' confint(est_pois)
#'
#' # confidence interval with "clustered" VCOV (w.r.t. the Origin factor)
#' confint(est_pois, se = "cluster")
#'
#'
confint.fixest = function(object, parm, level = 0.95, vcov, se, cluster,
              ssc = NULL, coef.col = FALSE, ...){

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = c("parm", "level", "se", "cluster"),
            valid_args = c("forceCovariance", "keepBounded"))
  }

  # Control
  check_arg(level, "numeric scalar gt{0.5} lt{1}")

  dots = list(...)

  IS_INTERNAL = isTRUE(dots$internal)
  if(IS_INTERNAL){
    if(length(object$coefficients) == 0){
      return(NULL)
    }
  }

  # The proper SE
  sum_object = summary(object, vcov = vcov, se = se, cluster = cluster, ssc = ssc, ...)
  coeftable = sum_object$coeftable

  se_all = coeftable[, 2]
  coef_all = coeftable[, 1]

  # the parameters for which we should compute the confint
  all_params = rownames(coeftable)

  # Ensuring names are all right
  names(se_all) = names(coef_all) = all_params

  if(missing(parm)){
    parm_use = all_params
  } else if(is.numeric(parm)){
    if(any(parm %% 1 != 0)){
      stop("If the argument 'parm' is numeric, it must be integers.")
    }

    parm_use = unique(na.omit(all_params[parm]))
    if(length(parm_use) == 0){

      if(IS_INTERNAL){
        return(NULL)
      }

      stop("There are ", length(all_params), " coefficients, the argument 'parm' does not correspond to any of them.")
    }
  } else if(is.character(parm)){
    parm_pblm = setdiff(parm, all_params)
    if(length(parm_pblm) > 0 && !IS_INTERNAL){
      stop("some parameters of 'parm' have no estimated coefficient: ", paste0(parm_pblm, collapse=", "), ".")
    }

    parm_use = intersect(parm, all_params)

    if(IS_INTERNAL && length(parm_use) == 0){
      return(NULL)
    }
  }

  # multiplicative factor
  fact = fixest_CI_factor(object, level, sum_object$cov.scaled)

  # The confints
  # Note that for glm models, there is no profiling
  lower_bound = coef_all[parm_use] + fact[1] * se_all[parm_use]
  upper_bound = coef_all[parm_use] + fact[2] * se_all[parm_use]

  val = (1 - level) / 2
  bound_names = paste0(round(100*c(val, 1-val), 1), " %")

  if(coef.col){
    res = data.frame(coefficient = parm_use, lower_bound, upper_bound, row.names = NULL)
    names(res)[-1] = bound_names
  } else {
    res = data.frame(lower_bound, upper_bound, row.names = parm_use)
    names(res) = bound_names
  }

  attr(res, "type") = attr(se_all, "type")

  res
}

#' Updates a `fixest` estimation
#'
#' Updates and re-estimates a `fixest` model (estimated with [`femlm`], [`feols`] or [`feglm`]). This function updates the formulas and use previous starting values to estimate a new `fixest` model. The data is obtained from the original `call`.
#'
#' @method update fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param fml.update Changes to be made to the original argument `fml`. See more information on [`update.formula`][stats::update.formula]. You can add/withdraw both variables and fixed-effects. E.g. `. ~ . + x2 | . + z2` would add the variable `x2` and the cluster `z2` to the former estimation.
#' @param nframes (Advanced users.) Defaults to 1. Number of frames up the stack where to perform the evaluation of the updated call. By default, this is the parent frame.
#' @param evaluate Logical, default is `TRUE`. If `FALSE`, only the updated call is returned.
#' @param ... Other arguments to be passed to the functions [`femlm`], [`feols`] or [`feglm`].
#'
#' @return
#' It returns a `fixest` object (see details in [`femlm`], [`feols`] or [`feglm`]).
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`predict.fixest`], [`summary.fixest`], [`vcov.fixest`], [`fixef.fixest`].
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Example using trade data
#' data(trade)
#'
#' # main estimation
#' est_pois = fepois(Euros ~ log(dist_km) | Origin + Destination, trade)
#'
#' # we add the variable log(Year)
#' est_2 = update(est_pois, . ~ . + log(Year))
#'
#' # we add another fixed-effect: "Product"
#' est_3 = update(est_2, . ~ . | . + Product)
#'
#' # we remove the fixed-effect "Origin" and the variable log(dist_km)
#' est_4 = update(est_3, . ~ . - log(dist_km) | . - Origin)
#'
#' # Quick look at the 4 estimations
#' etable(est_pois, est_2, est_3, est_4)
#'
update.fixest = function(object, fml.update, nframes = 1, evaluate = TRUE, ...){
  # Update method
  # fml.update: update the formula
  # If 1) SAME DATA and 2) SAME dep.var, then we make initialisation


  if(missing(fml.update)){
    fml.update = . ~ .
  } else {
    check_arg(fml.update, "formula")
  }

  check_arg(evaluate, "logical scalar")

  if(isTRUE(object$is_fit)){
    stop("update method not available for fixest estimations obtained from fit methods.")
  }

  if(!isScalar(nframes) || nframes < 1 || nframes %% 1 != 0){
    stop("Argument 'nframes' must be a single integer greater than, or equal to, 1.")
  }

  call_new = match.call()
  dots = list(...)

  dot_names = names(dots)
  if("fixef" %in% dot_names){
    stop("Argument 'fixef' is not accepted in the 'update.fixest' method. Please make modifications to fixed-effects directly in the argument 'fml.update'. (E.g. .~.|.+v5 to add variable v5 as a fixed-effect.)")
  }

  if(any(dot_names == "")){
    call_new_names = names(call_new)
    problems = call_new[call_new_names == ""][-1]
    stop("In 'update.fixest' the arguments of '...' are passed to the function ", object$method, ", and must be named. Currently there are un-named arguments (e.g. '", deparse_long(problems[[1]]), "').")
  }

  #
  # I) Linear formula update
  #

  fml_old = object$fml
  fml_linear = update(fml_old, fml_split(fml.update, 1))

  # Family information
  if(!is.null(dots$family)){
    if(object$method_type == "feols"){
      stop("'family' is not an argument of function feols().")
    } else if(object$method %in% c("femlm", "feNmlm", "fepois", "fenegbin")){
      family_new = match.arg(dots$family, c("poisson", "negbin", "gaussian", "logit"))
    }
  }

  #
  # II) fixed-effects updates
  #

  fml_fixef = NULL

  updt_fml_parts = fml_split(fml.update, raw = TRUE)
  n_parts = length(updt_fml_parts)

  if(n_parts > 2 + (object$method_type == "feols")){
    stop("The update formula cannot have more than ", 2 + (object$method_type == "feols"), " parts for the method ", object$method, ".")
  }

  is_fe = n_parts > 1 && !is_fml_inside(updt_fml_parts[[2]])

  fixef_vars = object$fixef_vars

  if(is_fe){

    fixef_old = object$fml_all$fixef

    # I use it as text to catch the var1^var2 FEs (update does not work)
    if(is.null(fixef_old)){
      fixef_old_text = "~ 1"
    } else {
      fixef_old_text = deparse_long(fixef_old)
    }

    fixef_new_fml = fml_maker(updt_fml_parts[[2]])
    fixef_new_text = deparse_long(fixef_new_fml)

    if(fixef_new_text == "~."){
      # nothing happens
      fixef_new = fixef_old

    } else if(fixef_new_text %in% c("~0", "~1")){
      fixef_new = ~1

    } else if(grepl("\\^", fixef_old_text) || grepl("\\^", fixef_new_text)){
      # we update manually.... dammmit
      # Note that what follows does not work ONLY when you have number^var or number^number
      # and both cases don't make much sense -- I need not control for them
      fml_text_old = gsub("\\^", "__666__", fixef_old_text)
      fml_text_new = gsub("\\^", "__666__", fixef_new_text)

      fixef_new_wip = update(as.formula(fml_text_old), as.formula(fml_text_new))

      fixef_new = as.formula(gsub("__666__", "^", fixef_new_wip))
    } else {
      fixef_new = update(fixef_old, fixef_new_fml)
    }

    if(length(all.vars(fixef_new)) > 0){
      # means there is a fixed-effect
      fml_fixef = fixef_new
    }

  } else if(!is.null(fixef_vars)){
    # the formula updated:
    fml_fixef = object$fml_all$fixef

  }

  #
  # III) IV updates
  #

  if(n_parts > 2 || (n_parts == 2 && !is_fe)){

    iv_new_fml = fml_maker(updt_fml_parts[[n_parts]])

    if(!is_fml_inside(iv_new_fml)){
      stop("The third part of the update formula in 'feols' must be a formula.")
    }

    iv_old = object$fml_all$iv

    if(is.null(iv_old)){
      fml_iv = iv_new_fml

    } else {
      fml_iv = update(iv_old, iv_new_fml)
    }

  } else {
    fml_iv = object$fml_all$iv
  }


  fml_new = merge_fml(fml_linear, fml_fixef, fml_iv)


  #
  # The call
  #

  call_old = object$call

  # we drop the argument fixef from old call (now it's in the fml_new)
  call_old$fixef = NULL

  # We also drop the arguments for multiple estimations:
  call_old$split = call_old$fsplit = NULL

  # new call: call_clear
  call_clear = call_old
  for(arg in setdiff(names(call_new)[-1], c("fml.update", "nframes", "evaluate", "object"))){
    call_clear[[arg]] = call_new[[arg]]
  }

  call_clear$fml = as.call(fml_new)

  if(!evaluate) return(call_clear)

  res = eval(call_clear, parent.frame(nframes))

  res
}


#' Extract the formula of a `fixest` fit
#'
#' This function extracts the formula from a `fixest` estimation (obtained with [`femlm`], [`feols`] or [`feglm`]). If the estimation was done with fixed-effects, they are added in the formula after a pipe (\dQuote{|}). If the estimation was done with a non linear in parameters part, then this will be added in the formula in between `I()`.
#'
#'
#' @param x An object of class `fixest`. Typically the result of a [`femlm`], [`feols`] or [`feglm`] estimation.
#' @param type A character scalar. Default is `type = "full"` which gives back a formula containing the linear part of the model along with the fixed-effects (if any) and the IV part (if any). If `type = "linear"` then only the linear formula is returned. If `type = "NL"` then only the non linear in parameters part is returned.
#' @param ... Not currently used.
#'
#' @return
#' It returns a formula.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`model.matrix.fixest`], [`update.fixest`], [`summary.fixest`], [`vcov.fixest`].
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # formula with the fixed-effect variable
#' formula(res)
#'
#' # linear part without the fixed-effects
#' formula(res, "linear")
#'
#'
formula.fixest = function(x, type = c("full", "linear", "iv", "NL"), ...){
  # Extract the formula from the object
  # we add the clusters in the formula if needed

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = "type")
  }

  if(isTRUE(x$is_fit)){
    stop("formula method not available for fixest estimations obtained from fit methods.")
  }

  check_set_arg(type, "match")

  if(type == "linear"){
    return(x$fml)

  } else if(type == "NL"){

    if(!x$method == "feNmlm"){
      stop("type = 'NL' is not valid for a ", x$method, " estimation.")
    }

    NL.fml = x$NL.fml
    if(is.null(NL.fml)){
      stop("There was no nonlinear part estimated, option type = 'NL' cannot be used.")
    }

    return(NL.fml)

  } else if(type == "iv"){
    if(is.null(x$fml_all$iv)){
      stop("type = 'iv' is only available for feols estimations with IV.")
    }
  }

  # Shall I add LHS ~ RHS + NL(NL fml) | fe | iv ???
  res = merge_fml(x$fml_all$linear, x$fml_all$fixef, x$fml_all$iv)

  res
}


#' Design matrix of a `fixest` object
#'
#' This function creates the left-hand-side or the right-hand-side(s) of a [`femlm`], [`feols`] or [`feglm`] estimation.
#'
#' @method model.matrix fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param data If missing (default) then the original data is obtained by evaluating the `call`. Otherwise, it should be a `data.frame`.
#' @param type Character vector or one sided formula, default is "rhs". Contains the type of matrix/data.frame to be returned. Possible values are: "lhs", "rhs", "fixef", "iv.rhs1" (1st stage RHS), "iv.rhs2" (2nd stage RHS), "iv.endo" (endogenous vars.), "iv.exo" (exogenous vars), "iv.inst" (instruments).
#' @param na.rm Default is `TRUE`. Should observations with NAs be removed from the matrix?
#' @param subset Logical or character vector. Default is `FALSE`. If `TRUE`, then the matrix created will be restricted only to the variables contained in the argument `data`, which can then contain a subset of the variables used in the estimation. If a character vector, then only the variables matching the elements of the vector via regular expressions will be created.
#' @param as.matrix Logical scalar, default is `FALSE`. Whether to coerce the result to a matrix.
#' @param as.df Logical scalar, default is `FALSE`. Whether to coerce the result to a data.frame.
#' @param collin.rm Logical scalar, default is `TRUE`. Whether to remove variables that were found to be collinear during the estimation. Beware: it does not perform a collinearity check.
#' @param ... Not currently used.
#'
#' @return
#' It returns either a vector, a matrix or a data.frame. It returns a vector for the dependent variable ("lhs"), a data.frame for the fixed-effects ("fixef") and a matrix for any other type.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`formula.fixest`], [`update.fixest`], [`summary.fixest`], [`vcov.fixest`].
#'
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' est = feols(y ~ poly(x1, 2) + x2, base)
#' head(model.matrix(est))
#'
#' # Illustration of subset
#'
#' # subset => character vector
#' head(model.matrix(est, subset = "x1"))
#'
#' # subset => TRUE, only works with data argument!!
#' head(model.matrix(est, data = base[, "x1", drop = FALSE], subset = TRUE))
#'
#'
#'
model.matrix.fixest = function(object, data, type = "rhs", na.rm = TRUE, subset = FALSE,
                 as.matrix = FALSE, as.df = FALSE, collin.rm = TRUE, ...){
  # We evaluate the formula with the past call
  # type: lhs, rhs, fixef, iv.endo, iv.inst, iv.rhs1, iv.rhs2
  # if fixef => return a DF

  # Checking the arguments
  if(is_user_level_call()){
    validate_dots(suggest_args = c("data", "type"))
  }

  # We allow type to be used in the location of data if data is missing
  if(!missing(data) && missing(type)){
    sc = sys.call()
    if(!"data" %in% names(sc)){
      if(!is.null(data) && (is.character(data) || "formula" %in% class(data))){
        # data is in fact the type
        type = data
        data = NULL
      }
    }
  }


  type = check_set_types(type, c("lhs", "rhs", "fixef", "iv.endo", "iv.inst", "iv.exo", "iv.rhs1", "iv.rhs2"))

  if(isTRUE(object$is_fit)){
    stop("model.matrix method not available for fixest estimations obtained from fit methods.")
  }

  if(any(grepl("^iv", type)) && !isTRUE(object$iv)){
    stop("The type", enumerate_items(grep("^iv", type, value = TRUE), "s.is"), " only valid for IV estimations.")
  }

  check_arg(subset, "logical scalar | character vector no na")

  check_set_arg(as.matrix, as.df, collin.rm, "logical scalar")

  # The formulas
  fml_full = formula(object, type = "full")
  fml_linear = formula(object, type = "linear")

  # Evaluation with the data
  original_data = FALSE
  if(missnull(data)){
    original_data = TRUE

    data = fetch_data(object, "To apply 'model.matrix.fixest', ")

  }

  # control of the data
  if(is.matrix(data)){
    if(is.null(colnames(data))){
      stop("If argument 'data' is to be a matrix, its columns must be named.")
    }
    data = as.data.frame(data)
  }

  if(!"data.frame" %in% class(data)){
    stop("The argument 'data' must be a data.frame or a matrix.")
  }

  data = as.data.frame(data)

  # Panel setup
  panel__meta__info = set_panel_meta_info(object, data)

  res = list()

  if("lhs" %in% type){
    lhs = list()

    namesLHS = all.vars(fml_linear[[2]])
    if(length(pblm <- setdiff(namesLHS, names(data)))){
      stop("In 'model.matrix', to create the LHS, the variable", enumerate_items(pblm, "s.is.quote"), " not in the data set.")
    }

    lhs_text = deparse_long(fml_linear[[2]])
    lhs[[lhs_text]] = eval(fml_linear[[2]], data)

    res[["lhs"]] = as.data.frame(lhs)
  }

  if("rhs" %in% type && !isTRUE(object$onlyFixef)){
    # we kick out the intercept if there is presence of fixed-effects
    fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))

    fml = fml_linear
    if(isTRUE(object$iv)){
      fml_iv = object$fml_all$iv
      fml = .xpd(..lhs ~ ..endo + ..rhs, ..lhs = fml[[2]], ..endo = fml_iv[[2]], ..rhs = fml[[3]])
    }

    linear.mat = error_sender(fixest_model_matrix_extra(
      object = object, newdata = data, original_data = original_data,
      fml = fml, fake_intercept = fake_intercept,
      subset = subset),
      "In 'model.matrix', the RHS could not be evaluated: ")

    if(collin.rm){
      qui = which(colnames(linear.mat) %in% object$collin.var)
      if(length(qui) == ncol(linear.mat)){
        linear.mat = NULL
      } else if(length(qui) > 0){
        linear.mat =  linear.mat[, -qui, drop = FALSE]
      }

      coefs = object$coefficients
      if(length(coefs) == ncol(linear.mat) && any(colnames(linear.mat) != names(coefs))){
        # we reorder the matrix
        # This can happen in multiple estimations, where we respect the
        # order of the user

        if(all(names(coefs) %in% colnames(linear.mat))){
          linear.mat = linear.mat[, names(coefs), drop = FALSE]
        }
      }
    }

    res[["rhs"]] = linear.mat
  }

  if("fixef" %in% type){

    if(is.null(object$fixef_vars)){
      stop("In model.matrix, the type 'fixef' is only valid for models with fixed-effects. This estimation does not contain fixed-effects.")
    }

    fixef_terms_full = fixef_terms(object$fml_all$fixef)
    fixef_terms = fixef_terms_full$fml_terms

    fixef_df = error_sender(prepare_df(fixef_terms_full$fe_vars, data, fastCombine = FALSE),
                "In 'model.matrix', problem evaluating the fixed-effects part of the formula:\n")

    isSlope = any(fixef_terms_full$slope_flag != 0)
    if(isSlope){
      slope_df = error_sender(prepare_df(fixef_terms_full$slope_vars, data),
                  "In 'model.matrix', problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n")

      fixef_df = cbind(fixef_df, slope_df)
    }

    res[["fixef"]] = fixef_df
  }

  if("iv.endo" %in% type){
    fml = object$iv_endo_fml

    endo.mat = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = TRUE), "In 'model.matrix', the endogenous variables could not be evaluated: ")

    if(collin.rm){
      qui = which(colnames(endo.mat) %in% object$collin.var)
      if(length(qui) == ncol(endo.mat)){
        endo.mat = NULL
      } else if(length(qui) > 0){
        endo.mat =  endo.mat[, -qui, drop = FALSE]
      }
    }

    res[["iv.endo"]] = endo.mat
  }

  if("iv.inst" %in% type){
    fml = object$fml_all$iv

    inst.mat = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = TRUE), "In 'model.matrix', the instruments could not be evaluated: ")

    if(collin.rm){
      qui = which(colnames(inst.mat) %in% object$collin.var)
      if(length(qui) == ncol(inst.mat)){
        inst.mat = NULL
      } else if(length(qui) > 0){
        inst.mat =  inst.mat[, -qui, drop = FALSE]
      }
    }

    res[["iv.inst"]] = inst.mat
  }

  if("iv.exo" %in% type){

    fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))
    fml = object$fml_all$linear

    exo.mat = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = fake_intercept), "In 'model.matrix', the instruments could not be evaluated: ")

    if(is.atomic(exo.mat) && length(exo.mat) == 1){
      # This is the intercept only
      # Two cases:
      is_int = attr(terms(fml), "intercept")
      if(is_int && is.null(object$fixef_vars)){
        # Valid intercept
        exo.mat = matrix(1, nrow(data))
      } else {
        # should be NULL
        exo.mat = NULL
      }
    } else if(collin.rm){
      qui = which(colnames(exo.mat) %in% object$collin.var)
      if(length(qui) == ncol(exo.mat)){
        exo.mat = NULL
      } else if(length(qui) > 0){
        exo.mat =  exo.mat[, -qui, drop = FALSE]
      }
    }

    res[["iv.exo"]] = exo.mat
  }

  if("iv.rhs1" %in% type){
    # First stage

    if(!isTRUE(object$iv)){
      stop("In model.matrix, the type 'iv.rhs1' is only valid for IV models. This estimation is no IV.")
    }

    fml = object$fml
    if(object$iv_stage == 2){
      fml_iv = object$fml_all$iv
      fml = .xpd(..lhs ~ ..inst + ..rhs, ..lhs = fml[[2]], ..inst = fml_iv[[3]], ..rhs = fml[[3]])
    }

    fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))
    # iv_rhs1 = error_sender(fixest_model_matrix(fml, data, fake_intercept = fake_intercept),
    #                        "In 'model.matrix', the RHS of the 1st stage could not be evaluated: ")
    iv_rhs1 = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = fake_intercept, subset = subset), "In 'model.matrix', the RHS of the 1st stage could not be evaluated: ")

    if(collin.rm){
      qui = which(colnames(iv_rhs1) %in% object$collin.var)
      if(length(qui) == ncol(iv_rhs1)){
        iv_rhs1 = NULL
      } else if(length(qui) > 0){
        iv_rhs1 =  iv_rhs1[, -qui, drop = FALSE]
      }
    }

    res[["iv.rhs1"]] = iv_rhs1
  }

  if("iv.rhs2" %in% type){
    # Second stage

    if(!isTRUE(object$iv)){
      stop("In model.matrix, the type 'iv.rhs2' is only valid for second stage IV models. This estimation is not even IV.")
    }

    if(!object$iv_stage == 2){
      stop("In model.matrix, the type 'iv.rhs2' is only valid for second stage IV models. This estimation is the first stage.")
    }

    # I) we get the fit
    stage_1 = object$iv_first_stage

    fit_vars = c()
    for(i in seq_along(stage_1)){
      fit_vars[i] = v = paste0("fit_", names(stage_1)[i])
      data[[v]] = predict(stage_1[[i]], newdata = data, sample = "original")
    }

    # II) we create the variables

    fml = object$fml
    fml = .xpd(..lhs ~ ..fit + ..rhs, ..lhs = fml[[2]], ..fit = fit_vars, ..rhs = fml[[3]])

    fake_intercept = !is.null(object$fixef_vars) && !(!is.null(object$slope_flag) && all(object$slope_flag < 0))
    # iv_rhs2 = error_sender(fixest_model_matrix(fml, data, fake_intercept = fake_intercept),
    #                        "In 'model.matrix', the RHS of the 2nd stage could not be evaluated: ")
    iv_rhs2 = error_sender(fixest_model_matrix_extra(object = object, newdata = data, original_data = original_data, fml = fml, fake_intercept = fake_intercept, subset = subset), "In 'model.matrix', the RHS of the 2nd stage could not be evaluated: ")

    if(collin.rm){
      qui = which(colnames(iv_rhs2) %in% object$collin.var)
      if(length(qui) == ncol(iv_rhs2)){
        iv_rhs2 = NULL
      } else if(length(qui) > 0){
        iv_rhs2 =  iv_rhs2[, -qui, drop = FALSE]
      }
    }

    res[["iv.rhs2"]] = iv_rhs2
  }

  # Formatting res
  if(length(res) == 0){
    return(NULL)
  } else if(length(type) > 1){
    res = res[type]
    res = do.call(cbind, unname(res))
  } else {
    res = res[[1]]
  }

  #
  # Removing obs if needed
  #

  check_0 = FALSE
  if(original_data){

    if(na.rm == FALSE){
      # We do nothing. Or shall I add NA values for obs not
      # included in the estimation?
      if(FALSE && length(object$obs_selection) > 0){

        # we reconstruct the full vector of obs
        # and we fill with NA
        obs_id = 1:nrow(data)
        for(i in seq_along(object$obs_selection)){
          obs_id = select_obs(obs_id, object$obs_selection[[i]])
        }

        res[!1:nrow(res) %in% obs_id, ] = NA

      }

    } else {
      for(i in seq_along(object$obs_selection)){
        check_0 = TRUE
        res = select_obs(res, object$obs_selection[[i]])
      }
    }



    na.rm = FALSE
  }

  if(na.rm){

    if(is.numeric(res) || all(sapply(res, is.numeric))){
      info = cpp_which_na_inf(res, nthreads = 1)
    } else {
      info = list(any_na_inf = anyNA(res))
      if(info$any_na_inf) info$is_na_inf = !complete.cases(res)
    }

    if(info$any_na_inf){
      check_0 = TRUE
      isNA_L = info$is_na_inf

      if(sum(isNA_L) == nrow(res)){
        warning("All observations contain NA values.")
        return(res[-which(isNA_L), , drop = FALSE])
      }

      res = select_obs(res, -which(isNA_L))
    }
  }


  if(as.matrix){
    res = as.matrix(res)
  } else if(as.df){
    res = as.data.frame(res)
  } else if(identical(type, "lhs")){
    res = res[[1]]
  }

  if(check_0 && !"fixef" %in% type){
    only_0 = cpppar_check_only_0(base::as.matrix(res), nthreads = 1)
    if(all(only_0 == 1)){
      stop("After removing NAs, not a single explanatory variable is different from 0.")

    } else if(any(only_0 == 1)){
      # At that point it must be either a matrix or a DF
      # (can't be a vector)
      res = res[, only_0 == 0, drop = FALSE]
    }
  }

  res
}


#' Extract the terms
#'
#' This function extracts the terms of a `fixest` estimation, excluding the fixed-effects part.
#'
#' @param x A `fixest` object. Obtained using the functions [`femlm`], [`feols`] or [`feglm`].
#' @param ... Not currently used.
#'
#' @return
#' An object of class `c("terms", "formula")` which contains the terms representation of a symbolic model.
#'
#'
#' @examples
#'
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = feols(Sepal.Length ~ Sepal.Width*Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # Terms of the linear part
#' terms(res)
#'
#'
terms.fixest = function(x, ...){
  terms(formula(x, type = "linear"))
}


#' Extracts the weights from a `fixest` object
#'
#' Simply extracts the weights used to estimate a `fixest` model.
#'
#' @param object A `fixest` object.
#' @param ... Not currently used.
#'
#' @return
#' Returns a vector of the same length as the number of observations in the original data set. Ignored observations due to NA or perfect fit are re-introduced and their weights set to NA.
#'
#' @seealso
#' [`feols`], [`fepois`][fixest::feglm], [`feglm`], [`fenegbin`][fixest::femlm], [`feNmlm`].
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width, iris, weights = ~as.integer(Sepal.Length) - 3.99)
#' weights(est)
#'
weights.fixest = function(object, ...){
  w = object[["weights"]]

  # To comply with stats default
  if(is.null(w)) return(NULL)

  w = fill_with_na(w, object)

  w
}



#' Residual standard deviation of `fixest` estimations
#'
#' Extract the estimated standard deviation of the errors from `fixest` estimations.
#'
#' @inheritParams weights.fixest
#'
#' @return
#' Returns a numeric scalar.
#'
#' @seealso
#' [`feols`], [`fepois`][fixest::feglm], [`feglm`], [`fenegbin`][fixest::femlm], [`feNmlm`].
#'
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width, iris)
#' sigma(est)
#'
#'
#'
sigma.fixest = function(object, ...){
  sqrt(deviance(object) / (object$nobs - object$nparams))
}


#' Extracts the deviance of a fixest estimation
#'
#' Returns the deviance from a `fixest` estimation.
#'
#' @inheritParams weights.fixest
#'
#' @return
#' Returns a numeric scalar equal to the deviance.
#'
#' @seealso
#' [`feols`], [`fepois`][fixest::feglm], [`feglm`], [`fenegbin`][fixest::femlm], [`feNmlm`].
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width, iris)
#' deviance(est)
#'
#' est_pois = fepois(Petal.Length ~ Petal.Width, iris)
#' deviance(est_pois)
#'
deviance.fixest = function(object, ...){

  if(isTRUE(object$lean)){
    # LATER: recompute it
    stop("The method 'deviance.fixest' cannot be applied to 'lean' fixest objects. Please re-estimate with 'lean = FALSE'.")
  }

  method = object$method
  family = object$family
  r = object$residuals
  w = object[["weights"]]
  if(is.null(w)) w = rep(1, length(r))

  if(is.null(r) && !method %in% c("fepois", "feglm")){
    stop("The method 'deviance.fixest' cannot be applied to a 'lean' summary. Please apply it to the estimation object directly.")
  }

  if(method %in% c("feols", "feols.fit") || (method %in% c("femlm", "feNmlm") && family == "gaussian")){
    res = sum(w * r**2)

  } else if(method %in% c("fepois", "feglm")){
    res = object$deviance

  } else {
    mu = object$fitted.values
    theta = ifelse(family == "negbin", object$theta, 1)

    # dev.resids function
    if(family == "poisson"){
      dev.resids = poisson()$dev.resids

    } else if(family == "logit"){
      dev.resids = binomial()$dev.resids

    } else if(family == "negbin"){
      dev.resids = function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))

    }

    y = r + mu

    res = sum(dev.resids(y, mu, w))
  }

  res
}



#' Hat values for `fixest` objects
#'
#' Computes the hat values for [`feols`] or [`feglm`] estimations. Only works when there are no fixed-effects.
#'
#' @param model A fixest object. For instance from feols or feglm.
#' @param ... Not currently used.
#'
#' @details
#' Hat values are not available for [`fenegbin`][fixest::femlm], [`femlm`] and [`feNmlm`] estimations.
#'
#' When there are fixed-effects, the hat values of the reduced form are different from the hat values of the full model. And we cannot get costlessly the hat values of the full model from the reduced form. It would require to reestimate the model with the fixed-effects as regular variables.
#'
#' @return
#' Returns a vector of the same length as the number of observations used in the estimation.
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Width, iris)
#' head(hatvalues(est))
#'
#'
hatvalues.fixest = function(model, ...){
  # Only works for feglm/feols objects + no fixed-effects
  # When there are fixed-effects the hatvalues of the reduced form is different from
  #  the hatvalues of the full model. And we cannot get costlessly the hatvalues of the full
  #  model from the reduced form. => we need to reestimate the model with the FEs as
  #  regular variables.

  if(isTRUE(model$lean)){
    # LATER: recompute it
    stop("The method 'hatvalues.fixest' cannot be applied to 'lean' fixest objects. Please re-estimate with 'lean = FALSE'.")
  }

  if(is_user_level_call()){
    validate_dots()
  }

  method = model$method_type
  family = model$family

  msg = "hatvalues.fixest: 'hatvalues' is not implemented for estimations with fixed-effects."

  # An error is in fact nicer than a message + NA return due to the interplay with sandwich
  if(!is.null(model$fixef_id)){
    stop(msg)
  }

  if(method == "feols"){
    X = model.matrix(model)

    res = cpp_diag_XUtX(X, model$cov.iid / model$sigma2)

  } else if(method == "feglm"){
    XW = model.matrix(model) * sqrt(model$irls_weights)
    res = cpp_diag_XUtX(XW, model$cov.iid)

  } else {
    stop("'hatvalues' is currently not implemented for function ", method, ".")
  }

  res
}

####
#### sandwich ####
####


#' Extracts the scores from a fixest estimation
#'
#' Extracts the scores from a fixest estimation.
#'
#' @param x A `fixest` object, obtained for instance from [`feols`].
#' @param ... Not currently used.
#'
#' @return
#' Returns a matrix of the same number of rows as the number of observations used for the estimation, and the same number of columns as there were variables.
#'
#' @examples
#'
#' data(iris)
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Width, iris)
#' head(estfun(est))
#'
estfun.fixest = function(x, ...){
  # 'scores' is an object always contained in fixest estimations

  if(isTRUE(x$lean)){
    # LATER: recompute it
    stop("The method 'estfun.fixest' cannot be applied to 'lean' fixest objects. Please re-estimate with 'lean = FALSE'.")
  }

  x$scores
}


#' Functions exported from \pkg{sandwich} to implement \pkg{fixest} methods
#'
#' The package \pkg{fixest} does not use `estfun` or `bread` from \pkg{sandwich}, but these methods have been implemented to allow users to leverage the variances from \pkg{sandwich}.
#'
#' * Here is the help from package \pkg{sandwich}: [`estfun`][sandwich::estfun] and [`bread`][sandwich::bread]. The help from package \pkg{fixest} is here: [`estfun.fixest`] and [`bread.fixest`].
#'
#'
#' @name sandwich_reexported
#' @keywords internal
NULL

#' @rdname sandwich_reexported
#' @name estfun
NULL

#' @rdname sandwich_reexported
#' @name bread
NULL


#' Extracts the bread matrix from fixest objects
#'
#' Extracts the bread matrix from fixest objects to be used to compute sandwich variance-covariance matrices.
#'
#' @param x A `fixest` object, obtained for instance from [`feols`].
#' @param ... Not currently used.
#'
#' @return
#' Returns a matrix of the same dimension as the number of variables used in the estimation.
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Width, iris)
#' bread(est)
#'
bread.fixest = function(x, ...){

  if(is_user_level_call()){
    validate_dots()
  }

  method = x$method_type
  family = x$family

  if(method == "feols"){

    res = x$cov.iid / x$sigma2 * x$nobs

  } else if(method == "feglm"){

    res = x$cov.iid * x$nobs

  } else {
    stop("'bread' is not currently implemented for function ", method, ".")
  }

  res
}




































































