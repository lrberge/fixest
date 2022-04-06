
#' A print facility for \code{fixest} objects.
#'
#' This function is very similar to usual \code{summary} functions as it provides the table of coefficients along with other information on the fit of the estimation. The type of output can be customized by the user (using function \code{setFixest_print}).
#'
#' @method print fixest
#'
#' @param x A \code{fixest} object. Obtained using the methods \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param n Integer, number of coefficients to display. By default, only the first 8 coefficients are displayed if \code{x} does not come from \code{\link[fixest]{summary.fixest}}.
#' @param type Either \code{"table"} (default) to display the coefficients table or \code{"coef"} to display only the coefficients.
#' @param fitstat A formula or a character vector representing which fit statistic to display. The types must be valid types of the function \code{\link[fixest]{fitstat}}. The default fit statistics depend on the type of estimation (OLS, GLM, IV, with/without fixed-effect). Providing the argument \code{fitstat} overrides the default fit statistics, you can however use the point "." to summon them back. Ex 1: \code{fitstat = ~ . + ll} adds the log-likelihood to the default values. Ex 2: \code{fitstat = ~ ll + pr2} only displays the log-likelihood and the pseudo-R2.
#' @param ... Other arguments to be passed to \code{\link[fixest]{vcov.fixest}}.
#'
#' @details
#'  It is possible to set the default values for the arguments \code{type} and \code{fitstat} by using the function \code{setFixest_print}.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
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
	fromSummary = isTRUE(x$summary)

	if(!missnull(fitstat)){
	    fitstat = fitstat_validate(fitstat, TRUE)
	}

	# User options
	set_defaults("fixest_print")

	# if NOT from summary, we consider the argument 'type'
	if(!fromSummary){
	    # checking argument type
	    check_arg_plus(type, "match(coef, table)")

	    if(type == "coef"){
	        print(coef(x))
	        return(invisible())
	    }
	}

	isNegbin = x$method == "fenegbin" || (x$method_type == "feNmlm" && x$family=="negbin")

	x = summary(x, fromPrint = TRUE, ...)

	check_arg(n, "integer scalar GE{1}")

	msgRemaining = ""
	nb_coef = length(coef(x)) - isNegbin
	if(missing(n) && is.null(x$n_print)){
		if(fromSummary && !isTRUE(x$summary_from_fit)){
			n = Inf
		} else {
			if(nb_coef <= 10){
				n = 10
			} else {
				n = 8
				msgRemaining = paste0("... ", nb_coef - n, " coefficients remaining (display them with summary() or use argument n)\n")
			}
		}

	} else {
	    if(!is.null(x$n_print)) n = x$n_print

	    if(n < nb_coef){
	        msgRemaining = paste0("... ", nb_coef - n, " coefficients remaining\n")
	    }
	}

	# We also add the collinearity message
	collinearity_msg = ""
	if(!is.null(x$collin.var)){
	    n_collin = length(x$collin.var)
	    collinearity_msg = paste0("... ", n_collin, " variable", plural(n_collin, "s.was"), " removed because of collinearity (", enumerate_items(x$collin.var, nmax = 3), ifelse(n_collin > 3, " [full set in $collin.var]", ""), ")\n")
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

#' Summary of a \code{fixest} object. Computes different types of standard errors.
#'
#' This function is similar to \code{print.fixest}. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
#'
#' @inheritParams feNmlm
#' @inheritParams aggregate.fixest
#'
#' @method summary fixest
#' @param vcov Versatile argument to specify the VCOV. In general, it is either a character scalar equal to a VCOV type, either a formula of the form: \code{vcov_type ~ variables}. The VCOV types implemented are: "iid", "hetero" (or "HC1"), "cluster", "twoway", "NW" (or "newey_west"), "DK" (or "driscoll_kraay"), and "conley". It also accepts object from \code{\link[fixest]{vcov_cluster}}, \code{\link[fixest:vcov_hac]{vcov_NW}}, \code{\link[fixest:vcov_hac]{NW}}, \code{\link[fixest:vcov_hac]{vcov_DK}}, \code{\link[fixest:vcov_hac]{DK}}, \code{\link[fixest]{vcov_conley}} and \code{\link[fixest:vcov_conley]{conley}}. It also accepts covariance matrices computed externally. Finally it accepts functions to compute the covariances. See the `vcov` documentation in the \href{https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1}{vignette}.
#' @param se Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: \code{se = "cluster"}, otherwise \code{se = "iid"}. Note that this argument is deprecated, you should use \code{vcov} instead.
#' @param cluster Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over \code{var1} and \code{var2} contained in the data.frame \code{base} used for the estimation. All the following \code{cluster} arguments are valid and do the same thing: \code{cluster = base[, c("var1", "var2")]}, \code{cluster = c("var1", "var2")}, \code{cluster = ~var1+var2}. If the two variables were used as fixed-effects in the estimation, you can leave it blank with \code{vcov = "twoway"} (assuming \code{var1} [resp. \code{var2}] was the 1st [res. 2nd] fixed-effect). You can interact two variables using \code{^} with the following syntax: \code{cluster = ~var1^var2} or \code{cluster = "var1^var2"}.
#' @param stage Can be equal to \code{2} (default), \code{1}, \code{1:2} or \code{2:1}. Only used if the object is an IV estimation: defines the stage to which \code{summary} should be applied. If \code{stage = 1} and there are multiple endogenous regressors or if \code{stage} is of length 2, then an object of class \code{fixest_multi} is returned.
#' @param object A \code{fixest} object. Obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param ssc An object of class \code{ssc.type} obtained with the function \code{\link[fixest]{ssc}}. Represents how the degree of freedom correction should be done.You must use the function \code{\link[fixest]{ssc}} for this argument. The arguments and defaults of the function \code{\link[fixest]{ssc}} are: \code{adj = TRUE}, \code{fixef.K="nested"}, \code{cluster.adj = TRUE}, \code{cluster.df = "min"}, \code{t.df = "min"}, \code{fixef.force_exact=FALSE)}. See the help of the function \code{\link[fixest]{ssc}} for details.
#' @param .vcov A user provided covariance matrix or a function computing this matrix. If a matrix, it must be a square matrix of the same number of rows as the number of variables estimated. If a function, it must return the previously mentioned matrix.
#' @param lean Logical, default is \code{FALSE}. Used to reduce the (memory) size of the summary object. If \code{TRUE}, then all objects of length N (the number of observations) are removed from the result. Note that some \code{fixest} methods may consequently not work when applied to the summary.
#' @param forceCovariance (Advanced users.) Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option to force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded (Advanced users -- \code{feNmlm} with non-linear part and bounded coefficients only.) Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param n Integer, default is 1000. Number of coefficients to display when the print method is used.
#' @param ... Only used if the argument \code{.vocv} is provided and is a function: extra arguments to be passed to that function.
#'
#' @section Compatibility with \pkg{sandwich} package:
#' The VCOVs from \code{sandwich} can be used with \code{feols}, \code{feglm} and \code{fepois} estimations. If you want to have a \code{sandwich} VCOV when using \code{summary.fixest}, you can use the argument \code{vcov} to specify the VCOV function to use (see examples).
#' Note that if you do so and you use a formula in the \code{cluster} argument, an innocuous warning can pop up if you used several non-numeric fixed-effects in the estimation (this is due to the function \code{\link[stats]{expand.model.frame}} used in \code{sandwich}).
#'
#' @return
#' It returns a \code{fixest} object with:
#' \item{cov.scaled}{The new variance-covariance matrix (computed according to the argument \code{se}).}
#' \item{se}{The new standard-errors (computed according to the argument \code{se}).}
#' \item{coeftable}{The table of coefficients with the new standard errors.}
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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
                          se = NULL, keepBounded = FALSE, n = 1000,
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
                                                     n = n, nthreads = nthreads, iv = TRUE)

                    stage_names[length(stage_names) + 1] = paste0("First stage: ", names(object$iv_first_stage)[i])
                }

	        } else {
	            # We keep the information on clustering => matters for wald tests of 1st stage
	            my_res = summary(object, vcov = vcov, ssc = ssc, lean = lean,
	                             forceCovariance = forceCovariance, n = n, nthreads = nthreads, iv = TRUE)

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

	    index = list("iv" = length(res))
	    all_names = list("iv" = stage_names)
	    res_multi = setup_multi(index, all_names, res)
	    attr(res_multi, "print_request") = "long"

	    return(res_multi)
	}


	# The new VCOV
	if(skip_vcov){
	    vcov = object$cov.scaled

	} else {
	    vcov = vcov.fixest(object, vcov = vcov, ssc = ssc, forceCovariance = forceCovariance,
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
            data = attr(my_res, "data")
            for(j in seq_along(data)){
                res[[length(res) + 1]] = data[[j]]
            }
        } else {
            res[[length(res) + 1]] = my_res
        }
    }

    # We return a simple list
    class(res) = NULL

    res
}

#' Obtain various statistics from an estimation
#'
#' Set of functions to directly extract some commonly used statistics, like the p-value or the table of coefficients, from estimations. This was first implemented for \code{fixest} estimations, but has some support for other models.
#'
#' @inheritParams etable
#'
#' @param object An estimation. For example obtained from \code{\link[fixest]{feols}}.
#' @param se [Fixest specific.] Character scalar. Which kind of standard error should be computed: \dQuote{iid}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are fixed-effects in the estimation: \code{se = "cluster"}, otherwise \code{se = "iid"}. Note that this argument is not needed if the argument \code{cluster} is present.
#' @param cluster [Fixest specific.] Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over \code{var1} and \code{var2} contained in the data.frame \code{base} used for the estimation. All the following \code{cluster} arguments are valid and do the same thing: \code{cluster = base[, c("var1, "var2")]}, \code{cluster = c("var1, "var2")}, \code{cluster = ~var1+var2}. If the two variables were used as clusters in the estimation, you could further use \code{cluster = 1:2} or leave it blank with \code{se = "twoway"} (assuming \code{var1} [resp. \code{var2}] was the 1st [res. 2nd] cluster).
#' @param ... Other arguments to be passed to \code{summary}.
#'
#' @details
#' This set of functions is primarily constructed for \code{fixest} estimations. Although it can work for non-\code{fixest} estimations, support for exotic estimation procedures that do not report standardized coefficient tables is highly limited.
#'
#' @return
#' Returns a table of coefficients, with in rows the variables and four columns: the estimate, the standard-error, the t-statistic and the p-value.
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
#'
#'
coeftable = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){
    # We don't explicitly refer to the other arguments

    check_arg(keep, drop, order, "NULL character vector no na")

    IS_FIXEST = "fixest" %in% class(object)
    # IS_FIXEST_MULI = "fixest_multi" %in% class(object)
    # IS_FIXEST_LIST = "fixest_multi" %in% class(object)

    if(IS_FIXEST){
        if(!isTRUE(object$summary) || !(all_missing(vcov, ssc, cluster) && ...length() > 0)){
            object = summary(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)
        }
    } else {
        if(!any(grepl("summary", class(object)))){
            # Known issue:
            # - if vcov/.../order are also args of summary, that will be a problem

            args = list(object = object)

            if(...length() > 0){
                args = append(args, list(...))
            }

            object = do.call("summary", args)
        }
    }

    # Let's find out the coefficients table
    if(IS_FIXEST){
        res = object$coeftable
    } else {
        list_mat = object[sapply(object, is.matrix)]

        ok = FALSE
        for(i in seq_along(list_mat)){
            mat = list_mat[[i]]
            if(!is.null(colnames(mat)) && any(grepl("(?i)(estimate|value|Pr\\()", colnames(mat)))){
                ok = TRUE
                res = mat
            }
        }

        if(ok == FALSE){
            stop("Sorry, the coeffficients table could not be extracted.")
        }

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

#' @describeIn coeftable Extracts the p-value of an estimation
pvalue = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){

    check_arg(keep, drop, order, "NULL character vector no na")

    mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4), sorry.")
    }

    res = mat[, 4]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = names(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(numeric(0))
        }

        res = res[r_names]
    }

    res
}

#' @describeIn coeftable Extracts the t-statistics of an estimation
tstat = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){

    check_arg(keep, drop, order, "NULL character vector no na")

    mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4), sorry.")
    }

    res = mat[, 3]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = names(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(numeric(0))
        }

        res = res[r_names]
    }

    res
}

#' @describeIn coeftable Extracts the standard-error of an estimation
se = function(object, vcov = NULL, ssc = NULL, cluster = NULL, keep, drop, order, ...){

    check_arg(keep, drop, order, "NULL character vector no na")

    if(is.matrix(object) && nrow(object) == ncol(object)){
        # special case => object is a VCOV matrix and NOT an estimation
        res = sqrt(diag(object))

    } else {

        mat = coeftable(object, vcov = vcov, ssc = ssc, cluster = cluster, ...)

        if(ncol(mat) != 4){
            stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4), sorry.")
        }

        res = mat[, 2]
        if(is.null(names(res))) {
            names(res) = rownames(mat)
        }
    }

    if(!missnull(keep) || !missnull(drop) || !missnull(order)){
        r_names = names(res)
        r_names = keep_apply(r_names, keep)
        r_names = drop_apply(r_names, drop)
        r_names = order_apply(r_names, order)

        if(length(r_names) == 0){
            return(numeric(0))
        }

        res = res[r_names]
    }

    res
}

#' Summary method for fixed-effects coefficients
#'
#' This function summarizes the main characteristics of the fixed-effects coefficients. It shows the number of fixed-effects that have been set as references and the first elements of the fixed-effects.
#'
#' @method summary fixest.fixef
#'
#' @param object An object returned by the function \code{\link[fixest]{fixef.fixest}}.
#' @param n Positive integer, defaults to 5. The \code{n} first fixed-effects for each fixed-effect dimension are reported.
#' @param ... Not currently used.
#'
#' @return
#' It prints the number of fixed-effect coefficients per fixed-effect dimension, as well as the number of fixed-effects used as references for each dimension, and the mean and variance of the fixed-effect coefficients. Finally, it reports the first 5 (arg. \code{n}) elements of each fixed-effect.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{femlm}}, \code{\link[fixest]{fixef.fixest}}, \code{\link[fixest]{plot.fixest.fixef}}.
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


#' Extract the Fixed-Effects from a \code{fixest} estimation.
#'
#' This function retrieves the fixed effects from a \code{fixest} estimation. It is useful only when there are one or more fixed-effect dimensions.
#'
#' @inheritParams feNmlm
#'
#' @param object A \code{fixest} estimation (e.g. obtained using \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}).
#' @param notes Logical. Whether to display a note when the fixed-effects coefficients are not regular.
#' @param sorted Logical, default is \code{TRUE}. Whether to order the fixed-effects by their names. If \code{FALSE}, then the order used in the demeaning algorithm is used.
#'
#' @details
#' If the fixed-effect coefficients not regular, then several reference points need to be set, leading to the coefficients to be NOT interpretable. If this is the case, then a warning is raised.
#'
#' @return
#' A list containing the vectors of the fixed effects.
#'
#' If there is more than 1 fixed-effect, then the attribute \dQuote{references} is created. This is a vector of length the number of fixed-effects, each element contains the number of coefficients set as references. By construction, the elements of the first fixed-effect dimension are never set as references. In the presence of regular fixed-effects, there should be Q-1 references (with Q the number of fixed-effects).
#'
#' @seealso
#' \code{\link[fixest]{plot.fixest.fixef}}. See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effect coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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
#' The package \pkg{fixest} uses the \code{fixef} method from \pkg{nlme}. Unfortunately, re-exporting this method is required in order not to attach package \pkg{nlme}.
#'
#' \itemize{
#' \item Here is the help from package \pkg{nlme}: \code{\link[nlme:fixed.effects]{fixef}}. The help from package \pkg{fixest} is here: \code{\link[fixest]{fixef.fixest}}.
#' }
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
#' This function plots the 5 fixed-effects with the highest and lowest values, for each of the fixed-effect dimension. It takes as an argument the fixed-effects obtained from the function \code{\link{fixef.fixest}} after an estimation using \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @method plot fixest.fixef
#'
#' @param x An object obtained from the function \code{\link{fixef.fixest}}.
#' @param n The number of fixed-effects to be drawn. Defaults to 5.
#' @param ... Not currently used.
#'
#' Note that the fixed-effect coefficients might NOT be interpretable. This function is useful only for fully regular panels.
#'
#' If the data are not regular in the fixed-effect coefficients, this means that several \sQuote{reference points} are set to obtain the fixed-effects, thereby impeding their interpretation. In this case a warning is raised.
#'
#' @seealso
#' \code{\link[fixest]{fixef.fixest}} to extract clouster coefficients. See also the main estimation function \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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


#' Collinearity diagnostics for \code{fixest} objects
#'
#' In some occasions, the optimization algorithm of \code{\link[fixest]{femlm}} may fail to converge, or the variance-covariance matrix may not be available. The most common reason of why this happens is colllinearity among variables. This function helps to find out which set of variables is problematic.
#'
#'
#' @param x A \code{fixest} object obtained from, e.g. functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param verbose An integer. If higher than or equal to 1, then a note is prompted at each step of the algorithm. By default \code{verbose = 0} for small problems and to 1 for large problems.
#'
#' @details
#' This function tests: 1) collinearity with the fixed-effect variables, 2) perfect multi-collinearity between the variables, 4) perfect multi-collinearity between several variables and the fixed-effects, and 4) identification issues when there are non-linear in parameters parts.
#'
#' @return
#' It returns a text message with the identified diagnostics.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Creating an example data base:
#' set.seed(1)
#' fe_1 = sample(3, 100, TRUE)
#' fe_2 = sample(20, 100, TRUE)
#' x = rnorm(100, fe_1)**2
#' y = rnorm(100, fe_2)**2
#' z = rnorm(100, 3)**2
#' dep = rpois(100, x*y*z)
#' base = data.frame(fe_1, fe_2, x, y, z, dep)
#'
#' # creating collinearity problems:
#' base$v1 = base$v2 = base$v3 = base$v4 = 0
#' base$v1[base$fe_1 == 1] = 1
#' base$v2[base$fe_1 == 2] = 1
#' base$v3[base$fe_1 == 3] = 1
#' base$v4[base$fe_2 == 1] = 1
#'
#' # Estimations:
#'
#' # Collinearity with the fixed-effects:
#' res_1 = femlm(dep ~ log(x) + v1 + v2 + v4 | fe_1 + fe_2, base)
#' collinearity(res_1)
#'
#' # => collinearity with the first fixed-effect identified, we drop v1 and v2
#' res_1bis = femlm(dep ~ log(x) + v4 | fe_1 + fe_2, base)
#' collinearity(res_1bis)
#'
#' # Multi-Collinearity:
#' res_2 =  femlm(dep ~ log(x) + v1 + v2 + v3 + v4, base)
#' collinearity(res_2)
#'
#'
collinearity = function(x, verbose){
	# x: fixest estimation

    # stop("Sorry, it does not work. A new version will hopefully come soon.")

	if(!inherits(x, "fixest")){
		stop("Argument 'x' must be a fixest object.")
	}

	# I) (linear) collinearity with fixed-effects
	# II) (linear) multi collinearity
	# III) (non-linear) overidentification

	# flags
    isFixef = !is.null(x$fixef_vars)
    isFE = FALSE
    isSlope = FALSE
    if(isFixef){
        isSlope = !is.null(x$fixef_terms)
        if(isSlope){
            fixef_terms = x$fixef_terms
            slope_flag = grepl("\\[", fixef_terms)
            isFE = any(!slope_flag)
        } else {
            isFE = TRUE
        }
    }

	linear_fml = fml_split(formula(x, "linear"), 2, split.lhs = TRUE)

	rhs_fml = fml_split(formula(x, "linear"), 1)
	linear.varnames = all.vars(rhs_fml[[3]])

	isLinear = length(linear.varnames) > 0

	NL_fml = x$NL.fml
	isNL = !is.null(NL_fml)
	coef = x$coefficients

	# Getting the data
	data = fetch_data(x, "To apply function 'collinearity', ")

	if(is.matrix(data)){
	    data = as.data.frame(data)
	} else {
	    class(data) = "data.frame"
	}

	if(isFE){
		linear_fml = update(linear_fml, ~ . + 1)
	}

	# Panel setup
	panel__meta__info = set_panel_meta_info(x, data)

	if(isLinear || isFixef || "(Intercept)" %in% names(coef)){
		# linear.matrix = model.matrix(linear_fml, data)
		linear.matrix = fixest_model_matrix(rhs_fml, data)
	}

	for(i in seq_along(x$obs_selection)){
	    linear.matrix = linear.matrix[x$obs_selection[[i]], , drop = FALSE]
	}

	if(isLinear){
	    # We do that to drop interaction variables that should not be there any more
	    # if factors with only NA values
	    varkeep = intersect(names(x$coefficients), colnames(linear.matrix))
	    if(length(varkeep) < ncol(linear.matrix)){
	        linear.matrix = linear.matrix[, varkeep, drop = FALSE]
	    }
	}


	if(isLinear){
	    # We center and scale to have something comparable across data sets
	    first_row = linear.matrix[1, ]
	    linear.matrix = scale(linear.matrix, center = FALSE, scale = TRUE)
	    # The constants => we set it to 1
	    constant_id = apply(linear.matrix, 2, anyNA)
	    constant_value = first_row[constant_id]
	    linear.matrix[is.na(linear.matrix)] = 1
	    linear_mat_noIntercept = linear.matrix[, -1, drop = FALSE]
	}

	n_obs = nrow(data)
	Q = length(x$fixef_id)

	# information display
	ccat = function(...) NULL
	if(missing(verbose)){
		verbose = FALSE
		# verbose only if task is long
		if(100 * nrow(linear.matrix) * (ncol(linear.matrix)**2 * Q**2) >= 1e9){
			verbose = TRUE
		}
	}
	if(verbose) ccat = cat

	#
	# 0) Data preparation (FE/slopes)
	#

	# Variables for FE / Slopes
	if(isSlope){
	    fixef_terms = x$fixef_terms
	    terms_full = extract_fe_slope(fixef_terms)
	    fixef_vars = terms_full$fixef_vars
	    slope_fe = terms_full$slope_fe
	    fe_all = terms_full$fe_all
	    slope_vars = terms_full$slope_vars
	    slope_terms = terms_full$slope_terms

	    # dictionary mapping fe var names to the ids of id_dummies_vect
	    dict_fe = 1:Q
	    names(dict_fe) = x$fixef_vars
	    slope_vars_unik = unique(slope_vars)

	    if(any(!slope_vars_unik %in% names(data))){
	        var_pblm = setdiff(slope_vars_unik, names(data))
	        stop("To check collinearity, we need to fetch some variables in the original dataset (", deparse_long(x$call$data), "). However, the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset any more.")
	    }

	    slope_var_list = list()
	    for(i in 1:length(slope_vars_unik)){
	        variable = all.vars(str2lang(slope_vars_unik[i]))

	        # as.numeric => we'll use cpp so required
	        svar = as.numeric(as.vector(eval(str2lang(slope_vars_unik[i]), data)))
	        if(length(svar) == 1) svar = rep(svar, n_obs)

	        slope_var_list[[slope_vars_unik[i]]] = svar
	    }

	    Q_slope = length(slope_fe)
	    slope_fe_id = x$fixef_id[slope_fe]
	    slope_var_all = slope_var_list[slope_vars]

	    Q_fe = length(fixef_vars)
	    fe_id = x$fixef_id[fixef_vars]

	} else if(isFE){
	    Q_fe = length(x$fixef_id)
	    fe_id = x$fixef_id
	}


	#
	# I) collinearity with clusters
	#

	ccat("Checking Collinearity: ")

	#
	# Variables are constant or equal to 0
	#

	if(isLinear){

	    if(isFE){
	        if(any(constant_id[-1])){
	            # var_problem = colnames(linear_mat_noIntercept)[is_const]
	            var_problem = colnames(linear_mat_noIntercept)[constant_id[-1]]
	            message = paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant, thus collinear with the fixed-effects.")
	            print(message)
	            return(invisible(message))
	        }

	    } else {

	        isIntercept = attr(terms(x$fml),"intercept")

	        if(isIntercept && any(constant_id[-1])){
	            var_problem = colnames(linear_mat_noIntercept)[constant_id[-1]]
	            message = paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant, thus collinear with the intercept.")
	            print(message)
	            return(invisible(message))
	        }

	        if(any(constant_value == 0)){
	            # constant and equal to 0
    	        var_problem = colnames(linear.matrix)[first_row == 0 & constant_id]
    	        message = paste0("Variable", enumerate_items(var_problem, "s.is.quote"), " constant and equal to 0.")

    	        print(message)
    	        return(invisible(message))
	        }
	    }
	}

	if(isFE && isLinear){
		ccat("simple with fixed-effects:")
		# We project each variable onto the cluster subspace

		cluster = fe_id
		Q_fe = length(cluster)
		for(q in 1:Q_fe){
			ccat(".")
			dum = cluster[[q]]
			k = max(dum)

			value = cpp_tapply_sum(Q = k, x = linear_mat_noIntercept, dum = dum)
			nb_per_cluster = cpp_table(Q = k, dum = dum)

			# residuals of the linear projection on the cluster space
			residuals = linear_mat_noIntercept - (value/nb_per_cluster)[dum, ]

			max_residuals = apply(abs(residuals), 2, max)

			if(any(max_residuals < 1e-6)){
				ccat("\n")
				varnames = colnames(linear_mat_noIntercept)
				collin_var = varnames[max_residuals < 1e-6]
				message = paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with fixed-effects '", names(cluster)[q], "'.")

				print(message)
				return(invisible(message))

			}

		}
		ccat("OK")
	}

	if(isSlope && isLinear){
	    ccat("simple with variables with varying slopes:")

	    for(q in 1:Q_slope){
	        ccat(".")
	        dum = fe_id[[q]]
	        my_var = slope_var_all[[q]]
	        k = max(dum)
	        value = cpp_tapply_sum(Q = k, x = linear_mat_noIntercept*my_var, dum = dum)

	        denom = as.vector(tapply(my_var**2, dum, sum))

	        # residuals of the linear projection on the cluster space
	        residuals = linear_mat_noIntercept - (value/denom)[dum, ]*my_var

	        max_residuals = apply(abs(residuals), 2, max)

	        if(any(max_residuals < 1e-6)){
	            ccat("\n")
	            varnames = colnames(linear_mat_noIntercept)
	            collin_var = varnames[max_residuals < 1e-6]

	            message = paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with variable with varying slope '", slope_vars[q], "' (on '", slope_fe[q], "').")

	            print(message)
	            return(invisible(message))

	        }

	    }
	    ccat("OK")
	}

	#
	# II) perfect multicollinearity
	#

	name2change = grepl("[^[:alnum:]\\._]", colnames(linear.matrix))
	dict_name = colnames(linear.matrix)[!grepl("(Intercept)", colnames(linear.matrix))]
	if(any(name2change)){
		linbase = as.data.frame(linear.matrix)
		names(linbase)[name2change] = gsub("[^[:alnum:]\\._]", "_", names(linbase)[name2change])
	} else {
		linbase = as.data.frame(linear.matrix)
	}
	linearVars = names(linbase)[!grepl("Intercept", names(linbase))]
	names(dict_name) = linearVars
	dict_name["_Intercept_"] = "(Intercept)"

	# for multicol without cluster
	mat_base = as.matrix(linbase)

	# we add possibly missing variables
	varmiss = setdiff(all.vars(linear_fml), names(linbase))
	for(v in varmiss) linbase[[v]] = data[[v]]

	if(isLinear && length(linearVars) >= 2){
		ccat(ifelse(Q >= 1, ", ", " "), "multiple:")
		for(v in linearVars){
			ccat(".")

			i = which(colnames(mat_base) == v)
			res = ols_fit(y = mat_base[, i], X = mat_base[, -i, drop = FALSE], w = 1, collin.tol = 1e-10,  nthreads = 1)

			max_residuals = max(abs(res$residuals))

			if(max_residuals < 1e-4){
				ccat("\n")
				coef_lm = res$coefficients
				vars = colnames(mat_base)[-i]
				collin_var = vars[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
				message = paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ".")

				print(message)
				return(invisible(message))
			}
		}
		ccat("OK")
	}

	#
	# II.b) perfect multicollinearity + fixed-effects
	#

	# linearVars = setdiff(colnames(linear.matrix), "(Intercept)")
	if(isLinear && length(linearVars) >= 2 && isFixef){
		ccat(", multiple with cluster")

		dum_names = names(x$fixef_id)
		n_clust = length(dum_names)
		new_dum_names = paste0("__DUM_", 1:n_clust)
		for(i in 1:n_clust){
			linbase[[paste0("__DUM_", i)]] = x$fixef_id[[i]]
		}

		for(v in linearVars){
			ccat(".")
			fml2estimate = as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))

			for(id_cluster in 1:n_clust){

				res = feols(fml2estimate, linbase, fixef = new_dum_names[id_cluster], warn = FALSE)

				max_residuals = max(abs(res$residuals))

				if(max_residuals < 1e-4){
					ccat("\n")
					coef_lm = coef(res)
					collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
					message = paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ", together with the fixed-effects ", dum_names[id_cluster], ".")

					print(message)
					return(invisible(message))
				}
			}
		}
		ccat("OK")
	}

	#
	# III) NL problem
	#

	if(isNL){
		NL_vars = all.vars(NL_fml)
		varNotHere = setdiff(NL_vars, c(names(coef), names(data)))
		if(length(varNotHere) > 0){
			stop("Some variables used to estimate the model (in the non-linear formula) are missing from the original data: ", paste0(varNotHere, collapse = ", "), ".")
		}

		var2send = intersect(NL_vars, names(data))
		env = new.env()
		for(var in var2send){
			assign(var, data[[var]], env)
		}
	}

	if(isNL && (length(coef) >= 2 || isFixef)){
		ccat(", in non-linear term:")
		if(isFixef){
			# we add the constant
			coef["CONSTANT"] = 1
			data$CONSTANT = 1
			linear.matrix = model.matrix(update(linear_fml, ~.-1+CONSTANT), data)
		}

		coef_name = names(coef)

		NL_coef = setdiff(all.vars(NL_fml), names(data))
		# L_coef = setdiff(coef_name, NL_coef)
		L_coef = colnames(linear.matrix)

		#
		# We compute mu:
		#

		mu = 0
		if(length(L_coef) > 0){
			mu = linear.matrix %*% coef[L_coef]
		}

		for(iter_coef in NL_coef){
			assign(iter_coef, coef[iter_coef], env)
		}

		# Evaluation of the NL part
		value_NL = eval(NL_fml[[2]], env)
		mu = mu + value_NL
		data$mu = mu

		if(var(mu) == 0){
			message = "Variance of the NL part is 0."
			print(message)
			return(invisible(message))
		}

		#
		# The loop
		#

		for(var in coef_name){
			ccat(".")
			# we modify the coef of var, then we fix it:
			# if we can find a set of other parameters that give the same fit,
			# then there is an identification issue

			if(var %in% L_coef){
				# if linear coef: we modify the fml and add an offset
				if(length(L_coef) == 1){
					fml = mu ~ 0
				} else {
					fml = as.formula(paste0("mu ~ 0+", paste0(setdiff(L_coef, var), collapse = "+")))
				}

				# offset with the new value
				offset = as.formula(paste0("~", coef[var] + 1, "*", var))

				res = feNmlm(fml, data = data, family = "gaussian", NL.fml = NL_fml, NL.start = as.list(coef), offset = offset)

			} else {
				# NL case:
				# we modify both fml and NL.fml

				if(isFixef){
					fml = update(x$fml, mu ~ . - 1 + CONSTANT)
				} else {
					fml = update(x$fml, mu ~ .)
				}

				# the new formula with the new variable
				NL_fml_char = as.character(NL_fml)[2]
				NL_fml_new = as.formula(paste0("~", gsub(paste0("\\b", var, "\\b"), 1 + coef[var], NL_fml_char)))

				if(length(NL_coef) == 1){
					# there is no more parameter to estimate => offset
					res = femlm(fml, data = data, family = "gaussian", offset = NL_fml_new)
				} else {
					res = feNmlm(fml, data = data, family = "gaussian", NL.fml = NL_fml_new, NL.start = as.list(coef))
				}

			}

			sum_resids = sum(abs(res$residuals))
			if(sum_resids < 1e-4){
				coef_esti = coef(res)
				coef_diff = abs(coef_esti - coef[names(coef_esti)])
				collin_var = names(coef_diff)[coef_diff > 1e-3]
				message = paste0("Coefficients ", show_vars_limited_width(c(var, collin_var)), " are not uniquely identifed.")

				print(message)
				return(invisible(message))
			}

		}
		ccat("OK")
	}

	ccat("\n")
	message = "No visible collinearity problem. (Doesn't mean there's none!)"
	print(message)
	return(invisible(message))

}



#' Treated and control sample descriptives
#'
#' This function shows the means and standard-deviations of several variables conditional on whether they are from the treated or the control group. The groups can further be split according to a pre/post variable. Results can be seamlessly be exported to Latex.
#'
#'
#' @param fml Either a formula of the type \code{var1 + ... + varN ~ treat} or \code{var1 + ... + varN ~ treat | post}. Either a data.frame/matrix containing all the variables for which the means are to be computed (they must be numeric of course). Both the treatment and the post variables must contain only exactly two values. You can use a point to select all the variables of the data set: \code{. ~ treat}.
#' @param base A data base containing all the variables in the formula \code{fml}.
#' @param treat_var Only if argument \code{fml} is *not* a formula. The vector identifying the treated and the control observations (the vector can be of any type but must contain only two possible values). Must be of the same length as the data.
#' @param post_var Only if argument \code{fml} is *not* a formula. The vector identifying the periods (pre/post) of the observations (the vector can be of any type but must contain only two possible values). The first value (in the sorted sense) of the vector is taken as the pre period. Must be of the same length as the data.
#' @param treat_first Which value of the 'treatment' vector should appear on the left? By default the max value appears first (e.g. if the treatment variable is a 0/1 vector, 1 appears first).
#' @param tex Should the result be displayed in Latex? Default is \code{FALSE}. Automatically set to \code{TRUE} if the table is to be saved in a file using the argument \code{file}.
#' @param treat_dict A character vector of length two. What are the names of the treated and the control? This should be a dictionary: e.g. \code{c("1"="Treated", "0" = "Control")}.
#' @param dict A named character vector. A dictionary between the variables names and an alias. For instance \code{dict=c("x"="Inflation Rate")} would replace the variable name \code{x} by \dQuote{Inflation Rate}.
#' @param file A file path. If given, the table is written in Latex into this file.
#' @param replace Default is \code{TRUE}, which means that when the table is exported, the existing file is not erased.
#' @param title Character string giving the Latex title of the table. (Only if exported.)
#' @param label Character string giving the Latex label of the table. (Only if exported.)
#' @param raw Logical, default is \code{FALSE}. If \code{TRUE}, it returns the information without formatting.
#' @param indiv Either the variable name of individual identifiers, a one sided formula, or a vector. If the data is that of a panel, this can be used to track the number of individuals per group.
#' @param prepostnames Only if there is a 'post' variable. The names of the pre and post periods to be displayed in Latex. Default is \code{c("Before", "After")}.
#' @param diff.inv Logical, default to \code{FALSE}. Whether to inverse the difference.
#'
#' @details
#' By default, when the user tries to apply this function to nun-numeric variables, an error is raised. The exception is when the all variables are selected with the dot (like in \code{. ~ treat}. In this case, non-numeric variables are automatically omitted (with a message).
#'
#' NAs are removed automatically: if the data contains NAs an information message will be prompted. First all observations containing NAs relating to the treatment or post variables are removed. Then if there are still NAs for the variables, they are excluded separately for each variable, and a new message detailing the NA breakup is prompted.
#'
#' @return
#' It returns a data.frame or a Latex table with the conditional means and statistical differences between the groups.
#'
#' @examples
#'
#' # Playing around with the DiD data
#' data(base_did)
#'
#' # means of treat/control
#' did_means(y+x1+period~treat, base_did)
#'
#' # same but inverting the difference
#' did_means(y+x1+period~treat, base_did, diff.inv = TRUE)
#'
#' # now treat/control, before/after
#' did_means(y+x1+period~treat|post, base_did)
#'
#' # same but with a new line giving the number of unique "indiv" for each case
#' did_means(y+x1+period~treat|post, base_did, indiv = "id")
#'
#' # same but with the treat case "0" coming first
#' did_means(y+x1+period~treat|post, base_did, indiv = ~id, treat_first = 0)
#'
#' # Selecting all the variables with "."
#' did_means(.~treat|post, base_did, indiv = "id")
#'
#'
did_means = function(fml, base, treat_var, post_var, tex = FALSE, treat_dict,
                     dict = getFixest_dict(), file, replace = FALSE, title,
                     label, raw = FALSE, indiv, treat_first, prepostnames = c("Before", "After"),
                     diff.inv = FALSE){
    # x is a data.frame
    # treat_vector is a list of IDs
    # treat_first: the treated-value to appear first

    # Ajouter protection quand un des groupes vaut 0!!!
    # ie un des treat-post n'a pas d'observation

    fml_in = fml

    #
    # Nber of units of observations (I do it before the main data)
    #

    IS_INDIV = FALSE
    if(!missing(indiv)){
        IS_INDIV = TRUE
        if("formula" %in% class(indiv)){

            check_arg(indiv, "os formula")
            indiv_varname = attr(terms(indiv), "term.labels")

            if(missing(base)){
                stop("To use 'indiv' as a formula, you must provide the argument 'base'.")
            }

            if(length(indiv_varname) > 1){
                stop("The argument 'indiv' must refer to only one variable.")
            }

            if(!all(all.vars(indiv) %in% names(base))){
                pblm = setdiff(all.vars(indiv), names(base))
                stop("In argument 'indiv': the variable", enumerate_items(pblm, "s.is")," not in the data set.")
            }

            indiv_var = try(eval(str2lang(indiv_varname), base), silent = TRUE)
            if("try-error" %in% class(indiv_var)){
                stop("Evaluation of 'indiv' raises and error:\n", indiv_var)
            }
        } else if(length(indiv) == 1 && is.character(indiv)){
            indiv_varname = indiv

            if(missing(base) || !indiv %in% names(base)){
                stop("To use 'indiv' as a 'character string' you must provide the argument 'base'.")
            } else {
                indiv_var = base[[indiv]]
            }
        } else {
            if(!missing(base) && length(indiv) != NROW(base)){
                stop("The length of 'indiv' must be the same as the data.")
            } else if(missing(base) && length(indiv) != NROW(fml)){
                stop("The length of 'indiv' must be the same as the data.")
            }
            indiv_var = indiv
        }
    } else {
        indiv_var = NULL
    }

    #
    # Extracting the information
    #

    usePost = FALSE
    if("formula" %in% class(fml_in)){
        if(missing(base) || !is.data.frame(base)){
            stop("If you provide a formula, a data.frame must be given in argument 'base'.")
        }

        # Creation of x and the condition
        if(!length(fml_in) == 3){
            stop("The formula must be of the type 'var1 + ... + varN ~ treat' or 'var1 + ... + varN ~ treat | post'.")
        }

        fml_parts = fml_split(fml_in)

        fml = fml_parts[[1]]
        pipe = NULL
        if(length(fml_parts) > 1){
            pipe = fml_parts[[2]]
        }

        #
        # Treat & post
        #

        if(!(all(all.vars(fml[[3]]) %in% names(base)))){
            pblm = setdiff(all.vars(fml[[3]]), names(base))
            stop("In the evaluation of the treatment variable: ", enumerate_items(pblm, "is"), " not in the data set.")
        }

        treat_var = try(eval(fml[[3]], base), silent = TRUE)
        if("try-error" %in% class(treat_var)){
            stop("Evaluation of the 'treatment' variable raises and error: \n", treat_var)
        }

        if(!is.null(pipe)){
            if(!(all(all.vars(pipe) %in% names(base)))){
                pblm = setdiff(all.vars(pipe), names(base))
                stop("In the evaluation of the 'post' variable: ", enumerate_items(pblm, "is"), " not in the data set.")
            }

            post_var = try(eval(pipe[[2]], base), silent = TRUE)
            if("try-error" %in% class(post_var)){
                stop("Evaluation of the 'post' variable raises and error: \n", treat_var)
            }

            usePost = TRUE
        }

        #
        # The variables
        #

        # special case: the point
        if(deparse(fml[[2]])[1] == "."){
            all_vars = setdiff(names(base), deparse(fml[[3]])[1])
            if(usePost) all_vars = setdiff(all_vars, as.character(pipe))
            if(IS_INDIV) all_vars = setdiff(all_vars, indiv_varname)

            if("data.table" %in% class(base)){
                mat_vars = as.data.frame(base[, all_vars, with = FALSE])
            } else {
                mat_vars = base[, all_vars, FALSE]
            }

            # we drop non-numeric info + note
            base_small = head(mat_vars, 10)
            is_num = sapply(base_small, function(x) is.numeric(x) || is.logical(x))
            pblm = all_vars[!is_num]

            if(all(!is_num)){
                stop("Function cannot be performed: not any numeric variable.")
            } else if(length(pblm) > 0){
                message("NOTE: The variable", enumerate_items(pblm, "s.is"), " removed because they are non-numeric.")
            }

            mat_vars = as.matrix(mat_vars[, is_num, FALSE])
        } else {
            # we swap the formula to use model.frame
            fml_x = as.formula(paste0("1~", deparse_long(fml[[2]]), "-1"))

            # Variables control:
            all_vars = all.vars(fml_x)
            pblm = setdiff(all_vars, names(base))
            if(length(pblm) > 0){
                stop("The variable", enumerate_items(pblm, "s.is"), " not in the data set (", deparse(substitute(base)), ").")
            }

            # Evaluation control
            base_small = head(base, 10)
            var2eval = attr(terms(fml_x), "term.labels")
            var2eval = gsub(":", "*", var2eval)
            for(i in seq_along(var2eval)){
                var = var2eval[i]
                x_small = try(eval(parse(text=var), base_small), silent = TRUE)
                if("try-error" %in% class(x_small)){
                    stop("Evaluation of the variable '", var, "' raises and error:\n", x_small)
                }

                if(!is.numeric(x_small) && !is.logical(x_small)){
                    stop("The variable ", var, " is not numeric, please provide only numeric variables.")
                }
            }

            mat_vars = prepare_matrix(fml_x, base)
        }

        if(is.logical(mat_vars)){
            mat_vars = as.numeric(mat_vars)
        }

        # other info
        x_name = colnames(mat_vars)
    } else {
        if(missing(treat_var)){
            stop("If argument 'fml' is not a formula, you must provide the argument 'treat_var'.")
        } else {
            mat_vars = fml_in

            if(NROW(mat_vars) != length(treat_var)){
                stop("The arguments 'x' and 'treat_var' must be of the same length.")
            }

            if(!missing(post_var) && !is.null(post_var)){
                if(NROW(mat_vars) != length(post_var)){
                    stop("If provided, the argument 'post_var' must be of the same length of 'x'.")
                }
                usePost = TRUE
            }

            # We exclude non numeric variables
            if("data.frame" %in% class(mat_vars)){
                is_num = sapply(mat_vars, function(x) is.numeric(x) || is.logical(x))
                if(any(!is_num)){
                    pblm = names(mat_vars)[!is_num]
                    stop("This function works only for numeric variables. Variable", enumerate_items(pblm, "s.is"), " not numeric.")
                }
                mat_vars = as.matrix(mat_vars)
            } else if(!is.matrix(mat_vars)){
                stop("If not a formula, argument 'fml' must be a data.frame or a matrix. Currently it is of class ", class(mat_vars)[1], ".")
            } else if(!is.numeric(mat_vars) && !is.logical(mat_vars)){
                stop("If not a formula, argument 'fml' must be a data.frame or a matrix with numeric variables. Currenlty its is not numeric.")
            }

            if(is.logical(mat_vars)) mat_vars = as.numeric(mat_vars)

            # other info
            x_name = colnames(mat_vars)
        }
    }

    #
    # NA control
    #

    # Treat & post
    ANY_NA = FALSE
    if(anyNA(treat_var) || (usePost && anyNA(post_var))){
        ANY_NA = TRUE
        qui_NA = is.na(treat_var)
        if(usePost) qui_NA = qui_NA | is.na(post_var)

        if(all(qui_NA)){
            msg = "treatment variable."
            if(usePost && anyNA(post_var)){
                if(anyNA(treat_var)){
                    msg = "'treatment' and 'post' variables."
                } else {
                    msg = "'post' variable."
                }
            }
            stop("All observation contain NA values for the ", msg)
        }

        mat_vars = mat_vars[!qui_NA, ]
        treat_var = treat_var[!qui_NA]
        if(usePost) post_var = post_var[!qui_NA]

        if(IS_INDIV) indiv_var = indiv_var[!qui_NA]

        msg = ifelse(usePost, "/post variables.", " variable.")
        message("NOTE: ", sum(qui_NA), " observations removed because of NA values in the treatment", msg)

    }



    # Treatment / post controls
    if(length(unique(treat_var)) != 2){
        n = length(unique(treat_var))
        msg = ifelse(n == 1, "only one value.", paste0(n, " values."))
        stop("This function supports only 2 conditional values for the treatment variable. Currently, it contains ", msg)
    }

    if(usePost && length(unique(post_var)) != 2){
        n = length(unique(post_var))
        msg = ifelse(n == 1, "only one value.", paste0(n, " values."))
        stop("This function supports only 2 conditional values for the 'post' variable. Currently, it contains ", msg)
    }

    if(usePost){
        check_arg(prepostnames, "character vector len(2) no na")
    }

    treat_cases = sort(unique(treat_var), decreasing = TRUE)
    if(!missing(treat_dict) && !is.null(treat_dict)){

        if(!isVector(treat_dict) || is.null(names(treat_dict))){
            stop("The argument 'treat_dict' must be a named character vector.")
        }

        pblm = setdiff(treat_cases, names(treat_dict))
        if(length(pblm) > 0){
            stop("The value", enumerate_items(pblm, "s.is"), " not in 'treat_dict'.")
        }
    } else {
        treat_dict = paste0("cond: ", treat_cases)
        names(treat_dict) = treat_cases
    }

    #
    # The main function
    #

    if(!missing(treat_first) && !is.null(treat_first)){
        if(!treat_first %in% treat_cases){
            stop("Argument 'treat_first' must be an element of the treated variable.")
        } else if(treat_first != treat_cases[1]){
            treat_cases = rev(treat_cases)
        }
    }

    if(!usePost){
        # the simple case
        res = .meanDiffTable(mat_vars, treat_var, treat_cases, indiv_var, raw = raw, diff.inv = diff.inv)

        if(any(attr(res, "na"))){
            nb_na = attr(res, "na")
            var_with_na = x_name[nb_na > 0]
            message("NA values were omitted for the variable", enumerate_items(paste0(var_with_na, " (", nb_na[nb_na > 0], ")"), "s"), ".")
        }

        attr(res, "na") = NULL


    } else {
        # the more complex case
        post_values = sort(unique(post_var))

        qui_pre = which(post_var == post_values[1])
        res_pre = .meanDiffTable(mat_vars[qui_pre, ], treat_var[qui_pre], treat_cases, indiv_var[qui_pre], raw = raw, diff.inv = diff.inv)

        qui_post = which(post_var == post_values[2])
        res_post = .meanDiffTable(mat_vars[qui_post, ], treat_var[qui_post], treat_cases, indiv_var[qui_post], raw = raw, diff.inv = diff.inv)

        nb_na = attr(res_pre, "na") + attr(res_post, "na")

        res = cbind(res_pre, res_post[, -1])

        if(any(nb_na)){
            var_with_na = x_name[nb_na > 0]
            message("NA values were omitted for the variable", enumerate_items(paste0(var_with_na, " (", nb_na[nb_na > 0], ")"), "s"), ".")
        }

        attr(res, "na") = NULL
    }

    #
    # Tex exportation
    #

    if(tex || !missing(file)){
        # we sink!
        if(!missing(file)){
            sink(file = file, append = (replace == FALSE))
            on.exit(sink())
        }

        if(!missing(title)){
            # if there is the title, then we create a "table" environment
            myTitle = title
            if(!missing(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
            cat("\\begin{table}[htbp]\\centering\n\\caption{",  myTitle, "}\n")
        }

        #
        # we write the tex code
        #

        treat_text = paste0(treat_dict[as.character(treat_cases)], collapse = " & ")

        if(!usePost){
            # The simple case
            cat("\\begin{tabular}{lcccc}\n")
            cat(" &  &  &  & \\tabularnewline\n\\hline\n\\hline\n")

            cat(" & ", treat_text,  "\\tabularnewline\n")

            # the variables "line"
            cat("Variables & mean (s.e.) & mean (s.e.) & Diff. & t-stat \\tabularnewline\n")
            cat("\\hline\n")
        } else {
            # The case with pre and post
            cat("\\begin{tabular}{lcccccccc}\n")
            cat(" &  &  &  &  &  &  &  & \\tabularnewline\n\\hline\n\\hline\n")

            # the treat post line
            cat(" & ", "\\multicolumn{4}{c}{", prepostnames[1], "} & \\multicolumn{4}{c}{", prepostnames[2], "}\\tabularnewline\n")

            # the 'treated' line
            cat(" & ", treat_text, " &  &  & ", treat_text, "\\tabularnewline\n")

            # the variables "line"
            cat("Variables & mean (s.e.) & mean (s.e.) & Diff. & t-stat & mean (s.e.) & mean (s.e.) & Diff. & t-stat \\tabularnewline\n")
            cat("\\hline\n")
        }


        # Changing the names of the variables
        if(!missnull(dict)){
            qui = which(res$vars %in% names(dict))
            for(i in qui) res$vars[i] = escape_latex(dict[res$vars[i]])
        }

        # The data
        for(i in 1:nrow(res)){
            cat(paste0(res[i,], collapse = " & "), "\\tabularnewline\n")
        }

        cat("\\hline\n\\hline\n &  &  &  & \\tabularnewline\n")
        cat("\\end{tabular}\n")

        if(!missing(title)) cat("\\end{table}\n")

        return(invisible(res))
    }

    res
}

.meanDiffTable = function(mat_vars, treat_var, treat_cases, indiv_var, raw, diff.inv = FALSE){
    # This function is the workhorse of did_means

    x_name = colnames(mat_vars)

    treat_01 = 1 * (treat_var == treat_cases[2])

    nthreads = getFixest_nthreads()
    info = cpppar_cond_means(mat_vars, treat_01, nthreads)

    means = info$means
    sds = info$sd
    ns = info$n
    n_01 = info$n_01

    name_1 = as.character(treat_cases[1])
    name_2 = as.character(treat_cases[2])

    means_1 = means[, 1]
    means_2 = means[, 2]

    sds_1 = sds[, 1]
    sds_2 = sds[, 2]

    n_1 = ns[, 1]
    n_2 = ns[, 2]

    format_1 = paste0(mysignif(means_1, 2), " (", mysignif(sds_1, 2), ")")
    format_2 = paste0(mysignif(means_2, 2), " (", mysignif(sds_2, 2), ")")

    if(diff.inv){
        D = (means_2 - means_1)
    } else {
        D = (means_1 - means_2)
    }

    sd_D = sqrt(sds_1**2/n_1 + sds_2**2/n_2)

    t_stat = signif(D/sd_D, 3)

    #
    # Formatting
    #

    if(raw){
        res = data.frame(vars = x_name, stringsAsFactors = FALSE)
        res[[paste0("Mean: ", treat_cases[1])]] = c(means_1, recursive = TRUE)
        res[[paste0("Mean: ", treat_cases[2])]] = c(means_2, recursive = TRUE)
        res[[paste0("se: ", treat_cases[1])]] = c(sds_1, recursive = TRUE)
        res[[paste0("se: ", treat_cases[2])]] = c(sds_2, recursive = TRUE)
        res[["Difference"]] = c(D, recursive=TRUE)
        res[["t-stat"]] = c(as.character(t_stat), recursive=TRUE)
    } else {
        res = data.frame(vars = x_name, stringsAsFactors = FALSE)
        res[[paste0("cond: ", treat_cases[1])]] = format_1
        res[[paste0("cond: ", treat_cases[2])]] = format_2
        res[["Difference"]] = c(mysignif(D, 3), recursive = TRUE)
        res[["t-stat"]] = c(as.character(t_stat), recursive = TRUE)

        res[nrow(res) + 1, ] = c("Observations", addCommas(n_01[1]), addCommas(n_01[2]), "", "")
        # Le nombre d'individus
        if(!missing(indiv_var) && !is.null(indiv_var)){
            nb_indiv =  tapply(indiv_var, treat_var, function(x) length(unique(x)))
            res[nrow(res) + 1, ] = c("# Individuals", addCommas(nb_indiv[name_1]), addCommas(nb_indiv[name_2]), "", "")
        }
    }

    attr(res, "na") = info$na

    res
}



#' Create, or interact variables with, factors
#'
#' Treat a variable as a factor, or interacts a variable with a factor. Values to be dropped/kept from the factor can be easily set. Note that to interact fixed-effects, this function should not be used: instead use directly the syntax \code{fe1^fe2}.
#'
#'
#' @inheritParams bin
#'
#' @param factor_var  A vector (of any type) that will be treated as a factor. You can set references (i.e. exclude values for which to create dummies) with the \code{ref} argument.
#' @param var A variable of the same length as \code{factor_var}. This variable will be interacted with the factor in \code{factor_var}. It can be numeric or factor-like. To force a numeric variable to be treated as a factor, you can add the \code{i.} prefix to a variable name. For instance take a numeric variable \code{x_num}: \code{i(x_fact, x_num)} will treat \code{x_num} as numeric while \code{i(x_fact, i.x_num)} will treat \code{x_num} as a factor (it's a shortcut to \code{as.factor(x_num)}).
#' @param ref A vector of values to be taken as references from \code{factor_var}. Can also be a logical: if \code{TRUE}, then the first value of \code{factor_var} will be removed. If \code{ref} is a character vector, partial matching is applied to values; use "@" as the first character to enable regular expression matching. See examples.
#' @param keep A vector of values to be kept from \code{factor_var} (all others are dropped). By default they should be values from \code{factor_var} and if \code{keep} is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param ref2 A vector of values to be dropped from \code{var}. By default they should be values from \code{var} and if \code{ref2} is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param keep2 A vector of values to be kept from \code{var} (all others are dropped). By default they should be values from \code{var} and if \code{keep2} is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param bin2 A list or vector defining the binning of the second variable. See help for the argument \code{bin} for details (or look at the help of the function \code{\link[fixest]{bin}}). You can use \code{.()} for \code{list()}.
#' @param ... Not currently used.
#'
#' @details
#' To interact fixed-effects, this function should not be used: instead use directly the syntax \code{fe1^fe2} in the fixed-effects part of the formula. Please see the details and examples in the help page of \code{\link[fixest]{feols}}.
#'
#' @return
#' It returns a matrix with number of rows the length of \code{factor_var}. If there is no interacted variable or it is interacted with a numeric variable, the number of columns is equal to the number of cases contained in \code{factor_var} minus the reference(s). If the interacted variable is a factor, the number of columns is the number of combined cases between \code{factor_var} and \code{var}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest:coefplot]{iplot}} to plot interactions or factors created with \code{i()}, \code{\link[fixest]{feols}} for OLS estimation with multiple fixed-effects.
#'
#' See the function \code{\link[fixest]{bin}} for binning variables.
#'
#' @examples
#'
#' #
#' # Simple illustration
#' #
#'
#' x = rep(letters[1:4], 3)[1:10]
#' y = rep(1:4, c(1, 2, 3, 4))
#'
#' # interaction
#' data.frame(x, y, i(x, y, ref = TRUE))
#'
#' # without interaction
#' data.frame(x, i(x, "b"))
#'
#' # you can interact factors too
#' z = rep(c("e", "f", "g"), c(5, 3, 2))
#' data.frame(x, z, i(x, z))
#'
#' # to force a numeric variable to be treated as a factor: use i.
#' data.frame(x, y, i(x, i.y))
#'
#' # Binning
#' data.frame(x, i(x, bin = list(ab = c("a", "b"))))
#'
#' # Same as before but using .() for list() and a regular expression
#' # note that to trigger a regex, you need to use an @ first
#' data.frame(x, i(x, bin = .(ab = "@a|b")))
#'
#' #
#' # In fixest estimations
#' #
#'
#' data(base_did)
#' # We interact the variable 'period' with the variable 'treat'
#' est_did = feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
#'
#' # => plot only interactions with iplot
#' iplot(est_did)
#'
#' # Using i() for factors
#' est_bis = feols(y ~ x1 + i(period, keep = 3:6) + i(period, treat, 5) | id, base_did)
#'
#' # we plot the second set of variables created with i()
#' # => we need to use keep (otherwise only the first one is represented)
#' coefplot(est_bis, keep = "trea")
#'
#' # => special treatment in etable
#' etable(est_bis, dict = c("6" = "six"))
#'
#' #
#' # Interact two factors
#' #
#'
#' # We use the i. prefix to consider week as a factor
#' data(airquality)
#' aq = airquality
#' aq$week = aq$Day %/% 7 + 1
#'
#' # Interacting Month and week:
#' res_2F = feols(Ozone ~ Solar.R + i(Month, i.week), aq)
#'
#' # Same but dropping the 5th Month and 1st week
#' res_2F_bis = feols(Ozone ~ Solar.R + i(Month, i.week, ref = 5, ref2 = 1), aq)
#'
#' etable(res_2F, res_2F_bis)
#'
#' #
#' # Binning
#' #
#'
#' data(airquality)
#'
#' feols(Ozone ~ i(Month, bin = "bin::2"), airquality)
#'
#' feols(Ozone ~ i(Month, bin = list(summer = 7:9)), airquality)
#'
#'
#'
i = function(factor_var, var, ref, keep, bin, ref2, keep2, bin2, ...){
    # Used to create interactions

    # Later: binning (bin = 1:3 // bin = list("a" = "[abc]")). Default name is bin name (eg "1:3")

    # gt = function(x) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s\n", sep = "")
    # t0 = proc.time()

    validate_dots(valid_args = c("f2", "f_name", "ref_special", "sparse"))

    dots = list(...)
    is_sparse = isTRUE(dots$sparse)

    mc = match.call()

    # Finding if it's a call from fixest
    FROM_FIXEST = is_fixest_call()

    # General checks
    check_arg(factor_var, "mbt vector")

    # NOTA:
    # the user can use the prefix "i." to tell the algorithm to consider the
    # variable var as a factor. This requires a non standard evaluation
    #

    var_name = NULL
    if(!is.null(mc$f2)){
        # should be only used internally

        var_name = deparse_long(mc$f2)
        var = dots$f2
        check_value(var, "vector", .arg_name = "f2")
    }

    # General information on the factor variable
    if(!is.null(dots$f_name)){
        f_name = dots$f_name
    } else {
        f_name = deparse_long(mc$factor_var)
    }
    f = factor_var # renaming for clarity

    # checking var
    IS_INTER_NUMERIC = IS_INTER_FACTOR = FALSE
    if(!missing(var)){
        # Evaluation of var (with possibly, user prepending with i.)
        if(is.null(var_name)){
            var_name = deparse_long(mc$var)
            if(grepl("^i\\.", var_name)){
                var_name = gsub("^i\\.", "", var_name)
                var = str2lang(var_name)
                check_value_plus(var, "evalset vector", .data = parent.frame(), .arg_name = "var")
                IS_INTER_FACTOR = TRUE
            } else {
                check_value(var, "vector", .arg_name = "var")
            }
        } else {
            # f2 was used
            IS_INTER_FACTOR = TRUE
        }

        # is it an interaction with a factor or with a continuous variable?
        if(length(var) == length(f)){
            # Logical => numeric by default
            # otherwise, just setting the flags

            if(IS_INTER_FACTOR == FALSE){
                # If previous condition == TRUE, means was requested by the user

                if(is.logical(var)){
                    var = as.numeric(var)
                    IS_INTER_NUMERIC = TRUE
                } else if(is.numeric(var)){
                    IS_INTER_NUMERIC = TRUE
                } else {
                    IS_INTER_FACTOR = TRUE
                }
            }

        } else if(length(var) <= 2 && !"ref" %in% names(mc)){
            # ref is implicitly called via the location of var
            ref = var
            dots$ref_special = FALSE
        } else {

            if(grepl("^[[:alpha:]\\.][[:alnum:]\\._]*:[[:alpha:]\\.][[:alnum:]\\._]*$", var_name)){
                info = strsplit(var_name, ":")[[1]]
                stop("In i(): When 'var' is equal to a product, please use I(", info[1], "*", info[2], ") instead of ", var_name, ".")
            } else {
                stop("The arguments 'var' and 'f' must be of the same length (currently ", length(var), " vs ", length(f), ").")
            }
        }
    }

    IS_INTER = IS_INTER_NUMERIC || IS_INTER_FACTOR

    if(isTRUE(dots$ref_special) && (!IS_INTER || IS_INTER_FACTOR)){
        ref = TRUE
    }

    #
    # QUFing + NA
    #

    if(IS_INTER){
        is_na_all = is.na(f) | is.na(var)
    } else {
        is_na_all = is.na(f)
    }

    if(!missing(bin)){
        bin = error_sender(eval_dot(bin), arg_name = "bin")

        if(!is.null(bin)){
            f = bin_factor(bin, f, f_name)
        }
    }

    if(IS_INTER_FACTOR && !MISSNULL(bin2)){
        bin2 = error_sender(eval_dot(bin2), arg_name = "bin2")

        if(!is.null(bin2)){
            var = bin_factor(bin2, var, f_name)
        }
    }

    if(IS_INTER_FACTOR){
        info = to_integer(f, var, add_items = TRUE, items.list = TRUE, sorted = TRUE, multi.join = "__%%__")
    } else {
        info = to_integer(f, add_items = TRUE, items.list = TRUE, sorted = TRUE)
    }

    fe_num = info$x
    items = info$items

    if(!IS_INTER_NUMERIC){
        # neutral var in C code
        var = 1
    }

    if(IS_INTER_FACTOR){
        items_split = strsplit(items, "__%%__", fixed = TRUE)

        f_items = sapply(items_split, `[`, 1)
        var_items = sapply(items_split, `[`, 2)
    } else {
        f_items = items
    }

    check_arg(ref, "logical scalar | vector no na")

    check_arg(ref2, keep, keep2, "vector no na")

    no_rm = TRUE
    id_drop = c()
    if(!missing(ref)){
        if(isTRUE(ref)){
            # We always delete the first value
            # Que ce soit items ici est normal (et pas f_items)
            id_drop = which(items == items[1])
        } else {
            id_drop = items_to_drop(f_items, ref, "factor_var")
        }
        ref_id = id_drop
    }


    if(!missing(keep)){
        id_drop = c(id_drop, items_to_drop(f_items, keep, "factor_var", keep = TRUE))
    }

    if(IS_INTER_FACTOR){
        if(!missing(ref2)){
            id_drop = c(id_drop, items_to_drop(var_items, ref2, "var"))
        }

        if(!missing(keep2)){
            id_drop = c(id_drop, items_to_drop(var_items, keep2, "var", keep = TRUE))
        }
    }

    if(length(id_drop) > 0){
        id_drop = unique(sort(id_drop))
        if(length(id_drop) == length(items)) stop("All items from the interaction have been removed.")
        who_is_dropped = id_drop
        no_rm = FALSE
    } else {
        # -1 is neutral
        who_is_dropped = -1
    }

    # The column names

    if(length(id_drop) > 0){
        items_name = items[-id_drop]
        f_items = f_items[-id_drop]
        if(IS_INTER_FACTOR){
            var_items = var_items[-id_drop]
        }
    } else {
        items_name = items
    }

    if(FROM_FIXEST){
        # Pour avoir des jolis noms c'est un vrai gloubiboulga,
        # mais j'ai pas trouve plus simple...
        if(IS_INTER_FACTOR){
            col_names = paste0("__CLEAN__", f_name, "::", f_items, ":", var_name, "::", var_items)
            col_names_full = NULL

        } else if(IS_INTER_NUMERIC){
            col_names = paste0("__CLEAN__", f_name, "::", f_items, ":", var_name)
            col_names_full = paste0(f_name, "::", items, ":", var_name)

        } else {
            col_names = paste0("__CLEAN__", f_name, "::", items_name)
            col_names_full = paste0(f_name, "::", items)
        }
    } else {

        if(IS_INTER_FACTOR){
            items_name = gsub("__%%__", ":", items_name, fixed = TRUE)
        }

        col_names = items_name
    }

    if(is_sparse){
        # Internal call: we return the row ids + the values + the indexes + the names

        if(length(who_is_dropped) > 0){
            valid_row = !is_na_all & !fe_num %in% who_is_dropped
        } else {
            valid_row = !is_na_all
        }

        # we need to ensure the IDs go from 1 to # Unique
        fe_colid = to_integer(fe_num[valid_row], sorted = TRUE)

        values = if(length(var) == 1) rep(1, length(valid_row)) else var
        res = list(rowid = which(valid_row), values = values,
                   colid = fe_colid, col_names = col_names)

        class(res) = "i_sparse"

        return(res)
    }

    res = cpp_factor_matrix(fe_num, is_na_all, who_is_dropped, var, col_names)
    # res => matrix with...
    #  - NAs where appropriate
    #  - appropriate number of columns
    #  - interacted if needed
    #


    # We send the information on the reference
    if(FROM_FIXEST){
        is_GLOBAL = FALSE
        for(where in 1:min(8, sys.nframe())){
            if(exists("GLOBAL_fixest_mm_info", parent.frame(where))){
                GLOBAL_fixest_mm_info = get("GLOBAL_fixest_mm_info", parent.frame(where))
                is_GLOBAL = TRUE
                break
            }
        }

        if(is_GLOBAL && !IS_INTER_FACTOR){
            # NOTA:
            # if you change stuff here => change them also in sunab()

            info = list()
            info$coef_names_full = col_names_full
            info$items = items
            if(!missing(ref)){
                info$ref_id = ref_id
                info$ref = items[ref_id]
            }
            info$is_num = is.numeric(items)
            info$f_name = f_name
            info$is_inter_num = IS_INTER_NUMERIC
            if(IS_INTER_NUMERIC){
                info$var_name = var_name
            }
            info$is_inter_fact = IS_INTER_FACTOR

            GLOBAL_fixest_mm_info[[length(GLOBAL_fixest_mm_info) + 1]] = info
            # re assignment
            assign("GLOBAL_fixest_mm_info", GLOBAL_fixest_mm_info, parent.frame(where))
        }

    }

    res
}


i_ref = function(factor_var, var, ref, bin, keep, ref2, keep2, bin2){
    # To automatically add references when i(x) is used

    mc = match.call()

    mc[[1]] = as.name("i")

    if(!any(c("ref", "keep", "ref2", "keep2") %in% names(mc))){
        mc$ref_special = TRUE
    }

    return(deparse_long(mc))
}

i_noref = function(factor_var, var, ref, bin, keep, ref2, keep2, bin2){
    # Used only in predict => to create data without restriction

    mc = match.call()

    mc[[1]] = as.name("i")

    mc$ref = mc$keep = mc$ref2 = mc$keep2 = NULL

    return(deparse_long(mc))
}


#' Bins the values of a variable (typically a factor)
#'
#' Tool to easily group the values of a given variable.
#'
#' @param x A vector whose values have to be grouped. Can be of any type but must be atomic.
#' @param bin A list of values to be grouped, a vector, a formula, or the special values \code{"bin::digit"} or \code{"cut::values"}. To create a new value from old values, use \code{bin = list("new_value"=old_values)} with \code{old_values} a vector of existing values. You can use \code{.()} for \code{list()}.
#' It accepts regular expressions, but they must start with an \code{"@"}, like in \code{bin="@Aug|Dec"}. It accepts one-sided formulas which must contain the variable \code{x}, e.g. \code{bin=list("<2" = ~x < 2)}.
#' The names of the list are the new names. If the new name is missing, the first value matched becomes the new name. In the name, adding \code{"@d"}, with \code{d} a digit, will relocate the value in position \code{d}: useful to change the position of factors. Use \code{"@"} as first item to make subsequent items be located first in the factor.
#' Feeding in a vector is like using a list without name and only a single element. If the vector is numeric, you can use the special value \code{"bin::digit"} to group every \code{digit} element.
#' For example if \code{x} represents years, using \code{bin="bin::2"} creates bins of two years.
#' With any data, using \code{"!bin::digit"} groups every digit consecutive values starting from the first value.
#' Using \code{"!!bin::digit"} is the same but starting from the last value.
#' With numeric vectors you can: a) use \code{"cut::n"} to cut the vector into \code{n} equal parts, b) use \code{"cut::a]b["} to create the following bins: \code{[min, a]}, \code{]a, b[}, \code{[b, max]}.
#' The latter syntax is a sequence of number/quartile (q0 to q4)/percentile (p0 to p100) followed by an open or closed square bracket. You can add custom bin names by adding them in the character vector after \code{'cut::values'}. See details and examples. Dot square bracket expansion (see \code{\link[fixest]{dsb}}) is enabled.
#'
#' @section "Cutting" a numeric vector:
#'
#' Numeric vectors can be cut easily into: a) equal parts, b) user-specified bins.
#'
#' Use \code{"cut::n"} to cut the vector into \code{n} (roughly) equal parts. Percentiles are used to partition the data, hence some data distributions can lead to create less than \code{n} parts (for example if P0 is the same as P50).
#'
#' The user can specify custom bins with the following syntax: \code{"cut::a]b]c]"etc}. Here the numbers \code{a}, \code{b}, \code{c}, etc, are a sequence of increasing numbers, each followed by an open or closed square bracket. The numbers can be specified as either plain numbers (e.g. \code{"cut::5]12[32["}), quartiles (e.g. \code{"cut::q1]q3["}), or percentiles (e.g. \code{"cut::p10]p15]p90]"}). Values of different types can be mixed: \code{"cut::5]q2[p80["} is valid provided the median (\code{q2}) is indeed greater than \code{5}, otherwise an error is thrown.
#'
#' The square bracket right of each number tells whether the numbers should be included or excluded from the current bin. For example, say \code{x} ranges from 0 to 100, then \code{"cut::5]"} will create two  bins: one from 0 to 5 and a second from 6 to 100. With \code{"cut::5["} the bins would have been 0-4 and 5-100.
#'
#' A factor is returned. The labels report the min and max values in each bin.
#'
#' To have user-specified bin labels, just add them in the character vector following \code{'cut::values'}. You don't need to provide all of them, and \code{NA} values fall back to the default label. For example, \code{bin = c("cut::4", "Q1", NA, "Q3")} will modify only the first and third label that will be displayed as \code{"Q1"} and \code{"Q3"}.
#'
#' @section \code{bin} vs \code{ref}:
#'
#' The functions \code{\link[fixest]{bin}} and \code{\link[fixest]{ref}} are able to do the same thing, then why use one instead of the other? Here are the differences:
#'
#' \itemize{
#' \item{}{\code{ref} always returns a factor. This is in contrast with \code{bin} which returns, when possible, a vector of the same type as the vector in input.}
#' \item{}{\code{ref} always places the values modified in the first place of the factor levels. On the other hand, \code{bin} tries to not modify the ordering of the levels. It is possible to make \code{bin} mimic the behavior of \code{ref} by adding an \code{"@"} as the first element of the list in the argument \code{bin}.}
#'  \item{}{when a vector (and not a list) is given in input, \code{ref} will place each element of the vector in the first place of the factor levels. The behavior of \code{bin} is totally different, \code{bin} will transform all the values in the vector into a single value in \code{x} (i.e. it's binning).}
#' }
#'
#' @return
#' It returns a vector of the same length as \code{x}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' To re-factor variables: \code{\link[fixest]{ref}}.
#'
#' @examples
#'
#' data(airquality)
#' month_num = airquality$Month
#' table(month_num)
#'
#' # Grouping the first two values
#' table(bin(month_num, 5:6))
#'
#' # ... plus changing the name to '10'
#' table(bin(month_num, list("10" = 5:6)))
#'
#' # ... and grouping 7 to 9
#' table(bin(month_num, list("g1" = 5:6, "g2" = 7:9)))
#'
#' # Grouping every two months
#' table(bin(month_num, "bin::2"))
#'
#' # ... every 2 consecutive elements
#' table(bin(month_num, "!bin::2"))
#'
#' # ... idem starting from the last one
#' table(bin(month_num, "!!bin::2"))
#'
#' # Using .() for list():
#' table(bin(month_num, .("g1" = 5:6)))
#'
#'
#' #
#' # with non numeric data
#' #
#'
#' month_lab = c("may", "june", "july", "august", "september")
#' month_fact = factor(month_num, labels = month_lab)
#'
#' # Grouping the first two elements
#' table(bin(month_fact, c("may", "jun")))
#'
#' # ... using regex
#' table(bin(month_fact, "@may|jun"))
#'
#' # ...changing the name
#' table(bin(month_fact, list("spring" = "@may|jun")))
#'
#' # Grouping every 2 consecutive months
#' table(bin(month_fact, "!bin::2"))
#'
#' # ...idem but starting from the last
#' table(bin(month_fact, "!!bin::2"))
#'
#' # Relocating the months using "@d" in the name
#' table(bin(month_fact, .("@5" = "may", "@1 summer" = "@aug|jul")))
#'
#' # Putting "@" as first item means subsequent items will be placed first
#' table(bin(month_fact, .("@", "aug", "july")))
#'
#' #
#' # "Cutting" numeric data
#' #
#'
#' data(iris)
#' plen = iris$Petal.Length
#'
#' # 3 parts of (roughly) equal size
#' table(bin(plen, "cut::3"))
#'
#' # Three custom bins
#' table(bin(plen, "cut::2]5]"))
#'
#' # .. same, excluding 5 in the 2nd bin
#' table(bin(plen, "cut::2]5["))
#'
#' # Using quartiles
#' table(bin(plen, "cut::q1]q2]q3]"))
#'
#' # Using percentiles
#' table(bin(plen, "cut::p20]p50]p70]p90]"))
#'
#' # Mixing all
#' table(bin(plen, "cut::2[q2]p90]"))
#'
#' # NOTA:
#' # -> the labels always contain the min/max values in each bin
#'
#' # Custom labels can be provided, just give them in the char. vector
#' # NA values lead to the default label
#' table(bin(plen, c("cut::2[q2]p90]", "<2", "]2; Q2]", NA, ">90%")))
#'
#'
#'
#' #
#' # With a formula
#' #
#'
#' data(iris)
#' plen = iris$Petal.Length
#'
#' # We need to use "x"
#' table(bin(plen, list("< 2" = ~x < 2, ">= 2" = ~x >= 2)))
#'
#'
bin = function(x, bin){
    check_arg(x, "vector mbt")

    bin = error_sender(eval_dot(bin), arg_name = "bin")

    check_arg(bin, "list | vector mbt")

    varname = deparse(substitute(x))[1]
    bin_factor(bin, x, varname)
}

#' Refactors a variable
#'
#' Takes a variables of any types, transforms it into a factors, and modifies the values of the factors. Useful in estimations when you want to set some value of a vector as a reference.
#'
#' @inheritSection bin \code{bin} vs \code{ref}
#'
#' @param x A vector of any type (must be atomic though).
#' @param ref A vector or a list, or special binning values (explained later). If a vector, it must correspond to (partially matched) values of the vector \code{x}. The vector \code{x} which will be transformed into a factor and these values will be placed first in the levels. That's the main usage of this function. You can also bin on-the-fly the values of \code{x}, using the same syntax as the function \code{\link[fixest]{bin}}. Here's a description of what bin does: To create a new value from old values, use \code{bin = list("new_value"=old_values)} with \code{old_values} a vector of existing values. You can use \code{.()} for \code{list()}.
#' It accepts regular expressions, but they must start with an \code{"@"}, like in \code{bin="@Aug|Dec"}. It accepts one-sided formulas which must contain the variable \code{x}, e.g. \code{bin=list("<2" = ~x < 2)}.
#' The names of the list are the new names. If the new name is missing, the first value matched becomes the new name. In the name, adding \code{"@d"}, with \code{d} a digit, will relocate the value in position \code{d}: useful to change the position of factors.
#' If the vector \code{x} is numeric, you can use the special value \code{"bin::digit"} to group every \code{digit} element.
#' For example if \code{x} represents years, using \code{bin="bin::2"} creates bins of two years.
#' With any data, using \code{"!bin::digit"} groups every digit consecutive values starting from the first value.
#' Using \code{"!!bin::digit"} is the same but starting from the last value.
#' With numeric vectors you can: a) use \code{"cut::n"} to cut the vector into \code{n} equal parts, b) use \code{"cut::a]b["} to create the following bins: \code{[min, a]}, \code{]a, b[}, \code{[b, max]}.
#' The latter syntax is a sequence of number/quartile (q0 to q4)/percentile (p0 to p100) followed by an open or closed square bracket. You can add custom bin names by adding them in the character vector after \code{'cut::values'}. See details and examples. Dot square bracket expansion (see \code{\link[fixest]{dsb}}) is enabled.
#'
#' @return
#' It returns a factor of the same length as \code{x}, where levels have been modified according to the argument \code{ref}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' To bin the values of a vectors: \code{\link[fixest]{bin}}.
#'
#' @examples
#'
#' data(airquality)
#'
#' # A vector of months
#' month_num = airquality$Month
#' month_lab = c("may", "june", "july", "august", "september")
#' month_fact = factor(month_num, labels = month_lab)
#' table(month_num)
#' table(month_fact)
#'
#' #
#' # Main use
#' #
#'
#' # Without argument: equivalent to as.factor
#' ref(month_num)
#'
#' # Main usage: to set a level first:
#' # (Note that partial matching is enabled.)
#' table(ref(month_fact, "aug"))
#'
#' # You can rename the level on-the-fly
#' # (Northern hemisphere specific!)
#' table(ref(month_fact, .("Hot month"="aug",
#'                         "Late summer" = "sept")))
#'
#'
#' # Main use is in estimations:
#' a = feols(Petal.Width ~ Petal.Length + Species, iris)
#'
#' # We change the reference
#' b = feols(Petal.Width ~ Petal.Length + ref(Species, "vers"), iris)
#'
#' etable(a, b)
#'
#'
#' #
#' # Binning
#' #
#'
#' # You can also bin factor values on the fly
#' # Using @ first means a regular expression will be used to match the values.
#' # Note that the value created is placed first.
#' # To avoid that behavior => use the function "bin"
#' table(ref(month_fact, .(summer = "@jul|aug|sep")))
#'
#' # Please refer to the example in the bin help page for more example.
#' # The syntax is the same.
#'
#'
#' #
#' # Precise relocation
#' #
#'
#' # You can place a factor at the location you want
#' #  by adding "@digit" in the name first:
#' table(ref(month_num, .("@5"=5)))
#'
#' # Same with renaming
#' table(ref(month_num, .("@5 five"=5)))
#'
#'
ref = function(x, ref){
    check_arg(x, "vector mbt")

    if(missing(ref)){
        return(as.factor(x))
    }

    ref = error_sender(eval_dot(ref), arg_name = "ref")
    check_arg(ref, "list | vector mbt")

    varname = deparse(substitute(x))[1]

    IS_SPECIAL = FALSE
    if(!is.list(ref)){
        if(is.character(ref[1]) && grepl("^(cut|bin)", ref[1])){
            IS_SPECIAL = TRUE
            if(!is.numeric(x)){
                stop(.dsb("To use the special binning '.[ref[1]]' the variable ",
                          "'.[varname]' must be numeric. Currently this is not the case ",
                          "(it is of class .[3KO, C?class(x)] instead)."))
            }
        } else {
            ref = as.list(ref)
        }
    }

    if(!IS_SPECIAL && !is.factor(x)){
        x = as.factor(x)
    }

    if(!IS_SPECIAL && ref[[1]] != "@"){
        ref_new = list("@")
        index = 1:length(ref) + 1
        ref_new[index] = ref
        names(ref_new)[index] = names(ref)
        ref = ref_new
    }

    bin_factor(ref, x, varname)
}



#' Expands formula macros
#'
#' Create macros within formulas and expand them with character vectors or other formulas.
#'
#' @inheritParams setFixest_fml
#'
#' @param fml A formula containing macros variables. Each macro variable must start with two dots. The macro variables can be set globally using \code{setFixest_fml}, or can be defined in \code{...}. Special macros of the form \code{..("regex")} can be used to fetch, through a regular expression, variables directly in a character vector (or in column names) given in the argument \code{data}. square brackets have a special meaning: Values in them are evaluated and parsed accordingly. Example: \code{y~x.[1:2] + z.[i]} will lead to \code{y~x1+x2+z3} if \code{i==3}. See examples.
#' @param add Either a character scalar or a one-sided formula. The elements will be added to the right-hand-side of the formula, before any macro expansion is applied.
#' @param lhs If present then a formula will be constructed with \code{lhs} as the full left-hand-side. The value of \code{lhs} can be a one-sided formula, a call, or a character vector. Note that the macro variables wont be applied. You can use it in combination with the argument \code{rhs}. Note that if \code{fml} is not missing, its LHS will be replaced by \code{lhs}.
#' @param rhs If present, then a formula will be constructed with \code{rhs} as the full right-hand-side. The value of \code{rhs} can be a one-sided formula, a call, or a character vector. Note that the macro variables wont be applied. You can use it in combination with the argument \code{lhs}. Note that if \code{fml} is not missing, its RHS will be replaced by \code{rhs}.
#' @param data Either a character vector or a data.frame. This argument will only be used if a macro of the type \code{..("regex")} is used in the formula of the argument \code{fml}. If so, any variable name from \code{data} that matches the regular expression will be added to the formula.
#'
#' @details
#' In \code{xpd}, the default macro variables are taken from \code{getFixest_fml}. Any value in the \code{...} argument of \code{xpd} will replace these default values.
#'
#' The definitions of the macro variables will replace in verbatim the macro variables. Therefore, you can include multi-part formulas if you wish but then beware of the order of the macros variable in the formula. For example, using the \code{airquality} data, say you want to set as controls the variable \code{Temp} and \code{Day} fixed-effects, you can do \code{setFixest_fml(..ctrl = ~Temp | Day)}, but then \code{feols(Ozone ~ Wind + ..ctrl, airquality)} will be quite different from \code{feols(Ozone ~ ..ctrl + Wind, airquality)}, so beware!
#'
#' @section Dot square bracket operator in formulas:
#'
#' In a formula, the dot square bracket (DSB) operator can: i) create manifold variables at once, or ii) capture values from the current environment and put them verbatim in the formula.
#'
#' Say you want to include the variables \code{x1} to \code{x3} in your formula. You can use \code{xpd(y ~ x.[1:3])} and you'll get \code{y ~ x1 + x2 + x3}.
#'
#' To summon values from the environment, simply put the variable in square brackets. For example: \code{for(i in 1:3) xpd(y.[i] ~ x)} will create the formulas \code{y1 ~ x} to \code{y3 ~ x} depending on the value of \code{i}.
#'
#' You can include a full variable from the environment in the same way: \code{for(y in c("a", "b")) xpd(.[y] ~ x)} will create the two formulas \code{a ~ x} and \code{b ~ x}.
#'
#' The DSB can even be used within variable names, but then the variable must be nested in character form. For example \code{y ~ .["x.[1:2]_sq"]} will create \code{y ~ x1_sq + x2_sq}. Using the character form is important to avoid a formula parsing error. Double quotes must be used. Note that the character string that is nested will be parsed with the function \code{\link[fixest]{dsb}}, and thus it will return a vector.
#'
#' By default, the DSB operator expands vectors into sums. You can add a comma, like in \code{.[, x]}, to expand with commas--the content can then be used within functions. For instance: \code{c(x.[, 1:2])} will create \code{c(x1, x2)} (and \emph{not} \code{c(x1 + x2)}).
#'
#' In all \code{fixest} estimations, this special parsing is enabled, so you don't need to use \code{xpd}.
#'
#' You can even use multiple square brackets within a single variable, but then the use of nesting is required. For example, the following \code{xpd(y ~ .[".[letters[1:2]]_.[1:2]"])} will create \code{y ~ a_1 + b_2}. Remember that the nested character string is parsed with \code{\link[fixest]{dsb}}, which explains this behavior.
#'
#' @section Regular expressions:
#'
#' You can catch several variable names at once by using regular expressions. To use regular expressions, you need to enclose it in the dot-dot or the regex function: \code{..("regex")} or \code{regex("regex")}. For example, \code{regex("Sepal")} will catch both the variables \code{Sepal.Length} and \code{Sepal.Width} from the \code{iris} data set. In a \code{fixest} estimation, the variables names from which the regex will be applied come from the data set. If you use \code{xpd}, you need to provide either a data set or a vector of names in the argument \code{data}.
#'
#' By default the variables are aggregated with a sum. For example in a data set with the variables x1 to x10, \code{regex("x(1|2)"} will yield \code{x1 + x2 + x10}. You can instead ask for "comma" aggregation by using a comma first, just before the regular expression: \code{y ~ sw(regex(,"x(1|2)"))} would lead to \code{y ~ sw(x1, x2, x10)}.
#'
#' Note that the dot square bracket operator (DSB, see before) is applied before the regular expression is evaluated. This means that \code{regex("x.[3:4]_sq")} will lead, after evaluation of the DSB, to \code{regex("x3_sq|x4_sq")}. It is a handy way to insert range of numbers in a regular expression.
#'
#'
#' @return
#' It returns a formula where all macros have been expanded.
#'
#' @author
#' Laurent Berge
#'
#'
#' @seealso
#' \code{\link[fixest]{setFixest_fml}} to set formula macros, and \code{\link[fixest]{dsb}} to modify character strings with the DSB operator.
#'
#' @examples
#'
#' # Small examples with airquality data
#' data(airquality)
#' # we set two macro variables
#' setFixest_fml(..ctrl = ~ Temp + Day,
#'               ..ctrl_long = ~ poly(Temp, 2) + poly(Day, 2))
#'
#' # Using the macro in lm with xpd:
#' lm(xpd(Ozone ~ Wind + ..ctrl), airquality)
#' lm(xpd(Ozone ~ Wind + ..ctrl_long), airquality)
#'
#' # You can use the macros without xpd() in fixest estimations
#' a = feols(Ozone ~ Wind + ..ctrl, airquality)
#' b = feols(Ozone ~ Wind + ..ctrl_long, airquality)
#' etable(a, b, keep = "Int|Win")
#'
#'
#' # Using .[]
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' i = 2:3
#' z = "species"
#' lm(xpd(y ~ x.[2:3] + .[z]), base)
#'
#' # No xpd() needed in feols
#' feols(y ~ x.[2:3] + .[z], base)
#'
#' #
#' # You can use xpd for stepwise estimations
#' #
#'
#' # Note that for stepwise estimations in fixest, you can use
#' # the stepwise functions: sw, sw0, csw, csw0
#' # -> see help in feols or in the dedicated vignette
#'
#' # we want to look at the effect of x1 on y
#' # controlling for different variables
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # We first create a matrix with all possible combinations of variables
#' my_args = lapply(names(base)[-(1:2)], function(x) c("", x))
#' (all_combs = as.matrix(do.call("expand.grid", my_args)))
#'
#' res_all = list()
#' for(i in 1:nrow(all_combs)){
#'   res_all[[i]] = feols(xpd(y ~ x1 + ..v, ..v = all_combs[i, ]), base)
#' }
#'
#' etable(res_all)
#' coefplot(res_all, group = list(Species = "^^species"))
#'
#' #
#' # You can use macros to grep variables in your data set
#' #
#'
#' # Example 1: setting a macro variable globally
#'
#' data(longley)
#' setFixest_fml(..many_vars = grep("GNP|ployed", names(longley), value = TRUE))
#' feols(Armed.Forces ~ Population + ..many_vars, longley)
#'
#' # Example 2: using ..("regex") or regex("regex") to grep the variables "live"
#'
#' feols(Armed.Forces ~ Population + ..("GNP|ployed"), longley)
#'
#' # Example 3: same as Ex.2 but without using a fixest estimation
#'
#' # Here we need to use xpd():
#' lm(xpd(Armed.Forces ~ Population + regex("GNP|ployed"), data = longley), longley)
#'
#' # Stepwise estimation with regex: use a comma after the parenthesis
#' feols(Armed.Forces ~ Population + sw(regex(,"GNP|ployed")), longley)
#'
#' # Multiple LHS
#' etable(feols(..("GNP|ployed") ~ Population, longley))
#'
#'
#' #
#' # lhs and rhs arguments
#' #
#'
#' # to create a one sided formula from a character vector
#' vars = letters[1:5]
#' xpd(rhs = vars)
#'
#' # Alternatively, to replace the RHS
#' xpd(y ~ 1, rhs = vars)
#'
#' # To create a two sided formula
#' xpd(lhs = "y", rhs = vars)
#'
#' #
#' # argument 'add'
#' #
#'
#' xpd(~x1, add = ~ x2 + x3)
#'
#' # also works with character vectors
#' xpd(~x1, add = c("x2", "x3"))
#'
#' # only adds to the RHS
#' xpd(y ~ x, add = ~bon + jour)
#'
#' #
#' # Dot square bracket operator
#' #
#'
#' # You can create multiple variables at once
#' xpd(y ~ x.[1:5] + z.[2:3])
#'
#' # You can summon variables from the environment
#' var = "a"
#' xpd(y ~ x.[var])
#'
#' # ... the variables can be multiple
#' vars = LETTERS[1:3]
#' xpd(y ~ x.[vars])
#'
#' # You can have "complex" variable names but they must be nested in character form
#' xpd(y ~ .["x.[vars]_sq"])
#'
#' # DSB can be used within regular expressions
#' re = c("GNP", "Pop")
#' xpd(Unemployed ~ regex(".[re]"), data = longley)
#'
#' # => equivalent to regex("GNP|Pop")
#'
#' # Use .[,var] (NOTE THE COMMA!) to expand with commas
#' # !! can break the formula if missused
#' vars = c("wage", "unemp")
#' xpd(c(y.[,1:3]) ~ csw(.[,vars]))
#'
#'
#' # Example of use of .[] within a loop
#' res_all = list()
#' for(p in 1:3){
#'   res_all[[p]] = feols(Ozone ~ Wind + poly(Temp, .[p]), airquality)
#' }
#'
#' etable(res_all)
#'
#' # The former can be compactly estimated with:
#' res_compact = feols(Ozone ~ Wind + sw(.[, "poly(Temp, .[1:3])"]), airquality)
#'
#' etable(res_compact)
#'
#' # How does it work?
#' # 1)  .[, stuff] evaluates stuff and, if a vector, aggregates it with commas
#' #     Comma aggregation is done thanks to the comma placed after the square bracket
#' #     If .[stuff], then aggregation is with sums.
#' # 2) stuff is evaluated, and if it is a character string, it is evaluated with
#' # the function dsb which expands values in .[]
#' #
#' # Wrapping up:
#' # 2) evaluation of dsb("poly(Temp, .[1:3])") leads to the vector:
#' #    c("poly(Temp, 1)", "poly(Temp, 2)", "poly(Temp, 3)")
#' # 1) .[, c("poly(Temp, 1)", "poly(Temp, 2)", "poly(Temp, 3)")] leads to
#' #    poly(Temp, 1), poly(Temp, 2), poly(Temp, 3)
#' #
#' # Hence sw(.[, "poly(Temp, .[1:3])"]) becomes:
#' #       sw(poly(Temp, 1), poly(Temp, 2), poly(Temp, 3))
#'
#'
xpd = function(fml, ..., add = NULL, lhs, rhs, data = NULL){
    .xpd(fml = fml, ..., add = add, lhs = lhs, rhs = rhs, data = data, check = TRUE,
         macro = TRUE, frame = parent.frame())
}

.xpd = function(fml, ..., add = NULL, lhs, rhs, data = NULL, check = FALSE, macro = FALSE, frame = NULL){

    is_lhs = !missing(lhs)
    is_rhs = !missing(rhs)
    if(is_lhs || is_rhs){

        if(check) check_arg(fml, .type = "formula", .up = 1)

        # Direct formula creation
        if(missing(fml)){
            res = if(is_lhs) 1 ~ 1 else ~ 1
        } else {
            if(is_lhs && length(fml)[1] == 2){
                # Means: RHS formula and we add the LHS
                res = 1 ~ 1
                res[[3]] = fml[[2]]

            } else {
                res = fml
            }
        }

        if(is_lhs){
            res[[2]] = value2stringCall(lhs, call = TRUE, check = check)
        }

        if(is_rhs){
            res[[length(res)]] = value2stringCall(rhs, call = TRUE, check = check)
        }

        fml = res

        if(!macro && missnull(add)) return(fml)

        # NOTA:
        # if we allow for macro implementation ex post:
        # This entails a 50% performance drop in terms of speed.
        # Now, without macro variables, speed is at 30us while it was 20us before
        # so in .xpd => macro argument

    } else if(!missing(add)){

        if(check){
            check_arg(fml, .type = "NULL formula", .up = 1)
        }

        if(missnull(fml)){
            fml = ~ 1
            fml[[2]] = value2stringCall(add, call = TRUE, check = check)
            add = NULL
        }

    } else if(check){
        check_arg(fml, .type = "formula mbt", .up = 1)
    }


    if(!missnull(add)){
        # Direct formula manipulation is too complicated (and I want to avoid ugly parentheses)
        # by string it's easy
        fml_dp = deparse_long(fml)

        add_txt = value2stringCall(add, call = FALSE, check = check)
        add_txt = gsub("^~", "", add_txt)

        fml = as.formula(paste0(fml_dp, "+", add_txt), frame)
    }

    macros = parse_macros(..., from_xpd = TRUE, check = check, frame = frame)

    is_macro = length(macros) != 0
    is_data = !missnull(data)
    fml_funs = all.vars(fml, functions = TRUE)
    is_brackets = "[" %in% fml_funs
    is_regex = any(c("regex", "..") %in% fml_funs)
    if(!(is_macro || is_data || is_brackets || is_regex)) return(fml)

    fml_dp = NULL

    if(is_macro){
        # We allow recursion + auto-blocking at 5th depth
        # to allow crazy things from the user side (but not too much)
        max_depth = 5
        depth  = 0
        any_done = TRUE
        qui = which(names(macros) %in% all.vars(fml))
        while(length(qui) > 0 && any_done && depth < max_depth){
            depth = depth + 1
            fml_dp_origin = fml_dp = deparse_long(fml)
            # We need to add a lookafter assertion: otherwise if we have ..ctrl + ..ctrl_long, there will be a replacement in ..ctrl_long
            for(i in qui){
                fml_dp = gsub(paste0(escape_regex(names(macros)[i]), "(?=$|[^[:alnum:]_\\.])"), macros[[i]], fml_dp, perl = TRUE)
            }

            fml = as.formula(fml_dp, frame)

            qui = which(names(macros) %in% all.vars(fml))
            any_done = fml_dp_origin != fml_dp
        }

        if(depth == max_depth && length(qui) > 0){
            warning(dsb("In xpd, max recursivity has been reached. ",
                        "Please revise the recursive definition of your variables. ",
                        "It concerns the variable.[*s_, q, C?names(macros)[qui]]."))
        }
    }

    if(is_regex){
        # We expand only if data is provided (it means the user wants us to check)
        # if .[]: we expand inside the ..(".[1:3]") => ..("1|2|3")

        if(is.null(fml_dp)) fml_dp = deparse_long(fml)
        is_brackets = grepl(".[", fml_dp, fixed = TRUE)

        if(is_data || is_brackets){

            if(is.null(fml_dp)) fml_dp = deparse_long(fml)

            if(is_data){
                check_arg(data, "character vector no na | matrix | data.frame")

                if(is.matrix(data)){
                    data = colnames(data)
                    if(is.null(data)){
                        stop("The argument 'data' must contain variable names. It is currently a matrix without column names.")
                    }
                } else if(is.data.frame(data)){
                    data = names(data)
                }
            }

            fml_dp_split = strsplit(fml_dp, '(?!<[[:alnum:]._])(regex|\\.\\.)\\((?=[,"])',
                                    perl = TRUE)[[1]]

            res = fml_dp_split
            for(i in 2:length(res)){

                re = sub('"\\).*', "", res[i])

                re_width = nchar(re)

                is_comma = grepl("^,", re)
                re = gsub("^,? ?\"", "", re)

                re = dot_square_bracket(re, frame, regex = TRUE)

                if(is_data){
                    vars = grep(re, data, value = TRUE, perl = TRUE)
                    if(length(vars) == 0){
                        vars = "1"
                    }

                    coll = if(is_comma) ", " else " + "

                    res[i] = paste0(paste(vars, collapse = coll), substr(res[i], re_width + 3, nchar(res[i])))
                } else {
                    res[i] = paste0("..(\"", re, "\")", substr(res[i], re_width + 3, nchar(res[i])))
                }

            }

            fml_txt = paste(res, collapse = "")
            fml = error_sender(as.formula(fml_txt, frame),
                               "Expansion of variables in ..(\"regex\"): coercion of the following text to a formula led to an error.\n",
                               fit_screen(paste0("     TEXT: ", fml_txt)),
                               "\n  PROBLEM: see below", up = 1)
        }
    }

    if("[" %in% all.vars(fml, functions = TRUE)){
        fml_txt = deparse_long(fml)
        if(grepl(".[", fml_txt, fixed = TRUE)){
            fml_txt = dot_square_bracket(fml_txt, frame)

            # level 1 recursivity
            if(grepl(".[", fml_txt, fixed = TRUE)){
                fml_txt = dot_square_bracket(fml_txt, frame)
            }

            fml = error_sender(as.formula(fml_txt, frame),
                               "Dot square bracket operator: coercion of the following text to a formula led to an error.\n",
                               fit_screen(paste0("     TEXT: ", fml_txt)),
                               "\n  PROBLEM: see below", up = 1)
        }

    }

    fml
}


#' Fast transform of any type of vector(s) into an integer vector
#'
#' Tool to transform any type of vector, or even combination of vectors, into an integer vector ranging from 1 to the number of unique values. This actually creates an unique identifier vector.
#'
#' @param ... Vectors of any type, to be transformed in integer.
#' @param sorted Logical, default is \code{FALSE}. Whether the integer vector should make reference to sorted values?
#' @param add_items Logical, default is \code{FALSE}. Whether to add the unique values of the original vector(s). If requested, an attribute \code{items} is created containing the values (alternatively, they can appear in a list if \code{items.list=TRUE}).
#' @param items.list Logical, default is \code{FALSE}. Only used if \code{add_items=TRUE}. If \code{TRUE}, then a list of length 2 is returned with \code{x} the integer vector and \code{items} the vector of items.
#' @param multi.df Logical, default is \code{FALSE}. If \code{TRUE} then a data.frame listing the unique elements is returned in the form of a data.frame. Ignored if \code{add_items = FALSE}.
#' @param multi.join Character scalar used to join the items of multiple vectors. The default is \code{"_"}. Ignored if \code{add_items = FALSE}.
#' @param internal Logical, default is \code{FALSE}. For programming only. If this function is used within another function, setting \code{internal = TRUE} is needed to make the evaluation of \code{...} valid. End users of \code{to_integer} should not care.
#'
#'
#' @return
#' Reruns a vector of the same length as the input vectors.
#' If \code{add_items=TRUE} and \code{items.list=TRUE}, a list of two elements is returned: \code{x} being the integer vector and \code{items} being the unique values to which the values in \code{x} make reference.
#'
#' @examples
#'
#' x1 = iris$Species
#' x2 = as.integer(iris$Sepal.Length)
#'
#' # transforms the species vector into integers
#' to_integer(x1)
#'
#' # To obtain the "items":
#' to_integer(x1, add_items = TRUE)
#' # same but in list form
#' to_integer(x1, add_items = TRUE, items.list = TRUE)
#'
#' # transforms x2 into an integer vector from 1 to 4
#' to_integer(x2, add_items = TRUE)
#'
#' # To have the sorted items:
#' to_integer(x2, add_items = TRUE, sorted = TRUE)
#'
#' # The result can safely be used as an index
#' res = to_integer(x2, add_items = TRUE, sorted = TRUE, items.list = TRUE)
#' all(res$items[res$x] == x2)
#'
#'
#' #
#' # Multiple vectors
#' #
#'
#' to_integer(x1, x2, add_items = TRUE)
#'
#' # You can use multi.join to handle the join of the items:
#' to_integer(x1, x2, add_items = TRUE, multi.join = "; ")
#'
to_integer = function(..., sorted = FALSE, add_items = FALSE, items.list = FALSE,
                      multi.df = FALSE, multi.join = "_", internal = FALSE){

    if(!internal) check_arg(..., "vector mbt")
    check_arg(sorted, add_items, items.list, "logical scalar")
    check_arg(multi.join, "character scalar")

    dots = list(...)

    # Removing NAs
    Q = length(dots)
    n_all = lengths(dots)
    n = n_all[1]

    if(length(unique(n_all)) != 1) stop("All elements in '...' should be of the same length (current lenghts are ", enumerate_items(n_all), ").")

    is_na = is.na(dots[[1]])
    for(q in seq(from = 2, length.out = Q - 1)){
        is_na = is_na | is.na(dots[[q]])
    }

    ANY_NA = FALSE
    if(any(is_na)){
        ANY_NA = TRUE

        if(all(is_na)){
            message("NOTE: All values are NA.")
            res = rep(NA, n)
            if(add_items){
                if(items.list){
                    res = list(x = res, items = NA)
                } else {
                    attr(res, "items") = NA
                }
            }

            return(res)
        }

        for(q in 1:Q) dots[[q]] = dots[[q]][!is_na]
    }

    #
    # Creating the ID
    #

    if(Q == 1){
        if(sorted && !is.numeric(dots[[1]]) && !is.character(dots[[1]])){
            # general way => works for any type with a sort method
            f = dots[[1]]
            res_raw = quickUnclassFactor(f, addItem = TRUE, sorted = FALSE)
            obs_1st = cpp_get_first_item(res_raw$x, length(res_raw$items))
            f_unik = f[obs_1st]
            f_order = order(f_unik)
            x_new = order(f_order)[res_raw$x]
            if(add_items){
                items_new = f_unik[f_order]
                res = list(x = x_new, items = items_new)
            } else {
                res = x_new
            }

        } else {
            res = quickUnclassFactor(dots[[1]], addItem = add_items, sorted = sorted)
        }

    } else {

        QUF_raw = list()
        for(q in 1:Q){
            QUF_raw[[q]] = quickUnclassFactor(dots[[q]], sorted = FALSE, addItem = TRUE)
        }

        # Then we combine
        power = floor(1 + log10(sapply(QUF_raw, function(x) length(x$items))))

        is_large = sum(power) > 14
        if(is_large){
            # 15 Aug 2021, finally found a solution. It was so obvious with hindsight...
            QUF_raw_value = lapply(QUF_raw, `[[`, 1)
            order_index = do.call(order, QUF_raw_value)
            index = cpp_combine_clusters(QUF_raw_value, order_index)
        } else {
            # quicker, but limited by the precision of doubles
            index = QUF_raw[[1]]$x
            for(q in 2:Q){
                index = index + QUF_raw[[q]]$x*10**sum(power[1:(q-1)])
            }
        }

        res = quickUnclassFactor(index, addItem = add_items || sorted, sorted = sorted)

        if(add_items || sorted){
            # we re order appropriately
            # f prefix means factor

            obs_1st = cpp_get_first_item(res$x, length(res$items))

            f_all = list()
            for(q in 1:Q){
                f_all[[q]] = dots[[q]][obs_1st]
            }

            f_order = do.call("order", f_all)

            x_new = order(f_order)[res$x]

            if(multi.df){
                # Putting into a DF => we take care of names
                mc_dots = match.call(expand.dots = FALSE)[["..."]]
                n_dots = length(mc_dots)
                mc_dots_names = names(mc_dots)
                if(is.null(mc_dots_names)) mc_dots_names = character(n_dots)

                my_names = character(n_dots)
                for(q in 1:n_dots){
                    if(nchar(mc_dots_names[q]) > 0){
                        my_names[q] = mc_dots_names[q]
                    } else {
                        my_names[q] = deparse_long(mc_dots[[q]])
                    }
                }

                names(f_all) = my_names

                f_df = as.data.frame(f_all)
                items_new = f_df[f_order, , drop = FALSE]
                row.names(items_new) = 1:nrow(items_new)
            } else {
                # we "paste" them
                arg_list = f_all
                arg_list$sep = multi.join
                f_char = do.call("paste", arg_list)
                items_new = f_char[f_order]
            }

            if(add_items){
                res = list(x = x_new, items = items_new)
            } else {
                res = x_new
            }
        }
    }

    if(ANY_NA){
        if(is.list(res)){
            x = res$x
        } else {
            x = res
        }

        x_na = rep(NA, n)
        x_na[!is_na] = x

        if(is.list(res)){
            res$x = x_na
        } else {
            res = x_na
        }

    }

    if(add_items && isFALSE(items.list)){
        res_tmp = res$x
        attr(res_tmp, "items") = res$items
        res = res_tmp
    }

    res
}



#' Centers a set of variables around a set of factors
#'
#' User-level access to internal demeaning algorithm of \code{fixest}.
#'
#' @inheritSection feols Varying slopes
#'
#' @param X A matrix, vector, data.frame or a list OR a formula OR a \code{\link[fixest]{feols}} estimation. If equal to a formula, then the argument \code{data} is required, and it must be of the type: \code{x1 + x2 ~ f1 + fe2} with on the LHS the variables to be centered, and on the RHS the factors used for centering. Note that you can use variables with varying slopes with the syntax \code{fe[v1, v2]} (see details in \code{\link[fixest]{feols}}). If a \code{feols} estimation, all variables (LHS+RHS) are demeaned and then returned (only if it was estimated with fixed-effects). Otherwise, it must represent the data to be centered. Of course the number of observations of that data must be the same as the factors used for centering (argument \code{f}).
#' @param f A matrix, vector, data.frame or list. The factors used to center the variables in argument \code{X}. Matrices will be coerced using \code{as.data.frame}.
#' @param slope.vars A vector, matrix or list representing the variables with varying slopes. Matrices will be coerced using \code{as.data.frame}. Note that if this argument is used it MUST be in conjunction with the argument \code{slope.flag} that maps the factors to which the varying slopes are attached. See examples.
#' @param slope.flag An integer vector of the same length as the number of variables in \code{f} (the factors used for centering). It indicates for each factor the number of variables with varying slopes to which it is associated. Positive values mean that the raw factor should also be included in the centering, negative values that it should be excluded. Sorry it's complicated... but see the examples it may get clearer.
#' @param data A data.frame containing all variables in the argument \code{X}. Only used if \code{X} is a formula, in which case \code{data} is mandatory.
#' @param weights Vector, can be missing or NULL. If present, it must contain the same number of observations as in \code{X}.
#' @param nthreads Number of threads to be used. By default it is equal to \code{getFixest_nthreads()}.
#' @param notes Logical, whether to display a message when NA values are removed. By default it is equal to \code{getFixest_notes()}.
#' @param iter Number of iterations, default is 2000.
#' @param tol Stopping criterion of the algorithm. Default is \code{1e-6}. The algorithm stops when the maximum absolute increase in the coefficients values is lower than \code{tol}.
#' @param na.rm Logical, default is \code{TRUE}. If \code{TRUE} and the input data contains any NA value, then any observation with NA will be discarded leading to an output with less observations than the input. If \code{FALSE}, if NAs are present the output will also be filled with NAs for each NA observation in input.
#' @param as.matrix Logical, if \code{TRUE} a matrix is returned, if \code{FALSE} it will be a data.frame. The default depends on the input, if atomic then a matrix will be returned.
#' @param im_confident Logical, default is \code{FALSE}. FOR EXPERT USERS ONLY! This argument allows to skip some of the preprocessing of the arguments given in input. If \code{TRUE}, then \code{X} MUST be a numeric vector/matrix/list (not a formula!), \code{f} MUST be a list, \code{slope.vars} MUST be a list, \code{slope.vars} MUST be consistent with \code{slope.flag}, and \code{weights}, if given, MUST be numeric (not integer!). Further there MUST be not any NA value, and the number of observations of each element MUST be consistent. Non compliance to these rules may simply lead your R session to break.
#' @param fixef.reorder Logical, default is \code{TRUE}. Whether to reorder the fixed-effects by frequencies before feeding them into the algorithm. If \code{FALSE}, the original fixed-effects order provided by the user is maintained. In general, reordering leads to faster and more precise performance.
#' @param ... Not currently used.
#'
#' @return
#' It returns a data.frame of the same number of columns as the number of variables to be centered.
#'
#' If \code{na.rm = TRUE}, then the number of rows is equal to the number of rows in input minus the number of NA values (contained in \code{X}, \code{f}, \code{slope.vars} or \code{weights}). The default is to have an output of the same number of observations as the input (filled with NAs where appropriate).
#'
#' A matrix can be returned if \code{as.matrix = TRUE}.
#'
#' @examples
#'
#' # Illustration of the FWL theorem
#' data(trade)
#'
#' base = trade
#' base$ln_dist = log(base$dist_km)
#' base$ln_euros = log(base$Euros)
#'
#' # We center the two variables ln_dist and ln_euros
#' #  on the factors Origin and Destination
#' X_demean = demean(X = base[, c("ln_dist", "ln_euros")],
#'                   f = base[, c("Origin", "Destination")])
#' base[, c("ln_dist_dm", "ln_euros_dm")] = X_demean
#'
#' est = feols(ln_euros_dm ~ ln_dist_dm, base)
#' est_fe = feols(ln_euros ~ ln_dist | Origin + Destination, base)
#'
#' # The results are the same as if we used the two factors
#' # as fixed-effects
#' etable(est, est_fe, se = "st")
#'
#' #
#' # Variables with varying slopes
#' #
#'
#' # You can center on factors but also on variables with varying slopes
#'
#' # Let's have an illustration
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' #
#' # We center y and x1 on species and x2 * species
#'
#' # using a formula
#' base_dm = demean(y + x1 ~ species[x2], data = base)
#'
#' # using vectors
#' base_dm_bis = demean(X = base[, c("y", "x1")], f = base$species,
#'                      slope.vars = base$x2, slope.flag = 1)
#'
#' # Let's look at the equivalences
#' res_vs_1 = feols(y ~ x1 + species + x2:species, base)
#' res_vs_2 = feols(y ~ x1, base_dm)
#' res_vs_3 = feols(y ~ x1, base_dm_bis)
#'
#' # only the small sample adj. differ in the SEs
#' etable(res_vs_1, res_vs_2, res_vs_3, keep = "x1")
#'
#' #
#' # center on x2 * species and on another FE
#'
#' base$fe = rep(1:5, 10)
#'
#' # using a formula => double square brackets!
#' base_dm = demean(y + x1 ~ fe + species[[x2]], data = base)
#'
#' # using vectors => note slope.flag!
#' base_dm_bis = demean(X = base[, c("y", "x1")], f = base[, c("fe", "species")],
#'                      slope.vars = base$x2, slope.flag = c(0, -1))
#'
#' # Explanations slope.flag = c(0, -1):
#' # - the first 0: the first factor (fe) is associated to no variable
#' # - the "-1":
#' #    * |-1| = 1: the second factor (species) is associated to ONE variable
#' #    *   -1 < 0: the second factor should not be included as such
#'
#' # Let's look at the equivalences
#' res_vs_1 = feols(y ~ x1 + i(fe) + x2:species, base)
#' res_vs_2 = feols(y ~ x1, base_dm)
#' res_vs_3 = feols(y ~ x1, base_dm_bis)
#'
#' # only the small sample adj. differ in the SEs
#' etable(res_vs_1, res_vs_2, res_vs_3, keep = "x1")
#'
#'
#'
#'
demean = function(X, f, slope.vars, slope.flag, data, weights,
                  nthreads = getFixest_nthreads(), notes = getFixest_notes(),
                  iter = 2000, tol = 1e-6, fixef.reorder = TRUE, na.rm = TRUE,
                  as.matrix = is.atomic(X),
                  im_confident = FALSE, ...) {


    ANY_NA = FALSE
    # SK: to reassign class if X is data.frame, data.table or tibble. Optimally you would preserve all attributes,
    # but using attributes<- is slow on data.frames. What I did in collapse is export the SET_ATTRIB and DUPLICATE_ATTRIB
    # functions from the C-API to use them internally in R -> copy attributes without checks at 0 cost, even for large data.frames.

    # LB: next line is needed if input data is matrix and as.matrix is set to FALSE
    clx = NULL
    is_fixest = inherits(X, "fixest")
    if(lX <- is.list(X) && !is_fixest) {
        clx <- oldClass(X)
        # SK: oldClass is faster and returns NULL when the list is plain. class returns the implicit class "list".
        # This is the key to fast R code -> all data.frame methods are super slow and should duly be avoided in internal code.
        oldClass(X) <- NULL
    }
    # LB: The reassignment of attributes to data.frames is actually a very good idea, thanks Seb!

    # LB: To avoid delayed evaluation problems (due to new default is.atomic(X))
    as_matrix = as.matrix

    # internal argument
    fe_info = FALSE

    # Step 1: formatting the input
    if(!im_confident){

        check_arg(X, "numeric vmatrix | list | formula | class(fixest) mbt")
        check_arg(iter, "integer scalar GE{1}")
        check_arg(tol, "numeric scalar GT{0}")
        check_arg(notes, "logical scalar")

        validate_dots(valid_args = "fe_info", stop = TRUE)
        dots = list(...)
        fe_info = isTRUE(dots$fe_info)

        data_mbt = TRUE
        if(inherits(X, "fixest")){
            data_mbt = FALSE
            if(!identical(X$method, "feols")){
                stop("This function only works for 'feols' estimations (not for ", X$method, ").")
            }

            if(!fe_info && !is.null(X$y_demeaned)){
                return(cbind(X$y_demeaned, X$X_demeaned))
            }

            if(!"fixef_vars" %in% names(X)){
                stop("To apply demeaning, the estimation must contain fixed-effects, this is currently not the case.")
            }

            # Otherwise => we create data and X: formula
            data = fetch_data(X)
            fml_linear = X$fml_all$linear
            fml_fixef = X$fml_all$fixef

            fml_vars = .xpd(~ ..lhs + ..rhs,
                            ..lhs = fml_linear[[2]],
                            ..rhs = fml_linear[[3]])

            if(!is.null(X$fml_all$iv)){
                fml_iv = X$fml_all$iv
                fml_vars = .xpd(~ ..vars + ..iv_lhs + ..iv_rhs,
                                ..vars = fml_vars,
                                ..iv_lhs = fml_iv[[2]],
                                ..iv_rhs = fml_iv[[3]])
            }

            fml = .xpd(lhs = fml_vars, rhs = fml_fixef)

            mc = match.call()
            if(!"tol" %in% names(mc)) tol = X$fixef.tol
            if(!"iter" %in% names(mc)) iter = X$fixef.iter

            weights = X[["weights"]]

            X = fml

            as_matrix = TRUE
        }

        #
        # X
        #

        fe_done = FALSE
        if(is.call(X)) {

            if(data_mbt) check_arg(data, "data.frame mbt")
            check_arg(X, "ts formula var(data)", .data = data)

            # Extracting the information
            terms_fixef = fixef_terms(.xpd(rhs = X[[3L]]))
            # We add the intercept only for only slope models, otherwise this would be void since the result would be always 0
            X = fixest_model_matrix(.xpd(lhs = quote(y), rhs = X[[2L]]), data, fake_intercept = any(terms_fixef$slope_flag >= 0))
            var_names = dimnames(X)[[2]]

            lX = FALSE # Needed for the rest of the code to work

            # FE
            fe_done = TRUE
            f = unclass(error_sender(prepare_df(terms_fixef$fe_vars, data),
                                     "Error when evaluating the fixed-effects variables: "))

            isSlope = any(terms_fixef$slope_flag != 0)

            if(isSlope) {

                slope.vars = unclass(error_sender(prepare_df(terms_fixef$slope_vars, data),
                                                  "Error when evaluating the variable with varying slopes: "))

                if(anyDuplicated(terms_fixef$slope_vars)){
                    # we need it!!!!
                    slope.vars = slope.vars[terms_fixef$slope_vars]
                }

                slope.flag = terms_fixef$slope_flag

                is_numeric = vapply(`attributes<-`(slope.vars, NULL), is.numeric, TRUE)
                if(!all(is_numeric)){
                    stop("In the fixed-effects part of the formula, variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope.vars)[!is_numeric], "s.is.quote"), " not.")
                }
            } else {
                slope.vars = list(0)
                slope.flag = rep(0L, length(f))
            }

        } else if(lX) {
            var_names = names(X)
            if(!any(clx == "data.frame")) {
                n_all = lengths(X)
                if(!all(n_all == n_all[1L])) {
                    n_unik = unique(n_all)
                    stop("In argument 'X' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                }
            }

            # SK: If your using Rcpp and the input is NumericVector or NumericMatrix, you should get a C++ error for wrongly typed data, so I don't see the need for this
            # LB: I disagree, such error messages are poor and hard to debug for end users.

            is_numeric = vapply(`attributes<-`(X, NULL), is.numeric, TRUE)
            if(!all(is_numeric)){
                # Faster than any(!is_numeric) => LB: yet a few nano seconds saved... ;-)
                n_non_num = sum(!is_numeric)
                stop("All values in 'X' must be numeric. This is not the case for ", n_non_num, " variable", plural(n_non_num), ".")
            }

        } else {
            # Here: numeric matrix or vector
            var_names = dimnames(X)[[2L]]
        }

        # nobs: useful for checks
        nobs = if(is.list(X)) length(.subset2(X, 1L)) else NROW(X)

        #
        # f + slope.vars
        #

        if(!fe_done){
            check_arg(f, "vmatrix | list mbt")
            check_arg(slope.vars, "numeric vmatrix | list")
            check_arg(slope.flag, "integer vector no na")

            if(is.list(f)) {
                if(!is.data.frame(f)) {
                    n_all = lengths(f)
                    if(!all(n_all == n_all[1L])) {
                        n_unik = unique(n_all)
                        stop("In argument 'f' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                    }
                    # f = as.data.frame(f)
                    # SK: as.data.frame is super slow, especially on large data. You already checked the lengths, so it's ok
                    # LB: true
                } else {
                    oldClass(f) <- NULL
                }

            } else {
                # SK: as.data.frame on matrix is also very slow, you could do in C++ as in at collapse::mctl.. but I suppose most will pass lists of factors anyway..
                # LB: In the long run I'll deal with that but that's really low priority
                f = if(!is.array(f)) list(f) else unclass(as.data.frame(f))
            }

            is_pblm = vapply(`attributes<-`(f, NULL), function(x) !(is.numeric(x) || is.character(x)), TRUE)
            if(any(is_pblm)) {
                for(i in which(is_pblm)) f[[i]] = as.character(f[[i]])
            }

            if(length(f[[1L]]) != nobs){
                stop("The number of observations in 'X' and in 'f' don't match (", if(lX) length(X[[1L]]) else NROW(X), " vs ", length(f[[1L]]), ").")
            }

            # Now the slopes
            isSlope = FALSE
            if(!missnull(slope.vars) || !missnull(slope.flag)) {

                isSlope = TRUE

                if(missnull(slope.vars)) stop("If argument 'slope.flag' is provided, you must also provide the argument 'slope.vars'.")
                if(missnull(slope.flag)) stop("If argument 'slope.vars' is provided, you must also provide the argument 'slope.flag'.")

                if(is.list(slope.vars)){
                    if(!is.data.frame(slope.vars)){
                        n_all = lengths(slope.vars)
                        if(!all(n_all == n_all[1L])) {
                            n_unik = unique(n_all)
                            stop("In argument 'slope.vars' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                        }

                    } else {
                        oldClass(slope.vars) <- NULL
                    }

                    is_numeric = vapply(`attributes<-`(slope.vars, NULL), is.numeric, TRUE)
                    # SK: Much faster than !sapply(slope.vars, is.numeric)
                    # LB: 10us diff, only impedes readability IMO
                    if(!all(is_numeric)) stop("In the argument 'slope.vars', the variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope.vars)[!is_numeric], "s.is.quote"), " not.")

                } else {
                    # SK: Same as above.. Necessary to convert ??
                    # LB: The new C++ code accepts lists of vector/matrices or vector/matrices directly // So if we're here (meaning vector or matrix), that's fine
                }

                if(length(slope.flag) != length(f)) stop("The argument 'slope.flag' must be a vector representing, for each fixed-effect, the number of variables with varying slopes associated to it. Problem: the lengths of slope.flag and the fixed-effect differ: ", length(slope.flag), " vs ", length(f), ".")


                # SK: if(sum(abs(slope.flag)) != length(slope.vars))  : This means that slope flag can only be 0, 1 or -1. In the documentation you only talk about positive and negative values. The change I made reflects this.
                # LB: Cmon Sebastian why change sum(abs(slope.flag)) != length(slope.vars) into sum(slope.flag != 0L) != length(slope.vars)?  The time gains lie in a handful of nano second and you've just introduced a bug!

                n_sv = if(is.list(slope.vars)) length(slope.vars) else NCOL(slope.vars)
                if(sum(abs(slope.flag)) != n_sv) stop("The argument 'slope.flag' must be a vector representing, for each fixed-effect, the number of variables with varying slopes associated to it. Problem: currently the number of variables in 'slope.flag' differ from the number of variables in 'slope.vars': ", sum(abs(slope.flag)), " vs ", n_sv, ".")

            } else {
                slope.vars = list(0)
                slope.flag = rep(0L, length(f))
            }
        }


        ## weights
        check_arg_plus(weights, "NULL numeric conv vector len(value) GE{0}", .value = nobs)
        is_weight = !missing(weights) && !is.null(weights)

        ## nthreads (Also seems unnecessarily costly: 250 microseconds)
        # LB: I know... but it has to be checked
        if(!missing(nthreads)) nthreads = check_set_nthreads(nthreads)

        #
        # FORMATTING
        #

        # NAs
        is_NA = FALSE
        info_X = cpp_which_na_inf(X, nthreads)

        if(info_X$any_na_inf){
            is_NA = info_X$is_na_inf
            n_na_x = 1L
        } else n_na_x = 0L


        if(anyNA(f)){
            is_na_fe = !complete.cases(f)
            n_na_fe = 1L
            is_NA = is_NA | is_na_fe
        } else n_na_fe = 0L


        is_na_slope = 0L
        if(isSlope){
            info_slopes = cpp_which_na_inf(slope.vars, nthreads)

            if(info_slopes$any_na_inf){
                is_na_slope = info_slopes$is_na_inf
                is_NA = is_NA | is_na_slope
            }
        }


        is_na_weights = 0L
        if(is_weight){
            info_weights = cpp_which_na_inf(weights, nthreads)

            if(info_weights$any_na_inf){
                is_na_weights = info_weights$is_na_inf
                is_NA = is_NA | is_na_weights
            }
        }

        n_na = sum(is_NA)
        if(n_na && notes) {
            # Here is all the error message stuff now: Only evaluated if needed
            if(n_na_x) n_na_x = sum(info_X$is_na_inf)
            if(n_na_fe) n_na_fe = sum(is_na_fe)
            slopes_msg = if(isSlope) paste0(", slope vars: ", sum(is_na_slope)) else ""
            weight_msg = if(is_weight) paste0(", weights: ", sum(is_na_weights)) else ""

            if(n_na == length(is_NA)) stop("All observations contain NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, slopes_msg, weight_msg, ").")
            message("NOTE: ", fsignif(n_na), " observation", plural(n_na), " removed because of NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, slopes_msg, weight_msg, ").")
        }

        if(n_na) {
            ANY_NA = TRUE
            # SK: subsetting integer is faster than logical !!
            # LB: Indeed, just checked, quite a diff. Good to know!
            cc <- which(!is_NA)
            X = select_obs(X, cc)

            f = lapply(f, `[`, cc)
            if(isSlope) slope.vars = lapply(slope.vars, `[`, cc)
            if(is_weight) weights = weights[cc]
        }

        if(!is_weight) weights = 1

    } else {
        # Need this here, otherwise error:
        var_names = if(lX) names(X) else dimnames(X)[[2L]]
        if(missing(weights) || is.null(weights)) weights = 1
        if(missnull(slope.vars) || missnull(slope.flag)){
            slope.vars = list(0)
            slope.flag = rep(0L, length(unclass(f)))
            # SK: unclass gives extra speed if f is data.frame, and no cost if f is list.
        }
    }

    #
    # Unclassing fes
    #

    quf_info_all = cpppar_quf_table_sum(x = f, y = 0, do_sum_y = FALSE, rm_0 = FALSE,
                                        rm_1 = FALSE, rm_single = FALSE, do_refactor = FALSE,
                                        r_x_sizes = 0, obs2keep = 0, only_slope = slope.flag < 0L,
                                        nthreads = nthreads)

    # table/sum_y/sizes
    fixef_table = quf_info_all$table
    fixef_sizes = lengths(fixef_table)

    if(fixef.reorder){
        new_order = order(fixef_sizes, decreasing = TRUE)
        if(is.unsorted(new_order)){
            fixef_table = fixef_table[new_order]
            fixef_sizes = fixef_sizes[new_order]
            quf_info_all$quf = quf_info_all$quf[new_order]

            # We reorder slope.vars only if needed (because it is a pain)
            if(sum(slope.flag != 0) >= 2){

                # here slope.vars have 2 or more variables:
                # either matrix or data.frame/list
                if(!is.list(slope.vars)){
                    slope.vars = as.data.frame(slope.vars)
                }
                # we ensure it's a plain list
                slope.vars = unclass(slope.vars)

                n_fixef = length(slope.flag)
                slope_vars_list = vector("list", n_fixef)
                id_current = 0
                for(i in 1:n_fixef){
                    n_slopes = abs(slope.flag[i])

                    if(n_slopes == 0){
                        slope_vars_list[[i]] = 0
                    } else {
                        slope_vars_list[[i]] = slope.vars[1:n_slopes + id_current]
                        id_current = id_current + n_slopes
                    }
                }

                # Reordering => currently it's quite clumsy (but I HATE VSs!)
                slope.vars = vector("list", sum(abs(slope.flag)))
                slope_vars_list_reordered = slope_vars_list[new_order]
                id_current = 0
                for(i in 1:n_fixef){
                    vs = slope_vars_list_reordered[[i]]
                    if(is.list(vs)){
                        n_vs = length(vs)
                        for(j in 1:n_vs){
                            slope.vars[[id_current + j]] = vs[[j]]
                        }
                        id_current = id_current + n_vs
                    }
                }

                # slope.flag must be reordered afterwards
                slope.flag = slope.flag[new_order]
            }
        }
    }

    fixef_table_vector = unlist(fixef_table)
    if(!is.integer(slope.flag)) slope.flag = as.integer(slope.flag)

    #
    # The algorithm
    #

    # y => list, X => matrix
    if(as_matrix){
        y = 0
    } else {
        # Quirk => y returns a list only if NCOL(y) > 1 or is.list(y)
        y = if(lX) X else list(X)
        # SK: This does the same, a bit more frugal
        # LB: I tend to prefer late evaluations instead of lX which requires bookkeeping
        X = 0
    }

    vars_demean = cpp_demean(y, X, weights, iterMax = iter,
                              diffMax = tol, r_nb_id_Q = fixef_sizes,
                              fe_id_list = quf_info_all$quf, table_id_I = fixef_table_vector,
                              slope_flag_Q = slope.flag, slope_vars_list = slope.vars,
                              r_init = 0, nthreads = nthreads)

    # Internal call
    if(fe_info){

        res = list(y = vars_demean$y_demean, X = vars_demean$X_demean,
                   weights = weights, fixef_id_list = quf_info_all$quf,
                   slope_vars = slope.vars, slope_flag = slope.flag, varnames = colnames(X))

        return(res)

    }


    if(as_matrix) {
        if(ANY_NA && na.rm == FALSE) {
            K = NCOL(vars_demean$X_demean)
            res = matrix(NA_real_, length(is_NA), K)
            res[cc, ] = vars_demean$X_demean

        } else {
            res = vars_demean$X_demean
        }

        dimnames(res)[[2L]] = var_names

    } else {
        if(ANY_NA && na.rm == FALSE) {
            n = length(is_NA)
            # SK: One-line efficient solution to the task:
            res = lapply(vars_demean$y_demean, function(x) `[<-`(rep(NA_real_, n), cc, value = x))
            # LB: nice compactness

        } else {
            res = vars_demean$y_demean
        }

        names(res) = var_names
        # SK: This makes your function class-agnostic to lists, data.frame's, data.table's and tibbles (with the exception for any custom data.frame row.names which are not preserved, but the as.data.frame solution also did not do that)
        if(!lX || length(clx)) attr(res, "row.names") <- .set_row_names(length(res[[1L]]))
        # SK: Here !lX || length(clx) means add row.names if either X was a matrix or is classed that is a data.frame
        oldClass(res) <- if(lX) clx else "data.frame"
        # SK: If X is a plain list, and since oldClass returns NULL, this will assign a null class again.
    }
    res
}

#' Extracts the observations used for the estimation
#'
#' This function extracts the observations used in \code{fixest} estimation.
#'
#' @param x A \code{fixest} object.
#'
#' @return
#' It returns a simple vector of integers.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#' base$y[1:5] = NA
#'
#' # Split sample estimations
#' est_split = feols(y ~ x1, base, split = ~species)
#' (obs_setosa = obs(est_split$setosa))
#' (obs_versi = obs(est_split$versicolor))
#'
#' est_versi = feols(y ~ x1, base, subset = obs_versi)
#'
#' etable(est_split, est_versi)
#'
#'
#'
#'
obs = function(x){
    check_arg(x, "class(fixest)")

    if(isTRUE(x$lean)){
        stop("obs() does not work with models estimated with 'lean = TRUE'.")
    }

    id = 1:x$nobs_origin

    for(i in seq_along(x$obs_selection)){
        id = id[x$obs_selection[[i]]]
    }

    return(id)
}




#' Check the fixed-effects convergence of a \code{feols} estimation
#'
#' Checks the convergence of a \code{feols} estimation by computing the first-order conditions of all fixed-effects (all should be close to 0)
#'
#' @param x A \code{\link[fixest]{feols}} estimation that should contain fixed-effects.
#' @param object An object returned by \code{check_conv_feols}.
#' @param type Either "short" (default) or "detail". If "short", only the maximum absolute FOC are displayed, otherwise the 2 smallest and the 2 largest FOC are reported for each fixed-effect and each variable.
#' @param ... Not currently used.
#'
#' Note that this function first re-demeans the variables, thus possibly incurring some extra computation time.
#'
#' @return
#' It returns a list of \code{N} elements, \code{N} being the number of variables in the estimation (dependent variable + explanatory variables +, if IV, endogenous variables and instruments). For each variable, all the first-order conditions for each fixed-effect are returned.
#'
#' @examples
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' base$FE = rep(1:30, 5)
#'
#' # one estimation with fixed-effects + varying slopes
#' est = feols(y ~ x1 | species[x2] + FE[x3], base)
#'
#' # Checking the convergence
#' conv = check_conv_feols(est)
#'
#' # We can check that al values are close to 0
#' summary(conv)
#'
#' summary(conv, "detail")
#'
#'
#'
check_conv_feols = function(x){

    check_arg(x, "class(fixest) mbt")
    if(!identical(x$method, "feols")){
        stop("This function only works for 'feols' estimations (not for ", x$method, ").")
    }

    if(!"fixef_vars" %in% names(x)){
        message("This function only works with fixed-effects (which are currently not present).")
        return(NULL)
    }

    # fixef names for information
    new_order = x$fe.reorder
    fixef_vars = x$fixef_vars[new_order]
    if(!is.null(x$fixef_terms)){

        slope_variables = x$slope_variables_reordered
        slope_flag = x$slope_flag_reordered

        # We reconstruct the terms
        fixef_names = c()
        start = c(0, cumsum(abs(slope_flag)))
        for(i in seq_along(slope_flag)){
            sf = slope_flag[i]
            if(sf >= 0){
                fixef_names = c(fixef_names, fixef_vars[i])
            }

            if(abs(sf) > 0){
                fixef_names = c(fixef_names, paste0(fixef_vars[i], "[[", names(slope_variables)[start[i] + 1:abs(sf)], "]]"))
            }
        }
    } else {
        fixef_names = fixef_vars
    }


    info = demean(x, fe_info = TRUE)

    res = check_conv(y = info$y, X = info$X, fixef_id_list = info$fixef_id_list,
                     slope_flag = info$slope_flag, slope_vars = info$slope_vars,
                     weights = info$weights, full = TRUE, fixef_names = fixef_names)

    names(res) = info$varnames

    class(res) = "fixest_check_conv"

    res
}

#' @rdname check_conv_feols
summary.fixest_check_conv = function(object, type = "short", ...){

    check_arg_plus(type, "match(short, detail)")
    if(is_user_level_call()){
        validate_dots(suggest_args = "type")
    }

    if(type == "short"){
        info_max_abs = lapply(object, function(x) sapply(x, function(y) max(abs(y))))

        info_max_abs = do.call("rbind", info_max_abs)

        cat("Maximum absolute value of the first-order conditions:\n\n")
        print(info_max_abs)

    } else {

        extract_and_format = function(x){
            # x: list

            x_sorted = lapply(x, sort)

            is_short = sapply(x, function(y) length(y) <= 4)

            x_small = list()
            for(i in seq_along(x)){
                xi = x_sorted[[i]]
                if(is_short[i]){
                    x_small[[i]] = c(xi, rep(NA, 4 - length(xi)))
                } else {
                    x_small[[i]] = c(head(xi, 2), tail(xi, 2))
                }
            }

            x_mat = format(do.call("rbind", x_small), digits = 3)

            x_fmt = c()
            for(i in seq_along(x)){
                xi = x_mat[i, ]
                if(is_short[i]){
                    x_fmt[[i]] = gsub(", +NA", "", paste0(xi, collapse = ", "))
                } else {
                    x_fmt[[i]] = paste0(xi[1], ", ", xi[2], ", ..., ", xi[3], ", ", xi[4])
                }
            }

            x_names = sfill(names(x))

            res = paste0(x_names, ": ", x_fmt)

            res
        }

        info = lapply(object, extract_and_format)

        cat("Smallest and largest values of the first-order conditions:\n\n")
        var_names = sfill(names(object), right = TRUE)
        intro = sfill("|", n = nchar(var_names[1]) + 1)
        n_var = length(object)
        for(i in 1:n_var){
            cat(var_names[i], "\n")
            value = paste0(intro, info[[i]])
            cat(value, sep = "\n")
            cat(gsub(" ", "_", intro), "\n")

            if(i < n_var) cat("\n")
        }
    }
}




#' Prints the number of unique elements in a data set
#'
#' This utility tool displays the number of unique elements in one or multiple data.frames as well as their number of NA values.
#'
#' @param x A formula, with data set names on the LHS and variables on the RHS, like \code{data1 + data2 ~ var1 + var2}. The following special variables are admitted: \code{"."} to get default values, \code{".N"} for the number of observations, \code{".U"} for the number of unique rows, \code{".NA"} for the number of rows with at least one NA. Variables can be combined with \code{"^"}, e.g. \code{df~id^period}; use \code{id\%^\%period} to also include the terms on both sides. Note that using \code{:} and \code{*} is equivalent to \code{^} and \code{\%^\%}. Sub select with \code{id[cond]}, when doing so \code{id} is automatically included. Conditions can be chained, as in \code{id[cond1, cond2]}. Use \code{NA(x, y)} in conditions instead of \code{is.na(x) | is.na(y)}. Use the \code{!!} operator to have both a condition and its opposite. To compare the keys in two data sets, use \code{data1:data2}. If not a formula, \code{x} can be: a vector (displays the # of unique values); a \code{data.frame} (default values are displayed), or a "sum" of data sets like in \code{x = data1 + data2}, in that case it is equivalent to \code{data1 + data2 ~ .}.
#' @param ... Not currently used.
#'
#' @section Special values and functions:
#'
#' In the formula, you can use the following special values: \code{"."}, \code{".N"}, \code{".U"}, and \code{".NA"}.
#'
#' \itemize{
#'
#' \item{\code{"."}}{Accesses the default values. If there is only one data set and the data set is \emph{not} a \code{data.table}, then the default is to display the number of observations and the number of unique rows. If the data is a \code{data.table}, the number of unique items in the key(s) is displayed instead of the number of unique rows (if the table has keys of course). If there are two or more data sets, then the default is to display the unique items for: a) the variables common across all data sets, if there's less than 4, and b) if no variable is shown in a), the number of variables common across at least two data sets, provided there are less than 5. If the data sets are data tables, the keys are also displayed on top of the common variables. In any case, the number of observations is always displayed.}
#'
#' \item{\code{".N"}}{Displays the number of observations.}
#'
#' \item{\code{".U"}}{Displays the number of unique rows.}
#'
#' \item{\code{".NA"}}{Displays the number of rows with at least one NA.}
#'
#' }
#'
#' @section The \code{NA} function:
#'
#' The special function \code{NA} is an equivalent to \code{is.na} but can handle several variables. For instance, \code{NA(x, y)} is equivalent to \code{is.na(x) | is.na(y)}. You can add as many variables as you want as arguments. If no argument is provided, as in \code{NA()}, it is identical to having all the variables of the data set as argument.
#'
#' @section Combining variables:
#'
#' Use the "hat", \code{"^"}, operator to combine several variables. For example \code{id^period} will display the number of unique values of id x period combinations.
#'
#' Use the "super hat", \code{"\%^\%"}, operator to also include the terms on both sides. For example, instead of writing \code{id + period + id^period}, you can simply write \code{id\%^\%period}.
#'
#' Alternatively, you can use \code{:} for \code{^} and \code{*} for \code{\%^\%}.
#'
#' @section Sub-selections:
#'
#' To show the number of unique values for sub samples, simply use \code{[]}. For example, \code{id[x > 10]} will display the number of unique \code{id} for which \code{x > 10}.
#'
#' Simple square brackets lead to the inclusion of both the variable and its subset. For example \code{id[x > 10]} is equivalent to \code{id + id[x > 10]}. To include only the sub selection, use double square brackets, as in \code{id[[x > 10]]}.
#'
#' You can add multiple sub selections at once, only separate them with a comma. For example \code{id[x > 10, NA(y)]} is equivalent to \code{id[x > 10] + id[NA(y)]}.
#'
#' Use the double negative operator, i.e. \code{!!}, to include both a condition and its opposite at once. For example \code{id[!!x > 10]} is equivalent to \code{id[x > 10, !x > 10]}. Double negative operators can be chained, like in \code{id[!!cond1 & !!cond2]}, then the cardinal product of all double negatived conditions is returned.
#'
#' @return
#' It returns a vector containing the number of unique values per element. If several data sets were provided, a list is returned, as long as the number of data sets, each element being a vector of unique values.
#'
#' @examples
#'
#' data = base_did
#' data$x1.L1 = round(lag(x1~id+period, 1, data))
#'
#' # By default, just the formatted number of observations
#' n_unik(data)
#'
#' # Or the nber of unique elements of a vector
#' n_unik(data$id)
#'
#' # number of unique id values and id x period pairs
#' n_unik(data ~.N + id + id^period)
#'
#' # use the %^% operator to include the terms on the two sides at once
#' # => same as id*period
#' n_unik(data ~.N + id %^% period)
#'
#' # using sub selection with []
#' n_unik(data ~.N + period[!NA(x1.L1)])
#'
#' # to show only the sub selection: [[]]
#' n_unik(data ~.N + period[[!NA(x1.L1)]])
#'
#' # you can have multiple values in [],
#' # just separate them with a comma
#' n_unik(data ~.N + period[!NA(x1.L1), x1 > 7])
#'
#' # to have both a condition and its opposite,
#' # use the !! operator
#' n_unik(data ~.N[!!NA(x1.L1)])
#'
#' # the !! operator works within condition chains
#' n_unik(data ~.N[!!NA(x1.L1) & !!x1 > 7])
#'
#' # Conditions can be distributed
#' n_unik(data ~ (id + period)[x1 > 7])
#'
#' #
#' # Several data sets
#' #
#'
#' # Typical use case: merging
#' # Let's create two data sets and merge them
#'
#' data(base_did)
#' base_main = base_did
#' base_extra = sample_df(base_main[, c("id", "period")], 100)
#' base_extra$id[1:10] = 111:120
#' base_extra$period[11:20] = 11:20
#' base_extra$z = rnorm(100)
#'
#' # You can use db1:db2 to compare the common keys in two data sets
#'  n_unik(base_main:base_extra)
#'
#' tmp = merge(base_main, base_extra, all.x = TRUE, by = c("id", "period"))
#'
#' # You can show unique values for any variable, as before
#' n_unik(tmp + base_main + base_extra ~ id[!!NA(z)] + id^period)
#'
#'
#'
n_unik = function(x){
    # returns a vector with the nber of unique values
    # attr("na.info") => nber of NA values, vector

    if(missing(x)){
        stop("Argument 'x' must be provided. Problem: it is missing.")
    }

    # Non standard-evaluation
    x_dp = deparse_long(substitute(x))
    if(!grepl("~", x_dp, fixed = TRUE)){
        if(grepl("[+:*]", x_dp)){
            # Non standard-evaluation
            x = as.formula(paste0(x_dp, "~ ."))
        } else {
            check_arg(x, "data.frame | vector l0 | ts formula")
        }
    } else {
        check_arg(x, "ts formula")
    }

    nthreads = getFixest_nthreads()

    # If vector
    comp_pairs = list()
    if(is.vector(x)){

        x_name = gsub("^[[:alpha:]\\.][[:alnum:]\\._]*\\$", "", x_dp)

        na.info = 0
        if(anyNA(x)){
            who_NA = is.na(x)
            na.info = sum(who_NA)
            x = x[!who_NA]
        }

        res = len_unique(x, nthreads)
        names(res) = x_name

        attr(res, "na.info") = na.info

        class(res) = "vec_n_unik"
        return(res)

    } else if(is.data.frame(x)){

        n_x = 1
        x_all = list(x)
        fml = ~.

    } else if(inherits(x, "formula")){

        x_terms = terms(x[1:2])

        fact_mat = attr(x_terms, "factors")

        x_all_names = rownames(fact_mat)

        # We get the comparison pairs
        vars = colnames(fact_mat)
        for(pair in grep(":", vars, fixed = TRUE, value = TRUE)){
            dict = setNames(1:length(x_all_names), x_all_names)
            pair_split = strsplit(pair, ":", fixed = TRUE)[[1]]
            comb = combn(pair_split, 2)
            for(i in 1:ncol(comb)){
                comp_pairs[[length(comp_pairs) + 1]] = unname(dict[comb[, i]])
            }
        }

        if(length(comp_pairs) > 0){
            comp_pairs = unique(comp_pairs)
        }

        # Construction of the list  + sanity check
        n_x = length(x_all_names)
        x_all = vector("list", n_x)
        for(i in 1:n_x){
            x_all[[i]] = error_sender(eval(str2lang(x_all_names[i]), parent.frame()),
                                      "The left-hand-side of the formula must contain valid data.frames. Problem in the evaluation of '", x_all_names[i], "':")
            check_value(x_all[[i]], "data.frame",
                        .message = paste0("The value '", x_all_names[i], "' (in the LHS of the formula) must be a data.frame."))
        }

        fml = x[c(1, 3)]

    }

    # If DF + formula
    # ex: fml = ~ id^year + author_id[sw0(is.na(author_name), is.na(affil_name), year == min_year)]
    # fml = ~ id^sw0(fe)

    # check variable names
    naked_vars = origin_vars = all.vars(fml)
    valid_vars = c(".N", ".U", ".NA", ".", unique(unlist(lapply(x_all, names))))

    check_value_plus(naked_vars, "multi match", .choices = valid_vars,
                     .message = paste0("The formula must only use variables in the data set", plural_len(x_all), "."))

    if("." %in% naked_vars){
        # The default values

        IS_DT = requireNamespace("data.table", quietly = TRUE)

        dot_default = ".N"
        x_keys = c()
        x_common_vars = c()

        # keys
        if(IS_DT){
            for(I in 1:n_x){
                if(inherits(x_all[[I]], "data.table")){
                    x_keys = c(x_keys, data.table::key(x_all[[I]]))
                }
            }
        }

        # common variables
        if(n_x > 1){

            # Common to all
            x_names_current = names(x_all[[1]])
            for(i in 2:n_x){
                x_names_current = intersect(x_names_current, names(x_all[[i]]))
                if(length(x_names_current) == 0) break
            }

            if(length(x_names_current) %in% 1:4){
                # No more than 4 common variables by default
                x_common_vars = x_names_current

                if(length(comp_pairs) > 0 && length(x_names_current) <= 3){
                    x_common_vars = paste0(x_names_current, collapse = "*")
                }

            } else if(n_x > 2){
                # Common to at least 2 data sets
                for(i in 1:(n_x - 1)){
                    x_names_current = names(x_all[[i]])
                    for(j in (i + 1):n_x){
                        qui_common = x_names_current %in% names(x_all[[j]])
                        if(any(qui_common)){
                            x_common_vars = c(x_common_vars, x_names_current[qui_common])
                        }
                    }
                }

                x_common_vars = unique(x_common_vars)

                if(length(x_common_vars) > 5){
                    #  we keep max 5
                    x_common_vars = c()
                }
            }
        }

        if(length(x_keys) > 0 || length(x_common_vars) > 0){
            dot_default = unique(c(dot_default, x_keys, x_common_vars))
        }

        if(length(dot_default) == 1){
            dot_default = c(".N", ".U")
        }

        rhs_txt = as.character(fml)[[2]]
        rhs_txt = gsub("(^| )\\.( |$)", dsb(".['+'c? dot_default]"), rhs_txt)
        fml = .xpd(rhs = rhs_txt)
    }

    if(any(naked_vars != origin_vars)){
        # we complete the partial matching
        fml_txt = deparse_long(fml)

        var_diff = which(naked_vars != origin_vars)

        for(i in var_diff){
            re = paste0("(?!<[[:alnum:]._])", origin_vars[i], "(?!=[[:alnum:]._])")
            fml_txt = gsub(re, naked_vars[i], fml_txt, perl = TRUE)
        }

        fml = as.formula(fml_txt)
    }

    tm = terms_hat(fml, n_unik = TRUE)
    all_vars = attr(tm, "term.labels")
    var_final = c()

    # stepwise function
    for(var in all_vars){
        # var = all_vars[1]

        if(grepl("combine_clusters(_fast)?", var)){
            var = gsub("combine_clusters(_fast)?", "to_integer", var)
        }

        IS_HAT_SW = grepl("sw0?\\)\\(", var)
        if(IS_HAT_SW){
            var = paste0(gsub("(sw0?)\\)\\(", "\\1(", var), ")")
        }

        if(grepl("^[[:alpha:].][[:alnum:]._]*\\[", var, perl = TRUE) && !grepl("sw0?\\(", var)){
            # We use sw to make id[cond1, cond2] work
            # That's what's called path dependency!

            var_new = gsub("^([[:alpha:].][[:alnum:]._]*)\\[", "\\1[sw0(", var)
            if(grepl("sw0([", var_new, fixed = TRUE)){
                var_new = sub("sw0([", "sw(", str_trim(var_new, -1), fixed = TRUE)
            }
            var_new = sub("\\]$", ")]", var_new)

            var = var_new
        }

        if(is_fun_in_char(var, "sw0?")){

            info = extract_fun(var, "sw0?", err_msg = "The stepwise function can be used only once per variable.")

            sw_value = eval(str2lang(info$fun))

            qui_double = grepl("!!", trimws(sw_value))
            while(any(qui_double)){
                i = which(qui_double)[1]

                value = sw_value[i]
                value_split = strsplit(value, "!!")[[1]]

                n_split = length(value_split)
                n_s = 2 ** (n_split - 1)
                new_values = rep(value_split[1], n_s)

                prefixes = do.call(expand.grid, rep(list(c("", "!")), n_split - 1))

                for(j in 1:ncol(prefixes)){
                    new_values = paste0(new_values, prefixes[, j], value_split[j + 1])
                }

                sw_value = insert(sw_value[-i], new_values, i)
                qui_double = grepl("!!", sw_value)
            }

            var_char_new = paste0(info$before, sw_value, info$after)

            if(IS_HAT_SW){
                var_char_new[1] = gsub(", (,|\\))", "\\1", var_char_new[1])

                qui_nested = grepl("^(to_integer\\(.+){2,}", var_char_new)
                if(any(qui_nested)){
                    # later => add check problem if another function is used
                    var_char_new[qui_nested] = paste0("to_integer(", gsub("to_integer\\(|\\)", "", var_char_new[qui_nested]), ")")
                }
            }

            if(grepl("[]", var_char_new[1], fixed = TRUE)){
                var_char_new[1] = gsub("[]", "", var_char_new[1], fixed = TRUE)
            }

            var_final = c(var_final, var_char_new)
        } else {
            var_final = c(var_final, var)
        }
    }

    var_final = unique(var_final)

    var_final_names = var_final
    # we take care of id^year cases
    if(any(who_combine <- is_fun_in_char(var_final, "to_integer"))){

        for(i in which(who_combine)){
            info = extract_fun(var_final[i], "to_integer")

            # we use sw() to get the vars
            sw_char = gsub("^[^\\(]+\\(", "sw(", info$fun)
            sw_value = eval(str2lang(sw_char))

            hat_value = paste(sw_value, collapse = "^")

            var_final_names[i] = paste0(info$before, hat_value, info$after)
        }
    }

    # NA counting + unique counting

    res_all = list()
    for(I in 1:n_x){

        x = x_all[[I]]
        x_list = unclass(x)
        x_list[["NA_fun"]] = function(...) NA_fun(..., df = x)

        vars_legit = c(".N", ".U", ".NA", names(x))

        res = c()
        na.info = c()

        for(i in seq_along(var_final)){

            vf = var_final[i]
            vf_name = var_final_names[i]
            na_i = 0

            vf = gsub("((?<=[^[:alnum:]_\\.])|^)NA\\(", "NA_fun(", vf, perl = TRUE)

            vf_call = str2lang(vf)

            if(!all(all.vars(vf_call) %in% vars_legit)){
                res_i = NA_real_
                vf_name = ""

            } else if(grepl("^\\.(N|U)($|\\[)", vf)){


                if(vf %in% c(".N", ".N[]")){
                    res_i = nrow(x)
                    vf_name = "# Observations"

                } else if(vf %in% c(".U", ".U[]")){
                    res_i = nrow(unique(x))
                    vf_name = "# Unique rows"

                } else {
                    # Other methods are faster but are less general

                    vf_new = gsub("(^\\.(N|U)\\[)|(\\]$)", "", vf)
                    do_unik = grepl("^\\.U", vf)

                    val = eval(str2lang(vf_new), x_list)

                    # we want to drop the NAs for indices
                    if(anyNA(val)){
                        val = val[!is.na(val)]
                    }

                    if(length(val) == 0){
                        res_i = 0

                    } else if(do_unik){
                        res_i = NROW(unique(x[val, ]))

                    } else {
                        res_i = if(is.logical(val)) sum(val) else length(val)
                    }

                    msg = if(do_unik) "# Unique rows" else "# Obs."
                    vf_new = gsub("NA_fun", "NA", vf_new, fixed = TRUE)
                    vf_name = paste0(msg, " with ", vf_new)
                }

            } else if(grepl("^(\\.NA|NA_fun\\(\\))", vf)) {

                if(vf %in% c(".NA", ".NA[]", "NA_fun()", "NA_fun()[]")){
                    res_i = sum(!complete.cases(x))
                    vf_name = "# Rows with NAs"

                } else {

                    subselect = NULL
                    sub_text = ""
                    if(grepl("\\]$", vf)){
                        sub_text_split = strsplit(vf, "[", fixed = TRUE)[[1]]
                        sub_text = paste0(sub_text_split[-1], collapse = "[")
                        sub_text = gsub("\\]$", "", sub_text)

                        subselect = eval(str2lang(sub_text), x_list)
                    }

                    if(length(subselect) == 0){
                        res_i = 0
                    } else {
                        res_i = sum(!complete.cases(x[subselect, ]))
                    }

                    vf_name = dsb("# Rows with NAs, with .[sub_text]")
                }
            } else {
                val = eval(vf_call, x_list)

                if(anyNA(val)){
                    who_NA = is.na(val)
                    na_i = sum(who_NA)
                    val = val[!who_NA]
                }

                if(grepl("^NA_fun\\(", vf)){
                    res_i = sum(val)

                    # now the message
                    leftover = sub("^[^\\)]+\\)", "", vf)
                    begin = substr(vf, 1, nchar(vf) - nchar(leftover))

                    my_fun = list("NA_fun" = function(x) enumerate_items(sapply(sys.call()[-1], deparse_long)))
                    na_vars = eval(str2lang(begin), my_fun)
                    msg_start = paste0("# NAs in ", na_vars)

                    msg_end = ""
                    if(grepl("[", vf, fixed = TRUE)){
                        msg_end = paste0(" with ", gsub("^\\[|\\]$", "", leftover))
                    }

                    vf_name = paste0(msg_start, msg_end)

                } else {
                    res_i = len_unique(val, nthreads)
                }

            }

            res[vf_name] = res_i
            na.info[i] = na_i
        }

        attr(res, "na.info") = na.info
        class(res) = "vec_n_unik"

        if(n_x == 1){
            return(res)
        }

        res_all[[x_all_names[I]]] = res
    }


    # Specific data set comparisons
    info_pairs = list()
    for(pair in comp_pairs){

        i_x = pair[1]
        i_y = pair[2]

        x = x_all[[i_x]]
        x_list = unclass(x)
        x_list[["NA_fun"]] = function(...) NA_fun(..., df = x)

        y = x_all[[i_y]]
        y_list = unclass(y)
        y_list[["NA_fun"]] = function(...) NA_fun(..., df = y)

        vars_legit = intersect(names(x), names(y))

        # two last elements: id and common
        all_rows_id_common = c()
        row_temp = rep(NA_real_, n_x + 2)

        for(i in seq_along(var_final)){

            vf = var_final[i]
            vf = gsub("((?<=[^[:alnum:]_\\.])|^)NA\\(", "NA_fun(", vf, perl = TRUE)
            vf_call = str2lang(vf)

            if(!all(all.vars(vf_call) %in% vars_legit) || grepl("^NA_fun", vf)){
                next

            } else {

                if(grepl("to_integer", vf, fixed = TRUE)){
                    # specific scheme for combinations

                    if(grepl("NA_fun()", vf, fixed = TRUE)){
                        # we don't allow that
                        next
                    }

                    vars_keep = all.vars(vf_call)

                    xy_list = list()
                    for(v in vars_keep){
                        xy_list[[v]] = c(x_list[[v]], y_list[[v]])
                    }

                    xy_list[["NA_fun"]] = function(...) NA_fun(..., df = stop("Internal error."))


                    val_xy = eval(vf_call, xy_list)

                    nrow_x = length(x_list[[1]])
                    nrow_y = length(y_list[[1]])
                    val_x = val_xy[1:nrow_x]
                    val_y = val_xy[(nrow_x + 1):(nrow_x + nrow_y)]

                } else {
                    val_x = eval(vf_call, x_list)
                    val_y = eval(vf_call, y_list)
                }


                if(anyNA(val_x)){
                    val_x = val_x[!is.na(val_x)]
                }

                if(anyNA(val_y)){
                    val_y = val_y[!is.na(val_y)]
                }

                # first col is ID
                val_x_unik = unique(val_x)
                val_y_unik = unique(val_y)

                n_x_not_in_y = sum(!val_x_unik %in% val_y_unik)
                n_y_not_in_x = sum(!val_y_unik %in% val_x_unik)
                n_common = sum(val_x_unik %in% val_y_unik)

                row = row_temp
                row[i_x] = n_x_not_in_y
                row[i_y] = n_y_not_in_x

                row[n_x + 1] = i
                row[n_x + 2] = n_common

                all_rows_id_common = rbind(all_rows_id_common, row)

            }
        }

        if(length(all_rows_id_common) > 0){
            attr(all_rows_id_common, "data_id") = c(i_x, i_y)
            info_pairs[[length(info_pairs) + 1]] = all_rows_id_common
        }
    }

    if(length(info_pairs) > 0){
        attr(res_all, "info_pairs") = info_pairs
    }

    class(res_all) = "list_n_unik"

    return(res_all)
}

#' @rdname n_unik
print.vec_n_unik = function(x, ...){

    hash = "## "

    x_names = sfill(names(x))
    na.info = attr(x, "na.info")
    na.info_format = sfill(fsignif(na.info))
    x_value = sfill(fsignif(x))

    n = length(x)
    na_col = paste0("(# NAs: ", na.info_format, ")")
    na_col[na.info == 0] = ""

    res = paste0(hash, x_names, ": ", x_value, " ", na_col)
    cat(res, sep = "\n")
}

#' @rdname n_unik
print.list_n_unik = function(x, ...){


    # I can't use #> anymore!!! The auto complete by the new version
    # of Rstudio drives me nuts!!!
    hash = "## "

    x_all = x
    n_x = length(x_all)

    if(n_x == 1) return(x_all[[1]])

    info_pairs = attr(x, "info_pairs")
    IS_PAIR = !is.null(info_pairs)

    # First, the variable names

    all_names = names(x_all[[1]])
    qui_0 = nchar(all_names) == 0
    i = 2
    while(any(qui_0) && i <= n_x){
        names_i = names(x_all[[i]])
        all_names[qui_0] = names_i[qui_0]
        qui_0 = nchar(all_names) == 0
        i = i + 1
    }

    # id_vars: used in info_pairs
    id_vars = seq_along(all_names)

    if(all(qui_0)){
        stop("Not any valid information to display: please check that all your variables exist in all data sets.")

    } else if(any(qui_0)){
        warning("Some variables could not be evaluated in any data set: please check that all your variables exist in all data sets.")

        all_names = all_names[!qui_0]
        id_vars = id_vars[!qui_0]

        for(i in 1:n_x){
            x = x_all[[i]]
            na.info = attr(x, "na.info")

            x = x[!qui_0]
            na.info = na.info[!qui_0]
            attr(x, "na.info") = na.info

            x_all[[i]] = x
        }
    }


    x_names = sfill(all_names)

    # If ambiguous pairs: we add data set suffixes
    if(n_x > 2 && IS_PAIR){
        add_suffix = function(x, i) paste0(x, " (", letters[i], ")")
        var_width = max(nchar("Exclusive "), nchar(x_names[1]))
    } else {
        add_suffix = function(x, i) x
        var_width = nchar(x_names[1])
    }

    var_col = paste0(hash, x_names, ":")
    na_intro = paste0(hash, sfill(" ", var_width), "|NAs:")
    var_col = insert_in_between(var_col, na_intro)
    # we add the first row of the data sets names + format
    var_col = sfill(c(hash, var_col), right = TRUE)

    print_mat = var_col

    KEEP_NA_ROW = rep(FALSE, length(x_names))

    for(i in 1:n_x){

        data_name = add_suffix(names(x_all)[i], i)

        x = x_all[[i]]

        na.info = attr(x, "na.info")
        KEEP_NA_ROW = KEEP_NA_ROW | na.info > 0
        na.info_format = fsignif(na.info)
        x_value = fsignif(x)
        x_value[is.na(x)] = "--"
        na.info_format[is.na(x)] = "--"


        width = max(nchar(na.info_format), nchar(x_value))

        na.info_format = sfill(na.info_format, width)
        x_value = sfill(x_value, width)

        x_col = insert_in_between(x_value, na.info_format)
        x_col = sfill(c(data_name, x_col))

        print_mat = cbind(print_mat, x_col)
    }

    if(!any(KEEP_NA_ROW)){
        print_mat[, 1] = substr(print_mat[, 1], 1, nchar(print_mat[1, 1]) - 4)
    }

    keep = c(TRUE, insert_in_between(TRUE, KEEP_NA_ROW))

    print_mat = print_mat[keep, ]

    # PAIR information
    if(IS_PAIR){

        # identifiers used for insertion
        id_vars_all = c(0, insert_in_between(id_vars, 0))
        id_vars_all = id_vars_all[keep]

        # we add two columns: id and common
        print_mat = cbind(print_mat, id_vars_all, "")

        insert_row_after_id = function(mat, row, col_id){

            for(i in 1:nrow(row)){
                i_id_mat = max(which(mat[, col_id] == row[i, col_id]))
                tmp_before = mat[1:i_id_mat, ]

                if(i_id_mat < nrow(mat)){
                    tmp_after = mat[(i_id_mat + 1):nrow(mat), ]
                } else {
                    tmp_after = NULL
                }

                mat = rbind(tmp_before, row[i, ], tmp_after)
            }

            mat
        }

        for(pair in info_pairs){

            data_id = attr(pair, "data_id")
            pair = as.data.frame(pair)

            if(n_x == 2){
                # data_col = paste0(hash, sfill(" ", var_width, right = TRUE), "|Excl:")
                data_col = paste0(hash, sfill("# Exclusive ", var_width, right = TRUE), "|")
            } else {
                data_col = paste0(hash, sfill("# Exclusive ", var_width, right = TRUE),
                                  "|", paste0(letters[data_id], collapse = ":"))
            }

            # formatting the NAs
            for(i in 1:(ncol(pair) - 2)){
                pair[[i]] = fsignif(pair[[i]])
            }
            pair[is.na(pair)] = " "
            pair = as.matrix(pair)

            # adding the data col
            pair = cbind(data_col, pair)

            print_mat = insert_row_after_id(print_mat, pair, n_x + 2)

        }

        # Formatting
        print_mat = print_mat[, -(n_x + 2)]
        print_mat[, 1] = sfill(print_mat[, 1], right = TRUE)
        for(j in 2:(n_x + 1)){
            print_mat[, j] = sfill(print_mat[, j])
        }

        # Last column
        common = print_mat[, n_x + 2]
        is_num = grepl("\\d", common)
        common[is_num] = sfill(fsignif(as.numeric(common[is_num])))
        common[is_num] = paste0("# Common: ", common[is_num])
        print_mat[, n_x + 2] = common
    }

    vec2print = apply(print_mat, 1, paste, collapse = " ")

    cat(vec2print, sep = "\n")

}


#' Formatted object size
#'
#' Tools that returns a formatted object size, where the appropriate unit is automatically chosen.
#'
#' @param x Any R object.
#' @param ... Not currently used.
#'
#' @return
#' Returns a character scalar.
#'
#' @examples
#'
#' osize(iris)
#'
#' data(trade)
#' osize(trade)
#'
#'
osize = function(x){

    size = as.numeric(utils::object.size(x))
    n = log10(size)

    if (n < 3) {
        res = paste0(size, " Octets.")
    } else if (n < 6) {
        res = paste0(fsignif(size/1000), " Ko.")
    } else {
        res = paste0(fsignif(size/1e+06), " Mo.")
    }

    class(res) = "osize"
    res
}

#' @rdname osize
print.osize = function(x, ...){
    cat(x, "\n")
}


#' Randomly draws observations from a data set
#'
#' This function is useful to check a data set. It gives a random number of rows of the input data set.
#'
#' @param x A data set: either a vector, a matrix or a data frame.
#' @param n The number of random rows/elements to sample randomly.
#' @param previous Logical scalar. Whether the results of the previous draw should be returned.
#'
#' @return
#' A data base (resp vector) with \code{n} rows (resp elements).
#'
#' @examples
#'
#' sample_df(iris)
#'
#' sample_df(iris, previous = TRUE)
#'
sample_df = function(x, n = 10, previous = FALSE){

    check_arg(n, "integer scalar")
    check_arg(previous, "logical scalar")

    if(MISSNULL(x)){
        if(missing(x)) stop("The argument 'x' cannot be missing, please provide it.")
        if(is.null(x)) stop("The argument 'x' must be a vector, matrix or data.frame. Currently it is NULL.")
    }

    all_draws = getOption("fixest_sample_df")
    x_dp = deparse(substitute(x))

    make_draw = TRUE
    if(previous){
        if(x_dp %in% names(all_draws)){
            draw__ = all_draws[[x_dp]]
            make_draw = FALSE
        } else {
            warning("No previous draw found for this data. Making a new draw.")
        }
    }

    is_unidim = is.null(dim(x)) || inherits(x, "table")
    n_data = if(is_unidim) length(x) else nrow(x)


    n = min(n, n_data)

    if(make_draw){
        # complicated var name to avoid issue with data.table variables
        draw__ = sample(n_data, n)
    }

    # saving
    all_draws[[x_dp]] = draw__
    options(fixest_sample_df = all_draws)

    # returning
    if(is_unidim){
        return(x[draw__])
    } else {
        return(x[draw__, ])
    }
}

len_unique = function(x, nthreads = getFixest_nthreads()){

    if(!is.vector(x)){
        return(length(unique(x)))
    }

    ANY_NA = FALSE
    if(is.numeric(x)){
        info_na = cpp_which_na_inf(x, nthreads)
        if(info_na$any_na_inf){
            ANY_NA = TRUE
            x = x[!info_na$is_na_inf]
        }
    } else {
        if(anyNA(x)){
            ANY_NA = TRUE
            x = x[!is.na(x)]
        }
    }

    if(length(x) == 0){
        return(0)
    }

    x_quf = quickUnclassFactor(x, addItem = TRUE, sorted = FALSE)

    length(x_quf$items) + ANY_NA
}

#' Replicates \code{fixest} objects
#'
#' Simple function that replicates \code{fixest} objects while (optionally) computing different standard-errors. Useful mostly in combination with \code{\link[fixest]{etable}} or \code{\link[fixest]{coefplot}}.
#'
#' @param x Either a \code{fixest} object, either a list of \code{fixest} objects created with \code{.l()}.
#' @param times Integer vector giving the number of repetitions of the vector of elements. By default \code{times = 1}. It must be either of length 1, either of the same length as the argument \code{x}.
#' @param each Integer scalar indicating the repetition of each element. Default is 1.
#' @param vcov A list containing the types of standard-error to be computed, default is missing. If not missing, it must be of the same length as \code{times}, \code{each}, or the final vector. Note that if the arguments \code{times} and \code{each} are missing, then \code{times} becomes equal to the length of \code{vcov}. To see how to summon a VCOV, see the dedicated section in the \href{https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1}{vignette}.
#' @param ... In \code{.l()}: \code{fixest} objects. In \code{rep()}: not currently used.
#'
#' @details
#' To apply \code{rep.fixest} on a list of \code{fixest} objects, it is absolutely necessary to use \code{.l()} and not \code{list()}.
#'
#' @return
#' Returns a list of the appropriate length. Each element of the list is a \code{fixest} object.
#'
#' @examples
#'
#' # Let's show results with different standard-errors
#'
#' est = feols(Ozone ~ Solar.R + Wind + Temp, data = airquality)
#'
#' my_vcov = list(~ Month, ~ Day, ~ Day + Month)
#'
#' etable(rep(est, vcov = my_vcov))
#'
#' coefplot(rep(est, vcov = my_vcov), drop = "Int")
#'
#' #
#' # To rep multiple objects, you need to use .l()
#' #
#'
#' est_bis = feols(Ozone ~ Solar.R + Wind + Temp | Month, airquality)
#'
#' etable(rep(.l(est, est_bis), vcov = my_vcov))
#'
#' # using each
#' etable(rep(.l(est, est_bis), each = 3, vcov = my_vcov))
#'
#'
rep.fixest = function(x, times = 1, each = 1, vcov, ...){
    # each is applied first, then times
    # x can be either a list of fixest objects, either a fixest object

    check_arg(x, "class(fixest, fixest_list) mbt")
    check_arg(times, "integer scalar GE{1} | integer vector no na GE{0}")
    check_arg(each, "integer scalar GE{1} | logical scalar")
    check_arg(vcov, "class(list)")

    if(is_user_level_call()){
        validate_dots(suggest_args = c("times", "each"), stop = TRUE)
    }

    # Checking the arguments
    IS_LIST = FALSE
    if("fixest_list" %in% class(x)){
        IS_LIST = TRUE
        class(x) = "list"

        n = length(x)

    } else {
        n = 1
    }

    if(is.logical(each) && each == FALSE){
        stop("Argument 'each' cannot be equal to FALSE.")
    }

    IS_MULTI_VCOV = !missing(vcov)
    if(IS_MULTI_VCOV){
        n_vcov = length(vcov)

        if(times == 1 && each == 1){
            if(isTRUE(each)){
                each = n_vcov
            } else {
                times = n_vcov
            }

        }
    }

    res_int = rep(1:n, times = times, each = each)
    n_res = length(res_int)

    if(IS_MULTI_VCOV){
        # Checking and expanding

        vcov_mapping = 1:n_res
        if(times == 1){
            if(n_vcov != each && n_vcov != n_res){
                stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", each, " or of length ", n_res, ".")
            }

            if(n_vcov == each) vcov_mapping = rep(1:each, times = n)

        } else if(each == 1){
            if(n_vcov != times && n_vcov != n_res){
                stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", times, " or of length ", n_res, ".")
            }

            if(n_vcov == times) vcov_mapping = rep(1:n_vcov, each = n)

        } else {
            if(n_vcov != n_res){
                stop("In rep, the argument 'vcov' (currently of length ", n_vcov, ") must be a list either of length ", n_res, ".")
            }
        }
    }

    res = vector("list", length(res_int))

    if(IS_MULTI_VCOV){
        for(i in 1:n_res){
            if(IS_LIST){
                res[[i]] = summary(x[[res_int[i]]], vcov = vcov[[vcov_mapping[i]]])
            } else {
                res[[i]] = summary(x, vcov = vcov[[vcov_mapping[i]]])
            }
        }

    } else {
        for(i in unique(res_int)){
            if(IS_LIST){
                res[res_int == i] = x[[i]]
            } else {
                res[res_int == i] = list(x)
            }
        }
    }

    for(i in 1:n_res){
        res[[i]]$model_id = res_int[i]
    }

    res
}

#' @rdname rep.fixest
rep.fixest_list = function(x, times = 1, each = 1, vcov, ...){
    rep.fixest(x, times = times, each = each, vcov = vcov, ...)
}

#' @rdname rep.fixest
.l = function(...){

    check_arg(..., "mbt class(fixest) | list")

    dots = list(...)
    if(all(sapply(dots, function(x) "fixest" %in% class(x)))){
        class(dots) = "fixest_list"

        return(dots)
    }

    if(length(dots) == 1){

        if("fixest_multi" %in% class(dots[[1]])){
            res = attr(dots[[1]], "data")
            class(res) = "fixest_list"
            return(res)
        }

        if(all(sapply(dots[[1]], function(x) "fixest" %in% class(x)))){
            res = dots[[1]]
            class(res) = "fixest_list"
            return(res)
        }
    }

    res = list()
    for(i in seq_along(dots)){
        if("fixest" %in% class(dots[[i]])){
            res[[length(res) + 1]] = dots[[i]]
        } else {
            obj = dots[[i]]

            if(class(obj) %in% "fixest_multi"){
                data = attr(obj, "data")
                for(j in seq_along(data)){
                    res[[length(res) + 1]] = data[[j]]
                }

            } else {
                for(j in seq_along(obj)){
                    if(!"fixest" %in% class(obj[[j]])){
                        stop("In .l(...), each argument must be either a fixest object, or a list of fixest objects. Problem: The ", n_th(j), " element of the ", n_th(i), " argument (the latter being a list) is not a fixest object.")
                    }

                    res[[length(res) + 1]] = obj[[j]]
                }
            }
        }
    }

    class(res) = "fixest_list"
    res
}


#### ................. ####
#### Internal Funs     ####
####

as.character.formula = function(x, ...) as.character.default(x, ...)

parse_macros = function(..., reset = FALSE, from_xpd = FALSE, check = TRUE, frame = NULL){
    set_up(1)

    if(check){
        check_arg(..., "dotnames os formula | character vector no na | numeric scalar | class(call, name)", .message = paste0("Each element of '...' must be a one-sided formula, and the name of each argument must start with two dots (ex: ", ifelse(from_xpd, "xpd(fml, ..ctrl = ~ x5 + x6)", "setFixest_fml(..ctrl = ~ x5 + x6)"), ").\nAlternatively it can be a character vector of variable names, or a numeric scalar."))
    }

    # Original macros
    fml_macro = getOption("fixest_fml_macro")
    if(reset || is.null(fml_macro)){
        fml_macro = list()
    } else if(!is.list(fml_macro)){
        warn_up("The value of getOption(\"fixest_fml_macro\") wasn't legal, it has been reset.")
        fml_macro = list()
        options(fixest_fml_macro = list())
    }

    # We check the names
    if(...length() == 0) return(fml_macro)

    dots = list(...)

    # Setting the new macros
    if(check){
        qui_pblm = !grepl("^\\.\\.", names(dots))
        if(any(qui_pblm)){
            arg_name = names(dots)[qui_pblm][1]
            correct_name = gsub("^\\.*", "\\.\\.", arg_name)
            stop_up("Each argument name must start with two dots. Use '", correct_name, "' instead of '", arg_name, "'.")
        }
    }

    for(v in names(dots)){
        fml_macro[[v]] = value2stringCall(dots[[v]], check = check, frame = frame)
    }

    fml_macro
}

value2stringCall = function(value_raw, call = FALSE, check = FALSE, frame = NULL){

    if(any(c("call", "name") %in% class(value_raw))){
        res = if(call) value_raw else deparse_long(value_raw)

    } else if(inherits(value_raw, "formula")){
        res = if(call) value_raw[[2]] else as.character(value_raw)[[2]]

    } else {

        if(check){
            value_raw = grep("[[:alnum:]]", value_raw, value = TRUE)
            if(length(value_raw)){
                # We need to check that it leads to a valid formula => otherwise problems later

                for(i in seq_along(value_raw)){
                    value_raw[i] = dot_square_bracket(value_raw[i], frame)
                }

                value_raw = paste(value_raw, collapse = " + ")

                my_call = error_sender(str2lang(value_raw), "The value '", value_raw, "' does not lead to a valid formula: ", up = 2)
                res = if(call) my_call else value_raw

            } else {
                res = if(call) 1 else "1"
            }

        } else {
            value_raw = value_raw[nzchar(value_raw)]
            if(length(value_raw)){
                value_raw = paste(value_raw, collapse = "+")
                res = if(call) str2lang(value_raw) else value_raw
            } else {
                res = if(call) 1 else "1"
            }
        }
    }

    if(check && "[" %in% all.vars(res, functions = TRUE)){
        res_txt = dot_square_bracket(deparse_long(res), frame)
        res = as.formula(res_txt)
    }

    res
}

dot_square_bracket = function(x, frame = .GlobalEnv, regex = FALSE, text = FALSE,
                              forced_merge = FALSE, up = 0){
    # transforms "x.[i]" into x1 if i==1
    # z = "XX" ; x = ".[z] + x.[1:5] + y.[1:2]_t"
    # x = "x.[a] | fe[b] + m.[o]"
    # .[,stuff] => means aggregation is comma based
    # if text, we allow:
    # - before.["text", stuff]after AND .[stuff, "text"]
    #   * .["text", stuff] => collapses stuff with "text"
    #     ie paste0(before, paste0(stuff, collapse = "text"), after)
    #   * .["text", stuff] => collapses stuff WITH the previous string with "text",
    #     ie paste0(paste0(before, stuff, collapse = "text), after)
    # regex = FALSE ; frame = .GlobalEnv ; text = FALSE
    # x = 'c(.[,y]) ~ sw(x.[,1:3])'
    # name = c("Jon", "Juliet") ; x = "bonjour .[' and ', name]" ; text = TRUE
    # name = c("Jon", "Juliet") ; x = "bonjour .[name, ' and ']" ; text = TRUE

    if(!grepl(".[", x, fixed = TRUE)) return(x)

    x_split_open = strsplit(x, ".[", fixed = TRUE)[[1]]

    # comma aggregation
    is_comma = grepl("^,", x_split_open[-1])
    if(any(is_comma) && !text){
        x_split_open[-1] = gsub("^, *", "", x_split_open[-1])
    }

    # nesting: a ~ .["x.[1:2]_sq"] + b
    any_nested = any(grepl("^\"", x_split_open[-1])) && !text
    if(any_nested){

        i_open_quote = setdiff(which(grepl("^\"", x_split_open)), 1)
        i_close_quote = which(grepl("\"\\]", x_split_open))

        x_split_new = x_split_open[1]
        for(i in i_open_quote){
            j = i_close_quote[i_close_quote >= i]

            xi_new_all = x_split_open[i:j]
            xi_new_all = gsub("\\]", "__close__", xi_new_all)

            xi_new = paste0(xi_new_all, collapse = "__open__")
            xi_new = gsub("\"__close__", "\"]", xi_new)

            x_split_new[[length(x_split_new) + 1]] = xi_new
        }

        n = length(x_split_open)
        if(j < n){
            for(i in (j+1):n){
                x_split_new[[length(x_split_new) + 1]] = x_split_open[i]
            }
        }

        x_split_open = x_split_new
    }

    # Finding the elements left/right of the closing bracket
    # we take care of indexing within brackets (e.g. x1 + .[m[1]])
    b_open = gregexpr("[", x_split_open[-1], fixed = TRUE)
    b_close = gregexpr("]", x_split_open[-1], fixed = TRUE)

    x_split = character(length(x_split_open) * 2 - 1)
    n = length(x_split)
    x_split[1] = x_split_open[1]

    for(i in seq_along(b_open)){

        if(b_open[[i]][1] != -1){
            # means there was an open [
            n_extend = length(b_close[[i]]) - length(b_open[[i]])
            index_closing = b_close[[i]][which.max(b_close[[i]] < c(b_open[[i]], rep(Inf, n_extend)))]
        } else {
            index_closing = b_close[[i]][1]
        }

        x_SO = x_split_open[i + 1]
        x_split_close_left = substr(x_SO, 1, index_closing - 1)
        x_split_close_right = substr(x_SO, index_closing + 1, nchar(x_SO))

        j = 2 * i
        x_split[[j]] = x_split_close_left
        x_split[[j + 1]] = x_split_close_right
    }

    if(any_nested){
        x_split = gsub("__open__", ".[", x_split)
        x_split = gsub("__close__", "]", x_split)
    }

    if(text){
        # NA means no operation
        operator = rep(NA_character_, n)
        do_split = rep(FALSE, n)
    }

    res = as.list(x_split)
    for(i in (1:n)[(1:n) %% 2 == 0]){

        # catching the operation
        do_lang = TRUE
        if(text){

            x_txt = trimws(x_split[i])

            op = NULL
            is_agg = is_split = FALSE
            if(grepl("^,", x_txt)){
                # default value of the aggregation
                op = ''
                x_txt = str_trim(x_txt, 1)
                is_agg = TRUE
            } else if(grepl("^/", x_txt)){
                # default value of the split
                op = "@, *"
                x_txt = str_trim(x_txt, 1)
                is_split = TRUE
            }

            if(!is_split && !is_agg){
                info_op = regexpr("^(('[^']*')|(\"[^\"]*\"))[,/]", x_txt)
                if(info_op != -1){
                    arg_pos = attr(info_op, "match.length")
                    op = substr(x_txt, 2, arg_pos - 2)
                    arg = substr(x_txt, arg_pos, arg_pos)
                    x_txt = str_trim(x_txt, arg_pos)

                    is_agg = arg == ","
                    is_split = !is_agg
                }
            }

            if(is_split){
                do_lang = FALSE

                operator[i] = op
                do_split[i] = TRUE
                my_call = x_txt

            } else if(is_agg){
                do_lang = FALSE

                x_txt = paste0("dsb_check_set_agg(", x_txt, ")")
                my_list = list(dsb_check_set_agg = dsb_check_set_agg)
                my_call = error_sender(eval(str2lang(x_txt), my_list, frame),
                                       "Dot square bracket operator: Evaluation of '.[",
                                       x_split[i], "]' led to an error:",
                                       up = up + 1)

                operator[i] = op
            }
        }


        if(do_lang){
            my_call = error_sender(str2lang(x_split[i]),
                                   "Dot square bracket operator: Evaluation of '.[",
                                   x_split[i], "]' led to an error:",
                                   up = up + 1)
        }

        if(is.character(my_call) && grepl(".[", my_call, fixed = TRUE)){
            # Nested call
            value = .dsb(my_call, frame = frame)

        } else {
            # Informative error message
            value = error_sender(as.character(eval(my_call, frame)),
                                 "Dot square bracket operator: Evaluation of '.[",
                                 x_split[i], "]' led to an error:",
                                 up = up + 1)
        }


        res[[i]] = value
    }

    if(any(is_comma)){
        is_comma = insert_in_between(FALSE, is_comma)
    } else if(!text){
        is_comma = rep(FALSE, length(res))
    }

    if(max(lengths(res)) == 1) {

        if(text && any(do_split)){
            res_txt = res[[1]]
            for(i in 2:length(res)){

                if(do_split[i]){
                    res_txt = paste0(res_txt, str_split(res[[i]], operator[i])[[1]])
                } else {
                    res_txt = paste0(res_txt, res[[i]])
                }
            }

            res = res_txt

        } else {
            res = paste(res, collapse = "")
        }

    } else {
        # first value is NEVER a vector ("" is added automatically in split)
        res_txt = res[[1]]
        i = 2
        while(i <= n){

            if(length(res[[i]]) == 1){

                if(text && do_split[i]){
                    res_txt = paste0(res_txt, str_split(res[[i]], operator[i])[[1]])
                } else {
                    res_txt = paste0(res_txt, res[[i]])
                }

                i = i + 1

            } else if(text){
                if(is.na(operator[i])){
                    if(forced_merge){
                        res_txt = paste0(res_txt, paste0(res[[i]], collapse = "_MERGE_"))
                    } else {
                        res_txt = paste0(res_txt, res[[i]])
                    }

                } else {
                    # If we're here => must be an aggregation requested
                    # (if it was a split, it would be of length 1)
                    res_txt = paste0(res_txt, paste0(res[[i]], collapse = operator[i]))
                }

                i = i + 1

            } else if(regex) {
                after = if(i != n) res[[i + 1]] else ""
                res_txt = paste0(res_txt, res[[i]], after, collapse = "|")
                i = i + 2

            } else {
                before_no_var = gsub("[[:alnum:]_\\.]+$", "", res_txt)
                var_before = substr(res_txt, nchar(before_no_var) + 1, nchar(res_txt))

                coll = if(is_comma[i]) ", " else " + "

                if(i != n && grepl("[[:alnum:]_\\.]", substr(res[[i + 1]], 1, 1))){
                    after = res[[i + 1]]
                    after_no_var = gsub("^[[:alnum:]_\\.]+", "", after)
                    var_after = substr(after, 1, nchar(after) - nchar(after_no_var))

                    res_txt = paste0(before_no_var, paste0(var_before, res[[i]], var_after, collapse = coll), after_no_var)
                    i = i + 2
                } else {
                    res_txt = paste0(before_no_var, paste0(var_before, res[[i]], collapse = coll))
                    i = i + 1
                }
            }
        }

        res = res_txt
    }


    # Recursivity prevention: no recursivity if text = TRUE
    DSB_RECURSIVE = TRUE
    is_rec = exists("DSB_RECURSIVE", parent.frame(), inherits = FALSE)

    if(!text && !is_rec && grepl(".[", res, fixed = TRUE)){
        res = dot_square_bracket(res, frame, regex)
    }

    res
}


dsb_check_set_agg = function(...){
    # Here we should have only one element, everything has been cleaned up before

    if(...length() > 1){
        if(...length() == 2){
            stop_up("The operator .[] accepts only up to two elements, if so the first one MUST be a string literal (eg: .['text', y]). Problem: it is not a string literal.")
        }
        stop_up("You cannot have more than two elements in between .[]. Currently there are ", ...length() + 1, ".")
    }

    mc = match.call()

    return(mc[[2]])
}

print_coeftable = function(coeftable, lastLine = "", show_signif = TRUE){
    # Simple function that does as the function coeftable but handles special cases
    # => to take care of the case when the coefficient is bounded

    if(!is.data.frame(coeftable)){
        class(coeftable) = NULL
        ct = as.data.frame(coeftable)
    } else {
        ct = coeftable
    }

    signifCode = c("***"=0.001, "** "=0.01, "*  "=0.05, ".  "=0.1)

    pvalues = ct[, 4]

    stars = cut(pvalues, breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
    stars[is.na(stars)] = ""

    whoIsLow = !is.na(pvalues) & pvalues < 2.2e-16

    # Note that it's a bit different than format => I don't like xxe-yy numbers, very hard to read: you can't see large/small nbers at first sight
    for(i in 1:3){
        ct[, i] = decimalFormat(ct[, i])
    }

    ct[!whoIsLow, 4] = format(ct[!whoIsLow, 4], digits = 5)

    ct[whoIsLow, 4] = "< 2.2e-16"
    ct[is.na(ct[, 4]), 4] = "NA"

    ct[, 5] = stars
    names(ct)[5] = ""

    print(ct)

    cat(lastLine)

    if(show_signif){
        cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    }
}

plot_single_cluster = function(x, n=5, addExp = FALSE, fe_name, ...){
    # It plots the n first and last most notable FEs

    isSlope = grepl("\\[", fe_name)

    # we compare with the average of the coefficients
    if(isSlope == FALSE) x = sort(x) - mean(x)

    x_name = names(x)

    k = length(x)
    nb_show = min(k, 2*n+1)
    mid_value = floor(nb_show/2)
    xlim = c(1, nb_show)
    xlim = xlim + c(-1, 1) * diff(xlim)/30

    if(k <= 2*n+1){
        # no need of space to split the data
        x_values = 1:k
        y_values = x
        isSplit = FALSE
    } else {
        nb_show = nb_show - 1 # because we don't want to show the middle point
        x_values = c(1:n, (n+2):(2*n+1))
        y_values = c(head(x, n), tail(x, n))
        isSplit = TRUE
    }

    # very specific case where the axis confonds with the boxing
    ylim = range(y_values)
    if(pretty(range(y_values))[1] < min(y_values) || tail(pretty(range(y_values)), 1) > max(y_values)){
        ylim = range(y_values) + c(-1, 1) * diff(range(y_values))/30
    }

    plot(x_values, y_values, ann = FALSE, axes = FALSE, xlim=xlim, ylim = ylim, col = 0)

    # display
    box()
    y = axis(2)
    abline(h = y, lty = 4, col = "lightgrey")

    # name information & points
    points(x_values, y_values)
    text(head(x_values, mid_value), head(y_values, mid_value), head(x_name, mid_value), pos = 3)
    text(tail(x_values, nb_show-mid_value), tail(y_values, nb_show-mid_value), tail(x_name, nb_show-mid_value), pos = 1)

    if(isSlope){
        fe = gsub("\\[.+", "", fe_name)
        slope = gsub(".+\\[|\\].+", "", fe_name)
        title(xlab = paste0("Slope Coefficients (", fe, ")"), line = 1)
        title(main = slope)
    } else {
        title(xlab = "Centered Fixed-Effects", line = 1)
        title(main = fe_name)
    }


    if(isSplit){
        abline(v = c(n+0.75, n+1.25), lty = 2)

        axis(1, at = c(n+0.75, n+1.25), col = "grey", labels = NA, lwd.ticks = 0)
        axis(3, at = c(n+0.75, n+1.25), col = "grey", labels = NA, lwd.ticks = 0)
    }

    coord = par("usr")
    # axis(1, at = coord[1:2], labels = c("coef", "exp(coef)"), lwd = 0, lwd.ticks = 1)

    if(addExp && !isSlope){
        axis(1, at = coord[1], labels = "coef", lwd = 0, lwd.ticks = 1, line = -1)
        axis(1, at = coord[2], labels = "exp(coef)", lwd = 0, lwd.ticks = 1, line = -1)
        axis(4, at=y, labels = signif(exp(y), 2))
    }

}

prepare_matrix = function(fml, base, fake_intercept = FALSE){
    # This function is way faster than model.matrix but does not accept factors
    # The argument fml **MUST** not have factors!

    rhs = if(length(fml) == 3) fml[c(1, 3)] else fml

    t = terms(rhs, data = base)

    all_var_names = attr(t, "term.labels")

    # We take care of interactions: references can be multiple, then ':' is legal
    all_vars = cpp_colon_to_star(all_var_names)

    # Forming the call
    if(attr(t, "intercept") == 1 && !fake_intercept){
        n = nrow(base)
        all_vars_call = str2lang(paste0("list('(Intercept)' = rep(1, ", n, "), ", paste0(all_vars, collapse = ", "), ")"))
        all_var_names = c("(Intercept)", all_var_names)
    } else {
        all_vars_call = str2lang(paste0("list(", paste0(all_vars, collapse = ", "), ")"))
    }

    # evaluation
    data_list = eval(all_vars_call, base)

    # Handling the multi columns case (ex: bs(x1), splines)
    # NOTA: I need to add a check for i() because of 1 value interactions
    #       not caught by the != nber of obs
    qui_i = grepl("\\bi\\(", all_var_names)
    if(any(lengths(data_list) != nrow(base)) || any(qui_i)){

        all_n = as.vector(lengths(data_list) / nrow(base))

        qui_pblm = which(all_n %% 1 != 0)
        if(length(qui_pblm) > 0){
            what = data_list[[qui_pblm]]
            reason = ifelse(is.null(nrow(what)), paste0("of length ", length(what)), paste0("with ", nrow(what), " rows"))

            stop("Evaluation of ", all_var_names[qui_pblm], " returns an object ", reason, " while the data set has ", nrow(base)," rows.", call. = FALSE)
        }

        all_n_vector = rep(all_n, all_n)

        new_names = as.list(all_var_names)
        for(i in which(all_n > 1 | qui_i)){
            my_names = colnames(data_list[[i]])
            if(is.null(my_names)){
                my_names = 1:all_n[i]
            }
            new_names[[i]] = paste0(all_var_names[i], my_names)
        }

        all_var_names = unlist(new_names)
    }

    res = do.call("cbind", data_list)

    colnames(res) = all_var_names

    res
}


fixest_model_matrix = function(fml, data, fake_intercept = FALSE, i_noref = FALSE, mf = NULL){
    # This functions takes in the formula of the linear part and the
    # data
    # It reformulates the formula (ie with lags and interactions)
    # then either apply a model.matrix
    # either applies an evaluation (which can be faster)
    #
    # fake_intercept => whether to add the intercept, only to make sure
    #  the factors are well created

    # fml = ~a*b+c+i(x1)+Temp:i(x2)+i(x3)/Wind

    if(length(fml) == 3) fml = fml[c(1, 3)]

    # We need to check a^b otherwise error is thrown in terms()
    rhs_txt = deparse_long(fml[[2]])
    if(grepl("\\^[[:alpha:]]", rhs_txt)){
        stop("The special operator '^' can only be used in the fixed-effects part of the formula. Please use ':' instead.")
    }

    #
    # Evaluation
    #

    t_fml = terms(fml)
    tl = attr(t_fml, "term.labels")

    if(length(tl) == 0){

        if(fake_intercept){
            return(1)
        }

        res = matrix(1, nrow = nrow(data), ncol = 1, dimnames = list(NULL, "(Intercept)"))
        return(res)
    }

    # We check for calls to i()
    qui_i = grepl("(^|[^[:alnum:]_\\.])i\\(", tl)
    IS_I = any(qui_i)
    if(IS_I){
        # OMG... why do I always have to reinvent the wheel???
        is_intercept = fake_intercept || (attr(t_fml,"intercept") == 1)
        i_naked = which(is_naked_fun(tl[qui_i], "i"))

        if(i_noref){
            for(i in seq_along(i_naked)){
                j = i_naked[i]
                txt = gsub("(^|(?<=[^[:alnum:]\\._]))i\\(", "i_noref(", tl[qui_i][j], perl = TRUE)
                tl[qui_i][j] = eval(str2lang(txt))
            }
        } else {
            for(i in seq_along(i_naked)){
                if(!is_intercept && i == 1) next

                j = i_naked[i]
                txt = gsub("(^|(?<=[^[:alnum:]\\._]))i\\(", "i_ref(", tl[qui_i][j], perl = TRUE)
                tl[qui_i][j] = eval(str2lang(txt))
            }
        }

        fml_no_inter = .xpd(rhs = tl[!qui_i])

        if(!is_intercept) tl = c("-1", tl)
        fml = .xpd(rhs = tl)

    }

    # Are there factors NOT in i()? If so => model.matrix is used
    dataNames = names(data)

    if(!is.null(mf)){
        useModel.matrix = TRUE
    } else {
        useModel.matrix = FALSE
        if(IS_I){
            linear.varnames = all.vars(fml_no_inter[[2]])
            is_num = sapply(data[, dataNames %in% linear.varnames, FALSE], is.numeric)
            if(length(is_num) > 0 && (any(!is_num) || grepl("factor", deparse_long(fml_no_inter)))){
                useModel.matrix = TRUE
            }
        } else {
            linear.varnames = all.vars(fml[[2]])
            is_num = sapply(data[, dataNames %in% linear.varnames, FALSE], is.numeric)
            if(length(is_num) == 0 || any(!is_num) || grepl("factor", deparse_long(fml))){
                useModel.matrix = TRUE
            }
        }
    }

    if(useModel.matrix){
        # to catch the NAs, model.frame needs to be used....
        if(is.null(mf)){
            mf = stats::model.frame(fml, data, na.action = na.pass)
        }

        linear.mat = stats::model.matrix(fml, mf)

        if(fake_intercept){
            who_int = which("(Intercept)" %in% colnames(linear.mat))
            if(length(who_int) > 0){
                linear.mat = linear.mat[, -who_int, drop = FALSE]
            }
        }
    } else {
        linear.mat = prepare_matrix(fml, data, fake_intercept)
    }

    if(any(grepl("__CLEAN__", colnames(linear.mat), fixed = TRUE))){
        new_names = clean_interact_names(colnames(linear.mat))

        colnames(linear.mat) = new_names
    }

    if(is.integer(linear.mat)){
        linear.mat = 1 * linear.mat
    }

    linear.mat
}


fixest_model_matrix_extra = function(object, newdata, original_data, fml, fake_intercept = FALSE, i_noref = FALSE, subset = FALSE){
    # Only used within model.matrix and predict
    # Overlay of fixest_model_matrix to take care of special things, eg:
    # - poly
    # - subset
    # - ?

    #
    # subset
    #

    # subset => we allow the extraction of only some variables
    all_vars = all_vars_with_i_prefix(fml[[3]])
    if(isFALSE(subset)){

        if(!original_data && any(!all_vars %in% names(newdata))){
            pblm = setdiff(all_vars, names(newdata))
            stop("In 'model.matrix', the variable", enumerate_items(pblm, "is.s.quote"), " in the formula but not in the argument 'data'. Use 'subset = TRUE' to enable the creation of partial data.")
        }

    } else {
        vars_keep = names(newdata)

        if(is.character(subset)){
            # ex: subset = c("x1$", "x2$")

            vars_keep = keep_apply(vars_keep, subset)
            if(length(vars_keep) == 0){
                stop("The variables in 'subset' do not match any variable in the 'data'.")
            }

            if(isFALSE(keep_apply("(Intercept)", subset, logical = TRUE))){
                fake_intercept = TRUE
            }

        } else {
            # intercept always removed if subset = TRUE!!!
            # that's the point of subset.

            fake_intercept = TRUE
        }

        if(!all(all_vars %in% vars_keep)){

            terms_all = attr(terms(fml), "term.labels")
            # We first check pure variables (because non pure variables are slower to check)
            is_var = grepl("^[\\.[:alpha:]][[:alnum:]\\._]*$", terms_all)
            terms_drop = is_var & !terms_all %in% vars_keep

            for(i in which(!is_var)){
                if(any(!all_vars_with_i_prefix(str2lang(terms_all[i])) %in% vars_keep)){
                    terms_drop[i] = TRUE
                }
            }

            if(all(terms_drop)){
                stop("Due to the use of the argument 'subset', not a single variable is left.")
            }

            fml = .xpd(rhs = terms_all[!terms_drop])
        }
    }

    #
    # Extra functions that need raw data-evaluation + single valued factors
    #

    mf = NULL
    if(!original_data){

        # What I don't like here is:
        # to take care of only two particular cases, namely
        #  + functions using the original data
        #  + single value factors
        # we incur a lot of computing time to EVERYONE
        # That's really not good.
        # The problem is that it's not so easy to catch those cases 100% without error.
        # I shall find a solution at some point.

        # if lean = TRUE, we should be avoiding that
        # => I don't know of a solution yet...

        if(length(fml) == 3) fml = fml[c(1, 3)]

        # We apply model.frame to the original data
        data = fetch_data(object, "To apply 'model.matrix.fixest', ")

        panel__meta__info = set_panel_meta_info(object, data)

        mf = model.frame(fml, data, na.action = na.pass)

        rm(panel__meta__info) # needed, so that when the data is recreated for real

        t_mf = terms(mf)
        xlev = .getXlevels(t_mf, mf)

        if(!identical(attr(t_mf,"variables"), attr(t_mf,"predvars")) || length(xlev) > 0){
            mf = model.frame(t_mf, newdata, xlev = xlev, na.action = na.pass)
        } else {
            mf = NULL
        }
    }

    GLOBAL_fixest_mm_info = list()

    new_matrix = fixest_model_matrix(fml, newdata, fake_intercept, i_noref, mf = mf)

    if(length(GLOBAL_fixest_mm_info) > 0){
        attr(new_matrix, "model_matrix_info") = GLOBAL_fixest_mm_info
    }

    new_matrix
}





fixef_terms = function(fml, stepwise = FALSE, origin_type = "feols"){
    # separate all terms of fml into fixed effects and varying slopes
    # fml: one sided formula
    # can also be a vector of terms
    # Je fais tout ce tralala a cause de ce putain de terms() qui ne marche pas avec a^b!!!! fait chier!

    # fixef_terms(~dum_1[[x1]] + dum_2[x2])
    # fixef_terms(~dum_1[[x1]] + dum_2 + dum_3[x2, x3])

    if(!is.vector(fml)){
        # we need to take care of the ^ used to combine variables
        fml_char = as.character(fml[2])
        if(grepl("^", fml_char, fixed = TRUE)){
            fml_char_new = gsub("^", "%^%", fml_char, fixed = TRUE)
            fml = as.formula(paste0("~", fml_char_new))
        }

        if(!origin_type %in% c("feols", "feglm") && grepl("[", fml_char, fixed = TRUE)){
            stop("The use of varying slopes is available only for the functions feols, feglm or fepois.")
        }

        t = terms(fml)

        if(stepwise){
            sw_info = extract_stepwise(tms = t)
            return(sw_info)
        }

        my_vars = attr(t, "term.labels")

    } else {
        my_vars = fml
    }

    # Further error checking (misuse of varying slopes)
    if(any(grepl("]", my_vars, fixed = TRUE))){
        var2check = my_vars[grepl("]", my_vars, fixed = TRUE)]
        var2check_double = var2check[grepl("]]", var2check, fixed = TRUE)]
        var2check_single = var2check[!grepl("]]", var2check, fixed = TRUE)]
        if(length(var2check_double) > 0){
            qui = !grepl("\\]\\]$", var2check_double) | !grepl("\\[\\[", var2check_double) | lengths(strsplit(var2check_double, "\\[")) != 3
            if(any(qui)){
                item_pblm = var2check_double[qui]
                msg = paste0("Square bracket are special characters used **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm), ").")
                class(msg) = "try-error"
                return(msg)
            }
        }

        if(length(var2check_single) > 0){
            qui = !grepl("\\]$", var2check_single) | lengths(strsplit(var2check_single, "\\[")) != 2
            if(any(qui)){
                item_pblm = var2check_single[qui]
                msg = paste0("Square bracket are special characters used **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm), ").")
                class(msg) = "try-error"
                return(msg)
            }
        }
    }

    # And yet again some error checking => i() should NOT be used
    if(any(grepl("^i\\(", my_vars))){
        # We create an error instead of simply correcting the syntax => this is because the function i is
        # very different and should not be confused

        var_pblm = my_vars[grepl("^i\\(", my_vars)][1]

        get_new_var = function(var, f, f2, ...) match.call()

        what = eval(str2lang(gsub("^i", "get_new_var", var_pblm)))
        n_var = sum(c("var", "f", "f2") %in% names(what))
        msg = if(n_var == 1) "Using i() to create fixed-effects is not possible, use directly the variable." else paste0("To interact fixed-effects, use the syntax fe1^fe2 (in your case ", deparse(what[[2]]), "^", deparse(what[[3]]), ").")

        stop("The function i() should not be used in the fixed-effects part of the formula. ", msg)
    }

    # Internal function
    reformulate_varslope = function(dum, ..., add_dum = TRUE){
        mc = match.call()
        dum = deparse_long(mc[["dum"]])

        if(add_dum){
            res = dum
        } else {
            res = c()
        }

        for(i in 3:(length(mc) - !add_dum)){
            value = deparse_long(mc[[i]])
            res = c(res, paste0(dum, "[[", value, "]]"))
        }
        res
    }

    qui_slope = grepl("[", my_vars, fixed = TRUE)
    if(any(qui_slope)){
        new_terms = c()

        vars_slope = my_vars[qui_slope]
        for(i in seq_along(my_vars)){
            if(qui_slope[i]) {
                v = my_vars[i]
                v_mod = paste0("reformulate_varslope(", gsub("\\[+", ", ", v))
                add_dum = ifelse(grepl("\\[\\[", v), ", add_dum = FALSE", "")
                v_mod = gsub("\\]+", paste0(add_dum, ")"), v_mod)

                new_v = eval(str2lang(v_mod))
                new_terms = c(new_terms, new_v)
            } else {
                new_terms = c(new_terms, my_vars[i])
            }
        }

        my_vars = new_terms
    }

    # we return fml_terms, fe_vars, slope_vars, slope_flag

    # we need to unique it
    my_vars = unique(my_vars)

    # we put the ^ back
    my_vars = gsub(" %^% ", "^", my_vars, fixed = TRUE)

    res = list(fml_terms = my_vars)
    # OLD version
    fe_vars_all = gsub("\\[.+", "", my_vars)
    is_slope_all = grepl("\\[", my_vars)
    slope_vars_all = rep(NA_character_, length(my_vars))
    slope_vars_all[is_slope_all] = gsub(".+\\[|\\]", "", my_vars[is_slope_all])

    # New version following the reworking of demeaning
    fe_vars = unique(fe_vars_all)
    #fe_vars: unique of the fixed-effects variables
    slope_flag = rep(0, length(fe_vars))
    # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect

    slope_vars_list = vector("list", length(fe_vars))
    names(slope_vars_list) = fe_vars
    for(i in seq_along(fe_vars)){
        qui = which(fe_vars_all == fe_vars[i])
        if(any(is_slope_all[qui])){
            nb = sum(is_slope_all[qui])
            slope_flag[i] = nb * (1 - 2*(length(qui) == nb))
            slope_vars_list[[i]] = slope_vars_all[qui][is_slope_all[qui]]
        }
    }

    res$fe_vars = fe_vars
    res$slope_flag = as.integer(slope_flag)
    res$slope_vars_list = slope_vars_list
    res$slope_vars = unlist(slope_vars_list, use.names = FALSE)

    res
}

prepare_df = function(vars, base, fastCombine = NA){
    # vars: vector of variables to evaluate

    # we drop NAs and make it unique
    vars = unique(vars[!is.na(vars)])
    all_var_names = vars

    do_combine = !is.na(fastCombine)

    changeNames = FALSE
    if(do_combine && any(grepl("^", vars, fixed = TRUE))){
        # special indicator to combine factors
        # ^ is a special character: only used to combine variables!!!

        vars = fml_combine(vars, fastCombine, vars = TRUE)

        changeNames = TRUE
    }

    all_vars = cpp_colon_to_star(vars)

    if(all(all_vars %in% names(base))){
        res = base[, all_vars, drop = FALSE]
    } else {
        all_vars_call = str2lang(paste0("list(", paste0(all_vars, collapse = ", "), ")"))
        data_list = try(eval(all_vars_call, base))

        # if error: we send it back to the main function
        if("try-error" %in% class(data_list)){
            return(data_list)
        }

        names(data_list) = all_var_names
        data_list$stringsAsFactors = FALSE

        res = do.call("data.frame", data_list)

        if(changeNames){
            all_var_names = rename_hat(all_var_names)
        }

        names(res) = all_var_names
    }

    res
}


# fml_char = "x + y + u^factor(v1, v2) + x5"
fml_combine = function(fml_char, fastCombine, vars = FALSE){
    # function that transforms "hat" interactions into a proper function call:
    # Origin^Destination^Product + Year becomes ~combine_clusters(Origin, Destination, Product) + Year

    fun2combine = ifelse(fastCombine, "combine_clusters_fast", "combine_clusters")

    # we need to change ^ into %^% otherwise terms sends error
    labels = attr(terms(.xpd(rhs = gsub("\\^(?=[^0-9])", "%^%", fml_char, perl = TRUE))), "term.labels")

    # now we work this out
    for(i in seq_along(labels)){
        lab = labels[i]
        if(grepl("^", lab, fixed = TRUE)){
            lab_split = trimws(strsplit(lab, "%^%", fixed = TRUE)[[1]])
            if(grepl("(", lab, fixed = TRUE)){
                # we add some error control -- imperfect, but... it's enough
                lab_collapsed = gsub("\\([^\\)]+\\)", "", lab)
                if(length(lab_split) != length(strsplit(lab_collapsed, "%^%", fixed = TRUE)[[1]])){
                    msg = "Wrong formatting of the fixed-effects interactions. The '^' operator should not be within parentheses."
                    stop(msg)
                }
            }
            labels[i] = paste0(fun2combine, "(", paste0(lab_split, collapse = ", "), ")")
        }
    }

    if(vars){
        return(labels)
    }

    fml = .xpd(rhs = labels)

    fml
}

# x = c('combine_clusters(bin(fe1, "!bin::2"), fe2)', 'fe3')
rename_hat = function(x){

    qui = grepl("combine_clusters", x, fixed = TRUE)
    if(!any(qui)) return(x)

    for(i in which(qui)){
        xi_new = gsub("combine_clusters(_fast)?", "sw", x[i])
        sw_only = extract_fun(xi_new, "sw")
        sw_eval = eval(str2lang(sw_only$fun))
        new_var = paste0(sw_only$before, paste0(sw_eval, collapse = "^"), sw_only$after)
        x[i] = new_var
    }

    x
}

prepare_cluster_mat = function(fml, base, fastCombine){
    # prepares the data.frame of the cluster variables

    fml_char = as.character(fml[2])
    changeNames = FALSE
    if(grepl("^", fml_char, fixed = TRUE)){
        # special indicator to combine factors
        fml = fml_combine(fml_char, fastCombine)
        changeNames = TRUE
    }

    t = terms(fml, data = base)

    all_var_names = attr(t, "term.labels")
    all_vars = gsub(":", "*", all_var_names)

    if(all(all_vars %in% names(base))){
        res = base[, all_vars, drop = FALSE]
    } else {
        all_vars_call = str2lang(paste0("list(", paste0(all_vars, collapse = ", "), ")"))
        data_list = eval(all_vars_call, base)
        names(data_list) = all_var_names
        data_list$stringsAsFactors = FALSE

        res = do.call("data.frame", data_list)

        if(changeNames){
            all_var_names = rename_hat(all_var_names)
        }

        names(res) = all_var_names
    }

    res
}

combine_clusters_fast = function(...){
    # This functions creates a new cluster from several clusters
    # basically: paste(cluster1, cluster2, ... etc, sep = "__")
    # but it tries to be faster than that because paste is very slow on large datasets

    cluster = list(...)
    Q = length(cluster)

    # The problem is that the clusters can be of ANY type...
    # So too much of a pain to take care of them in c++
    # => when I have time I'll do it, but pfff...

    # Not super efficient, but that's life
    ANY_NA = FALSE
    if(any(who_NA <- sapply(cluster, anyNA))){
        ANY_NA = TRUE
        who_NA = which(who_NA)
        IS_NA = is.na(cluster[[who_NA[1]]])
        for(i in who_NA[-1]){
            IS_NA = IS_NA | is.na(cluster[[i]])
        }

        # Nice error message comes later
        if(all(IS_NA)) return(rep(NA, length(IS_NA)))

        # we recreate the clusters
        for(i in 1:Q){
            cluster[[i]] = cluster[[i]][!IS_NA]
        }
    }

    # First we unclass
    for(i in 1:Q){
        cluster[[i]] = quickUnclassFactor(cluster[[i]])
    }

    # Then we combine
    power = floor(1 + log10(sapply(cluster, max)))

    if(sum(power) > 14){
        order_index = do.call(order, cluster)
        index = cpp_combine_clusters(cluster, order_index)
    } else {
        # quicker, but limited by the precision of doubles
        index = cluster[[1]]
        for(q in 2:Q){
            index = index + cluster[[q]]*10**sum(power[1:(q-1)])
        }
    }

    if(ANY_NA){
        # we recreate the return vector with appropriate NAs
        res = rep(NA_real_, length(IS_NA))
        res[!IS_NA] = index
    } else {
        res = index
    }

    return(res)
}

combine_clusters = function(...){
    # This functions creates a new cluster from several clusters
    # basically: paste(cluster1, cluster2, ... etc, sep = "_")
    # No, it's mow much more fatser than that!!!! Thanks for my new algo
    # => the paste only applies to the unique number of items


    clusters = to_integer(..., add_items = TRUE, items.list = TRUE, internal = TRUE)

    res = clusters$items[clusters$x]

    return(res)
}

listDefault = function(x, variable, value){
    # This function puts 'value' into the element 'variable' of list 'x'
    # IF it does not already exists in 'x'

    x_name = deparse(substitute(x))

    if(is.null(x[[variable]])){
        x[[variable]] = value

        assign(x_name, x, envir = parent.frame(n = 1))
    }

}

hgrid = function(lty = 3, col = "darkgray", ymin = -Inf, ymax = Inf, ...){
    # simple function that draws an horizontal grid

    # Finding the coordinates
    y = axis(2, lwd=0, labels = NA)

    y = y[y > ymin & y < ymax]

    # now drawing the lines
    if(length(y) > 0){
        abline(h = y, col = col, lty = lty, ...)
    }
}

vgrid = function(lty = 3, col = "darkgray", xmin = -Inf, xmax = Inf, ...){
    # simple function that draws an horizontal grid

    # Finding the coordinates
    x = axis(1, lwd=0, labels = NA)

    x = x[x > xmin & x < xmax]

    # now drawing the lines
    if(length(x) > 0){
        abline(v = x, col = col, lty = lty, ...)
    }
}

shade_area = function(y1, y2, x, xmin, xmax, col="grey", ...){
    # fonction plus pratique que polygon
    # elle permet de griser une partie delimitee par
    # y1 et y2 pour chacune des valeurs de x
    # on doit avoir la meme longueur de y1,y2 et x
    # exemple:
    # a=curve(x**2,-5,5)
    # shade_area(a$y+1,a$y-1,a$x)
    # qqes parametres graphiques:
    # lwd / border (couleur du bord, peut etre NA) / lty

    n = length(x)
    stopifnot(length(y1)==n | length(y1)==1)
    stopifnot(length(y2)==n | length(y2)==1)

    if(length(y1)==1) y1 = rep(y1,n)
    if(length(y2)==1) y2 = rep(y2,n)

    if(missing(xmin)) xmin = min(x)
    if(missing(xmax)) xmax = max(x)

    ind = which(x>=xmin & x<=xmax)
    x1 = x[ind] ; x2 = x[rev(ind)]
    polygon(c(x1,x2), c(y1[ind], y2[rev(ind)]), col=col, ...)
}


check_set_nthreads = function(nthreads){
    # Simple function that checks that the nber of threads is valid
    set_up(1)

    check_value(nthreads, "integer scalar GE{0} | numeric scalar GT{0} LT{1}", .message = paste0("The argument 'nthreads' must be an integer lower or equal to the number of threads available (", max(cpp_get_nb_threads(), 1), "). It can be equal to 0 which means all threads. Alternatively, if equal to a number strictly between 0 and 1, it represents the fraction of all threads to be used."))

    max_threads = cpp_get_nb_threads()

    # # To add later
    # if(cpp_is_in_fork()) return(1)

    if(nthreads == 0){
        nthreads = max(max_threads, 1)

    } else if(nthreads < 1){
        nthreads = max(ceiling(max_threads * nthreads), 1)

    } else if(nthreads > 1){
        if(max_threads == 0){
            warn_up("OpenMP not detected: cannot use ", nthreads, " threads, single-threaded mode instead.")
            nthreads = 1
        } else if(nthreads > max_threads){
            warn_up("Asked for ", nthreads, " threads while the maximum is ", max_threads, ". Set to ", max_threads, " threads instead.")
            nthreads = max_threads
        }

    }

    nthreads
}

clean_interact_names = function(x){
    # GOD, why is it so complicated? Why?
    # Why do I have to spend so much time on that crap?
    #
    # x = c("i(Month, keep = 4:7)__CLEAN__Month::5", "i(Month):Wind__CLEAN__Month::5", "Temp:i(bonjour(test), drop = 5:12):Wind__CLEAN__bonjour::5", "i(a, pi(b))__CLEAN__a::o")
    #
    # Speed consideration:
    # for 2000 names: 10ms simple case
    #                 28ms complex case
    # price to pay is OK
    #

    res = x

    who2clean = grepl("__CLEAN__", x, fixed = TRUE)

    x2clean = x[who2clean]

    x_split = strsplit(x2clean, "(^|(?<=[^[:alnum:]\\._]))(i|sunab(_att)?)\\(|__CLEAN__", perl = TRUE)

    x_left = sapply(x_split, function(v) v[1])
    x_right = sapply(x_split, function(v) v[2])
    x_alias = sapply(x_split, function(v) v[3])

    x_right = gsub("^[^\\(\\)]+(\\(|\\))", "\\1", x_right)

    qui_ok = substr(x_right, 1, 1) == ")"

    x_right[qui_ok] = substr(x_right[qui_ok], 2, 500)

    if(any(!qui_ok)){

        clean_x_right_hard = function(letter_vec){
            open = 1 + cumsum(letter_vec == "(")
            close = cumsum(letter_vec == ")")
            letter_vec = letter_vec[-(1:which.max(close - open == 0))]
            paste(letter_vec, collapse = "")
        }

        x_right[!qui_ok] = sapply(strsplit(x_right[!qui_ok], ""), clean_x_right_hard)
    }

    res[who2clean] = paste0(x_left, x_alias, x_right)

    return(res)
}

is_naked_fun = function(x, fun_pattern){
    # Why is it always so complicated... There must be an easier way
    # x = c("i(x1)", "i(I(x3))", "i(x3, x4, TRUE, drop = c(1, 3:5))", "Temp:i(x2)", "i(x3):Wind")

    x_split = strsplit(x, paste0("(^|(?<=[^[:alnum:]\\._]))", fun_pattern, "\\("), perl = TRUE)

    left = sapply(x_split, function(x) x[1])
    right = sapply(x_split, function(x) x[2])

    right_naked = function(r){

        if(!grepl("(", r, fixed = TRUE) && grepl("\\)$", r)){
            return(TRUE)
        }

        letter_vec = strsplit(r, "")[[1]]
        open = 1 + cumsum(letter_vec == "(")
        close = cumsum(letter_vec == ")")
        which.max(close - open == 0) == length(letter_vec)
    }

    left_ok  = nchar(left) == 0
    right_ok = rep(FALSE, length(right))
    if(any(left_ok)){
        right_ok[left_ok] = sapply(right[left_ok], right_naked)
    }

    left_ok & right_ok
}

set_defaults = function(opts_name){

    opts = getOption(opts_name)
    if(is.null(opts) || length(opts) == 0){
        return(NULL)
    }

    sysOrigin = sys.parent()
    mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
    args_in = names(mc)

    args2set = setdiff(names(opts), args_in)

    for(v in args2set){
        assign(v, opts[[v]], parent.frame())
    }

}

fetch_data = function(x, prefix = "", suffix = ""){
    # x: fixest estimation
    # We try different strategies:
    # 1) using the environment where the estimation was done
    # 2) the "parent.frame()" defined as the frame on top of ALL fixest functions
    # 3) the global environment, if it wasn't in 1)

    # Maybe I should keep only 1) => is there a reason to add the others?

    # 1) safest
    # 2) less safe but OK => note ???
    # 3) kind of dangerous => warning() ???

    if(is.null(x$call$data)) return(NULL)

    # 1) Environment of the call

    data = NULL
    try(data <- eval(x$call$data, x$call_env), silent = TRUE)

    if(!is.null(data)){
        return(data)
    }

    # 2) First non fixest frame

    fixest_funs = ls(getNamespace("fixest"))

    i = 2
    sysOrigin = sys.parent(i)
    while(sysOrigin != 0 && as.character(sys.call(sysOrigin)[[1]]) %in% fixest_funs){
        i = i + 1
        sysOrigin = sys.parent(i)
    }

    if(sysOrigin != 0){
        # We try again...
        try(data <- eval(x$call$data, parent.frame(sysOrigin)), silent = TRUE)

        if(!is.null(data)){
            return(data)
        }
    }

    # 3) Global environment

    if(!identical(parent.env(x$call_env), .GlobalEnv)){
        # ...and again
        try(data <- eval(x$call$data, .GlobalEnv), silent = TRUE)

        if(!is.null(data)){
            return(data)
        }
    }

    # => Error message

    if(nchar(prefix) == 0){
        msg = "W"
    } else {
        s = ifelse(grepl(" $", prefix), "", " ")
        if(grepl("\\. *$", prefix)){
            msg = paste0(prefix, s, "W")
        } else {
            msg = paste0(prefix, s, "w")
        }
    }

    if(nchar(prefix) == 0){
        msg = "W"
    } else if(grepl("\\. *$", prefix)){
        msg = paste0(gsub(" +$", "", prefix), " W")
    } else {
        msg = paste0(gsub(prefix, " +$", ""), " w")
    }

    if(nchar(suffix) > 0){
       suffix = gsub("^ +", "", suffix)
    }

    stop_up(msg, "e fetch the data in the enviroment where the estimation was made, but the data does not seem to be there any more (btw it was ", charShorten(deparse(x$call$data)[1], 15), "). ", suffix)


}

is_large_object = function(x){

    if(is.list(x)){
        if(length(x[[1]]) > 10000){
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else if(length(x)[1] > 10000){
        return(TRUE)
    }

    FALSE
}

build_flags = function(mc, ..., call_env){
    # Returns a list of arguments
    # If the arguments are too large => we use calls instead and save the calling environment
    # Note that the arguments that will be passed in ... must NOT be missing // they must be initialized to NULL
    # All that stuff just to avoid too large objects....

    dots = list(...)

    args = names(dots)

    res = list()
    for(i in seq_along(dots)){

        x = dots[[i]]
        if(is.null(x)){
            # nothing // NULL arguments are still in 'dots', so we want to avoid them
        } else if(is_large_object(x)){
            res[[args[i]]] = mc[[args[i]]]
            if(is.null(res$call_env)){
                if(missing(call_env)){
                    res$call_env = new.env(parent = parent.frame(2))
                } else {
                    res$call_env = call_env
                }
            }
        } else {
            res[[args[i]]] = x
        }
    }

    res
}

assign_flags = function(flags, ...){
    # flags == arguments
    # We get the arguments used when the function was originally called
    # we assign these arguments to the previous frame

    dots = list(...)

    res = list()

    for(arg in setdiff(names(flags), "call_env")){

        if(is.null(dots[[arg]])){
            # we only assign values that were not provided by the user
            x = flags[[arg]]
            if(is.name(x) || is.call(x)){
                x = eval(x, flags$call_env)
            }
            assign(arg, x, parent.frame())
        }
    }
}

items_to_drop = function(items, x, varname, keep = FALSE, argname,
                         keep_first = FALSE, up = 1, no_error = FALSE,
                         valid_ref = FALSE){
    # selection of items
    # the selection depends on the type of x
    # always returns the IDs of the items to drop

    set_up(up)

    if(missing(argname)){
        argname = deparse(substitute(x))
    }

    ref = argname == "ref" && !valid_ref

    if(keep_first && (keep || ref)){
        stop("Internal error: keep_first should not be used with 'keep' or 'ref'.")
    }

    if(is.factor(items)){
        items = as.character(items)
    }

    if(is.factor(x)){
        x = as.character(x)
    }


    if(inherits(x, "formula")){

        if(!identical(all.vars(x), "x")){
            stop_up("In argument '", argname, "', the formula must contain a single variable name: 'x'. So far '", deparse_long(x), "' is not valid.")
        }

        if(length(x) > 2){
            stop_up("In argument '", argname, "', if a formula, it must be one-sided. Problem: '", deparse_long(x), "' is two-sided.")
        }

        is_here = error_sender(eval(x[[2]], list(x = items)),
                               "In argument '", argname, "', the evaluation of the formula led to an error:")
        if(length(is_here) != length(items)){
            stop_up("In argument '", argname, "', the evaluation of the formula must return a logical vector of the same length as 'x'. Problem: '", deparse_long(x), "' returns a vector of length ", length(is_here), " (expected: ", length(items), ").")
        }

        if(!is.logical(is_here)){
            stop_up("In argument '", argname, "', the evaluation of the formula must return a logical vector. Problem: '", deparse_long(x), "' is not logical (instead it is of class ", enumerate_items(class(is_here)), ").")
        }

        is_here = !is.na(is_here) & is_here

        if(!no_error && !any(is_here)){
            stop_up("In argument '", argname, "', the evaluation of the formula must match at least one value. Problem: '", deparse_long(x), "' does not match any.")
        }

        if(keep){
            id_drop = which(!is_here)
        } else {
            id_drop = which(is_here)
        }

    } else if(is.character(x)){
        all_x = c()
        for(i in seq_along(x)){

            my_x = x[i]
            if(grepl("^@", my_x)){
                # A) regex
                pattern = substr(my_x, 2, nchar(my_x))
                new_x = grep(pattern, items, value = TRUE)
                if(length(new_x) == 0 && !no_error){
                    # strong checking!
                    stop_up("In argument '", argname, "', the regular expression '", pattern, "' does not match any value of '", varname, "'.")
                }
                all_x = c(all_x, new_x)
            } else {
                # B) partial matching
                if(no_error){
                    qui_ok = pmatch(tolower(my_x), tolower(items), duplicates.ok = TRUE)
                    my_x = items[qui_ok]
                    my_x = my_x[!is.na(my_x)]
                } else {
                    check_value_plus(my_x, "match", .choices = unique(items), .message = paste0("The argument '", argname, "' should contain values of the variable '", varname, "'."))
                }

                all_x = c(all_x, my_x)
            }

            if(keep_first && i == 1){
                all_x_first = all_x
            }

        }

        if(ref){
            if(length(all_x) == 1){
                id_drop = which(items %in% all_x)
            } else {
                id_drop = c(which(items %in% all_x[1]), which(items %in% all_x[-1]))
            }
        } else if(keep){
            id_drop = which(!items %in% all_x)
        } else {
            id_drop = which(items %in% all_x)

            if(keep_first){
                id_drop = unique(c(which(items %in% all_x_first), id_drop))
            }
        }

    } else {
        # exact matching
        if(keep){
            id_drop = which(!items %in% x)
        } else {

            if(ref){
                if(length(x) == 1){
                    id_drop = which(items %in% x)
                } else {
                    id_drop = c(which(items %in% x[1]), which(items %in% x[-1]))
                }
            } else {
                id_drop = which(items %in% x)

                if(keep_first){

                    if(!x[1] %in% items && !no_error){
                        stop_up("In argument '", argname, "', the value ", x[1], " does not match any value of '", varname, "'.")
                    }

                    id_drop = unique(c(which(items %in% x[1]), id_drop))
                }

            }

            if(length(id_drop) == 0 && !no_error){
                stop_up("In argument '", argname, "', the value", plural_len(x, "s.don't"), " match any value of '", varname, "'.")
            }
        }

    }

    id_drop
}



bin_factor = function(bin, x, varname, no_error = FALSE){
    # x = base_did$period
    # bin = list("9" = c(9, 2, 3), "7" = "@7|8", 5:6)
    # varname = "period" ; argname = "bin"

    set_up(1)

    argname = deparse(substitute(bin))
    check_arg(bin, "list | vector | os formula")

    #
    # DSB expansion
    #

    bin_dp = deparse_long(bin)
    if(grepl(".[", bin_dp, fixed = TRUE)){

        bin_names = names(bin)
        if(!is.null(bin_names)){
            qui = which(grepl(".[", bin_names, fixed = TRUE))
            for(i in qui){
                bin_names[i] = error_sender(dsb(bin_names[i], frame = parent.frame(2), nest = FALSE,
                                                vectorize = TRUE, collapse = ""),
                                            dsb("Error when binning: the name (.[bin_names[i]]) expanded",
                                                " with '.[]' led to an error:"))

            }
        }
        names(bin) = bin_names

        qui = which(sapply(bin, function(x) is.character(x) && any(grepl(".[", x, fixed = TRUE))))
        for(i in qui){
            value = as.list(bin[[i]])
            qui_bis = which(grepl(".[", value, fixed = TRUE))
            for(j in qui_bis){
                value[[j]] = error_sender(.dsb(value[[j]], frame = parent.frame(2), nest = FALSE),
                                           dsb("Error when binning: the name (.[value[[j]]]) expanded",
                                               " with '.[]' led to an error:"))
            }
            bin[[i]] = unlist(value)
        }

    }

    #
    # cut:: => special treatment, we short circuit
    #

    if(!is.list(bin) && (is.character(bin) && any(grepl("^cut::", bin)))){

        if(!grepl("^cut::", bin[1])){
            stop_up("To use the special binning 'cut::values', it must be in the first position of the vector, possibly followed by custom names. Currently 'cut::values' is in position ", which(grepl("^cut::", bin))[1], ".")
        }

        n_names = length(bin) - 1
        if(n_names > 1){
            my_cut = trimws(sub("^cut::", "", bin[1]))
            if(grepl("^[[:digit:]]+$", my_cut)){
                n_cuts = as.numeric(my_cut)
            } else {
                n_cuts = length(strsplit(my_cut, "\\[|\\]")[[1]]) + 1
            }

            if(n_names > n_cuts){
                stop_up("In the special binning 'cut::values', the number of custom names is greater than the number of bins (", n_names, " vs ", n_cuts, ").")
            }
        }

        if(!is.numeric(x)){
            stop_up("To use the special binning 'cut::values', the variable '", varname, "' must be numeric. Currently this is not the case (it is of class ", enumerate_items(class(x)), " instead).")
        }

        return(cut_vector(x, bin))
    }

    x_int = to_integer(x, add_items = TRUE, items.list = TRUE, sorted = TRUE)
    x_items = x_int$items

    do_factor = FALSE
    if(is.factor(x_items)){
        do_factor = TRUE
        x_items = as.character(x_items)
    }

    if(!is.list(bin)){
        if(is.character(bin) && any(grepl("^!?!?bin::", bin))){

            if(length(bin) > 1){
                stop_up("To use the special binning 'bin::digit', the argument '", argname, "' must be of length 1. Currently it is of length ", length(bin), ".")
            }

            d = gsub("^!?!?bin::", "", bin)
            if(any(grepl("[^[:digit:]]", d))){
                bin_type = gsub("^(!?!?bin).*", "\\1", bin)
                stop_up("In the argument bin, the special binning must be of the form '", bin_type, "::digit'. Currently this is not the case for '", bin_type, "::", d, "'.")
            }
            d = as.numeric(d)

            consecutive = grepl("^!", bin)
            from_last = grepl("^!!", bin)

            if(!consecutive){
                if(!is.numeric(x_items)){
                    stop_up("To use the special binning 'bin::digit', the variable '", varname, "' must be numeric. Currently this is not the case (it is of class ", enumerate_items(class(x_items)), " instead).")
                }

                new_x = (x_items %/% d) * d
                new_x_unik = unique(new_x)

            } else {
                n = length(x_items)
                x_seq = if(from_last) n:1 else 1:n
                new_x = ((x_seq - 1) %/% d) * d + 1
                new_x_unik = unique(new_x)
            }

            bin = list()
            for(i in seq_along(new_x_unik)){
                bin[[i]] = x_items[new_x == new_x_unik[i]]
            }

        } else {
            bin = list(bin)
        }
    }

    # we catch "@" as first item
    if(identical(bin[[1]], "@") && length(bin) > 1){
        # this means that the user wants to specify
        # the binned values as first elements of the new factor
        bin[[1]] = NULL
        if(is.null(names(bin))){
            names(bin) = paste0("@", seq_along(bin), " ", names(bin))
        } else {
            qui_ok = !grepl("^@\\d+", names(bin))
            names(bin)[qui_ok] = paste0("@", seq(sum(qui_ok)), " ", names(bin)[qui_ok])
        }
    }

    x_map = x_range = seq_along(x_int$items)
    id_bin = list()
    for(i in seq_along(bin)){
        id_bin[[i]] = items_to_drop(x_items, bin[[i]], varname, up = 2, argname = argname,
                                    keep_first = TRUE, no_error = no_error, valid_ref = TRUE)
    }

    # sanity check
    if(any(table(unlist(id_bin)) > 1)){
        t_id = table(unlist(id_bin))
        n_max = max(t_id)
        pblm = x_items[as.numeric(names(t_id)[t_id == n_max])]
        stop_up("In 'bin', some values are binned in different bins, it's of course not allowed. The value '", pblm, "' is in ", n_max, " bins.")
    }

    # recreating the factor

    # case no_error = TRUE
    if(any(lengths(id_bin) == 0)){
        qui_ok = lengths(id_bin) > 0
        id_bin = id_bin[qui_ok]
        if(length(id_bin) == 0) return(x)
        bin = bin[qui_ok]
    }

    id_not_core = FALSE
    for(i in seq_along(bin)){
        id_bin_i = id_bin[[i]]
        id_core = id_bin_i[1]

        id_change = x_range %in% setdiff(id_bin_i, id_core)
        id_not_core = id_not_core | id_change

        x_map = x_map - cumsum(id_change)
    }

    # final pass
    for(i in seq_along(bin)){
        id_bin_i = id_bin[[i]]
        id_core = id_bin_i[1]
        id_change = x_range %in% id_bin_i

        x_map[id_change] = x_map[id_core]
    }

    # changing the item values
    x_items_new = x_items

    bin_names = names(bin)
    custom_loc = FALSE # custom location
    if(is.null(bin_names)){
        bin_names = character(length(bin))

    } else if(any(grepl("^@\\d", bin_names))){
        custom_loc = TRUE
        do_factor = TRUE
        x_items_new = as.character(x_items_new)

        qui = grepl("^@\\d", bin_names)
        bin_location = sub("^@(\\d+).*", "\\1", bin_names)
        bin_location[!qui] = NA
        bin_location = as.numeric(bin_location)

        bin_names = sub("^@\\d+ *", "", bin_names)
    }

    for(i in seq_along(bin)){

        b_name = bin_names[[i]]
        if(nchar(b_name) == 0){
            b_name = as.character(x_items_new[id_bin[[i]][1]])
            bin_names[[i]] = b_name # useful for new_loc
        }

        if(is.numeric(x_items_new)){
            b_name_num = tryCatch(as.numeric(b_name), warning = function(x) "not numeric")
            if(identical(b_name_num, "not numeric")){
                # we convert to character
                do_factor = TRUE
                x_items_new = as.character(x_items_new)
                b_name_num = b_name
            }

            x_items_new[id_bin[[i]]] = b_name_num
        } else {
            x_items_new[id_bin[[i]]] = b_name
        }
    }

    x_items_new = x_items_new[!id_not_core]

    if(custom_loc){
        n_bins = length(bin_names)
        new_loc = rep(NA_real_, length(x_items_new))
        for(i in 1:n_bins){
            b_name = bin_names[i]
            new_loc[x_items_new == b_name] = bin_location[i]
        }

        loc_no_na = new_loc[!is.na(new_loc)]
        while(anyDuplicated(loc_no_na)){
            # regularization
            is_dup = duplicated(loc_no_na)
            value = min(loc_no_na[is_dup])
            i_first = which.max(loc_no_na == value)
            loc_no_na[loc_no_na >= value] = loc_no_na[loc_no_na >= value] + 1
            loc_no_na[i_first] = value
        }

        n_max = length(x_items_new)
        if(any(loc_no_na > n_max)){
            # regularization

            n_lnona = length(loc_no_na)
            if(n_lnona == 1){
                loc_no_na = n_max
            } else {
                order_loc_no_na = order(loc_no_na)

                loc_sorted = loc_no_na[order_loc_no_na]
                loc_sorted[n_lnona] = n_max
                for(i in n_lnona:2) {
                    if(loc_sorted[i] <= loc_sorted[i - 1]){
                        loc_sorted[i - 1] = loc_sorted[i] - 1
                    }
                }

                loc_no_na = loc_sorted[order(order_loc_no_na)]
            }
        }

        new_loc[!is.na(new_loc)] = loc_no_na

        for(i in seq_along(new_loc)){
            if(!i %in% new_loc){
                j = which(is.na(new_loc))[1]
                new_loc[j] = i
            }
        }

        x_map = new_loc[x_map]
        x_items_new = x_items_new[order(new_loc)]
    }

    if(do_factor){
        x_items_new = factor(x_items_new, levels = x_items_new)
    }

    x_new = x_items_new[x_map[x_int$x]]

    return(x_new)
}


cut_vector = function(x, bin){
    # We cut a vector into pieces
    # **only numeric vectors**
    # cut::a]b[c] cuts the vectors into four slices: [-Inf, a]; ]a, b[; [b, c]; ]c, Inf]
    # a, b, c should be numeric values
    # they can be replaced with: pXX or qX, percentiles and quartiles
    # cut::n splits the data into n equal sizes

    set_up(2)

    bin_names_usr = bin[-1]
    bin = bin[1]

    # checking bin
    if(!is.character(bin) || !grepl("^cut::", bin)){
        stop("Internal bug: Argument 'bin' should be equal to 'cut::stg' over here -- this is not the case.")
    }

    # cleaning x
    n = length(x)
    ANY_NA = anyNA(x)
    IS_NA = FALSE
    if(ANY_NA){
        IS_NA = is.na(x)

        if(all(IS_NA)){
            return(x)
        }

        x = x[!IS_NA]
    }

    # we sort x (needed later)
    x_order = order(x)
    x_sorted = x[x_order]

    my_cut = gsub("^cut::", "", bin)

    if(is_numeric_in_char(my_cut)){
        my_cut = as.numeric(my_cut)
        pct = cumsum(rep(1/my_cut, max(my_cut - 1, 1)))
        cut_points = quantile(x_sorted, pct)

        bounds = rep("]", length(cut_points))

    } else {
        bounds = regmatches(my_cut, gregexpr("\\[|\\]", my_cut))[[1]]
        values = strsplit(my_cut, "\\[|\\]")[[1]]

        n_cuts = length(values)

        if(n_cuts != length(bounds)){
            stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles/percentiles. Problem: each number must be followed by an open (or closed) square bracket, this is currentlly not the case.")
        }

        cut_points = numeric(n_cuts)
        for(i in seq_along(values)){

            v = trimws(values[i])

            if(is_numeric_in_char(v)){
                cut_points[i] = as.numeric(v)

            } else if(grepl("^p|P|q|Q", v)){
                p = gsub("^.", "", v)
                if(!is_numeric_in_char(p)){
                    stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect.")
                }

                p = as.numeric(p)


                if(grepl("^q|Q", v)){
                    # we transform the quartile into a percentile
                    if(!p %in% c(0:4)){
                        stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect. \n  The quartile must be an integer between 0 and 4.")
                    }

                    p = c(0, 25, 50, 75, 100)[p + 1]
                }

                if(!p %in% 0:100){
                    stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect. \n  The percentile must be an integer between 0 and 100.")
                }

                cut_points[i] = quantile(x, p/100)

            } else {
                stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value '", v, "', in '", bin, "', is incorrect. This is not a percentile nor a number.")
            }
        }

        if(n_cuts > 1 && any(cut_points[1:(n_cuts-1)] > cut_points[1 + 1:(n_cuts-1)])){
            i_pblm = which(cut_points[1:(n_cuts-1)] > cut_points[1 + 1:(n_cuts-1)])[1]
            stop_up("In 'bin', the format should be 'cut::a]b]' with 'a', 'b', etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The values 'a', 'b', etc should be increasing, but the ", n_th(i_pblm), " value (", cut_points[i_pblm], ") is larger than the ", n_th(i_pblm + 1), " (", cut_points[i_pblm + 1], ").")
        }

    }

    bin_names = character(length(cut_points) + 1)
    for(i in seq_along(bin_names_usr)){
        bi = bin_names_usr[i]
        if(!is.na(bi)){
            bin_names[i] = bi
        }
    }

    is_included = 1 * (bounds == "]")
    x_cut = cpp_cut(x_sorted, cut_points, is_included)

    x_int = x_cut$x_int
    isnt_empty = x_cut$isnt_empty
    value_min = x_cut$value_min
    value_max = x_cut$value_max
    is_int = x_cut$is_int

    n_bins = sum(isnt_empty)

    if(any(isnt_empty == 0)){
        i_empty = which(isnt_empty == 0)
        if(getFixest_notes()){
            message("When binning: in '", bin, "', the ", enumerate_items(n_th(i_empty)), " bin", plural_len(i_empty, "s.is"), " empty.")
        }

        x_int = cumsum(isnt_empty)[x_int]
        value_min = value_min[-i_empty]
        value_max = value_max[-i_empty]
        bin_names = bin_names[-i_empty]
    }

    # creating the labels
    labels = character(n_bins)
    if(is_int){
        # format A-B

        # we don't want to add a comma to the years! But to large numbers: yes
        fmt_fun = if(value_max[n_bins] > 9999) dreamerr::fsignif else function(x) as.character(x)

        for(i in 1:n_bins){
            if(value_min[i] == value_max[i]){
                labels[i] = fmt_fun(value_min[i])
            } else {
                labels[i] = paste0(fmt_fun(value_min[i]), "-", fmt_fun(value_max[i]))
            }
        }

    } else {
        # format [A, B]
        # we always write in inclusion, too hard to do the exclusion mentally
        d_min = NULL
        for(i in 1:n_bins){
            vmin = value_min[i]
            vmax = value_max[i]

            l10_diff = if(vmax == vmin) 0 else log10(vmax - vmin)
            # d: nber of digits
            d = if(l10_diff >= 0) 1 else ceiling(abs(l10_diff))
            if(!is.null(d_min) && d < d_min) d = d_min

            # we check consistency problems
            d_min = NULL
            vmax_fmt = sprintf("%.*f", d, vmax)
            if(i != n_bins){
                vmin_next_fmt = sprintf("%.*f", d, value_min[i + 1])
                if(vmax_fmt == vmin_next_fmt){
                    l10_diff = log10(value_min[i + 1] - vmax)
                    d = ceiling(abs(l10_diff))
                    d_min = d
                }
                vmax_fmt = cpp_add_commas(vmax, d, FALSE)
            }
            vmin_fmt = cpp_add_commas(vmin, d, FALSE)

            labels[i] = paste0("[", vmin_fmt, "; ", vmax_fmt, "]")
        }
    }

    for(i in seq_along(bin_names)){
        bi = bin_names[i]
        if(nchar(bi) > 0){
            labels[i] = bi
        }
    }


    # Te result
    if(ANY_NA){
        res = rep(NA_real_, n)
        res[!IS_NA] = x_int[order(x_order)]
    } else {
        res = x_int[order(x_order)]
    }

    res = factor(res, labels = labels)

    return(res)
}


fixest_pvalue = function(x, zvalue, vcov){
    # compute the pvalue for a fixest estimation

    if(use_t_distr(x)){

        if(missing(vcov)) {
            stop("Internal error (=bug): the argument 'vcov' should not be missing in fixest_pvalue().")

        } else if(is.null(attr(vcov, "dof.K"))){
            stop("Internal error (=bug): the attribute 'dof.K' from 'vcov' should not be NULL.")
        }

        # df.t is always an attribute of the vcov
        df.t = attr(vcov, "df.t")
        if(is.null(df.t)){
            df.t = max(nobs(x) - attr(vcov, "dof.K"), 1)
        }

        pvalue = 2*pt(-abs(zvalue), df.t)

    } else {
        pvalue = 2*pnorm(-abs(zvalue))
    }

    pvalue
}

fixest_CI_factor = function(x, level, vcov){

    val = (1 - level) / 2
    val = c(val, 1 - val)

    if(use_t_distr(x)){

        if(missing(vcov)) {
            stop("Internal error (=bug): the argument 'vcov' should not be missing in fixest_CI_factor().")

        } else if(is.null(attr(vcov, "dof.K"))){
            stop("Internal error (=bug): the attribute 'dof.K' from 'vcov' should not be NULL.")
        }

        # df.t is always an attribute of the vcov
        df.t = attr(vcov, "df.t")
        if(is.null(df.t)){
            df.t = max(nobs(x) - attr(vcov,"dof.K"), 1)
        }

        fact = qt(val, df.t)

    } else {
        fact = qnorm(val)
    }

    fact
}


#### ................. ####
#### Small Utilities ####
####

escape_regex = function(x){
    # escape special characters in regular expressions => to make it as "fixed"

    res = gsub("((?<=[^\\\\])|(?<=^))(\\$|\\.|\\+|\\(|\\)|\\[|\\]|\\?|\\^)", "\\\\\\2", x, perl = TRUE)
    res
}

.addCommas = function(x){

	if(!is.finite(x)) return(as.character(x))

    cpp_add_commas(x)
}

addCommas = function(x){
	sapply(x, .addCommas)
}

decimalFormat = function(x){

    who_valid = which(!is.na(x) & is.numeric(x))
    if(length(who_valid) == 0) return(x)

    res = x

    x_valid = x[who_valid]
    xPower = log10(abs(x_valid))

    if(min(xPower) > 0){
        pow_round = max(1, 6 - ceiling(xPower))
    } else if(min(xPower) < -5){
        pow_round = ceiling(abs(min(xPower))) + 2
    } else {
        pow_round = 6
    }

    res[who_valid] = round(x_valid, pow_round)

    res
}

numberFormat_single = function(x, type = "normal"){
	# For numbers higher than 1e9 => we apply a specific formatting
	# idem for numbers lower than 1e-4

    if(is.character(x)) return(x)

    if(is.na(x)){
        return(NA)
    }

	if(x == 0) return("0")

	exponent = floor(log10(abs(x)))

	if(-4 < exponent && exponent < 9){
		if(exponent > 0){
			return(addCommas(x))
		} else {
			return(as.character(decimalFormat(x)))
		}

	}

	left_value = round(x*10**-exponent, 3)

	if(type == "latex"){
		res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
	} else {
		res = paste0(left_value, "e", ifelse(exponent > 0, "+", ""), exponent)
	}

	res
}

numberFormatLatex = function(x){
	sapply(x, numberFormat_single, type = "latex")
}

numberFormatNormal = function(x){
	sapply(x, numberFormat_single)
}

mysignif = function (x, d = 2, r = 1){

    .mysignif = function(x, d, r) {
        if (is.na(x)) {
            return(NA)
        }

        if(abs(x) >= 10^(d - 1)){
            return(round(x, r))
        } else {
            return(signif(x, d))
        }
    }
    sapply(x, .mysignif, d = d, r = r)
}

format_nber_single = function(x, digits, round = FALSE, pow_above = 10, pow_below = -5, tex = FALSE){
    # Attention aux nombres ronds => pas chiffre apres la virgule!

    if(is.na(x) || !is.numeric(x)){
        return(x)
    }

    if(round){
        # super simple
        return(cpp_add_commas(x, digits))

    } else {
        whole = (x %% 1) == 0
        if(abs(x) >= 10^(digits - 1)){
            x = round(x, 1)
        } else {
            x = signif(x, digits)
        }
    }

    exponent = floor(log10(abs(x)))

    if(exponent >= pow_above || exponent <= pow_below){
        left_value = round(x*10**-exponent, min(digits, 2))

        if(tex){
            res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
        } else {
            res = paste0(left_value, "e", ifelse(exponent > 0, "+", ""), exponent)
        }

    } else if(exponent >= 0){
        r = max(-(exponent + 1 - digits), 1)
        res = cpp_add_commas(x, r, whole)

    } else if(abs(x) > 10**(-digits)){
        res = sprintf("%.*f", digits, x)

    } else {
        res = sprintf("%.*f", abs(exponent), x)
    }

    res
}

format_number = function(x, digits = 4, round = FALSE, pow_above = 10, pow_below = -5, tex = FALSE){
    sapply(x, format_nber_single, digits = digits, round = round, pow_above = pow_above, pow_below = pow_below, tex = tex)
}

index_2D_to_1D = function(i, j, n_j) 1 + (i - 1) * n_j + j - 1

par_fit = function(my_par, id){
    # simple function that extends the plot parameters
    my_par = rep(my_par, ceiling(max(id) / length(my_par)))
    my_par[id]
}

dict_apply = function(x, dict = NULL){

    check_arg(dict, "NULL named character vector no na", .message = "The argument 'dict' must be a dictionnary, ie a named vector (eg dict=c(\"old_name\"=\"new_name\")")

    if(missing(dict) || length(dict) == 0){
        return(x)
    }

    # We make the dictionary names space agnostic, adds a handful of us only
    if(any(grepl(" ", x, fixed = TRUE))){
        x_tmp = gsub(" ", "", x, fixed = TRUE)
    } else {
        x_tmp = x
    }

    if(any(grepl(" ", names(dict), fixed = TRUE))){
        names(dict) = gsub(" ", "", names(dict), fixed = TRUE)
        if(anyDuplicated(names(dict))){
            dup = duplicated(names(dict))
            stop("The dictionary 'dict' contains several items with the same names, it concerns ",
                 enumerate_items(names(dict)[dup]), " (note that spaces are ignored).")
        }
    }

    who = x_tmp %in% names(dict)
    x[who] = dict[as.character(x_tmp[who])]
    x
}

keep_apply = function(x, keep = NULL, logical = FALSE){

    if(missing(keep) || length(keep) == 0){
        if(logical){
            return(rep(TRUE, length(x)))
        } else {
            return(x)
        }
    }

    check_arg(keep, "character vector no na", .message = "The arg. 'keep' must be a vector of regular expressions (see help(regex)).")

    res = x

    qui_keep = rep(FALSE, length(x))
    for(var2keep in keep){
        if(grepl("^!%", var2keep)) var2keep = gsub("^!%", "%!", var2keep)

        vect2check = res
        if(grepl("^%", var2keep)){
            var2keep = gsub("^%", "", var2keep)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In keep, the special character '%' cannot be used here.")
                vect2check = res
            }
        }

        if(grepl("^!", var2keep)){
            qui_keep = qui_keep | !grepl(substr(var2keep, 2, nchar(var2keep)), vect2check)
        } else {
            qui_keep = qui_keep | grepl(var2keep, vect2check)
        }
    }

    if(logical) return(qui_keep)

    res[qui_keep]
}

drop_apply = function(x, drop = NULL){

    if(missing(drop) || length(drop) == 0){
        return(x)
    }

    check_arg(drop, "character vector no na", .message = "The arg. 'drop' must be a vector of regular expressions (see help(regex)). ")

    res = x

    for(var2drop in drop){
        if(grepl("^!%", var2drop)) var2drop = gsub("^!%", "%!", var2drop)

        vect2check = res
        if(grepl("^%", var2drop)){
            var2drop = gsub("^%", "", var2drop)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In drop, the special character '%' cannot be used here.")
                vect2check = res
            }
        }

        if(grepl("^!", var2drop)){
            res = res[grepl(substr(var2drop, 2, nchar(var2drop)), vect2check)]
        } else {
            res = res[!grepl(var2drop, vect2check)]
        }
    }

    res
}

order_apply = function(x, order = NULL){

    if(missing(order) || length(order) == 0){
        return(x)
    }

    check_arg(order, "character vector no na", .message = "The arg. 'order' must be a vector of regular expressions (see help(regex)). ")

    res = x

    for(var2order in rev(order)){
        if(grepl("^!%", var2order)) var2order = gsub("^!%", "%!", var2order)

        vect2check = res
        if(grepl("^%", var2order)){
            var2order = gsub("^%", "", var2order)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In order, the special character '%' cannot be used here.")
                vect2check = res
            }
        }

        if(grepl("^!", var2order)){
            who = !grepl(substr(var2order, 2, nchar(var2order)), vect2check)
            res = c(res[who], res[!who])
        } else {
            who = grepl(var2order, vect2check)
            res = c(res[who], res[!who])
        }
    }

    res
}

charShorten = function(x, width){
	# transforms characters => they can't go beyond a certain width
	# two dots are added to suggest longer character
	# charShorten("bonjour", 5) => "bon.."
	n = nchar(x)

	if(n > width && width > 3){
		res = substr(x, 1, width - 2)
		res = paste0(res, "..")
	} else {
		res = x
	}

	res
}


show_vars_limited_width = function(charVect, nbChars = 60, addS = FALSE){
	# There are different cases

	n = length(charVect)

	if(n == 1){
		return(charVect)
	}

	text_s = ifelse(addS, "s ", "")

	nb_char_vect = nchar(charVect)
	sumChar = cumsum(nb_char_vect) + (0:(n-1))*2 + 3 + 1

	if(max(sumChar) < nbChars){
		text = paste0(text_s, paste0(charVect[-n], collapse = ", "), " and ", charVect[n])
		return(text)
	}

	qui = max(which.max(sumChar > nbChars - 8) - 1, 1)

	nb_left = n - qui

	text = paste0(text_s, paste0(charVect[1:qui], collapse = ", "), " and ", nb_left, " other", ifelse(nb_left>1, "s", ""), ".")

	return(text)
}


char2num = function(x, addItem = FALSE){
	# we transform the data to numeric => faster analysis

	# special case
	qui = which(x == "")
	if(length(qui) > 0){
		x[qui] = "xxEMPTYxx"
	}

	x_unik = unique(x)
	dict = 1:length(x_unik)
	names(dict) = x_unik
	x_num = dict[x]

	names(x_num) = NULL

	if(addItem){
		res = list(x = x_num, items = x_unik)
		return(res)
	} else {
		return(x_num)
	}

}

quickUnclassFactor = function(x, addItem = FALSE, sorted = FALSE){
	# does as unclass(as.factor(x))
	# but waaaaay quicker

    not_num = !is.numeric(x)
    is_char_convert = not_num && !is.character(x)

    if(is_char_convert){
        res = cpp_quf_gnl(as.character(x))
    } else {
        res = cpp_quf_gnl(x)
    }

    if(sorted){

        if(not_num){
            items = x[res$x_unik]
        } else {
            items = res$x_unik
        }

        x = res$x_uf

        new_order = order(items)
        order_new_order = order(new_order)
        x_uf = order_new_order[x]

        if(addItem){
            if(is_char_convert){
                res = list(x = x_uf, items = as.character(items[new_order]))
            } else {
                res = list(x = x_uf, items = items[new_order])
            }

            return(res)
        } else {
            return(x_uf)
        }
    }

    if(addItem){

        if(not_num){
            if(is_char_convert){
                items = as.character(x[res$x_unik])
            } else {
                items = x[res$x_unik]
            }

            res = list(x = res$x_uf, items = items)
        } else {
            names(res) = c("x", "items")
        }

        return(res)
    } else {
        return(res$x_uf)
    }
}

missnull = function(x) missing(x) || is.null(x)

MISSNULL = function(x){
    # same a missnull but we also check evaluation

    if(missing(x)) return(TRUE)

    # we check evaluation and nullity
    error_sender(x, up = 1, arg_name = deparse(substitute(x)))

    is.null(x)
}

isScalar = function(x, int = FALSE) {
    if(length(x) == 1L && is.numeric(x) && is.finite(x)){
        if(int){
            return(x %% 1 == 0)
        } else {
            return(TRUE)
        }
    } else {
        return(FALSE)
    }
}

isLogical = function(x) length(x) == 1L && is.logical(x) && !is.na(x)

isSingleChar = function(x) length(x) == 1L && is.character(x) && !is.na(x)

fml2varnames = function(fml, combine_fun = FALSE){
	# This function transforms a one sided formula into a
	# character vector for each variable

	# In theory, I could just use terms.formula to extract the variable.
	# But I can't!!!! Because of this damn ^ argument.
	# I need to apply a trick

    # combine_fun: whether to add a call to combine_clusters_fast

    # Only the ^ users "pay the price"

    if("^" %in% all.vars(fml, functions = TRUE)){
        # new algo using fml_combine, more robust
        fml_char = as.character(fml)[2]
        all_var_names = fml_combine(fml_char, TRUE, vars = TRUE)
        if(!combine_fun){
            all_var_names = rename_hat(all_var_names)
        }

    } else {
        t = terms(fml)
        all_var_names = attr(t, "term.labels")
        # => maybe I should be more cautious below???? Sometimes ':' really means stuff like 1:5
        all_var_names = gsub(":", "*", all_var_names) # for very special cases
    }


	all_var_names
}

msg_na_inf = function(any_na, any_inf){
    if(any_na && any_inf){
        res = "NA and infinite values"
    } else if(any_na){
        res = "NA values"
    } else {
        res = "infinite values"
    }

    res
}


which_na_inf = function(x, nthreads){
    # returns the same elements as cppar_which_na_inf_mat

    is_num = sapply(x, is.numeric)

    if(all(is_num)){
        res = cpppar_which_na_inf_mat(as.matrix(x), nthreads)
    } else {

        # numeric
        if(any(is_num)){
            res = cpppar_which_na_inf_mat(as.matrix(x[, is_num, drop = FALSE]), nthreads)
        } else {
            res = list(any_na_inf = FALSE, any_na = FALSE, any_inf = FALSE, is_na_inf = FALSE)
        }

        # non numeric
        is_na = rowSums(is.na(x[, !is_num, drop = FALSE])) > 0

        res$any_na = res$any_na || any(is_na)
        res$is_na_inf = res$is_na_inf || res$any_na
        res$is_na_inf = res$is_na_inf | is_na
    }

    return(res)
}

extract_fe_slope = function(t){
    # input: x$fixef_terms which contains the slopes and the fixef

    fixef_vars = t[!grepl("\\[", t)]
    slope_terms = t[grepl("\\[", t)]
    slope_vars = gsub(".+\\[|\\]", "", slope_terms)
    slope_fe = gsub("\\[.+", "", slope_terms)
    fe_all = gsub("\\[.+", "", t)

    list(fixef_vars=fixef_vars, slope_vars=slope_vars, slope_fe=slope_fe, slope_terms=slope_terms, fe_all=fe_all)
}

deparse_long = function(x){
    dep_x = deparse(x, width.cutoff = 500)
    if(length(dep_x) == 1){
        return(dep_x)
    } else {
        return(paste(gsub("^ +", "", dep_x), collapse = ""))
    }
}

isVector = function(x){
    # it seems that when you subselect in data.table
    # sometimes it does not yield a vector
    # so i cannot use is.vector to check the consistency

    if(is.vector(x)){
        return(TRUE)
    } else {
        if(is.null(dim(x)) && !is.list(x)){
            return(TRUE)
        }
    }
    return(FALSE)
}


fill_with_na = function(x, object){
    if(is.null(object$obs_selection)){
        return(x)
    }

    res = rep(NA, object$nobs_origin)
    qui = 1:object$nobs_origin

    for(i in seq_along(object$obs_selection)){
        qui = qui[object$obs_selection[[i]]]
    }

    res[qui] = x

    return(res)
}

trim_obs_removed = function(x, object){
    if(is.null(object$obs_selection)){
        return(x)
    }

    for(i in seq_along(object$obs_selection)){
        x = x[object$obs_selection[[i]]]
    }

    x
}

is_operator = function(x, op) if(length(x) <= 1) FALSE else x[[1]] == op

fml_breaker = function(fml, op){
    res = list()
    k = 1
    while(is_operator(fml, op)){
        res[[k]] = fml[[3]]
        k = k + 1
        fml = fml[[2]]
    }
    res[[k]] = fml

    res
}

fml_maker = function(lhs, rhs){

    while(is_operator(lhs, "(")){
        lhs = lhs[[2]]
    }

    if(missing(rhs)){
        if(is_operator(lhs, "~")){
            return(lhs)
        }
        res = ~ .
        res[[2]] = lhs
    } else {

        while(is_operator(rhs, "(")){
            rhs = rhs[[2]]
        }

        res = . ~ .
        res[[2]] = lhs
        res[[3]] = rhs
    }

    res
}

fml_split_internal = function(fml, split.lhs = FALSE){

    fml_split_tilde = fml_breaker(fml, "~")
    k = length(fml_split_tilde)

    # NOTA in fml_breaker: to avoid copies, the order of elements returned is reversed

    # currently res is the LHS
    res = list(fml_split_tilde[[k]])

    if(k == 2){
        rhs = fml_breaker(fml_split_tilde[[1]], "|")
        l = length(rhs)

    } else if(k == 3){
        rhs  = fml_breaker(fml_split_tilde[[2]], "|")
        l = length(rhs)
        rhs_right = fml_breaker(fml_split_tilde[[1]], "|")

        if(length(rhs_right) > 1){
            stop_up("Problem in the formula: the formula in the RHS (expressing the IVs) cannot be multipart.")
        }

        # The rightmost element of the RHS is in position 1!!!!
        iv_fml = fml_maker(rhs[[1]], rhs_right[[1]])

        rhs[[1]] = iv_fml

    } else {
        # This is an error
        stop_up("Problem in the formula: you cannot have more than one RHS part containing a formula.")
    }

    if(!split.lhs){
        new_fml = fml_maker(res[[1]], rhs[[l]])

        res[[1]] = new_fml
        if(l == 1) return(res)
        res[2:l] = rhs[(l - 1):1]

    } else {
        res[1 + (1:l)] = rhs[l:1]
    }

    res
}


fml_split = function(fml, i, split.lhs = FALSE, text = FALSE, raw = FALSE){
    # I had to create that fonction to cope with the following formula:
    #
    #                       y ~ x1 | fe1 | u ~ z
    #
    # For Formula to work well, one would need to write insstead: y | x1 | fe1 | (u ~ z)
    # that set of parentheses are superfluous

    my_split = fml_split_internal(fml, split.lhs)

    if(raw){
        return(my_split)
    } else if(text){

        if(!missing(i)){
            return(deparse_long(my_split[[i]]))
        } else {
            return(sapply(my_split, deparse_long))
        }

    } else if(!missing(i)) {
        return(fml_maker(my_split[[i]]))

    } else {
        res = lapply(my_split, fml_maker)

        return(res)
    }

}

error_sender = function(expr, ..., clean, up = 0, arg_name){
    res = tryCatch(expr, error = function(e) structure(list(conditionCall(e), conditionMessage(e)), class = "try-error"))

    if("try-error" %in% class(res)){
        set_up(1 + up)
        msg = paste0(..., collapse = "")

        if(nchar(msg) == 0){
            if(missing(arg_name)){
                arg_name = deparse(substitute(expr))
            }
            msg = paste0("Argument '", arg_name, "' could not be evaluated: ")
            stop_up(msg, res[[2]])

        } else {

            # The expression that is evaluated is absolutely non informative for the user
            call_non_informative = deparse(substitute(expr), 100)[1]
            call_error = deparse(res[[1]], 100)[1]

            if(call_error == call_non_informative || call_error == "NULL" ||
               grepl("^(doTry|eval)", call_error)){
                call_error = ""

            } else {
                call_error = paste0("In ", call_error, ": ")
            }

            err = res[[2]]

            if(grepl("^in eval\\(str[^:]+:\n", err)){
                err = sub("^in eval\\(str[^:]+:\n", "", err)
            }

            if(!missing(clean)){

                if(grepl(" => ", clean)){
                    clean_split = strsplit(clean, " => ")[[1]]
                    from = clean_split[1]
                    to = clean_split[2]
                } else {
                    from = clean
                    to = ""
                }

                stop_up(msg, "\n  ", call_error, gsub(from, to, err))
            } else {
                stop_up(msg, "\n  ", call_error, err)
            }
        }
    }

    res
}

is_fml_inside = function(fml){
    # we remove parentheses first

    while(is_operator(fml, "(")){
        fml = fml[[2]]
    }

    is_operator(fml, "~")
}


merge_fml = function(fml_linear, fml_fixef = NULL, fml_iv = NULL){

    is_fe = length(fml_fixef) > 0
    is_iv = length(fml_iv) > 0

    if(!is_fe && !is_iv){
        res = fml_linear
    } else {
        fml_all = deparse_long(fml_linear)

        if(is_fe){
            # we add parentheses if necessary
            if(is_operator(fml_fixef[[2]], "|")){
                fml_all[[2]] = paste0("(", as.character(fml_fixef)[2], ")")
            } else {
                fml_all[[2]] = as.character(fml_fixef)[2]
            }
        }

        if(is_iv) fml_all[[length(fml_all) + 1]] = deparse_long(fml_iv)

       res = as.formula(paste(fml_all, collapse = "|"))
    }

    res
}


fixest_fml_rewriter = function(fml){
    # Currently performs the following
    # - expands lags
    # - protects powers: x^3 => I(x^3)
    #
    # fml = sw(f(y, 1:2)) ~ x1 + l(x2, 1:2) + x2^2 | fe1 | y ~ z::e + g^3
    # fml = y ~ 1 | id + period | l(x_endo, -1:1) ~ l(x_exo, -1:1)

    fml_text = deparse_long(fml)

    isPanel = grepl("(^|[^\\._[:alnum:]])(f|d|l)\\(", fml_text)
    isPower = grepl("^", fml_text, fixed = TRUE)

    if(isPanel){
        # We rewrite term-wise

        fml_parts = fml_split(fml, raw = TRUE)
        n_parts = length(fml_parts)

        #
        # LHS
        #

        # only panel: no power (bc no need), no interact

        # We tolerate multiple LHS and expansion
        lhs_text = fml_split(fml, 1, text = TRUE, split.lhs = TRUE)

        if(grepl("^(c|c?sw0?|list)\\(", lhs_text)){
            lhs_text2eval = gsub("^(c|c?sw0?|list)\\(", "sw(", lhs_text)
            lhs_names = eval(str2lang(lhs_text2eval))
        } else {
            lhs_names = lhs_text
        }

        lhs_all = error_sender(expand_lags_internal(lhs_names),
                               "Problem in the formula regarding lag/leads: ", clean = "__expand")

        if(length(lhs_all) > 1){
            lhs_fml = paste("c(", paste(lhs_all, collapse = ","), ")")
        } else {
            lhs_fml = lhs_all
        }

        lhs_text = lhs_fml

        #
        # RHS
        #

        # power + panel + interact

        if(isPower){
            # rhs actually also contains the LHS
            rhs_text = deparse_long(fml_parts[[1]])
            rhs_text = gsub("(?<!I\\()(\\b(\\.[[:alpha:]]|[[:alpha:]])[[:alnum:]\\._]*\\^[[:digit:]]+)", "I(\\1)", rhs_text, perl = TRUE)

            if(grepl("\\^[[:alpha:]]", rhs_text)){
                stop_up("The operator '^' between variables can be used only in the fixed-effects part of the formula. Otherwise, please use ':' instead.")
            }

            fml_rhs = as.formula(rhs_text)
        } else {
            fml_rhs = fml_maker(fml_parts[[1]])
        }

        rhs_terms = get_vars(fml_rhs)

        if(length(rhs_terms) == 0){
            rhs_text = "1"
        } else {
            rhs_terms = error_sender(expand_lags_internal(rhs_terms),
                                     "Problem in the formula regarding lag/leads: ", clean = "__expand")

            if(attr(terms(fml_rhs), "intercept") == 0){
                rhs_terms = c("-1", rhs_terms)
            }

            rhs_text = paste(rhs_terms, collapse = "+")
        }

        fml_linear = as.formula(paste0(lhs_text, "~", rhs_text))

        #
        # FE + IV
        #

        fml_fixef = fml_iv = NULL

        if(n_parts > 1){

            #
            # FE
            #

            # Only isPanel (although odd....)

            is_fe = !is_fml_inside(fml_parts[[2]])
            if(is_fe){

                if(identical(fml_parts[[2]], 0)){
                    fml_fixef = NULL

                } else {

                    fml_fixef = fml_maker(fml_parts[[2]])
                    fml_fixef_text = deparse_long(fml_fixef)

                    if(grepl("(l|d|f)\\(", fml_fixef_text)){
                        # We need to make changes
                        # 1st: for terms to work => we change ^ if present (sigh)

                        do_sub = grepl("^", fml_fixef_text, fixed = TRUE)

                        if(do_sub){
                            fml_fixef = as.formula(gsub("^", "%^%", fml_fixef_text, fixed = TRUE))
                        }

                        fixef_terms = attr(terms(fml_fixef), "term.labels")
                        fixef_text = error_sender(expand_lags_internal(fixef_terms),
                                                  "Problem in the formula regarding lag/leads: ", clean = "__expand")

                        if(do_sub){
                            fixef_text = gsub(" %^% ", "^", fixef_text, fixed = TRUE)
                        }

                        fml_fixef = as.formula(paste("~", paste(fixef_text, collapse = "+")))

                    }
                }
            }

            #
            # IV
            #

            if(n_parts == 3 || !is_fe){
                fml_iv = fml_maker(fml_parts[[n_parts]])

                fml_endo = .xpd(lhs = ~y, rhs = fml_iv[[2]])
                fml_inst = .xpd(lhs = ~y, rhs = fml_iv[[3]])

                endo_lag_expand = fixest_fml_rewriter(fml_endo)$fml
                inst_lag_expand = fixest_fml_rewriter(fml_inst)$fml

                fml_iv = .xpd(lhs = endo_lag_expand[[3]], rhs = inst_lag_expand[[3]])
            }
        }

        fml_new = merge_fml(fml_linear, fml_fixef, fml_iv)

    } else if(isPower){
        # This is damn costly.... but I have to deal with corner cases....
        # I think I should do that at the C level, this will be much faster I guess

        # We only take care of the RHS (we don't care about the LHS)
        no_lhs_text = gsub("^[^~]+~", "", fml_text)
        no_lhs_text = gsub("(?<!I\\()(\\b(\\.[[:alpha:]]|[[:alpha:]])[[:alnum:]\\._]*\\^[[:digit:]]+)", "I(\\1)", no_lhs_text, perl = TRUE)

        if(grepl("\\^[[:alpha:]]", no_lhs_text)){
            # We check if there is one ^ specifically in the RHS or in the IV part
            fml_parts = fml_split(fml, raw = TRUE)
            n_parts = length(fml_parts)

            rhs_txt = fml_split(fml, i = 1, text = TRUE)

            if(grepl("\\^[[:alpha:]]", rhs_txt)){
                stop_up("The operator '^' between variables can be used only in the fixed-effects part of the formula. Otherwise, please use ':' instead.")
            }

            if(is_fml_inside(fml_parts[[n_parts]])){
                iv_txt = fml_split(fml, i = n_parts, text = TRUE)

                if(grepl("\\^[[:alpha:]]", iv_txt)){
                    stop_up("The operator '^' between variables can be used only in the fixed-effects part of the formula. Otherwise, please use ':' instead.")
                }
            }

        }

        fml_new = as.formula(paste0(gsub("~.+", "", fml_text), "~", no_lhs_text))

    } else {
        res = list(fml = fml, isPanel = FALSE)
        return(res)
    }

    res = list(fml = fml_new, isPanel = isPanel)

    return(res)
}


check_set_types = function(x, types, msg){
    arg_name = deparse(substitute(x))
    check_arg(x, "os formula | character vector no na", .arg_name = arg_name, .up = 1)

    if("formula" %in% class(x)){
        x = attr(terms(x), "term.labels")
    }

    check_value_plus(x, "multi match", .choices = types, .arg_name = arg_name, .up = 1)

    x
}

check_set_digits = function(digits, up = 1){
    # Note that the argument name can be either digits or digits.stats

    set_up(up)
    check_value(digits, "integer scalar GE{1} | character scalar", .arg_name = deparse(substitute(digits)))

    if(is.character(digits)){
        d_type = substr(digits, 1, 1)
        d_value = substr(digits, 2, nchar(digits))

        # Control
        if(!d_type %in% c("r", "s")){
            arg = deparse(substitute(digits))
            stop_up("The argument '", arg, "' must start with 'r' (for round) or 's' (for significant). Currently it starts with '", d_type,"' which is not valid.\nExample of valid use: digits = 'r3'.")
        }

        round = d_type == "r"

        if(!grepl("^[0-9]$", d_value)){
            arg = deparse(substitute(digits))
            stop_up("The argument '", arg, "' must be equal to the character 'r' or 's' followed by a single digit. Currently '", digits,"' is not valid.\nExample of valid use: digits = '", d_type, "3'.")
        }

        digits = as.numeric(d_value)

    } else {
        round = FALSE
    }

    list(digits = digits, round = round)
}

get_vars = function(x){
    attr(terms(x), "term.labels")
}

mat_posdef_fix = function(X, tol = 1e-10){
    # X must be a symmetric matrix
    # We don't check it

    if(any(diag(X) < tol)){
        e = eigen(X)
        dm = dimnames(X)
        X = tcrossprod(e$vectors %*% diag(pmax(e$values, tol), nrow(X)), e$vectors)
        dimnames(X) = dm
    }

    return(X)
}


is_fixest_call = function(){
    nf = sys.nframe()

    if(nf < 5) return(FALSE)

    last_calls = sapply(tail(sys.calls(), 13), function(x) deparse(x)[1])

    any(grepl("fixest", last_calls[-length(last_calls)]))
}

all_vars_with_i_prefix = function(fml){
    # fml = a ~ x1^x2 + i(x3, i.x4) + x5*i(x6, I(i.x7))

    vars = all.vars(fml)
    if(any(grepl("^i\\..+", vars))){
        fml_dp = deparse_long(fml)
        # for i. to work it MUST be the second argument (var)
        # valid cases:
        # - i(x1, i.x2)
        # - i(var = i.x2, x1)
        # - i(x1, i.I(x7))
        # Not valid:
        # - i(x1, I(i.x7))
        #
        # This means that in the parsed formula, i. is always preceded by a space and
        # either a "," or a "="

        # Maybe later: add nice error messages reminding how to use i.

        qui_i = which(grepl("^i\\..+", vars))
        i_vars = vars[qui_i]
        for(i in seq_along(i_vars)){
            fml_split = strsplit(fml_dp, i_vars[i], fixed = TRUE)[[1]]
            n = length(fml_split) - 1
            for(j in 1:n){
                part = fml_split[j]
                if(grepl("(,|var =) *$", part)){
                    part = gsub("\\([^\\)]+\\)", "", part)
                    if(grepl("i\\(", part)){
                        # OK!
                        ii = qui_i[i]
                        vars[ii] = substr(vars[ii], 3, nchar(vars[ii]))
                    }
                }
            }
        }
    }

    vars
}


# function to normalize character vectors into variable names
as_varname = function(x){
    # "x1" => "x1"
    # "(x1)" => "`(x1)`"


    qui_pblm = grepl("[^[:alnum:]\\._]", x)
    if(any(qui_pblm)){
        x[qui_pblm] = paste0("`", x[qui_pblm], "`")
        x[qui_pblm] = gsub("``", "`", x[qui_pblm])
    }

    x
}

update_file = function(path, text){

    if(file.exists(path)){
        text_1 = readLines(path)

        text_clean_1 = unlist(strsplit(text_1, "\n"))
        text_clean_2 = unlist(strsplit(text, "\n"))

        text_clean_1 = text_clean_1[grepl("[[:alnum:][:punct:]]", text_clean_1)]
        text_clean_2 = text_clean_2[grepl("[[:alnum:][:punct:]]", text_clean_2)]

        do_write = length(text_clean_1) != length(text_clean_2) || any(text_clean_1 != text_clean_2)

    } else {
        do_write = TRUE
    }

    # We write only if the text is different
    if(do_write){
        message("writing '", path, "'", appendLF = FALSE)
        f = file(path, "w", encoding = "utf-8")
        writeLines(text, f)
        close(f)
        message(".")
    }
}


fetch_arg_deparse = function(arg){
    # Utility to find out what the user originally provided as argument
    # Only useful for functions that are deeply nested (typically vcov)

    sc = rev(sys.calls())

    # initialization of the name
    arg_name = deparse_long(sc[[2]][[arg]])

    n = length(sc)
    if(n > 2){
        for(i in 2:length(sc)){
            mc = sc[[i]]
            if(grepl("[tT]ry", mc[[1]])) next
            if(!arg %in% names(mc)) break
            arg_name = deparse_long(mc[[arg]])
        }
    }

    arg_name
}



# Simple utilities to avoid writing too much
# adds a small overhead but that's OK (about 1us which is OK)
any_missing = function(x1, x2, x3, x4, x5, x6){
    n = length(sys.call()) - 1
    switch(n,
           "1" = missing(x1),
           "2" = missing(x1) || missing(x2),
           "3" = missing(x1) || missing(x2) || missing(x3),
           "4" = missing(x1) || missing(x2) || missing(x3) || missing(x4),
           "5" = missing(x1) || missing(x2) || missing(x3) || missing(x4) || missing(x5),
           "6" = missing(x1) || missing(x2) || missing(x3) || missing(x4) || missing(x5) || missing(x6))
}

all_missing = function(x1, x2, x3, x4, x5, x6){
    n = length(sys.call()) - 1
    switch(n,
           "1" = missing(x1),
           "2" = missing(x1) && missing(x2),
           "3" = missing(x1) && missing(x2) && missing(x3),
           "4" = missing(x1) && missing(x2) && missing(x3) && missing(x4),
           "5" = missing(x1) && missing(x2) && missing(x3) && missing(x4) && missing(x5),
           "6" = missing(x1) && missing(x2) && missing(x3) && missing(x4) && missing(x5) && missing(x6))
}

use_t_distr = function(x){
    # whether to use the t-distribution or the normal
    # x: fixest estimation
    x$method %in% "feols" || (x$method %in% "feglm" && !x$family$family %in% c("poisson", "binomial"))
}


is_fun_in_char = function(fml_char, fun){
    extract_fun(fml_char, fun, bool = TRUE)
}

extract_fun = function(fml_char, fun, err_msg = NULL, bool = FALSE, drop = TRUE){
    # A formula, in character form // can be a vector
    # fun: a function name in regex form
    #
    # returns a list:
    # before
    # fun
    # after
    #
    # fml_char = "to_integer(id, sw0(fe))" ;  fun = "sw0?"

    only_one = TRUE

    regex = paste0(c("(?<=^)", "(?<=[^[:alnum:]\\._])"),
                   fun, "\\(", collapse = "|")

    is_there = grepl(regex, fml_char, perl = TRUE)

    if(bool) return(is_there)

    res = list()

    n = length(fml_char)
    for(i in 1:n){
        # fml_i = fml_char[1]

        fml_i = fml_char[i]
        if(is_there[i]){

            fml_split = strsplit(fml_i, regex, perl = TRUE)[[1]]

            if(only_one && length(fml_split) > 2){
                if(is.null(err_msg)){
                    stop_up("Only one function '", fun, "' can be used at a time.")
                } else {
                    stop_up(err_msg)
                }
            }

            # we need to add a last character otherwise => problems
            fml_value = paste0(fml_split[2], " ")
            fml_right_split = strsplit(fml_value, ")", fixed = TRUE)[[1]]
            n_open = pmax(lengths(strsplit(fml_right_split, "(", fixed = TRUE)) - 1, 0)

            n_parts = length(fml_right_split)
            fun_closed = which((1 + cumsum(n_open) - 1:n_parts) == 0)

            fun_name = substr(fml_i, nchar(fml_split[1]) + 1, nchar(fml_i))
            fun_name = gsub("\\(.+", "", fun_name)

            fun_char = paste0(fun_name, "(", paste(fml_right_split[1:fun_closed], collapse = ")"), ")")

            fml_right_rest = trimws(paste0(fml_right_split[-(1:fun_closed)], collapse = ")"))

            res[[i]] = list(before = fml_split[1],
                            fun = fun_char,
                            after = fml_right_rest)


        } else {
            res[[i]] = list(before = fml_i,
                            fun = "",
                            after = "")
        }
    }

    if(n == 1 && drop) res = res[[1]]

    return(res)
}

is_numeric_in_char = function(x){
    res = tryCatch(as.numeric(x), warning = function(x) "not numeric")
    !identical(res, "not numeric")
}

insert_in_between = function(x, y){
    n_x = length(x)
    n_y = length(y)

    if(n_y == 1) y = rep(y, n_x)
    if(n_x == 1) x = rep(x, n_y)
    n_x = length(x)

    res = rep(x, each = 2)
    res[2 * 1:n_x] = y

    res
}

str_trim = function(x, n_first = 0, n_last = 0){
    # positive values: first
    # negative values: last

    if(is.character(n_first)){
        n_first = nchar(n_first)
    } else if(n_first < 0){
        n_last = -n_first
        n_first = 0
    }

    res = substr(x, 1 + n_first, nchar(x) - n_last)
}

str_split = function(x, split){
    # Simple wrapper

    fixed = TRUE
    perl = FALSE
    if(grepl("@", split, fixed = TRUE)){
        if(grepl("^\\\\@", split)){
            split = str_trim(split, 1)
        } else if(grepl("^@", split)){
            split = str_trim(split, 1)
            fixed = FALSE
            perl = TRUE
        }
    }

    strsplit(x, split, fixed = fixed, perl = perl)
}

NA_fun = function(..., df){
    dots = list(...)
    n_dots = length(dots)

    if(n_dots == 0){
        res = !complete.cases(df)
        return(res)
    }

    res = is.na(dots[[1]])

    if(n_dots == 1){
        return(res)
    }

    for(i in 2:n_dots){
        res = res | is.na(dots[[i]])
    }
    res
}


eval_dot = function(x, up = 1){
    # Note that the right use of up is ESSENTIAL
    # if must refer to the parent frame from which the main call has been made
    # Otherwise an error will be thrown of "." not existing

    x_dp = deparse(substitute(x), 300)

    sysOrigin = sys.parent(up)
    mc = match.call(sys.function(sysOrigin), sys.call(sysOrigin))

    if(!x_dp %in% names(mc)){
        return(x)
    }

    my_list = list("." = list)

    eval(mc[[x_dp]], my_list, parent.frame(up + 1))
}


is_user_level_call = function(){
    length(sys.calls()) <= 2
}


is_calling_fun = function(pattern){
    sc_all = sys.calls()
    n_sc = length(sc_all)
    if(n_sc > 2){

        if(grepl(".fixest", sc_all[[n_sc - 1]][[1]], fixed = TRUE)){
            if(n_sc == 3){
                return(FALSE)
            }

            sc = sc_all[[n_sc - 3]]
        } else {
            sc = sc_all[[n_sc - 2]]
        }

        return(grepl(pattern, as.character(sc[[1]])[1]))
    }

    FALSE
}


#### ................. ####
#### Additional Methods ####
####

# Here we add common statistical functions

#' Extracts the number of observations form a \code{fixest} object
#'
#' This function simply extracts the number of observations form a \code{fixest} object, obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams summary.fixest
#'
#' @param ... Not currently used.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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
#' This function computes the AIC (Aikake's, an information criterion) from a \code{fixest} estimation.
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Optionally, more fitted objects.
#' @param k A numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC (i.e. \code{AIC=-2*LL+k*nparams}).
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2\times LogLikelihood + k\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this criterion on \code{\link[stats]{AIC}}.
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Other statictics methods: \code{\link[fixest]{BIC.fixest}}, \code{\link[fixest]{logLik.fixest}}, \code{\link[fixest]{nobs.fixest}}.
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
#' This function computes the BIC (Bayesian information criterion) from a \code{fixest} estimation.
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
#' You can have more information on this criterion on \code{\link[stats]{AIC}}.
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Other statistics functions: \code{\link[fixest]{AIC.fixest}}, \code{\link[fixest]{logLik.fixest}}.
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
#' This function extracts the log-likelihood from a \code{fixest} estimation.
#'
#' @inheritParams nobs.fixest
#'
#' @param ... Not currently used.
#'
#' @details
#' This function extracts the log-likelihood based on the model fit. You can have more information on the likelihoods in the details of the function \code{\link[fixest]{femlm}}.
#'
#' @return
#' It returns a numeric scalar.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Other statistics functions: \code{\link[fixest]{AIC.fixest}}, \code{\link[fixest]{BIC.fixest}}.
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

#' Extracts the coefficients from a \code{fixest} estimation
#'
#' This function extracts the coefficients obtained from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams etable
#'
#' @param agg Logical scalar, default is \code{TRUE}. If the coefficients of the estimation have been aggregated, whether to report the aggregated coefficients. If \code{FALSE}, the raw coefficients will be returned.
#' @param ... Not currently used.
#'
#' @details
#' The coefficients are the ones that have been found to maximize the log-likelihood of the specified model. More information can be found on the models from the estimations help pages: \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' Note that if the model has been estimated with fixed-effects, to obtain the fixed-effect coefficients, you need to use the function \code{\link[fixest]{fixef.fixest}}.
#'
#' @return
#' This function returns a named numeric vector.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{confint.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{etable}}, \code{\link[fixest]{fixef.fixest}}.
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
coef.fixest = coefficients.fixest = function(object, keep, drop, order, agg = TRUE, ...){

    check_arg(keep, drop, order, "NULL character vector no na")
    check_arg(agg, "logical scalar")

    if(isTRUE(object$is_agg) && agg){
        res = object$coeftable[, 1]
        names(res) = rownames(object$coeftable)
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


#' Extracts fitted values from a \code{fixest} fit
#'
#' This function extracts the fitted values from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. The fitted values that are returned are the \emph{expected predictor}.
#'
#' @inheritParams nobs.fixest
#'
#' @param type Character either equal to \code{"response"} (default) or \code{"link"}. If \code{type="response"}, then the output is at the level of the response variable, i.e. it is the expected predictor \eqn{E(Y|X)}. If \code{"link"}, then the output is at the level of the explanatory variables, i.e. the linear predictor \eqn{X\cdot \beta}.
#' @param na.rm Logical, default is \code{TRUE}. If \code{FALSE} the number of observation returned will be the number of observations in the original data set, otherwise it will be the number of observations used in the estimation.
#' @param ... Not currently used.
#'
#' @details
#' This function returns the \emph{expected predictor} of a \code{fixest} fit. The likelihood functions are detailed in \code{\link[fixest]{femlm}} help page.
#'
#' @return
#' It returns a numeric vector of length the number of observations used to estimate the model.
#'
#' If \code{type = "response"}, the value returned is the expected predictor, i.e. the expected value of the dependent variable for the fitted model: \eqn{E(Y|X)}.
#' If \code{type = "link"}, the value returned is the linear predictor of the fitted model, that is \eqn{X\cdot \beta} (remind that \eqn{E(Y|X) = f(X\cdot \beta)}).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{resid.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
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

#' Extracts residuals from a \code{fixest} object
#'
#' This function extracts residuals from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#'
#' @param type A character scalar, either \code{"response"} (default), \code{"deviance"}, \code{"pearson"}, or \code{"working"}. Note that the \code{"working"} corresponds to the residuals from the weighted least square and only applies to \code{\link[fixest]{feglm}} models.
#' @param na.rm Logical, default is \code{TRUE}. Whether to remove the observations with NAs from the original data set. If \code{FALSE}, then the vector returned is always of the same length as the original data set.
#' @param ... Not currently used.
#'
#'
#' @return
#' It returns a numeric vector of the length the number of observations used for the estimation (if \code{na.rm = TRUE}) or of the length of the original data set (if \code{na.rm = FALSE}).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{fitted.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
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

    check_arg_plus(type, "match")
    check_arg_plus(na.rm, "logical scalar")

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

#' Predict method for \code{fixest} fits
#'
#' This function obtains prediction from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams fitted.fixest
#' @inheritParams summary.fixest
#'
#' @param newdata A data.frame containing the variables used to make the prediction. If not provided, the fitted expected (or linear if \code{type = "link"}) predictors are returned.
#' @param sample Either "estimation" (default) or "original". This argument is only used when arg. 'newdata' is missing, and is ignored otherwise. If equal to "estimation", the vector returned matches the sample used for the estimation. If equal to "original", it matches the original data set (the observations not used for the estimation being filled with NAs).
#' @param se.fit Logical, default is \code{FALSE}. If \code{TRUE}, the standard-error of the predicted value is computed and returned in a column named \code{se.fit}. This feature is only available for OLS models not containing fixed-effects.
#' @param interval Either "none" (default), "confidence" or "prediction". What type of confidence interval to compute. Note that this feature is only available for OLS models not containing fixed-effects (GLM/ML models are not covered).
#' @param level A numeric scalar in between 0.5 and 1, defaults to 0.95. Only used when the argument 'interval' is requested, it corresponds to the width of the confidence interval.
#' @param fixef Logical scalar, default is \code{FALSE}. If \code{TRUE}, a data.frame is returned, with each column representing the fixed-effects coefficients for each observation in \code{newdata} -- with as many columns as fixed-effects. Note that when there are variables with varying slopes, the slope coefficients are returned (i.e. they are not multiplied by the variable).
#' @param vs.coef Logical scalar, default is \code{FALSE}. Only used when \code{fixef = TRUE} and when variables with varying slopes are present. If \code{TRUE}, the coefficients of the variables with varying slopes are returned instead of the coefficient multiplied by the value of the variables (default).
#' @param ... Not currently used.
#'
#'
#' @return
#' It returns a numeric vector of length equal to the number of observations in argument \code{newdata}.
#' If \code{newdata} is missing, it returns a vector of the same length as the estimation sample, except if \code{sample = "original"}, in which case the length of the vector will match the one of the original data set (which can, but also cannot, be the estimation sample).
#' If \code{fixef = TRUE}, a \code{data.frame} is returned.
#' If \code{se.fit = TRUE} or \code{interval != "none"}, the object returned is a data.frame with the following columns: \code{fit}, \code{se.fit}, and, if CIs are requested, \code{ci_low} and \code{ci_high}.
#'
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{update.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
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
	check_arg_plus(type, sample, "match")
	check_arg(fixef, vs.coef, "logical scalar")
	check_arg(se.fit, "logical scalar")
	check_arg(level, "numeric scalar GT{.50} LT{1}")
	check_arg_plus(interval, "match(none, confidence, prediction)")
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
		matrix_linear = error_sender(fixest_model_matrix_extra(object = object, newdata = newdata, original_data = FALSE, fml = rhs_fml, i_noref = TRUE), "Error when creating the linear matrix: ")

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


#' Confidence interval for parameters estimated with \code{fixest}
#'
#' This function computes the confidence interval of parameter estimates obtained from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams vcov.fixest
#'
#' @param parm The parameters for which to compute the confidence interval (either an integer vector OR a character vector with the parameter name). If missing, all parameters are used.
#' @param level The confidence level. Default is 0.95.
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
confint.fixest = function(object, parm, level = 0.95, vcov, se, cluster, ssc = NULL, ...){

    # Checking the arguments
    if(is_user_level_call()){
        validate_dots(suggest_args = c("parm", "level", "se", "cluster"),
                      valid_args = c("forceCovariance", "keepBounded"))
    }

	# Control
	if(!is.numeric(level) || !length(level) == 1 || level >= 1 || level <= .50){
		stop("The argument 'level' must be a numeric scalar greater than 0.50 and strictly lower than 1.")
	}

	# The proper SE
	sum_object = summary(object, vcov = vcov, se = se, cluster = cluster, ssc = ssc, ...)

	se_all = sum_object$coeftable[, 2]
	coef_all = sum_object$coeftable[, 1]

	# the parameters for which we should compute the confint
	all_params = names(coef_all)

	if(missing(parm)){
	    parm_use = all_params
	} else if(is.numeric(parm)){
	    if(any(parm %% 1 != 0)){
	        stop("If the argument 'parm' is numeric, it must be integers.")
	    }

	    parm_use = unique(na.omit(all_params[parm]))
	    if(length(parm_use) == 0){
	        stop("There are ", length(all_params), " coefficients, the argument 'parm' does not correspond to any of them.")
	    }
	} else if(is.character(parm)){
	    parm_pblm = setdiff(parm, all_params)
	    if(length(parm_pblm) > 0){
	        stop("some parameters of 'parm' have no estimated coefficient: ", paste0(parm_pblm, collapse=", "), ".")
	    }

	    parm_use = intersect(parm, all_params)
	}

	# multiplicative factor
	fact = fixest_CI_factor(object, level, sum_object$cov.scaled)

	# The confints
	# Note that for glm models, there is no profiling
	lower_bound = coef_all[parm_use] + fact[1] * se_all[parm_use]
	upper_bound = coef_all[parm_use] + fact[2] * se_all[parm_use]

	res = data.frame(lower_bound, upper_bound, row.names = parm_use)

	val = (1 - level) / 2
	names(res) = paste0(round(100*c(val, 1-val), 1), " %")

	attr(res, "type") = attr(se_all, "type")

	res
}

#' Updates a \code{fixest} estimation
#'
#' Updates and re-estimates a \code{fixest} model (estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). This function updates the formulas and use previous starting values to estimate a new \code{fixest} model. The data is obtained from the original \code{call}.
#'
#' @method update fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param fml.update Changes to be made to the original argument \code{fml}. See more information on \code{\link[stats]{update.formula}}. You can add/withdraw both variables and fixed-effects. E.g. \code{. ~ . + x2 | . + z2} would add the variable \code{x2} and the cluster \code{z2} to the former estimation.
#' @param nframes (Advanced users.) Defaults to 1. Number of frames up the stack where to perform the evaluation of the updated call. By default, this is the parent frame.
#' @param evaluate Logical, default is \code{TRUE}. If \code{FALSE}, only the updated call is returned.
#' @param ... Other arguments to be passed to the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @return
#' It returns a \code{fixest} object (see details in \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}).
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{fixef.fixest}}.
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


#' Extract the formula of a \code{fixest} fit
#'
#' This function extracts the formula from a \code{fixest} estimation (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). If the estimation was done with fixed-effects, they are added in the formula after a pipe (\dQuote{|}). If the estimation was done with a non linear in parameters part, then this will be added in the formula in between \code{I()}.
#'
#'
#' @param x An object of class \code{fixest}. Typically the result of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation.
#' @param type A character scalar. Default is \code{type = "full"} which gives back a formula containing the linear part of the model along with the fixed-effects (if any) and the IV part (if any). If \code{type = "linear"} then only the linear formula is returned. If \code{type = "NL"} then only the non linear in parameters part is returned.
#' @param ... Not currently used.
#'
#' @return
#' It returns a formula.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{model.matrix.fixest}}, \code{\link[fixest]{update.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}.
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

	check_arg_plus(type, "match")

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


#' Design matrix of a \code{fixest} object
#'
#' This function creates the left-hand-side or the right-hand-side(s) of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation.
#'
#' @method model.matrix fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param data If missing (default) then the original data is obtained by evaluating the \code{call}. Otherwise, it should be a \code{data.frame}.
#' @param type Character vector or one sided formula, default is "rhs". Contains the type of matrix/data.frame to be returned. Possible values are: "lhs", "rhs", "fixef", "iv.rhs1" (1st stage RHS), "iv.rhs2" (2nd stage RHS), "iv.endo" (endogenous vars.), "iv.exo" (exogenous vars), "iv.inst" (instruments).
#' @param na.rm Default is \code{TRUE}. Should observations with NAs be removed from the matrix?
#' @param subset Logical or character vector. Default is \code{FALSE}. If \code{TRUE}, then the matrix created will be restricted only to the variables contained in the argument \code{data}, which can then contain a subset of the variables used in the estimation. If a character vector, then only the variables matching the elements of the vector via regular expressions will be created.
#' @param as.matrix Logical scalar, default is \code{FALSE}. Whether to coerce the result to a matrix.
#' @param as.df Logical scalar, default is \code{FALSE}. Whether to coerce the result to a data.frame.
#' @param collin.rm Logical scalar, default is \code{TRUE}. Whether to remove variables that were found to be collinear during the estimation. Beware: it does not perform a collinearity check.
#' @param ... Not currently used.
#'
#' @return
#' It returns either a vector, a matrix or a data.frame. It returns a vector for the dependent variable ("lhs"), a data.frame for the fixed-effects ("fixef") and a matrix for any other type.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{formula.fixest}}, \code{\link[fixest]{update.fixest}}, \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}.
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

    check_arg_plus(as.matrix, as.df, collin.rm, "logical scalar")

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
#' This function extracts the terms of a \code{fixest} estimation, excluding the fixed-effects part.
#'
#' @param x A \code{fixest} object. Obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param ... Not currently used.
#'
#' @return
#' An object of class \code{c("terms", "formula")} which contains the terms representation of a symbolic model.
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

#### ............... ####
#### Setters/Getters ####
####

#' Sets/gets whether to display notes in \code{fixest} estimation functions
#'
#' Sets/gets the default values of whether notes (informing for NA and observations removed) should be displayed in \code{fixest} estimation functions.
#'
#' @param x A logical. If \code{FALSE}, then notes are permanently removed.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Change default with
#' setFixest_notes(FALSE)
#' feols(Ozone ~ Solar.R, airquality)
#'
#' # Back to default which is TRUE
#' setFixest_notes(TRUE)
#' feols(Ozone ~ Solar.R, airquality)
#'
setFixest_notes = function(x){
    check_arg(x, "mbt logical scalar")

	options("fixest_notes" = x)
}

#' @rdname setFixest_notes
getFixest_notes = function(){

    x = getOption("fixest_notes")
    if(length(x) != 1 || !is.logical(x) || is.na(x)){
        stop("The value of getOption(\"fixest_notes\") is currently not legal. Please use function setFixest_notes to set it to an appropriate value. ")
    }

    x
}

#' Sets/gets the number of threads to use in \code{fixest} functions
#'
#' Sets/gets the default number of threads to used in \code{fixest} estimation functions. The default is the maximum number of threads minus two.
#'
#'
#'
#' @param nthreads The number of threads. Can be: a) an integer lower than, or equal to, the maximum number of threads; b) 0: meaning all available threads will be used; c) a number strictly between 0 and 1 which represents the fraction of all threads to use. If missing, the default is to use 50\% of all threads.
#' @param save Either a logical or equal to \code{"reset"}. Default is \code{FALSE}. If \code{TRUE} then the value is set permanently at the project level, this means that if you restart R, you will still obtain the previously saved defaults. This is done by writing in the \code{".Renviron"} file, located in the project's working directory, hence we must have write permission there for this to work, and only works with Rstudio. If equal to "reset", the default at the project level is erased. Since there is writing in a file involved, permission is asked to the user.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' # Gets the current number of threads
#' getFixest_nthreads()
#'
#' # To set multi-threading off:
#' setFixest_nthreads(1)
#'
#' # To set it back to default:
#' setFixest_nthreads()
#'
#'
setFixest_nthreads = function(nthreads, save = FALSE){
	# By default, we use only 50% of threads (never use all)

    max_CRAN = as.numeric(Sys.getenv("OMP_THREAD_LIMIT"))
    max_CRAN[is.na(max_CRAN)] = 1000

	max_threads = min(cpp_get_nb_threads(), 1000, max_CRAN) # we cap at 1k nthreads

	check_arg_plus(save, "logical scalar | match(reset)")

	do_reset = identical(save, "reset")

	if(missing(nthreads) || is.null(nthreads)){
	    # We first get the default from the environment variable
	    # If it is missing => 50% of all threads

	    # 0.5 => 50% of all available threads (usually equiv to the nber of procs)

	    nthreads_default = renvir_get("fixest_nthreads")

	    if(!do_reset && !is.null(nthreads_default)){
	        if(!isScalar(nthreads_default) || nthreads_default < 0){
	            warning("The variable setting the number of threads in the .Renviron file is corrupted. It's value has been reset.")
	            renvir_update("fixest_nthreads", NULL)
	            nthreads_default = 0.5
	        }

	    } else {
	        nthreads_default = 0.5
	    }

	    nthreads = check_set_nthreads(nthreads_default)

	}

	nthreads = check_set_nthreads(nthreads)

	if(do_reset){
	    renvir_update("fixest_nthreads", NULL)
	} else if(save){
	    renvir_update("fixest_nthreads", nthreads)
	}

	options("fixest_nthreads" = nthreads)

	invisible()
}

#' @rdname setFixest_nthreads
getFixest_nthreads = function(){

    x = getOption("fixest_nthreads")
    if(length(x) != 1 || !is.numeric(x) || is.na(x) || x %% 1 != 0 || x < 0){
        stop("The value of getOption(\"fixest_nthreads\") is currently not legal. Please use function setFixest_nthreads to set it to an appropriate value. ")
    }

    x
}

#' Sets/gets the dictionary relabeling the variables
#'
#' Sets/gets the default dictionary used in the function \code{\link[fixest]{etable}}, \code{\link[fixest]{did_means}} and \code{\link[fixest]{coefplot}}. The dictionaries are used to relabel variables (usually towards a fancier, more explicit formatting) when exporting them into a Latex table or displaying in graphs. By setting the dictionary with \code{setFixest_dict}, you can avoid providing the argument \code{dict}.
#'
#'
#' @param dict A named character vector. E.g. to change my variable named "a" and "b" to (resp.) "$log(a)$" and "$bonus^3$", then use \code{dict = c(a="$log(a)$", b3="$bonus^3$")}. This dictionary is used in Latex tables or in graphs by the function \code{\link[fixest]{coefplot}}. If you want to separate Latex rendering from rendering in graphs, use an ampersand first to make the variable specific to \code{coefplot}.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' data(trade)
#' est = feols(log(Euros) ~ log(dist_km)|Origin+Destination+Product, trade)
#' # we export the result & rename some variables
#' esttex(est, dict = c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#'
#' # If you export many tables, it can be more convenient to use setFixest_dict:
#' setFixest_dict(c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#' esttex(est) # variables are properly relabeled
#'
setFixest_dict = function(dict){

	if(missing(dict) || is.null(dict)){
		options("fixest_dict" = NULL)
		return(invisible())
	}

	#
	# Controls
	#

	if(!is.character(dict) || !isVector(dict)){
		stop("Argument 'dict' must be a character vector.")
	}

	if(anyNA(dict)){
		stop("Argument 'dict' must be a character vector without NAs.")
	}

	# Formatting the names
	dict_names = names(dict)
	if(is.null(dict_names)){
		stop("Argument 'dict', the dictionary, must be a named vector. Currently it has no names.")
	}

	dict_names = gsub(" +", "", dict_names)
	td = table(dict_names)
	if(any(td > 1)){
		qui = which(dict_names %in% names(td)[td > 1])
		name_dup = unique(names(dict)[qui])
		stop("Argument 'dict' contains duplicated names: ", enumerate_items(name_dup))
	}

	options("fixest_dict" = dict)
}

#' @rdname setFixest_dict
getFixest_dict = function(){

    x = getOption("fixest_dict")
    if(length(x) > 0){
        if(!is.character(x) || !isVector(x) || anyNA(x)){
            stop("The value of getOption(\"fixest_dict\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value. ")
        }
    }

    x
}


#' @rdname print.fixest
setFixest_print = function(type = "table", fitstat = NULL){

    check_arg_plus(type, "match(coef, table)")

    if(!missnull(fitstat)){
        fitstat = fitstat_validate(fitstat, TRUE)
    }

    # Getting the existing defaults
    opts = getOption("fixest_print")

    if(is.null(opts)){
        opts = list()
    } else if(!is.list(opts)){
        warning("Wrong formatting of option 'fixest_print', all options are reset.")
        opts = list()
    }

    # Saving the default values
    mc = match.call()
    args_default = names(mc)[-1]

    # NOTA: we don't allow delayed evaluation => all arguments must have hard values
    for(v in args_default){
        opts[[v]] = eval(as.name(v))
    }

    options(fixest_print = opts)
}


#' @rdname print.fixest
getFixest_print = function(){

    x = getOption("fixest_print")
    if(!(is.null(x) || is.list(x))){
        stop("The value of getOption(\"fixest_print\") is currently not legal. Please use function setFixest_print to set it to an appropriate value. ")
    }

    x
}



#' Sets/gets formula macros
#'
#' You can set formula macros globally with \code{setFixest_fml}. These macros can then be used in \code{fixest} estimations or when using the function \code{\link[fixest:setFixest_fml]{xpd}}.
#'
#' @inherit xpd examples
#'
#' @param ... Definition of the macro variables. Each argument name corresponds to the name of the macro variable. It is required that each macro variable name starts with two dots (e.g. \code{..ctrl}). The value of each argument must be a one-sided formula or a character vector, it is the definition of the macro variable. Example of a valid call: \code{setFixest_fml(..ctrl = ~ var1 + var2)}. In the function \code{xpd}, the default macro variables are taken from \code{getFixest_fml}, any variable in \code{...} will replace these values. You can enclose values in \code{.[]}, if so they will be evaluated from the current environment. For example \code{..ctrl = ~ x.[1:2] + .[z]} will lead to \code{~x1 + x2 + var} if \code{z} is equal to \code{"var"}.
#' @param reset A logical scalar, defaults to \code{FALSE}. If \code{TRUE}, all macro variables are first reset (i.e. deleted).
#'
#' @details
#' In \code{xpd}, the default macro variables are taken from \code{getFixest_fml}. Any value in the \code{...} argument of \code{xpd} will replace these default values.
#'
#' The definitions of the macro variables will replace in verbatim the macro variables. Therefore, you can include multipart formulas if you wish but then beware of the order the the macros variable in the formula. For example, using the airquality data, say you want to set as controls the variable \code{Temp} and \code{Day} fixed-effects, you can do \code{setFixest_fml(..ctrl = ~Temp | Day)}, but then \code{feols(Ozone ~ Wind + ..ctrl, airquality)} will be quite different from \code{feols(Ozone ~ ..ctrl + Wind, airquality)}, so beware!
#'
#' @return
#' The function \code{getFixest_fml()} returns a list of character strings, the names corresponding to the macro variable names, the character strings corresponding to their definition.
#'
#' @seealso
#' \code{\link[fixest]{xpd}} to make use of formula macros.
#'
#'
#'
setFixest_fml = function(..., reset = FALSE){

    check_arg(reset, get, "logical scalar")

    fml_macro = parse_macros(..., reset = reset, frame = parent.frame())

    options("fixest_fml_macro" = fml_macro)

}

#' @rdname setFixest_fml
getFixest_fml = function(){
    fml_macro = getOption("fixest_fml_macro")
    if(is.null(fml_macro)){
        options("fixest_fml_macro" = list())
        fml_macro = list()
    } else if(!is.list(fml_macro)){
        options("fixest_fml_macro" = list())
        warning("The value of getOption(\"fixest_fml_macro\") is not legal, it has been reset. Please use only 'setFixest_fml' to set formula macro variables.")
        fml_macro = list()
    }

    fml_macro
}



#' Default arguments for fixest estimations
#'
#' This function sets globally the default arguments of fixest estimations.
#'
#' @inheritParams feols
#' @inheritParams feNmlm
#' @inheritParams feglm
#'
#' @param reset Logical, default to \code{FALSE}. Whether to reset all values.
#'
#' @return
#' The function \code{getFixest_estimation} returns the currently set global defaults.
#'
#' @examples
#'
#' #
#' # Example: removing singletons is FALSE by default
#' #
#'
#' # => changing this default
#'
#' # Let's create data with singletons
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#' base$fe_singletons = as.character(base$species)
#' base$fe_singletons[1:5] = letters[1:5]
#'
#' res          = feols(y ~ x1 + x2 | fe_singletons, base)
#' res_noSingle = feols(y ~ x1 + x2 | fe_singletons, base, fixef.rm = "single")
#'
#' # New defaults
#' setFixest_estimation(fixef.rm = "single")
#' res_newDefault = feols(y ~ x1 + x2 | fe_singletons, base)
#'
#' etable(res, res_noSingle, res_newDefault)
#'
#' # Resetting the defaults
#' setFixest_estimation(reset = TRUE)
#'
#'
#'
setFixest_estimation = function(data = NULL, panel.id = NULL, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10, lean = FALSE, verbose = 0, warn = TRUE, combine.quick = NULL, demeaned = FALSE, mem.clean = FALSE, glm.iter = 25, glm.tol = 1e-8, reset = FALSE){

    check_arg_plus(fixef.rm, "match(none, perfect, singleton, both)")
    check_arg(fixef.tol, collin.tol, glm.tol, "numeric scalar GT{0}")
    check_arg(fixef.iter, glm.iter, "integer scalar GE{1}")
    check_arg(verbose, "integer scalar GE{0}")
    check_arg(lean, warn, demeaned, mem.clean, reset, "logical scalar")
    check_arg(combine.quick, "NULL logical scalar")
    check_arg(panel.id, "NULL character vector len(,2) no na | os formula")

    if(!missing(data)){
        check_arg(data, "NULL data.frame | matrix")
        if(!is.null(data)){
            data = deparse_long(substitute(data))
            class(data) = "default_data"
        }
    }

    # Getting the existing defaults
    opts = getOption("fixest_estimation")

    if(reset || is.null(opts)){
        opts = list()
    } else if(!is.list(opts)){
        warning("Wrong formatting of option 'fixest_estimation', all options are reset.")
        opts = list()
    }

    # Saving the default values
    mc = match.call()
    args_default = setdiff(names(mc)[-1], "reset")

    # NOTA: we don't allow delayed evaluation => all arguments must have hard values
    for(v in args_default){
        opts[[v]] = eval(as.name(v))
    }

    options(fixest_estimation = opts)

}

#' @rdname setFixest_estimation
getFixest_estimation = function(){
    getOption("fixest_estimation")
}

#### ............... ####
#### Startup         ####
####



#' Permanently removes the fixest package startup message
#'
#' Package startup messages can be very annoying, although sometimes they can be necessary. Use this function to prevent \code{fixest}'s package startup message from popping when loading. This will be specific to your current project.
#'
#' @param x Logical, no default. If \code{FALSE}, the package startup message is removed.
#'
#' @details
#' Note that this function is introduced to cope with the first \code{fixest} startup message (in version 0.9.0).
#'
#' This function works only with R >= 4.0.0. There are no startup messages for R < 4.0.0.
#'
fixest_startup_msg = function(x){

    check_arg(x, "logical scalar mbt")

    config_update("fixest_startup_msg", x)

}

initialize_startup_msg = function(startup_msg){
    # When new versions of the package are installed => we reset the display of the startup message
    # we need to keep track of the versions for which this default has been set

    # NOTA:
    # - the variable fixest_version is written when the user uses fixest_startup_msg()
    # - if this function returns TRUE, then it forces the msg to pop

    # NOTA:
    # - one problem is that I check the version using a local variable
    # specific to a project.
    # - this means that when one creates a new project, the message will necessarily pop!
    # - so I MUST turn off the message for newly created projects.
    # otherwise it would be so annoying.
    # - still => this is a problem if the person uses fixest for the first time
    # -> the project can be deemed old, while in fact fixest was never used
    # so startup messages weren't necessary (bc it would break nothing in the existing code)
    # -> new way: I look at the R files to check whether fixest is used:
    # - if TRUE: startup message
    # - if FALSE: nothing
    # - that's quite costly, but should happen only the very first time the package is attached

    # Note that we must return the value of 'fixest_startup_msg' since these are
    # updated only at session restart (and hence are not directly accessible)

    # message("fixest_startup_msg")

    if(getRversion() < "4.0.0"){
        # No startup message for version < 4.0
        # because there's no way to monitor the messages
        return(FALSE)
    }

    if(is_Rmarkdown()){
        # Never in Rmarkdown: too ugly
        return(FALSE)
    }

    if(is.null(find_project_path())){
        return(FALSE)
    }

    # message("getting version")

    previous_version = config_get("fixest_version")
    is_corrupt_version = !is.null(previous_version) && !is_pkg_version(previous_version)

    # message("version is ", previous_version)

    if(is.null(previous_version)){
        # compatibility with previous versions
        # message("trying to get version from renviron")
        previous_version = renvir_get("fixest_version")
    }

    current_version = fixest_version()

    if(!is_pkg_version(current_version)){
        # If we're here, it's a bug: this should NEVER happen
        return(FALSE)
    }

    if(!is_pkg_version(previous_version)){
        # We first update the version
        # message("updating the version")
        config_update("fixest_version", current_version)

        # message("Is fixest used? ", is_fixest_used())

        # Is it a new project? Or was fixest simply never used before?
        if(!is_corrupt_version && is_fixest_used()){
            # => message
            # Since I register versions since 0.9.0, this means that the
            # version of fixest used was anterior => all msgs should pop

            config_update("fixest_startup_msg", TRUE)
            return(TRUE)
        } else {
            # fixest was never used or the version was corrupt
            # => we don't show any message since it will not break any existing code
            config_update("fixest_startup_msg", FALSE)
            return(FALSE)
        }

        # message("updating done ")

    } else if(!identical(previous_version, current_version)){

        if(version2num(current_version) < version2num(previous_version)){
            # Can happen in projects shared in the cloud
            # In that case, we don't touch the startup message

            msg = paste0("The current project used 'fixest' version ", previous_version, ", but the current version is only ", current_version, ". Maybe update the package?")
            packageStartupMessage(fit_screen(msg, 1))

        } else {

            # A) we update the version
            config_update("fixest_version", current_version)

            # B) we reset the value of fixest_startup_msg
            #    only if the previous_version is anterior to the version that introduced the
            #    message (means the message SHOULD pop since it would be the first time)

            max_version_msg = names(startup_msg)[1]

            if(version2num(previous_version) < version2num(max_version_msg)){
                # You force a startup message even if it was turned off in a previous version

                # use case:
                # - v0.9.0: startup message, user uses fixest_startup_msg(FALSE)
                # - v0.10.0: new breaking changes, you want to inform the user even if he had set
                # fixest_startup_msg(FALSE) in v0.9.0
                #

                config_update("fixest_startup_msg", previous_version)
                return(previous_version)

            } else {
                # The previous version is already posterior to the last message
                # => no startup message any more

                config_update("fixest_startup_msg", FALSE)
                return(FALSE)
            }
        }
    }

    # If null, we'll get the value thanks to renvir_get("fixest_startup_msg")
    # but in some instances, it may be corrupt, so we fix it
    res = config_get("fixest_startup_msg")
    if(is.null(res)){
        # corrupt situation (can occur in dev)
        config_update("fixest_startup_msg", FALSE)
        return(FALSE)
    }

    return(res)
}

version2num = function(x){
    sum(as.numeric(strsplit(x, "\\.")[[1]]) * c(1e6, 1e3, 1))
}

fixest_version = function(){
    as.character(packageVersion("fixest"))
}

is_pkg_version = function(x){
    length(x) == 1 && is.character(x) && length(strsplit(x, "\\.")[[1]]) == 3
}

is_fixest_used = function(){
    # To return TRUE:
    # - fixest in the files
    # - + file saved > 7 days
    #
    # - if fixest but file saved < 7 days, very likely a new project

    # Only level 1 recursivity
    files = list.files(pattern = "\\.(r|R)$")
    dirs = c("./", list.dirs(recursive = FALSE))
    sub_files = unlist(lapply(dirs, list.files, pattern = "\\.(r|R)$", full.names = TRUE))
    file_extra = if(file.exists(".Rprofile")) ".Rprofile" else NULL

    files = c(files, file_extra, sub_files)
    files = files[!dir.exists(files)]

    if(length(files) == 0) return(FALSE)

    big_text = lapply(files, readLines, warn = FALSE)

    # we get the files that have fixest in them
    id_fixest = which(sapply(big_text, function(x) any(grepl("fixest", x, fixed = TRUE))))

    fixest_files = files[id_fixest]
    if(length(fixest_files) == 0) return(FALSE)

    now = Sys.time()

    for(f in fixest_files){
        f_created = file.mtime(f)
        if("POSIXt" %in% class(f_created)){
            d = as.numeric(difftime(now, f_created, units = "days"))
            if(d > 7){
                return(TRUE)
            }
        }
    }

    return(FALSE)
}

renvir_get = function(key){
    # Get the values of envir variables
    # we also evaluate them

    value_raw = Sys.getenv(key)

    if(value_raw == ""){
        return(NULL)
    }

    # Any default value should be able to be evaluated "as such"
    value_clean = gsub("__%%;;", "\n", value_raw)
    value_clean = gsub("&quot;", '"', value_clean)
    value_clean = gsub("&apos;", "'", value_clean)

    value = eval(str2lang(value_clean))

    return(value)
}

find_project_path = function(force = FALSE){
    # finds the root directory
    # we just look up the search path to find the root
    # Only works for Rstudio projects!

    past_path = "init"
    path = normalizePath(".", "/")

    is_found = FALSE
    i = 1
    nmax = 10
    while(past_path != path && i <= nmax){
        i = i + 1
        if(length(list.files(path, pattern = "Rproj$")) > 0){
            is_found = TRUE
            break
        } else {
            past_path = path
            path = dirname(path)
        }
    }

    proj_path = NULL
    if(is_found){
        proj_path = path
    }

    if(force && is.null(proj_path)){
        proj_path = normalizePath(".", "/")
    }

    proj_path
}

find_renviron = function(path = NULL){
    # Simply attaches .Renviron to the project path

    if(is.null(path)){
        proj_path = find_project_path()
        if(is.null(proj_path)) return(NULL)
    } else {
        if(!dir.exists(path)){
            if(file.exists(path)){
                path = dirname(path)
            } else {
                stop_up("The path provided in 'save' does not exist.", .up = 2)
            }
        }

        proj_path = path
    }

    file.path(proj_path, ".Renviron")
}

renvir_update = function(key, value){
    # Updates the .Renviron file
    # asks permission to the user => avoids messing up their workspace!
    # I was thinking to add an argument path, given by the user... but in fact no
    # the .Renviron works only at the Rstudio project level so making the user think
    # that giving a path for saving would help is misleading, since the .Renviron from
    # that path very likely wouldn't be loaded

    check_arg(key, "character scalar mbt")
    check_arg(value, "NULL mbt")

    renv_path = find_renviron()

    if(is.null(renv_path)){
        message("The 'save' feature only works with Rstudio projects. The root directory of the Rstudio project could not be found: settings cannot be saved at the project level, sorry.")
        return(NULL)
    }

    message("To save the settings at the project level 'fixest' needs to update the '.Renviron' file, currently located at:\n\n ", renv_path, "\n\n If the path indeed leads to your current project, do you give persmission? ")

    consent = readline("ok/y/yes to consent:")
    consent = tolower(trimws(consent))

    if(!consent %in% c("ok", "y", "ye", "yes")){
        message("aborting save")
        return(NULL)
    }

    if(file.exists(renv_path)){
        file = file(renv_path, "r", encoding = "UTF-8")

        renvir_raw = readLines(file)

        close(file)
    } else {
        renvir_raw = ""
    }

    all_keys = trimws(gsub("=.*", "", renvir_raw))

    do_write = TRUE
    if(is.null(value)){

        line_to_drop = all_keys == key
        if(any(line_to_drop)){
            renvir_raw = renvir_raw[!line_to_drop]
        } else {
            do_write = TRUE
        }

    } else {

        # we need to do some extra legwork... => sys env don't do quotes
        value_text = paste(deparse(value, width.cutoff = 500), collapse = "\n")
        value_text = gsub("\n", "__%%;;", value_text)
        value_text = gsub("\"", "&quot;", value_text)
        value_text = gsub("'", "&apos;", value_text)

        key_line = all_keys == key
        renvir_raw = c(renvir_raw[!key_line], paste0(key, " = ", value_text))
    }

    if(do_write){
        file = file(renv_path, "w", encoding = "UTF-8")

        renvir_raw = writeLines(renvir_raw, file)

        close(file)
    }


}

find_config_path = function(){

    if(getRversion() < "4.0.0"){
        return(NULL)
    }

    dir = tools::R_user_dir("fixest", "config")

    # We create the directory if needed
    if(!dir.exists(dir)){
        dir.create(dir, recursive = TRUE)
    }

    dir = normalizePath(dir, "/")

    file.path(dir, "fixest_config.csv")
}


config_update = function(key, value){

    if(getRversion() < "4.0.0"){
        return(NULL)
    }

    path = find_config_path()
    proj = find_project_path(force = TRUE)

    if(file.exists(path)){
        data = read.csv(path)
    } else {
        data = data.frame(proj = proj, fixest_version = fixest_version(), stringsAsFactors = FALSE)
    }

    if(!key %in% names(data)){
        data[[key]] = NA_character_
    }

    if(!proj %in% data$proj){
        row = data[1, , drop = FALSE]
        for(i in 1:ncol(row)) row[1, i] = NA
        row[1, 1] = proj
        data = rbind(data, row)
    }

    i = which(data$proj %in% proj)

    data[["fixest_version"]][i] = fixest_version()

    if(is.null(value)) value = "NULL"
    data[[key]][i] = as.character(value)

    write.csv(data, path, row.names = FALSE)
}

config_get = function(key){

    path = find_config_path()

    if(is.null(path) || !file.exists(path)){
        return(NULL)
    }

    data = read.csv(path)

    proj = find_project_path(force = TRUE)

    if(!proj %in% data$proj){
        return(NULL)
    }

    i = which(data$proj %in% proj)

    value = data[[key]][i]

    if(is.character(value) && value %in% c("NULL", "TRUE", "FALSE")){
        value = str2lang(value)
    }

    value
}


#### .................. ####
#### DOCUMENTATION DATA ####
####



#' Trade data sample
#'
#' This data reports trade information between countries of the European Union (EU15).
#'
#' @usage
#' data(trade)
#'
#' @format
#' \code{trade} is a data frame with 38,325 observations and 6 variables named \code{Destination}, \code{Origin}, \code{Product}, \code{Year}, \code{dist_km} and \code{Euros}.
#'
#' \itemize{
#' \item{Origin: 2-digits codes of the countries of origin of the trade flow.}
#' \item{Destination: 2-digits codes of the countries of destination of the trade flow.}
#' \item{Products: Number representing the product categories (from 1 to 20).}
#' \item{Year: Years from 2007 to 2016}
#' \item{dist_km: Geographic distance in km between the centers of the countries of origin and destination.}
#' \item{Euros: The total amount in euros of the trade flow for the specific year/product category/origin-destination country pair.}
#'
#' }
#'
#' @source
#' This data has been extrated from Eurostat on October 2017.
#'
#'
#'
"trade"

































































