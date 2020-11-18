
#' A print facility for \code{fixest} objects.
#'
#' This function is very similar to usual \code{summary} functions as it provides the table of coefficients along with other information on the fit of the estimation. The type of output is customizable by the user (using function \code{\link[fixest]{setFixest_print.type}}).
#'
#' @method print fixest
#'
#' @param x A \code{fixest} object. Obtained using the methods \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param n Integer, number of coefficients to display. By default, only the first 8 coefficients are displayed if \code{x} does not come from \code{\link[fixest]{summary.fixest}}.
#' @param type Either \code{"table"} (default) to display the coefficients table or \code{"coef"} to display only the coefficients. By default the value is \code{getFixest_print.type()} which can be permanently set with \code{\link[fixest]{setFixest_print.type}}).
#' @param ... Other arguments to be passed to \code{\link[fixest]{vcov.fixest}}.
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
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # displaying the results
#' #  (by default SEs are clustered if FEs are used)
#' print(est_pois)
#'
#' # By default the coefficient table is displayed.
#' #  If the user wished to display only the coefficents, use option type:
#' print(est_pois, type = "coef")
#'
#' # To permanently display coef. only, use setFixest_print.type:
#' setFixest_print.type("coef")
#' est_pois
#' # back to default:
#' setFixest_print.type("table")
#'
#'
print.fixest <- function(x, n, type = getFixest_print.type(), ...){

    # checking the arguments
    validate_dots(suggest_args = c("n", "type", "se", "cluster"),
                  valid_args = c("se", "cluster", "dof", "forceCovariance", "keepBounded"))

    # The objects from the estimation and the summary are identical, except regarding the vcov
	fromSummary = "cov.scaled" %in% names(x)

	# if NOT from summary, we consider the argument 'type'
	if(!fromSummary){
	    # checking argument type
	    check_arg_plus(type, "match(coef, table)")

	    if(type == "coef"){
	        print(coef(x))
	        return(invisible())
	    }
	}

	isNegbin = x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && x$family=="negbin")

	x = summary(x, fromPrint = TRUE, ...)

	check_arg(n, "integer scalar GE{1}")

	msgRemaining = ""
	nb_coef = length(coef(x)) - isNegbin
	if(missing(n) && is.null(x$n_print)){
		if(fromSummary){
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
	}

	if(isFALSE(x$convStatus)){
	    last_warn = getOption("fixest_last_warning")
	    if(is.null(last_warn) || (proc.time() - last_warn)[3] > 1){
	        if(x$method %in% c("femlm", "feNmlm", "fenegbin")){
	            warning("The optimization algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
	        } else if(x$method == "feols"){
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

	if(x$method %in% c("femlm", "feNmlm")){
		family_format = c(poisson="Poisson", negbin="Negative Binomial", logit="Logit", gaussian="Gaussian")
		msg = ifelse(is.null(x$call$NL.fml), "", "Non-linear ")
		half_line = paste0(msg, "ML estimation, family = ", family_format[x$family])
	} else if(x$method == "feglm") {
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

	cat(half_line, ", Dep. Var.: ", as.character(x$fml)[[2]], "\n", sep="")
	# cat(msg, "ML estimation, family = ", family_format[x$family], ", Dep. Var.: ", as.character(x$fml)[[2]], "\n", sep="")
	cat("Observations:", addCommas(x$nobs), "\n")
	if(!is.null(x$fixef_terms)){
	    terms_full = extract_fe_slope(x$fixef_terms)
	    fixef_vars = terms_full$fixef_vars

	    if(length(fixef_vars) > 0){
	        cat("Fixed-effects: ", paste0(fixef_vars, ": ", addCommas(x$fixef_sizes[fixef_vars]), collapse=",  "), "\n", sep="")
	    }

	    cat("Varying slopes: ", paste0(terms_full$slope_vars, " (", terms_full$slope_fe, ": ", addCommas(x$fixef_sizes[terms_full$slope_fe]), ")", collapse=",  "), "\n", sep="")

	} else {
	    if(!is.null(x$fixef_sizes)) cat("Fixed-effects: ", paste0(x$fixef_vars, ": ", addCommas(x$fixef_sizes), collapse=",  "), "\n", sep="")
	}


	if(is.null(x$onlyFixef)){

		cat("Standard-errors:", se.type, "\n")

	    last_line = paste0(msgRemaining, collinearity_msg)

		# The matrix of coefficients
		if(isNegbin){
			if(nrow(coeftable) == 2){
				new_table = coeftable[1, , FALSE]
			} else {
				new_table = coeftable[-nrow(coeftable), ]
			}

			myPrintCoefTable(head(new_table, n), lastLine = last_line)

			theta = coeftable[".theta", 1]
			noDispInfo = ifelse(theta > 1000, "(theta >> 0, no sign of overdispersion, you may consider a Poisson model)", "")
			cat("Over-dispersion parameter: theta =", theta, noDispInfo, "\n")
		} else {
			myPrintCoefTable(head(coeftable, n), lastLine = last_line)
		}
	}

	if(x$method == "feols"){
	    ll_format = numberFormatNormal(logLik(x))
	    cat("Log-likelihood:", ll_format,                    "  Adj. R2:", round(r2(x, "ar2"), 5), "\n")
	    if(!is.null(x$fixef_sizes) && is.null(x$onlyFixef)){
	        cat(sprintf("% *s", 16 + nchar(ll_format), " "), "R2-Within:", round(r2(x, "wr2"), 5), "\n")
	    }
	} else {
		bic_ll = formatBicLL(BIC(x), x$loglik)
		cat("Log-likelihood:", bic_ll$ll,  "Adj. Pseudo-R2:", round(x$pseudo_r2, 5), "\n")
		cat("           BIC:", bic_ll$bic, "  Squared Cor.:", round(x$sq.cor, 5), "\n")
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
#' @method summary fixest
#'
#' @param se Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: \code{se = "cluster"}, otherwise \code{se = "standard"}. Note that this argument can be implicitly deduced from the argument \code{cluster}.
#' @param cluster Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over \code{var1} and \code{var2} contained in the data.frame \code{base} used for the estimation. All the following \code{cluster} arguments are valid and do the same thing: \code{cluster = base[, c("var1", "var2")]}, \code{cluster = c("var1", "var2")}, \code{cluster = ~var1+var2}. If the two variables were used as clusters in the estimation, you could further use \code{cluster = 1:2} or leave it blank with \code{se = "twoway"} (assuming \code{var1} [resp. \code{var2}] was the 1st [res. 2nd] cluster). You can interact two variables using \code{^} with the following syntax: \code{cluster = ~var1^var2} or \code{cluster = "var1^var2"}.
#' @param object A \code{fixest} object. Obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param dof An object of class \code{dof.type} obtained with the function \code{\link[fixest]{dof}}. Represents how the degree of freedom correction should be done.You must use the function \code{\link[fixest]{dof}} for this argument. The arguments and defaults of the function \code{\link[fixest]{dof}} are: \code{adj = TRUE}, \code{fixef.K="nested"}, \code{cluster.adj = TRUE}, \code{cluster.df = "conventional"}, \code{t.df = "conventional"}, \code{fixef.force_exact=FALSE)}. See the help of the function \code{\link[fixest]{dof}} for details.
#' @param .vcov A user provided covariance matrix or a function computing this matrix. If a matrix, it must be a square matrix of the same number of rows as the number of variables estimated. If a function, it must return the previsouly mentioned matrix.
#' @param lean Logical, default is \code{FALSE}. Used to reduce the (memory) size of the summary object. If \code{TRUE}, then all objects of length N (the number of observations) are removed from the result. Note that some \code{fixest} methods may consequently not work when applied to the summary.
#' @param forceCovariance (Advanced users.) Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option to force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded (Advanced users -- \code{feNmlm} with non-linear part and bounded coefficients only.) Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param n Integer, default is missing (means Inf). Number of coefficients to display when the print method is used.
#' @param ... Only used if the argument \code{.vocv} is provided and is a function: extra arguments to be passed to that function.
#'
#' @section Compatibility with \pkg{sandwich} package:
#' The VCOVs from \code{sandwich} can be used with \code{feols}, \code{feglm} and \code{fepois} estimations. If you want to have a \code{sandwich} VCOV when using \code{summary.fixest}, you can use the argument \code{.vcov} to specify the VCOV function to use (see examples).
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
#' sum_standard = summary(est_pois, se = "standard")
#' sum_hetero   = summary(est_pois, se = "hetero")
#' sum_oneway   = summary(est_pois, se = "cluster")
#' sum_twoway   = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' etable(sum_standard, sum_hetero, sum_oneway, sum_twoway, sum_threeway)
#'
#' # Alternative ways to cluster the SE:
#'
#' # two-way clustering: Destination and Product
#' # (Note that arg. se = "twoway" is implicitly deduced from the argument cluster)
#' summary(est_pois, cluster = c("Destination", "Product"))
#' summary(est_pois, cluster = trade[, c("Destination", "Product")])
#' summary(est_pois, cluster = list(trade$Destination, trade$Product))
#' summary(est_pois, cluster = ~Destination+Product)
#' # Since Destination and Product are used as fixed-effects, you can also use:
#' summary(est_pois, cluster = 2:3)
#'
#' # You can interact the clustering variables "live" using the var1 ^ var2 syntax.
#'
#' summary(est_pois, cluster = "Destination^Product")
#' summary(est_pois, cluster = ~Destination^Product)
#' # Equivalent to
#' summary(est_pois, cluster = paste(trade$Destination, trade$Product))
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
summary.fixest <- function(object, se, cluster, dof = getFixest_dof(), .vcov, lean = FALSE, forceCovariance = FALSE, keepBounded = FALSE, n, ...){
	# computes the clustered SD and returns the modified vcov and coeftable

	if(!is.null(object$onlyFixef)){
		# means that the estimation is done without variables
		return(object)
	}

	dots = list(...)

	# If cov.scaled exists => means that it has already been computed
	if(!is.null(object$cov.scaled) && "fromPrint" %in% names(dots)) return(object)

	# Checking arguments in ...
	if(!any(c("fromPrint", "nframes_up") %in% names(match.call()))){
	    # condition means NOT internal call => thus client call
	    if(missing(.vcov) || !is.function(.vcov)){
	        validate_dots(suggest_args = c("se", "cluster", "dof"))
	    }
	}

	check_arg(lean, "logical scalar")
	check_arg(n, "integer scalar GE{1}")
	if(!missing(n)){
	    object$n_print = n
	}

	if(is.null(dots$nframes_up)){
		nframes_up = 1 + !is.null(dots$fromPrint)
	} else {
		nframes_up = dots$nframes_up + 1
	}

	# The new VCOV
	if(!missnull(.vcov)){
	    n_coef = length(object$coefficients)
	    check_arg(.vcov, "square numeric matrix nrow(value) | function", .value = n_coef)

	    vcov_name = "Custom"
	    if(is.function(.vcov)){
	        arg_names = formalArgs(.vcov)
	        # we contruct the call
	        dots = list(...)

	        if(".vcov_args" %in% names(dots)){
	            # internal call
	            vcov_name = dots$vcov_name
	            dots = dots$.vcov_args
	        } else {
	            mc = match.call()
	            vcov_name = deparse_long(mc$.vcov)
	        }

	        vcov_name = gsub("sandwich::", "", vcov_name, fixed = TRUE)

	        # We shouldn't have a prior on the name of the first argument
	        dots[[arg_names[1]]] = as.name("object")
	        if("cluster" %in% arg_names && !missing(cluster)){
	            dots[["cluster"]] = as.name("cluster")
	        }

	        vcov = do.call(.vcov, dots)

	        check_value(vcov, "square numeric matrix nrow(value)", .value = n_coef,
	                    .message = paste0("If argument '.vcov' is to be a function, it should return a square numeric matrix of the same dimension as the number of coefficients (here ", n_coef, ")."))

	    } else {
	        # square matrix
	        vcov = .vcov
	    }

	    # We add the type of the matrix
	    attr(vcov, "type") = vcov_name

	    warn_ignore = c()
	    if(!missnull(se)) warn_ignore = "se"
	    if(length(warn_ignore) > 0){
	        warning("Since argument '.vcov' is provided, the argument", enumerate_items(warn_ignore, "s.quote.is"), " ignored.")
	    }

	} else {
	    vcov = vcov(object, se=se, cluster=cluster, dof=dof, forceCovariance = forceCovariance, keepBounded = keepBounded, nframes_up = nframes_up, attr = TRUE)
	}


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
	zvalue <- object$coefficients/se
	if(object$method %in% "feols" || (object$method %in% "feglm" && !object$family$family %in% c("poisson", "binomial"))){

	    # I have renamed t.df into G
	    t.df = attr(vcov, "G")

	    if(!is.null(t.df)){
	        pvalue <- 2*pt(-abs(zvalue), max(t.df - 1, 1))
	    } else {
	        pvalue <- 2*pt(-abs(zvalue), max(object$nobs - object$nparams, 1))
	    }

	} else {
	    pvalue <- 2*pnorm(-abs(zvalue))
	}


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
	coeftable[, 2] = se_format
	coeftable[, 3] = zvalue
	coeftable[, 4] = pvalue

	attr(coeftable, "type") = attr(se, "type") = attr(vcov, "type")

	object$cov.scaled = vcov
	object$coeftable = coeftable
	object$se = se

	if(lean){
	    object[c("fixef_id", "residuals", "fitted.values", "scores", "sumFE", "slope_variables_reordered")] = NULL
	}

	return(object)
}

#' @rdname summary.fixest
"summ"

summ = function(object, se, cluster, dof = getFixest_dof(), forceCovariance = FALSE, keepBounded = FALSE, ...){

    # we reiterate the call
    mc = match.call()
    mc[[1]] = as.name("summary")

    eval(mc, parent.frame())
}


#' R2s of \code{fixest} models
#'
#' Reports different R2s for \code{fixest} estimations (e.g. \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}).
#'
#' @param x A \code{fixest} object, e.g. obtained with function \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}.
#' @param type A character vector representing the R2 to compute. The R2 codes are of the form: "wapr2" with letters "w" (within), "a" (adjusted) and "p" (pseudo) possibly missing. E.g. to get the regular R2: use \code{type = "r2"}, the within adjusted R2: use \code{type = "war2"}, the pseudo R2: use \code{type = "pr2"}, etc. Use \code{"sq.cor"} for the squared correlation. By default, all R2s are computed.
#' @param full_names Logical scalar, default is \code{FALSE}. If \code{TRUE} then names of the vector in output will have full names instead of keywords (e.g. \code{Squared Correlation} instead of \code{sq.cor}, etc).
#'
#' @details
#' For R2s with no theoretical justification, like e.g. regular R2s for maximum likelihood models -- or within R2s for models without fixed-effects, NA is returned. The single measure to possibly compare all kinds of models is the squared correlation between the dependent variable and the expected predictor.
#'
#' The pseudo-R2 is also returned in the OLS case, it corresponds to the pseudo-R2 of the equivalent GLM model with a Gaussian family.
#'
#' For the adjusted within-R2s, the adjustment factor is \code{(n - nb_fe) / (n - nb_fe - K)} with \code{n} the number of observations, \code{nb_fe} the number of fixed-effects and \code{K} the number of variables.
#'
#' @return
#' Returns a named vector.
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
#' est = feols(log(Euros) ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Squared correlation:
#' r2(est, "sq.cor")
#'
#' # "regular" r2:
#' r2(est, "r2")
#'
#' # pseudo r2 (equivalent to GLM with Gaussian family)
#' r2(est, "pr2")
#'
#' # adjusted within r2
#' r2(est, "war2")
#'
#' # all four at once
#' r2(est, c("sq.cor", "r2", "pr2", "war2"))
#'
#' # same with full names instead of codes
#' r2(est, c("sq.cor", "r2", "pr2", "war2"), full_names = TRUE)
#'
r2 = function(x, type = "all", full_names = FALSE){
	# p: pseudo
	# w: within
	# a: adjusted

    check_arg(full_names, "logical scalar")

	if(!"fixest" %in% class(x)){
		stop("Only 'fixest' objects are supported.")
	}

	check_arg(type, "character vector no na", .message = "Argument 'type' must be a character vector (e.g. type = c(\"sq.cor\", \"r2\", \"pr2\")). (a: adjused, p: pseudo, w: within.)")

	# type_allowed next => ("count", "acount") ?
    dict_names = c("sq.cor" = "Squared Correlation", "r2" = "R2", "ar2" = "Adjusted R2", "pr2" = "Pseudo R2", "apr2" = "Adjusted Pseudo R2", "wr2" = "Within R2", "war2" = "Adjusted Within R2", "wpr2" = "Within Pseudo R2", "wapr2" = "Adjusted Within Pseudo R2")
	type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2")
	types_alias = c(par2 = "apr2", awr2 = "war2", pwr2 = "wpr2", wpar2 = "wapr2", pwar2 = "wapr2", pawr2 = "wapr2", apwr2 = "wapr2", awpr2 = "wapr2")
	if("all" %in% type){
		type_all = type_allowed
	} else {
		type_all = tolower(type)
		pblm = setdiff(type_all, c(type_allowed, names(types_alias)))
		if(length(pblm) > 0){
			stop("The r2 type", enumerate_items(pblm, "quote.s.is"), " not valid.")
		}
	}

	type_all = dict_apply(type_all, types_alias)

	is_ols = x$method == "feols"
	isFixef = "fixef_vars" %in% names(x)
	n = nobs(x)

	res = rep(NA, length(type_all))
	for(i in seq_along(type_all)){
		myType = type_all[i]

		if(myType == "sq.cor"){
		    res[i] = x$sq.cor
		    next
		}

		if(!grepl("p", myType) && !is_ols){
		    # non pseudo R2 not valid for non ols
		    next
		}

		if(grepl("w", myType) && !isFixef){
		    # within R2 not valid for models without FE
		    next
		}

		adj = grepl("a", myType)
		pseudo = grepl("p", myType)
		within = grepl("w", myType) && isFixef
		ifNullNA = function(x) ifelse(is.null(x), NA, x)
		if(within && isFixef){

		    if(isTRUE(x$onlyFixef)) next

		    if(is.null(x$ssr_fe_only) && !is.null(x$fixef_vars)){
		        # This is the case of feglm where there were no fe_only model estimated
		        # => we need to compute the FE model first

		        # 2019-11-26: now self contained call (no need for outer frame evaluation)
		        # res_fe = update(x, .~1|., glm.tol = 1e-2, fixef.tol = 1e-3, nframes = 2)

		        # constructing the data
		        newdata = cbind(data.frame(y = x$y), as.data.frame(x$fixef_id))
		        if(!is.null(x$fixef_terms)){
		            newdata = cbind(newdata, as.data.frame(x$slope_variables))
		        }
		        # Fe/slope only formula
		        new_fml = merge_fml(y ~ 1, x$fml_all$fixef)

		        # The fact that weights = x[["weights"]] is on purpose -- don't touch it
		        # same for offset = x[["offset"]] -- otherwise error is thrown in fixest_env
		        # x$family$family is also normal
		        res_fe = feglm(fml = new_fml, data = newdata, glm.tol = 1e-2, fixef.tol = 1e-3, family = x$family$family, weights = x[["weights"]], offset = x[["offset"]])

		        x$ssr_fe_only = cpp_ssq(res_fe$residuals)
		        x$ll_fe_only = logLik(res_fe)
		    }

			# within
			df_k = ifelse(adj, length(coef(x)), 0)
			if(pseudo){
				ll_fe_only = ifNullNA(x$ll_fe_only)
				ll = logLik(x)
				res[i] = 1 - (ll - df_k) / ll_fe_only
			} else {
				ssr_fe_only = ifNullNA(x$ssr_fe_only)
				nb_fe = x$nparams - length(coef(x))
				res[i] = 1 - x$ssr / ssr_fe_only * (n - nb_fe) / (n - nb_fe - df_k)
			}
		} else {
			df_k = ifelse(adj, x$nparams, 1 - pseudo)
			if(pseudo){
				ll_null = x$ll_null
				ll = logLik(x)
				if(adj){
					res[i] = 1 - (ll - x$nparams + 1) / ll_null
				} else {
					res[i] = 1 - ll / ll_null
				}
			} else {
				ssr_null = x$ssr_null
				df.intercept = 1 * (isFixef || any(grepl("(Intercept)", names(x$coefficients), fixed = TRUE)))
				res[i] = 1 - x$ssr / ssr_null * (n - df.intercept) / (n - df_k)
			}

		}
	}

	if(full_names){
	    names(res) = dict_apply(type_all, dict_names)
	} else {
	    names(res) = type_all
	}


	res
}


#' Obtain various statistics from an estimation
#'
#' Set of functions to directly extract some commonly used statistics, like the p-value or the table of coefficients, from estimations. This was first implemented for \code{fixest} estimations, but has some support for other models.
#'
#' @param object An estimation. For example obtained from \code{\link[fixest]{feols}}.
#' @param se [Fixest specific.] Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{hetero}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: \code{se = "cluster"}, otherwise \code{se = "standard"}. Note that this argument can be implicitly deduced from the argument \code{cluster}.
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
#' #  and using ctable(), the alias to coeftable()
#'
#' ctable(est, cluster = ~Origin + Product)
#' se(est, cluster = ~Origin + Product)
#' pvalue(est, cluster = ~Origin + Product)
#' tstat(est, cluster = ~Origin + Product)
#'
#' # Or you can cluster only once:
#' est_sum = summary(est, cluster = ~Origin + Product)
#' ctable(est_sum)
#' se(est_sum)
#' tstat(est_sum)
#' pvalue(est_sum)
#'
#'
#'
coeftable = ctable = function(object, se, cluster, ...){
    # We don't explicitely refer to the other arguments

    # We make the same call to summary if necessary
    mc = match.call()

    IS_FIXEST = "fixest" %in% class(object)

    if(!any(grepl("summary", class(object))) && (!IS_FIXEST || any(!names(mc) %in% c("", "object")) || !"cov.scaled" %in% names(object))){
        # We call summary
        mc[[1]] = as.name("summary")
        object = eval(mc, parent.frame())
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
            stop("No coefficient table found. Was the 'object' really an estimation?")
        }

    }

    res
}

#' @rdname coeftable
"ctable"

#' @describeIn coeftable Extracts the p-value of an estimation
pvalue = function(object, se, cluster, ...){

    mc = match.call()
    mc[[1]] = as.name("coeftable")

    mat = eval(mc, parent.frame())

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4). You can investigate the problem using function ctable().")
    }

    res = mat[, 4]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
    }
    res
}

#' @describeIn coeftable Extracts the t-statistics of an estimation
tstat = function(object, se, cluster, ...){

    mc = match.call()
    mc[[1]] = as.name("coeftable")

    mat = eval(mc, parent.frame())

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4). You can investigate the problem using function ctable().")
    }

    res = mat[, 3]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
    }
    res
}

#' @describeIn coeftable Extracts the standard-error of an estimation
se = function(object, se, cluster, ...){

    mc = match.call()
    mc[[1]] = as.name("coeftable")

    mat = eval(mc, parent.frame())

    if(ncol(mat) != 4){
        stop("No appropriate coefficient table found (number of columns is ", ncol(mat), " instead of 4). You can investigate the problem using function ctable().")
    }

    res = mat[, 2]
    if(is.null(names(res))) {
        names(res) = rownames(mat)
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
summary.fixest.fixef = function(object, n=5, ...){
	# This function shows some generic information on the fixed-effect coefficients

    # checking arguments in dots
    validate_dots(suggest_args = "n")

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
		if(!isRegular){
			cat("NOTE: The fixed-effects are NOT regular, so cannot be straighforwardly interpreted.\n")
		}
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
fixef.fixest = function(object, notes = getFixest_notes(), sorted = TRUE, ...){
	# object is a fixest object
	# This function retrieves the dummies

    check_arg(notes, sorted, "logical scalar")

    # Checking the arguments
    validate_dots(valid_args = "fixef.tol")

    dots = list(...)
    fixef.tol = dots$fixef.tol
    check_value_plus(fixef.tol, "NULL{1e-5} numeric scalar GT{0} LT{1}")

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


	    table_id_I = as.integer(unlist(lapply(fe_id_list, table), use.names = FALSE))

	    S_demean <- cpp_demean(y = S, X_raw = 0, n_vars_X = 0, r_weights = 0, iterMax = 1000L,
	                           diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
	                           fe_id_list = fe_id_list, table_id_I = table_id_I,
	                           slope_flag_Q = slope_flag, slope_vars_list = slope_variables,
	                           r_init = 0, nthreads = 1L, save_fixef = TRUE)

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

	        dumMat <- matrix(unlist(my_dum, use.names = FALSE), N, Q_fe) - 1
	        orderCluster <- matrix(unlist(lapply(my_dum, order), use.names = FALSE), N, Q_fe) - 1

	        nbCluster = sapply(my_dum, max)

	        fixef_values_tmp <- cpp_get_fe_gnl(Q_fe, N, rep(1, N), dumMat, nbCluster, orderCluster)

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

		dumMat <- matrix(unlist(id_dummies_vect), N, Q) - 1
		orderCluster <- matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

		nbCluster = sapply(fixef_id, max)

		fixef_values <- cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

		# the information on the references
		nb_ref = fixef_values[[Q+1]]
		fixef_values[[Q+1]] = NULL
	}

	# now saving & adding information
	all_clust = list()
	Q_all = ifelse(isSlope, length(fixef_terms), Q)
	for(i in 1:Q_all){
	    # We put it inthe right order, if requested
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
			message("NOTE: The fixed-effects are not regular, they cannot be straightforwardly interpreted.")
		}
	}

	# Family information
	attr(all_clust, "exponential") = FALSE
	if(object$method == "femlm" && object$family %in% c("poisson", "negbin")){
		attr(all_clust, "exponential") = TRUE
	} else if(object$method == "feglm" && object$family$link == "log"){
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
    validate_dots(suggest_args = "n")

	Q = length(x)

	mfrow = as.character(c(11, 12, 22, 22, 32, 32, 33, 33))

	fixef_names = names(x)
	slope_flag = grepl("\\[", fixef_names)

	if(Q > 1 && sum(attr(x, "references")[!slope_flag]) > sum(!slope_flag)-1){
		warning("The fixed-effects are not regular, they cannot be straightforwardly interpreted.", call. = FALSE)
	}

	# modification par:
	opar <- par(no.readonly =TRUE)
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

	if(class(x) != "fixest"){
		stop("Argument 'x' must be a fixest object.")
	}

	# I) (linear) collinearity with fixed-effects
	# II) (linear) multi collinearity
	# III) (non-linear) overidentification

	# flags
    isFixef = !is.null(x$fixef_id)
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
	if(grepl("[^:]::[^:]", deparse_long(rhs_fml[[3]]))){
	    new_fml = expand_interactions(rhs_fml)
	    linear.varnames = all.vars(new_fml[[3]])
	} else {
	    linear.varnames = all.vars(rhs_fml[[3]])
	}
	isLinear = length(linear.varnames) > 0

	NL_fml = x$NL.fml
	isNL = !is.null(NL_fml)
	coef = x$coefficients

	# Getting the data
	data = NULL
	dataName = x$call$data
	try(data <- eval(dataName, parent.frame()))

	if(is.null(data)){
		stop("To apply function 'collinearity', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", deparse_long(dataName), ").")
	}

	if(is.matrix(data)){
	    data = as.data.frame(data)
	} else {
	    class(data) = "data.frame"
	}

	if(isFE){
		linear_fml = update(linear_fml, ~ . + 1)
	}

	# Panel setup
	if(check_lag(linear_fml)){
	    if(!is.null(x$panel.info)){
	        if(is.null(attr(data, "panel_info"))){
	            # We try to recreate the panel
	            if(any(!names(x$panel.info) %in% c("", "data", "panel.id"))){
	                # This was NOT a standard panel creation
	                stop("The original data set was a fixest_panel, now it isn't any more. Please restore the original data to a panel to perform the collinearity check. NOTA: the original call to panel was:\n", deparse_long(x$panel.info))
	            } else {
	                panel__meta__info = panel_setup(data, x$panel.id, from_fixest = TRUE)
	            }
	        } else {
	            panel__meta__info = attr(data, "panel_info")
	        }
	    } else {
	        panel__meta__info = panel_setup(data, x$panel.id, from_fixest = TRUE)
	    }
	}

	if(isLinear || isFixef || "(Intercept)" %in% names(coef)){
		# linear.matrix = model.matrix(linear_fml, data)
		linear.matrix = fixest_model_matrix(rhs_fml, data)
	}

	if(!is.null(x$obsRemoved)){
	    linear.matrix = linear.matrix[-x$obsRemoved, , drop = FALSE]
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
	        stop("To check collinearity, we need to fetch some variables in the original dataset (", deparse_long(dataName), "). However, the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset any more.")
	    }

	    slope_var_list = list()
	    for(i in 1:length(slope_vars_unik)){
	        variable = all.vars(parse(text = slope_vars_unik[i]))

	        # as.numeric => we'll use cpp so required
	        svar = as.numeric(as.vector(eval(parse(text = slope_vars_unik[i]), data)))
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
		ccat("simple with cluster:")
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
				collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
				message = paste0("Variable '", dict_name[v], "' is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ".")

				print(message)
				return(invisible(message))
			}
		}
		ccat("OK")
	}

	#
	# II.b) perfect multicollinearity + cluster
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
#' @param fml Either a formula of the type \code{var1 + ... + var[N] ~ treat} or \code{var1 + ... + var[N] ~ treat | post}. Either a data.frame/matrix containing all the variables for which the means are to be computed (they must be numeric of course). Both the treatment and the post variables must contain only exactly two values. You can use a point to select all the variables of the data set: \code{. ~ treat}.
#' @param base A data base containing all the variables in the formula \code{fml}.
#' @param treat_var Only if argument \code{fml} is *not* a formula. The vector identifying the treated and the control observations (the vector can be of any type but must contain only two possible values). Must be of the same length as the data.
#' @param post_var Only if argument \code{fml} is *not* a formula. The vector identifying the periods (pre/post) of the observations (the vector can be of any type but must contain only two possible values). The first value (in the sorted sense) of the vector is taken as the pre period. Must be of the same length as the data.
#' @param treat_first Which value of the 'treatment' vector should appear on the left? By default the max value appears first (e.g. if the treatment variable is a 0/1 vector, 1 appears first).
#' @param tex Should the result be displayed in Latex? Default is \code{FALSE}. Automatically set to \code{TRUE} if the table is to be saved in a file using the argument \code{file}.
#' @param treat_dict A character vector of length two. What are the names of the treated and the control? This should be a dictionnary: e.g. \code{c("1"="Treated", "0" = "Control")}.
#' @param dict A named character vector. A dictionnary between the variables names and an alias. For instance \code{dict=c("x"="Inflation Rate")} would replace the variable name \code{x} by \dQuote{Inflation Rate}.
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
did_means = function(fml, base, treat_var, post_var, tex = FALSE, treat_dict, dict = getFixest_dict(), file, replace = FALSE, title, label, raw = FALSE, indiv, treat_first, prepostnames = c("Before", "After"), diff.inv = FALSE){
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

            indiv_var = try(eval(parse(text = indiv_varname), base), silent = TRUE)
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
            stop("The formula must be of the type 'var1 + ... + var[N] ~ treat' or 'var1 + ... + var[N] ~ treat | post'.")
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

            post_var = try(eval(pipe, base), silent = TRUE)
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
            for(i in qui) res$vars[i] = escape_latex(dict[res$vars[i]], up = 1)
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
#' Treat a variable as a factor, or interacts a variable with another treated as a factor. Values to be dropped/kept from the factor can be easily set.
#'
#' @param var A vector to be interacted with \code{f}. If the other argument \code{f} is missing, then this vector will be treated as the argument \code{f}.
#' @param f A vector (of any type) that will be treated as a factor. Must be of the same length as \code{var} if \code{var} is not missing.
#' @param ref A single value that belongs to the interacted variable (\code{f}). Can be missing.
#' @param drop A vector of values that belongs to the factor variable (\code{f}). If provided, all values from \code{f} that match \code{drop} will be removed.
#' @param keep A vector of values that belongs to the factor variable (\code{f}). If provided, only the values from \code{f} that match \code{keep} will be kept.
#'
#' @return
#' It returns a matrix with number of rows the length of \code{var}. The number of columns is equal to the number of cases contained in \code{f} minus the reference.
#'
#' @section Shorthand in \code{fixest} estimations:
#' In \code{fixest} estimations, instead of using \code{i(var, f, ref)}, you can instead use the following writing \code{var::f(ref)}. Note that this way of doing interactions is not endorsed any more and will likely be deprecated in the future.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{coefplot}} to plot interactions, \code{\link[fixest]{feols}} for OLS estimation with multiple fixed-effects.
#'
#' @examples
#'
#' #
#' # Simple illustration
#' #
#'
#' x = 1:10
#' y = rep(1:4, 3)[1:10]
#'
#' # interaction
#' cbind(x, y, i(x, y, 1))
#'
#' # without interaction
#' cbind(x, y, i(y, ref = 1))
#'
#' # You can interact factors too
#' z = rep(c("a", "b", "c"), c(5, 3, 2))
#' data.frame(z, y, i(z, y))
#'
#' #
#' # In fixest estimations
#' #
#'
#' data(base_did)
#' # We interact the variable 'period' with the variable 'treat'
#' est_did = feols(y ~ x1 + i(treat, period, 5) | id + period, base_did)
#'
#' # => special treatment in coefplot
#' coefplot(est_did)
#'
#' # Using i() for factors
#' est_bis = feols(y ~ x1 + i(period, keep = 3:6) + i(treat, period, 5) | id, base_did)
#'
#' coefplot(est_bis, only.inter = FALSE)
#'
#' # => special treatment in etable
#' etable(est_bis, dict = c("6" = "six"))
#'
#'
i = interact = function(var, f, ref, drop, keep){
    # Used to create interactions

    # gt = function(x) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s\n", sep = "")
    # t0 = proc.time()

    mc = match.call()

    if(missing(var) && missing(f)){
        stop("You must provide at least one of the two arguments: 'var' or 'f'.")
    }

    # Finding if it's a call from fixest
    FROM_FIXEST = sys.nframe() > 5 && any(sapply(sys.calls(), function(x) any(grepl("fixest", x[[1]], fixed = TRUE))))

    # General checks
    check_arg(var, "vector")
    check_arg(f, "vector")

    IS_INTER = TRUE
    if(missing(f) || missing(var)){
        # => i() is used to create a factor
        IS_INTER = FALSE

        if(missing(f)){
            f = var
            fe_name = deparse_long(mc$var)
        } else {
            fe_name = deparse_long(mc$f)
        }

    } else {
        var_name = deparse_long(mc$var)
        fe_name = deparse_long(mc$f)

        if(length(var) != length(f)){
            if(grepl("^[[:alpha:]\\.][[:alnum:]\\._]*:[[:alpha:]\\.][[:alnum:]\\._]*$", var_name)){
                info = strsplit(var_name, ":")[[1]]
                stop("In interact(): When 'var' is equal to a product, please use I(", info[1], "*", info[2], ") instead of ", var_name, ".")
            } else {
                stop("The arguments 'var' and 'f' must be of the same length (currently ", length(var), " vs ", length(f), ").")
            }
        }
    }

    # The NAs + recreation of f if necessary
    IS_FACTOR_INTER = FALSE
    if(IS_INTER){

        is_na_all = is.na(var) | is.na(f)

        if(!is.numeric(var)){
            IS_FACTOR_INTER = TRUE
            # It's an interaction between factors
            f_new = rep(NA_character_, length(f))
            f_new[!is_na_all] = paste0(var[!is_na_all], "__%%__", f[!is_na_all])
            f = f_new
            IS_INTER = FALSE
        }

    } else {
        is_na_all = is.na(f)
    }

    if(!IS_INTER){
        # neutral var in C code
        var = 1
    }

    # QUFing

    info = to_integer(f, add_items = TRUE, items.list = TRUE, sorted = TRUE)
    fe_num = info$x
    items = info$items

    check_arg(ref, "logical scalar | charin", .choices = items, .message = paste0("Argument 'ref' must be a single element of the variable '", fe_name, "'."))
    check_arg(drop, keep, "multi charin", .choices = items, .message = paste0("Argument '__ARG__' must consist of elements of the variable '", fe_name, "'."))

    no_rm = TRUE
    any_ref = FALSE
    id_drop = c()
    if(!missing(ref)){
        if(is.logical(ref)){
            if(ref == TRUE){
                # We always delete the first value
                any_ref = TRUE
                id_drop = ref = which(items == items[1])
            }
        } else {
            any_ref = TRUE
            id_drop = ref = which(items %in% ref)
        }
    }

    if(!missing(drop)){
        id_drop = c(id_drop, which(items %in% drop))
    }

    if(!missing(keep)){
        id_drop = c(id_drop, which(!items %in% keep))
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
    } else {
        items_name = items
    }

    if(FROM_FIXEST){
        # Pour avoir des jolis noms c'est un vrai gloubiboulga,
        # mais j'ai pas trouve plus simple...
        if(IS_FACTOR_INTER){
            name_split = strsplit(items_name, "__%%__", fixed = TRUE)
            var_items = sapply(name_split, function(x) x[1])
            f_items = sapply(name_split, function(x) x[2])
            col_names = paste0("__CLEAN__", var_name, "::", var_items, ":", fe_name, "::", f_items)

        } else if(IS_INTER){
            col_names = paste0("__CLEAN__", var_name, ":", fe_name, "::", items_name)

        } else {
            col_names = paste0("__CLEAN__", fe_name, "::", items_name)
        }
    } else {

        if(IS_FACTOR_INTER){
            name_split = strsplit(items_name, "__%%__", fixed = TRUE)
            items_name = sapply(name_split, paste, collapse = ":")
        }

        col_names = items_name
    }

    res = cpp_factor_matrix(fe_num, is_na_all, who_is_dropped, var, col_names)
    # res => matrix with...
    #  - NAs where appropriate
    #  - appropriate number of columns
    #  - interacted if needed
    #


    # We send the information on the reference
    if(FROM_FIXEST && IS_INTER){
        opt = getOption("fixest_interaction_ref")
        if(is.null(opt)){

            opt = list(fe_type = class(f))

            qui_drop = (1:length(items)) %in% id_drop

            if(any_ref){
                qui_drop[ref] = FALSE
                opt$items = items[!qui_drop]

                is_ref = rep(FALSE, length(opt$items))
                is_ref[sum(!qui_drop[1:ref])] = TRUE

                opt$is_ref = is_ref
            } else {
                opt$items = items[!qui_drop]
                opt$is_ref = rep(FALSE, length(opt$items))
            }

            opt$prefix = paste0(var_name, ":", fe_name)

            options("fixest_interaction_ref" = opt)
        }
    }

    res
}

#' @rdname i
"interact"

i_ref = function(var, f, ref, drop, keep){

    mc = match.call()

    mc[[1]] = as.name("i")

    if(!all(c("var", "f") %in% names(mc)) && !any(c("ref", "drop", "keep") %in% names(mc))){
        mc$ref = TRUE
    }

    return(deparse_long(mc))
}

i_noref = function(var, f, ref, drop, keep){

    mc = match.call()

    mc[[1]] = as.name("i")

    mc$ref = mc$drop = mc$keep = NULL

    return(deparse_long(mc))
}


#' @rdname setFixest_fml
xpd = function(fml, ...){
    check_arg(fml, .type = "formula mbt")

    macros = parse_macros(..., from_xpd = TRUE)

    if(length(macros) == 0) return(fml)

    qui = which(names(macros) %in% all.vars(fml))
    if(length(qui) > 0){
        fml_dp = deparse_long(fml)
        # We need to add a lookafter assertion: otherwise if we have ..ctrl + .. ctrl_long, there will be a replacement in ..ctrl_long
        for(i in qui){
            fml_dp = gsub(paste0(escape_regex(names(macros)[i]), "(?=$|[^[:alnum:]_\\.])"), macros[[i]], fml_dp, perl = TRUE)
        }
        fml = as.formula(fml_dp)
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
#' @param multi.join Logical, or character, scalar, defaults to \code{FALSE}. Only used if multiple vectors are to be transformed into integers. If \code{multi.join} is not \code{FALSE}, then the values of the different vectors will be collated using \code{\link[base]{paste}} with \code{collapse=multi.join}.
#'
#' @details
#' If multiple vectors have to be combined and \code{add_items=TRUE}, to have user readable values in the items, you should add the argument \code{multi.join} so that the values of the vectors are combined in a "user-readable" way. Note that in the latter case, the algorithm is much much slower.
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
#' # by default, the two vector are fast combined, and items are meaningless
#' to_integer(x1, x2, add_items = TRUE)
#'
#' # You can use multi.join to have human-readable values for the items:
#' to_integer(x1, x2, add_items = TRUE, multi.join = TRUE)
#'
#' to_integer(x1, x2, add_items = TRUE, multi.join = "; ")
#'
to_integer = function(..., sorted = FALSE, add_items = FALSE, items.list = FALSE, multi.join = FALSE){

    check_arg(..., "vector mbt")
    check_arg(sorted, add_items, items.list, "logical scalar")
    check_arg(multi.join, "scalar(logical, character)")

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

    do_join = FALSE
    if(!isFALSE(multi.join) && Q > 1){
        do_join = TRUE
        if(isTRUE(multi.join)) multi.join = "_"
    }

    if(Q == 1){
        if(sorted && is.factor(dots[[1]])){
            # Special treatment for factors => we keep their order
            f = dots[[1]][drop = TRUE]
            res = quickUnclassFactor(unclass(f), addItem = add_items, sorted = sorted)
            if(add_items){
                res$items = levels(f)[res$items]
            }

        } else {
            res = quickUnclassFactor(dots[[1]], addItem = add_items, sorted = sorted)
        }

    } else {

        if(do_join){
            myDots = dots
            myDots$sep = multi.join
            index = do.call("paste", myDots)

        } else {

            for(q in 1:Q){
                dots[[q]] = quickUnclassFactor(dots[[q]], sorted = sorted)
            }

            # Then we combine
            power = floor(1 + log10(sapply(dots, max)))

            if(sum(power) > 14){
                myDots = dots
                myDots$sep = multi.join
                index = do.call("paste", myDots)
            } else {
                # quicker, but limited by the precision of doubles
                index = dots[[1]]
                for(q in 2:Q){
                    index = index + dots[[q]]*10**sum(power[1:(q-1)])
                }
            }
        }

        res = quickUnclassFactor(index, addItem = add_items, sorted = sorted)
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
#' @param X A matrix, vector or a list. The vectors to be centered. There must be the same number of observations as in the factors used for centering (argument \code{f}).
#' @param f A matrix, vector or list. The factors used to center the variables in argument \code{X}.
#' @param weights Vector, can be missing or NULL. If present, it must contain the same number of observations as in \code{X}.
#' @param nthreads Number of threads to be used. By default it is equal to \code{getFixest_nthreads()}.
#' @param notes Logical, whether to display a message when NA values are removed. By default it is equal to \code{getFixest_notes()}.
#' @param iter Number of iterations, default is 2000.
#' @param tol Stopping criterion of the algorithm. Default is \code{1e-6}. The algorithm stops when the maximum absolute increase in the coefficients values is lower than \code{tol}.
#' @param im_confident Logical, default is \code{FALSE}. FOR EXPERT USERS ONLY! This argument allows to skip some of the preprocessing of the arguments given in input. If \code{TRUE}, then \code{X} MUST be a numeric matrix (not integer, numeric), \code{f} MUST be a list and \code{weights}, if given, MUST be numeric (not integer!). Further the three MUST NOT contain any NA values and MUST have the same number of observations. Non compliance to these rules may simply lead your R session to break.
#'
#' @return
#' It returns a matrix of the same number of columns as \code{X} in input. The number of rows is equal to the number of rows of \code{X} minus the number of NA values (contained in \code{X}, \code{f} or \code{weights}).
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
#'
#'
#'
demean = function(X, f, weights, nthreads = getFixest_nthreads(), notes = getFixest_notes(), iter = 2000, tol = 1e-6, im_confident = FALSE){
    # LATER: add slopes

    # Step 1: formatting the input
    if(!im_confident){
        check_arg(X, "numeric vmatrix | list mbt")
        check_arg(f, "vmatrix | list mbt")
        check_arg(iter, "integer scalar GE{1}")
        check_arg(tol, "numeric scalar GT{0}")
        check_arg(notes, "logical scalar")

        ## X
        if(is.list(X)){
            var_names = names(X)

            if(!is.data.frame(X)){
                n_all = lengths(X)
                if(any(diff(n_all) != 0)){
                    n_unik = unique(n_all)
                    stop("In argument 'X' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                }
                X = as.data.frame(X)
            }
            is_numeric = sapply(X, is.numeric)
            if(any(!is_numeric)){
                n_non_num = sum(is_numeric)
                stop("All values in 'X' must be numeric. This is not the case for ", n_non_num, " variable", plural(n_non_num), ".")
            }
            X = as.matrix(X)
        } else {
            if(!is.matrix(X)) X = as.matrix(X)
            var_names = colnames(X)
        }

        if(is.null(var_names)) var_names = paste0("V", 1:ncol(X))
        if(is.integer(X)) X = X * 1

        ## f
        if(is.list(f)){
            fe_names = names(f)

            if(!is.data.frame(f)){
                n_all = lengths(f)
                if(any(diff(n_all) != 0)){
                    n_unik = unique(n_all)
                    stop("In argument 'f' all elements of the list must be of the same length. This is currently not the case (ex: one is of length ", n_unik[1], " while another is of length ", n_unik[2], ").")
                }
                f = as.data.frame(f)
            }
        } else {
            fe_names = colnames(f)
            f = as.data.frame(f)
        }

        if(is.null(fe_names)) fe_names = paste0("Factor_", 1:length(f))

        is_pblm = sapply(f, function(x) !is.numeric(x) || !is.character(x))
        for(i in which(is_pblm)){
            f[[i]] = as.character(f[[i]])
        }

        if(nrow(f) != nrow(X)) stop("The number of observations in 'X' and in 'f' don't match (", nrow(X), " vs ", nrow(f), ").")


        ## weights
        check_arg_plus(weights, "NULL numeric conv vector len(data) GE{0}", .data = X)
        is_weight = !missing(weights) && !is.null(weights)

        ## nthreads
        nthreads = check_set_nthreads(nthreads)

        #
        # FORMATTING
        #

        # NAs
        is_NA = FALSE
        info_X = cpppar_which_na_inf_mat(X, nthreads)
        n_na_x = 0
        if(info_X$any_na_inf){
            is_NA = is_NA | info_X$is_na_inf
            n_na_x = sum(info_X$is_na_inf)
        }

        n_na_fe = 0
        if(anyNA(f)){
            is_na_fe = rowSums(is.na(f)) > 0
            n_na_fe = sum(is_na_fe)
            is_NA = is_NA | is_na_fe
        }

        weight_msg = ""
        if(is_weight){
            is_na_weights = is.na(weights)
            weight_msg = paste0(", weights: ", sum(is_na_weights))
            is_NA = is_NA | is_na_weights
        }

        n_na = sum(is_NA)
        if(n_na > 0 && notes){
            if(all(is_NA)) stop("All observations contain NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, weight_msg, ").")

            message("NOTE: ", signif_plus(n_na), " observation", plural(n_na), " removed because of NA values (Breakup: X: ", n_na_x, ", f: ", n_na_fe, weight_msg, ").")
        }

        if(n_na > 0){
            X = X[!is_NA, , drop = FALSE]
            f = f[!is_NA, , drop = FALSE]
            if(is_weight) weights = weights[!is_NA]
        }

        if(!is_weight) weights = 1

    } else {
        if(missing(weights) || is.null(weights)) weights = 1
    }

    #
    # Unclassing fes
    #

    Q = length(f)
    n = length(f[[1]])
    quf_info_all = cpppar_quf_table_sum(x = f, y = 0, do_sum_y = FALSE, rm_0 = FALSE, rm_1 = FALSE, rm_single = FALSE, do_refactor = FALSE, r_x_sizes = 0, obs2keep = 0, only_slope = rep(FALSE, Q), nthreads = nthreads)

    # table/sum_y/sizes
    fixef_table = quf_info_all$table
    fixef_sizes = lengths(fixef_table)
    fixef_table_vector = as.integer(unlist(fixef_table))

    #
    # The algorithm
    #

    vars_demean <- cpp_demean(y = 0, X, ncol(X), weights, iterMax = iter,
                              diffMax = tol, r_nb_id_Q = fixef_sizes,
                              fe_id_list = quf_info_all$quf, table_id_I = fixef_table_vector,
                              slope_flag_Q = rep(FALSE, Q), slope_vars_list = list(0),
                              r_init = 0, nthreads = nthreads)

    X_demean = vars_demean$X_demean
    colnames(X_demean) = var_names
    X_demean
}



#' Computes fit statistics of fixest objects
#'
#' Computes various fit statistics for \code{fixest} estimations.
#'
#' @param x A \code{fixest} estimation.
#' @param type Character vector. The type of fit statistic to be computed. So far, \code{"G"}, the 'working' number of observations, and \code{'theta'}, the over-dispersion parameter in negative binomial estimation, are available.
#' @param as.list Logical, default is \code{FALSE}. Only used when one statistic is to be computed (i.e. arg. \code{type} is of length 1). If \code{TRUE}, then a list whose name is the \code{type} is returned containing the statistic.
#' @param ... Other elements to be passed to other methods and may be used to compute the statistics (for example you can pass on arguments to compute the VCOV when using \code{type = "G"}).
#'
#' @return
#' If two or more types are to be computed, then a list is returned whose names correstpond to the argument \code{type}. If only one type is to be computed, then by default the result is simplified to a vector (or list, it depends on the type) containing the statistics; to return a list instead, \code{as.list=TRUE} needs to be used.
#'
#' @examples
#'
#' data(trade)
#' gravity = feols(log(Euros) ~ log(dist_km) | Destination + Origin, trade)
#'
#' # Extracting the 'working' number of observations used to compute the pvalues
#' fitstat(gravity, "G")
#'
#' # Idem, but when coimputing two-way SEs
#' fitstat(gravity, "G", se = "standard")
#'
#'
fitstat = function(x, type, as.list = FALSE, ...){

    # Later:
    # f => f stat
    #   * f.df / f.pval / f.stat
    # wf => within f stat
    #   * wf.df / wf.pval / wf.stat

    dots = list(...)
    valid_types = c("G", "theta")
    if(isTRUE(dots$give_types)){

        tex_alias = c("g" = "Size of 'working' sample", "theta" = "Over-dispersion")
        R_alias   = c("g" = "G", "theta" = "Over-dispersion")

        # type_alias = c("f" = "f.stat")

        res = list(types = tolower(valid_types), tex_alias = tex_alias, R_alias = R_alias)
        return(res)
    }

    check_arg(x, "class(fixest) mbt")
    check_arg_plus(type, "mbt multi match", .choices = valid_types)
    check_arg(as.list, "logical scalar")

    res_all = list()
    type_all = type

    for(i in seq_along(type_all)){
        type = type_all[i]

        if(type == "G"){
            fromSummary = "cov.scaled" %in% names(x)
            if(!fromSummary){
                my_vcov = vcov(x, attr = TRUE, ...)
            } else {
                my_vcov = x$cov.scaled
            }

            G = attr(my_vcov, "G")
            if(is.null(G)) G = x$nobs

            res_all[[type]] = G

        } else if(type == "theta"){
            isNegbin = x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && x$family=="negbin")
            if(isNegbin){
                theta = coef(x)[".theta"]
                names(theta) = "Overdispersion"
            } else {
                theta = NA
            }

            res_all[[type]] = theta
        }

    }

    if(length(type_all) == 1 && !as.list){
        res_all = res_all[[1]]
    }

    res_all
}

#### ................. ####
#### Internal Funs     ####
####


parse_macros = function(..., reset = FALSE, from_xpd = FALSE){
    set_up(1)

    check_arg(..., "dotnames os formula | character vector no na | numeric scalar | class(call, name)", .message = paste0("Each element of '...' must be a one-sided formula, and the name of each argument must start with two dots (ex: ", ifelse(from_xpd, "xpd(fml, ..ctrl = ~ x5 + x6)", "setFixest_fml(..ctrl = ~ x5 + x6)"), ").\nAlternatively it can be a character vector of variable names, or a numeric scalar."))

    # We require os formulas instead of character strings because:
    # 1) I find it more handy
    # 2) there is a parsing check from R

    # Original macros
    fml_macro = getOption("fixest_fml_macro")
    if(is.null(fml_macro)){
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
    qui_pblm = !grepl("^\\.\\.", names(dots))
    if(any(qui_pblm)){
        arg_name = names(dots)[qui_pblm][1]
        correct_name = gsub("^\\.*", "\\.\\.", arg_name)
        stop_up("Each argument name must start with two dots. Use '", correct_name, "' instead of '", arg_name, "'.")
    }

    for(v in names(dots)){
        fml_raw = dots[[v]]

        if(any(c("call", "name") %in% class(fml_raw))){
            fml_raw = deparse_long(fml_raw)
        }

        if(!"formula" %in% class(fml_raw)){
            fml_raw = fml_raw[grepl("[[:alnum:]]", fml_raw)]

            if(length(fml_raw) > 0){
                fml_str = paste("~", paste(fml_raw, collapse = " + "))
                fml = try(as.formula(fml_str))
                if("try-error" %in% class(fml)){
                    stop("The variable '", v, "' does not lead to a valid formula (i.e. ", fml_str, " is not valid).")
                }
            } else {
                fml = ~1
            }

        } else {
            fml = fml_raw
        }

        fml_macro[[v]] = deparse_long(fml[[2]])
    }

    fml_macro
}

# style_name = "fixef"
# keywords = c("title", "prefix", "suffix")
parse_style = function(x, keywords){
    # x is a character scalar, example x = "title:Variables:;below:_"
    # style_name: name of the style, only to format the error message
    # keywords: the vector of accepted keywords

    style_name = gsub(".+\\$", "", deparse(substitute(x)))

    set_up(1)

    res = setNames(as.list(rep("", length(keywords))), keywords)
    if(grepl("^ *$", x)) return(res)

    x_split = strsplit(x, ";")[[1]]
    if(!all(grepl(":", x_split))){
        stop_up("In argument 'style', the styles must be of the form 'keyword1:value1; keyword2:value2' etc. A colon is currently missing in the '", style_name, "' style (i.e. '", x, "' is not valid).")
    }

    kws = gsub("^ +|:.*", "", x_split)

    check_value(kws, "multi charin", .choices = keywords, .message = paste0("In argument 'style', the keywords of '", style_name, "' must be equal to ", enumerate_items(keywords, "quote.or"), "."))

    values = as.list(gsub("^[^:]+:", "", x_split))

    for(i in seq_along(kws)){
        res[[kws[i]]] = escape_latex(values[i], up = 3, TRUE)
    }

    res
}

myPrintCoefTable = function(coeftable, lastLine = "", show_signif = TRUE){
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

    for(i in 1:4){
        ct[, i] = decimalFormat(ct[, i])
    }

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


vcovClust <- function (cluster, myBread, scores, dof=FALSE, do.unclass=TRUE){
    # Internal function: no need for controls, they come beforehand
    # - cluster: the vector of dummies
    # - myBread: original vcov
    # - scores
    # Note: if length(unique(cluster)) == n (i.e. White correction), then the dof are such that vcovClust is equivalent to vcovHC(res, type="HC1")
    # Source: http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
    #         Cameron & Miller -- A Practitioner’s Guide to Cluster-Robust Inference

    n <- NROW(scores)

    # Control for cluster type
    if(do.unclass){
        cluster <- quickUnclassFactor(cluster)
    }

    Q <- max(cluster)
    RightScores = cpp_tapply_sum(Q, scores, cluster)

    # Finite sample correction:
    if(dof){
        dof_value  <- Q / (Q - 1)
    } else {
        dof_value = 1
    }

    # return(crossprod(RightScores%*%myBread) * dof_value)

    # Maybe I'll later add parallelism, but I don't want to add other arguments to the functions
    # calling vcov... and there are many
    # I don't want either to apply "hidden" parallelism, even using getOption
    # We'll see....
    xy = cpppar_matprod(RightScores, myBread, 1)
    res = cpppar_crossprod(xy, 1, 1) * dof_value
    res
}

prepare_matrix = function(fml, base, fake_intercept = FALSE){
    # This function is way faster than model.matrix but does not accept factors
    # The argument fml **MUST** not have factors!

    rhs = fml[c(1,3)]
    t = terms(rhs, data = base)

    all_var_names = attr(t, "term.labels")

    # We take care of interactions: references can be multiple, then ':' is legal
    all_vars = gsub(":", "*", all_var_names)

    if(any(qui_inter <- (grepl("^i(nteract)?\\(", all_var_names) & grepl(":", all_var_names, fixed = TRUE)))){
        # beware of ":" in drop/keep!!!

        for(arg in c("drop", "keep")){
            if(any(qui_ref <- grepl(paste0(arg, " =[^\\)]+:"), all_var_names[qui_inter]))){
                var_inter_ref = all_var_names[qui_inter][qui_ref]
                var_inter_ref_split = strsplit(var_inter_ref, paste0(arg, " = "))
                fun2apply = function(x) paste(gsub(":", "*", x[1]), x[2], sep = paste0(arg, " = "))
                new_var_inter_ref = sapply(var_inter_ref_split, fun2apply)
                all_vars[qui_inter][qui_ref] = new_var_inter_ref
            }
        }

    }

    # Forming the call
    if(attr(t, "intercept") == 1 && !fake_intercept){
        n = nrow(base)
        all_vars_call = parse(text = paste0("list('(Intercept)' = rep(1, ", n, "), ", paste0(all_vars, collapse = ", "), ")"))
        all_var_names = c("(Intercept)", all_var_names)
    } else {
        all_vars_call = parse(text = paste0("list(", paste0(all_vars, collapse = ", "), ")"))
    }

    # evaluation
    data_list <- eval(all_vars_call, base)

    # Handling the multi columns case (ex: bs(x1), splines)
    if(any(lengths(data_list) != nrow(base))){

        all_n = as.vector(lengths(data_list) / nrow(base))

        qui_pblm = which(all_n %% 1 != 0)
        if(length(qui_pblm) > 0){
            what = data_list[[qui_pblm]]
            reason = ifelse(is.null(nrow(what)), paste0("of length ", length(what)), paste0("with ", nrow(what), " rows"))

            stop("Evaluation of ", all_var_names[qui_pblm], " returns an object ", reason, " while the data set has ", nrow(base)," rows.", call. = FALSE)
        }

        all_n_vector = rep(all_n, all_n)

        new_names = as.list(all_var_names)
        for(i in which(all_n > 1)){
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


fixest_model_matrix = function(fml, data, fake_intercept = FALSE, i_noref = FALSE){
    # This functions takes in the formula of the linear part and the
    # data
    # It reformulates the formula (ie with lags and interactions)
    # then either apply a model.matrix
    # either applies an evaluation (which can be faster)
    #
    # fake_intercept => whether to add the intercept, only to make sure
    #  the factors are well created

    # fml = ~a*b+c+i(x1)+Temp:i(x2)+i(x3)/Wind

    # Modify the formula to add interactions
    if(grepl("::", deparse_long(fml[[3]]), fixed = TRUE)){
        fml = expand_interactions(fml)
    }

    #
    # Evaluation
    #

    t_fml = terms(fml)
    tl = attr(t_fml, "term.labels")

    if(length(tl) == 0){

        if(fake_intercept){
            return("NOT_LINEAR")
        }

        res = matrix(1, nrow = nrow(data), ncol = 1, dimnames = list(NULL, "(Intercept)"))
        return(res)
    }

    # We check for calls to i()
    qui_inter <- grepl("(^|[^[:alnum:]_\\.])i(nteract)?\\(", tl)
    IS_INTER = any(qui_inter)
    if(IS_INTER){
        # OMG... why do I always have to reinvent the wheel???
        is_intercept = fake_intercept || (attr(t_fml,"intercept") == 1)
        i_naked = which(is_naked_fun(tl[qui_inter], "i(nteract)?"))

        if(i_noref){
            for(i in seq_along(i_naked)){
                j = i_naked[i]
                txt = gsub("(^|(?<=[^[:alnum:]\\._]))i(nteract)?\\(", "i_noref(", tl[qui_inter][j], perl = TRUE)
                tl[qui_inter][j] = eval(parse(text = txt))
            }
        } else {
            for(i in seq_along(i_naked)){
                if(!is_intercept && i == 1) next

                j = i_naked[i]
                txt = gsub("(^|(?<=[^[:alnum:]\\._]))i(nteract)?\\(", "i_ref(", tl[qui_inter][j], perl = TRUE)
                tl[qui_inter][j] = eval(parse(text = txt))
            }
        }


        fml_no_inter = as.formula(paste0("y ~ ", paste(c(1, tl[!qui_inter]), collapse = "+")))
        fml = as.formula(paste0("y ~ ", paste(tl, collapse = "+")))

    }

    # Are there factors NOT in i()? If so => model.matrix is used
    dataNames = names(data)

    if(IS_INTER){
        linear.varnames = all.vars(fml_no_inter[[3]])
        is_num = sapply(data[, dataNames %in% linear.varnames, FALSE], is.numeric)
        if(length(is_num) > 0 && (any(!is_num) || grepl("factor", deparse_long(fml_no_inter)))){
            useModel.matrix = TRUE
        } else {
            useModel.matrix = FALSE
        }
    } else {
        linear.varnames = all.vars(fml[[3]])
        is_num = sapply(data[, dataNames %in% linear.varnames, FALSE], is.numeric)
        if(length(is_num) == 0 || any(!is_num) || grepl("factor", deparse_long(fml))){
            useModel.matrix = TRUE
        } else {
            useModel.matrix = FALSE
        }
    }

    if(useModel.matrix){
        # to catch the NAs, model.frame needs to be used....
        linear.mat = stats::model.matrix(fml, stats::model.frame(fml, data, na.action=na.pass))
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

fixef_terms = function(fml){
    # separate all terms of fml into fixed effects and varying slopes
    # fml: one sided formula
    # can also be a vector of terms

    # fixef_terms(~dum_1[[x1]] + dum_2[x2])
    # fixef_terms(~dum_1[[x1]] + dum_2 + dum_3[x2, x3])

    if(is.vector(fml)){
        fml = as.formula(paste0("~", paste0(fml, collapse = "+")))
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

    # we need to take care of the ^ used to combine variables
    fml_char = as.character(fml[2])
    if(grepl("^", fml_char, fixed = TRUE)){
        fml_char_new = gsub("^", "_impossible_var_name_", fml_char, fixed = TRUE)
        fml = as.formula(paste0("~", fml_char_new))
    }


    t = try(terms(fml))

    # if error, we send it back
    if("try-error" %in% class(t)){
        t = gsub("_impossible_var_name_", "^", t, fixed = TRUE)
        return(t)
    }

    my_vars = attr(t, "term.labels")

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

                new_v = eval(parse(text = v_mod))
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
    my_vars = gsub("_impossible_var_name_", "^", my_vars, fixed = TRUE)

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

        fun2combine = ifelse(fastCombine, "combine_clusters_fast", "combine_clusters")

        vars_new = gsub("([[:alpha:]\\.][[:alnum:]_\\.]*(\\^[[:alpha:]\\.][[:alnum:]_\\.]*)+)",
                        paste0(fun2combine, "(\\1)"), vars)

        vars_new = gsub("\\^", ", ", vars_new)

        vars = vars_new
        changeNames = TRUE
    }

    all_vars = gsub(":", "*", vars)

    if(all(all_vars %in% names(base))){
        res = base[, all_vars, drop = FALSE]
    } else {
        all_vars_call = parse(text = paste0("list(", paste0(all_vars, collapse = ", "), ")"))
        data_list <- try(eval(all_vars_call, base))

        # if error: we send it back to the main function
        if("try-error" %in% class(data_list)){
            return(data_list)
        }

        names(data_list) = all_var_names
        data_list$stringsAsFactors = FALSE

        res = do.call("data.frame", data_list)

        if(changeNames){
            qui = grepl("combine_clusters", all_var_names, fixed = TRUE)
            new_names = gsub("combine_clusters(_fast)?\\(|\\)", "", all_var_names[qui])
            new_names = gsub(", ?", "^", new_names)
            all_var_names[qui] = new_names
        }

        names(res) = all_var_names
    }

    res
}

prepare_cluster_mat = function(fml, base, fastCombine){
    # prepares the data.frame of the cluster variables

    fml_char = as.character(fml[2])
    changeNames = FALSE
    if(grepl("^", fml_char, fixed = TRUE)){
        # special indicator to combine factors

        fun2combine = ifelse(fastCombine, "combine_clusters_fast", "combine_clusters")

        fml_char_new = gsub("([[:alpha:]\\.][[:alnum:]_\\.]*(\\^[[:alpha:]\\.][[:alnum:]_\\.]*)+)",
                            paste0(fun2combine, "(\\1)"),
                            fml_char)

        fml_char_new = gsub("([[:alpha:]\\.][[:alnum:]_\\.]*)\\^([[:alpha:]\\.][[:alnum:]_\\.]*)", "\\1, \\2", fml_char_new)
        fml = as.formula(paste0("~", fml_char_new))
        changeNames = TRUE
    }

    t = terms(fml, data = base)

    all_var_names = attr(t, "term.labels")
    all_vars = gsub(":", "*", all_var_names)

    if(all(all_vars %in% names(base))){
        res = base[, all_vars, drop = FALSE]
    } else {
        all_vars_call = parse(text = paste0("list(", paste0(all_vars, collapse = ", "), ")"))
        data_list <- eval(all_vars_call, base)
        names(data_list) = all_var_names
        data_list$stringsAsFactors = FALSE

        res = do.call("data.frame", data_list)

        if(changeNames){
            qui = grepl("combine_clusters", all_var_names)
            new_names = gsub("combine_clusters(_fast)?\\(|\\)", "", all_var_names[qui])
            new_names = gsub(", ?", "^", new_names)
            all_var_names[qui] = new_names
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
        myDots = cluster
        myDots$sep = "_"
        index = do.call("paste", myDots)
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

    cluster = list(...)
    Q = length(cluster)

    # See comments in combine_cluster_fast
    ANY_NA = FALSE
    if(any(who_NA <- sapply(cluster, anyNA))){
        ANY_NA = TRUE
        who_NA = which(who_NA)
        IS_NA = is.na(cluster[[who_NA[1]]])
        for(i in who_NA[-1]){
            IS_NA = IS_NA | is.na(cluster[[i]])
        }

        if(all(IS_NA)) return(rep(NA, length(IS_NA)))

        # we recreate the clusters
        for(i in 1:Q){
            cluster[[i]] = cluster[[i]][!IS_NA]
        }
    }

    # We just paste
    myDots = cluster
    myDots$sep = "_"
    index = do.call("paste", myDots)

    if(ANY_NA){
        # we recreate the return vector with appropriate NAs
        res = rep(NA_character_, length(IS_NA))
        res[!IS_NA] = index
    } else {
        res = index
    }

    return(res)
}

expand_interactions_internal = function(x){
    # x == terms

    terms_all_list = as.list(x)
    qui = which(grepl("[^:]::[^:]", x))
    for(i in qui){
        my_term = x[i]
        terms_split = strsplit(x[i], "(?<=[^:])::(?=[^:])", perl = TRUE)[[1]]

        if(grepl("\\(", terms_split[2])){
            if(!grepl("^[\\.]?[[:alnum:]\\._]+\\(", terms_split[2])){
                stop("Problem in ", x[i], ": the format should be var::fe. See details.")
            }

            my_call = gsub("^[\\.]?[[:alnum:]\\._]+\\(", "interact_control(", terms_split[2])
            args = try(eval(parse(text = my_call)), silent = TRUE)
            fe_name = gsub("^([\\.]?[[:alnum:]\\._]+)\\(.+", "\\1", terms_split[2])
            if("try-error" %in% class(args)){
                stop("Problem in the interaction of the formula: Error in ", terms_split[1], "::", fe_name, gsub(".+interact_control", "", args))
            }

            new_term = paste0("interact(", terms_split[1], ", ", fe_name, ", ", paste(args, collapse = ", "), ")")

        } else {
            new_term = paste0("interact(", terms_split[1], ", ", terms_split[2], ")")
        }

        terms_all_list[[i]] = new_term

    }

    x = unlist(terms_all_list)

    return(x)
}

expand_interactions = function(fml){
    # The formula is simple (should contain only the RHS)

    fml_linear = fml_split(fml, 1)

    x = attr(terms(fml_linear), "term.labels")

    if (!any(grepl("[^:]::[^:]", x))){
        return(fml)

    } else {
        x = expand_interactions_internal(x)
    }

    lhs_fml = deparse_long(fml_linear[[2]])
    rhs_fml = paste(x, collapse = "+")

    as.formula(paste0(lhs_fml, "~", rhs_fml))
}

interact_control = function(ref, drop, keep){
    # Internal call
    # used to control the call to interact is valid

    counter = getOption("fixest_deprec_interact")
    if(is.null(counter)){
        options("fixest_deprec_interact" = TRUE)
        .Deprecated(msg = "Interactions with the syntax var::fe is deprecated and will disappear in 1 year from 12/11/2020. Please use the function i() instead.")
    }

    mc = match.call()

    res = c()
    if("ref" %in% names(mc)){
        res = paste0("ref = ", deparse_long(mc$ref))
    }

    if("drop" %in% names(mc)){
        res = paste0("drop = ", deparse_long(mc$drop))
    }

    if("keep" %in% names(mc)){
        res = paste0("keep = ", deparse_long(mc$keep))
    }

    res
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

shade_area <- function(y1, y2, x, xmin, xmax, col="grey", ...){
    # fonction plus pratique que polygon
    # elle permet de griser une partie delimitee par
    # y1 et y2 pour chacune des valeurs de x
    # on doit avoir la meme longueur de y1,y2 et x
    # exemple:
    # a=curve(x**2,-5,5)
    # shade_area(a$y+1,a$y-1,a$x)
    # qqes parametres graphiques:
    # lwd / border (couleur du bord, peut etre NA) / lty

    n <- length(x)
    stopifnot(length(y1)==n | length(y1)==1)
    stopifnot(length(y2)==n | length(y2)==1)

    if(length(y1)==1) y1 <- rep(y1,n)
    if(length(y2)==1) y2 <- rep(y2,n)

    if(missing(xmin)) xmin <- min(x)
    if(missing(xmax)) xmax <- max(x)

    ind <- which(x>=xmin & x<=xmax)
    x1 <- x[ind] ; x2 <- x[rev(ind)]
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

    x_split = strsplit(x2clean, "(^|(?<=[^[:alnum:]\\._]))i(nteract)?\\(|__CLEAN__", perl = TRUE)

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
    # x = c("i(x1)", "i(I(x3))", "interact(x3, x4, TRUE, drop = c(1, 3:5))", "Temp:i(x2)", "i(x3):Wind")

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



#### ................. ####
#### Small Utilities ####
####

escape_all = function(x){
    # we escape all
    res = gsub("((?<=[^\\\\])|(?<=^))(\\$|_|%|&|\\^|#)", "\\\\\\2", x, perl = TRUE)
    res
}

escape_latex = function(x_all, up = 0, noArg = FALSE){
    # This is super tricky to escape properly!
    # We do NOT escape within equations

    x_name = deparse(substitute(x_all))

    res = c()

    for(index in seq_along(x_all)){
        x = x_all[index]

        # 1) finding out equations, ie non escaped dollar signs
        dollars = gregexpr("((?<=[^\\\\])|(?<=^))\\$", x, perl = TRUE)[[1]]

        is_eq = FALSE
        if(length(dollars) > 1){
            is_eq = TRUE
            if(length(dollars) %% 2 != 0){
                my_arg = "T"
                if(!noArg){
                    my_arg = paste0("In argument '", x_name, "', t")
                }
                stop_up(up = up, my_arg, "here are ", length(dollars), " dollar signs in the following character string:\n", x, "\nIt will raise a Latex error (which '$' means equation? which means dollar-sign?): if you want to use a regular dollar sign, please escape it like that: \\\\$.")
            }
        }

        # 2) Escaping but conditionnally on not being in an equation
        if(is_eq){
            # Finding out the equations
            all_items = strsplit(paste0(x, " "), "((?<=[^\\\\])|(?<=^))\\$", perl = TRUE)[[1]]
            for(i in seq_along(all_items)){
                if(i %% 2 == 1){
                    all_items[i] = escape_all(all_items[i])
                }
            }

            res[index] = gsub(" $", "", paste(all_items, collapse = "$"))
        } else {
            res[index] = escape_all(x)
        }
    }

    if(!is.null(names(x_all))){
        names(res) = names(x_all)
    }

    res
}

escape_regex = function(x){
    # escape special characters in regular expressions => to make it as "fixed"

    res = gsub("((?<=[^\\\\])|(?<=^))(\\$|\\.|\\+|\\(|\\)|\\[|\\]|\\?|\\^)", "\\\\\\2", x, perl = TRUE)
    res
}

formatBicLL = function(bic, ll){
	# entry: bic and ll

	bic = numberFormatNormal(bic)
	ll = numberFormatNormal(ll)

	bic_split = strsplit(bic, "\\.")[[1]]
	ll_split = strsplit(ll, "\\.")[[1]]

	n = max(nchar(bic_split[1]), nchar(ll_split[1]))

	bic_new = sprintf("% *s.%s", n, bic_split[1], ifelse(length(bic_split) == 2, bic_split[2], 0))
	ll_new = sprintf("% *s.%s", n, ll_split[1], ifelse(length(ll_split) == 2, ll_split[2], 0))

	sep = "  "
	myWidth = max(nchar(c(bic_new, ll_new))) + length(sep) + 1

	bic_format = paste0(bic_new, sprintf("% *s", myWidth - nchar(bic_new), sep))
	ll_format = paste0(ll_new, sprintf("% *s", myWidth - nchar(ll_new), sep))

	list(bic = bic_format, ll = ll_format)
}

addCommas_single = function(x){

	if (!is.finite(x)) return(as.character(x))

	s = sign(x)
	x = abs(x)
	decimal = x - floor(x)
	if (decimal > 0){
		dec_string = substr(decimal, 2, 4)
	} else {
		dec_string = ""
	}

	entier = sprintf("%.0f", floor(x))
	quoi = rev(strsplit(entier, "")[[1]])
	n = length(quoi)
	sol = c()
	for (i in 1:n) {
		sol = c(sol, quoi[i])
		if (i%%3 == 0 && i != n) sol = c(sol, ",")
	}
	res = paste0(ifelse(s == -1, "-", ""), paste0(rev(sol), collapse = ""),
					 dec_string)
	res
}

addCommas = function(x){
	sapply(x, addCommas_single)
}

myRound_single = function(x, digits=5){
	# There can be non numeric values...
	# we give away the non numeric ones and round the others

	if(is.na(x)){
		return(NA)
	}

	if(is.numeric(x)){
		res = round(x, digits)
	} else {

		if(!grepl("[[:digit:]]", x)){
			# means it is a character
			res = x
		} else {
			res = round(as.numeric(x), digits)
		}
	}

	res
}

myRound = function(x, digits=5){
	sapply(x, myRound_single, digits = digits)
}

decimalFormat = function(x){

	decimalFormat_single = function(x){
		# for very small numbers: format 5.2e-08

		if(is.na(x) || !is.numeric(x)) return(x)

		xPower = log10(abs(x))

		if(xPower < -5){
			res = signif(x, 3)
		} else if(xPower < 0){
			res = round(x, 6)
		} else {
			res = round(x, max(1, 5 - ceiling(xPower)))
		}

		res
	}

	sapply(x, decimalFormat_single)
}

coefFormatLatex = function(x, digits = 4, power = 5){

	coefFormatLatex_single = function(x, digits, power){
		# format decimals: 5.3 10**-7 instead of 0.00000053
		# format large numbers 6356516.12464 => 6356516.1

		nbSignif = 3

		if(is.na(x)) return(x)

		if(!is.numeric(x)){
			if(grepl("[^[:digit:]e\\.-]", x)){
				return(x)
			} else {
				x = as.numeric(x)
			}
		}

		exponent = floor(log10(abs(x)))

		if(exponent >= 0){
			# return(sprintf("%.*f", max(1, digits - abs(exponent)), x))
			return(mysignif(x, digits))
		}

		if(abs(exponent) >= power){
			left_value = round(x*10**-exponent, 3)
			res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
		} else if(abs(x) > 10**(-digits)){
			res = sprintf("%.*f", digits, x)
		} else {
			res = sprintf("%.*f", abs(exponent), x)
		}

		res
	}

	sapply(x, coefFormatLatex_single, digits = digits, power = power)
}

numberFormat_single = function(x, type = "normal"){
	# For numbers higher than 1e9 => we apply a specific formatting
	# idem for numbers lower than 1e-4

    if(is.na(x)){
        return(NA)
        # return(ifelse(type == "normal", "--", ""))
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

    mysignif_single = function(x, d, r) {
        if (is.na(x)) {
            return(NA)
        }

        if (abs(x) >= 10^(d - 1)){
            return(round(x, r))
        } else {
            return(signif(x, d))
        }
    }
    sapply(x, mysignif_single, d = d, r = r)
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

    who = x %in% names(dict)
    x[who] = dict[as.character(x[who])]
    x
}

keep_apply = function(x, keep = NULL){

    if(missing(keep) || length(keep) == 0){
        return(x)
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

	if(!is.numeric(x)){
		# level and unclass is much slower
		x = as.character(x)
	}

    res = cpp_quf_gnl(x)

    if(sorted){

        if(is.character(x)){
            items = x[res$x_unik]
        } else {
            items = res$x_unik
        }

        x = res$x_uf

        new_order = order(items)
        order_new_order = order(new_order)
        x_uf = order_new_order[x]

        if(addItem){
            res = list(x = x_uf, items = items[new_order])
            return(res)
        } else {
            return(x_uf)
        }
    }

    if(addItem){

        if(is.character(x)){
            items = x[res$x_unik]
            res = list(x = res$x_uf, items = items)
        } else {
            names(res) = c("x", "items")
        }

        return(res)
    } else {
        return(res$x_uf)
    }
}

missnull = function(x){
	if(missing(x) || is.null(x)){
		return(TRUE)
	} else {
		return(FALSE)
	}
}

isScalar = function(x, int = FALSE){
	if(length(x) == 1 && is.numeric(x) && is.finite(x)){
	    if(int && !(x %% 1 == 0)) return(FALSE)
		return(TRUE)
	} else {
		return(FALSE)
	}
}

isLogical = function(x){
	if(length(x) == 1 && is.logical(x) && !is.na(x)){
		return(TRUE)
	} else {
		return(FALSE)
	}
}

isSingleChar = function(x){
	if(length(x) == 1 && is.character(x) && !is.na(x)){
		return(TRUE)
	} else {
		return(FALSE)
	}
}

fml2varnames = function(fml){
	# This function trandforms a one sided formula into a
	# character vector for each variable

	fml_char = gsub(" +", "", as.character(fml)[2])

	# In theory, I could just use terms.formula to extract the variable.
	# But I can't!!!! Because of this damn ^ argument.
	# I need to apply a trick
	if(grepl("[[:alpha:]\\.][[:alnum:]\\._]*\\^", fml_char)){
		fml_char_new = gsub("([[:alpha:]\\.][[:alnum:]\\._]*)\\^", "\\1_xXx_", fml_char)
		fml = as.formula(paste("~", fml_char_new))
	}

	t = terms(fml)
	all_var_names = attr(t, "term.labels")
	all_var_names = gsub("_xXx_", "^", all_var_names)
	all_var_names = gsub(":", "*", all_var_names) # for very special cases

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

format_se_type = function(x, width, by = FALSE){
    # we make 'nice' se types

    if(!grepl("\\(", x)){
        # means not clustered
        if(nchar(x) <= width) return(x)

        # We reduce each word to 3 letters (if needed)
        x_split = c("$", strsplit(x, "")[[1]]) # we add a non-letter flag, marking the beginning
        x_split_new = x_split
        end_word = length(x_split)
        non_letter_flag = grepl("[^[:alpha:]]", x_split) * (1:end_word)
        letter_flag = grepl("[[:alpha:]]", x_split) * (1:end_word)
        while(TRUE){
            start_word = which.max(non_letter_flag[1:end_word]) + 1
            # we truncate
            word_length = end_word - start_word + 1
            slack = length(x_split_new) - (width + 1)
            letters_to_rm = min(word_length - 4, slack)
            if(letters_to_rm > 0){
                i_max = end_word - letters_to_rm
                x_split_new = x_split_new[-((i_max+1):end_word)]
                x_split_new[i_max] = "."
            }

            lf = letter_flag[1:(start_word - 1)]
            if(all(lf == 0)) break

            # new end_word
            end_word = which.max(lf)
        }

        return(paste(x_split_new[-1], collapse = ""))
    } else if(x == "NA (not-available)"){
        return("not available")
    }

    # Now the FEs
    all_fe = gsub(".+\\((.+)\\)", "\\1", x)

    all_fe_split = gsub(" ", "", strsplit(all_fe, "&")[[1]])
    n_fe = length(all_fe_split)
    n_char = nchar(all_fe_split)

    if(n_fe == 1 && !grepl("\\^", all_fe_split[1])){
        if(by){
            se_formatted = paste0("by: ", all_fe_split[1])
        } else {
            se_formatted = paste0("1-way: ", all_fe_split[1])
        }

        if(nchar(se_formatted) > width){
            se_formatted = paste0(substr(se_formatted, 1, width - 2), "..")
        }
        return(se_formatted)
    }

    nb = ifelse(by, 3, 6)

    if(width < nb + sum(n_char) + (n_fe-1) * 3){
        qui = n_char > 5
        for(i in which(qui)){
            if(grepl("\\^", all_fe_split[i])){
                single_split = strsplit(all_fe_split[i], "\\^")[[1]]
                qui_bis = nchar(single_split) > 4
                single_split[qui_bis] = paste0(substr(single_split[qui_bis], 1, 3), ".")
                all_fe_split[i] = paste(single_split, collapse = "^")
            } else {
                all_fe_split[i] = paste0(substr(all_fe_split[i], 1, 4), ".")
            }
        }
    }

    if(by){
        se_formatted = paste0("by: ", paste(all_fe_split, collapse = " & "))
    } else {
        se_formatted = paste0(n_fe, "-way: ", paste(all_fe_split, collapse = " & "))
    }


    # if still too large, we trim right
    if(nchar(se_formatted) > width){
        # se_formatted = gsub("-way: ", "way: ", se_formatted)
        se_formatted = paste0(substr(se_formatted, 1, width - 2), "..")
    }

    se_formatted
}

format_se_type_latex = function(x, dict = c(), inline = FALSE){
    # we make 'nice' se types

    if(!grepl("\\(", x)){
        # means not clustered
        # we escape all
        return(escape_all(x))
    }

    # Now the FEs
    all_fe = gsub(".+\\((.+)\\)", "\\1", x)

    all_fe_split = gsub(" ", "", strsplit(all_fe, "&")[[1]])
    n_fe = length(all_fe_split)

    # Renaming the FEs

    all_fe_format = c()
    for(i in 1:length(all_fe_split)){
        fe = all_fe_split[i]

        if(fe %in% names(dict)){
            all_fe_format[i] = dict[fe]
        } else if(grepl("\\^", fe)){
            fe_split = strsplit(fe, "\\^")[[1]]
            who = fe_split %in% names(dict)
            fe_split[who] = dict[fe_split[who]]
            all_fe_format[i] = paste(fe_split, collapse = "$\\times$")
        } else {
            all_fe_format[i] = fe
        }
    }

    fe_format = paste(all_fe_format, collapse = " \\& ")

    # Full string
    nb = c("One", "Two", "Three", "Four")
    nway = paste0(nb[n_fe], "-way")

    if(inline){
        # se_formatted = paste0(nway, ": ", fe_format)
        se_formatted = fe_format
    } else {
        se_formatted = paste0(nway, " (", fe_format, ")")
    }
    # se_formatted = fe_format

    escape_latex(se_formatted)
}

tex_star = function(x){
    qui = nchar(x) > 0
    x[qui] = paste0("$^{", x[qui], "}$")
    x
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
    if(is.null(object$obsRemoved) && is.null(object$obs_selection)){
        return(x)
    }

    res = rep(NA, object$nobs_origin)
    qui = 1:object$nobs_origin

    if(!is.null(object$obsRemoved)){
        qui = qui[-object$obsRemoved]
    }

    for(i in seq_along(object$obs_selection)){
        qui = qui[object$obs_selection[[i]]]
    }

    res[qui] = x

    return(res)
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

error_sender = function(expr, ..., clean, up = 0){
    res = tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "try-error"))

    if("try-error" %in% class(res)){
        set_up(1 + up)
        msg = paste(..., collapse = "")
        if(!missing(clean)){
            stop_up(msg, gsub(clean, "", res))
        } else {
            stop_up(msg, res)
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
    # - expands interactions with :: (note that this will be deprecated)
    # - protects powers: x^3 => I(x^3)
    #
    # fml = sw(f(y, 1:2)) ~ x1 + l(x2, 1:2) + x2^2 | fe1 | y ~ z::e + g^3

    fml_text = deparse_long(fml)

    isPanel = grepl("(^|[^\\._[:alnum:]])(f|d|l)\\(", fml_text)
    isPower = grepl("^", fml_text, fixed = TRUE)
    isInteract = grepl("[^:]::[^:]", fml_text)

    if(isPanel || isInteract){
        # We rewrite term-wise

        fml_parts = fml_split(fml, raw = TRUE)
        n_parts = length(fml_parts)

        #
        # LHS
        #

        # only panel: no power (bc no need), no interact

        # We tolerate multiple LHS and expansion
        lhs_text = fml_split(fml, 1, text = TRUE, split.lhs = TRUE)
        if(isPanel){

            if(grepl("^(c|(c?(stepwise|sw)0?)|list)\\(", lhs_text)){
                lhs_text2eval = gsub("^(c|(c?(stepwise|sw)0?|list))\\(", "stepwise(", lhs_text)
                lhs_names = eval(parse(text = lhs_text2eval))
            } else {
                lhs_names = lhs_text
            }

            lhs_all = error_sender(expand_lags_internal(lhs_names),
                                   "Problem in the formula regarding lag/leads: ", clean = "__expand")

            if(length(lhs_all) > 1){
                lhs_fml = paste("c(", paste(lhs_all, collapse = "+"), ")")
            } else {
                lhs_fml = lhs_all
            }

            lhs_text = lhs_fml
        }

        #
        # RHS
        #

        # power + panel + interact

        if(isPower){
            # rhs actually also contains the LHS
            rhs_text = deparse_long(fml_parts[[1]])
            rhs_text = gsub("([\\.[:alpha:]][[:alnum:]\\._]*\\^[[:digit:]]+)", "I(\\1)", rhs_text)
            fml_rhs = as.formula(rhs_text)
        } else {
            fml_rhs = fml_maker(fml_parts[[1]])
        }

        rhs_terms = attr(terms(fml_rhs), "term.labels")

        if(isPanel){
            rhs_terms = error_sender(expand_lags_internal(rhs_terms),
                                     "Problem in the formula regarding lag/leads: ", clean = "__expand")
        }

        if(isInteract){
            rhs_terms = error_sender(expand_interactions_internal(rhs_terms),
                                     "Error in the interaction in the RHS of the formula: ")
        }

        rhs_text = paste(rhs_terms, collapse = "+")

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

                } else if(isPanel){

                    fml_fixef = fml_maker(fml_parts[[2]])
                    fml_fixef_text = deparse_long(fml_fixef)

                    if(grepl("(l|d|f)\\(", fml_fixef_text)){
                        # We need to make changes
                        # 1st: for terms to work => we change ^ if present (sigh)

                        do_sub = grepl("^", fml_fixef_text, fixed = TRUE)

                        if(do_sub){
                            fml_fixef = as.formula(gsub("^", "__impossible_var__", fml_fixef_text, fixed = TRUE))
                        }

                        fixef_terms = attr(terms(fml_fixef), "term.labels")
                        fixef_text = error_sender(expand_lags_internal(fixef_terms),
                                                  "Problem in the formula regarding lag/leads: ", clean = "__expand")

                        if(do_sub){
                            fixef_text = gsub("__impossible_var__", "^", fixef_text, fixed = TRUE)
                        }

                        fml_fixef = as.formula(paste("~", paste(fixef_text, collapse = "+")))

                    }
                } else {
                    fml_fixef = fml_maker(fml_parts[[2]])
                }
            }

            #
            # IV
            #

            if(n_parts == 3 || !is_fe){
                fml_iv = fml_maker(fml_parts[[n_parts]])

                fml_iv = fixest_fml_rewriter(fml_iv)$fml
            }
        }

        fml_new = merge_fml(fml_linear, fml_fixef, fml_iv)

    } else if(isPower){
        # It's faster not to call terms
        fml_text = gsub("([\\.[:alpha:]][[:alnum:]\\._]*\\^[[:digit:]]+)", "I(\\1)", fml_text)
        fml_new = as.formula(fml_text)

    } else {
        res = list(fml = fml, isPanel = FALSE)
        return(res)
    }

    res = list(fml = fml_new, isPanel = isPanel)

    return(res)
}

#### ................. ####
#### Aditional Methods ####
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

	if(object$method == "feols"){
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
#'
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
coef.fixest = coefficients.fixest = function(object, ...){
	object$coefficients
}

#' @rdname coef.fixest
"coefficients.fixest"


#' Extracts fitted values from a \code{fixest} fit
#'
#' This function extracts the fitted values from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. The fitted values that are returned are the \emph{expected predictor}.
#'
#' @inheritParams nobs.fixest
#'
#' @param type Character either equal to \code{"response"} (default) or \code{"link"}. If \code{type="response"}, then the output is at the level of the response variable, i.e. it is the expected predictor \eqn{E(Y|X)}. If \code{"link"}, then the output is at the level of the explanatory variables, i.e. the linear predictor \eqn{X\cdot \beta}.
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
fitted.fixest = fitted.values.fixest = function(object, type = c("response", "link"), ...){

    # Checking the arguments
    validate_dots(suggest_args = "type")

	type = match.arg(type)

	if(type == "response" || object$method == "feols"){
		res = object$fitted.values
	} else if(!is.null(object$mu)){
		res = object$mu
	} else if(object$method == "femlm"){
		family = object$family
		famFuns = switch(family,
							  poisson = ml_poisson(),
							  negbin = ml_negbin(),
							  logit = ml_logit(),
							  gaussian = ml_gaussian())

		res = famFuns$linearFromExpected(object$fitted.values)
	} else {
		res = object$family$linkfun(object$fitted.values)
	}

	# Nota: obs can be removed: either because of NA, either because perfect fit
	# Shall I put perfect fit as NA since they're out of the estimation???
	# Still pondering...
	# Actually adding them means a lot of work to ensure consitency (also in predict...)
	res = fill_with_na(res, object)

	res
}

#' @rdname fitted.fixest
#' @method fitted.values fixest
"fitted.values.fixest"

#' Extracts residuals from a \code{fixest} object
#'
#' This function extracts residuals from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#'
#' @param type A character scalar, either \code{"response"} (default), \code{"deviance"}, \code{"pearson"}, or \code{"working"}. Note that the \code{"working"} corresponds to the residuals from the weighted least square and only applies to \code{\link[fixest]{feglm}} models.
#' @param ... Not currently used.
#'
#' @details
#' The residuals returned are the difference between the dependent variable and the expected predictor.
#'
#' @return
#' It returns a numeric vector of the length the number of observations used for the estimation.
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
resid.fixest = residuals.fixest = function(object, type = c("response", "deviance", "pearson", "working"), ...){

    check_arg_plus(type, "match")

    method = object$method
    family = object$family

    r = object$residuals
    w = object[["weights"]]

    if(is.null(r)){
        stop("The method 'resid.fixest' cannot be applied to a 'lean' summary. Please apply it to the estimation object directly.")
    }

    if(method == "feols" || (method %in% c("feNmlm", "femlm") && family == "gaussian")){

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

    res = fill_with_na(res, object)

    res
}

#' @rdname resid.fixest
"residuals.fixest"

#' Predict method for \code{fixest} fits
#'
#' This function obtains prediction from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#' @inheritParams fitted.fixest
#'
#' @param newdata A data.frame containing the variables used to make the prediction. If not provided, the fitted expected (or linear if \code{type = "link"}) predictors are returned.
#' @param ... Not currently used.
#'
#' @return
#' It returns a numeric vector of length equal to the number of observations in argument \code{newdata}.
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
#' res = femlm(Sepal.Length ~ Petal.Length | Species, iris)
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
predict.fixest = function(object, newdata, type = c("response", "link"), ...){

    # Checking the arguments
    validate_dots(suggest_args = "type")

	# Controls
	type = match.arg(type)

	# if newdata is missing
	if(missing(newdata)){
		if(type == "response" || object$method == "feols"){
			res = object$fitted.values
		} else if(object$method == "femlm") {
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

	    res = fill_with_na(res, object)

		return(res)
	}

	if(!is.matrix(newdata) && !"data.frame" %in% class(newdata)){
		stop("Argument 'newdata' must be a data.frame.")
	}

	# we ensure it really is a clean data.frame
	newdata = as.data.frame(newdata)

	# We deconstruct it in four steps:
	# 1) cluster
	# 2) linear
	# 3) non-linear
	# 4) offset

	# + step 0: panel setup

	n = nrow(newdata)

	# message for the user NOT to use newdata if its the same data set
	# The user is not stupid: so just once
	mc = match.call()
	if(n == object$nobs_origin && deparse_long(object$call$data) == deparse_long(mc$newdata) && getFixest_notes()){
	    dont_warn = getOption("fixest_predict_dont_warn")
	    if(!isTRUE(dont_warn)){
	        message("NOTE: It looks like the data in 'newdata' is the same as the one used to run the regression. If so, to predict() on the existing data, you can leave the argument 'newdata' as missing, this is faster.")
	        options("fixest_predict_dont_warn" = TRUE)
	    }
	}

	# NOTA 2019-11-26: I'm pondering whether to include NA-related messages
	# (would it be useful???)


	# STEP 0: panel setup

	fml = object$fml
	if(check_lag(fml)){
	    if(!is.null(object$panel.info)){
	        if(is.null(attr(newdata, "panel_info"))){
                # We try to recreate the panel
	            if(any(!names(object$panel.info) %in% c("", "data", "panel.id"))){
	                # This was NOT a standard panel creation
	                stop("The estimation contained lags/leads and the original data was a 'fixest_panel' while the new data is not. Please set the new data as a panel first with the function panel(). NOTA: the original call to panel was:\n", deparse_long(object$panel.info))
	            } else {
	                panel__meta__info = panel_setup(newdata, object$panel.id, from_fixest = TRUE)
	            }
	        } else {
	            panel__meta__info = attr(newdata, "panel_info")
	        }
	    } else {
	        panel__meta__info = panel_setup(newdata, object$panel.id, from_fixest = TRUE)
	    }
	}

	# 1) Cluster

	# init cluster values
	value_cluster = 0

	fixef_vars = object$fixef_vars
	if(!is.null(fixef_vars)){

		n_cluster = length(fixef_vars)

		# Extraction of the FEs
		id_cluster = list()
		for(i in 1:n_cluster){
			# checking if the variable is in the newdata
		    fe_var = fixef_vars[i]
			variable = all.vars(parse(text = fe_var))
			isNotHere = !variable %in% names(newdata)
			if(any(isNotHere)){
				stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a fixed-effect variable).")
			}

			# The values taken by the FE variable
			fixef_values_possible = attr(object$fixef_id[[i]], "fixef_names")

			# Checking if ^ is present
			if(grepl("\\^", fe_var)){
			    # If fastCombine was used => we're screwed, impossible to recover

			    if(!is.character(fixef_values_possible)){
			        stop("You cannot use predict() based on the initial regression since the fixed-effect '", variable, "' was combined using an algorithm dropping the FE values (but fast). Please re-run the regression using the argument 'combine.quick=FALSE'.")
			    }

			    fe_var_new = gsub("([[:alpha:]_\\.][[:alnum:]_\\.]*(\\^[[:alpha:]_\\.][[:alnum:]_\\.]*)+)",
			                    "combine_clusters(\\1)", fe_var)

			    fe_var = gsub("\\^", ", ", fe_var_new)
			}

			# Obtaining the vector of clusters
			cluster_current = eval(parse(text = fe_var), newdata)

			cluster_current_num = unclass(factor(cluster_current, levels = fixef_values_possible))
			id_cluster[[i]] = cluster_current_num
		}

		names(id_cluster) = fixef_vars

		# Value of the cluster coefficients
		cluster_coef = fixef(object, sorted = FALSE)

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
		        variable = all.vars(parse(text = slope_vars_unik[i]))
		        isNotHere = !variable %in% names(newdata)
		        if(any(isNotHere)){
		            stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a variable with varying slope).")
		        }

		        slope_var_list[[slope_vars_unik[i]]] = eval(parse(text = slope_vars_unik[i]), newdata)
		    }

		    # Adding the FE values
		    for(var in fixef_vars){
		        cluster_current_num = id_cluster[[var]]
		        cluster_coef_current = cluster_coef[[var]]

		        value_cluster = value_cluster + cluster_coef_current[cluster_current_num]
		    }

		    # Adding the slopes
		    for(i in seq_along(slope_vars)){

		        cluster_current_num = id_cluster[[slope_fe[i]]]
		        cluster_coef_current = cluster_coef[[slope_terms[i]]]

		        value_cluster = value_cluster + cluster_coef_current[cluster_current_num] * slope_var_list[[slope_vars[i]]]
		    }


		} else {
		    # Adding only FEs
		    for(i in 1:n_cluster){
		        cluster_current_num = id_cluster[[i]]
		        cluster_coef_current = cluster_coef[[i]]

		        value_cluster = value_cluster + cluster_coef_current[cluster_current_num]
		    }
		}

		# dropping names
		value_cluster = as.vector(value_cluster)
	}

	# 2) Linear values

	coef = object$coefficients

	value_linear = 0
	rhs_fml = fml_split(fml, 1)
	if(grepl("[^:]::[^:]", deparse_long(rhs_fml[[3]]))){
	    new_fml = expand_interactions(rhs_fml)
	    linear.varnames = all.vars(new_fml[[3]])
	} else {
	    linear.varnames = all.vars(rhs_fml[[3]])
	}

	if(length(linear.varnames) > 0){
		# Checking all variables are there
		varNotHere = setdiff(linear.varnames, names(newdata))
		if(length(varNotHere) > 0){
			stop("The variable", enumerate_items(varNotHere, "s.quote"), " used to estimate the model (in fml) ", ifsingle(varNotHere, "is", "are"), " missing in the data.frame given by the argument 'newdata'.")
		}

		# We check if it's a panel or not (if so, we need to create it...)

		# we create the matrix
		matrix_linear = try(fixest_model_matrix(rhs_fml, newdata, i_noref = TRUE), silent = TRUE)
		if("try-error" %in% class(matrix_linear)){
		    stop("Error when creating the linear matrix: ", matrix_linear)
		}

		keep = intersect(names(coef), colnames(matrix_linear))
		value_linear = value_linear + as.vector(matrix_linear[, keep, drop = FALSE] %*% coef[keep])
	}

	# 3) Non linear terms

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

	# 4) offset value

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

	value_predicted = value_cluster + value_linear + value_NL + value_offset

	if(type == "link" || object$method == "feols"){
		res = value_predicted
	} else if(object$method == "femlm") {
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

	res
}

#' Extracts the variance/covariance of a \code{femlm} fit
#'
#' This function extracts the variance-covariance of estimated parameters from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams summary.fixest
#' @inheritParams nobs.fixest
#'
#' @param attr Logical, defaults to \code{FALSE}. Whether to include the attributes describing how the VCOV was computed.
#' @param ... Other arguments to be passed to \code{\link[fixest]{summary.fixest}}.
#'
#' The computation of the VCOV matrix is first done in \code{\link[fixest]{summary.fixest}}.
#'
#' @details
#' For an explanation on how the standard-errors are computed and what is the exact meaning of the arguments, please have a look at the dedicated vignette: \href{https://cran.r-project.org/package=fixest/vignettes/standard_errors.html}{On standard-errors}.
#'
#' @return
#' It returns a \eqn{N\times N} square matrix where \eqn{N} is the number of variables of the fitted model.
#' This matrix has an attribute \dQuote{type} specifying how this variance/covariance matrix has been computed (i.e. if it was created using a heteroskedascity-robust correction, or if it was clustered along a specific factor, etc).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{confint.fixest}}, \code{\link[fixest]{resid.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{fixef.fixest}}.
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
#' # By default, in the presence of FEs
#' # the VCOV is clustered along the first FE
#' vcov(est_pois)
#'
#' # Heteroskedasticity-robust VCOV
#' vcov(est_pois, se = "hetero")
#'
#' # "clustered" VCOV (with respect to the Product factor)
#' vcov(est_pois, se = "cluster", cluster = trade$Product)
#' # another way to make the same request:
#' # note that previously arg. se was optional since deduced from arg. cluster
#' vcov(est_pois, cluster = "Product")
#' # yet another way:
#' vcov(est_pois, cluster = ~Product)
#'
#' # Another estimation without fixed-effects:
#' est_pois_simple = femlm(Euros ~ log(dist_km) + log(Year), trade)
#'
#' # We can still get the clustered VCOV,
#' # but we need to give the argument cluster:
#' vcov(est_pois_simple, cluster = ~Product)
#'
#'
vcov.fixest = function(object, se, cluster, dof = getFixest_dof(), attr = FALSE, forceCovariance = FALSE, keepBounded = FALSE, ...){
	# computes the clustered vcov

    check_arg(attr, "logical scalar")

    if(!"nframes_up" %in% names(match.call())){
        # condition means client call
        validate_dots(suggest_args = c("se", "cluster", "dof"))
    }

	if(!is.null(object$onlyFixef)){
		# means that the estimation is done without variables
		stop("No explanatory variable was used: vcov() cannot be applied.")
	}

	# Default behavior se:
	suffix = ""
	if(missnull(se)){

		if(missing(cluster)){

		    se_default = getFixest_se()

		    n_fe = length(object$fixef_id)

		    if(n_fe == 0){
		        se = se_default$no_FE
		    } else if(n_fe == 1){
		        se = se_default$one_FE
		    } else {
		        se = se_default$two_FE
		    }

		} else {
			if("formula" %in% class(cluster)){
				# we just find the nway clustering and do only minor control

			    cluster = formula(cluster) # regularization to check it

				if(length(cluster) != 2){
					stop("If argument cluster is to be a formula, it must be one sided: e.g. ~fe_1+fe_2.")
				}

				all_vars = fml2varnames(cluster)
				doEval = TRUE
				nway = length(all_vars)

				if(nway > 4){
					stop("From argument 'cluster', ", nway, "-way clustering is deduced. However, only up to fourway clustering is supported.")
				}

			} else if(length(cluster) <= 4){
				nway = length(cluster)

			} else if(length(cluster) %in% c(object$nobs, object$nobs_origin)){
				nway = 1

			} else {
				stop("The length of argument cluster (", length(cluster), ") is invalid, no clustering can be deduced. Alternatively, you could use a formula.")
			}

			se = switch(nway, "1"="cluster", "2"="twoway", "3"="threeway", "4"="fourway")

			suffix = paste0(" [note: ", nway, "-way clustering deduced by cluster length.]")
		}

	}

	# Argument se
	if(!length(se) == 1) {
	    stop("Argument 'se' must be of length 1.")
	}

	if(isScalar(se) && se %in% 1:4){
	    # we allow for integer values
	    se = c("cluster", "twoway", "threeway", "fourway")[se]
	} else if(!is.character(se)){
		stop("Argument 'se' must be a character scalar equal to: 'standard', 'hetero', 'cluster', 'twoway', 'threeway' or 'fourway'.")
    }

	se.val = NULL
	try(se.val <- match.arg(se, c("standard", "white", "hetero", "cluster", "twoway", "threeway", "fourway", "1", "2", "3", "4")), silent = TRUE)
	if(is.null(se.val)){
		stop("Invalid argument 'se'. It should be equal to one of 'standard', 'hetero', 'cluster', 'twoway', 'threeway' or 'fourway'.")
	}

	if(se.val == "white") se.val = "hetero"

	dots = list(...)
	if(is.null(dots$nframes_up)){
		nframes_up = 1
	} else {
		nframes_up = dots$nframes_up + 1
	}

	# DoF related
	check_arg(dof, "class(dof.type)", .message = "The argument 'dof.type' must be an object created by the function dof().")

    dof.fixef.K = dof$fixef.K
    dof.adj = dof$adj
    is_exact = dof$fixef.force_exact
    is_cluster = dof$cluster.adj
    is_cluster_min = dof$cluster.df == "min"
    is_t_min = dof$t.df == "min"

	#
	# non-linear: handling bounded parameters
	#

	# We handle the bounded parameters:
	isBounded = object$isBounded
	if(is.null(isBounded)){
		isBounded = rep(FALSE, length(object$coefficients))
	}

	if(any(isBounded)){
		if(keepBounded){
			# we treat the bounded parameters as regular variables
			myScore = object$scores
			object$cov.unscaled = solve(object$hessian)
		} else {
			myScore = object$scores[, -which(isBounded), drop = FALSE]
		}
	} else {
		myScore = object$scores
	}


	#
	# Core function
	#

	n = object$nobs
	n_fe = n_fe_ok = length(object$fixef_id)

	# we adjust the fixef sizes to account for slopes
	fixef_sizes_ok = object$fixef_sizes
	isSlope = FALSE
	if(!is.null(object$fixef_terms)){
	    isSlope = TRUE
	    # The drop the fixef_sizes for only slopes
	    fixef_sizes_ok[object$slope_flag < 0] = 0
	    n_fe_ok = sum(fixef_sizes_ok > 0)
	}

	# How do we choose K? => argument dof

	if(dof.fixef.K == "none"){
	    # we do it with "minus" because of only slopes
	    K = object$nparams
	    if(n_fe_ok > 0){
	        K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1))
	    }
	} else if(dof.fixef.K == "full" || se.val %in% c("standard", "hetero")){
	    K = object$nparams
	    if(is_exact && n_fe >= 2 && n_fe_ok >= 1){
	        fe = fixef(object, notes = FALSE)
	        K = K + (n_fe_ok - 1) - sum(attr(fe, "references"))
	    }
	} else {
	    # nested
	    # we delay the adjustment
	    K = object$nparams
	}

	if(object$method == "feols"){
		if(se.val != "standard"){
			VCOV_raw = object$cov.unscaled / object$sigma2
		} else {
			VCOV_raw = object$cov.unscaled / ((n - 1) / (n - object$nparams))
		}
	} else {
		VCOV_raw = object$cov.unscaled
	}


	# Small sample adjustment
	correction.dof = ifelse(dof.adj, (n - 1) / (n - K), 1)


	# information on the variable used for the clustering
	type_info = ""

	is_nested = c()
	if(anyNA(VCOV_raw)){

		if(!forceCovariance){
		    last_warn = getOption("fixest_last_warning")
		    if(is.null(last_warn) || (proc.time() - last_warn)[3] > 1){
		        warning("Standard errors are NA because of likely presence of collinearity. Use function collinearity() to detect collinearity problems.", call. = FALSE)
		    }

		    attr(VCOV_raw, "type") = "NA (not-available)"

			return(VCOV_raw)
		} else {
			# VCOV_raw_forced = MASS::ginv(object$hessian)
			# if(anyNA(VCOV_raw_forced)) {
			# 	stop("The covariance matrix could not be 'forced'.")
			# }
		    info_inv = cpp_cholesky(object$hessian)
		    if(!is.null(info_inv$all_removed)){
		        # Means all variables are collinear! => can happen when using FEs
		        stop("All variables have virtually no effect on the dependent variable. Covariance is not defined.")
		    }

		    VCOV_raw_forced = info_inv$XtX_inv
		    if(any(info_inv$id_excl)){
		        n_collin = sum(info_inv$all_removed)
		        message("NOTE: ", n_letter(n_collin), " variable", plural(n_collin, "s.has"), " been found to be singular.")

		        VCOV_raw_forced = cpp_mat_reconstruct(VCOV_raw_forced, info_inv$id_excl)
		        VCOV_raw_forced[, info_inv$id_excl] = NA
		        VCOV_raw_forced[info_inv$id_excl, ] = NA
		    }

			object$cov.unscaled = VCOV_raw_forced
			return(vcov(object, se=se.val, cluster=cluster, dof=dof))
		}

	} else if(se.val == "standard"){

		vcov = VCOV_raw * correction.dof

	} else if(se.val == "hetero"){

	    # we make a n/(n-1) adjustment to match vcovHC(type = "HC1")
		vcov = crossprod(myScore %*% VCOV_raw) * correction.dof * ifelse(is_cluster, n/(n-1), 1)

	} else {
		# Clustered SD!
		nway = switch(se.val, cluster=1, twoway=2, threeway=3, fourway=4)

		# MISC

		# used twice later:
		msgRemoved = ""
		if(!is.null(object$call$na.rm) && object$call$na.rm){
			msgRemoved = " (additionnaly from the observations removed from the original estimation)"
		}

		#
		# Controls
		#

		# Controlling the clusters
		do.unclass = TRUE
		if(missing(cluster) || is.null(cluster)){

			if(is.null(object$fixef_id)){
				stop("To display clustered standard errors, you must provide the argument 'cluster'", ifelse(!missing(cluster), ", currently it is equal to NULL", ""), ".")

			} else if(length(object$fixef_id) < nway) {
				stop(nway, "-way clustering is asked for but the estimation was not performed with ", nway, " or more fixed-effects: You must provide the argument 'cluster' with ", nway, " clusters.")

			} else {
				cluster = object$fixef_id[1:nway]

				type_info = paste0(" (", paste0(object$fixef_vars[1:nway], collapse = " & "), ")")

				is_nested = 1:nway

				# in that specific case, there is no need of doing unclass.factor because already done
				do.unclass = FALSE
			}

		} else {

			#
			# Handle formulas
			#

			doEval = FALSE
			if("formula" %in% class(cluster)){

			    cluster = formula(cluster) # regularization to check it

				if(length(cluster) != 2){
					stop("If argument cluster is to be a formula, it must be one sided: e.g. ~dum_1+dum_2.")
				}

				all_var_names = fml2varnames(cluster)

				if(length(all_var_names) != nway){
					stop("Asked for ", nway, "-way clustering but evaluating argument cluster leads to ", length(all_var_names), " clusters (", enumerate_items(all_var_names), "). Please provide exactly ", nway, " clusters.")
				}

				cluster = all_var_names # Now a regular character vector

				doEval = TRUE
			}

			if(length(cluster) == nway && is.character(cluster)){

				if(all(cluster %in% object$fixef_vars)){
					# cluster == names of clusters used in the estimation
					type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")

					is_nested = which(names(object$fixef_id) %in% cluster)

					cluster = object$fixef_id[cluster]

					do.unclass = FALSE

				} else {
					cluster = gsub(" *", "", cluster)
					if(!doEval){
						is_ok = grepl("^[\\.[:alpha:]][[:alnum:]_\\.]*(\\^[\\.[:alpha:]][[:alnum:]_\\.]*)*$", cluster)
						if(any(!is_ok)){
							stop("In argument cluster, only variable names and the '^' operator are accepted. The expression", enumerate_items(cluster[!is_ok], "s.is"), " not valid.\nAlternatively, you can use a list of vectors.", suffix)
						}

					}

					cluster_fml = as.formula(paste0("~", paste0(cluster, collapse = " + ")))
					all_vars = all.vars(cluster_fml)

					if(all(all_vars %in% object$fixef_vars) || all(cluster %in% object$fixef_vars)){
						# Means dum_1^dum_2 with dum_1 and dum_2 used as clusters

						cluster_names = cluster
						type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")

						cluster = list()
						for(i in 1:nway){
							cname = cluster_names[i]
							if(cname %in% object$fixef_vars){
								cluster[[i]] = object$fixef_id[[cname]]
							} else {
								# combination
								if(grepl("^", cname, fixed = TRUE)){
									value_text = gsub("\\^", ", ", cname)
									value_text = paste0("combine_clusters_fast(", value_text, ")")
								}

								value_call = parse(text = value_text)
								value = eval(value_call, object$fixef_id)
								cluster[[i]] = value
							}
						}

					} else {
						# we try to get the variable from the base used in the estimation
						var2fetch = setdiff(all_vars, object$fixef_vars)

						# evaluation
						data = NULL
						try(data <- eval(object$call$data, parent.frame(nframes_up)), silent = TRUE)
						dataName = object$call$data

						if(is.null(data)){
							# We try to provide an informative error message
							stop("Cannot apply ", nway, "-way clustering with current 'cluster' argument. Variable", enumerate_items(var2fetch, "s.is.past"), " not used as fixed-effects in the estimation. We tried to fetch ", ifelse(length(var2fetch) == 1, "this variable", "these variables"), " in the original database in the parent.frame -- but the data doesn't seem to be there anymore (btw it was ", deparse_long(dataName), "). Alternatively, use a list of vectors.", suffix)
						}

						data = as.data.frame(data)

						# we check the variables are there
						# we use all_vars and not var2fetch: safer to catch all variables (case clustvar^datavar)
						if(any(!all_vars %in% names(data))){
							var_pblm = setdiff(all_vars, names(data))
							stop("In argument 'cluster', the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset. Alternatively, use a list of vectors.", suffix)
						}

						# we check length consistency
						if(NROW(data) != object$nobs_origin){
							stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(var2fetch, "s"), " in the original dataset (", deparse_long(dataName), "), yet the dataset doesn't have the same number of observations as was used in the estimation (", NROW(data), " instead of ", object$nobs_origin, ").", suffix)
						}

						if(length(object$obsRemoved) > 0){
							data = data[-object$obsRemoved, all_vars, drop = FALSE]
						} else {
							data = data[, all_vars, drop = FALSE]
						}

						for(i in seq_along(object$obs_selection)){
						    data = data[object$obs_selection[[i]], , drop = FALSE]
						}

						# Final check: NAs
						if(anyNA(data)){
							varsNA = sapply(data, anyNA)
							varsNA = names(varsNA)[varsNA]
							stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(varsNA, "s"), " in the original dataset (", deparse_long(dataName), "). But ", ifsingle(varsNA, "this variable", "these variables"), " contain", ifsingle(varsNA, "s", ""), " NAs", msgRemoved, ". This is not allowed.", suffix)
						}

						# We create the cluster list
						cluster_names = cluster
						type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")

						cluster = list()
						for(i in 1:nway){
							cname = cluster_names[i]
							if(cname %in% object$fixef_vars){
								cluster[[i]] = object$fixef_id[[cname]]

							} else if(cname %in% names(data)){
							    # data is already of the right size
								cluster[[i]] = data[[cname]]

							} else {
								# combination
								if(grepl("^", cname, fixed = TRUE)){
									value_text = gsub("\\^", ", ", cname)
									value_text = paste0("combine_clusters_fast(", value_text, ")")
								} else {
								    value_text = cname
								}

								value_call = parse(text = value_text)
								cluster[[i]] = eval(value_call, data)
							}
						}

					}

				}
			} else if(length(cluster) == nway && is.numeric(cluster)){
			    # You can use a number to tell which cluster to use

			    if(length(object$fixef_vars) == 0){
			        stop("You can use an integer in the argument 'cluster' only when there have been fixed-effects in the estimation. Currenlty this is not the case. Alternatively, arg. 'cluster' can be a formula, a vector of variables or a list of vectors.")
			    }

			    if(!all(cluster %% 1 == 0) || any(cluster < 1 | cluster > 4)){
			        msg = ifelse(!all(cluster %% 1 == 0), "it is not made of integers", "it contains values different from 1 to 4")
			        stop("Argument 'cluster' can be a numeric vector, if so it must have integer values between 1 and 4 (currently ", msg, ").")
			    }

			    if(length(object$fixef_vars) < max(cluster)){
			        nb_name = c("1st", "2nd", "3rd", "4th")
			        stop("In argument 'cluster', it is requested to cluster along the ", nb_name[max(cluster)], " fixed-effect, however the estimation was done with only ", length(object$fixef_vars), " fixed-effects. Alternatively, arg. 'cluster' can be a formula, a vector of variables or a list of vectors.")
			    }

			    # Eventually, it should be all right by now
			    type_info = paste0(" (", paste0(object$fixef_vars[cluster], collapse = " & "), ")")
			    cluster = object$fixef_id[cluster]

			} else if(nway == 1){
				if(!is.list(cluster) && (isVector(cluster) || is.factor(cluster))){
					cluster = list(cluster)

				} else if(! (is.list(cluster) && length(cluster) == 1)){
					stop("For one way clustering, the argument 'cluster' must be either the name of a cluster variable (e.g. \"dum_1\"), a vector (e.g. data$dum_1), a list containing the vector of clusters (e.g. list(data$dum_1)), or a one-sided formula (e.g. ~dum_1). Currently the class of cluster is ", enumerate_items(class(cluster)), ".", suffix)

				}
			} else if(length(cluster) != nway){

				msgList = "a list of vectors"
				if(is.list(cluster)) msgList = "a vector of variables names"
				stop(nway, "-way clustering is asked for, but the length of argument 'cluster' is ", length(cluster), " (it should be ", nway, "). Alternatively, you can use ", msgList, " or a one-sided formula.", suffix)

			} else if(!is.list(cluster)){
				stop("For ", nway, "-way clustering, the argument 'cluster' must be either a vector of cluster variables (e.g. c(\"", paste0("dum_", 1:nway, collapse = "\", \""), "\")), a list containing the vector of clusters (e.g. data[, c(\"", paste0("dum_", 1:nway, collapse = "\", \""), "\")]), or a one-sided formula (e.g. ~", paste0("dum_", 1:nway, collapse = "+"), "). Currently the class of cluster is: ", enumerate_items(class(cluster)), ".", suffix)
			}

			cluster = as.list(cluster)
		}

		# now we check the lengths:
		n_per_cluster = sapply(cluster, length)
		if(!all(diff(n_per_cluster) == 0)){
			stop("The vectors of the argument 'cluster' must be of the same length. Currently the lengths are: ", enumerate_items(n_per_cluster), ".")
		}

		# Either they are of the same length of the data
		if(n_per_cluster[1] != object$nobs){
			# Then two cases: either the user introduces the original data and it is OK
			if(n_per_cluster[1] == object$nobs_origin){
				# We modify the clusters
				for(i in 1:nway){
				    # first we take care of obsRemoved
				    if(!is.null(object$obsRemoved)){
				        cluster[[i]] = cluster[[i]][-object$obsRemoved]
				    }

				    for(j in seq_along(object$obs_selection)){
				        cluster[[i]] = cluster[[i]][object$obs_selection[[j]]]
				    }

				}
			} else {
				# If this is not the case: there is a problem
				stop("The length of the clusters (", n_per_cluster[1], ") does not match the number of observations in the estimation (", object$nobs, ").")
			}
		}

		# final NA check
		varsNA = sapply(cluster, anyNA)
		if(any(varsNA)){
			varsNA = which(varsNA)
			nb_name = c("1st", "2nd", "3rd", "4th")
			stop("In argument cluster, the ", enumerate_items(nb_name[varsNA]), " cluster variable", ifsingle(varsNA, " contains", "s contain"), " NAs", msgRemoved, ". This is not allowed.")
		}


		#
		# Calculus
		#

		# initialisation
		vcov = VCOV_raw * 0

		if(do.unclass){
			for(i in 1:nway){
				cluster[[i]] = quickUnclassFactor(cluster[[i]])
			}
		}

		# We recompute K
		if(dof.adj && dof.fixef.K == "nested" && n_fe_ok >= 1){
            if(do.unclass){
                # we need to find out which is nested
                is_nested = which(cpp_check_nested(object$fixef_id, cluster, object$fixef_sizes, n = n) == 1)
            } else {
                # no need to compute is_nested,
                # we created it earlier
            }

		    if(length(is_nested) == n_fe){
		        K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1))
		    } else {
		        if(is_exact && n_fe >= 2){
		            fe = fixef(object, notes = FALSE)
		            nb_ref = attr(fe, "references")

		            # Slopes are a pain in the neck!!!
		            if(length(is_nested) > 1){
	                    id_nested = intersect(names(nb_ref), names(object$fixef_id)[is_nested])
	                    nb_ref[id_nested] = object$fixef_sizes[id_nested]
		            }

		            total_refs = sum(nb_ref)

		            K = K - total_refs
		        } else {
		            K = K - (sum(fixef_sizes_ok[is_nested]) - sum(fixef_sizes_ok[is_nested] > 0))
		        }
		    }

		    K = max(K, length(object$coefficients) + 1)

		    correction.dof = (n - 1) / (n - K)
		}

		for(i in 1:nway){

			myComb = combn(nway, i)

			power = floor(1 + log10(sapply(cluster, max)))

			for(j in 1:ncol(myComb)){

				if(i == 1){
					index = cluster[[myComb[, j]]]
				} else if(i > 1){

					vars = myComb[, j]

					if(sum(power[vars]) > 14){
						myDots = cluster[vars]
						myDots$sep = "_"
						index = do.call("paste", myDots)
					} else {
						# quicker, but limited by the precision of integers
						index = cluster[[vars[1]]]
						for(k in 2:length(vars)){
							index = index + cluster[[vars[k]]]*10**sum(power[vars[1:(k-1)]])
						}
					}

					index = quickUnclassFactor(index)

				}

			    # When cluster.df == "min" => no dof here but later
				vcov = vcov + (-1)**(i+1) * vcovClust(index, VCOV_raw, myScore, dof = is_cluster && !is_cluster_min, do.unclass=FALSE)
			}
		}

		G_min = NULL
		if(is_cluster && is_cluster_min){
		    G_min = min(sapply(cluster, max))
		    correction.dof = correction.dof * G_min / (G_min - 1)
		}

		vcov = vcov * correction.dof

		if(is_t_min){
		    if(is.null(G_min)) G_min = min(sapply(cluster, max))

		    if(attr) base::attr(vcov, "G") = G_min
		}

	}

	if(any(diag(vcov) < 0)){
	    # We 'fix' it
	    e = eigen(vcov)
	    vcov = tcrossprod(e$vectors %*% diag(pmax(e$values, 1e-8)), e$vectors)
	    message("Variance contained negative values in the diagonal and was 'fixed' (a la Cameron, Gelbach & Miller 2011).")
	}

	sd.dict = c("standard" = "Standard", "hetero"="Heteroskedasticity-robust", "cluster"="Clustered", "twoway"="Two-way", "threeway"="Three-way", "fourway"="Four-way")

	if(attr){
	    base::attr(vcov, "type") = paste0(as.vector(sd.dict[se.val]), type_info)
	    base::attr(vcov, "dof.type") = paste0("dof(adj = ", dof.adj, ", fixef.K = '", dof.fixef.K, "', cluster.adj = ", is_cluster, ", cluster.df = '", dof$cluster.df, "', t.df = '", dof$t.df, "', fixef.force_exact = ", is_exact, ")")
	    base::attr(vcov, "dof.K") = K
	}



	vcov
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
confint.fixest = function(object, parm, level = 0.95, se, cluster, dof = getFixest_dof(), ...){

    # Checking the arguments
    validate_dots(suggest_args = c("parm", "level", "se", "cluster"),
                  valid_args = c("exact_dof", "forceCovariance", "keepBounded"))

	# Control
	if(!is.numeric(level) || !length(level) == 1 || level >= 1 || level <= .50){
		stop("The argument 'level' must be a numeric scalar greater than 0.50 and strictly lower than 1.")
	}

	# the parameters for which we should compute the confint
	all_params = names(object$coefficients)

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

	# The proper SE
	sum_object = summary(object, se = se, cluster = cluster, dof = dof, nframes_up = 1, ...)

	se_all = sum_object$se
	coef_all = object$coefficients

	# multiplicative factor
	val = (1 - level) / 2
	fact <- abs(qnorm(val))

	# The confints
	lower_bound = coef_all[parm_use] - fact * se_all[parm_use]
	upper_bound = coef_all[parm_use] + fact * se_all[parm_use]

	res = data.frame(lower_bound, upper_bound, row.names = parm_use)
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
#' est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, trade)
#'
#' # we add the variable log(Year)
#' est_2 <- update(est_pois, . ~ . + log(Year))
#'
#' # we add another fixed-effect: "Product"
#' est_3 <- update(est_2, . ~ . | . + Product)
#'
#' # we remove the fixed-effect "Origin" and the variable log(dist_km)
#' est_4 <- update(est_3, . ~ . - log(dist_km) | . - Origin)
#'
#' # Quick look at the 4 estimations
#' esttable(est_pois, est_2, est_3, est_4)
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

    if(isTRUE(object$fromFit)){
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
	    if(object$method == "feols"){
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

	if(n_parts > 2 + (object$method == "feols")){
	    stop("The update formula cannot have more than ", 2 + (object$method == "feols"), " parts for the method ", object$method, ".")
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
#' @param type A character scalar. Default is \code{type = "full"} which gives back a formula containing the linear part of the model along with the fixed-effects (if any) and the non-linear in parameters part (if any). If \code{type = "linear"} then only the linear formula is returned. If \code{type = "NL"} then only the non linear in parameters part is returned.
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
    validate_dots(suggest_args = "type")

    if(isTRUE(x$fromFit)){
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


#' Design matrix of a \code{femlm} model
#'
#' This function creates a design matrix of the linear part of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation. Note that it is only the linear part. The fixed-effects variables (which can be considered as factors) are excluded from the matrix.
#'
#' @method model.matrix fixest
#'
#' @inheritParams nobs.fixest
#'
#' @param data If missing (default) then the original data is obtained by evaluating the \code{call}. Otherwise, it should be a \code{data.frame}.
#' @param na.rm Default is \code{TRUE}. Should observations with NAs be removed from the matrix?
#' @param ... Not currently used.
#'
#' @return
#' It returns a design matrix.
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
#' # simple estimation on iris data, using "Species" fixed-effects
#' res = femlm(Sepal.Length ~ Sepal.Width*Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' head(model.matrix(res))
#'
#'
#'
model.matrix.fixest = function(object, data, na.rm = TRUE, ...){
	# We evaluate the formula with the past call

    # Checking the arguments
    validate_dots(suggest_args = "data")

    if(isTRUE(object$fromFit)){
        stop("model.matrix method not available for fixest estimations obtained from fit methods.")
    }

	# I) we obtain the right formula
	fml = object$fml

	# we kick out the intercept if there is presence of clusters
	if(attr(terms(fml), "intercept") == 1 && !is.null(object$fixef_vars)){
		fml = update(fml, . ~ . - 1)
	}

	# II) evaluation with the data
	original_data = FALSE
	if(missing(data)){
	    original_data = TRUE

		call_old = object$call

		data = NULL
		try(data <- eval(object$call$data, parent.frame()), silent = TRUE)

		if(is.null(data)){
			dataName = deparse_long(object$call$data)
			stop("To apply 'model.matrix.fixest', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", dataName, ").")
		}

	}

	# control of the data
	if(is.matrix(data)){
		if(is.null(colnames(data))){
			stop("If argument data is to be a matrix, its columns must be named.")
		}
		data = as.data.frame(data)
	}
	# The conversion of the data (due to data.table)
	if(!"data.frame" %in% class(data)){
		stop("The argument 'data' must be a data.frame or a matrix.")
	}

	data = as.data.frame(data)

	# Panel setup
	if(check_lag(fml)){
	    if(!is.null(object$panel.info)){
	        if(is.null(attr(data, "panel_info"))){
	            # We try to recreate the panel
	            if(any(!names(object$panel.info) %in% c("", "data", "panel.id"))){
	                # This was NOT a standard panel creation
	                stop("The original data set was a fixest_panel, now it isn't any more. Please restore the original data to a panel to perform model.matrix. NOTA: the original call to panel was:\n", deparse_long(object$panel.info))
	            } else {
	                panel__meta__info = panel_setup(data, object$panel.id, from_fixest = TRUE)
	            }
	        } else {
	            panel__meta__info = attr(data, "panel_info")
	        }
	    } else {
	        panel__meta__info = panel_setup(data, object$panel.id, from_fixest = TRUE)
	    }
	}

	linear.mat = fixest_model_matrix(fml, data)

	check_0 = FALSE
	if(original_data){
	    if(!is.null(object$obsRemoved)){
	        check_0 = TRUE
	        linear.mat = linear.mat[-object$obsRemoved, , drop = FALSE]
	    }

	    for(i in seq_along(object$obs_selection)){
	        check_0 = TRUE
	        linear.mat = linear.mat[object$obs_selection[[i]], , drop = FALSE]
	    }

	} else if(na.rm){
	    check_0 = TRUE
	    info = cpppar_which_na_inf_mat(linear.mat, nthreads = 1)

	    if(info$any_na_inf){
	        isNA_L = info$is_na_inf

	        if(sum(isNA_L) == nrow(linear.mat)){
	            warning("All observations contain NA values.")
	            return(linear.mat[-which(isNA_L), , drop = FALSE])
	        }

	        linear.mat = linear.mat[-which(isNA_L), , drop = FALSE]
	    }
	}

	if(check_0){
	    only_0 = cpppar_check_only_0(linear.mat, nthreads = 1)
	    if(all(only_0 == 1)){
	        stop("After removing NAs, not a single explanatory variable is different from 0.")

	    } else if(any(only_0 == 1)){
	        linear.mat = linear.mat[, only_0 == 0, drop = FALSE]
	    }
	}

    linear.mat
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



#' Replicates fixest objects
#'
#' Simple function that replicates fixest objects while (optionally) computing different standard-errors. Useful mostly in combination with \code{\link[fixest]{etable}} or \code{\link[fixest]{coefplot}}.
#'
#' @param x Either a fixest object, either a list of fixest objects created with \code{.l()}.
#' @param times Integer vector giving the number of repetitions of the vector of elements. By default \code{times = 1}. It must be either of length 1, either of the same length as the argument \code{x}.
#' @param each Integer scalar indicating the repetition of each element. Default is 1.
#' @param cluster A list containing the types of standard-error to be computed, default is missing. If not missing, it must be of the same length as \code{times}, \code{each}, or the final vector. Note that if the arguments \code{times} and \code{each} are missing, then \code{times} becomes equal to the length of \code{cluster}. (Note that \code{cluster} accepts the character values \code{"standard"} or \code{"hetero"} to compute non-clustered SEs.)
#' @param ... In \code{.l()}: \code{fixest} objects. In \code{rep()}: not currently used.
#'
#' @details
#' To apply \code{rep.fixest} on a list of fixest objects, it is absolutely necessary to use \code{.l()} and not \code{list()}.
#'
#' @return
#' Returns a list of the appropriate length. Each element of the list is a fixest object.
#'
#' @examples
#'
#' # Let's show results with different standard-errors
#'
#' est = feols(Ozone ~ Solar.R + Wind + Temp, data = airquality)
#'
#' my_cluster = list("Month", "Day", ~ Day + Month)
#'
#' etable(rep(est, cluster = my_cluster))
#'
#' coefplot(rep(est, cluster = my_cluster), drop = "Int")
#'
#' #
#' # To rep multiple objects, you need to use .l()
#' #
#'
#' est_bis = feols(Ozone ~ Solar.R + Wind + Temp | Month, airquality)
#'
#' etable(rep(.l(est, est_bis), cluster = my_cluster))
#'
#' # using each
#' etable(rep(.l(est, est_bis), each = 3, cluster = my_cluster))
#'
#'
rep.fixest = function(x, times = 1, each = 1, cluster, ...){
    # each is applied first, then times
    # x can be either a list of fixest objects, either a fixest object

    check_arg(x, "class(fixest, fixest_list) mbt")
    check_arg(times, "integer scalar GE{1} | integer vector no na GE{0}")
    check_arg(each, "integer scalar GE{1} | logical scalar")
    check_arg(cluster, "class(list)")

    validate_dots(suggest_args = c("times", "each"), stop = TRUE)

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

    IS_MULTI_CLUST = !missing(cluster)
    if(IS_MULTI_CLUST){
        n_clu = length(cluster)

        if(times == 1 && each == 1){
            if(isTRUE(each)){
                each = n_clu
            } else {
                times = n_clu
            }

        }
    }

    res_int = rep(1:n, times = times, each = each)
    n_res = length(res_int)

    if(IS_MULTI_CLUST){
        # Checking and expanding

        cluster_mapping = 1:n_res
        if(times == 1){
            if(n_clu != each && n_clu != n_res){
                stop("In rep, the argument 'cluster' (currently of length ", n_clu, ") must be a list either of length ", each, " or of length ", n_res, ".")
            }

            if(n_clu == each) cluster_mapping = rep(1:each, times = n)

        } else if(each == 1){
            if(n_clu != times && n_clu != n_res){
                stop("In rep, the argument 'cluster' (currently of length ", n_clu, ") must be a list either of length ", times, " or of length ", n_res, ".")
            }

            if(n_clu == times) cluster_mapping = rep(1:n_clu, each = n)

        } else {
            if(n_clu != n_res){
                stop("In rep, the argument 'cluster' (currently of length ", n_clu, ") must be a list either of length ", n_res, ".")
            }
        }

        se_all = vector("list", n_clu)
        for(m in 1:n_clu){
            if(identical(cluster[[m]], "standard")){
                se_all[[m]] = "standard"
            } else if(identical(cluster[[m]], "hetero")){
                se_all[[m]] = "hetero"
            }
        }

    }

    res = vector("list", length(res_int))

    if(IS_MULTI_CLUST){
        for(i in 1:n_res){
            if(IS_LIST){
                res[[i]] = summary(x[[res_int[i]]], se = se_all[[cluster_mapping[i]]], cluster = cluster[[cluster_mapping[i]]])
            } else {
                res[[i]] = summary(x, se = se_all[[cluster_mapping[i]]], cluster = cluster[[cluster_mapping[i]]])
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
rep.fixest_list = function(x, times = 1, each = 1, cluster, ...){
    rep.fixest(x, times = times, each = each, cluster = cluster, ...)
}

#' @rdname rep.fixest
.l = function(...){

    check_arg(..., "mbt class(fixest)")

    res = list(...)
    class(res) = "fixest_list"
    res
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
#'
#' # Back to default which is TRUE
#' getFixest_notes()
#'
setFixest_notes = function(x){

	if(missing(x) || length(x) != 1 || !is.logical(x) || is.na(x)){
		stop("Argument 'x' must be equal to TRUE or FALSE.")
	}

	options("fixest_notes" = x)
}

#' @rdname setFixest_notes
"getFixest_notes"

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
setFixest_nthreads = function(nthreads){
	# By default, we leave 2 nthreads (never use all)

    max_CRAN = as.numeric(Sys.getenv("OMP_THREAD_LIMIT"))
    max_CRAN[is.na(max_CRAN)] = 1000

	max_threads = min(cpp_get_nb_threads(), 1000, max_CRAN) # we cap at 1k nthreads

	if(missing(nthreads) || is.null(nthreads)){
		# New default => 50% of all available threads (usually equiv to the nber of procs)
		nthreads = check_set_nthreads(0.5)
	}

	nthreads = check_set_nthreads(nthreads)

	options("fixest_nthreads" = nthreads)

	invisible()
}

#' @rdname setFixest_nthreads
"getFixest_nthreads"

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
"getFixest_dict"

getFixest_dict = function(){

    x = getOption("fixest_dict")
    if(length(x) > 0){
        if(!is.character(x) || !isVector(x) || anyNA(x)){
            stop("The value of getOption(\"fixest_dict\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value. ")
        }
    }

    x
}

#' Sets/gets what \code{print} does to \code{fixest} estimations
#'
#' Sets/gets the default behavior of the \code{print} method for non-summary \code{fixest} estimations. Default is to display the coefficients table but it can be changed to displaying only the coefficients.
#'
#' @param x Either \code{"table"} or \code{"coef"}.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' res = feols(Sepal.Length ~ Sepal.Width + Petal.Length, iris)
#' # default is coef. table:
#' res
#' # can be changed to only the coefficients:
#' print(res, type = "coef")
#' setFixest_print.type("coef")
#' res # only the coefs
#'
#' # back to default
#' setFixest_print.type("table")
#'
setFixest_print.type = function(x){

    if(length(x) != 1 || !is.character(x) || is.na(x)){
        stop("Argument 'x' must be equal to 'table' or 'coef'.")
    }

    type = try(match.arg(x, c("coef", "table")), silent = TRUE)
    if("try-error" %in% class(type)){
        stop("Argument 'x' must be equal to 'table' or 'coef'.")
    }

    options("fixest_print.type" = type)
}

#' @rdname setFixest_print.type
"getFixest_print.type"

getFixest_print.type = function(){

    x = getOption("fixest_print.type")
    if(length(x) != 1 || !is.character(x) || is.na(x) || !x %in% c("coef", "table")){
        stop("The value of getOption(\"fixest_print.type\") is currently not legal. Please use function setFixest_print.type to set it to an appropriate value. ")
    }

    x
}



#' Type of degree of freedom in fixest summary
#'
#' Provides how the degrees of freedom should be calculated in \code{\link[fixest]{vcov.fixest}}/\code{\link[fixest]{summary.fixest}}.
#'
#' @param adj Logical scalar, defaults to \code{TRUE}. Whether to apply a small sample adjustment of the form \code{(n - 1) / (n - K)}, with \code{K} the number of estimated parameters. If \code{FALSE}, then no adjustment is made.
#' @param fixef.K Character scalar equal to \code{"nested"} (default), \code{"none"} or \code{"full"}. In the small sample adjustment, how to account for the fixed-effects parameters. If \code{"none"}, the fixed-effects parameters are discarded, meaning the number of parameters (\code{K}) is only equal to the number of variables. If \code{"full"}, then the number of parameters is equal to the number of variables plus the number of fixed-effects. Finally, if \code{"nested"}, then the number of parameters is equal to the number of variables plus the number of fixed-effects that *are not* nested in the clusters used to cluster the standard-errors.
#' @param fixef.force_exact Logical, default is \code{FALSE}. If there are 2 or more fixed-effects, these fixed-effects they can be irregular, meaning they can provide the same information. If so, the "real" number of parameters should be lower than the total number of fixed-effects. If \code{fixef.force_exact = TRUE}, then \code{\link[fixest]{fixef.fixest}} is first run to determine the exact number of parameters among the fixed-effects. Mostly, panels of the type individual-firm require \code{fixef.force_exact = TRUE} (but it adds computational costs).
#' @param cluster.adj Logical scalar, default is \code{TRUE}. How to make the small sample correction when clustering the standard-errors? If \code{TRUE} a \code{G/(G-1)} correction is performed with \code{G} the number of cluster values.
#' @param cluster.df Either "conventional" or "min" (default). Only relevant when the variance-covariance matrix is two-way clustered (or higher). It governs how the small sample adjustment for the clusters is to be performed. [Sorry for the jargon that follows.] By default a unique adjustment is made, of the form G_min/(G_min-1) with G_min the smallest G_i. If \code{cluster.df="conventional"} then the i-th "sandwich" matrix is adjusted with G_i/(G_i-1) with G_i the number of unique clusters.
#' @param t.df Either "conventional" or "min" (default). Only relevant when the variance-covariance matrix is clustered. It governs how the p-values should be computed. By default, the degrees of freedom of the Student t distribution is equal to the minimum size of the clusters with which the VCOV has been clustered. If \code{t.df="conventional"}, then the degrees of freedom of the Student t distribution is equal to the number of observations minus the number of estimated variables.
#'
#' @details
#'
#' The following vignette: \href{https://cran.r-project.org/package=fixest/vignettes/standard_errors.html}{On standard-errors}, describes in details how the standard-errors are computed in \code{fixest} and how you can replicate standard-errors from other software.
#'
#' @return
#' It returns a \code{dof.type} object.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}
#'
#' @examples
#'
#' #
#' # Equivalence with lm/glm standard-errors
#' #
#'
#' # LM
#' # In the absence of fixed-effects,
#' # by default, the standard-errors are computed in the same way
#'
#' res = feols(Petal.Length ~ Petal.Width + Species, iris)
#' res_lm = lm(Petal.Length ~ Petal.Width + Species, iris)
#' vcov(res) / vcov(res_lm)
#'
#' # GLM
#' # By default, there is no small sample adjustment in glm, as opposed to feglm.
#' # To get the same SEs, we need to use dof(adj = FALSE)
#'
#' res_pois = fepois(round(Petal.Length) ~ Petal.Width + Species, iris)
#' res_glm = glm(round(Petal.Length) ~ Petal.Width + Species, iris, family = poisson())
#' vcov(res_pois, dof = dof(adj = FALSE)) / vcov(res_glm)
#'
#' # Same example with the Gamma
#' res_gamma = feglm(round(Petal.Length) ~ Petal.Width + Species, iris, family = Gamma())
#' res_glm_gamma = glm(round(Petal.Length) ~ Petal.Width + Species, iris, family = Gamma())
#' vcov(res_gamma, dof = dof(adj = FALSE)) / vcov(res_glm_gamma)
#'
#' #
#' # Fixed-effects corrections
#' #
#'
#' # We create "irregular" FEs
#' base = data.frame(x = rnorm(10))
#' base$y = base$x + rnorm(10)
#' base$fe1 = rep(1:3, c(4, 3, 3))
#' base$fe2 = rep(1:5, each = 2)
#'
#' est = feols(y ~ x | fe1 + fe2, base)
#'
#' # fe1: 3 FEs
#' # fe2: 5 FEs
#'
#' #
#' # Clustered standard-errors: by fe1
#' #
#'
#' # Default: fixef.K = "nested"
#' #  => adjustment K = 1 + 5 (i.e. x + fe2)
#' summary(est)
#' attributes(vcov(est, attr = TRUE))[c("dof.type", "dof.K")]
#'
#'
#' # fixef.K = FALSE
#' #  => adjustment K = 1 (i.e. only x)
#' summary(est, dof = dof(fixef.K = "none"))
#' attr(vcov(est, dof = dof(fixef.K = "none"), attr = TRUE), "dof.K")
#'
#'
#' # fixef.K = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 1 (i.e. x + fe1 + fe2 - 1 restriction)
#' summary(est, dof = dof(fixef.K = "full"))
#' attr(vcov(est, dof = dof(fixef.K = "full"), attr = TRUE), "dof.K")
#'
#'
#' # fixef.K = TRUE & fixef.force_exact = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 2 (i.e. x + fe1 + fe2 - 2 restrictions)
#' summary(est, dof = dof(fixef.K = "full", fixef.force_exact = TRUE))
#' attr(vcov(est, dof = dof(fixef.K = "full", fixef.force_exact = TRUE), attr = TRUE), "dof.K")
#'
#' # There are two restrictions:
#' attr(fixef(est), "references")
#'
#' #
#' # To permanently set the default dof:
#' #
#'
#' # eg no small sample adjustment:
#' setFixest_dof(dof(adj = FALSE))
#'
#' # Factory default
#' setFixest_dof()
#'
dof = function(adj = TRUE, fixef.K = "nested", cluster.adj = TRUE, cluster.df = "min", t.df = "min", fixef.force_exact = FALSE){

    check_arg_plus(adj, "loose logical scalar conv")
    check_arg_plus(fixef.K, "match(none, full, nested)")
    check_arg_plus(cluster.df, "match(conventional, min)")
    check_arg_plus(t.df, "match(conventional, min)")
    check_arg(fixef.force_exact, cluster.adj, "logical scalar")

    res = list(adj = adj, fixef.K = fixef.K, cluster.adj = cluster.adj, cluster.df = cluster.df, t.df = t.df, fixef.force_exact = fixef.force_exact)
    class(res) = "dof.type"
    res
}


#' @rdname dof
#'
#' @param dof.type An object of class \code{dof.type} obtained with the function \code{\link[fixest]{dof}}.
setFixest_dof = function(dof.type = dof()){

    if(!"dof.type" %in% class(dof.type)){
        stop("The argument 'dof.type' must be an object created by the function dof().")
    }

    options("fixest_dof" = dof.type)
}

#' @rdname dof
"getFixest_dof"

getFixest_dof = function(){

    dof.type = getOption("fixest_dof")
    if(!"dof.type" %in% class(dof.type)){
        stop("The value of getOption(\"fixest_dof\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value.")
    }

    dof.type
}


#' Sets the default type of standard errors to be used
#'
#' This functions defines or extracts the default type of standard-errors to computed in \code{fixest} \code{\link[fixest:summary.fixest]{summary}}, and \code{\link[fixest:vcov.fixest]{vcov}}.
#'
#' @param no_FE Character scalar equal to either: \code{"standard"} (default), or \code{"hetero"}. The type of standard-errors to use by default for estimations without fixed-effects.
#' @param one_FE Character scalar equal to either: \code{"standard"}, \code{"hetero"}, or \code{"cluster"} (default). The type of standard-errors to use by default for estimations with \emph{one} fixed-effect.
#' @param two_FE Character scalar equal to either: \code{"standard"}, \code{"hetero"}, \code{"cluster"}, or \code{"twoway"} (default). The type of standard-errors to use by default for estimations with \emph{two or more} fixed-effects.
#'
#' @return
#' The function \code{getFixest_se()} returns a list with three elements containing the default for estimations i) wihtout, ii) with one, or iii) with two or more fixed-effects.
#'
#' @examples
#'
#' # By default:
#' # - no fixed-effect (FE): standard
#' # - one or more FEs: cluster
#'
#' data(base_did)
#' est_no_FE  = feols(y ~ x1, base_did)
#' est_one_FE = feols(y ~ x1 | id, base_did)
#' est_two_FE = feols(y ~ x1 | id + period, base_did)
#'
#' etable(est_no_FE, est_one_FE, est_two_FE)
#'
#' # Changing the default standard-errors
#' setFixest_se(no_FE = "hetero", one_FE = "standard", two_FE = "twoway")
#' etable(est_no_FE, est_one_FE, est_two_FE)
#'
#' # Reseting the defaults
#' setFixest_se()
#'
#'
setFixest_se = function(no_FE = "standard", one_FE = "cluster", two_FE = "cluster"){

    check_arg_plus(no_FE,  "match(standard, hetero)")
    check_arg_plus(one_FE, "match(standard, hetero, cluster)")
    check_arg_plus(two_FE, "match(standard, hetero, cluster, twoway)")

    options(fixest_se = list(no_FE = no_FE, one_FE = one_FE, two_FE = two_FE))


}

#' @rdname setFixest_se
getFixest_se = function(){

    se_default = getOption("fixest_se")

    if(is.null(se_default)){
        se_default = list(no_FE = "standard", one_FE = "cluster", two_FE = "cluster")
        options(fixest_se = se_default)
        return(se_default)
    }

    if(!is.list(se_default) || !all(unlist(se_default) %in% c("standard", "hetero", "cluster", "twoway"))){
        stop("The value of getOption(\"se_default\") is currently not legal. Please use function setFixest_se to set it to an appropriate value.")
    }

    se_default
}



#' Sets/gets and expands formula macros
#'
#' You can set formula macros globally with \code{setFixest_fml}. These macros can then be used in \code{fixest} estimations or when using the function \code{\link[fixest:setFixest_fml]{xpd}}.
#'
#' @param fml A formula containing macros variables. The macro variables can be set globally using \code{setFixest_fml}, or can be defined in \code{...}.
#' @param ... Definition of the macro variables. Each argument name corresponds to the name of the macro variable. It is required that each macro variable name starts with two dots (e.g. \code{..ctrl}). The value of each argument must be a one-sided formula or a character vector, it is the definition of the macro variable. Example of a valid call: \code{setFixest_fml(..ctrl = ~ var1 + var2)}. In the function \code{xpd}, the default macro variables are taken from \code{getFixest_fml}, any variable in \code{...} will replace these values.
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
#' a <- feols(Ozone ~ Wind + ..ctrl, airquality)
#' b <- feols(Ozone ~ Wind + ..ctrl_long, airquality)
#' etable(a, b, keep = "Int|Win")
#'
#' #
#' # You can use xpd for stepwise estimations
#' #
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
#' data(longley)
#' setFixest_fml(..many_vars = grep("GNP|ployed", names(longley), value = TRUE))
#' feols(Armed.Forces ~Population + ..many_vars, longley)
#'
#' #
#' # You can also put numbers in macros
#' #
#'
#' res_all = list()
#' for(p in 1:3){
#'   res_all[[p]] = feols(xpd(Ozone ~ Wind + poly(Temp, ..p), ..p = p), airquality)
#' }
#'
#' etable(res_all)
#'
#'
setFixest_fml = function(..., reset = FALSE){

    check_arg(reset, "logical scalar")

    fml_macro = parse_macros(..., reset = reset)

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


#' Sample data for difference in difference
#'
#' This data has been generated to illustrate the use of difference in difference functions in package \pkg{fixest}. This is a balanced panel of 104 individuals and 10 periods. About half the individuals are treated, the treatment having a positive effect on the dependent variable \code{y} after the 5th period. The effect of the treatment on \code{y} is gradual.
#'
#' @usage
#' data(base_did)
#'
#' @format
#' \code{base_did} is a data frame with 1,040 observations and 6 variables named \code{y}, \code{x1}, \code{id}, \code{period}, \code{post} and \code{treat}.
#'
#' \itemize{
#' \item{y: The dependent variable affected by the treatment.}
#' \item{x1: An explanatory variable.}
#' \item{id: Identifier of the individual.}
#' \item{period: From 1 to 10}
#' \item{post: Indicator taking value 1 if the period is strictly greater than 5, 0 otherwise.}
#' \item{treat: Indicator taking value 1 if the individual is treated, 0 otherwise.}
#'
#' }
#'
#' @source
#' This data has been generated from \pkg{R}.
#'
#' @seealso
#' The DiD functions of the package \pkg{fixest} are \code{\link[fixest]{did_estimate_yearly_effects}} and \code{\link[fixest]{did_plot_yearly_effects}}.
#'
#'
#'
"base_did"































































