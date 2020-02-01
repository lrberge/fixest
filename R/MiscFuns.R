
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
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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
#' \donttest{
#' # To permanently display coef. only, use setFixest_print.type:
#' setFixest_print.type("coef")
#' est_pois
#' # back to default:
#' setFixest_print.type("table")
#' }
#'
#'
print.fixest <- function(x, n, type = getFixest_print.type(), ...){


    # checking the arguments
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  dots_args = c("se", "cluster", "dof", "exact_dof", "forceCovariance", "keepBounded"),
                                  suggest_args = c("n", "type", "se", "cluster"))
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

    # The objects from the estimation and the summary are identical, except regarding the vcov
	fromSummary = "cov.scaled" %in% names(x)

	# if NOT from summary, we consider the argument 'type'
	if(!fromSummary){
	    # checking argument type
	    if(length(type) != 1 || !is.character(type) || is.na(type)){
	        stop("Argument 'type' must be equal to 'table' or 'coef'.")
	    }
	    type = try(match.arg(type, c("coef", "table")), silent = TRUE)
	    if("try-error" %in% class(type)){
	        stop("Argument 'type' must be equal to 'table' or 'coef'.")
	    }

	    if(type == "coef"){
	        print(coef(x))
	        return(invisible())
	    }
	}

	isNegbin = x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && x$family=="negbin")

	x = summary(x, fromPrint = TRUE, ...)

	msgRemaining = ""
	nb_coef = length(coef(x)) - isNegbin
	if(missing(n)){
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

	} else if(!length(n) == 1 || !is.numeric(n) || n<=0){
		stop("Argument 'n' must be a single positive integer.")
	} else if(n < nb_coef){
		msgRemaining = paste0("... ", nb_coef - n, " coefficients remaining\n")
	}

	if(isFALSE(x$convStatus)){
	    last_warn = getOption("fixest_last_warning")
	    if(is.null(last_warn) || (proc.time() - last_warn)[3] > 1){
	        warning("The optimization algorithm did not converge, the results are not reliable. (", x$message, ")", call. = FALSE)
	    }

	}

	coeftable = x$coeftable

	# The type of SE
	se.type = attr(coeftable, "type")

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
	        cat("Fixed-effects: ", paste0(fixef_vars, ": ", x$fixef_sizes[fixef_vars], collapse=",  "), "\n", sep="")
	    }

	    cat("Varying slopes: ", paste0(terms_full$slope_vars, " (", terms_full$slope_fe, ": ", x$fixef_sizes[terms_full$slope_fe], ")", collapse=",  "), "\n", sep="")

	} else {
	    if(!is.null(x$fixef_sizes)) cat("Fixed-effects: ", paste0(x$fixef_vars, ": ", x$fixef_sizes, collapse=",  "), "\n", sep="")
	}


	if(is.null(x$onlyFixef)){

		cat("Standard-errors:", se.type, "\n")

		# The matrix of coefficients
		if(isNegbin){
			if(nrow(coeftable) == 2){
				new_table = coeftable[1, , FALSE]
			} else {
				new_table = coeftable[-nrow(coeftable), ]
			}

			myPrintCoefTable(head(new_table, n), lastLine = msgRemaining)

			theta = coeftable[".theta", 1]
			noDispInfo = ifelse(theta > 1000, "(theta >> 0, no sign of overdispersion, you may consider a Poisson model)", "")
			cat("Over-dispersion parameter: theta =", theta, noDispInfo, "\n")
		} else {
			myPrintCoefTable(head(coeftable, n), lastLine = msgRemaining)
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


	if(!is.null(x$convStatus) && !x$convStatus && is.null(x$onlyFixef)){
		cat("# Evaluations:", x$iterations, "--", x$message, "\n")
	}

}

##

#' Summary of a \code{fixest} object. Computes different types of standard errors.
#'
#' This function is similar to \code{print.fixest}. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
#'
#' @method summary fixest
#'
#' @param se Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{White}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: \code{se = "cluster"}, otherwise \code{se = "standard"}. Note that this argument can be implicitly deduced from the argument \code{cluster}.
#' @param cluster Tells how to cluster the standard-errors (if clustering is requested). Can be either a list of vectors, a character vector of variable names, a formula or an integer vector. Assume we want to perform 2-way clustering over \code{var1} and \code{var2} contained in the data.frame \code{base} used for the estimation. All the following \code{cluster} arguments are valid and do the same thing: \code{cluster = base[, c("var1, "var2")]}, \code{cluster = c("var1, "var2")}, \code{cluster = ~var1+var2}. If the two variables were used as clusters in the estimation, you could further use \code{cluster = 1:2} or leave it blank with \code{se = "twoway"} (assuming \code{var1} [resp. \code{var2}] was the 1st [res. 2nd] cluster).
#' @param object A \code{fixest} object. Obtained using the functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#' @param dof An object of class \code{dof.type} obtained with the function \code{\link[fixest]{dof}}. Represent how the degree of freedom correction should be done. Defaults to \code{dof(fixef="nested", exact=FALSE, cluster = TRUE)}. See the help of the function \code{\link[fixest]{dof}} for details.
#' @param forceCovariance (Advanced users.) Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option to force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded (Advanced users -- feNmlm with non-linear part and bounded coefficients only.) Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param ... Not currently used.
#'
#' @return
#' It returns a \code{fixest} object with:
#' \item{cov.scaled}{The new variance-covariance matrix (computed according to the argument \code{se}).}
#' \item{se}{The new standard-errors (computed according to the argument \code{se}).}
#' \item{coeftable}{The table of coefficients with the new standard errors.}
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 cluster effects)
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Comparing different types of standard errors
#' sum_white    = summary(est_pois, se = "white")
#' sum_oneway   = summary(est_pois, se = "cluster")
#' sum_twoway   = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' esttable(sum_white, sum_oneway, sum_twoway, sum_threeway)
#'
#' # Alternative ways to cluster the SE:
#' \donttest{
#' # two-way clustering: Destination and Product
#' # (Note that arg. se = "twoway" is implicitly deduced from the argument cluster)
#' summary(est_pois, cluster = c("Destination", "Product"))
#' summary(est_pois, cluster = trade[, c("Destination", "Product")])
#' summary(est_pois, cluster = list(trade$Destination, trade$Product))
#' summary(est_pois, cluster = ~Destination+Product)
#' # Since Destination and Product are used as fixed-effects, you can also use:
#' summary(est_pois, cluster = 2:3)
#' }
#'
#'
summary.fixest <- function(object, se, cluster, dof = getFixest_dof(), forceCovariance = FALSE, keepBounded = FALSE, ...){
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
	    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
	                                  suggest_args = c("se", "cluster", "dof"))
	    if(any_invalid){
	        warning(attr(any_invalid, "msg"))
	    }
	}


	if(is.null(dots$nframes_up)){
		nframes_up = 1 + !is.null(dots$fromPrint)
	} else {
		nframes_up = dots$nframes_up + 1
	}

	# The new VCOV
	vcov = vcov(object, se=se, cluster=cluster, dof=dof, forceCovariance = forceCovariance, keepBounded = keepBounded, nframes_up = nframes_up)

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
	pvalue <- 2*pt(-abs(zvalue), object$nobs - object$nparams)

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


#' Estimations table (export the results of multiples estimations to a DF or to Latex)
#'
#' Aggregates the results of multiple estimations and displays them in the form of either a Latex table or a \code{data.frame}.
#'
#' @inheritParams summary.fixest
#'
#' @param ... Used to capture different \code{fixest} estimation objects (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). Note that any other type of element is discarded. Note that you can give a list of \code{fixest} objects.
#' @param digits Integer, default is 4. The number of digits to be displayed.
#' @param tex Logical: whether the results should be a data.frame or a Latex table. By default, this argument is \code{TRUE} if the argument \code{file} (used for exportation) is not missing; it is equal to \code{FALSE} otherwise.
#' @param fitstat A character vector or a one sided formula. A vector listing which fit statistics to display. The valid types are 'll', 'aic', 'bic' and r2 types like 'r2', 'pr2', 'war2', etc (see all valid types in \code{\link[fixest]{r2}}). The default value depends on the models to display. Example of use: \code{fitstat=c('sq.cor', 'ar2', 'war2')}, or \code{fitstat=~sq.cor+ar2+war2} using a formula.
#' @param title (Tex only.) Character scalar. The title of the Latex table.
#' @param sdBelow (Tex only.) Logical, default is \code{TRUE}. Should the standard-errors be displayed below the coefficients?
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be discarded. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{drop = "x[[:digit:]]{2}"}. If the first character is an exclamation mark, the effect is reversed (e.g. drop = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is dropped). See details.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions. This argument is applied post aliasing (see argument \code{dict}). Example: you have the following variables: \code{month1} to \code{month6}, then \code{x1} to \code{x5}, then \code{year1} to \code{year6}. If you want to display first the x's, then the years, then the months you could use: \code{order = c("x", "year")}. If the first character is an exclamation mark, the effect is reversed (e.g. order = "!Intercept" means: every variable that does not contain \dQuote{Intercept} goes first).  See details.
#' @param dict (Tex only.) A named character vector. It changes the original variable names to the ones contained in the \code{dict}. E.g. to change the variables named \code{a} and \code{b3} to (resp.) \dQuote{$log(a)$} and to \dQuote{$bonus^3$}, use \code{dict=c(a="$log(a)$",b3="$bonus^3$")}. By default it is equal to \code{getFixest_dict()}, a default dictionary which can be set with \code{\link[fixest]{setFixest_dict}}.
#' @param file A character scalar. If provided, the Latex (or data frame) table will be saved in a file whose path is \code{file}. If you provide this argument, then a Latex table will be exported, to export a regular \code{data.frame}, use argument \code{tex = FALSE}.
#' @param replace Logical, default is \code{FALSE}. Only used if option \code{file} is used. Should the exported table be written in a new file that replaces any existing file?
#' @param convergence Logical, default is missing. Should the convergence state of the algorithm be displayed? By default, convergence information is displayed if at least one model did not converge.
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)} for a Latex table and \code{c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)} for a data.frame (to conform with R's default).
#' @param label (Tex only.) Character scalar. The label of the Latex table.
#' @param subtitles Character vector of the same length as the number of models to be displayed. If provided, subtitles are added underneath the dependent variable name.
#' @param fixef_sizes (Tex only.) Logical, default is \code{FALSE}. If \code{TRUE} and fixed-effects were used in the models, then the number of "individuals" per fixed-effect dimension is also displayed.
#' @param yesNoFixef (Tex only.) A character vector of length 1 or 2. Default is \code{c("Yes", "No")}. This is the message displayed when a given cluster is (or is not) included in a regression. If \code{yesNoFixef} is of length 1, then the second element is the empty string.
#' @param family Logical, default is missing. Whether to display the families of the models. By default this line is displayed when at least two models are from different families.
#' @param keepFactors Logical, default is \code{TRUE}. If \code{FALSE}, then factor variables are displayed as fixed-effects and no coefficient is shown.
#' @param powerBelow (Tex only.) Integer, default is -5. A coefficient whose value is below \code{10**(powerBelow+1)} is written with a power in Latex. For example \code{0.0000456} would be written \code{4.56$\\times 10^{-5}$} by default. Setting \code{powerBelow = -6} would lead to \code{0.00004} in Latex.
#' @param interaction.combine (Tex only.) Character scalar, defaults to \code{" $\\times$ "}. When the estimation contains interactions, then the variables names (after aliasing) are combined with this argument. For example: if \code{dict = c(x1="Wind", x2="Rain")} and you have the following interaction \code{x1:x2}, then it will be renamed (by default) \code{Wind $\\times$ Rain} -- using \code{interaction.combine = "*"} would lead to \code{Wind*Rain}.
#' @param depvar (Data frame only.) Logical, default is missing. Whether a first line containing the dependent variables should be shown. By default, the dependent variables are shown only if they differ across models or if the argumen \code{file} is not missing.
#'
#' @details
#' The function \code{esttex} is equivalent to the function \code{etable} with argument \code{tex = TRUE}.
#'
#' The function \code{esttable} is equivalent to the function \code{etable} with argument \code{tex = FALSE}.
#'
#' @section Arguments drop and order:
#' The arguments \code{drop} and \code{order} use regular expressions. If you are aware of regular expressions, I urge you to learn it, since it is an extremely powerful way to manipulate character strings (and it exists across most programming languages).
#'
#' For example drop = "Wind" would drop any variable whose name contains "Wind". Note that variables such as "Temp:Wind" or "StrongWind" do contain "Wind", so would be dropped. To drop only the variable named "Wind", you need to use \code{drop = "^Wind$"} (with "^" meaning beginning, resp. "$" meaning end, of the string => this is the language of regular expressions).
#'
#' Although you can combine several regular expressions in a single character string using pipes, \code{drop} also accepts a vector of regular expressions.
#'
#' You can use the special character "!" (exclamation mark) to reverse the effect of the regular expression (this feature is specific to this fonction). For example \code{drop = "!Wind"} would drop any variable that does not contain "Wind".
#'
#' The argument \code{order} takes in a vector of regular expressions, the order will follow the elments of this vector. The vector gives a list of priorities, on the left the elements with highest priority. For example, order = c("Wind", "!Inter", "!Temp") would give highest priorities to the variables containing "Wind" (which would then appear first), second highest priority is the variables not containing "Inter", last, with lowest priority, the variables not containing "Temp". If you had the following variables: (Intercept), Temp:Wind, Wind, Temp you would end up with the following order: Wind, Temp:Wind, Temp, (Intercept).
#'
#' @return
#' If \code{tex = TRUE}, the lines composing the Latex table are returned invisibly while the table is directly prompted on the console.
#'
#' If \code{tex = FALSE}, the data.frame is directly returned. If the argument \code{file} is not missing, the \code{data.frame} is returned invisibly.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' aq = airquality
#' aq$Month = factor(aq$Month)
#'
#' est1 = feols(Ozone ~ Month / Wind + Temp, data = aq)
#' est2 = feols(Ozone ~ Wind + Temp | Month, data = aq)
#'
#' # Displaying the two results in a single table
#' etable(est1, est2)
#'
#' # drop: dropping the non interaction terms (see regexp help)
#' etable(est1, est2, drop = "^[[:alnum:]]+$")
#' # drop interactions
#' etable(est1, est2, drop = ":")
#' # keep only interactions ("!" reverses the effect)
#' etable(est1, est2, drop = "!:")
#'
#' # order: Wind variable first, intercept last
#' etable(est1, est2, order = c("Wind", "Month"))
#' etable(est1, est2, order = c("^Wind", "Wind", "Month"))
#' # Interactions, then Intercept, last ("!" reverses the effect)
#' etable(est1, est2, order = c("!Int", "!:"))
#'
#' # signifCode
#' etable(est1, est2, signifCode = c(" A"=0.01, " B"=0.05, " C"=0.1,
#'                                   " D"=0.15, " F"=1))
#'
#' # fitstat
#' etable(est1, est2, fitstat = ~r2+ar2+apr2+war2)
#'
#' # Adding a dictionnary (Tex only)
#' dict = c(Month5="May", Month6="Jun", Month7="Jul", Month8="Aug", Month9="Sep")
#' etable(est1, est2, dict = dict, tex = TRUE)
#'
etable = function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = getFixest_dof(), cluster, digits=4, tex, fitstat, title, sdBelow = TRUE, drop, order, dict = getFixest_dict(), file, replace=FALSE, convergence, signifCode, label, subtitles, fixef_sizes = FALSE, yesNoFixef = c("Yes", "No"), keepFactors = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ ", depvar){

    #
    # Checking the arguments
    #

    # Need to check for the presence of the se
    useSummary = TRUE
    if(missing(se) && missing(cluster)){
        useSummary = FALSE
    }

    if(!missing(se)){
        se = match.arg(se)
    } else {
        se = NULL
    }

    # check the dictionnary
    if(!is.null(dict) && (!is.character(dict) || is.null(names(dict)))){
        stop("The argument 'dict' must be a named character vector.")
    }

    # The depvar
    if(missing(depvar) && !missing(file)){
        depvar = TRUE
    }


    check_arg(digits, "singleIntegerGE1")
    check_arg(title, "singleCharacter")
    check_arg(sdBelow, "singleLogical")
    check_arg(drop, "characterVector", "The arg. 'drop' must be a vector of regular expressions (see help(regex)). REASON")
    check_arg(order, "characterVector", "The arg. 'order' must be a vector of regular expressions (see help(regex)). REASON")
    check_arg(file, "singleCharacter")
    check_arg(replace, "singleLogical")
    check_arg(convergence, "singleLogical")
    check_arg(signifCode, "numericVectorGE0LE1")
    check_arg(label, "singleCharacter")
    check_arg(subtitles, "characterVector")
    check_arg(fixef_sizes, "singleLogical")
    check_arg(yesNoFixef, "characterVector")
    check_arg(keepFactors, "singleLogical")
    check_arg(family, "singleLogical")
    check_arg(powerBelow, "singleIntegerLE-1")
    check_arg(interaction.combine, "singleCharacter")
    check_arg(tex, "singleLogical")
    check_arg(depvar, "singleLogical")

    if(missing(tex)){
        if(!missing(file)) {
            tex = TRUE
        } else {
            tex = FALSE
        }
    }

    # The signif codes
    if(missing(signifCode)){
        if(tex){
            signifCode = c("***"=0.01, "**"=0.05, "*"=0.10)
        } else {
            signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)
        }
    }

    # to get the model names
    dots_call = match.call(expand.dots = FALSE)[["..."]]

    info = results2formattedList(..., se=se, dof=dof, fitstat=fitstat, cluster=cluster, digits=digits, sdBelow=sdBelow, signifCode=signifCode, title=title, subtitles=subtitles, yesNoFixef=yesNoFixef, keepFactors=keepFactors, isTex = tex, useSummary=useSummary, dots_call=dots_call, powerBelow=powerBelow, dict=dict, interaction.combine=interaction.combine, convergence=convergence, show_family=family, drop=drop, file=file, order=order, label=label, fixef_sizes=fixef_sizes, show_depvar=depvar)

    if(tex){
        res = etable_internal_latex(info)
    } else {
        res = etable_internal_df(info)
    }

    if(!missnull(file)){
        sink(file = file, append = !replace)
        on.exit(sink())

        if(tex){
            cat(res, sep = "")
        } else {
            print(res)
        }

        return(invisible(res))
    } else {
        if(tex){
            cat(res, sep = "")
            return(invisible(res))
        } else {
            return(res)
        }
    }
}


#' @describeIn etable Exports the results of multiple \code{fixest} estimations in a Latex table.
esttex <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = getFixest_dof(), cluster, digits=4, fitstat, title, sdBelow=TRUE, drop, order, dict = getFixest_dict(), file, replace=FALSE, convergence, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, subtitles, fixef_sizes = FALSE, yesNoFixef = c("Yes", "No"), keepFactors = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ "){
	# drop: a vector of regular expressions
	# order: a vector of regular expressions
	# dict: a 'named' vector
	# file: a character string

    useSummary = TRUE
	if(missing(se) && missing(cluster)){
		useSummary = FALSE
	}

    if(!missing(se)){
		se = match.arg(se)
    } else {
	    se = NULL
    }

    # check the dictionnary
    if(!is.null(dict) && (!is.character(dict) || is.null(names(dict)))){
        stop("The argument 'dict' must be a named character vector.")
    }

    check_arg(interaction.combine, "singleCharacter")

	# to get the model names
	dots_call = match.call(expand.dots = FALSE)[["..."]]

	info = results2formattedList(..., se=se, dof=dof, fitstat=fitstat, cluster=cluster, digits=digits, sdBelow=sdBelow, signifCode=signifCode, subtitles=subtitles, yesNoFixef=yesNoFixef, keepFactors=keepFactors, isTex = TRUE, useSummary=useSummary, dots_call=dots_call, powerBelow=powerBelow, dict=dict, interaction.combine=interaction.combine)

    res = etable_internal_latex(info)

	if(!missnull(file)){
	    sink(file = file, append = !replace)
	    on.exit(sink())
	}

    cat(res, sep = "")
    return(invisible(res))

}

#' @describeIn etable Facility to display the results of multiple \code{fixest} estimations.
esttable <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = getFixest_dof(), cluster, depvar, drop, order, digits=4, fitstat, convergence, signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10), subtitles, keepFactors = FALSE, family){

	# Need to check for the presence of the se
    useSummary = TRUE
    if(missing(se) && missing(cluster)){
        useSummary = FALSE
    }

    if(!missing(se)){
        se = match.arg(se)
    } else {
        se = NULL
    }

	# to get the model names
	dots_call = match.call(expand.dots = FALSE)[["..."]]

	info = results2formattedList(..., se=se, dof = dof, cluster=cluster, digits=digits, signifCode=signifCode, subtitles=subtitles, keepFactors=keepFactors, useSummary=useSummary, dots_call=dots_call, fitstat=fitstat, yesNoFixef = c("Yes", "No"))

	res = etable_internal_df(info)

	return(res)
}




#' R2s of \code{fixest} models
#'
#' Reports different R2s for \code{fixest} estimations (e.g. \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}).
#'
#' @param x A \code{fixest} object, e.g. obtained with function \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}.
#' @param type A character vector representing the R2 to compute. The R2 codes are of the form: "wapr2" with letters "w" (within), "a" (adjusted) and "p" (pseudo) possibly missing. E.g. to get the regular R2: use \code{type = "r2"}, the within adjusted R2: use \code{type = "war2"}, the pseudo R2: use \code{type = "pr2"}, etc. Use \code{"sq.cor"} for the squared correlation. By default, all R2s are computed.
#'
#' @details
#' For R2s with no theoretical justification, like e.g. regular R2s for maximum likelihood models -- or within R2s for models without fixed-effects, NA is returned. The single measure to possibly compare all kinds of models is the squared correlation between the dependent variable and the expected predictor.
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
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Squared correlation:
#' r2(est_pois, "sq.cor")
#'
#' # "regular" r2:
#' r2(est_pois, "r2")
#'
#' # pseudo r2
#' r2(est_pois, "pr2")
#'
#' # within adjusted r2
#' r2(est_pois, "war2")
#'
#' # all four at once
#' r2(est_pois, c("sq.cor", "r2", "pr2", "war2"))
#'
r2 = function(x, type = "all"){
	# p: pseudo
	# w: within
	# a: adjusted
    # NOTA: wr2 not supported for feglm because fe_model incurs too much computational cost

	if(!"fixest" %in% class(x)){
		stop("Only 'fixest' objects are supported.")
	}

	if(!is.character(type) || !isVector(type)){
		stop("Argument 'type' must be a character vector (e.g. type = c(\"sq.cor\", \"r2\", \"pr2\")). (a: adjused, p: pseudo, w: within.)")
	}

	# type_allowed next => ("count", "acount") ?
    dict_names = c("sq.cor" = "Squared Correlation", "r2" = "R2", "ar2" = "Adjusted R2", "pr2" = "Pseudo R2", "apr2" = "Adjusted Pseudo R2", "wr2" = "Within R2", "war2" = "Within Adjusted R2", "wpr2" = "Within Pseudo R2", "wapr2" = "Within Adjusted Pseudo R2")
	type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2")
	if("all" %in% type){
		type_all = type_allowed
	} else {
		type_all = tolower(unique(type))
		pblm = setdiff(type_all, type_allowed)
		if(length(pblm) > 0){
			stop("The r2 type", enumerate_items(pblm, "s.is"), " not valid.")
		}
	}

	is_ols = x$method == "feols"
	isCluster = "fixef_vars" %in% names(x)
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

		if(grepl("w", myType) && !isCluster){
		    # within R2 not valid for models without FE
		    next
		}

		adj = grepl("a", myType)
		pseudo = grepl("p", myType)
		within = grepl("w", myType) && isCluster
		ifNullNA = function(x) ifelse(is.null(x), NA, x)
		if(within && isCluster){

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
		        new_fml = as.formula(paste0("y~1|", gsub(".+\\|", "", as.character(x$fml_full)[3])))

		        # The fact that weights = x[["weights"]] is on purpose -- don't touch it
		        # same for offset = x$offset -- otherwise error is thrown in fixest_env
		        res_fe = feglm(fml = new_fml, data = newdata, glm.tol = 1e-2, fixef.tol = 1e-3, family = x$family$family, weights = x[["weights"]], offset = x$offset)

		        x$ssr_fe_only = cpp_ssq(resid(res_fe))
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
				res[i] = 1 - cpp_ssq(resid(x)) / ssr_fe_only * (n - nb_fe) / (n - nb_fe - df_k)
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
				res[i] = 1 - drop(crossprod(resid(x))) / ssr_null * (n - 1) / (n - df_k)
			}

		}
	}

	names(res) = type_all

	res
}


#' Obtain various statistics from an estimation
#'
#' Set of functions to directly extract some commonly used statistics, like the p-value or the table of coefficients, from estimations. This was first implemented for \code{fixest} estimations, but has some support for other models.
#'
#' @param object An estimation. For example obtained from \code{\link[fixest]{feols}}.
#' @param se [Fixest specific.] Character scalar. Which kind of standard error should be computed: \dQuote{standard}, \dQuote{White}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}? By default if there are clusters in the estimation: \code{se = "cluster"}, otherwise \code{se = "standard"}. Note that this argument can be implicitly deduced from the argument \code{cluster}.
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

#' Summary method for cluster coefficients
#'
#' This function summarizes the main characteristics of the cluster coefficients. It shows the number of fixed-effects that have been set as references and the first elements of the fixed-effects.
#'
#' @method summary fixest.fixef
#'
#' @param object An object returned by the function \code{\link[fixest]{fixef.fixest}}.
#' @param n Positive integer, defaults to 5. The \code{n} first fixed-effects for each cluster are reported.
#' @param ... Not currently used.
#'
#' @return
#' It prints the number of fixed-effect coefficients per cluster, as well as the number of fixed-effects used as references for each cluster, and the mean and variance of the cluster coefficients. Finally it reports the first 5 elements of each cluster.
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
#' # printing some summary information on the cluster coefficients:
#' summary(fe_trade)
#'
#'
summary.fixest.fixef = function(object, n=5, ...){
	# This function shows some generic information on the clusters

    # checking arguments in dots
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = "n")
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
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
		mean_per_cluster = var_per_cluster = c()
		for(i in 1:Q){
			mean_per_cluster[i] = as.character(signif(mean(object[[i]]), 3))
			var_per_cluster[i] = as.character(signif(var(object[[i]]), 3))
		}
		res = as.data.frame(rbind(nb_per_cluster, nb_ref, mean_per_cluster, var_per_cluster))

		row_1 = paste0("Number of ", switch(info, "11" = "fixed-effects/slopes", "10"="fixed-effects", "1"="slopes"))

		rownames(res) = c(row_1, "Number of references", "Mean", "Variance")

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

	# We print the first 5 elements of each cluster
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
#' This function retrieves the fixed effects from a \code{fixest} estimation. It is useful only when there are one or more clusters.
#'
#' @inheritParams feNmlm
#'
#' @param object A \code{fixest} estimation (e.g. obtained using \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}).
#' @param notes Logical. Whether to display a note when the fixed-effects coefficients are not regular.
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
#' \code{\link[fixest]{plot.fixest.fixef}}. See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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
#' # The fixed-effects of the first cluster:
#' head(fe_trade$Origin)
#'
#' # Summary information:
#' summary(fe_trade)
#'
#' # Plotting them:
#' plot(fe_trade)
#'
fixef.fixest = function(object, notes = getFixest_notes(), ...){
	# object is a fixest object
	# This function retrieves the dummies

    # Checking the arguments
    any_invalid = check_dots_args(match.call(expand.dots = FALSE))
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
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

	isSlope = FALSE
	if(!is.null(object$fixef_terms)){
	    isSlope = TRUE
	    # This is an estimation with slopes
	    # we apply another method => we use the demeaning function

	    fixef_terms = object$fixef_terms
	    terms_full = extract_fe_slope(fixef_terms)
	    fixef_vars = terms_full$fixef_vars
	    slope_fe = terms_full$slope_fe
	    fe_all = terms_full$fe_all
	    slope_vars = terms_full$slope_vars
	    slope_terms = terms_full$slope_terms

	    # dictionary mapping fe var names to the ids of id_dummies_vect
	    dict_fe = 1:Q
	    names(dict_fe) = object$fixef_vars

	    # STEP 1: getting the variables with varying slopes (not saved in object to save space)
	    # We extract the slope variables
	    slope_vars_unik = unique(slope_vars)

	    # 2019-11-26: now the slope variables are directly in the objects
	    # # evaluation
	    # data = NULL
	    # try(data <- eval(object$call$data, parent.frame()), silent = TRUE)
	    # dataName = object$call$data
	    #
	    # if(is.null(data)){
	    #     # We try to provide an informative error message
	    #     stop("To get the coefficients for the variables with varying slopes, we need to fetch these variables (i.e. ", enumerate_items(slope_vars_unik), ") in the original dataset in the parent.frame -- but the data doesn't seem to be there anymore (btw it was ", deparse_long(dataName), ")")
	    # }
	    #
	    # data = as.data.frame(data)
	    #
	    # # we check the variables are there
	    # slope_vars_unik_raw = unique(sapply(slope_vars_unik, function(x) all.vars(parse(text = x))))
	    #
	    # if(any(!slope_vars_unik_raw %in% names(data))){
	    #     var_pblm = setdiff(slope_vars_unik_raw, names(data))
	    #     stop("To get the coefficients for the variables with varying slopes, we need to fetch these variables in the original dataset (", deparse_long(dataName), "). However, the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset any more.")
	    # }
	    #
	    # # we check length consistency
	    # if(NROW(data) != (object$nobs + length(object$obsRemoved))){
	    #     stop("To get the coefficients for the variables with varying slopes, we need to fetch these variables in the original dataset (", deparse_long(dataName), "), yet the dataset doesn't have the same number of observations as was used in the estimation (", NROW(data), " instead of ", object$nobs + length(object$obsRemoved), ").")
	    # }
	    #
	    # if(length(object$obsRemoved)){
	    #     data = data[-object$obsRemoved, slope_vars_unik_raw, drop = FALSE]
	    # } else {
	    #     data = data[, slope_vars_unik_raw, drop = FALSE]
	    # }

	    slope_var_list = list()
	    for(i in 1:length(slope_vars_unik)){

	        # # as.numeric => we'll use cpp so required
	        # svar = as.numeric(as.vector(eval(parse(text = slope_vars_unik[i]), data)))
	        # if(length(svar) == 1) svar = rep(svar, N)
	        #
	        # slope_var_list[[slope_vars_unik[i]]] = svar
	        slope_var_list[[slope_vars_unik[i]]] = object$slope_variables[[slope_vars_unik[i]]]
	    }

	    # STEP 2: demeaning

	    # This way to proceed is "dirty", I shall do some cleanup when I have more time
	    fixef_sizes = object$fixef_sizes

	    nb_cluster_all = as.integer(unlist(fixef_sizes[dict_fe[fe_all]]))

	    dum_vector = as.integer(unlist(id_dummies_vect[dict_fe[fe_all]])) - 1L

	    slope_flag = as.integer(grepl("\\[", fixef_terms))

	    what = slope_var_list[slope_vars]
	    names(what) = NULL
	    slope_vars_vector = unlist(what)

	    fixef_table_list = list()
	    for(i in 1:Q) fixef_table_list[[i]] = cpp_table(fixef_sizes[i], id_dummies_vect[[i]])
	    fixef_table_vector = as.integer(unlist(fixef_table_list[dict_fe[fe_all]]))

	    S_demean <- cpp_demean(y = S, X_raw = 0, r_weights = 0, iterMax = 1000L,
	                              diffMax = 1e-4, nb_cluster_all = nb_cluster_all,
	                              dum_vector = dum_vector, tableCluster_vector = fixef_table_vector,
	                              slope_flag = slope_flag, slope_vars = slope_vars_vector,
	                              r_init = 0, checkWeight = 1L, nthreads = 1L, save_fixef = TRUE)

	    fixef_coef = S_demean$fixef_coef

	    fixef_values = list()
	    cum_sizes = cumsum(fixef_sizes[fe_all])
	    start = 1 + c(0, cum_sizes)
	    for(i in seq_along(fixef_terms)){
	        fixef_values[[i]] = fixef_coef[start[i]:cum_sizes[i]]
	    }

	    #
	    # Now the referenes
	    #

	    # FE references
	    Q_fe = length(fixef_vars)
	    if(Q_fe >= 2){

	        my_dum = id_dummies_vect[dict_fe[fixef_vars]]

	        dumMat <- matrix(unlist(my_dum), N, Q_fe) - 1
	        orderCluster <- matrix(unlist(lapply(my_dum, order)), N, Q_fe) - 1

	        nbCluster = sapply(my_dum, max)

	        fixef_values_tmp <- cpp_get_fe_gnl(Q_fe, N, rep(1, N), dumMat, nbCluster, orderCluster)

	        # the information on the references
	        nb_ref_fe = fixef_values_tmp[[Q_fe+1]]
	    } else {
	        nb_ref_fe = integer(Q_fe)
	    }

	    # Slope references (if associated FE + constant)

        Q_slope = length(slope_fe)
        nb_ref_slope = integer(Q_slope)
        for(i in 1:Q_slope){

            my_dum = id_dummies_vect[[dict_fe[slope_fe[i]]]]
            my_order = order(my_dum)
            var_sorted = slope_var_list[[slope_vars[i]]][my_order]

            # if no associated FE => we check only 0 values
            if(!slope_fe[i] %in% fixef_vars){
                nb_ref_slope[i] = cpp_constant_dum(fixef_sizes[slope_fe[i]], var_sorted, my_dum[my_order], only_0 = TRUE)
            } else {
                nb_ref_slope[i] = cpp_constant_dum(fixef_sizes[slope_fe[i]], var_sorted, my_dum[my_order])
            }
        }

        nb_ref = integer(length(fixef_terms))
        nb_ref[slope_flag] = nb_ref_slope
        nb_ref[!slope_flag] = nb_ref_fe

        # we recreate that to avoid conditioning on isSlope later
        fixef_id = fixef_id[dict_fe[fe_all]]
        fixef_names = fixef_terms

    } else if(Q == 1){
		# This is the simplest case
		id = id_dummies_vect[[1]]

		myOrder = order(id)
		myDiff = c(1, diff(id[myOrder]))

		select = myOrder[myDiff == 1]

		fixef_values = list(S[select])

		# There are no references => no need to set nb_ref
	} else if(FALSE && Q == 2) {
		# specific method that is faster but specific to the case of 2 FE
	    # NOTA:
	    # - the next method is even faster
	    # - but this current method has the benefit to provide the components
		dum_i = id_dummies_vect[[1]] - 1L
		dum_j = id_dummies_vect[[2]] - 1L

		order_ij = order(dum_i, dum_j)
		i_sorted_index_j = dum_j[order_ij]

		order_ji = order(dum_j, dum_i)
		j_sorted_index_i = dum_i[order_ji]

		i_sorted_sumFE = S[order_ij]
		j_sorted_sumFE = S[order_ji]

		fe <- cpp_get_fe_2(clusterSize = object$fixef_sizes, i_sorted_index_j = i_sorted_index_j, i_sorted_sumFE = i_sorted_sumFE, j_sorted_index_i = j_sorted_index_i, j_sorted_sumFE = j_sorted_sumFE, r_cumtable_i = cumsum(table(dum_i)), r_cumtable_j = cumsum(table(dum_j)))

		# cpp_get_fe_2 returns a matrix (we will update cpp_get_fe_gnl to return a matrix later)
		# so we put it into a list (for now)

		fixef_values = list()
		fixef_values[[1]] = fe[fe[, 1] == 1, 3]
		fixef_values[[2]] = fe[fe[, 1] == 2, 3]

		# the number of references needed
		nb_ref = c(0, sum(fe[, 4]))

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
		cv = fixef_values[[i]]
		names(cv) = attr(fixef_id[[i]], "fixef_names")
		all_clust[[fixef_names[i]]] = cv
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
#' \item Here is the help from package \pkg{nlme}: \code{\link[nlme]{fixef}}. The help from package \pkg{fixest} is here: \code{\link[fixest]{fixef.fixest}}.
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
#' This function plots the 5 fixed-effects with the highest and lowest values, for each of the clusters. It takes as an argument the fixed-effects obtained from the function \code{\link{fixef.fixest}} after an estimation using \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @method plot fixest.fixef
#'
#' @param x An object obtained from the function \code{\link{fixef.fixest}}.
#' @param n The number of fixed-effects to be drawn. Defaults to 5.
#' @param ... Not currently used.
#'
#' Note that the fixed-effect coefficients might NOT be interpretable. This function is useful only for fully regular panels.
#'
#' If the data are not regular in the cluster coefficients, this means that several \sQuote{reference points} are set to obtain the fixed-effects, thereby impeding their interpretation. In this case a warning is raised.
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
    any_invalid = check_dots_args(match.call(expand.dots = FALSE), suggest_args = "n")
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

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



#' Finds observations to be removed from ML estimation with factors/clusters
#'
#' For Poisson, Negative Binomial or Logit estimations with fixed-effects, when the dependent variable is only equal to 0 (or 1 for Logit) for one cluster value this leads to a perfect fit for that cluster value by setting its associated cluster coefficient to \code{-Inf}. Thus these observations need to be removed before estimation. This function gives the observations to be removed. Note that by default the function \code{\link[fixest]{femlm}} or \code{\link[fixest]{feglm}} drops them before performing the estimation.
#'
#' @param fml A formula containing the dependent variable and the clusters. It can be of the type: \code{y ~ cluster_1 + cluster_2} or \code{y ~ x1 | cluster_1 + cluster_1} (in which case variables before the pipe are ignored).
#' @param data A data.frame containing the variables in the formula.
#' @param family Character scalar: either \dQuote{poisson} (default), \dQuote{negbin} or \dQuote{logit}.
#'
#' @return
#' It returns an integer vector of observations to be removed. If no observations are to be removed, an empty integer vector is returned. In both cases, it is of class \code{fixest.obs2remove}.
#' The vector has an attribute \code{cluster} which is a list giving the IDs of the clusters that have been removed, for each cluster dimension.
#'
#' @examples
#'
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' (x = obs2remove(v6 ~ Species, base))
#' attr(x, "cluster")
#'
#' # The two results are identical:
#' res_1 = femlm(v6 ~ Petal.Width | Species, base)
#' # => note + obsRemoved is created
#'
#' res_2 = femlm(v6 ~ Petal.Width | Species, base[-x, ])
#' # => no note because observations are removed before
#'
#' esttable(res_1, res_2)
#'
#' all(res_1$obsRemoved == x)
#'
obs2remove = function(fml, data, family = c("poisson", "negbin", "logit")){
	# in the formula, the clusters must be there:
	# either y ~ cluster_1 + cluster_2
	# either y ~ x1 + x2 | cluster_1 + cluster_2

	#
	# CONTROLS
	#

	# FAMILY

	family = match.arg(family)

	# FML

	if(!"formula" %in% class(fml) || length(fml) != 3){
		stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | cluster_1 + cluster_1' or of the type 'y ~ cluster_1 + cluster_2'.")
	}

	FML = Formula::Formula(fml)
	n_rhs = length(FML)[2]

	if(n_rhs > 2){
		stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | cluster_1 + cluster_1' or of the type 'y ~ cluster_1 + cluster_2'.")
	}

	# DATA

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
	if("data.table" %in% class(data)){
		# this is a local change only
		class(data) = "data.frame"
	}

	dataNames = names(data)

	# Extracting the variables
	vars_left = all.vars(formula(FML, lhs=1, rhs=0))
	cluster_fml = formula(FML, lhs=0, rhs=n_rhs)
	vars_clusters = all.vars(cluster_fml)

	if(length(left_missing <- setdiff(vars_left, dataNames)) > 0){
		stop("Left hand side could not be evaluated, following variables are missing from the data: ", paste0(left_missing, collapse = ", "), ".")
	}

	if(length(right_missing <- setdiff(vars_clusters, dataNames)) > 0){
		stop("The clsuters could not be evaluated, following variables are missing from the data: ", paste0(right_missing, collapse = ", "), ".")
	}

	# Evaluation variables
	lhs = as.vector(eval(fml[[2]], data))
	cluster_mat = model.frame(cluster_fml, data)
	cluster_name = names(cluster_mat)

	#
	# -- CORE --
	#

	Q = length(cluster_name)
	dummyOmises = list()
	obs2remove = c()
	for(q in 1:Q){

		dum_raw = cluster_mat[, q]

		# thisNames = getItems(dum_raw)
		# dum = quickUnclassFactor(dum_raw)
		dum_all = quickUnclassFactor(dum_raw, addItem = TRUE)
		dum = dum_all$x
		thisNames = dum_all$items
		k = length(thisNames)

		# We delete "all zero" outcome
		sum_y_clust = cpp_tapply_vsum(k, lhs, dum)
		n_perClust = cpp_table(k, dum)

		if(family %in% c("poisson", "negbin")){
			qui = which(sum_y_clust == 0)
		} else if(family == "logit"){
			qui = which(sum_y_clust == 0 | sum_y_clust == n_perClust)
		}

		if(length(qui > 0)){
			# We first delete the data:
			dummyOmises[[q]] = thisNames[qui]
			obs2remove = unique(c(obs2remove, which(dum %in% qui)))
		} else {
			dummyOmises[[q]] = character(0)
		}
	}

	names(dummyOmises) = cluster_name

	if(length(obs2remove) == 0){
		print("No observation to be removed.")
		obs2remove = integer(0)
		class(obs2remove) = "fixest.obs2remove"
		return(invisible(obs2remove))
	}

	class(obs2remove) = "fixest.obs2remove"
	attr(obs2remove, "family") = family
	attr(obs2remove, "cluster") = dummyOmises

	return(obs2remove)
}



#' Summary method for fixest.obs2remove objects
#'
#' This function synthesizes the information of function \code{\link[fixest]{obs2remove}}. It reports the number of observations to be removed as well as the number of clusters removed per cluster dimension.
#'
#' @method summary fixest.obs2remove
#'
#' @param object A \code{fixest.obs2remove} object obtained from function \code{\link[fixest]{obs2remove}}.
#' @param ... Not currently used.
#'
#'
#' @examples
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' x = obs2remove(v6 ~ Species, base)
#' summary(x)
#'
summary.fixest.obs2remove = function(object, ...){

	if(length(object) == 0){
		print("No observation to be removed.")
	} else {
		cat(length(object), " observations removed because of only zero", ifelse(attr(object, "family") == "logit", ", or only one,", ""), " outcomes.\n", sep = "")
		cluster = attr(object, "cluster")
		cat("# clusters removed: ", paste0(names(cluster), ": ", lengths(cluster), collapse = ", "), ".", sep = "")
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
#' This function tests: 1) collinearity with the cluster variables, 2) perfect multi-collinearity between the variables, 4) perfect multi-collinearity between several variables and the clusters, and 4) identification issues when there are non-linear in parameters parts.
#'
#' @return
#' It returns a text message with the identified diagnostics.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' \donttest{
#' # Creating an example data base:
#' cluster_1 = sample(3, 100, TRUE)
#' cluster_2 = sample(20, 100, TRUE)
#' x = rnorm(100, cluster_1)**2
#' y = rnorm(100, cluster_2)**2
#' z = rnorm(100, 3)**2
#' dep = rpois(100, x*y*z)
#' base = data.frame(cluster_1, cluster_2, x, y, z, dep)
#'
#' # creating collinearity problems:
#' base$v1 = base$v2 = base$v3 = base$v4 = 0
#' base$v1[base$cluster_1 == 1] = 1
#' base$v2[base$cluster_1 == 2] = 1
#' base$v3[base$cluster_1 == 3] = 1
#' base$v4[base$cluster_2 == 1] = 1
#'
#' # Estimations:
#'
#' # Collinearity with the cluster variables:
#' res_1 = femlm(dep ~ log(x) + v1 + v2 + v4 | cluster_1 + cluster_2, base)
#' collinearity(res_1)
#' # => collinearity with cluster identified, we drop v1 and v2
#' res_1bis = femlm(dep ~ log(x) + v4 | cluster_1 + cluster_2, base)
#' collinearity(res_1bis)
#'
#' # Multi-Collinearity:
#' res_2 =  femlm(dep ~ log(x) + v1 + v2 + v3 + v4, base)
#' collinearity(res_2)
#'
#' # In non-linear part:
#' res_3 = feNmlm(dep ~ log(z), base, NL.fml = ~log(a*x + b*y),
#'               NL.start = list(a=1, b=1), lower = list(a=0, b=0))
#' collinearity(res_3)
#' }
#'
collinearity = function(x, verbose){
	# x: fixest estimation

	if(class(x) != "fixest"){
		stop("Argument 'x' must be a fixest object.")
	}

	# I) (linear) collinearity with clusters
	# II) (linear) multi collinearity
	# III) (non-linear) overidentification

	# flags
    isCluster = !is.null(x$fixef_id)
    isFE = FALSE
    isSlope = FALSE
    if(isCluster){
        isSlope = !is.null(x$fixef_terms)
        if(isSlope){
            fixef_terms = x$fixef_terms
            slope_flag = grepl("\\[", fixef_terms)
            isFE = any(!slope_flag)
        } else {
            isFE = TRUE
        }
    }

	linear_fml = formula(Formula(formula(x, "linear")), lhs=0, rhs=1)

	rhs_fml = formula(Formula(formula(x, "linear")), lhs = 1, rhs = 1)
	if(grepl("[^:]::[^:]", deparse_long(rhs_fml[[3]]))){
	    new_fml = interact_fml(rhs_fml)
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

	if(isLinear || isCluster || "(Intercept)" %in% names(coef)){
		# linear.matrix = model.matrix(linear_fml, data)
		linear.matrix = fixest_model_matrix(rhs_fml, data)
	}

	if(!is.null(x$obsRemoved)){
	    linear.matrix = linear.matrix[-x$obsRemoved, ]
	    # We do that to drop interaction variables that should not be there any more
	    # if factors with only NA values
	    varkeep = intersect(names(x$coefficients), colnames(linear.matrix))
	    linear.matrix = linear.matrix[, varkeep]
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

	# for multicol without cluster
	mat_base = as.matrix(linbase)

	# we add possibly missing variables
	varmiss = setdiff(all.vars(linear_fml), names(linbase))
	for(v in varmiss) linbase[[v]] = data[[v]]

	if(isLinear && length(linearVars) >= 2){
		ccat(ifelse(Q >= 1, ", ", " "), "multiple:")
		for(v in linearVars){
			ccat(".")
			# fml2estimate = as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))
			# res = lm(fml2estimate, linbase)
			# sum_resid = sum(abs(resid(res)))

			i = which(colnames(mat_base) == v)
			res = ols_fit(y = mat_base[, i], X = mat_base[, -i, drop = FALSE], w = 1, nthreads = 1)

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
	if(isLinear && length(linearVars) >= 2 && isCluster){
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

				max_residuals = max(abs(resid(res)))

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

	if(isNL && (length(coef) >= 2 || isCluster)){
		ccat(", in non-linear term:")
		if(isCluster){
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

				if(isCluster){
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

			sum_resids = sum(abs(resid(res)))
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


#' Plots confidence intervals
#'
#' This function plots the results of estimations (coefficients and confidence intervals). It is flexible and handle interactions in a special way.
#'
#' @inheritParams etable
#'
#' @param object Can be either: i) an estimation object (obtained for example from \code{\link[fixest]{feols}}, ii) a matrix of coefficients table, or iii) a numeric vector of the point estimates -- the latter requiring the extra arguments \code{sd} or \code{ci_low} and \code{ci_top}.
#' @param sd The standard errors of the estimates. It may be missing.
#' @param ci_low If \code{sd} is not provided, the lower bound of the confidence interval. For each estimate.
#' @param ci_top If \code{sd} is not provided, the upper bound of the confidence interval. For each estimate.
#' @param x The value of the x-axis. If missing, the names of the argument \code{estimate} are used.
#' @param x.shift Shifts the confidence intervals bars to the left or right, depending on the value of \code{x.shift}. Default is 0.
#' @param ci.width The width of the extremities of the confidence intervals. Default is \code{0.1}.
#' @param ci_level Scalar between 0 and 1: the level of the CI. By default it is equal to 0.95.
#' @param add Default is \code{FALSE}, if the intervals are to be added to an existing graph. Note that if it is the case, then the argument \code{x} MUST be numeric.
#' @param pt.pch The patch of the coefficient estimates. Default is 20 (circle).
#' @param cex Numeric, default is \code{par("cex")}. Expansion factor for the points
#' @param pt.cex The size of the coefficient estimates. Default is the other argument \code{cex}.
#' @param col The color of the points and the confidence intervals. Default is 1 ("black"). Note that you can set the colors separately for each of them with \code{pt.col} and \code{ci.col}.
#' @param pt.col The color of the coefficient estimate. Default is equal to the other argument \code{col}.
#' @param ci.col The color of the confidence intervals. Default is equal to the other argument \code{col}.
#' @param lwd General liwe with. Default is par("lwd").
#' @param ci.lwd The line width of the confidence intervals. Default is equal to the other argument \code{lwd}.
#' @param ci.lty The line type of the confidence intervals. Default is 1.
#' @param grid Logical, default is \code{TRUE}. Whether a grid should be displayed. You can set the display of the grid with the argument \code{grid.par}.
#' @param grid.par List. Parameters of the grid. The default values are: \code{lty = 3} and \code{col = "gray"}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. You also have two additional arguments: use \code{horiz = FALSE} to disable the horizontal lines, and use \code{vert = FALSE} to disable the vertical lines. Eg: \code{grid.par = list(vert = FALSE, col = "red", lwd = 2)}.
#' @param zero Logical, default is \code{TRUE}. Whether the 0-line should be emphasized. You can set the parameters of that line with the argument \code{zero.par}.
#' @param zero.par List. Parameters of the zero-line. The default values are \code{col = "black"} and \code{lwd = 1}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. Example: \code{zero.par = list(col = "darkblue", lwd = 3)}.
#' @param join Logical, default depends on the situation. If \code{TRUE}, then the coefficient estimates are joined with a line. By default, it is equal to \code{TRUE} only if: i) interactions are plotted, ii) the x values are numeric and iii) a reference is found.
#' @param join.par List. Parameters of the line joining the cofficients. The default values are: \code{col = pt.col} and \code{lwd = lwd}. You can add any graphical parameter that will be passed to \code{\link[graphics]{lines}}. Eg: \code{join.par = list(lty = 2)}.
#' @param ref.line Logical, default depends on the situation. It is \code{TRUE} only if: i) interactions are plotted, ii) the x values are numeric and iii) a reference is found. If \code{TRUE}, then a vertical line is drawn at the level of the reference value. You can set the parameters of this line with the argument \code{ref.line.par}.
#' @param ref.line.par List. Parameters of the vertical line on the reference. The default values are: \code{col = "black"} and \code{lty = 2}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. Eg: \code{ref.line.par = list(lty = 1, lwd = 3)}.
#' @param xlim.add A numeric vector of length 1 or 2. It represents an extension factor of xlim, in percentage. Eg: \code{xlim.add = c(0, 0.5)} extends \code{xlim} of 50\% on the right. If of lentgh 1, positive values represent the right, and negative values the left (Eg: \code{xlim.add = -0.5} is equivalent to \code{xlim.add = c(0.5, 0)}).
#' @param ylim.add A numeric vector of length 1 or 2. It represents an extension factor of ylim, in percentage. Eg: \code{ylim.add = c(0, 0.5)} extends \code{ylim} of 50\% on the top. If of lentgh 1, positive values represent the top, and negative values the bottom (Eg: \code{ylim.add = -0.5} is equivalent to \code{ylim.add = c(0.5, 0)}).
#' @param only.params Logical, default is \code{FALSE}. If \code{TRUE} no graphic is displayed, only the values of \code{x} and \code{y} used in the plot are returned.
#' @param only.inter Logical, default is \code{TRUE}. If an interaction of the type of \code{var::fe} (see \code{\link[fixest]{feols}} help for details) is found, then only these interactions are plotted. If \code{FALSE}, then interactions are treated as regular coefficients.
#' @param ... Other arguments to be passed to \code{summary}, if \code{object} is an estimation, and/or to the function \code{plot} or \code{lines} (if \code{add = TRUE}).
#'
#' @seealso
#' See \code{\link[fixest]{setFixest_coefplot}} to set the default values of \code{coefplot}, and the estimation functions: e.g. \code{\link[fixest]{feols}}, \code{\link[fixest]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest]{fenegbin}}.
#'
#' @section Setting custom default values:
#' The function \code{coefplot} dispose of many arguments to parametrize the plots. Most of these arguments can be set once an for all using the function \code{\link[fixest]{setFixest_coefplot}}. See Example 3 below for a demonstration.
#'
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' #
#' # Example 1: Stacking two sets of results on the same graph
#' #
#'
#' # Estimation on Iris data with one fixed-effect (Species)
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Length +
#'             Sepal.Width | Species, iris)
#'
#'
#' # First graph with clustered standard-errors
#' # (the default when fixed-effects are present)
#' coefplot(est, x.shift = -.2)
#'
#' # 'x.shift' was used to shift the coefficients on the left.
#'
#'
#' # Second set of results: this time with
#' #  standard-errors that are not clustered.
#' coefplot(est, se = "standard", x.shift = .2,
#'          add = TRUE, col = 2, ci.lty = 2, pch=15)
#'
#'  # Note that we used 'se', an argument that will
#'  #  be passed to summary.fixest
#'
#' legend("topright", col = 1:2, pch = 20, lwd = 1, lty = 1:2,
#'        legend = c("Clustered", "Standard"), title = "Standard-Errors")
#'
#'
#' #
#' # Example 2: Interactions
#' #
#'
#'
#' # Now we estimate and plot the "yearly" treatment effects
#'
#' data(base_did)
#' base_inter = base_did
#'
#' # We interact the variable 'period' with the variable 'treat'
#' est_did = feols(y ~ x1 + i(treat, period, 5) | id+period, base_inter)
#'
#' # You could have written the following formula instead:
#' # y ~ x1 + treat::period(5) | id+period
#'
#' # In the estimation, the variable treat is interacted
#' #  with each value of period but 5, set as a reference
#'
#' # When estimations contain interactions, as before,
#' #  the default behavior of coefplot changes,
#' #  it now only plots interactions:
#' coefplot(est_did)
#'
#' # We can see that the graph is different from before:
#' #  - only interactions are shown,
#' #  - the reference is present,
#' #  - the estimates are joined.
#' # => this is fully flexible
#'
#' coefplot(est_did, ref.line = FALSE, join = FALSE)
#'
#' # Now to display all coefficients, use 'only.inter'
#' coefplot(est_did, only.inter = FALSE)
#'
#' #
#' # What if the interacted variable is not numeric?
#'
#' # Let's create a "month" variable
#' all_months = c("aug", "sept", "oct", "nov", "dec", "jan",
#'                "feb", "mar", "apr", "may", "jun", "jul")
#' base_inter$period_month = all_months[base_inter$period]
#'
#' # The new estimation
#' est = feols(y ~ x1 + i(treat, period_month, "oct") | id+period, base_inter)
#' # Since 'period_month' of type character, coefplot sorts it
#' coefplot(est)
#'
#' # To respect a plotting order, use a factor
#' base_inter$month_factor = factor(base_inter$period_month, levels = all_months)
#' est = feols(y ~ x1 + i(treat, month_factor, "oct") | id+period, base_inter)
#' coefplot(est)
#'
#'
#' #
#' # Example 3: Setting defaults
#' #
#'
#' # coefplot has many arguments, which makes it highly flexible.
#' # If you don't like the default style of coefplot. No worries,
#' # you can set *your* default by using the function
#' # setFixest_coefplot()
#'
#' dict = c("Petal.Length"="Length (Petal)", "Petal.Width"="Width (Petal)",
#'          "Sepal.Length"="Length (Sepal)", "Sepal.Width"="Width (Sepal)")
#'
#' setFixest_coefplot(ci.col = 2, pt.col = "darkblue", ci.lwd = 3,
#'                    pt.cex = 2, pt.pch = 15, ci.width = 0, dict = dict)
#'
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Length +
#'                 Sepal.Width | Species, iris)
#'
#' # Tadaaa!
#' coefplot(est)
#'
#' # To reset to the default settings:
#' setFixest_coefplot()
#' coefplot(est)
#'
#'
coefplot = function(object, sd, ci_low, ci_top, x, x.shift = 0, dict, drop, order, ci.width=0.1, ci_level = 0.95, add = FALSE, pt.pch = 20, cex = par("cex"), pt.cex = cex, col = 1, pt.col = col, ci.col = col, lwd = par("lwd"), ci.lwd = lwd, ci.lty = 1, grid = TRUE, grid.par = list(lty=3, col = "gray"), zero = TRUE, zero.par = list(col="black", lwd=1), join = FALSE, join.par = list(col = pt.col, lwd=lwd), ref.line, ref.line.par = list(col = "black", lty = 2), xlim.add, ylim.add, only.params = FALSE, only.inter = TRUE, ...){
	# creation de barres d'erreur
	# 1 segment vertical de la taille de l'IC
	# deux barres horizontales, a chaque bornes

    # I don't allow the different styles any more (maybe I can change that in the future)
    if(missing(dict)) dict = c()
    # style = c("bar", "interval", "tube")
    # style = match.arg(style)
    style = "bar"

    # We get the default values
    opts = getOption("fixest_coefplot")

    if(length(opts) >= 1){
        if(!is.list(opts)){
            warning("The default values of coefplot are ill-formed and therefore reset. Use only setFixest_coefplot for setting the default values.")
            opts = list()
            options("fixest_coefplot" = opts)
        } else {
            mc = match.call()
            arg2set = setdiff(names(opts), names(mc))
            for(arg in arg2set){
                assign(arg, opts[[arg]])
            }
        }
    }

    #
    # STEP 1 => getting the coefficient table + the CI
    #

    # This is NOT trivial because argument '...' refers to either summary (and we
    # have to find out which), either to plot or lines
    #

    # Object can either be:
    #  a) an estimation
    #  b) a coefficient table (ie a matrix)
    #  c) a vector of coefficients

    dots = list(...)

    if(is.list(object)){
        sum_exists = FALSE
        for(c_name in class(object)){
            if(exists(paste0("summary.", c_name), mode = "function")){
                sum_exists = TRUE
                break
            }
        }

        if(!sum_exists) stop("There is no summary method for objects of class ", c_name, ". 'coefplot' applies summary to the object to extract the coeftable. Maybe add directly the coeftable in object instead?")

        fun_name = paste0("summary.", c_name)
        args_name_sum = names(formals(fun_name))
        args_sum = intersect(names(dots), args_name_sum)

        # we kick out the summary arguments from dots
        dots[args_sum] = NULL

        # We reconstruct a call to coeftable
        mc_coeftable = match.call(expand.dots = TRUE)
        mc_coeftable[[1]] = as.name("coeftable")
        mc_coeftable[setdiff(names(mc_coeftable), c(args_sum, "object", ""))] = NULL

        mat = eval(mc_coeftable, parent.frame())

        sd = mat[, 2]
        estimate = mat[, 1]

        names(estimate) = rownames(mat)

        if("fml" %in% names(object)){
            depvar = gsub(" ", "", as.character(object$fml)[[2]])
            if(depvar %in% names(dict)) depvar = dict[depvar]
            listDefault(dots, "main", paste0("Effect on ", depvar))
        }

    } else if(is.matrix(object)){
        # object is a matrix containing the coefs and SEs

        m_names = tolower(colnames(object))
        if(ncol(object) == 4 || (grepl("estimate", m_names[1]) && grepl("std\\.? error", m_names[1]))){
            sd = object[, 2]
            estimate = object[, 1]

            names(estimate) = rownames(object)

        } else {
            stop("Argument 'object' is a matrix but it should contain 4 columns (the two first ones should be reporting the estimate and the standard-error). Either provide an appropriate matrix or give directly the vector of estimated coefficients in arg. estimate.")
        }
    } else if(length(object[1]) > 1 || !is.null(dim(object)) || !is.numeric(object)){
        stop("Argument 'object' must be either: i) an estimation object, ii) a matrix of coefficients table, or iii) a numeric vector of the point estimates. Currently it is neither of the three.")
    } else {
        # it's a numeric vector
        estimate = object
    }

	n <- length(estimate)

	if(missing(ci.lty)){
		ci.lty = ifelse(style == "bar", 1, 2)
	}

	if(missing(sd)){
		if(missing(ci_low) || missing(ci_top)) stop("If 'sd' is not provided, you must provide the arguments 'ci_low' and 'ci_top'.")

		ci025 = ci_low
		ci975 = ci_top

	} else {
		if(!missing(ci_low) || !missing(ci_top)) warning("Since 'sd' is provided, arguments 'ci_low' or 'ci_top' are ignored.")

		# We compue the CI
		nb = abs(qnorm((1-ci_level)/2))
		ci975 = estimate + nb*sd
		ci025 = estimate - nb*sd
	}

	#
	# STEP 2 => Plotting
	#

	# INTERACTIONS
	ref_id = NA
	if(only.inter && !is.null(names(estimate))){
    	all_vars = names(estimate)
    	if(any(grepl("::", all_vars))){

    	    is_info = FALSE
    	    if("fixest" %in% class(object)){
    	        is_info = TRUE
    	        interaction.info = object$interaction.info
    	        is_ref = interaction.info$is_ref
    	        items = interaction.info$items
    	        is_num = interaction.info$fe_type %in% c("numeric", "integer")
    	    }

    	    # We retrict only to interactions
    	    root_interaction = all_vars[grepl("::", all_vars)]
    	    # We keep only the first one !
    	    root_interaction = unique(gsub("::.+", "", root_interaction))[1]

    	    inter_keep = grepl(root_interaction, all_vars, fixed = TRUE)
    	    my_inter = all_vars[inter_keep]
    	    estimate = estimate[inter_keep]
    	    ci975 = ci975[inter_keep]
    	    ci025 = ci025[inter_keep]

    	    # We extract the name of the variables
    	    fe_name = gsub(".+:", "", root_interaction)
    	    if(fe_name %in% names(dict)) fe_name = dict[fe_name]
    	    listDefault(dots, "xlab", fe_name)

    	    var_name = gsub(":[[:alnum:]\\._]+", "", root_interaction)
    	    if(var_name %in% names(dict)) var_name = dict[var_name]
    	    listDefault(dots, "sub", paste0("(Interacted with ", var_name, ")"))

    	    # We construct the x-axis
    	    inter_values = gsub(".+::", "", my_inter)
    	    names(estimate) = inter_values

    	    inter_values_num = tryCatch(as.numeric(inter_values), warning = function(x) x)
    	    if(is_info){

    	        if(length(inter_values) != sum(!is_ref)){
    	            stop("Internal error regarding the lengths of vectors of coefficients.")
    	        }

    	        if(any(is_ref)){
    	            ref_id = which(is_ref)
    	        }

    	        my_values = my_ci_low = my_ci_high = rep(NA, length(is_ref))
    	        names(my_values) = names(my_ci_low) = names(my_ci_high) = items

    	        my_values[inter_values] = estimate
    	        my_ci_high[inter_values] = ci975
    	        my_ci_low[inter_values] = ci025

    	        qui = which(is.na(my_values))
    	        my_values[qui] = 0
    	        my_ci_high[qui] = my_values[qui]
    	        my_ci_low[qui] = my_values[qui]

    	        estimate = my_values
    	        ci975 = my_ci_high
    	        ci025 = my_ci_low

    	        if(is_num && missing(x)){
    	            names(estimate) = NULL
    	            x = items
    	            if(missing(join)) join = TRUE
    	            if(missing(ref.line)) ref.line = TRUE
    	        }


    	    } else if(is.numeric(inter_values_num)) {

    	        # We check these are "real" numbers and not just "codes"
    	        all_steps = diff(sort(inter_values_num))
    	        ts = table(all_steps)
    	        step_mode = as.numeric(names(ts)[which.max(ts)])
    	        all_steps_rescaled = all_steps / step_mode

    	        if(any(all_steps_rescaled < 10) && missing(x)){
    	            names(estimate) = NULL
    	            x = inter_values_num
    	            if(missing(join)) join = TRUE
    	        }
    	    }

	        n = length(estimate)
    	}
	}

	# We set the default after the interactions (which is the decision variable)
	if(missing(join)) join = FALSE
	if(missing(ref.line)) ref.line = FALSE

	# setting the names of the estimate
	if(!missing(x)){
	    if(length(x) != n) stop("Argument 'x' must have the same length as the number of coefficients (currently ", length(x), " vs ", n, ").")

	    if(!is.numeric(x)){
	        names(estimate) = x
	    }
	}

	# order/drop/dict
	if(!is.null(names(estimate))){

	    if(missnull(dict)){
	        dict = c()
	    } else {
	        check_arg(dict, "characterVector")
	        if(is.null(names(dict))){
                stop("Argument 'dict' must be a named character vector. Currently it has no names.")
	        }
	    }

	    # dict
	    all_vars = names(estimate)

	    who = all_vars %in% names(dict)
	    all_vars[who] = dict[all_vars[who]]
	    names(estimate) = all_vars

	    my_order = 1:n
	    names(my_order) = all_vars

	    # dropping some coefs
	    all_vars = drop_apply(all_vars, drop)

	    if(length(all_vars) == 0){
	        stop("Argument 'drop' has removed all variables!")
	    }

	    # ordering the coefs
	    all_vars = order_apply(all_vars, order)

	    qui = my_order[all_vars]

	    estimate = estimate[qui]
	    ci975 = ci975[qui]
	    ci025 = ci025[qui]
	    if(!missing(x)) x = x[qui]
	    n = length(qui)
	}

	# we create x_labels, x_value & x_at
	if(!missing(x) && is.numeric(x)){
	    my_xlim = range(c(x + x.shift, x - x.shift))

	    x_at = NULL
	    x_value = x + x.shift

	    if(is.null(names(estimate))){
	        x_labels = NULL
	    } else {
	        x_labels = names(estimate)
	    }

	} else {
	    x_at = 1:n
	    x_value = 1:n + x.shift

		if(is.null(names(estimate))){
		    x_labels = 1:n
		} else {
		    x_labels = names(estimate)
		}

		my_xlim = range(c(1:n + x.shift, 1:n - x.shift)) + c(-0.5, +0.5)
	}

	all_plot_args = unique(c(names(par()), names(formals(plot.default))))
	pblm = setdiff(names(dots), all_plot_args)
	if(length(pblm) > 0){
	    warning("The following argument", ifsingle(pblm, " is not a", "s are not"), " plotting argument", ifsingle(pblm, " and is", "s and are"), " therefore ignored: ", enumerate_items(pblm), ".")
	    dots[pblm] = NULL
	}

	# preparation of the do.call
	dots$col = col
	listDefault(dots, "xlab", "Variable")
	ylab = paste0("Estimate and ", ifelse(missing(sd), "", paste0(floor(ci_level*100), "% ")), "Conf. Int.")
	listDefault(dots, "ylab", ylab)

	# The limits

	# xlim
	if(!missnull(xlim.add)){
	    if("xlim" %in% names(dots)){
	        message("Since argument 'xlim' is provided, argument 'xlim.add' is ignored.")
	    } else {
	        if((!is.numeric(xlim.add) || !length(xlim.add) %in% 1:2)){
	            stop("Argument 'xlim.add' must be a numeric vector of length 1 or 2. It represents an extension factor of xlim, in percentage. (Eg: xlim.add = c(0, 0.5) extends xlim of 50% on the right.) If of lentgh 1, positive values represent the right, and negative values the left (Eg: xlim.add = -0.5 is equivalent to xlim.add = c(0.5, 0)).")
	        }

	        if(length(xlim.add) == 1){
	            if(xlim.add > 0) {
	                xlim.add = c(0, xlim.add)
	            } else {
	                xlim.add = c(xlim.add, 0)
	            }
	        }

	        x_width = diff(my_xlim)
	        my_xlim = my_xlim + xlim.add * x_width
	    }
	}
	listDefault(dots, "xlim", my_xlim)

	# ylim
	my_ylim = range(c(ci025, ci975))

	if(!missnull(ylim.add)){
	    if("ylim" %in% names(dots)){
	        message("Since argument 'ylim' is provided, argument 'ylim.add' is ignored.")
	    } else {
	        if((!length(ylim.add) %in% 1:2 || !is.numeric(ylim.add))){
	            stop("Argument 'ylim.add' must be a numeric vector of length 1 or 2. It represents an extension factor of ylim, in percentage. (Eg: ylim.add = c(0, 0.5) extends ylim of 50% on the top.) If of lentgh 1, positive values represent the top, and negative values the bottom (Eg: ylim.add = -0.5 is equivalent to ylim.add = c(0.5, 0)).")
	        }

	        if(length(ylim.add) == 1){
	            if(ylim.add > 0) {
	                ylim.add = c(0, ylim.add)
	            } else {
	                ylim.add = c(ylim.add, 0)
	            }
	        }

	        y_width = diff(my_ylim)
	        my_ylim = my_ylim + ylim.add * y_width
	    }
	}

	listDefault(dots, "ylim", my_ylim)

    dots$x = x_value

	dots$y = estimate

	if(style == "interval"){
		dots$type = "o"
		listDefault(dots, "lwd", 2)
	} else if(style == "tube"){
		dots$type = "n"
		listDefault(dots, "lwd", 2)
	} else {
		dots$type = "p"
	}

	if(only.params){
	    return(list(x = x_value, y = estimate, ylim = my_ylim))
	}

	if(!add){

	    dots$axes = FALSE

	    # Nude graph
	    first.par = dots
	    first.par$type = "n"
	    do.call("plot", first.par)

	    if(grid){
	        listDefault(grid.par, "col", "gray")
	        listDefault(grid.par, "lty", 3)
	        listDefault(grid.par, "vert", TRUE)
	        listDefault(grid.par, "horiz", TRUE)

	        vert = grid.par$vert
	        horiz = grid.par$horiz
	        grid.par$vert = grid.par$horiz = NULL

	        if(horiz){
	            do.call("hgrid", grid.par)
	        }

	        if(vert){
	            do.call("vgrid", grid.par)
	        }
	    }

	    if(zero){
	        listDefault(zero.par, "lwd", 1)
	        listDefault(zero.par, "col", "black")
	        zero.par$h = 0
	        do.call("abline", zero.par)
	    }

	    if(ref.line){

	        if(is.na(ref_id) && !"v" %in% names(ref.line.par)){
	            stop("You can use the argument 'ref.line' only when interactions are provided and a reference is found. You can still draw vertical lines by using 'v' in argument 'ref.line.par'. Example: ref.line.par=list(v = ", round(x_value[floor(length(x_value)/2)]), ", col=2).")
	        }

	        listDefault(ref.line.par, "v", x_value[ref_id])
	        listDefault(ref.line.par, "lty", 2)
	        do.call("abline", ref.line.par)
	    }

	    box()
	    axis(2)
	    axis(1, at = x_at, labels = x_labels)

	}

	if(style == "bar" && join){
		# We join the dots

	    listDefault(join.par, "lwd", lwd)
	    listDefault(join.par, "col", pt.col)

	    join.par$x = dots$x
	    join.par$y = dots$y

		do.call("lines", join.par)
	}


	if(style == "interval"){
		# the "tube"
		lines(x_value, ci025, lty = ci.lty, lwd = ci.lwd, col = ci.col)
		lines(x_value, ci975, lty = ci.lty, lwd = ci.lwd, col = ci.col)
	} else if(style == "tube"){
		# Here we use shade area
		shade_area(ci025, ci975, x_value, col = "lightgrey", lty=0)

		# dots$axes = NULL
		# dots$type = "o"
		# do.call("lines", dots)
	} else {
		for(i in 1:n){
		    x = x_value

			# a) barre verticale
			segments(x0=x[i], y0=ci025[i], x1=x[i], y1=ci975[i], lwd = ci.lwd, col = ci.col, lty = ci.lty)

			# Formatting the bar width

			if(length(ci.width) > 1){
			    stop("The argument 'ci.width' must be of length 1.")
			}

			if(is.character(ci.width)){

			    width_nb = tryCatch(as.numeric(gsub("%", "", ci.width)), warning = function(x) x)
			    if(!is.numeric(width_nb)){
			        stop("The value of 'ci.width' is not valid. It should be equal either to a number, either to a percentage (e.g. ci.width=\"3%\").")
			    }

			    if(grepl("%", ci.width)){
			        total_width = diff(par("usr")[1:2])
			        ci.width = total_width * width_nb / 100
			    } else {
			        ci.width = width_nb
			    }
			}



			# b) toppings
			# Only if not reference
			if(ci975[i] != ci025[i]){
			    #  i) ci975
			    segments(x0=x[i]-ci.width, y0=ci975[i], x1=x[i]+ci.width, y1=ci975[i], lwd = ci.lwd, col = ci.col, lty = ci.lty)
			    #  ii) ci025
			    segments(x0=x[i]-ci.width, y0=ci025[i], x1=x[i]+ci.width, y1=ci025[i], lwd = ci.lwd, col = ci.col, lty = ci.lty)
			}

		}
	}

	# Last the points
	if(!add){
	    # now the points or lines
	    if(dots$type != "n"){
	        point.par = dots[c("x", "y", "type", "cex", "col", "pch", "lty", "lwd")]
	        point.par$pch = pt.pch
	        point.par$cex = pt.cex
	        point.par$col = pt.col
	        point.par = point.par[lengths(point.par) > 0]
	        do.call("lines", point.par)
	    }
	} else {
	    dots$pch = pt.pch
	    dots$cex = pt.cex
	    dots$col = pt.col
	    do.call("lines", dots)
	}

	res = list(x = x_value, y = estimate, ylim = my_ylim)
	return(invisible(res))
}


#' Treated and control sample descriptives
#'
#' This function shows the means and standard-deviations of several variables conditional on whether they are from the treated or the control group. The groups can further be split according to a pre/post variable. Results can be seamlessly be exported to Latex.
#'
#' @inheritParams esttex
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

            check_arg(indiv, "osf")
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

        fml = extract_pipe(fml_in)$fml
        pipe = extract_pipe(fml_in)$pipe

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
        check_arg(prepostnames, "characterVector", "Argument 'prepostnames' must be a character vector of length 2. REASON")
        if(length(prepostnames) != 2){
            stop("Argument 'prepostnames' must be a character vector of length 2. It is currenlty of length ", length(prepostnames), ".")
        }
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
            for(i in qui) res$vars[i] = escape_latex(dict[res$vars[i]], depth = 1)
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



#' Interact variables with factors
#'
#' Interacts a variable with another treated as a factor, and sets a reference
#'
#' @param var A vector.
#' @param fe A vector (of any type). Must be of the same length as \code{var}.
#' @param ref A single value that belongs to the interacted variable (\code{fe}). Can be missing.
#' @param confirm Logical, default is \code{FALSE}. If the factor variable has over 100 cases, you need to use \code{confirm = TRUE} to carry on.
#'
#' @return
#' It returns a matrix with number of rows the length of \code{var}. The number of columns is equal to the number of cases contained in \code{fe} minus the reference.
#'
#' @section Shorthand in \code{fixest} estimations:
#' In \code{fixest} estimations, instead of using \code{i(var, fe, ref)}, you can instead use the following writing \code{var::fe(ref)}.
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
#' x = rnorm(10)
#' y = rep(1:4, 3)[1:10]
#'
#' cbind(x, y, i(x, y, 1))
#'
#' #
#' # In fixest estimations
#' #
#'
#' # NOTA: in fixest estimations, i(var, fe, ref) is equivalent to var::fe(ref)
#'
#' data(base_did)
#' # We interact the variable 'period' with the variable 'treat'
#' est_did = feols(y ~ x1 + i(treat, period, 5) | id+period, base_did)
#'
#' # You could have used the following formula instead:
#' # y ~ x1 + treat::period(5) | id+period
#'
i = interact = function(var, fe, ref, confirm = FALSE){
    # Used to create interactions

    mc = match.call()

    var_name = deparse_long(mc$var)
    fe_name = deparse_long(mc$fe)

    if(length(var) != length(fe)){
        stop("The arguments 'var' and 'fe' must be of the same length (currently ", length(var), " vs ", length(fe), ").")
    }

    # The NAs
    is_na_fe = is.na(fe)

    if(is.factor(fe)){
        # we respect the fact that fe is a factor => we will keep its ordering
        is_na_fe = is.na(fe)
        fe_no_na = fe[!is_na_fe, drop = TRUE]
        items = levels(fe_no_na)
        fe_num = rep(NA, length(fe))
        fe_num[!is_na_fe] = as.vector(unclass(fe_no_na))
    } else {
        if(any(is_na_fe)){
            quf = quf_sorted(fe[!is_na_fe], addItem = TRUE)
            fe_num = rep(NA, length(fe))
            fe_num[!is_na_fe] = quf$x
        } else {
            quf = quf_sorted(fe, addItem = TRUE)
            fe_num = quf$x
        }

        items = quf$items
    }

    noRef = FALSE
    if(!missing(ref)){
        # Controls

        if(length(ref) > 1) stop("The argument 'ref' must be of length 1 (currenlty it is of length ", length(ref), ").")
        if(is.na(ref)) stop("The argument 'ref' cannot be NA.")

        if(!ref %in% items){
            stop("Argument 'ref' is not an element of the variable ", fe_name, ".")
        }

        ref = which(items == ref)
    } else {
        noRef = TRUE
        ref = 1
    }

    if(length(items) > 100 && !confirm){
        stop("You are interacting ", var_name, " with a variable containing over 100 different values (exactly ", length(items), "). To proceed please add the argument 'confirm=TRUE'. Note that if you do not need the standard-errors, it is much faster to include the interactions in the fixed-effects part of the formula using ", fe_name, "[[", var_name, "]]. See details on how to add varying slopes.")
    }

    res = model.matrix(~ -1 + fe_num, model.frame(~ -1 + fe_num, data.frame(fe_num = factor(fe_num)), na.action = na.pass))

    if(noRef){
        res = res * var
        colnames(res) = paste0(var_name, ":", fe_name, "::", items)
    } else {
        res = res[, -ref] * var
        colnames(res) = paste0(var_name, ":", fe_name, "::", items[-ref])
    }

    # We send the information on the reference
    opt = getOption("fixest_interaction_ref")
    if(is.null(opt)){
        is_ref = rep(FALSE, length(items))
        if(noRef == FALSE){
            is_ref[ref] = TRUE
        }

        opt = list(is_ref = is_ref, items = items, fe_type = class(fe))

        options("fixest_interaction_ref" = opt)
    }

    res
}

#' @rdname i
"interact"

#### ................. ####
#### Internal Funs     ####
####

results2formattedList = function(..., se, dof = getFixest_dof(), cluster, digits=4, fitstat, sdBelow=TRUE, dict = NULL, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, subtitles, title, yesNoFixef = c("Yes", "No"), keepFactors = FALSE, isTex = FALSE, useSummary, dots_call, powerBelow, interaction.combine, convergence, show_family, drop, order, file, fixef_sizes = FALSE, show_depvar=FALSE){
    # This function is the core of the functions esttable and esttex

    # for error handling => refers to the right function
    my_call = deparse(sys.calls()[[sys.nframe()-1]])[1] # call can have svl lines
    nmax = 40
    if(nchar(my_call) > nmax) my_call = paste0(substr(my_call, 1, nmax-1), "...")
    my_call = paste0("error in ", my_call, ":\n")


    check_arg(signifCode, "numericVectorGE0LE1", call_depth = 1)
    if(is.null(names(signifCode))){
        stop(my_call, "The argument 'signifCode' must be a NAMED vector. It currently has no names.", call. = FALSE)
    }
    signifCode = sort(signifCode)

    check_arg(yesNoFixef, "characterVector", message = "The argument 'yesNoFixef' must be a character vector of length 1 or 2. REASON", call_depth = 1)
    if(length(yesNoFixef) == 1){
        yesNoFixef = c(yesNoFixef, "")
    }

    if(length(yesNoFixef) != 2){
        stop(my_call, "The argument 'yesNoFixef' must be of length 2.", call. = FALSE)
    }

    # at the moment: only fixest allowed
    allowed_types = "fixest"

    # We get all the models
    dots <- list(...)

    # formatting the names of the models
    dots_names = names(dots_call)
    if(!is.null(dots_names)){

        for(i in 1:length(dots_call)){
            if(dots_names[i] != ""){
                dots_call[[i]] = dots_names[i]
            } else {
                dots_call[[i]] = deparse_long(dots_call[[i]])
            }
        }
    }

    n = length(dots)

    if(n == 0) stop(my_call, "Not any estimation as argument.", call. = FALSE)

    all_models = list()
    model_names = list()
    k = 1
    for(i in 1:n){
        di = dots[[i]]

        if(any(allowed_types %in% class(di))){
            all_models[[k]] = di
            if(any(class(dots_call[[i]]) %in% c("call", "name"))){
                model_names[[k]] = deparse_long(dots_call[[i]])
            } else {
                model_names[[k]] = as.character(dots_call[[i]])
            }

            k = k+1
        } else if(length(class(di))==1 && class(di)=="list"){
            # we get into this list to get the fixest objects
            types = sapply(di, class)
            qui = which(types %in% allowed_types)
            for(m in qui){
                all_models[[k]] = di[[m]]

                # handling names
                if(n > 1){
                    if(is.null(names(di)[m]) || names(di)[m]==""){
                        model_names[[k]] = paste0(dots_call[[i]], "[[", m, "]]")
                    } else {
                        model_names[[k]] = paste0(dots_call[[i]], "$", names(di)[m])
                    }
                } else {
                    model_names[[k]] = as.character(names(di)[m])
                }

                k = k+1
            }
        }

    }

    if(length(all_models)==0) stop(my_call, "Not any proper model (fixest) as argument!", call. = FALSE)

    n_models <- length(all_models)

    # formatting the names (if needed)
    alternative_names = paste0("model ", 1:n_models)
    who2replace = sapply(model_names, function(x) length(x) == 0 || x == "")
    model_names[who2replace] = alternative_names[who2replace]

    # we keep track of the SEs
    se_type_list = list()

    check_interaction_reorder = FALSE
    var_list <- var_reorder_list <- coef_list <- coef_below <- sd_below <- list()
    depvar_list <- obs_list <- fitstat_list <- list()
    r2_list <- aic_list <- bic_list <- loglik_list <- convergence_list <- list()
    sqCor_list = family_list = theta_list = list()

    # To take care of factors
    fe_names = c()
    is_fe = vector(mode = "list", n_models)
    nb_fe = vector(mode = "list", n_models) # the number of items per factor

    slope_names = c()
    slope_flag_list = vector(mode = "list", n_models)

    # if there are subtitles
    if(!missing(subtitles)){
        if(length(subtitles) != n_models){
            stop(my_call, "If argument 'subtitles' is provided, it must be of the same length as the number of models. Current lengths: ", length(subtitles), " vs ", n_models, " models.", call. = FALSE)
        } else {
            isSubtitles = TRUE
        }

        if(isTex){
            subtitles = escape_latex(subtitles, depth = 2)
        }

    } else {
        subtitles = NULL
        isSubtitles = FALSE
    }

    if(!is.null(dict) && isTex){
        dict = escape_latex(dict, depth = 2)
    }

    #
    # fitstat: which R2 to display?
    #

    if(missing(fitstat)){
        # Default values:
        #   - if all OLS: typical R2
        #   - if any non-OLS: pseudo R2 + squared cor.
        is_ols = sapply(all_models, function(x) deparse(x$call[[1]]) == "feols")

        if(all(is_ols)){
            if(any(sapply(all_models, function(x) "fixef_vars" %in% names(x)))){
                # means any FE model
                fitstat = c("r2", "wr2")
            } else {
                fitstat = c("r2", "ar2")
            }
        } else {
            fitstat = c("sq.cor", "pr2", "bic")
        }


    } else if(isFALSE(fitstat) || (length(fitstat) == 1 && fitstat == "")){
        fitstat = NULL
    } else if("formula" %in% class(fitstat)){
        check_arg(fitstat, "osf", "Argument 'fitstat' must be a one sided formula (or a character vector) containing 'aic', 'bic', 'll', or valid r2 types names (see function r2). REASON", call_depth = 1)
        fitstat = attr(terms(fitstat), "term.labels")
    } else {
        check_arg(fitstat, "characterVector", "Argument 'fitstat' must be a character vector (or a one sided formula) containing 'aic', 'bic', 'll', or valid r2 types names (see function r2). REASON", call_depth = 1)
    }

    # checking the types
    fitstat_type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2", "ll", "aic", "bic")
    fitstat = unique(fitstat)

    pblm = setdiff(fitstat, fitstat_type_allowed)
    if(length(pblm) > 0){
        stop(my_call, "Argument 'fitstat' must be a character vector (or a one sided formula) containing 'aic', 'bic', 'll', or valid r2 types names. ", enumerate_items(pblm, "is.quote"), " not valid (see function r2).", call. = FALSE)
    }

    fitstat_dict_tex = c("sq.cor"="Squared Correlation", r2="R$^2$", ar2="Adjusted R$^2$", pr2="Pseudo R$^2$", apr2="Adjusted Pseudo R$^2$", wr2="Within R$^2$", war2="Within Adjusted R$^2$", wpr2="Within Pseudo R$^2$", wapr2="Whithin Adjusted Pseudo R$^2$", aic = "AIC", bic = "BIC", ll = "Log-Likelihood")

    fitstat_dict_R = c("sq.cor"="Squared Corr.", r2="R2", ar2="Adjusted R2", pr2="Pseudo R2", apr2="Adj. Pseudo R2", wr2="Within R2", war2="Within Adj. R2", wpr2="Within Pseudo R2", wapr2="Whithin Adj. Pseudo R2", aic = "AIC", bic = "BIC", ll = "Log-Likelihood")

    fitstat_dict = fitstat_dict_R
    if(isTex) fitstat_dict = fitstat_dict_tex

    # end: fitstat

    for(m in 1:n_models){

        # If se or cluster is provided, we use summary
        if(useSummary){
            x = summary(all_models[[m]], se=se, cluster, dof = dof, nframes_up = 2)
        } else {
            # What do we do if se not provided?
            # we apply summary only to the ones that are not summaries
            x = all_models[[m]]
            if(!"cov.scaled" %in% names(x)){
                # not a summary => we apply summary to trigger default behavior
                x = summary(x, dof = dof)
            }

        }
        se_type_list[[m]] = attr(x$se, "type")

        # family
        family = x$family
        if(x$method %in% c("femlm", "feNmlm")){
            fam = switch(family, poisson = "Poisson", negbin = "Neg. Bin.", gaussian = "Gaussian", logit = "Logit")
        } else if(x$method %in% c("feglm", "feglm.fit")){
            if(family$family == "poisson" && family$link == "log"){
                fam = "Poisson"
            } else if(family$family == "binomial" && family$link == "logit"){
                fam = "Logit"
            } else if(family$family == "binomial" && family$link == "probit"){
                fam = "Probit"
            } else {
                # we try to give the greatest details ()
                fam = paste0(family$family, '("', family$link, '")')
            }
        } else if(x$method %in% c("feols", "feols.fit")){
            fam = "OLS"
        } else if(x$method == "fepois"){
            fam = "Poisson"
        } else if(x$method == "fenegbin"){
            fam = "Neg. Bin."
        }
        family_list[[m]] = fam


        # Negbin parameter
        theta = all_models[[m]]$theta
        theta_list[[m]] = ifelse(is.null(theta), "", numberFormatNormal(theta))

        # variable dependante:
        depvar <- gsub(" ", "", as.character(x$fml)[[2]])

        a <- x$coeftable
        if(!is.data.frame(a)){
            class(a) <- NULL
            a = as.data.frame(a)
        }

        # We drop the .theta coefficient
        if(x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && family == "negbin")){
            quiTheta = rownames(a) == ".theta"
            a = a[!quiTheta, ]
        }

        #
        # START: Formatting of the factors / FEs / Slopes
        #

        # on enleve les facteurs des variables a garder
        if(!keepFactors){
            fact = rownames(a)
            qui_drop = grepl("factor(", fact, fixed = TRUE)
            a = a[!qui_drop, , FALSE]
            b = fact[qui_drop]
            c = sapply(b, function(x) strsplit(x, "factor(", fixed=TRUE)[[1]][2])
            d = sapply(c, function(x) strsplit(x, ")", fixed=TRUE)[[1]][1])
            factor_var = unique(d)

            # Now the number of items per factor
            if(length(factor_var) == 0){
                nbItems = character(0)
            } else {
                nbItems = addCommas(sapply(factor_var, function(x) 1+sum(grepl(x, b))))
            }
        } else {
            factor_var = c()
            nbItems = character(0)
        }

        # now the normal FEs
        if(!is.null(x$fixef_terms)){
            terms_full = extract_fe_slope(x$fixef_terms)
            fixef_vars = terms_full$fixef_vars

            factor_var = c(factor_var, fixef_vars, recursive=TRUE)

            new_items = addCommas(as.vector(x$fixef_sizes[fixef_vars]))
            names(new_items) = fixef_vars

            nbItems = c(nbItems, new_items)
        } else if(!is.null(x$fixef_vars)){
            factor_var = c(factor_var, x$fixef_vars, recursive=TRUE)

            new_items = addCommas(as.vector(x$fixef_sizes))
            names(new_items) = names(x$fixef_sizes)

            nbItems = c(nbItems, new_items)
        }

        nb_fe[[m]] = nbItems

        # Formatting

        lFactor = rep(yesNoFixef[1], length(factor_var))
        names(lFactor) = factor_var
        is_fe[[m]] = lFactor

        fe_names = unique(c(fe_names, factor_var, recursive=TRUE))

        #
        # SLOPES
        #

        if(!is.null(x$fixef_terms)){
            terms_full = extract_fe_slope(x$fixef_terms)
            slope_fe = terms_full$slope_fe
            slope_vars = terms_full$slope_vars

            # we change the names right away
            slope_fe_name = slope_fe
            slope_vars_name = slope_vars
            if(!is.null(dict)){
                qui = which(slope_fe %in% names(dict))
                if(length(qui) > 0){
                    slope_fe_name[qui] = dict[slope_fe[qui]]
                }

                qui = which(slope_vars %in% names(dict))
                if(length(qui) > 0){
                    slope_vars_name[qui] = dict[slope_vars[qui]]
                }
            }

            slope_var_full = paste0(slope_vars_name, " (", slope_fe_name, ")")

        } else {
            slope_var_full = c()
        }

        slope_flag = rep(yesNoFixef[1], length(slope_var_full))
        names(slope_flag) = slope_var_full
        slope_flag_list[[m]] = slope_flag

        slope_names = unique(c(slope_names, slope_var_full, recursive = TRUE))

        #
        #   END: FE/slope formatting
        #

        # on enleve les espaces dans les noms de variables
        var <- c(gsub(" ", "", row.names(a)))
        # renaming => Tex only
        if(isTex){
            qui = var %in% names(dict)
            var[qui] = dict[var[qui]]
            tv = table(var)
            if(any(tv > 1)){
                value_pblm = names(tv)[tv > 1][1]
                var_pblm = c(gsub(" ", "", row.names(a)))[var == value_pblm]
                stop("Problematic value for argument 'dict': The variables ", enumerate_items(var_pblm, "quote"), " have all the same alias ('", value_pblm, "') in the same estimation. This makes no sense, please provide a separate alias for each.")
            }

            # if there are still interactions, we rename them
            new_var = var
            var_left = var[!qui]
            if(length(var_left) > 0 && any(grepl(":", var_left))){
                check_interaction_reorder = TRUE


                qui_inter = grepl(":", var_left)
                inter = strsplit(var_left[qui_inter], "(?<=[^:]):(?=[^:])", perl = TRUE)

                fun_rename = function(x){
                    # We put the factors on the right
                    qui_factor = grepl("::", x)
                    if(any(qui_factor)){
                        res = x[base::order(qui_factor)]
                        res = gsub("::", " $=$ ", res)
                    } else {
                        res = x
                    }

                    who = res %in% names(dict)

                    res[who] = dict[res[who]]
                    paste0(res, collapse = interaction.combine)
                }

                inter_named = sapply(inter, fun_rename)
                new_inter = sapply(inter, function(x) fun_rename(sort(x)))

                var[!qui][qui_inter] = new_inter
                new_var[!qui][qui_inter] = new_inter
            }

            var_reorder_list[[m]] <- new_var
        } else {
            # We reorder the interaction terms alphabetically
            new_var = var
            qui = grepl(":", new_var)
            if(any(qui)){
                check_interaction_reorder = TRUE
                inter = strsplit(new_var[qui], ":")
                new_inter = sapply(inter, function(x) paste0(sort(x), collapse = ":"))
                new_var[qui] = new_inter
            }
            var_reorder_list[[m]] <- new_var
        }

        if(isTex){
            coef = coefFormatLatex(a[, 1], digits = digits, power = abs(powerBelow))
            se_value = coefFormatLatex(a[, 2], digits = digits, power = abs(powerBelow))
        } else {
            coef = as.character(round(a[, 1], digits))
            se_value = as.character(myRound(a[, 2], digits))
        }

        if(isTex){
            pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(tex_star(names(signifCode)), ""))
        } else {
            pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
        }

        # If the coefficient is bounded, we supress the 'stars'
        isBounded = grepl("bounded", se_value)
        if(any(isBounded)){
            pval[isBounded] = ""
        }

        structured_coef = c(paste0(coef, pval, " (", se_value, ")"))

        # saving the infos
        var_list[[m]] <- var
        names(structured_coef) <- var
        coef_list[[m]] <- structured_coef
        if(sdBelow){
            cb = c(paste0(coef, pval))
            sb = c(paste0("(", se_value, ")"))
            names(cb) = names(sb) = var
            coef_below[[m]] = cb
            sd_below[[m]] = sb
        }

        # La depvar
        depvar_list[[m]] <- depvar

        #
        #  Fit statistics
        #

        # Pseudo-R2 // AIC // BIC // N
        n <- nobs(x)
        obs_list[[m]] <- n
        convergence_list[[m]] = ifelse(is.null(x$convStatus), TRUE, x$convStatus)

        K <- x$nparams

        if(length(fitstat) == 0){
            fitstat_list[[m]] = NA
        } else {
            fistat_format = list()

            fun_format = ifelse(isTex, numberFormatLatex, numberFormatNormal)

            if("aic" %in% fitstat) fistat_format[["aic"]] = fun_format(round(AIC(x), 3))
            if("bic" %in% fitstat) fistat_format[["bic"]] = fun_format(round(BIC(x), 3))
            if("ll" %in% fitstat) fistat_format[["ll"]] = fun_format(logLik(x))

            # regular r2s
            r2_type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2")
            if(any(fitstat %in% r2_type_allowed)){
                all_r2 = r2(x, intersect(fitstat, r2_type_allowed))
                for(r2_val in names(all_r2)){
                    nb = 5 - (r2_val == "sq.cor") * 2
                    fistat_format[[r2_val]] = round(as.vector(all_r2[r2_val]), nb)
                }
            }

            fitstat_list[[m]] = fistat_format[fitstat]
        }

    }

    if(check_interaction_reorder){
        if(length(unique(unlist(var_reorder_list))) < length(unique(unlist(var_list)))){
            var_list = var_reorder_list
            for(m in 1:length(var_list)){
                names(coef_list[[m]]) <- var_list[[m]]
            }
        }
    }


    if(length(fitstat) > 0){
        attr(fitstat_list, "format_names") = fitstat_dict[fitstat]
    }

    if(isTex){
        if(missing(title)){
            title = "no title"
        } else {
            title = escape_latex(title, 2)
        }
    } else {
        if(missing(title)){
            title = NULL
        }
    }


    if((missing(convergence) && any(convergence_list == FALSE)) || (!missing(convergence) && convergence)){
        convergence = TRUE
    } else {
        convergence = FALSE
    }

    if((!missing(show_family) && show_family) || (missing(show_family) && length(unique(family_list)) > 1)){
        family = TRUE
    } else {
        family = FALSE
    }

    if((missing(show_depvar) && length(unique(unlist(depvar_list))) > 1) || (!missing(show_depvar) && show_depvar)){
        depvar = TRUE
    } else {
        depvar = FALSE
    }

    if(missing(drop)) drop = NULL
    if(missing(order)) order = NULL
    if(missing(file)) file = NULL
    if(missing(label)) label = NULL

    res = list(se_type_list=se_type_list, var_list=var_list, coef_list=coef_list, coef_below=coef_below, sd_below=sd_below, depvar_list=depvar_list, obs_list=obs_list, convergence_list=convergence_list, fe_names=fe_names, is_fe=is_fe, nb_fe=nb_fe, slope_flag_list = slope_flag_list, slope_names=slope_names, useSummary=useSummary, model_names=model_names, family_list=family_list, theta_list=theta_list, fitstat_list=fitstat_list, subtitles=subtitles, isSubtitles=isSubtitles, title=title, convergence=convergence, family=family, drop=drop, order=order, file=file, label=label, sdBelow=sdBelow, signifCode=signifCode, fixef_sizes=fixef_sizes, depvar=depvar, useSummary=useSummary, dict=dict, yesNoFixef=yesNoFixef)

    return(res)
}

etable_internal_latex = function(info){
    # Internal function to display the latex table

    n_models = length(info$depvar_list)
    # Getting the information
    se_type_list = info$se_type_list
    var_list = info$var_list
    coef_list = info$coef_list
    coef_below = info$coef_below
    sd_below = info$sd_below
    depvar_list = info$depvar_list
    obs_list = info$obs_list
    convergence_list = info$convergence_list
    fe_names = info$fe_names
    is_fe = info$is_fe
    nb_fe = info$nb_fe
    slope_names = info$slope_names
    slope_flag_list = info$slope_flag_list
    family_list = info$family_list
    theta_list = info$theta_list
    fitstat_list = info$fitstat_list
    subtitles = info$subtitles
    isSubtitles = info$isSubtitles
    title = info$title
    label = info$label
    drop = info$drop
    order = info$order
    file = info$file
    family = info$family
    convergence = info$convergence
    sdBelow = info$sdBelow
    signifCode = info$signifCode
    fixef_sizes = info$fixef_sizes
    dict = info$dict
    yesNoFixef = info$yesNoFixef

    #
    # prompting the infos gathered
    #

    # Starting the table
    myTitle = title
    if(!is.null(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
    start_table = paste0("\\begin{table}[htbp]\\centering\n\\caption{",  myTitle, "}\n")
    end_table = "\\end{table}"

    # intro and outro Latex tabular
    myAmpLine = paste0(paste0(rep(" ", length(depvar_list)+1), collapse="&"), "\\tabularnewline\n")
    intro_latex <- paste0("\\begin{tabular}{l", paste0(rep("c", n_models), collapse=""), "}\n",
                          myAmpLine,
                          "\\hline\n",
                          "\\hline\n")

    outro_latex <- "\\end{tabular}\n"

    # 1st lines => dep vars
    depvar_list = c(depvar_list, recursive = TRUE)

    qui = depvar_list %in% names(dict)
    who = depvar_list[qui]
    depvar_list[qui] = dict[who]
    depvar_list = escape_latex(depvar_list, depth = 2)

    # We write the dependent variables properly, with multicolumn when necessary
    # to do that, we count the number of occurences of each variable (& we respect the order provided by the user)
    nb_multi = 1
    names_multi = depvar_list[1]

    if(n_models > 1){
        k = 1
        old_dep = depvar_list[1]
        for(i in 2:length(depvar_list)){
            if(depvar_list[i] == old_dep){
                nb_multi[k] = nb_multi[k] + 1
            } else {
                k = k + 1
                nb_multi[k] = 1
                names_multi[k] = old_dep = depvar_list[i]
            }
        }
    }

    # now the proper format
    first_line <- "Dependent Variables:"
    if(length(nb_multi) == 1) first_line = "Dependent Variable:"
    for(i in 1:length(nb_multi)){
        if(nb_multi[i] == 1){
            # no multi column
            first_line = paste0(first_line, "&", names_multi[i])
        } else {
            first_line = paste0(first_line, "&\\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
        }
    }
    first_line = paste0(first_line, "\\\\\n")

    # Model line
    model_line = paste0("Model:&", paste0("(", 1:n_models, ")", collapse = "&"), "\\\\\n")

    # a simple line with only "variables" written in the first cell
    variable_line = "\\hline\n\\emph{Variables}\\tabularnewline\n"


    # Coefficients,  the tricky part
    coef_lines <- list()
    all_vars <- unique(c(var_list, recursive=TRUE))

    # dropping some coefs
    all_vars = drop_apply(all_vars, drop, depth = 2)

    # ordering the coefs
    all_vars = order_apply(all_vars, order, depth = 2)

    # changing the names of the coefs
    aliasVars = all_vars
    names(aliasVars) = all_vars

    qui = all_vars %in% names(dict)
    who = aliasVars[qui]
    aliasVars[qui] = dict[who]
    aliasVars = escape_latex(aliasVars, depth = 2)

    coef_mat <- all_vars
    for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
    coef_mat[is.na(coef_mat)] <- "  "
    if(sdBelow){
        coef_lines = c()
        for(v in all_vars){
            myCoef = mySd= myLine = c()
            for(m in 1:n_models){
                myCoef = c(myCoef, coef_below[[m]][v])
                mySd = c(mySd, sd_below[[m]][v])
            }

            myCoef[is.na(myCoef)] = "  "
            mySd[is.na(mySd)] = "  "
            myCoef = paste0(aliasVars[v], "&", paste0(myCoef, collapse="&"))
            mySd = paste0("  &", paste0(mySd, collapse="&"))
            myLines = paste0(myCoef, "\\\\\n", mySd, "\\\\\n")
            coef_lines = c(coef_lines, myLines)
        }
        coef_lines = paste0(coef_lines, collapse="")
    } else {
        coef_lines = paste0(paste0(apply(coef_mat, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")
    }

    # Fixed-effects (if needed)
    if(length(fe_names) > 0){
        dumIntro = paste0("\\hline\n\\emph{Fixed-Effects}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")

        for(m in 1:n_models) {
            quoi = is_fe[[m]][fe_names]
            quoi[is.na(quoi)] = yesNoFixef[2]
            is_fe[[m]] = quoi

            # We do the same for the number of items
            quoi = nb_fe[[m]][fe_names]
            quoi[is.na(quoi)] = "--"
            nb_fe[[m]] = quoi
        }

        all_fe = matrix(c(is_fe, recursive = TRUE), nrow = length(fe_names))

        # We change the names of the FEs
        for(i in seq_along(fe_names)){
            fe = fe_names[i]

            if(fe %in% names(dict)){
                fe_names[i] = dict[fe]
            } else if(grepl("\\^", fe)){
                fe_split = strsplit(fe, "\\^")[[1]]
                who = fe_split %in% names(dict)
                fe_split[who] = dict[fe_split[who]]
                fe_names[i] = paste(fe_split, collapse = "$\\times$")
            } else {
                fe_names[i] = fe
            }
        }
        fe_names = escape_latex(fe_names)

        all_fe = cbind(fe_names, all_fe)
        factor_lines <- paste0(paste0(apply(all_fe, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

        # For the number of items
        all_nb_Factors = matrix(c(nb_fe, recursive=TRUE), nrow = length(fe_names))
        fe_names_nbItems = paste0("# ", fe_names)
        all_nb_Factors = cbind(fe_names_nbItems, all_nb_Factors)
        nb_factor_lines <- paste0(paste0(apply(all_nb_Factors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")


    } else {
        factor_lines = NULL
        dumIntro = NULL
    }

    # Slopes (if needed)
    if(length(slope_names) > 0){
        slope_intro = paste0("\\hline\n\\emph{Varying Slopes}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")

        # reformat the yes/no slope
        for(m in 1:n_models) {
            quoi = slope_flag_list[[m]][slope_names]
            quoi[is.na(quoi)] = yesNoFixef[2]
            slope_flag_list[[m]] = quoi
        }

        # Changing the slope names
        qui = slope_names %in% names(dict)
        who = slope_names[qui]
        slope_names[qui] = dict[who]
        slope_names = escape_latex(slope_names)

        # Matrix with yes/no information
        all_slopes = matrix(c(slope_flag_list, recursive = TRUE), nrow = length(slope_names))
        all_slopes = cbind(slope_names, all_slopes)
        slope_lines <- paste0(paste0(apply(all_slopes, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

    } else {
        slope_intro = NULL
        slope_lines = NULL
    }

    # Subtitles
    if(isSubtitles){
        info_subtitles = paste0("  & ", paste(subtitles, collapse="&"), "\\\\\n")
    } else {
        info_subtitles = ""
    }

    # Convergence information
    info_convergence = ""
    if(convergence){
        info_convergence = paste0("Convergence &", paste(convergence_list, collapse="&"), "\\\\\n")
    }

    info_theta <- paste0("Overdispersion& ", paste(theta_list, collapse="&"), "\\\\\n")

    # information on family
    if(family){
        info_family <- paste0("Family& ", paste(family_list, collapse="&"), "\\\\\n")
    } else {
        info_family = ""
    }


    # The standard errors
    isUniqueSD = length(unique(unlist(se_type_list))) == 1
    nb_col = length(obs_list) + 1
    sd_intro = paste0("\\multicolumn{", nb_col, "}{l}{\\emph{")
    if(isUniqueSD){
        my_se = unique(unlist(se_type_list)) # it comes from summary
        # every model has the same type of SE
        if(my_se == "Standard") my_se = "Normal"
        if(my_se == "White") my_se = "White-corrected"

        # Now we modify the names of the clusters if needed
        my_se = format_se_type_latex(my_se, dict)

        info_SD = paste0("\\hline\n\\hline\n", sd_intro, my_se, " standard-errors in parentheses.}}\\\\\n")
        info_SD = paste0(info_SD, sd_intro, "Signif Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
        info_muli_se = ""
    } else {
        all_se_type = sapply(se_type_list, format_se_type_latex, dict = dict, inline = TRUE)
        info_muli_se = paste0("Standard-Error type& ", paste(all_se_type, collapse = "&"), "\\\\\n")

        info_SD = paste0("\\hline\n\\hline\n", sd_intro, "Signif Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
    }

    # Information on number of items
    supplemental_info = ""

    if(!fixef_sizes) nb_factor_lines = ""
    if(all(theta_list == "")) info_theta = ""

    #
    # Fit statistics
    #

    fit_info = paste0("\\hline\n\\emph{Fit statistics}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
    fit_info = paste0(fit_info, "Observations& ", paste(addCommas(obs_list), collapse="&"), "\\\\\n")
    fit_info = paste0(fit_info, nb_factor_lines, info_convergence, info_muli_se)
    if(!identical(fitstat_list, NA)){

        fit_names = attr(fitstat_list, "format_names")
        nb = length(fit_names)
        for(fit_id in 1:nb){
            fit = sapply(fitstat_list, function(x) x[[fit_id]])
            fit[is.na(fit)] = "--"
            fit_info = paste0(fit_info, fit_names[fit_id], " & ", paste0(fit, collapse = "&"), "\\\\\n")
        }
    }

    res = c(supplemental_info, start_table, intro_latex, first_line, info_subtitles, model_line,
    info_family, variable_line, coef_lines, info_theta, dumIntro, factor_lines,
    slope_intro, slope_lines, fit_info, info_SD, outro_latex, end_table)

    res = res[nchar(res) > 0]

    return(res)
}

etable_internal_df = function(info){

    n_models = length(info$depvar_list)
    # Getting the information
    se_type_list = info$se_type_list
    var_list = info$var_list
    coef_list = info$coef_list
    coef_below = info$coef_below
    sd_below = info$sd_below
    depvar_list = info$depvar_list
    obs_list = info$obs_list
    convergence_list = info$convergence_list
    fe_names = info$fe_names
    is_fe = info$is_fe
    nb_fe = info$nb_fe
    slope_names = info$slope_names
    slope_flag_list = info$slope_flag_list
    family_list = info$family_list
    theta_list = info$theta_list
    fitstat_list = info$fitstat_list
    title = info$title
    label = info$label
    drop = info$drop
    order = info$order
    file = info$file
    family = info$family
    convergence = info$convergence
    sdBelow = info$sdBelow
    signifCode = info$signifCode
    fixef_sizes = info$fixef_sizes
    depvar = info$depvar
    useSummary = info$useSummary
    model_names = info$model_names

    # naming differences
    titles = info$subtitles
    isTitles = info$isSubtitles

    # The coefficients

    all_vars <- unique(c(var_list, recursive=TRUE))

    # dropping some coefs
    all_vars = drop_apply(all_vars, drop, depth = 2)

    # ordering the coefs
    all_vars = order_apply(all_vars, order, depth = 2)

    coef_mat <- all_vars
    for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
    coef_mat[is.na(coef_mat)] <- "  "
    res = coef_mat

    if("Neg. Bin." %in% family_list){
        theta_line = c("Overdispersion:", unlist(theta_list))
        res = rbind(res, theta_line)
    }

    # The line with the dependent variable => defined here to get the width
    preamble = c()
    dep_width = 0
    if(depvar){
        preamble = rbind(c("Dependent Var.:", unlist(depvar_list)), preamble)
        dep_width = nchar(as.vector(preamble))
    }

    # Used to draw a line
    myLine = "______________________________________"
    longueur = apply(res, 2, function(x) max(nchar(as.character(x))))
    longueur = pmax(dep_width, longueur)
    theLine = sapply(longueur, function(x) sprintf("%.*s", x, myLine))
    theLine[1] = sprintf("%.*s", max(nchar(theLine[1]), 19), myLine)

    # The FEs
    if(length(fe_names)>0){

        for(m in 1:n_models) {
            quoi = is_fe[[m]][fe_names]
            quoi[is.na(quoi)] = "No"
            is_fe[[m]] = quoi
        }
        all_fe = matrix(c(is_fe, recursive=TRUE), nrow = length(fe_names))
        all_fe = cbind(fe_names, all_fe)
        factor_lines <- paste0(paste0(apply(all_fe, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

        myLine = "-------------------------------"

        res = rbind(res, c("Fixed-Effects:", sprintf("%.*s", longueur[-1], myLine)))
        factmat = matrix(c(strsplit(strsplit(factor_lines, "\n")[[1]], "&"), recursive = TRUE), ncol=n_models+1, byrow=TRUE)
        factmat[, ncol(factmat)]=gsub("\\", "", factmat[, ncol(factmat)], fixed = TRUE)
        res = rbind(res, factmat)
    }

    # The slopes
    if(length(slope_names) > 0){
        # reformatting the yes/no
        for(m in 1:n_models) {
            quoi = slope_flag_list[[m]][slope_names]
            quoi[is.na(quoi)] = "No"
            slope_flag_list[[m]] = quoi
        }
        all_slopes = matrix(c(slope_flag_list, recursive=TRUE), nrow = length(slope_names))
        all_slopes = cbind(slope_names, all_slopes)
        slope_lines <- paste0(paste0(apply(all_slopes, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

        myLine = "-------------------------------"

        res = rbind(res, c("Varying Slopes:", sprintf("%.*s", longueur[-1], myLine)))
        slope_mat = matrix(c(strsplit(strsplit(slope_lines, "\n")[[1]], "&"), recursive = TRUE), ncol=n_models+1, byrow = TRUE)
        slope_mat[, ncol(slope_mat)] = gsub("\\", "", slope_mat[, ncol(slope_mat)], fixed = TRUE)
        res = rbind(res, slope_mat)

    }


    # preamble created before because used to set the width
    if(length(preamble) > 0){
        # preamble = rbind(preamble, c("  ", theLine[-1]))
        preamble = rbind(preamble, rep("   ", length(theLine)))
        res <- rbind(preamble, res)
    }

    res <- rbind(res, theLine)

    # the line with the families
    if(family){
        # preamble = rbind(c("Family:", unlist(family_list)), preamble)
        res = rbind(res, c("Family", unlist(family_list)))
    }

    res <- rbind(res, c("Observations", addCommas(obs_list)))
    if(!useSummary || !any(grepl("\\(", unlist(se_type_list)))){
        se_type_format = c()
        for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]])
        res <- rbind(res, c("S.E. type", c(se_type_format, recursive = TRUE)))
    } else {
        main_type = gsub(" \\(.+", "", se_type_list[[1]])
        se_type_format = c()
        for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]], by = TRUE)
        res <- rbind(res, c(paste0("SE type: ", main_type), c(se_type_format, recursive = TRUE)))
    }

    # convergence status
    if(convergence){
        res <- rbind(res, c("Convergence", c(convergence_list, recursive = TRUE)))
    }

    #
    # Fit statistics
    #

    if(!identical(fitstat_list, NA)){

        fit_names = attr(fitstat_list, "format_names")
        nb = length(fit_names)
        for(fit_id in 1:nb){
            fit = sapply(fitstat_list, function(x) x[[fit_id]])
            fit[is.na(fit)] = "--"
            res <- rbind(res, c(fit_names[fit_id], fit))
        }
    }

    # if titles
    if(isTitles){
        modelNames = titles
    } else {
        # modelNames = paste0("model ", 1:n_models)
        modelNames = model_names
    }

    # we shorten the model names to fit the width
    for(m in 1:n_models) modelNames[m] = charShorten(modelNames[[m]], longueur[[1+m]])

    res <- as.data.frame(res)
    names(res) <- c("variables", modelNames)
    row.names(res) = res$variables
    res$variables = NULL

    # We rename theta when NB is used
    quiTheta = which(row.names(res) == ".theta")
    row.names(res)[quiTheta] = "Dispersion Parameter"

    if(!is.null(file)){
        sink(file = file, append = !replace)
        on.exit(sink())

        print(res)

        return(invisible(res))
    } else {
        return(res)
    }

}

myPrintCoefTable = function(coeftable, lastLine = ""){
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

    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

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


vcovClust <- function (cluster, myBread, scores, dof=FALSE, K, do.unclass=TRUE){
    # Internal function: no need for controls, they come beforehand
    # - cluster: the vector of dummies
    # - myBread: original vcov
    # - scores
    # Note: if length(unique(cluster)) == n (i.e. White correction), then the dof are such that vcovClust is equivalent to vcovHC(res, type="HC1")
    # Source: http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
    #         Cameron & Miller -- A Practitioner’s Guide to Cluster-Robust Inference

    n <- NROW(scores)
    if(missing(K)) K <- NCOL(scores)

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

prepare_matrix = function(fml, base){
    # This function is way faster than model.matrix but does not accept factors
    # The argument fml **MUST** not have factors!

    rhs = fml[c(1,3)]
    t = terms(rhs, data = base)

    all_var_names = attr(t, "term.labels")
    all_vars = gsub(":", "*", all_var_names)

    # Forming the call
    if(attr(t, "intercept") == 1){
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


fixest_model_matrix = function(fml, data){
    # This functions takes in the formula of the linear part and the
    # data
    # It reformulates the formula (ie with lags and interactions)
    # then either apply a model.matrix
    # either applies an evaluation (which can be faster)

    # Modify the formula to add interactions
    if(grepl("::", deparse_long(fml[[3]]))){
        fml = interact_fml(fml)
    }

    # Evaluation

    # we look at whether there are factor-like variables to be evaluated
    # if there is factors => model.matrix
    dataNames = names(data)
    linear.varnames = all.vars(fml[[3]])
    types = sapply(data[, dataNames %in% linear.varnames, FALSE], class)
    if(length(types) == 0 || grepl("factor", deparse_long(fml)) || any(types %in% c("character", "factor"))){
        useModel.matrix = TRUE
    } else {
        useModel.matrix = FALSE
    }

    if(useModel.matrix){
        # to catch the NAs, model.frame needs to be used....
        linear.mat = stats::model.matrix(fml, stats::model.frame(fml, data, na.action=na.pass))
    } else {
        linear.mat = prepare_matrix(fml, data)
    }

    if(any(grepl("^(i|interact)\\(", colnames(linear.mat)))){
        # the following is needed (later in fixest_env => means there are factors)
        useModel.matrix = TRUE

        # we change the names
        new_names = colnames(linear.mat)
        all_terms = attr(terms(fml), "term.labels")
        terms_inter = all_terms[grepl("^(i|interact)\\(", all_terms)]

        for(pattern in terms_inter){
            new_names = gsub(pattern, "", new_names, fixed = TRUE)
        }

        colnames(linear.mat) = new_names
    }

    attr(linear.mat, "useModel.matrix") = useModel.matrix

    linear.mat
}

terms_fixef = function(fml){
    # separate all terms of fml into fixed effects ans varying slopes
    # fml: one sided formula
    # can also be a vector of terms

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
    if(grepl("\\^", fml_char)){
        fml_char_new = gsub("\\^", "_impossible_var_name_", fml_char)
        fml = as.formula(paste0("~", fml_char_new))
    }


    t = try(terms(fml))

    # if error, we send it back
    if("try-error" %in% class(t)){
        t = gsub("_impossible_var_name_", "^", t)
        return(t)
    }

    my_vars = attr(t, "term.labels")

    # Further error checking (misuse of varying slopes)
    if(any(grepl("\\]", my_vars))){
        var2check = my_vars[grepl("\\]", my_vars)]
        var2check_double = var2check[grepl("\\]\\]", var2check)]
        var2check_single = var2check[!grepl("\\]\\]", var2check)]
        if(length(var2check_double) > 0){
            qui = !grepl("\\]\\]$", var2check_double) | !grepl("\\[\\[", var2check_double) | lengths(strsplit(var2check_double, "\\[")) != 3
            if(any(qui)){
                item_pblm = var2check_double[qui]
                msg = paste0("Square bracket are special characters use **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm), ").")
                class(msg) = "try-error"
                return(msg)
            }
        }

        if(length(var2check_single) > 0){
            qui = !grepl("\\]$", var2check_single) | lengths(strsplit(var2check_single, "\\[")) != 2
            if(any(qui)){
                item_pblm = var2check_single[qui]
                msg = paste0("Square bracket are special characters use **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm), ").")
                class(msg) = "try-error"
                return(msg)
            }
        }

    }

    qui_slope = grepl("\\[", my_vars)
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
    my_vars = gsub("_impossible_var_name_", "^", my_vars)

    res = list(fml_terms = my_vars, fe_vars = gsub("\\[.+", "", my_vars))
    res$slope_flag = grepl("\\[", my_vars)
    res$slope_vars = rep(NA, length(my_vars))
    res$slope_vars[res$slope_flag] = gsub(".+\\[|\\]", "", my_vars[res$slope_flag])

    res
}

only_slope = function(fe_terms){
    # returns a logical vector of whether a FE is only
    # related to a slope

    new_terms = terms_fixef(fe_terms)
    fixef_vars = unique(new_terms$fe_vars)
    res = tapply(!new_terms$slope_flag, new_terms$fe_vars, sum)[fixef_vars] == 0
    res
}

prepare_df = function(vars, base, fastCombine = NA){
    # vars: vector of variables to evaluate

    # we drop NAs and make it unique
    vars = unique(vars[!is.na(vars)])
    all_var_names = vars

    do_combine = !is.na(fastCombine)

    changeNames = FALSE
    if(do_combine && any(grepl("\\^", vars))){
        # special indicator to combine factors
        # ^ is a special character: only used to combine variables!!!

        fun2combine = ifelse(fastCombine, "combine_clusters_fast", "combine_clusters")

        vars_new = gsub("([[:alpha:]_\\.][[:alnum:]_\\.]*(\\^[[:alpha:]_\\.][[:alnum:]_\\.]*)+)",
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
            qui = grepl("combine_clusters", all_var_names)
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

    return(index)
}

combine_clusters = function(...){
    # This functions creates a new cluster from several clusters
    # basically: paste(cluster1, cluster2, ... etc, sep = "_")

    cluster = list(...)
    Q = length(cluster)

    # We just paste
    myDots = cluster
    myDots$sep = "_"
    index = do.call("paste", myDots)

    return(index)
}



interact_fml = function(fml){
    # The formula is simple (should contain only the RHS)
    # I use Formula for robustness

    x_FML = Formula(fml)

    x = attr(terms(formula(x_FML, lhs = 1, rhs = 1)), "term.labels")

    if (!any(grepl("[^:]::[^:]", x))){
        return(fml)
    } else {

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
    }

    lhs_fml = deparse_long(x_FML[[2]])
    rhs_fml = paste(x, collapse = "+")

    as.formula(paste0(lhs_fml, "~", rhs_fml))
}

interact_control = function(ref, confirm = FALSE){
    # Internal call
    # used to contral the call to interact is valid
    check_arg(confirm, "singleLogical")
    if(length(ref) > 1) stop("Argument 'ref' must be of length 1.")
    mc = match.call()

    res = c()
    if("ref" %in% names(mc)) res = paste0("ref = ", deparse_long(mc$ref))
    if("confirm" %in% names(mc)) res = c(res, paste0("confirm = ", confirm))

    res
}

add2fml <- function(fml, x){
    #
    stopifnot(is.character(x))

    my_call = parse(text = paste0("update(fml, .~.+", paste(x, collapse="+"), ")"))
    res  = eval(my_call)

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

#### ................. ####
#### Small Utilities ####
####

escape_all = function(x){
    # we escape all
    res = gsub("((?<=[^\\\\])|(?<=^))(\\$|_|%|&|\\^)", "\\\\\\2", x, perl = TRUE)
    res
}

escape_latex = function(x_all, depth = 0, noArg = FALSE){
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
                stop_depth(depth = depth, my_arg, "here are ", length(dollars), " dollar signs in the following character string:\n", x, "\nIt will raise a Latex error (which '$' means equation? which means dollar-sign?): if you want to use a regular dollar sign, please escape it like that: \\\\$.")
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

		if(exponent > 0){
			return(sprintf("%.*f", max(1, digits - abs(exponent)), x))
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

    if(is.na(x)) return("NA")
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

drop_apply = function(x, drop = NULL, depth = 1){

    if(missing(drop) || length(drop) == 0){
        return(x)
    }

    check_arg(drop, "characterVector", call_depth = depth, message = "The arg. 'drop' must be a vector of regular expressions (see help(regex)). REASON")

    res = x

    for(var2drop in drop){
        if(grepl("^!", var2drop)){
            res = res[grepl(substr(var2drop, 2, nchar(var2drop)), res)]
        } else {
            res = res[!grepl(var2drop, res)]
        }
    }

    res
}

order_apply = function(x, order = NULL, depth = 1){

    if(missing(order) || length(order) == 0){
        return(x)
    }

    check_arg(order, "characterVector", call_depth = depth, message = "The arg. 'order' must be a vector of regular expressions (see help(regex)). REASON")

    res = x

    for(var2order in rev(order)){
        if(grepl("^!", var2order)){
            who = !grepl(substr(var2order, 2, nchar(var2order)), res)
            res = c(res[who], res[!who])
        } else {
            who = grepl(var2order, res)
            res = c(res[who], res[!who])
        }
    }

    res
}

n_times = function(n){

    if(n <= 4){
        res = switch(n, "1"="once", "2"="twice", "3"="three times", "4"="four times")
    } else {
        res = paste0(n, " times")
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

quf_sorted = function(x, addItem = FALSE){
    # Same as QUF, but items are sorted

    quoi = quickUnclassFactor(x, TRUE)

    new_order = order(quoi$items)
    order_new_order = order(new_order)
    x_uf = order_new_order[quoi$x]

    if(addItem){
        res = list(x = x_uf, items = quoi$items[new_order])
        return(res)
    } else {
        return(x_uf)
    }
}

quickUnclassFactor = function(x, addItem = FALSE){
	# does as unclass(as.factor(x))
	# but waaaaay quicker

	if(!is.numeric(x)){
		# level and unclass is much slower
		x = as.character(x)
	}

    res = cpp_quf_gnl(x)

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
        return(x)
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
        se_formatted = paste0(nway, ": ", fe_format)
    } else {
        se_formatted = paste0(nway, " (", fe_format, ")")
    }

    se_formatted
}

tex_star = function(x){
    qui = nchar(x) > 0
    x[qui] = paste0("$^{", x[qui], "}$")
    x
}


extract_pipe = function(fml){
    # We extract the elements after the pipe

    FML = Formula::Formula(fml)
    n_fml = length(FML)
    n_rhs = n_fml[2]

    if(n_rhs == 1){
        fml_new = formula(FML, lhs = n_fml[1], rhs = 1)
        pipe = NULL
    } else if(n_rhs == 2){
        fml_new = formula(FML, lhs = 1, rhs = 1)
        pipe = as.expression(formula(FML, lhs = 0, rhs = 2)[[2]])
    } else {
        stop("fml must be at *most* a two part formula (currently it is ", n_rhs, " parts).")
    }

    list(fml=fml_new, pipe=pipe)
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
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
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

	all_BIC = c(-2*logLik(object) + 2*object$nparams*log(nobs(object)), otherBIC)

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
		resid = residuals(object)
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
#' Note that if the model has been estimated with clusters, to obtain the cluster coefficients, you need to use the function \code{\link[fixest]{fixef.fixest}}.
#'
#' @return
#' This function returns a named numeric vector.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{confint.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{esttable}}, \code{\link[fixest]{esttex}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # the coefficients of the variables:
#' coef(res)
#'
#' # the cluster coefficients:
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
#' # simple estimation on iris data, clustering by "Species"
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
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = "type")
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

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
	if(!is.null(object$obsRemoved)){
	    # we add NA values
	    tmp = rep(NA, object$nobs + length(object$obsRemoved))
	    tmp[-object$obsRemoved] = res
	    res = tmp
	}

	res
}

#' @rdname fitted.fixest
"fitted.values.fixest"

#' Extracts residuals from a \code{fixest} object
#'
#' This function extracts residuals from a fitted model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams nobs.fixest
#'
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
#' # simple estimation on iris data, clustering by "Species"
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we plot the residuals
#' plot(resid(res_poisson))
#'
resid.fixest = residuals.fixest = function(object, ...){
	object$residuals
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
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = "type")
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

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

	    if(!is.null(object$obsRemoved)){
	        # we add NA values
	        tmp = rep(NA, object$nobs + length(object$obsRemoved))
	        tmp[-object$obsRemoved] = res
	        res = tmp
	    }

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
	if(n == (object$nobs + length(object$obsRemoved)) && deparse_long(object$call$data) == deparse_long(mc$newdata) && getFixest_notes()){
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

		# Extraction of the clusters
		id_cluster = list()
		for(i in 1:n_cluster){
			# checking if the variable is in the newdata
		    fe_var = fixef_vars[i]
			variable = all.vars(parse(text = fe_var))
			isNotHere = !variable %in% names(newdata)
			if(any(isNotHere)){
				stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a cluster variable).")
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

			# Obtaining the unclassed vector of clusters
			cluster_current = eval(parse(text = fe_var), newdata)

			cluster_current_num = unclass(factor(cluster_current, levels = fixef_values_possible))
			id_cluster[[i]] = cluster_current_num
		}

		names(id_cluster) = fixef_vars

		# Value of the cluster coefficients
		cluster_coef = fixef(object)

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
	rhs_fml = formula(Formula(fml), lhs = 1, rhs = 1)
	if(grepl("[^:]::[^:]", deparse_long(rhs_fml[[3]]))){
	    new_fml = interact_fml(rhs_fml)
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
		# matrix_linear = stats::model.matrix(rhs_fml, stats::model.frame(rhs_fml, newdata, na.action=na.pass))
		matrix_linear = try(fixest_model_matrix(rhs_fml, newdata), silent = TRUE)
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
#' @param ... Other arguments to be passed to \code{\link[fixest]{summary.fixest}}.
#'
#' The computation of the VCOV matrix is first done in \code{\link[fixest]{summary.fixest}}.
#'
#' @return
#' It returns a \eqn{N\times N} square matrix where \eqn{N} is the number of variables of the fitted model.
#' This matrix has an attribute \dQuote{type} specifying how this variance/covariance matrix has been computed (i.e. was it created using White correction, or was it clustered along a specific factor, etc).
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
#' # "white" VCOV
#' vcov(est_pois, se = "white")
#'
#' # "clustered" VCOV (with respect to the Product factor)
#' vcov(est_pois, se = "cluster", cluster = trade$Product)
#' # another way to make the same request:
#' # note that previously arg. se was optional since deduced from arg. cluster
#' vcov(est_pois, cluster = "Product")
#' # yet another way:
#' vcov(est_pois, cluster = ~Product)
#'
#' # Another estimation without cluster:
#' est_pois_simple = femlm(Euros ~ log(dist_km) + log(Year), trade)
#'
#' # We can still get the clustered VCOV,
#' # but we need to give the argument cluster:
#' vcov(est_pois_simple, cluster = ~Product)
#'
#'
vcov.fixest = function(object, se, cluster, dof = getFixest_dof(), forceCovariance = FALSE, keepBounded = FALSE, ...){
	# computes the clustered vcov

    if(!"nframes_up" %in% names(match.call())){
        # condition means client call
        any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                      suggest_args = c("se", "cluster", "dof"))
        if(any_invalid){
            warning(attr(any_invalid, "msg"))
        }
    }

	if(!is.null(object$onlyFixef)){
		# means that the estimation is done without variables
		stop("No explanatory variable was used: vcov() cannot be applied.")
	}

	# Default behavior se:
	suffix = ""
	if(missnull(se)){
		if(missing(cluster)){
			if("fixef_vars" %in% names(object)){
				se = "cluster"
			} else {
				se = "standard"
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

			} else if(length(cluster) %in% (nobs(object) + 0:length(object$obsRemoved))){
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
		stop("Argument 'se' must be a character scalar equal to: 'standard', 'white', 'cluster', 'twoway', 'threeway' or 'fourway'.")
    }

	se.val = NULL
	try(se.val <- match.arg(se, c("standard", "white", "cluster", "twoway", "threeway", "fourway", "1", "2", "3", "4")), silent = TRUE)
	if(is.null(se.val)){
		stop("Invalid argument 'se'. It should be equal to one of 'standard', 'white', 'cluster', 'twoway', 'threeway' or 'fourway'.")
	}

	dots = list()
	if(is.null(dots$nframes_up)){
		nframes_up = 1
	} else {
		nframes_up = dots$nframes_up + 1
	}

	# check the dof
	if(!"dof.type" %in% class(dof)){
	    stop("The argument 'dof.type' must be an object created by the function dof().")
	} else {
	    dof.fixef = dof$fixef
	    is_exact = dof$exact
	    is_cluster = dof$cluster
	}

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
	isSlope = FALSE
	if(!is.null(object$fixef_terms)){
	    isSlope = TRUE
	    fixef_sizes_ok = object$fixef_sizes

	    size_to_drop = only_slope(object$fixef_terms)[names(fixef_sizes_ok)]

	    fixef_sizes_ok[size_to_drop] = 0
	    n_fe_ok = sum(fixef_sizes_ok > 0)
	} else {
	    fixef_sizes_ok = object$fixef_sizes
	}

	# How do we choose K? => argument dof

	if(dof.fixef == "false"){
	    # we do it with "minus" because of only slopes
	    K = object$nparams
	    if(n_fe_ok > 0){
	        K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1))
	    }
	} else if(dof.fixef == "true" || se.val %in% c("standard", "white")){
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
			VCOV_raw = object$cov.unscaled / (n / (n - object$nparams))
		}
	} else {
		VCOV_raw = object$cov.unscaled
	}


	# information on the variable used for the clustering
	type_info = ""

	is_nested = c()
	if(anyNA(VCOV_raw)){

		if(!forceCovariance){
		    last_warn = getOption("fixest_last_warning")
		    if(is.null(last_warn) || (proc.time() - last_warn)[3] > 1){
		        warning("Standard errors are NA because of likely presence of collinearity. Use function collinearity() to detect collinearity problems.", call. = FALSE)
		    }

			return(VCOV_raw)
		} else {
			VCOV_raw_forced = MASS::ginv(object$hessian)
			if(anyNA(VCOV_raw_forced)) {
				stop("The covariance matrix could not be 'forced'.")
			}

			object$cov.unscaled = VCOV_raw_forced
			return(vcov(object, se=se.val, cluster=cluster, dof=dof))
		}

	} else if(se.val == "standard"){

	    correction.dof = (n - 1) / (n - K)
		vcov = VCOV_raw * correction.dof

	} else if(se.val == "white"){

	    correction.dof = (n - 1) / (n - K)
		vcov = crossprod(myScore %*% VCOV_raw) * correction.dof

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
						if(NROW(data) != (object$nobs + length(object$obsRemoved))){
							stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(var2fetch, "s"), " in the original dataset (", deparse_long(dataName), "), yet the dataset doesn't have the same number of observations as was used in the estimation (", NROW(data), " instead of ", object$nobs + length(object$obsRemoved), ").", suffix)
						}

						if(length(object$obsRemoved)){
							data = data[-object$obsRemoved, all_vars, drop = FALSE]
						} else {
							data = data[, all_vars, drop = FALSE]
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
								value = data[[cname]]

								if(length(object$obsRemoved) > 0){
									value = value
								}

								cluster[[i]] = value
							} else {
								# combination
								if(grepl("^", cname, fixed = TRUE)){
									value_text = gsub("\\^", ", ", cname)
									value_text = paste0("combine_clusters_fast(", value_text, ")")
								}

								value_call = parse(text = value_text)
								value = eval(value_call, data)

								if(length(object$obsRemoved) > 0){
									value = value
								}

								cluster[[i]] = value
							}
						}

					}

				}
			} else if(length(cluster) == nway && is.numeric(cluster)){
			    # You can use a number to tell which cluster to use

			    if(length(object$fixef_vars) == 0){
			        stop("You can use an integer in the argument 'cluster' only when there have been fixed-effects in the estimation. Currenlty this is not the case. Alternatively, arg. 'cluster' can be a formula, a vector of variables or a list of vectors.")
			    }

			    if(!all(cluster %% 1 == 0) || any(cluster < 1 || cluster > 4)){
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
			if(n_per_cluster[1] == (object$nobs + length(object$obsRemoved))){
				# We modify the clusters
				for(i in 1:nway) cluster[[i]] = cluster[[i]][-object$obsRemoved]
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
		if(dof.fixef == "nested" && n_fe_ok >= 1){
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

				vcov = vcov + (-1)**(i+1) * vcovClust(index, VCOV_raw, myScore, dof = is_cluster, K, do.unclass=FALSE)
				vcov = vcov * ((n - 1) / (n - K))
			}
		}
	}

	if(any(diag(vcov)<0)){
		warning("Some variances are negative (likely problem of collinearity).")
	}

	sd.dict = c("standard" = "Standard", "white"="White", "cluster"="Clustered", "twoway"="Two-way", "threeway"="Three-way", "fourway"="Four-way")

	attr(vcov, "type") = paste0(as.vector(sd.dict[se.val]), type_info)
	attr(vcov, "dof.type") = paste0("dof(fixef = \"", dof.fixef, "\", exact = ", is_exact, ", cluster = ", is_cluster, ")")
	attr(vcov, "dof.K") = K

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
#' # We estimate the effect of distance on trade (with 3 cluster effects)
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
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  dots_args = c("exact_dof", "forceCovariance", "keepBounded"),
                                  suggest_args = c("parm", "level", "se", "cluster"))
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

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
#' @param fml.update Changes to be made to the original argument \code{fml}. See more information on \code{\link[stats]{update.formula}}. You can add/withdraw both variables and clusters. E.g. \code{. ~ . + x2 | . + z2} would add the variable \code{x2} and the cluster \code{z2} to the former estimation.
#' @param nframes (Advanced users.) Defaults to 1. Number of frames up the stack where to perform the evaluation of the updated call. By default, this is the parent frame.
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
#' # we add another cluster: "Product"
#' est_3 <- update(est_2, . ~ . | . + Product)
#'
#' # we remove the cluster "Origin" and the variable log(dist_km)
#' est_4 <- update(est_3, . ~ . - log(dist_km) | . - Origin)
#'
#' # Quick look at the 4 estimations
#' esttable(est_pois, est_2, est_3, est_4)
#'
update.fixest = function(object, fml.update, nframes = 1, ...){
	# Update method
	# fml.update: update the formula
	# If 1) SAME DATA and 2) SAME dep.var, then we make initialisation


	if(missing(fml.update)){
		fml.update = . ~ .
	} else {
		if(!"formula" %in% class(fml.update)){
			stop("The argument 'fml.update' is required.")
		}
	}

    if(isTRUE(object$fromFit)){
        stop("update method not available for fixest estimations obtained from fit methods.")
    }

    if(!isScalar(nframes) || nframes < 1 || nframes %% 1 != 0){
        stop("Argument 'nframes' must be a single integer greater than, or equal to, 1.")
    }

	FML = Formula(fml.update)
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
	fml = update(fml_old, formula(FML, lhs = 1, rhs = 1))
	fml_char = as.character(fml)

	useInit = TRUE
	if(fml[[2]] != fml_old[[2]]){
		# means different dependent variables
		# 	=> initialisation with past parameters is useless
		useInit = FALSE
	}

	# Family information
	if(!is.null(dots$family)){
		if(object$method %in% c("femlm", "feNmlm", "fepois", "fenegbin")){
			family_new = match.arg(dots$family, c("poisson", "negbin", "gaussian", "logit"))
			if(family_new != object$family){
				# if different families: initialisation is useless
				useInit = FALSE
			}
		} else if(object$method == "feglm"){
			useInit = FALSE # if the user uses argument family it means it's a different one
		} else {
			stop("'family' is not an argument of function feols().")
		}
	}

	if(!is.null(dots$na.rm) && dots$na.rm && !(!is.null(object$call$na.rm) && object$call$na.rm)){
		# Too complicated to initialize with na.rm
		# I would have to make tons of controls, and it would work
		# only in some cases...
		useInit = FALSE
	}

	if(!is.null(object$fixef_terms)){
	    # In case of slopes => no init
	    useInit = FALSE
	}

	#
	# II) evaluation data
	#

	# We find out if it is the same data
	#	=> only to find out if init is useful
	if(!is.null(dots$data)){
		useInit = FALSE
	}

	if(useInit){
		# we test if we can use the initialisation of the parameters

		# We are here ONLY if the data needs to be evaluated
		data = NULL
		try(data <- eval(object$call$data, parent.frame(nframes)), silent = TRUE)

		if(is.null(data) || is.function(data)){
			dataName = object$call$data
			stop("To apply 'update.fixest', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was '", deparse_long(dataName), "').")
		} else if(!is.data.frame(data)){
		    stop("To apply 'update.fixest', we fetch the original database in the parent.frame -- but currently the object '", deparse_long(object$call$data), "' is not a data.frame.")
		}

		# if same data => we apply init
		n_old = object$nobs + length(object$obsRemoved)
		n_new = nrow(data)
		if(n_old != n_new){
			useInit = FALSE
		}
	}

	#
	# III) cluster updates
	#

	fixef_vars = object$fixef_vars
	sumFE_init = NULL
	from_update = FALSE
	if(length(FML)[2] > 1){
		# modification of the clusters
	    if(!is.null(object$fixef_terms)){
	        fixef_old = as.formula(paste0("~", paste0(c(1, object$fixef_terms), collapse = "+")))
	    } else {
	        fixef_old = as.formula(paste0("~", paste0(c(1, fixef_vars), collapse = "+")))
	    }

	    # I use it as text to catch the var1^var2 FEs (update does not work)
	    fixef_old_text = deparse_long(fixef_old)
	    fixef_new_text = deparse_long(formula(FML, lhs = 0, rhs = 2))

	    if(fixef_new_text == "~."){
	        # nothing happens
	        fixef_new = fixef_old
	    } else if(fixef_new_text == "~1"){
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
	        fixef_new = update(fixef_old, formula(FML, lhs = 0, rhs = 2))
	    }

		if(useInit){
			# Only if we use the init => the starting cluster values
			isThere = sapply(fixef_vars, function(x) grepl(x, as.character(fixef_new)[2]))
			if(is.null(object$obsRemoved)){
				# ONLY when there is no cluster removed (otherwise computationaly too complex to be worth)
				if(all(isThere)){
					# we just have to put the old
					sumFE_init = object$sumFE
				} else if(any(isThere)){
					# we use the dummies only for the ones that are there
					my_fe = fixef(object)
					sumFE_init = 0
					for(i in which(isThere)){
						sumFE_init = sumFE_init + my_fe[[i]][object$fixef_id[[i]]]
					}
				}
			}
		}

		if(length(all.vars(fixef_new)) > 0){
			# means there is a cluster
			fml_new = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", as.character(fixef_new)[2]))
		} else {
			# there is no cluster
			fml_new = fml
		}

	} else if(!is.null(fixef_vars)){
		# Means we keep the same clusters
		from_update = TRUE

		# the starting value:
		sumFE_init = object$sumFE

		# the formula updated:
		fml_new = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", paste0(fixef_vars, collapse = "+")))

	} else {
		# there is no cluster in the initial model
		fml_new = fml
	}


	#
	# The call
	#

	call_old = object$call

	# we drop the argument fixef from old call (now it's in the fml_new)
	call_old$fixef = NULL

	# new call: call_clear
	call_clear = call_old
	for(arg in setdiff(names(call_new)[-1], c("fml.update", "nframes"))){
		call_clear[[arg]] = call_new[[arg]]
	}

	call_clear$fml = as.call(fml_new)

	if(useInit){
		# we can use the initialisation of parameters
		if(object$method %in% c("femlm", "feNmlm", "fenegbin")){
			if(object$family == "negbin"){
				if(is.null(dots$theta.init)){
					theta.init = object$theta.init
					call_clear$theta.init = theta.init
				}
			}
		}

		call_clear$from_update = from_update
		call_clear$sumFE_init = as.name("sumFE_init")
	}

	# The variable "sumFE_init" must be evaluated here!
	res = eval(call_clear, list(sumFE_init=sumFE_init), parent.frame(nframes))

	res
}


#' Extract the formula of a \code{fixest} fit
#'
#' This function extracts the formula from a \code{fixest} estimation (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). If the estimation was done with fixed-effects, they are added in the formula after a pipe (\dQuote{|}). If the estimation was done with a non linear in parameters part, then this will be added in the formula in between \code{I()}.
#'
#' @inheritParams nobs.fixest
#'
#' @param x An object of class \code{fixest}. Typically the result of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation.
#' @param type A character scalar. Default is \code{type = "full"} which gives back a formula containing the linear part of the model along with the clusters (if any) and the non-linear in parameters part (if any). If \code{type = "linear"} then only the linear formula is returned. If \code{type = "NL"} then only the non linear in parameters part is returned.
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
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # formula with the cluster variable
#' formula(res)
#' # linear part without the cluster variable
#' formula(res, "linear")
#'
#'
formula.fixest = function(x, type = c("full", "linear", "NL"), ...){
	# Extract the formula from the object
	# we add the clusters in the formula if needed

    # Checking the arguments
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = "type")
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

    if(isTRUE(x$fromFit)){
        stop("formula method not available for fixest estimations obtained from fit methods.")
    }

	type = match.arg(type)

	if(type == "linear"){
		return(x$fml)
	} else if(type == "NL"){

		if(x$method != "femlm"){
			stop("type = 'NL' is not valid for a ", x$method, " estimation.")
		}

		NL.fml = x$NL.fml
		if(is.null(NL.fml)){
			stop("There was no nonlinear part estimated, option type = 'NL' cannot be used.")
		}

		return(NL.fml)
	}

	if(is.null(x$NL.fml)){
	    if(is.null(x$fml_full)){
	        res = x$fml
	    } else {
	        res = x$fml_full
	    }
	} else if(is.null(x$fml_full)){
	    fml_char = deparse_long(x$fml)
	    nl_char = as.character(x$NL.fml)
	    res = as.formula(paste(fml_char, "+ I(", nl_char[2], ")"))
	} else {
	    fml_split = strsplit(deparse_long(x$fml_full), "\\|")[[1]]
	    nl_char = as.character(x$NL.fml)
	    res = as.formula(paste(fml_split[1], "+ I(", nl_char[2], ") |", fml_split[2]))
	}

	res
}


#' Design matrix of a \code{femlm} model
#'
#' This function creates a design matrix of the linear part of a \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimation. Note that it is only the linear part and the cluster variables (which can be considered as factors) are excluded from the matrix.
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
#' # simple estimation on iris data, clustering by "Species"
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
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = "data")
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

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
	if(original_data && !is.null(object$obsRemoved)){
	    check_0 = TRUE
	    linear.mat = linear.mat[-object$obsRemoved, , drop = FALSE]
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
	    only_0 = cpppar_check_only_0(linear.mat, nrow(linear.mat), nthreads = 1)
	    if(all(only_0 == 1)){
	        stop("After removing NAs, not a single explanatory variable is different from 0.")
	    } else if(any(only_0 == 1)){
	        linear.mat = linear.mat[, only_0 == 0, drop = FALSE]
	    }
	}

    linear.mat
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
#' \donttest{
#' # Default is TRUE
#' getFixest_notes()
#' # Change default with
#' setFixest_notes(FALSE)
#' }
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
#' @param nthreads An integer strictly greater than one and lower than the maximum number of threads (if OpenMP is available). If missing, the default is the maximum number of threads minus two.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' \donttest{
#' # Gets the current number of threads
#' getFixest_nthreads()
#' # To set multi-threading off:
#' setFixest_nthreads(1)
#' # To set it back to default:
#' setFixest_nthreads()
#' }
#'
#'
setFixest_nthreads = function(nthreads){
	# By default, we leave 2 nthreads (never use all)

	max_threads = min(get_nb_threads(), 1000) # we cap at 1k nthreads

	if(missing(nthreads) || is.null(nthreads)){
		nthreads = max(1, max_threads - 2)
	}

	if(length(nthreads) != 1 || !is.numeric(nthreads) || is.na(nthreads) || nthreads %% 1 != 0 || nthreads < 0){
		stop("Argument 'nthreads' must be equal to a positive integer.")
	}

	if(nthreads > 1){
		if(max_threads == 0){
			warning("OpenMP not detected: cannot use ", nthreads, " threads, single-threaded mode instead.")
			nthreads = 1
		} else if(nthreads > max_threads){
			warning("Asked for ", nthreads, " threads while the maximum is ", max_threads, ". Set to ", max_threads, " threads instead.")
			nthreads = max_threads
		}
	}

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

#' Sets/gets the dictionary used in \code{esttex}
#'
#' Sets/gets the default dictionary used in the function \code{\link[fixest]{esttex}}. The dictionaries are used to relabel variables (usually towards a fancier, more explicit formatting) when exporting them into a Latex table. By setting the dictionary with \code{setFixest_dict}, you can avoid providing the argument \code{dict} in function \code{\link[fixest]{esttex}}.
#'
#'
#' @param dict A named character vector. E.g. to change my variable named "a" and "b" to (resp.) "$log(a)$" and "$bonus^3$", then use \code{dict = c(a="$log(a)$", b3="$bonus^3$")}.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' \donttest{
#' data(trade)
#' est = feols(log(Euros) ~ log(dist_km)|Origin+Destination+Product, trade)
#' # we export the result & rename some variables
#' esttex(est, dict = c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#' # If you export many tables, it can be more convenient to use setFixest_dict:
#' setFixest_dict(c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#' esttex(est) # variables are properly relabeled
#' }
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

#' Sets/gets whether to remove NA/Inf values from \code{fixest} estimations
#'
#' Sets/gets the default policy of NA/Inf behavior in \code{fixest} estimations. By default, NA/Inf values are removed (and a note is displayed). If you prefer a no NA policy, just set \code{setFixest_na_inf.rm(FALSE)}.
#'
#' @param x A Logical.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' \donttest{
#' base = iris
#' base[1, 1] = NA
#' # default: NAs removed
#' res = feols(Sepal.Length ~ Sepal.Width, base)
#' # no tolerance: estimation fails
#' res = feols(Sepal.Length ~ Sepal.Width, base, na_inf.rm = FALSE)
#' # to set no tolerance as default:
#' setFixest_na_inf.rm(FALSE)
#' }
#'
setFixest_na_inf.rm = function(x){

    if(length(x) != 1 || !is.logical(x) || is.na(x)){
        stop("Argument 'x' must be either TRUE or FALSE.")
    }

    options("fixest_na_inf.rm" = x)
}

#' @rdname setFixest_na_inf.rm
"getFixest_na_inf.rm"

getFixest_na_inf.rm = function(){

    x = getOption("fixest_na_inf.rm")
    if(length(x) != 1 || !is.logical(x) || is.na(x)){
        stop("The value of getOption(\"fixest_na_inf.rm\") is currently not legal. Please use function setFixest_na_inf.rm to set it to an appropriate value. ")
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
#' \donttest{
#' res = feols(Sepal.Length ~ Sepal.Width + Petal.Length, base)
#' # default is coef. table:
#' res
#' # can be changed to only the coefficients:
#' print(res, type = "coef")
#' setFixest_print.type("coef")
#' res # only the coefs
#' }
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


#' Sets the defaults of coefplot
#'
#' You can set the default values of most arguments of \code{\link[fixest]{coefplot}} with this function.
#'
#' @inheritParams coefplot
#'
#' @param reset Logical, default is \code{TRUE}. If \code{TRUE}, then the arguments that *are not* set during the call are reset to their "factory"-default values. If \code{FALSE}, on the other hand, arguments that have already been modified are not changed.
#'
#' @return
#' Doesn't return anything.
#'
#' @seealso
#' \code{\link[fixest]{coefplot}}
#'
#' @examples
#'
#' # coefplot has many arguments, which makes it highly flexible.
#' # If you don't like the default style of coefplot. No worries,
#' # you can set *your* default by using the function
#' # setFixest_coefplot()
#'
#' # Estimation
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Length +
#'                 Sepal.Width | Species, iris)
#'
#' # Plot with default style
#' coefplot(est)
#'
#' # Now we permanently change some arguments
#' dict = c("Petal.Length"="Length (Petal)", "Petal.Width"="Width (Petal)",
#'          "Sepal.Length"="Length (Sepal)", "Sepal.Width"="Width (Sepal)")
#'
#' setFixest_coefplot(ci.col = 2, pt.col = "darkblue", ci.lwd = 3,
#'                    pt.cex = 2, pt.pch = 15, ci.width = 0, dict = dict)
#'
#' # Tadaaa!
#' coefplot(est)
#'
#' # To reset to the default settings:
#' setFixest_coefplot()
#' coefplot(est)
#'
setFixest_coefplot = function(dict, ci.width=0.1, ci_level = 0.95, pt.pch = 20, cex = par("cex"), pt.cex = cex, col = 1, pt.col = col, ci.col = col, lwd = par("lwd"), ci.lwd = lwd, ci.lty, grid = TRUE, grid.par = list(lty=3, col = "gray"), zero = TRUE, zero.par = list(col="black", lwd=1), join = FALSE, join.par = list(lwd=lwd), ref.line, ref.line.par = list(col = "black", lty = 2), only.inter = TRUE, reset = TRUE){

    #
    # Controls
    #

    check_arg(ci.width, "singleNumericGE0")
    check_arg(ci_level, "singleNumericGE0LT1")
    check_arg(lwd, "singleNumericGE0")
    check_arg(ci.lwd, "singleNumericGE0")
    check_arg(grid, "singleLogical")
    if(is.null(grid.par)){
        grid.par = list()
    } else if(!is.list(grid.par) ) {
        stop("Argument grid.par must be a list (even empty).")
    }
    check_arg(zero, "singleLogical")
    if(is.null(zero.par)){
        zero.par = list()
    } else if(!is.list(zero.par) ) {
        stop("Argument zero.par must be a list (even empty).")
    }
    check_arg(join, "singleLogical")
    if(is.null(join.par)){
        join.par = list()
    } else if(!is.list(join.par) ) {
        stop("Argument join.par must be a list (even empty).")
    }
    check_arg(ref.line, "singleLogical")
    if(is.null(ref.line.par)){
        ref.line.par = list()
    } else if(!is.list(ref.line.par) ) {
        stop("Argument ref.line.par must be a list (even empty).")
    }
    check_arg(only.inter, "singleLogical")
    check_arg(reset, "singleLogical")

    #
    # Code
    #

    if(reset){
        opts = list()
    } else {
        opts = getOption("fixest_coefplot")
        if(is.null(opts)){
            opts = list()
        } else if(!is.list(opts)){
            warning("Wrong format of getOption('fixest_coefplot'), the options of coefplot are reset.")
            opts = list()
        }
    }

    mc = match.call()

    all_args = setdiff(names(mc), c("", "reset"))

    for(arg in all_args){
        opts[[arg]] = eval(mc[[arg]])
    }

    options("fixest_coefplot" = opts)
}


#' Type of degree of freedom in fixest summary
#'
#' Provides how the degrees of freedom should be calculated in \code{\link[fixest]{vcov.fixest}}/\code{\link[fixest]{summary.fixest}}.
#'
#' @param fixef How to account for the fixed-effects parameters, defaults to \code{"nested"}. If \code{FALSE} or \code{"no"}, fixed-effects parameters are discarded, meaning the number of parameters is only equal to the number of variables. If \code{TRUE} or \code{yes}, then the number of parameters is equal to the number of variables plus the number of fixed-effects. Finally, if \code{nested}, then the number of parameters is equal to the number of variables an the number of fixed-effects that *are not* nested in the clusters used to cluster the standard-errors.
#' @param exact Logical, default is \code{FALSE}. If there are 2 or more fixed-effects, these fixed-effects they can be irregular, meaning they can provide the same information. If so, the "real" number of parameters should be lower than the total number of fixed-effects. If \code{exact = TRUE}, then \code{\link[fixest]{fixef.fixest}} is first run to determine the exact number of parameters among the fixed-effects. Mostly, panels of the type individual-firm-year require \code{exact = TRUE} (but it adds computational costs).
#' @param cluster Logical, default is \code{TRUE}. How to make the small sample correction when clustering the standard-errors? If \code{TRUE} a \code{G/(G-1)} correction is performed with \code{G} the number of cluster values.
#'
#' @return
#' It returns a \code{dof.type} object.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{setFixest_dof}}
#'
#' @examples
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
#' # Default: fixef = "nested"
#' #  => adjustment K = 1 + 5 (i.e. x + fe2)
#' summary(est)
#' attributes(vcov(est))[c("dof.type", "dof.K")]
#'
#'
#' # fixef = FALSE
#' #  => adjustment K = 1 (i.e. only x)
#' summary(est, dof = dof(fixef=FALSE))
#' attr(vcov(est, dof = dof(fixef=FALSE)), "dof.K")
#'
#'
#' # fixef = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 1 (i.e. x + fe1 + fe2 - 1 restriction)
#' summary(est, dof = dof(fixef=TRUE))
#' attr(vcov(est, dof = dof(fixef=TRUE)), "dof.K")
#'
#'
#' # fixef = TRUE & exact = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 2 (i.e. x + fe1 + fe2 - 2 restrictions)
#' summary(est, dof = dof(fixef=TRUE, exact = TRUE))
#' attr(vcov(est, dof = dof(fixef=TRUE, exact = TRUE)), "dof.K")
#'
#' # There are two restrictions:
#' attr(fixef(est), "references")
#'
#'
dof = function(fixef = "nested", exact = FALSE, cluster = TRUE){

    if(isLogical(fixef)){
        fixef = as.character(fixef)
    }

    check_arg(fixef, "singleCharacter", "Argument 'fixef' must be equal to TRUE, 'yes', FALSE, 'no', or 'nested' (default). REASON")

    check_arg(exact, "singleLogical")
    check_arg(cluster, "singleLogical")

    fixef = try(match.arg(tolower(fixef), c("no", "false", "yes", "true", "nested")), silent = TRUE)
    if("try-error" %in% class(fixef)){
        stop("Argument 'fixef' must be equal to TRUE, 'yes', FALSE, 'no', or 'nested' (default)")
    }

    if(fixef == "TRUE") fixef = "yes"
    if(fixef == "FALSE") fixef = "no"

    res = list(fixef = fixef, exact = exact, cluster = cluster)
    class(res) = "dof.type"
    res
}


#' Sets the default type of DoF correction in summary/vcov.fixest
#'
#' Sets or gets the  default type of DOF correction in \code{\link[fixest]{summary.fixest}} and \code{\link[fixest]{vcov.fixest}}
#'
#' @param dof.type An object of class \code{dof.type} obtained with the function \code{\link[fixest]{dof}}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{dof}}
#'
#' @return
#' The function \code{getFixest_dof} returns a \code{dof.type} object.
#'
#' @examples
#'
#' \dontrun{
#'
#' # If you never want DoF correction when computing the vcov
#' # of fixest object:
#'
#' setFixest_dof(dof(fixef = FALSE, cluster = FALSE))
#'
#' }
#'
#'
setFixest_dof = function(dof.type = dof()){

    if(!"dof.type" %in% class(dof.type)){
        stop("The argument 'dof.type' must be an object created by the function dof().")
    }

    options("fixest_dof" = dof.type)
}

#' @rdname setFixest_dof
"getFixest_dof"

getFixest_dof = function(){

    dof.type = getOption("fixest_dof")
    if(!"dof.type" %in% class(dof.type)){
        stop("The value of getOption(\"fixest_dof\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value.")
    }

    dof.type
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































































