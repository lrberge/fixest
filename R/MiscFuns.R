
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
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
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
		warning("The optimization algorithm did not converge, the results are not reliable. (", x$message, ")")
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
#' @param dof Logical, default is \code{TRUE}. Should there be a degree of freedom correction to the standard errors of the coefficients?
#' @param exact_dof Logical, default is \code{FALSE}. In case there were 2+ clusters in the estimation, it computes the exact number of degrees of freedom (this is not needed in case of balanced panels).
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
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
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
summary.fixest <- function(object, se, cluster, dof = TRUE, exact_dof = FALSE, forceCovariance = FALSE, keepBounded = FALSE, ...){
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
	pvalue <- 2*pnorm(-abs(zvalue))

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

#' Facility to export the results of multiple \code{fixest} estimations in a Latex table.
#'
#' This function aggregates the results of multiple estimations and display them in the form of  one Latex table whose row names are the variables and the columns contain the coefficients and standard-errors.
#'
#' @inheritParams summary.fixest
#'
#' @param ... Used to capture different \code{fixest} objects (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). Note that any other type of element is discarded. Note that you can give a list of \code{fixest} objects.
#' @param digits Integer, default is 4. The number of digits to be displayed.
#' @param pseudo Logical, default is \code{TRUE}. Should the pseudo R2 be displayed?
#' @param title Character scalar. The title of the Latex table.
#' @param sdBelow Logical, default is \code{TRUE}. Should the standard-errors be displayed below the coefficients?
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a regular expression (see \code{\link[base]{regex}} help for more info). There can be more than one regular expression. Each variable satisfying the regular expression will be discarded.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a regular expression (see \code{\link[base]{regex}} help for more info). There can be more than one regular expression. The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions.
#' @param dict A named character vector. It changes the original variable names to the ones contained in the \code{dict}. E.g. to change the variables named \code{a} and \code{b3} to (resp.) \dQuote{$log(a)$} and to \dQuote{$bonus^3$}, use \code{dict=c(a="$log(a)$",b3="$bonus^3$")}. By default it is equal to \code{getFixest_dict()}, a default dictionary can be set with \code{\link[fixest]{setFixest_dict}}.
#' @param file A character scalar. If provided, the Latex table will be saved in a file whose path is \code{file}.
#' @param replace Logical, default is \code{FALSE}. Only used if option \code{file} is used. Should the Latex table be written in a new file that replaces any existing file?
#' @param convergence Logical, default is missing. Should the convergence state of the algorithm be displayed? By default, convergence information is displayed if at least one model did not converge.
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)}.
#' @param label Character scalar. The label of the Latex table.
#' @param aic Logical, default is \code{FALSE}. Should the AIC be displayed?
#' @param sqCor Logical, default is \code{FALSE}. Should the squared correlation be displayed?
#' @param subtitles Character vector of the same length as the number of models to be displayed. If provided, subtitles are added underneath the dependent variable name.
#' @param fixef_sizes Logical, default is \code{FALSE}. If \code{TRUE} and fixed-effects were used in the models, then the number "individuals" per fixed-effect dimension is also displayed.
#' @param keepFactors Logical, default is \code{TRUE}. If \code{FALSE}, then factor variables are displayed as fixed-effects and no coefficient is shown.
#' @param bic Logical, default is \code{TRUE}.Should the BIC be reported?
#' @param loglik Logical, default is \code{TRUE}. Should the log-likelihood be reported?
#' @param yesNoFixef A character vector of length 2. Default is \code{c("Yes", "No")}. This is the message displayed when a given cluster is (or is not) included in a regression.
#' @param family A logical, default is missing. Whether to display the families of the models. By default this line is displayed when at least two models are from different families.
#' @param powerBelow Integer, default is -5. A coefficient whose value is below \code{10**(powerBelow+1)} is written with a power in Latex. For example \code{0.0000456} would be written \code{4.56$\\times 10^{-5}$} by default. Setting \code{powerBelow = -6} would lead to \code{0.00004} in Latex.
#'
#' @return
#' There is nothing returned, the result is only displayed on the console or saved in a file.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#'# two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'              Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' # We export the three results in one Latex table,
#' # with clustered standard-errors:
#' esttex(res1, res2, se = "cluster")
#'
#' # Changing the names & significance codes
#' esttex(res1, res2, dict = c(Sepal.Length = "The sepal length", Sepal.Width = "SW"),
#'         signifCode = c("**" = 0.1, "*" = 0.2, "n.s."=1))
#'
esttex <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = TRUE, cluster, digits=4, pseudo=TRUE, title, sdBelow=TRUE, drop, order, dict = getFixest_dict(), file, replace=FALSE, convergence, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, aic=FALSE, sqCor = FALSE, subtitles, fixef_sizes = FALSE, bic = TRUE, loglik = FALSE, yesNoFixef = c("Yes", "No"), keepFactors = TRUE, family, powerBelow = -5){
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

	# to get the model names
	dots_call = match.call(expand.dots = FALSE)[["..."]]

	info = results2formattedList(..., se=se, dof=dof, cluster=cluster, digits=digits, sdBelow=sdBelow, signifCode=signifCode, subtitles=subtitles, yesNoFixef=yesNoFixef, keepFactors=keepFactors, isTex = TRUE, useSummary=useSummary, dots_call=dots_call, powerBelow=powerBelow, dict=dict)

	n_models = length(info$depvar_list)
	# Getting the information
	se_type_list = info$se_type_list
	var_list = info$var_list
	coef_list = info$coef_list
	coef_below = info$coef_below
	sd_below = info$sd_below
	depvar_list = info$depvar_list
	obs_list = info$obs_list
	r2_list = info$r2_list
	aic_list = info$aic_list
	bic_list = info$bic_list
	loglik_list = info$loglik_list
	convergence_list = info$convergence_list
	sqCor_list = info$sqCor_list
	factorNames = info$factorNames
	isFactor = info$isFactor
	nbFactor = info$nbFactor
	slope_names = info$slope_names
	slope_flag_list = info$slope_flag_list
	family_list = info$family_list
	theta_list = info$theta_list

	if(!missing(subtitles)){
		isSubtitles = TRUE
	} else {
		isSubtitles = FALSE
	}

	#
	# prompting the infos gathered
	#

	# Starting the table
	myTitle = ifelse(!missing(title), title, "no title")
	if(!missing(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
	start_table = paste0("\\begin{table}[htbp]\\centering\n\\caption{",  .cleanPCT(myTitle), "}\n")
	end_table = "\\end{table}"

	# intro and outro Latex tabular
	myAmpLine = paste0(paste0(rep(" ", length(depvar_list)+1), collapse="&"), "\\tabularnewline\n")
	intro_latex <- paste0("\\begin{tabular}{l", paste0(rep("c", n_models), collapse=""), "}\n",
								 myAmpLine,
								 "\\hline\n",
								 "\\hline\n")

	outro_latex <- "\\end{tabular}\n"

	# 1st lines => dep vars
	# first_line <- paste0("Variables&", paste0(depvar_list, collapse="&"), "\\\\\n\\hline\n\\hline\n")
	depvar_list = c(depvar_list, recursive = TRUE)
	if(!missing(dict)){
		if(!is.character(dict)|| is.null(names(dict))) stop("the arg. 'dict' must be a named character vector.")
		depvar_list = c(depvar_list, recursive = TRUE)
		qui = which(depvar_list%in%names(dict))
		who = depvar_list[qui]
		depvar_list[qui] = dict[who]
	}

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
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regex).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop, all_vars)]
	}

	# ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regex).")
		for(var2order in rev(order)){
			who = grepl(var2order, all_vars)
			all_vars = c(all_vars[who], all_vars[!who])
		}
	}

	# changing the names of the coefs
	aliasVars = all_vars
	names(aliasVars) = all_vars

	if(!missing(dict)){
		if(!is.character(dict)|| is.null(names(dict))){
			stop("the arg. 'dict' must be a named character vector.")
		}
	}

	qui = all_vars %in% names(dict)
	who = aliasVars[qui]
	aliasVars[qui] = .cleanPCT(dict[who])

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

	# Factors (if needed)
	if(length(factorNames) > 0){
		dumIntro = paste0("\\hline\n\\emph{Fixed-Effects}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")

		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = yesNoFixef[2]
			isFactor[[m]] = quoi

			# We do the same for the number of items
			quoi = nbFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "--"
			nbFactor[[m]] = quoi
		}

		allFactors = matrix(c(isFactor, recursive=TRUE), nrow = length(factorNames))
		# We change the names of the factors
		qui = which(factorNames %in% names(dict))
		if(length(qui) > 0) {
			factorNames[qui] = dict[factorNames[qui]]
		}

		allFactors = cbind(factorNames, allFactors)
		factor_lines <- paste0(paste0(apply(allFactors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

		# For the number of items
		all_nb_Factors = matrix(c(nbFactor, recursive=TRUE), nrow = length(factorNames))
		factorNames_nbItems = paste0("# ", factorNames)
		all_nb_Factors = cbind(factorNames_nbItems, all_nb_Factors)
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

	# Fit statistics
	fit_part <- paste0("\\hline\n\\emph{Fit statistics}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
	# Misc
	info_aic <- paste0("AIC & ", paste(numberFormatLatex(aic_list), collapse="&"), "\\\\\n")
	info_loglik <- paste0("Log-Likelihood & ", paste(numberFormatLatex(loglik_list), collapse="&"), "\\\\\n")
	info_bic <- paste0("BIC & ", paste(numberFormatLatex(bic_list), collapse="&"), "\\\\\n")

	info_obs <- paste0("Observations& ", paste(addCommas(obs_list), collapse="&"), "\\\\\n")
	info_r2 <- paste0("Adj-pseudo $R^2$ &", paste(r2_list, collapse="&"), "\\\\\n")
	info_sqCor <- paste0("$R^2$ &", paste(sqCor_list, collapse="&"), "\\\\\n")

	# Convergence information
	info_convergence = ""
	if((missing(convergence) && any(convergence_list == FALSE)) || (!missing(convergence) && convergence)){
		info_convergence = paste0("Convergence &", paste(convergence_list, collapse="&"), "\\\\\n")
	}

	info_theta <- paste0("Overdispersion& ", paste(theta_list, collapse="&"), "\\\\\n")

	# information on family
	if((!missing(family) && family) || (missing(family) && length(unique(family_list)) > 1)){
		info_family <- paste0("Family& ", paste(family_list, collapse="&"), "\\\\\n")
	} else {
		info_family = ""
	}


	# The standard errors
	isUniqueSD = length(unique(unlist(se_type_list))) == 1
	if(isUniqueSD){
		my_se = unique(unlist(se_type_list)) # it comes from summary
		# every model has the same type of SE
		if(my_se == "Standard") my_se = "Normal"
		if(my_se == "White") my_se = "White-corrected"

		# Now we modify the names of the clusters if needed
		if(grepl("\\(", my_se)){
			# we extract the clusters
			se_cluster = strsplit(gsub("(^.+\\()|(\\))", "", my_se), " & ")[[1]]
			qui = se_cluster %in% names(dict)
			se_cluster[qui] = dict[se_cluster[qui]]
			new_se = gsub("\\(.+", "", my_se)
			my_se = paste0(new_se, "(", paste0(se_cluster, collapse = " & "), ")")
		}

		nb_col = length(obs_list) + 1
		info_SD = paste0("\\hline\n\\hline\n\\multicolumn{", nb_col, "}{l}{\\emph{", my_se, " standard-errors in parenthesis. Signif Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
		info_muli_se = ""
	} else {
		info_muli_se = paste0("Standard-Error type& ", paste(se_type_list, collapse="&"), "\\\\\n")
		info_SD = "\\hline\n\\hline\n\\\\\n"
	}

	# Information on number of items

	supplemental_info = "\\global\\long\\def\\sym#1{\\ifmmode^{#1}\\else\\ensuremath{^{#1}}\\fi}\n"

	if(!pseudo) info_r2 <- ""
	if(!sqCor) info_sqCor <- ""
	if(!aic) info_aic = ""
	if(!bic) info_bic = ""
	if(!loglik) info_loglik = ""
	if(!fixef_sizes) nb_factor_lines = ""
	if(all(theta_list == "")) info_theta = ""

	if(!missing(file)) sink(file = file, append = !replace)

	cat(paste0(supplemental_info,
				  start_table,
				  intro_latex,
				  first_line,
				  info_subtitles,
				  model_line,
				  info_family,
				  variable_line,
				  coef_lines,
				  info_theta,
				  dumIntro,
				  factor_lines,
				  slope_intro,
				  slope_lines,
				  fit_part,
				  info_obs,
				  nb_factor_lines,
				  info_convergence,
				  info_muli_se,
				  info_r2,
				  info_sqCor,
				  info_aic,
				  info_loglik,
				  info_bic,
				  info_SD,
				  outro_latex,
				  end_table))

	if(!missing(file)) sink()

}

#' Facility to display the results of multiple \code{fixest} estimations.
#'
#' This function aggregates the results of multiple estimations and display them in the form of only one table whose row names are the variables and the columns contain the coefficients and standard-errors.
#'
#' @inheritParams esttex
#' @inheritParams summary.fixest
#'
#' @param depvar Logical, default is missing. Whether a first line containing the dependent variables should be shown. By default, the dependent variables are shown only if they differ across models.
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)}.
#' @param titles A character vector. The length must match the number of models.
#'
#' @return
#' Returns a data.frame containing the formatted results.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' # estimation without clusters
#' res2 = update(res1, . ~ Sepal.Width | 0)
#'
#' # We export the two results in one Latex table:
#' esttable(res1, res2)
#'
#' # With clustered standard-errors + showing the dependent variable
#' esttable(res1, res2, se = "cluster", cluster = iris$Species, depvar = TRUE)
#'
#' # Changing the model names + the order of the variables
#' # + dropping the intercept.
#' esttable(model_1 = res1, res2,
#'           order = c("Width", "Petal"), drop = "Int",
#'           signifCode = c("**" = 0, "*" = 0.2, "n.s."=1))
#'
#'
#'
esttable <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = TRUE, cluster, depvar, drop, order, digits=4, pseudo=TRUE, convergence, signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10), titles, keepFactors = FALSE, family, bic = TRUE, loglik = FALSE){

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

	info = results2formattedList(..., se=se, dof = dof, cluster=cluster, digits=digits, signifCode=signifCode, titles=titles, keepFactors=keepFactors, useSummary=useSummary, dots_call=dots_call)

	n_models = length(info$depvar_list)
	# Getting the information
	se_type_list = info$se_type_list
	var_list = info$var_list
	coef_list = info$coef_list
	coef_below = info$coef_below
	sd_below = info$sd_below
	depvar_list = info$depvar_list
	obs_list = info$obs_list
	r2_list = info$r2_list
	aic_list = info$aic_list
	bic_list = info$bic_list
	loglik_list = info$loglik_list
	convergence_list = info$convergence_list
	sqCor_list = info$sqCor_list
	factorNames = info$factorNames
	isFactor = info$isFactor
	nbFactor = info$nbFactor
	slope_names = info$slope_names
	slope_flag_list = info$slope_flag_list
	useSummary = info$useSummary
	depvar_list = info$depvar_list
	model_names = info$model_names
	family_list = info$family_list
	theta_list = info$theta_list


	if(!missing(titles)){
		isTitles = TRUE
	} else {
	    isTitles = FALSE
	}

	# The coefficients

	all_vars <- unique(c(var_list, recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regex).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop, all_vars)]
	}

	# ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regex).")
		for(var2order in rev(order)){
			who = grepl(var2order, all_vars)
			all_vars = c(all_vars[who], all_vars[!who])
		}
	}

	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	res = coef_mat

	if("Neg. Bin." %in% family_list){
		theta_line = c("Overdispersion:", unlist(theta_list))
		res = rbind(res, theta_line)
	}

	# Used to draw a line
	myLine = "______________________________________"
	longueur = apply(res, 2, function(x) max(nchar(as.character(x))))
	theLine = sapply(longueur, function(x) sprintf("%.*s", x, myLine))
	theLine[1] = sprintf("%.*s", max(nchar(theLine[1]), 19), myLine)

	# The FEs
	if(length(factorNames)>0){

		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "No"
			isFactor[[m]] = quoi
		}
		allFactors = matrix(c(isFactor, recursive=TRUE), nrow = length(factorNames))
		allFactors = cbind(factorNames, allFactors)
		factor_lines <- paste0(paste0(apply(allFactors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

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

	# The line with the dependent variable
	preamble = c()
	if((missing(depvar) && length(unique(unlist(depvar_list))) > 1) || (!missing(depvar) && depvar)){
		preamble = rbind(c("Dependent Var.:", unlist(depvar_list)), preamble)
	}

	if(length(preamble) > 0){
		# preamble = rbind(preamble, c("  ", theLine[-1]))
		preamble = rbind(preamble, rep("   ", length(theLine)))
		res <- rbind(preamble, res)
	}

	res <- rbind(res, theLine)

	# the line with the families
	if((missing(family) && length(unique(unlist(family_list))) > 1) || (!missing(family) && family)){
		# preamble = rbind(c("Family:", unlist(family_list)), preamble)
		res = rbind(res, c("Family:", unlist(family_list)))
	}

	res <- rbind(res, c("Observations", addCommas(obs_list)))
	if(!useSummary){
		se_type_format = c()
		for(m in 1:n_models) se_type_format[m] = charShorten(se_type_list[[m]], longueur[[1+m]])
		res <- rbind(res, c("S.E. type", c(se_type_format, recursive = TRUE)))
	}

	# convergence status
	if((missing(convergence) && any(convergence_list == FALSE)) || (!missing(convergence) && convergence)){
		res <- rbind(res, c("Convergence", c(convergence_list, recursive = TRUE)))
	}

	res <- rbind(res, c("Squared-Corr.", c(sqCor_list, recursive = TRUE)))
	res <- rbind(res, c("Adj-pseudo R2", c(r2_list, recursive = TRUE)))
	# res <- rbind(res, c("AIC", c(aic_list, recursive = TRUE)))
	if(loglik) res <- rbind(res, c("Log-Likelihood", numberFormatNormal(loglik_list)))
	if(bic) res <- rbind(res, c("BIC", numberFormatNormal(bic_list)))

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

	return(res)
}




#' R2s of \code{fixest} models
#'
#' Reports different R2s for \code{fixest} estimations (e.g. \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}).
#'
#' @param x A \code{fixest} object, e.g. obtained with function \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}.
#' @param type A character vector representing the R2 to compute. The R2 codes are of the form: "wapr2" with letters "w" (within), "a" (adjusted) and "p" (pseudo) possibly missing. E.g. to get the regular R2: use \code{type = "r2"}, the within adjusted R2: use \code{type = "war2"}, the pseudo R2: use \code{type = "pr2"}, etc. Use \code{"sq.cor"} for the squared correlation. By default, all R2s are computed.
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
	type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2")
	if("all" %in% type){
		type_all = type_allowed
	} else {
		type_all = tolower(unique(type))
		pblm = setdiff(type_all, type_allowed)
		if(length(pblm) > 0){
			stop("The r2 type", enumerate_items(pblm, addS = TRUE), " not valid.")
		}
	}

	isCluster = "fixef_vars" %in% names(x)
	n = nobs(x)

	res = c()
	for(i in seq_along(type_all)){
		myType = type_all[i]

		if(myType == "sq.cor"){
		    res[i] = x$sq.cor
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

		        res_fe = update(x, .~1|., glm.tol = 1e-2, fixef.tol = 1e-3, nframes = 2)
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

#' Lags a variable using a formula
#'
#' Lags a variable using panel id+time identifiers in a formula.
#'
#'
#' @param x A formula of the type \code{var ~ id + time} where \code{var} is the variable to be lagged, \code{id} is a variable representing the panel id, and \code{time} is the time variable of the panel.
#' @param k An integer giving the number of lags. For leads, just use a negative number.
#' @param data Optional, the data.frame in which to evaluate the formula.
#' @param time.step The method to compute the lags. Can be equal to: \code{"unitary"} (default), \code{"consecutive"} or to a number. If \code{"unitary"}, then the largest common divisor between consecutive time periods is used (typically if the time variable represents years, it will be 1). This method can apply only to integer (or convertible to integer) variables. If \code{"consecutive"}, then the time variable can be of any type: two successive time periods represent a lag of 1. Finally, if the time variable is numeric, you can provide your own numeric time step.
#' @param fill How to fill the observations without defined lead/lag values. Default is \code{NA}.
#' @param duplicate.method If several observations have the same id and time values, then the notion of lag is not defined for them. If \code{duplicate.method = "none"} (default) and duplicate values are found, this leads to an error. You can use \code{duplicate.method = "first"} so that the first occurrence of identical id/time observations will be used as lag.
#' @param ... Not currently used.
#'
#' @return
#' It returns a vector of the same type and length as the variable to be lagged in the formula.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#' # simple example with an unbalanced panel
#' base = data.frame(id = rep(1:2, each = 4),
#'                   time = c(1, 2, 3, 4, 1, 4, 6, 9), x = 1:8)
#'
#' lag(x~id+time,  1, base) # lag 1
#' lag(x~id+time, -1, base) # lead 1
#'
#' lag(x~id+time, 2, base, fill = 0)
#'
#' # with time.step = "consecutive"
#' lag(x~id+time, 1, base, time.step = "cons")
#' # => works for indiv. 2 because 9 (resp. 6) is consecutive to 6 (resp. 4)
#' # mostly useful when the time variable is not a number:
#' # e.g. c("1991q1", "1991q2", "1991q3") etc
#'
#' # with duplicates
#' base_dup = data.frame(id = rep(1:2, each = 4),
#'                       time = c(1, 1, 1, 2, 1, 2, 2, 3), x = 1:8)
#'
#' # by default: error
#' \donttest{
#' lag(x~id+time, 1, base_dup)
#' }
#' # with duplicate.method = "first"
#' lag(x~id+time, 1, base_dup, duplicate.method = "first")
#'
#' \donttest{
#' # You can create variables without specifying the data within data.table:
#' library(data.table)
#' base = data.table(id = rep(1:2, each = 3), year = 1990 + rep(1:3, 2), x = 1:6)
#' base[, x.l1 := lag(x~id+year, 1)]
#' }
#'
#'
lag.formula = function(x, k, data, time.step = "unitary", fill = NA, duplicate.method = c("none", "first"), ...){
    # Arguments:
    # time.step: default: "consecutive", other option: "unitary" (where you find the most common step and use it -- if the data is numeric), other option: a number, of course the time must be numeric

    # LATER:
    # - add duplicate.method = "sum" // although it is not super useful in my opinion (but maybe other users disagree)
    # - several lags => matrix? I donno...

    # Checking arguments in ...
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = c("k", "data", "time.step", "fill"))
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }


    # Controls
    duplicate.method = match.arg(duplicate.method)

    check_arg(k, "integerScalar", "Argument 'k' must be a single integer.")

    if(!missing(data)){
        DATA_MISSING = FALSE
        if(is.matrix(data)){
            if(is.null(colnames(data))){
                stop("The variables of 'data' have no name (data is a matrix without column names).")
            }
            data = as.data.frame(data)
        } else if(!"data.frame" %in%class(data)){
            stop("Argument 'data' must be a data.frame.")
        } else if("data.table" %in% class(data)){
            data = as.data.frame(data)
        }
        existing_vars = names(data)
    } else {
        DATA_MISSING = TRUE
        existing_vars = ls(parent.frame())
    }

    vars = all.vars(x)
    qui_pblm = setdiff(vars, existing_vars)
    if(length(qui_pblm) > 0){
        stop("In the formula the variable", enumerate_items(qui_pblm, addS=TRUE), " not in the ", ifelse(DATA_MISSING, "environment.", "data."))
    }

    if(!isScalar(k, int = TRUE)) stop("Argument 'k' must be a single integer.")

    # LHS
    fml = x
    if(DATA_MISSING){
        value = eval(fml[[2]], parent.frame())
    } else {
        value = eval(fml[[2]], data)
    }

    # checking argument fill
    if(length(fill) != 1){
        stop("The length of argument 'fill' must be exaclty 1. Its current length is ", length(fill), ".")
    } else if(!is.na(fill)){
        if(is.list(fill)){
            stop("Argument fill must be a 'scalar', currenlty it's a list!")
        }

        if(is.numeric(value) && !is.numeric(fill)){
            stop("Argument 'fill' must be of the same type as ", deparse_long(fml[[2]]), ", which is numeric.")
        }

        if(!is.numeric(value) && is.numeric(fill)){
            stop("Argument 'fill' must be of the same type as ", deparse_long(fml[[2]]), ", which is not numeric while 'fill' is.")
        }
        # I could go further in checking but it's enough
    }

    # The id/time
    tm = terms(fml)
    var_id_time = attr(tm, "term.labels")
    if(length(var_id_time) != 2){
        stop("The formula must contain exactly two variables in the right hand side (currently there ", ifsingle(var_id_time, "is", "are"), length(var_id_time), ").")
    }

    if(DATA_MISSING){
        id = eval(parse(text = var_id_time[1]), parent.frame())
        time = eval(parse(text = var_id_time[2]), parent.frame())
    } else {
        id = eval(parse(text = var_id_time[1]), data)
        time = eval(parse(text = var_id_time[2]), data)
    }

    # copy for later (in duplicates)
    id_origin = id
    time_origin = time

    is_na = is.na(id) | is.na(time)
    na_flag = FALSE
    if(any(is_na)){
        na_flag = TRUE
        id = id[!is_na]
        time = time[!is_na]
    }

    # time.step
    if(length(time.step) != 1 || (!is.numeric(time.step) && !is.character(time.step))){
        stop("Argument time.step must be equal to 'unitary', 'consecutive' or to a number.")
    } else if(is.character(time.step)){
        ts = try(match.arg(time.step, c("unitary", "consecutive")), silent = TRUE)
        if("try-error" %in% class(ts)){
            stop("Argument time.step must be equal to 'unitary', 'consecutive' or to a number.")
        }
        time.step = ts
    }

    # unitary time.step: conversion to numeric before quf
    if(time.step == "unitary") {
        time_conversion = FALSE
        if(!is.numeric(time)){
            time_conversion = TRUE

            time_new = tryCatch(as.numeric(time), warning = function(x) x)

            if(!is.numeric(time_new)){
                stop("To use the 'unitary' time.step, the time variable must be numeric or at least convertible to numeric. So far the conversion has failed (time variable's class is currently ", enumerate_items(class(time), endVerb = FALSE), ").")
            }

            time = time_new
        }

        if(any(time %% 1 != 0)){
            stop("To use the 'unitary' time.step, the time variable", ifelse(time_conversion, " (which has been converted to numeric)", ""), " must be made of integers. So far this is not the case. Alternatively, you can give a number in time.step.")
        }

    } else if(time.step != "consecutive"){
        if(!is.numeric(time)){
            stop("If 'time.step' is a number, then the time variable must also be a number (this is not the cases: its class is currently ", enumerate_items(class(time), endVerb = FALSE), ").")
        }

        # if(any(time %% time.step != 0)){
        #     pblm = unique(head(time[time %% time.step != 0], 3))
        #     stop("If 'time.step' is a number, then it must be an exact divisor of all the values in the time variable. This is currently not the case: ", time.step, " is not a divisor of ", enumerate_items(pblm, or = TRUE, endVerb = FALSE), ".")
        # }
    }

    # Computation quf
    id = quickUnclassFactor(id)
    time_full = quickUnclassFactor(time, addItem = TRUE)

    #
    # WIP: define this unitary time step!!!! not straightforward at all!
    # what to do with the ones thatdon't fit the unit???
    # example: time = c(1, 3, 6, 11) => smallest unit is 2, but it does not divide the others

    # Releveling the time ID depending on the time.step
    if(time.step == "consecutive"){
        time = time_full$x
    } else if(time.step == "unitary"){
        time_unik = time_full$items
        all_steps = unique(diff(time_unik))
        my_step = cpp_pgcd(unique(all_steps))

        # we rescale time_unik
        time_unik_new = (time_unik - min(time_unik)) / my_step
        time = time_unik_new[time_full$x]


    } else {
        time_unik = time_full$items

        # consistency check
        all_steps = diff(time_unik)
        if(any(all_steps %% time.step != 0)){
            obs_pblm = which(all_steps %% time.step != 0)

            stop("If 'time.step' is a number, then it must be an exact divisor of all the difference between two consecutive time periods. This is currently not the case: ", time.step, " is not a divisor of ", all_steps[obs_pblm], " (the difference btw ", time_unik[obs_pblm + 1], " and ", , ").")
        }

        # we rescale time_unik // checks done beforehand
        time_unik_new = (time_unik - min(time_unik)) / time.step
        time = time_unik_new[time_full$x]
    }

    order_it = order(id, time)
    order_inv = order(order_it)

    id_sorted = id[order_it]
    time_sorted = time[order_it]

    # We check for duplicate rows => lag not defined for them
    if(duplicate.method == "none"){
        dup_info = cpp_find_duplicates(id_sorted, time_sorted)
        if(dup_info$n_dup > 0){

            if(na_flag){
                obs_ok = which(!is_na)
            } else {
                obs_ok = 1:length(id)
            }

            obs_pblm = obs_ok[order_inv[dup_info$obs]]
            id_dup = id_origin[obs_pblm]
            time_dup = time_origin[obs_pblm]

            stop("The panel identifiers contain duplicate values: this is not allowed since lag/leads are not defined for them. For example (id, time) = (", id_dup, ", ", time_dup, ") appears ", n_times(dup_info$n_dup), ". Please provide data without duplicates -- or you can also use duplicate.method = 'first'.")
        }
    }

    # we get the observation id!
    obs_lagged = cpp_lag_obs(id = id_sorted, time = time_sorted, nlag = k)

    # the lagged value
    value_lagged = value[obs_lagged]
    if(!is.na(fill)){
        qui_na = is.na(obs_lagged)
        value_lagged[qui_na] = fill
    }

    if(na_flag == FALSE){
        res = value_lagged[order_inv]
    } else{
        res = rep(NA, length(value))
        res[!is_na] = value_lagged[order_inv]
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
#' \code{\link[fixest]{plot.fixest.fixef}}. See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
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

	    # evaluation
	    data = NULL
	    try(data <- eval(object$call$data, parent.frame()), silent = TRUE)
	    dataName = object$call$data

	    if(is.null(data)){
	        # We try to provide an informative error message
	        stop("To get the coefficients for the variables with varying slopes, we need to fetch these variables (i.e. ", enumerate_items(slope_vars_unik, endVerb = FALSE), ") in the original dataset in the parent.frame -- but the data doesn't seem to be there anymore (btw it was ", deparse_long(dataName), ")")
	    }

	    data = as.data.frame(data)

	    # we check the variables are there
	    if(any(!slope_vars_unik %in% names(data))){
	        var_pblm = setdiff(slope_vars_unik, names(data))
	        stop("To get the coefficients for the variables with varying slopes, we need to fetch these variables in the original dataset (", deparse_long(dataName), "). However, the variable", enumerate_items(var_pblm, addS = TRUE), " not present in the original dataset any more.")
	    }

	    # we check length consistency
	    if(NROW(data) != (object$nobs + length(object$obsRemoved))){
	        stop("To get the coefficients for the variables with varying slopes, we need to fetch these variables in the original dataset (", deparse_long(dataName), "), yet the dataset doesn't have the same number of observations as was used in the estimation (", NROW(data), " instead of ", object$nobs + length(object$obsRemoved), ").")
	    }

	    if(length(object$obsRemoved)){
	        data = data[-object$obsRemoved, slope_vars_unik, drop = FALSE]
	    } else {
	        data = data[, slope_vars_unik, drop = FALSE]
	    }

	    slope_var_list = list()
	    for(i in 1:length(slope_vars_unik)){
	        variable = all.vars(parse(text = slope_vars_unik[i]))

	        # as.numeric => we'll use cpp so required
	        svar = as.numeric(as.vector(eval(parse(text = slope_vars_unik[i]), data)))
	        if(length(svar) == 1) svar = rep(svar, N)

	        slope_var_list[[slope_vars_unik[i]]] = svar
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

	        fixef_values <- cpp_get_fe_gnl(Q_fe, N, rep(1, N), dumMat, nbCluster, orderCluster)

	        # the information on the references
	        nb_ref_fe = fixef_values[[Q+1]]
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
#' \code{\link[fixest]{fixef.fixest}} to extract clouster coefficients. See also the main estimation function \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
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

		thisNames = getItems(dum_raw)
		dum = quickUnclassFactor(dum_raw)
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
	isLinear = length(all.vars(linear_fml)) > 0
	NL_fml = x$NL.fml
	isNL = !is.null(NL_fml)
	coef = x$coefficients

	# Getting the data
	data = NULL
	try(data <- eval(x$call$data, parent.frame()))

	if(is.null(data)){
		dataName = x$call$data
		stop("To apply function 'collinearity', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", deparse_long(dataName), ").")
	}

	if(!is.null(x$obsRemoved)){
		data = data[-x$obsRemoved, ]
	}

	if(isFE){
		linear_fml = update(linear_fml, ~ . + 1)
	}

	if(isLinear || isCluster || "(Intercept)" %in% names(coef)){
		linear.matrix = model.matrix(linear_fml, data)
	}

	if(isLinear && isCluster){
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
	        stop("To check collinearity, we need to fetch some variables in the original dataset (", deparse_long(dataName), "). However, the variable", enumerate_items(var_pblm, addS = TRUE), " not present in the original dataset any more.")
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
	        # Constant: not OK with FE
	        K = ncol(linear_mat_noIntercept)
	        is_const = logical(K)
	        for(k in 1:K){
	            is_const[k] = cpp_isConstant(linear_mat_noIntercept[, k])
	        }

	        if(any(is_const)){
	            var_problem = colnames(linear_mat_noIntercept)[is_const]
	            message = paste0("Variable", enumerate_items(var_problem, addS = TRUE), " constant, thus collinear with the fixed-effects.")
	            print(message)
	            return(invisible(message))
	        }

	    } else {
	        # constant and equal to 0
	        sum_all = colSums(abs(linear.matrix))
	        if(any(sum_all == 0)){
	            var_problem = colnames(linear.matrix)[sum_all == 0]
	            message = paste0("Variable", enumerate_items(var_problem, addS = TRUE), " constant and equal to 0.")

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

			sum_residuals = colSums(abs(residuals))

			if(any(sum_residuals < 1e-4)){
				ccat("\n")
				varnames = colnames(linear_mat_noIntercept)
				collin_var = varnames[sum_residuals < 1e-4]
				if(length(collin_var) == 1){
					message = paste0("Variable '", collin_var, "' is collinear with fixed-effects ", names(cluster)[q], ".")
				} else {
					message = paste0("Variables ", show_vars_limited_width(collin_var), " are collinear with fixed-effects '", names(cluster)[q], "'.")
				}

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

	        sum_residuals = colSums(abs(residuals))

	        if(any(sum_residuals < 1e-4)){
	            ccat("\n")
	            varnames = colnames(linear_mat_noIntercept)
	            collin_var = varnames[sum_residuals < 1e-4]

	            message = paste0("Variable", enumerate_items(collin_var, addS = TRUE), " collinear with variable with varying slope ", slope_vars[q], " (on ", slope_fe[q], ").")

	            print(message)
	            return(invisible(message))

	        }

	    }
	    ccat("OK")
	}

	#
	# II) perfect multicollinearity
	#

	name2change = grepl("\\)[[:alnum:]]", colnames(linear.matrix))
	if(any(name2change)){
		linbase = as.data.frame(linear.matrix)
		names(linbase)[name2change] = gsub("(\\(|\\))", "_", names(linbase)[name2change])
		linearVars = setdiff(names(linbase), "(Intercept)")
	} else {
		linbase = as.data.frame(linear.matrix)
		linearVars = setdiff(colnames(linear.matrix), "(Intercept)")
	}

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
			sum_resid = sum(abs(res$residuals))

			if(sum_resid < 1e-4){
				ccat("\n")
				coef_lm = res$coefficients
				collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
				message = paste0("Variable '", v, "' is collinear with variable", ifelse(length(collin_var)>=2, "s", ""), ": ", paste0(collin_var, collapse = ", "), ".")

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
			# data[[paste0("__DUM_", i)]] = x$fixef_id[[i]]
			linbase[[paste0("__DUM_", i)]] = x$fixef_id[[i]]
		}

		for(v in linearVars){
			ccat(".")
			fml2estimate = as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))

			for(id_cluster in 1:n_clust){

				res = feols(fml2estimate, linbase, cluster = new_dum_names[id_cluster], warn = FALSE)

				sum_resid = sum(abs(resid(res)))
				if(sum_resid < 1e-4){
					ccat("\n")
					coef_lm = coef(res)
					collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
					message = paste0("Variable '", v, "' is collinear with variable", ifelse(length(collin_var)>=2, "s", ""), ": ", paste0(collin_var, collapse = ", "), " and the fixed-effects ", dum_names[id_cluster], ".")

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


#' Estimates yearly treatment effects
#'
#' This facility helps to estimate yearly treatment effects in a difference-in-difference setup without having to manually make the estimation. It is made as general as possible such that non-\code{fixest} estimation functions can also be used.
#'
#' @param fml A formula containing the variables WITHOUT the yearly treatment effects (which will be added by this function).
#' @param data A \code{data.frame} containing all the variables.
#' @param treat_time Either a character vector of length two containing the name of the treatment variable and the name of the time variable (e.g. \code{c("treat", "year")}). Either a one-sided formula containing the treatment and the time (e.g. \code{~treat~year}).
#' @param reference The time period of reference. It should be a numeric scalar. The treatment will not be included for this time period so that it serves as reference.
#' @param returnData Logical, default is \code{FALSE}. If \code{TRUE}, then the original database with the yearly treatment variables is returned.
#' @param ... Other arguments to be passed to \code{estfun}, the estimation function.
#' @param estfun The estimation function. Defaults to \code{\link[fixest]{feols}}.
#'
#' @return
#' It returns an estimation object. In case of \code{fixest} estimations, it will return a \code{fixest} object.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{did_plot_yearly_effects}}, \code{\link[fixest]{errbar}}.
#'
#' @examples
#'
#' # Sample data illustrating the DiD
#' data(base_did)
#'
#' # Estimation of yearly effect (they are automatically added)
#' est = did_estimate_yearly_effects(y ~ x1 + treat + post, base_did,
#'                                   treat_time = ~treat+period, reference = 5)
#'
#' # Now we plot the results
#' did_plot_yearly_effects(est)
#'
#' # Now with fixed-effects:
#' est_fe = did_estimate_yearly_effects(y ~ x1 | id + period, base_did,
#'                                      treat_time = ~treat+period, reference = 5)
#' did_plot_yearly_effects(est_fe)
#'
#' # you can change the type of SE to be plotted:
#' did_plot_yearly_effects(est_fe, se = "cluster") # default
#' did_plot_yearly_effects(est_fe, se = "standard")
#'
#'
did_estimate_yearly_effects = function(fml, data, treat_time, reference, returnData = FALSE, ..., estfun = feols){
	# fml formula to estimate, if no variable: include the constant
	# data: a data frame
	# treat_time: the names of the treatment and time variables
	# reference: the time period to take as a reference
	# ... other arguments to be passed to femlm
	# estfun: the function to estimate the model

	# Controls and creation of the estimation data
	if("data.table" %in% class(data)){
		data_small = as.data.frame(data)
	} else if(is.data.frame(data)){
		data_small = data
	} else {
		stop("Argument 'data' must be a data.frame!")
	}

    # Formula
    if(missing(fml)) stop("You must provide argument 'fml' (currently it is missing).")
	if(length(fml) != 3) stop("The argument 'fml' must be a two sided formula, e.g. y~x1+x2.")

	# the treat and time variables
    if(missing(treat_time)) stop("You must provide argument 'treat_time' (currently it is missing).")
    if(is.character(treat_time)){
        if(length(treat_time) != 2){
            stop("Argument treat_time must be of length 2 (currently it is of length ", length(treat_time), ").")
        }

        qui_pblm = setdiff(treat_time, names(data))
        if(length(qui_pblm) != 0){
            stop("In argument treat_time, the variable", enumerate_items(qui_pblm, addS = TRUE), " not in the data.")
        }

        treat = data_small[[treat_time[1]]]
        time = data_small[[treat_time[2]]]

        treat_var = treat_time[1]
        time_var = treat_time[2]

    } else if("formula" %in% class(treat_time)){

        treat_time = formula(treat_time)

        if(length(treat_time) != 2) stop("In argument treat_time, the formula must be two sided, e.g. ~treat+year.")

        all_vars = all.vars(treat_time)
        qui_pblm = setdiff(all_vars, names(data))
        if(length(qui_pblm) != 0){
            stop("In argument treat_time, the variable", enumerate_items(qui_pblm, addS = TRUE), " not in the data.")
        }

        t = terms(treat_time, data = data)
        all_var_full = attr(t, "term.labels")
        if(length(all_var_full) != 2) stop("In argument treat_time, the formula must contain exactly two variables (currently it contains ", length(all_var_full), ").")

        treat = eval(parse(text = all_var_full[1]), data_small)
        time = eval(parse(text = all_var_full[2]), data_small)

        treat_var = all_var_full[1]
        time_var = all_var_full[2]

    } else {
        stop("Argument treat_time must be either a character vector, either a one-sided formula.")
    }


	if(!is.numeric(treat) || any(!treat %in% c(0, 1))){
	    obs = head(which(!treat %in% c(0, 1)), 3)
		stop("The treatment variable must be 0/1, with 1 repersenting the treated. The variable that you gave, ", treat_var, ", is not (e.g. observation", enumerate_items(obs, addS = TRUE, endVerb = FALSE), ".")
	}

	all_periods = sort(unique(time[!is.na(time)]))

	if(missing(reference)) stop("You must provide argument 'reference' (currently it is missing).")
	if(!length(reference) == 1 || !isVector(reference)) stop("Argument reference must be a numeric scalar.")
	if(!reference %in% all_periods){
		stop("The 'reference' must be a value of the time variable (currenlty ", reference, " does not belong to the time variable [", time_var, "]).")
	}

	# creating the yearly treatment effects variables
	select_periods = all_periods[all_periods != reference]
	treat_periods = gsub(" |-", "_", paste0("treat_", select_periods))
	for(i in seq_along(select_periods)){
		data_small[[treat_periods[i]]] = treat * (time == select_periods[i])
	}

	# creating the formula + estimation
	FML = Formula::Formula(fml)
	fml_fe = add2fml(FML, treat_periods)

	res = estfun(fml_fe, data_small, ...)

	res$reference = reference
	res$all_periods = select_periods
	res$time_variable = time_var
	attr(res, "isYearlyTreatmentEstimate") = TRUE

	# we tweak a bit the call
	current_call = match.call()
	est_call = res$call
	if(is.null(est_call)){
	    res$call = current_call
	} else {
	    est_call[["data"]] = current_call[["data"]]
	    res$call = est_call
	}

	if(returnData){
		res$data = data_small
	}

	res
}

#' Plots the results of yearly treatment effects estimation
#'
#' This function plots the results of a \code{\link[fixest]{did_estimate_yearly_effects}} estimation.
#'
#' @inheritParams errbar
#'
#' @param object An object returned by the function \code{\link[fixest]{did_estimate_yearly_effects}}.
#' @param ... Any other argument to be passed to \code{summary} or to \code{plot}.
#' @param style One of \code{"interval"} (default), \code{"bar"} or \code{"tube"}. The style of the plot.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{did_estimate_yearly_effects}}, \code{\link[fixest]{errbar}}.
#'
#' @examples
#'
#' # Sample data illustrating the DiD
#' data(base_did)
#'
#' # Estimation of yearly effect (they are automatically added)
#' est = did_estimate_yearly_effects(y ~ x1 + treat + post, base_did,
#'                                   treat_time = ~treat+period, reference = 5)
#'
#' # Now we plot the results
#' did_plot_yearly_effects(est)
#'
#' # Now with fixed-effects:
#' est_fe = did_estimate_yearly_effects(y ~ x1 | id + period, base_did,
#'                                      treat_time = ~treat+period, reference = 5)
#' did_plot_yearly_effects(est_fe)
#'
#' # you can change the type of SE to be plotted:
#' did_plot_yearly_effects(est_fe, se = "cluster") # default
#' did_plot_yearly_effects(est_fe, se = "standard")
#'
#'
did_plot_yearly_effects = function(object, x.shift = 0, w = 0.1, ci_level = 0.95, style = c("bar", "interval", "tube"), add = FALSE, col = 1, bar.col = col, bar.lwd = par("lwd"), bar.lty, grid = TRUE, grid.par = list(lty=1), bar.join = TRUE, ...){
	# object: a result from the estimate_yearly_effects function
	# ... anything to be passed to summary & to plot

    #-------------------#
    # Next in line:
    # when add = TRUE:
    #   * change the y-axis of the new value + add right axis
    #-------------------#

	# Controls
	style = match.arg(style)
	isYearlyTreatmentEstimate = attr(object, "isYearlyTreatmentEstimate")
	if(!isTRUE(isYearlyTreatmentEstimate)){
		stop("The argument 'object' must come from the function 'did_estimate_yearly_effects'.")
	}

	dots = list(...)

	# the coefficients
	if("fixest" %in% class(object)){
	    # more robust to future updates of the package
	    args_name_sum = names(formals("summary.fixest"))
	    args_sum = intersect(names(dots), args_name_sum)
	    if(length(args_sum) == 0){
	        coef_yearly = summary(object, nframes_up = 1)$coeftable
	    } else {
	        dots_sum = dots[args_sum]
	        dots_sum$object = object
	        dots_sum$nframes_up = 1
	        dots[args_sum] = NULL

	        # Rstudio supa slow when error in the following code
	        # (of course only when called from within did_plot_yearly_effects)
	        # => I MUST use try
	        coef_yearly = try(do.call("summary.fixest", dots_sum)$coeftable, silent = TRUE)
	        if("try-error" %in% class(coef_yearly)) stop(coef_yearly)
	        # I haven't found the cause yet
	    }

	} else {

	    sum_exists = FALSE
	    for(c_name in class(object)){
	        if(exists(paste0("summary.", c_name), mode = "function")){
	            sum_exists = TRUE
	            break
	        }
	    }

	    if(!sum_exists) stop("There is no summary method for objects of class ", c_name, ". did_plot_yearly_effects cannot be performed since it applies summary to the object (from which it extracts the coeftable).")

	    fun_name = paste0("summary.", c_name)
	    args_name_sum = names(formals(fun_name))
	    args_sum = intersect(names(dots), args_name_sum)
	    if(length(args_sum) == 0){
	        coef_yearly_sum = summary(object)
	    } else {
	        dots_sum = dots[args_sum]
	        dots_sum$object = object
	        dots[args_sum] = NULL
	        coef_yearly_sum = do.call(fun_name, dots_sum)
	    }

	    if("coeftable" %in% names(coef_yearly_sum)){
	        coef_yearly = coef_yearly_sum$coeftable
	        element = "coeftable"
	    } else if("coefficients" %in% names(coef_yearly_sum)){
	        coef_yearly = coef_yearly_sum$coefficients
	        element = "coefficients"
	    } else {
	        stop("The summary method does not return a coeftable or coefficients object, needed to plot the results.")
	    }

	    if(is.matrix(coef_yearly)){
	        ct_names = colnames(coef_yearly)
	    } else if("data.frame" %in% class(coef_yearly)){
	        ct_names = names(coef_yearly)
	    } else {
	        stop("The element ", element, " from the summary is not a matrix nor a data.frame => did_plot_yearly_effects cannot be performed.")
	    }
	}

	var_names = rownames(coef_yearly)

	var_select = which(grepl("treat_", var_names))
	coef_yearly = coef_yearly[var_select, ]

	# the periods used
	select_periods = object$all_periods
	# the variables name
	time_variable = object$time_variable

	# We add the reference point
	base2show = data.frame(estimate = coef_yearly[, "Estimate"], sd = coef_yearly[, "Std. Error"], x = select_periods)
	base2show = rbind(base2show, data.frame(estimate = 0, sd = 0, x = object$reference))
	base2show = base2show[order(base2show$x), ]

	# Arguments for the errbar and plot functions
	dots$estimate = base2show$estimate
	dots$sd = base2show$sd
	dots$x = base2show$x
	if(is.null(dots$xlab)) dots$xlab = time_variable
	if(is.null(dots$ylab)) dots$ylab = "Estimate of Yearly Treatment"

	# Now the plot
	mc = match.call(expand.dots = FALSE)
	for(v in setdiff(names(mc), c("", "object", "...", "bar.join"))){
	    dots[[v]] = mc[[v]]
	}

	dots$bar.join = bar.join # default different from the one of errbar

	do.call("errbar", dots)

	reference = object$reference
	abline(v=reference, lty = 2, col = "gray")

}


#' Plots confidence intervals
#'
#' This function draws confidence intervals in a graph.
#'
#' @param estimate Numeric vector. The point estimates.
#' @param sd The standard errors of the estimates. It may be missing.
#' @param ci_low If \code{sd} is not provided, the lower bound of the confidence interval. For each estimate.
#' @param ci_top If \code{sd} is not provided, the upper bound of the confidence interval. For each estimate.
#' @param x The value of the x-axis. If missing, the names of the argument \code{estimate} is used.
#' @param x.shift Shifts the confidence intervals bars to the left or right, depending on the value of \code{x.shift}. Default is 0.
#' @param w The width of the confidence intervals.
#' @param ci_level Scalar between 0 and 1: the level of the CI. By default it is equal to 0.95.
#' @param style If \dQuote{interval}: it plots a confidence interval. If \dQuote{bar}, it plots simply error bars. If \dQuote{tube}: as interval, but with a grayed area.
#' @param add Default is \code{FALSE}, if the intervals are to be added to an existing graph. Note that if it is the case, then the argument \code{x} MUST be numeric.
#' @param col Color of the point estimate and of the line joining them (if \code{style = "interval"}).
#' @param bar.col Color of the bars of the confidence interval. Defaults to \code{col}.
#' @param bar.lwd Line width of the confidence intervals, defaults to \code{1}.
#' @param bar.lty Line type of the confidence intervals, defaults to \code{1} for \code{style = "bar"} and to \code{2} for \code{style = "interval"}.
#' @param grid Whether to add an horizontal grid. Default is \code{TRUE}.
#' @param grid.par Graphical parameters used when plotting the grid in the background. Default is \code{list(lty=1)}.
#' @param bar.join Logical, default is \code{FALSE}. Whether to join the dots when \code{style = "bar"}.
#' @param only.params Logical, default is \code{FALSE}. If \code{TRUE}: no graph is plotted, only the \code{ylim} is returned. Useful to stack estimates from different estimations in the same graph.
#' @param ... Other arguments to be passed to the function \code{plot} or \code{lines} (if \code{add = TRUE}).
#'
#' @seealso
#' \code{\link[fixest]{did_estimate_yearly_effects}}, \code{\link[fixest]{did_plot_yearly_effects}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#' a = rnorm(100)
#' b = 0.5*a + rnorm(100)
#' c = -0.5*b + rnorm(100)
#'
#' est = summary(lm(a ~ c + b))
#'
#' errbar(est$coefficients, x.shift = -.2)
#'
#' errbar(est$coefficients, , x.shift = .2, add = TRUE, col = 2, bar.lty = 2, pch=15)
#'
errbar <- function(estimate, sd, ci_low, ci_top, x, x.shift = 0, w=0.1, ci_level = 0.95, style = c("bar", "interval", "tube"), add = FALSE, col = 1, bar.col = col, bar.lwd = par("lwd"), bar.lty, grid = TRUE, grid.par = list(lty=1), bar.join = FALSE, only.params = FALSE, ...){
	# creation de barres d'erreur
	# 1 segment vertical de la taille de l'IC
	# deux barres horizontales, a chaque bornes

    # We extract the point estimate and the SEs if estimate is a matrix/data.frame
    if(is.matrix(estimate)){
        m_names = tolower(colnames(estimate))
        if(m_names[1] == "estimate" & m_names[2] == "std. error"){
            sd = estimate[, 2]
            estimate = estimate[, 1]
        } else {
            stop("Argument estimate is a matrix but does not contain the columns 'Estimate' and 'Std. Error'. Either provide an appropriate matrix or give directly the vector of estimated coefficients in arg. estimate.")
        }
    }

	n <- length(estimate)

	style = match.arg(style)

	if(missing(bar.lty)){
		bar.lty = ifelse(style == "bar", 1, 2)
	}

	if(missing(sd)){
		if(missing(ci_low) | missing(ci_top)) stop("If 'sd' is not provided, you must provide the arguments 'ci_low' and 'ci_top'.")

		ci025 = ci_low
		ci975 = ci_top

	} else {
		if(!missing(ci_low) | !missing(ci_top)) warning("Since 'sd' is provided, arguments 'ci_low' or 'ci_top' are ignored.")

		# We compue the CI
		nb = abs(qnorm((1-ci_level)/2))
		ci975 = estimate + nb*sd
		ci025 = estimate - nb*sd
	}

	# if(add){
	# 	if(missing(x) || !is.numeric(x)){
	# 		stop("When the argument 'add' is used, you MUST provide a numeric argument for 'x'.")
	# 	}
	# }

	# we create x_labels, x_value & x_at
	if(missing(x)){
	    x_at = 1:n
	    x_value = 1:n + x.shift

		if(is.null(names(estimate))){
		    # x = factor(1:n + x.shift, labels = 1:n)
		    x_labels = 1:n
		} else {
		    # x = factor(1:n + x.shift, labels = names(estimate))
		    x_labels = names(estimate)
		}

		my_xlim = range(c(1:n + x.shift, 1:n - x.shift)) + c(-0.5, +0.5)
	} else{
		if(length(x) != n) stop("Argument 'x' must have the same length as 'estimate'.")

	    if(is.numeric(x)){
	        my_xlim = range(c(x + x.shift, x - x.shift))

	        x_at = NULL
	        x_value = x + x.shift
	        x_labels = NULL

	        # x = x + x.shift
	    } else {
	        # x = factor(1:n + x.shift, labels = x)
	        x_at = 1:n
	        x_value = 1:n + x.shift
	        x_labels = x

	        my_xlim = range(c(1:n + x.shift, 1:n - x.shift)) + c(-0.5, +0.5)
	    }
	}

	dots <- list(...)

	# preparation of the do.call
	dots$col = col
	listDefault(dots, "pch", 20)
	listDefault(dots, "xlab", "Variable")
	listDefault(dots, "ylab", "Estimate")
	listDefault(dots, "xlim", my_xlim)
	my_ylim = range(c(ci025, ci975))
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
	    return(list(ylim = my_ylim))
	}

	if(!add){

	    dots$axes = FALSE
	    if(grid){
	        first.par = dots
	        first.par$type = "n"
	        do.call("plot", first.par)

	        do.call("hgrid", grid.par)
	        # we emphasize 0:
	        grid.par$col = "black"
	        grid.par$h = 0
	        listDefault(grid.par, "lwd", 1.5)
	        do.call("abline", grid.par)

	        # now the points or lines
	        if(dots$type != "n"){
	            second.par = dots[c("x", "y", "type", "cex", "col", "pch", "lty", "lwd")]
	            second.par = second.par[lengths(second.par) > 0]
	            do.call("lines", second.par)
	        }
	    } else {
	        do.call("plot", dots)
	    }

	    box()
	    axis(2)
	    axis(1, at = x_at, labels = x_labels)

	} else {
		do.call("lines", dots)
	}

	if(style == "bar" && bar.join){
		# We join the dots
		third.par = dots[c("x", "y", "col", "lty", "lwd")]
		third.par = third.par[lengths(third.par) > 0]
		do.call("lines", third.par)
	}


	if(style == "interval"){
		# the "tube"
		lines(x_value, ci025, lty = bar.lty, lwd = bar.lwd, col = bar.col)
		lines(x_value, ci975, lty = bar.lty, lwd = bar.lwd, col = bar.col)
	} else if(style == "tube"){
		# Here we use shade area
		shade_area(ci025, ci975, x_value, col = "lightgrey", lty=0)

		dots$axes = NULL
		dots$type = "o"
		do.call("lines", dots)
	} else {
		for(i in 1:n){
			# if(is.factor(x)) x = unclass(x)
		    x = x_value

			# a) barre verticale
			segments(x0=x[i], y0=ci025[i], x1=x[i], y1=ci975[i], lwd = bar.lwd, col = bar.col, lty = bar.lty)

			# b)toppings
			#  i) ci975
			segments(x0=x[i]-w, y0=ci975[i], x1=x[i]+w, y1=ci975[i], lwd = bar.lwd, col = bar.col, lty = bar.lty)
			#  ii) ci025
			segments(x0=x[i]-w, y0=ci025[i], x1=x[i]+w, y1=ci025[i], lwd = bar.lwd, col = bar.col, lty = bar.lty)
		}
	}

	res = list(ylim = my_ylim)
	return(invisible(res))
}



#### ................. ####
#### Internal Funs ####
####

results2formattedList = function(..., se, dof = FALSE, cluster, digits=4, pseudo=TRUE, sdBelow=TRUE, dict = NULL, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, subtitles, titles, yesNoFixef = c("Yes", "No"), keepFactors = FALSE, isTex = FALSE, useSummary, dots_call, powerBelow){
    # This function is the core of the functions esttable and esttex

    signifCode = sort(signifCode)
    if(any(signifCode<0) | any(signifCode>1)) stop("The argument 'signifCode' must lie between 0 and 1.")

    if(length(yesNoFixef) != 2) stop("The argument 'yesNoFixef' must be of length 2.")

    # To take care of old verions:
    allowed_types = c("fixest", "femlm", "feNmlm")

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

    if(length(all_models)==0) stop("Not any proper model (fixest) as argument!")

    n_models <- length(all_models)

    # formatting the names (if needed)
    alternative_names = paste0("model ", 1:n_models)
    who2replace = sapply(model_names, function(x) length(x) == 0 || x == "")
    model_names[who2replace] = alternative_names[who2replace]

    # we keep track of the SEs
    se_type_list = list()

    var_list <- coef_list <- coef_below <- sd_below <- list()
    depvar_list <- obs_list <- list()
    r2_list <- aic_list <- bic_list <- loglik_list <- convergence_list <- list()
    sqCor_list = family_list = theta_list = list()

    # To take care of factors
    factorNames = c()
    isFactor = vector(mode = "list", n_models)
    nbFactor = vector(mode = "list", n_models) # the number of items per factor

    slope_names = c()
    slope_flag_list = vector(mode = "list", n_models)

    # if there are subtitles
    if(!missing(subtitles)){
        if(length(subtitles) != n_models){
            stop("If argument 'subtitles' is provided, it must be of the same length as the number of models.")
        } else {
            isSubtitles = TRUE
        }
    } else {
        isSubtitles = FALSE
    }

    # if there are titles
    if(!missing(titles)){
        if(length(titles) != n_models){
            stop("If argument 'titles' is provided, it must be of the same length as the number of models.")
        } else {
            isTitles = TRUE
        }
    } else {
        isTitles = FALSE
    }

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

        nbFactor[[m]] = nbItems

        # Formatting

        lFactor = rep(yesNoFixef[1], length(factor_var))
        names(lFactor) = factor_var
        isFactor[[m]] = lFactor

        factorNames = unique(c(factorNames, factor_var, recursive=TRUE))

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

        if(isTex){
            coef = coefFormatLatex(a[, 1], digits = digits, power = abs(powerBelow))
            se_value = coefFormatLatex(a[, 2], digits = digits, power = abs(powerBelow))
        } else {
            coef = as.character(round(a[, 1], digits))
            se_value = as.character(myRound(a[, 2], digits))
        }

        if(isTex){
            pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(paste0("\\sym{",names(signifCode),"}"), ""))
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

        # statistics
        # Pseudo-R2 // AIC // BIC // N
        n <- nobs(x)
        obs_list[[m]] <- n
        convergence_list[[m]] = ifelse(is.null(x$convStatus), TRUE, x$convStatus)

        K <- x$nparams
        ll <- logLik(x)
        bic_list[[m]] <- round(BIC(x), 3)
        aic_list[[m]] <- round(AIC(x), 3)
        loglik_list[[m]] <- round(logLik(x), 3)
        r2_list[[m]] <- round(as.vector(r2(x, "apr2")), 5)
        sqCor_list[[m]] <- round(as.vector(r2(x, "sq.cor")), 3)

    }

    res = list(se_type_list=se_type_list, var_list=var_list, coef_list=coef_list, coef_below=coef_below, sd_below=sd_below, depvar_list=depvar_list, obs_list=obs_list, r2_list=r2_list, aic_list=aic_list, bic_list=bic_list, loglik_list=loglik_list, convergence_list=convergence_list, sqCor_list=sqCor_list, factorNames=factorNames, isFactor=isFactor, nbFactor=nbFactor, slope_flag_list = slope_flag_list, slope_names=slope_names, useSummary=useSummary, model_names=model_names, family_list=family_list, theta_list=theta_list)

    return(res)
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

    # browser()

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
        dof_value  <- Q / (Q - 1) * (n - 1) / (n - K)
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

    if(attr(t, "intercept") == 1){
        n = nrow(base)
        all_vars_call = parse(text = paste0("list('(Intercept)' = rep(1, ", n, "), ", paste0(all_vars, collapse = ", "), ")"))
        data_list <- eval(all_vars_call, base)
        names(data_list)[-1] = all_var_names
    } else {
        all_vars_call = parse(text = paste0("list(", paste0(all_vars, collapse = ", "), ")"))
        data_list <- eval(all_vars_call, base)
        names(data_list) = all_var_names
    }

    res = do.call("cbind", data_list)


    res
}

terms_fixef = function(fml){
    # separate all terms of fml into fixed effects ans varying slopes

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
                msg = paste0("Square bracket are special characters use **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm, endVerb = FALSE), ").")
                class(msg) = "try-error"
                return(msg)
            }
        }

        if(length(var2check_single) > 0){
            qui = !grepl("\\]$", var2check_single) | lengths(strsplit(var2check_single, "\\[")) != 2
            if(any(qui)){
                item_pblm = var2check_single[qui]
                msg = paste0("Square bracket are special characters use **only** to designate varying slopes (see help). They are currenlty misused (it concerns ", enumerate_items(item_pblm, endVerb = FALSE), ").")
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

prepare_df = function(vars, base, fastCombine = NA){
    # vars: vector of variables to evaluate

    # we drop NAs an make it unique
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

.cleanPCT = function(x){
	# changes % into \% => to escape that character in Latex
	gsub("%", "\\%", x, fixed = TRUE)
	gsub("\\\\%", "\\%", x, fixed = TRUE) # if the user escaped: not done twice
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

enumerate_items = function (x, endVerb = "is", addS = FALSE, past = FALSE, or = FALSE){
	# function that enumerates items and add verbs
	endVerb = match.arg(as.character(endVerb), c("is", "has", "no", "contain", "FALSE"))
	if(endVerb == "FALSE") endVerb = "no"
	n = length(x)

	if(past){
		endWord = switch(endVerb, is = ifelse(n == 1, " was", " were"), no = "", contain = "contained", has="had")
	} else {
		endWord = switch(endVerb, is = ifelse(n == 1, " is", " are"), no = "", contain = ifelse(n == 1, " contains", " contain"), has = ifelse(n == 1, " has", " have"))
	}

	if (addS) {
		startWord = ifelse(n == 1, " ", "s ")
	} else {
		startWord = ""
	}

	if (n == 1) {
		res = paste0(startWord, x, endWord)
	} else {
	    and_or = ifelse(or, " or ", " and ")
		res = paste0(startWord, paste0(x[-n], collapse = ", "), and_or, x[n], endWord)
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

quickUnclassFactor = function(x, addItem = FALSE){
	# does as unclass(as.factor(x))
	# but waaaaay quicker

	if(!is.numeric(x)){
		# level and unclass are slower than applying char2num (about 2 times)
		x = as.character(x)
	}

	if(is.character(x)){
		res = char2num(x, addItem)
		return(res)
	}

	myOrder = order(x)
	x_sorted = x[myOrder]
	x_quf_sorted = cpp_unclassFactor(x_sorted)
	x_quf = x_quf_sorted[order(myOrder)]

	if(addItem){
		res = list(x = x_quf, items = cpp_unik(x_sorted, tail(x_quf_sorted, 1)))
		return(res)
	} else {
		return(x_quf)
	}
}


getItems = function(x){
	# to get the unique elements of x before quickunclassfactor
	# needs to be done because differs depending on the type of x

	if(is.character(x)){
		res = unique(x)
	} else if(is.factor(x)){
		res = levels(unique(x)[, drop=TRUE])
	} else {
		res = sort(unique(x))
	}

	return(res)
}

isVector = function(x){
	# it seems that when you subselect in data.table
	# sometimes it does not yield a vector
	# so i cannot use is.vecyor to check the consistency

	if(is.vector(x)){
		return(TRUE)
	} else {
		# if(class(x) %in% c("integer", "numeric", "character", "factor", "Date", "DateTime") && is.null(dim(x))){
		if(is.null(dim(x)) && !is.list(x)){
			return(TRUE)
		}
	}
	return(FALSE)
}

ifsingle = function(x, yes, no){
	if(length(x) == 1){
		return(yes)
	} else {
		return(no)
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

check_dots_args = function(mc, dots_args = c(), suggest_args = c()){
    # Function to catch the arguments passing in ...
    # we suggest some principal arguments

    fun_name = as.character(mc[[1]])

    args = names(mc$...)
    args = args[nchar(args) > 0]

    args_invalid = setdiff(args, dots_args)
    res = FALSE
    if(length(args_invalid) > 0){
        res = TRUE
        suggest_info = setdiff(suggest_args, names(mc))
        suggest = ""
        if(length(suggest_info) == 1){
            if(length(suggest_args) == 1){
                suggest = paste0(" (fyi, its main argument is ", suggest_info, ".)")
            } else {
                suggest = paste0(" (fyi, another of its main arguments is ", suggest_info, ".)")
            }
        } else if(length(suggest_info) >= 2){
            suggest = paste0(" (fyi, some of its main arguments are ", enumerate_items(suggest_info, endVerb = "no"), ".)")
        }

        msg = paste0(enumerate_items(args_invalid), " not ", ifsingle(args_invalid, "a valid argument", "valid arguments"), " for function ", fun_name, ".", suggest)
        attr(res, "msg") = msg
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

check_arg = function(x, type, message, mustBeThere = TRUE){
    # function that makes it easy to check arguments:
    #  provides precise and meaningful error messages

    type = tolower(type)

    stop_now = function(...){
        # message is a global

        reason = paste0(...)

        print(reason)

        # The original call
        my_call = deparse(sys.calls()[[sys.nframe()-2]])[1] # call can have svl lines
        nmax = 40
        if(nchar(my_call) > nmax) my_call = paste0(substr(my_call, 1, nmax-1), "...")

        # The formatted message
        msg_split = strsplit(message, " ?REASON ?")[[1]]

        msg_new = c(msg_split[1], reason, msg_split[-1])
        msg_new = paste(msg_new, collapse = " ")

        stop("error in ", my_call, ":\n", msg_new, call. = FALSE)
    }

    if(missing(x)){
        if(mustBeThere){
            stop_now("But it is missing.")
        } else {
            return(NULL)
        }
    }

    isSingle = FALSE
    if(grepl("single|scalar", type)){
        isSingle = TRUE
        if(length(x) == 0){
            stop_now("But it is of length 0.")
        } else if(length(x) != 1){
            stop_now("But it is of length ", length(x), ".")
        }
    }

    if(grepl("character", type) && !is.character(x)){
        stop_now("But it is not of type character.")
    }

    if(grepl("logical", type) && !is.logical(x)){
        stop_now("But it is not logical.")
    }

    if(grepl("numeric|integer", type) && !is.numeric(x)){
        stop_now("But it is not numeric.")
    }

    if(!grepl("naok", type) && anyNA(x)){
        if(isSingle){
            stop_now("But it is equal to NA.")
        } else {
            stop_now("But it contains NAs.")
        }
    }

    if(grepl("integer", type) && !all(x %% 1 == 0)){
        stop_now("But it is not integer (although numeric).")
    }

    # Greater than, lower than
    myExtract = function(expr, trim=2){
        # type is global
        start = gregexpr(expr, type)[[1]] + trim
        length = attr(start, "match.length") - trim
        res = substr(type, start, start + length - 1)
        as.numeric(res)
    }

    if(grepl(expr <- "ge[[:digit:]]+", type)){
        n = myExtract(expr)
        if( any(x < n) ) stop_now("But it is lower than ", n, ".")
    }

    if(grepl(expr <- "gt[[:digit:]]+", type)){
        n = myExtract(expr)
        if( any(x == n) ) stop_now("But it is equal to ", n, " (while it should be *striclty* greater).")
        if( any(x < n) ) stop_now("But it is lower than ", n, ".")
    }

    if(grepl(expr <- "le[[:digit:]]+", type)){
        n = myExtract(expr)
        if( !any(x > n) ) stop_now("But it is greater than ", n, ".")
    }

    if(grepl(expr <- "lt[[:digit:]]+", type)){
        n = myExtract(expr)
        if( any(x == n) ) stop_now("But it is equal to ", n, " (while it should be *striclty* lower).")
        if( any(x > n) ) stop_now("But it is greater than ", n, ".")
    }

}

# Avoids the problem of multiple lines deparse
deparse_long = function(x){
    dep_x = deparse(x)
    if(length(dep_x) == 1){
        return(dep_x)
    } else {
        return(paste(gsub("^ +", "", dep_x), collapse = ""))
    }
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
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the cluster coefficients, and the functions \code{\link[fixest]{esttable}} and \code{\link[fixest]{esttex}} to visualize the results of multiple estimations.
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

	n = nrow(newdata)

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
			variable = all.vars(parse(text = fixef_vars[i]))
			isNotHere = !variable %in% names(newdata)
			if(any(isNotHere)){
				stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a cluster variable).")
			}

			# Obtaining the unclassed vector of clusters
			cluster_current = eval(parse(text = fixef_vars[i]), newdata)
			cluster_current_unik = unique(cluster_current)

			fixef_values_possible = attr(object$fixef_id[[i]],"fixef_names")
			valueNotHere = setdiff(cluster_current_unik, fixef_values_possible)
			if(length(valueNotHere) > 0){
				stop("The fixed-effect value ", valueNotHere[1], " (fixed-effect ", fixef_vars[i], ") was not used in the initial estimation, prediction cannot be done for observations with that value. Prediction can be done only for fixed-effect values present in the main estimation.")
			}

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

	}

	# 2) Linear values

	coef = object$coefficients

	value_linear = 0
	rhs_fml = formula(Formula(object$fml), lhs = 0, rhs = 1)
	if(length(all.vars(rhs_fml)) > 0){
		# Checking all variables are there
		varNotHere = setdiff(all.vars(rhs_fml), names(newdata))
		if(length(varNotHere) > 0){
			stop("Some variables used to estimate the model (in fml) are missing from argument 'newdata': ", paste0(varNotHere, collapse = ", "), ".")
		}

		# we create the matrix
		matrix_linear = model.matrix(rhs_fml, newdata)

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
			stop("Some variables used to estimate the model (in the non-linear formula) are missing from argument 'newdata': ", paste0(varNotHere, collapse = ", "), ".")
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
					stop("In the offset, the variable", enumerate_items(varNotHere, addS = TRUE), " not present in 'newdata'.")
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
vcov.fixest = function(object, se, cluster, dof = TRUE, exact_dof = FALSE, forceCovariance = FALSE, keepBounded = FALSE, ...){
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
			myScore = object$score
			object$cov.unscaled = solve(object$hessian)
		} else {
			myScore = object$score[, -which(isBounded), drop = FALSE]
		}
	} else {
		myScore = object$score
	}


	#
	# Core function
	#

	n = object$nobs
	K = object$nparams

	if(exact_dof == TRUE){
	    dof = TRUE
		if(length(object$fixef_id) >= 2){
			fe = fixef(object)
			K = length(object$coefficients) + sum(object$fixef_sizes) - sum(attr(fe, "references"))
		} else {
		    exact_dof = FALSE
		    warning("The argument 'exact_dof' is useful only for estimations with 2 or plus fixed-effects.", call. = TRUE, immediate. = TRUE)
		}
	}

	correction.dof = n / (n - K*dof)

	if(object$method == "feols"){
		if(se.val != "standard"){
			VCOV_raw = object$cov.unscaled / object$sigma2
		} else {
			VCOV_raw = object$cov.unscaled / (n / (n - K))
		}
	} else {
		VCOV_raw = object$cov.unscaled
	}


	# information on the variable used for the clustering
	type_info = ""

	if(anyNA(VCOV_raw)){

		if(!forceCovariance){
			# warning("Standard errors are NA because of likely presence of collinearity. You can use option 'forceCovariance' to try to force the computation of the vcov matrix (to see what's wrong).", call. = FALSE)
			warning("Standard errors are NA because of likely presence of collinearity. Use function collinearity() to detect collinearity problems.", call. = FALSE)
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

		vcov = VCOV_raw * correction.dof

	} else if(se.val == "white"){

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
					stop("Asked for ", nway, "-way clustering but evaluating argument cluster leads to ", length(all_var_names), " clusters (", enumerate_items(all_var_names, endVerb = "no"), "). Please provide exactly ", nway, " clusters.")
				}

				cluster = all_var_names # Now a regular character vector

				doEval = TRUE
			}

			if(length(cluster) == nway && is.character(cluster)){

				if(all(cluster %in% object$fixef_vars)){
					# cluster == names of clusters used in the estimation
					type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")
					cluster = object$fixef_id[cluster]
					do.unclass = FALSE
				} else {
					cluster = gsub(" *", "", cluster)
					if(!doEval){
						is_ok = grepl("^[\\.[:alpha:]][[:alnum:]_\\.]*(\\^[\\.[:alpha:]][[:alnum:]_\\.]*)*$", cluster)
						if(any(!is_ok)){
							stop("In argument cluster, only variable names and the '^' operator are accepted. The expression", enumerate_items(cluster[!is_ok], addS = TRUE), " not valid.\nAlternatively, you can use a list of vectors.", suffix)
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
							stop("Cannot apply ", nway, "-way clustering with current 'cluster' argument. Variable", enumerate_items(var2fetch, past = TRUE, addS = TRUE), " not used as fixed-effects in the estimation. We tried to fetch ", ifelse(length(var2fetch) == 1, "this variable", "these variables"), " in the original database in the parent.frame -- but the data doesn't seem to be there anymore (btw it was ", deparse_long(dataName), "). Alternatively, use a list of vectors.", suffix)
						}

						data = as.data.frame(data)

						# we check the variables are there
						# we use all_vars and not var2fetch: safer to catch all variables (case clustvar^datavar)
						if(any(!all_vars %in% names(data))){
							var_pblm = setdiff(all_vars, names(data))
							stop("In argument 'cluster', the variable", enumerate_items(var_pblm, addS = TRUE), " not present in the original dataset. Alternatively, use a list of vectors.", suffix)
						}

						# we check length consistency
						if(NROW(data) != (object$nobs + length(object$obsRemoved))){
							stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(var2fetch, addS = TRUE, endVerb = "no"), " in the original dataset (", deparse_long(dataName), "), yet the dataset doesn't have the same number of observations as was used in the estimation (", NROW(data), " instead of ", object$nobs + length(object$obsRemoved), ").", suffix)
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
							stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(varsNA, addS = TRUE, endVerb = "no"), " in the original dataset (", deparse_long(dataName), "). But ", ifsingle(varsNA, "this variable", "these variables"), " contain", ifsingle(varsNA, "s", ""), " NAs", msgRemoved, ". This is not allowed.", suffix)
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
					stop("For one way clustering, the argument 'cluster' must be either the name of a cluster variable (e.g. \"dum_1\"), a vector (e.g. data$dum_1), a list containing the vector of clusters (e.g. list(data$dum_1)), or a one-sided formula (e.g. ~dum_1). Currently the class of cluster is ", enumerate_items(class(cluster), endVerb = "no"), ".", suffix)

				}
			} else if(length(cluster) != nway){

				msgList = "a list of vectors"
				if(is.list(cluster)) msgList = "a vector of variables names"
				stop(nway, "-way clustering is asked for, but the length of argument 'cluster' is ", length(cluster), " (it should be ", nway, "). Alternatively, you can use ", msgList, " or a one-sided formula.", suffix)

			} else if(!is.list(cluster)){
				stop("For ", nway, "-way clustering, the argument 'cluster' must be either a vector of cluster variables (e.g. c(\"", paste0("dum_", 1:nway, collapse = "\", \""), "\")), a list containing the vector of clusters (e.g. data[, c(\"", paste0("dum_", 1:nway, collapse = "\", \""), "\")]), or a one-sided formula (e.g. ~", paste0("dum_", 1:nway, collapse = "+"), "). Currently the class of cluster is: ", enumerate_items(class(cluster), endVerb = "no"), ".", suffix)
			}

			cluster = as.list(cluster)
		}

		# now we check the lengths:
		n_per_cluster = sapply(cluster, length)
		if(!all(diff(n_per_cluster) == 0)){
			stop("The vectors of the argument 'cluster' must be of the same length. Currently the lengths are: ", enumerate_items(n_per_cluster, endVerb = "no"), ".")
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
			stop("In argument cluster, the ", enumerate_items(nb_name[varsNA], endVerb = "no"), " cluster variable", ifsingle(varsNA, " contains", "s contain"), " NAs", msgRemoved, ". This is not allowed.")
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

				vcov = vcov + (-1)**(i+1) * vcovClust(index, VCOV_raw, myScore, dof, K, do.unclass=FALSE)

			}
		}
	}

	if(any(diag(vcov)<0)){
		warning("Some variances are negative (likely problem of collinearity).")
	}

	sd.dict = c("standard" = "Standard", "white"="White", "cluster"="Clustered", "twoway"="Two-way", "threeway"="Three-way", "fourway"="Four-way")
	dof_info = ""
	if(exact_dof) dof_info = " [exact dof corr.]"
	if(dof == FALSE) dof_info = " [no dof corr.]"
	attr(vcov, "type") = paste0(as.vector(sd.dict[se.val]), type_info, dof_info)

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
confint.fixest = function(object, parm, level = 0.95, se, cluster, dof = TRUE, ...){

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

		if(is.null(data)){
			dataName = object$call$data
			stop("To apply 'update.fixest', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", deparse_long(dataName), ").")
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

		fixef_new = update(fixef_old, formula(FML, lhs = 0, rhs = 2))

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
model.matrix.fixest = function(object, data, ...){
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
	if(missing(data)){
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

	res = model.matrix(fml, data)

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
#' Sets/gets the default dictionnay used in the function \code{\link[fixest]{esttex}}. The dictionaries are used to relabel variables (usually towards a fancier, more explicit formatting) when exporting them into a Latex table. By setting the dictionary with \code{setFixest_dict}, you can avoid providing the argument \code{dict} in function \code{\link[fixest]{esttex}}.
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
		name_dup = names(dict)[qui]
		stop("Argument 'dict' contains duplicated names: ", enumerate_items(name_dup, endVerb = "no"))
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































































