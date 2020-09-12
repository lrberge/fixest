#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Apr 23 16:41:47 2019
# Purpose: All estimation functions
#----------------------------------------------#



#' Fixed-effects OLS estimation
#'
#' Estimates OLS with any number of fixed-effects.
#'
#' @inheritParams femlm
#'
#' @param fml A formula representing the relation to be estimated. For example: \code{fml = z~x+y}. To include fixed-effects, insert them in this formula using a pipe: e.g. \code{fml = z~x+y | fe_1+fe_2}. You can combine two fixed-effects with \code{^}: e.g. \code{fml = z~x+y|fe_1^fe_2}, see details. You can also use variables with varying slopes using square brackets: e.g. in \code{fml = z~y|fe_1[x] + fe_2} the variable \code{x} will have one coefficient for each value of \code{fe_1} -- if you use varying slopes, please have a look at the details section (can't describe it all here).
#' @param weights A formula or a numeric vector. Each observation can be weighted, the weights must be greater than 0. If equal to a formula, it should be of one-sided: for example \code{~ var_weight}.
#' @param verbose Integer. Higher values give more information. In particular, it can detail the number of iterations in the demeaning algoritmh (the first number is the left-hand-side, the other numbers are the right-hand-side variables).
#' @param demeaned Logical, default is \code{FALSE}. Only used in the presence of fixed-effects: should the centered variables be returned? If \code{TRUE},  it creates the items \code{y_demeaned} and \code{X_demeaned}.
#' @param notes Logical. By default, two notes are displayed: when NAs are removed (to show additional information) and when some observations are removed because of collinearity. To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#'
#' @details
#' The method used to demean each variable along the fixed-effects is based on Berge (2018), since this is the same problem to solve as for the Gaussian case in a ML setup.
#'
#' @section Combining the fixed-effects:
#' You can combine two variables to make it a new fixed-effect using \code{^}. The syntax is as follows: \code{fe_1^fe_2}. Here you created a new variable which is the combination of the two variables fe_1 and fe_2. This is identical to doing \code{paste0(fe_1, "_", fe_2)} but more convenient.
#'
#' Note that pasting is a costly operation, especially for large data sets. Thus, the internal algorithm uses a numerical trick which is fast, but the drawback is that the identity of each observation is lost (i.e. they are now equal to a meaningless number instead of being equal to \code{paste0(fe_1, "_", fe_2)}). These \dQuote{identities} are useful only if you're interested in the value of the fixed-effects (that you can extract with \code{\link[fixest]{fixef.fixest}}). If you're only interested in coefficients of the variables, it doesn't matter. Anyway, you can use \code{combine.quick = FALSE} to tell the internal algorithm to use \code{paste} instead of the numerical trick. By default, the numerical trick is performed only for large data sets.
#'
#' @section Varying slopes:
#' You can add variables with varying slopes in the fixed-effect part of the formula. The syntax is as follows: fixef_var[var1, var2]. Here the variables var1 and var2 will be with varying slopes (one slope per value in fixef_var) and the fixed-effect fixef_var will also be added.
#'
#' To add only the variables with varying slopes and not the fixed-effect, use double square brackets: fixef_var[[var1, var2]].
#'
#' In other words:
#' \itemize{
#'   \item fixef_var[var1, var2] is equivalent to fixef_var + fixef_var[[var1]] + fixef_var[[var2]]
#'   \item fixef_var[[var1, var2]] is equivalent to fixef_var[[var1]] + fixef_var[[var2]]
#' }
#'
#' @section Lagging variables:
#'
#' To use leads/lags of variables in the estimation, you can: i) either provide the argument \code{panel.id}, ii) either set your data set as a panel with the function \code{\link[fixest]{panel}}. Doing either of the two will give you acceess to the lagging functions \code{\link[fixest]{l}} and \code{\link[fixest:l]{f}}.
#'
#' You can provide several leads/lags at once: e.g. if your formula is equal to \code{f(y) ~ l(x, -1:1)}, it means that the dependent variable is equal to the lead of \code{y}, and you will have as explanatory variables the lead of \code{x1}, \code{x1} and the lag of \code{x1}. See the examples in function \code{\link[fixest]{l}} for more details.
#'
#' @section Interactions:
#'
#' You can interact a numeric variable with a "factor-like" variable by using \code{interact(var, fe, ref)}, where \code{fe} is the variable to be interacted with and the argument \code{ref} is a value of \code{fe} taken as a reference (optional). Instead of using the function \code{\link[fixest:i]{interact}}, you can use the alias \code{i(var, fe, ref)} or even the highly specific syntax \code{var::fe(ref)}.
#'
#' It is important to note that *if you do not care about the standard-errors of the interactions*, then you can add interactions in the fixed-effects part of the formula (using the syntax fe[[var]], as explained in the section \dQuote{Varying slopes}).
#'
#' Using this specific way to create interactions leads to a different display of the interacted values in \code{\link[fixest]{etable}} and offers a special representation of the interacted coefficients in the function \code{\link[fixest]{coefplot}}. See examples.
#'
#' The function \code{\link[fixest:i]{interact}} has in fact more arguments, please see details in its associated help page.
#'
#' @section On standard-errors:
#'
#' Standard-errors can be computed in different ways, you can use the arguments \code{se} and \code{dof} in \code{\link[fixest]{summary.fixest}} to define how to compute them. By default, in the presence of fixed-effects, standard-errors are automatically clustered.
#'
#' The following vignette: \href{https://cran.r-project.org/package=fixest/vignettes/standard_errors.html}{On standard-errors} describes in details how the standard-errors are computed in \code{fixest} and how you can replicate standard-errors from other software.
#'
#' You can use the functions \code{\link[fixest]{setFixest_se}} and \code{\link[fixest:dof]{setFixest_dof}} to permanently set the way the standard-errors are computed.
#'
#'
#' @return
#' A \code{fixest} object.
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_full}{(When relevant.) The "full" formula containing the linear part and the fixed-effects.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{multicol}{Logical, if multicollinearity was found.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with the intercept).}
#' \item{ssr_fe_only}{Sum of the squared residuals of the model estimated with fixed-effects only.}
#' \item{ll_null}{The log-likelihood of the null model (containing only with the intercept).}
#' \item{ll_fe_only}{The log-likelihood of the model estimated with fixed-effects only.}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{fitted.values}{The fitted values.}
#' \item{linear.predictors}{The linear predictors.}
#' \item{residuals}{The residuals (y minus the fitted values).}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{cov.unscaled}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects coefficients for each observation.}
#' \item{offset}{(When relevant.) The offset formula.}
#' \item{weights}{(When relevant.) The weights formula.}
#' \item{obsRemoved}{(When relevant.) Vector of observations that were removed because of NA values.}
#' \item{collin.var}{(When relevant.) Vector containing the variables removed because of collinearity.}
#' \item{collin.coef}{(When relevant.) Vector of coefficients, where the values of the variables removed because of collinearity are NA.}
#' \item{y_demeaned}{Only when \code{demeaned = TRUE}: the centered dependent variable.}
#' \item{X_demeaned}{Only when \code{demeaned = TRUE}: the centered explanatory variable.}
#'
#'
#' @seealso
#' See also \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations. For plotting coefficients: see \code{\link[fixest]{coefplot}}.
#'
#' And other estimation methods: \code{\link[fixest]{femlm}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#' @examples
#'
#' #
#' # Just one set of fixed-effects:
#' #
#'
#' res = feols(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#' summary(res)
#'
#' #
#' # Varying slopes:
#' #
#'
#' res = feols(Sepal.Length ~ Petal.Length | Species[Sepal.Width], iris)
#' summary(res)
#'
#' #
#' # Combining the FEs:
#' #
#'
#' base = iris
#' base$fe_2 = rep(1:10, 15)
#' res_comb = feols(Sepal.Length ~ Petal.Length | Species^fe_2, base)
#' summary(res_comb)
#' fixef(res_comb)[[1]]
#'
#' #
#' # Using leads/lags:
#' #
#'
#' data(base_did)
#' # We need to set up the panel with the arg. panel.id
#' est1 = feols(y~l(x1, 0:1), base_did, panel.id = ~id+period)
#' est2 = feols(f(y)~l(x1, -1:1), base_did, panel.id = ~id+period)
#' etable(est1, est2, order = "f", drop="Int")
#'
#' #
#' # Using interactions:
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
#' # Now we can plot the result of the interaction with coefplot
#' coefplot(est_did)
#' # You have many more example in coefplot help
#'
#'
feols = function(fml, data, weights, offset, panel.id, fixef, fixef.tol = 1e-6, fixef.iter = 10000,
                 na_inf.rm = getFixest_na_inf.rm(), nthreads = getFixest_nthreads(),
                 verbose = 0, warn = TRUE, notes = getFixest_notes(), combine.quick,
                 demeaned = FALSE, mem.clean = FALSE, only.env = FALSE, env, ...){

	dots = list(...)

	# 1st: is the call coming from feglm?
	fromGLM = FALSE
	if("X" %in% names(dots)){
		fromGLM = TRUE
		# env is provided by feglm
		X = dots$X
		y = as.vector(dots$y)
		init = dots$means
		correct_0w = dots$correct_0w
	} else {
		time_start = proc.time()

		# we use fixest_env for appropriate controls and data handling
		if(missing(env)){
		    env = try(fixest_env(fml = fml, data = data, weights = weights, offset = offset, panel.id = panel.id, fixef = fixef, fixef.tol = fixef.tol, fixef.iter = fixef.iter, na_inf.rm = na_inf.rm, nthreads = nthreads, verbose = verbose, warn = warn, notes = notes, combine.quick = combine.quick, demeaned = demeaned, mem.clean = mem.clean, origin = "feols", mc_origin = match.call(), ...), silent = TRUE)
		} else if(r <- !is.environment(env) || !isTRUE(env$fixest_env)) {
		    stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
		}

		if("try-error" %in% class(env)){
			stop(format_error_msg(env, "feols"))
		}

		check_arg(only.env, "logical scalar")
		if(only.env){
		    return(env)
		}

		y = get("lhs", env)
		X = get("linear.mat", env)
		nthreads = get("nthreads", env)
		init = 0

		# offset
		offset = get("offset.value", env)
		isOffset = length(offset) > 1
		if(isOffset){
			y = y - offset
		}

		# weights
		weights = get("weights.value", env)
		isWeight = length(weights) > 1
		correct_0w = FALSE

		mem.clean = get("mem.clean", env)

		verbose = get("verbose", env)
		if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")
	}

	onlyFixef = length(X) == 1

	if(fromGLM){
		res = list(coefficients = NA)
	} else {
		res = get("res", env)
	}

	isFixef = get("isFixef", env)
	if(!isFixef){
		# No Fixed-effects
		y_demean = y
		X_demean = X
		res$means = 0
	} else {
		time_demean = proc.time()

		# Number of nthreads
		nthreads = min(nthreads, ncol(X) + 1 - onlyFixef)

		# fixef information
		fixef_sizes = get("fixef_sizes", env)
		fixef_table_vector = get("fixef_table_vector", env)
		fixef_id_list = get("fixef_id_list", env)

		slope_flag = get("slope_flag", env)
		slope_vars = get("slope_variables", env)

		if(mem.clean){
		    # we can't really rm many variables... but gc can be enough
		    # cpp_demean is the most mem intensive bit
		    gc()
		}

		vars_demean <- cpp_demean(y, X, weights, iterMax = fixef.iter,
		                          diffMax = fixef.tol, nb_cluster_all = fixef_sizes,
		                          dum_list = fixef_id_list, tableCluster_vector = fixef_table_vector,
		                          slope_flag = slope_flag, slope_vars = slope_vars,
		                          r_init = init, checkWeight = fromGLM, nthreads = nthreads)

		y_demean = vars_demean$y_demean
		X_demean = vars_demean$X_demean
		res$iterations = vars_demean$iterations
		if(fromGLM){
			res$means = vars_demean$means
		}

		if(mem.clean){
		    rm(vars_demean)
		}

		if(any(slope_flag > 0) && any(res$iterations > 300)){
		    # Maybe we have a convergence problem
		    # This is poorly coded, but it's a temporary fix
		    opt_fe = check_conv(y_demean, X_demean, fixef_id_list, slope_flag, slope_vars, weights)

		    # This is a bit too rough a check but it should catch the most problematic cases
		    if(any(opt_fe > 1e-4)){
		        msg = "There seems to be a convergence problem as regards the variables with varying slopes (in the RHS of your formula). The precision of the estimates may not be great. This is a known issue and there is work underway to solve it.\n As a workaround, you can use the variables with varying slopes as regular variables using the function interact (see ?interact)."

		        res$convStatus = FALSE

		        res$message = paste0("tol: ", signif_plus(fixef.tol), " iter: ", max(res$iterations))

		        if(fromGLM){
		            res$warn_varying_slope = msg
		        } else {
		            warning(msg)
		        }
		    }
		} else if(any(res$iterations >= fixef.iter)){
		    msg = paste0("Demeaning algorithm: Absence of convergence after reaching the maximum number of iterations (", fixef.iter, ").")

		    res$convStatus = FALSE
		    res$message = paste0("Maximum of ", fixef.iter, " iterations reached.")

		    if(fromGLM){
		        res$warn_varying_slope = msg
		    } else {
		        warning(msg)
		    }
		}

		if(verbose >= 1){
			if(length(fixef_sizes) > 1){
				cat("Demeaning in ", (proc.time() - time_demean)[3], "s (iter: ", paste0(c(tail(res$iterations, 1), res$iterations[-length(res$iterations)]), collapse = ", "), ")\n", sep="")
			} else {
				cat("Demeaning in ", (proc.time() - time_demean)[3], "s\n", sep="")
			}
		}
	}

	time_esti = proc.time()

	#
	# Estimation
	#

	if(mem.clean){
	    gc()
	}

	if(!onlyFixef){

	    est = ols_fit(y_demean, X_demean, weights, correct_0w, nthreads)

	    # Corner case: not any relevant variable
	    if(!is.null(est$all_removed)){
	        all_vars = colnames(X)

	        if(isFixef){
	            stop_up(ifsingle(all_vars, "The only variable ", "All variables"), enumerate_items(all_vars, "quote.is", nmax = 3), " collinear with the fixed effects. In such circumstances, the estimation is void.", up = fromGLM)
	        } else {
	            stop_up(ifsingle(all_vars, "The only variable ", "All variables"), enumerate_items(all_vars, "quote.is", nmax = 3), " virtually constant and equal to 0. In such circumstances, the estimation is void.", up = fromGLM)
	        }
	    }

		# Formatting the result
	    coef = est$coefficients
	    names(coef) = colnames(X)[!est$is_excluded]
		res$coefficients = coef
		# Additional stuff
		res$residuals = est$residuals
		res$multicol = est$multicol
		if(fromGLM) res$is_excluded = est$is_excluded

		if(demeaned){
		    res$y_demeaned = y_demean
		    res$X_demeaned = X_demean
		    colnames(res$X_demeaned) = colnames(X)
		}

	} else {
		res$residuals = y_demean
		res$coefficients = coef = NULL
		res$onlyFixef = TRUE
		res$multicol = FALSE

		if(demeaned){
		    res$y_demeaned = y_demean
		}
	}

	time_post = proc.time()
	if(verbose >= 1 && (time_post - time_esti)[3] > 0.05){
		cat("Estimation in ", (time_post - time_esti)[3], "s\n", sep="")
	}

	if(mem.clean){
	    gc()
	}

	if(fromGLM){
		res$fitted.values = y - res$residuals
		if(!onlyFixef){
			res$X_demean = X_demean
		}

		return(res)
	}

	#
	# Post processing
	#

	# Collinearity message
	collin.adj = 0
	if(res$multicol){
	    var_collinear = colnames(X)[est$is_excluded]
	    if(notes){
	        message(ifsingle(var_collinear, "The variable ", "Variables "), enumerate_items(var_collinear, "quote.has", nmax = 3), " been removed because of collinearity (see $collin.var).")
	    }

	    res$collin.var = var_collinear

	    # full set of coeffficients with NAs
	    collin.coef = setNames(rep(NA, ncol(X)), colnames(X))
	    collin.coef[!est$is_excluded] = res$coefficients
	    res$collin.coef = collin.coef

	    if(isFixef){
	        X = X[, !est$is_excluded, drop = FALSE]
	    }
	    X_demean = X_demean[, !est$is_excluded, drop = FALSE]

	    collin.adj = sum(est$is_excluded)
	}

	n = length(y)
	res$nparams = res$nparams - collin.adj
	df_k = res$nparams
	res$nobs = n

	if(isWeight) res$weights = weights

	if(onlyFixef){
		res$fitted.values = res$sumFE = y - res$residuals
	} else {

		# X_beta / fitted / sumFE
		if(isFixef){
			x_beta = cpppar_xbeta(X, coef, nthreads)
			res$sumFE = y - x_beta - res$residuals
			res$fitted.values = x_beta + res$sumFE
		} else {
			res$fitted.values = est$fitted.values
		}

		if(isOffset){
			res$fitted.values = res$fitted.values + offset
		}

		#
		# score + hessian + vcov
		if(isWeight){
			res$scores = (res$residuals * weights) * X_demean
		} else {
			res$scores = res$residuals * X_demean
		}

		res$hessian = est$xwx

		if(isWeight){
			res$sigma2 = cpp_ssq(res$residuals * sqrt(weights)) / (length(y) - df_k)
		} else {
			res$sigma2 = cpp_ssq(res$residuals) / (length(y) - df_k)
		}

		res$cov.unscaled = est$xwx_inv * res$sigma2

		rownames(res$cov.unscaled) = colnames(res$cov.unscaled) = names(coef)

		# se
		se = diag(res$cov.unscaled)
		se[se < 0] = NA
		se = sqrt(se)

		# coeftable
		zvalue <- coef/se
		pvalue <- 2*pt(-abs(zvalue), max(n - df_k, 1))

		coeftable <- data.frame("Estimate"=coef, "Std. Error"=se, "t value"=zvalue, "Pr(>|t|)"=pvalue)
		names(coeftable) <- c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
		row.names(coeftable) <- names(coef)

		attr(se, "type") = attr(coeftable, "type") = "Standard"
		res$coeftable = coeftable
		res$se = se
	}

	# fit stats
	if(!cpp_isConstant(res$fitted.values)){
	    res$sq.cor =  stats::cor(y, res$fitted.values)**2
	} else {
	    res$sq.cor = NA
	}

	res$ssr_null = cpp_ssr_null(y)
	sigma_null = sqrt(res$ssr_null/n)
	res$ll_null = -1/2/sigma_null^2*res$ssr_null - n*log(sigma_null) - n*log(2*pi)/2

	# fixef info
	if(isFixef){
		# For the within R2
		if(!onlyFixef){
			res$ssr_fe_only = cpp_ssq(y_demean)
			sigma = sqrt(res$ssr_fe_only/n)
			res$ll_fe_only = -1/2/sigma^2*res$ssr_fe_only - n*log(sigma) - n*log(2*pi)/2
		}
	}

	# other
	res$fml = get("fml", env)
	res$call = match.call()

	if(verbose >= 3) cat("Post-processing in ", (time_post - time_post)[3], "s\n", sep="")

	class(res) = "fixest"

	res
}

ols_fit = function(y, X, w, correct_0w = FALSE, nthreads){
    # No control here -- done before

    info_products = cpp_sparse_products(X, w, y, correct_0w, nthreads)
    xwx = info_products$XtX
    xwy = info_products$Xty

    multicol = FALSE
    info_inv = cpp_cholesky(xwx)

    if(!is.null(info_inv$all_removed)){
        # Means all variables are collinear! => can happen when using FEs
        return(list(all_removed = TRUE))
    }

    xwx_inv = info_inv$XtX_inv
    is_excluded = info_inv$id_excl

    multicol = any(is_excluded)

    if(multicol){
        beta = as.vector(xwx_inv %*% xwy[!is_excluded])
        fitted.values = cpppar_xbeta(X[, !is_excluded, drop = FALSE], beta, nthreads)
    } else {
        # avoids copies
        beta = as.vector(xwx_inv %*% xwy)
        fitted.values = cpppar_xbeta(X, beta, nthreads)
    }


    residuals = y - fitted.values

    res = list(xwx = xwx, coefficients = beta, fitted.values = fitted.values, xwx_inv = xwx_inv, multicol = multicol, residuals = residuals, is_excluded = is_excluded)

    res
}



check_conv = function(y, X, fixef_id_list, slope_flag, slope_vars, weights){
    # VERY SLOW!!!!
    # IF THIS FUNCTION LASTS => TO BE PORTED TO C++

    # y, X => variables that were demeaned

    # For each variable: we compute the optimal FE coefficient
    # it should be 0 if the algorithm converged

    Q = length(slope_flag)

    nobs = length(y)
    if(length(X) == 1){
        K = 1
    } else {
        K = NCOL(X) + 1
    }

    res = matrix(NA, K, Q)

    for(k in 1:K){
        if(k == 1){
            x = y
        } else {
            x = X[, k - 1]
        }

        for(q in 1:Q){
            fixef_id = fixef_id_list[[q]]

            if(slope_flag[q]){
                index_var = 1:nobs + (cumsum(slope_flag)[q] - 1) * nobs
                var = slope_vars[index_var]

                num = tapply(weights * x * var, fixef_id, sum)
                denom = tapply(weights * var^2, fixef_id, sum)
                res[k, q] = max(abs(num/denom))

            } else {
                res[k, q] = max(abs(tapply(weights * x, fixef_id, mean)))
            }
        }
    }

    res
}


#' Fixed-effects GLM estimations
#'
#' Estimates GLM models with any number of fixed-effects.
#'
#' @inheritParams feols
#' @inheritParams femlm
#' @inheritSection feols Combining the fixed-effects
#' @inheritSection feols Varying slopes
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#'
#' @param family Family to be used for the estimation. Defaults to \code{poisson()}. See \code{\link[stats]{family}} for details of family functions.
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1 (e.g. \code{start = 0}), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients). Default is missing.
#' @param etastart Numeric vector of the same length as the data. Starting values for the linear predictor. Default is missing.
#' @param mustart Numeric vector of the same length as the data. Starting values for the vector of means. Default is missing.
#' @param fixef.tol Precision used to obtain the fixed-effects. Defaults to \code{1e-6}. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations.
#' @param glm.iter Number of iterations of the glm algorithm. Default is 25.
#' @param glm.tol Tolerance level for the glm algorithm. Default is \code{1e-8}.
#' @param y Numeric vector of the dependent variable.
#' @param X Numeric matrix of the regressors.
#' @param fixef_mat Matrix/data.frame of the fixed-effects.
#' @param verbose Integer. Higher values give more information. In particular, it can detail the number of iterations in the demeaning algoritmh (the first number is the left-hand-side, the other numbers are the right-hand-side variables). It can also detail the step-halving algorithm.
#' @param notes Logical. By default, three notes are displayed: when NAs are removed, when some fixed-effects are removed because of only 0 (or 0/1) outcomes, or when a variable is dropped because of collinearity. To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#'
#' @details
#' The core of the GLM are the weighted OLS estimations. These estimations are performed with \code{\link[fixest]{feols}}. The method used to demean each variable along the fixed-effects is based on Berge (2018), since this is the same problem to solve as for the Gaussian case in a ML setup.
#'
#' @return
#' A \code{fixest} object.
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_full}{(When relevant.) The "full" formula containing the linear part and the fixed-effects.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{y}{(When relevant.) The dependent variable (used to compute the within-R2 when fixed-effects are present).}
#' \item{convStatus}{Logical, convergence status of the IRWLS algorithm.}
#' \item{irls_weights}{The weights of the last iteration of the IRWLS algorithm.}
#' \item{obsRemoved}{(When relevant.) Vector of observations that were removed because of NA values or because of only 0/1 outcome within a fixed-effect (depends on the family though).}
#' \item{fixef_removed}{(When relevant.) In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{deviance}{Deviance of the fitted model.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with the intercept).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{linear.predictors}{The linear predictors.}
#' \item{residuals}{The residuals (y minus the fitted values).}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{cov.unscaled}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects coefficients for each observation.}
#' \item{offset}{(When relevant.) The offset formula.}
#' \item{weights}{(When relevant.) The weights formula.}
#' \item{collin.var}{(When relevant.) Vector containing the variables removed because of collinearity.}
#' \item{collin.coef}{(When relevant.) Vector of coefficients, where the values of the variables removed because of collinearity are NA.}
#'
#'
#'
#'
#' @seealso
#' See also \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#' And other estimation methods: \code{\link[fixest]{feols}}, \code{\link[fixest]{femlm}}, \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#'
#' @examples
#'
#' # Default is a poisson model
#' res = feglm(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#'
#' # You could also use fepois
#' res_pois = fepois(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#'
#' # With the fit method:
#' res_fit = feglm.fit(iris$Sepal.Length, iris[, 2:3], iris$Species)
#'
#' # All results are identical:
#' etable(res, res_pois, res_fit)
#'
#'
feglm = function(fml, data, family = "poisson", offset, weights, panel.id, start = NULL,
                 etastart = NULL, mustart = NULL, fixef,
                 fixef.tol = 1e-6, fixef.iter = 10000, glm.iter = 25, glm.tol = 1e-8,
                 na_inf.rm = getFixest_na_inf.rm(), nthreads = getFixest_nthreads(),
                 warn = TRUE, notes = getFixest_notes(), verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    if(missing(weights)) weights = NULL

    time_start = proc.time()

    if(missing(env)){
        env = try(fixest_env(fml=fml, data=data, family = family, offset = offset, weights = weights, panel.id = panel.id, linear.start = start, etastart=etastart, mustart=mustart, fixef = fixef, fixef.tol=fixef.tol, fixef.iter=fixef.iter, glm.iter = glm.iter, glm.tol = glm.tol, na_inf.rm = na_inf.rm, nthreads = nthreads, warn=warn, notes=notes, verbose = verbose, combine.quick = combine.quick, mem.clean = mem.clean, origin = "feglm", mc_origin = match.call(), ...), silent = TRUE)

    } else {
        if(!is.environment(env)) stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not an environment.")
        if(is.null(env$fixest_env)) stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not a 'fixest' environment.")
    }

    if("try-error" %in% class(env)){
        mc = match.call()
        origin = ifelse(is.null(mc$origin), "feglm", mc$origin)
        stop(format_error_msg(env, origin))
    }

    check_arg(only.env, "logical scalar")
    if(only.env){
        return(env)
    }

    verbose = get("verbose", env)
    if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

    # workhorse is feglm.fit (OK if error msg leads to feglm.fit [clear enough])
    res = feglm.fit(env = env)

    res
}



#' @rdname feglm
feglm.fit = function(y, X, fixef_mat, family = "poisson", offset, weights, start = NULL,
                     etastart = NULL, mustart = NULL, fixef.tol = 1e-6, fixef.iter = 10000,
                     glm.iter = 25, glm.tol = 1e-8, na_inf.rm = getFixest_na_inf.rm(),
                     nthreads = getFixest_nthreads(), warn = TRUE, notes = getFixest_notes(), mem.clean = FALSE,
                     verbose = 0, only.env = FALSE, env, ...){

    dots = list(...)

    lean = isTRUE(dots$lean)
    means = 1
    if(!missing(env)){
        # This is an internal call from the function feglm
        # no need to further check the arguments
        # we extract them from the env

        if(r <- !is.environment(env) || !isTRUE(env$fixest_env)) {
            stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
        }

        # main variables
        if(missing(y)) y = get("lhs", env)
        if(missing(X)) X = get("linear.mat", env)
        if(!missing(fixef_mat) && is.null(fixef_mat)){
            assign("isFixef", FALSE, env)
        }

        if(missing(offset)) offset = get("offset.value", env)
        if(missing(weights)) weights = get("weights.value", env)

        # other params
        if(missing(fixef.tol)) fixef.tol = get("fixef.tol", env)
        if(missing(fixef.iter)) fixef.iter = get("fixef.iter", env)
        if(missing(glm.iter)) glm.iter = get("glm.iter", env)
        if(missing(glm.tol)) glm.tol = get("glm.tol", env)
        if(missing(warn)) warn = get("warn", env)
        if(missing(verbose)) verbose = get("verbose", env)

        # starting point of the fixed-effects
        if(!is.null(dots$means)) means = dots$means

        # init
        init.type = get("init.type", env)
        starting_values = get("starting_values", env)
        if(lean){
            # Call within here => either null model or fe only
            init.type = "default"
            if(!is.null(etastart)){
                init.type = "eta"
                starting_values = etastart
            }
        }

    } else {

        if(missing(weights)) weights = NULL

        time_start = proc.time()

        env = try(fixest_env(y = y, X = X, fixef_mat = fixef_mat, family = family, na_inf.rm = na_inf.rm, nthreads = nthreads, offset = offset, weights = weights, linear.start = start, etastart=etastart, mustart=mustart, fixef.tol = fixef.tol, fixef.iter = fixef.iter, glm.iter = glm.iter, glm.tol = glm.tol, notes=notes, mem.clean = mem.clean, warn=warn, verbose = verbose, origin = "feglm.fit", mc_origin = match.call(), ...), silent = TRUE)

        if("try-error" %in% class(env)){
            stop(format_error_msg(env, "feglm.fit"))
        }

        check_arg(only.env, "logical scalar")
        if(only.env){
            return(env)
        }

        verbose = get("verbose", env)
        if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

        # y/X
        y = get("lhs", env)
        X = get("linear.mat", env)

        # offset
        offset = get("offset.value", env)

        # weights
        weights = get("weights.value", env)

        # init
        init.type = get("init.type", env)
        starting_values = get("starting_values", env)
    }

    # Setup:
    family = get("family_funs", env)
    isFixef = get("isFixef", env)
    nthreads = get("nthreads", env)
    isWeight = length(weights) > 1
    isOffset = length(offset) > 1
    nobs <- length(y)
    onlyFixef = length(X) == 1

    # the preformatted results
    res = get("res", env)

    # glm functions:
    variance = family$variance
    linkfun = family$linkfun
    linkinv = family$linkinv
    sum_dev.resids = family$sum_dev.resids
    valideta = family$valideta
    validmu = family$validmu
    mu.eta = family$mu.eta
    family_equiv = family$family_equiv

    #
    # Init
    #

    if(init.type == "mu"){

        mu = starting_values

        if(!valideta(mu)){
            stop("In 'mustart' the values provided are not valid.")
        }

        eta = linkfun(mu)

    } else if(init.type == "eta"){

        eta = starting_values

        if(!valideta(eta)){
            stop("In 'etastart' the values provided are not valid.")
        }

        mu = linkinv(eta)

    } else if(init.type == "coef"){
        # If there are fixed-effects we MUST first compute the FE model with starting values as offset
        #   otherwise we are too far away from the solution and starting values may lead to divergence
        #   (hence step halving would be required)
        # This means that initializing with coefficients incurs large computational costs
        #   with fixed-effects

        start = get("start", env)
        offset_fe = offset + cpppar_xbeta(X, start, nthreads)

        if(isFixef){
            mustart = 0
            eval(family$initialize)
            eta = linkfun(mustart)

            # just a rough estimate (=> high tol values) [no benefit in high precision]
            model_fe = try(feglm.fit(X = 0, etastart = eta, offset = offset_fe, glm.tol = 1e-2, fixef.tol = 1e-2, env = env, lean = TRUE))

            if("try-error" %in% class(model_fe)){
                stop("Estimation failed during initialization when getting the fixed-effects, maybe change the values of 'start'? \n", model_fe)
            }

            eta = model_fe$linear.predictors
            mu = model_fe$fitted.values
            devold = model_fe$deviance
        } else {
            eta = offset_fe
            mu = linkinv(eta)
            devold = sum_dev.resids(y, mu, eta, wt = weights)
        }

        wols_old = list(fitted.values = eta - offset)

    } else {
        mustart = 0
        eval(family$initialize)
        eta = linkfun(mustart)
        mu = linkinv(eta)

        # NOTA: FE only => ADDS LOTS OF COMPUTATIONAL COSTS without convergence benefit
    }

    if(init.type != "coef"){
        # starting deviance with constant equal to 1e-5
        # this is important for getting in step halving early (when deviance goes awry right from the start)
        devold = sum_dev.resids(y, rep(linkinv(1e-5), nobs), rep(1e-5, nobs), wt = weights)
        wols_old = list(fitted.values = rep(1e-5, nobs))
    }

    if(!validmu(mu) || !valideta(eta)){
        stop("Current starting values are not valid.")
    }

    assign("nb_sh", 0, env)
    on.exit(warn_step_halving(env))

    if((init.type == "coef" && verbose >= 1) || verbose >= 4) {
        cat("Deviance at initializat.  = ", numberFormatNormal(devold), "\n", sep = "")
    }

    #
    # The main loop
    #

    wols = list(means = 1)
    conv = FALSE
    warning_msg = div_message = ""
    for (iter in 1:glm.iter) {

        mu.eta.val = mu.eta(mu, eta)
        var_mu = variance(mu)

        # controls
        any_pblm_mu = cpp_any_na_null(var_mu)
        if(any_pblm_mu){
            if (anyNA(var_mu)){
                stop("NAs in V(mu), at iteration ", iter, ".")
            } else if (any(var_mu == 0)){
                stop("0s in V(mu), at iteration ", iter, ".")
            }
        }

        if(anyNA(mu.eta.val)){
            stop("NAs in d(mu)/d(eta), at iteration ", iter, ".")
        }

        if(isOffset){
            z = (eta - offset) + (y - mu)/mu.eta.val
        } else {
            z = eta + (y - mu)/mu.eta.val
        }

        w = as.vector(weights * mu.eta.val**2 / var_mu)

        is_0w = w == 0
        any_0w = any(is_0w)
        if(any_0w && all(is_0w)){
            warning_msg = paste0("No informative observation at iteration ", iter, ".")
            div_message = "No informative observation."
            break
        }

        wols = feols(y = z, X = X, weights = w, means = wols$means, correct_0w = any_0w, env = env, fixef.tol = fixef.tol * 10**(iter==1), fixef.iter = fixef.iter, nthreads = nthreads, mem.clean = mem.clean, verbose = verbose - 1)

        # In theory OLS estimation is guaranteed to exist
        # yet, NA coef may happen with non-infinite very large values of z/w (e.g. values > 1e100)
        if(anyNA(wols$coefficients)){
            if(iter == 1){
                stop("Weighted-OLS returns NA coefficients at first iteration, step halving cannot be performed. Try other starting values?")
            }

            warning_msg = paste0("Divergence at iteration ", iter, ": ", msg, ". Weighted-OLS returns NA coefficients. Last evaluated coefficients with finite deviance are returned for information purposes.")
            div_message = "Weighted-OLS returned NA coefficients."
            wols = wols_old
            break
        }

        eta = wols$fitted.values
        if(isOffset){
            eta = eta + offset
        }

        mu = linkinv(eta)

        dev = sum_dev.resids(y, mu, eta, wt = weights)
        dev_evol = dev - devold

        if(verbose >= 1) cat("Iteration: ", sprintf("%02i", iter), " -- Deviance = ", numberFormatNormal(dev), " -- Evol. = ", dev_evol, "\n", sep = "")

        #
        # STEP HALVING
        #

        if(!is.finite(dev) || dev_evol > 0 || !valideta(eta) || !validmu(mu)){

            if(!is.finite(dev)){
                # we report step-halving but only for non-finite deviances
                # other situations are OK (it just happens)
                nb_sh = get("nb_sh", env)
                assign("nb_sh", nb_sh + 1, env)
            }

            eta_new = wols$fitted.values
            eta_old = wols_old$fitted.values

            iter_sh = 0
            do_exit = FALSE
            while(!is.finite(dev) || dev_evol > 0 || !valideta(eta_new) || !validmu(mu)){

                if(iter == 1 && (is.finite(dev) && valideta(eta_new) && validmu(mu)) && iter_sh >= 2){
                    # BEWARE FIRST ITERATION:
                    # at first iteration, the deviance can be higher than the init, and SH may not help
                    # we need to make sure we get out of SH before it's messed up
                    break
                } else if(iter_sh == glm.iter){

                    # if first iteration => means algo did not find viable solution
                    if(iter == 1){
                        stop("Algorithm failed at first iteration. Step-halving could not find a valid set of parameters.")
                    }

                    # Problem only if the deviance is non-finite or eta/mu not valid
                    # Otherwise, it means that we're at a maximum

                    if(!is.finite(dev) || !valideta(eta_new) || !validmu(mu)){
                        # message
                        msg = ifelse(!is.finite(dev), "non-finite deviance", "no valid eta/mu")

                        warning_msg = paste0("Divergence at iteration ", iter, ": ", msg, ". Step halving: no valid correction found. Last evaluated coefficients with finite deviance are returned for information purposes.")
                        div_message = paste0(msg, " despite step-halving")
                        wols = wols_old
                        do_exit = TRUE
                    }

                    break
                }

                iter_sh = iter_sh + 1
                eta_new = (eta_old + eta_new) / 2

                mu = linkinv(eta_new + offset)
                dev = sum_dev.resids(y, mu, eta_new + offset, wt = weights)
                dev_evol = dev - devold

                if(verbose >= 3) cat("Step-halving: iter =", iter_sh, "-- dev:", numberFormatNormal(dev), "-- evol:", numberFormatNormal(dev_evol), "\n")
            }

            if(do_exit) break

            # it worked: update
            eta = eta_new + offset
            wols$fitted.values = eta_new
            # NOTA: we must NOT end with a step halving => we need a proper weighted-ols estimation
            # we force the algorithm to continue
            dev_evol = Inf

            if(verbose >= 2){
                cat("Step-halving: new deviance = ", numberFormatNormal(dev), "\n", sep = "")
            }

        }

        if(abs(dev_evol)/(0.1 + abs(dev)) < glm.tol){
            conv = TRUE
            break
        } else {
            devold = dev
            wols_old = wols
        }
    }

    # Convergence flag
    if(!conv){
        if(iter == glm.iter){
            warning_msg = paste0("Absence of convergence: Maximum number of iterations reached (", glm.iter, "). Final deviance: ", numberFormatNormal(dev), ".")
            div_message = "no convergence: Maximum number of iterations reached"
        }

        res$convStatus = FALSE
        res$message = div_message
    } else {
        res$convStatus = TRUE
    }

    #
    # post processing
    #

    # repeated warning message later on
    # if(wols$multicol){
    #     warning_msg = paste(warning_msg, "Presence of collinearity. Use function collinearity() to pinpoint the problems.")
    # }

    # Collinearity message
    collin.adj = 0
    if(wols$multicol){
        var_collinear = colnames(X)[wols$is_excluded]
        if(notes) message(ifsingle(var_collinear, "The variable ", "Variables "), enumerate_items(var_collinear, "quote.has"), " been removed because of collinearity (see $collin.var).")

        res$collin.var = var_collinear

        # full set of coeffficients with NAs
        collin.coef = setNames(rep(NA, ncol(X)), colnames(X))
        collin.coef[!wols$is_excluded] = wols$coefficients
        res$collin.coef = collin.coef

        wols$X_demean = wols$X_demean[, !wols$is_excluded, drop = FALSE]
        X = X[, !wols$is_excluded, drop = FALSE]

        collin.adj = sum(wols$is_excluded)
    }


    res$irls_weights = w # weights from the iteratively reweighted least square

    res$coefficients = coef = wols$coefficients

    if(!is.null(wols$warn_varying_slope)){
        warning(wols$warn_varying_slope)
    }

    res$linear.predictors = wols$fitted.values
    if(isOffset){
        res$linear.predictors = res$linear.predictors + offset
    }
    res$fitted.values = linkinv(res$linear.predictors)
    res$residuals = y - res$fitted.values

    if(onlyFixef) res$onlyFixef = onlyFixef

    # dispersion + scores
    if(family$family %in% c("poisson", "binomial")){
        res$dispersion = 1
    } else {
        weighted_resids = wols$residuals * res$irls_weights
        # res$dispersion = sum(weighted_resids ** 2) / sum(res$irls_weights)
        # I use the second line to fit GLM's
        res$dispersion = sum(weighted_resids * wols$residuals) / (res$nobs - res$nparams)
    }

    res$working_residuals = wols$residuals

    if(!onlyFixef && !lean){
        # score + hessian + vcov

        # dispersion + scores
        if(family$family %in% c("poisson", "binomial")){
            res$scores = (wols$residuals * res$irls_weights) * wols$X_demean
            res$hessian = cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads)
        } else {
            res$scores = (weighted_resids / res$dispersion) * wols$X_demean
            res$hessian = cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads) / res$dispersion
        }

        # cov:
        # var <- NULL
        # try(var <- solve(res$hessian), silent = TRUE)
        # if(is.null(var) || wols$multicol){
        #     if(is.null(var)){
        #         warning_msg = paste(warning_msg, "Covariance not defined, presence of collinearity. Use function collinearity() to pinpoint the problems.")
        #         res$cov.unscaled = res$hessian * NA
        #     } else {
        #         warning_msg = paste(warning_msg, "Presence of collinearity in the IRLS stage. Use function collinearity() to pinpoint the problems.")
        #         res$cov.unscaled = var
        #     }
        # } else {
        #     res$cov.unscaled = var
        # }
        info_inv = cpp_cholesky(res$hessian, nthreads)
        if(!is.null(info_inv$all_removed)){
            # This should not occur, but I prefer to be safe
            stop("Not any single variable with a positive variance was found after the weighted-OLS stage. (If possible, could you send a replicable example to fixest's author? He's curious about when that actually happens, since in theory it should never happen.)")
        }

        var = info_inv$XtX_inv
        is_excluded = info_inv$id_excl

        if(any(is_excluded)){
            # There should be no remaining collinearity
            warning_msg = paste(warning_msg, "Residual collinearity was found after the weighted-OLS stage. The covariance is not defined. (This should not happen. If possible, could you send a replicable example to fixest's author? He's curious about when that actually happen.)")
            var = matrix(NA, length(is_excluded), length(is_excluded))
        }
        res$cov.unscaled = var

        rownames(res$cov.unscaled) = colnames(res$cov.unscaled) = names(coef)

        # se
        se = diag(res$cov.unscaled)
        se[se < 0] = NA
        se = sqrt(se)

        # coeftable
        zvalue <- coef/se
        use_t = !family$family %in% c("poisson", "binomial")
        if(use_t){
            pvalue <- 2*pt(-abs(zvalue), max(res$nobs - res$nparams, 1))
            ctable_names = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
        } else {
            pvalue <- 2*pnorm(-abs(zvalue))
            ctable_names = c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
        }

        coeftable <- data.frame("Estimate"=coef, "Std. Error"=se, "z value"=zvalue, "Pr(>|z|)"=pvalue)
        names(coeftable) <- ctable_names
        row.names(coeftable) <- names(coef)

        attr(se, "type") = attr(coeftable, "type") = "Standard"
        res$coeftable = coeftable
        res$se = se
    }

    if(nchar(warning_msg) > 0){
        if(warn){
            warning(warning_msg, call. = FALSE)
            options("fixest_last_warning" = proc.time())
        }
    }

    n = length(y)
    res$nobs = n

    res$nparams = res$nparams - collin.adj
    df_k = res$nparams

    # r2s
    if(!cpp_isConstant(res$fitted.values)){
        res$sq.cor = stats::cor(y, res$fitted.values)**2
    } else {
        res$sq.cor = NA
    }

    # deviance
    res$deviance = dev

    # simpler form for poisson
    if(family_equiv == "poisson"){
        if(isWeight){
            res$loglik = sum( (y * eta - mu - cpppar_lgamma(y + 1, nthreads)) * weights)
        } else {
            lfact = get("lfactorial", env)
            res$loglik = sum(y * eta - mu) - lfact
        }
    } else {
        res$loglik = family$aic(y = y, n = rep.int(1, n), mu = res$fitted.values, wt = weights, dev = dev) / -2
    }

    if(lean){
        return(res)
    }

    # The pseudo_r2
    if(family_equiv %in% c("poisson", "logit")){
        ll_null = env$model0$loglik
        fitted_null = linkinv(env$model0$constant)
    } else {
        if(verbose >= 1) cat("Null model:\n")
        model_null = feglm.fit(X = matrix(1, nrow = n, ncol = 1), fixef_mat = NULL, env = env, lean = TRUE)
        ll_null = model_null$loglik
        fitted_null = model_null$fitted.values
    }
    res$ll_null = ll_null
    res$pseudo_r2 = 1 - (res$loglik - df_k)/(ll_null - 1)
    res$ssr_null = cpp_ssq(y - fitted_null)


    # fixef info
    if(isFixef){
        if(onlyFixef){
            res$sumFE = res$linear.predictors
        } else {
            res$sumFE = res$linear.predictors - cpppar_xbeta(X, res$coefficients, nthreads)
        }

        if(isOffset){
            res$sumFE = res$sumFE - offset
        }

    }

    # other
    res$iterations = iter
    res$family = family
    class(res) = "fixest"

    return(res)
}


#' Fixed-effects maximum likelihood model
#'
#' This function estimates maximum likelihood models with any number of fixed-effects.
#'
#' @inheritParams feNmlm
#' @inherit feNmlm return details
#' @inheritSection feols Combining the fixed-effects
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#'
#' @param fml A formula representing the relation to be estimated. For example: \code{fml = z~x+y}. To include fixed-effects, you can 1) either insert them in this formula using a pipe (e.g. \code{fml = z~x+y|fixef_1+fixef_2}), or 2) either use the argument \code{fixef}.
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1 (e.g. \code{start = 0}, the default), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients).
#'
#' @details
#' Note that the functions \code{\link[fixest]{feglm}} and \code{\link[fixest]{femlm}} provide the same results when using the same families but differ in that the latter is a direct maximum likelihood optimization (so the two can really have different convergence rates).
#'
#' @return
#' A \code{fixest} object.
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_full}{(When relevant.) The "full" formula containing the linear part and the fixed-effects.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{convStatus}{Logical, convergence status.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{obsRemoved}{(When relevant.) In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the row numbers of the observations that were removed. Also reports the NA observations that were removed.}
#' \item{fixef_removed}{(When relevant.) In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{coefficients}{The named vector of estimated coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The log-likelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{ll_fe_only}{Log-likelihood of the model with only the fixed-effects.}
#' \item{ssr_null}{Sum of the squared residuals of the null model (containing only with the intercept).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{residuals}{The residuals (y minus the fitted values).}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{cov.unscaled}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects coefficients for each observation.}
#' \item{offset}{(When relevant.) The offset formula.}
#' \item{weights}{(When relevant.) The weights formula.}
#'
#'
#' @seealso
#' See also \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#' And other estimation methods: \code{\link[fixest]{feols}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest]{feNmlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#' On the unconditionnal Negative Binomial model:
#'
#' Allison, Paul D and Waterman, Richard P, 2002, "Fixed-Effects Negative Binomial Regression Models", Sociological Methodology 32(1) pp. 247--265
#'
#' @examples
#'
#' #
#' # Linear examples
#' #
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 fixed-effects
#' # 1) Poisson estimation
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # 2) Log-Log Gaussian estimation (with same FEs)
#' est_gaus = update(est_pois, log(Euros+1) ~ ., family="gaussian")
#'
#' # Comparison of the results using the function esttable
#' esttable(est_pois, est_gaus)
#' # Now using two way clustered standard-errors
#' esttable(est_pois, est_gaus, se = "twoway")
#'
#' # Comparing different types of standard errors
#' sum_white    = summary(est_pois, se = "white")
#' sum_oneway   = summary(est_pois, se = "cluster")
#' sum_twoway   = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' esttable(sum_white, sum_oneway, sum_twoway, sum_threeway)
#'
#'
#'
#'
femlm <- function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), start = 0, fixef,
						offset, panel.id, na_inf.rm = getFixest_na_inf.rm(),
						fixef.tol = 1e-5, fixef.iter = 10000,
						nthreads = getFixest_nthreads(), verbose = 0, warn = TRUE,
						notes = getFixest_notes(), theta.init, combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

	# This is just an alias

	res = try(feNmlm(fml=fml, data=data, family=family, fixef=fixef, offset=offset, panel.id = panel.id, start = start, na_inf.rm=na_inf.rm, fixef.tol=fixef.tol, fixef.iter=fixef.iter, nthreads=nthreads, verbose=verbose, warn=warn, notes=notes, theta.init = theta.init, combine.quick = combine.quick, mem.clean = mem.clean, origin="femlm", mc_origin_bis=match.call(), only.env=only.env, env=env, ...), silent = TRUE)

	if("try-error" %in% class(res)){
		stop(format_error_msg(res, "femlm"))
	}

	return(res)
}

#' @rdname femlm
fenegbin = function(fml, data, theta.init, start = 0, fixef, offset, panel.id,
                    na_inf.rm = getFixest_na_inf.rm(), fixef.tol = 1e-5,
                    fixef.iter = 10000, nthreads = getFixest_nthreads(),
                    verbose = 0, warn = TRUE, notes = getFixest_notes(), combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    # We control for the problematic argument family
    if("family" %in% names(match.call())){
        stop("Function fenegbin does not accept the argument 'family'.")
    }

    # This is just an alias

    res = try(feNmlm(fml = fml, data=data, family = "negbin", theta.init = theta.init, start = start, fixef = fixef, offset = offset, panel.id = panel.id, na_inf.rm = na_inf.rm, fixef.tol = fixef.tol, fixef.iter = fixef.iter, nthreads = nthreads, verbose = verbose, warn = warn, notes = notes, combine.quick = combine.quick, mem.clean = mem.clean, origin = "fenegbin", mc_origin_bis = match.call(), only.env=only.env, env=env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "fenegbin"))
    }

    return(res)
}

#' @rdname feglm
fepois = function(fml, data, offset, weights, panel.id, start = NULL, etastart = NULL, mustart = NULL,
                  fixef, fixef.tol = 1e-6, fixef.iter = 10000, glm.iter = 25, glm.tol = 1e-8,
                  na_inf.rm = getFixest_na_inf.rm(), nthreads = getFixest_nthreads(),
                  warn = TRUE, notes = getFixest_notes(), verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    # We control for the problematic argument family
    if("family" %in% names(match.call())){
        stop("Function fepois does not accept the argument 'family'.")
    }

    # This is just an alias

    res = try(feglm(fml = fml, data = data, family = "poisson", offset = offset, weights = weights, panel.id = panel.id, start = start, etastart = etastart, mustart = mustart, fixef = fixef, fixef.tol = fixef.tol, fixef.iter = fixef.iter, glm.iter = glm.iter, glm.tol = glm.tol, na_inf.rm = na_inf.rm, nthreads = nthreads, warn = warn, notes = notes, verbose = verbose, combine.quick = combine.quick, mem.clean = mem.clean, origin_bis = "fepois", mc_origin_bis = match.call(), only.env=only.env, env=env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "fepois"))
    }

    return(res)
}



#' Fixed effects nonlinear maximum likelihood models
#'
#' This function estimates maximum likelihood models (e.g., Poisson or Logit) with non-linear in parameters right-hand-sides and is efficient to handle any number of fixed effects. If you do not use non-linear in parameters right-hand-side, use \code{\link[fixest]{femlm}} or \code{\link[fixest]{feglm}} instead (design is simpler).
#'
#' @inheritParams panel
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#'
#' @param fml A formula. This formula gives the linear formula to be estimated (it is similar to a \code{lm} formula), for example: \code{fml = z~x+y}. To include fixed-effects variables, you can 1) either insert them in this formula using a pipe (e.g. \code{fml = z~x+y|fixef_1+fixef_2}), or 2) either use the argument \code{fixef}. To include a non-linear in parameters element, you must use the argment \code{NL.fml}.
#' @param start Starting values for the coefficients in the linear part (for the non-linear part, use NL.start). Can be: i) a numeric of length 1 (e.g. \code{start = 0}, the default), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients).
#' @param NL.fml A formula. If provided, this formula represents the non-linear part of the right hand side (RHS). Note that contrary to the \code{fml} argument, the coefficients must explicitly appear in this formula. For instance, it can be \code{~a*log(b*x + c*x^3)}, where \code{a}, \code{b}, and \code{c} are the coefficients to be estimated. Note that only the RHS of the formula is to be provided, and NOT the left hand side.
#' @param data A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this \code{data.frame} names. Can also be a matrix.
#' @param family Character scalar. It should provide the family. The possible values are "poisson" (Poisson model with log-link, the default), "negbin" (Negative Binomial model with log-link), "logit" (LOGIT model with log-link), "gaussian" (Gaussian model).
#' @param fixef Character vector. The names of variables to be used as fixed-effects. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier). Note that the recommended way to include fixed-effects is to insert them directly in the formula.
#' @param na_inf.rm Logical, default is \code{TRUE}. If the variables necessary for the estimation contain NA/Infs and \code{na_inf.rm = TRUE}, then all observations containing NA are removed prior to estimation and a note is displayed detailing the number of observations removed. Otherwise, an error is raised.
#' @param NL.start (For NL models only) A list of starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{NL.start=list(a=1,b=5,c=0)}. Though, there is an exception: if all parameters are to be given the same starting value, you can use a numeric scalar.
#' @param lower (For NL models only) A list. The lower bound for each of the non-linear parameters that requires one. Example: \code{lower=list(b=0,c=0)}. Beware, if the estimated parameter is at his lower bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param upper (For NL models only) A list. The upper bound for each of the non-linear parameters that requires one. Example: \code{upper=list(a=10,c=50)}. Beware, if the estimated parameter is at his upper bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param NL.start.init (For NL models only) Numeric scalar. If the argument \code{NL.start} is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to \code{NL.start.init}.
#' @param offset A formula or a numeric vector. An offset can be added to the estimation. If equal to a formula, it should be of the form (for example) \code{~0.5*x**2}. This offset is linearly added to the elements of the main formula 'fml'.
#' @param jacobian.method (For NL models only) Character scalar. Provides the method used to numerically compute the Jacobian of the non-linear part. Can be either \code{"simple"} or \code{"Richardson"}. Default is \code{"simple"}. See the help of \code{\link[numDeriv]{jacobian}} for more information.
#' @param useHessian Logical. Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.
#' @param hessian.args List of arguments to be passed to function \code{\link[numDeriv]{genD}}. Defaults is missing. Only used with the presence of \code{NL.fml}.
#' @param opt.control List of elements to be passed to the optimization method \code{\link[stats]{nlminb}}. See the help page of \code{\link[stats]{nlminb}} for more information.
#' @param nthreads Integer: Number of nthreads to be used (accelerates the algorithm via the use of openMP routines). The default is to use the total number of nthreads available minus two. You can set permanently the number of threads used within this package using the function \code{\link[fixest]{setFixest_nthreads}}.
#' @param verbose Integer, default is 0. It represents the level of information that should be reported during the optimisation process. If \code{verbose=0}: nothing is reported. If \code{verbose=1}: the value of the coefficients and the likelihood are reported. If \code{verbose=2}: \code{1} + information on the computing time of the null model, the fixed-effects coefficients and the hessian are reported.
#' @param theta.init Positive numeric scalar. The starting value of the dispersion parameter if \code{family="negbin"}. By default, the algorithm uses as a starting value the theta obtained from the model with only the intercept.
#' @param fixef.tol Precision used to obtain the fixed-effects. Defaults to \code{1e-5}. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations. Argument \code{fixef.tol} cannot be lower than \code{10000*.Machine$double.eps}. Note that this parameter is dynamically controlled by the algorithm.
#' @param fixef.iter Maximum number of iterations in fixed-effects algorithm (only in use for 2+ fixed-effects). Default is 10000.
#' @param deriv.iter Maximum number of iterations in the algorithm to obtain the derivative of the fixed-effects (only in use for 2+ fixed-effects). Default is 1000.
#' @param deriv.tol Precision used to obtain the fixed-effects derivatives. Defaults to \code{1e-4}. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations. Argument \code{deriv.tol} cannot be lower than \code{10000*.Machine$double.eps}.
#' @param warn Logical, default is \code{TRUE}. Whether warnings should be displayed (concerns warnings relating to convergence state).
#' @param notes Logical. By default, two notes are displayed: when NAs are removed (to show additional information) and when some observations are removed because of only 0 (or 0/1) outcomes in a fixed-effect setup (in Poisson/Neg. Bin./Logit models). To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#' @param combine.quick Logical. When you combine different variables to transform them into a single fixed-effects you can do e.g. \code{y ~ x | paste(var1, var2)}. The algorithm provides a shorthand to do the same operation: \code{y ~ x | var1^var2}. Because pasting variables is a costly operation, the internal algorithm may use a numerical trick to hasten the process. The cost of doing so is that you lose the labels. If you are interested in getting the value of the fixed-effects coefficients after the estimation, you should use \code{combine.quick = FALSE}. By default it is equal to \code{FALSE} if the number of observations is lower than 50,000, and to \code{TRUE} otherwise.
#' @param only.env (Advanced users.) Logical, default is \code{FALSE}. If \code{TRUE}, then only the environment used to make the estimation is returned.
#' @param mem.clean Logical, default is \code{FALSE}. Only to be used if the data set is large compared to the available RAM. If \code{TRUE} then intermediary objects are removed as much as possible and \code{\link[base]{gc}} is run before each substantial C++ section in the internal code to avoid memory issues.
#' @param env (Advanced users.) A \code{fixest} environment created by a \code{fixest} estimation with \code{only.env = TRUE}. Default is missing. If provided, the data from this environment will be used to perform the estimation.
#' @param ... Not currently used.
#'
#' @details
#' This function estimates maximum likelihood models where the conditional expectations are as follows:
#'
#' Gaussian likelihood:
#' \deqn{E(Y|X)=X\beta}{E(Y|X) = X*beta}
#' Poisson and Negative Binomial likelihoods:
#' \deqn{E(Y|X)=\exp(X\beta)}{E(Y|X) = exp(X*beta)}
#' where in the Negative Binomial there is the parameter \eqn{\theta}{theta} used to model the variance as \eqn{\mu+\mu^2/\theta}{mu+mu^2/theta}, with \eqn{\mu}{mu} the conditional expectation.
#' Logit likelihood:
#' \deqn{E(Y|X)=\frac{\exp(X\beta)}{1+\exp(X\beta)}}{E(Y|X) = exp(X*beta) / (1 + exp(X*beta))}
#'
#' When there are one or more fixed-effects, the conditional expectation can be written as:
#' \deqn{E(Y|X) = h(X\beta+\sum_{k}\sum_{m}\gamma_{m}^{k}\times C_{im}^{k}),}
#' where \eqn{h(.)} is the function corresponding to the likelihood function as shown before. \eqn{C^k} is the matrix associated to fixed-effect dimension \eqn{k} such that \eqn{C^k_{im}} is equal to 1 if observation \eqn{i} is of category \eqn{m} in the fixed-effect dimension \eqn{k} and 0 otherwise.
#'
#' When there are non linear in parameters functions, we can schematically split the set of regressors in two:
#' \deqn{f(X,\beta)=X^1\beta^1 + g(X^2,\beta^2)}
#' with first a linear term and then a non linear part expressed by the function g. That is, we add a non-linear term to the linear terms (which are \eqn{X*beta} and the fixed-effects coefficients). It is always better (more efficient) to put into the argument \code{NL.fml} only the non-linear in parameter terms, and add all linear terms in the \code{fml} argument.
#'
#' To estimate only a non-linear formula without even the intercept, you must exclude the intercept from the linear formula by using, e.g., \code{fml = z~0}.
#'
#' The over-dispersion parameter of the Negative Binomial family, theta, is capped at 10,000. If theta reaches this high value, it means that there is no overdispersion.
#'
#' @return
#' A \code{fixest} object.
#' \item{coefficients}{The named vector of coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{nobs}{The number of observations.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{call}{The call.}
#' \item{fml}{The linear formula of the call.}
#' \item{ll_null}{Log-likelihood of the null model (i.e. with the intercept only).}
#' \item{pseudo_r2}{The adjusted pseudo R2.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{sq.cor}{Squared correlation between the dependent variable and the expected predictor (i.e. fitted.values) obtained by the estimation.}
#' \item{hessian}{The Hessian of the parameters.}
#' \item{fitted.values}{The fitted values are the expected value of the dependent variable for the fitted model: that is \eqn{E(Y|X)}.}
#' \item{cov.unscaled}{The variance-covariance matrix of the parameters.}
#' \item{se}{The standard-error of the parameters.}
#' \item{scores}{The matrix of the scores (first derivative for each observation).}
#' \item{family}{The ML family that was used for the estimation.}
#' \item{residuals}{The difference between the dependent variable and the expected predictor.}
#' \item{sumFE}{The sum of the fixed-effects for each observation.}
#' \item{offset}{The offset formula.}
#' \item{NL.fml}{The nonlinear formula of the call.}
#' \item{bounds}{Whether the coefficients were upper or lower bounded. -- This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
#' \item{isBounded}{The logical vector that gives for each coefficient whether it was bounded or not. This can only be the case when a non-linear formula is included and the arguments 'lower' or 'upper' are provided.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{obsRemoved}{In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the row numbers of the observations that were removed. Also reports the NA observations that were removed.}
#' \item{fixef_removed}{In the case there were fixed-effects and some observations were removed because of only 0/1 outcome within a fixed-effect, it gives the list (for each fixed-effect dimension) of the fixed-effect identifiers that were removed.}
#' \item{theta}{In the case of a negative binomial estimation: the overdispersion parameter.}
#'
#'  @seealso
#' See also \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients, and the function \code{\link[fixest]{etable}} to visualize the results of multiple estimations.
#'
#' And other estimation methods: \code{\link[fixest]{feols}}, \code{\link[fixest]{femlm}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest:femlm]{fenegbin}}.
#'
#' @author
#' Laurent Berge
#'
#' @references
#'
#' Berge, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#' For models with multiple fixed-effects:
#'
#' Gaure, Simen, 2013, "OLS with multiple high dimensional category variables", Computational Statistics & Data Analysis 66 pp. 8--18
#'
#' On the unconditionnal Negative Binomial model:
#'
#' Allison, Paul D and Waterman, Richard P, 2002, "Fixed-Effects Negative Binomial Regression Models", Sociological Methodology 32(1) pp. 247--265
#'
#' @examples
#'
#' # This section covers only non-linear in parameters examples
#' # For linear relationships: use femlm or feglm instead
#'
#' # Generating data for a simple example
#' set.seed(1)
#' n = 100
#' x = rnorm(n, 1, 5)**2
#' y = rnorm(n, -1, 5)**2
#' z1 = rpois(n, x*y) + rpois(n, 2)
#' base = data.frame(x, y, z1)
#'
#' # Estimating a 'linear' relation:
#' est1_L = femlm(z1 ~ log(x) + log(y), base)
#' # Estimating the same 'linear' relation using a 'non-linear' call
#' est1_NL = feNmlm(z1 ~ 1, base, NL.fml = ~a*log(x)+b*log(y), NL.start = list(a=0, b=0))
#' # we compare the estimates with the function esttable (they are identical)
#' esttable(est1_L, est1_NL)
#'
#' # Now generating a non-linear relation (E(z2) = x + y + 1):
#' z2 = rpois(n, x + y) + rpois(n, 1)
#' base$z2 = z2
#'
#' # Estimation using this non-linear form
#' est2_NL = feNmlm(z2~0, base, NL.fml = ~log(a*x + b*y),
#'                NL.start = list(a=1, b=2), lower = list(a=0, b=0))
#' # we can't estimate this relation linearily
#' # => closest we can do:
#' est2_L = femlm(z2~log(x)+log(y), base)
#'
#' # Difference between the two models:
#' esttable(est2_L, est2_NL)
#'
#' # Plotting the fits:
#' plot(x, z2, pch = 18)
#' points(x, fitted(est2_L), col = 2, pch = 1)
#' points(x, fitted(est2_NL), col = 4, pch = 2)
#'
#'
feNmlm = function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml, fixef, na_inf.rm = getFixest_na_inf.rm(), NL.start, lower, upper, NL.start.init, offset, panel.id, start = 0, jacobian.method="simple", useHessian = TRUE, hessian.args = NULL, opt.control = list(), nthreads = getFixest_nthreads(), verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000, deriv.tol = 1e-4, deriv.iter = 1000, warn = TRUE, notes = getFixest_notes(), combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

	time_start = proc.time()

	if(missing(env)){
	    env = try(fixest_env(fml=fml, data=data, family=family, NL.fml=NL.fml, fixef=fixef, na_inf.rm=na_inf.rm, NL.start=NL.start, lower=lower, upper=upper, NL.start.init=NL.start.init, offset=offset, panel.id=panel.id, linear.start=start, jacobian.method=jacobian.method, useHessian=useHessian, opt.control=opt.control, nthreads=nthreads, verbose=verbose, theta.init=theta.init, fixef.tol=fixef.tol, fixef.iter=fixef.iter, deriv.iter=deriv.iter, warn=warn, notes=notes, combine.quick=combine.quick, mem.clean = mem.clean, mc_origin=match.call(), computeModel0=TRUE, ...), silent = TRUE)

	} else if(r <- !is.environment(env) || !isTRUE(env$fixest_env)) {
	    stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
	}

	check_arg(only.env, "logical scalar")
	if(only.env){
	    return(env)
	}


	if("try-error" %in% class(env)){
	    mc = match.call()
	    origin = ifelse(is.null(mc$origin), "feNmlm", mc$origin)
		stop(format_error_msg(env, origin))
	}

	verbose = get("verbose", env)
	if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

	# Objects needed for optimization + misc
	start = get("start", env)
	lower = get("lower", env)
	upper = get("upper", env)
	gradient = get("gradient", env)
	hessian = get("hessian", env)
	family = get("family", env)
	isLinear = get("isLinear", env)
	isNonLinear = get("isNL", env)
	opt.control = get("opt.control", env)

	lhs = get("lhs", env)

	family = get("family", env)
	famFuns = get("famFuns", env)
	model0 = get("model0", env)
	params = get("params", env)

	isFixef = get("isFixef", env)
	onlyFixef = get("onlyFixef", env)

	# the result
	res = get("res", env)

	# NO VARIABLE -- ONLY FIXED-EFFECTS
	if(onlyFixef){
		if(family == "negbin"){
			stop("To estimate the negative binomial model, you need at least one variable. (The estimation of the model with only the fixed-effects is not implemented.)")
		}

		res = femlm_only_clusters(env)
		res$call = match.call()
		res$onlyFixef = TRUE

		return(res)
	}

	# warnings => to avoid accumulation, but should appear even if the user stops the algorithm
	on.exit(warn_fixef_iter(env))

	#
	# Maximizing the likelihood
	#

	opt <- try(stats::nlminb(start=start, objective=femlm_ll, env=env, lower=lower, upper=upper, gradient=gradient, hessian=hessian, control=opt.control), silent = TRUE)

	if("try-error" %in% class(opt)){
		# We return the coefficients (can be interesting for debugging)
		iter = get("iter", env)
		origin = get("origin", env)
		warning_msg = paste0("[", origin, "] Optimization failed at iteration ", iter, ". Reason: ", gsub("^[^\n]+\n *(.+\n)", "\\1", opt))
		if(!"coef_evaluated" %in% names(env)){
			# big problem right from the start
			stop(warning_msg)
		} else {
			coef = get("coef_evaluated", env)
			warning(warning_msg, " Last evaluated coefficients returned.", call. = FALSE)
			return(coef)
		}

	} else {
		convStatus = TRUE
		warning_msg = ""
		if(!opt$message %in% c("X-convergence (3)", "relative convergence (4)", "both X-convergence and relative convergence (5)")){
			warning_msg = " The optimization algorithm did not converge, the results are not reliable."
			convStatus = FALSE
		}

		coef <- opt$par
	}



	# The Hessian
	hessian = femlm_hessian(coef, env = env)
	# we add the names of the non linear variables in the hessian
	if(isNonLinear || family == "negbin"){
		dimnames(hessian) = list(params, params)
	}

	# we create the Hessian without the bounded parameters
	hessian_noBounded = hessian

	# Handling the bounds
	if(!isNonLinear){
		NL.fml = NULL
		bounds = NULL
		isBounded = NULL
	} else {
	    nonlinear.params = get("nonlinear.params", env)
		# we report the bounds & if the estimated parameters are bounded
		upper_bound = upper[nonlinear.params]
		lower_bound = lower[nonlinear.params]

		# 1: are the estimated parameters at their bounds?
		coef_NL = coef[nonlinear.params]
		isBounded = rep(FALSE, length(params))
		isBounded[1:length(coef_NL)] = (coef_NL == lower_bound) | (coef_NL == upper_bound)

		# 2: we save the bounds
		upper_bound_small = upper_bound[is.finite(upper_bound)]
		lower_bound_small = lower_bound[is.finite(lower_bound)]
		bounds = list()
		if(length(upper_bound_small) > 0) bounds$upper = upper_bound_small
		if(length(lower_bound_small) > 0) bounds$lower = lower_bound_small
		if(length(bounds) == 0){
			bounds = NULL
		}

		# 3: we update the Hessian (basically, we drop the bounded element)
		if(any(isBounded)){
			hessian_noBounded = hessian[-which(isBounded), -which(isBounded), drop = FALSE]

			boundText = ifelse(coef_NL == upper_bound, "Upper bounded", "Lower bounded")[isBounded]

			attr(isBounded, "type") = boundText
		}

	}

	# Variance

	var <- NULL
	try(var <- solve(hessian_noBounded), silent = TRUE)
	if(is.null(var)){
		warning_msg = paste(warning_msg, "The information matrix is singular: presence of collinearity. Use function collinearity() to pinpoint the problems.")
		var = hessian_noBounded * NA
		se = diag(var)
	} else {
		se = diag(var)
		se[se < 0] = NA
		se = sqrt(se)
	}

	# Warning message
	if(nchar(warning_msg) > 0){
		if(warn){
		    warning("[femlm]:", warning_msg, call. = FALSE)
		    options("fixest_last_warning" = proc.time())
		}
	}

	# To handle the bounded coefficient, we set its SE to NA
	if(any(isBounded)){
		se = se[params]
		names(se) = params
	}

	zvalue <- coef/se
	pvalue <- 2*pnorm(-abs(zvalue))

	# We add the information on the bound for the se & update the var to drop the bounded vars
	se_format = se
	if(any(isBounded)){
		se_format[!isBounded] = decimalFormat(se_format[!isBounded])
		se_format[isBounded] = boundText
	}

	coeftable <- data.frame("Estimate"=coef, "Std. Error"=se_format, "z value"=zvalue, "Pr(>|z|)"=pvalue, stringsAsFactors = FALSE)
	names(coeftable) <- c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
	row.names(coeftable) <- params

	attr(se, "type") = attr(coeftable, "type") = "Standard"

	mu_both = get_mu(coef, env, final = TRUE)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# calcul pseudo r2
	loglik <- -opt$objective # moins car la fonction minimise
	ll_null <- model0$loglik

	# dummies are constrained, they don't have full dof (cause you need to take one value off for unicity)
	# this is an approximation, in some cases there can be more than one ref. But good approx.
	nparams = res$nparams
	pseudo_r2 <- 1 - (loglik - nparams + 1) / ll_null

	# Calcul residus
	expected.predictor = famFuns$expected.predictor(mu, exp_mu, env)
	residuals = lhs - expected.predictor

	# calcul squared corr
	if(cpp_isConstant(expected.predictor)){
		sq.cor = NA
	} else {
		sq.cor = stats::cor(lhs, expected.predictor)**2
	}
	ssr_null = cpp_ssr_null(lhs)

	# The scores
	scores = femlm_scores(coef, env)
	if(isNonLinear){
		# we add the names of the non linear params in the score
		colnames(scores) = params
	}

	n = length(lhs)

	# Saving
	res$coefficients = coef
	res$coeftable = coeftable
	res$loglik = loglik
	res$iterations = opt$iterations
	res$ll_null = ll_null
	res$ssr_null = ssr_null
	res$pseudo_r2 = pseudo_r2
	res$message = opt$message
	res$convStatus = convStatus
	res$sq.cor = sq.cor
	res$fitted.values = expected.predictor
	res$hessian = hessian
	res$cov.unscaled = var
	res$se = se
	res$scores = scores
	res$family = family
	res$residuals = residuals

	# The value of mu (if cannot be recovered from fitted())
	if(family == "logit"){
		qui_01 = expected.predictor %in% c(0, 1)
		if(any(qui_01)){
			res$mu = mu
		}
	} else if(family %in% c("poisson", "negbin")){
		qui_0 = expected.predictor == 0
		if(any(qui_0)){
			res$mu = mu
		}
	}


	if(!is.null(bounds)){
		res$bounds = bounds
		res$isBounded = isBounded
	}

	# Fixed-effects
	if(isFixef){

		useExp_fixefCoef = family %in% c("poisson")

		sumFE = attr(mu, "sumFE")
		if(useExp_fixefCoef){
			sumFE = rpar_log(sumFE, env)
		}

		res$sumFE = sumFE

		# The LL and SSR with FE only
		if("ll_fe_only" %in% names(env)){
			res$ll_fe_only = get("ll_fe_only", env)
			res$ssr_fe_only = get("ssr_fe_only", env)
		} else {
			# we need to compute it

			# indicator of whether we compute the exp(mu)
			useExp = family %in% c("poisson", "logit", "negbin")

			# mu, using the offset
			if(!is.null(res$offset)){
				mu_noDum = res$offset
			} else {
				mu_noDum = 0
			}

			if(length(mu_noDum) == 1) mu_noDum = rep(mu_noDum, n)

			exp_mu_noDum = NULL
			if(useExp_fixefCoef){
				exp_mu_noDum = rpar_exp(mu_noDum, env)
			}

			assign("fixef.tol", 1e-4, env) # no need of supa precision
			dummies = getDummies(mu_noDum, exp_mu_noDum, env, coef)

			exp_mu = NULL
			if(useExp_fixefCoef){
				# despite being called mu, it is in fact exp(mu)!!!
				exp_mu = exp_mu_noDum*dummies
				mu = rpar_log(exp_mu, env)
			} else {
				mu = mu_noDum + dummies
				if(useExp){
					exp_mu = rpar_exp(mu, env)
				}
			}

			res$ll_fe_only = famFuns$ll(lhs, mu, exp_mu, env, coef)
			ep = famFuns$expected.predictor(mu, exp_mu, env)
			res$ssr_fe_only = cpp_ssq(lhs - ep)

		}
	}

	if(family == "negbin"){
		theta = coef[".theta"]
		res$theta = theta

		if(notes && theta > 1000){
			message("Very high value of theta (", theta, "). There is no sign of overdisperion, you may consider a Poisson model.")
		}

	}

	class(res) <- "fixest"

	if(verbose > 0){
		cat("\n")
	}

	return(res)
}


####
#### Delayed Warnings ####
####


warn_fixef_iter = function(env){
	# Show warnings related to the nber of times the maximum of iterations was reached

	# For fixed-effect
	fixef.iter = get("fixef.iter", env)
	fixef.iter.limit_reached = get("fixef.iter.limit_reached", env)
	origin = get("origin", env)
	warn = get("warn", env)

	if(!warn) return(invisible(NULL))

	goWarning = FALSE
	warning_msg = ""
	if(fixef.iter.limit_reached > 0){
		goWarning = TRUE
		warning_msg = paste0(origin, ": [Getting the fixed-effects] iteration limit reached (", fixef.iter, ").", ifelse(fixef.iter.limit_reached > 1, paste0(" (", fixef.iter.limit_reached, " times.)"), " (Once.)"))
	}

	# For the fixed-effect derivatives
	deriv.iter = get("deriv.iter", env)
	deriv.iter.limit_reached = get("deriv.iter.limit_reached", env)

	if(deriv.iter.limit_reached > 0){
		prefix = ifelse(goWarning, paste0("\n", sprintf("% *s", nchar(origin) + 2, " ")), paste0(origin, ": "))
		warning_msg = paste0(warning_msg, prefix, "[Getting fixed-effects derivatives] iteration limit reached (", deriv.iter, ").", ifelse(deriv.iter.limit_reached > 1, paste0(" (", deriv.iter.limit_reached, " times.)"), " (Once.)"))
		goWarning = TRUE
	}

	if(goWarning){
		warning(warning_msg, call. = FALSE, immediate. = TRUE)
	}

}


warn_step_halving = function(env){

    nb_sh = get("nb_sh", env)
    warn = get("warn", env)

    if(!warn) return(invisible(NULL))

    if(nb_sh > 0){
        warning("feglm: Step halving due to non-finite deviance (", ifelse(nb_sh > 1, paste0(nb_sh, " times"), "once"), ").", call. = FALSE, immediate. = TRUE)
    }

}


format_error_msg = function(x, origin){
    # Simple formatting of the error msg

    # LATER:
    # - for object not found: provide a better error msg by calling the name of the missing
    #   argument => likely I'll need an match.call argument

    x = gsub("\n+$", "", x)

    if(grepl("^Error (in|:|: in) (fe|fixest)[^\n]+\n", x)){
        res = gsub("^Error (in|:|: in) (fe|fixest)[^\n]+\n *(.+)", "\\3", x)
    } else if(grepl("[Oo]bject '.+' not found", x) || grepl("memory|cannot allocate", x)) {
        res = x
    } else {
       res = paste0(x, "\nThis error was unforeseen by the author of the function ", origin, ". If you think your call to the function is legitimate, could you report?")
    }
    res
}




























