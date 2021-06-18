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
#' @param fml A formula representing the relation to be estimated. For example: \code{fml = z~x+y}. To include fixed-effects, insert them in this formula using a pipe: e.g. \code{fml = z~x+y | fe_1+fe_2}. You can combine two fixed-effects with \code{^}: e.g. \code{fml = z~x+y|fe_1^fe_2}, see details. You can also use variables with varying slopes using square brackets: e.g. in \code{fml = z~y|fe_1[x] + fe_2}, see details. To add IVs, insert the endogenous vars./instruments after a pipe, like in \code{y ~ x | c(x_endo1, x_endo2) ~ x_inst1 + x_inst2}. Note that it should always be the last element, see details. Multiple estimations can be performed at once: for multiple dep. vars, wrap them in \code{c()}: ex \code{c(y1, y2)}. For multiple indep. vars, use the stepwise functions: ex \code{x1 + csw(x2, x3)}. The formula \code{fml = c(y1, y2) ~ x1 + cw0(x2, x3)} leads to 6 estimation, see details.
#' @param weights A formula or a numeric vector. Each observation can be weighted, the weights must be greater than 0. If equal to a formula, it should be one-sided: for example \code{~ var_weight}.
#' @param verbose Integer. Higher values give more information. In particular, it can detail the number of iterations in the demeaning algorithm (the first number is the left-hand-side, the other numbers are the right-hand-side variables).
#' @param demeaned Logical, default is \code{FALSE}. Only used in the presence of fixed-effects: should the centered variables be returned? If \code{TRUE}, it creates the items \code{y_demeaned} and \code{X_demeaned}.
#' @param notes Logical. By default, two notes are displayed: when NAs are removed (to show additional information) and when some observations are removed because of collinearity. To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#' @param collin.tol Numeric scalar, default is \code{1e-10}. Threshold deciding when variables should be considered collinear and subsequently removed from the estimation. Higher values means more variables will be removed (if there is presence of collinearity). One signal of presence of collinearity is t-stats that are extremely low (for instance when t-stats < 1e-3).
#' @param y Numeric vector/matrix/data.frame of the dependent variable(s). Multiple dependent variables will return a \code{fixest_multi} object.
#' @param X Numeric matrix of the regressors.
#' @param fixef_df Matrix/data.frame of the fixed-effects.
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
#' In general, for convergence reasons, it is recommended to always add the fixed-effect and avoid using only the variable with varying slope (i.e. use single square brackets).
#'
#' @section Lagging variables:
#'
#' To use leads/lags of variables in the estimation, you can: i) either provide the argument \code{panel.id}, ii) either set your data set as a panel with the function \code{\link[fixest]{panel}}. Doing either of the two will give you acceess to the lagging functions \code{\link[fixest]{l}},  \code{\link[fixest:l]{f}} and \code{\link[fixest:l]{d}}.
#'
#' You can provide several leads/lags/differences at once: e.g. if your formula is equal to \code{f(y) ~ l(x, -1:1)}, it means that the dependent variable is equal to the lead of \code{y}, and you will have as explanatory variables the lead of \code{x1}, \code{x1} and the lag of \code{x1}. See the examples in function \code{\link[fixest]{l}} for more details.
#'
#' @section Interactions:
#'
#' You can interact a numeric variable with a "factor-like" variable by using \code{i(factor_var, continuous_var, ref)}, where \code{continuous_var} will be interacted with each value of \code{factor_var} and the argument \code{ref} is a value of \code{factor_var} taken as a reference (optional).
#'
#' Using this specific way to create interactions leads to a different display of the interacted values in \code{\link[fixest]{etable}} and offers a special representation of the interacted coefficients in the function \code{\link[fixest]{coefplot}}. See examples.
#'
#'  It is important to note that *if you do not care about the standard-errors of the interactions*, then you can add interactions in the fixed-effects part of the formula, it will be incomparably faster (using the syntax \code{factor_var[continuous_var]}, as explained in the section \dQuote{Varying slopes}).
#'
#' The function \code{\link[fixest:i]{i}} has in fact more arguments, please see details in its associated help page.
#'
#' @section On standard-errors:
#'
#' Standard-errors can be computed in different ways, you can use the arguments \code{se} and \code{dof} in \code{\link[fixest]{summary.fixest}} to define how to compute them. By default, in the presence of fixed-effects, standard-errors are automatically clustered.
#'
#' The following vignette: \href{https://lrberge.github.io/fixest/articles/standard_errors.html}{On standard-errors} describes in details how the standard-errors are computed in \code{fixest} and how you can replicate standard-errors from other software.
#'
#' You can use the functions \code{\link[fixest]{setFixest_se}} and \code{\link[fixest:dof]{setFixest_dof}} to permanently set the way the standard-errors are computed.
#'
#' @section Instrumental variables:
#'
#' To estimate two stage least square regressions, insert the relationship between the endogenous regressor(s) and the instruments in a formula, after a pipe.
#'
#' For example, \code{fml = y ~ x1 | x_endo ~ x_inst} will use the variables \code{x1} and \code{x_inst} in the first stage to explain \code{x_endo}. Then will use the fitted value of \code{x_endo} (which will be named \code{fit_x_endo}) and \code{x1} to explain \code{y}.
#' To include several endogenous regressors, just use "+", like in: \code{fml = y ~ x1 | x_endo1 + x_end2 ~ x_inst1 + x_inst2}.
#'
#' Of course you can still add the fixed-effects, but the IV formula must always come last, like in \code{fml = y ~ x1 | fe1 + fe2 | x_endo ~ x_inst}.
#'
#' If you want to estimate a model without exogenous variables, use \code{"1"} as a placeholder: e.g. \code{fml = y ~ 1 | x_endo + x_inst}.
#'
#' By default, the second stage regression is returned. You can access the first stage(s) regressions either directly in the slot \code{iv_first_stage} (not recommended), or using the argument \code{stage = 1} from the function \code{\link[fixest]{summary.fixest}}. For example \code{summary(iv_est, stage = 1)} will give the first stage(s). Note that using summary you can display both the second and first stages at the same time using, e.g., \code{stage = 1:2} (using \code{2:1} would reverse the order).
#'
#'
#' @section Multiple estimations:
#'
#' Multiple estimations can be performed at once, they just have to be specified in the formula. Multiple estimations yield a \code{fixest_multi} object which is \sQuote{kind of} a list of all the results but includes specific methods to access the results in a handy way. Please have a look at the dedicated vignette: \href{https://lrberge.github.io/fixest/articles/multiple_estimations.html}{Multiple estimations}.
#'
#' To include multiple dependent variables, wrap them in \code{c()} (\code{list()} also works). For instance \code{fml = c(y1, y2) ~ x1} would estimate the model \code{fml = y1 ~ x1} and then the model \code{fml = y2 ~ x1}.
#'
#' To include multiple independent variables, you need to use the stepwise functions. There are 4 stepwise functions: \code{sw}, \code{sw0}, \code{csw}, \code{csw0}. Of course \code{sw} stands for stepwise, and \code{csw} for cumulative stepwise. Let's explain that.
#' Assume you have the following formula: \code{fml = y ~ x1 + sw(x2, x3)}. The stepwise function \code{sw} will estimate the following two models: \code{y ~ x1 + x2} and \code{y ~ x1 + x3}. That is, each element in \code{sw()} is sequentially, and separately, added to the formula. Would have you used \code{sw0} in lieu of \code{sw}, then the model \code{y ~ x1} would also have been estimated. The \code{0} in the name means that the model without any stepwise element also needs to be estimated.
#' Finally, the prefix \code{c} means cumulative: each stepwise element is added to the next. That is, \code{fml = y ~ x1 + csw(x2, x3)} would lead to the following models \code{y ~ x1 + x2} and \code{y ~ x1 + x2 + x3}. The \code{0} has the same meaning and would also lead to the model without the stepwise elements to be estimated: in other words, \code{fml = y ~ x1 + csw0(x2, x3)} leads to the following three models: \code{y ~ x1}, \code{y ~ x1 + x2} and \code{y ~ x1 + x2 + x3}.
#'
#' Multiple independent variables can be combined with multiple dependent variables, as in \code{fml = c(y1, y2) ~ cw(x1, x2, x3)} which would lead to 6 estimations. Multiple estimations can also be combined to split samples (with the arguments \code{split}, \code{fsplit}).
#'
#' You can also add fixed-effects in a stepwise fashion. Note that you cannot perform stepwise estimations on the IV part of the formula (\code{feols} only).
#'
#' A note on performance. The feature of multiple estimations has been highly optimized for \code{feols}, in particular in the presence of fixed-effects. It is faster to estimate multiple models using the formula rather than with a loop. For non-\code{feols} models using the formula is roughly similar to using a loop performance-wise.
#'
#'
#' @return
#' A \code{fixest} object. Note that \code{fixest} objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{resid.fixest}}, etc) or functions (like for instance \code{\link[fixest]{fitstat}} to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then depending on the cases: \code{fixef}: the fixed-effects, \code{iv}: the IV part of the formula.}
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
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
#' \item{collin.var}{(When relevant.) Vector containing the variables removed because of collinearity.}
#' \item{collin.coef}{(When relevant.) Vector of coefficients, where the values of the variables removed because of collinearity are NA.}
#' \item{collin.min_norm}{The minimal diagonal value of the Cholesky decomposition. Small values indicate possible presence collinearity.}
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
#' # Basic estimation
#' #
#'
#' res = feols(Sepal.Length ~ Sepal.Width + Petal.Length, iris)
#' # You can specify clustered standard-errors in summary:
#' summary(res, cluster = ~Species)
#'
#' #
#' # Just one set of fixed-effects:
#' #
#'
#' res = feols(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#' # By default, the SEs are clustered according to the first fixed-effect
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
#' est1 = feols(y ~ l(x1, 0:1), base_did, panel.id = ~id+period)
#' est2 = feols(f(y) ~ l(x1, -1:1), base_did, panel.id = ~id+period)
#' etable(est1, est2, order = "f", drop="Int")
#'
#' #
#' # Using interactions:
#' #
#'
#' data(base_did)
#' # We interact the variable 'period' with the variable 'treat'
#' est_did = feols(y ~ x1 + i(period, treat, 5) | id+period, base_did)
#'
#' # Now we can plot the result of the interaction with coefplot
#' coefplot(est_did)
#' # You have many more example in coefplot help
#'
#' #
#' # Instrumental variables
#' #
#'
#' # To estimate Two stage least squares,
#' # insert a formula describing the endo. vars./instr. relation after a pipe:
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "fe1")
#' base$x_inst1 = 0.2 * base$x1 + 0.7 * base$x2 + rpois(150, 2)
#' base$x_inst2 = 0.2 * base$x2 + 0.7 * base$x3 + rpois(150, 3)
#' base$x_endo1 = 0.5 * base$y + 0.5 * base$x3 + rnorm(150, sd = 2)
#' base$x_endo2 = 1.5 * base$y + 0.5 * base$x3 + 3 * base$x_inst1 + rnorm(150, sd = 5)
#'
#' # Using 2 controls, 1 endogenous var. and 1 instrument
#' res_iv = feols(y ~ x1 + x2 | x_endo1 ~ x_inst1, base)
#'
#' # The second stage is the default
#' summary(res_iv)
#'
#' # To show the first stage:
#' summary(res_iv, stage = 1)
#'
#' # To show both the first and second stages:
#' summary(res_iv, stage = 1:2)
#'
#' # Adding a fixed-effect => IV formula always last!
#' res_iv_fe = feols(y ~ x1 + x2 | fe1 | x_endo1 ~ x_inst1, base)
#'
#' # With two endogenous regressors
#' res_iv2 = feols(y ~ x1 + x2 | x_endo1 + x_endo2 ~ x_inst1 + x_inst2, base)
#'
#' # Now there's two first stages => a fixest_multi object is returned
#' sum_res_iv2 = summary(res_iv2, stage = 1)
#'
#' # You can navigate through it by subsetting:
#' sum_res_iv2[iv = 1]
#'
#' # The stage argument also works in etable:
#' etable(res_iv, res_iv_fe, res_iv2, order = "endo")
#'
#' etable(res_iv, res_iv_fe, res_iv2, stage = 1:2, order = c("endo", "inst"),
#'        group = list(control = "!endo|inst"))
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult = feols(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
#'
#' # We can display the results for the first lhs:
#' etable(est_mult[lhs = 1])
#'
#' # And now the second (access can be made by name)
#' etable(est_mult[lhs = "Solar.R"])
#'
#' # Now we focus on the two last right hand sides
#' # (note that .N can be used to specify the last item)
#' etable(est_mult[rhs = 2:.N])
#'
#' # Combining with split
#' est_split = feols(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'                   airquality, split = ~ Month)
#'
#' # You can display everything at once with the print method
#' est_split
#'
#' # Different way of displaying the results with "compact"
#' summary(est_split, "compact")
#'
#' # You can still select which sample/LHS/RHS to display
#' est_split[sample = 1:2, lhs = 1, rhs = 1]
#'
feols = function(fml, data, weights, offset, subset, split, fsplit, cluster, se, dof, panel.id, fixef, fixef.rm = "none", fixef.tol = 1e-6,
                 fixef.iter = 10000, collin.tol = 1e-10, nthreads = getFixest_nthreads(), lean = FALSE, verbose = 0, warn = TRUE,
                 notes = getFixest_notes(), combine.quick, demeaned = FALSE, mem.clean = FALSE, only.env = FALSE, env, ...){

	dots = list(...)

	# 1st: is the call coming from feglm?
	fromGLM = FALSE
	skip_fixef = FALSE
	if("X" %in% names(dots)){
		fromGLM = TRUE
		# env is provided by feglm
		X = dots$X
		y = as.vector(dots$y)
		init = dots$means
		correct_0w = dots$correct_0w

		if(verbose){
		    time_start = proc.time()
		    gt = function(x, nl = TRUE) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s", ifelse(nl, "\n", ""), sep = "")
		    t0 = proc.time()
		}

	} else {
	    time_start = proc.time()
	    gt = function(x, nl = TRUE) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s", ifelse(nl, "\n", ""), sep = "")
	    t0 = proc.time()

		# we use fixest_env for appropriate controls and data handling

		if(missing(env)){
		    set_defaults("fixest_estimation")
		    call_env = new.env(parent = parent.frame())

		    env = try(fixest_env(fml = fml, data = data, weights = weights, offset = offset, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, panel.id = panel.id, fixef = fixef, fixef.rm = fixef.rm, fixef.tol = fixef.tol, fixef.iter = fixef.iter, collin.tol = collin.tol, nthreads = nthreads, lean = lean, verbose = verbose, warn = warn, notes = notes, combine.quick = combine.quick, demeaned = demeaned, mem.clean = mem.clean, origin = "feols", mc_origin = match.call(), call_env = call_env, ...), silent = TRUE)
		} else if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
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

		# demeaned variables
		if(!is.null(dots$X_demean)){
		    skip_fixef = TRUE
		    X_demean = dots$X_demean
		    y_demean = dots$y_demean
		}

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

		demeaned = get("demeaned", env)

		verbose = get("verbose", env)
		if(verbose >= 2) gt("Setup")
	}

	isFixef = get("isFixef", env)

	# Used to solve with the reduced model
	xwx = dots$xwx
	xwy = dots$xwy

	#
	# Split ####
	#

	do_split = get("do_split", env)
	if(do_split){

	    res = multi_split(env, feols)

	    return(res)
	}

	#
	# Multi fixef ####
	#

	do_multi_fixef = get("do_multi_fixef", env)
	if(do_multi_fixef){

	    res = multi_fixef(env, feols)

	    return(res)
	}

	#
	# Multi LHS and RHS ####
	#


	do_multi_lhs = get("do_multi_lhs", env)
	do_multi_rhs = get("do_multi_rhs", env)
	if(do_multi_lhs || do_multi_rhs){
	    assign("do_multi_lhs", FALSE, env)
	    assign("do_multi_rhs", FALSE, env)
	    do_iv = get("do_iv", env)

	    fml = get("fml", env)
	    lhs_names = get("lhs_names", env)
	    lhs = y

	    if(do_multi_lhs){
	        # We find out which LHS have the same NA patterns => saves a lot of computation

	        n_lhs = length(lhs)
	        lhs_group_is_na = list()
	        lhs_group_id = c()
	        lhs_group_n_na = c()
	        for(i in 1:n_lhs){
	            is_na_current = !is.finite(lhs[[i]])
	            n_na_current = sum(is_na_current)

	            if(i == 1){
	                lhs_group_id = 1
	                lhs_group_is_na[[1]] = is_na_current
	                lhs_group_n_na[1] = n_na_current

	            } else {
	                qui = which(lhs_group_n_na == n_na_current)

	                if(length(qui) > 0){

	                    if(n_na_current == 0){
	                        # no need to check the pattern
	                        lhs_group_id[i] = lhs_group_id[qui[1]]
	                        next
	                    }

	                    for(j in qui){
	                        if(all(is_na_current == lhs_group_is_na[[j]])){
	                            lhs_group_id[i] = lhs_group_id[j]
	                            next
	                        }
	                    }
	                }

	                # if here => new group because couldn't be matched
                    id = max(lhs_group_id) + 1
                    lhs_group_id[i] = id
                    lhs_group_is_na[[id]] = is_na_current
                    lhs_group_n_na[id] = n_na_current
	            }
	        }

	        # we make groups
	        lhs_group = list()
	        for(i in 1:max(lhs_group_id)){
	            lhs_group[[i]] = which(lhs_group_id == i)
	        }

	    } else if(do_multi_lhs == FALSE){
	        lhs_group_is_na = list(FALSE)
	        lhs_group_n_na = 0
	        lhs_group = list(1)
	        lhs = list(lhs) # I really abuse R shallow copy system...
	        names(lhs) = deparse_long(fml[[2]])
	    }

	    if(do_multi_rhs){
	        rhs_info_stepwise = get("rhs_info_stepwise", env)
	        multi_rhs_fml_full = rhs_info_stepwise$fml_all_full
	        multi_rhs_fml_sw = rhs_info_stepwise$fml_all_sw
	        multi_rhs_cumul = rhs_info_stepwise$is_cumul

	        linear_core = get("linear_core", env)
	        rhs = get("rhs_sw", env)

	        # Two schemes:
	        #  - if cumulative: we take advantage of it => both in demeaning and in estimation
	        #  - if regular stepwise => only in demeaning
	        # => of course this is dependent on the pattern of NAs
	        #

	        n_core_left = ifelse(length(linear_core$left) == 1, 0, ncol(linear_core$left))
	        n_core_right = ifelse(length(linear_core$right) == 1, 0, ncol(linear_core$right))

	        # rnc: running number of columns
	        rnc = n_core_left
	        if(rnc == 0){
	            col_start = integer(0)
	        } else {
	            col_start = 1:rnc
	        }

	        rhs_group_is_na = list()
	        rhs_group_id = c()
	        rhs_group_n_na = c()
	        rhs_n_vars = c()
	        rhs_col_id = list()
	        any_na_rhs = FALSE
	        for(i in seq_along(multi_rhs_fml_sw)){
	            # We evaluate the extra data and check the NA pattern

	            my_fml = multi_rhs_fml_sw[[i]]

	            if(i == 1 && (multi_rhs_cumul || identical(my_fml[[3]], 1))){
	                # That case is already in the main linear.mat => no NA
	                rhs_group_id = 1
	                rhs_group_is_na[[1]] = FALSE
	                rhs_group_n_na[1] = 0
	                rhs_n_vars[1] = 0
	                rhs[[1]] = 0
	                if(rnc == 0){
	                    rhs_col_id[[1]] = integer(0)
	                } else {
	                    rhs_col_id[[1]] = 1:rnc
	                }

	                next
	            }

	            rhs_current = rhs[[i]]
	            rhs_n_vars[i] = ncol(rhs_current)
	            info = cpppar_which_na_inf_mat(rhs_current, nthreads)

	            is_na_current = info$is_na_inf
	            if(multi_rhs_cumul && any_na_rhs){
	                # we cumulate the NAs
	                is_na_current = is_na_current | rhs_group_is_na[[rhs_group_id[i - 1]]]
	                info$any_na_inf = any(is_na_current)
	            }

	            n_na_current = 0
	            if(info$any_na_inf){
	                any_na_rhs = TRUE
	                n_na_current = sum(is_na_current)
	            } else {
	                # NULL would lead to problems down the road
	                is_na_current = FALSE
	            }

	            if(i == 1){
	                rhs_group_id = 1
	                rhs_group_is_na[[1]] = is_na_current
	                rhs_group_n_na[1] = n_na_current

	            } else {
	                qui = which(rhs_group_n_na == n_na_current)

	                if(length(qui) > 0){

	                    if(n_na_current == 0){
	                        # no need to check the pattern
	                        rhs_group_id[i] = rhs_group_id[qui[1]]
	                        next
	                    }

	                    go_next = FALSE
	                    for(j in qui){
	                        if(all(is_na_current == rhs_group_is_na[[j]])){
	                            rhs_group_id[i] = rhs_group_id[j]
	                            go_next = TRUE
	                            break
	                        }
	                    }

	                    if(go_next) next
	                }

	                # if here => new group because couldn't be matched
	                id = max(rhs_group_id) + 1
	                rhs_group_id[i] = id
	                rhs_group_is_na[[id]] = is_na_current
	                rhs_group_n_na[id] = n_na_current

	            }
	        }

	        # we make groups
	        rhs_group = list()
	        for(i in 1:max(rhs_group_id)){
	            rhs_group[[i]] = which(rhs_group_id == i)
	        }



	        # Finding the right column IDs to select

	        rhs_group_n_vars = rep(0, length(rhs_group)) # To get the total nber of cols per group

	        for(i in seq_along(multi_rhs_fml_sw)){

    	        if(multi_rhs_cumul){
    	            rnc = rnc + rhs_n_vars[i]
    	            if(rnc == 0){
    	                rhs_col_id[[i]] = integer(0)
    	            } else {
    	                rhs_col_id[[i]] = 1:rnc
    	            }
    	        } else {
    	            id = rhs_group_id[i]
    	            rhs_col_id[[i]] = c(col_start, seq(rnc + rhs_group_n_vars[id] + 1, length.out = rhs_n_vars[i]))
    	            rhs_group_n_vars[id] = rhs_group_n_vars[id] + rhs_n_vars[i]
    	        }
	        }

	        if(n_core_right > 0){
	            # We adjust
	            if(multi_rhs_cumul){
	                for(i in seq_along(multi_rhs_fml_sw)){
	                    id = rhs_group_id[i]
	                    gmax = max(rhs_group[[id]])
	                    rhs_col_id[[i]] = c(rhs_col_id[[i]], n_core_left + sum(rhs_n_vars[1:gmax]) + 1:n_core_right)
	                }
	            } else {
	                for(i in seq_along(multi_rhs_fml_sw)){
	                    id = rhs_group_id[i]
	                    rhs_col_id[[i]] = c(rhs_col_id[[i]], n_core_left + rhs_group_n_vars[id] + 1:n_core_right)
	                }
	            }
	        }

	    } else if(do_multi_rhs == FALSE){

	        multi_rhs_fml_full = list(.xpd(rhs = fml[[3]]))
	        multi_rhs_cumul = FALSE
	        rhs_group_is_na = list(FALSE)
	        rhs_group_n_na = 0
	        rhs_n_vars = 0
	        rhs_group = list(1)
	        rhs = list(0)
	        rhs_col_id = list(1:NCOL(X))
	        linear_core = list(left = X, right = 1)
	    }

	    isLinear_right = length(linear_core$right) > 1

	    isLinear = length(linear_core$left) > 1 || isLinear_right

	    n_lhs = length(lhs)
	    n_rhs = length(rhs)
	    res = vector("list", n_lhs * n_rhs)

	    rhs_names = sapply(multi_rhs_fml_full, function(x) as.character(x)[[2]])

	    for(i in seq_along(lhs_group)){
	        for(j in seq_along(rhs_group)){

	            # NA removal
	            no_na = FALSE
	            if(lhs_group_n_na[i] > 0){
	                if(rhs_group_n_na[j] > 0){
	                    is_na_current = lhs_group_is_na[[i]] | rhs_group_is_na[[j]]
	                } else {
	                    is_na_current = lhs_group_is_na[[i]]
	                }
	            } else if(rhs_group_n_na[j] > 0){
	                is_na_current = rhs_group_is_na[[j]]
	            } else {
	                no_na = TRUE
	            }

	            # Here it depends on whether there are FEs or not, whether it's cumul or not
	            my_lhs = lhs[lhs_group[[i]]]
	            if(isLinear){
	                my_rhs = linear_core[1]
	                if(multi_rhs_cumul){
	                    gmax = max(rhs_group[[j]])
	                    my_rhs[1 + (1:gmax)] = rhs[1:gmax]
	                } else {
	                    for(u in rhs_group[[j]]){
	                        if(length(rhs[[u]]) > 1){
	                            my_rhs[[length(my_rhs) + 1]] = rhs[[u]]
	                        }
	                    }
	                }

	                if(isLinear_right){
	                    my_rhs[[length(my_rhs) + 1]] = linear_core$right
	                }

	            } else{
	                rhs_len = lengths(rhs)
	                if(multi_rhs_cumul){
	                    gmax = max(rhs_group[[j]])
	                    my_rhs = rhs[rhs_len > 1 & seq_along(rhs) <= gmax]
	                } else {
	                    my_rhs = rhs[rhs_len > 1 & seq_along(rhs) %in% rhs_group[[j]]]
	                }

	                if(isLinear_right){
	                    my_rhs[[length(my_rhs) + 1]] = linear_core$right
	                }
	            }

	            len_all = lengths(my_rhs)
	            if(any(len_all == 1)){
	                my_rhs = my_rhs[len_all > 1]
	            }

	            if(!no_na){
	                # NA removal
	                for(u in seq_along(my_lhs)){
	                    my_lhs[[u]] = my_lhs[[u]][!is_na_current]
	                }

	                for(u in seq_along(my_rhs)){
	                    if(length(my_rhs[[u]]) > 1) my_rhs[[u]] = my_rhs[[u]][!is_na_current, , drop = FALSE]
	                }

	                my_env = reshape_env(env, obs2keep = which(!is_na_current), assign_lhs = FALSE, assign_rhs = FALSE)

	            } else {
	                my_env = reshape_env(env)
	            }

	            isLinear_current = TRUE
	            if(length(my_rhs) == 0){
	                X_all = 0
	                isLinear_current = FALSE
	            } else {
	                X_all = do.call("cbind", my_rhs)
	            }

	            if(do_iv){
	                # We need to GET them => they have been modified in my_env
	                iv_lhs = get("iv_lhs", my_env)
	                iv.mat = get("iv.mat", my_env)
	                n_inst = ncol(iv.mat)
	            }

	            if(isFixef){
	                # We batch demean

	                n_vars_X = ifelse(is.null(ncol(X_all)), 0, ncol(X_all))

	                # fixef information
	                fixef_sizes = get("fixef_sizes", my_env)
	                fixef_table_vector = get("fixef_table_vector", my_env)
	                fixef_id_list = get("fixef_id_list", my_env)

	                slope_flag = get("slope_flag", my_env)
	                slope_vars = get("slope_variables", my_env)

	                if(mem.clean) gc()

	                vars_demean = cpp_demean(my_lhs, X_all, weights, iterMax = fixef.iter,
	                                         diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
	                                         fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
	                                         slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
	                                         r_init = init, nthreads = nthreads)

	                X_demean = vars_demean$X_demean
	                y_demean = vars_demean$y_demean

	                if(do_iv){

	                    iv_vars_demean = cpp_demean(iv_lhs, iv.mat, weights, iterMax = fixef.iter,
	                                                diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
	                                                fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
	                                                slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
	                                                r_init = init, nthreads = nthreads)

	                    iv.mat_demean = iv_vars_demean$X_demean
	                    iv_lhs_demean = iv_vars_demean$y_demean
	                }

	            }

	            # We precompute the solution
	            if(do_iv){

	                if(isFixef){
	                    iv_products = cpp_iv_products(X = X_demean, y = y_demean,
	                                                  Z = iv.mat_demean, u = iv_lhs_demean,
	                                                  w = weights, nthreads = nthreads)
	                } else {
	                    iv_products = cpp_iv_products(X = X_all, y = my_lhs, Z = iv.mat,
	                                                  u = iv_lhs, w = weights, nthreads = nthreads)
	                }

	            } else {
	                if(isFixef){
	                    my_products = cpp_sparse_products(X_demean, weights, y_demean, nthreads = nthreads)
	                } else {
	                    my_products = cpp_sparse_products(X_all, weights, my_lhs, nthreads = nthreads)
	                }

	                xwx = my_products$XtX
	                xwy = my_products$Xty
	            }

	            for(ii in seq_along(my_lhs)){
	                i_lhs = lhs_group[[i]][ii]

	                for(jj in rhs_group[[j]]){

	                    qui = rhs_col_id[[jj]]
	                    if(isLinear_current){
	                        my_X = X_all[, qui, drop = FALSE]
	                    } else {
	                        my_X = 0
	                    }

	                    my_fml = .xpd(lhs = lhs_names[i_lhs], rhs = multi_rhs_fml_full[[jj]])
	                    current_env = reshape_env(my_env, lhs = my_lhs[[ii]], rhs = my_X, fml_linear = my_fml)

	                    if(do_iv){
	                        if(isLinear_current){
	                            qui_iv = c(1:n_inst, n_inst + qui)

	                            XtX = iv_products$XtX[qui, qui, drop = FALSE]
	                            Xty = iv_products$Xty[[ii]][qui]
	                        } else {
	                            qui_iv = 1:n_inst

	                            XtX = matrix(0, 1, 1)
	                            Xty = matrix(0, 1, 1)
	                        }

	                        my_iv_products = list(XtX = XtX,
	                                              Xty = Xty,
	                                              ZXtZX = iv_products$ZXtZX[qui_iv, qui_iv, drop = FALSE],
	                                              ZXtu = lapply(iv_products$ZXtu, function(x) x[qui_iv]))

	                        if(isFixef){
	                            my_res = feols(env = current_env, iv_products = my_iv_products,
	                                           X_demean = X_demean[ , qui, drop = FALSE],
	                                           y_demean = y_demean[[ii]],
	                                           iv.mat_demean = iv.mat_demean, iv_lhs_demean = iv_lhs_demean)
	                        } else {
	                            my_res = feols(env = current_env, iv_products = my_iv_products)
	                        }

	                    } else {
	                        if(isFixef){
	                            my_res = feols(env = current_env, xwx = xwx[qui, qui, drop = FALSE], xwy = xwy[[ii]][qui],
	                                           X_demean = X_demean[ , qui, drop = FALSE],
	                                           y_demean = y_demean[[ii]])
	                        } else {
	                            my_res = feols(env = current_env, xwx = xwx[qui, qui, drop = FALSE], xwy = xwy[[ii]][qui])
	                        }
	                    }


	                    res[[index_2D_to_1D(i_lhs, jj, n_rhs)]] = my_res
	                }
	            }

	        }
	    }

	    # Meta information for fixest_multi

	    index = list(lhs = n_lhs, rhs = n_rhs)
	    all_names = list(lhs = lhs_names, rhs = rhs_names)

	    # result
	    res_multi = setup_multi(index, all_names, res)

	    return(res_multi)
	}


	#
	# IV ####
	#

	do_iv = get("do_iv", env)
	if(do_iv){
        assign("do_iv", FALSE, env)
	    assign("verbose", 0, env)

	    # Loaded already
	    # y: lhs
	    # X: linear.mat

	    iv_lhs = get("iv_lhs", env)
	    iv_lhs_names = get("iv_lhs_names", env)
	    iv.mat = get("iv.mat", env) # we enforce (before) at least one variable in iv.mat
	    K = ncol(iv.mat)
	    n_endo = length(iv_lhs)
	    lean = get("lean", env)

	    # Simple check that the function is not misused
	    pblm = intersect(iv_lhs_names, colnames(X))
	    if(length(pblm) > 0){
	        any_exo = length(setdiff(colnames(X), iv_lhs_names)) > 0
	        msg = if(any_exo) "" else " If there is no exogenous variable, just use '1' in the first part of the formula."
	        stop("Endogenous variables should not be used as exogenous regressors. The variable", enumerate_items(pblm, "s.quote.were"), " found in the first part of the multipart formula: ", ifsingle(pblm, "it", "they"), " should not be there.", msg)
	    }

	    if(isFixef){
	        # we batch demean first

	        n_vars_X = ifelse(is.null(ncol(X)), 0, ncol(X))

	        if(mem.clean) gc()

	        if(!is.null(dots$iv_products)){
	            # means this is a call from multiple LHS/RHS

	            X_demean = dots$X_demean
	            y_demean = dots$y_demean

	            iv.mat_demean = dots$iv.mat_demean
	            iv_lhs_demean = dots$iv_lhs_demean

	            iv_products = dots$iv_products

	        } else {
	            # fixef information
	            fixef_sizes = get("fixef_sizes", env)
	            fixef_table_vector = get("fixef_table_vector", env)
	            fixef_id_list = get("fixef_id_list", env)

	            slope_flag = get("slope_flag", env)
	            slope_vars = get("slope_variables", env)

	            vars_demean = cpp_demean(y, X, weights, iterMax = fixef.iter,
	                                     diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
	                                     fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
	                                     slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
	                                     r_init = init, nthreads = nthreads)

	            iv_vars_demean = cpp_demean(iv_lhs, iv.mat, weights, iterMax = fixef.iter,
	                                        diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
	                                        fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
	                                        slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
	                                        r_init = init, nthreads = nthreads)

	            X_demean = vars_demean$X_demean
	            y_demean = vars_demean$y_demean

	            iv.mat_demean = iv_vars_demean$X_demean
	            iv_lhs_demean = iv_vars_demean$y_demean

	            # We precompute the solution
	            iv_products = cpp_iv_products(X = X_demean, y = y_demean, Z = iv.mat_demean,
	                                          u = iv_lhs_demean, w = weights, nthreads = nthreads)

	        }

	        if(n_vars_X == 0){
	            ZX_demean = iv.mat_demean
	            ZX = iv.mat
	        } else {
	            ZX_demean = cbind(iv.mat_demean, X_demean)
	            ZX = cbind(iv.mat, X)
	        }

	        # First stage(s)

	        ZXtZX = iv_products$ZXtZX
	        ZXtu  = iv_products$ZXtu

	        res_first_stage = list()

	        for(i in 1:n_endo){
	            current_env = reshape_env(env, lhs = iv_lhs[[i]], rhs = ZX, fml_iv_endo = iv_lhs_names[i])

	            my_res = feols(env = current_env, xwx = ZXtZX, xwy = ZXtu[[i]],
	                           X_demean = ZX_demean, y_demean = iv_lhs_demean[[i]],
	                           add_fitted_demean = TRUE, iv_call = TRUE, notes = FALSE)

	            # For the F-stats
	            if(n_vars_X == 0){
	                my_res$ssr_no_inst = cpp_ssq(iv_lhs_demean[[i]], weights)
	            } else {
	                fit_no_inst = ols_fit(iv_lhs_demean[[i]], X_demean, w = weights, correct_0w = FALSE,
	                                      collin.tol = collin.tol, nthreads = nthreads,
	                                      xwx = iv_products$XtX, xwy = ZXtu[[i]][-(1:K)])
	                my_res$ssr_no_inst = cpp_ssq(fit_no_inst$residuals, weights)
	            }

	            my_res$iv_stage = 1
	            my_res$iv_inst_names_xpd = colnames(iv.mat)

	            res_first_stage[[iv_lhs_names[i]]] = my_res
	        }

	        if(verbose >= 2) gt("1st stage(s)")

	        # Second stage

	        if(n_endo == 1){
	            res_FS = res_first_stage[[1]]
	            U = as.matrix(res_FS$fitted.values)
	            U_demean = as.matrix(res_FS$fitted.values_demean)
	        } else {
	            U_list = list()
	            U_dm_list = list()
	            for(i in 1:n_endo){
	                res_FS = res_first_stage[[i]]
	                U_list[[i]] = res_FS$fitted.values
	                U_dm_list[[i]] = res_FS$fitted.values_demean
	            }

	            U = do.call("cbind", U_list)
	            U_demean = do.call("cbind", U_dm_list)
	        }

	        colnames(U) = colnames(U_demean) = paste0("fit_", iv_lhs_names)

            if(n_vars_X == 0){
                UX = as.matrix(U)
                UX_demean = as.matrix(U_demean)
            } else {
                UX = cbind(U, X)
                UX_demean = cbind(U_demean, X_demean)
            }

	        XtX = iv_products$XtX
	        Xty = iv_products$Xty
	        iv_prod_second = cpp_iv_product_completion(XtX = XtX, Xty = Xty, X = X_demean,
	                                                   y = y_demean, U = U_demean, w = weights, nthreads = nthreads)

	        UXtUX = iv_prod_second$UXtUX
	        UXty  = iv_prod_second$UXty

	        resid_s1 = lapply(res_first_stage, function(x) x$residuals)

	        current_env = reshape_env(env, rhs = UX)
	        res_second_stage = feols(env = current_env, xwx = UXtUX, xwy = UXty,
	                                 X_demean = UX_demean, y_demean = y_demean,
	                                 resid_1st_stage = resid_s1, iv_call = TRUE, notes = FALSE)

	        # For the F-stats
	        if(n_vars_X == 0){
	            res_second_stage$ssr_no_endo = cpp_ssq(y_demean, weights)
	        } else {
	            fit_no_endo = ols_fit(y_demean, X_demean, w = weights, correct_0w = FALSE,
	                                  collin.tol = collin.tol, nthreads = nthreads,
	                                  xwx = XtX, xwy = Xty)
	            res_second_stage$ssr_no_endo = cpp_ssq(fit_no_endo$residuals, weights)
	        }

	    } else {
	        # fixef == FALSE

	        # We precompute the solution
	        if(!is.null(dots$iv_products)){
	            # means this is a call from multiple LHS/RHS
	            iv_products = dots$iv_products

	        } else {
	            iv_products = cpp_iv_products(X = X, y = y, Z = iv.mat,
	                                          u = iv_lhs, w = weights, nthreads = nthreads)

	        }

	        if(verbose >= 2) gt("IV products")

	        ZX = cbind(iv.mat, X)

	        # First stage(s)

	        ZXtZX = iv_products$ZXtZX
	        ZXtu  = iv_products$ZXtu

	        # Let's put the intercept first => I know it's not really elegant, but that's life
	        is_int = "(Intercept)" %in% colnames(X)
	        if(is_int){
	            nz = ncol(iv.mat)
	            nzx = ncol(ZX)
	            qui = c(nz + 1, (1:nzx)[-(nz + 1)])

	            ZX = ZX[, qui, drop = FALSE]
	            ZXtZX = ZXtZX[qui, qui, drop = FALSE]
	            for(i in seq_along(ZXtu)){
	                ZXtu[[i]] = ZXtu[[i]][qui]
	            }

	        }

	        res_first_stage = list()

	        for(i in 1:n_endo){
	            current_env = reshape_env(env, lhs = iv_lhs[[i]], rhs = ZX, fml_iv_endo = iv_lhs_names[i])
	            my_res = feols(env = current_env, xwx = ZXtZX, xwy = ZXtu[[i]],
	                           iv_call = TRUE, notes = FALSE)

	            # For the F-stats
	            fit_no_inst = ols_fit(iv_lhs[[i]], X, w = weights, correct_0w = FALSE,
	                                  collin.tol = collin.tol, nthreads = nthreads,
	                                  xwx = ZXtZX[-(1:K + is_int), -(1:K + is_int), drop = FALSE],
	                                  xwy = ZXtu[[i]][-(1:K + is_int)])
	            my_res$ssr_no_inst = cpp_ssq(fit_no_inst$residuals, weights)

	            my_res$iv_stage = 1
	            my_res$iv_inst_names_xpd = colnames(iv.mat)

	            res_first_stage[[iv_lhs_names[i]]] = my_res
	        }

	        if(verbose >= 2) gt("1st stage(s)")

	        # Second stage

	        if(n_endo == 1){
	            res_FS = res_first_stage[[1]]
	            U = as.matrix(res_FS$fitted.values)
	        } else {
	            U_list = list()
	            U_dm_list = list()
	            for(i in 1:n_endo){
	                res_FS = res_first_stage[[i]]
	                U_list[[i]] = res_FS$fitted.values
	            }

	            U = do.call("cbind", U_list)
	        }

	        colnames(U) = paste0("fit_", iv_lhs_names)

            UX = cbind(U, X)

	        XtX = iv_products$XtX
	        Xty = iv_products$Xty
	        iv_prod_second = cpp_iv_product_completion(XtX = XtX, Xty = Xty, X = X,
	                                                   y = y, U = U, w = weights, nthreads = nthreads)

	        UXtUX = iv_prod_second$UXtUX
	        UXty  = iv_prod_second$UXty

	        if(is_int){
	            nu = ncol(U)
	            nux = ncol(UX)
	            qui = c(nu + 1, (1:nux)[-(nu + 1)])

	            UX = UX[, qui, drop = FALSE]
	            UXtUX = UXtUX[qui, qui, drop = FALSE]
	            UXty = UXty[qui]
	        }

	        resid_s1 = lapply(res_first_stage, function(x) x$residuals)

	        current_env = reshape_env(env, rhs = UX)
	        res_second_stage = feols(env = current_env, xwx = UXtUX, xwy = UXty,
	                                 resid_1st_stage = resid_s1, iv_call = TRUE, notes = FALSE)

	        # For the F-stats
	        fit_no_endo = ols_fit(y, X, w = weights, correct_0w = FALSE,
	                              collin.tol = collin.tol, nthreads = nthreads,
	                              xwx = XtX, xwy = Xty)
	        res_second_stage$ssr_no_endo = cpp_ssq(fit_no_endo$residuals, weights)
	    }

	    if(verbose >= 2) gt("2nd stage")

	    #
	    # Wu-Hausman endogeneity test
	    #

	    # Current limitation => only standard vcov => later add argument (which would yield the full est.)?
	    # The problem of the full est. is that it takes memory very likely needlessly

	    if(isFixef){
	        ENDO_demean = do.call(cbind, iv_lhs_demean)
	        iv_prod_wh = cpp_iv_product_completion(XtX = UXtUX, Xty = UXty,
	                                               X = UX_demean, y = y_demean, U = ENDO_demean,
	                                               w = weights, nthreads = nthreads)

	        RHS_wh = cbind(ENDO_demean, UX_demean)
	        fit_wh = ols_fit(y_demean, RHS_wh, w = weights, correct_0w = FALSE, collin.tol = collin.tol,
	                         nthreads = nthreads, xwx = iv_prod_wh$UXtUX, xwy = iv_prod_wh$UXty)
	    } else {
	        ENDO = do.call(cbind, iv_lhs)
	        iv_prod_wh = cpp_iv_product_completion(XtX = UXtUX, Xty = UXty,
	                                               X = UX, y = y, U = ENDO,
	                                               w = weights, nthreads = nthreads)

	        RHS_wh = cbind(ENDO, UX)
	        fit_wh = ols_fit(y, RHS_wh, w = weights, correct_0w = FALSE, collin.tol = collin.tol,
	                         nthreads = nthreads, xwx = iv_prod_wh$UXtUX, xwy = iv_prod_wh$UXty)
	    }

	    df1 = n_endo
	    df2 = length(y) - (res_second_stage$nparams + df1)
	    if(any(fit_wh$is_excluded)){
	        stat = p = NA
	    } else {
	        qui = df1 + 1:df1 + ("(Intercept)" %in% names(res_second_stage$coefficients))
	        my_coef = fit_wh$coefficients[qui]
	        vcov_wh = fit_wh$xwx_inv[qui, qui] * cpp_ssq(fit_wh$residuals, weights) / df2
	        stat = drop(my_coef %*% solve(vcov_wh) %*% my_coef) / df1
	        p = pf(stat, df1, df2, lower.tail = FALSE)
	    }

	    res_second_stage$iv_wh = list(stat = stat, p = p, df1 = df1, df2 = df2)

	    #
	    # Sargan
	    #

	    if(n_endo < ncol(iv.mat)){
	        df = ncol(iv.mat) - n_endo
	        resid_2nd = res_second_stage$residuals

	        if(isFixef){
	            xwy = cpppar_xwy(ZX_demean, resid_2nd, weights, nthreads)
	            fit_sargan = ols_fit(resid_2nd, ZX_demean, w = weights, correct_0w = FALSE, collin.tol = collin.tol,
	                                 nthreads = nthreads, xwx = ZXtZX, xwy = xwy)
	        } else {
	            xwy = cpppar_xwy(ZX, resid_2nd, weights, nthreads)
	            fit_sargan = ols_fit(resid_2nd, ZX, w = weights, correct_0w = FALSE, collin.tol = collin.tol,
	                                 nthreads = nthreads, xwx = ZXtZX, xwy = xwy)
	        }

	        r = fit_sargan$residuals
	        stat = length(r) * (1 - cpp_ssq(r, weights) / cpp_ssr_null(resid_2nd))
	        p = pchisq(stat, df, lower.tail = FALSE)
	        res_second_stage$iv_sargan = list(stat = stat, p = p, df = df)
	    }

	    # extra information
	    res_second_stage$iv_inst_names_xpd = res_first_stage[[1]]$iv_inst_names_xpd
	    res_second_stage$iv_endo_names_fit = paste0("fit_", res_second_stage$iv_endo_names)

	    # Collinearity message

	    collin.vars = c(res_second_stage$collin.var, res_first_stage[[1]]$collin.var)
	    res_second_stage$collin.var = unique(collin.vars)
	    if(notes && length(collin.vars) > 0){
	        coll.endo = intersect(collin.vars, res_second_stage$iv_endo_names_fit)
	        coll.inst = intersect(collin.vars, res_second_stage$iv_inst_names_xpd)
	        coll.exo = setdiff(collin.vars, c(coll.endo, coll.inst))

	        n_c = length(collin.vars)
	        n_c_endo = length(coll.endo)
	        n_c_inst = length(coll.inst)
	        n_c_exo = length(coll.exo)

	        msg_endo = msg_exo = msg_inst = NULL
	        if(n_c_endo > 0){
	            msg_endo = paste0("The endogenous regressor", plural(n_c_endo), " ", enumerate_items(coll.endo, "quote", nmax = 3))
	        } else if(n_c_inst > 0){
	            msg_inst = paste0("the instrument", plural(n_c_inst), " ", enumerate_items(coll.inst, "quote", nmax = 3))
	        } else if(n_c_exo > 0){
	            msg_exo = paste0("the exogenous variable", plural(n_c_exo), " ", enumerate_items(coll.exo, "quote", nmax = 3))
	        }

	        msg = enumerate_items(c(msg_endo, msg_inst, msg_exo))
	        msg = gsub("^t", "T", msg)
	        message(msg, " ", plural(n_c, "has"), " been removed because of collinearity (see $collin.var).")

	    }

	    # if lean = TRUE: we clean the IV residuals (which were needed so far)
	    if(lean){
	        for(i in 1:n_endo){
	            res_first_stage[[i]]$residuals = NULL
	            res_first_stage[[i]]$fitted.values = NULL
	            res_first_stage[[i]]$fitted.values_demean = NULL
	        }

	        res_second_stage$residuals = NULL
	        res_second_stage$fitted.values = NULL
	        res_second_stage$fitted.values_demean = NULL
	    }

	    res_second_stage$iv_first_stage = res_first_stage

	    # meta info
	    res_second_stage$iv_stage = 2

	    return(res_second_stage)

	}


	#
	# Regular estimation ####
	#


	onlyFixef = length(X) == 1

	if(fromGLM){
		res = list(coefficients = NA)
	} else {
		res = get("res", env)
	}

	if(skip_fixef){
	    # Variables were already demeaned

	} else if(!isFixef){
		# No Fixed-effects
		y_demean = y
		X_demean = X
		res$means = 0
	} else {
		time_demean = proc.time()

		# Number of nthreads
		n_vars_X = ifelse(is.null(ncol(X)), 0, ncol(X))

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

		vars_demean = cpp_demean(y, X, weights, iterMax = fixef.iter,
		                            diffMax = fixef.tol, r_nb_id_Q = fixef_sizes,
		                            fe_id_list = fixef_id_list, table_id_I = fixef_table_vector,
		                            slope_flag_Q = slope_flag, slope_vars_list = slope_vars,
		                            r_init = init, nthreads = nthreads)

		y_demean = vars_demean$y_demean
		if(onlyFixef){
		    X_demean = matrix(1, nrow = length(y_demean))
		} else {
		    X_demean = vars_demean$X_demean
		}

		res$iterations = vars_demean$iterations
		if(fromGLM){
			res$means = vars_demean$means
		}

		if(mem.clean){
		    rm(vars_demean)
		}

		if(any(abs(slope_flag) > 0) && any(res$iterations > 300)){
		    # Maybe we have a convergence problem
		    # This is poorly coded, but it's a temporary fix
		    opt_fe = check_conv(y_demean, X_demean, fixef_id_list, slope_flag, slope_vars, weights)

		    # This is a bit too rough a check but it should catch the most problematic cases
		    if(any(opt_fe > 1e-4)){
		        msg = "There seems to be a convergence problem due to the presence of variables with varying slopes. The precision of the estimates may not be great."
		        if(any(slope_flag < 0)){
		            sugg = "This convergence problem mostly arises when there are varying slopes without their associated fixed-effect, as is the case in your estimation. Why not try to include the fixed-effect (i.e. use '[' instead of '[[')?"
		        } else {
		            sugg = "As a workaround, and if there are not too many slopes, you can use the variables with varying slopes as regular variables using the function interact (see ?interact)."
		        }
		        msg = paste(msg, sugg)

		        res$convStatus = FALSE

		        res$message = paste0("tol: ", signif_plus(fixef.tol), ", iter: ", max(res$iterations))

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
				gt("Demeaning", FALSE)
			    cat(" (iter: ", paste0(c(tail(res$iterations, 1), res$iterations[-length(res$iterations)]), collapse = ", "), ")\n", sep="")
			} else {
				gt("Demeaning")
			}
		}
	}

	#
	# Estimation
	#

	if(mem.clean){
	    gc()
	}

	if(!onlyFixef){

	    est = ols_fit(y_demean, X_demean, weights, correct_0w, collin.tol, nthreads, xwx, xwy)

	    if(mem.clean){
	        gc()
	    }

	    # Corner case: not any relevant variable
	    if(!is.null(est$all_removed)){
	        all_vars = colnames(X)

	        IN_MULTI = get("IN_MULTI", env)

	        if(isFixef){
	            msg = paste0(ifsingle(all_vars, "The only variable ", "All variables"), enumerate_items(all_vars, "quote.is", nmax = 3), " collinear with the fixed effects. In such circumstances, the estimation is void.")
	        } else {
	            msg =  paste0(ifsingle(all_vars, "The only variable ", "All variables"), enumerate_items(all_vars, "quote.is", nmax = 3), " virtually constant and equal to 0. In such circumstances, the estimation is void.")
	        }

	        if(IN_MULTI || !warn){

	            if(warn) warning(msg)

	            return(fixest_NA_results(env))

	        } else {
	            stop_up(msg, up = fromGLM)
	        }


	    }

		# Formatting the result
	    coef = est$coefficients
	    names(coef) = colnames(X)[!est$is_excluded]
		res$coefficients = coef
		# Additional stuff
		res$residuals = est$residuals
		res$multicol = est$multicol
		res$collin.min_norm = est$collin.min_norm
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
	if(verbose >= 1){
		gt("Estimation")
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

	#
	# IV correction
	#

	if(!is.null(dots$resid_1st_stage)){
	    # We correct the residual
	    is_int = "(Intercept)" %in% names(res$coefficients)
	    resid_new = cpp_iv_resid(res$residuals, res$coefficients, dots$resid_1st_stage, is_int, nthreads)
	    res$iv_residuals = res$residuals
	    res$residuals = resid_new
	}

	#
	# Hessian, score, etc
	#

	if(onlyFixef){
		res$fitted.values = res$sumFE = y - res$residuals
	} else {

	    if(mem.clean){
	        gc()
	    }

		# X_beta / fitted / sumFE
		if(isFixef){
			x_beta = cpppar_xbeta(X, coef, nthreads)
			res$sumFE = y - x_beta - res$residuals
			res$fitted.values = x_beta + res$sumFE
			if(isTRUE(dots$add_fitted_demean)){
			    res$fitted.values_demean = est$fitted.values
			}
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

		if(mem.clean){
		    gc()
		}

		res$sigma2 = cpp_ssq(res$residuals, weights) / (length(y) - df_k)

		res$cov.unscaled = est$xwx_inv * res$sigma2

		rownames(res$cov.unscaled) = colnames(res$cov.unscaled) = names(coef)

		# se
		se = diag(res$cov.unscaled)
		se[se < 0] = NA
		se = sqrt(se)

		# coeftable
		zvalue = coef/se
		pvalue = 2*pt(-abs(zvalue), max(n - df_k, 1))

		coeftable = data.frame("Estimate"=coef, "Std. Error"=se, "t value"=zvalue, "Pr(>|t|)"=pvalue)
		names(coeftable) = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
		row.names(coeftable) = names(coef)

		attr(se, "type") = attr(coeftable, "type") = "Standard"
		res$coeftable = coeftable
		res$se = se
	}

	# fit stats
	if(!cpp_isConstant(res$fitted.values)){
	    res$sq.cor = stats::cor(y, res$fitted.values)**2
	} else {
	    res$sq.cor = NA
	}

	if(mem.clean){
	    gc()
	}

	res$ssr_null = cpp_ssr_null(y, weights)
	res$ssr = cpp_ssq(res$residuals, weights)
	sigma_null = sqrt(res$ssr_null / ifelse(isWeight, sum(weights), n))
	res$ll_null = -1/2/sigma_null^2*res$ssr_null - (log(sigma_null) + log(2*pi)/2) * ifelse(isWeight, sum(weights), n)

	# fixef info
	if(isFixef){
		# For the within R2
		if(!onlyFixef){
			res$ssr_fe_only = cpp_ssq(y_demean, weights)
			sigma = sqrt(res$ssr_fe_only / ifelse(isWeight, sum(weights), n))
			res$ll_fe_only = -1/2/sigma^2*res$ssr_fe_only - (log(sigma) + log(2*pi)/2) * ifelse(isWeight, sum(weights), n)
		}
	}

	if(verbose >= 3) gt("Post-processing")

	class(res) = "fixest"

	do_summary = get("do_summary", env)
	if(do_summary){
	    se = get("se", env)
	    cluster = get("cluster", env)
	    lean = get("lean", env)
	    dof = get("dof", env)
	    agg = get("agg", env)
	    summary_flags = get("summary_flags", env)

	    # If lean = TRUE, 1st stage residuals are still needed for the 2nd stage
	    if(isTRUE(dots$iv_call) && lean){
	        r = res$residuals
	        fv = res$fitted.values
	        fvd = res$fitted.values_demean
	    }

	    res = summary(res, se = se, cluster = cluster, agg = agg, dof = dof, lean = lean, summary_flags = summary_flags)

	    if(isTRUE(dots$iv_call) && lean){
	        res$residuals = r
	        res$fitted.values = fv
	        res$fitted.values_demean = fvd
	    }
	}

	res
}

ols_fit = function(y, X, w, correct_0w = FALSE, collin.tol, nthreads, xwx = NULL, xwy = NULL){
    # No control here -- done before

    if(is.null(xwx)){
        info_products = cpp_sparse_products(X, w, y, correct_0w, nthreads)
        xwx = info_products$XtX
        xwy = info_products$Xty
    }

    multicol = FALSE
    info_inv = cpp_cholesky(xwx, collin.tol, nthreads)

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

    res = list(xwx = xwx, coefficients = beta, fitted.values = fitted.values, xwx_inv = xwx_inv, multicol = multicol, residuals = residuals, is_excluded = is_excluded, collin.min_norm = info_inv$min_norm)

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

    res = list()

    for(k in 1:K){
        if(k == 1){
            x = y
        } else {
            x = X[, k - 1]
        }

        res_tmp = c()
        index_slope = 1
        for(q in 1:Q){
            fixef_id = fixef_id_list[[q]]

            if(slope_flag[q] >= 0){
                res_tmp = c(res_tmp, max(abs(tapply(weights * x, fixef_id, mean))))
            }

            n_slopes = abs(slope_flag[q])
            if(n_slopes > 0){
                for(i in 1:n_slopes){
                    var = slope_vars[[index_slope]]

                    num = tapply(weights * x * var, fixef_id, sum)
                    denom = tapply(weights * var^2, fixef_id, sum)
                    res_tmp = c(res_tmp, max(abs(num/denom)))

                    index_slope = index_slope + 1
                }
            }
        }

        res[[k]] = res_tmp
    }

    res = do.call("rbind", res)

    res
}


#' @rdname feols
feols.fit = function(y, X, fixef_df, offset, split, fsplit, cluster, se, dof, weights, subset,
                     fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000,
                     collin.tol = 1e-10, nthreads = getFixest_nthreads(), lean = FALSE, warn = TRUE,
                     notes = getFixest_notes(), mem.clean = FALSE, verbose = 0, only.env = FALSE, env, ...){


    if(missing(weights)) weights = NULL

    time_start = proc.time()

    if(missing(env)){
        set_defaults("fixest_estimation")
        call_env = new.env(parent = parent.frame())

        env = try(fixest_env(y = y, X = X, fixef_df = fixef_df, offset = offset, weights = weights, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, fixef.rm = fixef.rm, fixef.tol=fixef.tol, fixef.iter=fixef.iter, collin.tol = collin.tol, nthreads = nthreads, lean = lean, warn=warn, notes=notes, verbose = verbose, mem.clean = mem.clean, origin = "feols.fit", mc_origin = match.call(), call_env = call_env, ...), silent = TRUE)

    } else if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)){
        stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
    }

    if("try-error" %in% class(env)){
        mc = match.call()
        origin = ifelse(is.null(mc$origin), "feols.fit", mc$origin)
        stop(format_error_msg(env, origin))
    }

    check_arg(only.env, "logical scalar")
    if(only.env){
        return(env)
    }

    verbose = get("verbose", env)
    if(verbose >= 2) cat("Setup in ", (proc.time() - time_start)[3], "s\n", sep="")

    # workhorse is feols (OK if error msg leads to feols [clear enough])
    res = feols(env = env)

    res

}


#' Fixed-effects GLM estimations
#'
#' Estimates GLM models with any number of fixed-effects.
#'
#' @inheritParams femlm
#' @inheritParams feols
#' @inheritSection feols Combining the fixed-effects
#' @inheritSection feols Varying slopes
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#' @inheritSection feols Multiple estimations
#'
#' @param family Family to be used for the estimation. Defaults to \code{gaussian()}. See \code{\link[stats]{family}} for details of family functions.
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1 (e.g. \code{start = 0}), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients). Default is missing.
#' @param etastart Numeric vector of the same length as the data. Starting values for the linear predictor. Default is missing.
#' @param mustart Numeric vector of the same length as the data. Starting values for the vector of means. Default is missing.
#' @param fixef.tol Precision used to obtain the fixed-effects. Defaults to \code{1e-6}. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations.
#' @param glm.iter Number of iterations of the glm algorithm. Default is 25.
#' @param glm.tol Tolerance level for the glm algorithm. Default is \code{1e-8}.
#' @param verbose Integer. Higher values give more information. In particular, it can detail the number of iterations in the demeaning algoritmh (the first number is the left-hand-side, the other numbers are the right-hand-side variables). It can also detail the step-halving algorithm.
#' @param notes Logical. By default, three notes are displayed: when NAs are removed, when some fixed-effects are removed because of only 0 (or 0/1) outcomes, or when a variable is dropped because of collinearity. To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#'
#' @details
#' The core of the GLM are the weighted OLS estimations. These estimations are performed with \code{\link[fixest]{feols}}. The method used to demean each variable along the fixed-effects is based on Berge (2018), since this is the same problem to solve as for the Gaussian case in a ML setup.
#'
#' @return
#' A \code{fixest} object. Note that \code{fixest} objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{resid.fixest}}, etc) or functions (like for instance \code{\link[fixest]{fitstat}} to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then, if relevant: \code{fixef}: the fixed-effects.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{y}{(When relevant.) The dependent variable (used to compute the within-R2 when fixed-effects are present).}
#' \item{convStatus}{Logical, convergence status of the IRWLS algorithm.}
#' \item{irls_weights}{The weights of the last iteration of the IRWLS algorithm.}
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
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
#' # Poisson estimation
#' res = feglm(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris, "poisson")
#'
#' # You could also use fepois
#' res_pois = fepois(Sepal.Length ~ Sepal.Width + Petal.Length | Species, iris)
#'
#' # With the fit method:
#' res_fit = feglm.fit(iris$Sepal.Length, iris[, 2:3], iris$Species, "poisson")
#'
#' # All results are identical:
#' etable(res, res_pois, res_fit)
#'
#' # Note that you have many more examples in feols
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult = fepois(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
#'
#' # We can display the results for the first lhs:
#' etable(est_mult[lhs = 1])
#'
#' # And now the second (access can be made by name)
#' etable(est_mult[lhs = "Solar.R"])
#'
#' # Now we focus on the two last right hand sides
#' # (note that .N can be used to specify the last item)
#' etable(est_mult[rhs = 2:.N])
#'
#' # Combining with split
#' est_split = fepois(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'                   airquality, split = ~ Month)
#'
#' # You can display everything at once with the print method
#' est_split
#'
#' # Different way of displaying the results with "compact"
#' summary(est_split, "compact")
#'
#' # You can still select which sample/LHS/RHS to display
#' est_split[sample = 1:2, lhs = 1, rhs = 1]
#'
#'
feglm = function(fml, data, family = "gaussian", offset, weights, subset, split, fsplit, cluster, se, dof, panel.id, start = NULL,
                 etastart = NULL, mustart = NULL, fixef, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                 glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(), lean = FALSE,
                 warn = TRUE, notes = getFixest_notes(), verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    if(missing(weights)) weights = NULL

    time_start = proc.time()

    if(missing(env)){
        set_defaults("fixest_estimation")
        call_env = new.env(parent = parent.frame())

        env = try(fixest_env(fml=fml, data=data, family = family, offset = offset, weights = weights, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, panel.id = panel.id, linear.start = start, etastart=etastart, mustart=mustart, fixef = fixef, fixef.rm = fixef.rm, fixef.tol=fixef.tol, fixef.iter=fixef.iter, collin.tol = collin.tol, glm.iter = glm.iter, glm.tol = glm.tol, nthreads = nthreads, lean = lean, warn=warn, notes=notes, verbose = verbose, combine.quick = combine.quick, mem.clean = mem.clean, origin = "feglm", mc_origin = match.call(), call_env = call_env, ...), silent = TRUE)

    } else if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)){
        stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
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
feglm.fit = function(y, X, fixef_df, family = "gaussian", offset, split, fsplit, cluster, se, dof, weights, subset, start = NULL,
                     etastart = NULL, mustart = NULL, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000,
                     collin.tol = 1e-10, glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(), lean = FALSE, warn = TRUE,
                     notes = getFixest_notes(), mem.clean = FALSE, verbose = 0, only.env = FALSE, env, ...){

    dots = list(...)

    lean_internal = isTRUE(dots$lean_internal)
    means = 1
    if(!missing(env)){
        # This is an internal call from the function feglm
        # no need to further check the arguments
        # we extract them from the env

        if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
            stop("Argument 'env' must be an environment created by a fixest estimation. Currently it is not ", ifelse(r, "an", "a 'fixest'"), " environment.")
        }

        # main variables
        if(missing(y)) y = get("lhs", env)
        if(missing(X)) X = get("linear.mat", env)
        if(!missing(fixef_df) && is.null(fixef_df)){
            assign("isFixef", FALSE, env)
        }

        if(missing(offset)) offset = get("offset.value", env)
        if(missing(weights)) weights = get("weights.value", env)

        # other params
        if(missing(fixef.tol)) fixef.tol = get("fixef.tol", env)
        if(missing(fixef.iter)) fixef.iter = get("fixef.iter", env)
        if(missing(collin.tol)) collin.tol = get("collin.tol", env)
        if(missing(glm.iter)) glm.iter = get("glm.iter", env)
        if(missing(glm.tol)) glm.tol = get("glm.tol", env)
        if(missing(warn)) warn = get("warn", env)
        if(missing(verbose)) verbose = get("verbose", env)

        # starting point of the fixed-effects
        if(!is.null(dots$means)) means = dots$means

        # init
        init.type = get("init.type", env)
        starting_values = get("starting_values", env)
        if(lean_internal){
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

        set_defaults("fixest_estimation")
        call_env = new.env(parent = parent.frame())

        env = try(fixest_env(y = y, X = X, fixef_df = fixef_df, family = family, nthreads = nthreads, lean = lean, offset = offset, weights = weights, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, linear.start = start, etastart=etastart, mustart=mustart, fixef.rm = fixef.rm, fixef.tol = fixef.tol, fixef.iter = fixef.iter, collin.tol = collin.tol, glm.iter = glm.iter, glm.tol = glm.tol, notes=notes, mem.clean = mem.clean, warn=warn, verbose = verbose, origin = "feglm.fit", mc_origin = match.call(), call_env = call_env, ...), silent = TRUE)

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


    #
    # Split ####
    #

    do_split = get("do_split", env)
    if(do_split){

        res = multi_split(env, feglm.fit)

        return(res)
    }

    #
    # Multi fixef ####
    #

    do_multi_fixef = get("do_multi_fixef", env)
    if(do_multi_fixef){

        res = multi_fixef(env, feglm.fit)

        return(res)
    }

    #
    # Multi LHS and RHS ####
    #

    do_multi_lhs = get("do_multi_lhs", env)
    do_multi_rhs = get("do_multi_rhs", env)
    if(do_multi_lhs || do_multi_rhs){

        res = multi_LHS_RHS(env, feglm.fit)

        return(res)
    }

    #
    # Regular estimation ####
    #

    # Setup:
    family = get("family_funs", env)
    isFixef = get("isFixef", env)
    nthreads = get("nthreads", env)
    isWeight = length(weights) > 1
    isOffset = length(offset) > 1
    nobs = length(y)
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

    if(mem.clean){
        gc()
    }

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
            model_fe = try(feglm.fit(X = 0, etastart = eta, offset = offset_fe, glm.tol = 1e-2, fixef.tol = 1e-2, env = env, lean_internal = TRUE))

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

    if(mem.clean){
        gc()
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

    wols_means = 1
    conv = FALSE
    warning_msg = div_message = ""
    for (iter in 1:glm.iter) {

        if(mem.clean){
            gc()
        }

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

        if(mem.clean && iter > 1){
            rm(wols)
            gc()
        }

        wols = feols(y = z, X = X, weights = w, means = wols_means, correct_0w = any_0w, env = env, fixef.tol = fixef.tol * 10**(iter==1), fixef.iter = fixef.iter, collin.tol = collin.tol, nthreads = nthreads, mem.clean = mem.clean, verbose = verbose - 1)

        if(isTRUE(wols$NA_model)){
            return(wols)
        }

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
        } else {
            wols_means = wols$means
        }

        eta = wols$fitted.values
        if(isOffset){
            eta = eta + offset
        }

        if(mem.clean){
            gc()
        }

        mu = linkinv(eta)

        dev = sum_dev.resids(y, mu, eta, wt = weights)
        dev_evol = dev - devold

        if(verbose >= 1) cat("Iteration: ", sprintf("%02i", iter), " -- Deviance = ", numberFormatNormal(dev), " -- Evol. = ", dev_evol, "\n", sep = "")

        #
        # STEP HALVING
        #

        no_SH = is.finite(dev) && (abs(dev_evol) < glm.tol || abs(dev_evol)/(0.1 + abs(dev)) < glm.tol)
        if(no_SH == FALSE && (!is.finite(dev) || dev_evol > 0 || !valideta(eta) || !validmu(mu))){

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

                if(mem.clean){
                    gc()
                }

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

    res$nparams = res$nparams - collin.adj

    res$irls_weights = w # weights from the iteratively reweighted least square

    res$coefficients = coef = wols$coefficients
    res$collin.min_norm = wols$collin.min_norm

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
        res$dispersion = sum(weighted_resids * wols$residuals) / max(res$nobs - res$nparams, 1)
    }

    res$working_residuals = wols$residuals

    if(!onlyFixef && !lean_internal){
        # score + hessian + vcov

        if(mem.clean){
            gc()
        }

        # dispersion + scores
        if(family$family %in% c("poisson", "binomial")){
            res$scores = (wols$residuals * res$irls_weights) * wols$X_demean
            res$hessian = cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads)
        } else {
            res$scores = (weighted_resids / res$dispersion) * wols$X_demean
            res$hessian = cpppar_crossprod(wols$X_demean, res$irls_weights, nthreads) / res$dispersion
        }

        info_inv = cpp_cholesky(res$hessian, collin.tol, nthreads)
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
        zvalue = coef/se
        use_t = !family$family %in% c("poisson", "binomial")
        if(use_t){
            pvalue = 2*pt(-abs(zvalue), max(res$nobs - res$nparams, 1))
            ctable_names = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
        } else {
            pvalue = 2*pnorm(-abs(zvalue))
            ctable_names = c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
        }

        coeftable = data.frame("Estimate"=coef, "Std. Error"=se, "z value"=zvalue, "Pr(>|z|)"=pvalue)
        names(coeftable) = ctable_names
        row.names(coeftable) = names(coef)

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

            if(mem.clean){
                gc()
            }

            res$loglik = sum( (y * eta - mu - cpppar_lgamma(y + 1, nthreads)) * weights)
        } else {
            # lfact is later used in model0 and is costly to compute
            lfact = sum(rpar_lgamma(y + 1, env))
            assign("lfactorial", lfact, env)
            res$loglik = sum(y * eta - mu) - lfact
        }
    } else {
        res$loglik = family$aic(y = y, n = rep.int(1, n), mu = res$fitted.values, wt = weights, dev = dev) / -2
    }

    if(lean_internal){
        return(res)
    }

    # The pseudo_r2
    if(family_equiv %in% c("poisson", "logit")){
        model0 = get_model_null(env, theta.init = NULL)
        ll_null = model0$loglik
        fitted_null = linkinv(model0$constant)

    } else {
        if(verbose >= 1) cat("Null model:\n")

        if(mem.clean){
            gc()
        }

        model_null = feglm.fit(X = matrix(1, nrow = n, ncol = 1), fixef_df = NULL, env = env, lean_internal = TRUE)
        ll_null = model_null$loglik
        fitted_null = model_null$fitted.values
    }
    res$ll_null = ll_null
    res$pseudo_r2 = 1 - (res$loglik - df_k)/(ll_null - 1)

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

    do_summary = get("do_summary", env)
    if(do_summary){
        se = get("se", env)
        cluster = get("cluster", env)
        lean = get("lean", env)
        dof = get("dof", env)
        agg = get("agg", env)
        summary_flags = get("summary_flags", env)

        # To compute the RMSE and lean = TRUE
        if(lean) res$ssr = cpp_ssq(res$residuals, weights)

        res = summary(res, se = se, cluster = cluster, agg = agg, dof = dof, lean = lean, summary_flags = summary_flags)
    }

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
#' @inheritSection feols Multiple estimations
#'
#' @param fml A formula representing the relation to be estimated. For example: \code{fml = z~x+y}. To include fixed-effects, insert them in this formula using a pipe: e.g. \code{fml = z~x+y|fixef_1+fixef_2}. Multiple estimations can be performed at once: for multiple dep. vars, wrap them in \code{c()}: ex \code{c(y1, y2)}. For multiple indep. vars, use the stepwise functions: ex \code{x1 + csw(x2, x3)}. The formula \code{fml = c(y1, y2) ~ x1 + cw0(x2, x3)} leads to 6 estimation, see details.
#' @param start Starting values for the coefficients. Can be: i) a numeric of length 1 (e.g. \code{start = 0}, the default), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients).
#'
#' @details
#' Note that the functions \code{\link[fixest]{feglm}} and \code{\link[fixest]{femlm}} provide the same results when using the same families but differ in that the latter is a direct maximum likelihood optimization (so the two can really have different convergence rates).
#'
#' @return
#' A \code{fixest} object. Note that \code{fixest} objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{resid.fixest}}, etc) or functions (like for instance \code{\link[fixest]{fitstat}} to access any fit statistic).
#' \item{nobs}{The number of observations.}
#' \item{fml}{The linear formula of the call.}
#' \item{call}{The call of the function.}
#' \item{method}{The method used to estimate the model.}
#' \item{family}{The family used to estimate the model.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then, if relevant: \code{fixef}: the fixed-effects; \code{NL}: the non linear part of the formula.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{fixef_vars}{The names of each fixed-effect dimension.}
#' \item{fixef_id}{The list (of length the number of fixed-effects) of the fixed-effects identifiers for each observation.}
#' \item{fixef_sizes}{The size of each fixed-effect (i.e. the number of unique identifierfor each fixed-effect dimension).}
#' \item{convStatus}{Logical, convergence status.}
#' \item{message}{The convergence message from the optimization procedures.}
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
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
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 fixed-effects
#' # 1) Poisson estimation
#' est_pois = femlm(Euros ~ log(dist_km) | Origin + Destination + Product, trade)
#'
#' # 2) Log-Log Gaussian estimation (with same FEs)
#' est_gaus = update(est_pois, log(Euros+1) ~ ., family = "gaussian")
#'
#' # Comparison of the results using the function etable
#' etable(est_pois, est_gaus)
#' # Now using two way clustered standard-errors
#' etable(est_pois, est_gaus, se = "twoway")
#'
#' # Comparing different types of standard errors
#' sum_hetero   = summary(est_pois, se = "hetero")
#' sum_oneway   = summary(est_pois, se = "cluster")
#' sum_twoway   = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' etable(sum_hetero, sum_oneway, sum_twoway, sum_threeway)
#'
#'
#' #
#' # Multiple estimations:
#' #
#'
#' # 6 estimations
#' est_mult = femlm(c(Ozone, Solar.R) ~ Wind + Temp + csw0(Wind:Temp, Day), airquality)
#'
#' # We can display the results for the first lhs:
#' etable(est_mult[lhs = 1])
#'
#' # And now the second (access can be made by name)
#' etable(est_mult[lhs = "Solar.R"])
#'
#' # Now we focus on the two last right hand sides
#' # (note that .N can be used to specify the last item)
#' etable(est_mult[rhs = 2:.N])
#'
#' # Combining with split
#' est_split = fepois(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'                   airquality, split = ~ Month)
#'
#' # You can display everything at once with the print method
#' est_split
#'
#' # Different way of displaying the results with "compact"
#' summary(est_split, "compact")
#'
#' # You can still select which sample/LHS/RHS to display
#' est_split[sample = 1:2, lhs = 1, rhs = 1]
#'
#'
#'
#'
femlm = function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), start = 0, fixef, fixef.rm = "perfect",
						offset, subset, split, fsplit, cluster, se, dof, panel.id, fixef.tol = 1e-5, fixef.iter = 10000,
						nthreads = getFixest_nthreads(), lean = FALSE, verbose = 0, warn = TRUE,
						notes = getFixest_notes(), theta.init, combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

	# This is just an alias

    call_env_bis = new.env(parent = parent.frame())

	res = try(feNmlm(fml = fml, data = data, family = family, fixef = fixef, fixef.rm = fixef.rm, offset = offset, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, panel.id = panel.id, start = start, fixef.tol=fixef.tol, fixef.iter=fixef.iter, nthreads=nthreads, lean = lean, verbose=verbose, warn=warn, notes=notes, theta.init = theta.init, combine.quick = combine.quick, mem.clean = mem.clean, origin = "femlm", mc_origin_bis = match.call(), call_env_bis = call_env_bis, only.env = only.env, env = env, ...), silent = TRUE)

	if("try-error" %in% class(res)){
		stop(format_error_msg(res, "femlm"))
	}

	return(res)
}

#' @rdname femlm
fenegbin = function(fml, data, theta.init, start = 0, fixef, fixef.rm = "perfect", offset, subset, split, fsplit, cluster, se, dof, panel.id,
                    fixef.tol = 1e-5, fixef.iter = 10000, nthreads = getFixest_nthreads(), lean = FALSE,
                    verbose = 0, warn = TRUE, notes = getFixest_notes(), combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    # We control for the problematic argument family
    if("family" %in% names(match.call())){
        stop("Function fenegbin does not accept the argument 'family'.")
    }

    # This is just an alias
    call_env_bis = new.env(parent = parent.frame())

    res = try(feNmlm(fml = fml, data=data, family = "negbin", theta.init = theta.init, start = start, fixef = fixef, fixef.rm = fixef.rm, offset = offset, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, panel.id = panel.id, fixef.tol = fixef.tol, fixef.iter = fixef.iter, nthreads = nthreads, lean = lean, verbose = verbose, warn = warn, notes = notes, combine.quick = combine.quick, mem.clean = mem.clean, origin = "fenegbin", mc_origin_bis = match.call(), call_env_bis = call_env_bis, only.env = only.env, env = env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "fenegbin"))
    }

    return(res)
}

#' @rdname feglm
fepois = function(fml, data, offset, weights, subset, split, fsplit, cluster, se, dof, panel.id,
                  start = NULL, etastart = NULL, mustart = NULL,
                  fixef, fixef.rm = "perfect", fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                  glm.iter = 25, glm.tol = 1e-8, nthreads = getFixest_nthreads(), lean = FALSE, warn = TRUE, notes = getFixest_notes(),
                  verbose = 0, combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

    # We control for the problematic argument family
    if("family" %in% names(match.call())){
        stop("Function fepois does not accept the argument 'family'.")
    }

    # This is just an alias
    call_env_bis = new.env(parent = parent.frame())

    res = try(feglm(fml = fml, data = data, family = "poisson", offset = offset, weights = weights, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, panel.id = panel.id, start = start, etastart = etastart, mustart = mustart, fixef = fixef, fixef.rm = fixef.rm, fixef.tol = fixef.tol, fixef.iter = fixef.iter, collin.tol = collin.tol, glm.iter = glm.iter, glm.tol = glm.tol, nthreads = nthreads, lean = lean, warn = warn, notes = notes, verbose = verbose, combine.quick = combine.quick, mem.clean = mem.clean, origin_bis = "fepois", mc_origin_bis = match.call(), call_env_bis = call_env_bis, only.env=only.env, env=env, ...), silent = TRUE)

    if("try-error" %in% class(res)){
        stop(format_error_msg(res, "fepois"))
    }

    return(res)
}



#' Fixed effects nonlinear maximum likelihood models
#'
#' This function estimates maximum likelihood models (e.g., Poisson or Logit) with non-linear in parameters right-hand-sides and is efficient to handle any number of fixed effects. If you do not use non-linear in parameters right-hand-side, use \code{\link[fixest]{femlm}} or \code{\link[fixest]{feglm}} instead (their design is simpler).
#'
#' @inheritParams summary.fixest
#' @inheritParams panel
#' @inheritSection feols Lagging variables
#' @inheritSection feols Interactions
#' @inheritSection feols On standard-errors
#' @inheritSection feols Multiple estimations
#'
#' @param fml A formula. This formula gives the linear formula to be estimated (it is similar to a \code{lm} formula), for example: \code{fml = z~x+y}. To include fixed-effects variables, insert them in this formula using a pipe (e.g. \code{fml = z~x+y|fixef_1+fixef_2}). To include a non-linear in parameters element, you must use the argment \code{NL.fml}. Multiple estimations can be performed at once: for multiple dep. vars, wrap them in \code{c()}: ex \code{c(y1, y2)}. For multiple indep. vars, use the stepwise functions: ex \code{x1 + csw(x2, x3)}. This leads to 6 estimation \code{fml = c(y1, y2) ~ x1 + cw0(x2, x3)}. See details.
#' @param start Starting values for the coefficients in the linear part (for the non-linear part, use NL.start). Can be: i) a numeric of length 1 (e.g. \code{start = 0}, the default), ii) a numeric vector of the exact same length as the number of variables, or iii) a named vector of any length (the names will be used to initialize the appropriate coefficients).
#' @param NL.fml A formula. If provided, this formula represents the non-linear part of the right hand side (RHS). Note that contrary to the \code{fml} argument, the coefficients must explicitly appear in this formula. For instance, it can be \code{~a*log(b*x + c*x^3)}, where \code{a}, \code{b}, and \code{c} are the coefficients to be estimated. Note that only the RHS of the formula is to be provided, and NOT the left hand side.
#' @param split A one sided formula representing a variable (eg \code{split = ~var}) or a vector. If provided, the sample is split according to the variable and one estimation is performed for each value of that variable. If you also want to include the estimation for the full sample, use the argument \code{fsplit} instead.
#' @param fsplit A one sided formula representing a variable (eg \code{split = ~var}) or a vector. If provided, the sample is split according to the variable and one estimation is performed for each value of that variable. This argument is the same as split but also includes the full sample as the first estimation.
#' @param data A data.frame containing the necessary variables to run the model. The variables of the non-linear right hand side of the formula are identified with this \code{data.frame} names. Can also be a matrix.
#' @param family Character scalar. It should provide the family. The possible values are "poisson" (Poisson model with log-link, the default), "negbin" (Negative Binomial model with log-link), "logit" (LOGIT model with log-link), "gaussian" (Gaussian model).
#' @param fixef Character vector. The names of variables to be used as fixed-effects. These variables should contain the identifier of each observation (e.g., think of it as a panel identifier). Note that the recommended way to include fixed-effects is to insert them directly in the formula.
#' @param subset A vector (logical or numeric) or a one-sided formula. If provided, then the estimation will be performed only on the observations defined by this argument.
#' @param NL.start (For NL models only) A list of starting values for the non-linear parameters. ALL the parameters are to be named and given a staring value. Example: \code{NL.start=list(a=1,b=5,c=0)}. Though, there is an exception: if all parameters are to be given the same starting value, you can use a numeric scalar.
#' @param lower (For NL models only) A list. The lower bound for each of the non-linear parameters that requires one. Example: \code{lower=list(b=0,c=0)}. Beware, if the estimated parameter is at his lower bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param upper (For NL models only) A list. The upper bound for each of the non-linear parameters that requires one. Example: \code{upper=list(a=10,c=50)}. Beware, if the estimated parameter is at his upper bound, then asymptotic theory cannot be applied and the standard-error of the parameter cannot be estimated because the gradient will not be null. In other words, when at its upper/lower bound, the parameter is considered as 'fixed'.
#' @param NL.start.init (For NL models only) Numeric scalar. If the argument \code{NL.start} is not provided, or only partially filled (i.e. there remain non-linear parameters with no starting value), then the starting value of all remaining non-linear parameters is set to \code{NL.start.init}.
#' @param offset A formula or a numeric vector. An offset can be added to the estimation. If equal to a formula, it should be of the form (for example) \code{~0.5*x**2}. This offset is linearly added to the elements of the main formula 'fml'.
#' @param jacobian.method (For NL models only) Character scalar. Provides the method used to numerically compute the Jacobian of the non-linear part. Can be either \code{"simple"} or \code{"Richardson"}. Default is \code{"simple"}. See the help of \code{\link[numDeriv]{jacobian}} for more information.
#' @param useHessian Logical. Should the Hessian be computed in the optimization stage? Default is \code{TRUE}.
#' @param hessian.args List of arguments to be passed to function \code{\link[numDeriv]{genD}}. Defaults is missing. Only used with the presence of \code{NL.fml}.
#' @param opt.control List of elements to be passed to the optimization method \code{\link[stats]{nlminb}}. See the help page of \code{\link[stats]{nlminb}} for more information.
#' @param nthreads The number of threads. Can be: a) an integer lower than, or equal to, the maximum number of threads; b) 0: meaning all available threads will be used; c) a number strictly between 0 and 1 which represents the fraction of all threads to use. The default is to use 50\% of all threads. You can set permanently the number of threads used within this package using the function \code{\link[fixest]{setFixest_nthreads}}.
#' @param verbose Integer, default is 0. It represents the level of information that should be reported during the optimisation process. If \code{verbose=0}: nothing is reported. If \code{verbose=1}: the value of the coefficients and the likelihood are reported. If \code{verbose=2}: \code{1} + information on the computing time of the null model, the fixed-effects coefficients and the hessian are reported.
#' @param theta.init Positive numeric scalar. The starting value of the dispersion parameter if \code{family="negbin"}. By default, the algorithm uses as a starting value the theta obtained from the model with only the intercept.
#' @param fixef.rm Can be equal to "perfect" (default), "singleton", "both" or "none". Controls which observations are to be removed. If "perfect", then observations having a fixed-effect with perfect fit (e.g. only 0 outcomes in Poisson estimations) will be removed. If "singleton", all observations for which a fixed-effect appears only once will be removed. The meaning of "both" and "none" is direct.
#' @param fixef.tol Precision used to obtain the fixed-effects. Defaults to \code{1e-5}. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations. Argument \code{fixef.tol} cannot be lower than \code{10000*.Machine$double.eps}. Note that this parameter is dynamically controlled by the algorithm.
#' @param fixef.iter Maximum number of iterations in fixed-effects algorithm (only in use for 2+ fixed-effects). Default is 10000.
#' @param deriv.iter Maximum number of iterations in the algorithm to obtain the derivative of the fixed-effects (only in use for 2+ fixed-effects). Default is 1000.
#' @param deriv.tol Precision used to obtain the fixed-effects derivatives. Defaults to \code{1e-4}. It corresponds to the maximum absolute difference allowed between two coefficients of successive iterations. Argument \code{deriv.tol} cannot be lower than \code{10000*.Machine$double.eps}.
#' @param warn Logical, default is \code{TRUE}. Whether warnings should be displayed (concerns warnings relating to convergence state).
#' @param notes Logical. By default, two notes are displayed: when NAs are removed (to show additional information) and when some observations are removed because of only 0 (or 0/1) outcomes in a fixed-effect setup (in Poisson/Neg. Bin./Logit models). To avoid displaying these messages, you can set \code{notes = FALSE}. You can remove these messages permanently by using \code{setFixest_notes(FALSE)}.
#' @param combine.quick Logical. When you combine different variables to transform them into a single fixed-effects you can do e.g. \code{y ~ x | paste(var1, var2)}. The algorithm provides a shorthand to do the same operation: \code{y ~ x | var1^var2}. Because pasting variables is a costly operation, the internal algorithm may use a numerical trick to hasten the process. The cost of doing so is that you lose the labels. If you are interested in getting the value of the fixed-effects coefficients after the estimation, you should use \code{combine.quick = FALSE}. By default it is equal to \code{FALSE} if the number of observations is lower than 50,000, and to \code{TRUE} otherwise.
#' @param only.env (Advanced users.) Logical, default is \code{FALSE}. If \code{TRUE}, then only the environment used to make the estimation is returned.
#' @param mem.clean Logical, default is \code{FALSE}. Only to be used if the data set is large compared to the available RAM. If \code{TRUE} then intermediary objects are removed as much as possible and \code{\link[base]{gc}} is run before each substantial C++ section in the internal code to avoid memory issues.
#' @param lean Logical, default is \code{FALSE}. If \code{TRUE} then all large objects are removed from the returned result: this will save memory but will block the possibility to use many methods. It is recommended to use the arguments \code{se} or \code{cluster} to obtain the appropriate standard-errors at estimation time, since obtaining different SEs won't be possible afterwards.
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
#' A \code{fixest} object. Note that \code{fixest} objects contain many elements and most of them are for internal use, they are presented here only for information. To access them, it is safer to use the user-level methods (e.g. \code{\link[fixest]{vcov.fixest}}, \code{\link[fixest]{resid.fixest}}, etc) or functions (like for instance \code{\link[fixest]{fitstat}} to access any fit statistic).
#' \item{coefficients}{The named vector of coefficients.}
#' \item{coeftable}{The table of the coefficients with their standard errors, z-values and p-values.}
#' \item{loglik}{The loglikelihood.}
#' \item{iterations}{Number of iterations of the algorithm.}
#' \item{nobs}{The number of observations.}
#' \item{nparams}{The number of parameters of the model.}
#' \item{call}{The call.}
#' \item{fml}{The linear formula of the call.}
#' \item{fml_all}{A list containing different parts of the formula. Always contain the linear formula. Then, if relevant: \code{fixef}: the fixed-effects; \code{NL}: the non linear part of the formula.}
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
#' \item{obs_selection}{(When relevant.) List containing vectors of integers. It represents the sequential selection of observation vis a vis the original data set.}
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
#' etable(est1_L, est1_NL)
#'
#' # Now generating a non-linear relation (E(z2) = x + y + 1):
#' z2 = rpois(n, x + y) + rpois(n, 1)
#' base$z2 = z2
#'
#' # Estimation using this non-linear form
#' est2_NL = feNmlm(z2 ~ 0, base, NL.fml = ~log(a*x + b*y),
#'                NL.start = 2, lower = list(a=0, b=0))
#' # we can't estimate this relation linearily
#' # => closest we can do:
#' est2_L = femlm(z2 ~ log(x) + log(y), base)
#'
#' # Difference between the two models:
#' etable(est2_L, est2_NL)
#'
#' # Plotting the fits:
#' plot(x, z2, pch = 18)
#' points(x, fitted(est2_L), col = 2, pch = 1)
#' points(x, fitted(est2_NL), col = 4, pch = 2)
#'
#'
feNmlm = function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml, fixef, fixef.rm = "perfect", NL.start, lower, upper, NL.start.init, offset, subset, split, fsplit, cluster, se, dof, panel.id, start = 0, jacobian.method="simple", useHessian = TRUE, hessian.args = NULL, opt.control = list(), nthreads = getFixest_nthreads(), lean = FALSE, verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000, deriv.tol = 1e-4, deriv.iter = 1000, warn = TRUE, notes = getFixest_notes(), combine.quick, mem.clean = FALSE, only.env = FALSE, env, ...){

	time_start = proc.time()

	if(missing(env)){

	    set_defaults("fixest_estimation")
	    call_env = new.env(parent = parent.frame())

	    env = try(fixest_env(fml = fml, data = data, family = family, NL.fml = NL.fml, fixef = fixef, fixef.rm = fixef.rm, NL.start = NL.start, lower = lower, upper = upper, NL.start.init = NL.start.init, offset = offset, subset = subset, split = split, fsplit = fsplit, cluster = cluster, se = se, dof = dof, panel.id = panel.id, linear.start = start, jacobian.method = jacobian.method, useHessian = useHessian, opt.control = opt.control, nthreads = nthreads, lean = lean, verbose = verbose, theta.init = theta.init, fixef.tol = fixef.tol, fixef.iter = fixef.iter, deriv.iter = deriv.iter, warn = warn, notes = notes, combine.quick = combine.quick, mem.clean = mem.clean, mc_origin = match.call(), call_env = call_env, computeModel0 = TRUE, ...), silent = TRUE)

	} else if((r <- !is.environment(env)) || !isTRUE(env$fixest_env)) {
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


	#
	# Split ####
	#

	do_split = get("do_split", env)
	if(do_split){

	    res = multi_split(env, feNmlm)

	    return(res)
	}

	#
	# Multi fixef ####
	#

	do_multi_fixef = get("do_multi_fixef", env)
	if(do_multi_fixef){

	    res = multi_fixef(env, feNmlm)

	    return(res)
	}

	#
	# Multi LHS and RHS ####
	#

	do_multi_lhs = get("do_multi_lhs", env)
	do_multi_rhs = get("do_multi_rhs", env)
	if(do_multi_lhs || do_multi_rhs){

	    res = multi_LHS_RHS(env, feNmlm)

	    return(res)

	}

	#
	# Regular estimation ####
	#


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
	params = get("params", env)

	isFixef = get("isFixef", env)
	onlyFixef = !isLinear && !isNonLinear && isFixef

	#
	# Model 0 + theta init
	#

	theta.init = get("theta.init", env)
	model0 = get_model_null(env, theta.init)

	# For the negative binomial:
	if(family == "negbin"){
	    theta.init = get("theta.init", env)
	    if(is.null(theta.init)){
	        theta.init = model0$theta
	    }

	    params = c(params, ".theta")
	    start = c(start, theta.init)
	    names(start) = params
	    upper = c(upper, 10000)
	    lower = c(lower, 1e-3)

	    assign("params", params, env)
	}

	assign("model0", model0, env)

	# the result
	res = get("res", env)

	# NO VARIABLE -- ONLY FIXED-EFFECTS
	if(onlyFixef){
		if(family == "negbin"){
			stop("To estimate the negative binomial model, you need at least one variable. (The estimation of the model with only the fixed-effects is not implemented.)")
		}

		res = femlm_only_clusters(env)
		res$onlyFixef = TRUE

		return(res)
	}

	# warnings => to avoid accumulation, but should appear even if the user stops the algorithm
	on.exit(warn_fixef_iter(env))

	#
	# Maximizing the likelihood
	#

	opt = try(stats::nlminb(start=start, objective=femlm_ll, env=env, lower=lower, upper=upper, gradient=gradient, hessian=hessian, control=opt.control), silent = TRUE)

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

		coef = opt$par
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

	var = NULL
	try(var <- solve(hessian_noBounded), silent = TRUE)
	if(is.null(var)){
		warning_msg = paste(warning_msg, "The information matrix is singular: presence of collinearity.")
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

	zvalue = coef/se
	pvalue = 2*pnorm(-abs(zvalue))

	# We add the information on the bound for the se & update the var to drop the bounded vars
	se_format = se
	if(any(isBounded)){
		se_format[!isBounded] = decimalFormat(se_format[!isBounded])
		se_format[isBounded] = boundText
	}

	coeftable = data.frame("Estimate"=coef, "Std. Error"=se_format, "z value"=zvalue, "Pr(>|z|)"=pvalue, stringsAsFactors = FALSE)
	names(coeftable) = c("Estimate", "Std. Error", "z value",  "Pr(>|z|)")
	row.names(coeftable) = params

	attr(se, "type") = attr(coeftable, "type") = "Standard"

	mu_both = get_mu(coef, env, final = TRUE)
	mu = mu_both$mu
	exp_mu = mu_both$exp_mu

	# calcul pseudo r2
	loglik = -opt$objective # moins car la fonction minimise
	ll_null = model0$loglik

	# dummies are constrained, they don't have full dof (cause you need to take one value off for unicity)
	# this is an approximation, in some cases there can be more than one ref. But good approx.
	nparams = res$nparams
	pseudo_r2 = 1 - (loglik - nparams + 1) / ll_null

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
			message("Very high value of theta (", theta, "). There is no sign of overdispersion, you may consider a Poisson model.")
		}

	}

	class(res) = "fixest"

	if(verbose > 0){
		cat("\n")
	}

	do_summary = get("do_summary", env)
	if(do_summary){
	    se = get("se", env)
	    cluster = get("cluster", env)
	    lean = get("lean", env)
	    dof = get("dof", env)
	    agg = get("agg", env)
	    summary_flags = get("summary_flags", env)

	    # To compute the RMSE and lean = TRUE
	    if(lean) res$ssr = cpp_ssq(res$residuals)

	    res = summary(res, se = se, cluster = cluster, dof = dof, agg = agg, lean = lean, summary_flags = summary_flags)
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
    #   argument => likely I'll need a match.call argument

    x = gsub("\n+$", "", x)

    if(grepl("^Error (in|:|: in) (fe|fixest|fun)[^\n]+\n", x)){
        res = gsub("^Error (in|:|: in) (fe|fixest|fun)[^\n]+\n *(.+)", "\\3", x)
    } else if(grepl("[Oo]bject '.+' not found", x) || grepl("memory|cannot allocate", x)) {
        res = x
    } else {
       res = paste0(x, "\nThis error was unforeseen by the author of the function ", origin, ". If you think your call to the function is legitimate, could you report?")
    }
    res
}


####
#### Multiple estimation tools ####
####



multi_split = function(env, fun){
    split = get("split", env)
    split.full = get("split.full", env)
    split.items = get("split.items", env)
    split.name = get("split.name", env)

    assign("do_split", FALSE, env)

    res_all = list()
    n_split = length(split.items)
    index = NULL
    all_names = NULL
    is_multi = FALSE
    for(i in 0:n_split){
        if(i == 0){
            if(split.full){
                my_env = reshape_env(env)
                my_res = fun(env = my_env)
            } else {
                next
            }
        } else {
            my_res = fun(env = reshape_env(env, obs2keep = which(split == i)))
        }

        res_all[[length(res_all) + 1]] = my_res
    }

    if(split.full){
        split.items = c("Full sample", split.items)
    }

    index = list(sample = length(res_all))
    all_names = list(sample = split.items, split.name = split.name)

    # result
    res_multi = setup_multi(index, all_names, res_all)

    return(res_multi)
}



multi_LHS_RHS = function(env, fun){
    do_multi_lhs = get("do_multi_lhs", env)
    do_multi_rhs = get("do_multi_rhs", env)

    assign("do_multi_lhs", FALSE, env)
    assign("do_multi_rhs", FALSE, env)

    nthreads = get("nthreads", env)

    # IMPORTANT NOTE:
    # contrary to feols, the preprocessing is only a small fraction of the
    # computing time in ML models
    # Therefore we don't need to optimize processing as hard as in FEOLS
    # because the gains are only marginal

    fml = get("fml", env)

    # LHS
    lhs_names = get("lhs_names", env)
    lhs = get("lhs", env)

    if(do_multi_lhs == FALSE){
        lhs = list(lhs)
    }

    # RHS
    if(do_multi_rhs){
        rhs_info_stepwise = get("rhs_info_stepwise", env)
        multi_rhs_fml_full = rhs_info_stepwise$fml_all_full
        multi_rhs_fml_sw = rhs_info_stepwise$fml_all_sw
        multi_rhs_cumul = rhs_info_stepwise$is_cumul

        linear_core = get("linear_core", env)
        rhs_sw = get("rhs_sw", env)

    } else {
        multi_rhs_fml_full = list(.xpd(rhs = fml[[3]]))
        multi_rhs_cumul = FALSE
        linear.mat = get("linear.mat", env)
        linear_core = list(left = linear.mat, right = 1)
        rhs_sw = list(1)
    }

    isLinear_left = length(linear_core$left) > 1
    isLinear_right = length(linear_core$right) > 1

    n_lhs = length(lhs)
    n_rhs = length(rhs_sw)
    res = vector("list", n_lhs * n_rhs)

    rhs_names = sapply(multi_rhs_fml_full, function(x) as.character(x)[[2]])

    for(i in seq_along(lhs)){
        for(j in seq_along(rhs_sw)){
            # reshaping the env => taking care of the NAs

            # Forming the RHS
            my_rhs = linear_core[1]

            if(multi_rhs_cumul){
                my_rhs[1 + 1:j] = rhs_sw[1:j]
            } else {
                my_rhs[2] = rhs_sw[j]
            }

            if(isLinear_right){
                my_rhs[[length(my_rhs) + 1]] = linear_core$right
            }

            n_all = lengths(my_rhs)
            if(any(n_all == 1)){
                my_rhs = my_rhs[n_all > 1]
            }

            if(length(my_rhs) == 0){
                my_rhs = 1
            } else {
                my_rhs = do.call("cbind", my_rhs)
            }

            if(length(my_rhs) == 1){
                is_na_current = !is.finite(lhs[[i]])
            } else {
                is_na_current = !is.finite(lhs[[i]]) | cpppar_which_na_inf_mat(my_rhs, nthreads)$is_na_inf
            }

            my_fml = .xpd(lhs = lhs_names[i], rhs = multi_rhs_fml_full[[j]])

            if(any(is_na_current)){
                my_env = reshape_env(env, which(!is_na_current), lhs = lhs[[i]], rhs = my_rhs, fml_linear = my_fml)
            } else {
                # We still need to check the RHS (only 0/1)
                my_env = reshape_env(env, lhs = lhs[[i]], rhs = my_rhs, fml_linear = my_fml, check_lhs = TRUE)
            }

            my_res = fun(env = my_env)

            res[[index_2D_to_1D(i, j, n_rhs)]] = my_res
        }
    }

    # Meta information for fixest_multi

    index = list(lhs = n_lhs, rhs = n_rhs)
    all_names = list(lhs = lhs_names, rhs = rhs_names)

    # result
    res_multi = setup_multi(index, all_names, res)

    return(res_multi)
}


multi_fixef = function(env, estfun){
    # Honestly had I known it was so painful, I wouldn't have done it...
    assign("do_multi_fixef", FALSE, env)

    multi_fixef_fml_full = get("multi_fixef_fml_full", env)
    combine.quick = get("combine.quick", env)
    fixef.rm = get("fixef.rm", env)
    family = get("family", env)
    origin_type = get("origin_type", env)
    nthreads = get("nthreads", env)

    data = get("data", env)

    n_fixef = length(multi_fixef_fml_full)

    data_results = list()
    for(i in 1:n_fixef){

        fml_fixef = multi_fixef_fml_full[[i]]

        if(length(all.vars(fml_fixef)) > 0){

            #
            # Evaluation of the fixed-effects
            #

            fixef_terms_full = fixef_terms(fml_fixef)

            # fixef_terms_full computed in the formula section
            fixef_terms = fixef_terms_full$fml_terms

            # FEs
            fixef_df = error_sender(prepare_df(fixef_terms_full$fe_vars, data, combine.quick),
                                     "Problem evaluating the fixed-effects part of the formula:\n")

            fixef_vars = names(fixef_df)

            # Slopes
            isSlope = any(fixef_terms_full$slope_flag != 0)
            slope_vars_list = list(0)
            if(isSlope){

                slope_df = error_sender(prepare_df(fixef_terms_full$slope_vars, data),
                                         "Problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n")

                slope_flag = fixef_terms_full$slope_flag
                slope_vars = fixef_terms_full$slope_vars
                slope_vars_list = fixef_terms_full$slope_vars_list

                # Further controls
                not_numeric = !sapply(slope_df, is.numeric)
                if(any(not_numeric)){
                    stop("In the fixed-effects part of the formula (i.e. in ", as.character(fml_fixef[2]), "), variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope_df)[not_numeric], "s.is.quote"), " not.")
                }

                # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect
                onlySlope = all(slope_flag < 0)

            }

            # fml update
            fml_fixef = .xpd(rhs = fixef_terms)

            #
            # NA
            #

            for(j in seq_along(fixef_df)){
                if(!is.numeric(fixef_df[[j]]) && !is.character(fixef_df[[j]])){
                    fixef_df[[j]] = as.character(fixef_df[[j]])
                }
            }

            is_NA = !complete.cases(fixef_df)

            if(isSlope){
                # Convert to double
                who_not_double = which(sapply(slope_df, is.integer))
                for(j in who_not_double){
                    slope_df[[j]] = as.numeric(slope_df[[j]])
                }

                info = cpppar_which_na_inf_df(slope_df, nthreads)
                if(info$any_na_inf){
                    is_NA = is_NA | info$is_na_inf
                }
            }

            if(any(is_NA)){
                # Remember that isFixef is FALSE so far => so we only change the reg vars
                my_env = reshape_env(env = env, obs2keep = which(!is_NA))

                # NA removal in fixef
                fixef_df = fixef_df[!is_NA, , drop = FALSE]

                if(isSlope){
                    slope_df = slope_df[!is_NA, , drop = FALSE]
                }
            } else {
                my_env = new.env(parent = env)
            }

            # We remove the linear part if needed


            if(get("do_multi_rhs", env)){
                linear_core = get("linear_core", my_env)
                if("(Intercept)" %in% colnames(linear_core$left)){
                    int_col = which("(Intercept)" %in% colnames(linear_core$left))
                    if(ncol(linear_core$left) == 1){
                        linear_core$left = 1
                    } else {
                        linear_core$left = linear_core$left[, -int_col, drop = FALSE]
                    }
                    assign("linear_core", linear_core, my_env)
                }
            } else {
                linear.mat = get("linear.mat", my_env)
                if("(Intercept)" %in% colnames(linear.mat)){
                    int_col = which("(Intercept)" %in% colnames(linear.mat))
                    if(ncol(linear.mat) == 1){
                        assign("linear.mat", 1, my_env)
                    } else {
                        assign("linear.mat", linear.mat[, -int_col, drop = FALSE], my_env)
                    }
                }
            }

            # We assign the fixed-effects
            lhs = get("lhs", my_env)

            # We delay the computation by using isSplit = TRUE and split.full = FALSE
            # Real QUF will be done in the last reshape env
            info_fe = setup_fixef(fixef_df = fixef_df, lhs = lhs, fixef_vars = fixef_vars, fixef.rm = fixef.rm, family = family, isSplit = TRUE, split.full = FALSE, origin_type = origin_type, isSlope = isSlope, slope_flag = slope_flag, slope_df = slope_df, slope_vars_list = slope_vars_list, nthreads = nthreads)

            fixef_id        = info_fe$fixef_id
            fixef_names     = info_fe$fixef_names
            fixef_sizes     = info_fe$fixef_sizes
            fixef_table     = info_fe$fixef_table
            sum_y_all       = info_fe$sum_y_all
            lhs             = info_fe$lhs

            obs2remove      = info_fe$obs2remove
            fixef_removed   = info_fe$fixef_removed
            message_fixef   = info_fe$message_fixef

            slope_variables = info_fe$slope_variables
            slope_flag      = info_fe$slope_flag

            fixef_id_res    = info_fe$fixef_id_res
            fixef_sizes_res = info_fe$fixef_sizes_res
            new_order       = info_fe$new_order

            assign("isFixef", TRUE, my_env)
            assign("new_order_original", new_order, my_env)
            assign("fixef_names", fixef_names, my_env)
            assign("fixef_vars", fixef_vars, my_env)

            assign_fixef_env(env, family, origin_type, fixef_id, fixef_sizes, fixef_table, sum_y_all, slope_flag, slope_variables, slope_vars_list)

            #
            # Formatting the fixef stuff from res
            #

            # fml & fixef_vars => other stuff will be taken care of in reshape
            res = get("res", my_env)
            res$fml_all$fixef = fml_fixef
            res$fixef_vars = fixef_vars
            if(isSlope){
                res$fixef_terms = fixef_terms
            }
            assign("res", res, my_env)

            #
            # Last reshape
            #

            my_env_est = reshape_env(my_env, assign_fixef = TRUE)

        } else {
            # No fixed-effect // new.env is indispensable => otherwise multi RHS/LHS not possible
            my_env_est = reshape_env(env)
        }

        data_results[[i]] = estfun(env = my_env_est)
    }

    index = list(fixef = n_fixef)
    fixef_names = sapply(multi_fixef_fml_full, function(x) as.character(x)[[2]])
    all_names = list(fixef = fixef_names)

    res_multi = setup_multi(index, all_names, data_results)

    if("lhs" %in% names(attr(res_multi, "meta")$index)){
        res_multi = res_multi[lhs = TRUE]
    }

    return(res_multi)

}





















