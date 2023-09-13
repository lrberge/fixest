#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Apr 23 15:35:53 2022
# ~: A few user-level + many internal funs
#----------------------------------------------#



#' Collinearity diagnostics for `fixest` objects
#'
#' In some occasions, the optimization algorithm of [`femlm`] may fail to converge, or the variance-covariance matrix may not be available. The most common reason of why this happens is collinearity among variables. This function helps to find out which set of variables is problematic.
#'
#'
#' @param x A `fixest` object obtained from, e.g. functions [`femlm`], [`feols`] or [`feglm`].
#' @param verbose An integer. If higher than or equal to 1, then a note is prompted at each step of the algorithm. By default `verbose = 0` for small problems and to 1 for large problems.
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
		stop("Argument `x` must be a fixest object.")
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
	data = fetch_data(x, "To apply function `collinearity`, ")

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
				message = paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with fixed-effects `", names(cluster)[q], "`.")

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

	            message = paste0("Variable", enumerate_items(collin_var, "s.is.quote"), " collinear with variable with varying slope `", slope_vars[q], "` (on `", slope_fe[q], "`).")

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
				message = paste0("Variable `", dict_name[v], "` is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ".")

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
					message = paste0("Variable `", dict_name[v], "` is collinear with variable", enumerate_items(dict_name[collin_var], "s.quote"), ", together with the fixed-effects ", dum_names[id_cluster], ".")

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
#' @param fml Either a formula of the type `var1 + ... + varN ~ treat` or `var1 + ... + varN ~ treat | post`. Either a data.frame/matrix containing all the variables for which the means are to be computed (they must be numeric of course). Both the treatment and the post variables must contain only exactly two values. You can use a point to select all the variables of the data set: `. ~ treat`.
#' @param base A data base containing all the variables in the formula `fml`.
#' @param treat_var Only if argument `fml` is *not* a formula. The vector identifying the treated and the control observations (the vector can be of any type but must contain only two possible values). Must be of the same length as the data.
#' @param post_var Only if argument `fml` is *not* a formula. The vector identifying the periods (pre/post) of the observations (the vector can be of any type but must contain only two possible values). The first value (in the sorted sense) of the vector is taken as the pre period. Must be of the same length as the data.
#' @param treat_first Which value of the 'treatment' vector should appear on the left? By default the max value appears first (e.g. if the treatment variable is a 0/1 vector, 1 appears first).
#' @param tex Should the result be displayed in Latex? Default is `FALSE`. Automatically set to `TRUE` if the table is to be saved in a file using the argument `file`.
#' @param treat_dict A character vector of length two. What are the names of the treated and the control? This should be a dictionary: e.g. `c("1"="Treated", "0" = "Control")`.
#' @param dict A named character vector. A dictionary between the variables names and an alias. For instance `dict=c("x"="Inflation Rate")` would replace the variable name `x` by \dQuote{Inflation Rate}.
#' @param file A file path. If given, the table is written in Latex into this file.
#' @param replace Default is `TRUE`, which means that when the table is exported, the existing file is not erased.
#' @param title Character string giving the Latex title of the table. (Only if exported.)
#' @param label Character string giving the Latex label of the table. (Only if exported.)
#' @param raw Logical, default is `FALSE`. If `TRUE`, it returns the information without formatting.
#' @param indiv Either the variable name of individual identifiers, a one sided formula, or a vector. If the data is that of a panel, this can be used to track the number of individuals per group.
#' @param prepostnames Only if there is a 'post' variable. The names of the pre and post periods to be displayed in Latex. Default is `c("Before", "After")`.
#' @param diff.inv Logical, default to `FALSE`. Whether to inverse the difference.
#'
#' @details
#' By default, when the user tries to apply this function to nun-numeric variables, an error is raised. The exception is when the all variables are selected with the dot (like in `. ~ treat`. In this case, non-numeric variables are automatically omitted (with a message).
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
                stop("To use `indiv` as a formula, you must provide the argument `base`.")
            }

            if(length(indiv_varname) > 1){
                stop("The argument `indiv` must refer to only one variable.")
            }

            if(!all(all.vars(indiv) %in% names(base))){
                pblm = setdiff(all.vars(indiv), names(base))
                stop("In argument `indiv`: the variable", enumerate_items(pblm, "s.is")," not in the data set.")
            }

            indiv_var = try(eval(str2lang(indiv_varname), base), silent = TRUE)
            if("try-error" %in% class(indiv_var)){
                stop("Evaluation of `indiv` raises and error:\n", indiv_var)
            }
        } else if(length(indiv) == 1 && is.character(indiv)){
            indiv_varname = indiv

            if(missing(base) || !indiv %in% names(base)){
                stop("To use `indiv` as a `character string` you must provide the argument `base`.")
            } else {
                indiv_var = base[[indiv]]
            }
        } else {
            if(!missing(base) && length(indiv) != NROW(base)){
                stop("The length of `indiv` must be the same as the data.")
            } else if(missing(base) && length(indiv) != NROW(fml)){
                stop("The length of `indiv` must be the same as the data.")
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
            stop("If you provide a formula, a data.frame must be given in argument `base`.")
        }

        # Creation of x and the condition
        if(!length(fml_in) == 3){
            stop("The formula must be of the type `var1 + ... + varN ~ treat` or `var1 + ... + varN ~ treat | post`.")
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
            stop("Evaluation of the `treatment` variable raises and error: \n", treat_var)
        }

        if(!is.null(pipe)){
            if(!(all(all.vars(pipe) %in% names(base)))){
                pblm = setdiff(all.vars(pipe), names(base))
                stop("In the evaluation of the `post` variable: ", enumerate_items(pblm, "is"), " not in the data set.")
            }

            post_var = try(eval(pipe[[2]], base), silent = TRUE)
            if("try-error" %in% class(post_var)){
                stop("Evaluation of the `post` variable raises and error: \n", treat_var)
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
                stop("The variable", enumerate_items(pblm, "s.is"), " not in the data set (", deparse_short(substitute(base)), ").")
            }

            # Evaluation control
            base_small = head(base, 10)
            var2eval = attr(terms(fml_x), "term.labels")
            var2eval = gsub(":", "*", var2eval)
            for(i in seq_along(var2eval)){
                var = var2eval[i]
                x_small = try(eval(parse(text=var), base_small), silent = TRUE)
                if("try-error" %in% class(x_small)){
                    stop("Evaluation of the variable `", var, "` raises and error:\n", x_small)
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
            stop("If argument `fml` is not a formula, you must provide the argument `treat_var`.")
        } else {
            mat_vars = fml_in

            if(NROW(mat_vars) != length(treat_var)){
                stop("The arguments `x` and `treat_var` must be of the same length.")
            }

            if(!missing(post_var) && !is.null(post_var)){
                if(NROW(mat_vars) != length(post_var)){
                    stop("If provided, the argument `post_var` must be of the same length of `x`.")
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
                stop("If not a formula, argument `fml` must be a data.frame or a matrix. Currently it is of class ", class(mat_vars)[1], ".")
            } else if(!is.numeric(mat_vars) && !is.logical(mat_vars)){
                stop("If not a formula, argument `fml` must be a data.frame or a matrix with numeric variables. Currenlty its is not numeric.")
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
                    msg = "`treatment` and `post` variables."
                } else {
                    msg = "`post` variable."
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
        stop("This function supports only 2 conditional values for the `post` variable. Currently, it contains ", msg)
    }

    if(usePost){
        check_arg(prepostnames, "character vector len(2) no na")
    }

    treat_cases = sort(unique(treat_var), decreasing = TRUE)
    if(!missing(treat_dict) && !is.null(treat_dict)){

        if(!isVector(treat_dict) || is.null(names(treat_dict))){
            stop("The argument `treat_dict` must be a named character vector.")
        }

        pblm = setdiff(treat_cases, names(treat_dict))
        if(length(pblm) > 0){
            stop("The value", enumerate_items(pblm, "s.is"), " not in `treat_dict`.")
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
            stop("Argument `treat_first` must be an element of the treated variable.")
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
#' Treat a variable as a factor, or interacts a variable with a factor. Values to be dropped/kept from the factor can be easily set. Note that to interact fixed-effects, this function should not be used: instead use directly the syntax `fe1^fe2`.
#'
#'
#' @inheritParams bin
#'
#' @param factor_var  A vector (of any type) that will be treated as a factor. You can set references (i.e. exclude values for which to create dummies) with the `ref` argument.
#' @param var A variable of the same length as `factor_var`. This variable will be interacted with the factor in `factor_var`. It can be numeric or factor-like. To force a numeric variable to be treated as a factor, you can add the `i.` prefix to a variable name. For instance take a numeric variable `x_num`: `i(x_fact, x_num)` will treat `x_num` as numeric while `i(x_fact, i.x_num)` will treat `x_num` as a factor (it's a shortcut to `as.factor(x_num)`).
#' @param ref A vector of values to be taken as references from `factor_var`. Can also be a logical: if `TRUE`, then the first value of `factor_var` will be removed. If `ref` is a character vector, partial matching is applied to values; use "@" as the first character to enable regular expression matching. See examples.
#' @param keep A vector of values to be kept from `factor_var` (all others are dropped). By default they should be values from `factor_var` and if `keep` is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param ref2 A vector of values to be dropped from `var`. By default they should be values from `var` and if `ref2` is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param keep2 A vector of values to be kept from `var` (all others are dropped). By default they should be values from `var` and if `keep2` is a character vector partial matching is applied. Use "@" as the first character to enable regular expression matching instead.
#' @param bin2 A list or vector defining the binning of the second variable. See help for the argument `bin` for details (or look at the help of the function [`bin`]). You can use `.()` for `list()`.
#' @param ... Not currently used.
#'
#' @details
#' To interact fixed-effects, this function should not be used: instead use directly the syntax `fe1^fe2` in the fixed-effects part of the formula. Please see the details and examples in the help page of [`feols`].
#'
#' @return
#' It returns a matrix with number of rows the length of `factor_var`. If there is no interacted variable or it is interacted with a numeric variable, the number of columns is equal to the number of cases contained in `factor_var` minus the reference(s). If the interacted variable is a factor, the number of columns is the number of combined cases between `factor_var` and `var`.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' [`iplot`][fixest::coefplot] to plot interactions or factors created with `i()`, [`feols`] for OLS estimation with multiple fixed-effects.
#'
#' See the function [`bin`] for binning variables.
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
                stop("In i(): When `var` is equal to a product, please use I(", info[1], "*", info[2], ") instead of ", var_name, ".")
            } else {
                stop("The arguments `var` and `f` must be of the same length (currently ", length(var), " vs ", length(f), ").")
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

    NO_ERROR = FALSE
    if(is_calling_fun("fixest_model_matrix_extra", full_search = TRUE, full_name = TRUE)){
        NO_ERROR = TRUE
    }

    no_rm = TRUE
    id_drop = c()
    if(!missing(ref)){
        if(isTRUE(ref)){
            # We always delete the first value
            # Que ce soit items ici est normal (et pas f_items)
            id_drop = which(items == items[1])
        } else {
            id_drop = items_to_drop(f_items, ref, "factor_var", no_error = NO_ERROR)
        }
        ref_id = id_drop
    }


    if(!missing(keep)){
        id_drop = c(id_drop, items_to_drop(f_items, keep, "factor_var", keep = TRUE, no_error = NO_ERROR))
    }

    if(IS_INTER_FACTOR){
        if(!missing(ref2)){
            id_drop = c(id_drop, items_to_drop(var_items, ref2, "var", no_error = NO_ERROR))
        }

        if(!missing(keep2)){
            id_drop = c(id_drop, items_to_drop(var_items, keep2, "var", keep = TRUE, no_error = NO_ERROR))
        }
    }

    if(length(id_drop) > 0){
        id_drop = unique(sort(id_drop))
        if(length(id_drop) == length(items)){
            if(FROM_FIXEST) {
                # we return something neutral in an estimation
                return(rep(0, length(f)))
            }

            stop("All items from the interaction have been removed.")
        }
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





#' Expands formula macros
#'
#' Create macros within formulas and expand them with character vectors or other formulas.
#'
#' @inheritParams setFixest_fml
#'
#' @param fml A formula containing macros variables. Each macro variable must start with two dots. The macro variables can be set globally using `setFixest_fml`, or can be defined in `...`. Special macros of the form `..("regex")` can be used to fetch, through a regular expression, variables directly in a character vector (or in column names) given in the argument `data` (note that the algorithm tries to "guess" the argument data when nested in function calls \[see example\]). You can negate the regex by starting with a `"!"`. Square brackets have a special meaning: Values in them are evaluated and parsed accordingly. Example: `y~x.[1:2] + z.[i]` will lead to `y~x1+x2+z3` if `i==3`. You can trigger the auto-completion of variables by using the `'..'` suffix, like in `y ~ x..` which would include `x1` and `x2`, etc. See examples.
#' @param add Either a character scalar or a one-sided formula. The elements will be added to the right-hand-side of the formula, before any macro expansion is applied.
#' @param lhs If present then a formula will be constructed with `lhs` as the full left-hand-side. The value of `lhs` can be a one-sided formula, a call, or a character vector. Note that the macro variables wont be applied. You can use it in combination with the argument `rhs`. Note that if `fml` is not missing, its LHS will be replaced by `lhs`.
#' @param rhs If present, then a formula will be constructed with `rhs` as the full right-hand-side. The value of `rhs` can be a one-sided formula, a call, or a character vector. Note that the macro variables wont be applied. You can use it in combination with the argument `lhs`. Note that if `fml` is not missing, its RHS will be replaced by `rhs`.
#' @param data Either a character vector or a data.frame. This argument will only be used if a macro of the type `..("regex")` is used in the formula of the argument `fml`. If so, any variable name from `data` that matches the regular expression will be added to the formula.
#' @param frame The environment containing the values to be expanded with the dot square bracket operator. Default is `parent.frame()`.
#'
#' @details
#' In `xpd`, the default macro variables are taken from `getFixest_fml`. Any value in the `...` argument of `xpd` will replace these default values.
#'
#' The definitions of the macro variables will replace in verbatim the macro variables. Therefore, you can include multi-part formulas if you wish but then beware of the order of the macros variable in the formula. For example, using the `airquality` data, say you want to set as controls the variable `Temp` and `Day` fixed-effects, you can do `setFixest_fml(..ctrl = ~Temp | Day)`, but then `feols(Ozone ~ Wind + ..ctrl, airquality)` will be quite different from `feols(Ozone ~ ..ctrl + Wind, airquality)`, so beware!
#'
#'
#' @section Dot square bracket operator in formulas:
#'
#' In a formula, the dot square bracket (DSB) operator can: i) create manifold variables at once, or ii) capture values from the current environment and put them verbatim in the formula.
#'
#' Say you want to include the variables `x1` to `x3` in your formula. You can use `xpd(y ~ x.[1:3])` and you'll get `y ~ x1 + x2 + x3`.
#'
#' To summon values from the environment, simply put the variable in square brackets. For example: `for(i in 1:3) xpd(y.[i] ~ x)` will create the formulas `y1 ~ x` to `y3 ~ x` depending on the value of `i`.
#'
#' You can include a full variable from the environment in the same way: `for(y in c("a", "b")) xpd(.[y] ~ x)` will create the two formulas `a ~ x` and `b ~ x`.
#'
#' The DSB can even be used within variable names, but then the variable must be nested in character form. For example `y ~ .["x.[1:2]_sq"]` will create `y ~ x1_sq + x2_sq`. Using the character form is important to avoid a formula parsing error. Double quotes must be used. Note that the character string that is nested will be parsed with the function [`dsb`], and thus it will return a vector.
#'
#' By default, the DSB operator expands vectors into sums. You can add a comma, like in `.[, x]`, to expand with commas--the content can then be used within functions. For instance: `c(x.[, 1:2])` will create `c(x1, x2)` (and *not* `c(x1 + x2)`).
#'
#' In all `fixest` estimations, this special parsing is enabled, so you don't need to use `xpd`.
#'
#' One-sided formulas can be expanded with the DSB operator: let `x = ~sepal + petal`, then `xpd(y ~ .[x])` leads to `color ~ sepal + petal`.
#'
#' You can even use multiple square brackets within a single variable, but then the use of nesting is required. For example, the following `xpd(y ~ .[".[letters[1:2]]_.[1:2]"])` will create `y ~ a_1 + b_2`. Remember that the nested character string is parsed with [`dsb`], which explains this behavior.
#'
#' When the element to be expanded i) is equal to the empty string or, ii) is of length 0, it is replaced with a neutral element, namely `1`. For example, `x = "" ; xpd(y ~ .[x])` leads to `y ~ 1`.
#'
#' @section Regular expressions:
#'
#' You can catch several variable names at once by using regular expressions. To use regular expressions, you need to enclose it in the dot-dot or the regex function: `..("regex")` or `regex("regex")`. For example, `regex("Sepal")` will catch both the variables `Sepal.Length` and `Sepal.Width` from the `iris` data set. In a `fixest` estimation, the variables names from which the regex will be applied come from the data set. If you use `xpd`, you need to provide either a data set or a vector of names in the argument `data`.
#'
#' By default the variables are aggregated with a sum. For example in a data set with the variables x1 to x10, `regex("x(1|2)"` will yield `x1 + x2 + x10`. You can instead ask for "comma" aggregation by using a comma first, just before the regular expression: `y ~ sw(regex(,"x(1|2)"))` would lead to `y ~ sw(x1, x2, x10)`.
#'
#' Note that the dot square bracket operator (DSB, see before) is applied before the regular expression is evaluated. This means that `regex("x.[3:4]_sq")` will lead, after evaluation of the DSB, to `regex("x3_sq|x4_sq")`. It is a handy way to insert range of numbers in a regular expression.
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
#' [`setFixest_fml`] to set formula macros, and [`dsb`] to modify character strings with the DSB operator.
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
#' # Auto completion with '..' suffix
#' #
#'
#' # You can trigger variables autocompletion with the '..' suffix
#' # You need to provide the argument data
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' xpd(y ~ x.., data = base)
#'
#' # In fixest estimations, this is automatically taken care of
#' feols(y ~ x.., data = base)
#'
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
#' # The basic use id to add variables in the formula
#' x = c("x1", "x2")
#' xpd(y ~ .[x])
#'
#' # Alternatively, one-sided formulas can be used and their content will be inserted verbatim
#' x = ~x1 + x2
#' xpd(y ~ .[x])
#'
#' # You can create multiple variables at once
#' xpd(y ~ x.[1:5] + z.[2:3])
#'
#' # You can summon variables from the environment to complete variables names
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
#' #
#' # In non-fixest functions: guessing the data allows to use regex
#' #
#'
#' # When used in non-fixest functions, the algorithm tries to "guess" the data
#' # so that ..("regex") can be directly evaluated without passing the argument 'data'
#' data(longley)
#' lm(xpd(Armed.Forces ~ Population + ..("GNP|ployed")), longley)
#'
#' # same for the auto completion with '..'
#' lm(xpd(Armed.Forces ~ Population + GN..), longley)
#'
#'
xpd = function(fml, ..., add = NULL, lhs, rhs, data = NULL, frame = parent.frame()){

    if(MISSNULL(data)){
        # We "guess" the data
        sc = sys.calls()
        n_sc = length(sc)
        if(n_sc > 1){
            mc = tryCatch(match.call(definition = sys.function(n_sc - 1), call = sys.call(n_sc - 1)), error = function(e) NULL)
            if("data" %in% names(mc)){
                data = tryCatch(eval(mc$data, parent.frame(2)), error = function(e) NULL)
            }
        }
    }

    .xpd(fml = fml, ..., add = add, lhs = lhs, rhs = rhs, data = data, check = TRUE,
         macro = TRUE, frame = frame)
}

.xpd = function(fml, ..., add = NULL, lhs, rhs, data = NULL, check = FALSE, macro = FALSE, frame = .GlobalEnv){

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

        attr(fml, ".Environment") = frame

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

            attr(fml, ".Environment") = frame
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

    data_vars = NULL
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
                        stop("The argument `data` must contain variables names. It is currently a matrix without column names.")
                    }
                } else if(is.data.frame(data)){
                    data = names(data)
                }

                data_vars = data
            }

            fml_dp_split = strsplit(fml_dp, '(?<![[:alnum:]._])(regex|\\.\\.)\\((?=[,"])',
                                    perl = TRUE)[[1]]

            res = fml_dp_split
            for(i in 2:length(res)){

                re = sub('"\\).*', "", res[i])

                re_width = nchar(re)

                is_comma = grepl("^,", re)
                re = gsub("^,? ?\"", "", re)

                re = dot_square_bracket(re, frame, regex = TRUE)

                if(is_data){
                    if(substr(re, 1, 1) == "!"){
                        # negation
                        re = str_trim(re, 1)
                        vars = data[!grepl(re, data, perl = TRUE)]
                    } else {
                        vars = grep(re, data, value = TRUE, perl = TRUE)
                    }

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
                               fit_screen(paste0("...TEXT: `", fml_txt, "`")),
                               "\nPROBLEM: see below", up = 1)
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
                               fit_screen(paste0("...TEXT: `", fml_txt, "`")),
                               "\nPROBLEM: see below", up = 1)
        }

    }

    if(is_data){
        var_to_complete = grep("[[:alnum:]]\\.\\.$", all.vars(fml), value = TRUE)
        n_var = length(var_to_complete)
        if(n_var > 0){

            if(is.null(data_vars)){
                check_arg(data, "character vector no na | matrix | data.frame")

                if(is.matrix(data)){
                    data = colnames(data)
                    if(is.null(data)){
                        stop("The argument `data` must contain variables names. It is currently a matrix without column names.")
                    }
                } else if(is.data.frame(data)){
                    data = names(data)
                }

                data_vars = data
            }

            fml_txt = deparse_long(fml)
            for(i in 1:n_var){
                var = var_to_complete[i]
                vars_filled = data_vars[startsWith(data_vars, str_trim(var, n_last = 2))]
                if(length(vars_filled) == 0) next

                vars_new = paste0(vars_filled, collapse = " + ")

                pattern = dsb("(?<![[:alnum:]._])\\Q.[var]\\E(?![[:alnum:]._])")

                fml_txt = gsub(pattern, vars_new, fml_txt, perl = TRUE)
            }

            fml = error_sender(as.formula(fml_txt, frame),
                               "Expansion of variables ending with `..` did not work. Coercion of the following text to a formula led to an error.\n",
                               fit_screen(paste0("...TEXT: `", fml_txt, "`")),
                               "\nPROBLEM: see below", up = 1)
        }
    }

    fml
}


#' Centers a set of variables around a set of factors
#'
#' User-level access to internal demeaning algorithm of `fixest`.
#'
#' @inheritSection feols Varying slopes
#'
#' @param X A matrix, vector, data.frame or a list OR a formula OR a [`feols`] estimation. If equal to a formula, then the argument `data` is required, and it must be of the type: `x1 + x2 ~ f1 + fe2` with on the LHS the variables to be centered, and on the RHS the factors used for centering. Note that you can use variables with varying slopes with the syntax `fe[v1, v2]` (see details in [`feols`]). If a `feols` estimation, all variables (LHS+RHS) are demeaned and then returned (only if it was estimated with fixed-effects). Otherwise, it must represent the data to be centered. Of course the number of observations of that data must be the same as the factors used for centering (argument `f`).
#' @param f A matrix, vector, data.frame or list. The factors used to center the variables in argument `X`. Matrices will be coerced using `as.data.frame`.
#' @param slope.vars A vector, matrix or list representing the variables with varying slopes. Matrices will be coerced using `as.data.frame`. Note that if this argument is used it MUST be in conjunction with the argument `slope.flag` that maps the factors to which the varying slopes are attached. See examples.
#' @param slope.flag An integer vector of the same length as the number of variables in `f` (the factors used for centering). It indicates for each factor the number of variables with varying slopes to which it is associated. Positive values mean that the raw factor should also be included in the centering, negative values that it should be excluded. Sorry it's complicated... but see the examples it may get clearer.
#' @param data A data.frame containing all variables in the argument `X`. Only used if `X` is a formula, in which case `data` is mandatory.
#' @param weights Vector, can be missing or NULL. If present, it must contain the same number of observations as in `X`.
#' @param nthreads Number of threads to be used. By default it is equal to `getFixest_nthreads()`.
#' @param notes Logical, whether to display a message when NA values are removed. By default it is equal to `getFixest_notes()`.
#' @param iter Number of iterations, default is 2000.
#' @param tol Stopping criterion of the algorithm. Default is `1e-6`. The algorithm stops when the maximum absolute increase in the coefficients values is lower than `tol`.
#' @param na.rm Logical, default is `TRUE`. If `TRUE` and the input data contains any NA value, then any observation with NA will be discarded leading to an output with less observations than the input. If `FALSE`, if NAs are present the output will also be filled with NAs for each NA observation in input.
#' @param as.matrix Logical, if `TRUE` a matrix is returned, if `FALSE` it will be a data.frame. The default depends on the input, if atomic then a matrix will be returned.
#' @param im_confident Logical, default is `FALSE`. FOR EXPERT USERS ONLY! This argument allows to skip some of the preprocessing of the arguments given in input. If `TRUE`, then `X` MUST be a numeric vector/matrix/list (not a formula!), `f` MUST be a list, `slope.vars` MUST be a list, `slope.vars` MUST be consistent with `slope.flag`, and `weights`, if given, MUST be numeric (not integer!). Further there MUST be not any NA value, and the number of observations of each element MUST be consistent. Non compliance to these rules may simply lead your R session to break.
#' @param fixef.reorder Logical, default is `TRUE`. Whether to reorder the fixed-effects by frequencies before feeding them into the algorithm. If `FALSE`, the original fixed-effects order provided by the user is maintained. In general, reordering leads to faster and more precise performance.
#' @param ... Not currently used.
#'
#' @return
#' It returns a data.frame of the same number of columns as the number of variables to be centered.
#'
#' If `na.rm = TRUE`, then the number of rows is equal to the number of rows in input minus the number of NA values (contained in `X`, `f`, `slope.vars` or `weights`). The default is to have an output of the same number of observations as the input (filled with NAs where appropriate).
#'
#' A matrix can be returned if `as.matrix = TRUE`.
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
#' This function extracts the observations used in `fixest` estimation.
#'
#' @param x A `fixest` object.
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
#' (obs_setosa = obs(est_split[[1]]))
#' (obs_versi = obs(est_split[sample = "versi", drop = TRUE]))
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




#' Check the fixed-effects convergence of a `feols` estimation
#'
#' Checks the convergence of a `feols` estimation by computing the first-order conditions of all fixed-effects (all should be close to 0)
#'
#' @param x A [`feols`] estimation that should contain fixed-effects.
#' @param object An object returned by `check_conv_feols`.
#' @param type Either "short" (default) or "detail". If "short", only the maximum absolute FOC are displayed, otherwise the 2 smallest and the 2 largest FOC are reported for each fixed-effect and each variable.
#' @param ... Not currently used.
#'
#' Note that this function first re-demeans the variables, thus possibly incurring some extra computation time.
#'
#' @return
#' It returns a list of `N` elements, `N` being the number of variables in the estimation (dependent variable + explanatory variables +, if IV, endogenous variables and instruments). For each variable, all the first-order conditions for each fixed-effect are returned.
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






#' Replicates `fixest` objects
#'
#' Simple function that replicates `fixest` objects while (optionally) computing different standard-errors. Useful mostly in combination with [`etable`] or [`coefplot`].
#'
#' @param x Either a `fixest` object, either a list of `fixest` objects created with `.l()`.
#' @param times Integer vector giving the number of repetitions of the vector of elements. By default `times = 1`. It must be either of length 1, either of the same length as the argument `x`.
#' @param each Integer scalar indicating the repetition of each element. Default is 1.
#' @param vcov A list containing the types of standard-error to be computed, default is missing. If not missing, it must be of the same length as `times`, `each`, or the final vector. Note that if the arguments `times` and `each` are missing, then `times` becomes equal to the length of `vcov`. To see how to summon a VCOV, see the dedicated section in the [vignette](https://lrberge.github.io/fixest/articles/fixest_walkthrough.html#the-vcov-argument-1).
#' @param ... In `.l()`: `fixest` objects. In `rep()`: not currently used.
#'
#' @details
#' To apply `rep.fixest` on a list of `fixest` objects, it is absolutely necessary to use `.l()` and not `list()`.
#'
#' @return
#' Returns a list of the appropriate length. Each element of the list is a `fixest` object.
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
            res = dots[[1]]
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

            if("fixest_multi" %in% class(obj)){
                for(j in seq_along(obj)){
                    res[[length(res) + 1]] = obj[[j]]
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
                                       "Dot square bracket operator: Evaluation of `.[",
                                       x_split[i], "]` led to an error:",
                                       up = up + 1)

                operator[i] = op
            }
        }


        if(do_lang){
            my_call = error_sender(str2lang(x_split[i]),
                                   "Dot square bracket operator: Evaluation of `.[",
                                   x_split[i], "]` led to an error:",
                                   up = up + 1)
        }

        if(is.character(my_call) && grepl(".[", my_call, fixed = TRUE)){
            # Nested call
            value = .dsb(my_call, frame = frame)

        } else {
            # Informative error message
            value = error_sender(as.character(eval(my_call, frame)),
                                 "Dot square bracket operator: Evaluation of `.[",
                                 x_split[i], "]` led to an error:",
                                 up = up + 1)

            if(length(value) == 2 && value[1] == "~"){
                value = char_to_vars(value[2])
            }
        }

        if(length(value) == 0 || (is.character(value) && length(value) == 1 && nchar(value) == 0)){
            # Neutral element in formulas
            value = "1"
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
        stop("The special operator `^` can only be used in the fixed-effects part of the formula. Please use `:` instead.")
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
        } else {
            # => predict, newdata
            # case i() + anything that requires evaluation based on the raw data (poly, factor, etc)
            # is there any drawback?
            fml = formula(mf)
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
            stop("In `model.matrix`, the variable", enumerate_items(pblm, "is.s.quote"), " in the formula but not in the argument `data`. Use `subset = TRUE` to enable the creation of partial data.")
        }

    } else {
        vars_keep = names(newdata)

        if(is.character(subset)){
            # ex: subset = c("x1$", "x2$")

            vars_keep = keep_apply(vars_keep, subset)
            if(length(vars_keep) == 0){
                stop("The variables in `subset` do not match any variable in the `data`.")
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
                stop("Due to the use of the argument `subset`, not a single variable is left.")
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
        data = fetch_data(object, "To apply `model.matrix.fixest`, ")

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

    I_IGNORE_ERRORS = TRUE

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
                    msg = "Wrong formatting of the fixed-effects interactions. The `^` operator should not be within parentheses."
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

    check_value(nthreads, "integer scalar GE{0} | numeric scalar GT{0} LT{1}", .message = paste0("The argument `nthreads` must be an integer lower or equal to the number of threads available (", max(cpp_get_nb_threads(), 1), "). It can be equal to 0 which means all threads. Alternatively, if equal to a number strictly between 0 and 1, it represents the fraction of all threads to be used."))

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
        stop("Internal error: keep_first should not be used with `keep` or `ref`.")
    }

    if(is.factor(items)){
        items = as.character(items)
    }

    if(is.factor(x)){
        x = as.character(x)
    }


    if(inherits(x, "formula")){

        if(!identical(all.vars(x), "x")){
            stop_up("In argument `", argname, "`, the formula must contain a single variable name: `x`. So far `", deparse_long(x), "` is not valid.")
        }

        if(length(x) > 2){
            stop_up("In argument `", argname, "`, if a formula, it must be one-sided. Problem: `", deparse_long(x), "` is two-sided.")
        }

        is_here = error_sender(eval(x[[2]], list(x = items)),
                               "In argument `", argname, "`, the evaluation of the formula led to an error:")
        if(length(is_here) != length(items)){
            stop_up("In argument `", argname, "`, the evaluation of the formula must return a logical vector of the same length as `x`. Problem: `", deparse_long(x), "` returns a vector of length ", length(is_here), " (expected: ", length(items), ").")
        }

        if(!is.logical(is_here)){
            stop_up("In argument `", argname, "`, the evaluation of the formula must return a logical vector. Problem: `", deparse_long(x), "` is not logical (instead it is of class ", enumerate_items(class(is_here)), ").")
        }

        is_here = !is.na(is_here) & is_here

        if(!no_error && !any(is_here)){
            stop_up("In argument `", argname, "`, the evaluation of the formula must match at least one value. Problem: `", deparse_long(x), "` does not match any.")
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
                    stop_up("In argument `", argname, "`, the regular expression `", pattern, "` does not match any value of `", varname, "`.")
                }
                all_x = c(all_x, new_x)
            } else {
                # B) partial matching
                if(no_error){
                    qui_ok = pmatch(tolower(my_x), tolower(items), duplicates.ok = TRUE)
                    my_x = items[qui_ok]
                    my_x = my_x[!is.na(my_x)]
                } else {
                    check_value_plus(my_x, "match", .choices = unique(items), .message = paste0("The argument `", argname, "` should contain values of the variable `", varname, "`."))
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
                        stop_up("In argument `", argname, "`, the value ", x[1], " does not match any value of `", varname, "`.")
                    }

                    id_drop = unique(c(which(items %in% x[1]), id_drop))
                }

            }

            if(length(id_drop) == 0 && !no_error){
                stop_up("In argument `", argname, "`, the value", plural_len(x, "s.don't"), " match any value of `", varname, "`.")
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
                                                " with `.[]` led to an error:"))

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
                                               " with `.[]` led to an error:"))
            }
            bin[[i]] = unlist(value)
        }

    }

    #
    # cut:: => special treatment, we short circuit
    #

    if(!is.list(bin) && (is.character(bin) && any(grepl("^cut::", bin)))){

        if(!grepl("^cut::", bin[1])){
            stop_up("To use the special binning `cut::values`, it must be in the first position of the vector, possibly followed by custom names. Currently `cut::values` is in position ", which(grepl("^cut::", bin))[1], ".")
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
                stop_up("In the special binning `cut::values`, the number of custom names is greater than the number of bins (", n_names, " vs ", n_cuts, ").")
            }
        }

        if(!is.numeric(x)){
            stop_up("To use the special binning `cut::values`, the variable `", varname, "` must be numeric. Currently this is not the case (it is of class ", enumerate_items(class(x)), " instead).")
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
                stop_up("To use the special binning `bin::digit`, the argument `", argname, "` must be of length 1. Currently it is of length ", length(bin), ".")
            }

            d = gsub("^!?!?bin::", "", bin)
            if(any(grepl("[^[:digit:]]", d))){
                bin_type = gsub("^(!?!?bin).*", "\\1", bin)
                stop_up("In the argument bin, the special binning must be of the form `", bin_type, "::digit`. Currently this is not the case for `", bin_type, "::", d, "`.")
            }
            d = as.numeric(d)

            consecutive = grepl("^!", bin)
            from_last = grepl("^!!", bin)

            if(!consecutive){
                if(!is.numeric(x_items)){
                    stop_up("To use the special binning `bin::digit`, the variable `", varname, "` must be numeric. Currently this is not the case (it is of class ", enumerate_items(class(x_items)), " instead).")
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
        stop_up("In `bin`, some values are binned in different bins, it's of course not allowed. The value `", pblm, "` is in ", n_max, " bins.")
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
        stop("Internal bug: Argument `bin` should be equal to `cut::stg` over here -- this is not the case.")
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
            stop_up("In `bin`, the format should be `cut::a]b]` with `a`, `b`, etc, numbers or quartiles/percentiles. Problem: each number must be followed by an open (or closed) square bracket, this is currentlly not the case.")
        }

        cut_points = numeric(n_cuts)
        for(i in seq_along(values)){

            v = trimws(values[i])

            if(is_numeric_in_char(v)){
                cut_points[i] = as.numeric(v)

            } else if(grepl("^p|P|q|Q", v)){
                p = gsub("^.", "", v)
                if(!is_numeric_in_char(p)){
                    stop_up("In `bin`, the format should be `cut::a]b]` with `a`, `b`, etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value `", v, "`, in `", bin, "`, is incorrect.")
                }

                p = as.numeric(p)


                if(grepl("^q|Q", v)){
                    # we transform the quartile into a percentile
                    if(!p %in% c(0:4)){
                        stop_up("In `bin`, the format should be `cut::a]b]` with `a`, `b`, etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value `", v, "`, in `", bin, "`, is incorrect. \n  The quartile must be an integer between 0 and 4.")
                    }

                    p = c(0, 25, 50, 75, 100)[p + 1]
                }

                if(!p %in% 0:100){
                    stop_up("In `bin`, the format should be `cut::a]b]` with `a`, `b`, etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value `", v, "`, in `", bin, "`, is incorrect. \n  The percentile must be an integer between 0 and 100.")
                }

                cut_points[i] = quantile(x, p/100)

            } else {
                stop_up("In `bin`, the format should be `cut::a]b]` with `a`, `b`, etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The value `", v, "`, in `", bin, "`, is incorrect. This is not a percentile nor a number.")
            }
        }

        if(n_cuts > 1 && any(cut_points[1:(n_cuts-1)] > cut_points[1 + 1:(n_cuts-1)])){
            i_pblm = which(cut_points[1:(n_cuts-1)] > cut_points[1 + 1:(n_cuts-1)])[1]
            stop_up("In `bin`, the format should be `cut::a]b]` with `a`, `b`, etc, numbers or quartiles (resp. percentiles) of the form qX (resp. pX) with X a number. \n  The values `a`, `b`, etc should be increasing, but the ", n_th(i_pblm), " value (", cut_points[i_pblm], ") is larger than the ", n_th(i_pblm + 1), " (", cut_points[i_pblm + 1], ").")
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
            message("When binning: in `", bin, "`, the ", enumerate_items(n_th(i_empty)), " bin", plural_len(i_empty, "s.is"), " empty.")
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
            stop("Internal error (=bug): the argument `vcov` should not be missing in fixest_pvalue().")

        } else if(is.null(attr(vcov, "dof.K"))){
            stop("Internal error (=bug): the attribute `dof.K` from `vcov` should not be NULL.")
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

fixest_CI_factor = function(x, level, vcov = NULL, df.t = NULL){

    val = (1 - level) / 2
    val = c(val, 1 - val)

    if(use_t_distr(x)){
        
        if(missing(vcov) && missing(df.t)){
            stop("Internal error (=bug): the arguments `vcov` and `df.t` ",
                 "should not be both missing in fixest_CI_factor().")

        } 
        
        if(!missing(df.t)){
            if(!is.numeric(df.t) && !length(df.t) == 1){
                stop("Internal error (=bug): the arguments `df.t` ",
                     "should be numeric in fixest_CI_factor().")
            }
        } else {
            if(is.null(attr(vcov, "dof.K"))){
                stop("Internal error (=bug): the attribute `dof.K` from `vcov` should not be NULL.")
            }

            # df.t is always an attribute of the vcov
            df.t = attr(vcov, "df.t")
            if(is.null(df.t)){
                df.t = max(nobs(x) - attr(vcov,"dof.K"), 1)
            }
        }
        
        fact = qt(val, df.t)

    } else {
        fact = qnorm(val)
    }

    fact
}


#### ................. ####
#### Small Utilities   ####
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

    check_arg(dict, "NULL named character vector no na", .message = "The argument `dict` must be a dictionnary, ie a named vector (eg dict=c(\"old_name\"=\"new_name\")")

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
            stop("The dictionary `dict` contains several items with the same names, it concerns ",
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

    check_arg(keep, "character vector no na", .message = "The arg. `keep` must be a vector of regular expressions (see help(regex)).")

    res = x

    qui_keep = rep(FALSE, length(x))
    for(var2keep in keep){
        if(grepl("^!%", var2keep)) var2keep = gsub("^!%", "%!", var2keep)

        vect2check = res
        if(grepl("^%", var2keep)){
            var2keep = gsub("^%", "", var2keep)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In keep, the special character `%` cannot be used here.")
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

    check_arg(drop, "character vector no na", .message = "The arg. `drop` must be a vector of regular expressions (see help(regex)). ")

    res = x

    for(var2drop in drop){
        if(grepl("^!%", var2drop)) var2drop = gsub("^!%", "%!", var2drop)

        vect2check = res
        if(grepl("^%", var2drop)){
            var2drop = gsub("^%", "", var2drop)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In drop, the special character `%` cannot be used here.")
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

    check_arg(order, "character vector no na", .message = "The arg. `order` must be a vector of regular expressions (see help(regex)). ")

    res = x

    for(var2order in rev(order)){
        if(grepl("^!%", var2order)) var2order = gsub("^!%", "%!", var2order)

        vect2check = res
        if(grepl("^%", var2order)){
            var2order = gsub("^%", "", var2order)
            vect2check = names(res)
            if(is.null(vect2check)){
                # warning("In order, the special character `%` cannot be used here.")
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

charShorten = function(x, width, keep.digits = FALSE){
	# transforms characters => they can't go beyond a certain width
	# two dots are added to suggest longer character
	# charShorten("bonjour", 5) => "bon.."
	n = nchar(x)

	if(n > width && width > 3){
	    if(keep.digits){
	        trailing_digits = dsb("'\\d+$'X ? x")
	        n_d = nchar(trailing_digits)
	        if(length(n_d) == 1 && n_d > 0){
	            res = dsb(".[`max(width - n_d - 2, 3)`k ? x]...[trailing_digits]")
	            return(res)
	        }
	    }

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
        # => maybe I should be more cautious below???? Sometimes `:` really means stuff like 1:5
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


deparse_short = function(x){
    x_dp = deparse(x)
    if(length(x_dp) > 1){
        x_dp = paste0(x_dp, "...")
    }

    x_dp
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
            msg = paste0("Argument `", arg_name, "` could not be evaluated: ")
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

       res = as.formula(paste(fml_all, collapse = "|"), .GlobalEnv)
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
                stop_up("The operator `^` between variables can be used only in the fixed-effects part of the formula. Otherwise, please use `:` instead.")
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
                stop_up("The operator `^` between variables can be used only in the fixed-effects part of the formula. Otherwise, please use `:` instead.")
            }

            if(is_fml_inside(fml_parts[[n_parts]])){
                iv_txt = fml_split(fml, i = n_parts, text = TRUE)

                if(grepl("\\^[[:alpha:]]", iv_txt)){
                    stop_up("The operator `^` between variables can be used only in the fixed-effects part of the formula. Otherwise, please use `:` instead.")
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
            stop_up("The argument `", arg, "` must start with `r` (for round) or `s` (for significant). Currently it starts with `", d_type,"` which is not valid.\nExample of valid use: digits = `r3`.")
        }

        round = d_type == "r"

        if(!grepl("^[0-9]$", d_value)){
            arg = deparse(substitute(digits))
            stop_up("The argument `", arg, "` must be equal to the character `r` or `s` followed by a single digit. Currently `", digits,"` is not valid.\nExample of valid use: digits = `", d_type, "3`.")
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
                    stop_up("Only one function `", fun, "` can be used at a time.")
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

    substr(x, 1 + n_first, nchar(x) - n_last)
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


is_calling_fun = function(pattern, full_search = FALSE, full_name = FALSE){
    sc_all = sys.calls()
    n_sc = length(sc_all)
    if(n_sc > 2){

        if(full_search){
            fun_all = sapply(tail(sc_all, 13), function(x) deparse(x)[1])

            if(full_name){
                pattern = .dsb("^.[pattern]\\(")
            }

            res = any(grepl(pattern, fun_all))
        } else {
            if(grepl(".fixest", sc_all[[n_sc - 1]][[1]], fixed = TRUE)){
                if(n_sc == 3){
                    return(FALSE)
                }

                sc = sc_all[[n_sc - 3]]
            } else {
                sc = sc_all[[n_sc - 2]]
            }

            fun_name = deparse(sc)[1]
            if(full_name){
                pattern = .dsb("^.[pattern]\\(")
            }

            res = grepl(pattern, fun_name)
        }

        return(res)
    }

    FALSE
}

# x = ~a + (b + c) + d*k
extract_vars_simple = function(x){
    # We don't apply terms!
    # => m * (a + b) is not expanded
    # a bit slow

    check_arg(x, "formula")

    res = c()
    x = x[[length(x)]]

    while(is_operator(x, "+")){
        res = c(res, deparse_long(x[[3]]))
        x = x[[2]]
    }

    res = c(res, deparse_long(x))

    rev(res)
}

# x = "a + (b + c) + d*k"
char_to_vars = function(x){

    check_arg(x, "character scalar")

    if(!grepl("(", x, fixed = TRUE)){
        res = trimws(strsplit(x, "+", fixed = TRUE)[[1]])
        return(res)
    }

    res = c()
    x = str2lang(x)

    while(is_operator(x, "+")){
        res = c(res, deparse_long(x[[3]]))
        x = x[[2]]
    }

    res = c(res, deparse_long(x))

    rev(res)
}

not_too_many_messages = function(key){
    # we don't proc in less than 2 seconds
    # avoid ugly looping issues
    all_times = getOption("fixest_all_timings")
    old_time = all_times[[key]]
    if(is.null(old_time) || as.numeric(Sys.time() - old_time) > 2){
        if(is.null(all_times)){
            all_times = list()
        }
        all_times[[key]] = Sys.time()
        options(fixest_all_timings = all_times)
        return(TRUE)
    }
    
    FALSE
}


#### ............... ####
#### Setters/Getters ####
####

#' Sets/gets whether to display notes in `fixest` estimation functions
#'
#' Sets/gets the default values of whether notes (informing for NA and observations removed) should be displayed in `fixest` estimation functions.
#'
#' @param x A logical. If `FALSE`, then notes are permanently removed.
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

#' Sets/gets the number of threads to use in `fixest` functions
#'
#' Sets/gets the default number of threads to used in `fixest` estimation functions. The default is the maximum number of threads minus two.
#'
#'
#'
#' @param nthreads The number of threads. Can be: a) an integer lower than, or equal to, the maximum number of threads; b) 0: meaning all available threads will be used; c) a number strictly between 0 and 1 which represents the fraction of all threads to use. If missing, the default is to use 50% of all threads.
#' @param save Either a logical or equal to `"reset"`. Default is `FALSE`. If `TRUE` then the value is set permanently at the project level, this means that if you restart R, you will still obtain the previously saved defaults. This is done by writing in the `".Renviron"` file, located in the project's working directory, hence we must have write permission there for this to work, and only works with Rstudio. If equal to "reset", the default at the project level is erased. Since there is writing in a file involved, permission is asked to the user.
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



#' Transforms a character string into a dictionary
#'
#' Transforms a single character string containing a dictionary in a textual format into a proper dictionary, that is a named character vector
#'
#' @param x A character scalar of the form `"variable 1: definition \n variable 2: definition"` etc. Each line of this character must contain at most one definition with, on the left the variable name, and on the right its definition. The separation between the variable and its definition must be a colon followed with a single space (i.e. ": "). You can stack definitions within a single line by making use of a semi colon: `"var1: def; var2: def"`. White spaces on the left and right are ignored. You can add commented lines with a `"#"`. Non-empty, non-commented lines that don't have the proper format witll raise an error.
#'
#' @details
#' This function is mostly used in combination with [`setFixest_dict`] to set the dictionary to be used in the function [`etable`].
#'
#' @return
#' It returns a named character vector.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' [`etable`], [`setFixest_dict`]
#'
#' @examples
#'
#' x = "# Main vars
#'      mpg: Miles per gallon
#'      hp: Horsepower
#'
#'      # Categorical variables
#'      cyl: Number of cylinders; vs: Engine"
#'
#' as.dict(x)
#'
#'
#'
as.dict = function(x){
    check_arg(x, "character scalar")

    text_split = strsplit(x, ";|\n")[[1]]

    text_clean = trimws(text_split)
    text_clean = text_clean[nchar(text_clean) > 0 & !grepl("^#", text_clean)]

    if(any(!grepl(": ", text_clean))){
        stop("All definitions must be of the form 'variable: definition', the colon and space are indispensible. Currently the following line is not valid:\n", text_clean[!grepl(": ", text_clean)][1])
    }

    dict_names = sapply(strsplit(text_clean, ": ", fixed = TRUE), `[[`, 1)
    if(anyDuplicated(dict_names)){
        tb = sort(table(dict_names), decreasing = TRUE)[1]
        stop("The names contain duplicate entries: this is not allowed. The value ",
             names(tb), " appears ", dreamerr::n_times(tb), ".")
    }

    values = substr(text_clean, nchar(dict_names) + 3, nchar(text_clean))

    setNames(values, dict_names)
}

#' Sets/gets the dictionary relabeling the variables
#'
#' Sets/gets the default dictionary used in the function [`etable`], [`did_means`] and [`coefplot`]. The dictionaries are used to relabel variables (usually towards a fancier, more explicit formatting) when exporting them into a Latex table or displaying in graphs. By setting the dictionary with `setFixest_dict`, you can avoid providing the argument `dict`.
#'
#'
#' @param dict A named character vector or a character scalar. E.g. to change my variable named "a" and "b" to (resp.) "$log(a)$" and "$bonus^3$", then use `dict = c(a="$log(a)$", b3="$bonus^3$")`. Alternatively you can feed a character scalar containing the dictionary in the form `"variable 1: definition \n variable 2: definition"`. In that case the function [`as.dict`] will be applied to get a proper dictionary. This dictionary is used in Latex tables or in graphs by the function [`coefplot`]. If you want to separate Latex rendering from rendering in graphs, use an ampersand first to make the variable specific to `coefplot`.
#' @param ... You can add arguments of the form: `variable_name = "Definition"`. This is an alternative to using a named vector in the argument `dict`.
#' @param reset Logical, default is `FALSE`. If `TRUE`, then the dictionary is reset. Note that the default dictionary always relabels the variable "(Intercept)" in to "Constant". To overwrite it, you need to add "(Intercept)" explicitly in your dictionary.
#'
#' @details
#' By default the dictionary only grows. This means that successive calls with not erase the previous definitions unless the argument `reset` has been set to `TRUE`.
#'
#' The default dictionary is equivalent to having `setFixest_dict("(Intercept)" = "Constant")`. To change this default, you need to provide a new definition to `"(Intercept)"` explicitly.
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
#' etable(est, dict = c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#'
#' # If you export many tables, it can be more convenient to use setFixest_dict:
#' setFixest_dict(c("log(Euros)"="Euros (ln)", Origin="Country of Origin"))
#' etable(est) # variables are properly relabeled
#'
#' # The dictionary only 'grows'
#' # Here you get the previous two variables + the new one that are relabeled
#' # Btw you set the dictionary directly using the argument names:
#' setFixest_dict(Destination = "Country of Destination")
#' etable(est)
#'
#' # Another way to set a dictionary: with a character string:
#' # See the help page of as.dict
#' dict = "log(dist_km): Distance (ln); Product: Type of Good"
#' setFixest_dict(dict)
#' etable(est)
#'
#' # And now we reset:
#' setFixest_dict(reset = TRUE)
#' etable(est)
#'
setFixest_dict = function(dict = NULL, ..., reset = FALSE){

    check_arg(dict, "NULL named character vector no na | character scalar")

    if(is.null(names(dict)) && !is.null(dict)){
        dict = error_sender(as.dict(dict), "In setFixest_dict, problem when coercing the dictionay with as.dict().")
    }

    check_arg(..., "dotnames character scalar",
              .message = .dsb("In '...', each argument must be named. ",
                              "The argument name corresponds to the variable to be renamed while",
                              " the value must be a character scalar (how the variable should be renamed)."))

    dots = list(...)
    dict = as.list(dict)
    dict[names(dots)] = dots

    if(reset){
        core_dict = list("(Intercept)" = "Constant")
    } else {
        core_dict = getOption("fixest_dict")
        if(is.null(core_dict)) core_dict = list()
    }

    core_dict[names(dict)] = dict

	options("fixest_dict" = unlist(core_dict))
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
#' You can set formula macros globally with `setFixest_fml`. These macros can then be used in `fixest` estimations or when using the function [`xpd`][fixest::setFixest_fml].
#'
#' @inherit xpd examples
#'
#' @param ... Definition of the macro variables. Each argument name corresponds to the name of the macro variable. It is required that each macro variable name starts with two dots (e.g. `..ctrl`). The value of each argument must be a one-sided formula or a character vector, it is the definition of the macro variable. Example of a valid call: `setFixest_fml(..ctrl = ~ var1 + var2)`. In the function `xpd`, the default macro variables are taken from `getFixest_fml`, any variable in `...` will replace these values. You can enclose values in `.[]`, if so they will be evaluated from the current environment. For example `..ctrl = ~ x.[1:2] + .[z]` will lead to `~x1 + x2 + var` if `z` is equal to `"var"`.
#' @param reset A logical scalar, defaults to `FALSE`. If `TRUE`, all macro variables are first reset (i.e. deleted).
#'
#' @details
#' In `xpd`, the default macro variables are taken from `getFixest_fml`. Any value in the `...` argument of `xpd` will replace these default values.
#'
#' The definitions of the macro variables will replace in verbatim the macro variables. Therefore, you can include multipart formulas if you wish but then beware of the order the the macros variable in the formula. For example, using the airquality data, say you want to set as controls the variable `Temp` and `Day` fixed-effects, you can do `setFixest_fml(..ctrl = ~Temp | Day)`, but then `feols(Ozone ~ Wind + ..ctrl, airquality)` will be quite different from `feols(Ozone ~ ..ctrl + Wind, airquality)`, so beware!
#'
#' @return
#' The function `getFixest_fml()` returns a list of character strings, the names corresponding to the macro variable names, the character strings corresponding to their definition.
#'
#' @seealso
#' [`xpd`] to make use of formula macros.
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
#' @param reset Logical scalar, default is `FALSE`. Whether to reset all values.
#'
#' @return
#' The function `getFixest_estimation` returns the currently set global defaults.
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
setFixest_estimation = function(data = NULL, panel.id = NULL, fixef.rm = "perfect",
                                fixef.tol = 1e-6, fixef.iter = 10000, collin.tol = 1e-10,
                                lean = FALSE, verbose = 0, warn = TRUE, combine.quick = NULL,
                                demeaned = FALSE, mem.clean = FALSE, glm.iter = 25,
                                glm.tol = 1e-8, reset = FALSE){

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
#' `trade` is a data frame with 38,325 observations and 6 variables named `Destination`, `Origin`, `Product`, `Year`, `dist_km` and `Euros`.
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

































































