#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Apr 19 17:23:31 2021
# ~: Various fit statistics
#----------------------------------------------#


####
#### stats -- User level ####
####



#' Print method for fit statistics of fixest estimations
#'
#' Displays a brief summary of selected fit statistics from the function \code{\link[fixest]{fitstat}}.
#'
#' @inherit fitstat examples
#'
#' @param x An object resulting from the \code{\link[fixest]{fitstat}} function.
#' @param na.rm Logical, default is \code{FALSE}. If \code{TRUE}, the statistics that are missing are not displayed.
#' @param ... Not currently used.
#'
#'
print.fixest_fitstat = function(x, na.rm = FALSE, ...){

    dots = list(...)

    # glue = function(x) gsub("( +)(,)", "\\2\\1", paste(x[nchar(x) > 0], collapse = ", "))
    glue = function(x) paste(x[nchar(x) > 0], collapse = ", ")

    all_types = fitstat(give_types = TRUE)
    dict_type = all_types$R_alias

    if(na.rm){
        IS_NA = sapply(x, function(z) identical(z, NA))
        x = x[!IS_NA]
        if(length(x) == 0) return(invisible(NULL))
    }

    if(isTRUE(dots$group.solo)){
        # Later integration in print
        solo = lengths(x) == 1
        x_solo = unlist(x[solo])
        n_solo = length(x_solo)

        if(n_solo > 1){
            x_solo_names = names(x[solo])
            x_solo = unlist(x[solo])

            x = x[!solo]

            nr = ceiling(n_solo / 2)
            m = matrix(c(x_solo, NA)[1:(2*nr)], nr, 2, byrow = TRUE)
            m_alias = matrix(c(unlist(dict_type[x_solo_names]), "")[1:(2*nr)], nr, 2, byrow = TRUE)

            if(nr > 1 && anyNA(m[nr, ])){
                m[nr, ] = m[nr, 2:1]
                m_alias[nr, ] = m_alias[nr, 2:1]
            }

            cols_new = list()
            for(i in 1:2){
                col = m[,i]
                gt1 = abs(col) > 1 & !is.na(col)
                col_char = numberFormatNormal(col)
                col_char[gt1] = sfill(col_char[gt1], anchor = ".", na = "")
                col_char[!gt1] = sfill(col_char[!gt1], anchor = ".", na = "")
                col_char = sfill(col_char, right = TRUE)

                # The aliases
                col_char = paste0(sfill(m_alias[, i]), ": ", col_char)
                if(anyNA(col)){
                    col_char[is.na(col)] = sfill("", nchar(col_char[1]))
                }

                cols_new[[i]] = col_char
            }
            mm = do.call(cbind, cols_new)

            cat(apply(mm, 1, paste, collapse = "   "), sep = "\n")
        }

        if(length(x) == 0) return(invisible(NULL))
    }

    res = c()

    # Formatting the stats and pvalues
    is_stat = sapply(x, function(z) "stat" %in% names(z))
    is_pval = sapply(x, function(z) "p" %in% names(z))

    if(any(is_stat)){
        stat_all = sapply(x[is_stat], function(z) z$stat)
        stat_all = paste0("stat = ", sfill(sfill(numberFormatNormal(stat_all), anchor = "."), right = TRUE))
        for(i in seq_along(stat_all)){
            j = which(is_stat)[i]
            x[[j]]$stat = stat_all[i]
        }
    }

    if(any(is_pval)){
        pval_all = sapply(x[is_pval], function(z) z$p)
        low = pval_all < 2.2e-16 & !is.na(pval_all)
        pval_all[low] = 2.2e-16
        pval_all = paste0("p ", ifelse(low, "< ", "= "), sfill(numberFormatNormal(pval_all), right = TRUE))

        for(i in seq_along(pval_all)){
            j = which(is_pval)[i]
            x[[j]]$p = pval_all[i]
        }
    }

    # formatting for multiple endo regs
    qui_right = rep(FALSE, length(x))

    for(i in seq_along(x)){
        type = names(x)[i]

        v = x[[i]]

        if(grepl("::", type, fixed = TRUE)){
            dict = getFixest_dict()
            rename_fun = function(x) paste0(dict_type[x[1]], ", ", paste(sapply(x[-1], dict_apply, dict = dict), collapse = "::"))
            test_name = rename_fun(strsplit(type, "::")[[1]])
            qui_right[i] = TRUE
        } else {
            if(type %in% names(dict_type)){
                test_name = dict_type[type]
            } else if(!grepl(".", type, fixed = TRUE)) {
                warning("Current type '", type, "' has no alias.")
                test_name = type
            } else {
                type_split = strsplit(type, ".", fixed = TRUE)[[1]]
                test_name = paste0(dict_type[type_split[1]], ", ", dict_type[type_split[2]])
            }

        }

        if(length(v) == 1){
            # Basic display
            res[i] = paste0(test_name, "! ", numberFormatNormal(v))
        } else {

            stat_line = p_line = dof_line = vcov_line = ""
            if(!is.null(v$stat)) stat_line = v$stat
            if(!is.null(v$p)) p_line = v$p
            if(!is.null(v$df)) dof_line = paste0("on ", numberFormatNormal(v$df), " DoF")
            if(!is.null(v$df1)) dof_line = paste0("on ", numberFormatNormal(v$df1), " and ", numberFormatNormal(v$df2), " DoF")
            if(!is.null(v$vcov)) vcov_line = paste0("VCOV: ", v$vcov)

            res[i] = paste0(test_name, "! ", glue(c(stat_line, p_line, dof_line, vcov_line)), ".")
        }
    }

    if(any(qui_right)){
        res2fix = strsplit(res[qui_right], "!", fixed = TRUE)
        left = sapply(res2fix, function(x) x[1])
        right = sapply(res2fix, function(x) x[2])
        res_fixed = paste0(sfill(left, right = TRUE), "!", right)
        res[qui_right] = res_fixed
    }

    cat(gsub("!", ":", sfill(res, anchor = "!"), fixed = TRUE), sep = "\n")

}



#' Register custom fit statistics
#'
#' Enables the registration of custom fi statistics that can be easily summoned with the function \code{\link[fixest]{fitstat}}.
#'
#' @inherit fitstat examples
#' @inherit fitstat seealso
#'
#' @param type A character scalar giving the type-name.
#' @param fun A function to be applied to a \code{fixest} estimation. It must return either a scalar, either a list. Note that for the print method to work correctly, the names of the items of the list must be one of: \code{stat}, \code{p}, \code{df}, \code{df1}, \code{df2}, \code{vcov}. Only the print method is affected by this.
#' @param alias An alias to be used in lieu of the type name in the display methods (ie when used in the function \code{\link[fixest]{print.fixest_fitstat}} or \code{\link[fixest]{etable}}).
#'
#'
fitstat_register = function(type, fun, alias){

    check_arg(type, "character scalar mbt")
    check_arg(fun, "function mbt")
    check_arg(alias, "NULL character scalar")

    # We check the type is not conflicting
    existing_types = fitstat(give_types = TRUE)$types

    opts = getOption("fixest_fitstat_user")

    if(type %in% setdiff(existing_types, names(opts))){
        stop("The type name '", type, "' is the same as one built-in type. Please choose another one.")
    }

    if(missnull(alias)){
        alias = type
    }

    res = list(fun = fun, alias = alias)

    opts[[type]] = res

    options(fixest_fitstat_user = opts)

    invisible(NULL)
}

#' Computes fit statistics of fixest objects
#'
#' Computes various fit statistics for \code{fixest} estimations.
#'
#' @param x A \code{fixest} estimation.
#' @param type Character vector or one sided formula. The type of fit statistic to be computed. The classic ones are: n, rmse, r2, pr2, f, wald, ivf, ivwald. You have the full list in the details section or use \code{show_types = TRUE}. Further, you can register your own types with \code{\link[fixest]{fitstat_register}}.
#' @param simplify Logical, default is \code{FALSE}. By default a list is returned whose names are the selected types. If \code{simplify = TRUE} and only one type is selected, then the element is directly returned (ie will not be nested in a list).
#' @param verbose Logical, default is \code{TRUE}. If \code{TRUE}, an object of class \code{fixest_fitstat} is returned (so its associated print method will be triggered). If \code{FALSE} a simple list is returned instead.
#' @param show_types Logical, default is \code{FALSE}. If \code{TRUE}, only prompts all available types.
#' @param ... Other elements to be passed to other methods and may be used to compute the statistics (for example you can pass on arguments to compute the VCOV when using \code{type = "g"} or \code{type = "wald"}.).
#'
#' @section Registering your own types:
#'
#' You can register custom fit statistics with the function \code{fitstat_register}.
#'
#' @section Available types:
#'
#' The types are case sensitive, please use lower case only. The types available are:
#'
#' \itemize{
#' \item{\code{n}, \code{ll}, \code{aic}, \code{bic}, \code{rmse}: }{The number of observations, the log-likelihood, the AIC, the BIC and the root mean squared error, respectively.}
#' \item{\code{my}: }{Mean of the dependent variable.}
#' \item{\code{g}: }{The degrees of freedom used to compute the t-test (it influences the p-values of the coefficients). When the VCOV is clustered, this value is equal to the minimum cluster size, otherwise, it is equal to the sample size minus the number of variables.}
#' \item{\code{r2}, \code{ar2}, \code{wr2}, \code{awr2}, \code{pr2}, \code{apr2}, \code{wpr2}, \code{awpr2}: }{All r2 that can be obtained with the function \code{\link[fixest]{r2}}. The \code{a} stands for 'adjusted', the \code{w} for 'within' and the \code{p} for 'pseudo'. Note that the order of the letters \code{a}, \code{w} and \code{p} does not matter.}
#' \item{\code{theta}: }{The over-dispersion parameter in Negative Binomial models. Low values mean high overdispersion. }
#' \item{\code{f}, \code{wf}: }{The F-tests of nullity of the coefficients. The \code{w} stands for 'within'. These types return the following values: \code{stat}, \code{p}, \code{df1} and \code{df2}. If you want to display only one of these, use their name after a dot: e.g. \code{f.stat} will give the statistic of the F-test, or \code{wf.p} will give the p-values of the F-test on the projected model (i.e. projected onto the fixed-effects).}
#' \item{\code{wald}: }{Wald test of joint nullity of the coefficients. This test always excludes the intercept and the fixed-effects. These type returns the following values: \code{stat}, \code{p}, \code{df1}, \code{df2} and \code{vcov}. The element \code{vcov} reports the way the VCOV matrix was computed since it directly influences this statistic.}
#' \item{\code{ivf}, \code{ivf1}, \code{ivf2}, \code{ivfall}: }{These statistics are specific to IV estimations. They report either the IV F-test (namely the Cragg-Donald F statistic in the presence of only one endogenous regressor) of the first stage (\code{ivf} or \code{ivf1}), of the second stage (\code{ivf2}) or of both (\code{ivfall}). The F-test of the first stage is commonly named weak instrument test. The value of \code{ivfall} is only useful in \code{\link[fixest]{etable}} when both the 1st and 2nd stages are displayed (it leads to the 1st stage F-test(s) to be displayed on the 1st stage estimation(s), and the 2nd stage one on the 2nd stage estimation -- otherwise, \code{ivf1} would also be displayed on the 2nd stage estimation). These types return the following values: \code{stat}, \code{p}, \code{df1} and \code{df2}.}
#' \item{\code{ivwald}, \code{ivwald1}, \code{ivwald2}, \code{ivwaldall}: }{These statistics are specific to IV estimations. They report either the IV Wald-test of the first stage (\code{ivwald} or \code{ivwald1}), of the second stage (\code{ivwald2}) or of both (\code{ivwaldall}). The Wald-test of the first stage is commonly named weak instrument test. Note that if the estimation was done with a robust VCOV and there is only one endogenous regressor, this is equivalent to the Kleibergen-Paap statistic. The value of \code{ivwaldall} is only useful in \code{\link[fixest]{etable}} when both the 1st and 2nd stages are displayed (it leads to the 1st stage Wald-test(s) to be displayed on the 1st stage estimation(s), and the 2nd stage one on the 2nd stage estimation -- otherwise, \code{ivwald1} would also be displayed on the 2nd stage estimation). These types return the following values: \code{stat}, \code{p}, \code{df1}, \code{df2}, and \code{vcov}.}
#' \item{\code{cd}: }{The Cragg-Donald test for weak instruments.}
#' \item{\code{kpr}: }{The Kleibergen-Paap test for weak instruments.}
#' \item{\code{wh}: }{This statistic is specific to IV estimations. Wu-Hausman endogeneity test. H0 is the absence of endogeneity of the instrumented variables. It returns the following values: \code{stat}, \code{p}, \code{df1}, \code{df2}.}
#' \item{\code{sargan}: }{Sargan test of overidentifying restrictions. H0: the instruments are not correlated with the second stage residuals. It returns the following values: \code{stat}, \code{p}, \code{df}.}
#' \item{\code{lr}, \code{wlr}: }{Likelihood ratio and within likelihood ratio tests. It returns the following elements: \code{stat}, \code{p}, \code{df}. Concerning the within-LR test, note that, contrary to estimations with \code{femlm} or \code{feNmlm}, estimations with \code{feglm}/\code{fepois} need to estimate the model with fixed-effects only which may prove time-consuming (depending on your model). Bottom line, if you really need the within-LR and estimate a Poisson model, use \code{femlm} instead of \code{fepois} (the former uses direct ML maximization for which the only FEs model is a by product).}
#' }
#'
#'
#'
#'
#' @return
#' By default an object of class \code{fixest_fitstat} is returned. Using \code{verbose = FALSE} returns a simple a list. Finally, if only one type is selected, \code{simplify = TRUE} leads to the selected type to be returned.
#'
#' @examples
#'
#' data(trade)
#' gravity = feols(log(Euros) ~ log(dist_km) | Destination + Origin, trade)
#'
#' # Extracting the 'working' number of observations used to compute the pvalues
#' fitstat(gravity, "g", simplify = TRUE)
#'
#' # Some fit statistics
#' fitstat(gravity, ~ rmse + r2 + wald + wf)
#'
#' # You can use them in etable
#' etable(gravity, fitstat = ~ rmse + r2 + wald + wf)
#'
#' # For wald and wf, you could show the pvalue instead:
#' etable(gravity, fitstat = ~ rmse + r2 + wald.p + wf.p)
#'
#' # Now let's display some statistics that are not built-in
#' # => we use fitstat_register to create them
#'
#' # We need: a) type name, b) the function to be applied
#' #          c) (optional) an alias
#'
#' fitstat_register("tstand", function(x) tstat(x, se = "stand")[1], "t-stat (regular)")
#' fitstat_register("thc", function(x) tstat(x, se = "heter")[1], "t-stat (HC1)")
#' fitstat_register("t1w", function(x) tstat(x, se = "clus")[1], "t-stat (clustered)")
#' fitstat_register("t2w", function(x) tstat(x, se = "twow")[1], "t-stat (2-way)")
#'
#' # Now we can use these keywords in fitstat:
#' etable(gravity, fitstat = ~ . + tstand + thc + t1w + t2w)
#'
#' # Note that the custom stats we created are can easily lead
#' # to errors, but that's another story!
#'
#'
fitstat = function(x, type, simplify = FALSE, verbose = TRUE, show_types = FALSE, ...){

    r2_types = c("sq.cor", "cor2", "r2", "ar2", "pr2", "apr2", "par2", "wr2", "war2", "awr2", "wpr2", "pwr2", "wapr2", "wpar2", "awpr2", "apwr2", "pawr2", "pwar2")

    opts = getOption("fixest_fitstat_user")
    user_types = names(opts)

    dots = list(...)

    if(isTRUE(dots$give_types) || isTRUE(show_types)){

        # Compound types => types yielding several values

        # F-test etc
        comp_types = c("f", "wf", "ivf", "ivf1", "ivf2", "ivfall", "wald", "ivwald", "ivwald1", "ivwald2", "ivwaldall", "wh", "sargan", "lr", "wlr", "kpr", "cd")
        full_comp_types = paste(comp_types, rep(c("stat", "p"), each = length(comp_types)), sep = ".")

        comp_alias = c(f = "F-test", wf = "F-test (projected)", ivfall = "F-test (IV only)", ivf1 = "F-test (1st stage)", ivf2 = "F-test (2nd stage)", wald = "Wald (joint nullity)", ivwaldall = "Wald (IV only)", ivwald1 = "Wald (1st stage)", ivwald2 = "Wald (2nd stage)", wh = "Wu-Hausman", sargan = "Sargan", lr = "LR", wlr = "LR (within)", df1 = "DoF (first)", df2 = "DoF (second)", df = "DoF", kpr = "Kleibergen-Paap", cd = "Cragg-Donald")
        my_names = paste(names(comp_alias), rep(c("stat", "p"), each = length(comp_alias)), sep = ".")
        full_comp_alias = setNames(paste0(comp_alias, ", ", rep(c("stat.", "p-value"), each = length(comp_alias))), my_names)

        # "Regular" types
        valid_types = c("n", "ll", "aic", "bic", "rmse", "g", "my", "theta", r2_types, comp_types, full_comp_types, user_types)

        user_alias = sapply(opts, function(x) x$alias)

        tex_alias = c(n = "Observations", ll = "Log-Likelihood", aic = "AIC", bic = "BIC", my = "Dependent variable mean", g = "Size of the 'effective' sample", rmse = "RMSE", theta = "Over-dispersion", sq.cor = "Squared Correlation", cor2 = "Squared Correlation", r2="R$^2$", ar2="Adjusted R$^2$", pr2="Pseudo R$^2$", apr2="Adjusted Pseudo R$^2$", wr2="Within R$^2$", war2="Within Adjusted R$^2$", wpr2="Within Pseudo R$^2$", wapr2="Whithin Adjusted Pseudo R$^2$", comp_alias, full_comp_alias, user_alias)

        R_alias = c(n = "Observations", ll = "Log-Likelihood", aic = "AIC", bic = "BIC", my = "Dep. Var. mean", g = "G", rmse = "RMSE", theta = "Over-dispersion", sq.cor = "Squared Cor.", cor2 = "Squared Cor.", r2="R2", ar2="Adj. R2", pr2="Pseudo R2", apr2="Adj. Pseudo R2", wr2="Within R2", war2="Within Adj. R2", wpr2="Within Pseudo R2", wapr2="Whithin Adj. Pseudo R2", comp_alias, full_comp_alias, user_alias)

        # add r2 type alias
        type_alias = c(ivf = "ivf1", ivf.stat = "ivf1.stat", ivf.p = "ivf1.p", ivwald = "ivwald1", ivwald.stat = "ivwald1", ivwald.p = "ivwald1.p", par2 = "apr2", awr2 = "war2", pwr2 = "wpr2", wpar2 = "wapr2", pwar2 = "wapr2", pawr2 = "wapr2", apwr2 = "wapr2", awpr2 = "wapr2", sq.cor = "cor2")

        if(show_types){
            cat("Available types:\n", paste(valid_types, collapse = ", "), "\n")
            return(invisible(NULL))
        }

        res = list(types = valid_types, tex_alias = tex_alias, R_alias = R_alias, type_alias = type_alias)
        return(res)
    }

    check_arg(x, "class(fixest) mbt")
    check_arg_plus(type, "character vector no na | os formula")
    check_arg(simplify, verbose, "logical scalar")

    if("formula" %in% class(type)){
        type = gsub(" ", "", strsplit(deparse_long(type[[2]]), "+", fixed = TRUE)[[1]])
    }
    type = tolower(type)

    # To update
    if(!isTRUE(x$summary) || any(c("se", "cluster", "dof") %in% names(dots))){
        x = summary(x, ...)
    }

    IS_ETABLE = isTRUE(dots$etable)
    set_value = function(vec, value){
        if(length(vec) == 1) return(vec)

        if(value != ""){

            if(!value %in% names(vec)){
                stop_up("'", value, "' is not a valid component of the statistic '", root, "'. Only the following are valid: ", enumerate_items(names(vec), "quote"), ".")
            }

            return(vec[[value]])
        } else if(IS_ETABLE){
            return(vec[1])
        } else {
            return(vec)
        }
    }

    res_all = list()
    type_all = type

    for(i in seq_along(type_all)){
        type = type_all[i]

        # Big if
        if(type == "n"){
            res_all[[type]] = x$nobs

        } else if(type == "ll"){
            res_all[[type]] = logLik(x)

        } else if(type == "aic"){
            res_all[[type]] = AIC(x)

        } else if(type == "bic"){
            res_all[[type]] = BIC(x)

        } else if(type == "rmse"){
            res_all[[type]] = if(!is.null(x$ssr)) sqrt(x$ssr / x$nobs) else sqrt(cpp_ssq(resid(x)) / x$nobs)

        } else if(type == "g"){
            my_vcov = x$cov.scaled

            G = attr(my_vcov, "G")
            if(is.null(G)) G = x$nobs - x$nparams

            res_all[[type]] = G

        } else if(type == "my"){
            res_all[[type]] = mean(model.matrix(x, type = "lhs"))

        } else if(type == "theta"){
            isNegbin = x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && x$family == "negbin")
            if(isNegbin){
                theta = coef(x)[".theta"]
                names(theta) = "Overdispersion"
            } else {
                theta = NA
            }

            res_all[[type]] = theta

        } else if(type %in% r2_types){
            res_all[[type]] = r2(x, type)

        } else {

            #
            # Types for which root.value is allowed
            #

            if(grepl(".", type, fixed = TRUE)){
                root = gsub("\\..+", "", type)
                value = gsub(".+\\.", "", type)
            } else {
                root = type
                value = ""
            }

            # We need to normalize the types
            from_to = c("ivwald" = "ivwald1", "ivf" = "ivf1")
            if(root %in% names(from_to)){
                type = gsub(root, from_to[root], type, fixed = TRUE)
                root = from_to[root]
            }

            if(root == "f"){
                if(!is.null(x$ssr)){
                    df1 = degrees_freedom(x, "k") - 1
                    df2 = degrees_freedom(x, "t")

                    if(isTRUE(x$iv) && x$iv_stage == 2){
                        # We need to compute the SSR
                        w = 1
                        if(!is.null(x$weights)) w = x$weights
                        ssr = cpp_ssq(x$iv_residuals, w)
                    } else {
                        ssr = x$ssr
                    }

                    stat = ((x$ssr_null - ssr) / df1) / (ssr / df2)
                    p = pf(stat, df1, df2, lower.tail = FALSE)
                    vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }
            } else if(root == "wf"){
                df1 = length(x$coefficients)

                if(!is.null(x$ssr_fe_only) && df1 > 0){
                    df2 = x$nobs - x$nparams
                    stat = ((x$ssr_fe_only - x$ssr) / df1) / (x$ssr / df2)
                    p = pf(stat, df1, df2, lower.tail = FALSE)
                    vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }


            } else if(root == "ivfall"){
                if(isTRUE(x$iv)){
                    if(x$iv_stage == 1){
                        df1 = degrees_freedom(x, vars = x$iv_inst_names_xpd)
                        df2 = degrees_freedom(x, "resid")

                        stat = ((x$ssr_no_inst - x$ssr) / df1) / (x$ssr / df2)
                        p = pf(stat, df1, df2, lower.tail = FALSE)
                        vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                        res_all[[type]] = set_value(vec, value)

                    } else {
                        # f stat for the second stage

                        df1 = degrees_freedom(x, vars = x$iv_endo_names_fit)
                        df2 = degrees_freedom(x, "resid")

                        w = 1
                        if(!is.null(x$weights)) w = x$weights
                        ssr = cpp_ssq(x$iv_residuals, w)

                        stat = ((x$ssr_no_endo - ssr) / df1) / (ssr / df2)
                        p = pf(stat, df1, df2, lower.tail = FALSE)
                        vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                        res_all[[type]] = set_value(vec, value)
                    }

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "ivf1"){

                if(isTRUE(x$iv)){
                    df1 = degrees_freedom(x, vars = x$iv_inst_names_xpd, stage = 1)
                    df2 = degrees_freedom(x, "resid", stage = 1)

                    if(x$iv_stage == 1){

                        stat = ((x$ssr_no_inst - x$ssr) / df1) / (x$ssr / df2)
                        p = pf(stat, df1, df2, lower.tail = FALSE)
                        vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                        res_all[[type]] = set_value(vec, value)

                    } else {
                        x_first = x$iv_first_stage

                        for(endo in names(x_first)){
                            stat = ((x_first[[endo]]$ssr_no_inst - x_first[[endo]]$ssr) / df1) / (x_first[[endo]]$ssr / df2)
                            p = pf(stat, df1, df2, lower.tail = FALSE)
                            vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                            res_all[[paste0(type, "::", endo)]] = set_value(vec, value)
                        }
                    }

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "ivf2"){
                if(isTRUE(x$iv) && x$iv_stage == 2){
                    # f stat for the second stage

                    df1 = degrees_freedom(x, vars = x$iv_endo_names_fit)
                    df2 = degrees_freedom(x, "resid")

                    w = 1
                    if(!is.null(x$weights)) w = x$weights
                    ssr = cpp_ssq(x$iv_residuals, w)

                    stat = ((x$ssr_no_endo - ssr) / df1) / (ssr / df2)
                    p = pf(stat, df1, df2, lower.tail = FALSE)
                    vec = list(stat = stat, p = p, df1 = df1, df2 = df2)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "kpr"){
                if(isTRUE(x$iv) && x$iv_stage == 2){
                    # The KP rank test is computed in a specific function

                    vec = kp_stat(x)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "cd"){
                if(isTRUE(x$iv) && x$iv_stage == 2){
                    # The KP rank test is computed in a specific function

                    vec = cd_stat(x)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "wald"){
                # Joint nullity of the coefficients
                # if FE => on the projected model only

                qui = !names(x$coefficients) %in% "(Intercept)"
                my_coef = x$coefficients[qui]
                df1 = length(my_coef)

                if(df1 > 0){
                    df2 = degrees_freedom(x, "resid")

                    # The VCOV is always full rank in here
                    stat = .wald(x, names(my_coef))
                    p = pf(stat, df1, df2, lower.tail = FALSE)
                    vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = attr(x$cov.scaled, "type"))
                    res_all[[type]] = set_value(vec, value)

                } else {
                    # Only fixef
                    res_all[[type]] = NA
                }

            } else if(root == "ivwaldall"){

                if(isTRUE(x$iv)){
                    if(x$iv_stage == 1){

                        df1 = degrees_freedom(x, vars = x$iv_inst_names_xpd)
                        df2 = degrees_freedom(x, "resid")

                        stat = .wald(x, x$iv_inst_names_xpd)

                        p = pf(stat, df1, df2, lower.tail = FALSE)
                        vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = attr(x$cov.scaled, "type"))
                        res_all[[type]] = set_value(vec, value)

                    } else {
                        # wald stat for the second stage

                        df1 = degrees_freedom(x, vars = x$iv_endo_names_fit)
                        df2 = degrees_freedom(x, "resid")

                        stat = .wald(x, x$iv_endo_names_fit)

                        p = pf(stat, df1, df2, lower.tail = FALSE)
                        vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = attr(x$cov.scaled, "type"))
                        res_all[[type]] = set_value(vec, value)
                    }

                } else {
                    res_all[[type]] = NA
                }

            } else if(root %in% "ivwald1"){

                if(isTRUE(x$iv)){
                    df1 = degrees_freedom(x, vars = x$iv_inst_names_xpd, stage = 1)
                    df2 = degrees_freedom(x, "resid", stage = 1)


                    if(x$iv_stage == 1){

                        stat = .wald(x, x$iv_inst_names_xpd)

                        p = pf(stat, df1, df2, lower.tail = FALSE)
                        vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = attr(x$cov.scaled, "type"))
                        res_all[[type]] = set_value(vec, value)

                    } else {
                        x_first = x$iv_first_stage
                        inst = x$iv_inst_names_xpd

                        for(endo in names(x_first)){

                            my_x_first = x_first[[endo]]

                            if(is.null(my_x_first$cov.scaled)){
                                # We compute the VCOV like for the second stage
                                if(is.null(x$se_info) || !is.null(dots$se) || !is.null(dots$cluster)){
                                    my_x_first = summary(my_x_first, ...)
                                } else {
                                    my_x_first = summary(my_x_first, se = x$se_info$se, cluster = x$se_info$cluster, dof = x$se_info$dof, ...)
                                }
                            }

                            stat = .wald(my_x_first, inst)

                            p = pf(stat, df1, df2, lower.tail = FALSE)
                            vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = attr(x$cov.scaled, "type"))
                            res_all[[paste0(type, "::", endo)]] = set_value(vec, value)
                        }
                    }

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "ivwald2"){
                if(isTRUE(x$iv) && x$iv_stage == 2){
                    # wald stat for the second stage

                    df1 = degrees_freedom(x, vars = x$iv_endo_names_fit)
                    df2 = degrees_freedom(x, "resid")

                    stat = .wald(x, x$iv_endo_names_fit)

                    p = pf(stat, df1, df2, lower.tail = FALSE)
                    vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = attr(x$cov.scaled, "type"))
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "wh"){
                if(isTRUE(x$iv) && x$iv_stage == 2){
                    # Wu Hausman stat for the second stage
                    vec = x$iv_wh
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "sargan"){
                if(!is.null(x$iv_sargan)){
                    # Wu Hausman stat for the second stage
                    vec = x$iv_sargan
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "lr"){
                if(!is.null(x$ll_null)){

                    stat = 2 * (x$loglik - x$ll_null)
                    df = x$nparams
                    p = pchisq(stat, df, lower.tail = FALSE)

                    vec = list(stat = stat, p = p, df = df)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root == "wlr"){
                if(!is.null(x$ll_fe_only) || (x$method_type == "feglm" && !is.null(x$fixef_id))){

                    if(x$method_type == "feglm"){
                        # estimation of the FE only model

                        newdata = cbind(data.frame(y = x$y), as.data.frame(x$fixef_id))
                        if(!is.null(x$fixef_terms)){
                            newdata = cbind(newdata, as.data.frame(x$slope_variables))
                        }
                        new_fml = merge_fml(y ~ 1, x$fml_all$fixef)
                        res_fe = feglm(fml = new_fml, data = newdata, glm.tol = 1e-2, fixef.tol = 1e-3, family = x$family$family, weights = x$weights, offset = x$offset)

                        ll_fe_only = logLik(res_fe)
                    } else {
                        ll_fe_only = x$ll_fe_only
                    }

                    stat = 2 * (x$loglik - ll_fe_only)
                    df = length(x$coefficients)
                    p = pchisq(stat, df, lower.tail = FALSE)

                    vec = list(stat = stat, p = p, df = df)
                    res_all[[type]] = set_value(vec, value)

                } else {
                    res_all[[type]] = NA
                }

            } else if(root %in% user_types){

                res_all[[type]] = set_value(opts[[root]]$fun(x), value)

            } else {
                stop("The type '", type, "' is not supported, see details of ?fitstat, or use fitstat(show_types = TRUE) to get the names of all supported tests.")
            }
        }

    }

    if(simplify){
        verbose = FALSE
    }

    if(verbose){
        class(res_all) = "fixest_fitstat"
    }

    if(length(type_all) == 1 && simplify){
        return(res_all[[1]])
    }

    res_all
}


#' Wald test of nullity of coefficients
#'
#' Wald test used to test the joint nullity of a set of coefficients.
#'
#' @inheritParams print.fixest
#' @inheritParams summary.fixest
#' @inheritParams etable
#'
#' @param print Logical, default is \code{TRUE}. If \code{TRUE}, then a verbose description of the test is prompted on the R console. Otherwise only a named vector containing the test statistics is returned.
#' @param ... Any other element to be passed to \code{\link[fixest]{summary.fixest}}.
#'
#' @details
#' The type of VCOV matrix plays a crucial role in this test. Use the arguments \code{se} and \code{cluster} to change the type of VCOV for the test.
#'
#' @return
#' A named vector containing the following elements is returned: \code{stat}, \code{p}, \code{df1}, and \code{df2}. They correspond to the test statistic, the p-value, the first and second degrees of freedoms.
#'
#' If no valud coefficient is found, the value \code{NA} is returned.
#'
#' @examples
#'
#' data(airquality)
#'
#' est = feols(Ozone ~ Solar.R + Wind + poly(Temp, 3), airquality)
#'
#' # Testing the joint nullity of the Temp polynomial
#' wald(est, "poly")
#'
#' # Same but with clustered SEs
#' wald(est, "poly", cluster = "Month")
#'
#' # Now: all vars but the polynomial and the intercept
#' wald(est, drop = "Inte|poly")
#'
#' #
#' # Toy example: testing pre-trends
#' #
#'
#' data(base_did)
#'
#' est_did = feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
#'
#' # The graph of the coefficients
#' coefplot(est_did)
#'
#' # The pre-trend test
#' wald(est_did, "period::[1234]$")
#'
#' # If "period::[1234]$" looks weird to you, check out
#' # regular expressions: e.g. see ?regex.
#' # Learn it, you won't regret it!
#'
#'
wald = function(x, keep = NULL, drop = NULL, print = TRUE, se, cluster, ...){
    # LATER:
    # - keep can be a list
    #   * list("fit_" = 1, "x5$")
    #   * regex = restriction. No "=" => 0

    check_arg(x, "class(fixest)")
    check_arg(keep, drop, "NULL character vector no na")

    if(isTRUE(x$onlyFixef)) return(NA)

    dots = list(...)
    if(!isTRUE(x$summary) || !missing(se) || !missing(cluster) || ...length() > 0){
        x = summary(x, se = se, cluster = cluster, ...)
    }

    if(missing(keep) && missing(drop)){
        drop = "\\(Intercept\\)"
    }

    coef_name = names(x$coefficients)
    coef_name = keep_apply(coef_name, keep)
    coef_name = drop_apply(coef_name, drop)

    if(length(coef_name) == 0) return(NA)

    qui = names(x$coefficients) %in% coef_name
    my_coef = x$coefficients[qui]
    df1 = length(my_coef)

    df2 = x$nobs - x$nparams

    # The VCOV is always full rank in here
    stat = drop(my_coef %*% solve(x$cov.scaled[qui, qui]) %*% my_coef) / df1
    p = pf(stat, df1, df2, lower.tail = FALSE)
    vcov = attr(x$cov.scaled, "type")
    vec = list(stat = stat, p = p, df1 = df1, df2 = df2, vcov = vcov)

    if(print){
        cat("Wald test, H0: ", ifsingle(coef_name, "", "joint "), "nullity of ", enumerate_items(coef_name), "\n", sep  ="")
        cat(" stat = ", numberFormatNormal(stat),
            ", p-value ", ifelse(p < 2.2e-16, "< 2.2e-16", paste0("= ", numberFormatNormal(p))),
            ", on ", numberFormatNormal(df1), " and ", numberFormatNormal(df2), " DoF,",
            "VCOV: ", vcov, ".", sep = "")

        return(invisible(vec))
    } else {
        return(vec)
    }
}

fitstat_validate = function(x, vector = FALSE){
    check_value(x, "NA | os formula | charin(FALSE) | character vector no na", .arg_name = "fitstat", .up = 1)

    if("formula" %in% class(x)){
        x = attr(terms(update(~ DEFAULT, x)), "term.labels")
    } else if (length(x) == 1 && (isFALSE(x) || is.na(x))){
        x = c()
    }

    # checking the types
    fitstat_fun_types = fitstat(give_types = TRUE)
    fitstat_type_allowed = fitstat_fun_types$types
    x = unique(x)
    type_alias = fitstat_fun_types$type_alias

    if(any(x %in% names(type_alias))){
        i = intersect(x, names(type_alias))
        x[x %in% i] = type_alias[x[x %in% i]]
    }

    pblm = setdiff(x, c(fitstat_type_allowed, "DEFAULT"))
    if(length(pblm) > 0){
        stop_up("In fitstat, argument 'x' must be a one sided formula (or a character vector) containing valid types from the function fitstat (see details in ?fitstat or use fitstat(show_types = TRUE)). The type", enumerate_items(pblm, "s.is.quote"), " not valid.")
    }

    if(length(x) == 0){
        x = NA
    } else if(vector){
        x = gsub("DEFAULT", ".", x, fixed = TRUE)
    } else {
        x = as.formula(paste("~", paste(gsub("DEFAULT", ".", x), collapse = "+")))
    }

    x
}


#' R2s of \code{fixest} models
#'
#' Reports different R2s for \code{fixest} estimations (e.g. \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}).
#'
#' @param x A \code{fixest} object, e.g. obtained with function \code{\link[fixest]{feglm}} or \code{\link[fixest]{feols}}.
#' @param type A character vector representing the R2 to compute. The R2 codes are of the form: "wapr2" with letters "w" (within), "a" (adjusted) and "p" (pseudo) possibly missing. E.g. to get the regular R2: use \code{type = "r2"}, the within adjusted R2: use \code{type = "war2"}, the pseudo R2: use \code{type = "pr2"}, etc. Use \code{"cor2"} for the squared correlation. By default, all R2s are computed.
#' @param full_names Logical scalar, default is \code{FALSE}. If \code{TRUE} then names of the vector in output will have full names instead of keywords (e.g. \code{Squared Correlation} instead of \code{cor2}, etc).
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
#' r2(est, "cor2")
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
#' r2(est, c("cor2", "r2", "pr2", "war2"))
#'
#' # same with full names instead of codes
#' r2(est, c("cor2", "r2", "pr2", "war2"), full_names = TRUE)
#'
r2 = function(x, type = "all", full_names = FALSE){
    # p: pseudo
    # w: within
    # a: adjusted

    check_arg(full_names, "logical scalar")

    if(!"fixest" %in% class(x)){
        stop("Only 'fixest' objects are supported.")
    }

    check_arg(type, "character vector no na", .message = "Argument 'type' must be a character vector (e.g. type = c(\"cor2\", \"r2\", \"pr2\")). (a: adjused, p: pseudo, w: within.)")

    # type_allowed next => ("count", "acount") ?
    dict_names = c("cor2" = "Squared Correlation", "r2" = "R2", "ar2" = "Adjusted R2", "pr2" = "Pseudo R2", "apr2" = "Adjusted Pseudo R2", "wr2" = "Within R2", "war2" = "Adjusted Within R2", "wpr2" = "Within Pseudo R2", "wapr2" = "Adjusted Within Pseudo R2")
    type_allowed = c("cor2", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2")
    types_alias = c(par2 = "apr2", awr2 = "war2", pwr2 = "wpr2", wpar2 = "wapr2", pwar2 = "wapr2", pawr2 = "wapr2", apwr2 = "wapr2", awpr2 = "wapr2", sq.cor = "cor2")
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

    is_ols = x$method_type == "feols"
    isFixef = "fixef_vars" %in% names(x)
    n = nobs(x)

    res = rep(NA, length(type_all))
    for(i in seq_along(type_all)){
        myType = type_all[i]

        if(myType == "cor2"){
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

                if(isTRUE(x$lean)){
                    stop("Within R2s are not available for 'lean' fixest objects. Please reestimate with 'lean = FALSE'.")
                }

                # constructing the data
                newdata = cbind(data.frame(y = x$y), as.data.frame(x$fixef_id))
                if(!is.null(x$fixef_terms)){
                    newdata = cbind(newdata, as.data.frame(x$slope_variables))
                }
                # Fe/slope only formula
                new_fml = merge_fml(y ~ 1, x$fml_all$fixef)

                # x$family$family is also normal
                res_fe = feglm(fml = new_fml, data = newdata, glm.tol = 1e-2, fixef.tol = 1e-3, family = x$family$family, weights = x$weights, offset = x$offset)

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


#' Gets the degrees of freedom of a \code{fixest} estimation
#'
#' Simple utility to extract the degrees of freedom from a \code{fixest} estimation.
#'
#' @inheritParams vcov.fixest
#'
#' @param x A \code{fixest} estimation.
#' @param type Character scalar, equal to "k", "resid", "t". If "k", then the number of regressors is returned. If "resid", then it is the "residuals degree of freedom", i.e. the number of observations minus the number of regressors. If "t", it is the degrees of freedom used in the t-test. Note that these values are affected by how the VCOV of \code{x} is computed, in particular when the VCOV is clustered.
#' @param vars A vector of variable names, of the regressors. This is optional. If provided, then \code{type} is set to 1 by default and the number of regressors contained in \code{vars} is returned. This is only useful in the presence of collinearity and we want a subset of the regressors only. (Mostly for internal use.)
#' @param stage Either 1 or 2. Only concerns IV regressions, which stage to look at.
#'
#'
#' @examples
#'
#' # First: an estimation
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#' est = feols(y ~ x1 + x2 | species, base)
#'
#' # "Normal" standard-errors (SE)
#' est_standard = summary(est, se = "st")
#'
#' # Clustered SEs
#' est_clustered = summary(est, se = "clu")
#'
#' # The different degrees of freedom
#'
#' # => different type 1 DoF (because of the clustering)
#' degrees_freedom(est_standard, type = "k")
#' degrees_freedom(est_clustered, type = "k") # fixed-effects are excluded
#'
#' # => different type 2 DoF (because of the clustering)
#' degrees_freedom(est_standard, type = "resid") # => equivalent to the df.residual from lm
#' degrees_freedom(est_clustered, type = "resid")
#'
#'
#'
degrees_freedom = function(x, type, vars = NULL, se = NULL, cluster = NULL, dof = NULL, stage = 2){
    check_arg(x, "class(fixest) mbt")
    check_arg_plus(type, "match(k, resid, t)")
    check_arg(stage, "integer scalar GE{1} LE{2}")
    check_arg(vars, "character vector no na")

    if(stage == 1 && isTRUE(x$iv) && x$iv_stage == 2){
        x = x$iv_first_stage[[1]]
    }

    if(!missnull(vars)){
        if(!missing(type) && type != "k"){
            warning("The argument 'type' is ignored when the argument 'vars' is present. Type 'k' is returned.")
        }

        vars_keep = intersect(vars, names(x$coefficients))
        return(length(vars_keep))
    }

    if(missing(type)){
        stop("The argument 'type' is required but is currently missing.")
    }

    if(!isTRUE(x$summary) || !missnull(se) || !missnull(cluster) || !missnull(dof)){
        x = summary(x, se = se, cluster = cluster, dof = dof)
    }

    vcov = x$cov.scaled

    if(is.null(vcov)){
        dof.K = x$nparams
        t.df = NULL
    } else {
        t.df = attr(vcov, "G")
        dof.K = attr(vcov, "dof.K")
    }


    if(type == "k"){
        res = dof.K
    } else if(type == "resid"){
        res = x$nobs - dof.K
    } else if(type == "t"){
        if(is.null(t.df)){
            res = nobs(x) - dof.K
        } else {
            res = t.df - 1
        }
    }

    res
}

####
#### Stats -- internal ####
####

.wald = function(x, var){
    # x: fixest estimation

    coef = x$coefficients

    vcov = x$cov.scaled
    if(is.null(vcov)){
        # => INTERNAL ERROR
        stop("INTERNAL ERROR: .wald should be applied only to objects with a VCOV already computed. Could you report the error to the maintainer of fixest?")
    }

    var_keep = intersect(var, names(coef))

    if(length(var_keep) == 0){
        # All vars removed bc of collin => stat = 0
        return(0)
    }

    # To handle errors (rare but can happen)
    if(any(diag(vcov) < 1e-10)){
        vcov = mat_posdef_fix(vcov)
    }

    chol_decomp = cpp_cholesky(vcov[var_keep, var_keep, drop = FALSE])
    vcov_inv = chol_decomp$XtX_inv

    if(isTRUE(chol_decomp$all_removed)){
        # Can happen for indicators + clustering at the same level
        return(1e6)
    }

    if(any(chol_decomp$id_excl)){
        var_keep = var_keep[!chol_decomp$id_excl]
    }

    my_coef = coef[var_keep]

    stat = drop(my_coef %*% vcov_inv %*% my_coef) / length(my_coef)

    stat
}


kp_stat = function(x){
    # internal function => x must be a fixest object
    #
    # The code here is a translation of the ranktest.jl function from the Vcov.jl package
    # from @matthieugomez (see https://github.com/matthieugomez/Vcov.jl)
    #


    if(!isTRUE(x$iv) || !x$iv_stage == 2) return(NA)

    # Necessary data

    X_proj = as.matrix(resid(summary(x, stage = 1)))

    # while the projection of X has already been computed, we need to do it
    # for Z => that's a pain in the neck
    # computational cost is much lower is I include the KP computation at
    # estimation time.... since here we do the job twice... we'll see

    Z = model.matrix(x, type = "iv.inst")
    Z_proj = proj_on_U(x, Z)

    k = n_endo = ncol(X_proj)
    l = n_inst = ncol(Z)

    # We assume l >= k
    q = min(k, l) - 1

    # Let's go

    if(n_endo == 1){
        PI = t(coef(summary(x, stage = 1)))
    } else {
        PI = coef(summary(x, stage = 1))
    }
    PI = PI[, colnames(PI) %in% x$iv_inst_names_xpd, drop = FALSE]

    Fmat = chol(crossprod(Z_proj))
    Gmat = chol(crossprod(X_proj))
    theta = Fmat %*% t(solve(t(Gmat)) %*% PI)
    # theta: n_inst x n_endo

    if(n_inst == n_endo){
        svd_decomp = svd(theta)
        u = svd_decomp$u
        vt = t(svd_decomp$v)
    } else {
        # we need full decomp => un optimized decomp
        svd_decomp = mat_svd(theta)
        u = svd_decomp$u
        vt = svd_decomp$vt
    }

    u_sub = u[k:l, k:l, drop = FALSE]
    vt_sub = vt[k, k]

    vt_k = vt[1:k, k, drop = FALSE]

    ssign = function(x) if(x == 0) 1 else sign(x)

    # There may be more sign problems here
    if(k == l){
        a_qq = ssign(u_sub[1]) * u[1:l, k:l, drop = FALSE]
        b_qq = ssign(vt_sub[1]) * t(vt_k)
    } else {
        a_qq = u[1:l, k:l] %*% (solve(u_sub) %*% mat_sqrt(u_sub %*% t(u_sub)))
        b_qq = mat_sqrt(vt_sub %*% t(vt_sub)) %*% (solve(t(vt_sub)) %*% t(vt_k))
    }

    # kronecker
    kronv = kronecker(b_qq, t(a_qq))
    lambda = kronv %*% c(theta)

    # There is need to compute the vcov specifically for this case
    # We do it the same way as it was for x

    if(identical(x$se_info$se, "standard")){
        vlab = chol(tcrossprod(kronv) / nrow(X_proj))

    } else {
        K = t(kronecker(Gmat, Fmat))

        my_scores = do.call(cbind, lapply(1:ncol(X_proj), function(i) Z_proj * X_proj[, i]))

        x_new = x
        x_new$scores = my_scores

        se = x$summary_flags$se
        cluster = x$summary_flags$cluster
        dof = x$summary_flags$dof

        meat = vcov(x_new, se = se, cluster = cluster, dof = dof, meat_only = TRUE)
        vhat = solve(K, t(solve(K, meat)))

        # DOF correction now
        n = nobs(x) - (x$se_info$se == "cluster")
        df_resid = degrees_freedom(x, "resid", stage = 1)
        vhat = vhat * n / df_resid

        vlab = kronv %*% vhat %*% t(kronv)
    }

    r_kp = t(lambda) %*% solve(vlab, lambda)
    # cat("KP r:\n") ; print(r_kp)

    # Now the results
    kp_df = n_inst - n_endo + 1
    kp_f = r_kp / n_inst
    kp_p = pchisq(kp_f, kp_df, lower.tail = FALSE)

    list(stat = kp_f, p = kp_p, df = kp_df)
}

cd_stat = function(x){
    # internal function
    # x: fixest object
    #

    if(!isTRUE(x$iv) || !x$iv_stage == 2) return(NA)

    # Necessary data
    X = model.matrix(x, type = "iv.endo")
    X_proj = proj_on_U(x, X)

    Z = model.matrix(x, type = "iv.inst")
    Z_proj = proj_on_U(x, Z)

    n = nobs(x)
    n_endo = ncol(X_proj)
    n_inst = ncol(Z_proj)

    # The stat
    if(FALSE){
        # => just use the canonical correlation, it's faster

        df_resid = degrees_freedom(x, resid, stage = 1)
        V = resid(summary(x, stage = 1))

        P_Z_proj = Z_proj %*% solve(crossprod(Z_proj)) %*% t(Z_proj)

        SIGMA_vv_hat = (t(X) %*% V) / df_resid

        inv_sqrt_SIGMA_vv_hat = solve(mat_sqrt(SIGMA_vv_hat))

        G = (t(inv_sqrt_SIGMA_vv_hat) %*% t(X_proj) %*% P_Z_proj %*% X_proj %*% inv_sqrt_SIGMA_vv_hat) / n_inst

        cd = min(eigen(G)$values)

    } else {
        cc = min(cancor(X_proj, Z_proj)$cor)
        cd = ((n - n_endo - n_inst - 1) / n_inst) / ((1 - cc^2) / cc^2)
    }

    cd
}


####
#### Misc internal funs ####
####

mat_sqrt = function(A){
    e = eigen(A)
    # e$vectors %*% diag(x = sqrt(e$values), nrow = nrow(A)) %*% t(e$vectors)
    ev = pmax(e$values, 1e-10)
    M = e$vectors %*% diag(x = ev**.25, nrow = nrow(A))
    tcrossprod(M)
}


mat_svd = function(A){
    # From https://rpubs.com/aaronsc32/singular-value-decomposition-r
    # https://math.stackexchange.com/questions/2359992/how-to-resolve-the-sign-issue-in-a-svd-problem

    ATA = crossprod(A)
    ATA.e = eigen(ATA)
    v = ATA.e$vectors

    AAT = tcrossprod(A)
    AAT.e = eigen(AAT)
    u = AAT.e$vectors
    r = sqrt(ATA.e$values)

    # we want the last diag element of vt to be positive
    vt = t(v)
    k = nrow(vt)
    if(nrow(u) == k && vt[k, k] < 0){
        vt[k, ] = - vt[k, ]
        u[, k] = - u[, k]
    }

    list(u = u, vt = vt, d = r)
}





proj_on_U = function(x, Z){
    # Projects some variables (either the endo or the inst) on the
    # exogenous variables (U)
    #
    # x: a fixest IV estimation
    #
    # I write Z in the arguments, but it can be any matrix (like X)
    #

    U = model.matrix(x, type = "iv.exo")
    is_U = !is.null(U)

    if(!is.null(x$fixef_vars)){
        # we need to center the variables
        # Of course, if there are varying slopes, that's a pain in the neck

        if(!is.null(x$slope_flag)){
            Z_dm = demean(Z, x$fixef_id[order(x$fixef_sizes, decreasing = TRUE)],
                          slope.vars = x$slope_variables_reordered,
                          slope.flag = x$slope_flag_reordered)

            if(is_U){
                U_dm = demean(U, x$fixef_id[order(x$fixef_sizes, decreasing = TRUE)],
                              slope.vars = x$slope_variables_reordered,
                              slope.flag = x$slope_flag_reordered)
            }

        } else {
            Z_dm = demean(Z, x$fixef_id)
            if(is_U){
                U_dm = demean(U, x$fixef_id)
            }
        }

        if(is_U){
            Z_proj = resid(feols.fit(Z_dm, U_dm))
        } else {
            Z_proj = Z_dm
        }

    } else if(is_U){
        Z_proj = resid(feols.fit(Z, U))
    } else {
        Z_proj = Z
    }

    Z_proj
}









































