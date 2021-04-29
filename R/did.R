#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Thu Apr 29 11:24:49 2021
# ~: DiD functions
#----------------------------------------------#






#' Sun and Abraham interactions
#'
#' User-level method to implement staggered difference-in-difference estimations a la Sun and Abraham (Journal of Econometrics, forthcoming).
#'
#' @param cohort A vector representing the cohort. It should represent the period at which the treatment has been received (and thus be fixed for each unit).
#' @param period A vector representing the period. It can be either a relative time period (with negative values representing the before the treatment and positive values after the treatment), or a regular time period. In the latter case, the relative time period will be created from the cohort information (which represents the time at which the treatment has been received).
#' @param ref.c A vector of references for the cohort. By default the never treated cohorts are taken as reference and the always treated are excluded from the estimation. You can add more references with this argument, which means that dummies will not be created for them (but they will remain in the estimation).
#' @param ref.p A vector of references for the (relative!) period. By default the first relative period (RP) and the 0 are taken as reference. You can instead use your own references (i.e. RPs for which dummies will not be created -- but these observations remain in the sample). Please note that you will need at least two references. You can use the special variables \code{.F} and \code{.L} to access the first and the last relative periods. Using these notations, the default is \code{c(.F, 0)}.
#'
#' @details
#' This function creates a matrix of \code{cohort x relative_period} interactions, and if used within a \code{fixest} estimation, the coefficients will automatically be aggregated to obtain the ATT for each relative period. In practice, the coefficients are aggregated with the \code{\link[fixest]{aggregate.fixest}} function whose argument \code{agg} is automatically set to the appropriate value.
#'
#' The SA method requires relative periods (negative/positive for before/after the treatment). Either the user can compute the RP (relative periods) by his/her own, either the RPs are computed on the fly from the periods and the cohorts (which then should represent the treatment period).
#'
#' The never treated, which are the cohorts displaying only negative RPs are used as references (i.e. no dummy will be constructed for them). On the other hand, the always treated are removed from the estimation, by means of adding NAs for each of their observations.
#'
#' If the RPs have to be constructed on the fly, any cohort that is not present in the period is considered as never treated. This means that if the period ranges from 1995 to 2005, \code{cohort = 1994} will be considered as never treated, although it should be considered as always treated: so be careful.
#'
#' If you construct your own relative periods, the controls cohorts should have only negative RPs.
#'
#' @return
#' If not used within a \code{fixest} estimation, this function will return a matrix of interacted coefficients.
#'
#' @examples
#'
#'
#'
#'
#'
sunab = function(cohort, period, ref.c = NULL, ref.p = c(.F, 0), att = FALSE, no_agg = FALSE){
    # LATER:
    # - add id or indiv argument, just to remove always treated
    # - add argument bin.p

    check_arg(cohort, "mbt vector")
    check_arg(period, "mbt vector len(data)", .data = cohort)
    check_arg(ref.c, "NULL vector no na")
    check_arg(att, no_agg, "logical scalar")

    period_name = deparse_long(substitute(period))
    period_name = gsub("^[[:alpha:]][[:alpha:]_\\.]*\\$", "", period_name)

    # Finding out what kind of data that is

    # NAness (big perf hit)
    n_origin = length(cohort)
    IS_NA = which(is.na(cohort) | is.na(period))
    ANY_NA = length(IS_NA) > 0
    if(ANY_NA){
        cohort = cohort[-IS_NA]
        period = period[-IS_NA]
    }

    n = length(cohort)

    # Case 1
    # cohort: typically the year of treatment
    # period: relative period
    # we don't need to do anything


    # Case 2
    # cohort: year of treatment
    # period: the year
    # => we compute the relative period, we exclude the never/always treated

    # We find out the never/always treated with:
    #  absence of variance in the relative time

    period_unik = unique(period)
    cohort_unik = unique(cohort)

    # CASE 1
    is_CASE_1 = FALSE
    if(is.numeric(period) && 0 %in% period_unik && min(period_unik) < 0 && max(period_unik) > 0){
        # Case 1 => we don't need to do anything
        is_CASE_1 = TRUE
    } else {

        # CASE 2 => construction of the relative period
        refs = setdiff(cohort_unik, period_unik)

        if(length(refs) == length(cohort_unik)){
            stop("Problem in the creation of the relative time periods. We expected the cohort to be the treated period, yet not a single 'cohort' value was found in 'period'.")
        }

        qui_keep = which(!cohort %in% refs)
        cohort_valid = cohort[qui_keep]
        period_valid = period[qui_keep]

        if(is.numeric(period_valid) || is.numeric(cohort_valid)){

            # simplest case
            if(!is.numeric(cohort_valid)) cohort_valid = as.numeric(cohort_valid)
            if(!is.numeric(period_valid)) period_valid = as.numeric(period_valid)

            rel_period = period_valid - cohort_valid
        } else {
            # difficult case
            sunik_period = sort(unique(period_valid))
            dict_period = seq_along(sunik_period)
            names(dict_period) = sunik_period

            period_valid = dict_period[as.character(period_valid)]
            cohort_valid = dict_period[as.character(cohort_valid)]

            rel_period = period_valid - cohort_valid
        }

        # I put a negative number so that they are not considered as always treated
        new_period = rep(-1, n)
        new_period[qui_keep] = rel_period

        period = new_period
    }

    # Now the reference expressed in relative periods
    check_arg_plus(ref.p, "evalset integer vector no na", .data = list(.F = min(period), .L = max(period)))
    if(missing(ref.p)) ref.p = c(min(period), 0)

    #
    #  we find out the never/always treated
    #

    cohort_int = quickUnclassFactor(cohort)
    c_order = order(cohort_int)
    info = cpp_find_never_always_treated(cohort_int[c_order], period[c_order])

    if(!is.null(ref.c)){
        qui_drop = which(cohort_int %in% info$ref | period %in% ref.p | cohort %in% ref.c)
    } else {
        qui_drop = which(cohort_int %in% info$ref | period %in% ref.p)
    }
    qui_NA = info$always_treated

    # All references have been removed => pure i() without ref
    cohort = cohort[-qui_drop]
    period = period[-qui_drop]

    res_raw = i(f = period, f2 = cohort, f_name = period_name)

    # We extend the matrix

    if(ANY_NA){
        # We DON'T remove the references from the estimation
        res = matrix(NA_real_, nrow = n_origin, ncol = ncol(res_raw), dimnames = list(NULL, colnames(res_raw)))
        res[-IS_NA, ][qui_drop, ] = 0
        res[-IS_NA, ][-qui_drop, ] = res_raw
        if(length(qui_NA) > 0){
            res[-IS_NA, ][qui_NA, ] = NA_real_
        }
    } else {
        res = matrix(0, nrow = n_origin, ncol = ncol(res_raw), dimnames = list(NULL, colnames(res_raw)))
        res[-qui_drop, ] = res_raw
        if(length(qui_NA) > 0){
            res[qui_NA, ] = NA_real_
        }
    }


    # We add the agg argument to GLOBAL_fixest_mm_info
    if(!no_agg){
        is_GLOBAL = FALSE
        for(where in 1:min(8, sys.nframe())){
            if(exists("GLOBAL_fixest_mm_info", parent.frame(where))){
                GLOBAL_fixest_mm_info = get("GLOBAL_fixest_mm_info", parent.frame(where))
                is_GLOBAL = TRUE
                break
            }
        }

        if(is_GLOBAL){
            if(att){
                agg = c("ATT" = paste0("\\Q", period_name, "\\E::[[:digit:]]+"))
            } else {
                agg = paste0("(\\Q", period_name, "\\E)::(-?[[:digit:]]+)")
            }

            GLOBAL_fixest_mm_info$agg = agg
            # re assignment
            assign("GLOBAL_fixest_mm_info", GLOBAL_fixest_mm_info, parent.frame(where))
        }
    }

    res
}

#' @rdname sunab
sunab_att = function(cohort, period, ref.c = NULL, ref.p = c(.F, 0)){
    sunab(cohort, period, ref.c, ref.p, att = TRUE)
}



#' Aggregates the values of DiD coefficients a la Sun and Abraham
#'
#' Simple tool that aggregates the value of CATT coefficients in staggered difference-in-difference setups (see details).
#'
#' @param x A \code{fixest} object.
#' @param agg A character scalar describing the variable names to be aggregated, it is pattern-based. All variables that match the pattern will be aggregated. It must be of the form \code{"(root)"}, the parentheses must be there and the resulting variable name will be \code{"root"}. You can add another root with parentheses: \code{"(root1)regex(root2)"}, in which case the resulting name is \code{"root1::root2"}. To name the resulting variable differently you can pass a named vector: \code{c("name" = "pattern")} or \code{c("name" = "pattern(root2)")}. It's a bit intricate sorry, please see the examples.
#' @param full Logical scalar, defaults to \code{FALSE}. If \code{TRUE}, then all coefficients are returned, not only the aggregated coefficients.
#' @param use_weights Logical, default is \code{TRUE}. If the estimation was weighted, whether the aggregation should take into account the weights. Basically if the weights reflected frequency it should be \code{TRUE}.
#' @param ... Arguments to be passed to \code{\link[fixest]{summary.fixest}}.
#'
#' @details
#' This is a function helping to replicate the estimator from Sun and Abraham (2020). You first need to perform an estimation with cohort and relative periods dummies (typically using the function \code{\link[fixest]{i}}), this leads to estimators of the cohort average treatment effect on the treated (CATT). Then you can use this function to retrieve the average treatment effect on each relative period, or for any other way you wish to aggregate the CATT.
#'
#' Note that contrary to the SA article, here the cohort share in the sample is considered to be a perfect measure for the cohort share in the population.
#'
#' @return
#' It returns a matrix representing a table of coefficients.
#'
#' @references
#' Liyang Sun and Sarah Abraham, forthcoming, "Estimating Dynamic Treatment Effects in Event Studies with Heterogeneous Treatment Effects". Journal of Econometrics.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' #
#' # DiD example
#' #
#'
#' # first we set up the data
#'
#' set.seed(1)
#' n_group = 20
#' n_per_group = 5
#'
#' id_i = paste0((1:n_group), ":", rep(1:n_per_group, each = n_group))
#' id_t = 1:10
#'
#' base = expand.grid(id = id_i, year = id_t)
#' base$group = as.numeric(gsub(":.+", "", base$id))
#'
#' base$year_treated = base$group
#' base$year_treated[base$group > 10] = 10000
#' base$treat_post = (base$year >= base$year_treated) * 1
#' base$time_to_treatment = pmax(base$year - base$year_treated, -1000)
#' base$treated = (base$year_treated < 10000) * 1
#'
#' # The effect of the treatment is cohort specific and increases with time
#' base$y_true = base$treat_post * (1 + 1 * base$time_to_treatment - 1 * base$group)
#' base$y = base$y_true + rnorm(nrow(base))
#'
#'
#' # The controls have a time_to_treatment equal to -1000
#'
#' # we drop the always treated
#' base = base[base$group > 1,]
#'
#' # Now we perform the estimation
#' res_naive = feols(y ~ i(treated, time_to_treatment,
#'                         ref = -1, drop = -1000) | id + year, base)
#'
#' res_cohort = feols(y ~ i(time_to_treatment, f2 = group,
#'                          drop = c(-1, -1000)) | id + year, base)
#'
#' coefplot(res_naive, ylim = c(-6, 8))
#' att_true = tapply(base$y_true, base$time_to_treatment, mean)[-1]
#' points(-9:8 + 0.15, att_true, pch = 15, col = 2)
#'
#' # The aggregate effect for each period
#' agg_coef = aggregate(res_cohort, "(ti.*nt)::(-?[[:digit:]]+)")
#' x = c(-9:-2, 0:8) + .35
#' points(x, agg_coef[, 1], pch = 17, col = 4)
#' ci_low = agg_coef[, 1] - 1.96 * agg_coef[, 2]
#' ci_up = agg_coef[, 1] + 1.96 * agg_coef[, 2]
#' segments(x0 = x, y0 = ci_low, x1 = x, y1 = ci_up, col = 4)
#'
#' legend("topleft", col = c(1, 2, 4), pch = c(20, 15, 17),
#'        legend = c("Naive", "True", "Sun & Abraham"))
#'
#'
#' # The ATT
#' aggregate(res_cohort, c("ATT" = "treatment::[^-]"))
#' mean(base[base$treat_post == 1, "y_true"])
#'
#' # With etable
#' etable(res_naive, res_cohort, agg = "(ti.*nt)::(-?[[:digit:]]+):gro")
#'
aggregate.fixest = function(x, agg, full = FALSE, use_weights = TRUE, ...){
    # Aggregates the value of coefficients

    check_arg(x, "class(fixest) mbt")
    check_arg(agg, "character scalar")
    check_arg(full, "logical scalar")
    # => later => extend it to more than one set of vars to agg

    dots = list(...)
    from_summary = isTRUE(dots$from_summary)

    is_name = !is.null(names(agg))

    if(!is_name && !grepl("(", agg, fixed = TRUE)){
        stop("Argument 'agg' must be a character in which the pattern to match must be in between parentheses. So far there are no parenthesis: please have a look at the examples.")
    }

    coef = x$coefficients
    cname = names(coef)

    qui = grepl(agg, cname, perl = TRUE)
    if(!any(qui)){
        if(from_summary){
            # We make it silent when aggregate is used in summary
            # => this way we can pool calls to agg even for models that don't have it
            # ==> useful in etable eg
            return(x$coeftable)
        } else {
            stop("The argument 'agg' does not match any variable.")
        }
    }

    if(!isTRUE(x$summary)){
        x = summary(x, ...)
    }

    cname_select = cname[qui]
    if(is_name){
        root = rep(names(agg), length(cname_select))
        val = gsub(paste0(".*", agg, ".*"), "\\1", cname_select, perl = TRUE)
    } else {
        root = gsub(paste0(".*", agg, ".*"), "\\1", cname_select, perl = TRUE)
        val = gsub(paste0(".*", agg, ".*"), "\\2", cname_select, perl = TRUE)
    }

    V = x$cov.scaled

    mm = model.matrix(x)

    name_df = unique(data.frame(root, val, stringsAsFactors = FALSE))

    c_all = c()
    se_all = c()
    for(i in 1:nrow(name_df)){

        r = name_df[i, 1]
        v = name_df[i, 2]
        v_names = cname_select[root == r & val == v]

        if(use_weights && !is.null(x$weights)){
            shares = colSums(x$weights * sign(mm[, v_names, drop = FALSE]))
        } else {
            shares = colSums(sign(mm[, v_names, drop = FALSE]))
        }

        shares = shares / sum(shares)

        # The coef
        c_value = sum(shares * coef[v_names])

        # The variance
        n = length(v_names)
        s1 = matrix(shares, n, n)
        s2 = matrix(shares, n, n, byrow = TRUE)

        var_value = sum(s1 * s2 * V[v_names, v_names])
        se_value = sqrt(var_value)

        c_all[length(c_all) + 1] = c_value
        se_all[length(se_all) + 1] = se_value
    }

    # th z & p values
    zvalue <- c_all/se_all
    if(x$method_type == "feols" || (x$method %in% "feglm" && !x$family$family %in% c("poisson", "binomial"))){

        # I have renamed t.df into G
        t.df = attr(vcov, "G")

        if(!is.null(t.df)){
            pvalue <- 2*pt(-abs(zvalue), max(t.df - 1, 1))
        } else {
            pvalue <- 2*pt(-abs(zvalue), max(x$nobs - x$nparams, 1))
        }

    } else {
        pvalue <- 2*pnorm(-abs(zvalue))
    }

    res = cbind(c_all, se_all, zvalue, pvalue)
    if(max(nchar(val)) == 0){
        rownames(res) = name_df[[1]]
    } else {
        rownames(res) = apply(name_df, 1, paste, collapse = "::")
    }

    colnames(res) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

    if(full){
        table_origin = x$coeftable
        i_min = min(which(qui)) - 1
        before = if(i_min > 0) table_origin[1:i_min, , drop = FALSE] else NULL

        i_after = (1:nrow(table_origin)) > i_min & !qui
        after = if(any(i_after)) table_origin[i_after, , drop = FALSE] else NULL

        res = rbind(before, res, after)

        attr(res, "type") = attr(table_origin, "type")
    }

    res
}











































