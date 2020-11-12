#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Jan 28 10:03:04 2020
# ~: Set of deprecated (near defunct) functions
#----------------------------------------------#

#
# I found a great tutorial for documenting deprecated functions
# https://www.r-bloggers.com/r-deprecate-functions-with-roxygen2-3/
#

## fixest-deprecated.r

#' @title Deprecated functions in package \pkg{fixest}.
#'
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("-deprecated")}.
#'
#' @name fixest-deprecated
#'
#' @keywords internal
NULL


#' @rdname fixest-deprecated
#' @name did_estimate_yearly_effects-deprecated
#' @section \code{did_estimate_yearly_effects}:
#' From \code{fixest} version 0.3.0 onwards, estimating yearly effects with a simple interaction of the type \code{var::year(ref)} in the formula is preferred (see details in \code{\link[fixest]{feols}}). You can then plot the yearly treatment with the function \code{\link[fixest]{coefplot}}. You have examples detailed in \code{\link[fixest]{coefplot}}.
NULL

#' @title Estimates yearly treatment effects
#'
#' @description This facility helps to estimate yearly treatment effects in a difference-in-difference setup without having to manually make the estimation. It is made as general as possible such that non-\code{fixest} estimation functions can also be used.
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
#' @details
#' From \code{fixest} version 0.3.0 onwards, estimating yearly effects with a simple interaction of the type \code{var::year(ref)} in the formula is preferred (see details in \code{\link[fixest]{feols}}). You can then plot the yearly treatment with the function \code{\link[fixest]{coefplot}}. You have examples detailed in \code{\link[fixest]{coefplot}}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
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

    counter = getOption("fixest_deprec_did_estimate_yearly")
    if(is.null(counter)){
        options("fixest_deprec_did_estimate_yearly" = TRUE)
        .Deprecated("coefplot", msg = "From fixest version 0.3.0 onwards, estimating yearly effects with a simple interaction of the type var::year(ref) in the formula is preferred (see details in feols). You can then plot the yearly treatment with the function coefplot(). You have examples detailed in coefplot.")
    }


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
            stop("In argument treat_time, the variable", enumerate_items(qui_pblm, "s"), " not in the data.")
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
            stop("In argument treat_time, the variable", enumerate_items(qui_pblm, "s"), " not in the data.")
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
        stop("The treatment variable must be 0/1, with 1 repersenting the treated. The variable that you gave, ", treat_var, ", is not (e.g. observation", enumerate_items(obs, "s"), ".")
    } else if(length(unique(treat)) == 1){
        stop("The treatment variable is equal to ", treat[1], " for all observations: there must be a problem!")
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

#' @rdname fixest-deprecated
#' @name did_plot_yearly_effects-deprecated
#' @section \code{did_plot_yearly_effects}:
#' From \code{fixest} version 0.3.0 onwards, estimating yearly effects with a simple interaction of the type \code{var::year(ref)} in the formula is preferred (see details in \code{\link[fixest]{feols}}). You can then plot the yearly treatment with the function \code{\link[fixest]{coefplot}}. You have examples detailed in \code{\link[fixest]{coefplot}}.
NULL

#' @title Plots the results of yearly treatment effects estimation
#'
#' @description This function plots the results of a \code{\link[fixest]{did_estimate_yearly_effects}} estimation.
#'
#' @inheritParams errbar
#'
#' @param object An object returned by the function \code{\link[fixest]{did_estimate_yearly_effects}}.
#' @param ... Any other argument to be passed to \code{summary} or to \code{plot}.
#' @param style One of \code{"bar"} (default), \code{"interval"} or \code{"tube"}. The style of the plot.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
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

    counter = getOption("fixest_deprec_did_plot_yearly_effects")
    if(is.null(counter)){
        options("fixest_deprec_did_plot_yearly_effects" = TRUE)
        .Deprecated("coefplot", msg = "From fixest version 0.3.0 onwards, estimating yearly effects with a simple interaction of the type var::year(ref) in the formula is preferred (see details in feols). You can then plot the yearly treatment with the function coefplot(). You have examples detailed in coefplot.")
    }

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
    if(is.null(dots$ylab)) dots$ylab = paste0("Estimate of Yearly Treatment on ", deparse(object$fml[[2]]))

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



#' @rdname fixest-deprecated
#' @name errbar-deprecated
#' @section \code{errbar}:
#' The function \code{errbar} is now deprecated and will be removed in the near future. Use \code{\link[fixest]{coefplot}} instead, which dispose of much more possibilities..
NULL

#' @title Plots confidence intervals
#'
#' @description This function draws confidence intervals in a graph.
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
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
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
#' errbar(est$coefficients, x.shift = .2, add = TRUE, col = 2, bar.lty = 2, pch=15)
#'
errbar <- function(estimate, sd, ci_low, ci_top, x, x.shift = 0, w=0.1, ci_level = 0.95, style = c("bar", "interval", "tube"), add = FALSE, col = 1, bar.col = col, bar.lwd = par("lwd"), bar.lty, grid = TRUE, grid.par = list(lty=1), bar.join = FALSE, only.params = FALSE, ...){

    sc = sys.calls()
    if(length(sc) < 3 || !'do.call("errbar", dots)' == deparse(rev(sys.calls())[[2]])){
        # We warn only for direct calls

        counter = getOption("fixest_deprec_errbar")
        if(is.null(counter)){
            options("fixest_deprec_errbar" = TRUE)
            .Deprecated("coefplot", msg = "The function errbar is now deprecated and will be removed in the near future. Use coefplot instead, which dispose of much more possibilities.")
        }
    }

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

#' @rdname fixest-deprecated
#' @name obs2remove-deprecated
#' @section \code{obs2remove}:
#' This function will be removed in 1 year from 12/11/2020.
NULL

#' Finds observations to be removed from ML estimation with fixed-effects
#'
#' For Poisson, Negative Binomial or Logit estimations with fixed-effects, when the dependent variable is only equal to 0 (or 1 for Logit) for one fixed-effect value this leads to a perfect fit for that fixed-effect value by setting its associated fixed-effect coefficient to \code{-Inf}. Thus these observations need to be removed before estimation. This function gives the observations to be removed. Note that by default the function \code{\link[fixest]{femlm}} or \code{\link[fixest]{feglm}} drops them before performing the estimation.
#'
#' @param fml A formula containing the dependent variable and the fixed-effects. It can be of the type: \code{y ~ fixef_1 + fixef_2} or \code{y ~ x1 | fixef_1 + fixef_1} (in which case variables before the pipe are ignored).
#' @param data A data.frame containing the variables in the formula.
#' @param family Character scalar: either \dQuote{poisson} (default), \dQuote{negbin} or \dQuote{logit}.
#'
#' @return
#' It returns an integer vector of observations to be removed. If no observations are to be removed, an empty integer vector is returned. In both cases, it is of class \code{fixest.obs2remove}.
#' The vector has an attribute \code{fixef} which is a list giving the IDs of the fixed-effects that have been removed, for each fixed-effect dimension.
#'
#' @seealso
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
#'
#' @examples
#'
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' (x = obs2remove(v6 ~ Species, base))
#' attr(x, "fixef")
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
    # in the formula, the fixed-effects must be there:
    # either y ~ fixef_1 + fixef_2
    # either y ~ x1 + x2 | fixef_1 + fixef_2

    counter = getOption("fixest_deprec_obs2remove")
    if(is.null(counter)){
        options("fixest_deprec_obs2remove" = TRUE)
        .Deprecated(msg = "This function is deprecated and will disappear in 1 year from 12/11/2020.")
    }

    #
    # CONTROLS
    #

    # FAMILY

    family = match.arg(family)

    # FML

    if(!"formula" %in% class(fml) || length(fml) != 3){
        stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | fixef_1 + fixef_2' or of the type 'y ~ fixef_1 + fixef_2'.")
    }

    FML = Formula::Formula(fml)
    n_rhs = length(FML)[2]

    if(n_rhs > 2){
        stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | fixef_1 + fixef_2' or of the type 'y ~ fixef_1 + fixef_2'.")
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
    attr(obs2remove, "fixef") = dummyOmises

    return(obs2remove)
}

#' @rdname fixest-deprecated
#' @name summary.fixest.obs2remove-deprecated
#' @section \code{summary.fixest.obs2remove}:
#' This function will be removed in 1 year from 12/11/2020.
NULL

#' Summary method for fixest.obs2remove objects
#'
#' This function synthesizes the information of function \code{\link[fixest]{obs2remove}}. It reports the number of observations to be removed as well as the number of fixed-effects removed per fixed-effect dimension.
#'
#' @method summary fixest.obs2remove
#'
#' @param object A \code{fixest.obs2remove} object obtained from function \code{\link[fixest]{obs2remove}}.
#' @param ... Not currently used.
#'
#' @seealso
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
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
        cluster = attr(object, "fixef")
        cat("# fixed-effects removed: ", paste0(names(cluster), ": ", lengths(cluster), collapse = ", "), ".", sep = "")
    }

}







