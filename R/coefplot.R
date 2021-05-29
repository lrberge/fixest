#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Feb 10 09:46:01 2020
# ~: Simple function to display results
# from multiple estimations
#----------------------------------------------#



#' Plots confidence intervals and point estimates
#'
#' This function plots the results of estimations (coefficients and confidence intervals). The function \code{iplot} restricts the output to variables created with \code{\link[fixest]{i}}, either interactions with factors or raw factors.
#'
#' @inheritParams etable
#' @inheritSection etable Arguments keep, drop and order
#'
#' @param object Can be either: i) an estimation object (obtained for example from \code{\link[fixest]{feols}}, ii) a list of estimation objects (several results will be plotted at once), iii) a matrix of coefficients table, iv) a numeric vector of the point estimates -- the latter requiring the extra arguments \code{sd} or \code{ci_low} and \code{ci_high}.
#' @param sd The standard errors of the estimates. It may be missing.
#' @param ci_low If \code{sd} is not provided, the lower bound of the confidence interval. For each estimate.
#' @param ci_high If \code{sd} is not provided, the upper bound of the confidence interval. For each estimate.
#' @param horiz A logical scalar, default is \code{FALSE}. Whether to display the confidence intervals horizontally instead of vertically.
#' @param x The value of the x-axis. If missing, the names of the argument \code{estimate} are used.
#' @param x.shift Shifts the confidence intervals bars to the left or right, depending on the value of \code{x.shift}. Default is 0.
#' @param ci.width The width of the extremities of the confidence intervals. Default is \code{0.1}.
#' @param ci_level Scalar between 0 and 1: the level of the CI. By default it is equal to 0.95.
#' @param add Default is \code{FALSE}, if the intervals are to be added to an existing graph. Note that if it is the case, then the argument \code{x} MUST be numeric.
#' @param pt.pch The patch of the coefficient estimates. Default is 1 (circle).
#' @param cex Numeric, default is 1. Expansion factor for the points
#' @param pt.cex The size of the coefficient estimates. Default is the other argument \code{cex}.
#' @param col The color of the points and the confidence intervals. Default is 1 ("black"). Note that you can set the colors separately for each of them with \code{pt.col} and \code{ci.col}.
#' @param pt.col The color of the coefficient estimates. Default is equal to the other argument \code{col}.
#' @param ci.col The color of the confidence intervals. Default is equal to the other argument \code{col}.
#' @param lwd General line with. Default is 1.
#' @param pt.lwd The line width of the coefficient estimates. Default is equal to the other argument \code{lwd}.
#' @param ci.lwd The line width of the confidence intervals. Default is equal to the other argument \code{lwd}.
#' @param ci.lty The line type of the confidence intervals. Default is 1.
#' @param grid Logical, default is \code{TRUE}. Whether a grid should be displayed. You can set the display of the grid with the argument \code{grid.par}.
#' @param grid.par List. Parameters of the grid. The default values are: \code{lty = 3} and \code{col = "gray"}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. You also have two additional arguments: use \code{horiz = FALSE} to disable the horizontal lines, and use \code{vert = FALSE} to disable the vertical lines. Eg: \code{grid.par = list(vert = FALSE, col = "red", lwd = 2)}.
#' @param zero Logical, default is \code{TRUE}. Whether the 0-line should be emphasized. You can set the parameters of that line with the argument \code{zero.par}.
#' @param zero.par List. Parameters of the zero-line. The default values are \code{col = "black"} and \code{lwd = 1}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. Example: \code{zero.par = list(col = "darkblue", lwd = 3)}.
#' @param pt.join Logical, default is \code{FALSE}. If \code{TRUE}, then the coefficient estimates are joined with a line.
#' @param pt.join.par List. Parameters of the line joining the coefficients. The default values are: \code{col = pt.col} and \code{lwd = lwd}. You can add any graphical parameter that will be passed to \code{\link[graphics]{lines}}. Eg: \code{pt.join.par = list(lty = 2)}.
#' @param ref Used to add points equal to 0 (typically to visualize reference points). Either: i) "auto" (default), ii) a character vector of length 1, iii) a list of length 1, iv) a named integer vector of length 1, or v) a numeric vector. By default, in \code{iplot}, if the argument \code{ref} has been used in the estimation, these references are automatically added. If ii), ie a character scalar, then that coefficient equal to zero is added as the first coefficient. If a list or a named integer vector of length 1, then the integer gives the position of the reference among the coefficients and the name gives the coefficient name. A non-named numeric value of \code{ref} only works if the x-axis is also numeric (which can happen in \code{iplot}).
#' @param ref.line Logical or numeric, default is "auto", whose behavior depends on the situation. It is \code{TRUE} only if: i) interactions are plotted, ii) the x values are numeric and iii) a reference is found. If \code{TRUE}, then a vertical line is drawn at the level of the reference value. Otherwise, if numeric a vertical line will be drawn at that specific value.
#' @param ref.line.par List. Parameters of the vertical line on the reference. The default values are: \code{col = "black"} and \code{lty = 2}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. Eg: \code{ref.line.par = list(lty = 1, lwd = 3)}.
#' @param xlim.add A numeric vector of length 1 or 2. It represents an extension factor of xlim, in percentage. Eg: \code{xlim.add = c(0, 0.5)} extends \code{xlim} of 50\% on the right. If of length 1, positive values represent the right, and negative values the left (Eg: \code{xlim.add = -0.5} is equivalent to \code{xlim.add = c(0.5, 0)}).
#' @param ylim.add A numeric vector of length 1 or 2. It represents an extension factor of ylim, in percentage. Eg: \code{ylim.add = c(0, 0.5)} extends \code{ylim} of 50\% on the top. If of length 1, positive values represent the top, and negative values the bottom (Eg: \code{ylim.add = -0.5} is equivalent to \code{ylim.add = c(0.5, 0)}).
#' @param only.params Logical, default is \code{FALSE}. If \code{TRUE} no graphic is displayed, only the values of \code{x} and \code{y} used in the plot are returned.
#' @param ... Other arguments to be passed to \code{summary}, if \code{object} is an estimation, and/or to the function \code{plot} or \code{lines} (if \code{add = TRUE}).
#' @param sep The distance between two estimates -- only when argument \code{object} is a list of estimation results.
#' @param as.multiple Logical: default is \code{FALSE}. Only when \code{object} is a single estimation result: whether each coefficient should have a different color, line type, etc. By default they all get the same style.
#' @param bg Background color for the plot. By default it is white.
#' @param ci.join Logical default to \code{FALSE}. Whether to join the extremities of the confidence intervals. If \code{TRUE}, then you can set the graphical parameters with the argument \code{ci.join.par}.
#' @param ci.join.par A list of parameters to be passed to \code{\link[graphics]{lines}}. Only used if \code{ci.join=TRUE}. By default it is equal to \code{list(lwd = lwd, col = col, lty = 2)}.
#' @param ci.fill Logical default to \code{FALSE}. Whether to fille the confidence intervals with a color. If \code{TRUE}, then you can set the graphical parameters with the argument \code{ci.fill.par}.
#' @param ci.fill.par A list of parameters to be passed to \code{\link[graphics]{polygon}}. Only used if \code{ci.fill=TRUE}. By default it is equal to \code{list(col = "lightgray", alpha = 0.5)}. Note that \code{alpha} is a special parameter that adds transparency to the color (ranges from 0 to 1).
#' @param group A list, default is missing. Each element of the list reports the coefficients to be grouped while the name of the element is the group name. Each element of the list can be either: i) a character vector of length 1, ii) of length 2, or ii) a numeric vector. If equal to: i) then it is interpreted as a pattern: all element fitting the regular expression will be grouped (note that you can use the special character "^^" to clean the beginning of the names, see example), if ii) it corrsponds to the first and last elements to be grouped, if iii) it corresponds to the coefficients numbers to be grouped. If equal to a character vector, you can use a percentage to tell the algorithm to look at the coefficients before aliasing (e.g. \code{"\%varname"}). Example of valid uses: \code{group=list(group_name=\"pattern\")}, \code{group=list(group_name=c(\"var_start\", \"var_end\"))}, \code{group=list(group_name=1:2))}. See details.
#' @param group.par A list of parameters controlling the display of the group. The parameters controlling the line are: \code{lwd}, \code{tcl} (length of the tick), \code{line.adj} (adjustment of the position, default is 0), \code{tick} (whether to add the ticks), \code{lwd.ticks}, \code{col.ticks}. Then the parameters controlling the text: \code{text.adj} (adjustment of the position, default is 0), \code{text.cex}, \code{text.font}, \code{text.col}.
#' @param pt.bg The background color of the point estimate (when the \code{pt.pch} is in 21 to 25). Defaults to NULL.
#' @param lab.cex The size of the labels of the coefficients. Default is missing. It is automatically set by an internal algorithm which can go as low as \code{lab.min.cex} (another argument).
#' @param lab.min.cex The minimum size of the coefficients labels, as set by the internal algorithm. Default is 0.85.
#' @param lab.max.mar The maximum size the left margin can take when trying to fit the coefficient labels into it (only when \code{horiz = TRUE}). This is used in the internal algorithm fitting the coefficient labels. Default is \code{0.25}.
#' @param lab.fit The method to fit the coefficient labels into the plotting region (only when \code{horiz = FALSE}). Can be \code{"auto"} (the default), \code{"simple"}, \code{"multi"} or \code{"tilted"}. If \code{"simple"}, then the classic axis is drawn. If \code{"multi"}, then the coefficient labels are fit horizontally across several lines, such that they don't collide. If \code{"tilted"}, then the labels are tilted. If \code{"auto"}, an automatic choice between the three is made.
#' @param main The title of the plot. Default is \code{"Effect on __depvar__"}. You can use the special variable \code{__depvar__} to set the title (useful when you set the plot default with \code{\link[fixest]{setFixest_coefplot}}).
#' @param value.lab The label to appear on the side of the coefficient values. If \code{horiz = FALSE}, the label appears in the y-axis. If \code{horiz = TRUE}, then it appears on the x-axis. The default is equal to \code{"Estimate and __ci__ Conf. Int."}, with \code{__ci__} a special variable giving the value of the confidence interval.
#' @param xlab The label of the x-axis, default is \code{NULL}. Note that if \code{horiz = TRUE}, it overrides the value of the argument \code{value.lab}.
#' @param ylab The label of the y-axis, default is \code{NULL}. Note that if \code{horiz = FALSE}, it overrides the value of the argument \code{value.lab}.
#' @param sub A subtitle, default is \code{NULL}.
#' @param style A character scalar giving the style of the plot to be used. You can set styles with the function \code{\link[fixest]{setFixest_coefplot}}, setting all the default values of the function. If missing, then it switches to either "default" or "iplot", depending on the calling function.
#'
#' @seealso
#' See \code{\link[fixest]{setFixest_coefplot}} to set the default values of \code{coefplot}, and the estimation functions: e.g. \code{\link[fixest]{feols}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:femlm]{fenegbin}}.
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
#' # Estimation results with clustered standard-errors
#' # (the default when fixed-effects are present)
#' est_clu = summary(est)
#' # Now with "regular" standard-errors
#' est_std = summary(est, se = "standard")
#'
#' # You can plot the two results at once
#' coefplot(list(est_clu, est_std))
#'
#'
#' # Alternatively, you can use the argument x.shift
#' # to do it sequentially:
#'
#' # First graph with clustered standard-errors
#' coefplot(est, x.shift = -.2)
#'
#' # 'x.shift' was used to shift the coefficients on the left.
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
#' est_did = feols(y ~ x1 + i(period, treat, 5) | id+period, base_inter)
#'
#' # In the estimation, the variable treat is interacted
#' #  with each value of period but 5, set as a reference
#'
#' # coefplot will show all the coefficients:
#' coefplot(est_did)
#'
#' # Note that the grouping of the coefficients is due to 'group = "auto"'
#'
#' # If you want to keep only the coefficients
#' # created with i() (ie the interactions), use iplot
#' iplot(est_did)
#'
#' # When estimations contain interactions, as before,
#' #  the default behavior of coefplot changes,
#' #  it now only plots interactions:
#' coefplot(est_did)
#'
#' # We can see that the graph is different from before:
#' #  - only interactions are shown,
#' #  - the reference is present,
#' # => this is fully flexible
#'
#' iplot(est_did, ref.line = FALSE, pt.join = TRUE)
#'
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
#' est = feols(y ~ x1 + i(period_month, treat, "oct") | id+period, base_inter)
#' # Since 'period_month' of type character, coefplot sorts it
#' iplot(est)
#'
#' # To respect a plotting order, use a factor
#' base_inter$month_factor = factor(base_inter$period_month, levels = all_months)
#' est = feols(y ~ x1 + i(month_factor, treat, "oct") | id+period, base_inter)
#' iplot(est)
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
#'                 Sepal.Width + i(Species), iris)
#'
#' # And that's it
#' coefplot(est)
#'
#' # You can set separate default values for iplot
#' setFixest_coefplot("iplot", pt.join = TRUE, pt.join.par = list(lwd = 2, lty = 2))
#' iplot(est)
#'
#' # To reset to the default settings:
#' setFixest_coefplot("all", reset = TRUE)
#' coefplot(est)
#'
#' #
#' # Example 4: group + cleaning
#' #
#'
#' # You can use the argument group to group variables
#' # You can further use the special character "^^" to clean
#' #  the beginning of the coef. name: particularly useful for factors
#'
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Length +
#'                 Sepal.Width + Species, iris)
#'
#' # No grouping:
#' coefplot(est)
#'
#' # now we group by Sepal and Species
#' coefplot(est, group = list(Sepal = "Sepal", Species = "Species"))
#'
#' # now we group + clean the beginning of the names using the special character ^^
#' coefplot(est, group = list(Sepal = "^^Sepal.", Species = "^^Species"))
#'
#'
coefplot = function(object, ..., style = NULL, sd, ci_low, ci_high, x, x.shift = 0, horiz = FALSE,
                    dict = getFixest_dict(), keep, drop, order, ci.width = "1%",
                    ci_level = 0.95, add = FALSE, pt.pch = c(20, 17, 15, 21, 24, 22), pt.bg = NULL, cex = 1,
                    pt.cex = cex, col = 1:8, pt.col = col, ci.col = col, lwd = 1, pt.lwd = lwd,
                    ci.lwd = lwd, ci.lty = 1, grid = TRUE, grid.par = list(lty=3, col = "gray"),
                    zero = TRUE, zero.par = list(col="black", lwd=1), pt.join = FALSE,
                    pt.join.par = list(col = pt.col, lwd=lwd), ci.join = FALSE,
                    ci.join.par = list(lwd = lwd, col = col, lty = 2), ci.fill = FALSE,
                    ci.fill.par = list(col = "lightgray", alpha = 0.5), ref = "auto",
                    ref.line = "auto", ref.line.par = list(col = "black", lty = 2), lab.cex,
                    lab.min.cex = 0.85, lab.max.mar = 0.25, lab.fit = "auto", xlim.add,
                    ylim.add, only.params = FALSE, sep, as.multiple = FALSE,
                    bg, group = "auto", group.par = list(lwd=2, line=3, tcl=0.75),
                    main = "Effect on __depvar__", value.lab = "Estimate and __ci__ Conf. Int.",
                    ylab = NULL, xlab = NULL, sub = NULL){

    # Set up the dictionary
    if(is.null(dict)){
        dict = c()
    } else if(any(grepl("^&", names(dict)))){
        # Speficic markuo to identify coefplot aliases
        dict_amp = dict[grepl("^&", names(dict))]
        names(dict_amp) = gsub("^&", "", names(dict_amp))
        dict[names(dict_amp)] = dict_amp
    }

    check_arg_plus(lab.fit, "match(auto, simple, multi, tilted)")

    dots = list(...)

    ylab_add_ci = missing(ci_low)
    is_iplot = isTRUE(dots$internal.only.i)

    #
    # Setting the default values ####
    #

    opts = getOption("fixest_coefplot")

    check_arg(style, "NULL character scalar")
    if(is.null(style)){
        if(is_iplot){
            style = "iplot"
        } else {
            style = "default"
        }
    }

    old_style = style
    style = try(match.arg(style, choices = names(opts)), silent = TRUE)
    if("try-error" %in% class(style)){
        warning("Unknow style '", old_style, "', using default style instead.")
        style = "default"
    }

    my_opt = opts[[style]]

    if(!is.list(my_opt)){
        warning("The default values of coefplot for style '", style, "' are ill-formed and therefore reset. Use only setFixest_coefplot for setting the default values.")
        setFixest_coefplot(style = style, reset = TRUE)
    } else if(length(my_opt) > 0){

        arg_2_add = setdiff(names(opts[["default"]]), names(my_opt))
        if(length(arg_2_add) > 0){
            my_opt[arg_2_add] = opts[["default"]][arg_2_add]
            # we reorder
            arg_list = names(formals(coefplot))
            my_fact = factor(names(my_opt), levels = arg_list)
            my_opt = my_opt[base::order(my_fact)]
        }

        if(is_iplot){
            sysOrigin = sys.parent()
            mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin))
        } else {
            mc = match.call()
        }

        arg2set = setdiff(names(my_opt), names(mc))
        for(arg in arg2set){
            my_arg = my_opt[[arg]]
            if(is.name(my_arg) || is.call(my_arg)){
                if(length(my_opt$extra_values) > 0){
                    assign(arg, eval(my_arg, my_opt$extra_values))
                } else {
                    assign(arg, eval(my_arg))
                }
            } else {
                assign(arg, my_arg)
            }
        }
    }

    #
    # Getting the parameters => function coefplot_prms
    #

    info = coefplot_prms(object = object, ..., sd = sd, ci_low = ci_low, ci_high = ci_high,
                         x = x, x.shift = x.shift, dict = dict, keep = keep, drop = drop,
                         order = order, ci_level = ci_level, ref = ref, only.i = is_iplot,
                         sep = sep, as.multiple = as.multiple)

    prms = info$prms
    AXIS_AS_NUM = info$num_axis
    x_at = info$at
    x_labels = info$labels
    x_labels_raw = info$x_labels_raw
    varlist = info$varlist
    dots_drop = info$dots_drop
    my_xlim = info$xlim
    suggest_ref_line = info$suggest_ref_line
    multiple_est = info$multiple_est

    if(only.params){
        return(list(prms = prms, is_iplot = is_iplot, at = x_at, labels = x_labels))
    }

    dots = dots[!names(dots) %in% dots_drop]

    ci_low = prms$ci_low
    ci_high = prms$ci_high
    x_value = prms$x

    if(horiz){
        # We just reverse the x/y
        tmp = prms$x
        prms$x = prms$y
        prms$y = tmp
    }

    check_arg(xlab, "NULL character vector")
    if(is.null(xlab) && is_iplot){
        xlab = "__i__"
    }

    #
    # Title ####
    #

    # xlab / main / ylab / sub

    check_arg(value.lab, main, "character scalar")
    check_arg(xlab, ylab, sub, "character scalar NULL")

    if(horiz){
        if(is.null(xlab)) xlab = value.lab
    } else {
        if(is.null(ylab)) ylab = value.lab
    }

    if(is.null(xlab)) xlab = ""
    if(is.null(ylab)) ylab = ""
    if(is.null(sub)) sub = ""

    main = replace_and_make_callable(main, varlist)
    ylab = replace_and_make_callable(ylab, varlist)
    xlab = replace_and_make_callable(xlab, varlist)
    sub = replace_and_make_callable(sub, varlist)

    main = expr_builder(main)
    ylab = expr_builder(ylab)
    xlab = expr_builder(xlab)
    sub = expr_builder(sub)

    dots$main = ""
    dots$ylab = ""
    dots$xlab = ""
    dots$sub = ""

    #
    # group = auto + Height ####
    #

    # The value of group = "auto" => renaming the labels
    if(identical(group, "auto") && is_iplot == FALSE){
        # we change the names of interactions
        qui = grepl(":", x_labels_raw)
        if(any(qui)){
            x_inter = gsub(":.+", "", x_labels_raw[qui])
            tx_inter = table(x_inter)
            qui_auto = names(tx_inter)[tx_inter >= 2]

            group = list()

            for(i in seq_along(qui_auto)){
                var_left = qui_auto[i]
                qui_select = substr(x_labels_raw, 1, nchar(var_left)) == var_left

                # we require to be next to each other
                if(!all(diff(which(qui_select)) == 1)) next

                x_select = x_labels_raw[qui_select]

                is_inter = TRUE
                n_trim = 0
                group_regex = NULL
                if(all(grepl("[^:]:[^:].*::", x_select))){
                    # treat:period::1
                    first_part = strsplit(x_select[1], "::")[[1]][1]
                    var_right = gsub(".+:", "", first_part)
                    n_max = nchar(first_part) + 2

                } else if(all(grepl("::.*[^:]:[^:]", x_select))){
                    # period::1:treat
                    # Note that we always put 'treat' on the left
                    second_part = strsplit(x_select[1], "::")[[1]][2]
                    var_right = var_left
                    var_left = gsub(".+:", "", second_part)
                    n_max = nchar(var_right) + 2
                    n_trim = nchar(var_left) + 1

                    group_regex = paste0("%", escape_regex(var_right), "::.+:", escape_regex(var_left))

                } else if(all(grepl("::", x_select))) {
                    is_inter = FALSE
                    # This is a fixest factor
                    n_max = nchar(var_left) + 2

                } else {
                    # we need to find out by ourselves...
                    n = min(nchar(x_select[1:2]))
                    x_split = strsplit(substr(x_select[1:2], 1, n), "")
                    start_diff = which.max(x_split[[1]] != x_split[[2]])

                    ok = FALSE
                    for(n_max in n:(nchar(var_left) + 2)){
                        if(all(grepl(substr(x_select[1], 1, n_max), x_select, fixed = TRUE))){
                            ok = TRUE
                            break
                        }
                    }

                    if(!ok) next

                    var_right = gsub(".+:", "", substr(x_select[1], 1, n_max))
                }

                if(is_inter){
                    v_name = dict_apply(c(var_left, var_right), dict)
                    group_name = replace_and_make_callable("__x__ %*% (__y__ == ldots)", list(x = v_name[1], y = v_name[2]), text_as_expr = TRUE)

                    if(is.null(group_regex)){
                        group[[group_name]] = escape_regex(paste0("%", var_left, ":", var_right))
                    } else {
                        group[[group_name]] = group_regex
                    }


                } else {
                    v_name = dict_apply(var_left, dict)
                    group_name = replace_and_make_callable("__x__", list(x = v_name), text_as_expr = TRUE)

                    group[[group_name]] = paste0("%^", escape_regex(var_left))

                }


                # We update the labels
                x_labels[qui_select] = substr(x_select, n_max + 1, nchar(x_select) - n_trim)
            }
        }
    }

    IS_GROUP = !identical(group, "auto") && !add && !missing(group) && !is.null(group) && length(group) > 0 && !is.null(x_labels)

    line_height = par("mai")[1] / par("mar")[1]

    if(IS_GROUP){
        # This means there are groups!
        # We compute the height of the group based on the parameters
        # passed in group.par
        # => it will be used to adjust the margins

        # we need: tcl + text.adj

        tcl = 0.75
        if(!is.null(group.par$tcl)){
            tcl = max(group.par$tcl)
        }

        line.adj = 0
        if(!is.null(group.par$line.adj)){
            line.adj = max(group.par$line.adj)
        }

        text.adj = 0
        if(!is.null(group.par$text.adj)){
            text.adj = max(group.par$text.adj)
        }

        text.cex = 1
        if(!is.null(group.par$text.cex)){
            text.cex = max(group.par$text.cex)
        }

        group.height = (0.5 + tcl + text.cex + text.adj + line.adj)

        # We perform basic checks on the group
        if(!is.list(group)) stop("Argument 'group' must be a list.")

        # We just take care of the special group + cleaning case
        # when group starts with ^^

        for(i in seq_along(group)){
            my_group = group[[i]]
            if(! (length(my_group) == 1 && is.character(my_group) && substr(my_group, 1, 2) == "^^") ) next

            if(grepl("^%", my_group)){
                qui = grepl(gsub("^%", "", my_group), x_labels_raw)
            } else {
                qui = grepl(my_group, x_labels)
            }

            if(!any(qui)){
                warning("In argument 'group', the pattern: \"", my_group, "\", did not match any coefficient name.")
                group[[i]] = NULL
                if(length(group) == 0) IS_GROUP = FALSE
                next
            }

            group[[i]] = range(which(qui))
            x_labels = gsub(my_group, "", x_labels)
        }


    } else {
        group.height = 0
    }


    #
    # The limits ####
    #

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

    # ylim
    my_ylim = range(c(ci_low, ci_high))

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


    #
    # Margins and cex setting ####
    #

    if(!missing(lab.cex)){
        lab.min.cex = lab.cex
    } else {
        lab.cex = 1
    }

    if(horiz){
        # we make the labels fit into the margin

        xlab.line = 3
        sub.line = 4

        # Algorithm:
        # - if lab.cex is provided:
        #  => we fit the labels into the margin // we extend the margin until it fits
        # - if lab.cex is NOT provided:
        #  => we extend the margin so that the label fit. If the margin exceeds lab.max.mar% of the plot space, we reduce the cex
        #   etc, until it fits. If it still does not fit => we extend the margin until it fits.
        #

        nlines = group.height + 2
        # The last 2 means 1 line on each side of the label

        if(lab.cex > lab.min.cex){
            # No algorithm
            lab.width_in = max(strwidth(x_labels, units = "in", cex = lab.cex))
        } else {
            # Algorithm

            max_mar_width_in = par("pin")[1] * lab.max.mar

            while(nlines * line_height + lab.width_in > max_mar_width_in && lab.cex > lab.min.cex){
                lab.cex = 0.95 * lab.cex
                lab.width_in = max(strwidth(x_labels, units = "in", cex = lab.cex))
            }

            # Final step => if we go too far, we etend the margin, even beyond max_mar_width_in
            if(lab.cex < lab.min.cex){
                lab.cex = lab.min.cex
                lab.width_in = max(strwidth(x_labels, units = "in", cex = lab.cex))
            }

        }

        group.baseline = 2 + lab.width_in / line_height

        nlines = nlines + lab.width_in / line_height

        ylab.line = nlines

        if(ylab != ""){
            nlines = nlines + ceiling(strheight(ylab, units = "in") / line_height)
        }

        total_width = nlines * line_height

        if(total_width > par("mai")[2]){
            new_mai = par("mai")
            new_mai[2] = total_width + 0.05
            op = par(mai = new_mai)
            on.exit(par(op))
        }

        my_xlim = rev(my_xlim)

    } else {
        # We adjust the margin only if there are groups and they
        # don't fit in the original margin
        # or if there is a xlab and groups at the same time

        ylab.line = 3

        LINE_MIN_TILTED = 0.25
        LINE_MAX_TILTED = 3

        if(lab.fit == "auto"){
            # We want to display ALL labels
            in_to_usr = diff(my_xlim) / par("pin")[1]

            w_all = strwidth(x_labels, cex = lab.cex, units = "in") * in_to_usr
            em = strwidth("M", units = "in") * in_to_usr
            is_collided = ((w_all[-1] + w_all[-length(w_all)]) / 2 + em) > 1

            while(any(is_collided) && lab.cex > lab.min.cex){
                lab.cex = 0.95 * lab.cex
                w_all = strwidth(x_labels, cex = lab.cex, units = "in") * in_to_usr
                is_collided = ((w_all[-1] + w_all[-length(w_all)]) / 2 + em) > 1
            }

            # switch to multi lines
            if(any(is_collided) || lab.cex < lab.min.cex){
                lab.info = xaxis_labels(at = x_at, labels = x_labels, line.max = 1, minCex = lab.min.cex, trunc = Inf, only.params = TRUE)

                if(lab.info$failed){
                    # switch to tilted
                    lab.fit = "tilted"
                    lab.info = xaxis_biased(at = x_at, labels = x_labels, line.min = LINE_MIN_TILTED, line.max = LINE_MAX_TILTED, cex = seq(lab.cex, lab.min.cex, length.out = 4), trunc = Inf, only.params = TRUE)
                    lab.cex = lab.info$cex
                    nlines = lab.info$height_line + LINE_MIN_TILTED

                } else {
                    lab.fit = "multi"
                    nlines = lab.info$height_line
                    lab.cex = lab.info$cex
                }

            } else {
                lab.fit = "simple"
                nlines = 2 * lab.cex
            }
        } else if(lab.fit == "multi"){
            lab.info = xaxis_labels(at = x_at, labels = x_labels, line.max = 1, minCex = lab.min.cex, trunc = Inf, only.params = TRUE)
            lab.cex = lab.info$cex
            nlines = lab.info$height_line

            if(length(unique(lab.info$line)) == 1){
                lab.fit = "simple"
                nlines = 2 * lab.cex
            }

        } else if(lab.fit == "tilted"){
            lab.info = xaxis_biased(at = x_at, labels = x_labels, line.min = LINE_MIN_TILTED, line.max = LINE_MAX_TILTED, cex = seq(lab.cex, lab.min.cex, length.out = 4), trunc = Inf, only.params = TRUE)
            lab.cex = lab.info$cex
            nlines = lab.info$height_line + LINE_MIN_TILTED

        } else if(lab.fit == "simple"){
            nlines = 2 * lab.cex
        }

        # browser()

        if(IS_GROUP){
            # tcl was set before
            group.baseline = nlines + 1 - 0.75 + tcl

            nlines = nlines + group.height

            xlab.line = nlines
        } else {
            xlab.line = 3
        }

        if(xlab != ""){
            xlab.line = xlab.line + ceiling(strheight(xlab, units = "in") / line_height) - 1
        }

        sub.line = xlab.line + 1

        if(sub != ""){
            nlines = sub.line + 1
        } else if(xlab != ""){
            nlines = sub.line
        }

        total_height = nlines * line_height

        if(total_height > par("mai")[1]){
            new_mai = par("mai")
            new_mai[1] = total_height + 0.05
            op = par(mai = new_mai)
            on.exit(par(op))
        }

    }

    #
    # Plot (start) ####
    #


    all_plot_args = unique(c(names(par()), names(formals(plot.default))))
    pblm = setdiff(names(dots), all_plot_args)
    if(length(pblm) > 0){
        # warning("The following argument", ifsingle(pblm, " is not a", "s are not"), " plotting argument", ifsingle(pblm, " and is", "s and are"), " therefore ignored: ", enumerate_items(pblm), ".")
        dots[pblm] = NULL
    }

    # preparation of the do.call
    dots$col = col

    if(horiz){
        listDefault(dots, "xlim", my_ylim)
        listDefault(dots, "ylim", my_xlim)
    } else {
        listDefault(dots, "xlim", my_xlim)
        listDefault(dots, "ylim", my_ylim)
    }

    dots$x = prms$x
    dots$y = prms$y
    dots$type = "p"

    #
    # Plot Call ####
    #

    if(!add){

        dots$axes = FALSE

        # Nude graph
        first.par = dots
        first.par$type = "n"
        do.call("plot", first.par)

        # The background
        if(!missing(bg)){
            dx = diff(my_xlim)
            dy = diff(my_ylim)
            if(horiz){
                rect(xleft = my_ylim[1] - dy, ybottom = my_xlim[1] - dx, xright = my_ylim[2] + dy, ytop = my_xlim[2] + dx, col = bg)
            } else {
                rect(xleft = my_xlim[1] - dx, ybottom = my_ylim[1] - dy, xright = my_xlim[2] + dx, ytop = my_ylim[2] + dy, col = bg)
            }

        }

        # The grid
        if(grid){
            listDefault(grid.par, "col", "gray")
            listDefault(grid.par, "lty", 3)
            listDefault(grid.par, "vert", TRUE)
            listDefault(grid.par, "horiz", TRUE)

            vert = grid.par$vert
            g.horiz = grid.par$horiz
            grid.par$vert = grid.par$horiz = NULL

            if(g.horiz){
                do.call("hgrid", grid.par)
            }

            if(vert){
                do.call("vgrid", grid.par)
            }
        }

        if(zero){
            listDefault(zero.par, "lwd", 1)
            listDefault(zero.par, "col", "black")

            if(horiz){
                zero.par$v = 0
            } else {
                zero.par$h = 0
            }

            do.call("abline", zero.par)
        }

        # Tha annotations

        if(main != ""){
            title(main = main)
        }

        if(xlab != ""){
            title(xlab = xlab, line = xlab.line)
        }

        if(ylab != ""){
            title(ylab = ylab, line = ylab.line)
        }

        if(sub != ""){
            title(sub = sub, line = sub.line)
        }

        # Reference line

        check_arg(ref.line, "logical scalar | numeric vector no na | charin(auto)")

        if(identical(ref.line, "auto")){
            ref.line = suggest_ref_line && length(unique(prms[prms$is_ref, "estimate_names_raw"])) == 1
        }

        is_ref_line = !isFALSE(ref.line)
        line_direction = if(horiz) "h" else "v"
        if(is.logical(ref.line)){
            if(ref.line) {
                ref_pblm = is.null(prms$is_ref) || !any(prms$is_ref)

                if(ref_pblm && !"v" %in% names(ref.line.par)){
                    warning("You can use the argument 'ref.line = TRUE' only when a 'natural' reference is found, or if you provided a reference with argument 'ref'. You can still draw vertical lines by using 'v' in argument 'ref.line.par'. Example: ref.line.par=list(v = ", round(x_value[floor(length(x_value)/2)]), ", col=2).")
                } else if(!ref_pblm){
                    where = tapply(prms[prms$is_ref, "x"], prms[prms$is_ref, "estimate_names_raw"], mean)
                    listDefault(ref.line.par, line_direction, where)
                }
            }

        } else {
            listDefault(ref.line.par, line_direction, ref.line)
        }

        if(is_ref_line){
            listDefault(ref.line.par, "lty", 2)
            do.call("abline", ref.line.par)
        }

        box()

        if(horiz){
            axis(1)

            if(AXIS_AS_NUM){
                axis(2, las = 1)
            } else {
                if(any(grepl("^&", x_labels))){
                    # means we call expression()
                    # drawback => expressions can overlap
                    qui = grepl("^&", x_labels)
                    if(any(!qui)){
                        axis(2, at = x_at[!qui], labels = x_labels[!qui], las = 1, cex = lab.cex)
                    }

                    for(i in which(qui)){
                        axis(2, at = x_at[i], labels = expr_builder(x_labels[i]), las = 1, cex = lab.cex)
                    }

                } else {
                    # easy case: only character
                    axis(2, at = x_at, labels = x_labels, las = 1, cex = lab.cex)
                }
            }
        } else {
            axis(2)

            # browser()

            if(AXIS_AS_NUM){
                axis(1)
            } else {
                if(lab.fit == "simple"){
                    if(any(grepl("^&", x_labels))){
                        # means we call expression()
                        # drawback => expressions can overlap
                        qui = grepl("^&", x_labels)
                        if(any(!qui)){
                            axis(1, at = x_at[!qui], labels = x_labels[!qui], cex.axis = lab.cex)
                        }

                        for(i in which(qui)){
                            axis(1, at = x_at[i], labels = expr_builder(x_labels[i]), cex.axis = lab.cex)
                        }

                    } else {
                        # easy case: only character
                        axis(1, at = x_at, labels = x_labels, cex.axis = lab.cex)
                    }
                } else if(lab.fit == "multi"){

                    axis(1, at = x_at, labels = NA, tcl=-0.25)

                    for(my_line in unique(lab.info$line)){
                        qui_line = my_line == lab.info$line
                        x_labels_current = x_labels[qui_line]
                        x_at_current = x_at[qui_line]

                        if(any(grepl("^&", x_labels_current))){
                            qui = grepl("^&", x_labels_current)
                            if(any(!qui)){
                                axis(1, at = x_at_current[!qui], labels = x_labels_current[!qui], cex.axis = lab.cex, line = my_line, lwd = 0)
                            }

                            for(i in which(qui)){
                                axis(1, at = x_at_current[i], labels = expr_builder(x_labels_current[i]), cex.axis = lab.cex, line = my_line, lwd = 0)
                            }

                        } else {
                            # easy case: only character
                            axis(1, at = x_at_current, labels = x_labels_current, cex.axis = lab.cex, line = my_line, lwd = 0)
                        }
                    }
                } else if(lab.fit == "tilted"){

                    axis(1, at = x_at, labels = NA, tcl = -0.25, lwd = 0, lwd.ticks = 1)

                    if(any(grepl("^&", x_labels))){
                        qui = grepl("^&", x_labels)
                        if(any(!qui)){
                            xaxis_biased(1, at = x_at[!qui], labels = x_labels[!qui], cex = lab.cex, angle = lab.info$angle, line.min = LINE_MIN_TILTED)
                        }

                        for(i in which(qui)){
                            xaxis_biased(1, at = x_at[i], labels = expr_builder(x_labels[i]), cex = lab.cex, angle = lab.info$angle, line.min = LINE_MIN_TILTED)
                        }

                    } else {
                        # easy case: only character
                        xaxis_biased(1, at = x_at, labels = x_labels, cex = lab.cex, angle = lab.info$angle, line.min = LINE_MIN_TILTED)
                    }
                }
            }


        }



    }

    n_id = length(unique(prms$id))

    #
    # pt.join ####
    #

    check_arg(pt.join, "logical scalar")

    if(pt.join){
        # We join the dots

        listDefault(pt.join.par, "lwd", lwd)
        listDefault(pt.join.par, "col", pt.col)

        if(multiple_est){
            for(i in 1:n_id){
                my_join.par = pt.join.par

                for(param in c("col", "lwd", "lty")){
                    if(!is.null(pt.join.par[[param]])){
                        my_join.par[[param]] = par_fit(pt.join.par[[param]], i)
                    }
                }

                my_join.par$x = prms$x[prms$id == i]
                my_join.par$y = prms$y[prms$id == i]

                do.call("lines", my_join.par)
            }
        } else {
            pt.join.par$x = prms$x
            pt.join.par$y = prms$y

            do.call("lines", pt.join.par)
        }
    }


    #
    # The confidence intervals ####
    #

    x = x_value

    if(!is.na(ci.lwd) && ci.lwd > 0){

        # a) barre verticale
        if(horiz){
            segments(x0=ci_low, y0=x, x1=ci_high, y1=x, lwd = par_fit(ci.lwd, prms$id), col = par_fit(ci.col, prms$id), lty = par_fit(ci.lty, prms$id))
        } else {
            segments(x0=x, y0=ci_low, x1=x, y1=ci_high, lwd = par_fit(ci.lwd, prms$id), col = par_fit(ci.col, prms$id), lty = par_fit(ci.lty, prms$id))
        }

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
                if(horiz){
                    total_width = diff(par("usr")[3:4])
                } else {
                    total_width = diff(par("usr")[1:2])
                }

                ci.width = total_width * width_nb / 100 / 2
            } else {
                ci.width = width_nb / 2
            }
        }

        # b) toppings
        # Only if not a reference
        qui = ci_high != ci_low

        if(horiz){
            #  i) ci_high
            segments(x0=ci_high[qui], y0=x[qui]-ci.width, x1=ci_high[qui], y1=x[qui]+ci.width, lwd = par_fit(ci.lwd, prms$id[qui]), col = par_fit(ci.col, prms$id[qui]), lty = par_fit(ci.lty, prms$id[qui]))
            #  ii) ci_low
            segments(x0=ci_low[qui], y0=x[qui]-ci.width, x1=ci_low[qui], y1=x[qui]+ci.width, lwd = par_fit(ci.lwd, prms$id[qui]), col = par_fit(ci.col, prms$id[qui]), lty = par_fit(ci.lty, prms$id[qui]))

        } else {
            #  i) ci_high
            segments(x0=x[qui]-ci.width, y0=ci_high[qui], x1=x[qui]+ci.width, y1=ci_high[qui], lwd = par_fit(ci.lwd, prms$id[qui]), col = par_fit(ci.col, prms$id[qui]), lty = par_fit(ci.lty, prms$id[qui]))
            #  ii) ci_low
            segments(x0=x[qui]-ci.width, y0=ci_low[qui], x1=x[qui]+ci.width, y1=ci_low[qui], lwd = par_fit(ci.lwd, prms$id[qui]), col = par_fit(ci.col, prms$id[qui]), lty = par_fit(ci.lty, prms$id[qui]))

        }

    }


    #
    # ci.join ####
    #

    if(ci.join){
        # We join the extremities of the conf. int.

        listDefault(ci.join.par, "lwd", lwd)
        listDefault(ci.join.par, "col", col)
        listDefault(ci.join.par, "lty", 2)

        if(multiple_est){

            for(i in 1:n_id){
                my_join.par = ci.join.par

                for(param in c("col", "lwd", "lty")){
                    if(!is.null(ci.join.par[[param]])){
                        my_join.par[[param]] = par_fit(ci.join.par[[param]], i)
                    }
                }

                if(horiz){
                    my_join.par$y = prms$y[prms$id == i]

                    my_join.par$x = prms$ci_high[prms$id == i]
                    do.call("lines", my_join.par)

                    my_join.par$x = prms$ci_low[prms$id == i]
                    do.call("lines", my_join.par)

                } else {
                    my_join.par$x = prms$x[prms$id == i]

                    my_join.par$y = prms$ci_high[prms$id == i]
                    do.call("lines", my_join.par)

                    my_join.par$y = prms$ci_low[prms$id == i]
                    do.call("lines", my_join.par)
                }

            }
        } else {

            if(horiz){
                ci.join.par$y = prms$y

                ci.join.par$x = prms$ci_high
                do.call("lines", ci.join.par)

                ci.join.par$x = prms$ci_low
                do.call("lines", ci.join.par)

            } else {
                ci.join.par$x = prms$x

                ci.join.par$y = prms$ci_high
                do.call("lines", ci.join.par)

                ci.join.par$y = prms$ci_low
                do.call("lines", ci.join.par)
            }

        }
    }

    #
    # ci.fill ####
    #

    if(ci.fill){
        # We join the extremities of the conf. int.

        listDefault(ci.fill.par, "col", rgb(0.5, 0.5, 0.5))
        listDefault(ci.fill.par, "alpha", 0.5)
        listDefault(ci.fill.par, "border", NA)

        if(multiple_est){
            # We create the appropriate colors
            my_cols = par_fit(ci.fill.par$col, 1:n_id)
            my_alpha = par_fit(ci.fill.par$alpha, 1:n_id)
            ci.fill.par$alpha = NULL

            for(i in 1:n_id){
                current_col = as.vector(col2rgb(my_cols[i], alpha = TRUE)) / 255
                if(current_col[4] == 1){
                    my_cols[i] = rgb(current_col[1], current_col[2], current_col[3], my_alpha[i])
                }
            }
            ci.fill.par$col = my_cols

            for(i in 1:n_id){
                my_join.par = ci.fill.par

                for(param in c("col", "lwd", "lty")){
                    if(!is.null(ci.fill.par[[param]])){
                        my_join.par[[param]] = par_fit(ci.fill.par[[param]], i)
                    }
                }

                if(horiz){
                    my_join.par$y = c(prms$y[prms$id == i], rev(prms$y[prms$id == i]))

                    my_join.par$x = c(prms$ci_high[prms$id == i], rev(prms$ci_low[prms$id == i]))
                    do.call("polygon", my_join.par)

                } else {
                    my_join.par$x = c(prms$x[prms$id == i], rev(prms$x[prms$id == i]))

                    my_join.par$y = c(prms$ci_high[prms$id == i], rev(prms$ci_low[prms$id == i]))
                    do.call("polygon", my_join.par)
                }

            }
        } else {
            # colors => adding alpha if needed
            my_alpha = ci.fill.par$alpha[1]
            ci.fill.par$alpha = NULL
            my_col = as.vector(col2rgb(ci.fill.par$col[1], alpha = TRUE)) / 255
            if(my_col[4] == 1){
                ci.fill.par$col = rgb(my_col[1], my_col[2], my_col[3], my_alpha)
            }

            if(horiz){
                ci.fill.par$y = c(prms$y, rev(prms$y))
                ci.fill.par$x = c(prms$ci_high, rev(prms$ci_low))
                do.call("polygon", ci.fill.par)
            } else {
                ci.fill.par$x = c(prms$x, rev(prms$x))
                ci.fill.par$y = c(prms$ci_high, rev(prms$ci_low))
                do.call("polygon", ci.fill.par)
            }

        }
    }

    #
    # Point estimates ####
    #

    if(!add){
        # now the points or lines
        if(dots$type != "n"){
            point.par = dots[c("x", "y", "type", "cex", "col", "pch", "lty", "lwd")]
            point.par$pch = par_fit(pt.pch, prms$id)
            point.par$cex = par_fit(pt.cex, prms$id)
            point.par$col = par_fit(pt.col, prms$id)
            point.par$lwd = par_fit(pt.lwd, prms$id)
            if(!is.null(pt.bg)) point.par$bg = par_fit(pt.bg, prms$id)
            point.par = point.par[lengths(point.par) > 0]
            do.call("points", point.par)
        }
    } else {
        dots$pch = par_fit(pt.pch, prms$id)
        dots$cex = par_fit(pt.cex, prms$id)
        dots$col = par_fit(pt.col, prms$id)
        dots$lwd = par_fit(pt.lwd, prms$id)
        if(!is.null(pt.bg)) point.par$bg = par_fit(pt.bg, prms$id)
        do.call("points", dots)
    }

    #
    # Group ####
    #

    if(IS_GROUP){

        if(!is.list(group)) stop("Argument 'group' must be a list.")

        axis_prms_fit = c("lwd", "tcl", "line", "tick", "lty", "col", "lwd.ticks", "col.ticks")

        # The line and text adjustments
        line.adj = group.par$line.adj
        if(is.null(line.adj)){
            line.adj = 0
        }
        group.par$line.adj = NULL

        text.adj = group.par$text.adj
        if(is.null(text.adj)){
            text.adj = 0
        }
        group.par$text.adj = NULL

        # axis
        listDefault(group.par, "lwd", 2)
        listDefault(group.par, "tcl", 0.75)
        group.par$line = group.baseline + line.adj

        axis_par = group.par[!grepl("^text", names(group.par))]
        axis_par$side = 1
        axis_par$labels = NA

        # label
        listDefault(group.par, "text.line", group.par$line - 0.5 + text.adj)
        listDefault(group.par, "text.cex", 1)
        listDefault(group.par, "text.font", 1)
        listDefault(group.par, "text.col", 1)

        for(i in seq_along(group)){
            group_name = names(group)[i]
            my_group = group[[i]]

            # formatting properly
            if(is.numeric(my_group)){
                g_start = min(my_group)
                g_end = max(my_group)

                if(any(my_group %% 1 != 0)){
                    stop("The elements of the argument 'group' must be integers (e.g. group=list(group_name=1:2)). Currently they are numeric but not integers. Alternatively, you could use coefficient names (see details).")
                }

                if(g_start < 1){
                    warning("The elements of the argument 'group' must be integers ranging from 1 to the number of coefficients. The value of ", g_start, " has been replaced by 1.")
                }

                if(g_end > nrow(prms)){
                    warning("The elements of the argument 'group' must be integers ranging from 1 to the number of coefficients (here ", g_end, "). The value of ", g_end, " has been replaced by ", g_end, ".")
                }

            } else {

                if(!is.character(my_group)){
                    stop("The elements of the argument 'group' must be either: i) a character string of length 1 or 2, or ii) integers. Currently it is not character nor numeric.\nExample of valid use: group=list(group_name=\"pattern\"), group=list(group_name=c(\"var_start\", \"var_end\")), group=list(group_name=1:2))")
                }

                if(!length(my_group) %in% 1:2){
                    stop("The elements of the argument 'group' must be either: i) a character string of length 1 or 2, or ii) integers. Currently it is a character of length ", length(my_group), ".\nExample of valid use: group=list(group_name=\"pattern\"), group=list(group_name=c(\"var_start\", \"var_end\")), group=list(group_name=1:2))")
                }

                if(length(my_group) == 1){
                    # This is a pattern
                    if(grepl("^%", my_group)){
                        qui = grepl(gsub("^%", "", my_group), x_labels_raw)
                    } else {
                        qui = grepl(my_group, x_labels)
                    }

                    if(!any(qui)){
                        warning("In argument 'group', the pattern: \"", my_group, "\", did not match any coefficient name.")
                        next
                    }

                } else {
                    # This is a pattern
                    check = c(FALSE, FALSE)
                    qui = rep(FALSE, length(x_labels))
                    for(i in 1:2){
                        val = my_group[i]
                        if(grepl("^%", val)){
                            val = gsub("^%", "", val)
                            check = val %in% x_labels_raw
                            qui = qui | x_labels_raw %in% val
                        } else {
                            check = val %in% x_labels
                            qui = qui | x_labels %in% val
                        }
                    }

                    check = my_group %in% x_labels
                    if(!all(check)){
                        warning("In argument 'group', the value", enumerate_items(my_group[check], "s.quote"), ", did not match any coefficient name.")
                        next
                    }

                }

                qui_nb = which(qui)

                g_start = min(qui_nb)
                g_end = max(qui_nb)
            }

            my_axis = axis_par
            my_axis$at = x_at[c(g_start, g_end)]

            for(p in intersect(names(my_axis), axis_prms_fit)){
                my_axis[[p]] = par_fit(my_axis[[p]], i)
            }

            side = 1 + horiz
            my_axis$side = side

            info = do.call("axis", my_axis)

            if(!is.null(group_name)){
                if(grepl("^&", group_name)){
                    group_name = eval(str2lang(gsub("^&", "", group_name)))
                }

                axis(side, mean(info), labels = group_name, tick = FALSE, line = par_fit(group.par$text.line, i), cex.axis = par_fit(group.par$text.cex, i), font = par_fit(group.par$text.font, i), col.axis = par_fit(group.par$text.col, i))
            }
        }
    }

    res = list(prms=prms, is_iplot = is_iplot, at = x_at, labels = x_labels)
    return(invisible(res))
}




coefplot_prms = function(object, ..., sd, ci_low, ci_high, x, x.shift = 0, dict, keep, drop, order, ci_level = 0.95, ref = "auto", only.i = TRUE, sep, as.multiple = FALSE){

    # get the default for:
    # dict, ci.level, ref


    dots = list(...)
    is_internal = isTRUE(dots$internal__)

    varlist = list(ci = paste0(ci_level * 100, "%"))
    dots_drop = c()

    suggest_ref_line = FALSE
    multiple_est = FALSE
    NO_NAMES = FALSE
    is_iplot = only.i
    AXIS_AS_NUM = FALSE
    if(is_internal == FALSE && ((is.list(object) && class(object)[1] == "list") || "fixest_multi" %in% class(object))){
        # This is a list of estimations

        if("fixest_multi" %in% class(object)){
            # WIP:
            # - add tags of the sample
            # - add legend
            object = attr(object, "data")
        }

        #
        # Multiple estimations ####
        #

        multiple_est = TRUE

        mc = match.call()
        mc$object = as.name("my__object__")
        mc$only.params = TRUE
        mc$internal__ = TRUE

        nb_est = length(object)

        rerun = FALSE
        first = TRUE

        while(first || rerun){
            first = FALSE

            res = varlist = list()
            all_inter = c()
            all_inter_root = c()
            dots_drop = c()
            num_axis = TRUE
            suggest_ref_line = TRUE
            for(i in 1:nb_est){
                # cat("Eval =", i, "\n")
                my__object__ = object[[i]]
                prms = try(eval(mc), silent = TRUE)

                if("try-error" %in% class(prms)){
                    stop("The ", n_th(i), " element of 'object' raises and error:\n", prms)
                }

                # Some meta variables
                varlist$ci = unique(prms$varlist$ci)
                varlist$depvar = unique(prms$varlist$depvar)
                varlist$var = unique(prms$varlist$var)
                varlist$fe = unique(prms$varlist$fe)
                varlist$i = unique(prms$varlist$i)

                dots_drop = unique(c(dots_drop, prms$dots_drop))
                suggest_ref_line = suggest_ref_line && prms$suggest_ref_line

                num_axis = num_axis && prms$num_axis

                df_prms = prms$prms
                df_prms$est_nb = i
                res[[i]] = df_prms
            }
        }

        AXIS_AS_NUM = num_axis

        all_estimates = do.call("rbind", res)

        # We respect the order provided by the user
        my_names_order = unique(all_estimates$estimate_names)
        my_order = 1:length(my_names_order)
        names(my_order) = my_names_order
        all_estimates$id_order = my_order[as.character(all_estimates$estimate_names)]
        all_estimates = all_estimates[base::order(all_estimates$id_order, all_estimates$est_nb), ]

        # we rescale -- beware of multiple est whose x is numeric!!!

        if(nb_est > 1){
            if(missing(sep)){
                all_sep = c(0.2, 0.2, 0.18, 0.16)
                if(length(all_sep) < nb_est - 1){
                    sep = 1 / (nb_est-1) * 0.7
                } else {
                    sep = all_sep[nb_est - 1]
                }

            } else {
                n_sep = length(sep)
                if(n_sep > 1){
                    if(n_sep < nb_est - 1){
                        sep = sep[n_sep]
                    } else {
                        sep = sep[nb_est - 1]
                    }
                }
            }

            if(AXIS_AS_NUM){
                all_estimates$x_new = all_estimates$x + ((all_estimates$est_nb - 1) / (nb_est - 1) - 0.5) * ((nb_est - 1) * sep)
            } else {
                all_estimates$x_new = all_estimates$id_order + ((all_estimates$est_nb - 1) / (nb_est - 1) - 0.5) * ((nb_est - 1) * sep)
            }

        } else {
            sep = 0
            all_estimates$x_new = all_estimates$id_order
        }

        # The coefficients

        estimate = all_estimates$y
        names(estimate) = all_estimates$estimate_names
        ci_low = all_estimates$ci_low
        ci_high = all_estimates$ci_high
        x = all_estimates$x_new

        prms = data.frame(estimate = estimate, ci_low = ci_low, ci_high = ci_high, x = x, id = all_estimates$est_nb, estimate_names = all_estimates$estimate_names, estimate_names_raw = all_estimates$estimate_names_raw)
        if(!is.null(all_estimates$is_ref)){
            prms$is_ref = all_estimates$is_ref
        }

    } else {

        #
        # Single estimation ####
        #


        if(is.list(object)){

            sum_exists = FALSE
            for(c_name in class(object)){
                if(exists(paste0("summary.", c_name), mode = "function")){
                    sum_exists = TRUE
                    break
                }
            }

            if(!sum_exists){
                # stop("There is no summary method for objects of class ", c_name, ". 'coefplot' applies summary to the object to extract the coeftable. Maybe add directly the coeftable in object instead?")
                mat = coeftable(object)
            } else {
                fun_name = paste0("summary.", c_name)
                args_name_sum = names(formals(fun_name))
                args_sum = intersect(names(dots), args_name_sum)

                dots_drop = args_sum

                # we kick out the summary arguments from dots
                dots[args_sum] = NULL

                # We reconstruct a call to coeftable
                mc_coeftable = match.call(expand.dots = TRUE)
                mc_coeftable[[1]] = as.name("coeftable")
                mc_coeftable[setdiff(names(mc_coeftable), c(args_sum, "object", ""))] = NULL

                mat = eval(mc_coeftable, parent.frame())
            }

            sd = mat[, 2]
            estimate = mat[, 1]

            names(estimate) = rownames(mat)

            if("fml" %in% names(object)){
                depvar = gsub(" ", "", as.character(object$fml)[[2]])
                if(depvar %in% names(dict)) depvar = dict[depvar]
                varlist$depvar = depvar
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

        n = length(estimate)

        if(missing(sd)){
            if(missing(ci_low) || missing(ci_high)) stop("If 'sd' is not provided, you must provide the arguments 'ci_low' and 'ci_high'.")

            varlist$ci = NULL
        } else {
            if(!missing(ci_low) || !missing(ci_high)) warning("Since 'sd' is provided, arguments 'ci_low' or 'ci_high' are ignored.")

            # We compute the CI
            nb = abs(qnorm((1 - ci_level)/2))
            ci_high = estimate + nb*sd
            ci_low = estimate - nb*sd
        }

        #
        # iplot ####
        #

        ref_id = NA
        xlab_suggest = NULL
        if(is_iplot){
            if(is.null(names(estimate))){
                stop("'iplot' must be used only with fixest objects containing variables created with i(). Currently it does not seem to be the case.")
            }

            all_vars = names(estimate)

            if(!any(grepl("::", all_vars))){
                stop("'iplot' must be used only with fixest objects containing variables created with i(). Currently it does not seem to be the case.")
            }

            # Four cases:
            # - factor_var::value
            # - factor_var::value:xnum
            # - xnum:factor_var::value
            # - factor_var::value:xfact::value

            # Restriction:
            # it only accepts "pure" i() variables

            # We can handle only case 1, 2, and 3
            # case 4 is too messy
            # case 4: multiple lines + legend ?

            # We take the first i() in the list
            # after having applied keep_apply

            all_vars = keep_apply(all_vars, keep)

            mm_info = object$model_matrix_info

            # Finding out which to display
            ok = FALSE
            for(i in seq_along(mm_info)){
                info = mm_info[[i]]
                if(isFALSE(info$is_inter_fact) && any(info$coef_names %in% all_vars)){
                    # That's the one
                    ok = TRUE
                    break
                }
            }

            if(!ok){
                # Now we look at the cause
                msg = if(!missnull(keep)) "reshape your 'keep' argument?" else "note that this function only works with i() variables (which should not be interacted with any other variable)."
                stop("No variable was selected: ", msg)
            }

            # "keep" here works differently => new arg. i.select?

            ANY_AUTO_REF = length(info$ref_id) > 0
            SHOW_REF = (identical(ref, "auto") || isTRUE(ref) || identical(ref, "all")) && ANY_AUTO_REF
            SHOW_REF_FIRST = SHOW_REF && !identical(ref, "all")

            # Global variables for fill_coef function
            names_coef = names(estimate)
            names_all = info$coef_names_full
            new_names = info$items

            is_rm = FALSE
            ID_rm = NULL
            if(ANY_AUTO_REF){
                if(SHOW_REF_FIRST){
                    ID_rm = info$ref_id[-1]
                } else if(SHOW_REF == FALSE){
                    ID_rm = info$ref_id
                }
                is_rm = length(ID_rm) > 0
            }

            if(is_rm){
                names_all = names_all[-ID_rm]
                new_names = new_names[-ID_rm]
            }

            fill_coef = function(coef){
                # we get the vector right
                res = rep(0, length(names_all))
                names(res) = names_all

                # we need it for CI
                names(coef) = names_coef
                inter_names = intersect(names_all, names_coef)
                res[inter_names] = coef[inter_names]

                names(res) = new_names
                res
            }

            estimate = fill_coef(estimate)
            ci_high = fill_coef(ci_high)
            ci_low = fill_coef(ci_low)
            estimate_names = new_names
            estimate_names_raw = names_all

            if(isTRUE(info$is_num) && missing(x)){
                AXIS_AS_NUM = TRUE
                names(estimate) = NULL
                x = info$items
                if(is_rm) x = x[-ID_rm]
            }

            # ref
            if(SHOW_REF){
                suggest_ref_line = length(info$ref_id) == 1 && info$is_inter_num
                is_ref = seq_along(estimate) == info$ref_id[1]

            } else {
                is_ref = rep(FALSE, length(estimate))
            }

            n = length(estimate)

            varlist$i = dict_apply(gsub("::.*", "", names_all[1]), dict)
        }

        # The DF of all the parameters
        prms = data.frame(estimate = estimate, ci_low = ci_low, ci_high = ci_high)
        if(is_iplot){
            prms$estimate_names = estimate_names
            prms$estimate_names_raw = estimate_names_raw
            prms$is_ref = is_ref
        } else if(!is.null(names(estimate))){
            prms$estimate_names = names(estimate)
            prms$estimate_names_raw = names(estimate)
        } else {
            NO_NAMES = TRUE
            prms$estimate_names = paste0("c", 1:nrow(prms))
            prms$estimate_names_raw = prms$estimate_names
        }

        # setting the names of the estimate
        if(!missing(x)){
            if(length(x) != n) stop("Argument 'x' must have the same length as the number of coefficients (currently ", length(x), " vs ", n, ").")

            if(!is.numeric(x)){
                names(estimate) = x
                prms$estimate_names = names(estimate)
            } else if(NO_NAMES) {
                AXIS_AS_NUM = TRUE
            }

            prms$x = x
        }

        # We add the reference
        if(!(identical(ref, "auto") || identical(ref, "all")) && length(ref) > 0 && !isFALSE(ref)){

            if(AXIS_AS_NUM){
                if(!is.numeric(ref) || length(ref) > 1){
                    check_arg(ref, "numeric scalar")
                } else {
                    names(ref) = "reference"
                }

            } else {
                if(is.null(names(ref))){
                    if(!is.character(ref) || length(ref) > 1){
                        check_arg(ref, "character scalar", .message = "Argument 'ref' must be either: a single character, either a list or a named integer vector of length 1 (The integer gives the position of the reference among the coefficients).")
                    } else {
                        refname = ref
                        ref = list()
                        ref[[refname]] = 1
                    }
                }

                ref = unlist(ref)

                if(!isScalar(ref, int = TRUE)){
                    reason = ifelse(length(ref) == 1, " an integer", " of length 1")
                    stop("Argument 'ref' must be either: a single character, either a list or a named integer vector of length 1. The integer gives the position of the reference among the coefficients. Currently this is not ", reason, ".")
                }
            }

            # we recreate the parameters
            n = nrow(prms)
            prms$is_ref = FALSE
            ref_row = data.frame(estimate = 0, ci_low = 0, ci_high = 0, estimate_names = names(ref), estimate_names_raw = names(ref), is_ref = TRUE)
            if(AXIS_AS_NUM){
                ref_row$x = unname(ref)
                prms = rbind(prms, ref_row)
                prms = prms[base::order(prms$x), ]
                x = prms$x

            } else {
                prms = rbind(prms, ref_row)
                if(ref > n) ref = n + 1
                ids = 1:n
                ids[ids >= ref] = ids[ids >= ref] + 1
                prms = prms[base::order(c(ids, ref)), ]
            }

        }
    }

    n = nrow(prms)

    #
    # order/drop/dict ####
    #

    if(!AXIS_AS_NUM && !is.null(prms$estimate_names)){

        # dict
        if(missnull(dict)){
            dict = c()
        } else {
            prms$estimate_names = dict_apply(prms$estimate_names, dict)
        }

        # dropping some coefs
        all_vars = unique(prms$estimate_names)

        if(!missing(keep) && length(keep) > 0){
            all_vars = keep_apply(all_vars, keep)

            if(length(all_vars) == 0){
                stop("Argument 'keep' has removed all variables!")
            }

            prms = prms[prms$estimate_names %in% all_vars,]
        }

        if(!missing(drop) && length(drop) > 0){
            all_vars = drop_apply(all_vars, drop)

            if(length(all_vars) == 0){
                stop("Argument 'drop' has removed all variables!")
            }

            prms = prms[prms$estimate_names %in% all_vars,]
        }

        # ordering the coefs
        if(!missing(order) && length(order) > 0){
            all_vars = order_apply(all_vars, order)

            my_order = 1:length(all_vars)
            names(my_order) = all_vars
            prms$id_order = my_order[prms$estimate_names]

            prms = prms[base::order(prms$id_order), ]
        }

        estimate = prms$estimate
        names(estimate) = prms$estimate_names
        ci_high = prms$ci_high
        ci_low = prms$ci_low
        if(!is.null(prms$x)) x = prms$x

        n = nrow(prms)
    }

    #
    # we create x_labels, x_value & x_at
    #

    # id: used for colors/lty etc
    if(!multiple_est){
        if(as.multiple){
            prms$id = 1:nrow(prms)
        } else {
            prms$id = 1
        }
    }

    if(multiple_est){
        # We don't allow x.shift

        my_xlim = range(x) + c(-1, 1) * ((nb_est - 1) * sep)
        x_value = x

        # We allow identical aliases to display properly (they can be identified with 'group')
        quoi = unique(prms[, c("estimate_names", "estimate_names_raw")])

        x_labels = quoi$estimate_names
        x_labels_raw = quoi$estimate_names_raw
        x_at = 1:length(x_labels)

    } else if(!missing(x) && is.numeric(x)){
        my_xlim = range(c(x + x.shift, x - x.shift))

        x_value = x + x.shift

        if(NO_NAMES){
            x_at = NULL
            x_labels = x_labels_raw = NULL
        } else {
            x_at = x
            x_labels = prms$estimate_names
            x_labels_raw = prms$estimate_names_raw
        }

    } else {
        x_at = 1:n
        x_value = 1:n + x.shift

        if(NO_NAMES){
            x_labels = x_labels_raw = 1:n
        } else {
            x_labels = prms$estimate_names
            x_labels_raw = prms$estimate_names_raw
        }

        my_xlim = range(c(1:n + x.shift, 1:n - x.shift)) + c(-0.5, +0.5)
    }

    prms$x = x_value
    prms$y = prms$estimate

    return(list(prms = prms, num_axis = AXIS_AS_NUM, at = x_at, labels = x_labels, x_labels_raw = x_labels_raw, varlist = varlist, dots_drop = dots_drop, xlim = my_xlim, suggest_ref_line = suggest_ref_line, multiple_est = multiple_est))
}



replace_and_make_callable = function(text, varlist, text_as_expr = FALSE){
    # used to make a character string callable and substitutes variables inside text with the variables
    # in varlist
    # ex: "Interacted with __var__" becomes "Interacted with x_beta"
    # or: "&paste(\"Interacted with \", x[beta])"

    if(length(text) > 1) stop("Internal problem: length of text should not be greater than 1.")

    text_split = strsplit(paste0(text, " "), "__")[[1]]

    if(length(text_split) < 3){
        # Nothing to be done!
        return(text)
    } else {
        # We need to replace the variables
        is_var = seq_along(text_split) %% 2 == 0
        my_variables = text_split[is_var]

        if(length(varlist) == 0 || any(!my_variables %in% names(varlist))){
            info = "No special variable is available for this estimation."
            if(length(varlist) > 0){
                info = paste0("In this estimation, the only special variable", enumerate_items(paste0("__", names(varlist), "__"), "s.is.start"), ". ")
            }

            # warning(info, enumerate_items(paste0("__", setdiff(my_variables, names(varlist)), "__"), "is"), " not valid, thus ignored.", call. = FALSE)

            return("")

            not_var = !my_variables %in% names(varlist)
            is_var[is_var][not_var] = FALSE
            my_variables = intersect(my_variables, names(varlist))
            if(length(my_variables) == 0){
                return(text)
            }
        }

        my_variables_values = varlist[my_variables]

        if(any(lengths(varlist[my_variables]) > 1)){
            qui = which(lengths(varlist[my_variables]) > 1)[1]
            n_val = lengths(varlist[my_variables])[qui]
            warning("The special variable __", my_variables[qui], "__ takes ", n_val, " values, only the first is used.", call. = FALSE)
            my_variables_values = sapply(my_variables_values, function(x) x[1])
        } else {
            my_variables_values = unlist(my_variables_values)
        }

        # we prepare the result (we drop the last space)
        n = length(text_split)
        text_new = text_split
        text_new[n] = gsub(" $", "", text_new[n])
        if(nchar(text_new[n]) == 0){
            text_new = text_new[-n]
            is_var = is_var[-n]
        }


        # Do the variables contain expressions?
        is_expr = grepl("^&", my_variables_values)
        if(any(is_expr)){

            expr_drop = function(x){
                if(grepl("^&", x)){
                    res = gsub("^&", "", x)
                    if(grepl("^expression\\(", res)){
                        res = gsub("(^expression\\()|(\\)$)", "", res)
                    } else if(grepl("^substitute\\(", res)){
                        res = deparse(eval(str2lang(res)))
                    }
                } else {
                    res = x
                }
                res
            }

            my_vars = sapply(my_variables_values, expr_drop)

            if(text_as_expr){
                text_new[is_var][is_expr] = my_vars[is_expr]

                if(all(is_expr)){
                    text_new = paste0("&expression(", paste(text_new, collapse = " "), ")")
                } else {
                    my_var_no_expr = paste0('"', my_vars, '"')[!is_expr]
                    new_names = paste0("x___", seq_along(my_var_no_expr))
                    text_new[is_var][!is_expr] = new_names
                    text_new = paste0("&substitute(", paste(text_new, collapse = " "), ", list(", paste0(new_names, "=", my_var_no_expr, collapse = ", "), "))")
                }
            } else {
                text_new = paste0('"', text_new, '"')
                text_new[is_var][is_expr] = my_vars[is_expr]

                if(all(is_expr)){
                    text_new = paste0("&expression(paste(", paste(text_new, collapse = ", "), "))")
                } else {
                    my_var_no_expr = paste0('"', my_vars, '"')[!is_expr]
                    new_names = paste0("x___", seq_along(my_var_no_expr))
                    text_new[is_var][!is_expr] = new_names
                    text_new = paste0("&substitute(paste(", paste(text_new, collapse = ", "), "), list(", paste0(new_names, "=", my_var_no_expr, collapse = ", "), "))")
                }

            }

            return(text_new)
        } else {
            # They don't contain expressions => fine, we just replace with the variables
            if(text_as_expr){
                my_vars = paste0('"', my_variables_values, '"')
                new_names = paste0("x___", seq_along(my_vars))
                text_new[is_var] = new_names
                text_new = paste0("&substitute(", paste(text_new, collapse = ""), ", list(", paste0(new_names, "=", my_vars, collapse = ", "), "))")
            } else {
                text_new[is_var] = my_variables_values
                text_new = paste(text_new, collapse = "")
            }

            return(text_new)
        }

    }
}

expr_builder = function(x){

    if(grepl("^&", x)){
        my_expr = gsub("^&", "", x)
        if(grepl("^(expression|substitute)\\(", my_expr)){
            # direct evaluation
            res = eval(str2lang(my_expr), parent.frame(2))
        } else {
            res = eval(str2lang(paste0("expression(", my_expr, ")")))
        }
    } else {
        res = x
    }

    res
}

####
#### iplot ####
####


gen_iplot = function(){
    # iplot has the same arguments as coefplot
    # we make all changes in coefplot
    # I automatically generate the iplot function, matching all coefplot arguments

    coefplot_args = formals(coefplot)

    arg_name = names(coefplot_args)
    arg_default = sapply(coefplot_args, deparse_long)

    #
    # iplot
    #

    qui_keep = !arg_name %in% c("object", "...")

    iplot_args = paste0(arg_name[qui_keep], " = ", arg_default[qui_keep], collapse = ", ")
    iplot_args = gsub(" = ,", ",", iplot_args)

    coefplot_call = paste0(arg_name[qui_keep], " = ", arg_name[qui_keep], collapse = ", ")

    iplot_fun = paste0("iplot = function(object, ..., ", iplot_args, "){\n\n",
                        "\tcoefplot(object = object, ..., ", coefplot_call, ", internal.only.i = TRUE)\n}")

    iplot_rox = "#' @describeIn coefplot Plots the coefficients generated with i()"

    # Writing the functions

    f = file("R/iplot.R", "w", encoding = "utf-8")

    intro = c("# Do not edit by hand\n# => iplot calls coefplot internally\n\n\n")

    s = "\n\n\n\n"
    text = c(intro, s, iplot_rox, iplot_fun, s)
    writeLines(text, f)
    close(f)

}



####
#### Default setting ####
####

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
setFixest_coefplot = function(style, horiz = FALSE, dict = getFixest_dict(), keep,
                              ci.width = "1%", ci_level = 0.95, pt.pch = 20, pt.bg = NULL,
                              cex = 1, pt.cex = cex, col = 1:8, pt.col = col, ci.col = col,
                              lwd = 1, pt.lwd = lwd, ci.lwd = lwd, ci.lty = 1, grid = TRUE,
                              grid.par = list(lty = 3, col = "gray"), zero = TRUE,
                              zero.par = list(col = "black", lwd = 1), pt.join = FALSE,
                              pt.join.par = list(col = pt.col, lwd = lwd), ci.join = FALSE,
                              ci.join.par = list(lwd = lwd, col = col, lty = 2), ci.fill = FALSE,
                              ci.fill.par = list(col = "lightgray", alpha = 0.5), ref.line = "auto",
                              ref.line.par = list(col = "black", lty = 2), lab.cex, lab.min.cex = 0.85,
                              lab.max.mar = 0.25, lab.fit = "auto", xlim.add, ylim.add, sep, bg,
                              group = "auto", group.par = list(lwd = 2, line = 3, tcl = 0.75),
                              main = "Effect on __depvar__", value.lab = "Estimate and __ci__ Conf. Int.",
                              ylab = NULL, xlab = NULL, sub = NULL, reset = FALSE){

    fm_cp = formals(coefplot)
    arg_list = names(fm_cp)
    # arg_no_default = c("object", "sd", "ci_low", "ci_high", "drop", "order", "ref", "add", "only.params", "as.multiple", "...", "x", "x.shift")
    # m = fm_cp[!names(fm_cp) %in% arg_no_default]
    # cat(gsub(" = ,", ",", paste0(names(m), " = ", sapply(m, deparse), collapse = ", ")))

    iplot_default = list()

    #
    # Controls
    #

    check_arg(style, "character scalar")
    check_arg(ci.width, "scalar(numeric, character) GE{0}")
    check_arg(ci_level, "numeric scalar GT{0} LT{1}")
    check_arg(lwd, ci.lwd, "numeric scalar GE{0}")
    check_arg(grid, zero, "logical scalar")

    check_arg_plus("L0 list NULL{list()}", grid.par, zero.par, pt.join.par, ref.line.par)

    check_arg(reset, "logical scalar")

    #
    # Code
    #

    if(missing(style)){
        style = "default"
    }

    if(style == "all" && reset){
        opts = list(default = list(), iplot = iplot_default)
        options("fixest_coefplot" = opts)
        return(invisible(NULL))
    }

    opts = getOption("fixest_coefplot")
    if(!is.list(opts)){
        warning("Wrong format of getOption('fixest_coefplot'), the options of coefplot are reset.")
        opts = list(default = list(), iplot = iplot_default)
    }

    if(reset){
        if(style == "iplot"){
            my_opt = iplot_default
        } else {
            my_opt = list()
        }

    } else {
        my_opt = opts[[style]]

        if(is.null(my_opt)){
            my_opt = list()
        } else if(!is.list(opts)){
            warning("Wrong format of getOption('fixest_coefplot') for style '", style, "', the options of coefplot for this style are reset.")
            my_opt = list()
        }
    }

    mc = match.call()

    all_args = setdiff(names(mc), c("", "reset"))
    mc = mc[all_args]

    # now we find which args for which we delay evaluation
    arg_var_default = lapply(fm_cp[!names(fm_cp) %in% all_args], all.vars)
    arg_var_default = arg_var_default[lengths(arg_var_default) > 0]
    default.eval = sapply(arg_var_default, function(x) any(x %in% all_args))
    if(any(default.eval)){
        default.arg = names(default.eval)[default.eval]
        mc[default.arg] = fm_cp[default.arg]
        # we re order
        my_fact = factor(names(mc), levels = arg_list)
        mc = mc[order(my_fact)]
        all_args = names(mc)
    }

    extra_values = my_opt$extra_values
    if(!is.list(extra_values)) extra_values = list()

    for(arg in all_args){
        my_arg = mc[[arg]]

        my_arg_vars = all.vars(my_arg)
        if(length(my_arg_vars) == 0 || !(any(my_arg_vars %in% setdiff(arg_list, "dict")))){
            if("par" %in% all.names(my_arg)){
                my_opt[[arg]] = my_arg
            } else {
                my_opt[[arg]] = eval(my_arg)
            }

        } else {
            my_extra_args = setdiff(my_arg_vars, arg_list)
            if(length(my_extra_args) > 0){
                for(v in my_extra_args) extra_values[[v]] = eval(str2lang(v), parent.frame())
            }
            # control that the order is proper
            arg2eval = intersect(my_arg_vars, arg_list)
            order_arg = which(arg_list == arg)
            order_arg2eval = sapply(arg2eval, function(x) which(arg_list == x))
            if(any(order_arg2eval > order_arg)){
                qui = which(order_arg2eval > order_arg)[1]
                stop("In argument ", arg, ": its value (", deparse(my_arg), ") is set with the argument ", arg2eval[qui], " which occurs after. If you set an argument with the value of another argument: you must  use arguments appearing before.")
            }
            my_opt[[arg]] = my_arg
        }
    }

    if(length(extra_values) > 0) my_opt$extra_values = extra_values

    opts[[style]] = my_opt

    options("fixest_coefplot" = opts)
}

#' @rdname setFixest_coefplot
getFixest_coefplot = function(){
    getOption("fixest_coefplot")
}































