#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Feb 10 09:46:01 2020
# ~: Simple function to display results
# from multiple estimations
#----------------------------------------------#



#' Plots confidence intervals
#'
#' This function plots the results of estimations (coefficients and confidence intervals). It is flexible and handles interactions in a special way.
#'
#' @inheritParams etable
#'
#' @param object Can be either: i) an estimation object (obtained for example from \code{\link[fixest]{feols}}, ii) a list of estimation objects (several results will be plotted at once), iii) a matrix of coefficients table, iv) a numeric vector of the point estimates -- the latter requiring the extra arguments \code{sd} or \code{ci_low} and \code{ci_high}.
#' @param sd The standard errors of the estimates. It may be missing.
#' @param ci_low If \code{sd} is not provided, the lower bound of the confidence interval. For each estimate.
#' @param ci_high If \code{sd} is not provided, the upper bound of the confidence interval. For each estimate.
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
#' @param pt.join Logical, default depends on the situation. If \code{TRUE}, then the coefficient estimates are joined with a line. By default, it is equal to \code{TRUE} only if: i) interactions are plotted, ii) the x values are numeric and iii) a reference is found.
#' @param pt.join.par List. Parameters of the line joining the cofficients. The default values are: \code{col = pt.col} and \code{lwd = lwd}. You can add any graphical parameter that will be passed to \code{\link[graphics]{lines}}. Eg: \code{pt.join.par = list(lty = 2)}.
#' @param ref.line Logical, default depends on the situation. It is \code{TRUE} only if: i) interactions are plotted, ii) the x values are numeric and iii) a reference is found. If \code{TRUE}, then a vertical line is drawn at the level of the reference value. You can set the parameters of this line with the argument \code{ref.line.par}.
#' @param ref.line.par List. Parameters of the vertical line on the reference. The default values are: \code{col = "black"} and \code{lty = 2}. You can add any graphical parameter that will be passed to \code{\link[graphics]{abline}}. Eg: \code{ref.line.par = list(lty = 1, lwd = 3)}.
#' @param xlim.add A numeric vector of length 1 or 2. It represents an extension factor of xlim, in percentage. Eg: \code{xlim.add = c(0, 0.5)} extends \code{xlim} of 50\% on the right. If of lentgh 1, positive values represent the right, and negative values the left (Eg: \code{xlim.add = -0.5} is equivalent to \code{xlim.add = c(0.5, 0)}).
#' @param ylim.add A numeric vector of length 1 or 2. It represents an extension factor of ylim, in percentage. Eg: \code{ylim.add = c(0, 0.5)} extends \code{ylim} of 50\% on the top. If of lentgh 1, positive values represent the top, and negative values the bottom (Eg: \code{ylim.add = -0.5} is equivalent to \code{ylim.add = c(0.5, 0)}).
#' @param only.params Logical, default is \code{FALSE}. If \code{TRUE} no graphic is displayed, only the values of \code{x} and \code{y} used in the plot are returned.
#' @param only.inter Logical, default is \code{TRUE}. If an interaction of the type of \code{var::fe} (see \code{\link[fixest]{feols}} help for details) is found, then only these interactions are plotted. If \code{FALSE}, then interactions are treated as regular coefficients.
#' @param ... Other arguments to be passed to \code{summary}, if \code{object} is an estimation, and/or to the function \code{plot} or \code{lines} (if \code{add = TRUE}).
#' @param sep The distance between two estimates -- only when argument \code{object} is a list of estimation results.
#' @param as.multiple Logical: default is \code{FALSE}. Only when \code{object} is a single estimation result: whether each coefficient should have a different color, line type, etc. By default they all get the same style.
#' @param bg Background color for the plot. By default it is white.
#' @param ci.join Logical default to \code{FALSE}. Whether to join the extremities of the confidence intervals. If \code{TRUE}, then you can set the graphical parameters with the argument \code{ci.join.par}.
#' @param ci.join.par A list of parameters to be passed to \code{\link[graphics]{lines}}. Only used if \code{ci.join=TRUE}. By default it is equal to \code{list(lwd = lwd, col = col, lty = 2)}.
#' @param ci.fill Logical default to \code{FALSE}. Whether to fille the confidence intervals with a color. If \code{TRUE}, then you can set the graphical parameters with the argument \code{ci.fill.par}.
#' @param ci.fill.par A list of parameters to be passed to \code{\link[graphics]{polygon}}. Only used if \code{ci.fill=TRUE}. By default it is equal to \code{list(col = "lightgray", alpha = 0.5)}. Note that \code{alpha} is a special parameter that adds transparency to the color (ranges from 0 to 1).
#' @param group A list, default is missing. Each element of the list reports the coefficients to be grouped while the name of the element is the group name. Each element of the list can be either: i) a character vector of length 1, ii) of length 2, or ii) a numeric vector. If equal to: i) then it is interpreted as a pattern: all element fitting the regular expression will be grouped, if ii) it corrsponds to the first and last elements to be grouped, if iii) it corresponds to the coefficients numbers to be grouped. If equal to a character vector, you can use a percentage to tell the algorithm to look at the coefficients before aliasing (e.g. \code{"%varname"}). Example of valid uses: \code{group=list(group_name=\"pattern\")}, \code{group=list(group_name=c(\"var_start\", \"var_end\"))}, \code{group=list(group_name=1:2))}. See details.
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
#' coefplot(est_did, ref.line = FALSE, pt.join = FALSE)
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
coefplot = function(object, ..., style, sd, ci_low, ci_high, x, x.shift = 0, dict, drop, order, ci.width="1%", ci_level = 0.95, add = FALSE, pt.pch = 20, cex = par("cex"), pt.cex = cex, col = 1:8, pt.col = col, ci.col = col, lwd = par("lwd"), ci.lwd = lwd, ci.lty = 1, grid = TRUE, grid.par = list(lty=3, col = "gray"), zero = TRUE, zero.par = list(col="black", lwd=1), pt.join = FALSE, pt.join.par = list(col = pt.col, lwd=lwd), ci.join = FALSE, ci.join.par = list(lwd = lwd, col = col, lty = 2), ci.fill = FALSE, ci.fill.par = list(col = "lightgray", alpha = 0.5), ref = "&auto", ref.line = "auto", ref.line.par = list(col = "black", lty = 2), xlim.add, ylim.add, only.params = FALSE, only.inter = TRUE, sep, as.multiple = FALSE, bg, group = "auto", group.par = list(lwd=2, line=3, tcl=0.75, lab.line=2.5, lab.cex = 1), main = "Effect on __depvar__", ylab = "Estimate and __ci__ Conf. Int.", xlab = "", sub = ""){

    if(missing(dict)) dict = c()

    dots = list(...)

    #
    # We get the default values
    #

    ylab_add_ci = missing(ci_low)

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
    # Getting the parameters => function coefplot_prms
    #

    info = coefplot_prms(object = object, ..., sd = sd, ci_low = ci_low, ci_high = ci_high, x = x, x.shift = x.shift, dict = dict, drop = drop, order = order, ci_level = ci_level, ref = ref, only.inter = only.inter, sep = sep, as.multiple = as.multiple)

    prms = info$prms
    is_interaction = info$is_interaction
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
        return(list(prms=prms, is_interaction = is_interaction, at = x_at, labels = x_labels))
    }

    dots = dots[!names(dots) %in% dots_drop]

    ci_low = prms$ci_low
    ci_high = prms$ci_high
    x_value = prms$x

    #
    # Title ####
    #

    # xlab / main / ylab / sub

    check_arg(xlab, "singleCharacter")
    check_arg(main, "singleCharacter")
    check_arg(ylab, "singleCharacter")
    check_arg(sub, "singleCharacter")

    main = replace_and_make_callable(main, varlist)
    ylab = replace_and_make_callable(ylab, varlist)
    xlab = replace_and_make_callable(xlab, varlist)
    sub = replace_and_make_callable(sub, varlist)

    dots$main = expr_builder(main)
    dots$ylab = expr_builder(ylab)
    dots$xlab = expr_builder(xlab)
    dots$sub = expr_builder(sub)

    #
    # group = auto ####
    #

    # The value of group = "auto" => renaming the labels
    if(identical(group, "auto") && is_interaction == FALSE){
        # we change the names of interactions
        qui = grepl(":", x_labels_raw)
        if(any(qui)){
            x_inter = gsub(":.+", "", x_labels_raw[qui])
            tx_inter = table(x_inter)
            qui_auto = names(tx_inter)[tx_inter >= 2]

            group = list()

            for(i in seq_along(qui_auto)){
                var_left = qui_auto[i]
                qui_select = grepl(paste0(var_left, ":"), x_labels_raw, fixed = TRUE)

                # we require to be next to each other
                if(!all(diff(which(qui_select)) == 1)) next

                x_select = x_labels_raw[qui_select]

                if(all(grepl(":.+::", x_select))){
                    # This is a fixest call
                    first_part = strsplit(x_select[1], "::")[[1]][1]
                    var_right = gsub(".+:", "", first_part)
                    n_max = nchar(first_part) + 2
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

                v_name = dict_apply(c(var_left, var_right), dict)
                # group_name = paste0("&substitute(x %*% (y == ldots), list(x = \"", v_name[1], "\", y = \"", v_name[2], "\"))")
                group_name = replace_and_make_callable("__x__ %*% (__y__ == ldots)", list(x = v_name[1], y = v_name[2]), text_as_expr = TRUE)

                group[[group_name]] = escape_regex(paste0("%", var_left, ":", var_right))

                # We update the labels
                x_labels[qui_select] = substr(x_select, n_max + 1, nchar(x_select))
            }
        }
    }

    #
    # Plot (start) ####
    #


    all_plot_args = unique(c(names(par()), names(formals(plot.default))))
    pblm = setdiff(names(dots), all_plot_args)
    if(length(pblm) > 0){
        warning("The following argument", ifsingle(pblm, " is not a", "s are not"), " plotting argument", ifsingle(pblm, " and is", "s and are"), " therefore ignored: ", enumerate_items(pblm), ".")
        dots[pblm] = NULL
    }

    # preparation of the do.call
    dots$col = col

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

    listDefault(dots, "ylim", my_ylim)

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

        # The backgound
        if(!missing(bg)){
            dx = diff(my_xlim)
            dy = diff(my_ylim)
            rect(xleft = my_xlim[1] - dx, ybottom = my_ylim[1] - dy, xright = my_xlim[2] + dx, ytop = my_ylim[2] + dy, col = bg)
        }

        # The grid
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

        # Reference line

        if(identical(ref.line, "auto")){
            ref.line = suggest_ref_line && length(unique(prms[prms$is_ref, "estimate_names_raw"])) == 1
        }

        if(!isLogical(ref.line)){
            stop("Argument 'ref.line' must be either a logical, either equal to 'auto'. Currently it is none of these.")
        } else if(ref.line){

            ref_pblm = is.null(prms$is_ref) || !any(prms$is_ref)

            if(ref_pblm && !"v" %in% names(ref.line.par)){
                warning("You can use the argument 'ref.line' only when interactions are provided and a reference is found, or if you provided a reference with argument 'ref'. You can still draw vertical lines by using 'v' in argument 'ref.line.par'. Example: ref.line.par=list(v = ", round(x_value[floor(length(x_value)/2)]), ", col=2).")
            } else {

                if(!ref_pblm){
                    where = tapply(prms[prms$is_ref, "x"], prms[prms$is_ref, "estimate_names_raw"], mean)
                }

                listDefault(ref.line.par, "v", where)
                listDefault(ref.line.par, "lty", 2)
                do.call("abline", ref.line.par)
            }
        }

        box()
        axis(2)

        if(AXIS_AS_NUM){
            axis(1)
        } else {
            if(any(grepl("^&", x_labels))){
                # means we call expression()
                # drawback => expressions can overlap
                qui = grepl("^&", x_labels)
                if(any(!qui)){
                    axis(1, at = x_at[!qui], labels = x_labels[!qui])
                }

                for(i in which(qui)){
                    # my_expr = gsub("^&", "", x_labels[i])
                    # if(grepl("^(expression|substitute)\\(", my_expr)){
                    #     # direct evaluation
                    #     my_lab = eval(parse(text = my_expr))
                    # } else {
                    #     my_lab = eval(parse(text = paste0("expression(", my_expr, ")")))
                    # }
                    axis(1, at = x_at[i], labels = expr_builder(x_labels[i]))
                }

            } else {
                # easy case: only character
                axis(1, at = x_at, labels = x_labels)
            }
        }

    }

    n_id = length(unique(prms$id))

    #
    # pt.join ####
    #

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
        segments(x0=x, y0=ci_low, x1=x, y1=ci_high, lwd = par_fit(ci.lwd, prms$id), col = par_fit(ci.col, prms$id), lty = par_fit(ci.lty, prms$id))

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
        # Only if not a reference
        qui = ci_high != ci_low
        #  i) ci_high
        segments(x0=x[qui]-ci.width, y0=ci_high[qui], x1=x[qui]+ci.width, y1=ci_high[qui], lwd = par_fit(ci.lwd, prms$id[qui]), col = par_fit(ci.col, prms$id[qui]), lty = par_fit(ci.lty, prms$id[qui]))
        #  ii) ci_low
        segments(x0=x[qui]-ci.width, y0=ci_low[qui], x1=x[qui]+ci.width, y1=ci_low[qui], lwd = par_fit(ci.lwd, prms$id[qui]), col = par_fit(ci.col, prms$id[qui]), lty = par_fit(ci.lty, prms$id[qui]))
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

                my_join.par$x = prms$x[prms$id == i]

                my_join.par$y = prms$ci_high[prms$id == i]
                do.call("lines", my_join.par)

                my_join.par$y = prms$ci_low[prms$id == i]
                do.call("lines", my_join.par)
            }
        } else {
            ci.join.par$x = prms$x

            ci.join.par$y = prms$ci_high
            do.call("lines", ci.join.par)

            ci.join.par$y = prms$ci_low
            do.call("lines", ci.join.par)
        }
    }

    #
    # ci.fill ####
    #

    if(ci.fill){
        # We join the extremities of the conf. int.

        listDefault(ci.fill.par, "col", rgb(0.5, 0.5, 0.5, 0.5))
        listDefault(ci.fill.par, "alpha", 0.5)

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

                my_join.par$x = c(prms$x[prms$id == i], rev(prms$x[prms$id == i]))

                my_join.par$y = c(prms$ci_high[prms$id == i], rev(prms$ci_low[prms$id == i]))
                do.call("polygon", my_join.par)
            }
        } else {
            # colors => adding alpha if needed
            my_alpha = ci.fill.par$alpha[1]
            ci.fill.par$alpha = NULL
            my_col = as.vector(col2rgb(ci.fill.par$col[1], alpha = TRUE)) / 255
            if(my_col[4] == 1){
                ci.fill.par$col = rgb(my_col[1], my_col[2], my_col[3], my_alpha)
            }

            ci.fill.par$x = c(prms$x, rev(prms$x))
            ci.fill.par$y = c(prms$ci_high, rev(prms$ci_low))
            do.call("polygon", ci.fill.par)
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
            point.par = point.par[lengths(point.par) > 0]
            do.call("lines", point.par)
        }
    } else {
        dots$pch = par_fit(pt.pch, prms$id)
        dots$cex = par_fit(pt.cex, prms$id)
        dots$col = par_fit(pt.col, prms$id)
        do.call("lines", dots)
    }

    #
    # Group ####
    #

    if(!add && !missing(group) && !is.null(group) && !is.null(x_labels)){

        if(!is.list(group)) stop("Argument 'group' must be a list.")

        axis_prms_fit = c("lwd", "tcl", "line", "tick", "lty", "col", "lwd.ticks", "col.ticks")


        # axis
        listDefault(group.par, "lwd", 2)
        listDefault(group.par, "tcl", 0.75)
        listDefault(group.par, "line", 3)

        axis_par = group.par[!grepl("^lab", names(group.par))]
        axis_par$side = 1
        axis_par$labels = NA

        # label
        listDefault(group.par, "lab.line", 2.5)
        listDefault(group.par, "lab.cex", 1)
        listDefault(group.par, "lab.font", 1)
        listDefault(group.par, "lab.col", 1)

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

            info = do.call("axis", my_axis)

            if(!is.null(group_name)){
                if(grepl("^&", group_name)){
                    group_name = eval(parse(text = gsub("^&", "", group_name)))
                }

                axis(1, mean(info), labels = group_name, tick = FALSE, line = par_fit(group.par$lab.line, i), cex.axis = par_fit(group.par$lab.cex, i), font = par_fit(group.par$lab.font, i), col.axis = par_fit(group.par$lab.col, i))
            }
        }
    }

    res = list(prms=prms, is_interaction = is_interaction, at = x_at, labels = x_labels)
    return(invisible(res))
}




coefplot_prms = function(object, ..., sd, ci_low, ci_high, x, x.shift = 0, dict, drop, order, ci_level = 0.95, ref = "&auto", only.inter = TRUE, sep, as.multiple = FALSE){

    # get the default for:
    # dict, ci.level, ref, only.inter


    dots = list(...)
    is_internal = dots$internal__
    if(is.null(is_internal)) is_internal = FALSE

    varlist = list(ci = paste0(ci_level * 100, "%"))
    dots_drop = c()

    suggest_ref_line = FALSE
    multiple_est = FALSE
    NO_NAMES = FALSE
    IS_INTER = FALSE
    AXIS_AS_NUM = FALSE
    if(is_internal == FALSE && is.list(object) && class(object)[1] == "list"){
        # This is a list of estimations

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

                # dealing with interactions
                if(prms$is_interaction){
                    if(i > 1 && any(!all_inter)){
                        rerun = TRUE
                        mc$only.inter = FALSE
                        break
                    } else {
                        all_inter[i] = TRUE
                        my_root = names(prms$is_interaction)
                        if(i == 1 || all(all_inter_root == my_root)){
                            all_inter_root[i] = my_root
                        } else {
                            rerun = TRUE
                            mc$only.inter = FALSE
                            break
                        }
                    }
                } else {
                    all_inter[i] = FALSE
                }

                # Some meta variables
                varlist$ci = unique(prms$varlist$ci)
                varlist$depvar = unique(prms$varlist$depvar)
                varlist$var = unique(prms$varlist$var)
                varlist$fe = unique(prms$varlist$fe)

                dots_drop = unique(c(dots_drop, prms$dots_drop))
                suggest_ref_line = suggest_ref_line && prms$suggest_ref_line

                num_axis = num_axis && prms$num_axis

                df_prms = prms$prms
                df_prms$est_nb = i
                res[[i]] = df_prms
            }
        }

        IS_INTER = all(all_inter)
        AXIS_AS_NUM = num_axis

        all_estimates = do.call("rbind", res)

        # We respect the order provided by the user
        my_names_order = unique(all_estimates$estimate_names)
        my_order = 1:length(my_names_order)
        names(my_order) = my_names_order
        all_estimates$id_order = my_order[all_estimates$estimate_names]
        all_estimates = all_estimates[base::order(all_estimates$id_order, all_estimates$est_nb), ]

        # we rescale
        if(missing(sep)){
            all_sep = c(0.2, 0.2, 0.18, 0.16)
            if(length(all_sep) < nb_est - 1){
                sep = 1 / (n-1) * 0.7
            }
            sep = all_sep[nb_est - 1]
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

        all_estimates$x_new = all_estimates$id_order + ((all_estimates$est_nb - 1) / (nb_est - 1) - 0.5) * ((nb_est - 1) * sep)

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

        n <- length(estimate)

        if(missing(sd)){
            if(missing(ci_low) || missing(ci_high)) stop("If 'sd' is not provided, you must provide the arguments 'ci_low' and 'ci_high'.")

            varlist$ci = NULL
        } else {
            if(!missing(ci_low) || !missing(ci_high)) warning("Since 'sd' is provided, arguments 'ci_low' or 'ci_high' are ignored.")

            # We compue the CI
            nb = abs(qnorm((1-ci_level)/2))
            ci_high = estimate + nb*sd
            ci_low = estimate - nb*sd
        }

        #
        # Interactions ####
        #

        ref_id = NA
        if(only.inter && !is.null(names(estimate))){
            all_vars = names(estimate)
            if(any(grepl("::", all_vars))){

                IS_INTER = TRUE

                is_info = FALSE
                if("fixest" %in% class(object)){
                    is_info = TRUE
                    interaction.info = object$interaction.info
                    is_ref = interaction.info$is_ref
                    items = interaction.info$items
                    is_num = is.numeric(interaction.info$fe_type)
                }

                # We retrict only to interactions
                root_interaction = all_vars[grepl("::", all_vars)]
                # We keep only the first one !
                root_interaction = unique(gsub("::.+", "", root_interaction))[1]

                names(IS_INTER) = root_interaction

                inter_keep = grepl(root_interaction, all_vars, fixed = TRUE)
                my_inter = estimate_names_raw = all_vars[inter_keep]
                estimate = estimate[inter_keep]
                ci_high = ci_high[inter_keep]
                ci_low = ci_low[inter_keep]

                if(!is_info){
                    is_ref = rep(FALSE, length(my_inter))
                }

                # We extract the name of the variables
                fe_name = gsub(".+:", "", root_interaction)
                if(fe_name %in% names(dict)) fe_name = dict[fe_name]
                varlist$fe = fe_name

                var_name = gsub(":[[:alnum:]\\._]+", "", root_interaction)
                if(var_name %in% names(dict)) var_name = dict[var_name]
                varlist$var = var_name

                # We construct the x-axis
                inter_values = estimate_names = gsub(".+::", "", my_inter)
                names(estimate) = inter_values

                inter_values_num = tryCatch(as.numeric(inter_values), warning = function(x) x)
                if(is_info){

                    if(identical(ref, "&auto")){
                        # We add the reference used in the estimation

                        if(length(inter_values) != sum(!is_ref)){
                            stop("Internal error regarding the lengths of vectors of coefficients.")
                        }

                        if(any(is_ref)){
                            ref_id = which(is_ref)
                        }

                        my_values = my_ci_low = my_ci_high = rep(NA, length(is_ref))
                        names(my_values) = names(my_ci_low) = names(my_ci_high) = items

                        my_values[inter_values] = estimate
                        my_ci_high[inter_values] = ci_high
                        my_ci_low[inter_values] = ci_low

                        qui = which(is.na(my_values))
                        my_values[qui] = 0
                        my_ci_high[qui] = my_values[qui]
                        my_ci_low[qui] = my_values[qui]

                        estimate = my_values
                        ci_high = my_ci_high
                        ci_low = my_ci_low
                        estimate_names = items
                        estimate_names_raw = paste0(root_interaction, "::", items)

                        # We suggest a reference
                        suggest_ref_line = any(interaction.info$fe_type %in% c("numeric", "integer", "factor"))
                    } else {
                        is_ref = rep(FALSE, length(estimate))
                    }

                    if(is_num && missing(x)){
                        AXIS_AS_NUM = TRUE
                        names(estimate) = NULL
                        x = items
                    }

                } else if(is.numeric(inter_values_num)) {

                    # We check these are "real" numbers and not just "codes"
                    all_steps = diff(sort(inter_values_num))
                    ts = table(all_steps)
                    step_mode = as.numeric(names(ts)[which.max(ts)])
                    all_steps_rescaled = all_steps / step_mode

                    if(any(all_steps_rescaled < 10) && missing(x)){
                        AXIS_AS_NUM = TRUE
                        x = inter_values_num
                    }
                }

                n = length(estimate)
            }
        }

        # The DF of all the parameters
        prms = data.frame(estimate = estimate, ci_low = ci_low, ci_high = ci_high)
        if(IS_INTER){
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
        if(!identical(ref, "&auto") && length(ref) > 0 && !isFALSE(ref)){

            if(is.null(names(ref))){
                if(!is.character(ref) || length(ref) > 1){
                    check_arg(ref, "singleCharacter", "Argument 'ref' must be either: a single character, either a list or a named integer vector. REASON")
                } else {
                    refname = ref
                    ref = list()
                    ref[[refname]] = 1
                }
            }

            ref = unlist(ref)
            if(!isScalar(ref, int = TRUE)){
                reason = ifelse(length(ref) == 1, " an integer", " of length 1")
                stop("Argument 'ref' must be either: a single character, either a list or a named integer vector. The integer gives the position of the reference among the coefficients. Currently this is not ", reason, ".")
            }

            # we recreate the parameters
            n = nrow(prms)
            prms$is_ref = FALSE
            ref_row = data.frame(estimate = 0, ci_low = 0, ci_high = 0, estimate_names = names(ref), estimate_names_raw = names(ref), is_ref = TRUE)
            prms = rbind(prms, ref_row)
            if(ref > n) ref = n + 1
            ids = 1:n
            ids[ids >= ref] = ids[ids >= ref] + 1
            prms = prms[base::order(c(ids, ref)), ]

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

    return(list(prms=prms, is_interaction = IS_INTER, num_axis = AXIS_AS_NUM, at = x_at, labels = x_labels, x_labels_raw = x_labels_raw, varlist=varlist, dots_drop=dots_drop, xlim = my_xlim, suggest_ref_line=suggest_ref_line, multiple_est=multiple_est))
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

            warning(info, enumerate_items(paste0("__", setdiff(my_variables, names(varlist)), "__"), "is"), " not valid, thus ignored.", call. = FALSE)

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
                        res = deparse(eval(parse(text = res)))
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
            res = eval(parse(text = my_expr), parent.frame(2))
        } else {
            res = eval(parse(text = paste0("expression(", my_expr, ")")))
        }
    } else {
        res = x
    }

    res
}





































