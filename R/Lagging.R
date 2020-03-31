#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sun Nov 24 09:18:19 2019
# ~: Lagging
#----------------------------------------------#

# -------------------------------------------------------------- #
# Contains all functions to ensure the easy lagging of variables
# within fixest estimations.
#
# -------------------------------------------------------------- #

# Make sure data.table knows we know we're using it
.datatable.aware = TRUE


panel_setup = function(data, panel.id, time.step = "unitary", duplicate.method = c("none", "first"), DATA_MISSING = FALSE, from_fixest = FALSE){
    # Function to setup the panel.
    # Used in lag.formula, panel, and fixest_env (with argument panel.id and panel.args)
    # DATA_MISSING: arg used in lag.formula


    # for error handling => refers to the right function
    my_call = deparse(sys.calls()[[sys.nframe()-1]])[1] # call can have svl lines
    nmax = 40
    if(nchar(my_call) > nmax) my_call = paste0(substr(my_call, 1, nmax-1), "...")
    my_call = paste0(my_call, ":\n ")


    duplicate.method = try(match.arg(duplicate.method, choices = c("none", "first")), silent = TRUE)
    if(class(duplicate.method) == "try-error"){
        stop(my_call, "The argument 'duplicate.method' must be one of: 'none' or 'first'.", call. = FALSE)
    }

    if("formula" %in% class(panel.id)){
        tm = terms_hat(panel.id)
        var_id_time = attr(tm, "term.labels")
        if(length(var_id_time) != 2){
            stop(my_call, "The formula of the argument 'panel.id' must contain exactly two variables in the right hand side (currently there ", ifsingle(var_id_time, "is ", "are "), length(var_id_time), ").", call. = FALSE)
        }

        all_vars = all.vars(tm)
    } else if(is.character(panel.id) && length(panel.id) == 2){
        all_vars = var_id_time = panel.id
    } else if(is.character(panel.id) && length(panel.id) == 1){
        var_id_time = gsub("(^ +| +$)", "", strsplit(panel.id, ",")[[1]])
        all_vars = var_id_time = var_id_time[nchar(var_id_time) > 0]
        if(length(var_id_time) != 2){
            stop(my_call, "The argument 'panel.id' must be either: i) a one sided formula (e.g. ~id+time), ii) a character vector of length 2 (e.g. c('id', 'time'), or iii) a character scalar of two variables separated by a comma (e.g. 'id,time'). Currently it is neither of the three." , call. = FALSE)
        }
    } else {
        stop(my_call, "The argument 'panel.id' must be either: i) a one sided formula (e.g. ~id+time), ii) a character vector of length 2 (e.g. c('id', 'time'), or iii) a character scalar of two variables separated by a comma (e.g. 'id,time'). Currently it is neither of the three." , call. = FALSE)
    }


    if(DATA_MISSING){
        if(!all(all_vars %in% ls(parent.frame(2)))){
            pblm = setdiff(all_vars, ls(parent.frame(2)))
            stop(my_call, "In the argument 'panel.id', the variable", enumerate_items(pblm, "s.is.past.quote"), " not found in the environment.", call. = FALSE)
        }
        id = eval(parse(text = var_id_time[1]), parent.frame(2))
        time = eval(parse(text = var_id_time[2]), parent.frame(2))
    } else {
        if(!all(all_vars %in% names(data))){
            pblm = setdiff(all_vars, names(data))
            stop(my_call, "In the argument 'panel.id', the variable", enumerate_items(pblm, "s.is.quote"), " not in the data set.", call. = FALSE)
        }
        id = eval(parse(text = var_id_time[1]), data)
        time = eval(parse(text = var_id_time[2]), data)
    }

    panel.id = as.formula(paste0("~", var_id_time[1], "+", var_id_time[2]))

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
        stop(my_call, "The argument 'time.step' must be equal to 'unitary', 'consecutive' or to a number.", call. = FALSE)
    } else if(is.character(time.step)){
        ts = try(match.arg(time.step, c("unitary", "consecutive", "within.consecutive")), silent = TRUE)
        if("try-error" %in% class(ts)){
            stop(my_call, "The argument 'time.step' must be one of 'unitary', 'consecutive' or 'within.consecutive', or a number representing the time step.", call. = FALSE)
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
                if(from_fixest){
                    stop(my_call, "The time variable must be numeric or at least convertible to numeric. So far the conversion has failed (time variable's class is currently ", enumerate_items(class(time)), "). Alternatively, you can have more options to set up the panel using the function panel().", call. = FALSE)
                } else {
                    stop(my_call, "To use the 'unitary' time.step, the time variable must be numeric or at least convertible to numeric. So far the conversion has failed (time variable's class is currently ", enumerate_items(class(time)), ").", call. = FALSE)
                }
            }

            time = time_new
        }

        if(any(time %% 1 != 0)){
            if(from_fixest){
                stop(my_call, "The time variable", ifelse(time_conversion, " (which has been converted to numeric)", ""), " must be made of integers. So far this is not the case. Alternatively, you can have more options to set up the panel using the function panel().", call. = FALSE)
            } else {
                stop(my_call, "To use the 'unitary' time.step, the time variable", ifelse(time_conversion, " (which has been converted to numeric)", ""), " must be made of integers. So far this is not the case. Alternatively, you can give a number in time.step.", call. = FALSE)
            }
        }

    } else if(!is.character(time.step)){
        if(!is.numeric(time)){
            stop(my_call, "If 'time.step' is a number, then the time variable must also be a number (this is not the case: its class is currently ", enumerate_items(class(time)), ").", call. = FALSE)
        }
    }

    # Computation quf
    id = quickUnclassFactor(id)
    time_full = quf_sorted(time, addItem = TRUE)

    #
    # WIP: define this unitary time step!!!! not straightforward at all!
    # what to do with the ones that don't fit the unit???
    # example: time = c(1, 3, 6, 11) => smallest unit is 2, but it does not divide the others

    # NOTA: "time" is not used as a CPP index: so we don't care that it goes from 1 to K
    #       or from 0 to K-1. We only care that the numbers are consecutive

    # Releveling the time ID depending on the time.step
    if(time.step %in% c("consecutive", "within.consecutive")){
        # for within.consecutive, we deal with it after sorting
        time = time_full$x
    } else if(time.step == "unitary"){
        time_unik = time_full$items
        all_steps = unique(diff(time_unik))
        my_step = cpp_pgcd(unique(all_steps))

        # we rescale time_unik
        time_unik_new = (time_unik - min(time_unik)) / my_step
        time = time_unik_new[time_full$x]

        if(my_step != 1){
            message("NOTE: unitary time step taken: ", my_step, ".")
        }


    } else {
        time_unik = time_full$items

        # consistency check
        all_steps = diff(time_unik)
        if(any(all_steps %% time.step != 0)){
            obs_pblm = which(all_steps %% time.step != 0)

            stop(my_call, "If 'time.step' is a number, then it must be an exact divisor of all the difference between two consecutive time periods. This is currently not the case: ", time.step, " is not a divisor of ", all_steps[obs_pblm][1], " (the difference btw the time periods ", time_unik[obs_pblm[1] + 1], " and ", time_unik[obs_pblm[1]], ").", call. = FALSE)
        }

        # we rescale time_unik // checks done beforehand
        time_unik_new = (time_unik - min(time_unik)) / time.step
        time = time_unik_new[time_full$x]
    }

    # Here time is always integer: we convert it if necessary (hasten 'order' calls)
    time = as.integer(time)

    order_it = order(id, time)
    order_inv = order(order_it)

    id_sorted = id[order_it]
    time_sorted = time[order_it]

    if(time.step == "within.consecutive"){
        time_sorted = 1:length(time_sorted)
    }

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

            stop(my_call, "The panel identifiers contain duplicate values: this is not allowed since lag/leads are not defined for them. For example (id, time) = (", id_dup, ", ", time_dup, ") appears ", n_times(dup_info$n_dup), ". Please provide data without duplicates -- or you can also use duplicate.method = 'first' (see Details).", call. = FALSE)
        }
    }

    res = list(order_it = order_it, order_inv = order_inv, id_sorted = id_sorted, time_sorted = time_sorted, na_flag = na_flag, panel.id = panel.id)
    if(na_flag) res$is_na = is_na
    res
}


#' @describeIn l Forwards a variable (inverse of lagging) in a \code{fixest} estimation
f = function(x, lead = 1, fill = NA){
    l(x, -lead, fill)
}


#' Lags a variable in a \code{fixest} estimation
#'
#' Produce lags or leads in the formulas of \code{fixest} estimations or when creating variables in a \code{\link[data.table]{data.table}}. The data must be set as a panel beforehand (either with the function \code{\link[fixest]{panel}} or with the argument \code{panel.id} in the estimation).
#'
#' @param x The variable.
#' @param lag A vector of integers giving the number of lags. Negative values lead to leads. This argument can be a vector when using it in fixest estimations. When creating variables in a \code{\link[data.table]{data.table}}, it **must** be of length one.
#' @param lead A vector of integers giving the number of leads. Negative values lead to lags. This argument can be a vector when using it in fixest estimations. When creating variables in a \code{\link[data.table]{data.table}}, it **must** be of length one.
#' @param fill A scalar, default is \code{NA}. How to fill the missing values due to the lag/lead? Note that in a \code{fixest} estimation, 'fill' must be numeric (not required when creating new variables).
#'
#' @return
#' These functions can only be used i) in a formula of a \code{fixest} estimation, or ii) when creating variables within a \code{fixest_panel} object (obtained with function \code{\link[fixest]{panel}}) which is alaos a \code{\link[data.table]{data.table}}.
#'
#' @seealso
#' The function \code{\link[fixest]{panel}} changes \code{data.frames} into a panel from which the functions \code{l} and \code{f} can be called. Otherwise you can set the panel 'live' during the estimation using the argument \code{panel.id} (see for example in the function \code{\link[fixest]{feols}}).
#'
#' @examples
#'
#' data(base_did)
#'
#' # Setting a data set as a panel...
#' pdat = panel(base_did, ~id+period)
#'
#' # ...then using the functions l and f
#' est1 = feols(y~l(x1, 0:1), pdat)
#' est2 = feols(f(y)~l(x1, -1:1), pdat)
#' est3 = feols(l(y)~l(x1, 0:3), pdat)
#' etable(est1, est2, est3, order = c("f", "^x"), drop="Int")
#'
#' # or using the argument panel.id
#' feols(f(y)~l(x1, -1:1), base_did, panel.id = ~id+period)
#'
#' # l() and f() can also be used within a data.table:
#' if(require("data.table")){
#'   pdat_dt = panel(as.data.table(base_did), ~id+period)
#'   # Now since pdat_dt is also a data.table
#'   #   you can create lags/leads directly
#'   pdat_dt[, x1_l1 := l(x1)]
#'   pdat_dt[, c("x1_l1_fill0", "y_f2") := .(l(x1, fill = 0), f(y, 2))]
#' }
#'
#'
#'
l = function(x, lag = 1, fill = NA){

    value = x

    sys_calls = rev(sys.calls())

    # To improve => you don't wanna check all frames, only the relevant ones
    from_fixest = FALSE
    for(where in 1:min(6, sys.nframe())){
        if(exists("panel__meta__info", parent.frame(where))){
            from_fixest = TRUE
            meta_info = get("panel__meta__info", parent.frame(where))
            break
        }
    }

    if(from_fixest == FALSE){
        # Using l/f within data.table

        fl_authorized = getOption("fixest_fl_authorized")

        if(fl_authorized){
            # Further control
            check_arg(lag, "singleInteger", "When creating lags (or leads) within a data.table with l() (or f()), the argument 'lag' must be a single integer. REASON")

            if(length(fill) != 1){
                stop("The length of argument 'fill' must be exaclty 1. Its current length is ", length(fill), ".")
            } else if(!is.na(fill)){
                if(is.list(fill)){
                    stop("Argument fill must be a 'scalar', currenlty it's a list!")
                }

                if(is.numeric(value) && !is.numeric(fill)){
                    mc = match.call()
                    stop("Argument 'fill' must be of the same type as ", deparse_long(mc$x), ", which is numeric.")
                }

                if(!is.numeric(value) && is.numeric(fill)){
                    mc = match.call()
                    stop("Argument 'fill' must be of the same type as ", deparse_long(mc$x), ", which is not numeric while 'fill' is.")
                }
                # I could go further in checking but it's enough
            }

            # we fetch the panel information
            qui  = which.max(grepl("^.\\[\\.data\\.table", sys_calls))

            if(length(qui) == 0){
                options(fixest_fl_authorized = FALSE)
                stop("Unknown error when trying to create a lag (or lead) with function l() (or f()) within data.table. This is a 'fixest' error, may be worth reporting if needed.")
            }

            var = gsub("^[^\\(]+\\(|,.+", "", sys_calls[qui])
            m = eval(parse(text = var), envir = sys.frames()[qui])

            if(!"fixest_panel" %in% class(m)){
                stop("You can use l() or f() only when the data set is of class 'fixest_panel', you can use function panel() to set it.")
            }

            meta_info = attr(m, "panel_info")
        } else {
            stop("Function l() (or f()) is only callable within 'fixest' estimations or within a variable creation with data.table (i.e. using ':=') where the data set is a 'fixest_panel' (obtained from panel()). Alternatively, you can use lag.formula().")
        }
    }

    # we get the observation id!

    if(length(meta_info$id_sorted) != length(meta_info$time_sorted)){
        stop("Internal error: lengths of the individuals and the time vectors are different.")
    }

    obs_lagged = cpp_lag_obs(id = meta_info$id_sorted, time = meta_info$time_sorted, nlag = lag)

    # the lagged value => beware of NAs in IDs!
    if(meta_info$na_flag){
        value_sorted = (value[!meta_info$is_na])[meta_info$order_it]
    } else {
        value_sorted = value[meta_info$order_it]
    }

    value_lagged = value_sorted[obs_lagged]

    if(!is.na(fill)){
        qui_na = is.na(obs_lagged)
        value_lagged[qui_na] = fill
    }

    if(meta_info$na_flag == FALSE){
        res = value_lagged[meta_info$order_inv]
    } else{
        res = rep(NA, length(value))
        res[!meta_info$is_na] = value_lagged[meta_info$order_inv]
    }

    res
}


l_expand = function(x, k=1, fill){

    mc = match.call()

    if(missing(x)) stop("Argument 'x' cannot be missing.")
    check_arg(k, "integerVector")
    check_arg(fill, "numericScalarNaok")

    mc_new = mc

    res = c()
    for(i in 1:length(k)){

        if(k[i] == 0){
            res[i] = deparse_long(mc_new$x)
        } else {
            if(k[i] < 0){
                mc_new[[1]] = as.name("f")
                mc_new$k = as.numeric(-k[i])
            } else {
                mc_new[[1]] = as.name("l")
                mc_new$k = as.numeric(k[i])
            }

            res[i] = gsub("x = |k = ", "", deparse_long(mc_new))
        }
    }

    res
}

# same as l_expand, but opposite sign
f_expand = function(x, k=1, fill){

    mc = match.call()

    if(missing(x)) stop("Argument 'x' cannot be missing.")
    check_arg(k, "integerVector")
    check_arg(fill, "numericScalarNaok")

    mc_new = mc

    res = c()
    for(i in 1:length(k)){

        if(k[i] == 0){
            res[i] = deparse_long(mc_new$x)
        } else {
            if(k[i] < 0){
                mc_new[[1]] = as.name("l")
                mc_new$k = as.numeric(-k[i])
            } else {
                mc_new[[1]] = as.name("f")
                mc_new$k = as.numeric(k[i])
            }

            res[i] = gsub("x = |k = ", "", deparse_long(mc_new))
        }
    }


    res
}


#' Lags a variable using a formula
#'
#' Lags a variable using panel id + time identifiers in a formula.
#'
#' @inheritParams panel
#'
#'
#' @param x A formula of the type \code{var ~ id + time} where \code{var} is the variable to be lagged, \code{id} is a variable representing the panel id, and \code{time} is the time variable of the panel.
#' @param k An integer giving the number of lags. Default is 1. For leads, just use a negative number.
#' @param data Optional, the data.frame in which to evaluate the formula. If not provided, variables will be fetched in the current environment.
#' @param fill Scalar. How to fill the observations without defined lead/lag values. Default is \code{NA}.
#' @param ... Not currently used.
#'
#' @return
#' It returns a vector of the same type and length as the variable to be lagged in the formula.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' Alternatively, the function \code{\link[fixest]{panel}} changes a \code{data.frame} into a panel from which the functions \code{l} and \code{f} (creating leads and lags) can be called. Otherwise you can set the panel 'live' during the estimation using the argument \code{panel.id} (see for example in the function \code{\link[fixest]{feols}}).
#'
#' @examples
#' # simple example with an unbalanced panel
#' base = data.frame(id = rep(1:2, each = 4),
#'                   time = c(1, 2, 3, 4, 1, 4, 6, 9), x = 1:8)
#'
#' base$lag1 = lag(x~id+time,  1, base) # lag 1
#' base$lead1 = lag(x~id+time, -1, base) # lead 1
#' base$lag2_fill0 = lag(x~id+time, 2, base, fill = 0)
#' # with time.step = "consecutive"
#' base$lag1_consecutive = lag(x~id+time, 1, base, time.step = "consecutive")
#' #   => works for indiv. 2 because 9 (resp. 6) is consecutive to 6 (resp. 4)
#' base$lag1_within.consecutive = lag(x~id+time, 1, base, time.step = "within")
#' #   => now two consecutive years within each indiv is one lag
#'
#' print(base)
#'
#' # Argument time.step = "consecutive" is
#' # mostly useful when the time variable is not a number:
#' # e.g. c("1991q1", "1991q2", "1991q3") etc
#'
#' # with duplicates
#' base_dup = data.frame(id = rep(1:2, each = 4),
#'                       time = c(1, 1, 1, 2, 1, 2, 2, 3), x = 1:8)
#'
#' # Error because of duplicate values for (id, time)
#' try(lag(x~id+time, 1, base_dup))
#'
#'
#' # Error is bypassed, lag corresponds to first occurence of (id, time)
#' lag(x~id+time, 1, base_dup, duplicate.method = "first")
#'
#'
#' # Playing with time steps
#' base = data.frame(id = rep(1:2, each = 4),
#'                   time = c(1, 2, 3, 4, 1, 4, 6, 9), x = 1:8)
#'
#' # time step: 0.5 (here equivalent to lag of 1)
#' lag(x~id+time, 2, base, time.step = 0.5)
#'
#' # Error: wrong time step
#' try(lag(x~id+time, 2, base, time.step = 7))
#'
#' # Adding NAs + unsorted IDs
#' base = data.frame(id = rep(1:2, each = 4),
#'                   time = c(4, NA, 3, 1, 2, NA, 1, 3), x = 1:8)
#'
#' base$lag1 = lag(x~id+time, 1, base)
#' base$lag1_within = lag(x~id+time, 1, base, time.step = "w")
#' base_bis = base[order(base$id, base$time),]
#'
#' print(base_bis)
#'
#' # You can create variables without specifying the data within data.table:
#' if(require("data.table")){
#'   base = data.table(id = rep(1:2, each = 3), year = 1990 + rep(1:3, 2), x = 1:6)
#'   base[, x.l1 := lag(x~id+year, 1)]
#' }
#'
#'
#'
lag.formula = function(x, k = 1, data, time.step = "unitary", fill = NA, duplicate.method = c("none", "first"), ...){
    # Arguments:
    # time.step: default: "consecutive", other option: "unitary" (where you find the most common step and use it -- if the data is numeric), other option: a number, of course the time must be numeric

    # LATER:
    # - add duplicate.method = "sum" // although it is not super useful in my opinion (but maybe other users disagree)
    # - several lags => matrix? I donno...

    # IMPORTANT TO BE DONE:
    # - check argument 'fill' in here (not in panel_setup)

    # Checking arguments in ...
    any_invalid = check_dots_args(match.call(expand.dots = FALSE),
                                  suggest_args = c("k", "data", "time.step", "fill"))
    if(any_invalid){
        warning(attr(any_invalid, "msg"))
    }

    # Controls
    duplicate.method = match.arg(duplicate.method)

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
        stop("In the formula the variable", enumerate_items(qui_pblm, "s.is.quote"), " not in the ", ifelse(DATA_MISSING, "environment. Add argument data?", "data."))
    }

    # checking k
    check_arg(k, "integerScalar", "Argument 'k' must be a single integer.")

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

    meta_info = panel_setup(data, panel.id = x, time.step = time.step, duplicate.method = duplicate.method, DATA_MISSING = DATA_MISSING)

    # we get the observation id!
    obs_lagged = cpp_lag_obs(id = meta_info$id_sorted, time = meta_info$time_sorted, nlag = k)

    # the lagged value => beware of NAs in IDs!
    if(meta_info$na_flag){
        value_sorted = (value[!meta_info$is_na])[meta_info$order_it]
    } else {
        value_sorted = value[meta_info$order_it]
    }

    value_lagged = value_sorted[obs_lagged]

    if(!is.na(fill)){
        qui_na = is.na(obs_lagged)
        value_lagged[qui_na] = fill
    }

    if(meta_info$na_flag == FALSE){
        res = value_lagged[meta_info$order_inv]
    } else{
        res = rep(NA, length(value))
        res[!meta_info$is_na] = value_lagged[meta_info$order_inv]
    }

    res
}



#' Constructs a \code{fixest} panel data base
#'
#' Constructs a \code{fixest} panel data base out of a data.frame which allows to use leads and lags in \code{fixest} estimations and to create new variables from leads and lags if the data.frame was also a \code{\link[data.table]{data.table}}.
#'
#' @param data A data.frame.
#' @param panel.id The panel identifiers. Can either be: i) a one sided formula (e.g. \code{panel.id~id+time}), ii) a character vector of length 2 (e.g. \code{panel.id=c('id', 'time')}, or iii) a character scalar of two variables separated by a comma (e.g. \code{panel.id='id,time'}). Note that you can combine variables with \code{^} only inside formulas (see the dedicated section in \code{\link[fixest]{feols}}).
#' @param time.step The method to compute the lags. Can be equal to: \code{"unitary"} (default), \code{"consecutive"}, \code{"within.consecutive"}, or to a number. If \code{"unitary"}, then the largest common divisor between consecutive time periods is used (typically if the time variable represents years, it will be 1). This method can apply only to integer (or convertible to integer) variables. If \code{"consecutive"}, then the time variable can be of any type: two successive time periods represent a lag of 1. If \code{"witihn.consecutive"} then **within a given id**, two successive time periods represent a lag of 1. Finally, if the time variable is numeric, you can provide your own numeric time step.
#' @param duplicate.method If several observations have the same id and time values, then the notion of lag is not defined for them. If \code{duplicate.method = "none"} (default) and duplicate values are found, this leads to an error. You can use \code{duplicate.method = "first"} so that the first occurrence of identical id/time observations will be used as lag.
#'
#' @details
#' This function allows you to use leads and lags in a \code{fixest} estimation without having to provide the argument \code{panel.id}. It also offers more options on how to set the panel (with the additional arguments 'time.step' and 'duplicate.method').
#'
#' When the initial data set was also a \code{data.table}, not all operations are supported and some may dissolve the \code{fixest_panel}. This is the case when creating subselections of the initial data with additional attributes (e.g. pdt[x>0, .(x, y, z)] would dissolve the \code{fixest_panel}, meaning only a data.table would be the result of the call).
#'
#' If the initial data set was also a \code{data.table}, then you can create new variables from lags and leads using the functions \code{\link[fixest]{l}}() and \code{\link[fixest]{f}}(). See the example.
#'
#'
#' @return
#' It returns a data base identical to the one given in input, but with an additional attribute: \dQuote{panel_info}. This attribute contains vectors used to efficiently create lags/leads of the data. When the data is subselected, some bookeeping is performed on the attribute \dQuote{panel_info}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' The estimation methods \code{\link[fixest]{feols}}, \code{\link[fixest]{fepois}} and \code{\link[fixest]{feglm}}.
#'
#' The functions \code{\link[fixest]{l}} and \code{\link[fixest]{f}} to create lags and leads within \code{}fixest_panel objects.
#'
#' @examples
#'
#' data(base_did)
#'
#' # Setting a data set as a panel...
#' pdat = panel(base_did, ~id+period)
#'
#' # ...then using the functions l and f
#' est1 = feols(y~l(x1, 0:1), pdat)
#' est2 = feols(f(y)~l(x1, -1:1), pdat)
#' est3 = feols(l(y)~l(x1, 0:3), pdat)
#' etable(est1, est2, est3, order = c("f", "^x"), drop="Int")
#'
#' # or using the argument panel.id
#' feols(f(y)~l(x1, -1:1), base_did, panel.id = ~id+period)
#'
#' # You can use panel.id in various ways:
#' pdat = panel(base_did, ~id+period)
#' # is identical to:
#' pdat = panel(base_did, c("id", "period"))
#' # and also to:
#' pdat = panel(base_did, "id,period")
#'
#' # l() and f() can also be used within a data.table:
#' if(require("data.table")){
#'   pdat_dt = panel(as.data.table(base_did), ~id+period)
#'   # Now since pdat_dt is also a data.table
#'   #   you can create lags/leads directly
#'   pdat_dt[, x1_l1 := l(x1)]
#'   pdat_dt[, c("x1_l1_fill0", "y_f2") := .(l(x1, fill = 0), f(y, 2))]
#' }
#'
#'
panel = function(data, panel.id, time.step = "unitary", duplicate.method = c("none", "first")){

    if(missing(data)){
        stop("You must provide the argument 'data'.")
    } else if(!"data.frame" %in% class(data)){
        stop("Argument 'data' must be a data.frame.")
    }

    mc = match.call()

    meta_info = panel_setup(data, panel.id = panel.id, time.step = time.step, duplicate.method = duplicate.method)
    meta_info$call = mc

    # R makes a shallow copy of data => need to do it differently with DT

    if("data.table" %in% class(data) && requireNamespace("data.table", quietly=TRUE)){
        res = data.table::copy(data)
        data.table::setattr(res, "panel_info", meta_info)
        data.table::setattr(res, "class", c("fixest_panel", "data.table", "data.frame"))
    } else {
        res = data
        attr(res, "panel_info") = meta_info
        class(res) = unique(c("fixest_panel", class(res)))
    }

    res
}


#' Dissolves a \code{fixest} panel
#'
#' Transforms a \code{fixest_panel} object into a regular data.frame.
#'
#' @param x A \code{fixest_panel} object (obtained from function \code{\link[fixest]{panel}}).
#'
#' @return
#' Returns a data set of the exact same dimension. Only the attribute 'panel_info' is erased.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' Alternatively, the function \code{\link[fixest]{panel}} changes a \code{data.frame} into a panel from which the functions \code{l} and \code{f} (creating leads and lags) can be called. Otherwise you can set the panel 'live' during the estimation using the argument \code{panel.id} (see for example in the function \code{\link[fixest]{feols}}).
#'
#' @examples
#'
#' data(base_did)
#'
#' # Setting a data set as a panel
#' pdat = panel(base_did, ~id+period)
#'
#' # ... allows you to use leads and lags in estimations
#' feols(y~l(x1, 0:1), pdat)
#'
#' # Now unpanel => returns the initial data set
#' class(pdat) ; dim(pdat)
#' new_base = unpanel(pdat)
#' class(new_base) ; dim(new_base)
#'
#'
#'
unpanel = function(x){

    if("data.table" %in% class(x) && requireNamespace("data.table", quietly=TRUE)){
        data.table::setattr(x, "panel_info", NULL)
        data.table::setattr(x, "class", setdiff(class(x), "fixest_panel"))

        return(invisible(x))
    } else {
        attr(x, "panel_info") = NULL
        class(x) = setdiff(class(x), "fixest_panel")
    }

    x
}


#' Method to subselect from a \code{fixest_panel}
#'
#' Subselection from a \code{fixest_panel} which has been created with the function \code{\link[fixest]{panel}}. Also allows to create lag/lead variables with functions \code{\link[fixest]{l}}()/\code{\link[fixest]{f}}() if the \code{fixest_panel} is also a \code{\link[data.table]{data.table}}.
#'
#' @param x A \code{fixest_panel} object, created with the function \code{\link[fixest]{panel}}.
#' @param i Row subselection. Allows \code{\link[data.table]{data.table}} style selection (provided the data is also a data.table).
#' @param j Variable selection. Allows \code{\link[data.table]{data.table}} style selection/variable creation (provided the data is also a data.table).
#' @param ... Other arguments to be passed to \code{[.data.frame} or \code{\link[data.table]{data.table}} (or whatever the class of the initial data).
#'
#' @details
#' If the original data was also a data.table, some calls to \code{[.fixest_panel} may dissolve the \code{fixest_panel} object and return a regular data.table. This is the case for subselections with additional arguments. If so, a note is displayed on the console.
#'
#' @return
#' It returns a \code{fixest_panel} data base, with the attributes allowing to create lags/leads properly bookkeeped.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' Alternatively, the function \code{\link[fixest]{panel}} changes a \code{data.frame} into a panel from which the functions \code{l} and \code{f} (creating leads and lags) can be called. Otherwise you can set the panel 'live' during the estimation using the argument \code{panel.id} (see for example in the function \code{\link[fixest]{feols}}).
#'
#' @examples
#'
#' data(base_did)
#'
#' # Creating a fixest_panel object
#' pdat = panel(base_did, ~id+period)
#'
#' # Subselections of fixest_panel objects bookkeeps the leads/lags engine
#' pdat_small = pdat[!pdat$period %in% c(2, 4), ]
#' a = feols(y~l(x1, 0:1), pdat_small)
#'
#' # we obtain the same results, had we created the lags "on the fly"
#' base_small = base_did[!base_did$period %in% c(2, 4), ]
#' b = feols(y~l(x1, 0:1), base_small, panel.id = ~id+period)
#' etable(a, b)
#'
#'
#' # Using data.table to create new lead/lag variables
#' if(require("data.table")){
#'   pdat_dt = panel(as.data.table(base_did), ~id+period)
#'
#'   # Variable creation
#'   pdat_dt[, x_l1 := l(x1)]
#'   pdat_dt[, c("x_l1", "x_f1_2") := .(l(x1), f(x1)**2)]
#'
#'   # Estimation on a subset of the data
#'   #  (the lead/lags work appropriately)
#'   feols(y~l(x1, 0:1), pdat_dt[!period %in% c(2, 4)])
#' }
#'
#'
"[.fixest_panel" = function(x, i, j, ...){
    # we need to perform proper bookkeeping

    info = attr(x, "panel_info")
    mc = match.call()

    IS_DT = FALSE
    if("data.table" %in% class(x) && requireNamespace("data.table", quietly=TRUE)){
        IS_DT = TRUE
        # data.table is really hard to handle....
        # Not very elegant... but that's life!

        # When modifications are too risky, I dissolve the

        mc_new = mc
        mc_new[[1]] = as.name('[')

        if(grepl(":=", deparse(mc$j)[1])){
            # Variable creation is OK

            jvalue = mc$j
            jtext = deparse_long(jvalue)
            # we check if f() or l() is used
            if(grepl("[^\\._[:alnum:]](l|f)\\(", jtext)){
                # We authorize it but only in 'simple' variable creation
                if(any(!names(mc) %in% c("", "x", "j"))){
                    # We don't allow i either
                    stop("When creating lags (resp. leads) with the function l() (resp. f()) within a data.table, only the argument 'j' is allowed.\nExample: 'data[, x_lag := l(x)]' is OK, while 'data[x>0, x_lag := l(x)]' is NOT ok.")
                }
                options(fixest_fl_authorized = TRUE)
                on.exit(options(fixest_fl_authorized = FALSE))
            }

            # I have to do it that way... really not elegant...
            mc_new[[1]] = as.name('[.data.table')
            eval(mc_new, asNamespace("data.table"), parent.frame())

            return(invisible(TRUE))
        } else {

            x_dt = data.table::copy(x)
            data.table::setattr(x_dt, "class", setdiff(class(x), "fixest_panel"))
            mc_new$x = as.name('x_dt')

            res = eval(mc_new)

            if(any(!names(mc) %in% c("", "x", "i"))) {
                # If any argument other than i is used => we dissolve the fixest panel
                message("NOTE: The fixest panel is dissolved.")
                return(res)
            } else {
                data.table::setattr(res, "class", c("fixest_panel", class(res)))
            }

        }
    } else {
        res = "[.data.frame"(x, i, j, ...)
    }

    if(!missing(i)){
        if(IS_DT){
            # data.table is quite a pain in the neck to handle...

            data.table::set(x_dt, j = "x__ID__x", value = 1:nrow(x_dt))
            select = eval(parse(text = paste0("x_dt[", deparse_long(mc_new$i), "]")))$x__ID__x

        } else {
            select = i
        }

        # we modify the indexes
        id = info$id_sorted[info$order_inv]
        time = info$time_sorted[info$order_inv]

        if(info$na_flag == FALSE){
            id = id[select]
            time = time[select]
        } else{
            # We don't forget to add the NAs!
            id_tmp = time_tmp = rep(NA, info$is_na)
            id_tmp[!info$is_na] = id
            time_tmp[!info$is_na] = time
            id = id_tmp[select]
            time = time_tmp[select]
        }

        is_na = is.na(id) | is.na(time)
        na_flag = FALSE
        if(any(is_na)){
            na_flag = TRUE
            id = id[!is_na]
            time = time[!is_na]
        }

        order_it = order(id, time)
        order_inv = order(order_it)

        new_info = list(order_it = order_it, order_inv=order_inv, id_sorted=id[order_it], time_sorted=time[order_it], na_flag = na_flag, panel.id = info$panel.id, call = info$call)
        if(na_flag) new_info$is_na = is_na
        attr(res, "panel_info") = new_info
    }

    # if(is.null(dim(res))){
    #     class(res) = "fixest_panel_vector"
    # }

    return(res)
}

terms_hat = function(fml, fastCombine = TRUE){

    fml_char = as.character(fml[length(fml)])
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

    t = terms(fml)

    t
}

check_lag = function(fml){
    fml_txt = deparse_long(fml)
    grepl("((^)|[^\\._[:alnum:]])(l|f)\\(", fml_txt)
}


reformulate_terms = function(x){
    # x: a vector of terms: x = attr(terms(fml), "term.labels")

    res = x

    if(any(grepl("\\b(l|f)\\(", x))){
        # Means lagging required

        term_all_expand = gsub("^l\\(", "l_expand(", x)
        term_all_expand = gsub("([^\\._[:alnum:]])l\\(", "\\1l_expand(", term_all_expand)

        term_all_expand = gsub("^f\\(", "f_expand(", term_all_expand)
        term_all_expand = gsub("([^\\._[:alnum:]])f\\(", "\\1f_expand(", term_all_expand)

        qui_expand = which(grepl("(l|f)_expand", term_all_expand))

        if(length(qui_expand) > 0){
            terms_all_list = as.list(x)

            # we select only the first function
            end_l_expand = function(x){
                x_split = strsplit(x, "")[[1]]
                open = x_split == "("
                close = x_split == ")"
                which.max(cumsum(open) > 0 & (cumsum(open) - cumsum(close)) == 0)
            }

            for(k in qui_expand){

                terms_split = strsplit(term_all_expand[k], "(f|l)_expand")[[1]]
                is_l = strsplit(term_all_expand[k], "_expand")[[1]]
                is_l = grepl("l$", is_l[-length(is_l)])

                slice = terms_split[1]
                for(i in 1:(length(terms_split) - 1)){

                    val = terms_split[i + 1]

                    end = end_l_expand(val)
                    if(end == nchar(val)){
                        quoi = eval(parse(text = paste0(ifelse(is_l[i], "l", "f"), "_expand", val)))
                    } else {
                        what = paste0(ifelse(is_l[i], "l", "f"), "_expand", substr(val, 1, end))
                        end_text = substr(val, end + 1, nchar(val))
                        quoi = paste0(eval(parse(text = what)), end_text)
                    }

                    slice = paste0(rep(slice, each = length(quoi)), quoi)

                }

                terms_all_list[[k]] = slice

            }

            res = unlist(terms_all_list)
        }

    }

    res
}


rewrite_fml = function(x){

    x_FML = Formula(x)

    # If the user wants to combine the LHS: so be it...
    lhs_fml = paste(reformulate_terms(deparse_long(x_FML[[2]])), collapse = "+")

    rhs_terms = attr(terms(formula(x_FML, lhs = 1, rhs = 1)), "term.labels")
    rhs_fml = paste(reformulate_terms(rhs_terms), collapse = "+")

    if(length(x_FML)[2] > 1){
        fixef_terms = attr(terms(formula(x_FML, lhs = 1, rhs = 2)), "term.labels")
        fixef_fml = paste("|", paste(reformulate_terms(fixef_terms), collapse = "+"))
    } else {
        fixef_fml = ""
    }


    as.formula(paste0(lhs_fml, "~", rhs_fml, fixef_fml))
}



