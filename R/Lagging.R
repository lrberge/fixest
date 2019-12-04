#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sun Nov 24 09:18:19 2019
# ~: Lagging
#----------------------------------------------#

# ------------------------------------------------------------- #
# Contains all functions to ensure the easy lagging of variable
# within fixest estimations.
#
# ------------------------------------------------------------- #


# Work in progress:
# - include the possibility to use var1^var2 -- OK
# - include NA handling in bookeeping
# - drop the methods of extraction ([[ and $)... not very useful, l/f should not be applied to vectors

panel_setup = function(data, panel.id, time.step = "unitary", duplicate.method = c("none", "first"), DATA_MISSING = FALSE, from_fixest = FALSE){
    # Function to setup the panel.
    # Used in lag.formula, panel, and fixest_env (with argument panel.id and panel.args)
    # DATA_MISSING: arg used in lag.formula

    # The id/time -- Beware when using var1^var2

    duplicate.method = match.arg(duplicate.method)

    tm = terms_hat(panel.id)
    var_id_time = attr(tm, "term.labels")
    if(length(var_id_time) != 2){
        stop("The formula of the panel IDs must contain exactly two variables in the right hand side (currently there ", ifsingle(var_id_time, "is ", "are "), length(var_id_time), ").")
    }

    if(DATA_MISSING){
        if(!all(all.vars(tm) %in% ls(parent.frame(2)))){
            pblm = setdiff(all.vars(tm), ls(parent.frame(2)))
            stop("In the panel IDs, the variable", enumerate_items(pblm, "s.is.past.quote"), " not found in the environment.")
        }
        id = eval(parse(text = var_id_time[1]), parent.frame(2))
        time = eval(parse(text = var_id_time[2]), parent.frame(2))
    } else {
        if(!all(all.vars(tm) %in% names(data))){
            pblm = setdiff(all.vars(tm), ls(parent.frame(2)))
            stop("In the panel IDs, the variable", enumerate_items(pblm, "s.is.quote"), " not in the data set.")
        }
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
                stop("To use the 'unitary' time.step, the time variable must be numeric or at least convertible to numeric. So far the conversion has failed (time variable's class is currently ", enumerate_items(class(time)), ").")
            }

            time = time_new
        }

        if(any(time %% 1 != 0)){
            stop("To use the 'unitary' time.step, the time variable", ifelse(time_conversion, " (which has been converted to numeric)", ""), " must be made of integers. So far this is not the case. Alternatively, you can give a number in time.step.")
        }

    } else if(time.step != "consecutive"){
        if(!is.numeric(time)){
            stop("If 'time.step' is a number, then the time variable must also be a number (this is not the case: its class is currently ", enumerate_items(class(time)), ").")
        }

        # if(any(time %% time.step != 0)){
        #     pblm = unique(head(time[time %% time.step != 0], 3))
        #     stop("If 'time.step' is a number, then it must be an exact divisor of all the values in the time variable. This is currently not the case: ", time.step, " is not a divisor of ", enumerate_items(pblm, or = TRUE, verb = FALSE), ".")
        # }
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
    if(time.step == "consecutive"){
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

            stop("If 'time.step' is a number, then it must be an exact divisor of all the difference between two consecutive time periods. This is currently not the case: ", time.step, " is not a divisor of ", all_steps[obs_pblm][1], " (the difference btw the time periods ", time_unik[obs_pblm[1] + 1], " and ", time_unik[obs_pblm[1]], ").")
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

            stop("The panel identifiers contain duplicate values: this is not allowed since lag/leads are not defined for them. For example (id, time) = (", id_dup, ", ", time_dup, ") appears ", n_times(dup_info$n_dup), ". Please provide data without duplicates -- or you can also use duplicate.method = 'first' (see Details).")
        }
    }

    res = list(order_it = order_it, order_inv = order_inv, id_sorted = id_sorted, time_sorted = time_sorted, na_flag = na_flag)
    if(na_flag) res$is_na = is_na
    res
}

# Add argument fill to f and l
f = function(x, lead = 1, fill = NA){
    l(x, -lead, fill)
}


l = function(x, lag = 1, fill = NA){

    value = x

    sys_calls = rev(sys.calls())

    # To improve => you don't wanna check all frames, only the relevant ones
    from_fixest = FALSE
    for(where in 1:5){
        if(exists("panel__meta__info", parent.frame(where))){
            from_fixest = TRUE
            meta_info = get("panel__meta__info", parent.frame(where))
        }
    }

    if(from_fixest == FALSE){
        if(any(grepl("^.\\[\\.data\\.table", sys_calls))){
            # this is a call from within data table => OK
            # Getting the meta info is a bit tricky

            qui  = which.max(grepl("^.\\[\\.data\\.table", sys_calls))
            var = gsub("^[^\\(]+\\(|,.+", "", sys_calls[qui])
            m = eval(parse(text = var), envir = sys.frames()[qui])

            if(!"fixest_panel" %in% class(m)){
                stop("You can use l() or f() only when the data set is of class 'fixest_panel', you can use function panel() to set it.")
            }

            meta_info = attr(m, "panel_info")
        } else {
            stop("Function l() is only callable within 'fixest' estimations or within a variable creation with data.table where the data set is a fixest_panel (obtained from panel()).")
        }
    }

    # we get the observation id!
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
        existing_vars = ls(parent.frame(2))
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



panel = function(data, panel.id, time.step = "unitary", duplicate.method = c("none", "first")){

    if(missing(data)){
        stop("You must provide the argument 'data'.")
    } else if(!"data.frame" %in% class(data)){
        stop("Argument 'data' must be a data.frame.")
    }

    duplicate.method = match.arg(duplicate.method)
    check_arg(panel.id, "osf")

    meta_info = panel_setup(data, panel.id = panel.id, time.step = time.step, duplicate.method = duplicate.method)

    # R makes a shallow copy of data => need to do it differently with DT

    if("data.table" %in% class(data) && isNamespaceLoaded("data.table")){
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

unpanel = function(x){

    if("data.table" %in% class(x) && isNamespaceLoaded("data.table")){
        data.table::setattr(x, "panel_info", NULL)
        data.table::setattr(x, "class", setdiff(class(x), "fixest_panel"))
    } else {
        attr(x, "panel_info") = NULL
        class(x) = setdiff(class(x), "fixest_panel")
    }

    x
}


# "$.fixest_panel" = function(x, name){
#     x[[name]]
# }
#
# "[[.fixest_panel" = function(x, name){
#     res = "[[.data.frame"(x, name)
#     if(!is.null(res)){
#         attr(res, "panel_info") = attr(x, "panel_info")
#     }
#
#     class(res) = "fixest_panel_vector"
#
#     res
# }

"[.fixest_panel" = function(x, i, j, ...){
    # we need to perform proper bookkeeping

    info = attr(x, "panel_info")
    mc = match.call()

    if("data.table" %in% class(x)){
        # data.table is really hard to handle....
        # Not very elegant... but that's life!

        mc_new = mc
        mc_new[[1]] = as.name('[')

        if(grepl(":=", deparse(mc))){
            # Pfff data.table is REALLY a pain in the neck
            # doesn't look 100% safe (say 95%)

            # Does NOT work:
            # x_name = deparse(mc$x)
            # eval(parse(text = paste0("setDT(", x_name, ")")), parent.frame())
            # res <- eval(mc_new, parent.frame())
            # eval(parse(text = paste0("setattr(", x_name, ", 'class', c('panel', class(", x_name, ")))")), parent.frame())


            # does not work:
            # return(NextMethod())

            mc_new[[1]] = as.name('[.data.table')
            # eval(mc_new, environment(data.table), parent.frame()) # Dont Work
            # eval(mc_new, as.environment("package:data.table"), parent.frame()) # Dont Work
            eval(mc_new, asNamespace("data.table"), parent.frame())
            # does not work:
            # eval(data.table:::"[.data.table"(x, i, j, ...), parent.frame())

            return(invisible(NULL))
        } else {
            x_dt = copy(x)
            setattr(x_dt, "class", setdiff(class(x), "panel"))
            mc_new$x = as.name('x_dt')
            res = eval(mc_new)
            setattr(res, "class", c("fixest_panel", class(res)))
        }

    } else {
        res = "[.data.frame"(x, i, j, ...)
    }

    if(!missing(i)){
        if("data.table" %in% class(x)){
            # data.table is quite a pain in the neck to handle...

            x_dt[, xxxIDxxx := 1:.N]
            mc_new$j = as.name("xxxIDxxx")
            select = eval(mc_new)

        } else {
            select = i
        }

        # we modify the indexes
        id = info$id_sorted[info$order_inv]
        time = info$time_sorted[info$order_inv]

        id = id[select]
        time = time[select]

        order_it = order(id, time)
        order_inv = order(order_it)
        attr(res, "panel_info") = list(order_it = order_it, order_inv=order_inv, id_sorted=id[order_it], time_sorted=time[order_it])
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



