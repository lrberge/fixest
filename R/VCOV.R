#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Thu Apr 01 10:25:53 2021
# ~: VCOV related stuff
#----------------------------------------------#


####
#### User-level ####
####



#' Computes the variance/covariance of a \code{fixest} object
#'
#' This function extracts the variance-covariance of estimated parameters from a model estimated with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}.
#'
#' @inheritParams summary.fixest
#' @inheritParams nobs.fixest
#'
#' @param attr Logical, defaults to \code{FALSE}. Whether to include the attributes describing how the VCOV was computed.
#' @param ... Other arguments to be passed to \code{\link[fixest]{summary.fixest}}.
#'
#' The computation of the VCOV matrix is first done in \code{\link[fixest]{summary.fixest}}.
#'
#' @details
#' For an explanation on how the standard-errors are computed and what is the exact meaning of the arguments, please have a look at the dedicated vignette: \href{https://lrberge.github.io/fixest/articles/standard_errors.html}{On standard-errors}.
#'
#' @seealso
#' You can also compute VCOVs with the following functions: \code{\link[fixest]{vcov_cluster}}, \code{\link[fixest]{vcov_hac}}, \code{\link[fixest]{vcov_conley}}.
#'
#' @return
#' It returns a \eqn{K\times K} square matrix where \eqn{K} is the number of variables of the fitted model.
#' If \code{attr = TRUE}, this matrix has an attribute \dQuote{type} specifying how this variance/covariance matrix has been computed.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{confint.fixest}}, \code{\link[fixest]{resid.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # Load panel data
#' data(base_did)
#'
#' # Simple estimation on a panel
#' est = feols(y ~ x1, base_did)
#'
#' # ======== #
#' # IID VCOV #
#' # ======== #
#'
#' # By default the VCOV assumes iid errors:
#' se(vcov(est))
#'
#' # You can make the call for an iid VCOV explicitly:
#' se(vcov(est, "iid"))
#'
#' #
#' # Heteroskedasticity-robust VCOV
#' #
#'
#' # By default the VCOV assumes iid errors:
#' se(vcov(est, "hetero"))
#'
#' # => note that it also accepts vcov = "White" and vcov = "HC1" as aliases.
#'
#' # =============== #
#' # Clustered VCOVs #
#' # =============== #
#'
#' # To cluster the VCOV, you can use a formula of the form cluster ~ var1 + var2 etc
#' # Let's cluster by the panel ID:
#' se(vcov(est, cluster ~ id))
#'
#' # Alternative ways:
#'
#' # -> cluster is implicitly assumed when a one-sided formula is provided
#' se(vcov(est, ~ id))
#'
#' # -> using the argument cluster instead of vcov
#' se(vcov(est, cluster = ~ id))
#'
#' # For two-/three- way clustering, just add more variables:
#' se(vcov(est, ~ id + period))
#'
#' # -------------------|
#' # Implicit deduction |
#' # -------------------|
#' # When the estimation contains FEs, the dimension on which to cluster
#' # is directly inferred from the FEs used in the estimation, so you don't need
#' # to explicitly add them.
#'
#' est_fe = feols(y ~ x1 | id + period, base_did)
#'
#' # Clustered along "id"
#' se(vcov(est_fe, "cluster"))
#'
#' # Clustered along "id" and "period"
#' se(vcov(est_fe, "twoway"))
#'
#'
#' # =========== #
#' # Panel VCOVs #
#' # =========== #
#'
#' # ---------------------|
#' # Newey West (NW) VCOV |
#' # ---------------------|
#' # To obtain NW VCOVs, use a formula of the form NW ~ id + period
#' se(vcov(est, NW ~ id + period))
#'
#' # If you want to change the lag:
#' se(vcov(est, NW(3) ~ id + period))
#'
#' # Alternative way:
#'
#' # -> using the vcov_NW function
#' se(vcov(est, vcov_NW(unit = "id", time = "period", lag = 3)))
#'
#' # -------------------------|
#' # Driscoll-Kraay (DK) VCOV |
#' # -------------------------|
#' # To obtain DK VCOVs, use a formula of the form DK ~ period
#'
#' se(vcov(est, DK ~ period))
#'
#' # If you want to change the lag:
#' se(vcov(est, DK(3) ~ period))
#'
#' # Alternative way:
#'
#' # -> using the vcov_DK function
#' se(vcov(est, vcov_DK(time = "period", lag = 3)))
#'
#' # -------------------|
#' # Implicit deduction |
#' # -------------------|
#' # When the estimation contains a panel identifier, you don't need
#' # to re-write them later on
#'
#' est_panel = feols(y ~ x1, base_did, panel.id = ~id + period)
#'
#' # Both methods, NM and DK, now work automatically
#' se(vcov(est_panel, "NW"))
#' se(vcov(est_panel, "DK"))
#'
#'
#' # =================================== #
#' # VCOVs robust to spatial correlation #
#' # =================================== #
#'
#' data(quakes)
#' est_geo = feols(depth ~ mag, quakes)
#'
#' # ------------|
#' # Conley VCOV |
#' # ------------|
#' # To obtain a Conley VCOV, use a formula of the form conley(cutoff) ~ lat + lon
#' # with lat/lon the latitude/longitude variable names in the data set
#' se(vcov(est_geo, conley(100) ~ lat + long))
#'
#' # Alternative way:
#'
#' # -> using the vcov_DK function
#' se(vcov(est_geo, vcov_conley(lat = "lat", lon = "long", cutoff = 100)))
#'
#' # -------------------|
#' # Implicit deduction |
#' # -------------------|
#' # By default the latitude and longitude are directly fetched in the data based
#' # on pattern matching. So you don't have to specify them.
#' # Furhter, an automatic cutoff is deduced by default.
#'
#' # The following works:
#' se(vcov(est_geo, "conley"))
#'
#'
#' # ======================== #
#' # Small Sample Corrections #
#' # ======================== #
#'
#' # You can change the way the small sample corrections are done with the argument ssc.
#' # The argument ssc must be created by the ssc function
#' se(vcov(est, ssc = ssc(adj = FALSE)))
#'
#' # You can add directly the call to ssc in the vcov formula.
#' # You need to add it like a variable:
#' se(vcov(est, iid ~ ssc(adj = FALSE)))
#' se(vcov(est, DK ~ period + ssc(adj = FALSE)))
#'
#'
#'
vcov.fixest = function(object, vcov = NULL, se = NULL, cluster, ssc = NULL, attr = FALSE, forceCovariance = FALSE,
                       keepBounded = FALSE, nthreads = getFixest_nthreads(), ...){
    # computes the clustered vcov

    check_arg(attr, "logical scalar")
    is_attr = attr

    dots = list(...)

    # START :: SECTION used only internally in fixest_env
    only_varnames = isTRUE(dots$only_varnames)
    data_names = dots$data_names
    if(only_varnames){
        # Used internally in fixest_env to find out which variable to keep
        # => we need panel.id, so we can remove the NAs from it if it is implicitly deduced to be used
        # => idem for fixef_vars
        object = list(panel.id = dots$panel.id, fixef_vars = dots$fixef_vars)
    }
    #   END


    if(isTRUE(object$NA_model)){
        # means that the estimation is done without any valid variable
        return(object$cov.scaled)
    }

    if("dof" %in% names(dots)){
        if(is.null(getOption("fixest_warn_dof_arg"))){
            warning("The argument 'dof' is deprecated. Please use 'ssc' instead.")
            options(fixest_warn_dof_arg = TRUE)
        }
        ssc = dots$dof
    }

    if(is_user_level_call()){
        if(!is_function_in_it(vcov)){
            validate_dots(suggest_args = c("vcov", "ssc"), valid_args = "dof")
        }
    }

    # All the available VCOVs
    all_vcov = getOption("fixest_vcov_builtin")
    all_vcov_names = unlist(lapply(all_vcov, `[[`, "name"))
    all_vcov_names = all_vcov_names[nchar(all_vcov_names) > 0]

    # We transform se and cluster into vcov
    vcov_vars = var_values_all = var_names_all = NULL
    vcov = oldargs_to_vcov(se, cluster, vcov)

    sandwich = !isFALSE(dots$sandwich) # I know it's weird

    if(!is.null(object$onlyFixef)){
        # means that the estimation is done without variables
        return(NULL)
    }

    # If it's a summary => we give the vcov directly without further computation! except if arguments are provided which would mean that the user wants the new vcov
    if(isTRUE(object$summary) && missnull(vcov) && missnull(ssc)){
        vcov = object$cov.scaled
        if(!is_attr) {
            all_attr = names(attributes(vcov))
            for(v in setdiff(all_attr, c("dim", "dimnames"))){
                attr(vcov, v) = NULL
            }
        }
        return(vcov)
    }


    # Default behavior vcov:
    suffix = ""
    if(missnull(vcov)){

        vcov_default = getFixest_vcov()

        if(!is.null(object$panel.id)){
            # Panel has precedence over FEs
            vcov = vcov_default$panel

        } else {
            n_fe = length(object$fixef_id)

            if(n_fe == 0){
                vcov = vcov_default$no_FE
            } else if(n_fe == 1){
                vcov = vcov_default$one_FE
            } else {
                vcov = vcov_default$two_FE
            }

            if(!is.null(object$fixef_sizes) && object$fixef_sizes[1] == 1){
                # Special case => cleaner output
                vcov = vcov_default$no_FE
            }
        }

    }


    if(isTRUE(object$lean)){
        # we can't compute the SE because scores are gone!
        # LATER: recompute the scores (costly but maybe only possibility for user?)
        #
        # so only IID VCOV is valid

        ok = FALSE
        if(is.character(vcov) && length(vcov) == 1){
            vcov_id = which(sapply(all_vcov, function(x) vcov %in% x$name))

            if(length(vcov_id) == 1){
                vcov_select = all_vcov[[vcov_id]]
                ok = identical(vcov_select$vcov_label, "IID")
            }
        }

        if(!ok) stop("VCOV of 'lean' fixest objects cannot be computed. Please re-estimate with 'lean = FALSE'.")
    }

    ####
    #### ... vcov parsing ####
    ####

    # Checking the value of vcov
    check_arg_plus(vcov, "match | formula | function | matrix | list len(1) | class(fixest_vcov_request)", .choices = all_vcov_names)

    user_vcov_name = NULL
    if(is.list(vcov) && !inherits(vcov, "fixest_vcov_request")){
        # We already ensured it was a list of length 1
        user_vcov_name = names(vcov)
        vcov = vcov[[1]]
    }

    if(is.function(vcov)){

        if(only_varnames){
            return(character(0))
        }

        # cleaning dots
        for(var in c("summary_flags")){
            dots[[var]] = NULL
        }

        # The name
        arg_list = dots
        if(".vcov_args" %in% names(dots)){
            # internal call
            vcov_name = dots$vcov_name
            arg_list = dots$.vcov_args
        } else {
            vcov_name = attr(vcov, "deparsed_arg")
            if(is.null(vcov_name)){
                vcov_name = fetch_arg_deparse("vcov")
            } else {
                # cleaning
                attr(vcov, "deparsed_arg") = NULL
            }

            # Getting the right arguments
            arg_list = catch_fun_args(vcov, arg_list, exclude_args = "vcov", keep_dots = TRUE)
        }

        if(!is.null(user_vcov_name)){
            vcov_name = user_vcov_name

        } else {
            vcov_name = gsub("sandwich::", "", vcov_name, fixed = TRUE)
            vcov_name = gsub("^function\\([^\\)]+\\) ", "", vcov_name)
            vcov_name = gsub("^vcov$", "Custom", vcov_name)
        }

        # We shouldn't have a prior on the name of the first argument
        arg_names = formalArgs(vcov)
        arg_list[[arg_names[1]]] = object

        vcov = do.call(vcov, arg_list)

        n_coef = length(object$coefficients)
        check_value(vcov, "square numeric matrix nrow(value)", .value = n_coef,
                    .message = paste0("If argument 'vcov' is to be a function, it should return a square numeric matrix of the same dimension as the number of coefficients (here ", n_coef, ")."))

        # We add the type of the matrix
        attr(vcov, "type") = vcov_name
        attr(vcov, "dof.K") = object$nparams

        return(vcov)
    }

    if(is.matrix(vcov)){

        if(only_varnames){
            stop("A custom VCOV matrix cannot be used directly in a fixest estimation.")
        }

        # Check that this makes sense
        n_coef = length(object$coefficients)
        check_value(vcov, "square matrix nrow(value)", .value = n_coef)

        attr(vcov, "type") = if(is.null(user_vcov_name)) "Custom" else user_vcov_name

        attr(vcov, "dof.K") = object$nparams

        return(vcov)
    }

    extra_args = NULL
    if(inherits(vcov, "fixest_vcov_request")){
        if(!is.null(vcov$ssc)) ssc = vcov$ssc
        var_names_all = vcov$var_names_all
        var_values_all = vcov$vcov_vars
        extra_args = vcov$extra_args
        vcov = vcov$vcov

    }

    if(inherits(vcov, "formula")){

        vcov_fml = vcov

        if(length(vcov_fml) == 2){
            vcov = ""
            vcov_vars = fml2varnames(vcov_fml, combine_fun = TRUE)
        } else {
            vcov = deparse_long(vcov_fml[[2]])

            is_extra = grepl("(", vcov, fixed = TRUE)
            if(is_extra){
                vcov = trimws(gsub("\\(.+", "", vcov))
            }

            check_arg_plus(vcov, "match", .message = "If a formula, the arg. 'vcov' must be of the form 'vcov_type ~ vars'. The vcov_type must be a supported VCOV type.", .choices = all_vcov_names)

            if(is_extra){
                new_req = eval(vcov_fml[[2]], environment(vcov_fml))
                extra_args = new_req$extra_args
            }

            vcov_vars = fml2varnames(vcov_fml[c(1, 3)], combine_fun = TRUE)
        }

        qui_ssc = grepl("ssc(", vcov_vars, fixed = TRUE)
        if(any(qui_ssc)){
            ssc_txt = vcov_vars[qui_ssc]
            ssc = eval(str2lang(ssc_txt))
            vcov_vars = vcov_vars[!qui_ssc]
        }

    }

    # Here vcov **must** be a character scalar

    vcov_id = which(sapply(all_vcov, function(x) vcov %in% x$name))

    if(length(vcov_id) != 1){
        stop("Unexpected problem in the selection of the VCOV. This is an internal error. Could you report?")
    }

    vcov_select = all_vcov[[vcov_id]]

    if(length(vcov_select$vars) == 0){
        var_names_all = character(0)

        if(only_varnames){
            return(var_names_all)
        }
    } else {

        if(!is.null(var_values_all)){
            # If we're here, this means that vcov_vars
            # has been provided via a vcov_request:
            #   - either via retro-compatibility (use of "cluster" with a vector or list)
            #   - either via feeding vectors into the user level version of the VCOV requested
            #
            # => so we don't have to do anything, checking will come later and is common with
            # the else{} of this condition

            is_int_all = list()

            if(only_varnames){
                return(character(0))
            }

            # We trim the observations if needed
            for(i in seq_along(var_values_all)){
                value = var_values_all[[i]]
                if(length(value) != object$nobs_origin){
                    stop("To compute the ", vcov_select$vcov_label, " VCOV, you need to provide variables with the same number of observations as in the original data set. Currently there are ", length(value), " instead of ", object$nobs_origin, ".")
                }

                var_values_all[[i]] = trim_obs_removed(value, object)
            }


        } else {

            patterns_split = strsplit(vcov_select$patterns, " ?\\+ ?")
            n_patterns = lengths(patterns_split)

            if(isTRUE(object$is_fit)){
                # No automatic deduction for fit methods,except for clustered SEs
                n_FE = length(object$fixef_vars)
                if(vcov_select$vcov_label %in% "Clustered"){
                    if(n_FE == 0){
                        stop("To compute a clustered VCOV from a '.fit' estimation, you need to provide the variables over which to cluster with the function vcov_cluster(). E.g. vcov = vcov_cluster(cluster = df) with df a data.frame containing the clusters. Or make a regular estimation, i.e. not from '.fit'.")
                    }
                    if(!n_FE %in% n_patterns){
                        stop("To compute a clustered VCOV from a '.fit' estimation, you need to provide the variables over which to cluster with the function vcov = vcov_cluster(cluster = df), with df a data.frmae containing the clusters. Or make a regular estimation, i.e. not from '.fit'.")
                    }
                } else {
                    stop("The estimations from '.fit' methods don't support ", vcov_select$vcov_label, " VCOVs. To use them, perform a regular estimation instead.")
                }
            }

            if(only_varnames == FALSE){
                data = fetch_data(object)
                data_names = names(data)
            }


            if(!length(vcov_vars) %in% n_patterns){
                fml_display = paste0("~", paste0(vcov_select$patterns[n_patterns != 0], collapse = "+"))
                stop("In the argument 'vcov', the number of variables in the RHS of the formula (", length(vcov_vars), ") is not valid. The formula should correspond to ", ifunit(n_patterns, "", "one of "), enumerate_items(fml_display, "or"), ".")
            }


            var_values_all = list()
            # => list of variables used to compute the VCOV
            var_names_all = c()
            # => same as var_values_all but the variable names

            is_int_all = list()
            # => flag to avoid reapplying to_integer

            pattern = patterns_split[[which(n_patterns == length(vcov_vars))]]

            # We find all the variable names and then evaluate them
            for(i in seq_along(vcov_select$vars)){
                vcov_var_name = names(vcov_select$vars)[i]
                vcov_var_value = vcov_select$vars[[i]]

                if(vcov_var_name %in% pattern){
                    # => provided by the user, we find which one it corresponds to
                    # based on the pattern

                    vname = vcov_vars[which(pattern == vcov_var_name)]

                    vname_all = all.vars(str2expression(vname))
                    if(!all(vname_all %in% data_names)){
                        pblm = setdiff(vname_all, data_names)
                        stop("The variable", enumerate_items(pblm, "s.quote"), " used to compute the VCOV ", plural_len(pblm, "is"), " not in the original data set. Only variables in the data set can be used.")
                    }

                    var_names_all[vcov_var_name] = rename_hat(vname)
                    if(only_varnames) {
                        var_names_all[vcov_var_name] = vname_all[1]
                        # To handle combined clusters
                        var_names_all = c(var_names_all, vname_all[-1])
                        next
                    }

                    value = eval(str2expression(vname), data)
                    var_values_all[[vcov_var_name]] = trim_obs_removed(value, object)

                } else {

                    ####
                    #### ... --- guessing the variables ####
                    ####

                    # not provided by the user: GUESSED!
                    # 3 types of guesses:
                    # - fixef (fixed-effects used in the estimation)
                    # - panel.id (panel identifiers used in the estimation)
                    # - regex (based on variable names)
                    #


                    guesses = vcov_var_value$guess_from
                    # err_msg, vector of elements of the form c("expect" = "pblm")
                    err_msg = c()
                    for(k in seq_along(guesses)){
                        type = names(guesses)[k]

                        if(type == "fixef"){

                            if(is.null(object$fixef_vars)){
                                msg = setNames(nm = "the fixed-effects",
                                               "there was no fixed-effect in the estimation")
                                err_msg = c(err_msg, msg)
                                next
                            }

                            if(length(object$fixef_vars) < guesses$fixef){
                                msg = setNames(nm = paste0("the ", n_th(guesses$fixef), " fixed-effect"),
                                               paste0("the estimation had only ", n_letter(length(object$fixef_vars)), " fixed-effect", plural_len(object$fixef_vars)))
                                err_msg = c(err_msg, msg)
                                next
                            }

                            var_names_all[vcov_var_name] = object$fixef_vars[guesses$fixef]
                            if(only_varnames) break

                            is_int_all[[vcov_var_name]] = TRUE

                            var_values_all[[vcov_var_name]] = object$fixef_id[[guesses$fixef]]
                            break

                        } else if(type == "panel.id"){

                            if(is.null(object$panel.id)){
                                msg = setNames(nm = "the 'panel.id' identifiers",
                                               "no 'panel.id' was set in this estimation")
                                err_msg = c(err_msg, msg)
                                next
                            }

                            vname = object$panel.id[guesses$panel.id]

                            vname_all = all.vars(str2expression(vname))
                            if(!all(vname_all %in% data_names)){
                                pblm = setdiff(vname_all, data_names)
                                msg = setNames(nm = "the 'panel.id' identifiers",
                                               paste0("the variable", enumerate_items(pblm, "s"), " set in 'panel.id' ", plural_len(pblm, "is"), " not in the data set"))
                                err_msg = c(err_msg, msg)
                                next
                            }

                            # ATTENTION:
                            # if the panel.id variable has NA values, and there are no lags used in the estimation
                            # (which would trigger obs removal), then we cannot use the panel.id
                            # for the VCOV reliably since the sample would differ from the sample used in the
                            # estimation

                            var_names_all[vcov_var_name] = vname
                            if(only_varnames) break

                            value = eval(str2expression(vname), data)
                            var_values_all[[vcov_var_name]] = trim_obs_removed(value, object)
                            break

                        } else if(type == "regex"){

                            ok = FALSE
                            msg = NULL
                            for(pat in guesses$regex){
                                var_id = which(grepl(pat, data_names))

                                if(length(var_id) == 0){
                                    if(is.null(msg)){
                                        msg = setNames(nm = "the variable names of the data set",
                                                       paste0("no match was found for ", vcov_var_value$label))
                                    }
                                } else if(length(var_id) == 1){
                                    ok = TRUE
                                    vname = data_names[var_id]
                                    break
                                } else {
                                    msg = setNames(nm = "the variable names of the data set",
                                                   paste0("several matches were found: ",
                                                          enumerate_items(data_names[var_id], "quote.enum i")))
                                }
                            }

                            if(ok){
                                var_names_all[vcov_var_name] = vname
                                if(only_varnames) break

                                var_values_all[[vcov_var_name]] = trim_obs_removed(data[[vname]], object)
                                break
                            } else {
                                err_msg = c(err_msg, msg)
                            }
                        }
                    }

                    if(is.null(var_names_all[vcov_var_name]) &&
                       length(err_msg) > 0 && !isTRUE(vcov_var_value$optional)){

                        if(length(err_msg) == 1){
                            expect = names(err_msg)
                            pblm = err_msg
                        } else {
                            expect = enumerate_items(names(err_msg), "enum a.or")
                            pblm = enumerate_items(err_msg, "enum a")
                        }

                        stop("To compute the ", vcov_select$vcov_label, " VCOV, we need a variable for the ",
                             vcov_var_value$label, ". Since you didn't provide it in the formula, we typically deduce it from ", expect,
                             ". Problem: ", pblm, ". Please provide it in the formula.")
                    }
                }
            }

            if(only_varnames){
                return(var_names_all)
            }
        }


        # Extra checking + putting into int if requested
        for(i in seq_along(var_values_all)){
            vcov_var_name = names(var_values_all)[i]
            vname = var_names_all[vcov_var_name]
            value = var_values_all[[i]]

            vcov_var_value = vcov_select$vars[[vcov_var_name]]

            if(!isTRUE(is_int_all[[vcov_var_name]]) && anyNA(value)){
                # First condition means: it is a fixed-effect used in the estimation => no need to check
                stop("The variable '", vname, "' used to estimate the VCOV (employed as ", vcov_var_value$label, ") has NA values which would lead to a sample used to compute the VCOV different from the sample used to estimate the parameters. This would lead to wrong inference.\nPossible solutions: i) ex ante prune them, ii) impute them, or iii) use the argument 'vcov' at estimation time.")
            }

            if(isTRUE(vcov_var_value$to_int)){
                # if to_int => we do not have to check the type since it can be applied to any type

                if(!isTRUE(is_int_all[[vcov_var_name]])){
                    var_values_all[[vcov_var_name]] = quickUnclassFactor(value)
                }

            } else if(!is.null(vcov_var_value$expected_type)){
                check_value(value, .type = vcov_var_value$expected_type,
                            .prefix = paste0("To compute the VCOV, the ", vcov_var_value$label, " (here '", vname, "') "))
            }
        }

        vcov_vars = var_values_all

    }


    # Checking the nber of threads
    if(!missing(nthreads)) nthreads = check_set_nthreads(nthreads)


    ####
    #### ... scores ####
    ####

    # We handle the bounded parameters:
    isBounded = object$isBounded
    if(is.null(isBounded)){
        isBounded = rep(FALSE, length(object$coefficients))
    }

    if(any(isBounded)){
        if(keepBounded){
            # we treat the bounded parameters as regular variables
            scores = object$scores
            object$cov.iid = solve(object$hessian)
        } else {
            scores = object$scores[, -which(isBounded), drop = FALSE]
        }
    } else {
        scores = object$scores
    }


    ####
    #### ... bread ####
    ####

    n = object$nobs

    if(object$method_type == "feols"){
        iid_names = all_vcov[[1]]$name
        if(vcov %in% iid_names){
            bread = object$cov.iid / ((n - 1) / (n - object$nparams))
        } else {
            bread = object$cov.iid / object$sigma2
        }
    } else {
        bread = object$cov.iid
    }

    if(anyNA(bread)){

        IS_NA_VCOV = FALSE

        if(!forceCovariance){
            IS_NA_VCOV = TRUE
        } else {
            info_inv = cpp_cholesky(object$hessian)
            if(!is.null(info_inv$all_removed)){
                # Means all variables are collinear! => can happen when using FEs
                IS_NA_VCOV = TRUE
            }
        }

        if(IS_NA_VCOV){
            attr(bread, "type") = "NA (not-available)"

            if(is_attr){
                attr(bread, "dof.K") = object$nparams
                attr(bread, "df.t") = NA
            }

            return(bread)
        }

        VCOV_raw_forced = info_inv$XtX_inv
        if(any(info_inv$id_excl)){
            n_collin = sum(info_inv$all_removed)
            message("NOTE: ", n_letter(n_collin), " variable", plural(n_collin, "s.has"), " been found to be singular.")

            VCOV_raw_forced = cpp_mat_reconstruct(VCOV_raw_forced, info_inv$id_excl)
            VCOV_raw_forced[, info_inv$id_excl] = NA
            VCOV_raw_forced[info_inv$id_excl, ] = NA
        }

        bread = VCOV_raw_forced

    }

    ####
    #### ... vcov no adj ####
    ####

    # we compute the vcov. The adjustment (which is a pain in the neck) will come after that
    # Here vcov is ALWAYS a character scalar

    # ssc related => we accept NULL
    # we check ssc since it can be used by the funs
    if(missnull(ssc)) ssc = getFixest_ssc()
    check_arg(ssc, "class(ssc.type)", .message = "The argument 'ssc' must be an object created by the function ssc().")

    fun_name = vcov_select$fun_name
    args = list(bread = bread, scores = scores, vars = vcov_vars, ssc = ssc,
                sandwich = sandwich, nthreads = nthreads,
                var_names_all = var_names_all)

    for(a in names(extra_args)){
        # I have to add this weird condition (because of the aliases)
        if(!is.null(extra_args[[a]])) args[[a]] = extra_args[[a]]
    }

    vcov_noAdj = do.call(fun_name, args)

    if(!sandwich){
        return(vcov_noAdj)
    }

    dimnames(vcov_noAdj) = dimnames(bread)

    ####
    #### ... ssc ####
    ####

    # ssc is a ssc object in here

    n_fe = n_fe_ok = length(object$fixef_id)

    # we adjust the fixef sizes to account for slopes
    fixef_sizes_ok = object$fixef_sizes
    isSlope = FALSE
    if(!is.null(object$fixef_terms)){
        isSlope = TRUE
        # The drop the fixef_sizes for only slopes
        fixef_sizes_ok[object$slope_flag < 0] = 0
        n_fe_ok = sum(fixef_sizes_ok > 0)
    }

    nested_vars = sapply(vcov_select$vars, function(x) isTRUE(x$rm_nested))
    any_nested_var = length(nested_vars) > 0 && length(vcov_vars) > 0 && any(nested_vars[names(vcov_vars)])

    if(ssc$fixef.K == "none"){
        # we do it with "minus" because of only slopes
        K = object$nparams
        if(n_fe_ok > 0){
            K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1))
        }
    } else if(ssc$fixef.K == "full" || !any_nested_var){
        K = object$nparams
        if(ssc$fixef.force_exact && n_fe >= 2 && n_fe_ok >= 1){
            fe = fixef(object, notes = FALSE)
            K = K + (n_fe_ok - 1) - sum(attr(fe, "references"))
        }
    } else {
        # nested
        # we delay the adjustment
        K = object$nparams
    }

    #
    # NESTING (== pain in the neck)
    #

    if(ssc$fixef.K == "nested" && n_fe_ok > 0 && any_nested_var){
        # OK, let's go checking....
        # We always try to minimize computation.
        # So we maximize deduction and apply computation only in last resort.

        nested_vcov_var_names = intersect(names(nested_vars[nested_vars]), names(vcov_vars))
        nested_var_names = var_names_all[nested_vcov_var_names]

        is_nested = rep(FALSE, n_fe)
        for(i in 1:n_fe){
            # We check for each FE
            fe_name = object$fixef_vars[i]

            if(fe_name %in% nested_var_names){
                is_nested[i] = TRUE
                next

            } else if(grepl("^", fe_name, fixed = TRUE)){
                # nesting of the FE in the parent term
                fe_name_split = strsplit(fe_name, "^", fixed = TRUE)[[1]]
                if(any(fe_name_split %in% nested_var_names)){
                    is_nested[i] = TRUE
                    next
                }
            }
        }


        # We check the remaining FEs, only if necessary
        if(!all(is_nested) && !all(nested_var_names %in% object$fixef_vars)){
            # Note that if all(nested_var_names %in% object$fixef_vars) == TRUE
            # Then all the variables used to compute the VCOV are part of the FEs
            # So the nesting would have been correctly spotted right before.
            # Only caveat: if the user has a^b in FE and includes the parent terms (a and b), we may forget some nested FEs
            # But then that's the user problem bc it's an erroneous specification

            vcov_vars_nesting = vcov_vars[nested_vcov_var_names]
            # We put the non integer to integer (needed)
            id_to_int = which(sapply(vcov_select$vars[nested_vcov_var_names], function(x) !isTRUE(x$to_int)))

            for(i in id_to_int){
                vcov_vars_nesting[[i]] = quickUnclassFactor(vcov_vars_nesting[[i]])
            }

            id2check = which(is_nested == FALSE)
            info = cpp_check_nested(object$fixef_id[id2check], vcov_vars_nesting, object$fixef_sizes[id2check], n = n)
            is_nested[id2check] = info == 1
        }


        if(sum(is_nested) == n_fe){
            # All FEs are removed, we add 1 for the intercept
            K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1)) + 1
        } else {
            if(ssc$fixef.force_exact && n_fe >= 2){
                fe = fixef(object, notes = FALSE)
                nb_ref = attr(fe, "references")

                # Slopes are a pain in the neck!!!
                if(sum(is_nested) > 1){
                    id_nested = intersect(names(nb_ref), names(object$fixef_id)[is_nested])
                    nb_ref[id_nested] = object$fixef_sizes[id_nested]
                }

                total_refs = sum(nb_ref)

                K = K - total_refs
            } else {
                K = K - (sum(fixef_sizes_ok[is_nested]) - sum(fixef_sizes_ok[is_nested] > 0))
            }
        }

        # below for consistency => should *not* be triggered
        K = max(K, length(object$coefficients) + 1)
    }

    # Small sample adjustment
    ss_adj = attr(vcov_noAdj, "ss_adj")
    if(!is.null(ss_adj)){
        if(is.function(ss_adj)){
            ss_adj = ss_adj(n = n, K = K)
        }
        attr(vcov_noAdj, "ss_adj") = NULL
    } else {
        ss_adj = ifelse(ssc$adj, (n - 1) / (n - K), 1)
    }

    vcov_mat = vcov_noAdj * ss_adj

    ####
    #### ... vcov attributes ####
    ####


    if(any(diag(vcov_mat) < 0)){
        # We 'fix' it
        all_attr = attributes(vcov_mat)
        vcov_mat = mat_posdef_fix(vcov_mat)
        attributes(vcov_mat) = all_attr
        message("Variance contained negative values in the diagonal and was 'fixed' (a la Cameron, Gelbach & Miller 2011).")
    }

    if(is_attr){

        min_cluster_size = attr(vcov_mat, "min_cluster_size")
        if(is.numeric(ssc$t.df)){
            attr(vcov_mat, "df.t") = ssc$t.df

        } else if(ssc$t.df == "min" && !is.null(min_cluster_size)){
            attr(vcov_mat, "df.t") = max(min_cluster_size - 1, 1)

        } else {
            attr(vcov_mat, "df.t") = max(n - K, 1)
        }

        if(is.null(attr(vcov_mat, "type", exact = TRUE))){
            type_info = attr(vcov_mat, "type_info")
            attr(vcov_mat, "type") = paste0(vcov_select$vcov_label, type_info)
            attr(vcov_mat, "type_info") = NULL
        }

        attr(vcov_mat, "ssc") = ssc
        attr(vcov_mat, "dof.K") = K
    } else {
        # We clean the attributes
        all_attr = names(attributes(vcov_mat))
        for(v in setdiff(all_attr, c("dim", "dimnames"))){
            attr(vcov_mat, v) = NULL
        }
    }

    vcov_mat
}




####
#### Small Sample Correction ####
####





#' Governs the small sample correction in \code{fixest} VCOVs
#'
#' Provides how the small sample correction should be calculated in \code{\link[fixest]{vcov.fixest}}/\code{\link[fixest]{summary.fixest}}.
#'
#' @param adj Logical scalar, defaults to \code{TRUE}. Whether to apply a small sample adjustment of the form \code{(n - 1) / (n - K)}, with \code{K} the number of estimated parameters. If \code{FALSE}, then no adjustment is made.
#' @param fixef.K Character scalar equal to \code{"nested"} (default), \code{"none"} or \code{"full"}. In the small sample adjustment, how to account for the fixed-effects parameters. If \code{"none"}, the fixed-effects parameters are discarded, meaning the number of parameters (\code{K}) is only equal to the number of variables. If \code{"full"}, then the number of parameters is equal to the number of variables plus the number of fixed-effects. Finally, if \code{"nested"}, then the number of parameters is equal to the number of variables plus the number of fixed-effects that *are not* nested in the clusters used to cluster the standard-errors.
#' @param fixef.force_exact Logical, default is \code{FALSE}. If there are 2 or more fixed-effects, these fixed-effects they can be irregular, meaning they can provide the same information. If so, the "real" number of parameters should be lower than the total number of fixed-effects. If \code{fixef.force_exact = TRUE}, then \code{\link[fixest]{fixef.fixest}} is first run to determine the exact number of parameters among the fixed-effects. Mostly, panels of the type individual-firm require \code{fixef.force_exact = TRUE} (but it adds computational costs).
#' @param cluster.adj Logical scalar, default is \code{TRUE}. How to make the small sample correction when clustering the standard-errors? If \code{TRUE} a \code{G/(G-1)} correction is performed with \code{G} the number of cluster values.
#' @param cluster.df Either "conventional" or "min" (default). Only relevant when the variance-covariance matrix is two-way clustered (or higher). It governs how the small sample adjustment for the clusters is to be performed. [Sorry for the jargon that follows.] By default a unique adjustment is made, of the form G_min/(G_min-1) with G_min the smallest G_i. If \code{cluster.df="conventional"} then the i-th "sandwich" matrix is adjusted with G_i/(G_i-1) with G_i the number of unique clusters.
#' @param t.df Either "conventional", "min" (default) or an integer scalar. Only relevant when the variance-covariance matrix is clustered. It governs how the p-values should be computed. By default, the degrees of freedom of the Student t distribution is equal to the minimum size of the clusters with which the VCOV has been clustered. If \code{t.df="conventional"}, then the degrees of freedom of the Student t distribution is equal to the number of observations minus the number of estimated variables. You can also pass a number to manually specify the DoF of the t-distribution.
#'
#' @details
#'
#' The following vignette: \href{https://lrberge.github.io/fixest/articles/standard_errors.html}{On standard-errors}, describes in details how the standard-errors are computed in \code{fixest} and how you can replicate standard-errors from other software.
#'
#' @return
#' It returns a \code{ssc.type} object.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{vcov.fixest}}
#'
#' @examples
#'
#' #
#' # Equivalence with lm/glm standard-errors
#' #
#'
#' # LM
#' # In the absence of fixed-effects,
#' # by default, the standard-errors are computed in the same way
#'
#' res = feols(Petal.Length ~ Petal.Width + Species, iris)
#' res_lm = lm(Petal.Length ~ Petal.Width + Species, iris)
#' vcov(res) / vcov(res_lm)
#'
#' # GLM
#' # By default, there is no small sample adjustment in glm, as opposed to feglm.
#' # To get the same SEs, we need to use ssc(adj = FALSE)
#'
#' res_pois = fepois(round(Petal.Length) ~ Petal.Width + Species, iris)
#' res_glm = glm(round(Petal.Length) ~ Petal.Width + Species, iris, family = poisson())
#' vcov(res_pois, ssc = ssc(adj = FALSE)) / vcov(res_glm)
#'
#' # Same example with the Gamma
#' res_gamma = feglm(round(Petal.Length) ~ Petal.Width + Species, iris, family = Gamma())
#' res_glm_gamma = glm(round(Petal.Length) ~ Petal.Width + Species, iris, family = Gamma())
#' vcov(res_gamma, ssc = ssc(adj = FALSE)) / vcov(res_glm_gamma)
#'
#' #
#' # Fixed-effects corrections
#' #
#'
#' # We create "irregular" FEs
#' base = data.frame(x = rnorm(10))
#' base$y = base$x + rnorm(10)
#' base$fe1 = rep(1:3, c(4, 3, 3))
#' base$fe2 = rep(1:5, each = 2)
#'
#' est = feols(y ~ x | fe1 + fe2, base)
#'
#' # fe1: 3 FEs
#' # fe2: 5 FEs
#'
#' #
#' # Clustered standard-errors: by fe1
#' #
#'
#' # Default: fixef.K = "nested"
#' #  => adjustment K = 1 + 5 (i.e. x + fe2)
#' summary(est)
#' attributes(vcov(est, attr = TRUE))[c("ssc", "dof.K")]
#'
#'
#' # fixef.K = FALSE
#' #  => adjustment K = 1 (i.e. only x)
#' summary(est, ssc = ssc(fixef.K = "none"))
#' attr(vcov(est, ssc = ssc(fixef.K = "none"), attr = TRUE), "dof.K")
#'
#'
#' # fixef.K = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 1 (i.e. x + fe1 + fe2 - 1 restriction)
#' summary(est, ssc = ssc(fixef.K = "full"))
#' attr(vcov(est, ssc = ssc(fixef.K = "full"), attr = TRUE), "dof.K")
#'
#'
#' # fixef.K = TRUE & fixef.force_exact = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 2 (i.e. x + fe1 + fe2 - 2 restrictions)
#' summary(est, ssc = ssc(fixef.K = "full", fixef.force_exact = TRUE))
#' attr(vcov(est, ssc = ssc(fixef.K = "full", fixef.force_exact = TRUE), attr = TRUE), "dof.K")
#'
#' # There are two restrictions:
#' attr(fixef(est), "references")
#'
#' #
#' # To permanently set the default ssc:
#' #
#'
#' # eg no small sample adjustment:
#' setFixest_ssc(ssc(adj = FALSE))
#'
#' # Factory default
#' setFixest_ssc()
#'
ssc = function(adj = TRUE, fixef.K = "nested", cluster.adj = TRUE, cluster.df = "min",
               t.df = "min", fixef.force_exact = FALSE){

    check_arg_plus(adj, "loose logical scalar conv")
    check_arg_plus(fixef.K, "match(none, full, nested)")
    check_arg_plus(cluster.df, "match(conventional, min)")
    check_arg_plus(t.df, "match(conventional, min) | integer scalar GT{0}")
    check_arg(fixef.force_exact, cluster.adj, "logical scalar")

    res = list(adj = adj, fixef.K = fixef.K, cluster.adj = cluster.adj, cluster.df = cluster.df,
               t.df = t.df, fixef.force_exact = fixef.force_exact)
    class(res) = "ssc.type"

    res
}

#' @describeIn ssc This function is deprecated and will be removed at some point (in 6 months from August 2021). Exactly the same as \code{ssc}.
dof = function(adj = TRUE, fixef.K = "nested", cluster.adj = TRUE, cluster.df = "min",
               t.df = "min", fixef.force_exact = FALSE){

    if(is.null(getOption("fixest_warn_dof"))){
        warning("The function 'dof' is deprecated. Please use function 'ssc' instead.")
        options(fixest_warn_dof = TRUE)
    }

    ssc(adj = adj, fixef.K = fixef.K, cluster.adj = cluster.adj, cluster.df = cluster.df,
        t.df = t.df, fixef.force_exact = fixef.force_exact)
}


####
#### User-level ####
####


#' Clustered VCOV
#'
#' Computes the clustered VCOV of \code{fixest} objects.
#'
#' @param x A \code{fixest} object.
#' @param cluster Either i) a character vector giving the names of the variables onto which to cluster, or ii) a formula giving those names, or iii) a vector/list/data.frame giving the hard values of the clusters. Note that in cases i) and ii) the variables are fetched directly in the data set used for the estimation.
#' @param ssc An object returned by the function \code{\link[fixest]{ssc}}. It specifies how to perform the small sample correction.
#'
#' @return
#' If the first argument is a \code{fixest} object, then a VCOV is returned (i.e. a symmetric matrix).
#'
#' If the first argument is not a \code{fixest} object, then a) implicitly the arguments are shifted to the left (i.e. \code{vcov_cluster(~var1 + var2)} is equivalent to \code{vcov_cluster(cluster = ~var1 + var2)}) and b) a VCOV-\emph{request} is returned and NOT a VCOV. That VCOV-request can then be used in the argument \code{vcov} of various \code{fixest} functions (e.g. \code{\link[fixest]{vcov.fixest}} or even in the estimation calls).
#'
#' @author
#' Laurent Berge
#'
#' @references
#' Cameron AC, Gelbach JB, Miller DL (2011). "Robust Inference with Multiway Clustering." \emph{Journal of Business & Economic Statistics}, 29(2), 238-249. doi:10.1198/jbes.2010.07136.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#' base$clu = rep(1:5, 30)
#'
#' est = feols(y ~ x1, base)
#'
#' # VCOV: using a formula giving the name of the clusters
#' vcov_cluster(est, ~species + clu)
#'
#' # works as well with a character vector
#' vcov_cluster(est, c("species", "clu"))
#'
#' # you can also combine the two with '^'
#' vcov_cluster(est, ~species^clu)
#'
#' #
#' # Using VCOV requests
#' #
#'
#' # per se: pretty useless...
#' vcov_cluster(~species)
#'
#' # ...but VCOV-requests can be used at estimation time:
#' # it may be more explicit than...
#' feols(y ~ x1, base, vcov = vcov_cluster("species"))
#'
#' # ...the equivalent, built-in way:
#' feols(y ~ x1, base, vcov = ~species)
#'
#' # The argument vcov does not accept hard values,
#' # so you can feed them with a VCOV-request:
#' feols(y ~ x1, base, vcov = vcov_cluster(rep(1:5, 30)))
#'
#'
vcov_cluster = function(x, cluster = NULL, ssc = NULL){
    # User-level function to compute clustered SEs
    # typically we only do checking and reshaping here

    # slide_args allows the implicit allocation of arguments
    # it makes semi-global changes => the values of the args here are modified
    slide_args(x, cluster = cluster, ssc = ssc)
    IS_REQUEST = is.null(x)

    check_value(ssc, "NULL class(ssc.type)", .message = "The argument 'ssc' must be an object created by the function ssc().")


    # We create the request

    use_request = FALSE
    vcov_vars = var_names_all = NULL
    if(is.null(cluster)){
        vcov = "cluster"
    } else if(inherits(cluster, "formula")){
        vcov = cluster
    } else if(is.character(cluster) && length(cluster) <= 4){
        vcov = .xpd(lhs = "cluster", rhs = cluster)
    } else {
        use_request = TRUE
        vcov = "cluster"
        vcov_vars = cluster
    }

    # The main task is to check the validity of the vcov_vars in input
    # and give proper names
    if(!is.null(vcov_vars)){

        sc = sys.call()
        if("cluster" %in% names(sc)){
            cl_name = fetch_arg_deparse("cluster")
        } else {
            cl_name = deparse_long(sc[[2]])
        }

        if(is.atomic(vcov_vars)){
            vcov_vars = list(cl1 = vcov_vars)
            var_names_all["cl1"] = gsub("^.+\\$", "", cl_name)

        } else if(!is.list(vcov_vars) || length(vcov_vars) > 4){
            # We send an error + a reminder of the accepted types
            check_value(vcov_vars, "os formula | character vector | list len(,4)", .arg_name = "cluster")

        } else {
            vcov_var_names = paste0("cl", 1:length(vcov_vars))

            if(!is.null(names(vcov_vars))){
                var_names_all[vcov_var_names] = names(vcov_vars)

            } else {

                len = length(vcov_var_names)
                new_name = ifunit(len, cl_name, paste0(cl_name, "..", 1:len))

                var_names_all[vcov_var_names] = new_name
            }

            # We ensure it's a plain list
            vcov_vars = unclass(vcov_vars)
            names(vcov_vars) = vcov_var_names

            # Some further checking
            if(!is.atomic(vcov_vars[[1]])){
                stop("If the argument 'cluster' it to be a list, it must be made of atomic vectors representing the vectors on which to cluster. Currently its first element is of class '", class(vcov_vars[[1]])[1], "'.")
            }

        }

    }

    if(use_request || !is.null(ssc)){
        vcov_request = list(vcov = vcov, vcov_vars = vcov_vars,
                            var_names_all = var_names_all, ssc = ssc)
        class(vcov_request) = "fixest_vcov_request"

    } else {
        # Everything can fit into a vcov formula
        vcov_request = vcov
    }

    if(IS_REQUEST){
        res = vcov_request
    } else {
        res = vcov(x, vcov = vcov_request)
    }

    res
}



#' HAC VCOVs
#'
#' Set of functions to compute the VCOVs robust to different forms correlation in panel or time series settings.
#'
#' There are currently three VCOV types: Newey-West applied to time series, Newey-West applied to a panel setting (when the argument 'unit' is not missing), and Driscoll-Kraay.
#'
#' The functions on this page without the prefix "vcov_" do not compute VCOVs directly but are meant to be used in the argument \code{vcov} of \code{fixest} functions (e.g. in \code{\link[fixest]{vcov.fixest}} or even in the estimation calls).
#'
#' Note that for Driscoll-Kraay VCOVs, to ensure its properties the number of periods should be long enough (a minimum of 20 periods or so).
#'
#' @inheritParams vcov_cluster
#'
#' @param unit A character scalar or a one sided formula giving the name of the variable representing the units of the panel.
#' @param time A character scalar or a one sided formula giving the name of the variable representing the time.
#' @param lag An integer scalar, default is \code{NULL}. If \code{NULL}, then the default lag is equal to \code{n_t^0.25} with \code{n_t} the number of time periods (as of Newey and West 1987) for panel Newey-West and Driscoll-Kraay. The default for the time series Newey-West is computed via \code{\link[sandwich:NeweyWest]{bwNeweyWest}} which implements the Newey and West 1994 method.
#'
#' @section Lag selection:
#'
#' The default lag selection depends on whether the VCOV applies to a panel or a time series.
#'
#' For panels, i.e. panel Newey-West or Driscoll-Kraay VCOV, the default lag is \code{n_t^0.25} with \code{n_t} the number of time periods. This is based on Newey and West 1987.
#'
#' For time series Newey-West, the default lag is found thanks to the \code{\link[sandwich:NeweyWest]{bwNeweyWest}} function from the \code{sandwich} package. It is based on Newey and West 1994.
#'
#'
#' @references
#' Newey WK, West KD (1987). "A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix." \emph{Econometrica}, 55(3), 703-708. doi:10.2307/1913610.
#'
#' Driscoll JC, Kraay AC (1998). "Consistent Covariance Matrix Estimation with Spatially Dependent Panel Data." \emph{The Review of Economics and Statistics}, 80(4), 549-560. doi:10.1162/003465398557825.
#'
#' Millo G (2017). "Robust Standard Error Estimators for Panel Models: A Unifying Approach" \emph{Journal of Statistical Software}, 82(3). doi:10.18637/jss.v082.i03.
#'
#' @return
#' If the first argument is a \code{fixest} object, then a VCOV is returned (i.e. a symmetric matrix).
#'
#' If the first argument is not a \code{fixest} object, then a) implicitly the arguments are shifted to the left (i.e. \code{vcov_DK(~year)} is equivalent to \code{vcov_DK(time = ~year)}) and b) a VCOV-\emph{request} is returned and NOT a VCOV. That VCOV-request can then be used in the argument \code{vcov} of various \code{fixest} functions (e.g. \code{\link[fixest]{vcov.fixest}} or even in the estimation calls).
#'
#' @examples
#'
#' data(base_did)
#'
#' #
#' # During the estimation
#' #
#'
#' # Panel Newey-West, lag = 2
#' feols(y ~ x1, base_did, NW(2) ~ id + period)
#'
#' # Driscoll-Kraay
#' feols(y ~ x1, base_did, DK ~ period)
#'
#' # If the estimation is made with a panel.id, the dimensions are
#' # automatically deduced:
#' est = feols(y ~ x1, base_did, "NW", panel.id = ~id + period)
#' est
#'
#' #
#' # Post estimation
#' #
#'
#' # If missing, the unit and time are automatically deduced from
#' # the panel.id used in the estimation
#' vcov_NW(est, lag = 2)
#'
#'
#' @name vcov_hac
NULL

#' @rdname vcov_hac
vcov_DK = function(x, time = NULL, lag = NULL, ssc = NULL){
    # unit and time MUST be variables of the data set
    # otherwise: too error prone + extremely complex to set up + very edge case, so not worth it
    #

    # slide_args allows the implicit allocation of arguments
    # it makes semi-global changes => the values of the args here are modified
    slide_args(x, time = time, lag = lag, ssc = ssc)
    IS_REQUEST = is.null(x)

    check_value(lag, "NULL integer scalar GE{1}")
    check_value(ssc, "NULL class(ssc.type)", .message = "The argument 'ssc' must be an object created by the function ssc().")

    # time
    check_value(time, "NULL character scalar | os formula")

    if(inherits(time, "formula")){
        time_txt = fml2varnames(time)
        check_value(time_txt, "character scalar", .message = "The argument 'time' must be composed of only one variable.")
    }

    # recreating the call
    vcov = .xpd(lhs = "dk", rhs = time)

    if(!is.null(lag) || !is.null(ssc)){
        extra_args = list(lag = lag)
        vcov_request = list(vcov = vcov, ssc = ssc, extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
    } else {
        # Everything can fit into a vcov formula
        vcov_request = vcov
    }

    if(IS_REQUEST){
        res = vcov_request
    } else {
        res = vcov(x, vcov = vcov_request)
    }

    res

}

#' @rdname vcov_hac
vcov_NW = function(x, unit = NULL, time = NULL, lag = NULL, ssc = NULL){
    # unit and time MUST be variables of the data set
    # otherwise: too error prone + extremely complex to set up + very edge case, so not worth it
    #

    # slide_args allows the implicit allocation of arguments
    # it makes semi-global changes => the values of the args here are modified
    slide_args(x, unit = unit, time = time, lag = lag, ssc = ssc)
    IS_REQUEST = is.null(x)

    check_value(lag, "NULL integer scalar GE{1}")
    check_value(ssc, "NULL class(ssc.type)", .message = "The argument 'ssc' must be an object created by the function ssc().")

    # unit
    check_value(unit, "NULL character scalar | os formula")

    if(inherits(unit, "formula")){
        unit_txt = fml2varnames(unit)
        check_value(unit_txt, "character scalar", .message = "The argument 'unit' must be composed of only one variable.")
    }

    # time
    check_value(time, "NULL character scalar | os formula")

    if(inherits(time, "formula")){
        time_txt = fml2varnames(time)
        check_value(time_txt, "character scalar", .message = "The argument 'time' must be composed of only one variable.")
    }

    # recreating the call
    vcov = .xpd(lhs = "nw", rhs = c(unit, time))

    if(!is.null(lag) || !is.null(ssc)){
        extra_args = list(lag = lag)
        vcov_request = list(vcov = vcov, ssc = ssc, extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
    } else {
        # Everything can fit into a vcov formula
        vcov_request = vcov
    }

    if(IS_REQUEST){
        res = vcov_request
    } else {
        res = vcov(x, vcov = vcov_request)
    }

    res

}


#' Conley VCOV
#'
#' Compute VCOVs robust to spatial correlation, a la Conley (1999).
#'
#' This function computes VCOVs that are robust to spatial correlations by assuming a correlation between the units that are at a geographic distance lower than a given cutoff.
#'
#' The kernel is uniform.
#'
#' If the cutoff is not provided, an estimation of it is given. This cutoff ensures that a minimum of units lie within it and is robust to sub-sampling. This automatic cutoff is only here for convenience, the most appropriate cutoff shall depend on the application and shall be provided by the user.
#'
#' The function \code{conley} does not compute VCOVs directly but is meant to be used in the argument \code{vcov} of \code{fixest} functions (e.g. in \code{\link[fixest]{vcov.fixest}} or even in the estimation calls).
#'
#' @inheritParams vcov_cluster
#'
#' @param lat A character scalar or a one sided formula giving the name of the variable representing the latitude. The latitude must lie in [-90, 90], [0, 180] or [-180, 0].
#' @param lon A character scalar or a one sided formula giving the name of the variable representing the longitude. The longitude must be in [-180, 180], [0, 360] or [-360, 0].
#' @param cutoff The distance cutoff, in km. You can express the cutoff in miles by writing the number in character form and adding "mi" as a suffix: cutoff = "100mi" would be 100 miles. If missing, a rule of thumb is used to deduce the cutoff.
#' @param pixel A positive numeric scalar, default is 0. If a positive number, the coordinates of each observation are pooled into \code{pixel} x \code{pixel} km squares. This lowers the precision but can (depending on the cases) greatly improve computational speed at a low precision cost. Note that if the \code{cutoff} was expressed in miles, then \code{pixel} will also be in miles.
#' @param distance How to compute the distance between points. It can be equal to "triangular" (default) or "spherical". The latter case corresponds to the great circle distance and is more precise than triangular but is a bit more intensive computationally.
#'
#' @return
#' If the first argument is a \code{fixest} object, then a VCOV is returned (i.e. a symmetric matrix).
#'
#' If the first argument is not a \code{fixest} object, then a) implicitly the arguments are shifted to the left (i.e. \code{vcov_conley("lat", "long")} is equivalent to \code{vcov_conley(lat = "lat", lon = "long")}) and b) a VCOV-\emph{request} is returned and NOT a VCOV. That VCOV-request can then be used in the argument \code{vcov} of various \code{fixest} functions (e.g. \code{\link[fixest]{vcov.fixest}} or even in the estimation calls).
#'
#' @references
#' Conley TG (1999). "GMM Estimation with Cross Sectional Dependence", \emph{Journal of Econometrics}, 92, 1-45.
#'
#' @examples
#'
#' data(quakes)
#'
#' # We use conley() in the vcov argument of the estimation
#' feols(depth ~ mag, quakes, conley(100))
#'
#' # Post estimation
#' est = feols(depth ~ mag, quakes)
#' vcov_conley(est, cutoff = 100)
#'
#'
#'
vcov_conley = function(x, lat = NULL, lon = NULL, cutoff = NULL, pixel = 0,
                       distance = "triangular", ssc = NULL){

    # slide_args allows the implicit allocation of arguments
    # it makes semi-global changes => the values of the args here are modified
    slide_args(x, lat = lat, lon = lon, cutoff = cutoff, pixel = pixel, distance = distance, ssc = ssc)
    IS_REQUEST = is.null(x)

    check_value(ssc, "NULL class(ssc.type)", .message = "The argument 'ssc' must be an object created by the function ssc().")

    # lat
    check_value(lat, "NULL character scalar | os formula")

    if(inherits(lat, "formula")){
        lat_txt = fml2varnames(lat)
        check_value(lat_txt, "character scalar", .message = "The argument 'lat' must be composed of only one variable.")
    }


    # lon
    check_value(lon, "NULL character scalar | os formula")

    if(inherits(lon, "formula")){
        lon_txt = fml2varnames(lon)
        check_value(lon_txt, "character scalar", .message = "The argument 'lon' must be composed of only one variable.")
    }

    # cutoff
    cutoff = check_set_cutoff(cutoff)

    # pixel
    check_value(pixel, "numeric scalar GE{0}",
                .prefix = "In vcov.fixest, Conley VCOV cannot be computed: the argument 'pixel'")

    # distance
    check_value_plus(distance, "match(triangular, spherical)")

    # recreating the call
    vcov = .xpd(lhs = "conley", rhs = c(lat, lon))

    extra_args = list(cutoff = cutoff, pixel = pixel, distance = distance)
    vcov_request = list(vcov = vcov, ssc = ssc, extra_args = extra_args)
    class(vcov_request) = "fixest_vcov_request"

    if(IS_REQUEST){
        res = vcov_request
    } else {
        res = vcov(x, vcov = vcov_request)
    }

    res
}


####
#### SETUP ####
####



vcov_setup = function(){

    #
    # iid
    #

    vcov_iid_setup = list(name = c("iid", "normal", "standard"), fun_name = "vcov_iid_internal", vcov_label = "IID")

    #
    # Heteroskedasticity robust
    #

    vcov_hetero_setup = list(name = c("hetero", "white", "hc1"), fun_name = "vcov_hetero_internal", vcov_label = "Heteroskedasticity-robust")

    #
    # cluster(s)
    #

    vcov_clust_setup = list(name = c("cluster", ""), fun_name = "vcov_cluster_internal", vcov_label = "Clustered")
    vcov_clust_setup$vars = list(cl1 = list(guess_from = list(panel.id = 1, fixef = 1), label = "clusters", to_int = TRUE, rm_nested = TRUE),
                                 cl2 = list(optional = TRUE, label = "second cluster", to_int = TRUE, rm_nested = TRUE),
                                 cl3 = list(optional = TRUE, label = "third cluster", to_int = TRUE, rm_nested = TRUE),
                                 cl4 = list(optional = TRUE, label = "fourth cluster", to_int = TRUE, rm_nested = TRUE))
    vcov_clust_setup$patterns = c("", "cl1", "cl1 + cl2", "cl1 + cl2 + cl3",
                                  "c1 + cl2 + cl3 + cl4")

    # Other keywords => direct two-/three-/four-way clustering
    cl1 = list(guess_from = list(fixef = 1), label = "first cluster", to_int = TRUE, rm_nested = TRUE)
    cl2 = list(guess_from = list(fixef = 2), label = "second cluster", to_int = TRUE, rm_nested = TRUE)
    cl3 = list(guess_from = list(fixef = 3), label = "third cluster", to_int = TRUE, rm_nested = TRUE)
    cl4 = list(guess_from = list(fixef = 4), label = "fourth cluster", to_int = TRUE, rm_nested = TRUE)

    vcov_twoway_setup = list(name = "twoway", fun_name = "vcov_cluster_internal", vcov_label = "Clustered")
    vcov_twoway_setup$vars = list(cl1 = cl1, cl2 = cl2)
    vcov_twoway_setup$patterns = c("", "cl1 + cl2")

    vcov_threeway_setup = list(name = "threeway", fun_name = "vcov_cluster_internal", vcov_label = "Clustered")
    vcov_threeway_setup$vars = list(cl1 = cl1, cl2 = cl2, cl3 = cl3)
    vcov_threeway_setup$patterns = c("", "cl1 + cl2 + cl3")

    vcov_fourway_setup = list(name = "fourway", fun_name = "vcov_cluster_internal", vcov_label = "Clustered")
    vcov_fourway_setup$vars = list(cl1 = cl1, cl2 = cl2, cl3 = cl3, cl4 = cl4)
    vcov_fourway_setup$patterns = c("", "cl1 + cl2 + cl3 + cl4")


    #
    # newey_west
    #

    # The variables
    unit = list(guess_from = list(panel.id = 1), label = "panel ID", to_int = TRUE, rm_nested = TRUE,
                optional = TRUE)
    time = list(guess_from = list(panel.id = 2), label = "time", rm_nested = TRUE)

    vcov_newey_west_setup = list(name = c("NW", "newey_west"),
                                 fun_name = "vcov_newey_west_internal",
                                 vcov_label = "Newey-West")

    vcov_newey_west_setup$vars = list(unit = unit, time = time)
    vcov_newey_west_setup$arg_main = "lag"
    vcov_newey_west_setup$rdname = "vcov_hac"
    vcov_newey_west_setup$patterns = c("", "time", "unit + time")


    #
    # driscoll_kraay
    #

    vcov_driscoll_kraay_setup = list(name = c("DK", "driscoll_kraay"),
                                     fun_name = "vcov_driscoll_kraay_internal",
                                     vcov_label = "Driscoll-Kraay")

    vcov_driscoll_kraay_setup$vars = list(time = time)
    vcov_driscoll_kraay_setup$arg_main = "lag"
    vcov_driscoll_kraay_setup$rdname = "vcov_hac"
    vcov_driscoll_kraay_setup$patterns = c("", "time")


    #
    # conley
    #

    # The variables
    lat = list(guess_from = list(regex = c("^lat(itude)?$", "^lat_.+")),
               label = "latitude",
               expected_type = "numeric vector",
               rm_nested = TRUE)

    lng = list(guess_from = list(regex = c("^(lon|lng|long|longitude)$", "^(lng|lon|long)_.+")),
               label = "longitude",
               expected_type = "numeric vector",
               rm_nested = TRUE)

    vcov_conley_setup = list(name = "conley", fun_name = "vcov_conley_internal", vcov_label = "Conley")
    vcov_conley_setup$vars = list(lat = lat, lng = lng)
    vcov_conley_setup$arg_main = c("cutoff", "pixel", "distance")
    vcov_conley_setup$rdname = "vcov_conley"
    vcov_conley_setup$patterns = c("", "lat + lng")


    #
    # conley hac
    #

    vcov_conley_hac_setup = list(name = c("conley_hac", "hac_conley"), fun_name = "vcov_conley_hac_internal", vcov_label = "Conley-HAC")
    # The variables (already defined earlier)
    unit_conleyHAC = unit
    unit_conleyHAC$optional = NULL
    vcov_conley_hac_setup$vars = list(lat = lat, lng = lng, unit = unit_conleyHAC, time = time)
    vcov_conley_hac_setup$arg_main = c("cutoff", "pixel", "distance", "lag")
    vcov_conley_hac_setup$patterns = c("", "lat + lng", "time", "lat + lng + unit + time")

    #
    # Saving all the vcov possibilities
    #

    all_vcov = list(vcov_iid_setup,
                    vcov_hetero_setup,
                    vcov_clust_setup,
                    vcov_twoway_setup,
                    vcov_threeway_setup,
                    vcov_fourway_setup,
                    vcov_newey_west_setup,
                    vcov_driscoll_kraay_setup,
                    vcov_conley_setup)
    # vcov_conley_hac_setup)

    options(fixest_vcov_builtin = all_vcov)

}



####
#### Internal ####
####




vcovClust = function (cluster, myBread, scores, adj = FALSE, do.unclass = TRUE, sandwich = TRUE, nthreads = 1){
    # Internal function: no need for controls, they come beforehand
    # - cluster: the vector of dummies
    # - myBread: original vcov
    # - scores
    # Note: if length(unique(cluster)) == n (i.e. White correction), then the adj are such that vcovClust is equivalent to vcovHC(res, type="HC1")
    # Source: http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
    #         Cameron & Miller -- A Practitioner's Guide to Cluster-Robust Inference

    n = NROW(scores)

    # Control for cluster type
    if(do.unclass){
        cluster = quickUnclassFactor(cluster)
    }

    Q = max(cluster)
    RightScores = cpp_tapply_sum(Q, scores, cluster)

    # Finite sample correction:
    if(adj){
        adj_value = Q / (Q - 1)
    } else {
        adj_value = 1
    }

    if(!sandwich){
        res = cpppar_crossprod(RightScores, 1, nthreads) * adj_value
        return(res)
    }

    xy = cpppar_matprod(RightScores, myBread, nthreads)
    res = cpppar_crossprod(xy, 1, nthreads) * adj_value
    res
}


vcov_iid_internal = function(bread, ...){
    bread
}

vcov_hetero_internal = function(bread, scores, sandwich, ssc, nthreads, ...){

    if(!sandwich){
        vcov_noAdj = cpppar_crossprod(scores, 1, nthreads)
    } else {
        # we make a n/(n-1) adjustment to match vcovHC(type = "HC1")

        n = nrow(scores)
        adj = ifelse(ssc$cluster.adj, n / (n-1), 1)

        vcov_noAdj = cpppar_crossprod(cpppar_matprod(scores, bread, nthreads), 1, nthreads) * adj
    }

    vcov_noAdj
}


vcov_cluster_internal = function(bread, scores, vars, ssc, sandwich, nthreads, var_names_all, ...){

    # aliasing to add (a bit of) clarity
    cluster = vars
    nway = length(cluster)

    # initialization
    vcov = bread * 0
    meat = 0

    for(i in 1:nway){

        myComb = combn(nway, i)

        power = floor(1 + log10(sapply(cluster, max)))

        for(j in 1:ncol(myComb)){

            if(i == 1){
                index = cluster[[myComb[, j]]]

            } else if(i > 1){

                vars = myComb[, j]

                if(sum(power[vars]) > 14){
                    my_clusters = cluster[vars]

                    order_index = do.call(order, my_clusters)
                    index = cpp_combine_clusters(my_clusters, order_index)
                } else {
                    # quicker, but limited by the precision of integers
                    index = cluster[[vars[1]]]
                    for(k in 2:length(vars)){
                        index = index + cluster[[vars[k]]]*10**sum(power[vars[1:(k-1)]])
                    }
                }

                index = quickUnclassFactor(index)

            }

            vcov = vcov + (-1)**(i+1) * vcovClust(index, bread, scores,
                                                  adj = ssc$cluster.adj && ssc$cluster.df == "conventional",
                                                  sandwich = sandwich, nthreads = nthreads)

        }
    }

    G_min = min(sapply(cluster, max))
    if(ssc$cluster.adj && ssc$cluster.df == "min"){
        vcov = vcov * G_min / (G_min - 1)
        attr(vcov, "G") = G_min
    }

    if(!sandwich){
        return(vcov)
    }

    # I know I duplicate the information, but they refer to two different things
    attr(vcov, "min_cluster_size") = G_min

    var_names_all = var_names_all[nchar(var_names_all) > 0]
    if(length(var_names_all) > 0){
        attr(vcov, "type_info") = paste0(" (", paste(var_names_all, collapse = " & "), ")")
    } else {
        attr(vcov, "type") = switch(nway,
                                    "1" = "Clustered",
                                    "2" = "Two-way",
                                    "3" = "Three-way",
                                    "4" = "Four-way")
    }

    vcov
}


vcov_newey_west_internal = function(bread, scores, vars, ssc, sandwich, nthreads, lag = NULL, ...){
    # Function that computes Newey-West VCOV

    # Setting up

    unit = vars$unit
    time = to_integer(vars$time, sorted = TRUE)

    n_time = max(time)

    # The bartlett weights
    set_weights = function(lag){
        res = seq(1, 0, by = -(1/(lag + 1)))
        # we must halve the first weight (we'll add the transpose internally)
        res[1] = 0.5
        res
    }

    is_panel = !is.null(unit)
    if(is_panel){

        # Lag: simple rule of thumb
        if(missnull(lag)){
            lag = floor(n_time^(1/4))
        }

        w = set_weights(lag)

        my_order = order(time, unit)
        time_ro = time[my_order]
        unit_ro = unit[my_order]

        dup = cpp_find_duplicates(unit_ro, time_ro)

        if(dup$n_dup > 0){
            stop("In vcov.fixest, Newey-West VCOV cannot be computed: there are (unit x time) duplicates. You have to sort that out first, eg by creating new units free of duplicates or dropping duplicates. Or you can use a Driscoll-Kraay VCOV for which duplicates does not matter.", call. = FALSE)
        }

        scores_ro = scores[my_order, , drop = FALSE]

        n_unit = max(unit_ro)

        meat = cpp_newey_west_panel(scores_ro, w, unit_ro, n_unit,
                                         time_ro, n_time, nthreads)

    } else {

        if(max(time) < length(time)){
            stop("In vcov.fixest, Newey-West VCOV cannot be computed: there are time duplicates. You may provide a panel identifier (if relevant) to sort that out. Or you can use a Driscoll-Kraay VCOV for which duplicates does not matter.", call. = FALSE)
        }

        my_order = order(time)
        time_ro = time[my_order]
        scores_ro = scores[my_order, , drop = FALSE]

        # Lag: the default relies on bwNeweyWest
        if(missnull(lag)){

            # Later => feed in directly the matrix when the new version of sandwich is ready

            # we drop the intercept (except if there's only the intercept)
            is_intercept = (colnames(scores_ro) == "(Intercept)") & (ncol(scores_ro) > 1)
            if(any(is_intercept)){
                x = structure(list(scores = scores_ro[, !is_intercept, drop = FALSE]), class = "fixest")
            } else {
                x = structure(list(scores = scores_ro), class = "fixest")
            }

            lag = sandwich::bwNeweyWest(x, weights = 1)
            lag = floor(lag)
        }

        w = set_weights(lag)

        meat = cpp_newey_west(scores_ro, w, nthreads)
    }

    if(sandwich){
        vcov = prepare_sandwich(bread, meat, nthreads)
    } else {
        vcov = meat
    }

    if(ssc$cluster.adj){
        vcov = vcov * n_time / (n_time - 1)
    }

    attr(vcov, "G") = n_time
    attr(vcov, "min_cluster_size") = n_time

    attr(vcov, "type_info") = paste0(" (L=", lag, ")")

    vcov
}


vcov_driscoll_kraay_internal = function(bread, scores, vars, ssc, sandwich, nthreads, lag = NULL, ...){
    # Function that computes Driscoll-Kraay VCOV

    # Setting up
    time = to_integer(vars$time, sorted = TRUE)

    n_time = max(time)

    # Lag: simple rule of thumb
    if(missnull(lag)){
        lag = floor(n_time^(1/4))
    }

    w = seq(1, 0, by = -(1/(lag + 1)))
    # we halve the first weight (since we add the transpose in the internal code)
    w[1] = 0.5

    my_order = order(time)
    time_ro = time[my_order]
    scores_ro = scores[my_order, , drop = FALSE]

    meat = cpp_driscoll_kraay(scores_ro, w, time_ro, n_time, nthreads)

    if(sandwich){
        vcov = prepare_sandwich(bread, meat, nthreads)
    } else {
        vcov = meat
    }

    if(ssc$cluster.adj){
        vcov = vcov * n_time / (n_time - 1)
    }

    attr(vcov, "G") = n_time
    attr(vcov, "min_cluster_size") = n_time

    attr(vcov, "type_info") = paste0(" (L=", lag, ")")

    vcov
}

vcov_conley_internal = function(bread, scores, vars, sandwich, nthreads,
                                cutoff = NULL, pixel = 0, distance = "triangular", ...){

    lon = vars$lng
    lat = vars$lat

    #
    # START :: checks

    # lon
    lon_range = range(lon)
    if(lon_range[1] < -360 || lon_range[2] > 360 || diff(lon_range) > 360){
        stop("In vcov.fixest, Conley VCOV cannot be computed: the longitude is outside the [-180, 180], [0, 360] or [-360, 0] range (current range is [", fsignif(lon_range[1]), ", ", fsignif(lon_range[2]), "]).", call. = FALSE)
    }

    # Later:
    # - I think that I really don't need to rescale but I should check first in the cpp code
    # that I don't use the -180/180 bounds
    if(lon_range[2] > 180){
        lon = lon - 180
    } else if(lon_range[1] < -180){
        lon = lon + 180
    }

    # lat
    lat_range = range(lat)
    if(lat_range[1] < -180 || lat_range[2] > 180 || diff(lat_range) > 180){
        stop("In vcov.fixest, Conley VCOV cannot be computed: the latitude is outside the [-90, 90], [0, 180] or [-180, 0] range (current range is [", fsignif(lat_range[1]), ", ", fsignif(lat_range[2]), "]).", call. = FALSE)
    }

    if(lat_range[2] > 90){
        # we normalize
        lat = lat - 90
    } else if(lat_range[1] < -90){
        lat = lat + 90
    }

    # cutoff
    cutoff = check_set_cutoff(cutoff)

    if(is.null(cutoff)){
        cutoff = cutoff_deduce(lat, lon)
    }

    # pixel
    check_value(pixel, "numeric scalar GE{0}",
                .prefix = "In vcov.fixest, Conley VCOV cannot be computed: the argument 'pixel'")

    # distance
    check_value_plus(distance, "match(triangular, spherical)")
    distance = switch(distance,
                      spherical = 1,
                      triangular = 2)

    #   END :: checks
    #

    metric = attr(cutoff, "metric")
    if(metric == "mi"){
        pixel = pixel * 1.60934
    }

    if(pixel > 0){
        pixel_lat = pixel / 111
        lat_cell = ((lat + 90) %/% pixel_lat) * pixel_lat - 90

        pixel_lon = pixel / (111 * cos(lat * pi / 180))
        lon_cell = ((lon + 180) %/% pixel_lon) * pixel_lon - 180
    } else {
        lat_cell = lat
        lon_cell = lon
    }

    deg2rad = function(x) x / 180 * pi

    id_full = to_integer(lat_cell, lon_cell, sorted = TRUE, add_items = TRUE, items.list = TRUE, multi.df = TRUE)

    lon_ro = deg2rad(id_full$items$lon_cell)
    lat_ro = deg2rad(id_full$items$lat_cell)

    n_cases = length(lon_ro)
    scores_ro = cpp_tapply_sum(n_cases, scores, id_full$x)

    meat = cpp_vcov_conley(scores_ro, lon_ro, lat_ro, distance = distance, cutoff = cutoff, nthreads = nthreads)

    if(sandwich){
        vcov = prepare_sandwich(bread, meat, nthreads)
    } else {
        vcov = meat
    }

    scale = if(metric == "km") 1 else 1 / 1.60934

    attr(vcov, "type_info") = paste0(" (", cutoff * scale, metric, ")")

    vcov
}

####
#### Utilities ####
####

oldargs_to_vcov = function(se, cluster, vcov, .vcov = NULL){
    # Transforms se + cluster + .vcov into a vcov call

    set_up(1)

    if(missnull(se) && missnull(cluster) && missnull(.vcov)){
        if(!missnull(vcov)){
            return(vcov)
        } else {
            return(NULL)
        }
    }

    if(!missnull(vcov)){

        if(is_function_in_it(vcov)){
            return(vcov)

        } else {
            id = c(!missnull(se), !missnull(cluster), !missnull(.vcov))
            msg = c("se", "cluster", ".vcov")[id]
            stop_up("You cannot use the argument 'vcov' in combination with ", enumerate_items(msg, "or.quote"), ".")
        }
    }

    if(!missnull(.vcov)){

        if(!is_function_in_it(.vcov)){
            if(!missnull(se) || !missnull(cluster)){
                id = c(!missnull(se), !missnull(cluster))
                msg = c("se", "cluster")[id]
                stop_up("You cannot use the argument '.vcov' in combination with ", enumerate_items(msg, "or.quote"), ".")
            }
        }

        check_arg(.vcov, "matrix | function")

        vcov = .vcov
        attr(vcov, "deparsed_arg") = fetch_arg_deparse(".vcov")

    } else {
        all_vcov = getOption("fixest_vcov_builtin")
        all_vcov_names = unlist(lapply(all_vcov, `[[`, "name"))
        all_vcov_names = all_vcov_names[nchar(all_vcov_names) > 0]

        if(missnull(se)) se = "cluster"
        check_value_plus(se, "match", .choices = all_vcov_names, .prefix = "Argument 'se' (which has been replaced by arg. 'vcov')")

        if(missnull(cluster)){
            vcov = se

        } else if(!se %in% c("cluster", "twoway", "threeway", "frouway")){
            stop_up("The VCOV requested with the argument se = '", se, "' is not compatible with the use of the argument 'cluster' (i.e. you have to choose!).")

        } else {
            if(inherits(cluster, "formula")){
                vcov = cluster
            } else if(is.character(cluster) && length(cluster) <= 4){
                vcov = .xpd(lhs = "cluster", rhs = cluster)
            } else {
                # We create a vcov request
                vcov = vcov_cluster(cluster = cluster)
            }
        }

        assign("se", NULL, parent.frame())
        assign("cluster", NULL, parent.frame())
    }

    vcov
}

catch_fun_args = function(fun, dots, exclude_args = NULL, erase_args = FALSE, keep_dots = FALSE){

    arg_names = formalArgs(fun)

    if(any(exclude_args %in% arg_names)){
        fname = deparse(substitute(fun))
        pblm = intersect(exclude_args, arg_names)
        stop_up("When the argument '", fname, "' is a function, this function cannot contain any argument named ", enumerate_items(pblm, "or.quote"), " since ", ifsingle(pblm, "this name is", "these names are"), " reserved.")
    }

    # we keep only the variables in formalArgs
    if(keep_dots){
        arg_list = dots
    } else {
        vars = intersect(names(dots), arg_names)
        arg_list = list()
        if(length(vars) > 0){
            arg_list[vars] = dots[vars]
        }
    }

    # variables from the main function call
    sysOrigin = sys.parent()
    mc = match.call(sys.function(sysOrigin), sys.call(sysOrigin), expand.dots = FALSE)
    args_previous = setdiff(names(mc)[-1], c("", "..."))

    for(var in intersect(args_previous, arg_names)){
        arg_list[[var]] = eval(str2expression(var), parent.frame())
        # NOTA: shall I put to null these arguments? Since they are now used in
        # this function... I don't know... I think it's a good idea but there can be problems.
        # Let's do it then! I'm just adding the argument erase_args.
        #
        # I think it's a good idea bc it's super easy to do it here while it's cumbersome to
        # do it in the calling fun

        if(erase_args) assign(var, NULL, parent.frame())
    }

    arg_list
}


slide_args = function(x, ...){
    # if the first argument is implicitly provided
    # slides the arguments by one slot
    #
    # we assign directly in the parent frame

    sysOrigin = sys.parent()
    sc = sys.calls()[[sysOrigin]]
    mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin))

    dots = list(...)

    if(!"x" %in% names(sc)){
        if(!missing(x)){
            if(inherits(x, "fixest")){
                # OK, regular x
            } else {
                # => implicit specification
                # sliding

                implicit_args = setdiff(names(mc), c(names(sc), ""))
                remaining_args = setdiff(names(dots), names(sc))

                for(i in seq_along(implicit_args)){
                    if(implicit_args[i] == "x"){
                        value = x
                    } else {
                        value = dots[[implicit_args[i]]]
                    }

                    assign(remaining_args[i], value, parent.frame())
                }

                assign("x", NULL, parent.frame())

            }
        } else {
            assign("x", NULL, parent.frame())
        }
    } else {
        if(missing(x)){
            assign("x", NULL, parent.frame())
        } else {
            check_arg(x, "class(fixest)", .up = 1)
        }
    }
}


prepare_sandwich = function(bread, meat, nthreads = 1){
    cpppar_matprod(cpppar_matprod(bread, meat, nthreads), bread, nthreads)
}

check_set_cutoff = function(cutoff){

    set_up(1)

    if(is.null(cutoff)){
        return(NULL)
        # stop("In vcov.fixest, Conley VCOV cannot be computed: the argument 'cutoff' (the cutoff distance in km) must be provided.", call. = FALSE)
    }

    check_value(cutoff, "numeric scalar GE{0} | character scalar",
                .message = "In vcov.fixest, Conley VCOV cannot be computed: the argument 'cutoff' must be a numeric for kilometers, or a number in character form with the 'mi' suffix for miles (ex: 100 means 100km, '100mi' means 100 miles).")

    metric = attr(cutoff, "metric")
    if(is.null(metric)) metric = "km"

    if(is.character(cutoff)){
        if(!grepl("^[[:digit:]]+ ?((m|mi|mil|mile|miles)|(k|km))?$", cutoff)){
            stop("In vcov.fixest, Conley VCOV cannot be computed: the argument 'cutoff', equal to '", cutoff,  "' is not valid. It must be either a numeric scalar or a number in character form with the suffix 'mi' (ex: 100 means 100km, '100mi' means 100 miles).", call. = FALSE)
        }

        if(grepl("k", cutoff)){
            cutoff = as.numeric(gsub("k.*", "", cutoff))
        } else if(grepl("m", cutoff)){
            cutoff = as.numeric(gsub("m.*", "", cutoff)) * 1.60934
            metric = "mi"
        } else {
            cutoff = as.numeric(cutoff)
        }
    }

    attr(cutoff, "metric") = metric

    cutoff
}



gen_vcov_aliases = function(){
    # only for vcov functions having one or more main argument

    all_vcov = getOption("fixest_vcov_builtin")

    fun_core = '
    #\' @rdname __RDNAME__
    __FUN__ = function(__ARGS__){
        extra_args = list(__EXTRA__)
        vcov_request = list(vcov = "__VCOV__", extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
        vcov_request
    }'


    text = c()
    all_fun_names = c()

    for(i in seq_along(all_vcov)){
        vcov_current = all_vcov[[i]]

        arg_main = vcov_current$arg_main
        if(!is.null(arg_main)){
            vcov_names_all = vcov_current$name

            vcov_name = vcov_names_all[1]
            rdname = vcov_current$rdname
            args = paste0(arg_main, " = NULL", collapse = ", ")
            extra = paste0(arg_main, " = ", arg_main, collapse = ", ")

            core = gsub("__RDNAME__", rdname, fun_core)
            core = gsub("__VCOV__", vcov_name, core)
            core = gsub("__ARGS__", args, core)
            core = gsub("__EXTRA__", extra, core)

            for(fname in vcov_names_all){
                text = c(text, gsub("__FUN__", fname, core, fixed = TRUE))
                all_fun_names = c(all_fun_names, fname)
            }
        }
    }

    # Writing the functions
    intro = c("# Do not edit by hand\n# => aliases some VCOV functions\n\n\n")

    text = c(intro, text)

    update_file("R/VCOV_aliases.R", text)

    # We also add the exports to the name space
    header = "# Auto-exports::vcov_aliases"
    fun_line = paste0("export(", paste0(all_fun_names, collapse = ", "), ")")

    NAMESPACE = readLines("NAMESPACE")

    if(header %in% NAMESPACE){
        # update
        qui = which(NAMESPACE == header)
        NAMESPACE[qui + 1] = fun_line
    } else {
        # creation
        NAMESPACE = c(NAMESPACE, "\n\n", header, fun_line, "\n\n")
    }

    update_file("NAMESPACE", NAMESPACE)


}


great_circle_distance = function (lat1, long1, lat2, long2){
    deg2rad = function(x) x / 180 * pi

    lat1 = deg2rad(lat1)
    lat2 = deg2rad(lat2)
    long1 = deg2rad(long1)
    long2 = deg2rad(long2)

    R = 6371
    delta.long = (long2 - long1)
    delta.lat = (lat2 - lat1)
    a = sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
    c = 2 * asin(pmin(1, sqrt(a)))
    d = R * c

    return(d)
}


cutoff_deduce = function(lat, lon){
    # To deduce the cutoff, we apply a basic rule of thumb

    quoi = unique(data.frame(lat, lon))

    order_latlon = with(quoi, order(lat, lon))
    order_lonlat = with(quoi, order(lon, lat))

    quoi_latlon = quoi[order_latlon, ]
    quoi_lonlat = quoi[order_lonlat, ]

    N = nrow(quoi)
    keep = 1:(N - 3)

    # To increase stability, we take the closest among the 3 closest lat or long
    d_latlon_1 = with(quoi_latlon, great_circle_distance(lat[-N], lon[-N], lat[-(1:1)], lon[-(1:1)]))
    d_latlon_2 = with(quoi_latlon, great_circle_distance(lat[-((N - 1):N)], lon[-((N - 1):N)], lat[-(1:2)], lon[-(1:2)]))
    d_latlon_3 = with(quoi_latlon, great_circle_distance(lat[-((N - 2):N)], lon[-((N - 2):N)], lat[-(1:3)], lon[-(1:3)]))

    d_lonlat_1 = with(quoi_lonlat, great_circle_distance(lat[-N], lon[-N], lat[-1], lon[-1]))
    d_lonlat_2 = with(quoi_lonlat, great_circle_distance(lat[-((N - 1):N)], lon[-((N - 1):N)], lat[-(1:2)], lon[-(1:2)]))
    d_lonlat_3 = with(quoi_lonlat, great_circle_distance(lat[-((N - 2):N)], lon[-((N - 2):N)], lat[-(1:3)], lon[-(1:3)]))

    d_latlon = pmin(d_latlon_1[keep], d_latlon_2[keep], d_latlon_3[keep])
    d_lonlat = pmin(d_lonlat_1[keep], d_lonlat_2[keep], d_lonlat_3[keep])

    # We take twice the median distance => corresponds to two units radius
    # roughly 50% of units have at least around 8 neighbors (=> we ensure a minimum)
    cutoff = (median(d_latlon) + median(d_lonlat))
    # we round it to make it nicer and more robust
    m = floor(log10(cutoff)) - 1
    if(cutoff > 20) m = max(m, 1)
    cutoff = (cutoff %/% 10**m) * 10**m

    attr(cutoff, "metric") = "km"

    cutoff
}

is_function_in_it = function(x){
    if(is.function(x)){
        return(TRUE)
    } else if(is.list(x) && !inherits(x, "fixest_vcov_request") && length(x) > 0 && is.function(x[[1]])){
        return(TRUE)
    }
    return(FALSE)
}

####
#### SET / GET ========== ####
####



#' @rdname ssc
#'
#' @param ssc.type An object of class \code{ssc.type} obtained with the function \code{\link[fixest]{ssc}}.
setFixest_ssc = function(ssc.type = ssc()){

    if(!"ssc.type" %in% class(ssc.type)){
        stop("The argument 'ssc' must be an object created by the function ssc().")
    }

    options("fixest_ssc" = ssc.type)
}

#' @rdname ssc
getFixest_ssc = function(){

    ssc = getOption("fixest_ssc")
    if(!"ssc.type" %in% class(ssc)){
        stop("The value of getOption(\"fixest_ssc\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value.")
    }

    ssc
}


#' Sets the default type of standard errors to be used
#'
#' This functions defines or extracts the default type of standard-errors to computed in \code{fixest} \code{\link[fixest:summary.fixest]{summary}}, and \code{\link[fixest:vcov.fixest]{vcov}}.
#'
#' @param no_FE Character scalar equal to either: \code{"iid"} (default), or \code{"hetero"}. The type of standard-errors to use by default for estimations without fixed-effects.
#' @param one_FE Character scalar equal to either: \code{"iid"}, \code{"hetero"}, or \code{"cluster"} (default). The type of standard-errors to use by default for estimations with \emph{one} fixed-effect.
#' @param two_FE Character scalar equal to either: \code{"iid"}, \code{"hetero"}, \code{"cluster"} (default), or \code{"twoway"}. The type of standard-errors to use by default for estimations with \emph{two or more} fixed-effects.
#' @param panel Character scalar equal to either: \code{"iid"}, \code{"hetero"}, \code{"cluster"} (default), or \code{"driscoll_kraaay"}. The type of standard-errors to use by default for estimations with the argument \code{panel.id} set up. Note that panel has precedence over the presence of fixed-effects.
#' @param all Character scalar equal to either: \code{"iid"}, or \code{"hetero"}. By default is is NULL. If provided, it sets all the SEs to that value.
#' @param reset Logical, default is \code{FALSE}. Whether to reset to the default values.
#'
#' @return
#' The function \code{getFixest_vcov()} returns a list with three elements containing the default for estimations i) without, ii) with one, or iii) with two or more fixed-effects.
#'
#' @examples
#'
#' # By default:
#' # - no fixed-effect (FE): standard
#' # - one or more FEs: cluster
#' # - panel: cluster on panel id
#'
#' data(base_did)
#' est_no_FE  = feols(y ~ x1, base_did)
#' est_one_FE = feols(y ~ x1 | id, base_did)
#' est_two_FE = feols(y ~ x1 | id + period, base_did)
#' est_panel = feols(y ~ x1 | id + period, base_did, panel.id = ~id + period)
#'
#' etable(est_no_FE, est_one_FE, est_two_FE)
#'
#' # Changing the default standard-errors
#' setFixest_vcov(no_FE = "hetero", one_FE = "iid",
#'                two_FE = "twoway", panel = "drisc")
#' etable(est_no_FE, est_one_FE, est_two_FE, est_panel)
#'
#' # Resetting the defaults
#' setFixest_vcov(reset = TRUE)
#'
#'
setFixest_vcov = function(no_FE = "iid", one_FE = "cluster", two_FE = "cluster",
                          panel = "cluster", all = NULL, reset = FALSE){

    # NOTE:
    # The default values should ALWAYS be working.
    # That's why I don't allow conley SEs nor NW SEs
    # => they can be not working at times

    check_arg_plus(no_FE,  "match(iid, hetero)")
    check_arg_plus(one_FE, "match(iid, hetero, cluster)")
    check_arg_plus(two_FE, "match(iid, hetero, cluster, twoway)")
    check_arg_plus(panel, "match(iid, hetero, cluster, DK, driscoll_kraay)")
    check_arg_plus(all,  "NULL match(iid, hetero)")
    check_arg_plus(reset, "logical scalar")

    opts = getOption("fixest_vcov_default")
    if(is.null(opts) || !is.list(opts) || reset){
        opts = list(no_FE = "iid", one_FE = "cluster", two_FE = "cluster", panel = "cluster")
    }

    if(!is.null(all)){
        opts$no_FE = opts$one_FE = opts$two_FE = opts$panel = all
    }


    args = intersect(c("no_FE", "one_FE", "two_FE", "panel"), names(match.call()))

    for(a in args){
        opts[[a]] = eval(as.name(a))
    }

    options(fixest_vcov_default = opts)

}

#' @rdname setFixest_vcov
getFixest_vcov = function(){

    vcov_default = getOption("fixest_vcov_default")

    if(is.null(vcov_default)){
        vcov_default = list(no_FE = "iid", one_FE = "cluster", two_FE = "cluster", panel = "cluster")
        options(fixest_vcov_default = vcov_default)
        return(vcov_default)
    }

    vcov_default
}






