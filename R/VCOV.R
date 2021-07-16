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
#' @return
#' It returns a \eqn{N\times N} square matrix where \eqn{N} is the number of variables of the fitted model.
#' This matrix has an attribute \dQuote{type} specifying how this variance/covariance matrix has been computed (i.e. if it was created using a heteroskedascity-robust correction, or if it was clustered along a specific factor, etc).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. \code{\link[fixest]{summary.fixest}}, \code{\link[fixest]{confint.fixest}}, \code{\link[fixest]{resid.fixest}}, \code{\link[fixest]{predict.fixest}}, \code{\link[fixest]{fixef.fixest}}.
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 fixed-effects)
#' est_pois = femlm(Euros ~ log(dist_km) + log(Year) | Origin + Destination +
#'                  Product, trade)
#'
#' # By default, in the presence of FEs
#' # the VCOV is clustered along the first FE
#' vcov(est_pois)
#'
#' # Heteroskedasticity-robust VCOV
#' vcov(est_pois, se = "hetero")
#'
#' # "clustered" VCOV (with respect to the Product factor)
#' vcov(est_pois, se = "cluster", cluster = trade$Product)
#' # another way to make the same request:
#' # note that previously arg. se was optional since deduced from arg. cluster
#' vcov(est_pois, cluster = "Product")
#' # yet another way:
#' vcov(est_pois, cluster = ~Product)
#'
#' # Another estimation without fixed-effects:
#' est_pois_simple = femlm(Euros ~ log(dist_km) + log(Year), trade)
#'
#' # We can still get the clustered VCOV,
#' # but we need to give the argument cluster:
#' vcov(est_pois_simple, cluster = ~Product)
#'
#'
vcov.fixest = function(object, vcov, se, cluster, dof = NULL, attr = FALSE, forceCovariance = FALSE,
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

    if(!any(c("keep_se_info", "only_varnames") %in% names(match.call(expand.dots = TRUE)))){
        # means NOT a client call
        validate_dots(suggest_args = c("vcov", "dof"), valid_args = c("meat_only"))
    }

    # We transform se and cluster into vcov
    vcov_vars = NULL
    if(!missnull(se) || !missnull(cluster)){
        if(!missnull(vcov)){
            id = c(!missnull(se), !missnull(cluster))
            msg = c("se", "cluster")[id]
            stop("You cannot use the argument 'vcov' in combination with ", enumerate_items(msg, "or.quote"), ".")
        }

        check_arg_plus(se, "NULL{'cluster'} match", .choices = c("standard", "iid", "white", "hc1", "hetero", "cluster", "twoway", "threeway", "fourway"), .message = "Argument 'se' (which has been replaced by arg. 'vcov') should be equal to one of 'iid', 'hetero', 'cluster', 'twoway', 'threeway' or 'fourway'.")

        se = switch(se,
                    white = "hetero",
                    hc1 = "hetero",
                    standard = "iid",
                    se)

        if(se %in% c("standard", "hetero") || missnull(cluster)){
            vcov = se
        } else {
            if(inherits(cluster, "formula")){
                vcov = cluster
            } else {
                vcov = "cluster"
                vcov_vars = cluster
            }
        }
    }

    meat_only = isTRUE(dots$meat_only)

    if(!is.null(object$onlyFixef)){
        # means that the estimation is done without variables
        return(NULL)
    }

    # If it's a summary => we give the vcov directly without further computation! except if arguments are provided which would mean that the user wants the new vcov
    if(isTRUE(object$summary) && missnull(vcov) && missnull(dof)){
        vcov = object$cov.scaled
        if(!is_attr) {
            attr(vcov, "se_info") = NULL
            attr(vcov, "dof.type") = NULL
            attr(vcov, "dof.K") = NULL
        }
        return(vcov)
    }


    # Default behavior vcov:
    suffix = ""
    if(missnull(vcov)){

        vcov_default = getFixest_vcov()

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


    if(isTRUE(object$lean) && !identical(vcov, "standard")){
        # we can't compute the SE because scores are gone!
        # LATER: recompute the scores (costly but maybe only possibility for user?)
        stop("VCOV of 'lean' fixest objects cannot be computed. Please re-estimate with 'lean = FALSE'.")
    }

    ####
    #### vcov parsing ####
    ####

    all_vcov = getOption("fixest_vcov_builtin")

    all_vcov_names = unlist(lapply(all_vcov, `[[`, "name"))
    all_vcov_names = all_vcov_names[nchar(all_vcov_names) > 0]

    # Checking the value of vcov
    check_arg_plus(vcov, "match | formula | function | matrix", .choices = all_vcov_names)

    # WIP => to implement (not that difficult)
    if(is.function(vcov)){
        stop("Argument 'vcov' does not accept functions yet.")

        return(vcov__)
    }

    if(is.matrix(vcov)){
        # Check that this makes sense
        nr = ncol(object$cov.iid)
        check_value(vcov, "square matrix nrow(value)", .value = nr)

        return(vcov)
    }

    if(inherits(vcov, "formula")){

        vcov_fml = vcov

        if(length(vcov_fml) == 2){
            vcov = ""
            vcov_vars = fml2varnames(vcov_fml, combine_fun = TRUE)
        } else {
            vcov = deparse_long(vcov_fml[[2]])

            check_arg_plus(vcov, "match(standard, iid, hetero, hc1, white, cluster, twoway, threeway, fourway, hac, conley, conley_hac)", .message = "If a formula, the arg. 'vcov' must be of the form 'vcov_type ~ vars'. The vcov_type must be a supported VCOV type.")

            vcov_vars = fml2varnames(vcov_fml[c(1, 3)], combine_fun = TRUE)

            qui_dof = grepl("dof(", vcov_vars, fixed = TRUE)
            if(any(qui_dof)){
                dof_txt = vcov_vars[qui_dof]
                dof = eval(str2lang(dof_txt))
                vcov_vars = vcov_vars[!qui_dof]
            }
        }

    }

    # Here vcov **must** be a character vector

    vcov_id = which(sapply(all_vcov, function(x) vcov %in% x$name))

    if(length(vcov_id) != 1){
        stop("Unexpected problem in the selection of the VCOV. This is an internal error. Could you report?")
    }

    vcov_select = all_vcov[[vcov_id]]

    if(length(vcov_select$vars) > 0){
        # We find out the vcov

        if(only_varnames == FALSE){
            data = fetch_data(object)
            data_names = names(data)
        }

        patterns_split = strsplit(vcov_select$patterns, " ?\\+ ?")
        n_patterns = lengths(patterns_split)
        if(!length(vcov_vars) %in% n_patterns){
            stop("In the argument 'vcov', the number of variables in the RHS of the formula (", length(vcov_vars), ") is not valid, it should correspond to ", ifsingle(n_patterns, "", "one of "), enumerate_items(vcov_select$patterns[n_patterns != 0], "quote.or"), ".")
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

            if(!vcov_var_name %in% pattern && isTRUE(vcov_var_value$optional)){
                next

            } else if(vcov_var_name %in% pattern){
                # => provided by the user, we find which one it corresponds to
                # based on the pattern

                vname = vcov_vars[which(pattern == vcov_var_name)]

                vname_all = all.vars(str2expression(vname))
                if(!all(vname_all %in% data_names)){
                    pblm = setdiff(vname_all, data_names)
                    stop("The variable", enumerate_items(pblm, "s"), " used to compute the VCOV ", plural_len(pblm, "is"), " not in the original data set. Only variables in the data set can be used.")
                }

                var_names_all[vcov_var_name] = vname
                if(only_varnames) next

                var_values_all[[vcov_var_name]] = trim_obs_removed(data[[vname]], object)

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
                            # stop("The variable(s) used to compute the VCOV are typically deduced from the fixed-effects used in the estimation, but this one had no fixed-effect! Hence, please provide the variable name(s) explicitly in the RHS of the 'vcov' formula.")
                        }

                        if(length(object$fixef_vars) < guesses$fixef){
                            msg = setNames(nm = paste0("the ", n_th(guesses$fixef), " fixed-effect"),
                                           paste0("the estimation had only ", length(object$fixef_vars), " fixed-effect", plural_len(object$fixef_vars)))
                            err_msg = c(err_msg, msg)
                            next
                            # stop("One variable used to compute the VCOV is typically deduced from the ", n_th(guesses$fixef), " fixed-effect used in the estimation, but this one had only ", length(object$fixef_vars), " fixed-effects! Hence, please provide all variable names explicitly in the RHS of the 'vcov' formula.")
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
                            # stop("One variable used to compute the VCOV is typically deduced from the 'panel.id' used in the estimation, but this one had no 'panel.id' specified! Hence, please provide the variable names explicitly in the RHS of the 'vcov' formula.")
                        }

                        vname = object$panel.id[guesses$fixef]

                        vname_all = all.vars(str2expression(vname))
                        if(!all(vname_all %in% data_names)){
                            pblm = setdiff(vname_all, data_names)
                            msg = setNames(nm = "the 'panel.id' identifiers",
                                           paste0("the variable", enumerate_items(pblm, "s"), " set in 'panel.id' ", plural_len(pblm, "is"), " not in the data set"))
                            err_msg = c(err_msg, msg)
                            next
                            # stop("The variable ", vname, " used to compute the VCOV is not in the original data set. Only variables in the data set can be used.")
                        }

                        # ATTENTION:
                        # if the panel.id variable has NA values, and there are no lags used in the estimation
                        # (which would trigger obs removal), then we cannot use the panel.id
                        # for the VCOV reliably since the sample would differ from the sample used in the
                        # estimation

                        var_names_all[vcov_var_name] = vname
                        if(only_varnames) break

                        var_values_all[[vcov_var_name]] = trim_obs_removed(data[[vname]], object)
                        break

                    } else if(type == "regex"){

                        ok = FALSE
                        msg = NULL
                        for(pat in guesses$regex){
                            var_id = which(grepl(pat, data_names))

                            if(length(var_id) == 0){
                                if(is.null(msg)){
                                    msg = setnames(nm = "the variable names of the data set",
                                                   paste0("no match was found for ", vcov_var_value$label))
                                }
                            } else if(length(var_id) == 1){
                                ok = TRUE
                                vname = data_names[var_id]
                                break
                            } else {
                                msg = setnames(nm = "the variable names of the data set",
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

                if(length(err_msg) > 0){

                    if(length(err_msg) == 1){
                        expect = names(err_msg)
                        pblm = err_msg
                    } else {
                        expect = enumerate_items(names(err_msg), "enum a.or")
                        pblm = enumerate_items(err_msg, "enum a")
                    }

                    stop("To compute the ", vcov_select$vcov_label, " VCOV, we need a variable for the ",
                         vcov_var_value$label, ". Since you didn't provide it in the formula, we typically deduce it from ", expect,
                         ". Problem: ", pblm, ".")
                }
            }
        }

        if(only_varnames){
            return(var_names_all)
        }

        # Extra checking + putting into int if requested
        for(i in seq_along(var_values_all)){
            vcov_var_name = names(var_values_all)[i]
            vname = var_names_all[vcov_var_name]
            value = var_values_all[[i]]

            vcov_var_value = vcov_select$vars[[vcov_var_name]]

            if(anyNA(value)){
                stop("The variable '", vname, "' used to estimate the VCOV (employed as ", vcov_var_value$label, ") has NA values which would lead to a sample used to compute the VCOV different from the sample used to estimate the parameters. This would lead to wrong inference.\nPossible solutions: i) ex ante prune them or ii) use the argument 'vcov' at estimation time.")
            }

            if(vcov_var_value$to_int){
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

    } else if(only_varnames){
        return(character(0))
    }


    # Checking the nber of threads
    if(!missing(nthreads)) nthreads = check_set_nthreads(nthreads)


    ####
    #### scores ####
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
    #### bread ####
    ####

    n = object$nobs

    if(object$method_type == "feols"){
        if(vcov != "iid"){
            bread = object$cov.iid / object$sigma2
        } else {
            bread = object$cov.iid / ((n - 1) / (n - object$nparams))
        }
    } else {
        bread = object$cov.iid
    }

    if(anyNA(bread)){

        if(!forceCovariance){
            last_warn = getOption("fixest_last_warning")
            if(is.null(last_warn) || (proc.time() - last_warn)[3] > 1){
                warning("Standard errors are NA because of likely presence of collinearity.", call. = FALSE)
            }

            attr(bread, "type") = "NA (not-available)"

            return(bread)
        } else {
            info_inv = cpp_cholesky(object$hessian)
            if(!is.null(info_inv$all_removed)){
                # Means all variables are collinear! => can happen when using FEs
                stop("All variables have virtually no effect on the dependent variable. Covariance is not defined.")
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

    }

    ####
    #### vcov no adj ####
    ####

    # we compute the vcov. The adjustment (which is a pain in the neck) will come after that
    # Here vcov is ALWAYS a character scalar

    fun_name = vcov_select$fun_name
    args = list(bread = bread, scores = scores, vars = vcov_vars)
    vcov_noAdj = do.call(fun_name, args)


    ####
    #### DoF adj ####
    ####

    # dof is a dof object in here

    # DoF related => we accept NULL
    check_arg_plus(dof, "NULL{getFixest_dof()} class(dof.type)", .message = "The argument 'dof.type' must be an object created by the function dof().")

    dof.fixef.K = dof$fixef.K
    dof.adj = dof$adj
    is_exact = dof$fixef.force_exact
    is_cluster = dof$cluster.adj
    is_cluster_min = dof$cluster.df == "min"
    is_t_min = dof$t.df == "min"

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

    # How do we choose K? => argument dof

    nested_vars = sapply(vcov_select$vars, function(x) isTRUE(x$rm_nested))
    any_nested_var = length(nested_vars) > 0 && length(vcov_vars > 0) && any(nested_vars[names(vcov_vars)])

    if(dof$fixef.K == "none"){
        # we do it with "minus" because of only slopes
        K = object$nparams
        if(n_fe_ok > 0){
            K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1))
        }
    } else if(dof$fixef.K == "full" || !any_nested_var){
        K = object$nparams
        if(dof$fixef.force_exact && n_fe >= 2 && n_fe_ok >= 1){
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

    if(dof$fixef.K == "nested" && n_fe_ok > 0 && any_nested_var){
        # OK, let's go checking....
        # We always try to minimize computation.
        # So we maximize deduction and apply computation only in last resort.

        nested_vcov_var_names = intersect(names(nested_vars[nested_vars]), names(vcov_vars))
        nested_var_names = var_names_all[nested_vcov_var_names]




        if(!any(grepl("^", object$fixef_vars, fixed = TRUE))){
            # simple cases only

            id_nested = which(object$fixef_vars %in% cluster)

        }


        id_nested = c()
        for(v_nested in nested_var_names){



        }


        if(all(nested_var_names %in% object$fixef_vars)){

            check_nested = FALSE
            # We do that to avoid checking nestedness later
            if(all(object$fixef_vars %in% cluster)){
                # everyone nested (also works for var1^var2)
                is_nested = 1:length(object$fixef_id)

            } else if(!any(grepl("^", object$fixef_vars, fixed = TRUE))){
                # simple cases only
                is_nested = which(names(object$fixef_id) %in% cluster)

            } else if(!any(grepl("^", cluster, fixed = TRUE))){
                # simple cases in cluster
                check_var_in_there = function(x){
                    # cluster is a global
                    if(x %in% cluster){
                        return(TRUE)

                    } else if(grepl("^", x, fixed = TRUE)){
                        x_split = strsplit(x, "^", fixed = TRUE)[[1]]
                        if(any(x_split %in% cluster)){
                            return(TRUE)
                        }
                    }

                    return(FALSE)
                }

                is_nested = which(sapply(names(object$fixef_id), check_var_in_there))
            } else {
                # too complex to apply tricks => we make real check
                check_nested = TRUE
            }

            cluster = object$fixef_id[cluster]

            do.unclass = FALSE

        }


    }

    # We recompute K depending on nesting
    if(dof$adj && dof$fixef.K == "nested" && n_fe_ok >= 1){

        if(check_nested){
            # we need to find out which is nested
            is_nested = which(cpp_check_nested(object$fixef_id, cluster, object$fixef_sizes, n = n) == 1)
        } else {
            # no need to compute is_nested,
            # we created it earlier
        }

        if(length(is_nested) == n_fe){
            # All FEs are removed, we add 1 for the intercept
            K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1)) + 1
        } else {
            if(is_exact && n_fe >= 2){
                fe = fixef(object, notes = FALSE)
                nb_ref = attr(fe, "references")

                # Slopes are a pain in the neck!!!
                if(length(is_nested) > 1){
                    id_nested = intersect(names(nb_ref), names(object$fixef_id)[is_nested])
                    nb_ref[id_nested] = object$fixef_sizes[id_nested]
                }

                total_refs = sum(nb_ref)

                K = K - total_refs
            } else {
                K = K - (sum(fixef_sizes_ok[is_nested]) - sum(fixef_sizes_ok[is_nested] > 0))
            }
        }

        # below for consistency => should not be triggered
        K = max(K, length(object$coefficients) + 1)

        correction.dof = (n - 1) / (n - K)
    }



    # Small sample adjustment
    correction.dof = ifelse(dof$adj, (n - 1) / (n - K), 1)





    #
    # Core function
    #


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

    # How do we choose K? => argument dof

    if(dof.fixef.K == "none"){
        # we do it with "minus" because of only slopes
        K = object$nparams
        if(n_fe_ok > 0){
            K = K - (sum(fixef_sizes_ok) - (n_fe_ok - 1))
        }
    } else if(dof.fixef.K == "full" || se.val %in% c("iid", "hetero")){
        K = object$nparams
        if(is_exact && n_fe >= 2 && n_fe_ok >= 1){
            fe = fixef(object, notes = FALSE)
            K = K + (n_fe_ok - 1) - sum(attr(fe, "references"))
        }
    } else {
        # nested
        # we delay the adjustment
        K = object$nparams
    }

    # Small sample adjustment
    correction.dof = ifelse(dof.adj, (n - 1) / (n - K), 1)


    # information on the variable used for the clustering
    type_info = ""

    is_nested = c()
    if(se.val == "iid"){

        vcov = bread * correction.dof

    } else if(se.val == "hetero"){

        # we make a n/(n-1) adjustment to match vcovHC(type = "HC1")
        if(meat_only){
            meat = cpppar_crossprod(scores, 1, nthreads)
            return(meat)

        } else {
            vcov = cpppar_crossprod(cpppar_matprod(scores, bread, nthreads), 1, nthreads) * correction.dof * ifelse(is_cluster, n/(n-1), 1)
        }

        dimnames(vcov) = dimnames(bread)

    } else {
        # Clustered SE!
        nway = switch(se.val, cluster=1, twoway=2, threeway=3, fourway=4)

        # MISC

        # used twice later:
        msgRemoved = ""
        if(!is.null(object$call$na.rm) && object$call$na.rm){
            msgRemoved = " (additionnaly from the observations removed from the original estimation)"
        }

        #
        # Controls
        #

        # Controlling the clusters
        do.unclass = check_nested = TRUE
        if(missing(cluster) || is.null(cluster)){

            if(is.null(object$fixef_id)){
                stop("To display clustered standard errors, you must provide the argument 'cluster'", ifelse(!missing(cluster), ", currently it is equal to NULL", ""), ".")

            } else if(length(object$fixef_id) < nway) {
                stop(nway, "-way clustering is asked for but the estimation was not performed with ", nway, " or more fixed-effects: You must provide the argument 'cluster' with ", nway, " clusters.")

            } else {
                cluster = object$fixef_id[1:nway]

                type_info = paste0(" (", paste0(object$fixef_vars[1:nway], collapse = " & "), ")")

                is_nested = 1:nway

                # in that specific case, there is no need of doing unclass.factor because already done
                do.unclass = FALSE
                check_nested = FALSE
            }

        } else {

            #
            # Handle formulas
            #

            doEval = FALSE
            if("formula" %in% class(cluster)){

                cluster = formula(cluster) # regularization to check it

                if(length(cluster) != 2){
                    stop("If argument cluster is to be a formula, it must be one sided: e.g. ~dum_1+dum_2.")
                }

                all_var_names = fml2varnames(cluster)

                if(length(all_var_names) != nway){
                    stop("Asked for ", nway, "-way clustering but evaluating argument cluster leads to ", length(all_var_names), " clusters (", enumerate_items(all_var_names), "). Please provide exactly ", nway, " clusters.")
                }

                cluster = all_var_names # Now a regular character vector

                doEval = TRUE
            }

            if(length(cluster) == nway && is.character(cluster)){

                if(all(cluster %in% object$fixef_vars)){
                    # cluster == names of clusters used in the estimation
                    type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")

                    check_nested = FALSE
                    # We do that to avoid checking nestedness later
                    if(all(object$fixef_vars %in% cluster)){
                        # everyone nested (also works for var1^var2)
                        is_nested = 1:length(object$fixef_id)

                    } else if(!any(grepl("^", object$fixef_vars, fixed = TRUE))){
                        # simple cases only
                        is_nested = which(names(object$fixef_id) %in% cluster)

                    } else if(!any(grepl("^", cluster, fixed = TRUE))){
                        # simple cases in cluster
                        check_var_in_there = function(x){
                            # cluster is a global
                            if(x %in% cluster){
                                return(TRUE)

                            } else if(grepl("^", x, fixed = TRUE)){
                                x_split = strsplit(x, "^", fixed = TRUE)[[1]]
                                if(any(x_split %in% cluster)){
                                    return(TRUE)
                                }
                            }

                            return(FALSE)
                        }

                        is_nested = which(sapply(names(object$fixef_id), check_var_in_there))
                    } else {
                        # too complex to apply tricks => we make real check
                        check_nested = TRUE
                    }

                    cluster = object$fixef_id[cluster]

                    do.unclass = FALSE

                } else {
                    cluster = gsub(" *", "", cluster)
                    if(!doEval){
                        is_ok = grepl("^[\\.[:alpha:]][[:alnum:]_\\.]*(\\^[\\.[:alpha:]][[:alnum:]_\\.]*)*$", cluster)
                        if(any(!is_ok)){
                            stop("In argument cluster, only variable names and the '^' operator are accepted. The expression", enumerate_items(cluster[!is_ok], "s.is"), " not valid.\nAlternatively, you can use a list of vectors.")
                        }

                    }

                    cluster_fml = as.formula(paste0("~", paste0(cluster, collapse = " + ")))
                    all_vars = all.vars(cluster_fml)

                    if(all(all_vars %in% object$fixef_vars) || all(cluster %in% object$fixef_vars)){
                        # Means dum_1^dum_2 with dum_1 and dum_2 used as clusters

                        cluster_names = cluster
                        type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")

                        cluster = list()
                        for(i in 1:nway){
                            cname = cluster_names[i]
                            if(cname %in% object$fixef_vars){
                                cluster[[i]] = object$fixef_id[[cname]]
                            } else {
                                # combination
                                if(grepl("^", cname, fixed = TRUE)){
                                    value_text = gsub("\\^", ", ", cname)
                                    value_text = paste0("combine_clusters_fast(", value_text, ")")
                                }

                                value_call = str2lang(value_text)
                                value = eval(value_call, object$fixef_id)
                                cluster[[i]] = value
                            }
                        }

                    } else {
                        # we try to get the variable from the base used in the estimation
                        var2fetch = setdiff(all_vars, object$fixef_vars)

                        # evaluation
                        data = fetch_data(object, paste0("Cannot apply ", nway, "-way clustering with current 'cluster' argument. Variable", enumerate_items(var2fetch, "s.is.past"), " not used as fixed-effects in the estimation so ", plural_len(var2fetch, "need"), " to be taken from the data. "), " Alternatively, use a list of vectors.")

                        data = as.data.frame(data)

                        # we check the variables are there
                        # we use all_vars and not var2fetch: safer to catch all variables (case clustvar^datavar)

                        if(any(!all_vars %in% names(data))){
                            var_pblm = setdiff(all_vars, names(data))
                            stop("In argument 'cluster', the variable", enumerate_items(var_pblm, "s.is"), " not present in the original dataset. Alternatively, use a list of vectors.")
                        }

                        # we check length consistency
                        if(NROW(data) != object$nobs_origin){
                            stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(var2fetch, "s"), " in the original dataset (", deparse_long(object$call$data), "), yet the dataset doesn't have the same number of observations as was used in the estimation (", NROW(data), " instead of ", object$nobs_origin, ").")
                        }

                        data = data[, all_vars, drop = FALSE]

                        for(i in seq_along(object$obs_selection)){
                            data = data[object$obs_selection[[i]], , drop = FALSE]
                        }

                        # Final check: NAs
                        if(anyNA(data)){
                            varsNA = sapply(data, anyNA)
                            varsNA = names(varsNA)[varsNA]
                            stop("To evaluate argument 'cluster', we fetched the variable", enumerate_items(varsNA, "s"), " in the original dataset (", deparse_long(object$call$data), "). But ", ifsingle(varsNA, "this variable", "these variables"), " contain", ifsingle(varsNA, "s", ""), " NAs", msgRemoved, ". This is not allowed since the sample used to compute the SEs would be different from the sample used in the estimation.")
                        }

                        # We create the cluster list
                        cluster_names = cluster
                        type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")

                        cluster = list()
                        for(i in 1:nway){
                            cname = cluster_names[i]
                            if(cname %in% object$fixef_vars){
                                cluster[[i]] = object$fixef_id[[cname]]

                            } else if(cname %in% names(data)){
                                # data is already of the right size
                                cluster[[i]] = data[[cname]]

                            } else {
                                # combination
                                if(grepl("^", cname, fixed = TRUE)){
                                    value_text = gsub("\\^", ", ", cname)
                                    value_text = paste0("combine_clusters_fast(", value_text, ")")
                                } else {
                                    value_text = cname
                                }

                                value_call = str2lang(value_text)
                                cluster[[i]] = eval(value_call, data)
                            }
                        }

                    }

                }
            } else if(length(cluster) == nway && is.numeric(cluster)){
                # You can use a number to tell which cluster to use

                if(length(object$fixef_vars) == 0){
                    stop("You can use an integer in the argument 'cluster' only when there have been fixed-effects in the estimation. Currenlty this is not the case. Alternatively, arg. 'cluster' can be a formula, a vector of variables or a list of vectors.")
                }

                if(!all(cluster %% 1 == 0) || any(cluster < 1 | cluster > 4)){
                    msg = ifelse(!all(cluster %% 1 == 0), "it is not made of integers", "it contains values different from 1 to 4")
                    stop("Argument 'cluster' can be a numeric vector, if so it must have integer values between 1 and 4 (currently ", msg, ").")
                }

                if(length(object$fixef_vars) < max(cluster)){
                    nb_name = c("1st", "2nd", "3rd", "4th")
                    stop("In argument 'cluster', it is requested to cluster along the ", nb_name[max(cluster)], " fixed-effect, however the estimation was done with only ", length(object$fixef_vars), " fixed-effects. Alternatively, arg. 'cluster' can be a formula, a vector of variables or a list of vectors.")
                }

                # Eventually, it should be all right by now
                type_info = paste0(" (", paste0(object$fixef_vars[cluster], collapse = " & "), ")")
                cluster = object$fixef_id[cluster]

            } else if(nway == 1){
                if(!is.list(cluster) && (isVector(cluster) || is.factor(cluster))){
                    cluster = list(cluster)

                } else if(! (is.list(cluster) && length(cluster) == 1)){
                    stop("For one way clustering, the argument 'cluster' must be either the name of a cluster variable (e.g. \"dum_1\"), a vector (e.g. data$dum_1), a list containing the vector of clusters (e.g. list(data$dum_1)), or a one-sided formula (e.g. ~dum_1). Currently the class of cluster is ", enumerate_items(class(cluster)), ".")

                } else if(!is.null(names(cluster))){
                    type_info = paste0(" (", names(cluster), ")")
                }

            } else if(length(cluster) != nway){

                msgList = "a list of vectors"
                if(is.list(cluster)) msgList = "a vector of variables names"
                stop(nway, "-way clustering is asked for, but the length of argument 'cluster' is ", length(cluster), " (it should be ", nway, "). Alternatively, you can use ", msgList, " or a one-sided formula.")

            } else if(!is.list(cluster)){
                stop("For ", nway, "-way clustering, the argument 'cluster' must be either a vector of cluster variables (e.g. c(\"", paste0("dum_", 1:nway, collapse = "\", \""), "\")), a list containing the vector of clusters (e.g. data[, c(\"", paste0("dum_", 1:nway, collapse = "\", \""), "\")]), or a one-sided formula (e.g. ~", paste0("dum_", 1:nway, collapse = "+"), "). Currently the class of cluster is: ", enumerate_items(class(cluster)), ".")

            } else if(!is.null(names(cluster))){
                type_info = paste0(" (", paste0(names(cluster), collapse = " & "), ")")
            }

            cluster = as.list(cluster)
        }

        # now we check the lengths:
        n_per_cluster = sapply(cluster, length)
        if(!all(diff(n_per_cluster) == 0)){
            stop("The vectors of the argument 'cluster' must be of the same length. Currently the lengths are: ", enumerate_items(n_per_cluster), ".")
        }

        # Either they are of the same length of the data
        if(n_per_cluster[1] != object$nobs){
            # Then two cases: either the user introduces the original data and it is OK
            if(n_per_cluster[1] == object$nobs_origin){
                # We modify the clusters
                for(i in 1:nway){
                    # first we take care of the observations removed

                    for(j in seq_along(object$obs_selection)){
                        cluster[[i]] = cluster[[i]][object$obs_selection[[j]]]
                    }

                }
            } else {
                # If this is not the case: there is a problem
                stop("The length of the clusters (", fsignif(n_per_cluster[1]), ") does not match the number of observations in the estimation (", fsignif(object$nobs), ").")
            }
        }

        # final NA check
        varsNA = sapply(cluster, anyNA)
        if(any(varsNA)){
            varsNA = which(varsNA)
            nb_name = c("1st", "2nd", "3rd", "4th")
            stop("In argument cluster, the ", enumerate_items(nb_name[varsNA]), " cluster variable", ifsingle(varsNA, " contains", "s contain"), " NAs", msgRemoved, ". This is not allowed.")
        }


        #
        # Calculus
        #

        # initialisation
        vcov = bread * 0
        meat = 0

        if(do.unclass){
            for(i in 1:nway){
                cluster[[i]] = quickUnclassFactor(cluster[[i]])
            }
        }



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

                # When cluster.df == "min" => no dof here but later
                if(meat_only){
                    meat = meat + (-1)**(i+1) * vcovClust(index, bread, scores, dof = is_cluster && !is_cluster_min, do.unclass=FALSE, meat_only = TRUE, nthreads = nthreads)
                } else {
                    vcov = vcov + (-1)**(i+1) * vcovClust(index, bread, scores, dof = is_cluster && !is_cluster_min, do.unclass=FALSE, nthreads = nthreads)
                }

            }
        }

        G_min = NULL
        if(is_cluster && is_cluster_min){
            G_min = min(sapply(cluster, max))
            correction.dof = correction.dof * G_min / (G_min - 1)
        }

        if(meat_only){
            if(!is.null(G_min)) meat = meat * G_min / (G_min - 1)
            return(meat)
        }

        vcov = vcov * correction.dof

        if(is_t_min){
            if(is.null(G_min)) G_min = min(sapply(cluster, max))

            if(is_attr) base::attr(vcov, "G") = G_min
        }

    }

    if(any(diag(vcov) < 0)){
        # We 'fix' it
        vcov = mat_posdef_fix(vcov)
        message("Variance contained negative values in the diagonal and was 'fixed' (a la Cameron, Gelbach & Miller 2011).")
    }

    sd.dict = c("iid" = "IID", "hetero"="Heteroskedasticity-robust", "cluster"="Clustered", "twoway"="Two-way", "threeway"="Three-way", "fourway"="Four-way")

    if(is_attr){
        base::attr(vcov, "type") = paste0(as.vector(sd.dict[se.val]), type_info)
        base::attr(vcov, "dof.type") = paste0("dof(adj = ", dof.adj, ", fixef.K = '", dof.fixef.K, "', cluster.adj = ", is_cluster, ", cluster.df = '", dof$cluster.df, "', t.df = '", dof$t.df, "', fixef.force_exact = ", is_exact, ")")
        base::attr(vcov, "dof.K") = K
    }

    if(isTRUE(dots$keep_se_info)){
        if(missing(cluster)) cluster = NULL
        attr(vcov, "se_info") = list(se = se.val, dof = dof, cluster = cluster)
    }

    vcov
}









#' Type of degree of freedom in fixest summary
#'
#' Provides how the degrees of freedom should be calculated in \code{\link[fixest]{vcov.fixest}}/\code{\link[fixest]{summary.fixest}}.
#'
#' @param adj Logical scalar, defaults to \code{TRUE}. Whether to apply a small sample adjustment of the form \code{(n - 1) / (n - K)}, with \code{K} the number of estimated parameters. If \code{FALSE}, then no adjustment is made.
#' @param fixef.K Character scalar equal to \code{"nested"} (default), \code{"none"} or \code{"full"}. In the small sample adjustment, how to account for the fixed-effects parameters. If \code{"none"}, the fixed-effects parameters are discarded, meaning the number of parameters (\code{K}) is only equal to the number of variables. If \code{"full"}, then the number of parameters is equal to the number of variables plus the number of fixed-effects. Finally, if \code{"nested"}, then the number of parameters is equal to the number of variables plus the number of fixed-effects that *are not* nested in the clusters used to cluster the standard-errors.
#' @param fixef.force_exact Logical, default is \code{FALSE}. If there are 2 or more fixed-effects, these fixed-effects they can be irregular, meaning they can provide the same information. If so, the "real" number of parameters should be lower than the total number of fixed-effects. If \code{fixef.force_exact = TRUE}, then \code{\link[fixest]{fixef.fixest}} is first run to determine the exact number of parameters among the fixed-effects. Mostly, panels of the type individual-firm require \code{fixef.force_exact = TRUE} (but it adds computational costs).
#' @param cluster.adj Logical scalar, default is \code{TRUE}. How to make the small sample correction when clustering the standard-errors? If \code{TRUE} a \code{G/(G-1)} correction is performed with \code{G} the number of cluster values.
#' @param cluster.df Either "conventional" or "min" (default). Only relevant when the variance-covariance matrix is two-way clustered (or higher). It governs how the small sample adjustment for the clusters is to be performed. [Sorry for the jargon that follows.] By default a unique adjustment is made, of the form G_min/(G_min-1) with G_min the smallest G_i. If \code{cluster.df="conventional"} then the i-th "sandwich" matrix is adjusted with G_i/(G_i-1) with G_i the number of unique clusters.
#' @param t.df Either "conventional" or "min" (default). Only relevant when the variance-covariance matrix is clustered. It governs how the p-values should be computed. By default, the degrees of freedom of the Student t distribution is equal to the minimum size of the clusters with which the VCOV has been clustered. If \code{t.df="conventional"}, then the degrees of freedom of the Student t distribution is equal to the number of observations minus the number of estimated variables.
#'
#' @details
#'
#' The following vignette: \href{https://lrberge.github.io/fixest/articles/standard_errors.html}{On standard-errors}, describes in details how the standard-errors are computed in \code{fixest} and how you can replicate standard-errors from other software.
#'
#' @return
#' It returns a \code{dof.type} object.
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
#' # To get the same SEs, we need to use dof(adj = FALSE)
#'
#' res_pois = fepois(round(Petal.Length) ~ Petal.Width + Species, iris)
#' res_glm = glm(round(Petal.Length) ~ Petal.Width + Species, iris, family = poisson())
#' vcov(res_pois, dof = dof(adj = FALSE)) / vcov(res_glm)
#'
#' # Same example with the Gamma
#' res_gamma = feglm(round(Petal.Length) ~ Petal.Width + Species, iris, family = Gamma())
#' res_glm_gamma = glm(round(Petal.Length) ~ Petal.Width + Species, iris, family = Gamma())
#' vcov(res_gamma, dof = dof(adj = FALSE)) / vcov(res_glm_gamma)
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
#' attributes(vcov(est, attr = TRUE))[c("dof.type", "dof.K")]
#'
#'
#' # fixef.K = FALSE
#' #  => adjustment K = 1 (i.e. only x)
#' summary(est, dof = dof(fixef.K = "none"))
#' attr(vcov(est, dof = dof(fixef.K = "none"), attr = TRUE), "dof.K")
#'
#'
#' # fixef.K = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 1 (i.e. x + fe1 + fe2 - 1 restriction)
#' summary(est, dof = dof(fixef.K = "full"))
#' attr(vcov(est, dof = dof(fixef.K = "full"), attr = TRUE), "dof.K")
#'
#'
#' # fixef.K = TRUE & fixef.force_exact = TRUE
#' #  => adjustment K = 1 + 3 + 5 - 2 (i.e. x + fe1 + fe2 - 2 restrictions)
#' summary(est, dof = dof(fixef.K = "full", fixef.force_exact = TRUE))
#' attr(vcov(est, dof = dof(fixef.K = "full", fixef.force_exact = TRUE), attr = TRUE), "dof.K")
#'
#' # There are two restrictions:
#' attr(fixef(est), "references")
#'
#' #
#' # To permanently set the default dof:
#' #
#'
#' # eg no small sample adjustment:
#' setFixest_dof(dof(adj = FALSE))
#'
#' # Factory default
#' setFixest_dof()
#'
dof = function(adj = TRUE, fixef.K = "nested", cluster.adj = TRUE, cluster.df = "min", t.df = "min", fixef.force_exact = FALSE){

    check_arg_plus(adj, "loose logical scalar conv")
    check_arg_plus(fixef.K, "match(none, full, nested)")
    check_arg_plus(cluster.df, "match(conventional, min)")
    check_arg_plus(t.df, "match(conventional, min)")
    check_arg(fixef.force_exact, cluster.adj, "logical scalar")

    res = list(adj = adj, fixef.K = fixef.K, cluster.adj = cluster.adj, cluster.df = cluster.df, t.df = t.df, fixef.force_exact = fixef.force_exact)
    class(res) = "dof.type"
    res
}



####
#### Internal ####
####




vcovClust = function (cluster, myBread, scores, dof=FALSE, do.unclass=TRUE, meat_only = FALSE, nthreads = 1){
    # Internal function: no need for controls, they come beforehand
    # - cluster: the vector of dummies
    # - myBread: original vcov
    # - scores
    # Note: if length(unique(cluster)) == n (i.e. White correction), then the dof are such that vcovClust is equivalent to vcovHC(res, type="HC1")
    # Source: http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf
    #         Cameron & Miller -- A Practitioner’s Guide to Cluster-Robust Inference

    n = NROW(scores)

    # Control for cluster type
    if(do.unclass){
        cluster = quickUnclassFactor(cluster)
    }

    Q = max(cluster)
    RightScores = cpp_tapply_sum(Q, scores, cluster)

    # Finite sample correction:
    if(dof){
        dof_value = Q / (Q - 1)
    } else {
        dof_value = 1
    }

    if(meat_only){
        res = cpppar_crossprod(RightScores, 1, nthreads) * dof_value
        return(res)
    }

    xy = cpppar_matprod(RightScores, myBread, nthreads)
    res = cpppar_crossprod(xy, 1, nthreads) * dof_value
    res
}


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

    vcov_clust_setup = list(name = c("cluster", ""), fun_name = "vcov_cluster_internal", vcov_label = "clustered")
    vcov_clust_setup$vars = list(cl1 = list(guess_from = list(fixef = 1), label = "clusters", to_int = TRUE, rm_nested = TRUE),
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

    vcov_twoway_setup = list(name = "twoway", fun_name = "vcov_cluster_internal", vcov_label = "2-way clustered")
    vcov_twoway_setup$vars = list(cl1, cl2)
    vcov_twoway_setup$patterns = c("", "cl1 + cl2")

    vcov_threeway_setup = list(name = "threeway", fun_name = "vcov_cluster_internal", vcov_label = "3-way clustered")
    vcov_threeway_setup$vars = list(cl1, cl2, cl3)
    vcov_threeway_setup$patterns = c("", "cl1 + cl2 + cl3")

    vcov_fourway_setup = list(name = "fourway", fun_name = "vcov_cluster_internal", vcov_label = "4-way clustered")
    vcov_fourway_setup$vars = list(cl1, cl2, cl3, cl4)
    vcov_fourway_setup$patterns = c("", "cl1 + cl2 + cl3 + cl4")


    #
    # hac (panel)
    #

    vcov_hac_setup = list(name = "hac", fun_name = "vcov_hac_internal", vcov_label = "HAC")
    # The variables
    id = list(guess_from = list(panel.id = 1), label = "panel ID", to_int = TRUE, rm_nested = TRUE)
    time = list(guess_from = list(panel.id = 2), label = "time", expected_type = "numeric vector")
    vcov_hac_setup$vars = list(id = id, time = time)
    vcov_hac_setup$arg_main = "lag"
    vcov_hac_setup$patterns = c("", "id + time")


    #
    # conley
    #

    vcov_conley_setup = list(name = "conley", fun_name = "vcov_conley_internal", vcov_label = "Conley")
    # The variables
    lat = list(guess_from = list(regex = c("^lat(itude)?$", "^lat_.+")),
               label = "latitude",
               expected_type = "numeric vector",
               rm_nested = TRUE)

    lng = list(guess_from = list(regex = c("^lng$", "^long?(itude)?$", "^(lng|lon|long)_.+")),
               label = "longitude",
               expected_type = "numeric vector",
               rm_nested = TRUE)

    vcov_conley_setup$vars = list(lat = lat, lng = lng)
    vcov_conley_setup$arg_main = "cutoff"
    vcov_conley_setup$patterns = c("", "lat + lng")


    #
    # conley hac
    #

    vcov_conley_hac_setup = list(name = c("conley_hac", "hac_conley"), fun_name = "vcov_conley_hac_internal", vcov_label = "Conley-HAC")
    # The variables (already defined earlier)
    vcov_conley_hac_setup$vars = list(lat = lat, lng = lng, time = time)
    vcov_conley_hac_setup$arg_main = ""
    vcov_conley_hac_setup$patterns = c("", "lat + lng", "time", "lat + lng + time")

    #
    # Saving all the vcov possibilities
    #

    all_vcov = list(vcov_iid_setup,
                    vcov_hetero_setup,
                    vcov_clust_setup,
                    vcov_twoway_setup,
                    vcov_threeway_setup,
                    vcov_fourway_setup,
                    vcov_hac_setup,
                    vcov_conley_setup,
                    vcov_conley_hac_setup)

    options(fixest_vcov_builtin = all_vcov)

}

vcov_iid_internal = function(bread, scores, vars, meat_only){
    bread
}

vcov_hetero_internal = function(bread, scores, vars, meat_only){

    # we make a n/(n-1) adjustment to match vcovHC(type = "HC1")
    if(meat_only){
        meat = cpppar_crossprod(scores, 1, nthreads)
        return(meat)

    } else {
        # vcov = cpppar_crossprod(cpppar_matprod(scores, bread, nthreads), 1, nthreads) * correction.dof * ifelse(is_cluster, n/(n-1), 1)
        vcov_noAdj = cpppar_crossprod(cpppar_matprod(scores, bread, nthreads), 1, nthreads)
    }

    vcov_noAdj
}

# vcov_xx_internal arguments
# - bread: typically the VCOV_raw
# - scores
# - vars: a *list* of variables all of the same length which matches the scores
#
vcov_cluster_internal = function(bread, scores, vars){

    nway = length(vars)



}

####
#### SET / GET ========== ####
####



#' @rdname dof
#'
#' @param dof.type An object of class \code{dof.type} obtained with the function \code{\link[fixest]{dof}}.
setFixest_dof = function(dof.type = dof()){

    if(!"dof.type" %in% class(dof.type)){
        stop("The argument 'dof.type' must be an object created by the function dof().")
    }

    options("fixest_dof" = dof.type)
}

#' @rdname dof
"getFixest_dof"

getFixest_dof = function(){

    dof.type = getOption("fixest_dof")
    if(!"dof.type" %in% class(dof.type)){
        stop("The value of getOption(\"fixest_dof\") is currently not legal. Please use function setFixest_dict to set it to an appropriate value.")
    }

    dof.type
}


#' Sets the default type of standard errors to be used
#'
#' This functions defines or extracts the default type of standard-errors to computed in \code{fixest} \code{\link[fixest:summary.fixest]{summary}}, and \code{\link[fixest:vcov.fixest]{vcov}}.
#'
#' @param no_FE Character scalar equal to either: \code{"standard"} (default), or \code{"hetero"}. The type of standard-errors to use by default for estimations without fixed-effects.
#' @param one_FE Character scalar equal to either: \code{"standard"}, \code{"hetero"}, or \code{"cluster"} (default). The type of standard-errors to use by default for estimations with \emph{one} fixed-effect.
#' @param two_FE Character scalar equal to either: \code{"standard"}, \code{"hetero"}, \code{"cluster"} (default), or \code{"twoway"}. The type of standard-errors to use by default for estimations with \emph{two or more} fixed-effects.
#' @param all Character scalar equal to either: \code{"standard"}, or \code{"hetero"}. By default is is NULL. If provided, it sets all the SEs to that value.
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
#'
#' data(base_did)
#' est_no_FE  = feols(y ~ x1, base_did)
#' est_one_FE = feols(y ~ x1 | id, base_did)
#' est_two_FE = feols(y ~ x1 | id + period, base_did)
#'
#' etable(est_no_FE, est_one_FE, est_two_FE)
#'
#' # Changing the default standard-errors
#' setFixest_vcov(no_FE = "hetero", one_FE = "standard", two_FE = "twoway")
#' etable(est_no_FE, est_one_FE, est_two_FE)
#'
#' # Resetting the defaults
#' setFixest_vcov(reset = TRUE)
#'
#'
setFixest_vcov = function(no_FE = "standard", one_FE = "cluster", two_FE = "cluster", panel = "cluster", all = NULL, reset = FALSE){

    check_arg_plus(no_FE,  "match(standard, hetero)")
    check_arg_plus(one_FE, "match(standard, hetero, cluster)")
    check_arg_plus(two_FE, "match(standard, hetero, cluster, twoway)")
    check_arg_plus(panel, "match(standard, hetero, cluster, twoway, hac)")
    check_arg_plus(all,  "NULL match(standard, hetero)")
    check_arg_plus(reset, "logical scalar")

    opts = getOption("fixest_vcov_default")
    if(is.null(opts) || !is.list(opts) || reset){
        opts = list(no_FE = "standard", one_FE = "cluster", two_FE = "cluster", panel = "cluster")
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
        vcov_default = list(no_FE = "standard", one_FE = "cluster", two_FE = "cluster", panel = "cluster")
        options(fixest_vcov_default = vcov_default)
        return(vcov_default)
    }

    vcov_default
}






