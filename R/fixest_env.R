#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Jun 10 23:46:47 2019
# Purpose: sets up the params + performs
#           all necessary checks
#----------------------------------------------#


fixest_env = function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml = NULL,
                       fixef, NL.start, lower, upper, NL.start.init,
                       offset = NULL, subset, split = NULL, fsplit = NULL, linear.start = 0, jacobian.method = "simple",
                       useHessian = TRUE, hessian.args = NULL, opt.control = list(),
                       cluster, se, dof,
                       y, X, fixef_df, panel.id, fixef.rm = "perfect",
                       nthreads = getFixest_nthreads(), lean = FALSE,
                       verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000, collin.tol = 1e-14,
                       deriv.iter = 5000, deriv.tol = 1e-4, glm.iter = 25, glm.tol = 1e-8,
                       etastart, mustart,
                       warn = TRUE, notes = getFixest_notes(), combine.quick, demeaned = FALSE,
                       origin_bis, origin = "feNmlm", mc_origin, mc_origin_bis, mc_origin_ter,
                       computeModel0 = FALSE, weights = NULL,
                       debug = FALSE, mem.clean = FALSE, call_env = NULL, call_env_bis, ...){

    # INTERNAL function:
    # the estimation functions need input data in the exact format without any mistake possible (bc of c++)
    # this function takes care of implementing all the checks needed and providing the proper formatting
    # it also prepares the list returned in the est. funs

    # Only an environment is returned
    env = new.env(parent = emptyenv())

    # The function of origin
    if(!missnull(origin_bis)){
        origin = origin_bis
    }

    # The match.call of origin
    if(!missnull(mc_origin_bis)){
        mc_origin = mc_origin_bis
    }

    if(!missnull(call_env_bis)){
        call_env = call_env_bis
    }

    # NOTE on argument "origin":
    # I cannot reliably get the origin from deparse(mc_origin[[1]]) because the user might do
    # myfun = femlm // myfun(stuff) => then myfun would be the origin // I NEED this argument....

    if(origin %in% c("fenegbin", "femlm")){
        origin_type = "feNmlm"
    } else if(origin == "fepois") {
        origin_type = "feglm"
    } else {
        origin_type = origin
    }

    isFit = FALSE
    if(grepl("\\.fit", origin_type)){
        origin_type = gsub("\\.fit", "", origin_type)
        isFit = TRUE
    }

    #
    # Arguments control
    main_args = c("fml", "data", "panel.id", "offset", "subset", "split", "fsplit", "cluster", "se", "dof", "fixef.rm", "fixef.tol", "fixef.iter", "fixef", "nthreads", "lean", "verbose", "warn", "notes", "combine.quick", "start", "only.env", "mem.clean")
    femlm_args = c("family", "theta.init", "linear.start", "opt.control", "deriv.tol", "deriv.iter")
    feNmlm_args = c("NL.fml", "NL.start", "lower", "upper", "NL.start.init", "jacobian.method", "useHessian", "hessian.args")
    feglm_args = c("family", "weights", "glm.iter", "glm.tol", "etastart", "mustart", "collin.tol")
    feols_args = c("weights", "demeaned", "collin.tol")
    internal_args = c("debug")

    deprec_old_new = c()

    common_args = c(main_args, internal_args, deprec_old_new)

    # FIT methods
    feglm.fit_args = c(setdiff(common_args, c("fml", "data")), feglm_args, "y", "X", "fixef_df")
    feols.fit_args = c(setdiff(common_args, c("fml", "data")), feols_args, "y", "X", "fixef_df")

    args = names(mc_origin)
    args = args[nchar(args) > 0]

    # Checking the args
    valid_args = switch(origin,
                        feols = c(common_args, feols_args),
                        feglm = c(common_args, feglm_args),
                        feglm.fit = feglm.fit_args,
                        feols.fit = feols.fit_args,
                        femlm = c(common_args, femlm_args),
                        fenegbin = c(common_args, femlm_args),
                        fepois = c(common_args, feglm_args),
                        feNmlm = c(common_args, femlm_args, feNmlm_args))

    # stop if nonlinear args in linear funs (bc very likely it's a misunderstanding)
    if(origin != "feNmlm" && any(feNmlm_args %in% args)){
        qui_pblm = intersect(feNmlm_args, args)
        stop("Argument", enumerate_items(qui_pblm, "s.is"), " not valid for function ", origin, ". These are arguments of function feNmlm only.")
    }

    args_invalid = setdiff(args, valid_args)
    if(length(args_invalid) > 0){
        if(warn) warning(substr(deparse(mc_origin)[1], 1, 15), "...: ", enumerate_items(args_invalid, "is"), " not ", ifsingle(args_invalid, "a valid argument", "valid arguments"), " for function ", origin, ".", call. = FALSE)
    }

    args_deprec = deprec_old_new[deprec_old_new %in% args]
    if(length(args_deprec) > 0){
        dots = list(...)
        for(i in seq_along(args_deprec)){
            # we assign the deprecated argument to the new one if not provided
            old = args_deprec[i]
            new = names(args_deprec)[i]
            warning("Argument '", old, "' is deprecated. Use '", new, "' instead.", immediate. = TRUE)
            if(!new %in% args){
                assign(new, dots[[old]])
            }
        }

    }

    #
    # Internal quirks
    #

    if(origin_type %in% "feols"){
        family = "gaussian"
    }

    if(debug){
        verbose = 100
        assign("verbose", 100, env)
    } else if(!isScalar(verbose)){
        stop("Argument verbose must be an integer scalar.")
    } else {
        assign("verbose", verbose, env)
    }

    if(origin_type == "feglm"){
        # starting strategy
        start_provided = c(!is.null(linear.start), !is.null(etastart), !is.null(mustart))
        names(start_provided) = c("start", "etastart", "mustart")
        if(is.null(linear.start)) linear.start = 0
    }

    #
    # Formatting + checks
    #

    if(debug) cat(" - Formatting + checks\n")

    # we check the family => only for femlm/feNmlm and feglm
    # PROBLEM: the argument family has the same name in femlm and feglm but different meanings
    if(origin_type == "feNmlm"){
        family_name = try(match.arg(family), silent = TRUE)
        if("try-error" %in% class(family_name)){
            #  then we try with deparse
            family_dep = deparse(mc_origin$family)

            if(length(family_dep) > 1){
                stop("Argument family must be equal to 'poisson', 'logit', 'negbin' or 'gaussian'.")
            }

            family_name = try(match.arg(family_dep, c("poisson", "negbin", "logit", "gaussian")), silent = TRUE)
            if("try-error" %in% class(family_name)){
                stop("Argument family must be equal to 'poisson', 'logit', 'negbin' or 'gaussian'.")
            }
        }
        family = family_name

    } else if(origin_type == "feglm"){
        # We construct the family function, with some more variables

        # Family handling
        if(is.character(family)){
            family = error_sender(get(family, mode = "function", envir = parent.frame(2)), "Problem in the argument family:\n")
        }

        if(is.function(family)) {
            family = family()
        }

        if(is.null(family$family)) {
            stop("Argument 'family': the family is not recognized.")
        }

        family_type = switch(family$family,
                             poisson = "poisson",
                             quasipoisson = "poisson",
                             binomial = "logit",
                             quasibinomial = "logit",
                             "gaussian")

        if(family_type == "poisson" && family$link != "log") family_type = "gaussian"

        family_equiv = "none"
        if(family$family == "poisson" && family$link == "log") family_equiv = "poisson"
        if(family$family == "binomial" && family$link == "logit") family_equiv = "logit"

        family$family_type = family_type
        family$family_equiv = family_equiv

        #
        # Custom functions (we will complete these functions later in this code, after lhs and nthreads are set)
        # they concern only logit and poisson (so far)
        #

        dev.resids = family$dev.resids
        family$sum_dev.resids = function(y, mu, eta, wt) sum(dev.resids(y, mu, wt))

        fun_mu.eta = family$mu.eta
        family$mu.eta = function(mu, eta) fun_mu.eta(eta)

        if(is.null(family$valideta)) family$valideta = function(...) TRUE
        if(is.null(family$validmu)) family$validmu = function(...) TRUE

        # QUIRKINESS => now family becomes the family name and the functions become family_funs
        family_funs = family
        family = family_type

        computeModel0 = family_equiv %in% c("poisson", "logit")
    }

    check_arg(demeaned, notes, warn, mem.clean, "logical scalar")

    check_arg_plus(fixef.rm, "match(singleton, perfect, both, none)")

    check_arg(collin.tol, "numeric scalar GT{0}")

    show_notes = notes
    notes = c()

    # flags for NA infinite vales => will be used in the message (good to discrimintae b/w NA and inf)
    ANY_INF = FALSE
    ANY_NA = FALSE
    anyNA_sample = FALSE
    isNA_sample = FALSE
    message_NA = ""

    # summary information
    do_summary = FALSE

    # Note on mem.clean
    # we remove the most intermediary objects as possible
    # we apply gc only from time to time since it is costly
    # flag gc triggered to monitor when gc was last triggered => avoids triggering it too often
    gc2trig = TRUE

    #
    # nthreads argument
    nthreads = check_set_nthreads(nthreads)

    # The family functions (for femlm only)
    famFuns = switch(family,
                     poisson = ml_poisson(),
                     negbin = ml_negbin(),
                     logit = ml_logit(),
                     gaussian = ml_gaussian())


    #
    # Formula handling ####
    #

    fml_no_xpd = NULL # will be returned if expansion is performed
    isPanel = FALSE
    do_iv = FALSE
    multi_fixef = FALSE
    fake_intercept = FALSE
    fml_fixef = fml_iv = NULL
    fixef_terms_full = NULL
    complete_vars = c()
    if(isFit){
        isFixef = !missnull(fixef_df)
    } else {

        #
        # The data

        if(missing(data)) stop("You must provide the argument 'data' (currently it is missing).")
        check_value(data, "matrix | data.frame", .arg_name = "data")

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

        #
        # The fml => controls + setup
        if(missing(fml)) stop("You must provide the argument 'fml' (currently it is missing).")
        check_arg(fml, "ts formula")

        fml = formula(fml) # we regularize the formula to check it

        # We apply expand for macros => we return fml_no_xpd
        if(length(getFixest_fml()) > 0 || ".." %in% all.vars(fml, functions = TRUE)){
            fml_no_xpd = fml
            fml = .xpd(fml, data = dataNames, macro = TRUE)
        }

        #
        # Checking the validity of the formula
        #

        fml_parts = fml_split(fml, raw = TRUE)
        n_parts = length(fml_parts)

        # checking iv
        if(n_parts > 2){
            if(!origin_type == "feols"){
                stop("The argument 'fml' cannot contain more than two parts separated by a pipe ('|'). IVs are available only for 'feols'.\nThe syntax is: DEP VAR ~ EXPL VARS | FIXED EFFECTS. (No IVs allowed.)")
            } else {
                if(n_parts > 3){
                    stop("In feols, the argument 'fml' cannot contain more than three parts separated by a pipe ('|').\nThe syntax is: DEP VAR ~ EXPL VARS | FIXED EFFECTS | IV FORMULA.")
                }

                if(!is_fml_inside(fml_parts[[3]])){
                    stop("in feols, the third part of the formula must be the IV formula, and it does not look like a formula right now.")
                }

                if(is_fml_inside(fml_parts[[2]])){
                    stop("In feols, problem in the formula: the RHS must contain at most one formula (the IV formula), currently there are two formulas.")
                }

                do_iv = TRUE
            }
        } else if(n_parts == 2){
            if(is_fml_inside(fml_parts[[2]])){
                # 2nd part a formula only allowed in feols
                if(origin_type == "feols"){
                    do_iv = TRUE
                } else {
                    stop("The RHS of the formula must represent the fixed-effects (and can't be equal to a formula: only feols supports this).")
                }
            }
        }

        # Checking the presence
        complete_vars = all_vars_with_i_prefix(fml)
        if(any(!complete_vars %in% c(dataNames, ".F", ".L"))){
            # Question: is the missing variable a scalar from the global environment?
            var_pblm = setdiff(complete_vars, dataNames)
            n_pblm = length(var_pblm)
            var_pblm_dp = character(n_pblm)
            type = c()
            ok = TRUE
            for(i in 1:n_pblm){
                var = var_pblm[i]
                if(exists(var, envir = call_env)){
                    if(length(var) < 5){
                        var_pblm_dp = deparse_long(eval(str2lang(var), call_env))
                        type[i] = "scalar"
                    } else {
                        ok = FALSE
                        type[i] = if(length(var) == nrow(data)) "var" else "other"
                    }
                } else {
                    ok = FALSE
                    type[i] = "other"
                }
            }

            if(ok){
                # we rewrite the formula with hard values
                fml_dp = deparse_long(fml)
                for(i in 1:n_pblm){
                    pattern = paste0("(^|[^[:alnum:]._])\\Q", var_pblm[i], "\\E($|[^[:alnum:]._])")
                    sub = paste0("\\1", var_pblm_dp[i], "\\2")
                    fml_dp = gsub(pattern, sub, fml_dp, perl = TRUE)
                }

                fml = str2lang(fml_dp)
                fml_parts = fml_split(fml, raw = TRUE)

            } else {

                var_pblm = var_pblm[type != "scalar"]
                type = type[type != "scalar"]

                # Error => we take the time to provide an informative error message
                LHS = all.vars(fml_parts[[1]][[2]])
                RHS = all.vars(fml_parts[[1]][[3]])

                msg_builder = function(var_pblm, type, qui){
                    msg = ""
                    if(any(type[qui] == "var")){
                        msg = " Note that fixest does not accept variables from the global enviroment, they must be in the data set"
                        extra = "."
                        if(!all(type[qui] == "var")){
                            extra = paste0(" (it concerns ", enumerate_items(var_pblm[qui][type[qui] == "var"]), ").")
                        }
                        msg = paste0(msg, extra)
                    }

                    msg
                }

                if(any(!LHS %in% dataNames)){
                    qui = which(var_pblm %in% LHS)
                    msg = msg_builder(var_pblm, type, qui)

                    stop("The variable", enumerate_items(var_pblm[qui], "s.is.quote"), " in the LHS of the formula but not in the data set.", msg)
                }

                if(any(!RHS %in% dataNames)){
                    qui = which(var_pblm %in% RHS)
                    msg = msg_builder(var_pblm, type, qui)

                    stop("The variable", enumerate_items(var_pblm[qui], "s.is.quote"), " in the RHS", ifunit(n_parts, "", " (first part)"), " of the formula but not in the data set.", msg)
                }

                part_2 = all.vars(fml_parts[[2]])
                if(any(!part_2 %in% dataNames)){
                    qui = which(var_pblm %in% part_2)
                    msg_end = msg_builder(var_pblm, type, qui)

                    msg = ifelse(is_fml_inside(fml_parts[[2]]), "IV", "fixed-effects")
                    stop("The variable", enumerate_items(var_pblm[qui], "s.is.quote"), " in the ", msg, " part of the formula but not in the data set.", msg_end)
                }

                part_3 = all.vars(fml_parts[[3]])
                qui = which(var_pblm %in% part_3)
                msg = msg_builder(var_pblm, type, qui)
                stop("The variable", enumerate_items(var_pblm[qui], "s.is.quote"), " in the IV part of the formula but not in the data set.", msg)
            }

        }


        #
        # ... Panel setup ####
        #

        if(debug) cat(" ---> Panel setup\n")

        info = fixest_fml_rewriter(fml)
        isPanel = info$isPanel
        fml = info$fml

        if(isPanel){
            panel.info = NULL
            if(!is.null(attr(data, "panel_info"))){
                if(!missnull(panel.id)){
                    warning("The argument 'panel.id' is provided but argument 'data' is already a 'fixest_panel' object. Thus the argument 'panel.id' is ignored.", immediate. = TRUE)
                }

                panel__meta__info = attr(data, "panel_info")
                panel.id = panel__meta__info$panel.id
                panel.info = panel__meta__info$call
            } else {
                # Later: automatic deduction using the first two clusters
                if(missnull(panel.id)){
                    stop("To use lag/leads (with l()/f()): either provide the argument 'panel.id' with the panel identifiers OR set your data as a panel with function panel().")
                }
                panel__meta__info = panel_setup(data, panel.id, from_fixest = TRUE)
            }
            class(data) = "data.frame"
        }

        fml_parts = fml_split(fml, raw = TRUE)

        # The different parts of the formula
        fml_fixef = fml_iv = NULL
        if(n_parts == 3){
            fml_fixef = fml_maker(fml_parts[[2]])
            fml_iv = fml_maker(fml_parts[[3]])
        } else if(n_parts == 2){
            if(do_iv){
                fml_iv = fml_maker(fml_parts[[2]])
            } else {
                fml_fixef = fml_maker(fml_parts[[2]])
            }
        }

        fml_linear = fml_maker(fml_parts[[1]])

        #
        # ... Fixed-effects ####
        #

        if(debug) cat(" ---> Fixed effects\n")

        # for clarity, arg fixef is transformed into fml_fixef
        if(!is.null(fml_fixef) && !missnull(fixef)){
            stop("To add fixed-effects: either include them in argument 'fml' using a pipe ('|'), either use the argument 'fixef'. You cannot use both!")

        } else if(!missing(fixef) && length(fixef) >= 1){
            # We check the argument => character vector
            check_arg_plus(fixef, "multi match", .choices = dataNames, .message = "Argument 'fixef', when provided, must be a character vector of variable names.")

            # we transform it into a formula
            fml_fixef = as.formula(paste0("~", paste(fixef, collapse = " + ")))
        }

        isFixef = !is.null(fml_fixef)

        # we check the fixed-effects
        if(isFixef){
            # => extraction of the stepwise information
            fixef_info_stepwise = error_sender(fixef_terms(fml_fixef, stepwise = TRUE, origin_type = origin_type),
                                         "The fixed-effects part of the formula (", charShorten(as.character(fml_fixef)[2], 15),
                                         ") is not valid:\n", clean = "_impossible_var_name_ => ^")

            if(fixef_info_stepwise$do_multi){
                fml_fixef = fml_tl = fixef_info_stepwise$fml
                multi_fixef = TRUE
                # We will compute the fixed-effects afterwards
                isFixef = FALSE

            } else {
                # The function already returns the terms, so no need to recall it
                fml_tl = fixef_info_stepwise$tl
                if(length(fml_tl) == 0) isFixef = FALSE
            }


            fixef_terms_full = error_sender(fixef_terms(fml_tl),
                                            "The fixed-effects part of the formula (", charShorten(as.character(fml_fixef)[2], 15),
                                            ") is not valid:\n", clean = "_impossible_var_name_ => ^")

            # if fixed-effects are provided, we make sure there is an
            # intercept so that factors can be handled properly
            #
            # Beware: only slope: we don't add it
            if(length(fixef_terms_full$slope_flag) == 0){
                # We can be here when SW and the first SW is no fixef
                # I redefine it just for sake of clarity
                fake_intercept = FALSE
            } else {
                fake_intercept = any(fixef_terms_full$slope_flag >= 0)
            }

        }

    }

    # We need to get the variables from the argument cluster and add them to fml_full
    # so that subset works properly
    if(!missnull(cluster)){

        check_arg(cluster, "os formula | character vector")

        if(is.character(cluster)){
            cluster = .xpd(rhs = cluster)
        }

        cluster_origin = cluster

        # Check
        clust_vars = all.vars(cluster)
        if(any(!clust_vars %in% dataNames)){
            stop("In argument 'cluster' the variable", enumerate_items(setdiff(clust_vars, dataNames), "s.is.quote"), " is not in the data set. It must be composed of variable names only.")
        }

        complete_vars = unique(c(complete_vars, clust_vars))
    }

    #
    # Data formatting ####
    #

    if(debug) cat(" - Data formattiing\n")


    #
    # ... subset ####
    #

    if(debug) cat(" ---> subset\n")

    obs_selection = list()
    multi_lhs = FALSE
    nobs_origin = NULL
    if(isFit){
        # We have to first eval y if isFit in order to check subset properly
        if(missing(y)){
            stop("You must provide argument 'y' when using ", origin, ".")
        }

        lhs = check_value_plus(y, "numeric vmatrix conv | data.frame")
        lhs_names = NULL
        if(is.data.frame(lhs)){

            lhs_names = names(lhs)

            if(ncol(lhs) > 1){
                multi_lhs = TRUE

                lhs_new = list()
                for(i in 1:ncol(lhs)){
                    lhs_current = lhs[[i]]
                    if(!is.numeric(lhs_current) || !is.logical(lhs_current)){
                        stop("The ", n_th(i), " column of argument 'y' is not numeric. Estimation cannot be done.")
                    }

                    if(is.logical(lhs_current) || is.integer(lhs_current)){
                        lhs_new[[i]] = as.numeric(lhs_current)
                    } else {
                        lhs_new[[i]] = lhs_current
                    }
                }

                lhs = lhs_new
                if(mem.clean) rm(lhs_new)

            } else {
                if(!is.numeric(lhs[[1]]) || !is.logical(lhs[[1]])){
                    stop("The argument 'y' is not numeric. Estimation cannot be done.")
                }

                lhs = lhs[[1]]
            }
        } else if(is.matrix(lhs)){

            lhs_names = colnames(lhs)

            if(ncol(lhs) > 1){
                multi_lhs = TRUE

                lhs_new = list()
                for(i in 1:ncol(lhs)){
                    lhs_new[[i]] = lhs[, i]
                }

                lhs = lhs_new
                if(mem.clean) rm(lhs_new)
            } else {
                lhs = as.vector(lhs)
            }
        }

        n_lhs = if(multi_lhs) length(lhs) else 1
        if(is.null(lhs_names)){

            lhs_names_raw = deparse_long(mc_origin[["y"]])

            if(n_lhs > 1){
                lhs_names = paste0(lhs_names_raw, 1:n_lhs)
            } else {
                lhs_names = lhs_names_raw
            }
        }

        # we reconstruct a formula
        fml_linear = .xpd(lhs = lhs_names, rhs = 1)

        nobs_origin = nobs = if(multi_lhs) length(lhs[[1]]) else length(lhs)
    }

    # delayed.subset only concerns isFit // subsetting must always occur before NA checking
    isSubset = FALSE
    delayed.subset = FALSE
    if(!missing(subset)){
        if(!is.null(subset)){

            isSubset = TRUE

            if("formula" %in% class(subset)){

                if(isFit){
                    stop("In ", origin, " the subset cannot be a formula. You must provide either an integer vector or a logical vector.")
                }

                check_value(subset, "os formula var(data)", .data = data)
                subset.value = subset[[2]]
                subset = check_value_plus(subset.value, "evalset integer vmatrix ncol(1) gt{0} | logical vector len(data)", .data = data, .prefix = "In argument 'subset', the expression")

            } else {

                if(isFit){
                    check_value(subset, "integer vmatrix ncol(1) gt{0} | logical vector len(data)", .data = lhs)
                } else {
                    check_value(subset, "integer vmatrix ncol(1) gt{0} | logical vector len(data)",
                                .prefix = "If not a formula, argument 'subset'", .data = data)
                }

            }

            if(is.matrix(subset)) subset = as.vector(subset)

            isNA_subset = is.na(subset)
            if(anyNA(isNA_subset)){
                if(all(isNA_subset)){
                    stop("In argument 'subset', all values are NA, estimation cannot be done.")
                }

                if(is.logical(subset)){
                    subset[isNA_subset] = FALSE
                } else {
                    subset = subset[!isNA_subset]
                }
            }

            if(is.logical(subset)){
                subset = which(subset)
            }

            if(isFit){
                delayed.subset = TRUE
                lhs = lhs[subset]
            } else {
                nobs_origin = NROW(data)

                # subsetting creates a deep copy. We avoid copying the entire data set.
                # Note that we don't need to check that complete_vars exists in the data set
                # that will be done in the dedicated sections.

                if(missnull(offset)) offset = NULL
                if(missnull(weights)) weights = NULL
                if(missnull(split)) split = NULL
                if(missnull(fsplit)) fsplit = NULL
                if(missnull(NL.fml)) NL.fml = NULL

                additional_vars = collect_vars(NL.fml, offset, weights, split, fsplit)
                complete_vars = c(complete_vars, additional_vars)

                complete_vars = intersect(unique(complete_vars), names(data))

                data = data[subset, complete_vars, drop = FALSE]
            }

            # We add subset to the obs selected
            obs_selection = list(subset = subset)

        } else if(!is.null(mc_origin$subset)){
            dp = deparse_long(mc_origin$subset)
            if((grepl("[[", dp, fixed = TRUE) || grepl("$", dp, fixed = TRUE)) && dp != 'x[["subset"]]'){
                # we avoid this behavior
                stop("Argument 'subset' (", dp, ") is evaluated to NULL. This is likely not what you want.")
            }
        }
    }


    #
    # ... The left hand side ####
    #

    if(debug) cat(" ---> LHS\n")

    # the LHS for isFit has been done just before subset

    # evaluation
    if(isFit == FALSE){

        # The LHS must contain only values in the DF
        namesLHS = all.vars(fml_linear[[2]])
        lhs_text = deparse_long(fml_linear[[2]])
        if(length(namesLHS) == 0){
            stop("The right hand side of the formula (", lhs_text, ") contains no variable!")

        }

        lhs_text2eval = gsub("^(c|(c?sw0?))\\(", "list(", lhs_text)

        lhs = error_sender(eval(str2lang(lhs_text2eval), data),
                           "Evaluation of the left-hand-side (equal to ", lhs_text, ") raises an error: \n")

        if(is.list(lhs)){
            # we check the consistency
            lhs_names = eval(str2lang(gsub("^list\\(", "sw(", lhs_text2eval)))

            n_all = lengths(lhs)
            if(!all(n_all == n_all[1])){
                i = which(n_all != n_all[1])[1]
                stop("The evaluation of the multiple left hand sides raises an error: the first LHS is of length ", n_all[1], " while the ", n_th(i), " is of length ", n_all[i], ". They should all be of the same length")
            }

            # individual check + conversions
            for(i in seq_along(lhs)){
                lhs[[i]] = check_value_plus(lhs[[i]], "numeric vector conv", .prefix = paste0("Problem in the ", n_th(i), " left hand side. It"))
            }

            if(length(lhs) == 1){
                lhs = lhs[[1]]
            } else {
                multi_lhs = TRUE
            }

        } else {
            lhs = check_value_plus(lhs, "numeric vmatrix ncol(1) conv", .prefix = "The left hand side")
            if(is.matrix(lhs)) lhs = as.vector(lhs)
            lhs_names = lhs_text
        }

        if(is.list(lhs)){
            nobs = length(lhs[[1]])
        } else {
            nobs = length(lhs)
        }
    }

    if(is.null(nobs_origin)) {
        nobs_origin = nobs
    }

    check_LHS_const = FALSE
    anyNA_y = FALSE
    msgNA_y = ""
    if(multi_lhs){

        if(family %in% c("poisson", "negbin")){
            for(i in seq_along(lhs)){
                if(any(lhs[[i]] < 0, na.rm = TRUE)){
                    stop("Presence of negative values in the ", n_th(i), " dependent variable: this is not allowed for the \"", family, "\" family.")
                }
            }
        }

        if(origin_type == "feNmlm" && family %in% "logit"){
            if(!all(lhs[[i]]==0 | lhs[[i]]==1, na.rm = TRUE)){
                stop("The ", n_th(i), " dependent variable has values different from 0 or 1.\nThis is not allowed for the \"logit\" family.")
            }

        }
    } else {

        lhs_clean = lhs # copy used for NA case
        info = cpppar_which_na_inf_vec(lhs, nthreads)
        if(info$any_na_inf){

            if(info$any_na) ANY_NA = TRUE
            if(info$any_inf) ANY_INF = TRUE

            anyNA_y = TRUE
            isNA_y = info$is_na_inf
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_y
            msgNA_y = paste0("LHS: ", numberFormatNormal(sum(isNA_y)))

            lhs_clean = lhs[!isNA_y]

            if(mem.clean){
                rm(isNA_y)
            }
        }

        if(mem.clean){
            rm(info)
        }

        # we check the var is not a constant
        if(cpp_isConstant(lhs_clean)){
            # We delay this evaluation
            check_LHS_const = TRUE
        }

        if(family %in% c("poisson", "negbin") && any(lhs_clean < 0)){
            stop("Negative values of the dependent variable are not allowed for the \"", family, "\" family.")
        }

        if(origin_type == "feNmlm" && family %in% "logit" && !all(lhs_clean==0 | lhs_clean==1)){
            stop("The dependent variable has values different from 0 or 1.\nThis is not allowed for the \"logit\" family.")
        }
    }


    #
    # ... linear part ####
    #

    if(debug) cat(" ---> Linear part\n")

    interaction.info = NULL
    agg = NULL
    # This variables will be set globally from within fixest_model_matrix!
    GLOBAL_fixest_mm_info = list()
    multi_rhs = FALSE
    if(isFit){

        isLinear = FALSE
        if(missing(X) || is.null(X)){
            if(!isFixef){
                stop("Argument X must be provided in the absence of fixed-effects.")
            } else {
                isLinear = FALSE
            }
        } else {
            isLinear = TRUE
        }

        if(isLinear){

            if("data.frame" %in% class(X)){
                X = as.matrix(X)
            }

            if(isVector(X)){
                if(length(X) == 1){
                    X = matrix(X, nrow = nobs, ncol = 1)
                } else {
                    X = matrix(X, ncol = 1)
                }

                colnames(X) = deparse(mc_origin[["X"]])[1]
            }

            if(!is.matrix(X) || !is.numeric(X)){
                stop("Argument X must be a matrix or a vector (in case of only one variable).")
            }

            if(nrow(X) != nobs){
                stop("The number of rows of X (", nrow(X), ") must be of the same dimension as y (", nobs, ").")
            }

            linear.mat = X

            # we coerce to matrix (avoids sparse matrix) + transform into double
            if(!is.matrix(linear.mat) || is.integer(linear.mat[1])){
                linear.mat = as.matrix(X) * 1
            }

            if(is.null(colnames(linear.mat))){
                colnames(linear.mat) = paste0("X", 1:ncol(linear.mat))
            }

            # The formula
            rhs = as_varname(colnames(linear.mat))
            fml_linear = .xpd(lhs = fml_linear[[2]], rhs = rhs)

            linear.varnames = NULL

            if(delayed.subset){
                linear.mat = linear.mat[subset, , drop = FALSE]
            }

        }

    } else {
        isLinear = FALSE

        linear.varnames = all_vars_with_i_prefix(fml_linear[[3]])

        fml_terms = terms(fml_linear)
        if(length(linear.varnames) > 0 || attr(fml_terms, "intercept") == 1){
            isLinear = TRUE
        }

        if(isLinear){

            #
            # Handling multiple RHS
            #

            rhs_info_stepwise = error_sender(extract_stepwise(fml_linear), "Problem in the RHS of the formula: ")

            multi_rhs = rhs_info_stepwise$do_multi


            #
            # We construct the linear matrix
            #

            if(multi_rhs){
                # We construct:
                # - left and right cores (vars that are always there)
                # - center

                isLinear = FALSE # => superseded by the multi_rhs mechanism

                fml_core_left = rhs_info_stepwise$fml_core_left
                fml_core_right = rhs_info_stepwise$fml_core_right
                fml_all_sw = rhs_info_stepwise$fml_all_sw

                linear_core = list()
                linear_core$left = error_sender(fixest_model_matrix(fml_core_left, data, fake_intercept),
                                                "Evaluation of the right-hand-side of the formula raises an error: ")

                linear_core$right = error_sender(fixest_model_matrix(fml_core_right, data, TRUE),
                                                "Evaluation of the right-hand-side of the formula raises an error: ")

                rhs_sw = list()
                for(i in seq_along(fml_all_sw)){
                    rhs_sw[[i]] = error_sender(fixest_model_matrix(fml_all_sw[[i]], data, TRUE),
                                               "Evaluation of the right-hand-side of the formula raises an error: ")
                }

            } else {
                # Regular, single RHS

                linear.mat = error_sender(fixest_model_matrix(fml_linear, data, fake_intercept),
                                          "Evaluation of the right-hand-side of the formula raises an error: ")

                if(identical(linear.mat, 1)){
                    isLinear = FALSE
                }

            }
        }

        if("sunab" %in% names(GLOBAL_fixest_mm_info)){
            agg = GLOBAL_fixest_mm_info$sunab$agg
            do_summary = TRUE
        }

    }

    if(check_LHS_const){
        # Estimation cannot be done if fixed-effects or the constant are present
        # Note that the fit stats get messed up in such estimation => but since it's super niche I don't correct for it

        if(isFixef){
            stop("The dependent variable is a constant. The estimation with fixed-effects cannot be done.")

        } else if(isLinear && (!isFit && attr(fml_terms, "intercept") == 1)){
            stop("The dependent variable is a constant. The estimation cannot be done. If you really want to carry on, please remove the intercept first.")
        }
    }

    # Further controls (includes na checking)
    msgNA_L = ""
    if(multi_rhs){

        anyNA_L = FALSE
        isNA_L = FALSE
        if(length(linear_core$left) > 1){
            info = cpppar_which_na_inf_mat(linear_core$left, nthreads)
            if(info$any_na_inf){
                anyNA_L = TRUE
                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE
                isNA_L = isNA_L | info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_L
            }
        }

        if(length(linear_core$right) > 1){
            info = cpppar_which_na_inf_mat(linear_core$right, nthreads)
            if(info$any_na_inf){
                anyNA_L = TRUE
                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE
                isNA_L = isNA_L | info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_L
            }
        }

        if(anyNA_L){
            msgNA_L = paste0("RHS: ", numberFormatNormal(sum(isNA_L)))
        }
    }

    if(isLinear){

        linear.params = colnames(linear.mat)
        anyNA_L = FALSE
        info = cpppar_which_na_inf_mat(linear.mat, nthreads)
        if(info$any_na_inf){
            anyNA_L = TRUE

            if(info$any_na) ANY_NA = TRUE
            if(info$any_inf) ANY_INF = TRUE

            isNA_L = info$is_na_inf
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_L
            msgNA_L = paste0("RHS: ", numberFormatNormal(sum(isNA_L)))

            if(mem.clean){
                rm(isNA_L)
            }

        }

        if(mem.clean){
            rm(info)
            gc()
            gc2trig = FALSE
        }

    } else if(multi_rhs && origin_type %in% "feNmlm"){
        # We need to assign linear.params
        v_core = unlist(lapply(linear_core, colnames))
        v_sw = unlist(lapply(rhs_sw, colnames))
        linear.varnames = linear.params = unique(c(v_core, v_sw))

    } else {
        linear.params = linear.start = linear.varnames = NULL
    }


    #
    # ... nonlinear part ####
    #

    if(debug) cat(" ---> Non-linear part\n")

    msgNA_NL = ""
    if(!missnull(NL.fml)){

        if(!"formula" %in% class(NL.fml)) stop("Argument 'NL.fml' must be a formula.")
        NL.fml = formula(NL.fml) # we regularize the formula

        isNonLinear = TRUE
        nl.call = NL.fml[[length(NL.fml)]]

        allnames = all.vars(nl.call)
        nonlinear.params = allnames[!allnames %in% dataNames]
        nonlinear.varnames = allnames[allnames %in% dataNames]

        if(length(nonlinear.params) == 0){
            warning("As there is no parameter to estimate in argument 'NL.fml', this argument is ignored. If you want to add an offset, use argument 'offset'.")
        }

        data_NL = data[, nonlinear.varnames, drop = FALSE]

        # Update 08-2020 => We now allow non-numeric variables in the NL part (they can be used as identifiers)
        anyNA_NL = FALSE
        if(any(quiNA <- sapply(data_NL, anyNA))){
            anyNA_NL = TRUE
            quiNA = which(quiNA)

            isNA_NL = is.na(data_NL[[quiNA[1]]])
            for(i in quiNA[-1]){
                isNA_NL = isNA_NL | is.na(data_NL[[quiNA[i]]])
            }

            ANY_NA = TRUE
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_NL
            msgNA_NL = paste0("NL: ", numberFormatNormal(sum(isNA_NL)))

            if(mem.clean){
                rm(isNA_NL)
                gc2trig = TRUE
            }

        }

    } else {
        isNonLinear = FALSE
        nl.call = 0
        allnames = nonlinear.params = nonlinear.varnames = character(0)
        NL.fml = NULL
    }

    params = c(nonlinear.params, linear.params)
    lparams = length(params)
    varnames = c(nonlinear.varnames, linear.varnames)

    # Attention les parametres non lineaires peuvent etre vides
    isNL = length(nonlinear.params) > 0

    #
    # ... Offset ####
    #

    if(debug) cat(" ---> offset\n")

    isOffset = FALSE
    offset.value = 0
    msgNA_offset = ""
    if(!missing(offset)){
        if(!is.null(offset)){
            isOffset = TRUE

            if("formula" %in% class(offset)){

                if(isFit){
                    stop("In ", origin, " the offset cannot be a formula. You must provide a numeric vector.")
                }

                check_value(offset, "os formula var(data)", .data = data)
                offset.value = offset[[2]]
                check_value_plus(offset.value, "evalset numeric vmatrix ncol(1) conv", .data = data, .prefix = "In argument 'offset', the expression")

            } else {

                check_value_plus(offset, "numeric vmatrix ncol(1) conv", .prefix = "If not a formula, argument 'offset'")

                if(length(offset) == 1){
                    offset.value = rep(offset, nobs)
                } else {
                    if(length(offset) != nobs_origin){
                        stop("The offset's length should be equal to the data's length (currently it's ", numberFormatNormal(length(offset)), " instead of ", numberFormatNormal(nobs_origin), ").")
                    } else {
                        offset.value = offset
                    }

                    if(length(obs_selection$subset) > 0){
                        offset.value = offset.value[obs_selection$subset]
                    }
                }
            }

            if(is.matrix(offset.value)) offset.value = as.vector(offset.value)

            if(delayed.subset){
                offset.value = offset.value[subset]
            }

            anyNA_offset = FALSE
            info = cpppar_which_na_inf_vec(offset.value, nthreads)
            if(info$any_na_inf){
                anyNA_offset = TRUE

                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE

                isNA_offset = info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_offset
                msgNA_offset = paste0("Offset: ", numberFormatNormal(sum(isNA_offset)))

                if(mem.clean){
                    rm(isNA_offset)
                }
            }

            if(mem.clean){
                rm(info)
                gc2trig = TRUE
            }

        } else if(is.null(offset)){
            # msg if it's not what the user wanted
            if(!is.null(mc_origin$offset)){
                dp = deparse_long(mc_origin$offset)

                if((grepl("[[", dp, fixed = TRUE) || grepl("$", dp, fixed = TRUE)) && dp != 'x$offset'){
                    # we avoid this behavior
                    stop("Argument 'offset' (", dp, ") is evaluated to NULL. This is likely not what you want.")
                }
            }
        }
    }

    #
    # ... Weights ####
    #

    if(debug) cat(" ---> weights\n")

    any0W = anyNA_W = FALSE
    msgNA_weight = message_0W = ""
    weights.value = 1
    isWeight = FALSE
    if(!missing(weights)){
        if(!is.null(weights)){
            isWeight = TRUE

            if("formula" %in% class(weights)){

                if(isFit){
                    stop("In ", origin, " the weights cannot be a formula. You must provide a numeric vector.")
                }

                check_value(weights, "os formula var(data)", .data = data)
                weights.value = weights[[2]]
                check_value_plus(weights.value, "evalset numeric vmatrix ncol(1) conv", .data = data, .prefix = "In argument 'weights', the expression")


            } else {

                check_value_plus(weights, "numeric vmatrix ncol(1) conv", .prefix = "If not a formula, argument 'weights'")

                if(length(weights) == 1){
                    if(weights == 1){
                        isWeight = FALSE
                        weights.value = 1
                    } else {
                        # No point in having all obs the same weight... yet...
                        weights.value = rep(weights, nobs)
                    }
                } else {
                    if(length(weights) != nobs_origin){
                        stop("The weights's length should be equal to the data's length (currently it's ", numberFormatNormal(length(weights)), " instead of ", numberFormatNormal(nobs_origin), ").")
                    } else {
                        weights.value = weights
                    }

                    if(length(obs_selection$subset) > 0){
                        weights.value = weights.value[obs_selection$subset]
                    }
                }
            }

            if(is.matrix(weights.value)) weights.value = as.vector(weights.value)

            if(delayed.subset){
                weights.value = weights.value[subset]
            }

            anyNA_weights = FALSE
            info = cpppar_which_na_inf_vec(weights.value, nthreads)
            if(info$any_na_inf){
                anyNA_weights = TRUE

                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE

                isNA_W = info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_W
                msgNA_weight = paste0("Weights: ", numberFormatNormal(sum(isNA_W)))

                if(mem.clean){
                    rm(isNA_W)
                }
            }

            if(mem.clean){
                rm(info)
                gc2trig = TRUE
            }

            # need to check NA if there are some (ie if NA in vector: any() gives TRUE if any TRUE, NA if no TRUE)
            anyNeg = any(weights.value < 0)
            if(!is.na(anyNeg) && anyNeg){
                stop("The vector of weights contains negative values, this is not allowed (e.g. obs. ", enumerate_items(head(which(weights.value < 0), 3)), ").")
            }

            # we remove 0 weights
            any0W = any(weights.value == 0)
            any0W = !is.na(any0W) && any0W
            if(any0W){
                is0W = weights.value == 0
                is0W = is0W & !is.na(is0W)
                message_0W = paste0(numberFormatNormal(sum(is0W)), " observation", plural(sum(is0W)), " removed because of 0-weight.")
                notes = c(notes, message_0W)
            }

        } else if(!is.null(mc_origin$weights)){
            dp = deparse_long(mc_origin$weights)
            if((grepl("[[", dp, fixed = TRUE) || grepl("$", dp, fixed = TRUE)) && dp != 'x$weights'){
                # we avoid this behavior
                stop("Argument 'weights' (", dp, ") is evaluated to NULL. This is likely not what you want.")
            }
        }
    }

    if(mem.clean && gc2trig){
        gc()
        gc2trig = FALSE
    }


    #
    # ... Split ####
    #

    if(debug) cat(" ---> split\n")

    isSplit = FALSE
    msgNA_split = ""
    split.full = FALSE
    if(!missnull(split) || !missnull(fsplit)){

        if(!missnull(fsplit)){
            if(!missnull(split)){
                stop("You cannot provide the two arguments 'split' and 'fsplit' at the same time.")
            }
            split.full = TRUE
            split = fsplit
        }

        if("formula" %in% class(split)){

            if(isFit){
                stop("In ", origin, " the split cannot be a formula. You must provide a numeric vector.")
            }

            check_value(split, "os formula var(data)", .data = data)
            split.value = split[[2]]
            split.name = as.character(split.value)
            split = check_value_plus(split.value, "evalset vector", .data = data, .prefix = "In argument 'split', the expression")

        } else {

            if(isFit){
                check_value(split, "vector len(value)", .value = nobs_origin)
            } else {
                check_value(split, "vector len(value) | character scalar", .value = nobs_origin, .prefix = "If not a formula, argument 'split'")
            }

            if(length(split) == 1){
                split.name = split
                check_value(split, "charin", .choices = dataNames, .message = "If equal to a character scalar, argument 'split' must be a variable name.")
                split = data[[split]]
            } else {
                split.name = gsub("^[[:alpha:]][[:alnum:]\\._]*\\$", "", deparse_long(mc_origin$split))

                if(length(obs_selection$subset) > 0){
                    split = split[obs_selection$subset]
                }
            }

        }

        if(delayed.subset){
            split = split[subset]
        }

        isSplit = TRUE
        # We will unclass split after NA removal

        anyNA_split = FALSE
        isNA_S = is.na(split)
        if(any(isNA_S)){
            anyNA_split = TRUE

            if(info$any_na) ANY_NA = TRUE

            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_S
            msgNA_split = paste0("split: ", numberFormatNormal(sum(isNA_S)))

            if(mem.clean){
                rm(isNA_S)
            }
        }

    }

    if(mem.clean && gc2trig){
        gc()
        gc2trig = FALSE
    }

    #
    # ... IV ####
    #

    if(debug) cat(" ---> IV\n")

    msgNA_iv = ""
    if(do_iv){
        # We create the IV LHS and the IV RHS

        iv_fml_parts = fml_split(fml_iv, split.lhs = TRUE, raw = TRUE)

        #
        # LHS
        #

        # Now LHS evaluated as a regular formula
        iv_endo_fml = .xpd(lhs = 1, rhs = iv_fml_parts[[1]])
        iv_lhs_mat = error_sender(fixest_model_matrix(iv_endo_fml, data, TRUE),
                                  "Evaluation of the left-hand-side of the IV part (equal to ", deparse_long(iv_fml_parts[[1]]), ") raises an error: \n")

        # Not great => to improve later
        iv_lhs = list()
        iv_lhs_names = colnames(iv_lhs_mat)
        for(i in 1:ncol(iv_lhs_mat)){
            iv_lhs[[iv_lhs_names[i]]] = iv_lhs_mat[, i]
        }

        # NAs
        n_NA_iv = c(0, 0)
        info = cpppar_which_na_inf_df(iv_lhs, nthreads)
        if(info$any_na_inf){

            if(info$any_na) ANY_NA = TRUE
            if(info$any_inf) ANY_INF = TRUE

            anyNA_iv = TRUE
            isNA_iv = info$is_na_inf
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_iv
            n_NA_iv[1] = sum(isNA_iv)
        }

        #
        # RHS
        #

        if(length(all.vars(iv_fml_parts[[2]])) == 0){
            stop("In the IV part, the RHS must contain at least one variable.")
        }

        iv.mat = error_sender(fixest_model_matrix(fml_iv, data, fake_intercept = TRUE),
                              "Problem in the IV part. Evaluation of the right-hand-side raises an error: ")

        inst_names = attr(terms(fml_iv), "term.labels")

        # NAs
        info = cpppar_which_na_inf_mat(iv.mat, nthreads)
        if(info$any_na_inf){

            if(info$any_na) ANY_NA = TRUE
            if(info$any_inf) ANY_INF = TRUE

            anyNA_iv = TRUE
            isNA_iv = info$is_na_inf
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_iv
            n_NA_iv[2] = sum(isNA_iv)
        }

        if(any(n_NA_iv > 0)){
            msgNA_iv = paste0("IV: ", numberFormatNormal(n_NA_iv[1]), "/", numberFormatNormal(n_NA_iv[2]))
        }

    }


    #
    # ... Cluster ####
    #

    if(debug) cat(" ---> cluster\n")

    msgNA_cluster = ""
    if(!missnull(cluster)){
        # cluster was checked already before subset => to get the right variables
        # Here cluster is a formula

        do_summary = TRUE

        cluster_terms = error_sender(fixef_terms(cluster), "Problem in the argument 'cluster':\n", clean = "_impossible_var_name_ => ^")
        if(any(cluster_terms$slope_flag != 0)){
            stop("You cannot use variables with varying slopes in the argument 'cluster'.")
        }

        cluster_terms = cluster_terms$fe_vars
        if(isFixef && all(cluster_terms %in% fixef_terms_full$fe_vars)){
                # we do nothing => the algo will use the FE ids
        } else {
            cluster = error_sender(prepare_df(cluster_terms, data, TRUE),
                                   "Problem evaluating the argument 'cluster':\n")

            # Type conversion
            for(i in seq_along(cluster)){
                if(!is.numeric(cluster[[i]]) && !is.character(cluster[[i]])){
                    cluster[[i]] = as.character(cluster[[i]])
                }
            }

            if(anyNA(cluster)){
                isNA_cluster = !complete.cases(cluster)

                ANY_NA = TRUE
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_cluster
                msgNA_cluster = paste0("cluster: ", numberFormatNormal(sum(isNA_cluster)))

                if(mem.clean){
                    rm(isNA_cluster)
                    gc2trig = TRUE
                }

            }
        }
    } else {
        # Cluster origin was first defined when we first check cluster, way up
        cluster_origin = cluster = NULL
    }

    if(!missnull(se)){
        do_summary = TRUE

        check_arg_plus(se, "match", .choices = c("standard", "white", "hetero", "cluster", "twoway", "threeway", "fourway", "1", "2", "3", "4"), .message = "Argument argument 'se' should be equal to one of 'standard', 'hetero', 'cluster', 'twoway', 'threeway' or 'fourway'.")

        # we check consistency
        if(isFixef){
            n_fe = length(fixef_terms_full$fe_vars)
        } else {
            n_fe = 0
        }

        n_clu = c(cluster = 1, twoway = 2, threeway = 3, fourway = 4)[se]
        if(missnull(cluster) && !is.na(n_clu) && n_clu > n_fe){
            stop("In argument 'se': ", n_letter(n_clu), "-way clustering cannot be done with the current estimation which has ", ifelse(n_fe == 0, "no fixed-effect.", paste0("only ", n_letter(n_fe), " fixed-effects.")), " Please provide the argument cluster.")
        }

    } else {
        se = NULL
    }

    if(!missnull(dof)){
        do_summary = TRUE
        if(!identical(class(dof), "dof.type")){
            stop("The argument 'dof' must be an object obtained from the function dof().")
        }
    } else {
        dof = NULL
    }

    check_arg(lean, "logical scalar")
    if(lean){
        do_summary = TRUE
    }


    #
    # Fixed-effects ####
    #

    if(debug) cat(" - Fixed-effects\n")

    # NOTA on stepwise FEs:
    # I wanted to evaluate all the FEs first, then send the evaluated stuff for later. Like in stepwise linear.
    # This is actually a bad idea because FEs are too complex to manipulate (damn SLOPES!!!).
    # This means that lags in the FEs + stepwise FEs will never be supported.

    isSlope = onlySlope = FALSE
    if(isFixef){
        # The main fixed-effects construction

        if(isFit){

            #
            # ... From fit ####
            #

            if(isVector(fixef_df)){
                fixef_df = data.frame(x = fixef_df, stringsAsFactors = FALSE)
                names(fixef_df) = deparse(mc_origin[["fixef_df"]])[1]

            } else if(is.list(fixef_df)){
                all_len = lengths(fixef_df)
                if(any(diff(all_len) != 0)){
                    stop("The lengths of the vectors in fixef_df differ (currently it is: ", enumerate_items(all_len), ").")
                }
                fixef_df = as.data.frame(fixef_df)
            }

            if(!is.matrix(fixef_df) && !"data.frame" %in% class(fixef_df)){
                stop("Argument fixef_df must be a vector, a matrix, a list or a data.frame (currently its class is ", enumerate_items(class(fixef_df)), ").")
            }

            if(is.matrix(fixef_df) && is.null(colnames(fixef_df))){
                colnames(fixef_df) = paste0("fixef_", 1:ncol(fixef_df))
                fixef_df = as.data.frame(fixef_df)
            }

            if(nrow(fixef_df) != nobs){
                stop("The number of observations of fixef_df (", nrow(fixef_df), ") must match the length of y (", nobs, ").")
            }

            fixef_vars = names(fixef_df)

            # The formula
            fml_fixef = .xpd(rhs = fixef_vars)

            if(delayed.subset){
                fixef_df = fixef_df[subset, , drop = FALSE]
            }

        } else {
            #
            # ... Regular ####
            #

            # Managing fast combine
            if(missnull(combine.quick)){
                if(NROW(data) > 5e4){
                    combine.quick = TRUE
                } else {
                    combine.quick = FALSE
                }
            }

            # fixef_terms_full computed in the formula section
            fixef_terms = fixef_terms_full$fml_terms

            # FEs
            fixef_df = error_sender(prepare_df(fixef_terms_full$fe_vars, data, combine.quick),
                                     "Problem evaluating the fixed-effects part of the formula:\n")

            fixef_vars = names(fixef_df)

            # Slopes
            isSlope = any(fixef_terms_full$slope_flag != 0)
            if(isSlope){

                slope_df = error_sender(prepare_df(fixef_terms_full$slope_vars, data),
                                         "Problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n")

                slope_flag = fixef_terms_full$slope_flag
                slope_vars = fixef_terms_full$slope_vars
                slope_vars_list = fixef_terms_full$slope_vars_list

                # Further controls
                not_numeric = !sapply(slope_df, is.numeric)
                if(any(not_numeric)){
                    stop("In the fixed-effects part of the formula (i.e. in ", as.character(fml_fixef[2]), "), variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope_df)[not_numeric], "s.is.quote"), " not.")
                }

                # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect
                onlySlope = all(slope_flag < 0)

            }

            # fml update
            fml_fixef = .xpd(rhs = fixef_terms)

        }

        #
        # ... NA handling ####
        #

        if(debug) cat(" ---> NA handling\n")

        # We change non-numeric to character (important for parallel qufing)
        is_not_num = sapply(fixef_df, function(x) !is.numeric(x))
        if(any(is_not_num)){
            for(i in which(is_not_num)){
                # we don't convert Dates to numeric because conversion can lead to NA in some cases!
                # e.g dates < 1970 with non-robust parsers
                if(!is.character(fixef_df[[i]])){
                    fixef_df[[i]] = as.character(fixef_df[[i]])
                }
            }
        }

        msgNA_fixef = ""
        if(anyNA(fixef_df)){
            isNA_fixef = !complete.cases(fixef_df)

            ANY_NA = TRUE
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_fixef
            msgNA_fixef = paste0("Fixed-effects: ", numberFormatNormal(sum(isNA_fixef)))

            if(mem.clean){
                rm(isNA_fixef)
                gc2trig = TRUE
            }

        }

        # NAs in slopes
        if(isSlope){

            if(mem.clean && gc2trig){
                gc()
                gc2trig = FALSE
            }

            # Convert to double => to remove in the future
            who_not_double = which(sapply(slope_df, is.integer))
            for(i in who_not_double){
                slope_df[[i]] = as.numeric(slope_df[[i]])
            }

            info = cpppar_which_na_inf_df(slope_df, nthreads)
            if(info$any_na_inf){

                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE

                isNA_slope = info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_slope
                msgNA_slope = paste0("Var. Slopes: ", numberFormatNormal(sum(isNA_slope)))

            } else {
                msgNA_slope = ""
            }

            if(mem.clean){
                rm(info)
                gc()
            }

        } else {
            msgNA_slope = ""
            fixef_terms = fixef_vars # both are identical if no slope
        }

        # Removing observations
        obs2remove_NA = c()
        if(anyNA_sample || any0W){
            # we remove all NAs obs + 0 weight obs
            gc2trig = TRUE

            details = c(msgNA_y, msgNA_L, msgNA_iv, msgNA_NL, msgNA_fixef, msgNA_slope,
                        msgNA_offset, msgNA_weight, msgNA_split, msgNA_cluster)
            msg_details = paste(details[nchar(details) > 0], collapse = ", ")

            nbNA = sum(isNA_sample)

            if(anyNA_sample){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                message_NA = paste0(numberFormatNormal(nbNA), " observation", plural(nbNA), " removed because of ", msg, " (", msg_details, ").")
                notes = c(notes, message_NA)
            }

            if(nbNA == nobs){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                stop("All observations contain ", msg, ". Estimation cannot be done. Breakup: ", msg_details, ".")
            }

            if(any0W){
                # 0 weight => like NAs
                isNA_sample = isNA_sample | is0W
                nbNA = sum(isNA_sample)

                if(nbNA == nobs){
                    if(anyNA_sample){
                        msg = msg_na_inf(ANY_NA, ANY_INF)
                        stop("All observations contain ", msg, " or are 0-weight. Estimation cannot be done. (0-weight: ", sum(is0W), ", breakup ", msg, ": ", msg_details, ")")
                    } else {
                        stop("All observations are 0-weight. Estimation cannot be done.")
                    }
                }
            }

            # we drop the NAs from the fixef matrix
            fixef_df = fixef_df[!isNA_sample, , drop = FALSE]
            obs2remove_NA = which(isNA_sample)

            if(isSlope){
                slope_df = slope_df[!isNA_sample, , drop = FALSE]
            }

            # we change the LHS variable
            if(is.list(lhs)){
                for(i in seq_along(lhs)){
                    lhs[[i]] = lhs[[i]][-obs2remove_NA]
                }
            } else {
                lhs = lhs[-obs2remove_NA]

                if(cpp_isConstant(lhs)){
                    # We absolutely need to control for that, otherwise, the code breaks later on

                    message(ifsingle(notes, "NOTE: ", "NOTES: "), paste(notes, collapse = "\n       "))
                    msg = "NAs"
                    if(any0W && !anyNA_sample) msg = "0-weight"
                    if(any0W && anyNA_sample) msg = "NAs and 0-weight"
                    stop("The dependent variable (after cleaning for ", msg, ") is a constant. Estimation cannot be done.")
                }
            }
        }

        #
        # ... QUF setup ####
        #

        if(debug) cat(" ---> QUF\n")

        info_fe = setup_fixef(fixef_df = fixef_df, lhs = lhs, fixef_vars = fixef_vars,
                              fixef.rm = fixef.rm, family = family, isSplit = isSplit,
                              split.full = split.full, origin_type = origin_type,
                              isSlope = isSlope, slope_flag = slope_flag, slope_df = slope_df,
                              slope_vars_list = slope_vars_list, nthreads = nthreads)

        Q               = info_fe$Q
        fixef_id        = info_fe$fixef_id
        fixef_names     = info_fe$fixef_names
        fixef_sizes     = info_fe$fixef_sizes
        fixef_table     = info_fe$fixef_table
        sum_y_all       = info_fe$sum_y_all
        lhs             = info_fe$lhs

        obs2remove      = info_fe$obs2remove
        fixef_removed   = info_fe$fixef_removed
        message_fixef   = info_fe$message_fixef

        slope_variables = info_fe$slope_variables
        slope_flag      = info_fe$slope_flag

        fixef_id_res    = info_fe$fixef_id_res
        fixef_sizes_res = info_fe$fixef_sizes_res
        new_order       = info_fe$new_order

        notes = c(notes, message_fixef)

        if(length(obs2remove_NA) > 0){
            # we update the value of obs2remove (will contain both NA and removed bc of outcomes)
            if(length(obs2remove) > 0){
                index_noNA = (1:(tail(obs2remove_NA, 1) + tail(obs2remove, 1)))[-obs2remove_NA]
                obs2remove_fixef = index_noNA[obs2remove]
            } else {
                obs2remove_fixef = c()
            }

            obs2remove = sort(c(obs2remove_NA, obs2remove_fixef))
        }

    } else {
        # There is no fixed-effect
        Q = 0

        # NA management is needed to create obs2remove
        if(anyNA_sample || any0W){

            nbNA = sum(isNA_sample)

            details = c(msgNA_y, msgNA_L, msgNA_iv, msgNA_NL, msgNA_offset,
                        msgNA_weight, msgNA_split, msgNA_cluster)
            msg_details = paste(details[nchar(details) > 0], collapse = ", ")

            if(anyNA_sample){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                message_NA = paste0(numberFormatNormal(nbNA), " observation", plural(nbNA), " removed because of ", msg, " (", msg_details, ").")
                notes = c(notes, message_NA)

                if(nbNA == nobs){
                    stop("All observations contain NAs. Estimation cannot be done. (Breakup: ", msg_details, ".)")
                }
            }

            if(any0W){
                # 0 weight => like NAs
                isNA_sample = isNA_sample | is0W
                nbNA = sum(isNA_sample)

                if(nbNA == nobs){
                    if(anyNA_sample){
                        stop("All observations are either NA or 0-weight. Estimation cannot be done. (0-weight: ", sum(is0W), ", breakup NA: ", msg_details, ".)")
                    } else {
                        stop("All observations are 0-weight. Estimation cannot be done.")
                    }
                }
            }

            # we drop the NAs from the fixef matrix
            obs2remove = which(isNA_sample)
        } else {
            obs2remove = c()
        }
    }

    # Messages
    if(show_notes && length(notes) > 0){
        message(ifsingle(notes, "NOTE: ", "NOTES: "), paste(notes, collapse = "\n       "))
    }

    #
    # NA removal ####
    #

    if(debug) cat(" - NA Removal\n")

    # NA & problem management
    if(length(obs2remove) > 0){
        # we kick out the problems (both NA related and fixef related)

        if(isLinear){
            # We drop only 0 variables (may happen for factors)
            linear.mat = select_obs(linear.mat, -obs2remove, nthreads)

            # useful for feNmlm
            linear.params = colnames(linear.mat)
            params = c(nonlinear.params, linear.params)
            lparams = length(params)
        }

        if(Q == 0){
            # if Q > 0: done already when managing the fixed-effects
            lhs = select_obs(lhs, -obs2remove)
        }

        if(do_iv){
            iv_lhs = select_obs(iv_lhs, -obs2remove)
            iv.mat = select_obs(iv.mat, -obs2remove, nthreads, "instrument")
        }

        if(isOffset){
            offset.value = offset.value[-obs2remove]
        }

        if(isWeight){
            weights.value = weights.value[-obs2remove]
        }

        if(isNonLinear){
            data_NL = data_NL[-obs2remove, , drop = FALSE]
        }

        if(isSplit){
            split = split[-obs2remove]
        }

        if(multi_rhs){
            linear_core = select_obs(linear_core, -obs2remove, nthreads)
            rhs_sw = select_obs(rhs_sw, -obs2remove, nthreads)
        }

        if(!is.null(cluster) && !"formula" %in% class(cluster)){
            cluster = select_obs(cluster, -obs2remove)
        }

    }

    # NEW NUMBER OF OBS AFTER NA/FE REMOVAL
    if(is.list(lhs)){
        nobs = length(lhs[[1]])
    } else {
        nobs = length(lhs)
    }

    # If presence of FEs => we exclude the intercept only if not all slopes
    if(Q > 0 && onlySlope == FALSE){
        # If there is a linear intercept, we withdraw it
        # We drop the intercept:

        if(isLinear && "(Intercept)" %in% colnames(linear.mat)){
            var2remove = which(colnames(linear.mat) == "(Intercept)")
            if(ncol(linear.mat) == length(var2remove)){
                isLinear = FALSE
                linear.params = NULL
                params = nonlinear.params
                lparams = length(params)
                varnames = nonlinear.varnames
            } else{
                linear.mat = linear.mat[, -var2remove, drop = FALSE]
                linear.params = colnames(linear.mat)
                params = c(nonlinear.params, linear.params)
                lparams = length(params)
                varnames = c(nonlinear.varnames, linear.varnames)
            }
        }

    }

    if(do_iv && "(Intercept)" %in% colnames(iv.mat)){
        var2remove = which(colnames(iv.mat) == "(Intercept)")
        iv.mat = iv.mat[, -var2remove, drop = FALSE]
    }


    #
    # Other setups ####
    #

    if(debug) cat(" - Other setups\n")

    # Starting values:
    # NOTA: within this function it is called linear.start while in the client functions, it is called start
    if(isLinear && !missing(linear.start)){
        # we format linear.start => end product is a either a 1-length vector or a named vector

        n_vars = ncol(linear.mat)
        n_linstart = length(linear.start)

        if(n_linstart == 0){
            stop("If 'start' is provided, its length should be equal to the number of variables (", n_vars, ") or equal to 1 (or a named-vector whose names are the variable names). Currently it is zero-length.")

        } else if(!is.numeric(linear.start)){
            stop("Argument 'start' must be a numeric vector (currently its class is ", class(linear.start)[1], ").")

        } else if(!isVector(linear.start)){
            stop("Argument 'start' must be a numeric vector (currently it is not a vector).")

        } else if(is.null(names(linear.start))){
            if(n_linstart == 1){
                # Nothing to do
            } else if(n_linstart == n_vars){
                # we name it
                names(linear.start) = linear.params
            } else {
                nota = ""
                if(isFixef && n_linstart == n_vars + 1) nota = ", note that there should be no intercept when fixed-effects are present"
                if(!isFixef && n_linstart == n_vars - 1) nota = ", maybe the intercept is missing?"

                stop("The length of 'start' should be ", n_vars," (it is currently equal to ", length(linear.start), nota, ").")
            }
        } else {
            # this is a named vector => if no variable is found: warning
            if(!any(names(linear.start) %in% linear.params)){
                warning("In argument 'start': no name found to match the variables.")
            }
        }
    }

    # Starting values feglm
    if(origin_type == "feglm"){

        assign("starting_values", NULL, env)

        if(sum(start_provided) >= 2){
            stop("You cannot provide more than one initialization method (currently there are ", enumerate_items(names(start_provided)[start_provided == 1]), ").")
        }

        # controls for start performed already
        if(start_provided["etastart"] || start_provided["mustart"]){

            # controls are identical for the two
            if(start_provided["mustart"]){
                starting_values = mustart
                arg_name = "mustart"
            } else {
                starting_values = etastart
                arg_name = "etastart"
            }


            # controls + avoid integers
            check_value_plus(starting_values, "numeric vector conv", .arg_name = arg_name)

            if(length(starting_values) == 1){
                starting_values = rep(starting_values, nobs)

            } else if(length(starting_values) == nobs){
                # Fine (although that's weird the user knows the final length ex ante)

            } else if(length(starting_values) == nobs){ # nobs => original nb of obs
                if(length(obs2remove) > 0){
                    starting_values = starting_values[-obs2remove]
                }

            } else {
                stop("The length '", arg_name, "' (", numberFormatNormal(length(starting_values)), ") differs from the length of the data (", numberFormatNormal(nobs), ").")
            }

            info = cpppar_which_na_inf_vec(starting_values, nthreads)
            if(info$any_na_inf){
                msg = msg_na_inf(info$any_na, info$any_inf)

                obs_origin = (1:nobs)
                if(length(obs2remove) > 0){
                    obs_origin = obs_origin[-obs2remove]
                }

                obs = head(obs_origin[which(info$is_na_inf)], 3)

                stop("The argument '", arg_name, "' contains ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide starting values without ", msg, ".")
                # This is not exactly true, because if the NAs are the observations removed, it works
                # but it's hard to clearly explain it concisely
            }

            assign("starting_values", starting_values, env)
        }

        assign("init.type", "default", env)
        if(any(start_provided)){
            assign("init.type", c("coef", "eta", "mu")[start_provided], env)
        }

        if(!isLinear && get("init.type", env) == "coef"){
            stop("You cannot initialize feglm with coefficients when there is no variable to be estimated (only the fixed-effects).")
        }

    }

    #
    # families of feglm
    #

    if(origin_type == "feglm"){
        # We optimize Poisson and Logit

        family_equiv = family_funs$family_equiv

        if(family_equiv == "poisson"){
            family_funs$linkfun = function(mu) cpppar_log(mu, nthreads)
            family_funs$linkinv = function(eta) cpppar_poisson_linkinv(eta, nthreads)

            if(multi_lhs == FALSE && isSplit == FALSE){
                # => we simply delay the assignment for multi_lhs == TRUE or split
                # It will be done in reshape_env
                y_pos = lhs[lhs > 0]
                qui_pos = lhs > 0
                if(isWeight){
                    constant = sum(weights.value[qui_pos] * y_pos * cpppar_log(y_pos, nthreads) - weights.value[qui_pos] * y_pos)
                    sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(wt[qui_pos] * y_pos * eta[qui_pos]) + sum(wt * mu))
                } else {
                    constant = sum(y_pos * cpppar_log(y_pos, nthreads) - y_pos)
                    sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(y_pos * eta[qui_pos]) + sum(mu))
                }

                family_funs$sum_dev.resids = sum_dev.resids
            }


            family_funs$mu.eta = function(mu, eta) mu
            family_funs$validmu = function(mu) cpppar_poisson_validmu(mu, nthreads)

        } else if(family_equiv == "logit"){
            family_funs$linkfun = function(mu) cpppar_logit_linkfun(mu, nthreads)
            family_funs$linkinv = function(eta) cpppar_logit_linkinv(eta, nthreads)
            if(isWeight){
                sum_dev.resids = function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, wt, nthreads))
            } else {
                sum_dev.resids = function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, 1, nthreads))
            }
            family_funs$sum_dev.resids = sum_dev.resids

            family_funs$mu.eta = function(mu, eta) cpppar_logit_mueta(eta, nthreads)
        }

        assign("family_funs", family_funs, env)

    }

    if(lparams == 0 && Q == 0 && !multi_fixef && !multi_rhs) stop("No parameter to be estimated.")

    check_arg(useHessian, "logical scalar")
    assign("hessian.args", hessian.args, env)

    if(origin == "feNmlm"){
        jacobian.method = check_arg_plus(jacobian.method, "match(simple, Richardson)")
    }

    #
    # Controls: The non linear part
    #

    if(isNL){
        if(missing(NL.start.init)){
            if(missing(NL.start)) stop("There must be starting values for NL parameters. Please use argument NL.start (or NL.start.init).")
            if(length(NL.start) == 1 && is.numeric(NL.start)){
                # NL.start is used for NL.start.init
                if(is.na(NL.start)) stop("Argument 'NL.start' is currently equal to NA. Please provide a valid starting value.")
                NL.start.init = NL.start
                NL.start <- list()
                NL.start[nonlinear.params] <- NL.start.init

            } else if(typeof(NL.start) != "list"){
                stop("NL.start must be a list.")
            } else if(any(!nonlinear.params %in% names(NL.start))){
                stop(paste("Some NL parameters have no starting values:\n", paste(nonlinear.params[!nonlinear.params %in% names(NL.start)], collapse=", "), ".", sep=""))
            } else {
                # we restrict NL.start to the nonlinear.params
                NL.start = NL.start[nonlinear.params]
            }

        } else {
            if(length(NL.start.init)>1) stop("NL.start.init must be a scalar.")
            if(!is.numeric(NL.start.init)) stop("NL.start.init must be numeric!")
            if(!is.finite(NL.start.init)) stop("Non-finite values as starting values, you must be kidding me...")

            if(missing(NL.start)){
                NL.start <- list()
                NL.start[nonlinear.params] <- NL.start.init
            } else {
                if(typeof(NL.start)!="list") stop("NL.start must be a list.")
                if(any(!names(NL.start) %in% params)) stop(paste("Some parameters in 'NL.start' are not in the formula:\n", paste(names(NL.start)[!names(NL.start) %in% params], collapse=", "), ".", sep=""))

                missing.params <- nonlinear.params[!nonlinear.params%in%names(NL.start)]
                NL.start[missing.params] <- NL.start.init
            }
        }
    } else {
        NL.start <- list()
    }

    #
    # The upper and lower limits (only NL + negbin)
    #

    if(!missnull(lower)){

        if(typeof(lower )!= "list"){
            stop("'lower' MUST be a list.")
        }

        lower[params[!params %in% names(lower)]] = -Inf
        lower = unlist(lower[params])

    }	else {
        lower = rep(-Inf, lparams)
        names(lower) = params
    }

    if(!missnull(upper)){

        if(typeof(upper) != "list"){
            stop("'upper' MUST be a list.")
        }

        upper[params[!params %in% names(upper)]] = Inf
        upper = unlist(upper[params])

    }	else {
        upper = rep(Inf, lparams)
        names(upper) = params
    }

    lower = c(lower)
    upper = c(upper)

    #
    # Controls: user defined gradient
    #

    #
    # PRECISION + controls ####
    #

    if(debug) cat(" - Precision setting + controls\n")

    # The main precision

    check_value(fixef.tol, "numeric scalar GT{(10000*.Machine$double.eps)}")
    check_arg(fixef.iter, "integer scalar GE{1}")

    if(origin_type == "feNmlm"){
        check_arg(deriv.tol, "numeric scalar GT{(10000*.Machine$double.eps)}")
        check_arg(deriv.iter, "integer scalar GE{1}")

        # Other: opt.control
        check_arg(opt.control, "list l0", .message = "Argument opt.control must be a list of controls for the function nlminb (see help for details).")

    }

    if(origin_type == "feglm"){
        check_arg(glm.tol, "numeric scalar gt{(1000*.Machine$double.eps)}")
        check_arg(glm.iter, "integer scalar GE{1}")

        assign("glm.iter", glm.iter, env)
        assign("glm.tol", glm.tol, env)
    }

    # other precisions
    NR.tol = fixef.tol/100

    # Initial checks are done
    nonlinear.params <- names(NL.start) #=> in the order the user wants

    # starting values of all parameters (initialized with the NL ones):
    start = NL.start

    # control for the linear start => we can provide coefficients
    # from past estimations. Coefficients that are not provided are set
    # to 0
    if(length(linear.start) > 1){
        what = linear.start[linear.params]
        what[is.na(what)] = 0
        linear.start = what
    }
    start[linear.params] = linear.start
    params = names(start)
    start = unlist(start)
    start = c(start)
    lparams = length(params)
    names(start) = params

    # The right order of upper and lower
    upper = upper[params]
    lower = lower[params]

    ####
    #### Sending to the env ####
    ####

    if(debug) cat(" - sending to the env\n")

    #
    # General elements
    #

    # Control
    assign("fixef.iter", fixef.iter, env)
    assign("fixef.iter.limit_reached", 0, env) # for warnings if max iter is reached
    assign("fixef.tol", fixef.tol, env)
    assign("collin.tol", collin.tol, env)
    assign("fixef.rm", fixef.rm, env)

    # Data
    if(isLinear){
        assign("linear.mat", linear.mat, env)
    } else {
        assign("linear.mat", 1, env)
    }
    assign("offset.value", offset.value, env)
    assign("weights.value", weights.value, env)
    assign("lhs", lhs, env)
    assign("lhs_names", lhs_names, env)
    assign("nobs", nobs, env)
    assign("nobs_origin", nobs_origin, env)

    # Other
    assign("family", family, env)
    assign("famFuns", famFuns, env)
    assign("fml", fml_linear, env)
    assign("origin", origin, env)
    assign("origin_type", origin_type, env)
    assign("warn", warn, env)
    assign("mem.clean", mem.clean, env)
    assign("nthreads", nthreads, env)
    assign("demeaned", demeaned, env)

    # Summary
    assign("do_summary", do_summary, env)
    assign("lean", lean, env)
    if(do_summary){
        assign("cluster", cluster, env)
        assign("se", se, env)
        assign("dof", dof, env)
        assign("agg", agg, env)

        assign("summary_flags", build_flags(mc_origin, se = se, cluster = cluster_origin, dof = dof, agg = agg), env)
    }

    # Multi

    assign("do_multi_lhs", multi_lhs, env)

    assign("do_multi_rhs", multi_rhs, env)
    if(multi_rhs){
        assign("rhs_info_stepwise", rhs_info_stepwise, env)
        assign("linear_core", linear_core, env)
        assign("rhs_sw", rhs_sw, env)

    }

    assign("do_multi_fixef", multi_fixef, env)
    if(multi_fixef){
        assign("multi_fixef_fml_full", fixef_info_stepwise$fml_all_full, env)

        if(missnull(combine.quick)){
            combine.quick = TRUE
        }
        assign("combine.quick", combine.quick, env)

        var_sw = unique(fixef_info_stepwise$sw_all_vars)

        if(length(var_sw) > 0){
            if(length(obs2remove) > 0){
                assign("data", data[-obs2remove, var_sw, drop = FALSE], env)
            } else {
                assign("data", data[, var_sw, drop = FALSE], env)
            }
        }
    }

    assign("IN_MULTI", multi_lhs || multi_rhs || multi_fixef || isSplit, env)


    # IV
    assign("do_iv", do_iv, env)
    if(do_iv){
        assign("iv_lhs", iv_lhs, env)
        assign("iv_lhs_names", iv_lhs_names, env)
        assign("iv.mat", iv.mat, env)
    }

    #
    # The MODEL0 => to get the init of the theta for the negbin
    #

    check_arg(theta.init, "numeric scalar gt{0}")
    if(missing(theta.init)){
        theta.init = NULL
    }
    assign("theta.init", theta.init, env)

    assign("start", start, env)
    assign("isNL", isNL, env)
    assign("linear.params", linear.params, env)
    assign("isLinear", isLinear, env)
    assign("nonlinear.params", nonlinear.params, env)
    assign("nonlinear.varnames", nonlinear.varnames, env)
    assign("params", params, env)

    #
    # The Fixed-effects
    #


    assign("isFixef", isFixef, env)
    if(isFixef){

        if(!isSlope){
            slope_vars_list = list(0)
        }

        assign("new_order_original", new_order, env)
        assign("fixef_names", fixef_names, env)
        assign("fixef_vars", fixef_vars, env)

        assign_fixef_env(env, family, origin_type, fixef_id, fixef_sizes, fixef_table,
                         sum_y_all, slope_flag, slope_variables, slope_vars_list)
    }


    #
    # Specific to femlm/feNmlm
    #

    if(origin_type == "feNmlm"){

        # basic NL
        envNL = new.env()
        if(isNL){
            data_NL = as.data.frame(data_NL)
            for(var in nonlinear.varnames) assign(var, data_NL[[var]], envNL)
            for(var in nonlinear.params) assign(var, start[var], envNL)
        }
        assign("envNL", envNL, env)
        assign("nl.call", nl.call, env)
        # NO user defined gradient: too complicated, not really efficient
        assign("isGradient", FALSE, env)
        assign("lower", lower, env)
        assign("upper", upper, env)

        # other
        assign("iter", 0, env)
        assign("jacobian.method", jacobian.method, env)
        # Pour gerer les valeurs de mu:
        assign("coefMu", list(), env)
        assign("valueMu", list(), env)
        assign("valueExpMu", list(), env)
        assign("wasUsed", TRUE, env)
        # Pour les valeurs de la Jacobienne non lineaire
        assign("JC_nbSave", 0, env)
        assign("JC_nbMaxSave", 1, env)
        assign("JC_savedCoef", list(), env)
        assign("JC_savedValue", list(), env)
        # PRECISION
        assign("NR.tol", NR.tol, env)
        assign("deriv.tol", deriv.tol, env)
        # ITERATIONS
        assign("deriv.iter", deriv.iter, env)
        assign("deriv.iter.limit_reached", 0, env) # for warnings if max iter is reached
        # OTHER
        assign("useAcc", TRUE, env)
        assign("warn_0_Hessian", FALSE, env)
        assign("warn_overfit_logit", FALSE, env)
        # Misc
        assign("opt.control", opt.control, env)

        # To monitor how the FEs are computed (if the problem is difficult or not)
        assign("firstIterCluster", 1e10, env) # the number of iterations in the first run
        assign("firstRunCluster", TRUE, env) # flag for first entry in get_dummies
        assign("iterCluster", 1e10, env) # the previous number of fixef iterations
        assign("evolutionLL", Inf, env) # diff b/w two successive LL
        assign("pastLL", 0, env)
        assign("iterLastPrecisionIncrease", 0, env) # the last iteration when precision was increased
        assign("nbLowIncrease", 0, env) # number of successive evaluations with very low LL increments
        assign("nbIterOne", 0, env) # nber of successive evaluations with only 1 iter to get the clusters
        assign("difficultConvergence", FALSE, env)

        # Same for derivatives
        assign("derivDifficultConvergence", FALSE, env)
        assign("firstRunDeriv", TRUE, env) # flag for first entry in derivative
        assign("accDeriv", TRUE, env) # Logical: flag for accelerating deriv
        assign("iterDeriv", 1e10, env) # number of iterations in the derivatives step

        # NL: On teste les valeurs initiales pour informer l'utilisateur

        if(isNL){

            mu = error_sender(eval(nl.call, envir = envNL), "The non-linear part (NL.fml) could not be evaluated.:\n")

            # No numeric vectors
            check_value(mu, "numeric vector", .message = "Evaluation of NL.fml should return a numeric vector.")

            # Handling NL.fml errors
            if(length(mu) != nrow(data_NL)){
                if(length(mu) == 1 && length(nonlinear.varnames) == 0){
                    stop("No variable from the data is in the 'NL.fml', there must be a problem.")
                } else {
                    stop("Evaluation of NL.fml leads to ", length(mu), " observations while there are ", nrow(data_NL), " observations in the data base. They should be of the same lenght.")
                }

            }

            if(anyNA(mu)){
                stop("Evaluating NL.fml leads to NA values: this is not supported. Maybe it's a problem with the starting values, maybe it's another problem.")
            }

        } else {
            mu = eval(nl.call, envir = envNL)
        }

        # On sauvegarde les valeurs de la partie non lineaire
        assign("nbMaxSave", 2, env) # nombre maximal de valeurs a sauvegarder
        assign("nbSave", 1, env)  # nombre de valeurs totales sauvegardees
        assign("savedCoef", list(start[nonlinear.params]), env)
        assign("savedValue", list(mu), env)

        # Mise en place du calcul du gradient
        gradient = femlm_gradient
        hessian <- NULL
        if(useHessian) hessian <- femlm_hessian
        assign("gradient", gradient, env)
        assign("hessian", hessian, env)
    }

    # split
    assign("do_split", isSplit, env)
    if(isSplit){
        split = to_integer(split, add_items = TRUE, sorted = TRUE, items.list = TRUE)
        split.items = as.character(split$items)
        split = split$x

        assign("split", split, env)
        assign("split.items", split.items, env)
        assign("split.name", split.name, env)
        assign("split.full", split.full, env)
    }

    # fixest tag
    assign("fixest_env", TRUE, env)

    #
    # Res ####
    #

    #
    # Preparation of the results list (avoids code repetition in estimation funs)
    #

    res = list(nobs=nobs, nobs_origin=nobs_origin, fml=fml_linear, call = mc_origin,
               call_env = call_env, method = origin, method_type = origin_type)

    fml_all = list()
    fml_all$linear = fml_linear
    fml_all$fixef = fml_fixef
    fml_all$NL = NL.fml
    fml_all$iv = fml_iv
    res$fml_all = fml_all

    if(!is.null(fml_no_xpd)) res$fml_no_xpd = fml_no_xpd

    if(isFit) res$fromFit = TRUE

    if(do_iv){
        res$iv = TRUE
        res$iv_inst_names = inst_names
        res$iv_n_inst = ncol(iv.mat)
        res$iv_endo_names = iv_lhs_names
        res$iv_endo_fml = iv_endo_fml
    }

    # nber of params
    K = length(params)
    if(isFixef){
        if(isSlope){
            K = K + sum(fixef_sizes * (1 + abs(slope_flag) - (slope_flag < 0)))
            # now the references
            if(any(slope_flag >= 0)){
                K = K + 1 - sum(slope_flag >= 0)
            }
        } else {
            K = K + sum(fixef_sizes - 1) + 1
        }
    }
    res$nparams = K

    # NL information
    if(isNL){
        res$NL.fml = NL.fml
    }

    # Fixed-effects information
    if(isFixef){

        res$fixef_vars = fixef_vars

        if(isSlope){
            res$fixef_terms = fixef_terms
            res$slope_flag = slope_flag[order(new_order)]
            res$slope_flag_reordered = slope_flag
            res$slope_variables_reordered = slope_variables
            res$fe.reorder = new_order
        }

        res$fixef_id = fixef_id_res
        res$fixef_sizes = fixef_sizes_res

    }

    # Observations removed (either NA or fixed-effects)

    if(length(obs2remove) > 0){
        obs_selection$obsRemoved = -obs2remove
        if(isFixef && any(lengths(fixef_removed) > 0)){
            res$fixef_removed = fixef_removed
        }
    }

    res$obs_selection = obs_selection

    # offset and weight
    if(isOffset){
        res$offset = offset.value
    }
    if(isWeight){
        res$weights = weights.value
    }

    # We save lhs in case of feglm, for later within-r2 evaluation
    if(origin_type == "feglm" && isFixef && !multi_lhs){
        # in multi_lhd => will be done in reshape_env
        res$y = lhs
    }

    # Panel information
    if(isPanel){
        res$panel.id = panel.id
        if(!is.null(panel.info)){
            res$panel.info = panel.info
        }
    }

    # Interaction information
    if(length(GLOBAL_fixest_mm_info) > 0){
        res$model_matrix_info = GLOBAL_fixest_mm_info
        if("sunab" %in% names(GLOBAL_fixest_mm_info)){
            res$is_sunab = TRUE
        }
    }

    assign("res", res, env)

    env
}


setup_fixef = function(fixef_df, lhs, fixef_vars, fixef.rm, family, isSplit, split.full = FALSE,
                       origin_type, isSlope, slope_flag, slope_df, slope_vars_list,
                       fixef_names_old = NULL, fixef_sizes = NULL, obs2keep = NULL, nthreads){

    isSplitNoFull = identical(fixef_names_old, "SPLIT_NO_FULL")
    isRefactor = !isSplitNoFull && !is.null(fixef_names_old)

    Q = length(fixef_vars) # terms: contains FEs + slopes

    rm_0 = !family == "gaussian"
    rm_1 = family == "logit"
    rm_single = fixef.rm %in% c("singleton", "both")
    do_sum_y = !origin_type %in% c("feols", "feglm")

    multi_lhs = is.list(lhs)

    if(multi_lhs){
        # we delay the removal of FEs to later, when called in reshape env
        rm_0 = rm_1 = do_sum_y = FALSE
    }

    only_slope = FALSE
    if(isSlope){
        # only slope: not any non-slope
        only_slope = slope_flag < 0

        # To prevent the case where one variable can have != varying slopes,
        # we need this condition
        # (super corner case....)
        my_vs_vars = unlist(slope_vars_list, use.names = FALSE)
        if(length(slope_df) == length(my_vs_vars)){
            # simple case
            slope_variables = as.list(slope_df)

        } else {
            slope_variables = list()
            for(i in seq_along(my_vs_vars)){
                slope_variables[[i]] = slope_df[[my_vs_vars[i]]]
            }
            # unfortunately, we need the names...
            names(slope_variables) = my_vs_vars
        }

        if(!is.null(obs2keep)){
            for(i in seq_along(slope_variables)){
                slope_variables[[i]] = slope_variables[[i]][obs2keep]
            }
        }
    } else {
        slope_flag = rep(0L, length(fixef_df))
        slope_variables = list(0)
    }

    if(isSplit && split.full == FALSE){
        # We don't do anything => it will be taken care of in the splits
        res = list(Q = Q, fixef_id = fixef_df, fixef_names = "SPLIT_NO_FULL",
                   sum_y_all = 0, fixef_sizes = 0, fixef_table = 0, obs2remove = NULL,
                   fixef_removed = NULL, message_fixef = NULL, lhs = lhs,
                   slope_variables = slope_variables, slope_flag = slope_flag, new_order = 1:Q)

        return(res)
    }

    do_keep = !is.null(obs2keep)

    if(isSplitNoFull && do_keep){
        # Here fixef_df is a DF or a list (outcome of previous calls) => that's why I need to use select_obs
        fixef_df = select_obs(fixef_df, obs2keep)
    }

    if(is.null(obs2keep)){
        obs2keep = 0
    }

    if(is.null(fixef_sizes)){
        fixef_sizes = 0
    }

#     quoi = list(fixef_df=fixef_df, lhs=lhs, do_sum_y=do_sum_y, rm_0=rm_0, rm_1=rm_1, rm_single=rm_single, only_slope=only_slope, nthreads=nthreads, isRefactor=isRefactor, fixef_sizes=fixef_sizes, obs2keep=obs2keep)
#     save(quoi, file = "../_PROBLEM/fepois crashes/problem_quf.Rdata")
# stop()

    quf_info_all = cpppar_quf_table_sum(x = fixef_df, y = lhs, do_sum_y = do_sum_y,
                                        rm_0 = rm_0, rm_1 = rm_1, rm_single = rm_single,
                                        only_slope = only_slope, nthreads = nthreads,
                                        do_refactor = isRefactor, r_x_sizes = fixef_sizes, obs2keep = obs2keep)

    fixef_id = quf_info_all$quf
    # names

    fixef_names = list()
    if(isRefactor == FALSE){
        is_string = sapply(fixef_df, is.character)
        for(i in 1:length(fixef_id)){
            if(is_string[i]){
                fixef_names[[i]] = fixef_df[[i]][quf_info_all$items[[i]]]
            } else {
                fixef_names[[i]] = quf_info_all$items[[i]]
            }
        }

    } else {
        # If here, we're in refactoring. The new fixef are integers
        for(i in 1:length(fixef_id)){
            fixef_names[[i]] = fixef_names_old[[i]][quf_info_all$items[[i]]]
        }
    }


    # table/sum_y/sizes
    fixef_table = quf_info_all$table
    sum_y_all = quf_info_all$sum_y
    fixef_sizes = lengths(fixef_table)

    # If observations have been removed:
    fixef_removed = list()
    obs2remove = message_fixef = c()

    if(!is.null(quf_info_all$obs_removed)){

        # which obs are removed
        obs2remove = which(quf_info_all$obs_removed)

        # update of the lhs
        # if multi_lhs, only reason we're here is because of rm_single, which is performed
        # only in main fixest_env
        if(!identical(lhs, list(0))){
            # list(0): should be unnecessary, when assign_lhs = FALSE
            # (singletons shouldn't be there any more, so no reason to be here with a list)
            lhs = select_obs(lhs, -obs2remove)
        }

        # update of the slope variables
        if(isSlope){
            for(i in seq_along(slope_variables)) slope_variables[[i]] = slope_variables[[i]][-obs2remove]
        }

        # Names of the FE removed
        for(i in 1:length(fixef_id)){
            if(is.character(fixef_df[[i]])){
                fixef_removed[[i]] = fixef_df[[i]][quf_info_all$fe_removed[[i]]]
            } else {
                fixef_removed[[i]] = quf_info_all$fe_removed[[i]]
            }
        }

        names(fixef_removed) = fixef_vars

        # Then the "Notes"
        nb_missing = lengths(fixef_removed)
        if(rm_0 == FALSE){
            n_single = sum(nb_missing)
            message_fixef = paste0(fsignif(n_single), " fixed-effect singleton", plural(n_single, "s.was"), " removed (", numberFormatNormal(length(obs2remove)), " observation", plural_len(obs2remove), ifelse(Q == 1, "", paste0(", breakup: ", paste0(nb_missing, collapse = "/"))), ").")
        } else {
            message_fixef = paste0(paste0(fsignif(nb_missing), collapse = "/"), " fixed-effect", plural(sum(nb_missing)), " (", numberFormatNormal(length(obs2remove)), " observation", plural_len(obs2remove), ") removed because of only ", ifelse(rm_1, "0 (or only 1)", "0"), " outcomes", ifelse(rm_single && !rm_1, " or singletons", ""), ".")
        }

    }

    if(multi_lhs == FALSE || family == "gaussian"){
        #
        # we save the fixed-effects IDs + sizes (to be returned in "res") [original order!]
        #

        fixef_id_res = list()
        for(i in 1:Q){
            dum = fixef_id[[i]]
            attr(dum, "fixef_names") = as.character(fixef_names[[i]])
            fixef_id_res[[fixef_vars[i]]] = dum
        }

        # The real size is equal to nb_coef * nb_slopes
        if(isSlope){
            fixef_sizes_real = fixef_sizes * (1 + abs(slope_flag) - (slope_flag < 0))
        } else {
            fixef_sizes_real = fixef_sizes
        }

        fixef_sizes_res = fixef_sizes
        names(fixef_sizes_res) = fixef_vars

        #
        # We re-order the fixed-effects
        #

        if(any(fixef_sizes_real != sort(fixef_sizes_real, decreasing = TRUE))){
            # FE with the most cases first (limits precision problems)

            new_order = order(fixef_sizes_real, decreasing = TRUE)

            fixef_sizes = fixef_sizes[new_order]
            fixef_id = fixef_id[new_order]
            sum_y_all = sum_y_all[new_order]
            fixef_table = fixef_table[new_order]

            if(isSlope){
                slope_variables = slope_variables[unlist(slope_vars_list[new_order], use.names = FALSE)]
                slope_flag = slope_flag[new_order]
            }

        } else {
            new_order = 1:Q
        }
    }

    # saving
    res = list(Q = Q, fixef_id = fixef_id, fixef_names = fixef_names, sum_y_all = sum_y_all,
               fixef_sizes = fixef_sizes, fixef_table = fixef_table, obs2remove = obs2remove,
               fixef_removed = fixef_removed, message_fixef = message_fixef, lhs = lhs)

    res$slope_variables = slope_variables
    res$slope_flag = slope_flag

    if(multi_lhs == FALSE || family == "gaussian"){
        res$fixef_id_res = fixef_id_res
        res$fixef_sizes_res = fixef_sizes_res
        res$new_order = new_order
    }

    return(res)
}





assign_fixef_env = function(env, family, origin_type, fixef_id, fixef_sizes, fixef_table,
                            sum_y_all, slope_flag, slope_variables, slope_vars_list){

    Q = length(fixef_id)

    assign("fixef_id_list", fixef_id, env)
    assign("fixef_sizes", fixef_sizes, env)
    assign("sum_y", sum_y_all, env)
    assign("fixef_table", fixef_table, env)
    assign("familyConv", family, env) # new family -- used in convergence, can be modified

    if(family %in% c("negbin")){
        # for the deriv_xi_other
        fixef_id.matrix_cpp = NULL
        if(Q > 1){
            fixef_id.matrix_cpp = matrix(unlist(fixef_id), ncol = Q) - 1
        }

        assign("fixef_id.matrix_cpp", fixef_id.matrix_cpp, env)
    }

    # New cpp functions
    # This is where we send the elements needed for convergence in cpp
    if(origin_type == "feNmlm"){
        assign("fixef_id_vector", as.integer(unlist(fixef_id) - 1), env)
    }

    assign("fixef_table_vector", as.integer(unlist(fixef_table)), env)
    assign("sum_y_vector", unlist(sum_y_all), env)

    if(origin_type == "feNmlm"){

        useExp_fixefCoef = family %in% c("poisson")

        nobs = length(fixef_id[[1]])

        if(useExp_fixefCoef){
            assign("saved_sumFE", rep(1, nobs), env)
        } else {
            assign("saved_sumFE", rep(0, nobs), env)
        }

        if(family == "gaussian" && Q >= 2){
            assign("fixef_invTable", 1/unlist(fixef_table), env)
        } else {
            assign("fixef_invTable", 1, env)
        }

        if(family %in% c("negbin", "logit")){
            fixef_order = list()
            for(i in 1:Q){
                fixef_order[[i]] = order(fixef_id[[i]]) - 1
            }
            assign("fixef_cumtable_vector", as.integer(unlist(lapply(fixef_table, cumsum))), env)
            assign("fixef_order_vector", as.integer(unlist(fixef_order)), env)
        } else {
            # we need to assign it anyway
            assign("fixef_cumtable_vector", 1L, env)
            assign("fixef_order_vector", 1L, env)
        }
    }

    # Slopes
    assign("slope_flag", slope_flag, env)
    assign("slope_variables", slope_variables, env)
    assign("slope_vars_list", slope_vars_list, env)
}



reshape_env = function(env, obs2keep = NULL, lhs = NULL, rhs = NULL, assign_lhs = TRUE,
                       assign_rhs = TRUE, fml_linear = NULL, fml_fixef = NULL, fml_iv_endo = NULL,
                       check_lhs = FALSE, assign_fixef = FALSE){
    # env: environment from an estimation
    # This functions reshapes the environment to perform the new estimation
    # either by selecting some observation (in split)
    # either by changing the depvar (leading to new NAs, etc)
    # obs2keep => must be a which()

    # assing_lhs/assign_rhs:
    # => avoids assigning inappropriate values in res
    # => present only in OLS!

    # Mainly:
    # - dropping the NAs from the new y
    # - reformatting of res
    # - refactoring the fixef and VSs
    # - dropping the factor values from i() that are only 0

    # gt = function(x) cat(sfill(x, 20), ": ", -(t0 - (t0<<-proc.time()))[3], "s\n", sep = "")
    # t0 = proc.time()

    new_env = new.env(parent = env)

    isFixef     = get("isFixef", env)
    origin_type = get("origin_type", env)
    family      = get("family", env)
    nthreads    = get("nthreads", env)

    save_lhs = !missnull(lhs)
    save_rhs = !missnull(rhs)

    if(assign_lhs == FALSE){
        lhs = list(0)
    } else if(is.null(lhs)){
        lhs = get("lhs", env)
    }

    lhs_done = FALSE

    res = get("res", new_env)

    #
    # fixef ####
    #

    if(isFixef && (length(obs2keep) > 0 || check_lhs || assign_fixef)){

        # gt("nothing")

        if(assign_lhs && length(obs2keep) > 0){
            lhs = select_obs(lhs, obs2keep)
        }

        # gt("fixef, dropping lhs")

        fixef_df       = get("fixef_id_list", env)
        fixef_vars      = get("fixef_vars", env)
        fixef.rm        = get("fixef.rm", env)
        fixef_names_old = get("fixef_names", env)
        fixef_sizes     = get("fixef_sizes", env)

        slope_flag      = get("slope_flag", env)
        slope_df        = get("slope_variables", env)
        slope_vars_list = get("slope_vars_list", env)
        isSlope         = any(slope_flag != 0)

        # gt("fixef, dropping")

        # We refactor the fixed effects => we may remove even more obs
        info_fe = setup_fixef(fixef_df = fixef_df, lhs = lhs, fixef_vars = fixef_vars,
                              fixef.rm = fixef.rm, family = family, isSplit = FALSE,
                              origin_type = origin_type, isSlope = isSlope, slope_flag = slope_flag,
                              slope_df = slope_df, slope_vars_list = slope_vars_list,
                              fixef_names_old = fixef_names_old, fixef_sizes = fixef_sizes,
                              obs2keep = obs2keep, nthreads = nthreads)

        # gt("fixef, recompute")

        Q               = info_fe$Q
        fixef_id        = info_fe$fixef_id
        fixef_names     = info_fe$fixef_names
        fixef_sizes     = info_fe$fixef_sizes
        fixef_table     = info_fe$fixef_table
        sum_y_all       = info_fe$sum_y_all
        lhs             = info_fe$lhs

        obs2remove      = info_fe$obs2remove
        fixef_removed   = info_fe$fixef_removed
        message_fixef   = info_fe$message_fixef

        slope_variables = info_fe$slope_variables
        slope_flag      = info_fe$slope_flag

        fixef_id_res    = info_fe$fixef_id_res
        fixef_sizes_res = info_fe$fixef_sizes_res
        new_order       = info_fe$new_order

        if(length(obs2remove) > 0){

            if(is.null(obs2keep)){
                obs2keep = 1:length(fixef_df[[1]])
            }

            obs2keep = obs2keep[-obs2remove]
        }

        lhs_done = TRUE

        assign_fixef_env(new_env, family, origin_type, fixef_id, fixef_sizes, fixef_table,
                         sum_y_all, slope_flag, slope_variables, slope_vars_list)

        assign_fixef = TRUE

    } else if(isFixef){
        slope_flag = get("slope_flag", env)
        isSlope = any(slope_flag != 0)
    }


    # This is very tedious
    if(!is.null(obs2keep) || assign_fixef){
        # Recreating the values

        #
        # RM obs = TRUE ####
        #

        if(!is.null(obs2keep)){
            # I set both values to FALSE since they're done here
            save_lhs = save_rhs = FALSE

            #
            # The left hand side
            #

            if(assign_lhs){
                if(lhs_done == FALSE){
                    # lhs: can be vector or list
                    lhs = select_obs(lhs, obs2keep)
                }
                assign("lhs", lhs, new_env)
            }

            #
            # The linear mat
            #

            if(assign_rhs){

                if(is.null(rhs)){
                    rhs = get("linear.mat", env)
                }

                isLinear = FALSE
                if(length(rhs) > 1){
                    isLinear = TRUE
                    rhs = select_obs(rhs, obs2keep, nthreads)

                    linear.params = colnames(rhs)
                    nonlinear.params = get("nonlinear.params", env)
                    params = c(nonlinear.params, linear.params)
                    assign("linear.params", linear.params, new_env)
                    assign("params", params, new_env)

                    # useful when feNmlm or feglm
                    if(origin_type %in% c("feNmlm", "feglm")){

                        start = get("start", env)
                        if(!is.null(start)){
                            new_start = setNames(start[params], params)
                            new_start[is.na(new_start)] = 0
                            assign("start", new_start, new_env)
                        }

                        if(origin_type == "feNmlm"){
                            upper = get("upper", env)
                            new_upper = setNames(upper[params], params)
                            new_upper[is.na(new_upper)] = Inf
                            assign("upper", new_upper, new_env)

                            lower = get("lower", env)
                            new_lower = setNames(lower[params], params)
                            new_lower[is.na(new_lower)] = -Inf
                            assign("lower", new_lower, new_env)
                        }
                    }

                    assign("linear.mat", rhs, new_env)
                }
                assign("isLinear", isLinear, new_env)
            }

            #
            # The IV
            #

            do_iv = get("do_iv", env)
            if(do_iv){
                # LHS
                iv_lhs = get("iv_lhs", env)
                for(i in seq_along(iv_lhs)){
                    iv_lhs[[i]] = iv_lhs[[i]][obs2keep]
                }
                assign("iv_lhs", iv_lhs, new_env)

                # RHS
                iv.mat = get("iv.mat", env)
                iv.mat = select_obs(iv.mat, obs2keep, nthreads, msg = "instrument")
                assign("iv.mat", iv.mat, new_env)

            }

            #
            # The cluster
            #

            do_summary = get("do_summary", env)
            if(do_summary){
                cluster = get("cluster", env)
                if(is.data.frame(cluster)){
                    cluster = cluster[obs2keep, , drop = FALSE]
                    assign("cluster", cluster, new_env)
                }
            }

            #
            # Stepwise RHS
            #

            if(get("do_multi_rhs", env)){
                linear_core = get("linear_core", env)
                rhs_sw = get("rhs_sw", env)

                assign("linear_core", select_obs(linear_core, obs2keep, nthreads), new_env)
                assign("rhs_sw", select_obs(rhs_sw, obs2keep, nthreads), new_env)
            }

            #
            # New data if stepwise estimation
            #

            # So far this is used only in stepwise FE
            if(get("do_multi_fixef", env)){
                data = get("data", env)
                assign("data", data[obs2keep, , drop = FALSE], new_env)
            }

            #
            # The non-linear part, the weight and the offset
            #

            offset.value = get("offset.value", env)
            isOffset = length(offset.value) > 1
            if(isOffset){
                assign("offset.value", offset.value[obs2keep], new_env)
            }

            weights.value = get("weights.value", env)
            isWeight = length(weights.value) > 1
            if(isWeight){
                # used in family$family_equiv
                weights.value = weights.value[obs2keep]
                assign("weights.value", weights.value, new_env)
            }

            isNL = get("isNL", env)
            if(isNL){
                envNL = get("envNL", env)
                envNL_new = new.env()
                for(var in names(envNL)){
                    x = get(var, envNL)
                    if(length(x) > 1){
                        assign(var, x[obs2keep], envNL_new)
                    } else {
                        assign(var, x, envNL_new)
                    }
                }
                assign("envNL", envNL_new, new_env)
            }


            # GLM specific stuff
            if(origin_type == "feglm"){

                family_funs = get("family_funs", env)
                if(!is.list(lhs) && family_funs$family_equiv == "poisson"){

                    y_pos = lhs[lhs > 0]
                    qui_pos = lhs > 0

                    if(isWeight){
                        constant = sum(weights.value[qui_pos] * y_pos * cpppar_log(y_pos, nthreads) - weights.value[qui_pos] * y_pos)
                        sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(wt[qui_pos] * y_pos * eta[qui_pos]) + sum(wt * mu))
                    } else {
                        constant = sum(y_pos * cpppar_log(y_pos, nthreads) - y_pos)
                        sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(y_pos * eta[qui_pos]) + sum(mu))
                    }

                    family_funs$sum_dev.resids = sum_dev.resids

                    assign("family_funs", family_funs, new_env)
                }

                starting_values = get("starting_values", env)
                if(!is.null(starting_values)){
                    assign("starting_values", starting_values[obs2keep], new_env)
                }

            }

            #
            # Reformatting "res"
            #

            # offset and weight
            if(isOffset){
                res$offset = offset.value
            }
            if(isWeight){
                res$weights = weights.value
            }

            # We save lhs in case of feglm, for later within-r2 evaluation
            if(origin_type == "feglm" && isFixef){
                res$y = lhs
            }

            # obs2keep
            res$nobs = length(obs2keep)
            assign("nobs", res$nobs, new_env)
            if(is.null(res$obs_selection)){
                res$obs_selection = list(obs2keep)
            } else {
                res$obs_selection[[length(res$obs_selection) + 1]] = obs2keep
            }

        }

        # Nber of parameters
        params = get("params", new_env)
        K = length(params)
        if(isFixef){
            fixef_sizes = get("fixef_sizes", new_env)
            slope_flag = get("slope_flag", new_env)
            isSlope = any(slope_flag != 0)
            if(isSlope){
                K = K + sum(fixef_sizes * (1 + abs(slope_flag) - (slope_flag < 0)))
                # now the references
                if(any(slope_flag >= 0)){
                    K = K + 1 - sum(slope_flag >= 0)
                }
            } else {
                K = K + sum(fixef_sizes - 1) + 1
            }
        }
        res$nparams = K

        # Fixed-effects information
        if(isFixef){
            if(isSlope){
                res$slope_flag = slope_flag[order(new_order)]
                res$slope_flag_reordered = slope_flag
                res$slope_variables_reordered = slope_variables
                res$fe.reorder = new_order
            }

            res$fixef_id = fixef_id_res
            res$fixef_sizes = fixef_sizes_res
        }

    }

    #
    # Save LHS/RHS ####
    #

    if(save_lhs){
        # Here lhs is ALWAYS a vector
        assign("lhs", lhs, new_env)

        if(origin_type == "feglm"){
            family_funs = get("family_funs", env)

            if(family_funs$family_equiv == "poisson"){

                y_pos = lhs[lhs > 0]
                qui_pos = lhs > 0

                weights.value = get("weights.value", env)
                isWeight = length(weights.value) > 1

                if(isWeight){
                    constant = sum(weights.value[qui_pos] * y_pos * cpppar_log(y_pos, nthreads) - weights.value[qui_pos] * y_pos)
                    sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(wt[qui_pos] * y_pos * eta[qui_pos]) + sum(wt * mu))
                } else {
                    constant = sum(y_pos * cpppar_log(y_pos, nthreads) - y_pos)
                    sum_dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(y_pos * eta[qui_pos]) + sum(mu))
                }

                family_funs$sum_dev.resids = sum_dev.resids

                assign("family_funs", family_funs, new_env)
            }
        }

        if(origin_type == "feglm" && isFixef){
            res$y = lhs
        }
    }

    if(save_rhs){
        assign("linear.mat", rhs, new_env)

        isLinear = FALSE
        if(length(rhs) > 1){
            isLinear = TRUE

            linear.params = colnames(rhs)
            nonlinear.params = get("nonlinear.params", env)
            params = c(nonlinear.params, linear.params)
            assign("linear.params", linear.params, new_env)
            assign("params", params, new_env)

            # useful when feNmlm or feglm
            if(origin_type %in% c("feNmlm", "feglm")){

                start = get("start", env)
                if(!is.null(start)){
                    new_start = setNames(start[params], params)
                    new_start[is.na(new_start)] = 0
                    assign("start", new_start, new_env)
                }

                if(origin_type == "feNmlm"){
                    upper = get("upper", env)
                    new_upper = setNames(upper[params], params)
                    new_upper[is.na(new_upper)] = Inf
                    assign("upper", new_upper, new_env)

                    lower = get("lower", env)
                    new_lower = setNames(lower[params], params)
                    new_lower[is.na(new_lower)] = -Inf
                    assign("lower", new_lower, new_env)
                }
            }
        }

        assign("isLinear", isLinear, new_env)

        # Finally the DoF
        params = get("params", new_env)
        K = length(params)
        if(isFixef){
            fixef_sizes = get("fixef_sizes", new_env)
            slope_flag = get("slope_flag", new_env)
            if(isSlope){
                K = K + sum(fixef_sizes * (1 + abs(slope_flag) - (slope_flag < 0)))
                # now the references
                if(any(slope_flag >= 0)){
                    K = K + 1 - sum(slope_flag >= 0)
                }
            } else {
                K = K + sum(fixef_sizes - 1) + 1
            }
        }
        res$nparams = K

        if(!isLinear && !isTRUE(res$iv)){
            res$onlyFixef = TRUE
        }
    }

    # We save the formulas
    if(!is.null(fml_linear)){
        res$fml = res$fml_all$linear = fml_linear
    }

    if(!is.null(fml_fixef)){
        res$fml_all$fixef = fml_fixef
    }

    if(!is.null(fml_iv_endo)){
        fml_linear = res$fml_all$linear
        fml_iv = res$fml_all$iv

        is_pblm = "try-error" %in% class(try(str2lang(fml_iv_endo), silent = TRUE))
        if(is_pblm){
            # Avoids a bug
            fml_iv_endo = paste0("`", fml_iv_endo, "`")
        }

        new_fml_linear = .xpd(..y ~ ..inst + ..rhs,
                             ..y = fml_iv_endo, ..rhs = fml_linear[[3]], ..inst = fml_iv[[3]])
        res$fml = res$fml_all$linear = new_fml_linear
        res$fml_all$iv = NULL
    }

    assign("fixest_env", TRUE, new_env)
    assign("res", res, new_env)

    return(new_env)
}


#' Stepwise estimation tools
#'
#' Functions to perform stepwise estimations.
#'
#' @param ... Represents formula variables to be added in a stepwise fashion to an estimation.
#'
#' @details
#'
#' To include multiple independent variables, you need to use the stepwise functions. There are 4 stepwise functions: sw, sw0, csw, csw0. Let's explain that.
#'
#' Assume you have the following formula: \code{fml = y ~ x1 + sw(x2, x3)}. The stepwise function \code{sw} will estimate the following two models: \code{y ~ x1 + x2} and \code{y ~ x1 + x3}. That is, each element in \code{sw()} is sequentially, and separately, added to the formula. Would have you used \code{sw0} in lieu of \code{sw}, then the model \code{y ~ x1} would also have been estimated. The \code{0} in the name implies that the model without any stepwise element will also be estimated.
#'
#' Finally, the prefix \code{c} means cumulative: each stepwise element is added to the next. That is, \code{fml = y ~ x1 + csw(x2, x3)} would lead to the following models \code{y ~ x1 + x2} and \code{y ~ x1 + x2 + x3}. The \code{0} has the same meaning and would also lead to the model without the stepwise elements to be estimated: in other words, \code{fml = y ~ x1 + csw0(x2, x3)} leads to the following three models: \code{y ~ x1}, \code{y ~ x1 + x2} and \code{y ~ x1 + x2 + x3}.
#'
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # Regular stepwise
#' feols(y ~ sw(x1, x2, x3), base)
#'
#' # Cumulative stepwise
#' feols(y ~ csw(x1, x2, x3), base)
#'
#' # Using the 0
#' feols(y ~ x1 + x2 + sw0(x3), base)
#'
#' @name stepwise
NULL

sw = csw = function(...){
    mc = match.call(expand.dots = TRUE)

    n = length(mc) - 1

    if(n == 0){
        return("")
    }

    res = vector("character", n)

    for(i in 1:n){
        res[[i]] = deparse_long(mc[[i + 1]])
    }

    if(grepl("^c", as.character(mc[[1]]))){
        attr(res, "is_cumul") = TRUE
    }

    res
}

csw0 = sw0 = function(...){
    mc = match.call(expand.dots = TRUE)

    n = length(mc) - 1

    if(n == 0){
        return("")
    }

    res = vector("character", n + 1)

    for(i in 1:n){
        res[[i + 1]] = deparse_long(mc[[i + 1]])
    }

    if(grepl("^c", as.character(mc[[1]]))){
        attr(res, "is_cumul") = TRUE
    }

    res
}

# To add a proper documentation:

#' @rdname stepwise
sw <- sw

#' @rdname stepwise
csw <- csw

#' @rdname stepwise
sw0 <- sw0

#' @rdname stepwise
csw0 <- csw0


extract_stepwise = function(fml, tms, all_vars = TRUE){
    # fml = y ~ sw(x1, x2)
    # fml = ~ fe1 + csw(fe2, fe3)

    is_fml = !missing(fml)
    if(is_fml){
        n_parts = length(fml)
        osf = n_parts == 2
        fml_txt = deparse_long(fml[[n_parts]])
        do_stepwise = grepl("(^|[^[:alnum:]_\\.])c?sw0?\\(", fml_txt)
        tl = attr(terms(fml), "term.labels")

    } else {
        osf = TRUE
        tl = attr(tms, "term.labels")
        tl = gsub("_impossible_var_name_", "^", tl, fixed = TRUE)
        do_stepwise = any(grepl("(^|[^[:alnum:]_\\.])c?sw0?\\(", tl))

    }

    if(do_stepwise){
        # We will partially create the linear matrix

        qui = grepl("(^|[^[:alnum:]_\\.])c?sw0?\\(", tl)
        if(sum(qui) >= 2){
            stop("You cannot use more than one stepwise function.")
        }

        if(!is_naked_fun(tl[qui], "c?sw0?")){
            stop("You cannot combine stepwise functions with any other element.")
        }

        tl_new = tl

        sw_text = tl[qui]
        sw_terms = eval(str2lang(sw_text))
        is_cumul = isTRUE(attr(sw_terms, "is_cumul"))

        if(length(sw_terms) == 1){

            if(nchar(sw_terms) == 0){
                tl_new = tl[!qui]
                if(length(tl_new) == 0){
                    tl_new = 1
                }

            } else {
                tl_new[qui] = sw_terms
            }

            # lhs_text: created in the LHS section
            if(osf){
                fml_new = .xpd(rhs = tl_new)
            } else {
                fml_new = .xpd(lhs = fml[[2]], rhs = tl_new)
            }

            if(is_fml){
                res = list(do_multi = FALSE, fml = fml_new)
            } else {
                res = list(do_multi = FALSE, tl = tl_new)
            }

            return(res)

        } else {

            tl_new[qui] = "..STEPWISE_VARIABLES"
            fml_raw = .xpd(rhs = tl_new)

            fml_all_full = list()
            fml_all_sw = list()
            for(i in seq_along(sw_terms)){
                if(nchar(sw_terms[i]) == 0 && length(tl_new) > 1){
                    # adding 1 when empty is "ugly"
                    fml_all_full[[i]] = .xpd(rhs = tl_new[!qui])
                    if(osf){
                        fml_all_sw[[i]] = ~ 1
                    } else {
                        fml_all_sw[[i]] = lhs ~ 1
                    }
                } else {
                    if(is_cumul){
                        fml_all_full[[i]] = .xpd(fml_raw, ..STEPWISE_VARIABLES = sw_terms[1:i])
                    } else {
                        fml_all_full[[i]] = .xpd(fml_raw, ..STEPWISE_VARIABLES = sw_terms[i])
                    }

                    if(osf){
                        fml_all_sw[[i]] = .xpd(rhs = sw_terms[i])
                    } else {
                        fml_all_sw[[i]] = .xpd(lhs = quote(lhs),  rhs = sw_terms[i])
                    }

                }
            }

            # New mechanism => we will evaluate everything at the first call for linear stepwise
            # we respect the order of the user. So we need the core left and right

            qui = which(qui)
            if(is_cumul){
                # We add the first terms to the main formula
                tl_left = tl[seq(1, length.out = qui)]
                tl_left[qui] = sw_terms[1]
            } else {
                tl_left = tl[seq(1, length.out = qui - 1)]
            }

            tl_right = tl[seq(qui + 1, length.out = length(tl) - qui)]

            if(length(tl_left) == 0) tl_left = ""
            if(length(tl_right) == 0) tl_right = ""

            if(osf){
                fml_core_left = .xpd(rhs = tl_left)
                fml_core_right = .xpd(rhs = tl_right)
            } else {
                fml_core_left = .xpd(lhs = quote(lhs), rhs = tl_left)
                fml_core_right = .xpd(lhs = quote(lhs), rhs = tl_right)
            }


            if(all_vars){
                if(is_fml){
                    sw_all_vars = all_vars_with_i_prefix(fml[[3]])
                } else {
                    # we need to recreate the formula because of impossible_var_name in tms
                    sw_all_vars = all_vars_with_i_prefix(.xpd(rhs = tl))
                }
            } else {
                sw_all_vars = all_vars_with_i_prefix(.xpd(rhs = sw_terms))
            }

            # Creating the current formula
            if(is_cumul){
                # Cumulative => first item in main linear matrix
                tl_new[qui] = sw_terms[1]
            } else {
                tl_new = tl[!qui]
                if(length(tl_new) == 0){
                    tl_new = 1
                }
            }

            # => this is only useful to deduce fake_intercept
            if(osf){
                fml_new = .xpd(rhs = tl_new)
            } else {
                fml_new = .xpd(lhs = quote(lhs), rhs = tl_new)
            }

        }
    } else if(is_fml) {
        res = list(do_multi = FALSE, fml = fml)
        return(res)
    } else {
        res = list(do_multi = FALSE, tl = tl)
        return(res)
    }

    res = list(do_multi = TRUE, fml = fml_new, fml_all_full = fml_all_full,
               fml_all_sw = fml_all_sw, is_cumul = is_cumul, fml_core_left = fml_core_left,
               fml_core_right = fml_core_right, sw_all_vars = sw_all_vars)

    return(res)
}



# a = list(1:5, 6:10)
# b = 1:5
# d = list(matrix(1:10, 5, 2), 1, matrix(1:5, 5, 1))
# select_obs(a, 1:2) ; select_obs(b, 1:2) ; select_obs(d, 1:2)
select_obs = function(x, index, nthreads = 1, msg = "explanatory variable"){
    # => selection of observations.
    # Since some objects can be of multiple types, this avoids code repetition and increases clarity.

    if(!is.list(x)){

        if(is.matrix(x)){
            x = x[index, , drop = FALSE]

            only_0 = cpppar_check_only_0(x, nthreads)
            if(all(only_0 == 1)){
                stop("After removing NAs (or perfect fit fixed-effects), not a single explanatory variable is different from 0.")

            } else if(any(only_0 == 1)){
                x = x[, only_0 == 0, drop = FALSE]
            }

        } else if(length(x) > 0){
            # We check the length because RHS = 1 means not linear (we don't want to subset on that)
            x = x[index]
        }

    } else if(is.data.frame(x)){
        x = x[index, , drop = FALSE]

    } else if(is.matrix(x[[1]]) || length(x[[1]]) == 1){
        # Means RHS
        for(i in seq_along(x)){
            if(length(x[[i]]) > 1){
                x[[i]] = x[[i]][index, , drop = FALSE]

                only_0 = cpppar_check_only_0(x[[i]], nthreads)
                if(all(only_0 == 1)){
                    stop("After removing NAs (or perfect fit fixed-effects), not a single ", msg, " is different from 0.")

                } else if(any(only_0 == 1)){
                    x[[i]] = x[[i]][, only_0 == 0, drop = FALSE]
                }
            }
        }

    } else {
        for(i in seq_along(x)){
            x[[i]] = x[[i]][index]
        }
    }

    x
}

collect_vars = function(...){
    # utility tool to collect the variables used in some arguments

    dots = list(...)
    vars = c()

    for(i in seq_along(dots)){
        di = dots[[i]]
        if("formula" %in% class(di)){
            vars = c(vars, all.vars(di))
        } else if(length(di) < 5 && is.character(di)){
            vars = c(vars, di)
        }
    }

    unique(vars)
}

fixest_NA_results = function(env){
    # Container for NA results
    # so far the non-linear part is not covered

    res = get("res", env)

    X = get("linear.mat", env)

    n = ncol(X)
    NA_values = rep(NA_real_, n)

    coef = se = NA_values
    names(coef) = names(se) = colnames(X)

    coeftable = data.frame("Estimate" = NA_values, "Std. Error" = NA_values, "t value" = NA_values, "Pr(>|t|)" = NA_values)
    names(coeftable) = c("Estimate", "Std. Error", "t value",  "Pr(>|t|)")
    row.names(coeftable) = names(coef)

    cov.scaled = matrix(NA, n, n)

    attr(se, "type") = attr(coeftable, "type") = attr(cov.scaled, "type") = "Void"
    res$coeftable = coeftable
    res$se = se
    res$cov.scaled = cov.scaled

    res$summary = TRUE

    # Fit stats
    res$nobs = nrow(X)
    res$ssr = res$ssr_null = res$ssr_fe_only = res$sigma2 = res$loglik = res$ll_null = res$ll_fe_only = res$pseudo_r2 = res$deviance = res$sq.cor = NA_real_

    res$NA_model = TRUE
    class(res) = "fixest"

    res
}


