#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Jun 10 23:46:47 2019
# Purpose: sets up the params + performs
#           all necessay checks
#----------------------------------------------#


fixest_env <- function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml,
                       fixef, NL.start, lower, upper, NL.start.init,
                       offset, subset, split = NULL, fsplit = NULL, linear.start = 0, jacobian.method = "simple",
                       useHessian = TRUE, hessian.args = NULL, opt.control = list(),
                       y, X, fixef_mat, panel.id, fixef.rm = "perfect",
                       nthreads = getFixest_nthreads(),
                       verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000, collin.tol = 1e-14,
                       deriv.iter = 5000, deriv.tol = 1e-4, glm.iter = 25, glm.tol = 1e-8,
                       etastart, mustart,
                       warn = TRUE, notes = getFixest_notes(), combine.quick, demeaned = FALSE,
                       origin_bis, origin = "feNmlm", mc_origin, mc_origin_bis, mc_origin_ter,
                       computeModel0 = FALSE, weights,
                       debug = FALSE, mem.clean = FALSE, ...){

    # INTERNAL function:
    # the estimation functions need input data in the exact format without any mistake possible (bc of c++)
    # this function takes care of implementing all the checks needed and providing the proper formatting
    # it also prepare the list returned in the est. funs

    # Only an environment is returned
    env <- new.env(parent = emptyenv())

    # the function of origin
    if(!missnull(origin_bis)){
        origin = origin_bis
    }

    # The match.call of origin
    if(!missnull(mc_origin_bis)){
        mc_origin = mc_origin_bis
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
    main_args = c("fml", "data", "panel.id", "offset", "subset", "split", "fsplit", "fixef.rm", "fixef.tol", "fixef.iter", "fixef", "nthreads", "verbose", "warn", "notes", "combine.quick", "start", "only.env", "mem.clean")
    femlm_args = c("family", "theta.init", "linear.start", "opt.control", "deriv.tol", "deriv.iter")
    feNmlm_args = c("NL.fml", "NL.start", "lower", "upper", "NL.start.init", "jacobian.method", "useHessian", "hessian.args")
    feglm_args = c("family", "weights", "glm.iter", "glm.tol", "etastart", "mustart", "collin.tol")
    feols_args = c("weights", "demeaned", "collin.tol")
    internal_args = c("debug")

    deprec_old_new = c()

    common_args = c(main_args, internal_args, deprec_old_new)

    # FIT methods
    feglm.fit_args = c(setdiff(common_args, c("fml", "data")), feglm_args, "y", "X", "fixef_mat")

    args = names(mc_origin)
    args = args[nchar(args) > 0]

    # Checking the args
    valid_args = switch(origin,
                        feols = c(common_args, feols_args),
                        feglm = c(common_args, feglm_args),
                        feglm.fit = feglm.fit_args,
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
    # IV estimation
    #

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
    # Formattting + checks
    #

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
            family <- get(family, mode = "function", envir = parent.frame(2))
        }

        if(is.function(family)) {
            family <- family()
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
    # Formatting data ####
    #

    fml_no_xpd = NULL # will be returned if expansion is performed
    isPanel = FALSE
    if(isFit){
        isFixef = !missnull(fixef_mat)
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
        check_value(fml, "ts formula", .arg_name = "fml")

        if(!"formula" %in% class(fml)) stop("The argument 'fml' must be a formula.")
        fml = formula(fml) # we regularize the formula to check it
        # if(length(fml) != 3) stop("The formula must be two sided: e.g. y~x1+x2, or y~x1+x2|fe1+fe2.")

        # We apply expand for macros => we return fml_no_xpd
        if(length(getFixest_fml()) > 0){
            fml_no_xpd = fml
            fml = xpd(fml)
        }

        #
        # ... Panel setup ####
        #

        if(check_lag(fml)){
            panel.info = NULL
            isPanel = TRUE
            if(!is.null(attr(data, "panel_info"))){
                if(!missing(panel.id)){
                    warning("The argument 'panel.id' is provided but argument 'data' is already a 'fixest_panel' object. Thus the argument 'panel.id' is ignored.", immediate. = TRUE)
                }

                panel__meta__info = attr(data, "panel_info")
                panel.id = panel__meta__info$panel.id
                panel.info = panel__meta__info$call
            } else {
                # Later: automatic deduction using the first two clusters
                if(missing(panel.id)){
                    stop("To use lag/leads (with l()/f()): either provide the argument 'panel.id' with the panel identifiers OR set your data as a panel with function panel().")
                }
                panel__meta__info = panel_setup(data, panel.id, from_fixest = TRUE)
            }
            class(data) = "data.frame"
        }

        FML = Formula::Formula(fml)
        n_rhs = length(FML)[2]

        if(n_rhs > 2){
            stop("The argument 'fml' cannot contain more than two parts separated by a pipe ('|').")
        }

        # We test there are no problems in the formula
        info = try(terms(formula(FML, lhs = 1, rhs = 1)), silent = TRUE)
        if("try-error" %in% class(info)){
            dp = deparse(formula(FML, lhs = 1, rhs = 1))
            stop("The formula: ", dp[1], ifsingle(dp, "", "..."), ", is not valid:\n", gsub("^[^\n]+\n", "", info))
        }

        # we check the FE part
        if(n_rhs == 2){
            info = terms_fixef(formula(FML, lhs = 0, rhs = 2))
            if("try-error" %in% class(info)){
                dp = deparse(formula(FML, lhs = 0, rhs = 2))
                stop("The fixed-effects part of the formula: ", dp[1], ifsingle(dp, "", "..."), ", is not valid:\n", gsub("^[^\n]+\n", "", info))
            }
        }

        if(isPanel){
            fml = try(rewrite_fml(fml), silent = TRUE)
            if("try-error" %in% class(fml)){
                stop("Problem in the formula regarding lag/leads: ", gsub("__expand", "", fml))
            }
            FML = Formula::Formula(fml)
        }

        # for clarity, arg fixef is transformed into fixef_vars
        fixef_vars = NULL
        isFixef = FALSE
        if(n_rhs == 2){
            isFixef = TRUE
            if(missnull(fixef)){
                fixef_vars = formula(FML, lhs = 0, rhs = 2)
                fml = formula(FML, lhs = 1, rhs = 1)
            } else {
                stop("To add fixed-effects: either include them in argument 'fml' using a pipe ('|'), either use the argument 'fixef'. You cannot use both!")
            }
        } else if(!missing(fixef) && length(fixef) >= 1){
            isFixef = TRUE

            # We check the argument => character vector

            check_arg_plus(fixef, "multi match", .choices = dataNames, .message = "Argument 'fixef', when provided, must be a character vector of variable names.")

            fixef_vars = fixef
        }

        # fml_full now does not contain the FEs, but will later contain them
        fml_full = fml
    }

    #
    # ... subset ####
    #

    nobs_origin = NULL
    if(isFit){
        # We have to first eval y if isFit in order to check subset properly
        if(missing(y)){
            stop("You must provide argument 'y' when using ", origin_type, ".fit.")
        }

        lhs = check_value_plus(y, "numeric vector conv")

        # we reconstruct a formula
        fml = as.formula(paste0(deparse_long(mc_origin[["y"]]), "~1"))

        nobs_origin = nobs = length(lhs)
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
                subset = check_value_plus(subset.value, "evalset integer vector gt{0} | logical vector len(data)", .data = data, .prefix = "In argument 'subset', the expression")

            } else {

                if(isFit){
                    check_value(subset, "integer vector gt{0} | logical vector len(data)", .data = lhs)
                } else {
                    check_value(subset, "integer vector gt{0} | logical vector len(data)",
                                .prefix = "If not a formula, argument 'subset'", .data = data)
                }
            }

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

            if(isFit){
                delayed.subset = TRUE
                lhs = lhs[subset]
            } else {
                nobs_origin = NROW(data)

                # subsetting creates a deep copy. We avoid copying the entire data set.
                var2keep = all.vars(fml)
                if(isFixef){
                    if(is.character(fixef_vars)){
                        var2keep = c(var2keep, fixef_vars)
                    } else {
                        var2keep = c(var2keep, all.vars(fixef_vars))
                    }
                }
                if(!missnull(NL.fml)){
                    check_value(NL.fml, "os formula")
                    var2keep = c(var2keep, all.vars(NL.fml))
                }
                var2keep = intersect(unique(var2keep), names(data))

                data = data[subset, var2keep, drop = FALSE]
            }


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

    # the LHS for isFit has been done just before subset

    # evaluation
    multi_lhs = FALSE
    if(isFit == FALSE){

        # The LHS must contain only values in the DF
        namesLHS = all.vars(fml[[2]])
        if(length(namesLHS) == 0){
            stop("The right hand side of the formula (", deparse_long(fml[[2]]), ") contains no variable!")

        } else if(!all(namesLHS %in% dataNames)){
            not_there = namesLHS[!namesLHS %in% dataNames]
            stop(ifsingle(not_there, "The v", "V"), "ariable", enumerate_items(not_there, "s.is"), " in the LHS of the formula but not in the dataset.")
        }

        lhs_text = gsub("^(c|(c?(stepwise|sw)0?))\\(", "list(", deparse_long(fml[[2]]))

        lhs = try(eval(parse(text = lhs_text), data), silent = TRUE)

        if("try-error" %in% class(lhs)){
            stop("Evaluation of the left-hand-side (equal to ", deparse_long(fml[[2]]), ") raises an error: \n", lhs)

        } else if(is.list(lhs)){
            # we check the consistency
            lhs_names = eval(parse(text = gsub("^list\\(", "stepwise(", lhs_text)))

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
            lhs = check_value_plus(lhs, "numeric vector conv", .prefix = "The left hand side")
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

    anyNA_y = FALSE
    if(multi_lhs){
        msgNA_y = "LHS: --"
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
        msgNA_y = "LHS: 0"
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
            stop("The dependent variable is a constant. Estimation cannot be done.")
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

    fixef_terms_full = NULL
    interaction.info = NULL
    fake_intercept = FALSE
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

            # we coerce to matrix (avoids sparse matrix) + transform into double
            linear.mat = as.matrix(X) * 1

            if(is.null(colnames(linear.mat))){
                colnames(linear.mat) = paste0("X", 1:ncol(linear.mat))
            }

            # The formula
            fml = update(fml, as.formula(paste0(".~", paste0(colnames(linear.mat), collapse = "+"))))
            fml_full = fml

            linear.varnames = NULL

            if(delayed.subset){
                linear.mat = linear.mat[subset, , drop = FALSE]
            }

        }

    } else {
        isLinear = FALSE
        options("fixest_interaction_ref" = NULL)

        rhs_txt = deparse_long(fml[[3]])
        if(grepl("[^:]::[^:]", rhs_txt)){

            new_fml = try(interact_fml(fml), silent = TRUE)
            if("try-error" %in% class(new_fml)){
                stop("Error in the right-hand-side of the formula: ", new_fml)
            }
            linear.varnames = all.vars(new_fml[[3]])
        } else {
            linear.varnames = all.vars(fml[[3]])
        }

        if(length(linear.varnames) > 0 || attr(terms(fml), "intercept") == 1){
            isLinear = TRUE
        }

        if(isLinear){

            if(!all(linear.varnames %in% dataNames)){
                not_there = setdiff(linear.varnames, dataNames)
                stop(ifsingle(not_there, "The v", "V"), "ariable", enumerate_items(not_there, "s.is"), " in the RHS of the formula but not in the dataset.")
            }

            if(isFixef){
                # if dummies are provided, we make sure there is an
                # intercept so that factors can be handled properly
                #
                # Beware: only slope: we don't add it

                if(is.character(fixef_vars)){
                    # There is no varying slope allowed
                    fake_intercept = TRUE
                } else {
                    # On en profite pour extraire les termes
                    fixef_terms_full = terms_fixef(fixef_vars)
                    if("try-error" %in% class(fixef_terms_full)){
                        stop("Problem extracting the terms of the fixed-effects part of the formula:\n", fixef_terms_full)
                    }
                    fake_intercept = any(fixef_terms_full$slope_flag >= 0)

                }
            }

            #
            # We construct the linear matrix
            #

            if(grepl("(^|[^[:alnum:]_\\.])c?(stepwise|sw)0?\\(", rhs_txt)){
                # We will partially create the linear matrix

                t_fml = terms(fml)
                tl = attr(t_fml, "term.labels")

                qui = grepl("(^|[^[:alnum:]_\\.])c?(stepwise|sw)0?\\(", tl)
                if(sum(qui) >= 2){
                    stop("You cannot use more than one stepwise function.")
                }

                if(!is_naked_fun(tl[qui], "c?(stepwise|sw)0?")){
                    stop("You cannot combine stepwise functions with any other element.")
                }

                tl_new = tl

                sw_text = tl[qui]
                sw_terms = eval(parse(text = sw_text))
                is_cumul = isTRUE(attr(sw_terms, "is_cumul"))

                lhs_txt = deparse_long(fml[[2]])

                if(length(sw_terms) == 1){

                    if(nchar(sw_terms) == 0){
                        tl_new = tl[!qui]
                        if(length(tl_new) == 0){
                            tl_new = 1
                        }

                    } else {
                        tl_new[qui] = sw_terms
                    }

                    fml = as.formula(paste0(lhs_txt, "~", paste(tl_new, collapse = " + ")))

                } else {
                    multi_rhs = TRUE
                    tl_new[qui] = "..STEPWISE_VARIABLES"
                    fml_raw = xpd(~..rhs, ..rhs = tl_new)

                    fml_rhs_all = list()
                    fml_all_sw = list()
                    for(i in seq_along(sw_terms)){
                        if(nchar(sw_terms[i]) == 0 && length(tl_new) > 1){
                            # adding 1 when empty is "ugly"
                            fml_rhs_all[[i]] = xpd(~..rhs, ..rhs = tl_new[!qui])
                            fml_all_sw[[i]] = lhs ~ 1
                        } else {
                            if(is_cumul){
                                fml_rhs_all[[i]] = xpd(fml_raw, ..STEPWISE_VARIABLES = sw_terms[1:i])
                            } else {
                                fml_rhs_all[[i]] = xpd(fml_raw, ..STEPWISE_VARIABLES = sw_terms[i])
                            }

                            fml_all_sw[[i]] = xpd(lhs ~ ..STEPWISE_VARIABLES, ..STEPWISE_VARIABLES = sw_terms[i])
                        }
                    }

                    sw_all_vars = all.vars(xpd(~ ..sw, ..sw = sw_terms))

                    # Creating the current formula
                    if(is_cumul){
                        # Cumulative => first item in main linear matrix
                        tl_new[qui] = sw_terms[1]
                        sw_terms = sw_terms[-1]
                    } else {
                        tl_new = tl[!qui]
                        if(length(tl_new) == 0){
                            tl_new = 1
                        }
                    }

                    # Plus tard quand j'aurais develope mes donnes en fixest_semi_sparse => evaluer tous les termes
                    fml = xpd(lhs ~ ..rhs, ..rhs = tl_new)
                    # fml = as.formula(paste0(lhs_txt, "~", paste(tl_new, collapse = " + ")))

                }

            }

            linear.mat = try(fixest_model_matrix(fml, data, fake_intercept), silent = TRUE)
            if("try-error" %in% class(linear.mat)){
                stop("Evaluation of the right-hand-side of the formula raises an error: ", linear.mat)
            }

            useModel.matrix = attr(linear.mat, "useModel.matrix")
            attr(linear.mat, "useModel.matrix") = NULL

            # Interaction information => if no interaction: NULL
            interaction.info = getOption("fixest_interaction_ref")


        }

    }


    # further controls (includes na checking)
    msgNA_L = ""
    if(isLinear){

        linear.params <- colnames(linear.mat)
        anyNA_L = FALSE
        info = cpppar_which_na_inf_mat(linear.mat, nthreads)
        if(info$any_na_inf){
            anyNA_L = TRUE

            if(info$any_na) ANY_NA = TRUE
            if(info$any_inf) ANY_INF = TRUE

            isNA_L = info$is_na_inf
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_L
            msgNA_L = paste0(", RHS: ", numberFormatNormal(sum(isNA_L)))

            if(mem.clean){
                rm(isNA_L)
            }

        } else {
            msgNA_L = ", RHS: 0"
        }

        if(mem.clean){
            rm(info)
            gc()
            gc2trig = FALSE
        }

    } 	else {
        linear.params <- linear.start <- linear.varnames <- NULL
        useModel.matrix = FALSE
    }


    #
    # ... nonlinear part ####
    #

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
            msgNA_NL = paste0(", NL: ", numberFormatNormal(sum(isNA_NL)))

            if(mem.clean){
                rm(isNA_NL)
                gc2trig = TRUE
            }

        } else {
            msgNA_NL = ", NL: 0"
        }

    } else {
        isNonLinear = FALSE
        nl.call = 0
        allnames = nonlinear.params = nonlinear.varnames = character(0)
    }

    params <- c(nonlinear.params, linear.params)
    lparams <- length(params)
    varnames <- c(nonlinear.varnames, linear.varnames)

    # Attention les parametres non lineaires peuvent etre vides
    if(length(nonlinear.params)==0) isNL = FALSE
    else isNL = TRUE

    #
    # ... Offset ####
    #

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
                check_value_plus(offset.value, "evalset numeric vector conv", .data = data, .prefix = "In argument 'offset', the expression")

            } else {

                check_value_plus(offset, "numeric vector conv", .prefix = "If not a formula, argument 'offset'")

                if(length(offset) == 1){
                    offset.value = rep(offset, nobs)
                } else if(length(offset) != nobs){
                    stop("The offset's length should be equal to the data's length (currently it's ", numberFormatNormal(length(offset)), " instead of ", numberFormatNormal(nobs), ").")
                } else {
                    offset.value = offset
                }
            }

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
                msgNA_offset = paste0(", Offset: ", numberFormatNormal(sum(isNA_offset)))

                if(mem.clean){
                    rm(isNA_offset)
                }
            } else {
                msgNA_offset = ", Offset: 0"
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
                check_value_plus(weights.value, "evalset numeric vector conv", .data = data, .prefix = "In argument 'weights', the expression")


            } else {

                check_value_plus(weights, "numeric vector conv", .prefix = "If not a formula, argument 'weights'")

                if(length(weights) == 1){
                    if(weights == 1){
                        isWeight = FALSE
                        weights.value = 1
                    } else {
                        # No point in having all obs the same weight... yet...
                        weights.value = rep(weights, nobs)
                    }
                } else if(length(weights) != nobs){
                    stop("The weights's length should be equal to the data's length (currently it's ", numberFormatNormal(length(weights)), " instead of ", numberFormatNormal(nobs), ").")
                } else {
                    weights.value = weights
                }
            }

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
                msgNA_weight = paste0(", Weights: ", numberFormatNormal(sum(isNA_W)))

                if(mem.clean){
                    rm(isNA_W)
                }
            } else {
                msgNA_weight = ", Weights: 0"
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
            if((grepl("[[", dp, fixed = TRUE) || grepl("$", dp, fixed = TRUE)) && dp != 'x[["weights"]]'){
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
            split = check_value_plus(split.value, "evalset vector", .data = data, .prefix = "In argument 'split', the expression")

        } else {

            if(isFit){
                check_value(split, "vector len(value)", .value = nobs_origin)
            } else {
                check_value(split, "vector len(data)", .data = data, .prefix = "If not a formula, argument 'split'")
            }
        }

        if(delayed.subset){
            split = split[subset]
        }

        isSplit = TRUE
        split = to_integer(split, add_items = TRUE, sorted = TRUE, items.list = TRUE)
        split.items = split$items
        split = split$x

        anyNA_split = FALSE
        isNA_S = is.na(split)
        if(any(isNA_S)){
            anyNA_split = TRUE

            if(info$any_na) ANY_NA = TRUE

            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_S
            msgNA_split = paste0(", split: ", numberFormatNormal(sum(isNA_S)))

            if(mem.clean){
                rm(isNA_S)
            }
        } else {
            msgNA_split = ", split: 0"
        }
    }

    if(mem.clean && gc2trig){
        gc()
        gc2trig = FALSE
    }


    #
    # Handling Fixed-effects ####
    #

    isSlope = onlySlope = FALSE
    if(isFixef){
        # The main fixed-effects construction

        if(isFit){

            #
            # ... From fit ####
            #

            if(isVector(fixef_mat)){
                fixef_mat = data.frame(x = fixef_mat, stringsAsFactors = FALSE)
                names(fixef_mat) = deparse(mc_origin[["fixef_mat"]])[1]

            } else if(is.list(fixef_mat)){
                all_len = lengths(fixef_mat)
                if(any(diff(all_len) != 0)){
                    stop("The lengths of the vectors in fixef_mat differ (currently it is: ", enumerate_items(all_len), ").")
                }
                fixef_mat = as.data.frame(fixef_mat)
            }

            if(!is.matrix(fixef_mat) && !"data.frame" %in% class(fixef_mat)){
                stop("Argument fixef_mat must be a vector, a matrix, a list or a data.frame (currently its class is ", enumerate_items(class(fixef_mat)), ").")
            }

            if(is.matrix(fixef_mat) && is.null(colnames(fixef_mat))){
                colnames(fixef_mat) = paste0("fixef_", 1:ncol(fixef_mat))
                fixef_mat = as.data.frame(fixef_mat)
            }

            if(nrow(fixef_mat) != nobs){
                stop("The number of observations of fixef_mat (", nrow(fixef_mat), ") must match the length of y (", nobs, ").")
            }

            fixef_vars = names(fixef_mat)

            # The formula
            # fml = update(fml, as.formula(paste0(".~.|", paste0(fixef_vars, collapse = "+"))))
            fml_char = as.character(fml)
            fml_full = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", paste0(fixef_vars, collapse = "+")))

            if(delayed.subset){
                fixef_mat = fixef_mat[subset, , drop = FALSE]
            }

        } else {
            #
            # ... Regular ####
            #

            # Regular way
            if(is.character(fixef_vars) && any(grepl("^", fixef_vars, fixed = TRUE))){
                # we make it a formula
                fixef_vars = as.formula(paste0("~", paste0(fixef_vars, collapse = "+")))
            }

            if(is.character(fixef_vars)){
                # fixef_vars are valid fixef names
                fixef_mat = data[, fixef_vars, drop = FALSE]

            } else if(!is.character(fixef_vars)){
                # means the fixef_vars is a formula
                fixef_fml = fixef_vars

                # we check that the fixef variables are indeed in the data
                fixef_vars_all = all.vars(fixef_fml)
                if(!all(fixef_vars_all %in% dataNames)){
                    var_problem = setdiff(fixef_vars_all, dataNames)
                    stop("The fixed-effects variable", enumerate_items(var_problem, "s.is"), " not in the data.")
                }

                # Managing fast combine
                if(missnull(combine.quick)){
                    if(NROW(data) > 5e4){
                        combine.quick = TRUE
                    } else {
                        combine.quick = FALSE
                    }
                }

                if(!is.null(fixef_terms_full)){
                    fixef_terms_full = terms_fixef(fixef_fml)
                    if("try-error" %in% class(fixef_terms_full)){
                        stop("Problem extracting the terms of the fixed-effects part of the formula:\n", fixef_terms_full)
                    }
                }


                fixef_terms = fixef_terms_full$fml_terms

                # FEs
                fixef_mat = prepare_df(fixef_terms_full$fe_vars, data, combine.quick)
                if("try-error" %in% class(fixef_mat)){
                    stop("Problem evaluating the fixed-effects part of the formula:\n", fixef_mat)
                }
                fixef_vars = names(fixef_mat)

                # Slopes
                isSlope = any(fixef_terms_full$slope_flag != 0)
                if(isSlope){

                    if(!origin_type %in% c("feols", "feglm")){
                        stop("The use of varying slopes is available only for the functions feols, feglm or fepois.")
                    }

                    slope_mat = prepare_df(fixef_terms_full$slope_vars, data)
                    if("try-error" %in% class(slope_mat)){
                        stop("Problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n", slope_mat)
                    }
                    slope_flag = fixef_terms_full$slope_flag
                    slope_vars = fixef_terms_full$slope_vars
                    slope_vars_list = fixef_terms_full$slope_vars_list

                    # Further controls
                    not_numeric = !sapply(slope_mat, is.numeric)
                    if(any(not_numeric)){
                        stop("In the fixed-effects part of the formula (i.e. in ", as.character(fixef_fml[2]), "), variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope_mat)[not_numeric], "s.is"), " not.")
                    }

                    # slope_flag: 0: no Varying slope // > 0: varying slope AND fixed-effect // < 0: varying slope WITHOUT fixed-effect
                    onlySlope = all(slope_flag < 0)

                }

                # fml update
                fml_char = as.character(fml)
                fml_full = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", paste0(fixef_terms, collapse = "+")))
            }
        }

        #
        # ... NA handling ####
        #

        # We change non-numeric to character (important for parallel qufing)
        is_not_num = sapply(fixef_mat, function(x) !is.numeric(x))
        if(any(is_not_num)){
            for(i in which(is_not_num)){
                if(!is.character(fixef_mat[[i]])){
                    fixef_mat[[i]] = as.character(fixef_mat[[i]])
                }
            }
        }

        if(anyNA(fixef_mat)){
            isNA_fixef = rowSums(is.na(fixef_mat)) > 0

            ANY_NA = TRUE
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_fixef
            msgNA_fixef = paste0(", Fixed-effects: ", numberFormatNormal(sum(isNA_fixef)))

            if(mem.clean){
                rm(isNA_fixef)
                gc2trig = TRUE
            }

        } else {
            msgNA_fixef = ", Fixed-effects: 0"
        }

        # NAs in slopes
        if(isSlope){

            if(mem.clean && gc2trig){
                gc()
                gc2trig = FALSE
            }

            # Convert to double
            who_not_double = which(sapply(slope_mat, is.integer))
            for(i in who_not_double){
                slope_mat[[i]] = as.numeric(slope_mat[[i]])
            }

            info = cpppar_which_na_inf_df(slope_mat, nthreads)
            if(info$any_na_inf){

                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE

                isNA_slope = info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_slope
                msgNA_slope = paste0(", Var. Slopes: ", numberFormatNormal(sum(isNA_slope)))

            } else {
                msgNA_slope = ", Var. Slopes: 0"
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

            nbNA = sum(isNA_sample)

            if(anyNA_sample){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                message_NA = paste0(numberFormatNormal(nbNA), " observation", plural(nbNA), " removed because of ", msg, " (", msgNA_y, msgNA_L, msgNA_NL, msgNA_fixef, msgNA_slope, msgNA_offset, msgNA_weight, msgNA_split, ").")
                notes = c(notes, message_NA)
            }

            if(nbNA == nobs){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                stop("All observations contain ", msg, ". Estimation cannot be done. (Breakup: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_fixef, msgNA_slope, msgNA_offset, msgNA_weight, msgNA_split, ".)")
            }

            if(any0W){
                # 0 weight => like NAs
                isNA_sample = isNA_sample | is0W
                nbNA = sum(isNA_sample)

                if(nbNA == nobs){
                    if(anyNA_sample){
                        msg = msg_na_inf(ANY_NA, ANY_INF)
                        stop("All observations contain ", msg, " or are 0-weight. Estimation cannot be done. (0-weight: ", sum(is0W), ", breakup ", msg, ": ", msgNA_y, msgNA_L, msgNA_NL, msgNA_fixef, msgNA_slope, msgNA_offset, msgNA_weight, msgNA_split, ")")
                    } else {
                        stop("All observations are 0-weight. Estimation cannot be done.")
                    }
                }
            }

            # we drop the NAs from the fixef matrix
            fixef_mat = fixef_mat[!isNA_sample, , drop = FALSE]
            obs2remove_NA = which(isNA_sample)
            index_noNA = (1:nobs)[!isNA_sample]

            if(isSlope){
                slope_mat = slope_mat[!isNA_sample, , drop = FALSE]
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

        info_fe = setup_fixef(fixef_mat = fixef_mat, lhs = lhs, fixef_vars = fixef_vars, fixef.rm = fixef.rm, family = family, isSplit = isSplit, split.full = split.full, origin_type = origin_type, isSlope = isSlope, slope_flag = slope_flag, slope_mat = slope_mat, slope_vars_list = slope_vars_list, nthreads = nthreads)

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

        # Q = length(fixef_vars) # terms: contains FEs + slopes
        #
        # fixef_removed = list()
        # obs2remove = c()
        #
        # rm_0 = !family == "gaussian"
        # rm_1 = family == "logit"
        # rm_single = fixef.rm %in% c("singleton", "both")
        # do_sum_y = !origin_type %in% c("feols", "feglm")
        #
        # if(isSplit){
        #     # we delay the removal of FEs
        #     rm_0 = rm_1 = FALSE
        # }
        #
        # only_slope = FALSE
        # if(isSlope){
        #     # only slope: not any non-slope
        #     only_slope = slope_flag < 0
        #
        #     # shallow copy
        #     slope_variables = as.list(slope_mat)
        # }
        #
        # if(mem.clean && gc2trig){
        #     gc()
        #     gc2trig = FALSE
        # }
        #
        # quf_info_all = cpppar_quf_table_sum(x = fixef_mat, y = lhs, do_sum_y = do_sum_y, rm_0 = rm_0, rm_1 = rm_1, rm_single = rm_single, only_slope = only_slope, nthreads = nthreads)
        #
        # fixef_id = quf_info_all$quf
        # # names
        # fixef_names = list()
        # is_string = sapply(fixef_mat, is.character)
        # for(i in 1:length(fixef_id)){
        #     if(is_string[i]){
        #         fixef_names[[i]] = fixef_mat[[i]][quf_info_all$items[[i]]]
        #     } else {
        #         fixef_names[[i]] = quf_info_all$items[[i]]
        #     }
        # }
        # # table/sum_y/sizes
        # fixef_table = quf_info_all$table
        # sum_y_all = quf_info_all$sum_y
        # fixef_sizes = lengths(fixef_table)
        #
        # # If observations have been removed:
        # if(!is.null(quf_info_all$obs_removed)){
        #
        #     # which obs are removed
        #     obs2remove = which(quf_info_all$obs_removed)
        #
        #     # update of the lhs
        #     lhs = lhs[-obs2remove]
        #
        #     # update of the slope variables
        #     if(isSlope){
        #         for(i in seq_along(slope_variables)) slope_variables[[i]] = slope_variables[[i]][-obs2remove]
        #     }
        #
        #     # Names of the FE removed
        #     for(i in 1:length(fixef_id)){
        #         if(is_string[i]){
        #             fixef_removed[[i]] = fixef_mat[[i]][quf_info_all$fe_removed[[i]]]
        #         } else {
        #             fixef_removed[[i]] = quf_info_all$fe_removed[[i]]
        #         }
        #     }
        #
        #     names(fixef_removed) = fixef_vars
        #
        #     # Then the "Notes"
        #     nb_missing = lengths(fixef_removed)
        #     message_fixef = paste0(paste0(nb_missing, collapse = "/"), " fixed-effect", plural(sum(nb_missing)), " (", numberFormatNormal(length(obs2remove)), " observation", plural_len(obs2remove), ") removed because of only ", ifelse(family=="logit", "zero (or only one)", "zero"), " outcomes.")
        #     notes = c(notes, message_fixef)
        # }

        # #
        # # we save the fixed-effects IDs + sizes (to be returned in "res") [original order!]
        # #
        #
        # fixef_id_res = list()
        # for(i in 1:Q){
        #     dum = fixef_id[[i]]
        #     attr(dum, "fixef_names") = as.character(fixef_names[[i]])
        #     fixef_id_res[[fixef_vars[i]]] = dum
        # }
        #
        #
        # # The real size is equal to nb_coef * nb_slopes
        # if(isSlope){
        #     fixef_sizes_real = fixef_sizes * (1 + abs(slope_flag) - (slope_flag < 0))
        # } else {
        #     fixef_sizes_real = fixef_sizes
        # }
        #
        # fixef_sizes_res = fixef_sizes
        # names(fixef_sizes_res) = fixef_vars
        #
        # #
        # # We re-order the fixed-effects
        # #
        #
        # if(any(fixef_sizes_real != sort(fixef_sizes_real, decreasing = TRUE))){
        #     # FE with the most cases first (limits precision problems)
        #
        #     new_order = order(fixef_sizes_real, decreasing = TRUE)
        #
        #     fixef_sizes = fixef_sizes[new_order]
        #     fixef_id = fixef_id[new_order]
        #     sum_y_all = sum_y_all[new_order]
        #     fixef_table = fixef_table[new_order]
        #
        #     if(isSlope){
        #         slope_variables = slope_variables[unlist(slope_vars_list[new_order], use.names = FALSE)]
        #         slope_flag = slope_flag[new_order]
        #     }
        #
        # } else {
        #     new_order = 1:Q
        # }


        notes = c(notes, message_fixef)

        if(length(obs2remove_NA) > 0){
            # we update the value of obs2remove (will contain both NA and removed bc of outcomes)
            if(length(obs2remove) > 0){
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

            if(anyNA_sample){
                message_NA = paste0(numberFormatNormal(nbNA), " observation", plural(nbNA), " removed because of NA values (", msgNA_y, msgNA_L, msgNA_NL, msgNA_offset, msgNA_weight, msgNA_split, ").")
                notes = c(notes, message_NA)

                if(nbNA == nobs){
                    stop("All observations contain NAs. Estimation cannot be done. (Breakup: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_offset, msgNA_weight, msgNA_split, ".)")
                }
            }

            if(any0W){
                # 0 weight => like NAs
                isNA_sample = isNA_sample | is0W
                nbNA = sum(isNA_sample)

                if(nbNA == nobs){
                    if(anyNA_sample){
                        stop("All observations are either NA or 0-weight. Estimation cannot be done. (0-weight: ", sum(is0W), ", breakup NA: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_offset, msgNA_weight, msgNA_split, ".)")
                    } else {
                        stop("All observations are 0-weight. Estimation cannot be done.")
                    }
                }
            }

            # note = ifelse((anyNA_sample + any0W) == 2, "NOTES: ", "NOTE: ")
            # if(notes) message(note, message_NA, ifelse(anyNA_sample && any0W, "\n       ", ""), message_0W)

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

    # NA & problem management
    if(length(obs2remove) > 0){
        # we kick out the problems (both NA related and fixef related)

        # We drop only 0 variables (may happen for factors)
        linear.mat = linear.mat[-obs2remove, , drop = FALSE]

        if(useModel.matrix){
            # There are factors => possibly some vars are only 0 now that NAs are removed

            only_0 = cpppar_check_only_0(linear.mat, nthreads)
            if(all(only_0 == 1)){
                stop("After removing NAs, not a single explanatory variable is different from 0.")
            } else if(any(only_0 == 1)){
                linear.mat = linear.mat[, only_0 == 0, drop = FALSE]

                # useful when feNmlm
                linear.params <- colnames(linear.mat)
                params <- c(nonlinear.params, linear.params)
                lparams <- length(params)
                varnames <- c(nonlinear.varnames, linear.varnames)
            }
        }

        if(Q == 0){
            # if Q > 0: done already when managing the fixed-effects
            if(is.list(lhs)){
                for(i in seq_along(lhs)){
                    lhs[[i]] = lhs[[i]][-obs2remove]
                }
            } else {
                lhs = lhs[-obs2remove]
            }

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
                params <- nonlinear.params
                lparams <- length(params)
                varnames <- nonlinear.varnames
            } else{
                linear.mat = linear.mat[, -var2remove, drop = FALSE]
                linear.params <- colnames(linear.mat)
                params <- c(nonlinear.params, linear.params)
                lparams <- length(params)
                varnames <- c(nonlinear.varnames, linear.varnames)
            }
        }
    }


    #
    # Other setups ####
    #

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

    if(lparams == 0 && Q == 0) stop("No parameter to be estimated.")

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

        lower[params[!params %in% names(lower)]] <- -Inf
        lower <- unlist(lower[params])

    }	else {
        lower <- rep(-Inf, lparams)
        names(lower) <- params
    }

    if(!missnull(upper)){

        if(typeof(upper) != "list"){
            stop("'upper' MUST be a list.")
        }

        upper[params[!params %in% names(upper)]] <- Inf
        upper <- unlist(upper[params])

    }	else {
        upper <- rep(Inf, lparams)
        names(upper) <- params
    }

    lower <- c(lower)
    upper <- c(upper)

    #
    # Controls: user defined gradient
    #

    #
    # PRECISION + controls ####
    #

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
    start[linear.params] <- linear.start
    params <- names(start)
    start <- unlist(start)
    start <- c(start)
    lparams <- length(params)
    names(start) <- params

    # The right order of upper and lower
    upper = upper[params]
    lower = lower[params]

    ####
    #### Sending to the env ####
    ####

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
    assign("fml", fml, env)
    assign("origin", origin, env)
    assign("origin_type", origin_type, env)
    assign("warn", warn, env)
    assign("mem.clean", mem.clean, env)
    assign("nthreads", nthreads, env)

    # Multi
    assign("do_multi_rhs", multi_rhs, env)
    if(multi_rhs){
        assign("fml_rhs_all", fml_rhs_all, env)
        assign("fml_all_sw", fml_all_sw, env)
        assign("is_cumul", is_cumul, env)
        assign("fake_intercept", fake_intercept, env)

        if(length(sw_all_vars) > 0){
            if(length(obs2remove) > 0){
                assign("data", data[-obs2remove, sw_all_vars, drop = FALSE], env)
            } else {
                assign("data", data[, sw_all_vars, drop = FALSE], env)
            }
        }

    }

    assign("do_multi_lhs", multi_lhs, env)

    #
    # The MODEL0 => to get the init of the theta for the negbin
    #

    check_arg(theta.init, "numeric scalar gt{0}")
    if(missing(theta.init)){
        theta.init = NULL
    }

    if(origin_type == "feNmlm" || computeModel0){
        model0 = get_model_null(env, theta.init)
        theta.init = model0$theta
    } else {
        model0 = NULL
    }

    # For the negative binomial:
    if(family == "negbin"){
        params = c(params, ".theta")
        start = c(start, theta.init)
        names(start) = params
        upper = c(upper, 10000)
        lower = c(lower, 1e-3)
    }

    assign("model0", model0, env)
    onlyFixef = !isLinear && !isNonLinear && Q > 0
    assign("onlyFixef", onlyFixef, env)
    assign("start", start, env)
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

        if(isSlope == FALSE){
            slope_flag = rep(0L, length(fixef_vars))
            slope_variables = list(0)
            slope_vars_list = list(0)
        }

        assign("new_order_original", new_order, env)
        assign("fixef_names", fixef_names, env)
        assign("fixef_vars", fixef_vars, env)

        assign_fixef_env(env, nobs, family, origin_type, fixef_id, fixef_sizes, fixef_table, sum_y_all, slope_flag, slope_variables, slope_vars_list)
    }


    #
    # Specific to femlm/feNmlm
    #

    if(origin_type == "feNmlm"){

        # basic NL
        envNL = new.env()
        assign("isNL", isNL, env)
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
            mu = NULL
            try(mu <- eval(nl.call, envir = envNL), silent = FALSE)
            if(is.null(mu)){
                # the non linear part could not be evaluated - ad hoc message
                stop("The non-linear part (NL.fml) could not be evaluated. There may be a problem in 'NL.fml'.")
            }

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
        # if(isLinear){
        #     start_coef_linear = unlist(start[linear.params])
        #     mu <- mu + cpppar_xbeta(linear.mat, start_coef_linear, nthreads)
        # }

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
        assign("split", split, env)
        assign("split.items", split.items, env)
        assign("split.full", split.full, env)
    }

    # fixest tag
    assign("fixest_env", TRUE, env)

    #
    # Res ####
    #

    #
    # Preparation of the results list (avoid code repetition in estimation funs)
    #

    res = list(nobs=nobs, nobs_origin=nobs_origin, fml=fml, call = mc_origin, method = origin)

    if(isFixef) res$fml_full = fml_full

    if(!is.null(fml_no_xpd)) res$fml_no_xpd = fml_no_xpd

    if(isFit) res$fromFit = TRUE

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
            original_order = order(new_order)
            res$slope_flag = slope_flag[order(new_order)]
            res$slope_flag_reordered = slope_flag
            res$slope_variables_reordered = slope_variables
            res$fe.reorder = new_order
        }

        res$fixef_id = fixef_id_res
        res$fixef_sizes = fixef_sizes_res

    }

    # Observations removed (either NA or fixed-effects)
    if(isSubset){
        # subset does not accept duplicate values eg c(1, 1, 1, 2)
        all_obs = 1:nobs_before_subset
        obs2keep = all_obs[subset]
        obs2remove_subset = all_obs[-obs2keep]
        obs2remove = sort(c(obs2remove_subset, obs2keep[obs2remove]))
    }

    if(length(obs2remove) > 0){
        res$obsRemoved = obs2remove
        if(isFixef && any(lengths(fixef_removed) > 0)){
            res$fixef_removed = fixef_removed
        }
    }

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
    if(!is.null(interaction.info)){
        res$interaction.info = interaction.info
    }

    assign("res", res, env)

    env
}


setup_fixef = function(fixef_mat, lhs, fixef_vars, fixef.rm, family, isSplit, split.full = FALSE, origin_type, isSlope, slope_flag, slope_mat, slope_vars_list, fixef_names_old = NULL, fixef_sizes = NULL, obs2keep = NULL, nthreads){

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

        # shallow copy
        slope_variables = as.list(slope_mat)

        if(!is.null(obs2keep)){
            for(i in seq_along(slope_variables)){
                slope_variables[[i]] = slope_variables[[i]][obs2keep]
            }
        }
    } else {
        slope_flag = rep(0L, length(fixef_mat))
        slope_variables = list(0)
    }

    if(isSplit && split.full == FALSE){
        # We don't do anything => it will be taken care of in the splits
        res = list(Q = Q, fixef_id = fixef_mat, fixef_names = "SPLIT_NO_FULL", sum_y_all = 0, fixef_sizes = 0, fixef_table = 0, obs2remove = NULL, fixef_removed = NULL, message_fixef = NULL, lhs = lhs, slope_variables = slope_variables, slope_flag = slope_flag, new_order = 1:Q)

        return(res)
    }

    do_keep = !is.null(obs2keep)

    if(isSplitNoFull && do_keep){
        # Here fixef_mat is a DF
        fixef_mat = fixef_mat[obs2keep, , drop = FALSE]
    }

    if(is.null(obs2keep)){
        obs2keep = 0
    }

    if(is.null(fixef_sizes)){
        fixef_sizes = 0
    }

    quf_info_all = cpppar_quf_table_sum(x = fixef_mat, y = lhs, do_sum_y = do_sum_y, rm_0 = rm_0, rm_1 = rm_1, rm_single = rm_single, only_slope = only_slope, nthreads = nthreads, do_refactor = isRefactor, r_x_sizes = fixef_sizes, obs2keep = obs2keep)


    fixef_id = quf_info_all$quf
    # names

    fixef_names = list()
    if(isRefactor == FALSE){
        is_string = sapply(fixef_mat, is.character)
        for(i in 1:length(fixef_id)){
            if(is_string[i]){
                fixef_names[[i]] = fixef_mat[[i]][quf_info_all$items[[i]]]
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

        # update of the lhs (only if not multi LHS)
        if(multi_lhs == FALSE) lhs = lhs[-obs2remove]

        # update of the slope variables
        if(isSlope){
            for(i in seq_along(slope_variables)) slope_variables[[i]] = slope_variables[[i]][-obs2remove]
        }

        # Names of the FE removed
        for(i in 1:length(fixef_id)){
            if(is_string[i]){
                fixef_removed[[i]] = fixef_mat[[i]][quf_info_all$fe_removed[[i]]]
            } else {
                fixef_removed[[i]] = quf_info_all$fe_removed[[i]]
            }
        }

        names(fixef_removed) = fixef_vars

        # Then the "Notes"
        nb_missing = lengths(fixef_removed)
        if(rm_0 == FALSE){
            n_single = sum(nb_missing)
            message_fixef = paste0(n_single, " fixed-effect singleton", plural(n_single, "s.was"), " removed (", numberFormatNormal(length(obs2remove)), " observation", plural_len(obs2remove), ifelse(Q == 1, "", paste0(", breakup: ", paste0(nb_missing, collapse = "/"))), ").")
        } else {
            message_fixef = paste0(paste0(nb_missing, collapse = "/"), " fixed-effect", plural(sum(nb_missing)), " (", numberFormatNormal(length(obs2remove)), " observation", plural_len(obs2remove), ") removed because of only ", ifelse(rm_1, "0 (or only 1)", "0"), " outcomes", ifelse(rm_single && !rm_1, " or singletons", ""), ".")
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
    res = list(Q = Q, fixef_id = fixef_id, fixef_names = fixef_names, sum_y_all = sum_y_all, fixef_sizes = fixef_sizes, fixef_table = fixef_table, obs2remove = obs2remove, fixef_removed = fixef_removed, message_fixef = message_fixef, lhs = lhs)

    res$slope_variables = slope_variables
    res$slope_flag = slope_flag

    if(multi_lhs == FALSE || family == "gaussian"){
        res$fixef_id_res = fixef_id_res
        res$fixef_sizes_res = fixef_sizes_res
        res$new_order = new_order
    }

    return(res)
}





assign_fixef_env = function(env, nobs, family, origin_type, fixef_id, fixef_sizes, fixef_table, sum_y_all, slope_flag, slope_variables, slope_vars_list){

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



reshape_env = function(env, obs2keep = NULL, lhs = NULL, rhs = NULL, assign_lhs = TRUE, assign_rhs = TRUE, fml = NULL){
    # env: environment from an estimation
    # This functions reshapes the environment to perform the neww estimation
    # either by selecting some observation (in split)
    # either by changing the depvar (leading to new NAs, etc)
    # obs2keep => must be a which()

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

    if(isFixef && length(obs2keep) > 0){

        # gt("nothing")

        if(assign_lhs){
            if(is.list(lhs)){
                # list means multiple LHS
                for(i in seq_along(lhs)){
                    lhs[[i]] = lhs[[i]][obs2keep]
                }
            } else {
                lhs = lhs[obs2keep]
            }
        }

        # gt("fixef, dropping lhs")

        fixef_mat       = get("fixef_id_list", env)
        fixef_vars      = get("fixef_vars", env)
        fixef.rm        = get("fixef.rm", env)
        fixef_names_old = get("fixef_names", env)
        fixef_sizes     = get("fixef_sizes", env)

        slope_flag      = get("slope_flag", env)
        slope_mat       = get("slope_variables", env)
        slope_vars_list = get("slope_vars_list", env)
        isSlope         = any(slope_flag != 0)

        # gt("fixef, dropping")

        # browser()
        # We refactor the fixed effects => we may remove even more obs
        info_fe = setup_fixef(fixef_mat = fixef_mat, lhs = lhs, fixef_vars = fixef_vars, fixef.rm = fixef.rm, family = family, isSplit = FALSE, origin_type = origin_type, isSlope = isSlope, slope_flag = slope_flag, slope_mat = slope_mat, slope_vars_list = slope_vars_list, fixef_names_old = fixef_names_old, fixef_sizes = fixef_sizes, obs2keep = obs2keep, nthreads = nthreads)

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

            if(assign_lhs && is.list(lhs)){
                if(is.list(lhs)){
                    # list means multiple LHS
                    for(i in seq_along(lhs)){
                        lhs[[i]] = lhs[[i]][-obs2remove]
                    }
                }
            }

            obs2keep = obs2keep[-obs2remove]
        }

        lhs_done = TRUE

        assign_fixef_env(new_env, nobs, family, origin_type, fixef_id, fixef_sizes, fixef_table, sum_y_all, slope_flag, slope_variables, slope_vars_list)

    } else if(isFixef){
        slope_flag = get("slope_flag", env)
        isSlope = any(slope_flag != 0)
    }

    # gt("fixef")

    # This is very tedious
    if(!is.null(obs2keep)){
        # Recreating the values

        #
        # The left hand side
        #

        if(assign_lhs){
            if(lhs_done == FALSE){
                if(is.list(lhs)){
                    # list means multiple LHS
                    for(i in seq_along(lhs)){
                        lhs[[i]] = lhs[[i]][obs2keep]
                    }
                } else {
                    lhs = lhs[obs2keep]
                }
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
                rhs = rhs[obs2keep, , drop = FALSE]

                only_0 = cpppar_check_only_0(rhs, nthreads)
                if(all(only_0 == 1)){
                    stop("After removing NAs, not a single explanatory variable is different from 0.")

                } else if(any(only_0 == 1)){
                    rhs = rhs[, only_0 == 0, drop = FALSE]

                    linear.params = colnames(rhs)
                    nonlinear.params = get("nonlinear.params", env)
                    params = c(nonlinear.params, linear.params)
                    assign("linear.params", linear.params, new_env)
                    assign("params", params, new_env)

                    # useful when feNmlm or feglm
                    if(origin_type %in% c("feNmlm", "feglm")){

                        start = get("start", env)
                        qui = names(start) %in% params
                        assign("start", start[qui], new_env)

                        if(origin_type == "feNmlm"){
                            upper = get("upper", env)
                            lower = get("lower", env)

                            qui = names(upper) %in% params
                            assign("upper", upper[qui], new_env)
                            assign("lower", lower[qui], new_env)
                        }
                    }

                }

                assign("linear.mat", rhs, new_env)
            }
            assign("isLinear", isLinear, new_env)
        }

        #
        # New data if stepwise estimation
        #

        if(!is.null(env$data) && isTRUE(env$do_multi_rhs)){
            data = get("data", env)
            assign("data", data[obs2keep, , drop = FALSE], new_env)
        }

        #
        # The non-linear part, the weight and the offset
        #

        isOffset = length(env$offset.value) > 1
        if(isOffset){
            offset.value = get("offset.value", env)
            assign("offset.value", offset.value[obs2keep], new_env)
        }

        isWeight = length(env$weights.value) > 1
        if(isWeight){
            weights.value = get("weights.value", env)
            assign("weights.value", weights.value[obs2keep], new_env)
        }

        if(isTRUE(env$isNL)){
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

            if(!is.list(lhs) && family$family_equiv == "poisson"){

                family_funs = get("family_funs", env)

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

        # gt("res, values")

        #
        # Reformatting "res"
        #

        res = get("res", env)

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

        # gt("res, fixef")

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
        if(is.null(res$obsKept)){
            res$obsKept = obs2keep
        } else {
            res$obsKept_bis = obs2keep
        }

        assign("res", res, new_env)

        # gt("res, obsRemoved")

    } else {

        if(save_lhs){
            # Here lhs is ALWAYS a vector
            assign("lhs", lhs, new_env)

            if(origin_type == "feglm" && family$family_equiv == "poisson"){

                family_funs = get("family_funs", env)

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
                    qui = names(start) %in% params
                    assign("start", start[qui], new_env)

                    if(origin_type == "feNmlm"){
                        upper = get("upper", env)
                        lower = get("lower", env)

                        qui = names(upper) %in% params
                        assign("upper", upper[qui], new_env)
                        assign("lower", lower[qui], new_env)
                    }
                }
            }

            assign("isLinear", isLinear, new_env)

            # Finally the DoF
            res = get("res", env)
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

            if(!isLinear){
                res$onlyFixef = TRUE
            }

            assign("res", res, new_env)
        }

    }

    # We save the formula
    if(!is.null(fml)){

        res = get("res", new_env)

        if(isFixef){
            fe = formula(Formula(res$fml_full), lhs = 0, rhs = 2)
            res$fml_full = as.formula(paste0(deparse_long(fml), "|", deparse_long(fe[[2]])))
        }
        res$fml = fml

        assign("res", res, new_env)
    }


    assign("fixest_env", TRUE, new_env)

    return(new_env)
}


cstepwise = csw = stepwise = sw = function(...){
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

cstepwise0 = csw0 = stepwise0 = sw0 = function(...){
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

# cstepwise = csw = function(...){
#     mc = match.call(expand.dots = TRUE)
#
#     n = length(mc) - 1
#
#     if(n == 0){
#         return("")
#     }
#
#     res = vector("character", n)
#     res[[1]] = deparse_long(mc[[2]])
#
#     res_add = vector("character", n - 1)
#
#     if(n >= 2){
#         for(i in 2:n){
#             value = deparse_long(mc[[i + 1]])
#             res_add[[i - 1]] = value
#             res[[i]] = paste0(res[[i - 1]], " + ", value)
#         }
#     }
#
#     attr(res, "additional_vars") = res_add
#
#     res
# }
#
# cstepwise0 = csw0 = function(...){
#     mc = match.call(expand.dots = TRUE)
#
#     n = length(mc) - 1
#
#     if(n == 0){
#         return("")
#     }
#
#     res = vector("character", n + 1)
#     res[[2]] = deparse_long(mc[[2]])
#
#     res_add = vector("character", n)
#
#     for(i in 2:n){
#         value = deparse_long(mc[[i + 1]])
#         res_add[[i - 1]] = value
#         res[[i + 1]] = paste0(res[[i]], " + ", value)
#     }
#
#     attr(res, "additional_vars") = res_add
#
#     res
# }


















