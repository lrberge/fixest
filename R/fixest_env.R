#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Jun 10 23:46:47 2019
# Purpose: sets up the params + performs
#           all necessay checks
#----------------------------------------------#


fixest_env <- function(fml, data, family=c("poisson", "negbin", "logit", "gaussian"), NL.fml,
                       fixef, na_inf.rm = getFixest_na_inf.rm(), NL.start, lower, upper, NL.start.init,
                       offset, linear.start = 0, jacobian.method = "simple",
                       useHessian = TRUE, hessian.args = NULL, opt.control = list(),
                       y, X, fixef_mat, panel.id,
                       nthreads = getFixest_nthreads(),
                       verbose = 0, theta.init, fixef.tol = 1e-5, fixef.iter = 10000,
                       deriv.iter = 5000, deriv.tol = 1e-4, glm.iter = 25, glm.tol = 1e-8,
                       etastart, mustart,
                       warn = TRUE, notes = getFixest_notes(), combine.quick,
                       origin_bis, origin = "feNmlm", mc_origin, mc_origin_bis, mc_origin_ter,
                       computeModel0 = FALSE, weights,
                       from_update = FALSE, object, sumFE_init, debug = FALSE, ...){

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
    main_args = c("fml", "data", "panel.id", "offset", "na_inf.rm", "fixef.tol", "fixef.iter", "fixef", "nthreads", "verbose", "warn", "notes", "combine.quick", "start")
    femlm_args = c("family", "theta.init", "linear.start", "opt.control", "deriv.tol", "deriv.iter")
    feNmlm_args = c("NL.fml", "NL.start", "lower", "upper", "NL.start.init", "jacobian.method", "useHessian", "hessian.args")
    feglm_args = c("family", "weights", "glm.iter", "glm.tol", "etastart", "mustart")
    feols_args = c("weights")
    internal_args = c("debug", "object", "from_update", "sumFE_init")

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
        if(warn) warning(substr(deparse(mc_origin)[1], 1, 15), "...: ", enumerate_items(args_invalid, "is"), " not ", ifsingle(args_invalid, "a valid argument", "valid arguments"), " for function ", origin, ".", call. = FALSE, immediate. = TRUE)
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

        fun_dev = family$dev.resids
        family$dev.resids = function(y, mu, eta, wt) sum(fun_dev(y, mu, wt))
        fun_mu.eta = family$mu.eta
        family$mu.eta = function(mu, eta) fun_mu.eta(eta)

        if(is.null(family$valideta)) family$valideta = function(...) TRUE
        if(is.null(family$validmu)) family$validmu = function(...) TRUE

        # QUIRKINESS => now family becomes the family name and the functions become family_funs
        family_funs = family
        family = family_type

        computeModel0 = family_equiv %in% c("poisson", "logit")
    }

    if(!isLogical(notes)){
        stop("Argument 'notes' must be a single logical.")
    }

    show_notes = notes
    notes = c()

    if(!isLogical(warn)){
        stop("Argument 'warn' must be a single logical.")
    }

    if(!isLogical(na_inf.rm)){
        stop("Argument 'na_inf.rm' must be a single logical.")
    }

    # flags for NA infinite vales => will be used in the message (good to discrimintae b/w NA and inf)
    ANY_INF = FALSE
    ANY_NA = FALSE
    anyNA_sample = FALSE
    isNA_sample = FALSE
    message_NA = ""

    #
    # nthreads argument
    if(!isScalar(nthreads) || (nthreads %% 1) != 0 || nthreads <= 0){
        stop("The argument 'nthreads' must be an integer greater or equal to 1 and lower than the number of threads available (", max(get_nb_threads(), 1), ").")
    }

    if(nthreads > 1){
        max_threads = get_nb_threads()
        if(max_threads == 0){
            warning("OpenMP not detected: cannot use ", nthreads, " threads, single-threaded mode instead.")
            nthreads = 1
        } else if(nthreads > max_threads){
            warning("Asked for ", nthreads, " threads while the maximum is ", max_threads, ". Set to ", max_threads, " threads instead.")
            nthreads = max_threads
        }
    }

    # The family functions (for femlm only)
    famFuns = switch(family,
                     poisson = ml_poisson(),
                     negbin = ml_negbin(),
                     logit = ml_logit(),
                     gaussian = ml_gaussian())


    #
    # Formatting data ####
    #

    isPanel = FALSE
    if(isFit){
        isFixef = !missnull(fixef_mat)
    } else {

        #
        # The data
        if(missing(data)) stop("You must provide the argument 'data' (currently it is missing).")
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
        if(!"formula" %in% class(fml)) stop("The argument 'fml' must be a formula.")
        fml = formula(fml) # we regularize the formula to check it
        if(length(fml) != 3) stop("The formula must be two sided: e.g. y~x1+x2, or y~x1+x2|fe1+fe2.")

        #
        # Panel setup
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
                stop("Problem in the formula: ", gsub("_expand", "", fml))
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

            if(!isVector(fixef) || !is.character(fixef)){
                stop("Argument 'fixef', when provided, must be a character vector of variable names.")
            }

            if(any(!fixef %in% dataNames)){
                var_problem = fixef[!fixef %in% dataNames]
                stop("The argument 'fixef' must be variable names! Variable", enumerate_items(var_problem, "s.is"), " not in the data.")

            }

            fixef_vars = fixef
        }

        # fml_full now does not contain the clusters, but will later contain them
        fml_full = fml
    }


    #
    # The left hand side
    #

    # evaluation
    if(isFit){

        if(missing(y)){
            stop("You must provide argument 'y' when using ", origin_type, ".fit.")
        }

        if(!is.numeric(y)){
            stop("Argument 'y' must be numeric.")
        }

        if(!is.null(dim(y)) || is.list(y)){
            stop("Argument y must be a vector.")
        }

        # for cpp (as double):
        lhs = as.numeric(as.vector(y))

        # we reconstruct a formula
        fml = as.formula(paste0(deparse_long(mc_origin[["y"]]), "~1"))

    } else {

        # The LHS must contain only values in the DF
        namesLHS = all.vars(fml[[2]])
        if(length(namesLHS) == 0){
            stop("The right hand side of the formula (", deparse_long(fml[[2]]), ") contains no variable!")
        } else if(!all(namesLHS %in% dataNames)){
            not_there = namesLHS[!namesLHS %in% dataNames]
            stop(ifsingle(not_there, "The v", "V"), "ariable", enumerate_items(not_there, "s.is"), " in the LHS of the formula but not in the dataset.")
        }

        lhs = try(eval(fml[[2]], data), silent = TRUE)

        if("try-error" %in% class(lhs)){
            stop("Evaluation of the left-hand-side (equal to ", deparse_long(fml[[2]]), ") raises an error: \n", lhs)
        } else if(is.logical(lhs) || is.integer(lhs)){
            lhs = as.numeric(as.vector(lhs))
        } else if(!is.numeric(lhs)){
            cls = class(lhs)
            stop("The left hand side (", deparse_long(fml[[2]]), ") is not numeric. The class", ifsingle(cls, "", "e"), enumerate_items(cls, "s.is"), " not supported.")
        } else {
            lhs = as.numeric(as.vector(lhs)) # complex
        }
    }

    lhs_clean = lhs # copy used for NA case
    nobs = length(lhs)

    anyNA_y = FALSE
    msgNA_y = "LHS: 0"
    info = cpppar_which_na_inf_vec(lhs, nthreads)
    if(info$any_na_inf){
        if(!na_inf.rm){
            msg = msg_na_inf(info$any_na, info$any_inf)
            obs = head(which(info$is_na_inf), 3)
            stop("The left hand side of the fomula contains ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide data without ", msg, " (or use na_inf.rm).")
        } else {
            # If na_inf.rm => we keep track of the NAs
            if(info$any_na) ANY_NA = TRUE
            if(info$any_inf) ANY_INF = TRUE

            anyNA_y = TRUE
            isNA_y = info$is_na_inf
            anyNA_sample = TRUE
            isNA_sample = isNA_sample | isNA_y
            msgNA_y = paste0("LHS: ", numberFormatNormal(sum(isNA_y)))
        }

        lhs_clean = lhs[!isNA_y]
    }

    # we check the var is not a constant
    if(cpp_isConstant(lhs_clean)){
        stop("The dependent variable is a constant. Estimation cannot be done.")
    }

    if(family %in% c("poisson", "negbin") && any(lhs_clean < 0)){
        stop("Negative values of the dependant variable are not allowed for the \"", family, "\" family.")
    }

    if(origin_type == "feNmlm" && family %in% "logit" && !all(lhs_clean==0 | lhs_clean==1)){
        stop("The dependent variable has values different from 0 or 1.\nThis is not allowed with the \"logit\" family.")
    }

    #
    # Controls and setting of the linear part:
    #

    interaction.info = NULL
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
        }

    } else {
        isLinear = FALSE
        options("fixest_interaction_ref" = NULL)

        if(grepl("[^:]::[^:]", deparse_long(fml[[3]]))){

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
            linear.fml = fml
        }

        if(isLinear){

            if(!all(linear.varnames %in% dataNames)){
                not_there = setdiff(linear.varnames, dataNames)
                stop(ifsingle(not_there, "The v", "V"), "ariable", enumerate_items(not_there, "s.is"), " in the RHS of the formula but not in the dataset.")
            }

            if(isFixef){
                # if dummies are provided, we make sure there is an
                # intercept so that factors can be handled properly
                linear.fml = update(linear.fml, ~.+1)
            }

            #
            # We construct the linear matrix
            #

            linear.mat = try(fixest_model_matrix(fml, data), silent = TRUE)
            if("try-error" %in% class(linear.mat)){
                stop("Evaluation of the right-hand-side of the formula raises an error: ", linear.mat)
            }

            useModel.matrix = attr(linear.mat, "useModel.matrix")

            # Interaction information => if no interaction: NULL
            interaction.info = getOption("fixest_interaction_ref")


            # # we look at whether there are factor-like variables to be evaluated
            # # if there is factors => model.matrix
            # types = sapply(data[, dataNames %in% linear.varnames, FALSE], class)
            # if(length(types) == 0 || grepl("factor", deparse_long(linear.fml)) || any(types %in% c("character", "factor"))){
            #     useModel.matrix = TRUE
            # } else {
            #     useModel.matrix = FALSE
            # }
            #
            # if(useModel.matrix){
            #     # linear.mat = stats::model.matrix(linear.fml, data)
            #     # to catch the NAs, model.frame needs to be used....
            #     linear.mat = try(stats::model.matrix(linear.fml, stats::model.frame(linear.fml, data, na.action=na.pass)), silent = TRUE)
            #     if("try-error" %in% class(linear.mat)){
            #         stop("Evaluation of the right-hand-side of the formula raises an error: \n", linear.mat)
            #     }
            # } else {
            #     linear.mat = try(prepare_matrix(linear.fml, data), silent = TRUE)
            #     if("try-error" %in% class(linear.mat)){
            #         stop("Evaluation of the right-hand-side of the formula raises an error: \n", linear.mat)
            #     }
            # }

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
            if(!na_inf.rm){
                # Default behavior: no NA tolerance
                quiNA = colSums(linear.mat[info$is_na_inf, , drop=FALSE]) > 0
                whoIsNA = linear.params[quiNA]
                msg = msg_na_inf(info$any_na, info$any_inf)
                obs = which(info$is_na_inf)

                stop("The right hand side contains ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide data without ", msg, " (or use na_inf.rm). FYI the variable", enumerate_items(whoIsNA, "s.is"), " concerned.")
            } else {
                # If na_inf.rm => we keep track of the NAs
                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE

                isNA_L = info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_L
                msgNA_L = paste0(", RHS: ", numberFormatNormal(sum(isNA_L)))
            }

        } else {
            msgNA_L = ", RHS: 0"
        }

    } 	else {
        linear.params <- linear.start <- linear.varnames <- NULL
        useModel.matrix = FALSE
    }


    #
    # The nonlinear part:
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
        qui_num = sapply(data_NL, is.numeric)
        if(any(!qui_num)){
            var_pblm = nonlinear.varnames[qui_num]
            stop("In NL.fml, the variable", enumerate_items(var_pblm, "s.is"), " not numeric. This is not allowed.")
        }
        data_NL = as.numeric(as.matrix(data_NL))

        # Control for NAs
        anyNA_NL = FALSE
        info = which_na_inf(data_NL, nthreads)
        if(info$any_na_inf){
            anyNA_NL = TRUE
            if(!na_inf.rm){
                quiNA = colSums(data_NL[info$is_na_inf, , drop=FALSE]) > 0
                whoIsNA = nonlinear.varnames[quiNA]
                msg = msg_na_inf(info$any_na, info$any_inf)
                obs = which(info$is_na_inf)

                stop("The non-linear part (NL.fml) contains ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide data without ", msg, " (or use na_inf.rm). FYI the variable", enumerate_items(whoIsNA, "s.is"), " concerned.")

            } else {
                # If na_inf.rm => we keep track of the NAs
                if(info$any_na) ANY_NA = TRUE
                if(info$any_inf) ANY_INF = TRUE

                isNA_NL = info$is_na_inf
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_NL
                msgNA_NL = paste0(", NL: ", numberFormatNormal(sum(isNA_NL)))
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
    # Offset
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

                offset = formula(offset) # regularization

                if(length(offset) != 2){
                    stop("Argument 'offset' must be a one-sided formula (e.g.): ~ 1+x^2 ; or a numeric vector.")
                }

                offset.call = offset[[length(offset)]]
                vars.offset = all.vars(offset.call)

                if(any(!vars.offset %in% dataNames)){
                    var_missing = setdiff(vars.offset, dataNames)
                    stop("In the argument 'offset': the variable", enumerate_items(var_missing, "s.is"), " not in the data.")
                }

                offset.value = try(eval(offset.call, data), silent = TRUE)
                if("try-error" %in% class(offset.value)){
                    stop("Evaluation of the offset (equal to ", deparse_long(offset.call), ") raises and error: \n", offset.value)
                }

            } else {
                if(!is.numeric(offset) || !isVector(offset)){
                    stop("The argument 'offset' must be either a one-sided formula (e.g. ~x1), either a numeric vector.")
                }

                if(length(offset) == 1){
                    offset.value = rep(offset, nobs)
                } else if(length(offset) != nobs){
                    stop("The offset's length should be equal to the data's length (currently it's ", numberFormatNormal(length(offset)), " instead of ", numberFormatNormal(nobs), ").")
                } else {
                    offset.value = offset
                }
            }

            anyNA_offset = FALSE
            info = cpppar_which_na_inf_vec(offset.value, nthreads)
            if(info$any_na_inf){
                anyNA_offset = TRUE
                if(na_inf.rm == FALSE){
                    msg = msg_na_inf(info$any_na, info$any_inf)
                    obs = head(which(info$is_na_inf), 3)
                    stop("The offset contains ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide data without ", msg, " (or use na_inf.rm).")
                } else {
                    if(info$any_na) ANY_NA = TRUE
                    if(info$any_inf) ANY_INF = TRUE

                    isNA_offset = info$is_na_inf
                    anyNA_sample = TRUE
                    isNA_sample = isNA_sample | isNA_offset
                    msgNA_offset = paste0(", Offset: ", numberFormatNormal(sum(isNA_offset)))
                }
            } else {
                msgNA_offset = ", Offset: 0"
            }
        } else if(is.null(offset)){
            # msg if it's not what the user wanted
            if(!is.null(mc_origin$offset) && deparse_long(mc_origin$offset) != "x$offset"){
                stop("Argument 'offset' (", deparse_long(mc_origin$offset), ") is evaluated to NULL. This is likely not what you want.")
            }
        }
    }

    #
    # Weights
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

                weights = formula(weights) # regularization

                if(length(weights) != 2){
                    stop("The argument weights must be either a one-sided formula (e.g. ~x1), either a numeric vector.")
                }

                weights.call = weights[[length(weights)]]
                vars.weights = all.vars(weights.call)

                if(any(!vars.weights %in% dataNames)){
                    var_missing = vars.weights[!vars.weights %in% dataNames]
                    stop("In the argument 'weights': the variable", enumerate_items(var_missing, "s.is"), " not in the data." )
                }

                weights.value = try(eval(weights.call, data), silent = TRUE)
                if("try-error" %in% class(weights.value)){
                    stop("Evaluation of the weights (equal to ", deparse_long(weights.call), ") raises and error: \n", weights.value)
                }

            } else {
                if(!is.numeric(weights) || !isVector(weights)){
                    stop("The argument weights must be either a one-sided formula (e.g. ~x1), either a numeric vector. Currently it is not ", ifelse(!is.numeric(weights), "even numeric", "a vector despite being numeric"), ".")
                }

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

            anyNA_weights = FALSE
            info = cpppar_which_na_inf_vec(weights.value, nthreads)
            if(info$any_na_inf){
                anyNA_weights = TRUE

                if(na_inf.rm == FALSE){
                    msg = msg_na_inf(info$any_na, info$any_inf)
                    obs = head(which(info$is_na_inf), 3)
                    stop("The weights contain ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide data without ", msg, " (or use na_inf.rm).")

                } else {
                    if(info$any_na) ANY_NA = TRUE
                    if(info$any_inf) ANY_INF = TRUE

                    isNA_W = info$is_na_inf
                    anyNA_sample = TRUE
                    isNA_sample = isNA_sample | isNA_W
                    msgNA_weight = paste0(", Weights: ", numberFormatNormal(sum(isNA_W)))
                }
            } else {
                msgNA_weight = ", Weights: 0"
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

        } else if(!is.null(mc_origin$weights) && deparse_long(mc_origin$weights) != 'x[["weights"]]'){
            # we avoid this behavior
            stop("Argument 'weights' (", deparse_long(mc_origin$weights), ") is evaluated to NULL. This is likely not what you want.")
        }
    }


    #
    # Handling Clusters ####
    #

    isSlope = onlySlope = FALSE
    if(from_update){
        # Fixed-effects information coming from the update method

        # means that there is no modification of past clusters

        # we retrieve past information
        fixef_id = object$fixef_id
        obs2remove = object$obsRemoved
        fixef_removed = object$fixef_removed
        fixef_vars = object$fixef_vars
        fixef_sizes = object$fixef_sizes
        Q = length(fixef_sizes)

        # We still need to recreate some objects though
        names(fixef_id) = NULL

        if(length(obs2remove) > 0){
            lhs = lhs[-obs2remove]
        }

        fixef_names = sum_y_all = fixef_table = list()
        for(i in 1:Q){
            k = fixef_sizes[i]
            dum = fixef_id[[i]]
            sum_y_all[[i]] = cpp_tapply_vsum(k, lhs, dum)

            fixef_table[[i]] = cpp_table(k, dum)
            fixef_names[[i]] = attr(fixef_id[[i]], "fixef_names")
        }

        #
        # saving the FEs IDs and size before reordering (note: here, case with slopes never occurs)
        #

        fixef_id_res = list()
        for(i in 1:Q){
            dum = fixef_id[[i]]
            attr(dum, "fixef_names") = as.character(fixef_names[[i]])
            fixef_id_res[[fixef_vars[i]]] = dum
        }

        fixef_sizes_res = fixef_sizes
        names(fixef_sizes_res) = fixef_vars

        #
        # Re-ordering the CC
        #

        if(any(fixef_sizes != sort(fixef_sizes, decreasing = TRUE))){
            new_order = order(fixef_sizes, decreasing = TRUE)

            fixef_sizes = fixef_sizes[new_order]
            fixef_id = fixef_id[new_order]
            sum_y_all = sum_y_all[new_order]
            fixef_table = fixef_table[new_order]
        }

        # The formula with the clusters
        fml_char = as.character(fml)
        fml_full = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", paste0(fixef_vars, collapse = "+")))

    } else if(isFixef){
        # The main fixed-effects construction

        if(isFit){

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

        } else {
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
                } else if(!isLogical(combine.quick)){
                    stop("Argument 'combine.quick' must be a single logical.")
                }

                # start: NEW

                fixef_terms_full = terms_fixef(fixef_fml)
                if("try-error" %in% class(fixef_terms_full)){
                    stop("Problem extracting the terms of the fixed-effects part of the formula:\n", fixef_terms_full)
                }

                fixef_terms = fixef_terms_full$fml_terms

                # FEs
                fixef_mat = prepare_df(fixef_terms_full$fe_vars, data, combine.quick)
                if("try-error" %in% class(fixef_mat)){
                    stop("Problem evaluating the fixed-effects part of the formula:\n", fixef_mat)
                }
                fixef_vars = names(fixef_mat)

                # Slopes
                isSlope = any(fixef_terms_full$slope_flag)
                if(isSlope){

                    if(!origin_type %in% c("feols", "feglm")){
                        stop("The use of varying slopes is available only for functions feols, feglm or fepois.")
                    }

                    slope_mat = prepare_df(fixef_terms_full$slope_vars, data)
                    if("try-error" %in% class(slope_mat)){
                        stop("Problem evaluating the variables with varying slopes in the fixed-effects part of the formula:\n", slope_mat)
                    }
                    slope_flag = fixef_terms_full$slope_flag
                    slope_fe = fixef_terms_full$fe_vars
                    slope_vars = fixef_terms_full$slope_vars

                    # Further controls
                    not_numeric = !sapply(slope_mat, is.numeric)
                    if(any(not_numeric)){
                        stop("In the fixed-effects part of the formula (i.e. in ", as.character(fixef_fml[2]), "), variables with varying slopes must be numeric. Currently variable", enumerate_items(names(slope_mat)[not_numeric], "s.is"), " not.")
                    }

                    onlySlope = all(slope_flag)

                }

                # fml update
                # fml = update(fml, as.formula(paste0(".~.|", paste0(fixef_terms, collapse = "+"))))
                fml_char = as.character(fml)
                fml_full = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", paste0(fixef_terms, collapse = "+")))

                #   end: NEW

                # OLD -- deprec:
                # fixef_mat = prepare_cluster_mat(fixef_fml, data, combine.quick)
                # we change fixef_vars to become a vector of characters
                # fixef_vars = names(fixef_mat)
            }
        }

        # # We change factors to character
        # isFactor = sapply(fixef_mat, is.factor)
        # if(any(isFactor)){
        #     for(i in which(isFactor)){
        #         fixef_mat[[i]] = as.character(fixef_mat[[i]])
        #     }
        # }

        # We change non-numeric to character (impotant for parallel qufing)
        is_not_num = sapply(fixef_mat, function(x) !is.numeric(x))
        if(any(is_not_num)){
            for(i in which(is_not_num)){
                fixef_mat[[i]] = as.character(fixef_mat[[i]])
            }
        }

        if(anyNA(fixef_mat)){
            isNA_cluster = rowSums(is.na(fixef_mat)) > 0

            if(!na_inf.rm){
                # Default behavior, NA not allowed
                var_problem = fixef_vars[sapply(fixef_mat, anyNA)]
                obs = head(which(isNA_cluster), 3)
                stop("The fixed-effects variables contain NA values. Please provide data without NA (or use na_inf.rm).", ifelse(ncol(fixef_mat) > 1, paste0(" FYI the clusters with NA are: ", paste0(var_problem, collapse = ", "), "."), ""), " (e.g. observation", enumerate_items(obs, "s"), ".)")

            } else {
                # If na_inf.rm => we keep track of the NAs
                anyNA_sample = TRUE
                isNA_sample = isNA_sample | isNA_cluster
                msgNA_cluster = paste0(", Clusters: ", numberFormatNormal(sum(isNA_cluster)))
            }
        } else {
            msgNA_cluster = ", Clusters: 0"
        }

        # NAs in slopes
        if(isSlope){

            mat_tmp = matrix(as.numeric(unlist(slope_mat)), nrow(slope_mat))
            info = cpppar_which_na_inf_mat(mat_tmp, nthreads)
            if(info$any_na_inf){

                if(!na_inf.rm){
                    # Default behavior: no NA tolerance
                    quiNA = colSums(slope_mat[info$is_na_inf, , drop=FALSE]) > 0
                    whoIsNA = names(slope_mat)[quiNA]
                    msg = msg_na_inf(info$any_na, info$any_inf)
                    obs = which(info$is_na_inf)

                    stop("In the fixed-effects part of the formula, variables with varying slopes contain ", msg, " (e.g. observation", enumerate_items(obs, "s"), "). Please provide data without ", msg, " (or use na_inf.rm). FYI the variable", enumerate_items(whoIsNA, "s.is"), " concerned.")
                } else {
                    # If na_inf.rm => we keep track of the NAs
                    if(info$any_na) ANY_NA = TRUE
                    if(info$any_inf) ANY_INF = TRUE

                    isNA_slope = info$is_na_inf
                    anyNA_sample = TRUE
                    isNA_sample = isNA_sample | isNA_slope
                    msgNA_slope = paste0(", Var. Slopes: ", numberFormatNormal(sum(isNA_slope)))
                }

            } else {
                msgNA_slope = ", Var. Slopes: 0"
            }
        } else {
            msgNA_slope = ""
            fixef_terms = fixef_vars # both are identical if no slope
        }

        # Removing observations
        obs2remove_NA = c()
        if(anyNA_sample || any0W){
            # we remove all NAs obs + 0 weight obs

            nbNA = sum(isNA_sample)

            if(anyNA_sample){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                message_NA = paste0(numberFormatNormal(nbNA), " observation", plural(nbNA), " removed because of ", msg, " (Breakup: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_cluster, msgNA_slope, msgNA_offset, msgNA_weight, ").")
                notes = c(notes, message_NA)
            }

            if(nbNA == nobs){
                msg = msg_na_inf(ANY_NA, ANY_INF)
                stop("All observations contain ", msg, ". Estimation cannot be done. (Breakup: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_cluster, msgNA_slope, msgNA_offset, msgNA_weight, ")")
            }

            if(any0W){
                # 0 weight => like NAs
                isNA_sample = isNA_sample | is0W
                nbNA = sum(isNA_sample)

                if(nbNA == nobs){
                    if(anyNA_sample){
                        msg = msg_na_inf(ANY_NA, ANY_INF)
                        stop("All observations contain ", msg, " or are 0-weight. Estimation cannot be done. (0-weight: ", sum(is0W), ", breakup ", msg, ": ", msgNA_y, msgNA_L, msgNA_NL, msgNA_cluster, msgNA_slope, msgNA_offset, msgNA_weight, ")")
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

        # QUF setup ####

        Q = length(fixef_terms) # terms: contains FEs + slopes



        fixef_removed = slope_variables = list()
        obs2remove = c()

        type = switch(family, gaussian = 0, logit = 2, 1)
        do_sum_y = !origin_type %in% c("feols", "feglm")

        only_slope = FALSE
        if(isSlope){
            # only slope: not any non-slope
            only_slope = as.vector(tapply(!slope_flag, slope_fe, sum)[fixef_vars]) == 0
        }

        quf_info_all = cpppar_quf_table_sum(x = fixef_mat, y = lhs, do_sum_y = do_sum_y, type = type, only_slope = only_slope, nthreads = nthreads)

        fixef_id = quf_info_all$quf
        # names
        fixef_names = list()
        is_string = sapply(fixef_mat, is.character)
        for(i in 1:length(fixef_id)){
            if(is_string[i]){
                fixef_names[[i]] = fixef_mat[[i]][quf_info_all$items[[i]]]
            } else {
                fixef_names[[i]] = quf_info_all$items[[i]]
            }
        }
        # table/sum_y/sizes
        fixef_table = quf_info_all$table
        sum_y_all = quf_info_all$sum_y
        fixef_sizes = lengths(fixef_table)

        # If observations have been removed:
        if(!is.null(quf_info_all$obs_removed)){

            # which obs are removed
            obs2remove = which(quf_info_all$obs_removed)

            # update of the lhs
            lhs = lhs[-obs2remove]

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
            message_cluster = paste0(paste0(nb_missing, collapse = "/"), " fixed-effect", plural(sum(nb_missing)), " (", numberFormatNormal(length(obs2remove)), " observation", plural_len(obs2remove), ") removed because of only ", ifelse(family=="logit", "zero (or only one)", "zero"), " outcomes.")
            notes = c(notes, message_cluster)
        }

        # If slopes: we need to recreate some values (quf/table/sum_y)
        if(isSlope){

            # The slope variables
            for(i in which(slope_flag)){
                slope_variables[[i]] = slope_mat[[slope_vars[i]]]
                if(length(slope_variables[[i]]) == 1){
                    # si l'utilisateur utilise une constante comme variable...
                    # il faut aussi controler pour quand il fait n'importe quoi...
                    slope_variables[[i]] = rep(slope_variables[[i]], length(lhs))
                }
            }

            dict = 1:length(fixef_vars)
            names(dict) = fixef_vars
            new_id = dict[slope_fe]

            fixef_id = fixef_id[new_id]
            fixef_names = fixef_names[new_id]
            sum_y_all = sum_y_all[new_id]
            fixef_table = fixef_table[new_id]
            fixef_sizes = lengths(fixef_table)
        }

        if(length(obs2remove_NA) > 0){
            # we update the value of obs2remove (will contain both NA and removed bc of outcomes)
            if(length(obs2remove) > 0){
                obs2remove_cluster = index_noNA[obs2remove]
            } else {
                obs2remove_cluster = c()
            }

            obs2remove = sort(c(obs2remove_NA, obs2remove_cluster))
        }

        #
        # we save the cluster IDs + sizes (to be returned in "res") [original order!]
        #

        if(isSlope){
            # we only have one ID per FE used (not 1 per slope)
            index = c()
            for(var in fixef_vars){
                index = c(index, which.max(slope_fe == var))
            }
        } else {
            index = 1:Q
        }


        fixef_id_res = list()
        for(i in index){
            dum = fixef_id[[i]]
            attr(dum, "fixef_names") = as.character(fixef_names[[i]])
            if(isSlope){
                fixef_id_res[[slope_fe[i]]] = dum
            } else {
                fixef_id_res[[fixef_vars[i]]] = dum
            }
        }

        fixef_sizes_res = fixef_sizes[index]
        names(fixef_sizes_res) = fixef_vars

        #
        # We re-order the clusters
        #

        IS_REORDER = FALSE
        if(any(fixef_sizes != sort(fixef_sizes, decreasing = TRUE))){
            IS_REORDER = TRUE
            # FE with the most cases first (limits precision problems)

            new_order = order(fixef_sizes, decreasing = TRUE)

            fixef_sizes = fixef_sizes[new_order]
            fixef_id = fixef_id[new_order]
            sum_y_all = sum_y_all[new_order]
            fixef_table = fixef_table[new_order]

            if(isSlope){
                slope_variables = slope_variables[new_order]
            }

        }

    } else {
        # There is no fixed-effect
        Q = 0

        # NA management is needed to create obs2remove
        if(anyNA_sample || any0W){

            nbNA = sum(isNA_sample)

            if(anyNA_sample){
                message_NA = paste0(numberFormatNormal(nbNA), " observation", plural(nbNA), " removed because of NA values (Breakup: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_offset, msgNA_weight, ").")
                notes = c(notes, message_NA)

                if(nbNA == nobs){
                    stop("All observations contain NAs. Estimation cannot be done. (Breakup: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_cluster, msgNA_offset, msgNA_weight, ")")
                }
            }

            if(any0W){
                # 0 weight => like NAs
                isNA_sample = isNA_sample | is0W
                nbNA = sum(isNA_sample)

                if(nbNA == nobs){
                    if(anyNA_sample){
                        stop("All observations are either NA or 0-weight. Estimation cannot be done. (0-weight: ", sum(is0W), ", breakup NA: ", msgNA_y, msgNA_L, msgNA_NL, msgNA_offset, msgNA_weight, ")")
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

        # if(isNL || (isLinear && useModel.matrix)){
        #     data = data[-obs2remove, , drop = FALSE]
        # }
        #
        # # We recreate the linear matrix and the LHS
        # if(isLinear) {
        #     if(useModel.matrix){
        #         # means there are factors
        #         # linear.mat = stats::model.matrix(linear.fml, data)
        #         linear.mat = fixest_model_matrix(linear.fml, data)
        #     } else {
        #         linear.mat = linear.mat[-obs2remove, , drop = FALSE]
        #     }
        # }

        # We drop only 0 variables (may happen for factors)
        linear.mat = linear.mat[-obs2remove, , drop = FALSE]

        if(useModel.matrix){
            # There are factors => possibly some vars are only 0 now that NAs are removed

            only_0 = cpppar_check_only_0(linear.mat, nrow(linear.mat), nthreads)
            if(all(only_0 == 1)){
                stop("After removing NAs, not a single explanatory variable is different from 0.")
            } else if(any(only_0 == 1)){
                linear.mat = linear.mat[, only_0 == 0, drop = FALSE]
            }
        }

        if(Q == 0){
            # if Q > 0: done already when managing the clusters
            lhs = lhs[-obs2remove]
        }

        if(isOffset){
            offset.value = offset.value[-obs2remove]
        }

        if(isWeight){
            weights.value = weights.value[-obs2remove]
        }

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


            # controls
            if(!is.numeric(starting_values) || !isVector(starting_values)){
                stop("Argument '", arg_name, "' must be a numeric vector.")
            }
            starting_values = as.numeric(as.vector(starting_values)) # avoids integers

            if(length(starting_values) == 1){
                starting_values = rep(starting_values, length(lhs))
            } else if(length(starting_values) == length(lhs)){
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

    # families of feglm
    if(origin_type == "feglm"){

        family_equiv = family_funs$family_equiv

        if(family_equiv == "poisson"){
            family_funs$linkfun = function(mu) cpppar_log(mu, nthreads)
            family_funs$linkinv = function(eta) cpppar_poisson_linkinv(eta, nthreads)
            y_pos = lhs[lhs > 0]
            qui_pos = lhs > 0
            if(isWeight){
                constant = sum(weights.value[qui_pos] * y_pos * cpppar_log(y_pos, nthreads) - weights.value[qui_pos] * y_pos)
                dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(wt[qui_pos] * y_pos * eta[qui_pos]) + sum(wt * mu))
            } else {
                constant = sum(y_pos * cpppar_log(y_pos, nthreads) - y_pos)
                dev.resids = function(y, mu, eta, wt) 2 * (constant - sum(y_pos * eta[qui_pos]) + sum(mu))
            }

            family_funs$dev.resids = dev.resids

            family_funs$mu.eta = function(mu, eta) mu
            family_funs$validmu = function(mu) cpppar_poisson_validmu(mu, nthreads)

        } else if(family_equiv == "logit"){
            family_funs$linkfun = function(mu) cpppar_logit_linkfun(mu, nthreads)
            family_funs$linkinv = function(eta) cpppar_logit_linkinv(eta, nthreads)
            if(isWeight){
                dev.resids = function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, wt, nthreads))
            } else {
                dev.resids = function(y, mu, eta, wt) sum(cpppar_logit_devresids(y, mu, 1, nthreads))
            }
            family_funs$dev.resids = dev.resids

            family_funs$mu.eta = function(mu, eta) cpppar_logit_mueta(eta, nthreads)
        }

        assign("family_funs", family_funs, env)

    }


    if(lparams == 0 && Q == 0) stop("No parameter to be estimated.")

    if(!isLogical(useHessian)){
        stop("'useHessian' must be a single 'logical'!")
    }
    assign("hessian.args", hessian.args, env)

    if(origin == "feNmlm"){

        if(!isSingleChar(jacobian.method)){
            stop("Argument 'jacobian.method' must be a character scalar equal to 'simple' or 'Richardson'.")
        } else {
            value = try(match.arg(jacobian.method, c("simple", "Richardson")))
            if("try-error" %in% class(value)){
                stop("Argument jacobian.method does not match 'simple' or 'Richardson' (currently equal to ", jacobian.method, ").")
            }
            jacobian.method = value
        }

    }

    #
    # Controls: The non linear part
    #

    if(isNL){
        if(missing(NL.start.init)){
            if(missing(NL.start)) stop("There must be starting values for NL parameters. Please use argument NL.start (or NL.start.init).")
            if(typeof(NL.start) != "list") stop("NL.start must be a list.")
            if(any(!nonlinear.params %in% names(NL.start))) stop(paste("Some NL parameters have no starting values:\n", paste(nonlinear.params[!nonlinear.params %in% names(NL.start)], collapse=", "), ".", sep=""))

            # we restrict NL.start to the nonlinear.params
            NL.start = NL.start[nonlinear.params]
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

    # OFFSET
    assign("offset.value", offset.value, env)

    # WEIGHTS
    assign("weights.value", weights.value, env)

    #
    # PRECISION + controls ####
    #

    # The main precision

    if(!isScalar(fixef.tol) || fixef.tol <= 0 || fixef.tol >1){
        stop("If provided, argument 'fixef.tol' must be a strictly positive scalar lower than 1.")
    } else if(fixef.tol < 10000*.Machine$double.eps){
        stop("Argument 'fixef.tol' cannot be lower than ", signif(10000*.Machine$double.eps))
    }

    if(!isScalar(fixef.iter) || fixef.iter < 1){
        stop("Argument fixef.iter must be an integer greater than 0.")
    }

    if(origin_type == "feNmlm"){
        if(!isScalar(deriv.iter) || deriv.iter < 1){
            stop("Argument deriv.iter must be an integer greater than 0.")
        }

        if(!isScalar(deriv.tol) || deriv.tol <= 0 || deriv.tol >1){
            stop("If provided, argument 'deriv.tol' must be a strictly positive scalar lower than 1.")
        } else if(deriv.tol < 10000*.Machine$double.eps){
            stop("Argument 'deriv.tol' cannot be lower than ", signif(10000*.Machine$double.eps))
        }

        # Other: opt.control
        if(!is.list(opt.control)){
            stop("Argument opt.control must be a list of controls for the function nlminb (see help for details).")
        }

    }

    if(origin_type == "feglm"){

        if(!isScalar(glm.iter) || glm.iter < 1){
            stop("Argument glm.iter must be an integer greater than, or equal to, 1.")
        }

        if(!isScalar(glm.tol) || glm.tol <= 0 || glm.tol > 1){
            stop("Argument 'glm.tol' must be a strictly positive scalar lower than 1.")
        } else if(glm.tol < 1000*.Machine$double.eps){
            stop("Argument 'glm.tol' cannot be lower than ", signif(1000*.Machine$double.eps))
        }

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

    #
    # The MODEL0 => to get the init of the theta for the negbin
    #

    # Bad location => rethink the design of the code
    assign("famFuns", famFuns, env)
    assign("family", family, env)
    assign("nobs", length(lhs), env)
    assign("lhs", lhs, env)
    assign("nthreads", nthreads, env)

    if(missing(theta.init)){
        theta.init = NULL
    } else {
        if(!isScalar(theta.init) || theta.init <= 0){
            stop("the argument 'theta.init' must be a strictly positive scalar.")
        }
    }

    if(origin_type == "feNmlm" || computeModel0){
        model0 <- get_model_null(env, theta.init)
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

    # On balance les donnees a utiliser dans un nouvel environnement
    if(isLinear) assign("linear.mat", linear.mat, env)

    ####
    #### Sending to the env ####
    ####

    useExp_clusterCoef = family %in% c("poisson")

    # The dummies
    assign("isFixef", isFixef, env)
    if(isFixef){
        assign("fixef_id", fixef_id, env)
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

        # the saved sumFE
        if(!missnull(sumFE_init)){
            # information on starting values coming from update method
            doExp = ifelse(useExp_clusterCoef, exp, I)

            # Means it's the full fixef properly given
            assign("saved_sumFE", doExp(sumFE_init), env)

        } else if(useExp_clusterCoef){
            assign("saved_sumFE", rep(1, length(lhs)), env)
        } else {
            assign("saved_sumFE", rep(0, length(lhs)), env)
        }

        # New cpp functions
        # This is where we send the elements needed for convergence in cpp
        assign("fixef_id_vector", as.integer(unlist(fixef_id) - 1), env)
        assign("fixef_table_vector", as.integer(unlist(fixef_table)), env)
        assign("sum_y_vector", unlist(sum_y_all), env)

        if(origin_type == "feNmlm"){
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

        if(isSlope){
            assign("slope_flag", as.integer(slope_flag), env)
            assign("slope_variables", as.numeric(unlist(slope_variables)), env)
        } else {
            assign("slope_flag", rep(0L, length(fixef_vars)), env)
            assign("slope_variables", 0, env)
        }

    }

    # basic NL
    envNL = new.env()
    assign("isNL", isNL, env)
    if(isNL){
        for(var in nonlinear.varnames) assign(var, data[[var]], envNL)
        for(var in nonlinear.params) assign(var, start[var], envNL)
    }
    assign("envNL", envNL, env)
    assign("nl.call", nl.call, env)
    # NO user defined gradient: too complicated, not really efficient
    assign("isGradient", FALSE, env)
    assign("lower", lower, env)
    assign("upper", upper, env)

    # other
    assign("lhs", lhs, env)
    assign("isLinear", isLinear, env)
    assign("linear.params", linear.params, env)
    assign("nonlinear.params", nonlinear.params, env)
    assign("params", params, env)
    assign("nobs", length(lhs), env)
    assign("jacobian.method", jacobian.method, env)
    assign("famFuns", famFuns, env)
    assign("family", family, env)
    assign("iter", 0, env)
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
    assign("fixef.tol", fixef.tol, env)
    assign("NR.tol", NR.tol, env)
    assign("deriv.tol", deriv.tol, env)
    # ITERATIONS
    assign("fixef.iter", fixef.iter, env)
    assign("deriv.iter", deriv.iter, env)
    assign("fixef.iter.limit_reached", 0, env) # for warnings if max iter is reached
    assign("deriv.iter.limit_reached", 0, env) # for warnings if max iter is reached
    # OTHER
    assign("useAcc", TRUE, env)
    assign("warn_0_Hessian", FALSE, env)
    assign("warn_overfit_logit", FALSE, env)
    # Misc
    assign("origin", origin, env)
    assign("warn", warn, env)
    assign("opt.control", opt.control, env)


    # To monitor how the FEs are computed (if the problem is difficult or not)
    assign("firstIterCluster", 1e10, env) # the number of iterations in the first run
    assign("firstRunCluster", TRUE, env) # flag for first enrty in get_dummies
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
        if(!isVector(mu) || !is.numeric(mu)){
            stop("Evaluation of NL.fml should return a numeric vector. (This is currently not the case.)")
        }

        # Handling NL.fml errors
        if(length(mu) != nrow(data)){
            stop("Evaluation of NL.fml leads to ", length(mu), " observations while there are ", nrow(data), " observations in the data base. They should be of the same lenght.")
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
    if(isLinear) {
        # mu <- mu + c(linear.mat%*%unlist(start[linear.params]))
        start_coef_linear = unlist(start[linear.params])
        mu <- mu + cpppar_xbeta(linear.mat, start_coef_linear, nthreads)
    }

    # Mise en place du calcul du gradient
    if(origin_type == "feNmlm"){
        gradient = femlm_gradient
        hessian <- NULL
        if(useHessian) hessian <- femlm_hessian
        assign("gradient", gradient, env)
        assign("hessian", hessian, env)
    }

    assign("fml", fml, env)
    assign("model0", model0, env)

    onlyFixef = !isLinear && !isNonLinear && Q > 0
    assign("onlyFixef", onlyFixef, env)
    if(onlyFixef){
        assign("linear.mat", 0, env)
    }

    assign("start", start, env)

    #
    # Res ####
    #

    #
    # Preparation of the results list (avoid code repetition in estimation funs)
    #

    res = list(nobs=length(lhs), fml=fml, call = mc_origin, method = origin)

    if(isFixef) res$fml_full = fml_full

    if(isFit) res$fromFit = TRUE

    # nber of params
    K = length(params)
    if(isFixef){
        if(isSlope){
            K = K + sum(fixef_sizes) - sum(slope_flag == FALSE) + any(slope_flag == FALSE)
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

            # We also add the variables => makes fixef self-contained
            # It looks odd but the user may ask for the same var to have
            # varying slopes in svl dimensions
            sv = list()
            slope_vars_unik = slope_vars[!is.na(slope_vars)]
            if(IS_REORDER) slope_variables = slope_variables[order(new_order)]
            for(var in slope_vars_unik){
                sv[[var]] = slope_variables[[which.max(slope_vars == var)]]
            }

            res$slope_variables = sv
        }

        res$fixef_id = fixef_id_res
        res$fixef_sizes = fixef_sizes_res

    }

    # Observations removed (either NA or clusters)
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
    if(origin_type == "feglm" && isFixef){
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


