#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Sep 22 15:30:09 2020
# ~: Code snippets to share
#----------------------------------------------#


####
#### Sparse model.matrix ####
####

# Here's a quick and dirty (hence limited) code to create a sparse model.matrix from
# a fixest estimation

# What does it do? It returns a list of two elements:
# - mat_RHS: the *sparse* version of the RHS (excluding the FEs)
# - mat_FE: the sparse matrix of the fixed-effects
# NOTE: both are Matrix matrices from library(Matrix)

# Limitations:
# - does not automatically remove references from factors
# - does not handle estimations containing variables with varying slopes
#

# Benefits (yes there are some!):
# - avoids copies of the data in the process of constructing the matrix
# - handles interactions between any number of factors/numeric variables
# - the non-sparse model matrix is never created, hence it's very fast
#

# Be sure to run the *three* functions sparse_model_matrix, mult_wrap and mult_sparse


# x: fixest estimation
sparse_model_matrix = function(x){

    require(Matrix)

    # Linear formula
    fml_lin = formula(x, "lin")

    data = fixest:::fetch_data(x)
    data = as.data.frame(data)
    data = data[obs(x), ]

    #
    # Step 1: Linear matrix
    #

    vars = attr(terms(fml_lin), "term.labels")

    if(length(vars) == 0){
        # Case only FEs
        mat = NULL
    } else {

        # Since we don't want to evaluate the factors,
        # the code is a bit intricate because we have to catch them before
        # any interaction takes place
        #
        # that's why I wrap interactions in a function (mult_sparse())
        #

        # Below, we evaluate all the variables in a "sparse" way

        vars_calls = lapply(vars, mult_wrap)

        n = length(vars)
        variables_list = vector("list", n)
        for(i in 1:n){
            variables_list[[i]] = eval(vars_calls[[i]], data)
        }

        # To create the sparse matrix, we need the indexes

        total_cols = 0
        running_cols = c(0)
        for(i in 1:n){
            xi = variables_list[[i]]
            if(inherits(xi, "sparse_var")){
                total_cols = total_cols + xi$n_cols
            } else {
                total_cols = total_cols + NCOL(xi)
            }
            running_cols[i + 1] = total_cols
        }

        # We just create a sparse matrix and fill it

        # 1) creating the indexes + names

        # NOTA: I use lists to avoid creating copies
        rowid = 1:nrow(data)
        id_all = values_all = names_all = vector("list", n)
        for(i in 1:n){
            xi = variables_list[[i]]
            if(inherits(xi, "sparse_var")){
                id_all[[i]] = cbind(xi$rowid, running_cols[i] + xi$colid)
                values_all[[i]] = xi$values
                names_all[[i]] = paste0(vars[[i]], "::", xi$col_names)
            } else if(NCOL(xi) == 1){
                id_all[[i]] = cbind(rowid, running_cols[i] + 1)
                values_all[[i]] = xi
                names_all[[i]] = vars[[i]]
            } else {
                colid = rep(1:NCOL(xi), each = nrow(data))
                id_all[[i]] = cbind(rep(rowid, NCOL(xi)), running_cols[i] + colid)
                values_all[[i]] = as.vector(xi)
                if(!is.null(colnames(xi))){
                    names_all[[i]] = paste0(vars[[i]], colnames(xi))
                } else {
                    names_all[[i]] = paste0(vars[[i]], 1:NCOL(xi))
                }
            }
        }

        id_mat = do.call(rbind, id_all)
        values_vec = unlist(values_all)
        names_vec = unlist(names_all)

        # 2) filling the matrix: one shot, no copies

        mat = Matrix(0, nrow(data), total_cols, dimnames = list(NULL, names_vec))
        mat[id_mat] = values_vec
    }

    #
    # Step 2: the fixed-effects
    #

    if(length(x$fixef_id) == 0){
        mat_FE = NULL
    } else {
        # Same process, but easier

        rowid = 1:nrow(data)
        total_cols = sum(x$fixef_sizes)
        running_cols = c(0, x$fixef_sizes)
        n_FE = length(x$fixef_sizes)
        id_all = names_all = vector("list", n_FE)
        for(i in 1:n_FE){
            xi = x$fixef_id[[i]]
            id_all[[i]] = cbind(rowid, running_cols[i] + xi)
            names_all[[i]] = paste0(names(x$fixef_id)[i], "::", attr(xi, "fixef_names"))
        }

        id_mat = do.call(rbind, id_all)
        names_vec = unlist(names_all)

        mat_FE = Matrix(0, nrow(data), total_cols, dimnames = list(NULL, names_vec))
        mat_FE[id_mat] = 1
    }

    res = list(mat_RHS = mat, mat_FE = mat_FE)

    res
}

# Internal: modifies the calls so that each variable/interaction is evaluated with mult_sparse
mult_wrap = function(x){
    # x: character string of a variable to be evaluated
    # ex: "x1" => mult_sparse(x1)
    #     "x1:factor(x2):x3" => mult_sparse(x3, factor(x2), x1)
    #
    # We also add the argument sparse to i()
    #     "x1:i(species, TRUE)" => mult_sparse(x1, i(species, TRUE, sparse = TRUE))

    x_call = str2lang(x)

    res = (~ mult_sparse())[[2]]

    if(length(x_call) == 1 || x_call[[1]] != ":"){
        res[[2]] = x_call

    } else {
        res[[2]] = x_call[[3]]
        tmp = x_call[[2]]

        while(length(tmp) == 3 && tmp[[1]] == ":"){
            res[[length(res) + 1]] = tmp[[3]]
            tmp = tmp[[2]]
        }

        res[[length(res) + 1]] = tmp
    }

    # We also add sparse to i() if found
    for(i in 2:length(res)){
        ri = res[[i]]
        if(length(ri) > 1 && ri[[1]] == "i"){
            ri[["sparse"]] = TRUE
            res[[i]] = ri
        }
    }

    if(length(res) > 2){
        # we restore the original order
        res[-1] = rev(res[-1])
    }

    return(res)
}

# Internal function to evaluate the variables (and interactions) in a sparse way
mult_sparse = function(...){
    # Only sparsifies factor variables
    # Takes care of interactions

    dots = list(...)
    n = length(dots)

    num_var = NULL
    factor_list = list()
    info_i = NULL
    is_i = is_factor = FALSE
    # You can't have interactions between i and factors, it's either

    for(i in 1:n){
        xi = dots[[i]]
        if(is.numeric(xi)){
            # We stack the product
            num_var = if(is.null(num_var)) xi else xi * num_var
        } else if(inherits(xi, "i_sparse")){
            is_i = TRUE
            info_i = xi
        } else {
            is_factor = TRUE
            factor_list[[length(factor_list) + 1]] = xi
        }
    }

    if(!is_i && !is_factor){
        return(num_var)
    }

    if(is_factor){
        factor_list$add_items = TRUE
        factor_list$items.list = TRUE

        fact_as_int = do.call(to_integer, factor_list)

        values = if(is.null(num_var)) rep(1, length(fact_as_int$x)) else num_var

        rowid = seq_along(values)
        res = list(rowid = rowid, colid = fact_as_int$x, values = values,
                   col_names = fact_as_int$items, n_cols = length(fact_as_int$items))
    } else {

        values = info_i$values
        if(!is.null(num_var)){
            num_var = num_var[info_i$rowid]
            values = values * num_var
        }

        res = list(rowid = info_i$rowid, colid = info_i$colid, values = values,
                   col_names = info_i$col_names, n_cols = length(info_i$col_names))
    }

    class(res) = "sparse_var"

    res
}

library(fixest)
est = feols(mpg ~ poly(hp, 2) + as.factor(gear)*as.factor(carb) + i(am) | cyl, mtcars)

sparse_model_matrix(est)


####
#### Marginal effects ####
####

# Quick and Dirty implementation of marginal effects
# Very limited but does the job. Easy to expand.

require(dreamerr) ; require(fixest)
meffect = function(x, at_means = TRUE, vcov, se, cluster, ...){
    # x: fixest object

    check_arg(x, "class(fixest) mbt")
    check_arg(at_means, "logical scalar")
    if(isFALSE(at_means)) stop("Sorry, so far only the marginal effects at means is available.")

    is_FE = FALSE
    if(!is.null(x$fixef_vars)){
        if(!x$method_type == "feglm"){
            stop("Models with fixed-effects are not yet supported for method ", x$method, ", use factors instead (if possible).")
        }
        is_FE = TRUE
    }

    if(!isTRUE(x$summary) || !missing(se) || !missing(cluster)){
        x = summary(x, vcov = vcov, se = se, cluster = cluster, ...)
    }

    coef = x$coefficients
    vars = names(coef)

    vcov = x$cov.scaled

    if(!is_FE){
        m = model.matrix(x)
        m_mean = colMeans(m[, vars])

    } else {
        # x => feglm
        if(x$family$family %in% c("poisson", "binomial")){
            m = x$scores / (x$working_residuals * x$irls_weights)
        } else {
            m = x$scores / (x$working_residuals * x$irls_weights / x$dispersion)
        }

        m_mean = colMeans(m)
    }

    mu_mean = mean(predict(x, type = "link"))

    # 1) the SE
    # 2) the ME with formatting

    # 1) The standard errors via the delta method (at means)
    if(x$method_type == "feols" || (x$method_type == "feNmlm" && x$family == "gaussian")){
        d_f = function(x) 1
        dd_f = function(x) 0
    } else if((x$method_type == "feNmlm" && x$family == "poisson") || (x$method_type == "feglm" && x$family$family == "poisson")){
        d_f = function(x) exp(x)
        dd_f = function(x) exp(x)
    } else if(x$method_type == "feglm" && x$family$family == "binomial" && x$family$link == "probit"){
        d_f = dnorm
        dd_f = function(x) 1/(sqrt(2*pi)) * -x * exp(-x**2/2)
    } else {
        stop("Current family is not supported.")
    }

    # jacobian
    jac = tcrossprod(coef, m_mean) * dd_f(mu_mean)
    diag(jac) = diag(jac) + d_f(mu_mean)
    vcov_new = jac %*% vcov %*% t(jac)

    sd = sqrt(diag(vcov_new))
    value = coef * d_f(mu_mean)
    z = value/sd

    res = data.frame(var = vars, dydx = value, se = sd, z = z, pvalue = 2 * pnorm(-abs(z)), CI_95_low = value - 1.96*sd, CI_95_high = value + 1.96*sd, row.names = seq_along(vars))

    res
}

#
# Example
#

# Note:
# - check with Stata's: margins, dydx(x1 x2) atmeans

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")

# I set the default SE to "iid" for comparability with Stata
setFixest_vcov(all = "iid")

res = fepois(y ~ x1 + x2, base)
meffect(res)

res = feglm(1*(y > 5) ~ x1 + x2, base, family = binomial("probit"))
meffect(res)

#
# FE example

# With fixed-effects as factor
res_lin = fepois(y ~ x1 + x2 + species, base)
meffect(res_lin)

# With fixed-effects absorbed
res_fe = fepois(y ~ x1 + x2 | species, base)
meffect(res_fe)




####
#### Imai & Song ####
####

# Equivalence between the "unit" wfe estimation and the same estimation in fixest
# A) function to compute the unit weights
# B) comparaison of the estimations
#
# other:
#  - AJPS 2019 paper: https://imai.fas.harvard.edu/research/files/FEmatch.pdf
#  - original code for unit weights (function GenWeightsUnit, line 1043): https://github.com/insongkim/wfe/blob/master/src/wfe.c

#
# A) weights
#

# Function to compute the weights (with full checks)

require(dreamerr) ; require(data.table)
IS_weights_unit = function(index, base){

    # ARGUMENT CHECK #

    check_arg(index, "mbt data.frame ncol(2) | os formula | character vector len(2)",
              .message = "The argument 'index' must contain the a) unit and b) the treatment variables. It can be either: i) a data.frame of two columns, ii) a one sided formula, or iii) a character vector of length 2.")

    if(!is.data.frame(index)){
        check_arg(base, "data.frame mbt", .message = "If the argument 'index' is not a data.frame, you must provide a data.frame in the argument 'base'.")

        if("formula" %in% class(index)){
            check_arg(index, "formula var(data)", .data = base)

            index = model.frame(index, base)

            check_arg(index, "data.frame ncol(2)", .message = "If a formula, then 'index' must lead to a data.frame with two variables.")

        } else {
            # flexible matching of names
            check_arg_plus(index, "multi match", .choices = names(base), .message = "If the argument 'index' is a character vector, it must match (at least partially) the names of the data.frame in argument 'base'.")

            index = as.data.frame(base)[, index]
        }

    }

    # Here 'index' is a data.frame with two variables
    res = as.data.table(index)
    names(res) = c("unit", "treatment")

    res[, n_treated := sum(treatment), by = unit]
    res[, n_not_treated := sum(1 - treatment), by = unit]
    res[, m_size := treatment * n_not_treated + (1 - treatment) * n_treated]
    res[, w := 1 + treatment * (m_size / n_treated) + (1 - treatment) * (m_size / n_not_treated)]
    res[is.na(w), w := 0]

    res$w
}


#
# B) Estimations
#

library(fixest) ; library(wfe)

# Data preparation

data(base_did)
set.seed(0)
base = data.table(base_did)
# we trim first/last periods of some guys to add variation in the weights (otherwise: all equal to 2 or 0)
id_first = sample(108, 20) ; id_last = sample(108, 20)
base = base[!(id %in% id_first & period <= 3) & !(id %in% id_last & period >= 7)]
base[, treat_post := treat*post]

# Estimations

system.time(res_wfe <- wfe(y ~ treat_post + x1, data = base, treat = "treat_post", unit.index = "id", time.index = "period", method = "unit",
                           qoi = "ate", hetero.se=TRUE, auto.se=TRUE))


system.time(res_feols <- feols(y ~ treat_post + x1 | id, base, weights = w <- IS_weights_unit(~ id + treat_post, base)))

# Comparing the weights
all.equal(w, summary(res_wfe)$W$W.it)

# Timings on my laptop:
#    wfe: 320ms
# fixest:  10ms

# Coefficients and SEs:
rbind(fixest = coef(res_feols), wfe = coef(res_wfe))
rbind(fixest = se(res_feols), wfe = se(res_wfe))
















































