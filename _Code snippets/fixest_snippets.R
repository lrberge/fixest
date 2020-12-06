#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Sep 22 15:30:09 2020
# ~: Code snippets to share
#----------------------------------------------#


####
#### Marginal effects ####
####

# Quick and Dirty implementation of marginal effects
# Very limited but does the job. Easy to expand.

require(dreamerr) ; require(fixest)
meffect = function(x, at_means = TRUE, se, cluster, ...){
    # x: fixest object

    check_arg(x, "class(fixest) mbt")
    check_arg(at_means, "logical scalar")
    if(isFALSE(at_means)) stop("Sorry, so far, only the marginal effects at means is available.")

    if(!isTRUE(x$summary) || !missing(se) || !missing(cluster)){
        x = summary(x, se = se, cluster = cluster, ...)
    }

    coef = x$coefficients
    vars = names(coef)

    vcov = x$cov.scaled

    m = model.matrix(x)
    m_mean = colMeans(m[, vars])
    mu_mean = as.vector(m_mean %*% coef)

    # 1) the SE
    # 2) the ME with formatting

    # 1) The standard errors via the delta method (at means)
    if(x$method_type == "feols" || (x$method_type == "feNmlm" && x$family == "gaussian")){
        d_f = function(x) 1
        dd_f = function(x) 0
    } else if((x$method_type == "feNmlm" && x$family == "gaussian") || (x$method_type == "feglm" && x$family$family == "poisson")){
        d_f = function(x) exp(x)
        dd_f = function(x) exp(x)
    } else if(x$method_type == "feglm" && x$family$family == "binomial" && x$family$link == "probit"){
        d_f = dnorm
        dd_f = function(x) 1/(sqrt(2*pi)) * -x * exp(-x**2/2)
    } else {
        stop("Current family is not supported.")
    }

    # jacobian from the formula
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

# I set the default SE to "standard" for comparability with Stata
setFixest_se(no_FE = "standard")

res = fepois(y ~ x1 + x2, base)
meffect(res)

res = feglm(1*(y > 5) ~ x1 + x2, base, family = binomial("probit"))
meffect(res)




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
















































