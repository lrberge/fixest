#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Mon Jun 01 08:21:26 2020
# ~: Secondary methods
#----------------------------------------------#


#' Extracts the weights from a fixest object
#'
#' Simply extracts the weights used to estimate a \code{fixest} model.
#'
#' @param object A \code{fixest} object.
#' @param ... Not currently used.
#'
#' @return
#' Returns a vector of the same length as the number of observations in the original data set. Ignored observations due to NA or perfect fit are re-introduced and their weights set to NA.
#'
#' @seealso
#' \code{\link[fixest]{feols}}, \code{\link[fixest]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width, iris, weights = ~as.integer(Sepal.Length) - 3.99)
#' weights(est)
#'
weights.fixest = function(object, ...){
    w = object[["weights"]]

    # To comply with stats default
    if(is.null(w)) return(NULL)

    w = fill_with_na(w, object$obsRemoved)

    w
}



#' Residual standard deviation of fixest estimations
#'
#' Extract the estimated standard deviation of the errors from \code{fixest} estimations.
#'
#' @inheritParams weights.fixest
#'
#' @return
#' Returns a numeric scalar.
#'
#' @seealso
#' \code{\link[fixest]{feols}}, \code{\link[fixest]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
#'
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width, iris)
#' sigma(est)
#'
#'
#'
sigma.fixest = function(object, ...){
    sqrt(deviance(object) / (object$nobs - object$nparams))
}


#' Extracts the deviance of a fixest estimation
#'
#' Returns the deviance from a \code{fixest} estimation.
#'
#' @inheritParams weights.fixest
#'
#' @return
#' Returns a numeric scalar equal to the deviance.
#'
#' @seealso
#' \code{\link[fixest]{feols}}, \code{\link[fixest]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width, iris)
#' deviance(est)
#'
#' est_pois = fepois(Petal.Length ~ Petal.Width, iris)
#' deviance(est_pois)
#'
deviance.fixest = function(object, ...){

    method = object$method
    family = object$family
    r = object$residuals
    w = object[["weights"]]
    if(is.null(w)) w = rep(1, length(r))

    if(method == "feols" || (method %in% c("femlm", "feNmlm") && family == "gaussian")){
        res = sum(w * r**2)

    } else if(method %in% c("fepois", "feglm")){
        res = object$deviance

    } else {
        mu = object$fitted.values
        theta = ifelse(family == "negbin", object$theta, 1)

        # dev.resids function
        if(family == "poisson"){
            dev.resids = poisson()$dev.resids

        } else if(family == "logit"){
            dev.resids = binomial()$dev.resids

        } else if(family == "negbin"){
            dev.resids = function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))

        }

        y = r + mu

        res = sum(dev.resids(y, mu, w))
    }

    res
}



hatvalues.fixest = function(model, ...){
    # Only works for feglm/feols objects
    # to integrate



}

# #' Extracts the scores from a fixest estimation
# #'
# #' Extracts the scores from a fixest estimation.
# #'
# #' @param x A \code{fixest} object, obtained for instance from \code{\link[fixest]{feols}}.
# #' @param ... Not currently used.
# #'
# #' @return
# #' Returns a matrix of the same number of rows as the number of observations used for the estimation, and the same number of columns as there were variables.
# #'
# #' @examples
# #'
# #' est = feols(Petal.Length ~ Petal.Width + Sepal.Width, iris)
# #' head(estfun(est))
# #'
# estfun.fixest = function(x, ...){
#     # The scores are an object always contained in fixest estimations
#     x$scores
# }












































































