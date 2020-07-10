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
#' \code{\link[fixest]{feols}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
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
#' \code{\link[fixest]{feols}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
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
#' \code{\link[fixest]{feols}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest]{feglm}}, \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{feNmlm}}.
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



#' Hat values for \code{fixest} objects
#'
#' Computes the hat values for \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}} estimations. Only works when there are no fixed-effects.
#'
#' @param model A fixest object. For instance from feols or feglm.
#' @param ... Not currently used.
#'
#' @details
#' Hat values are not available for \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{femlm}} and \code{\link[fixest]{feNmlm}} estimations.
#'
#' When there are fixed-effects, the hat values of the reduced form are different from the hat values of the full model. And we cannot get costlessly the hat values of the full model from the reduced form. It would require to reestimate the model with the fixed-effects as regular variables.
#'
#' @return
#' Returns a vector of the same length as the number of observations used in the estimation.
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Width, iris)
#' head(hatvalues(est))
#'
#'
hatvalues.fixest = function(model, ...){
    # Only works for feglm/feols objects + no fixed-effects
    # When there are fixed-effects the hatvalues of the reduced form is different from
    #  the hatvalues of the full model. And we cannot get costlessly the hatvalues of the full
    #  model from the reduced form. => we need to reestimate the model with the FEs as
    #  regular variables.

    validate_dots()

    method = model$method
    family = model$family

    if(method == "feols"){

        if(!is.null(model$fixef_id)){
            stop("'hatvalues' is not implemented when the estimation contains fixed-effects.")
        }

        X = model.matrix(model)

        res = cpp_diag_XUtX(X, model$cov.unscaled / model$sigma2)

    } else if(method %in% c("fepois", "feglm")){

        if(!is.null(model$fixef_id)){
            stop("'hatvalues' is not implemented when the estimation contains fixed-effects.")
        }

        XW = model.matrix(model) * sqrt(model$irls_weights)
        res = cpp_diag_XUtX(XW, model$cov.unscaled)

    } else {
        stop("'hatvalues' is not currently implemented for function ", method, ".")
    }

    res
}

#' Extracts the scores from a fixest estimation
#'
#' Extracts the scores from a fixest estimation.
#'
#' @param x A \code{fixest} object, obtained for instance from \code{\link[fixest]{feols}}.
#' @param ... Not currently used.
#'
#' @return
#' Returns a matrix of the same number of rows as the number of observations used for the estimation, and the same number of columns as there were variables.
#'
#' @examples
#'
#' est = feols(Petal.Length ~ Petal.Width + Sepal.Width, iris)
#' head(estfun(est))
#'
estfun.fixest = function(x, ...){
    # 'scores' is an object always contained in fixest estimations
    x$scores
}


#' Functions exported from \pkg{sandwich} to implement \pkg{fixest} methods
#'
#' The package \pkg{fixest} does not use \code{estfun} from \pkg{sandwich}, but this method has been implemented to allow users to leverage the variances from \pkg{sandwich}.
#'
#' \itemize{
#' \item Here is the help from package \pkg{sandwich}: \code{\link[sandwich:estfun]{estfun}}. The help from package \pkg{fixest} is here: \code{\link[fixest]{estfun.fixest}}.
#' }
#'
#'
#' @name estfun_reexported
#' @keywords internal
NULL

#' @rdname estfun_reexported
#' @name estfun
NULL


bread.fixest = function(x, ...){
    validate_dots()

    method = x$method
    family = x$family

    if(method == "feols"){

        res = x$cov.unscaled / x$sigma2 * x$nobs

    } else if(method %in% c("fepois", "feglm")){

        res = x$cov.unscaled * x$nobs

    } else {
        stop("'bread' is not currently implemented for function ", method, ".")
    }

    res
}






































































