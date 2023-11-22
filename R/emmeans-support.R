##### emmeans support
##### contributed by Russ Lenth russell-lenth[at]uiowa.edu

### CAUTION: For models with fixed effects, this code will yield the correct estimates
### of estimated marginal means -- but the wrong SEs because the intercept is among
### the fixed effects, hence not included in vcov(). CONTRASTS among EMMs will
### work fine because those give zero weight to the intercept.

### Another note: Only main model predictors are included in the reference grid.
### Not any fixed effects, not any instrumental variables

#' Support for emmeans package
#'
#' If \pkg{emmeans} is installed, its functionality is supported for \code{fixest}
#' or \code{fixest_multi} objects. Its reference grid is based on the main part
#' of the model, and does not include fixed effects or instrumental variables.
#' Note that any desired arguments to \code{vcov()} may be passed as optional
#' arguments in \code{emmeans::emmeans()} or \code{emmeans::ref_grid()}.
#'
#' @note
#' When fixed effects are present, estimated marginal means (EMMs) are estimated
#' correctly, provided equal weighting is used. However, the SEs of these EMMs
#' will be incorrect - often dramatically - because the estimated variance of
#' the intercept is not available. However, \emph{contrasts} among EMMs can be
#' estimated and tested with no issues, because these do not involve the
#' intercept.
#'
#' @examples
#' if(requireNamespace("emmeans") && requireNamespace(AER)) {
#'     data(Fatalities, package = "AER")
#'     Fatalities$frate <- with(Fatalities, fatal/pop * 10000)
#'     fat.mod = feols(frate ~ breath * jail * beertax | state + year, data = Fatalities)
#'     emm = emmeans(fat.mod, ~ breath*jail, cluster = ~ state + year)
#'     emm   ### SEs and CIs are incorrect
#'
#'     contrast(emm, "consec", by = "breath")   ### results are reliable
#' }
#' @name emmeans_support
NULL


### DO NOT export - this is done dynamically in .onLoad()
recover_data.fixest = function (object, ...) {
    fcall = object$call
    dat = recover_data(fcall, delete.response(terms(object)), object$na.action,
                       pwts = weights(object), ...)
    if(!is.character(dat)) { # i.e., we actually did recover data
        # if any fixed effects, make them a prior offset
        fe = try(fixef(object), silent = TRUE)
        if(!inherits(fe, "try-error")) {
            dat$.static.offset. = sum(sapply(fe, mean))
            attr(dat, "predictors") = c(attr(dat, "predictors"), ".static.offset.")
        }
    }
    dat
}

emm_basis.fixest = function(object, trms, xlev, grid, ...) {
    bhat = coef(object)
    nm = if(is.null(names(bhat))) row.names(bhat) else names(bhat)
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    assign = attr(X, "assign")
    X = X[, nm, drop = FALSE]
    bhat = as.numeric(bhat)
    V = vcov(object, ...) # we're using fixest's method as it already allows vcov option

    # fixest handles rank deficiencies by tossing columns. Plus it seems hard to
    # recover the needed portion of X after the fixed effects are orthogonalized out.
    # So at this point we are just ignoring non-estimability issues (and generally
    # bhat will NOT be aligned with linfct, creating mysterious errors downstream)
    nbasis = estimability::all.estble

    misc = list()
    if (!is.null(fam<- object$family)) {
        if(is.character(fam))
            fam = list(family = fam, link = "identity")
        misc = .std.link.labels(fam, misc)
        dfargs = list(df = Inf)
    }
    else
        dfargs = list(df = object$nobs - object$nparams)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc,
         model.matrix = .cmpMM(object$qr, assign = assign))
}

# multi methods ------------------------------------------------
recover_data.fixest_multi = function(object, which, ...) {
    if(missing(which))
        return("For 'fixest_multi' objects, you must use 'which' to specify which response to use.")
    recover_data.fixest(object[[which]], ...)
}

emm_basis.fixest_multi = function(object, trms, xlev, grid, which, ...)
    emm_basis.fixest(object[[which]], trms, xlev, grid, ...)

