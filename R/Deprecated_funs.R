#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Jan 28 10:03:04 2020
# ~: Set of deprecated (near defunct) functions
#----------------------------------------------#

#
# I found a great tutorial for documenting deprecated functions
# https://www.r-bloggers.com/r-deprecate-functions-with-roxygen2-3/
#

# REMOVED FUNCTIONS:
#
# Version 0.8.0 (December 2020):
#   - did_estimate_yearly_effects
#   - did_plot_yearly_effects
#   - errbar
#

## fixest-deprecated.r

#' @title Deprecated functions in package \pkg{fixest}.
#'
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("-deprecated")}.
#'
#' @name fixest-deprecated
#'
#' @keywords internal
NULL




#' @rdname fixest-deprecated
#' @name obs2remove-deprecated
#' @section \code{obs2remove}:
#' This function will be removed in 1 year from 12/11/2020.
NULL

#' Finds observations to be removed from ML estimation with fixed-effects
#'
#' For Poisson, Negative Binomial or Logit estimations with fixed-effects, when the dependent variable is only equal to 0 (or 1 for Logit) for one fixed-effect value this leads to a perfect fit for that fixed-effect value by setting its associated fixed-effect coefficient to \code{-Inf}. Thus these observations need to be removed before estimation. This function gives the observations to be removed. Note that by default the function \code{\link[fixest]{femlm}} or \code{\link[fixest]{feglm}} drops them before performing the estimation.
#'
#' @param fml A formula containing the dependent variable and the fixed-effects. It can be of the type: \code{y ~ fixef_1 + fixef_2} or \code{y ~ x1 | fixef_1 + fixef_1} (in which case variables before the pipe are ignored).
#' @param data A data.frame containing the variables in the formula.
#' @param family Character scalar: either \dQuote{poisson} (default), \dQuote{negbin} or \dQuote{logit}.
#'
#' @return
#' It returns an integer vector of observations to be removed. If no observations are to be removed, an empty integer vector is returned. In both cases, it is of class \code{fixest.obs2remove}.
#' The vector has an attribute \code{fixef} which is a list giving the IDs of the fixed-effects that have been removed, for each fixed-effect dimension.
#'
#' @seealso
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
#'
#' @examples
#'
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' (x = obs2remove(v6 ~ Species, base))
#' attr(x, "fixef")
#'
#' # The two results are identical:
#' res_1 = femlm(v6 ~ Petal.Width | Species, base)
#' # => note + obsRemoved is created
#'
#' res_2 = femlm(v6 ~ Petal.Width | Species, base[-x, ])
#' # => no note because observations are removed before
#'
#' esttable(res_1, res_2)
#'
#' all(res_1$obs_selection$obsRemoved == x)
#'
obs2remove = function(fml, data, family = c("poisson", "negbin", "logit")){
    # in the formula, the fixed-effects must be there:
    # either y ~ fixef_1 + fixef_2
    # either y ~ x1 + x2 | fixef_1 + fixef_2

    counter = getOption("fixest_deprec_obs2remove")
    if(is.null(counter)){
        options("fixest_deprec_obs2remove" = TRUE)
        .Deprecated(msg = "This function is deprecated and will disappear in 1 year from 12/11/2020.")
    }

    #
    # CONTROLS
    #

    # FAMILY

    family = match.arg(family)

    # FML

    if(!"formula" %in% class(fml) || length(fml) != 3){
        stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | fixef_1 + fixef_2' or of the type 'y ~ fixef_1 + fixef_2'.")
    }

    fml_parts = fml_split(fml, raw = TRUE)
    n_rhs = length(fml_parts)

    if(n_rhs > 2){
        stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | fixef_1 + fixef_2' or of the type 'y ~ fixef_1 + fixef_2'.")
    }

    # DATA

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

    # Extracting the variables
    vars_left = all.vars(fml_split(fml, 1, split.lhs = TRUE))
    cluster_fml = fml_split(fml, n_rhs)
    vars_clusters = all.vars(cluster_fml)

    if(length(left_missing <- setdiff(vars_left, dataNames)) > 0){
        stop("Left hand side could not be evaluated, following variables are missing from the data: ", paste0(left_missing, collapse = ", "), ".")
    }

    if(length(right_missing <- setdiff(vars_clusters, dataNames)) > 0){
        stop("The clsuters could not be evaluated, following variables are missing from the data: ", paste0(right_missing, collapse = ", "), ".")
    }

    # Evaluation variables
    lhs = as.vector(eval(fml[[2]], data))
    cluster_mat = model.frame(cluster_fml, data)
    cluster_name = names(cluster_mat)

    #
    # -- CORE --
    #

    Q = length(cluster_name)
    dummyOmises = list()
    obs2remove = c()
    for(q in 1:Q){

        dum_raw = cluster_mat[, q]

        # thisNames = getItems(dum_raw)
        # dum = quickUnclassFactor(dum_raw)
        dum_all = quickUnclassFactor(dum_raw, addItem = TRUE)
        dum = dum_all$x
        thisNames = dum_all$items
        k = length(thisNames)

        # We delete "all zero" outcome
        sum_y_clust = cpp_tapply_vsum(k, lhs, dum)
        n_perClust = cpp_table(k, dum)

        if(family %in% c("poisson", "negbin")){
            qui = which(sum_y_clust == 0)
        } else if(family == "logit"){
            qui = which(sum_y_clust == 0 | sum_y_clust == n_perClust)
        }

        if(length(qui > 0)){
            # We first delete the data:
            dummyOmises[[q]] = thisNames[qui]
            obs2remove = unique(c(obs2remove, which(dum %in% qui)))
        } else {
            dummyOmises[[q]] = character(0)
        }
    }

    names(dummyOmises) = cluster_name

    if(length(obs2remove) == 0){
        print("No observation to be removed.")
        obs2remove = integer(0)
        class(obs2remove) = "fixest.obs2remove"
        return(invisible(obs2remove))
    }

    class(obs2remove) = "fixest.obs2remove"
    attr(obs2remove, "family") = family
    attr(obs2remove, "fixef") = dummyOmises

    return(obs2remove)
}

#' @rdname fixest-deprecated
#' @name summary.fixest.obs2remove-deprecated
#' @section \code{summary.fixest.obs2remove}:
#' This function will be removed in 1 year from 12/11/2020.
NULL

#' Summary method for fixest.obs2remove objects
#'
#' This function synthesizes the information of function \code{\link[fixest]{obs2remove}}. It reports the number of observations to be removed as well as the number of fixed-effects removed per fixed-effect dimension.
#'
#' @method summary fixest.obs2remove
#'
#' @param object A \code{fixest.obs2remove} object obtained from function \code{\link[fixest]{obs2remove}}.
#' @param ... Not currently used.
#'
#' @seealso
#' \code{\link{fixest-deprecated}}
#'
#' @keywords
#' internal
#'
#' @examples
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' x = obs2remove(v6 ~ Species, base)
#' summary(x)
#'
summary.fixest.obs2remove = function(object, ...){

    if(length(object) == 0){
        print("No observation to be removed.")
    } else {
        cat(length(object), " observations removed because of only zero", ifelse(attr(object, "family") == "logit", ", or only one,", ""), " outcomes.\n", sep = "")
        cluster = attr(object, "fixef")
        cat("# fixed-effects removed: ", paste0(names(cluster), ": ", lengths(cluster), collapse = ", "), ".", sep = "")
    }

}







