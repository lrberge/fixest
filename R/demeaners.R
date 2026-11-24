#' Demeaning backend constructors
#'
#' `MapDemeaner()` selects fixest's native method of alternating projections.
#' `LsmrDemeaner()` selects the optional `withinr` LSMR backend.
#'
#' @param fixef_tol Convergence tolerance for the MAP backend.
#' @param fixef_maxiter Maximum number of fixed-effect iterations.
#' @param extraProj,iter_warmup,iter_projAfterAcc,iter_grandAcc Arguments passed
#'   to [demeaning_algo()] for the MAP backend.
#' @param fixef_atol,fixef_btol Absolute tolerances for the LSMR backend. The
#'   current `withinr` API accepts one tolerance, so fixest passes the larger
#'   of the two values.
#' @param preconditioner LSMR preconditioner selection. Character values
#'   `"auto"`, `"additive"`, `"off"`, and `"diagonal"` are supported. Built
#'   `withinr` preconditioners and `withinr::AdditiveSchwarz()` configurations
#'   are passed through to `withinr`.
#' @param local_size Optional local reorthogonalization window size, or `NULL`.
#'
#' @return A `fixest_demeaner` object accepted by `feols()` and `demean()`.
#'
#' @export
MapDemeaner = function(fixef_tol = 1e-6, fixef_maxiter = 10000,
                       extraProj = 0, iter_warmup = 15,
                       iter_projAfterAcc = 40, iter_grandAcc = 4){

  check_arg(fixef_tol, "numeric scalar GT{(10000*.Machine$double.eps)}")
  check_arg(fixef_maxiter, "integer scalar GE{1}")
  check_arg("integer scalar", extraProj, iter_warmup, iter_projAfterAcc, iter_grandAcc)

  res = list(
    fixef_tol = fixef_tol,
    fixef_maxiter = as.integer(fixef_maxiter),
    fixef.algo = demeaning_algo(
      extraProj = extraProj,
      iter_warmup = iter_warmup,
      iter_projAfterAcc = iter_projAfterAcc,
      iter_grandAcc = iter_grandAcc,
      internal = TRUE
    )
  )

  class(res) = c("fixest_demeaner_map", "fixest_demeaner")
  res
}

#' @rdname MapDemeaner
#' @export
LsmrDemeaner = function(fixef_atol = 1e-6, fixef_btol = 1e-6,
                        fixef_maxiter = 10000, preconditioner = "auto",
                        local_size = NULL){

  check_arg(fixef_atol, fixef_btol, "numeric scalar GT{0}")
  check_arg(fixef_maxiter, "integer scalar GE{1}")

  if(!is.null(local_size)){
    check_arg(local_size, "integer scalar GE{1}")
    local_size = as.integer(local_size)
  }

  preconditioner = fixest_check_lsmr_preconditioner(preconditioner)

  res = list(
    fixef_atol = fixef_atol,
    fixef_btol = fixef_btol,
    fixef_maxiter = as.integer(fixef_maxiter),
    preconditioner = preconditioner,
    local_size = local_size
  )

  class(res) = c("fixest_demeaner_lsmr", "fixest_demeaner")
  res
}

fixest_check_lsmr_preconditioner = function(preconditioner){

  if(is.null(preconditioner)){
    return(NULL)
  }

  if(is.character(preconditioner) && length(preconditioner) == 1L && !is.na(preconditioner)){
    return(preconditioner)
  }

  if(inherits(preconditioner, "within_preconditioner") ||
     inherits(preconditioner, "within_additive_schwarz")){
    return(preconditioner)
  }

  stop("Argument `preconditioner` must be a character scalar, NULL, ",
       "a withinr preconditioner, or a withinr AdditiveSchwarz configuration.")
}

fixest_check_demeaner = function(demeaner){

  if(is.null(demeaner)){
    return(invisible(NULL))
  }

  if(!inherits(demeaner, "fixest_demeaner")){
    stop("Argument `demeaner` must be NULL or created by MapDemeaner() ",
         "or LsmrDemeaner().")
  }

  invisible(NULL)
}

fixest_resolve_demeaner = function(demeaner, fixef.tol, fixef.iter, fixef.algo,
                                   explicit_demeaning_args = character()){

  if(!is.null(demeaner) && length(explicit_demeaning_args) > 0){
    pblm = paste0("`", explicit_demeaning_args, "`", collapse = ", ")
    stop("When `demeaner` is provided, do not also provide ",
         pblm, ". ",
         "Put these controls inside MapDemeaner() or LsmrDemeaner().")
  }

  if(is.null(fixef.algo)){
    fixef.algo = demeaning_algo(internal = TRUE)
  }

  if(is.null(demeaner)){
    demeaner = MapDemeaner(
      fixef_tol = fixef.tol,
      fixef_maxiter = fixef.iter,
      extraProj = fixef.algo$extraProj,
      iter_warmup = fixef.algo$iter_warmup,
      iter_projAfterAcc = fixef.algo$iter_projAfterAcc,
      iter_grandAcc = fixef.algo$iter_grandAcc
    )
  } else {
    fixest_check_demeaner(demeaner)
  }

  demeaner
}

fixest_demeaner_tol = function(demeaner){

  if(inherits(demeaner, "fixest_demeaner_lsmr")){
    return(max(demeaner$fixef_atol, demeaner$fixef_btol))
  }

  demeaner$fixef_tol
}

fixest_demeaner_iter = function(demeaner){
  demeaner$fixef_maxiter
}

fixest_demeaner_algo = function(demeaner){

  if(inherits(demeaner, "fixest_demeaner_map")){
    return(demeaner$fixef.algo)
  }

  demeaning_algo(internal = TRUE)
}

fixest_lsmr_preconditioner_arg = function(demeaner){

  preconditioner = demeaner$preconditioner

  if(is.character(preconditioner)){
    preconditioner_lower = tolower(preconditioner)

    if(preconditioner_lower == "auto"){
      return("additive")
    }

    if(preconditioner_lower %in% c("additive", "off", "diagonal")){
      return(preconditioner_lower)
    }

    warning("Unsupported LSMR preconditioner \"", preconditioner,
            "\"; using the withinr default.", call. = FALSE)
    return(NULL)
  }

  preconditioner
}

fixest_withinr_supports_lsmr_options = function(){
  requireNamespace("withinr", quietly = TRUE) &&
    exists("LsmrOptions", envir = asNamespace("withinr"), inherits = FALSE)
}

#' @export
print.fixest_demeaner_map = function(x, ...){
  cat("<fixest_demeaner: MAP>\n")
  cat("  fixef_tol: ", x$fixef_tol, "\n", sep = "")
  cat("  fixef_maxiter: ", x$fixef_maxiter, "\n", sep = "")
  invisible(x)
}

#' @export
print.fixest_demeaner_lsmr = function(x, ...){
  preconditioner = x$preconditioner
  if(is.character(preconditioner)){
    preconditioner = paste0('"', preconditioner, '"')
  } else if(is.null(preconditioner)){
    preconditioner = "NULL"
  } else {
    preconditioner = class(preconditioner)[1]
  }

  cat("<fixest_demeaner: LSMR>\n")
  cat("  fixef_atol: ", x$fixef_atol, "\n", sep = "")
  cat("  fixef_btol: ", x$fixef_btol, "\n", sep = "")
  cat("  fixef_maxiter: ", x$fixef_maxiter, "\n", sep = "")
  cat("  preconditioner: ", preconditioner, "\n", sep = "")
  cat("  local_size: ", if(is.null(x$local_size)) "NULL" else x$local_size, "\n", sep = "")
  invisible(x)
}
