\name{MapDemeaner}
\alias{MapDemeaner}
\alias{LsmrDemeaner}
\title{Demeaning backend constructors}
\usage{
MapDemeaner(
  fixef_tol = 1e-06,
  fixef_maxiter = 10000,
  extraProj = 0,
  iter_warmup = 15,
  iter_projAfterAcc = 40,
  iter_grandAcc = 4
)

LsmrDemeaner(
  fixef_atol = 1e-06,
  fixef_btol = 1e-06,
  fixef_maxiter = 10000,
  preconditioner = "auto",
  local_size = NULL
)
}
\arguments{
\item{fixef_tol}{Convergence tolerance for the MAP backend.}

\item{fixef_maxiter}{Maximum number of fixed-effect iterations.}

\item{extraProj, iter_warmup, iter_projAfterAcc, iter_grandAcc}{Arguments passed to \code{\link[=demeaning_algo]{demeaning_algo()}} for the MAP backend.}

\item{fixef_atol, fixef_btol}{Absolute tolerances for the LSMR backend. The current \code{withinr} API accepts one tolerance, so fixest passes the larger of the two values.}

\item{preconditioner}{LSMR preconditioner selection. Character values \code{"auto"}, \code{"additive"}, \code{"off"}, and \code{"diagonal"} are supported. Built \code{withinr} preconditioners and \code{withinr::AdditiveSchwarz()} configurations are passed through to \code{withinr}.}

\item{local_size}{Optional local reorthogonalization window size, or \code{NULL}.}
}
\value{
A \code{fixest_demeaner} object accepted by \code{feols()} and \code{demean()}.
}
\description{
\code{MapDemeaner()} selects fixest's native method of alternating projections.
\code{LsmrDemeaner()} selects the optional \code{withinr} LSMR backend.
}
