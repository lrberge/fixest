% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{coef.fixest}
\alias{coef.fixest}
\alias{coefficients.fixest}
\title{Extracts the coefficients from a \code{fixest} estimation}
\usage{
\method{coef}{fixest}(object, keep, drop, order, collin = FALSE, agg = TRUE, ...)

\method{coefficients}{fixest}(object, keep, drop, order, collin = FALSE, agg = TRUE, ...)
}
\arguments{
\item{object}{A \code{fixest} object. Obtained using the functions \code{\link{femlm}}, \code{\link{feols}} or \code{\link{feglm}}.}

\item{keep}{Character vector. This element is used to display only a subset of variables. This
should be a vector of regular expressions (see \code{\link[base:regex]{base::regex}} help for more info). Each
variable satisfying any of the regular expressions will be kept. This argument is applied post
aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display
only \code{x1} to \code{x9}, then you could use \code{keep = "x[[:digit:]]$"}. If the first character is an
exclamation mark, the effect is reversed (e.g. keep = "!Intercept" means: every variable that
does not contain \dQuote{Intercept} is kept). See details.}

\item{drop}{Character vector. This element is used if some variables are not to be displayed.
This should be a vector of regular expressions (see \code{\link[base:regex]{base::regex}} help for more info). Each
variable satisfying any of the regular expressions will be discarded. This argument is applied
post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to
display only \code{x1} to \code{x9}, then you could use \verb{drop = "x[[:digit:]]\{2\}}". If the first character
is an exclamation mark, the effect is reversed (e.g. drop = "!Intercept" means: every variable
that does not contain \dQuote{Intercept} is dropped). See details.}

\item{order}{Character vector. This element is used if the user wants the variables to be
ordered in a certain way. This should be a vector of regular expressions (see \code{\link[base:regex]{base::regex}}
help for more info). The variables satisfying the first regular expression will be placed first,
then the order follows the sequence of regular expressions. This argument is applied post
aliasing (see argument \code{dict}). Example: you have the following variables: \code{month1} to \code{month6},
then \code{x1} to \code{x5}, then \code{year1} to \code{year6}. If you want to display first the x's, then the
years, then the months you could use: \code{order = c("x", "year")}. If the first character is an
exclamation mark, the effect is reversed (e.g. order = "!Intercept" means: every variable that
does not contain \dQuote{Intercept} goes first).  See details.}

\item{collin}{Logical, default is \code{FALSE}. Whether the coefficients removed because of collinearity should be also returned as \code{NA}. It cannot be used when coefficients aggregation is also used.}

\item{agg}{Logical scalar, default is \code{TRUE}. If the coefficients of the estimation have been aggregated, whether to report the aggregated coefficients. If \code{FALSE}, the raw coefficients will be returned.}

\item{...}{Not currently used.}
}
\value{
This function returns a named numeric vector.
}
\description{
This function extracts the coefficients obtained from a model estimated with
\code{\link{femlm}}, \code{\link{feols}} or \code{\link{feglm}}.
}
\details{
The coefficients are the ones that have been found to maximize the log-likelihood of the specified model. More information can be found on the models from the estimations help pages: \code{\link{femlm}}, \code{\link{feols}} or \code{\link{feglm}}.

Note that if the model has been estimated with fixed-effects, to obtain the fixed-effect coefficients, you need to use the function \code{\link{fixef.fixest}}.
}
\examples{

# simple estimation on iris data, using "Species" fixed-effects
res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
            Petal.Width | Species, iris)

# the coefficients of the variables:
coef(res)

# the fixed-effects coefficients:
fixef(res)


}
\seealso{
See also the main estimation functions \code{\link{femlm}}, \code{\link{feols}} or \code{\link{feglm}}. \code{\link{summary.fixest}}, \code{\link{confint.fixest}}, \code{\link{vcov.fixest}}, \code{\link{etable}}, \code{\link{fixef.fixest}}.
}
\author{
Laurent Berge
}
