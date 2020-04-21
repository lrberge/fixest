\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{fixest}}

\section{Changes in version 0.xx.xx (2020-xx-xx)}{

    \subsection{Bug fixes}{
        \itemize{
            \item Bug with estimations with varying slopes if the fixed-effect relative to the slope is not in its decreasing order (thanks to Davide Proserpio).
            \item Bug when interacting two variables with the \code{var::fe} syntax with \code{confirm = TRUE} and no reference.
            \item Bug in \code{etable} when the standard-errors where \code{NA}.
        }
    }

    \subsection{User-visible changes}{
        \itemize{
            \item Interactions in \code{var::fe(ref)} now accept multiple references (i.e. \code{ref} can be a vector).
            \item In \code{etable}, the variable names of non-Latex output can now be changed.
        }
    }

    \subsection{Error-handling}{
        \itemize{
            \item Improved error-handling with dreamerr's functions.
        }
    }
}

\section{Changes in version 0.4.1 (2020-04-13)}{
  \subsection{Bug fixes}{
  \itemize{
    \item Major bug leading R to crash when using non-linear-in-parameters right-hand-sides in feNmlm. Only occured when some observations were removed from the data set (due to NAness or to perfect fit). [Thanks to @marissachilds, GH issue #17.]
    \item In the \code{collinearity} help pages: an example could lead to an error (due to random data generation). It has been removed.
    \item In \code{collinearity}, corrected the problem of display of the intercept in some situations.
    \item Defaults for the arguments \code{cex} and \code{lwd} in \code{coefplot} have been changed to 1 and 1 (instead of par("cex") and par("lwd")). Otherwise this led to the creation of \code{Rplots.pdf} in the working directory (thanks to Kurt Hornik).
    \item Corrected a typo in the article's title in the vignette.
    }
  }

  \subsection{Help}{
  \itemize{
    \item Rewriting of sections, correction of small mistakes (wrong argument names), dropping completely the 'cluster' terminology (meant for fixed-effects), addition of what is contained in fixest objects.
    }
  }

  \subsection{Other}{
  \itemize{
    \item Adding a README.md.
    \item Small corrections in the vignette.
    }
  }
}

\section{Changes in version 0.4.0 (2020-03-27)}{

    \subsection{User visible changes: Latex export}{
        \itemize{
            \item Better Latex special character escaping.
            \item New argument \code{fixef_sizes.simplify}, which provides the sizes of the fixed-effects in parentheses when there is no ambiguity.
            \item You can suppress the line with the significance codes with \code{signifCode = NA}.
            \item New argument \code{float} which decides whether to embed the table into a table environment. By default it is set to \code{TRUE} if a \code{title} or \code{label} is present.
            \item New argument \code{keep} to select the variables to keep in the table.
            \item New way to keep/drop/order variables with the special argument "\%". If you use "\%var", then it makes reference to the original variable name, not the aliased one (which is the default).
            \item Better rendering of significant digits.
        }
    }

    \subsection{User visible changes: coefplot}{
        \itemize{
            \item Argument \code{horiz}. The coefficients can now be displayed horizontally instead of vertically.
            \item The coefficient labels, when in the x-axis, can now be displayed in three different ways thanks to the new argument \code{lab.fit}: "simple", the classic axis, "multi", the labels appear across multiple lines to avoid collision, and "tilted" for tilted labels.
            \item The margins now automatically fit.
            \item Argument \code{style} allows you to set styles with the function \code{setFixest_coefplot}, you can then summon the style in \code{coefplot} with this argument.
            \item Use the ampersand to set dictionnary variables specific to \code{coefplot} (e.g. setFixest_dict(c(var1 = "a $\\times$ b", "&var1" = "&a\%*\%b")), the variable \code{var1} will be rendered differently in \code{etable} and in \code{coefplot}).
            \item Better display of groups (with the arguments \code{group} and \code{group.par}).
        }
    }

    \subsection{New methods}{
        \itemize{
            \item \code{terms.fixest} giving the terms of the estimation.
        }
    }

    \subsection{Other}{
        \itemize{
            \item All \code{donttest} sections were removed from help pages.
        }
    }

}

\section{Changes in version 0.3.1 (2020-02-09)}{

    \subsection{Major bug fix}{
        \itemize{
            \item[panel] Fixed faulty memory access when taking the lead of a variable.
        }
    }

    \subsection{Other bug fixes}{
        \itemize{
            \item[esttable/esttex] These two functions were replaced by the function \code{etable}. In the process, some of their arguments were "lost", this is now corrected.
            \item[etable] Better escaping of special characters.
            \item[estimations] Bug when particular non-numeric vectors were used in explanatory variables.
        }
    }

    \subsection{New features}{
        \itemize{
            \item[coefplot] The function \code{coefplot} now accepts lists of estimations.
        }
    }

}

\section{Changes in version 0.3.0 (2020-02-01)}{

    \subsection{New feature: Lagging}{
         \itemize{
            \item You can now add lags and leads in any \code{fixest} estimations. You only need to provide the panel identifiers with the new argument \code{panel.id}, then you're free to use the new functions \code{l()} for lags and \code{f()} for leads.
            \item You can also set up a panel data set using the function \code{panel} which allows you to use the lagging functions without having to provide the argument \code{panel.id}, and which dispose of more options for setting the panel.
        }
    }

    \subsection{New feature: Interactions}{
         \itemize{
            \item You can now add interactions in formulas with a new syntax: \code{var::fe(ref)}
            \item The command \code{var::fe(ref)} interacts the variable \code{var} with each value of \code{fe} and sets \code{ref} as a reference. Note that if you don't use the argument \code{ref}, the command \code{var::fe} is identical to \code{var:factor(fe)}.
            \item Using \code{var::fe(ref)} to write interactions opens up a special treatment of such variables in the exporting function \code{etable} and in the coefficient plotting function \code{coefplot}.
        }
    }

    \subsection{New feature: \code{coefplot}}{
         \itemize{
            \item You can plot coefficients and their associated confidence intervals with the function \code{coefplot}.
            \item \code{coefplot} dispose of many options, whose default values can be set with the function \code{setFixest_coefplot}.
            \item As for the function \code{etable}, you can easily rename/drop/order the coefficients.
            \item \code{coefplot} detects when interactions have been used and offers a special display for it.
        }
    }

    \subsection{New functions}{
         \itemize{
            \item[etable] Estimations table: new function to export the results of multiple estimations. Replaces the two functions \code{esttex} and \code{esttable} (the two functions still exist but they will be deprecated in the future).
            \item[Lagging] New functions related to lagging: \code{l}, \code{f}, \code{panel}, \code{unpanel} and \code{[.fixest_panel}.
            \item[Utilities] A set of small utility functions has been added. They allow to extract part a coefficient table or parts of it (like the t-statistics of the standard-error) from an estimation. These functions are \code{coeftable}, \code{ctable} (an alias to \code{coeftable}), \code{se}, \code{tstat} and \code{pvalue}.
            \item[coefplot] The functions \code{coefplot} and \code{setFixest_coefplot}.
            \item[dof] New function to set the type of degree of freedom adjustment when computing the variance-covariance matrix. You can permanently set the type of DoF adjustment with the new function setFixest_dof().
        }
    }

    \subsection{User visible changes}{
         \itemize{
            \item[all estimations] A key pre-processing step has been paralellized => algorithm faster in general and much faster for multi-FEs.
            \item[predict & fitted] Predict and fitted now returns vectors of the length equal to the one of original data.
            \item[standard-errors] New ways to compute the standard-errors have been implemented. In particular, now it account for the "nestedness" of the fixed-effects in the clusters by default. You can freely change how to compute the degrees of freedom correction with the function dof().
            \item[r2] Computation of the within-R2 for feglm models is now self-contained.
            \item[all estimations] New, more accurate, stopping criterion for 2+ fixed-effects.
            \item[feols] Estimations are slightly faster.
            \item[etable/esttex] When there are interactions, R may change the order of the interactions, making two interactions in two different estimations look different while they are in fact the same (e.g. x3:x2 and x2:x3). Now esstable automatically reorders the interactions when needed for comparison across estimations.
            \item[etable/esttable] The type of standard errors is now always shown.
            \item[etable/esttex] The aliases provided by 'dict' are also applied within interactions. For example: \code{dict=c(x1="Wind", x2="Rain")}, with an estimation with the following variables 'x1', 'x2', 'x1:x2' will lead to the following aliases in Latex 'Wind', 'Rain' and 'Wind $\\times$ Rain'.
            \item[etable/esttex] Interactions of similar values but of different order (e.g. x1:x2 and x2:x1) are reorderd to appear in the same lines.
            \item[etable/esttex] The i) type of standard errors and ii) the significance codes, are now displayed in two separate lines (otherwise the line would be too wide).
            \item[etable/esttex] Argument \code{yesNoFixef} can be of length one, defaulting the second element to the empty string.
            \item[etable/esttex] Escaping of Latex special characters is now much more robust.
        }
    }

    \subsection{Bug correction}{
         \itemize{
            \item[all estimations] Fixed: bug when functions in the formula returned matrices.
            \item[update] Fixed: error message when the data is missing.
            \item[feglm] Fixed: bug double estimation when family not equal to poisson or logit
            \item[feglm] Fixed: severe bug occurring for families not equal to poisson or logit
            \item[predict] Fixed: bug when the estimation contained combined FEs.
            \item[summary] Regarding small sample only: now Student t distribution is used instead of the Normal to compute the pvalue.
            \item[esttex] Different variables with the same aliases (given by the argument 'dict') now appear in the same row.
            \item[esttex] Arguments 'drop' and 'order' are now applied post aliasing (alias given by the argument 'dict').
            \item[esttex] But when exporting multi-way standard errors.
            \item[r2] Small bug regarding objects obtained with \code{did_estimate_yearly_effects}.
            \item[estimations] bug when using weights in \code{feglm}.
        }
    }
}

\section{Changes in version 0.2.1 (2019-11-22)}{

    \subsection{Major bug correction}{
        \itemize{
            \item lag.formula: Bug introduced from previous update which could lead to wrong results. Now fixed.
        }
    }


    \subsection{Major user visible changes}{
         \itemize{
            \item[All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are string vectors.
        }
    }



}

\section{Changes in version 0.2.0 (2019-11-19)}{

    \subsection{New function}{
         \itemize{
            \item[did_means] New function \code{did_means} to conveniently compare means of groups of observations (both treat/control and pre/post). Contains tools to easily export in Latex.
        }
    }

    \subsection{Major user visible changes}{
        \itemize{
            \item[All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are of type integer or double.
            \item[esttex, esttable] New argument 'fitstat' to select which fit statistic to display. The default adapts to the models. Old arguments (loglik, bic, aic, sq.cor) are dropped.
            \item[esttable] Significantly better rendering of SE types.
            \item[r2] Now NA is returned for R2s that have no theoretical justification (e.g. within R2 when no FEs, or 'regular' R2 for ML models).
        }
    }

    \subsection{Minor user visible changes}{
        \itemize{
            \item[did_plot_yearly_effects] Now the name of the dependent variable appears on the y-axis.
            \item[esttex] Usage of the \code{sym} macro in Latex is dropped.
        }
    }

    \subsection{Bug correction}{
        \itemize{
            \item[fixef.fixest] bug could appear when using varying slopes coefficients in some specific circumstances (when the slope FEs were different from the regular FEs).
            \item[fixef.fixest] bug when many regular FEs jointly with varying slopes.
            \item[fixef.fixest] regarding slope coefficients: now the algorithm also evaluates functions of variables.
            \item[esttable] Width of the "separating lines" now appropriately set for long dependent variable names.
            \item[esttex] Spelling mistake corrected.
            \item[estimations] Bug could occur for extremely small data sets (< 10 observations).
        }
    }

    \subsection{Error handling}{
        \itemize{
            \item[esttex, esttable] More informative error messages in functions esttex and esttable.
        }
    }
}

\section{Changes in version 0.1.2 (2019-10-04)}{
    \subsection{Major bug correction}{
        \itemize{
            \item lag.formula: When the data was not in a particular format, the results could be wrong. Now corrected.
        }
    }
}

\section{Changes in version 0.1.1 (2019-09-20)}{
    \subsection{Major bug correction}{
	    \itemize{
    	    \item feglm: bug when a) the deviance at initialization was higher than the deviance of the first iteration of the IRWLS and b) the step-halving was unable to find a lower deviance. This led the estimation to fail with an error although it should have been performed properly.
    	    \item did_estimate_yearly_effects: bug when the estimation involved periods with negative values
	    }
    }

    \subsection{Minor bug correction}{
        \itemize{
            \item esttex: bug regarding the number of digits of negative coefficients to be displayed
            \item esttex: now properly escaping the percentage and the underscore for exports in Latex
            \item esttex: bug when changing the names of the dependent variables using a dictionnary
            \item vcov: some warning messages were misleading
            \item update: bug update when using the argument nframes
            \item update: bug when updating the function fepois
        }
    }

    \subsection{Error handling}{
        \itemize{
            \item Better error messages for: did_estimate_yearly_effects, main estimation functions, setFixest_dict, fepois and fenegbin.
        }
    }
}

\section{Version 0.1.0 (2019-09-03)}{
	\subsection{First version}{
	    \itemize{
	    \item This package is an effort to create a family of fast and user-friendly functions to perform estimations with multiple fixed-effects (F.E.).

		\item Estimations with fixed-effects (or call it factor variables) is a staple in social science. Hence having a package gathering many methods with fast execution time is of prime importance. At the time of this version, this is the fastest existing method to perform F.E. estimations (often by orders of magnitude, compared to the most efficient alternative methods [both in R and Stata]). The underlying method to obtain the F.E. is based on Berge 2018, and the workhorse of the code is in c++ parallelized via OpenMP (btw thanks Rcpp for simplifying coders' life!).

		\item This package is the follow up of the (now deprecated) package \pkg{FENmlm} which performed fixed-effects estimations but for only four likelihood families. Package \pkg{fixest} completely supersedes \pkg{FENmlm} by extending the method to regular OLS and all GLM families, and adding new utility functions. Further, the design of the functions has been completely overhauled and extended towards much more user-friendliness. Massive effort has been put into providing a set of informative error messages to the user for quick debugging of her workflow (e.g. one of the functions contains over 100 different errors messages).
	    }
	}
}

