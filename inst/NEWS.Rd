\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{fixest}}

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

