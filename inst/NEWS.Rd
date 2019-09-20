\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{fixest}}

\section{Changes in version 0.1.1 (20 September 2019)}{
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

\section{Version 0.1.0 (3 September 2019)}{
	\subsection{First version}{
	    \itemize{
	    \item This package is an effort to create a family of fast and user-friendly functions to perform estimations with multiple fixed-effects (F.E.).

		\item Estimations with fixed-effects (or call it factor variables) is a staple in social science. Hence having a package gathering many methods with fast execution time is of prime importance. At the time of this version, this is the fastest existing method to perform F.E. estimations (often by orders of magnitude, compared to the most efficient alternative methods [both in R and Stata]). The underlying method to obtain the F.E. is based on Berge 2018, and the workhorse of the code is in c++ parallelized via OpenMP (btw thanks Rcpp for simplifying coders' life!).

		\item This package is the follow up of the (now deprecated) package \pkg{FENmlm} which performed fixed-effects estimations but for only four likelihood families. Package \pkg{fixest} completely supersedes \pkg{FENmlm} by extending the method to regular OLS and all GLM families, and adding new utility functions. Further, the design of the functions has been completely overhauled and extended towards much more user-friendliness. Massive effort has been put into providing a set of informative error messages to the user for quick debugging of her workflow (e.g. one of the functions contains over 100 different errors messages).
	    }
	}
}

