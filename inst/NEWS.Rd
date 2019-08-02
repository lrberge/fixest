\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{fixest}}

\section{Changes in version 0.1.0 (2019-XX-XX)}{
	\subsection{First version, youhou!}{
	    \itemize{
	    \item This package is an effort to create a family of fast and user-friendly functions to perform estimations with multiple fixed-effects (F.E.).

		\item Estimations with fixed-effects (or call it factor variables) is a staple in social science. Hence having a package gathering many methods with fast execution time is of prime importance. At the time of this version, this is the fastest existing method to perform F.E. estimations (often by orders of magnitude, compared to the most efficient alternative methods [both in R and Stata]). The underlying method to obtain the F.E. is based on Berge 2018, and the workhorse of the code is in c++ parallelized via OpenMP (btw thanks Rcpp for simplifying coders' life!).

		\item This package is the follow up of the (now deprecated) package \pkg{FENmlm} which performed fixed-effects estimations but for only four likelihood families. Package \pkg{fixest} completely supersedes \pkg{FENmlm} by extending the method to regular OLS and all GLM families, and adding new utility functions. Further, the design of the functions has been completely overhauled and extended towards much more user-friendliness. Massive effort has been put into providing a set of informative error messages to the user for quick debugging of her workflow (e.g. one of the functions contains over 100 different errors messages).
	    }
	}
}

