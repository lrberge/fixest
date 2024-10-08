% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_tools.R
\name{n_unik}
\alias{n_unik}
\alias{print.vec_n_unik}
\alias{print.list_n_unik}
\title{Prints the number of unique elements in a data set}
\usage{
n_unik(x)

\method{print}{vec_n_unik}(x, ...)

\method{print}{list_n_unik}(x, ...)
}
\arguments{
\item{x}{A formula, with data set names on the LHS and variables on the RHS,
like \code{data1 + data2 ~ var1 + var2}. The following special variables are
admitted: \code{"."} to get default values, \code{".N"} for the number of observations, \code{".U"}
for the number of unique rows, \code{".NA"} for the number of rows with at least one NA.
Variables can be combined with \code{"^"}, e.g. \code{df~id^period}; use \code{id\%^\%period} to also
include the terms on both sides. Note that using \code{:} and \code{*}
is equivalent to \code{^} and \verb{\%^\%}. Sub select with \code{id[cond]}, when doing so \code{id}
is automatically included. Conditions can be chained, as in \code{id[cond1, cond2]}.
Use \code{NA(x, y)} in conditions instead of \code{is.na(x) | is.na(y)}. Use the \verb{!!}
operator to have both a condition and its opposite. To compare the keys
in two data sets, use \code{data1:data2}. If not a formula, \code{x} can be: a vector
(displays the # of unique values); a \code{data.frame} (default values are displayed),
or a "sum" of data sets like in \code{x = data1 + data2}, in that case it is equivalent
to \code{data1 + data2 ~ .}.}

\item{...}{Not currently used.}
}
\value{
It returns a vector containing the number of unique values per element. If several
data sets were provided, a list is returned, as long as the number of data sets,
each element being a vector of unique values.
}
\description{
This utility tool displays the number of unique elements in one or multiple data.frames
as well as their number of NA values.
}
\section{Special values and functions}{


In the formula, you can use the following special values: \code{"."}, \code{".N"}, \code{".U"}, and \code{".NA"}.

\describe{

\item{\code{"."}}{Accesses the default values. If there is only one data set and the
data set is \emph{not} a \code{data.table}, then the default is to display the number of
observations and the number of unique rows. If the data is a \code{data.table}, the number
of unique items in the key(s) is displayed instead of the number of unique rows
(if the table has keys of course). If there are two or more data sets, then the
default is to display the unique items for: a) the variables common across all data sets,
if there's less than 4, and b) if no variable is shown in a), the number of variables
common across at least two data sets, provided there are less than 5. If the data sets are
data tables, the keys are also displayed on top of the common variables. In any case, the
number of observations is always displayed.}

\item{\code{".N"}}{Displays the number of observations.}

\item{\code{".U"}}{Displays the number of unique rows.}

\item{\code{".NA"}}{Displays the number of rows with at least one NA.}

}
}

\section{The \code{NA} function}{


The special function \code{NA} is an equivalent to \code{is.na} but can handle several variables.
For instance, \code{NA(x, y)} is equivalent to \code{is.na(x) | is.na(y)}. You can add as
many variables as you want as arguments. If no argument is provided, as in \code{NA()},
it is identical to having all the variables of the data set as argument.
}

\section{Combining variables}{


Use the "hat", \code{"^"}, operator to combine several variables. For example \code{id^period}
will display the number of unique values of id x period combinations.

Use the "super hat", \code{"\%^\%"}, operator to also include the terms on both sides.
For example, instead of writing \code{id + period + id^period}, you can simply write \code{id\%^\%period}.

Alternatively, you can use \code{:} for \code{^} and \code{*} for \verb{\%^\%}.
}

\section{Sub-selections}{


To show the number of unique values for sub samples, simply use \verb{[]}.
For example, \code{id[x > 10]} will display the number of unique \code{id} for which \code{x > 10}.

Simple square brackets lead to the inclusion of both the variable and its subset.
For example \code{id[x > 10]} is equivalent to \code{id + id[x > 10]}.
To include only the sub selection, use double square brackets, as in \code{id[[x > 10]]}.

You can add multiple sub selections at once, only separate them with a comma.
For example \code{id[x > 10, NA(y)]} is equivalent to \code{id[x > 10] + id[NA(y)]}.

Use the double negative operator, i.e. \verb{!!}, to include both a condition and
its opposite at once. For example \code{id[!!x > 10]} is equivalent to \code{id[x > 10, !x > 10]}.
Double negative operators can be chained, like in \code{id[!!cond1 & !!cond2]}, then the
cardinal product of all double negatived conditions is returned.
}

\examples{

data = base_did
data$x1.L1 = round(lag(x1~id+period, 1, data))

# By default, just the formatted number of observations
n_unik(data)

# Or the nber of unique elements of a vector
n_unik(data$id)

# number of unique id values and id x period pairs
n_unik(data ~.N + id + id^period)

# use the \%^\% operator to include the terms on the two sides at once
# => same as id*period
n_unik(data ~.N + id \%^\% period)

# using sub selection with []
n_unik(data ~.N + period[!NA(x1.L1)])

# to show only the sub selection: [[]]
n_unik(data ~.N + period[[!NA(x1.L1)]])

# you can have multiple values in [],
# just separate them with a comma
n_unik(data ~.N + period[!NA(x1.L1), x1 > 7])

# to have both a condition and its opposite,
# use the !! operator
n_unik(data ~.N[!!NA(x1.L1)])

# the !! operator works within condition chains
n_unik(data ~.N[!!NA(x1.L1) & !!x1 > 7])

# Conditions can be distributed
n_unik(data ~ (id + period)[x1 > 7])

#
# Several data sets
#

# Typical use case: merging
# Let's create two data sets and merge them

data(base_did)
base_main = base_did
base_extra = sample_df(base_main[, c("id", "period")], 100)
base_extra$id[1:10] = 111:120
base_extra$period[11:20] = 11:20
base_extra$z = rnorm(100)

# You can use db1:db2 to compare the common keys in two data sets
 n_unik(base_main:base_extra)

tmp = merge(base_main, base_extra, all.x = TRUE, by = c("id", "period"))

# You can show unique values for any variable, as before
n_unik(tmp + base_main + base_extra ~ id[!!NA(z)] + id^period)



}
\author{
Laurent Berge
}
