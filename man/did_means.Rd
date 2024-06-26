% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscfuns.R
\name{did_means}
\alias{did_means}
\title{Treated and control sample descriptives}
\usage{
did_means(
  fml,
  base,
  treat_var,
  post_var,
  tex = FALSE,
  treat_dict,
  dict = getFixest_dict(),
  file,
  replace = FALSE,
  title,
  label,
  raw = FALSE,
  indiv,
  treat_first,
  prepostnames = c("Before", "After"),
  diff.inv = FALSE
)
}
\arguments{
\item{fml}{Either a formula of the type \code{var1 + ... + varN ~ treat} or
\code{var1 + ... + varN ~ treat | post}. Either a data.frame/matrix containing
all the variables for which the means are to be computed (they must be numeric of course).
Both the treatment and the post variables must contain only exactly two values.
You can use a point to select all the variables of the data set: \code{. ~ treat}.}

\item{base}{A data base containing all the variables in the formula \code{fml}.}

\item{treat_var}{Only if argument \code{fml} is \emph{not} a formula. The vector identifying
the treated and the control observations (the vector can be of any type but must
contain only two possible values). Must be of the same length as the data.}

\item{post_var}{Only if argument \code{fml} is \emph{not} a formula. The vector identifying
the periods (pre/post) of the observations (the vector can be of any type but
must contain only two possible values). The first value (in the sorted sense)
of the vector is taken as the pre period. Must be of the same length as the data.}

\item{tex}{Should the result be displayed in Latex? Default is \code{FALSE}. Automatically
set to \code{TRUE} if the table is to be saved in a file using the argument \code{file}.}

\item{treat_dict}{A character vector of length two. What are the names of the treated
and the control? This should be a dictionary: e.g. \code{c("1"="Treated", "0" = "Control")}.}

\item{dict}{A named character vector. A dictionary between the variables names and an alias.
For instance \code{dict=c("x"="Inflation Rate")} would replace the variable name \code{x}
by \dQuote{Inflation Rate}.}

\item{file}{A file path. If given, the table is written in Latex into this file.}

\item{replace}{Default is \code{TRUE}, which means that when the table is exported, the
existing file is not erased.}

\item{title}{Character string giving the Latex title of the table. (Only if exported.)}

\item{label}{Character string giving the Latex label of the table. (Only if exported.)}

\item{raw}{Logical, default is \code{FALSE}. If \code{TRUE}, it returns the information without formatting.}

\item{indiv}{Either the variable name of individual identifiers, a one sided formula,
or a vector. If the data is that of a panel, this can be used to track the number
of individuals per group.}

\item{treat_first}{Which value of the 'treatment' vector should appear on the left?
By default the max value appears first (e.g. if the treatment variable is a 0/1 vector,
1 appears first).}

\item{prepostnames}{Only if there is a 'post' variable. The names of the pre and post
periods to be displayed in Latex. Default is \code{c("Before", "After")}.}

\item{diff.inv}{Logical, default to \code{FALSE}. Whether to inverse the difference.}
}
\value{
It returns a data.frame or a Latex table with the conditional means and statistical differences between the groups.
}
\description{
This function shows the means and standard-deviations of several variables conditional on whether they are from the treated or the control group. The groups can further be split according to a pre/post variable. Results can be seamlessly be exported to Latex.
}
\details{
By default, when the user tries to apply this function to nun-numeric variables, an error is raised. The exception is when the all variables are selected with the dot (like in \code{. ~ treat}. In this case, non-numeric variables are automatically omitted (with a message).

NAs are removed automatically: if the data contains NAs an information message will be prompted. First all observations containing NAs relating to the treatment or post variables are removed. Then if there are still NAs for the variables, they are excluded separately for each variable, and a new message detailing the NA breakup is prompted.
}
\examples{

# Playing around with the DiD data
data(base_did)

# means of treat/control
did_means(y+x1+period~treat, base_did)

# same but inverting the difference
did_means(y+x1+period~treat, base_did, diff.inv = TRUE)

# now treat/control, before/after
did_means(y+x1+period~treat|post, base_did)

# same but with a new line giving the number of unique "indiv" for each case
did_means(y+x1+period~treat|post, base_did, indiv = "id")

# same but with the treat case "0" coming first
did_means(y+x1+period~treat|post, base_did, indiv = ~id, treat_first = 0)

# Selecting all the variables with "."
did_means(.~treat|post, base_did, indiv = "id")


}
