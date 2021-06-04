#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Thu Jul 09 09:52:31 2020
# ~: etable
#----------------------------------------------#


# GEN ALIASES ####
# Don't forget to run gen_etable_aliases before cran submissions
# gen_etable_aliases()


#' Estimations table (export the results of multiples estimations to a DF or to Latex)
#'
#' Aggregates the results of multiple estimations and displays them in the form of either a Latex table or a \code{data.frame}.
#'
#' @inheritParams summary.fixest
#'
#' @param ... Used to capture different \code{fixest} estimation objects (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). Note that any other type of element is discarded. Note that you can give a list of \code{fixest} objects.
#' @param digits Integer or character scalar. Default is 4 and represents the number of significant digits to be displayed for the coefficients and standard-errors. To apply rounding instead of significance use, e.g., \code{digits = "r3"} which will round at the first 3 decimals. If character, it must be of the form \code{"rd"} or \code{"sd"} with \code{d} a digit (\code{r} is for round and \code{s} is for significance). For the number of digits for the fit statistics, use \code{digits.stats}. Note that when significance is used it does not exactly display the number of significant digits: see details for its exact meaning.
#' @param digits.stats Integer or character scalar. Default is 5 and represents the number of significant digits to be displayed for the fit statistics. To apply rounding instead of significance use, e.g., \code{digits = "r3"} which will round at the first 3 decimals. If character, it must be of the form \code{"rd"} or \code{"sd"} with \code{d} a digit (\code{r} is for round and \code{s} is for significance). Note that when significance is used it does not exactly display the number of significant digits: see details for its exact meaning.
#' @param tex Logical: whether the results should be a data.frame or a Latex table. By default, this argument is \code{TRUE} if the argument \code{file} (used for exportation) is not missing; it is equal to \code{FALSE} otherwise.
#' @param fitstat A character vector or a one sided formula (both with only lowercase letters). A vector listing which fit statistics to display. The valid types are 'n', 'll', 'aic', 'bic' and r2 types like 'r2', 'pr2', 'war2', etc (see all valid types in \code{\link[fixest]{r2}}). Also accepts valid types from the function \code{\link[fixest]{fitstat}}. The default value depends on the models to display. Example of use: \code{fitstat=c('n', 'cor2', 'ar2', 'war2')}, or \code{fitstat=~n+cor2+ar2+war2} using a formula. You can use the dot to refer to default values:\code{ ~ . + ll} would add the log-likelihood to the default fit statistics.
#' @param title (Tex only.) Character scalar. The title of the Latex table.
#' @param float (Tex only.) Logical. By default, if the argument \code{title} or \code{label} is provided, it is set to \code{TRUE}. Otherwise, it is set to \code{FALSE}.
#' @param sdBelow Logical or \code{NULL} (default). Should the standard-errors be displayed below the coefficients? If \code{NULL}, then this is \code{TRUE} for Latex and \code{FALSE} otherwise.
#' @param keep Character vector. This element is used to display only a subset of variables. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be kept. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{keep = "x[[:digit:]]$"}. If the first character is an exclamation mark, the effect is reversed (e.g. keep = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is kept). See details.
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be discarded. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{drop = "x[[:digit:]]{2}"}. If the first character is an exclamation mark, the effect is reversed (e.g. drop = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is dropped). See details.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions. This argument is applied post aliasing (see argument \code{dict}). Example: you have the following variables: \code{month1} to \code{month6}, then \code{x1} to \code{x5}, then \code{year1} to \code{year6}. If you want to display first the x's, then the years, then the months you could use: \code{order = c("x", "year")}. If the first character is an exclamation mark, the effect is reversed (e.g. order = "!Intercept" means: every variable that does not contain \dQuote{Intercept} goes first).  See details.
#' @param dict A named character vector or a logical scalar. It changes the original variable names to the ones contained in the \code{dict}ionary. E.g. to change the variables named \code{a} and \code{b3} to (resp.) \dQuote{$log(a)$} and to \dQuote{$bonus^3$}, use \code{dict=c(a="$log(a)$",b3="$bonus^3$")}. By default, it is equal to \code{getFixest_dict()}, a default dictionary which can be set with \code{\link[fixest]{setFixest_dict}}. You can use \code{dict = FALSE} to disable it..
#' @param file A character scalar. If provided, the Latex (or data frame) table will be saved in a file whose path is \code{file}. If you provide this argument, then a Latex table will be exported, to export a regular \code{data.frame}, use argument \code{tex = FALSE}.
#' @param replace Logical, default is \code{FALSE}. Only used if option \code{file} is used. Should the exported table be written in a new file that replaces any existing file?
#' @param convergence Logical, default is missing. Should the convergence state of the algorithm be displayed? By default, convergence information is displayed if at least one model did not converge.
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)} for a Latex table and \code{c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)} for a data.frame (to conform with R's default). To suppress the significance codes, use \code{signifCode=NA} or \code{signifCode=NULL}. Can also be equal to \code{"letters"}, then the default becomes \code{c("a"=0.01, "b"=0.05, "c"=0.10)}.
#' @param label (Tex only.) Character scalar. The label of the Latex table.
#' @param subtitles Character vector or list. The elements should be of length 1 or of the same length as the number of models. If a list, the names of the list will be displayed on the leftmost column. By default it is equal to \code{list("auto")} which means that if the object is a split sample estimation, the sample will be automatically added as a sub-title.
#' @param fixef_sizes (Tex only.) Logical, default is \code{FALSE}. If \code{TRUE} and fixed-effects were used in the models, then the number of "individuals" per fixed-effect dimension is also displayed.
#' @param fixef_sizes.simplify Logical, default is \code{TRUE}. Only used if \code{fixef_sizes = TRUE}. If \code{TRUE}, the fixed-effects sizes will be displayed in parentheses instead of in a separate line if there is no ambiguity (i.e. if the size is constant across models).
#' @param family Logical, default is missing. Whether to display the families of the models. By default this line is displayed when at least two models are from different families.
#' @param keepFactors Logical, default is \code{TRUE}. If \code{FALSE}, then factor variables are displayed as fixed-effects and no coefficient is shown.
#' @param powerBelow (Tex only.) Integer, default is -5. A coefficient whose value is below \code{10**(powerBelow+1)} is written with a power in Latex. For example \code{0.0000456} would be written \code{4.56$\\times 10^{-5}$} by default. Setting \code{powerBelow = -6} would lead to \code{0.00004} in Latex.
#' @param interaction.combine (Tex only.) Character scalar, defaults to \code{" $\\times$ "}. When the estimation contains interactions, then the variables names (after aliasing) are combined with this argument. For example: if \code{dict = c(x1="Wind", x2="Rain")} and you have the following interaction \code{x1:x2}, then it will be renamed (by default) \code{Wind $\\times$ Rain} -- using \code{interaction.combine = "*"} would lead to \code{Wind*Rain}.
#' @param depvar Logical, default is \code{TRUE}. Whether a first line containing the dependent variables should be shown.
#' @param coefstat One of \code{"se"} (default), \code{"tstat"} or \code{"confint"}. The statistic to report for each coefficient: the standard-error, the t-statistics or the confidence interval. You can adjust the confidence interval with the argument \code{ci}.
#' @param ci Level of the confidence interval, defaults to \code{0.95}. Only used if \code{coefstat = confint}.
#' @param style.tex An object created by the function \code{\link[fixest]{style.tex}}. It represents the style of the Latex table, see the documentation of \code{\link[fixest]{style.tex}}.
#' @param style.df An object created by the function \code{\link[fixest]{style.df}. }It represents the style of the data frame returned (if \code{tex = FALSE}), see the documentation of \code{\link[fixest]{style.df}}.
#' @param notes (Tex only.) Character vector. If provided, a \code{"notes"} section will be added at the end right after the end of the table, containing the text of this argument. Note that if it is a vector, it will be collapsed with new lines.
#' @param group A list. The list elements should be vectors of regular expressions. For each elements of this list: A new line in the table is created, all variables that are matched by the regular expressions are discarded (same effect as the argument \code{drop}) and \code{TRUE} or \code{FALSE} will appear in the model cell, depending on whether some of the previous variables were found in the model. Example: \code{group=list("Controls: personal traits"=c("gender", "height", "weight"))} will create an new line with \code{"Controls: personal traits"} in the leftmost cell, all three variables gender, height and weight are discarded, \code{TRUE} appearing in each model containing at least one of the three variables (the style of \code{TRUE}/\code{FALSE} is governed by the argument \code{yesNo}). You can control the placement of the new row by using 1 or 2 special characters at the start of the row name. The meaning of these special characters are: 1) \code{"^"}: 1st, \code{"_"}: last, row; 2) \code{"^"}: coef., \code{"-"}: fixed-effect, \code{"_"}: stats, section. For example: \code{group=list("_Controls"=stuff)} will place the line last in the stats sections, and using \code{group=list("^_Controls"=stuff)} will make the row appear first in the stats section. For details, see the dedicated section.
#' @param extraline A list or a one sided formula. The list elements should be either a logical scalar,a vector of the same length as the number of models, or a function. This argument can be many things, please have a look at the dedicated help section; a simplified description follows. For each elements of this list: A new line in the table is created, the list name being the row name and the vector being the content of the cells. Example: \code{extraline=list("Sub-sample"=c("<20 yo", "all", ">50 yo"))} will create an new line with \code{"Sub-sample"} in the leftmost cell, the vector filling the content of the cells for the three models. You can control the placement of the new row by using 1 or 2 special characters at the start of the row name. The meaning of these special characters are: 1) \code{"^"}: 1st, \code{"_"}: last, row; 2) \code{"^"}: coef., \code{"-"}: fixed-effect, \code{"_"}: stats, section. For example: \code{group=list("_Controls"=stuff)} will place the line last in the stats sections, and using \code{group=list("^_Controls"=stuff)} will make the row appear first in the stats section. For details, see the dedicated section.
#' @param fixef.group Logical scalar or list (default is \code{NULL}). If equal to \code{TRUE}, then all fixed-effects always appearing jointly in models will be grouped in one row. If a list, its elements must be character vectors of regular expressions and the list names will be the row names. For ex. \code{fixef.group=list("Dates fixed-effects"="Month|Day")} will remove the \code{"Month"} and \code{"Day"} fixed effects from the display and replace them with a single row named "Dates fixed-effects". You can monitor the placement of the new row with the special characters telling where to place the row within a section: \code{"^"} (first), or \code{"_"} (last); and in which section it should appear: \code{"^"} (coef.), \code{"-"} (fixed-effects), or \code{"_"} (stat.). These two special characters must appear first in the row names. Please see the dedicated section.
#' @param placement (Tex only.) Character string giving the position of the float in Latex. Default is "htbp". It must consist of only the characters 'h', 't', 'b', 'p', 'H' and '!'. Reminder: h: here; t: top; b: bottom; p: float page; H: definitely here; !: prevents Latex to look for other positions. Note that it can be equal to the empty string (and you'll get the default placement).
#' @param drop.section Character vector which can be of length 0 (i.e. equal to \code{NULL}). Can contain the values "fixef", "slopes" or "stats". It would drop, respectively, the fixed-effects section, the variables with varying slopes section or the fit statistics section.
#' @param reset (\code{setFixest_etable} only.) Logical, default is \code{FALSE}. If \code{TRUE}, this will reset all the default values that were already set by the user in previous calls.
#' @param .vcov A function to be used to compute the standard-errors of each fixest object. You can pass extra arguments to this function using the argument \code{.vcov_args}. See the example.
#' @param .vcov_args A list containing arguments to be passed to the function \code{.vcov}.
#' @param poly_dict Character vector, default is \code{c("", " square", " cube")}. When raw polynomials (\code{x^2}, etc) are used, the variables are automatically renamed and \code{poly_dict} rules the display of the power. For powers greater than the number of elements of the vector, the value displayed is \code{$^{pow}$} in Latex and \code{^ pow} in the R console.
#' @param postprocess.tex A function that will postprocess the character vector defining the latex table. Only when \code{tex = TRUE}. By default it is equal to \code{NULL}, meaning that there is no postprocessing. When \code{tex = FALSE}, see the argument \code{postprocess.df}. See details.
#' @param postprocess.df A function that will postprocess.tex the resulting data.frame. Only when \code{tex = FALSE}. By default it is equal to \code{NULL}, meaning that there is no postprocessing. When \code{tex = TRUE}, see the argument \code{postprocess.tex}.
#' @param fit_format Character scalar, default is \code{"__var__"}. Only used in the presence of IVs. By default the endogenous regressors are named \code{fit_varname} in the second stage. The format of the endogenous regressor to appear in the table is governed by \code{fit_format}. For instance, by default, the prefix \code{"fit_"} is removed, leading to only \code{varname} to appear. If \code{fit_format = "$\\\\hat{__var__}$"}, then \code{"$\\hat{varname}$"} will appear in the table.
#' @param coef.just (DF only.) Either \code{"."}, \code{"("}, \code{"l"}, \code{"c"} or \code{"r"}, default is \code{NULL}. How the coefficients should be justified. If \code{NULL} then they are right aligned if \code{sdBelow = FALSE} and aligned to the dot if \code{sdBelow = TRUE}. The keywords stand respectively for dot-, parenthesis-, left-, center- and right-aligned.
#'
#' @details
#' The function \code{esttex} is equivalent to the function \code{etable} with argument \code{tex = TRUE}.
#'
#' The function \code{esttable} is equivalent to the function \code{etable} with argument \code{tex = FALSE}.
#'
#' You can permanently change the way your table looks in Latex by using \code{setFixest_etable}. The following vignette gives an example as well as illustrates how to use the \code{style} and postprocessing functions: \href{https://lrberge.github.io/fixest/articles/exporting_tables.html}{Exporting estimation tables}.
#'
#' When the argument \code{postprocessing.tex} is not missing, two additional tags will be included in the character vector returned by \code{etable}: \code{"\%start:tab\\n"} and \code{"\%end:tab\\n"}. These can be used to identify the start and end of the tabular and are useful to insert code within the \code{table} environment.
#'
#' @section How does \code{digits} handle the number of decimals displayed?:
#'
#' The default display of decimals is the outcome of an algorithm. Let's take the example of \code{digits = 3} which "kind of" requires 3 significant digits to be displayed.
#'
#' For numbers greater than 1 (in absolute terms), their integral part is always displayed and the number of decimals shown is equal to \code{digits} minus the number of digits in the integral part. This means that \code{12.345} will be displayed as \code{12.3}. If the number of decimals should be 0, then a single decimal is displayed to suggest that the number is not whole. This means that \code{1234.56} will be displayed as \code{1234.5}. Note that if the number is whole, no decimals are shown.
#'
#' For numbers lower than 1 (in absolute terms), the number of decimals displayed is equal to \code{digits} except if there are only 0s in which case the first significant digit is shown. This means that \code{0.01234} will be displayed as \code{0.012} (first rule), and that 0.000123 will be displayed as \code{0.0001} (second rule).
#'
#' @section Arguments keep, drop and order:
#' The arguments \code{keep}, \code{drop} and \code{order} use regular expressions. If you are not aware of regular expressions, I urge you to learn it, since it is an extremely powerful way to manipulate character strings (and it exists across most programming languages).
#'
#' For example drop = "Wind" would drop any variable whose name contains "Wind". Note that variables such as "Temp:Wind" or "StrongWind" do contain "Wind", so would be dropped. To drop only the variable named "Wind", you need to use \code{drop = "^Wind$"} (with "^" meaning beginning, resp. "$" meaning end, of the string => this is the language of regular expressions).
#'
#' Although you can combine several regular expressions in a single character string using pipes, \code{drop} also accepts a vector of regular expressions.
#'
#' You can use the special character "!" (exclamation mark) to reverse the effect of the regular expression (this feature is specific to this function). For example \code{drop = "!Wind"} would drop any variable that does not contain "Wind".
#'
#' You can use the special character "\%" (percentage) to make reference to the original variable name instead of the aliased name. For example, you have a variable named \code{"Month6"}, and use a dictionary \code{dict = c(Month6="June")}. Thus the variable will be displayed as \code{"June"}. If you want to delete that variable, you can use either \code{drop="June"}, or \code{drop="\%Month6"} (which makes reference to its original name).
#'
#' The argument \code{order} takes in a vector of regular expressions, the order will follow the elements of this vector. The vector gives a list of priorities, on the left the elements with highest priority. For example, order = c("Wind", "!Inter", "!Temp") would give highest priorities to the variables containing "Wind" (which would then appear first), second highest priority is the variables not containing "Inter", last, with lowest priority, the variables not containing "Temp". If you had the following variables: (Intercept), Temp:Wind, Wind, Temp you would end up with the following order: Wind, Temp:Wind, Temp, (Intercept).
#'
#' @section The argument \code{extraline}:
#'
#' The argument \code{extraline} adds well... extra lines to the table. It accepts either a list, or a one-sided formula.
#'
#' If a one-sided formula: then the elements in the formula must represent either \code{extraline} macros, either fit statistics (i.e. valid types of the function \code{\link[fixest]{fitstat}}). One new line will be added for each element of the formula, they will be appended right after the coefficients. To register \code{extraline} macros, you must first register them in \code{\link[fixest]{extraline_register}}.
#'
#' If a list: then the elements must be either: a) a logical scalar, b) a vector of the same length as the number of models, c) a function to be applied to each model and which returns a scalar, or d) a one-sided formula of \code{extraline} macros (registered with \code{\link[fixest]{extraline_register}}) or valid \code{\link[fixest]{fitstat}} types.
#'
#' When a list, the elements of type a), b) and c) must have a name! Example: \code{extraline = list("Controls" = TRUE)} is good, \code{extraline = list(TRUE)} is bad. The name for element of type d) is optional. Example: \code{extraline = list("^^My F-stat" = ~f)} will add the row named \code{"My F-stat"} on top of the coefficients (for the placement, see the dedicated section), while \code{extraline = list(~f)} will add the F-stat row with its default name (\code{"F-test"}) at its default placement (below the coefficients).
#'
#' You can combine elements in the list, so that the following \code{extraline = list(~r2, "Controls" = TRUE), ~f+ivf} is valid.
#'
#' @section Controlling the placement of extra lines:
#'
#' The arguments \code{group}, \code{extraline} and \code{fixef.group} allow to add customized lines in the table. They can be defined via a list where the list name will be the row name. By default, the placement of the extra line is right after the coefficients (except for \code{fixef.group}, covered in the last paragraph). For instance, \code{group = list("Controls" = "x[[:digit:]]")} will create a line right after the coefficients telling which models contain the control variables.
#'
#' But the placement can be customized. The previous example (of the controls) will be used for illustration (the mechanism for \code{extraline} and \code{fixef.group} is identical).
#'
#' The row names accept 2 special characters at the very start. The first governs the placement of the new line within the section: it can be equal to \code{"^"}, meaning first line, or \code{"_"}, meaning last line. The second character tells in which section the line should appear: it can be equal to \code{"^"}, \code{"-"}, or \code{"_"}, meaning respectively the coefficients, the fixed-effects and the statistics section (which typically appear at the top, mid and bottom of the table).
#'
#' Let's have some examples. Using the previous example, writing \code{"^_Controls"} would place the new line in at the top of the statistics section. Writing \code{"_-Controls"} places it as the last row of the fixed-effects section; \code{"^^Controls"} at the top row of the coefficients section; etc...
#'
#' On top of that there are shortcuts to avoid writing the two special characters. If only one special character is found, it is assumed to reflect the section, unless it corresponds to the default section, a case where the character then reflects the position within the section. An example will make it clear. Writing \code{"^Controls"} would place the line at the top of the coefficients section (since the default placement would have been this section). \code{"_Controls"} would place it at the bottom of the statistics section.
#'
#' The placement in \code{fixef.group} is defined similarly, only the default placement is different. Its default placement is at the top of the fixed-effects section. This means that the only difference with \code{group} and \code{extraline} is when a single special character is used. Here using, e.g., \code{"^My FEs"} would place the row at the bottom of the coefficients section, since \code{"^"} would refer to the section (and not the row within the section).
#'
#' @return
#' If \code{tex = TRUE}, the lines composing the Latex table are returned invisibly while the table is directly prompted on the console.
#'
#' If \code{tex = FALSE}, the data.frame is directly returned. If the argument \code{file} is not missing, the \code{data.frame} is printed and returned invisibly.
#'
#' @seealso
#' See also the main estimation functions \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}. Use \code{\link[fixest]{summary.fixest}} to see the results with the appropriate standard-errors, \code{\link[fixest]{fixef.fixest}} to extract the fixed-effects coefficients.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' aq = airquality
#'
#' est1 = feols(Ozone ~ i(Month) / Wind + Temp, data = aq)
#' est2 = feols(Ozone ~ i(Month, Wind) + Temp | Month, data = aq)
#'
#' # Displaying the two results in a single table
#' etable(est1, est2)
#'
#' # keep/drop: keeping only interactions
#' etable(est1, est2, keep = " x ")
#' # or using drop  (see regexp help):
#' etable(est1, est2, drop = "^(Month|Temp|\\()")
#'
#' # keep/drop: dropping interactions
#' etable(est1, est2, drop = " x ")
#' # or using keep ("!" reverses the effect):
#' etable(est1, est2, keep = "! x ")
#'
#' # order: Wind variable first, intercept last (note the "!" to reverse the effect)
#' etable(est1, est2, order = c("Wind", "!Inter"))
#' # Month, then interactions, then the rest
#' etable(est1, est2, order = c("^Month", " x "))
#'
#' #
#' # dict
#' #
#'
#' # You can rename variables with dict = c(var1 = alias1, var2 = alias2, etc)
#' # You can also rename values taken by factors.
#' # Here's a full example:
#' dict = c(Temp = "Temperature", "Month::5"="May", "6"="Jun")
#' etable(est1, est2, dict = dict)
#' # Note the difference of treatment between Jun and May
#'
#' # Assume the following dictionary:
#' dict = c("Month::5"="May", "Month::6"="Jun", "Month::7"="Jul",
#'          "Month::8"="Aug", "Month::9"="Sep")
#'
#' # We would like to keep only the Months, but now the names are all changed...
#' # How to do?
#' # We can use the special character '%' to make reference to the original names.
#'
#' etable(est1, est2, dict = dict, keep = "%Month")
#'
#' #
#' # signifCode
#' #
#'
#' etable(est1, est2, signifCode = c(" A"=0.01, " B"=0.05, " C"=0.1, " D"=0.15, " F"=1))
#'
#' #
#' # Using the argument style to customize Latex exports
#' #
#'
#' # If you don't like the default layout of the table, no worries!
#' # You can modify many parameters with the argument style
#'
#' # To drop the headers before each section, use:
#' # Note that a space adds an extra line
#' style_noHeaders = style.tex(var.title = "", fixef.title = "", stats.title = " ")
#' etable(est1, est2, dict = dict, tex = TRUE, style.tex = style_noHeaders)
#'
#' # To change the lines of the table + dropping the table footer
#' style_lines = style.tex(line.top = "\\toprule", line.bottom = "\\bottomrule",
#'                     tablefoot = FALSE)
#' etable(est1, est2, dict = dict, tex = TRUE, style.tex = style_lines)
#'
#' # Or you have the predefined type "aer"
#' etable(est1, est2, dict = dict, tex = TRUE, style.tex = style.tex("aer"))
#'
#' #
#' # Group and extraline
#' #
#'
#' # Sometimes it's useful to group control variables into a single line
#' # You can achieve that with the group argument
#'
#' setFixest_fml(..ctrl = ~ poly(Wind, 2) + poly(Temp, 2))
#' est_c0 = feols(Ozone ~ Solar.R, data = aq)
#' est_c1 = feols(Ozone ~ Solar.R + ..ctrl, data = aq)
#' est_c2 = feols(Ozone ~ Solar.R + Solar.R^2 + ..ctrl, data = aq)
#'
#' etable(est_c0, est_c1, est_c2, group = list(Controls = "poly"))
#'
#' # 'group' here does the same as drop = "poly", but adds an extra line
#' # with TRUE/FALSE where the variables were found
#'
#' # 'extraline' adds an extra line, where you can add the value for each model
#' est_all  = feols(Ozone ~ Solar.R + Temp + Wind, data = aq)
#' est_sub1 = feols(Ozone ~ Solar.R + Temp + Wind, data = aq[aq$Month %in% 5:6, ])
#' est_sub2 = feols(Ozone ~ Solar.R + Temp + Wind, data = aq[aq$Month %in% 7:8, ])
#' est_sub3 = feols(Ozone ~ Solar.R + Temp + Wind, data = aq[aq$Month == 9, ])
#'
#' etable(est_all, est_sub1, est_sub2, est_sub3,
#'        extraline = list("Sub-sample" = c("All", "May-June", "Jul.-Aug.", "Sept.")))
#'
#' # You can monitor the placement of the new lines with two special characters
#' # at the beginning of the row name.
#' # 1) "^" or "_" which mean first or last line of the section
#' # 2) "^", "-" or "_" which mean the coefficients, the fixed-effects or the
#' # statistics section.
#' #
#' # Ex: starting with "^_" will place the line at the top of the stat. section
#' #     starting with "_-" will place the line at the bottom of the FEs section
#' #     etc.
#' #
#' # You can use a single character which will represent the section, unless
#' # it's
#'
#' # Examples
#' etable(est_c0, est_c1, est_c2, group = list("_Controls" = "poly"))
#' etable(est_all, est_sub1, est_sub2, est_sub3,
#'        extraline = list("^Sub-sample" = c("All", "May-June", "Jul.-Aug.", "Sept.")))
#' # Note that since the default placement is the coefficients section,
#' # a single "^" then refers to the position within the section. We end
#' # up at the top of the coefficients section.
#'
#' #
#' # fixef.group
#' #
#'
#' # You can group the fixed-effects line with fixef.group
#'
#' est_0fe = feols(Ozone ~ Solar.R + Temp + Wind, aq)
#' est_1fe = feols(Ozone ~ Solar.R + Temp + Wind | Month, aq)
#' est_2fe = feols(Ozone ~ Solar.R + Temp + Wind | Month + Day, aq)
#'
#' # A) automatic way => simply use fixef.group = TRUE
#'
#' etable(est_0fe, est_2fe, fixef.group = TRUE)
#'
#' # Note that when grouping would lead to inconsistencies across models,
#' # it is avoided
#'
#' etable(est_0fe, est_1fe, est_2fe, fixef.group = TRUE)
#'
#' # B) customized way => use a list
#'
#' etable(est_0fe, est_2fe, fixef.group = list("Dates" = "Month|Day"))
#'
#' # Note that when a user grouping would lead to inconsistencies,
#' # the term partial replaces yes/no and the fixed-effects are not removed.
#'
#' etable(est_0fe, est_1fe, est_2fe, fixef.group = list("Dates" = "Month|Day"))
#'
#' # Using customized placement => as with 'group' and 'extraline',
#' # the user can control the placement of the new line.
#' # See the previous 'group' examples and the dedicated section in the help.
#'
#' # On top of the coefficients:
#' etable(est_0fe, est_2fe, fixef.group = list("^^Dates" = "Month|Day"))
#'
#' # Last line of the statistics
#' etable(est_0fe, est_2fe, fixef.group = list("_Dates" = "Month|Day"))
#'
#'
#'
#' #
#' # Using custom functions to compute the standard errors
#' #
#'
#' # You can customize the way you compute the SEs with the argument .vcov
#' # Let's use some covariances from the sandwich package
#'
#' etable(est_c0, est_c1, est_c2, .vcov = sandwich::vcovHC)
#'
#' # To add extra arguments to vcovHC, you need to use .vcov_args
#' etable(est_c0, est_c1, est_c2, .vcov = sandwich::vcovHC, .vcov_args = list(type = "HC0"))
#'
#'
#' #
#' # Customize which fit statistic to display
#' #
#'
#' # You can change the fit statistics with the argument fitstat
#' # and you can rename them with the dictionnary
#' etable(est1, est2, fitstat = ~ r2 + n + G)
#'
#' # If you use a formula, '.' means the default:
#' etable(est1, est2, fitstat = ~ ll + .)
#'
#'
#' #
#' # Computing a different SE for each model
#' #
#'
#' est = feols(Ozone ~ Solar.R + Wind + Temp, data = aq)
#'
#' #
#' # Method 1: use summary
#'
#' s1 = summary(est, "standard")
#' s2 = summary(est, cluster = ~ Month)
#' s3 = summary(est, cluster = ~ Day)
#' s4 = summary(est, cluster = ~ Day + Month)
#'
#' etable(list(s1, s2, s3, s4))
#'
#' #
#' # Method 2: using a list in the argument 'cluster'
#'
#' est_bis = feols(Ozone ~ Solar.R + Wind + Temp | Month, data = aq)
#' etable(list(est, est_bis), cluster = list("standard", ~ Month))
#'
#' #
#' # Method 3: Using rep()
#'
#' etable(rep(est, cluster = list("standard", ~ Month)))
#'
#' # When using rep on 2 or more objects, you need to embed them in .l()
#' etable(rep(.l(est, est_bis), cluster = list("standard", ~ Month, ~ Day)))
#'
#' # Using each to order differently
#' etable(rep(.l(est, est_bis), each = 3, cluster = list("standard", ~ Month, ~ Day)))
#'
#'
etable = function(..., se = NULL, dof = NULL, cluster = NULL, stage = 2, agg = NULL,
                  .vcov, .vcov_args = NULL, digits = 4, digits.stats = 5, tex,
                  fitstat, title, coefstat = "se", ci = 0.95, sdBelow = NULL,
                  keep, drop, order, dict, file, replace = FALSE, convergence,
                  signifCode, label, float, subtitles = list("auto"), fixef_sizes = FALSE,
                  fixef_sizes.simplify = TRUE, keepFactors = TRUE, family, powerBelow = -5,
                  interaction.combine = " $\\times $ ", depvar = TRUE, style.tex = NULL,
                  style.df = NULL, notes = NULL, group = NULL, extraline = NULL,
                  fixef.group = NULL, placement = "htbp", drop.section = NULL,
                  poly_dict = c("", " square", " cube"), postprocess.tex = NULL,
                  postprocess.df = NULL, fit_format = "__var__", coef.just = NULL){

    #
    # Checking the arguments
    #

    # Need to check for the presence of the se
    useSummary = TRUE
    if(missnull(se) && missnull(cluster) && missing(.vcov) && missing(stage) && missnull(agg)){
        useSummary = FALSE
    }

    if(!missnull(se)){
        check_arg_plus(se, "match(standard, white, hetero, cluster, twoway, threeway, fourway)")
    } else {
        se = NULL
    }

    # The depvar
    if(missing(depvar) && !missing(file)){
        depvar = TRUE
    }

    check_arg(tex, "logical scalar")
    if(missing(tex)){
        if(!missing(file)) {
            tex = TRUE
        } else {
            tex = FALSE
        }
    }

    # The signif codes
    if(missing(signifCode)){
        if(tex){
            signifCode = c("***"=0.01, "**"=0.05, "*"=0.10)
        } else {
            signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)
        }
    }

    # Float or not
    check_arg(float, "logical scalar")
    if(missing(float)){
        if(!missing(title) || !missing(label)){
            float = TRUE
        } else {
            float = FALSE
        }
    } else if(float && (!missing(title) || !missing(label))) {
        what = c("title", "label")[c(!missing(title), !missing(label))]
        warning("Since float = TRUE, the argument", enumerate_items(what, "s.is"), " ignored", immediate. = TRUE, call. = FALSE)
    }

    if(missing(dict) && !tex){
        dict = TRUE
    }

    dots = error_sender(list(...), "Some elements in '...' could not be evaluated: ")

    if(".up" %in% names(dots)){
        # internal call from esttable/esttex
        .up = 2
        dots[[".up"]] = NULL
    } else {
        .up = 1
    }

    if(length(dots) == 0){
        stop("You must provide at least one element in '...'.")
    }

    # Getting the model names
    if(.up == 2){
        # it's pain in the necky
        sysOrigin = sys.parent()
        mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
        dots_call = mc[["..."]]
    } else {
        dots_call = match.call(expand.dots = FALSE)[["..."]]
    }


    #
    # postprocess.tex
    #

    # Note that I need to catch the arguments from the two pp functions
    # => so that the same call lead to valid evaluations

    opts = getOption("fixest_etable")
    pp_tex = pp_df = NULL

    check_arg(postprocess.tex, "NULL function arg(1,)")

    if(!is.null(postprocess.tex)){
        pp_tex = postprocess.tex
    } else if(!is.null(opts$postprocess.tex)){
        pp_tex = opts$postprocess.tex
    }

    check_arg(postprocess.df, "NULL function arg(1,)")
    if(!is.null(postprocess.df)){
        pp_df = postprocess.df
    } else if(!is.null(opts$postprocess.df)){
        pp_df = opts$postprocess.df
    }

    my_postprocess = pp_other = NULL
    if(tex){
        my_postprocess = pp_tex
        pp_other = pp_df
    } else {
        my_postprocess = pp_df
        pp_other = pp_tex
    }

    DO_POSTPROCESS = !is.null(my_postprocess)

    # Catching the arguments
    pp_args = list()
    if(!is.null(names(dots))){
        # Catching the arguments
        if(!is.null(my_postprocess)){
            fm = formalArgs(my_postprocess)
            qui = names(dots) %in% fm
            pp_args = dots[qui]
            dots = dots[!qui]
        }

        if(!is.null(pp_other)){
            fm = formalArgs(pp_other)
            qui = names(dots) %in% fm
            dots = dots[!qui]
        }

    }

    if(length(dots) == 0){
        stop("After cleaning the arguments to the postprocessing functions, there is no element left in '...'. Please provide at least one.")
    }

    for(i in seq_along(dots)){
        if(!any(c("fixest", "fixest_list", "fixest_multi") %in% class(dots[[i]])) && !is.list(dots[[i]])){
            msg = ""
            if(!is.null(names(dots))){
                msg = paste0(" (named '", names(dots)[i], "')")
            }
            stop("The ", n_th(i), " element of '...'", msg, " is not valid: it should be a fixest object or a list, it is neither.")
        }
    }

    info = results2formattedList(dots = dots, se=se, dof=dof, fitstat_all=fitstat,
                                 cluster=cluster, stage=stage, agg = agg, .vcov=.vcov,
                                 .vcov_args=.vcov_args, digits=digits, digits.stats=digits.stats,
                                 sdBelow=sdBelow, signifCode=signifCode, coefstat = coefstat,
                                 ci = ci, title=title, float=float, subtitles=subtitles,
                                 keepFactors=keepFactors, tex = tex, useSummary=useSummary,
                                 dots_call=dots_call, powerBelow=powerBelow, dict=dict,
                                 interaction.combine=interaction.combine, convergence=convergence,
                                 family=family, keep=keep, drop=drop, file=file, order=order,
                                 label=label, fixef_sizes=fixef_sizes,
                                 fixef_sizes.simplify=fixef_sizes.simplify,
                                 depvar=depvar, style.tex=style.tex, style.df=style.df,
                                 replace=replace, notes = notes, group = group, extraline=extraline,
                                 fixef.group=fixef.group, placement = placement,
                                 drop.section = drop.section, poly_dict = poly_dict,
                                 tex_tag = DO_POSTPROCESS, fit_format = fit_format,
                                 coef.just = coef.just, .up = .up)

    if(tex){
        res = etable_internal_latex(info)
    } else {
        res = etable_internal_df(info)
    }

    if(!missnull(file)){
        sink(file = file, append = !replace)
        on.exit(sink())

        if(DO_POSTPROCESS){
            # res = my_postprocess(res)
            pp_args_all = list(res)
            if(length(pp_args) > 0){
                pp_args_all[names(pp_args)] = pp_args
            }
            res = do.call(my_postprocess, pp_args_all)

            if(is.null(res)){
                return(invisible(NULL))
            }

            if(tex){
                res = res[!res %in% c("%start:tab\n", "%end:tab\n")]
            }
        }

        if(tex){
            cat(res, sep = "")
            # We add extra whitespaces => otherwise the Latex file is a bit cluttered
            cat("\n\n")
        } else {
            if(is.data.frame(res)){
                print(res)
            } else {
                cat(res)
            }
        }

        return(invisible(res))
    } else {

        if(DO_POSTPROCESS){

            pp_args_all = list(res)
            if(length(pp_args) > 0){
                pp_args_all[names(pp_args)] = pp_args
            }
            res = do.call(my_postprocess, pp_args_all)

            if(is.null(res)){
                return(invisible(NULL))
            }

            if(tex){
                res = res[!res %in% c("%start:tab\n", "%end:tab\n")]
            }
        }

        if(tex){
            cat(res, sep = "")
            return(invisible(res))
        } else {
            if(is.data.frame(res)){
                return(res)
            } else {
                cat(res)
                return(invisible(res))
            }
        }
    }
}

gen_etable_aliases = function(){
    # esttex and esttable are the functions that existed when the package was launched
    # Now the two have been merged into etable
    # I like it much better
    # I wanted to deprecate them, but maintainance with that function is very easy

    etable_args = formals(etable)

    arg_name = names(etable_args)
    arg_default = sapply(etable_args, deparse_long)

    #
    # esttable
    #

    qui_df = !arg_name %in% c("tex", "title", "label", "float", "style.tex", "notes", "placement", "postprocess.tex")

    esttable_args = paste0(arg_name[qui_df], " = ", arg_default[qui_df], collapse = ", ")
    esttable_args = gsub(" = ,", ",", esttable_args)

    etable_call = paste0(arg_name[qui_df][-1], " = ", arg_name[qui_df][-1], collapse = ", ")

    esttable_fun = paste0("esttable = function(", esttable_args, "){\n\n",
                          "\tetable(..., ", etable_call, ", tex = FALSE, .up = 2)\n}")

    esttable_rox = "#' @describeIn etable Exports the results of multiple \\code{fixest} estimations in a Latex table."

    #
    # esttex
    #

    qui_tex = !arg_name %in% c("tex", "style.df", "postprocess.df", "coef.just")

    esttex_args = paste0(arg_name[qui_tex], " = ", arg_default[qui_tex], collapse = ", ")
    esttex_args = gsub(" = ,", ",", esttex_args)

    etable_call = paste0(arg_name[qui_tex][-1], " = ", arg_name[qui_tex][-1], collapse = ", ")

    esttex_fun = paste0("esttex = function(", esttex_args, "){\n\n",
                        "\tetable(..., ", etable_call, ", tex = TRUE, .up = 2)\n}")

    esttex_rox = "#' @describeIn etable Exports the results of multiple \\code{fixest} estimations in a Latex table."



    # Writing the functions

    f = file("R/etable_aliases.R", "w", encoding = "utf-8")

    intro = c("# Do not edit by hand\n# => aliases to the function etable\n\n\n")

    s = "\n\n\n\n"
    text = c(intro, s, esttable_rox, esttable_fun, s, esttex_rox, esttex_fun, s)
    writeLines(text, f)
    close(f)


}

results2formattedList = function(dots, se, dof = getFixest_dof(), cluster, stage = 2, agg = NULL, .vcov, .vcov_args = NULL, digits = 4, digits.stats = 5, fitstat_all, sdBelow=NULL, dict, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), coefstat = "se", ci = 0.95, label, subtitles, title, float = FALSE, replace = FALSE, keepFactors = FALSE, tex = FALSE, useSummary, dots_call, powerBelow = -5, interaction.combine, convergence, family, drop, order, keep, file, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, depvar = FALSE, style.tex = NULL, style.df=NULL, notes = NULL, group = NULL, extraline=NULL, fixef.group = NULL, placement = "htbp", drop.section = NULL, poly_dict = c("", " square", " cube"), tex_tag = FALSE, fit_format = "__var__", coef.just = NULL, .up = 1){

    # This function is the core of the function etable

    set_up(.up)

    # Setting the default values (we take extra care for "style")
    if(tex){
        check_arg(style.tex, "NULL class(fixest_style_tex)")
        # The variable style will be changed via the defaults
        style_user = style.tex
    } else {
        check_arg(style.df, "NULL class(fixest_style_df)")
        # The variable style will be changed via the defaults
        style_user = style.df
    }


    #
    # Setting the default
    #

    opts = getOption("fixest_etable")
    if(length(opts) > 0){
        sysOrigin = sys.parent(.up)
        mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
        args_usr = setdiff(names(mc), c("style.tex", "style.df"))

        # Arguments for which the defaults should not be changed in etable, tex = FALSE
        if(!tex) args_usr = c(args_usr, "signifCode")

        # We modify only non-user provided arguments
        for(v in names(opts)){
            if(!v %in% args_usr){
                assign(v, opts[[v]])
            }
        }
    }

    if(tex){
        if(!"style.tex" %in% names(opts)){
            style = fixest::style.tex(main = "base")
        } else {
            style = style.tex
        }
    } else if(!tex){
        if(!"style.df" %in% names(opts)){
            style = fixest::style.df()
        } else {
            # We rename style.df into style
            style = style.df
        }
    }

    # Setting the style defaults
    # We modify the default set with setFixest_etable()
    if(length(style_user) > 0){
        style[names(style_user)] = style_user
    }


    #
    # Full control
    #

    check_arg(title, "character scalar")
    check_arg_plus(coefstat, "match(se, tstat, confint)")

    check_arg_plus(notes, "NULL{''} character vector no na")
    if(length(notes) > 1) notes = paste(notes, collapse = "\n")

    check_arg("logical scalar", replace, convergence, fixef_sizes, fixef_sizes.simplify, keepFactors, family, tex, depvar)
    check_arg("NULL logical scalar", sdBelow)

    isTex = tex
    if(missing(family)){
        show_family = NULL
    } else {
        show_family = family
    }

    if(is.null(sdBelow)){
        sdBelow = isTex
    }

    # start: coef.just
    check_arg(coef.just, "NULL charin", .choices = c(".", "(", "l", "c", "r"))

    if(!isTex){
        # this parameter is only used in DF
        if(is.null(coef.just) && sdBelow){
            coef.just = "."
        }

        if(coefstat == "confint" && identical(coef.just, ".")){
            coef.just = "c"
        }

        if(identical(coef.just, "r")){
            coef.just = NULL
        }
    }
    #   end: coef.just

    # digits argument + formatting
    digits_list = check_set_digits(digits, up = 2)
    digits = digits_list$digits
    round = digits_list$round

    fun_format = function(x) format_number(x, digits = digits, round = round, pow_below = powerBelow, tex = isTex)

    digits_list = check_set_digits(digits.stats, up = 2)
    digits.stats = digits_list$digits
    round.stats = digits_list$round

    fun_format_stats = function(x) format_number(x, digits = digits.stats, round = round.stats, pow_below = powerBelow, tex = isTex)

    # we rename the argument
    if(missing(depvar)){
        # Default is set in the end
        show_depvar = NULL
    } else {
        show_depvar = depvar
    }

    check_arg(keep, drop, order, "character vector no na NULL", .message = "The arg. '__ARG__' must be a vector of regular expressions (see help(regex)).")

    check_arg(file, label, interaction.combine, "character scalar")
    check_arg_plus(signifCode, "NULL NA | match(letters) | named numeric vector no na GE{0} LE{1}")
    check_arg_plus(subtitles, "NULL{list()} character vector no na | NA | list")

    check_arg(poly_dict, "character vector no na")

    check_arg(powerBelow, "integer scalar LE{-1}")

    check_arg(ci, "numeric scalar GT{0.5} LT{1}")

    check_arg(dict, "NULL logical scalar | named character vector no na")

    check_arg(placement, "character scalar")
    if(isTex){
        if(nchar(placement) > 0){
            p_split = strsplit(placement, "")[[1]]
            check_value(p_split, "strict multi charin(h, t, b, p, H, !)", .message = "Argument 'placement' must be a character string containing only the following characters: 'h', 't', 'b', 'p', 'H', and '!'.")
        }
    }

    check_arg_plus(drop.section, "NULL multi match(fixef, slopes, stats)")

    check_arg_plus(group, "NULL{list()} named list l0")
    check_arg_plus(extraline, "NULL{list()} list l0 | os formula")
    # we check it more in depth later

    check_arg(fixef.group, "NULL{list()} logical scalar | named list l0")
    if(isFALSE(fixef.group)){
        fixef.group = list()
    }

    IS_FIXEF_GROUP = length(fixef.group) > 0

    check_arg(fit_format, "character scalar")
    if(!grepl("__var__", fit_format, fixed = TRUE)){
        stop("The argument 'fit_format' should include the special name '__var__' that will be replaced by the variable name. So far it does not contain it.")
    }


    # yesNo => used later when setting the fixed-effect line
    yesNo = style$yesNo

    # default values for dict
    if(missing(dict)){
        if(isTex || !missing(file)){
            dict = getFixest_dict()
        } else {
            dict = NULL
        }
    } else if(isTRUE(dict)) {
        dict = getFixest_dict()
        if(is.null(dict)){
            dict = c("____" = "OH OH OH")
        }
    } else if(isFALSE(dict)) {
        dict = NULL
    }

    add_signif = TRUE
    if(identical(signifCode, "letters")){
        signifCode = c("a"=0.01, "b"=0.05, "c"=0.10)
    } else if(length(signifCode) == 0 || (length(signifCode) == 1 && is.na(signifCode))){
        add_signif = FALSE
    } else {
        signifCode = sort(signifCode)
    }

    # subtitles => must be a list
    # We get the automatic substitles, if split is used
    AUTO_SUBTITLES = FALSE
    if(is.list(subtitles)){
        if(length(subtitles) > 0){
            qui = sapply(subtitles, function(x) identical(x, "auto"))
            if(any(qui)){
                subtitles = subtitles[!qui]
                AUTO_SUBTITLES = TRUE
            }
        }
    } else if(anyNA(subtitles)){
        subtitles = list()
    } else {
        # It's a character vector
        if(identical(subtitles, "auto")){
            subtitles = list()
            AUTO_SUBTITLES = TRUE
        } else {
            subtitles = list(subtitles)
        }
    }

    # formatting the names of the models
    dots_names = names(dots_call)
    if(!is.null(dots_names)){

        for(i in 1:length(dots_call)){
            if(dots_names[i] != ""){
                dots_call[[i]] = dots_names[i]
            } else {
                dots_call[[i]] = deparse_long(dots_call[[i]])
            }
        }
    }

    n = length(dots)

    if(n == 0) stop_up("Not any estimation as argument.")

    all_models = list()
    model_names = list()
    auto_subtitles = list()
    model_id = NULL
    k = 1
    for(i in 1:n){
        di = dots[[i]]

        if("fixest_multi" %in% class(di)){
            meta = attr(di, "meta")
            di = attr(di, "data")

            if(AUTO_SUBTITLES){
                if(!is.null(meta$all_names[["sample"]])){
                    my_subtitles = list()

                    if(tex == FALSE){
                        # We need to rename to avoid FE problem duplicate row names
                        my_subtitles$title = paste0("Sample (", dict_apply(meta$all_names$split.name, dict), ")")
                    } else {
                        my_subtitles$title = dict_apply(meta$all_names$split.name, dict)
                    }

                    n_mod = length(di)

                    if(!"sample" %in% names(meta$index)){
                        my_subtitles$value = rep(meta$all_names$sample, n_mod)
                    } else {
                        my_subtitles$value = meta$all_names$sample[meta$tree[, "sample"]]
                    }

                    my_subtitles$index = seq(k, length.out = n_mod)

                    auto_subtitles[[length(auto_subtitles) + 1]] = my_subtitles

                }
            }

        }

        if("fixest" %in% class(di)){
            all_models[[k]] = di
            if(any(class(dots_call[[i]]) %in% c("call", "name"))){
                model_names[[k]] = deparse_long(dots_call[[i]])
            } else {
                model_names[[k]] = as.character(dots_call[[i]])
            }

            k = k+1
        } else if(any(c("list", "fixest_list") %in% class(di))){
            # we get into this list to get the fixest objects
            types = sapply(di, function(x) class(x)[1])
            qui = which(types == "fixest")
            for(m in qui){
                all_models[[k]] = di[[m]]

                # handling names
                if(n > 1){
                    if(is.null(names(di)[m]) || names(di)[m]==""){
                        model_names[[k]] = paste0(dots_call[[i]], "[[", m, "]]")
                    } else {
                        model_names[[k]] = paste0(dots_call[[i]], "$", names(di)[m])
                    }
                } else {
                    model_names[[k]] = as.character(names(di)[m])
                    id = di[[m]]$model_id
                    if(!is.null(id)){
                        model_id[k] = id
                    }
                }

                k = k+1
            }
        }

    }

    if(length(all_models)==0) stop_up("Not any 'fixest' model as argument!")

    n_models = length(all_models)

    auto_subtitles_clean = list()
    if(length(auto_subtitles) > 0){
        # We reconstruct the subtitles properly

        for(i in seq_along(auto_subtitles)){
            my_sub = auto_subtitles[[i]]
            if(my_sub$title %in% names(auto_subtitles_clean)){
                value = auto_subtitles_clean[[my_sub$title]]
            } else {
                value = character(n_models)
            }

            value[my_sub$index] = my_sub$value
            auto_subtitles_clean[[my_sub$title]] = value
        }
    }


    # Formatting the names (if needed)
    alternative_names = paste0("model ", 1:n_models)
    who2replace = sapply(model_names, function(x) length(x) == 0 || x == "")
    model_names[who2replace] = alternative_names[who2replace]
    if(length(model_id) == n_models){
        index = rep(0, max(model_id))
        sub_index = rep(1, n_models)

        for(i in 1:n_models){
            id = model_id[i]
            index[id] = index[id] + 1
            sub_index[i] = index[id]
        }

        model_names = paste0("model ", model_id, ".", sub_index)
    }

    #
    # ... summary ####
    #

    check_arg(stage, "integer vector no na len(,2) GE{1} LE{2}")

    if(!missing(.vcov)){
        sysOrigin = sys.parent(.up)
        mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
        vcov_name = deparse_long(mc$.vcov)

        # We check the arguments
        check_arg(.vcov, "function", .message = "The argument '.vcov' must be a function to compute the variance-covariance matrix.")
        check_arg_plus(.vcov_args, "NULL{list()} list", .message = "The argument '.vcov_args' must be a list of arguments to be passed to the function in '.vcov'.")
    }

    IS_MULTI_CLUST = FALSE
    if(!missing(cluster) && identical(class(cluster), "list")){
        IS_MULTI_CLUST = TRUE

        check_value(cluster, "list len(value)", .value = n_models, .message = "If 'cluster' is a list, it must be of the same length as the number of models.")

        if(!missnull(se)){
            message("Argument 'se' is ignored, to use 'standard' or 'hetero' standard-errors please place these keywords in the argument 'cluster'.")
        }

        se = vector("list", n_models)
        for(m in 1:n_models){
            if(identical(cluster[[m]], "standard")){
                se[[m]] = "standard"
            } else if(identical(cluster[[m]], "hetero")){
                se[[m]] = "hetero"
            }
        }
    } else if(!missnull(se) && length(se) > 1){
        IS_MULTI_CLUST = TRUE
        check_value(se, "character vector len(value)", .value = n_models, .message = "If 'se' is a vector, it must be of the same length as the number of models.")

        cluster = vector("list", n_models)
    }

    # If se or cluster is provided, we use summary
    # if is_mult, we'll have to unroll the results
    check_mult = any(qui_iv <- sapply(all_models, function(x) isTRUE(x$iv)))
    is_mult = FALSE
    if(check_mult){
        stage = unique(stage)
        # we check if we'll have multiple results
        if(length(stage) > 1){
            # Multiple for sure
            is_mult = TRUE
        } else if(1 %in% stage && any(sapply(all_models[qui_iv], function(x) length(x$iv_endo_names) > 1))){
            # Means we have multiple end, so svl first stages
            is_mult = TRUE
        }
    }

    # The following two objects are only used if is_mult
    all_models_bis = list()
    iv_times = c()
    # we'll use iv subtitles only if length(stage) > 1
    iv_sub = c()
    for(m in 1:n_models){
        if(useSummary){
            if(missing(.vcov)){
                if(IS_MULTI_CLUST){
                    x = summary(all_models[[m]], se = se[[m]], cluster = cluster[[m]], dof = dof, stage = stage, agg = agg)
                } else {
                    x = summary(all_models[[m]], se = se, cluster = cluster, dof = dof, stage = stage, agg = agg)
                }
            } else {
                x = summary(all_models[[m]], stage = stage, .vcov = .vcov, .vcov_args = .vcov_args, vcov_name = vcov_name, agg = agg)
            }

        } else {
            # What do we do if se not provided?
            # we apply summary only to the ones that are not summaries
            x = all_models[[m]]
            if(!"cov.scaled" %in% names(x)){
                # not a summary => we apply summary to trigger default behavior
                x = summary(x, dof = dof)
            }
        }

        if(is_mult){

            if("fixest_multi" %in% class(x)){
                data = attr(x, "data")
                for(i in seq_along(data)){
                    all_models_bis[[length(all_models_bis) + 1]] = data[[i]]
                }

                iv_times[m] = length(data)
            } else {
                all_models_bis[[length(all_models_bis) + 1]] = x
                iv_times[m] = 1
            }

        } else {
            all_models[[m]] = x
        }
    }

    if(is_mult){
        # We've got mutliple estimations from summary => we expand the subtitles
        all_models = all_models_bis

        model_names = rep(model_names, iv_times)
        i_max = rep(iv_times, iv_times)
        suffix = paste0(".", unlist(lapply(iv_times, seq)))
        suffix[i_max == 1] = ""
        model_names = paste0(model_names, suffix)

        if(length(auto_subtitles) > 0){
            # We expand the subtitles
            for(i in seq_along(auto_subtitles_clean)){
                auto_subtitles_clean[[i]] = rep(auto_subtitles_clean[[i]], iv_times)
            }
        }

        if(AUTO_SUBTITLES){
            my_stages = sapply(all_models, function(x) ifelse(is.null(x$iv_stage), 3, x$iv_stage))
            if(length(unique(my_stages)) > 1){
                dict_stage = c("First", "Second", " ")
                auto_subtitles_clean[["IV stages"]] = dict_stage[my_stages]
            }
        }

        n_models = length(all_models)
    }

    if(length(auto_subtitles_clean) > 0){
        for(i in seq_along(auto_subtitles_clean)){
            my_title = names(auto_subtitles_clean)[i]
            subtitles[[dict_apply(my_title, dict)]] = dict_apply(auto_subtitles_clean[[i]], dict)
        }
    }

    # if there are subtitles
    if(length(subtitles) == 0){
        isSubtitles = FALSE
    } else {

        n_all = lengths(subtitles)
        for(i in which(n_all == 1)){
            subtitles[[i]] = rep(subtitles[[i]], n_models)
        }

        if(any(!n_all %in% c(1, n_models))){
            stop_up("If argument 'subtitles' is provided, it must be of the same length as the number of models. Current lengths: ", setdiff(n_all, c(1, n_models))[1], " vs ", n_models, " models.")
        }

        isSubtitles = TRUE

        if(isTex){
            for(i in seq_along(subtitles)){
                subtitles[[i]] = escape_latex(subtitles[[i]], up = 2)
            }
        }

    }

    #
    # We check the group and extraline arguments
    #

    if(missing(drop)) drop = NULL
    for(i in seq_along(group)){
        check_value(group[[i]], "character vector", .message = "The elements of argument 'group' must be character vectors of regular expressions.")
        drop = unique(c(drop, group[[i]]))
    }

    #
    # ... extraline ####
    #

    if(inherits(extraline, "formula")){
        # If a formula => to summon registered stats
        extraline = extraline_extractor(extraline, tex = isTex)
    }


    el_new = list() # I need it to cope with list(~f+ivf+macro, "my vars" = TRUE)
    # => the first command will create several lines
    el_names = uniquify_names(names(extraline))
    for(i in seq_along(extraline)){
        check_value(extraline[[i]], "scalar | vector(character, numeric, logical) len(value) | function | os formula",
                    .message = paste0("The elements of argument 'extraline' must be vectors of length ", n_models, ", logical scalars, functions, or one-sided formulas."),
                    .value = n_models)

        el = extraline[[i]]
        if("formula" %in% class(el)){
            el_tmp = extraline_extractor(el, el_names[i], tex = isTex)
            for(k in seq_along(el_tmp)){
                el_new[[names(el_tmp)[k]]] = el_tmp[[k]]
            }
        } else {
            if(is.null(el_names) || nchar(el_names[i]) == 0){
                stop_up("The argument 'extraline' must have names that will correspond to the row names. This is not the case for the ", n_th(i), " element.")
            }

            if(!is.function(el) && length(el) < n_models){
                # we extend
                el = rep(el, n_models)
            }

            el_new[[el_names[i]]] = el
        }
    }

    extraline = el_new

    # Now we catch the functions + normalization of the names
    el_fun_id = NULL
    if(length(extraline) > 0){
        el_fun_id = which(sapply(extraline, is.function)) # set of ID such that el[id] is a function
        if(length(el_fun_id) > 0){
            el_origin = extraline
        }

        if(length(unique(names(extraline))) != length(extraline)){
            new_names = uniquify_names(names(extraline))
            names(extraline) = new_names
        }
    }


    # end: extraline

    # we keep track of the SEs
    se_type_list = list()

    check_interaction_reorder = FALSE
    var_list <- var_reorder_list <- coef_list <- coef_below <- sd_below <- list()
    depvar_list <- obs_list <- fitstat_list <- list()
    r2_list <- aic_list <- bic_list <- loglik_list <- convergence_list <- list()
    sqCor_list = family_list = list()

    # To take care of factors
    fe_names = c()
    is_fe = vector(mode = "list", n_models)
    nb_fe = vector(mode = "list", n_models) # the number of items per factor

    slope_names = c()
    slope_flag_list = vector(mode = "list", n_models)

    if(!is.null(dict) && isTex){
        dict = escape_latex(dict, up = 2)
    }

    #
    # ... fitstat ####
    #

    if("fitstat" %in% names(opts)){
        fitstat_all = opts$fitstat
    }

    if(missing(fitstat_all)){
        # => do default
        fitstat_all = "."

    } else if(isFALSE(fitstat_all) || (length(fitstat_all) == 1 && (is.na(fitstat_all) || fitstat_all == ""))){
        fitstat_all = NULL
        drop.section = c(drop.section, "stats")

    } else if("formula" %in% class(fitstat_all)){
        check_arg(fitstat_all, "os formula", .message = "Argument 'fitstat' must be a one sided formula (or a character vector) containing valid types from the function fitstat (see details in ?fitstat).")

        fitstat_all = gsub(" ", "", strsplit(deparse_long(fitstat_all[[2]]), "+", fixed = TRUE)[[1]])

    } else {
        check_arg(fitstat_all, "character vector no na", .message = "Argument 'fitstat' must be a one sided formula (or a character vector) containing valid types from the function fitstat (see details in ?fitstat).")

    }

    if("." %in% fitstat_all){
        # Default values:
        #   - if all OLS: typical R2
        #   - if any non-OLS: pseudo R2 + squared cor.
        is_ols = sapply(all_models, function(x) x$method_type == "feols")

        if(all(is_ols)){
            if(any(sapply(all_models, function(x) "fixef_vars" %in% names(x)))){
                # means any FE model
                fitstat_default = c("r2", "wr2")
            } else {
                fitstat_default = c("r2", "ar2")
            }
        } else {
            fitstat_default = c("cor2", "pr2", "bic")
        }

        fitstat_default = c("n", fitstat_default)

        if(any(sapply(all_models, function(x) !is.null(x$theta)))){
            fitstat_default = c(fitstat_default, "theta")
        }

        fitstat_default = setdiff(fitstat_default, fitstat_all)

        if(length(fitstat_default) > 0){
            i = which(fitstat_all == ".")[1]
            if(i == length(fitstat_all)){
                fitstat_all = c(fitstat_all[0:(i-1)], fitstat_default)
            } else {
                fitstat_all = c(fitstat_all[0:(i-1)], fitstat_default, fitstat_all[(i+1):length(fitstat_all)])
            }
        }

    }

    fitstat_all = tolower(fitstat_all)

    # checking the types
    fitstat_fun_types = fitstat(give_types = TRUE)
    fitstat_type_allowed = fitstat_fun_types$types
    fitstat_all = unique(fitstat_all)
    type_alias = fitstat_fun_types$type_alias

    if(any(fitstat_all %in% names(type_alias))){

        i = intersect(fitstat_all, names(type_alias))
        fitstat_all[fitstat_all %in% i] = type_alias[fitstat_all[fitstat_all %in% i]]
    }

    pblm = setdiff(fitstat_all, fitstat_type_allowed)
    if(length(pblm) > 0){
        stop_up("Argument 'fitstat' must be a one sided formula (or a character vector) containing valid types from the function fitstat (see details in ?fitstat or use fitstat(show_types = TRUE)). The type", enumerate_items(pblm, "s.is.quote"), " not valid.")
    }

    if(isTex){
        fitstat_dict = fitstat_fun_types$tex_alias
    } else {
        fitstat_dict = fitstat_fun_types$R_alias
    }

    # Renaming the stats with user-provided aliases
    if(length(dict) > 0){
        user_stat_name = intersect(names(fitstat_dict), names(dict))
        if(length(user_stat_name) > 0){
            fitstat_dict[user_stat_name] = dict[user_stat_name]
        }

        # Now with the aliases
        user_stat_name = intersect(names(type_alias), names(dict))
        if(length(user_stat_name) > 0){
            fitstat_dict[type_alias[user_stat_name]] = dict[user_stat_name]
        }
    }

    # end: fitstat

    #
    # ... Loop ####
    #

    reformat_fitstat = FALSE
    for(m in 1:n_models){
        x = all_models[[m]]
        se_type_list[[m]] = attr(x$se, "type")

        # family
        family = x$family
        if(x$method %in% c("femlm", "feNmlm")){
            fam = switch(family, poisson = "Poisson", negbin = "Neg. Bin.", gaussian = "Gaussian", logit = "Logit")
        } else if(x$method %in% c("feglm", "feglm.fit")){
            if(family$family == "poisson" && family$link == "log"){
                fam = "Poisson"
            } else if(family$family == "binomial" && family$link == "logit"){
                fam = "Logit"
            } else if(family$family == "binomial" && family$link == "probit"){
                fam = "Probit"
            } else {
                # we try to give the greatest details ()
                fam = paste0(family$family, '("', family$link, '")')
            }
        } else if(x$method %in% c("feols", "feols.fit")){
            fam = "OLS"
        } else if(x$method == "fepois"){
            fam = "Poisson"
        } else if(x$method == "fenegbin"){
            fam = "Neg. Bin."
        }
        family_list[[m]] = fam

        # variable dependante:
        depvar <- gsub(" ", "", as.character(x$fml)[[2]])

        a <- x$coeftable
        if(!is.data.frame(a)){
            class(a) <- NULL
            a = as.data.frame(a)
        }

        # We drop the .theta coefficient
        if(x$method == "fenegbin" || (x$method %in% c("femlm", "feNmlm") && family == "negbin")){
            quiTheta = rownames(a) == ".theta"
            a = a[!quiTheta, ]
        }

        #
        # START: Formatting of the factors / FEs / Slopes
        #

        # on enleve les facteurs des variables a garder
        if(!keepFactors){
            fact = rownames(a)
            qui_drop = grepl("factor(", fact, fixed = TRUE)
            a = a[!qui_drop, , FALSE]
            b = fact[qui_drop]
            c = sapply(b, function(x) strsplit(x, "factor(", fixed=TRUE)[[1]][2])
            d = sapply(c, function(x) strsplit(x, ")", fixed=TRUE)[[1]][1])
            factor_var = unique(d)

            # Now the number of items per factor
            if(length(factor_var) == 0){
                nbItems = character(0)
            } else {
                nbItems = addCommas(sapply(factor_var, function(x) 1+sum(grepl(x, b))))
            }
        } else {
            factor_var = c()
            nbItems = character(0)
        }

        #
        # Fixed-effects
        #

        if(!is.null(x$fixef_terms)){
            terms_full = extract_fe_slope(x$fixef_terms)
            fixef_vars = terms_full$fixef_vars

            factor_var = c(factor_var, fixef_vars, recursive=TRUE)

            new_items = addCommas(as.vector(x$fixef_sizes[fixef_vars]))
            names(new_items) = fixef_vars

            nbItems = c(nbItems, new_items)
        } else if(!is.null(x$fixef_vars)){
            factor_var = c(factor_var, x$fixef_vars, recursive=TRUE)

            new_items = addCommas(as.vector(x$fixef_sizes))
            names(new_items) = names(x$fixef_sizes)

            nbItems = c(nbItems, new_items)
        }

        nb_fe[[m]] = nbItems

        # Formatting

        lFactor = rep(yesNo[1], length(factor_var))
        names(lFactor) = factor_var
        is_fe[[m]] = lFactor

        fe_names = unique(c(fe_names, factor_var, recursive=TRUE))

        #
        # SLOPES
        #

        if(!is.null(x$fixef_terms)){
            terms_full = extract_fe_slope(x$fixef_terms)
            slope_fe = terms_full$slope_fe
            slope_vars = terms_full$slope_vars

            # we change the names right away
            slope_fe_name = slope_fe
            slope_vars_name = slope_vars
            if(!is.null(dict)){
                qui = which(slope_fe %in% names(dict))
                if(length(qui) > 0){
                    slope_fe_name[qui] = dict[slope_fe[qui]]
                }

                qui = which(slope_vars %in% names(dict))
                if(length(qui) > 0){
                    slope_vars_name[qui] = dict[slope_vars[qui]]
                }
            }

            if(nchar(style$slopes.format) > 0){
                slope_format = style$slopes.format
                slope_var_full = c()
                for(i in seq_along(slope_vars_name)){
                    slope_var_full[i] = gsub("__slope__", slope_fe_name[i], gsub("__var__", slope_vars_name[i], slope_format, fixed = TRUE), fixed = TRUE)
                }

            } else {
                slope_var_full = paste0(slope_vars_name, " (", slope_fe_name, ")")
            }


        } else {
            slope_var_full = c()
        }

        slope_flag = rep(yesNo[1], length(slope_var_full))
        names(slope_flag) = slope_var_full
        slope_flag_list[[m]] = slope_flag

        slope_names = unique(c(slope_names, slope_var_full, recursive = TRUE))

        #
        #   END: FE/slope formatting
        #

        # on enleve les espaces dans les noms de variables
        var <- var_origin <- c(gsub(" ", "", row.names(a)))
        # renaming => Tex only
        if(isTex || !is.null(dict)){
            qui = var %in% names(dict)
            var[qui] = dict[var[qui]]
            tv = table(var)
            if(any(tv > 1)){
                value_pblm = names(tv)[tv > 1][1]
                var_pblm = c(gsub(" ", "", row.names(a)))[var == value_pblm]
                stop("Problematic value for argument 'dict': The variables ", enumerate_items(var_pblm, "quote"), " have all the same alias ('", value_pblm, "') in the same estimation. This is not supported, please provide a separate alias for each.")
            }

            # if there are still interactions, we rename them
            new_var = var
            var_left = var[!qui]
            if(length(var_left) > 0 && any(grepl(":|\\*", var_left))){
                check_interaction_reorder = TRUE


                qui_inter = grepl("(?<=[^:\\*])(::|:|\\*)(?=[^:\\*])", var_left, perl = TRUE)
                inter = strsplit(var_left[qui_inter], "(?<=[^:\\*])(:|\\*)(?=[^:\\*])", perl = TRUE)

                fun_rename = function(x){
                    # We put the factors on the right

                    qui_factor = grepl("::", x)
                    if(any(qui_factor)){
                        res = x
                        x_split = strsplit(res, "::")

                        # There may be svl factors
                        for(i in which(qui_factor)){
                            if(res[i] %in% names(dict)){
                                res[i] = dict[res[i]]
                            } else {
                                value_split = dict_apply(x_split[[i]], dict)

                                if(isTex){
                                    res[i] = paste(value_split, collapse = " $=$ ")
                                } else {
                                    res[i] = paste(value_split, collapse = " = ")
                                }
                            }
                        }

                        # We put the factors on the right
                        res = res[base::order(qui_factor)]

                    } else {
                        res = x
                    }

                    who = res %in% names(dict)
                    res[who] = dict[res[who]]

                    if(isTex){
                        res = paste0(res, collapse = interaction.combine)
                    } else {
                        res = paste0(res, collapse = " x ")
                    }

                    return(res)
                }

                inter_named = sapply(inter, fun_rename)
                new_inter = sapply(inter, function(x) fun_rename(sort(x)))

                var[!qui][qui_inter] = inter_named
                new_var[!qui][qui_inter] = new_inter

            }

            # IV: fit_
            if(any(qui_iv <- grepl("^fit_", var_left))){
                # IVs
                iv_vars = gsub("^fit_", "", var_left[qui_iv])

                if(fit_format == "__var__"){
                    iv_vars = dict_apply(iv_vars, dict)

                } else {
                    for(i in seq_along(iv_vars)){
                        iv_vars = gsub("__var__", dict_apply(iv_vars[i], dict), fit_format)
                    }
                }

                var[!qui][qui_iv] = new_var[!qui][qui_iv] = iv_vars
            }

            # We take care of poly
            # DEPRECATED: this should not be done!!!!!!!!!!!
            # Poly is no regular polynomial but orthogonal polynomial!!!!!
            # This will imply difference in estimates with raw powers, we should not conflate them

            # We take care of raw powers
            if(any(grepl("^I\\([^\\^]+\\^[[:digit:]]+\\)", var_left))){

                # We clean only I(var^d)
                qui_pow = grepl("^I\\([^\\^]+\\^[[:digit:]]+\\)$", var_left)
                if(any(qui_pow)){
                    pow_var = gsub("^I\\(|\\^.+$", "", var_left[qui_pow])
                    pow_digit = as.numeric(gsub(".+\\^|\\)", "", var_left[qui_pow]))

                    pow_digit_clean = poly_dict[pow_digit]
                    if(isTex){
                        pow_digit_clean[is.na(pow_digit_clean)] = paste0("$^{", pow_digit[is.na(pow_digit_clean)], "}")
                    } else {
                        pow_digit_clean[is.na(pow_digit_clean)] = paste0(" ^ ", pow_digit[is.na(pow_digit_clean)])
                    }

                    pow_named = paste0(dict_apply(pow_var, dict), pow_digit_clean)

                    var[!qui][qui_pow] = new_var[!qui][qui_pow] = pow_named
                }
            }

            names(new_var) = names(var) = var_origin
            var_reorder_list[[m]] <- new_var

        } else {
            # We reorder the interaction terms alphabetically
            new_var = var
            qui = grepl("(?<=[^:]):(?=[^:])", new_var, perl = TRUE)
            if(any(qui)){
                check_interaction_reorder = TRUE
                inter = strsplit(new_var[qui], "(?<=[^:]):(?=[^:])", perl = TRUE)
                new_inter = sapply(inter, function(x) paste0(sort(x), collapse = ":"))
                new_var[qui] = new_inter
            }

            names(new_var) = names(var) = var_origin

            var_reorder_list[[m]] <- new_var
        }

        coef = fun_format(a[, 1])

        if(coefstat == "se"){
            se_value = fun_format(a[, 2])

        } else if(coefstat == "tstat"){
            se_value = fun_format(a[, 3])

        } else if(coefstat == "confint"){
            se_value = apply(confint(x, level = ci), 1, function(z) paste0("[", fun_format(z[1]), "; ", fun_format(z[2]), "]"))

        } else {
            stop("Wrong value for the argument 'coefstat': ", coefstat, " is not supported.")
        }

        if(add_signif){
            if(isTex){
                pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(tex_star(names(signifCode)), ""))
            } else {
                pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
            }
            pval[is.na(pval)] = ""
        } else {
            pval = rep("", nrow(a))
        }


        # If the coefficient is bounded, we suppress the 'stars'
        isBounded = grepl("bounded", se_value)
        if(any(isBounded)){
            pval[isBounded] = ""
        }

        if(coefstat != "confint"){
            structured_coef = c(paste0(coef, pval, " (", se_value, ")"))
        } else {
            structured_coef = c(paste0(coef, pval, " ", se_value))
        }

        # saving the infos
        var_list[[m]] <- var
        names(structured_coef) <- var
        coef_list[[m]] <- structured_coef
        if(sdBelow){
            cb = c(paste0(coef, pval))
            if(coefstat != "confint"){
                sb = c(paste0("(", se_value, ")"))
            } else {
                sb = se_value
            }

            names(cb) = names(sb) = var
            coef_below[[m]] = cb
            sd_below[[m]] = sb
        }

        # La depvar
        depvar_list[[m]] <- depvar

        #
        #  Fit statistics
        #

        n = nobs(x)
        obs_list[[m]] = n
        convergence_list[[m]] = ifelse(is.null(x$convStatus), TRUE, x$convStatus)

        K = x$nparams

        if(length(fitstat_all) == 0){
            fitstat_list[[m]] = NA
        } else {

            my_stats = lapply(fitstat(x, fitstat_all, etable = TRUE, verbose = FALSE), fun_format_stats)

            if(any(grepl("::", names(my_stats), fixed = TRUE))){
                reformat_fitstat = TRUE
            }

            fitstat_list[[m]] = my_stats
        }

        #
        # Extraline, when function
        #

        for(i in el_fun_id){
            f = el_origin[[i]]

            if(m == 1){
                extraline[[i]] = f(x)

            } else {
                extraline[[i]] = c(extraline[[i]], f(x))
            }
        }


    }

    if(check_interaction_reorder){
        if(length(unique(unlist(var_reorder_list))) < length(unique(unlist(var_list)))){
            var_list = var_reorder_list
            for(m in 1:length(var_list)){
                names(coef_list[[m]]) = var_list[[m]]
                if(sdBelow){
                    names(coef_below[[m]]) = var_list[[m]]
                    names(sd_below[[m]]) = var_list[[m]]
                }

            }
        }
    }

    #
    # Fitstat reformating
    #

    if(length(fitstat_all) > 0){

        if(reformat_fitstat){
            # PAIN IN THE NECK!

            all_names = unique(unlist(lapply(fitstat_list, names)))
            all_names_new = c()
            for(v in fitstat_all){
                all_names_new = c(all_names_new, all_names[gsub("::.+", "", all_names) == v])
            }

            fitstat_list_new = list()
            for(i in seq_along(fitstat_list)){
                my_list = fitstat_list[[i]]
                fitstat_list_new[[i]] = lapply(all_names_new, function(x) ifelse(is.null(my_list[[x]]), NA, my_list[[x]]))
            }

            # Now the aliases
            fitstat_dict_new = c()
            fun_rename = function(x) {
                if(!grepl("::", x, fixed = TRUE)) return(fitstat_dict[x])

                xx = strsplit(x, "::")[[1]]
                paste0(fitstat_dict[xx[1]], ", ", dict_apply(xx[2], dict))
            }

            for(v in all_names_new){
                fitstat_dict_new[v] = fun_rename(v)
            }

            fitstat_list = fitstat_list_new
            attr(fitstat_list, "format_names") = fitstat_dict_new

        } else {
            attr(fitstat_list, "format_names") = fitstat_dict[fitstat_all]
        }
    }

    if(isTex){
        if(missing(title)){
            title = "no title"
        } else {
            title = escape_latex(title, 2)
        }
    } else {
        if(missing(title)){
            title = NULL
        }
    }


    if((missing(convergence) && any(convergence_list == FALSE)) || (!missing(convergence) && convergence)){
        convergence = TRUE
    } else {
        convergence = FALSE
    }

    if(isTRUE(show_family) || (is.null(show_family) && length(unique(family_list)) > 1)){
        family = TRUE
    } else {
        family = FALSE
    }

    if((is.null(show_depvar) && length(unique(unlist(depvar_list))) > 1) || isTRUE(show_depvar)){
        depvar = TRUE
    } else {
        depvar = FALSE
    }

    if(missing(keep)) keep = NULL
    if(missing(order)) order = NULL
    if(missing(file)) file = NULL
    if(missing(label)) label = NULL

    if(IS_FIXEF_GROUP && length(fe_names) > 0){

        # Table of presence of fixed-effects
        is_fe_mat = c()
        for(i in 1:length(is_fe)){
            # The order of the FEs needs not be the same since is_FE is
            # sequentially constructed (hence the following construct)
            is_fe_mat = cbind(is_fe_mat, 1 * (!is.na(is_fe[[i]][fe_names])))
        }

        # row => FEs
        # column => models

        if(isTRUE(fixef.group) && length(fe_names) > 1){

            fe_groups = list()
            n_fe = nrow(is_fe_mat)
            g = 1 # => group
            i_left = 1:n_fe
            while(length(i_left) > 1){

                i = i_left[1]
                i_left = i_left[-1]
                any_done = FALSE
                for(j in i_left){

                    # we check if same pattern
                    if(all(is_fe_mat[i, ] == is_fe_mat[j, ])){
                        any_done = TRUE
                        if(length(fe_groups) < g){
                            # we initialize it
                            fe_groups[[g]] = c(i, j)
                        } else {
                            # we add j
                            fe_groups[[g]] = c(fe_groups[[g]], j)
                        }

                        i_left = setdiff(i_left, j)
                    }
                }

                if(any_done){
                    g = g + 1
                }
            }

            # We simply group the fixed-effects + we apply dict
            fe_names_origin = fe_names
            for(my_group in fe_groups){
                # We drop the old
                group_fe = fe_names_origin[my_group]
                fe_names = setdiff(fe_names, group_fe)

                # and add the new
                group_fe_new = rename_fe(group_fe, dict)
                new_fe = enumerate_items(group_fe_new)

                fe_names = c(new_fe, fe_names)

                # Now we update is_fe
                for(i in seq_along(is_fe)){
                    current_fe = is_fe[[i]]
                    if(all(group_fe %in% names(current_fe))){
                        current_fe[new_fe] = yesNo[1]
                        is_fe[[i]] = current_fe
                    }
                }
            }


        } else {
            # Custom grouping from the user

            # Behavior:
            # The values in the row:
            # - if all the row's FEs are in the model => TRUE
            # - if some of the row's FEs are in the model => partial
            # - if none is in the model => FALSE
            #
            # The regular fixed-effects:
            # - FEs of a group are removed if they are either always fully present or fully absent
            # - if an FE from a group is partially present at least once
            #   * it will stay as regular FE
            #   * this avoids messing up from the user side
            #

            fe_names_full = rename_fe(fe_names, dict)
            names(fe_names_full) = fe_names

            is_inconsistent = rep(FALSE, length(fixef.group))
            fixef.extralines = list()
            fe2remove = c()
            for(i in seq_along(fixef.group)){
                my_group = keep_apply(fe_names_full, fixef.group[[i]])

                if(length(my_group) == 0){
                    # not any FE found, even partial => we skip the row
                    # we skip the row
                    next
                }

                my_group_origin = names(my_group)
                is_there = c()
                for(m in 1:length(is_fe)){
                    # We check full and partial presence

                    fe_model = names(is_fe[[m]])
                    if(all(my_group_origin %in% fe_model)){
                        is_there[m] = TRUE

                    } else if(any(my_group_origin %in% fe_model)){
                        # partial presence
                        is_there[m] = NA
                        is_inconsistent[i] = TRUE

                    } else {
                        is_there[m] = FALSE
                    }

                }

                ok_removal = !any(is.na(is_there))
                if(ok_removal){
                    fe2remove = c(fe2remove, my_group_origin)
                }

                is_there = yesNo[2 - is_there]
                is_there[is.na(is_there)] = "partial"

                # Now we create the extra line
                el_name = names(fixef.group)[i]
                if(grepl("^(\\^|_|-)", el_name)){
                    # => the user specifies the placement => replaces the default

                    # If one ^ only => we replace it with _^ to override default of extraline
                    if(!grepl("^\\^(_|-|\\^)", el_name)){
                        el_name = paste0("_", el_name)
                    }

                } else {
                    # we add the default placement
                    el_name = paste0("^-", el_name)
                }

                fixef.extralines[[el_name]] = is_there

            }

            if(length(fixef.extralines) > 0){
                # We add it to extra lines
                extraline[names(fixef.extralines)] = fixef.extralines
            }

            if(length(fe2remove) > 0){
                fe_names = setdiff(fe_names, fe2remove)
            }

            # Warning for inconsistencies
            if(any(is_inconsistent)){
                i = which(is_inconsistent)[1]
                if(length(is_inconsistent) == 1){
                    msg = "the group leads to an inconsistent row (defined by "
                } else if(all(is_inconsistent)){
                    msg = "all groups lead to inconsistent rows (e.g. the one defined by "
                } else {
                    msg = paste0("some groups lead to inconsistent rows (the ", n_th(i), " one defined by ")
                }
                warn_up("In 'fixef.group', ", msg, deparse_long(fixef.group[[i]]), ").\nTo create inconsistent rows: use drop.section = 'fixef' combined with the arghument 'extraline'.")
            }

        }

    }

    res = list(se_type_list=se_type_list, var_list=var_list, coef_list=coef_list, coef_below=coef_below, sd_below=sd_below, depvar_list=depvar_list, obs_list=obs_list, convergence_list=convergence_list, fe_names=fe_names, is_fe=is_fe, nb_fe=nb_fe, slope_flag_list = slope_flag_list, slope_names=slope_names, useSummary=useSummary, model_names=model_names, family_list=family_list, fitstat_list=fitstat_list, subtitles=subtitles, isSubtitles=isSubtitles, title=title, convergence=convergence, family=family, keep=keep, drop=drop, order=order, file=file, label=label, sdBelow=sdBelow, signifCode=signifCode, fixef_sizes=fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, depvar=depvar, useSummary=useSummary, dict=dict, yesNo=yesNo, add_signif=add_signif, float=float, coefstat=coefstat, ci=ci, style=style, notes=notes, group=group, extraline=extraline, placement=placement, drop.section=drop.section, tex_tag=tex_tag, fun_format = fun_format, coef.just = coef.just)

    return(res)
}

etable_internal_latex = function(info){
    # Internal function to display the latex table

    n_models = length(info$depvar_list)
    # Getting the information
    se_type_list = info$se_type_list
    var_list = info$var_list
    coef_list = info$coef_list
    coef_below = info$coef_below
    sd_below = info$sd_below
    depvar_list = info$depvar_list
    depvar = info$depvar # flag for including the depvar
    obs_list = info$obs_list
    convergence_list = info$convergence_list
    fe_names = info$fe_names
    is_fe = info$is_fe
    nb_fe = info$nb_fe
    slope_names = info$slope_names
    slope_flag_list = info$slope_flag_list
    family_list = info$family_list
    fitstat_list = info$fitstat_list
    subtitles = info$subtitles
    isSubtitles = info$isSubtitles
    title = info$title
    label = info$label
    keep = info$keep
    drop = info$drop
    order = info$order
    file = info$file
    family = info$family
    convergence = info$convergence
    sdBelow = info$sdBelow
    signifCode = info$signifCode
    fixef_sizes = info$fixef_sizes
    fixef_sizes.simplify = info$fixef_sizes.simplify
    dict = info$dict
    yesNo = info$yesNo
    add_signif = info$add_signif
    float = info$float
    coefstat = info$coefstat
    ci = info$ci
    style = info$style
    notes = info$notes
    group = info$group
    extraline = info$extraline
    placement = info$placement
    drop.section = info$drop.section
    tex_tag = info$tex_tag
    fun_format = info$fun_format

    # Formatting the searating lines
    if(nchar(style$line.top) > 1) style$line.top = paste0(style$line.top, "\n")
    if(nchar(style$line.bottom) > 1) style$line.bottom = paste0(style$line.bottom, "\n")

    #
    # prompting the infos gathered
    #

    # Starting the table
    myTitle = title
    if(!is.null(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
    if(float){
        if(nchar(placement) > 0) placement = paste0("[", placement, "]")
        start_table = paste0("\\begin{table}", placement, "\n\\centering\n\\caption{",  myTitle, "}\n")
        end_table = "\\end{table}"
    } else {
        start_table = ""
        end_table = ""
    }


    # intro and outro Latex tabular
    # \begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lc}
    if(style$tabular == "normal"){
        intro_latex <- paste0("\\begin{tabular}{l", paste0(rep("c", n_models), collapse=""), "}\n", style$line.top)
        outro_latex <- "\\end{tabular}\n"
    } else if(style$tabular == "*"){
        intro_latex <- paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}l", paste0(rep("c", n_models), collapse=""), "}\n", style$line.top)
        outro_latex <- "\\end{tabular*}\n"
    } else if(style$tabular == "X"){
        intro_latex <- paste0("\\begin{tabularx}{\\textwidth}{", paste0(rep("X", n_models + 1), collapse = ""), "}\n", style$line.top)
        outro_latex <- "\\end{tabularx}\n"
    }

    # 1st lines => dep vars
    depvar_list = dict_apply(c(depvar_list, recursive = TRUE), dict)
    depvar_list = escape_latex(depvar_list, up = 2)

    # We write the dependent variables properly, with multicolumn when necessary
    # to do that, we count the number of occurrences of each variable (& we respect the order provided by the user)
    nb_multi = 1
    names_multi = depvar_list[1]

    if(n_models > 1){
        k = 1
        old_dep = depvar_list[1]
        for(i in 2:length(depvar_list)){
            if(depvar_list[i] == old_dep){
                nb_multi[k] = nb_multi[k] + 1
            } else {
                k = k + 1
                nb_multi[k] = 1
                names_multi[k] = old_dep = depvar_list[i]
            }
        }
    }

    # now the proper format
    first_line = escape_latex(style$depvar.title)
    if(length(nb_multi) == 1){
        first_line = gsub("(s)", "", first_line, fixed = TRUE)
    } else {
        first_line = gsub("(s)", "s", first_line, fixed = TRUE)
    }

    for(i in 1:length(nb_multi)){
        if(nb_multi[i] == 1){
            # no multi column
            first_line = paste0(first_line, "&", names_multi[i])
        } else {
            first_line = paste0(first_line, "&\\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
        }
    }
    first_line = paste0(first_line, "\\\\\n")

    if(!depvar) first_line = NULL

    # Model line
    if(nchar(style$model.format) > 0){
        m = style$model.format
        if(grepl("1", m, fixed = TRUE)){
            info_split = strsplit(m, "1", fixed = TRUE)[[1]]
            m_info = 1:n_models

        } else if(grepl("i", m, fixed = TRUE)){
            info_split = strsplit(m, "i", fixed = TRUE)[[1]]
            unit = "i.ii.iii.iv.v.vi.vii.viii.ix.x.xi.xii.xiii.xiv.xv.xvi.xvii.xviii.xix.xx"
            m_info = strsplit(unit, ".", fixed = TRUE)[[1]][1:n_models]

        } else if(grepl("I", m, fixed = TRUE)){
            info_split = strsplit(m, "I", fixed = TRUE)[[1]]
            unit = "I.II.III.IV.V.VI.VII.VIII.IX.X.XI.XII.XIII.XIV.XV.XVI.XVII.XVIII.XIX.XX"
            m_info = strsplit(unit, ".", fixed = TRUE)[[1]][1:n_models]

        } else if(grepl("a", m, fixed = TRUE)){
            info_split = strsplit(m, "a", fixed = TRUE)[[1]]
            m_info = letters[1:n_models]

        } else if(grepl("A", m, fixed = TRUE)){
            info_split = strsplit(m, "A", fixed = TRUE)[[1]]
            m_info = LETTERS[1:n_models]

        }

        right = ifelse(length(info_split) > 1, info_split[2], "")
        model_format = paste0(info_split[1], m_info, right)

    } else {
        model_format = paste0("(", 1:n_models, ")")
    }
    model_line = paste0(escape_latex(style$model.title), "&", paste0(model_format, collapse = " & "), "\\\\\n")

    # a simple line with only "variables" written in the first cell
    if(nchar(style$var.title) == 0){
        coef_title = ""
    } else if(style$var.title == "\\midrule"){
        coef_title = "\\midrule "
    } else {
        coef_title = paste0(escape_latex(style$var.title), "& ", paste(rep(" ", n_models), collapse = " & "), "\\\\\n")
    }

    # Coefficients, the tricky part
    coef_lines <- list()

    # we need to loop not to lose names
    all_vars = c()
    for(vars in var_list){
        all_vars = c(all_vars, vars[!vars %in% all_vars])
    }
    # all_vars <- unique(unlist(var_list))

    for(i in seq_along(group)){
        gi = group[[i]]
        present = sapply(var_list, function(x) any(keep_apply(x, gi, TRUE)))

        group[[i]] = present
    }

    # keeping some coefs
    all_vars = keep_apply(all_vars, keep)

    # dropping some coefs
    all_vars = drop_apply(all_vars, drop)

    # ordering the coefs
    all_vars = order_apply(all_vars, order)

    if(length(all_vars) == 0) stop_up("Not any variable was selected, please reframe your keep/drop arguments.")

    # changing the names of the coefs
    aliasVars = all_vars
    names(aliasVars) = all_vars

    qui = all_vars %in% names(dict)
    who = aliasVars[qui]
    aliasVars[qui] = dict[who]
    aliasVars = escape_latex(aliasVars, up = 2)

    coef_mat <- all_vars
    for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
    coef_mat[is.na(coef_mat)] <- "  "
    if(sdBelow){
        coef_lines = c()
        for(v in all_vars){
            myCoef = mySd = myLine = c()
            for(m in 1:n_models){
                myCoef = c(myCoef, coef_below[[m]][v])
                mySd = c(mySd, sd_below[[m]][v])
            }

            myCoef[is.na(myCoef)] = "  "
            mySd[is.na(mySd)] = "  "
            myCoef = paste0(aliasVars[v], " & ", paste0(myCoef, collapse = " & "))
            mySd = paste0("  &", paste0(mySd, collapse = " & "))
            myLines = paste0(myCoef, "\\\\\n", mySd, "\\\\\n")
            coef_lines = c(coef_lines, myLines)
        }
        coef_lines = paste0(coef_lines, collapse="")
    } else {
        coef_mat[, 1] = aliasVars
        coef_lines = paste0(paste0(apply(coef_mat, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")
    }

    #
    # Fixed-effects (if needed)
    #

    nb_FE_lines = ""
    if(length(fe_names) > 0 && !"fixef" %in% drop.section){

        if(nchar(style$fixef.title) == 0){
            fixef_title = ""
        } else if(style$fixef.title == "\\midrule"){
            fixef_title = "\\midrule "
        } else {
            fixef_title = paste0(escape_latex(style$fixef.title), "& ", paste(rep(" ", n_models), collapse = " & "), "\\\\\n")
        }

        # The number of FEs
        for(m in 1:n_models) {
            quoi = is_fe[[m]][fe_names]
            quoi[is.na(quoi)] = yesNo[2]
            is_fe[[m]] = quoi

            # We do the same for the number of items
            quoi = nb_fe[[m]][fe_names]
            quoi[is.na(quoi)] = "--"
            nb_fe[[m]] = quoi
        }

        all_fe = matrix(c(is_fe, recursive = TRUE), nrow = length(fe_names))

        # We change the names of the FEs
        fe_names = rename_fe(fe_names, dict)
        fe_names_raw = escape_latex(fe_names)

        fe_names = paste0(style$fixef.prefix, fe_names_raw, style$fixef.suffix)

        # The fixed-effects sizes

        if(fixef_sizes){
            # For the number of items
            all_nb_FEs = matrix(c(nb_fe, recursive=TRUE), nrow = length(fe_names))

            is_complex = rep(TRUE, nrow(all_nb_FEs))
            if(fixef_sizes.simplify){
                check_simplified = function(x) length(unique(x[grepl("[[:digit:]]", x)])) == 1
                is_simple = apply(all_nb_FEs, 1, check_simplified)
                for(i in which(is_simple)){
                    x = all_nb_FEs[i, ]
                    nb = x[grepl("[[:digit:]]", x)][1]
                    fe_names[i] = paste0(fe_names[i], " (", nb, ")")
                }
                is_complex = !is_simple
            }

            if(any(is_complex)){
                fe_names_nbItems = paste0(style$fixef_sizes.prefix, fe_names_raw[is_complex], style$fixef_sizes.suffix)
                all_nb_FEs = cbind(fe_names_nbItems, all_nb_FEs[is_complex, , drop = FALSE])
                nb_FE_lines <- paste0(paste0(apply(all_nb_FEs, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")
            }

        }

        all_fe = cbind(fe_names, all_fe)
        FE_lines <- paste0(paste0(apply(all_fe, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")

    } else {
        FE_lines = NULL
        fixef_title = NULL
    }

    #
    # Slopes (if needed)
    #

    if(length(slope_names) > 0 &&  !"slopes" %in% drop.section){

        if(nchar(style$slopes.title) == 0){
            slope_intro = ""
        } else if(style$slopes.title == "\\midrule"){
            slope_intro = "\\midrule "
        } else {
            slope_intro = paste0(escape_latex(style$slopes.title), "& ", paste(rep(" ", n_models), collapse = " & "), "\\\\\n")
        }

        # reformat the yes/no slope
        for(m in 1:n_models) {
            quoi = slope_flag_list[[m]][slope_names]
            quoi[is.na(quoi)] = yesNo[2]
            slope_flag_list[[m]] = quoi
        }

        # Changing the slope names
        qui = slope_names %in% names(dict)
        who = slope_names[qui]
        slope_names[qui] = dict[who]
        slope_names = escape_latex(slope_names)

        # Matrix with yes/no information
        all_slopes = matrix(c(slope_flag_list, recursive = TRUE), nrow = length(slope_names))
        all_slopes = cbind(slope_names, all_slopes)
        slope_lines <- paste0(paste0(apply(all_slopes, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")

    } else {
        slope_intro = NULL
        slope_lines = NULL
    }

    # Subtitles
    if(isSubtitles){
        # info_subtitles = paste0("  & ", paste(subtitles, collapse = " & "), "\\\\\n")
        n_sub = length(subtitles)
        sub_names = names(subtitles)
        if(is.null(sub_names)){
            sub_names = character(n_sub)
        }

        info_subtitles = character(n_sub)
        for(i in 1:n_sub){
            info_subtitles[[i]] = paste0(sub_names[i], "  & ", tex_multicol(subtitles[[i]]), "\\\\\n")
        }
        info_subtitles = paste(info_subtitles, collapse = "")

    } else {
        info_subtitles = ""
    }

    # Convergence information
    info_convergence = ""
    if(convergence){
        info_convergence = paste0("Convergence &", paste(convergence_list, collapse = " & "), "\\\\\n")
    }

    # information on family
    if(family){
        info_family <- paste0(" &  ", paste(family_list, collapse = " & "), "\\\\\n")
    } else {
        info_family = ""
    }


    # The standard errors => if tablefoot = TRUE
    info_SD = ""
    info_muli_se = ""

    # We go through this in order to get the multi se if needed
    isUniqueSD = length(unique(unlist(se_type_list))) == 1
    nb_col = length(obs_list) + 1
    sd_intro = paste0("\\multicolumn{", nb_col, "}{l}{\\emph{")

    if(isUniqueSD){
        my_se = unique(unlist(se_type_list)) # it comes from summary
        # every model has the same type of SE
        if(my_se == "Standard") my_se = "Normal"

        # Now we modify the names of the clusters if needed
        my_se = format_se_type_latex(my_se, dict)

        if(coefstat == "se"){
            coefstat_sentence = " standard-errors in parentheses"
        } else if(coefstat == "tstat"){
            coefstat_sentence = " co-variance matrix, t-stats in parentheses"
        } else {
            coefstat_sentence = paste0(" co-variance matrix, ", round(ci*100), "\\% confidence intervals in brackets")
        }
        info_SD = paste0(escape_latex(style$tablefoot.title), sd_intro, my_se, coefstat_sentence, "}}\\\\\n")

        if(add_signif){
            info_SD = paste0(info_SD, sd_intro, "Signif. Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
        }

        info_muli_se = ""

    } else {
        all_se_type = sapply(se_type_list, format_se_type_latex, dict = dict, inline = TRUE)

        if(coefstat == "se"){
            coefstat_sentence = "Standard-Errors"
        } else {
            coefstat_sentence = "Co-variance"
        }

        info_muli_se = paste0(coefstat_sentence, "& ", paste(all_se_type, collapse = "&"), "\\\\\n")

        if(add_signif){
            info_SD = paste0(escape_latex(style$tablefoot.title), sd_intro, "Signif. Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
        } else {
            myAmpLine = paste0(paste0(rep(" ", length(depvar_list)+1), collapse = " & "), "\\tabularnewline\n")
            info_SD = paste0(escape_latex(style$tablefoot.title), myAmpLine, "\\\\\n")
        }
    }

    if(style$tablefoot){
        if(!identical(style$tablefoot.value, "default")){
            value = style$tablefoot.value
            if(isUniqueSD){
                # my_se: computed right before
                value = gsub("__se_type__", my_se, value)
            }

            info_SD = paste0(escape_latex(style$tablefoot.title), paste(sd_intro, value, "}}\\\\\n", collapse = ""))
        }
    } else {
        info_SD = ""
    }


    # Information on number of items
    supplemental_info = ""

    #
    # Fit statistics
    #

    if(!"stats" %in% drop.section){
        if(nchar(style$stats.title) == 0){
            stat_title = ""
        } else if(style$stats.title == "\\midrule"){
            stat_title = "\\midrule "
        } else {
            stat_title = paste0(escape_latex(style$stats.title), "& ", paste(rep(" ", n_models), collapse = "&"), "\\\\\n")
        }

        stat_lines = paste0(nb_FE_lines, info_convergence, info_muli_se)

        if(!all(sapply(fitstat_list, function(x) all(is.na(x))))){

            fit_names = attr(fitstat_list, "format_names")
            nb = length(fit_names)
            for(fit_id in 1:nb){
                fit = sapply(fitstat_list, function(x) x[[fit_id]])
                if(all(is.na(fit))) next
                fit[is.na(fit)] = ""
                stat_lines = paste0(stat_lines, fit_names[fit_id], " & ", paste0(fit, collapse = "&"), "\\\\\n")
            }
        }
    } else {
        stat_title = stat_lines = ""
    }


    # Notes
    info_notes = ""
    if(nchar(notes) > 0){
        info_notes = paste0("\n", escape_latex(style$notes.title), notes, "\n")
    }

    #
    # Group
    #

    create.fixef_title = FALSE

    for(i in seq_along(group)){
        gi = group[[i]]
        gi_format = yesNo[2 - gi]

        gi_name = names(group)[i]

        gi_full = ""
        gi_where = "coef"

        gi_top = FALSE
        if(grepl("^(\\^|_|-)", gi_name)){
            # sec: section

            row = substr(gi_name, 1, 1)
            gi_name = substr(gi_name, 2, nchar(gi_name))
            if(row == "-"){
                # FE section
                sec = "fixef"

            } else if(grepl("^(\\^|_|-)", gi_name)){
                sec = substr(gi_name, 1, 1)
                gi_name = substr(gi_name, 2, nchar(gi_name))

            } else {
                if(row == "^"){
                    # implicit location
                    sec = "coef"
                } else {
                    sec = row
                    row = "_"
                }
            }

            gi_top = row == "^"
            gi_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat", sec)
        }

        gi_full = paste0(gi_full, gi_name, " & ", paste0(gi_format, collapse = " & "), "\\\\\n")

        if(gi_top){
            if(gi_where == "coef"){
                coef_lines = c(gi_full, coef_lines)
            } else if(gi_where == "fixef"){
                create.fixef_title = TRUE
                FE_lines = c(gi_full, FE_lines)
            } else {
                stat_lines = c(gi_full, stat_lines)
            }
        } else {
            if(gi_where == "coef"){
                coef_lines = c(coef_lines, gi_full)
            } else if(gi_where == "fixef"){
                create.fixef_title = TRUE
                FE_lines = c(FE_lines, gi_full)
            } else {
                stat_lines = c(stat_lines, gi_full)
            }
        }
    }

    #
    # Extra lines
    #

    for(i in seq_along(extraline)){
        el = extraline[[i]]

        # The format depends on the type
        if(is.logical(el)){
            if(length(el) == 1){
                el_format = rep(yesNo[2 - el], n_models)
            } else {
                el_format = yesNo[2 - el]
            }
        } else if(is.numeric(el)){
            el_format = fun_format(el)
        } else {
            el_format = el
        }

        el_format[is.na(el_format)] = ""

        el_name = names(extraline)[i]

        el_full = ""
        el_where = "coef"

        el_top = FALSE
        if(grepl("^(\\^|_|-)", el_name)){
            # sec: section

            row = substr(el_name, 1, 1)
            el_name = substr(el_name, 2, nchar(el_name))
            if(row == "-"){
                # FE section
                sec = "fixef"

            } else if(grepl("^(\\^|_|-)", el_name)){
                sec = substr(el_name, 1, 1)
                el_name = substr(el_name, 2, nchar(el_name))

            } else {
                if(row == "^"){
                    # implicit location
                    sec = "coef"
                } else {
                    sec = row
                    row = "_"
                }
            }

            el_top = row == "^"
            el_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat", sec)
        }

        el_full = paste0(el_full, el_name, " & ", paste0(el_format, collapse = " & "), "\\\\\n")

        if(el_top){
            if(el_where == "coef"){
                coef_lines = c(el_full, coef_lines)
            } else if(el_where == "fixef"){
                create.fixef_title = TRUE
                FE_lines = c(el_full, FE_lines)
            } else {
                stat_lines = c(el_full, stat_lines)
            }
        } else {
            if(el_where == "coef"){
                coef_lines = c(coef_lines, el_full)
            } else if(el_where == "fixef"){
                create.fixef_title = TRUE
                FE_lines = c(FE_lines, el_full)
            } else {
                stat_lines = c(stat_lines, el_full)
            }
        }

    }

    if(create.fixef_title && is.null(fixef_title) && !"fixef" %in% drop.section){
        if(nchar(style$fixef.title) == 0){
            fixef_title = ""
        } else if(style$fixef.title == "\\midrule"){
            fixef_title = "\\midrule "
        } else {
            fixef_title = paste0(escape_latex(style$fixef.title), "& ", paste(rep(" ", n_models), collapse = " & "), "\\\\\n")
        }
    }

    # Stacking var and stat
    coef_stack = c(coef_title, coef_lines)
    stat_stack = c(stat_title, stat_lines)
    fixef_stack = c(fixef_title, FE_lines, slope_intro, slope_lines)

    # Now we place the fixed-effects
    if(style$fixef.where == "var"){
        coef_stack = c(coef_stack, fixef_stack)
    } else {
        stat_stack = c(stat_stack, fixef_stack)
    }

    if(tex_tag){
        start_tag = "%start:tab\n"
        end_tag = "%end:tab\n"
    } else {
        start_tag = end_tag = ""
    }

    res = c(supplemental_info, start_table, start_tag, intro_latex, first_line, info_subtitles, model_line, info_family, coef_stack, stat_stack, info_SD, style$line.bottom, outro_latex, end_tag, info_notes, end_table)

    res = res[nchar(res) > 0]

    return(res)
}

etable_internal_df = function(info){

    n_models = length(info$depvar_list)
    # Getting the information
    se_type_list = info$se_type_list
    var_list = info$var_list
    coef_list = info$coef_list
    coef_below = info$coef_below
    sd_below = info$sd_below
    depvar_list = info$depvar_list
    convergence_list = info$convergence_list
    fe_names = info$fe_names
    is_fe = info$is_fe
    nb_fe = info$nb_fe
    slope_names = info$slope_names
    slope_flag_list = info$slope_flag_list
    family_list = info$family_list
    fitstat_list = info$fitstat_list
    title = info$title
    label = info$label
    keep = info$keep
    drop = info$drop
    order = info$order
    file = info$file
    family = info$family
    convergence = info$convergence
    sdBelow = info$sdBelow
    signifCode = info$signifCode
    fixef_sizes = info$fixef_sizes
    depvar = info$depvar
    useSummary = info$useSummary
    model_names = info$model_names
    coefstat = info$coefstat
    dict = info$dict
    group = info$group
    extraline = info$extraline
    style = info$style
    fun_format = info$fun_format
    drop.section = info$drop.section
    coef.just = info$coef.just

    # naming differences
    subtitles = info$subtitles
    isSubtitles = info$isSubtitles

    # The coefficients

    depvar_list = dict_apply(unlist(depvar_list), dict)

    # all_vars <- unique(c(var_list, recursive=TRUE))
    # we need to loop not to lose names
    all_vars = c()
    for(vars in var_list){
        all_vars = c(all_vars, vars[!vars %in% all_vars])
    }

    for(i in seq_along(group)){
        gi = group[[i]]

        present = sapply(var_list, function(x) any(keep_apply(x, gi, TRUE)))

        group[[i]] = present
    }

    # keeping some coefs
    all_vars = keep_apply(all_vars, keep)

    # dropping some coefs
    all_vars = drop_apply(all_vars, drop)

    # ordering the coefs
    all_vars = order_apply(all_vars, order)

    if(length(all_vars) == 0) stop_up("Not any variable was selected, please reframe your keep/drop arguments.")

    sdBelow = info$sdBelow
    if(sdBelow){
        coef_below = info$coef_below
        sd_below = info$sd_below
        coef_se_mat = c()
        for(v in all_vars){
            myCoef = mySd= myLine = c()
            for(m in 1:n_models){
                myCoef = c(myCoef, coef_below[[m]][v])
                mySd = c(mySd, sd_below[[m]][v])
            }

            myCoef[is.na(myCoef)] = "  "
            mySd[is.na(mySd)] = "  "

            coef_se_mat = rbind(coef_se_mat, myCoef, mySd)
        }

        # The tricky part: the row names!!!!
        # => this is a pain in the neck and currently the behavior is a bit odd
        # but so be it.
        all_names = rep(all_vars, each = 2)
        n_vars = length(all_vars)
        # We create batches of names
        empty_names = character(n_vars)
        my_batch = sprintf("% *s", 2:15, " ")

        if(n_vars <= 14){
            empty_names = my_batch[1:n_vars]
        } else {
            empty_names[1:14] = my_batch

            # This accounts for thousands of variables without bug
            all_chars_raw = all_chars = c(".", ":", "_", "*", ";", "~", "-", "=", letters)
            n_round = 1
            my_char = all_chars[1]
            i_char = 1
            i = 14
            while(i < n_vars){
                if(nchar(my_char) > 10){
                    i_char = i_char + 1

                    if(i_char > length(all_chars)){
                        # Specific case => many many vars
                        n_round = n_round + 1
                        my_list = list()
                        for(r in 1:n_round){
                            my_list[[r]] = all_chars_raw
                        }
                        quoi = do.call(expand.grid, my_list)
                        quoi = quoi[!apply(quoi, 1, function(x) length(unique(x)) == 1), ]
                        new_chars = apply(quoi, 1, paste, collapse = "")
                        new_chars = setdiff(new_chars, empty_names)
                        all_chars = new_chars

                        i_char = 1
                    }
                    my_char = all_chars[i_char]
                }

                my_batch = sprintf("%s% *s", my_char, 0:(15 - nchar(my_char)), "")
                n_batch = length(my_batch)
                n_max = min(n_vars - i, n_batch)
                empty_names[i + 1:n_max] = my_batch[1:n_max]

                i = i + n_max
                my_char = paste0(my_char, all_chars_raw[i_char])
            }
        }

        my_names = character(2 * n_vars)
        my_names[1 + 2 * 0:(n_vars - 1)] = all_vars
        my_names[2 + 2 * 0:(n_vars - 1)] = empty_names

        coef_mat = cbind(my_names, coef_se_mat)
    } else {
        coef_mat <- all_vars
        for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
        coef_mat[is.na(coef_mat)] <- "  "
    }

    if(!is.null(coef.just)){

        if(coef.just == "."){
            align_fun = function(x) sfill(sfill(x, anchor = "."), right = TRUE)
        } else if(coef.just == "("){
            align_fun = function(x) sfill(sfill(x, anchor = "("), right = TRUE)
        } else if(coef.just == "c"){
            align_fun = function(x) format(x, justify = "centre")
        } else if(coef.just == "l"){
            align_fun = function(x) sfill(x, anchor = "left")
        }

        new_mat = list(coef_mat[, 1])
        for(i in 2:ncol(coef_mat)){
            new_mat[[i]] = align_fun(coef_mat[, i])
        }

        coef_mat = do.call("cbind", new_mat)
    }

    res = coef_mat

    #
    # Group
    #

    before_stat = after_stat = c()
    before_fixef = after_fixef = c()

    for(i in seq_along(group)){
        gi = group[[i]]
        gi_format = style$yesNo[2 - gi]

        gi_name = names(group)[i]
        gi_where = "coef"

        gi_top = FALSE
        if(grepl("^(\\^|_|-)", gi_name)){
            # sec: section

            row = substr(gi_name, 1, 1)
            gi_name = substr(gi_name, 2, nchar(gi_name))
            if(row == "-"){
                # FE section
                sec = "fixef"

            } else if(grepl("^(\\^|_|-)", gi_name)){
                sec = substr(gi_name, 1, 1)
                gi_name = substr(gi_name, 2, nchar(gi_name))

            } else {
                if(row == "^"){
                    # implicit location
                    sec = "coef"
                } else {
                    sec = row
                    row = "_"
                }
            }

            gi_top = row == "^"
            gi_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat", sec)
        }

        my_line = c(gi_name, gi_format)

        if(gi_top){
            if(gi_where == "coef"){
                res = rbind(my_line, res)
            } else if(gi_where == "fixef"){
                before_fixef = rbind(before_fixef, my_line)
            } else {
                before_stat = rbind(before_stat, my_line)
            }
        } else {
            if(gi_where == "coef"){
                res = rbind(res, my_line)
            } else if(gi_where == "fixef"){
                after_fixef = rbind(after_fixef, my_line)
            } else {
                after_stat = rbind(after_stat, my_line)
            }
        }

    }

    #
    # Extra lines
    #

    for(i in seq_along(extraline)){
        el = extraline[[i]]
        yesNo = style$yesNo

        # The format depends on the type
        if(is.logical(el)){
            if(length(el) == 1){
                el_format = rep(yesNo[2 - el], n_models)
            } else {
                el_format = yesNo[2 - el]
            }
        } else if(is.numeric(el)){
            el_format = fun_format(el)
        } else {
            el_format = el
        }

        el_format[is.na(el_format)] = ""

        el_where = "coef"

        el_name = names(extraline)[i]
        el_top = FALSE

        if(grepl("^(\\^|_|-)", el_name)){
            # sec: section

            row = substr(el_name, 1, 1)
            el_name = substr(el_name, 2, nchar(el_name))
            if(row == "-"){
                # FE section
                sec = "fixef"

            } else if(grepl("^(\\^|_|-)", el_name)){
                sec = substr(el_name, 1, 1)
                el_name = substr(el_name, 2, nchar(el_name))

            } else {
                if(row == "^"){
                    # implicit location
                    sec = "coef"
                } else {
                    sec = row
                    row = "_"
                }
            }

            el_top = row == "^"
            el_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat", sec)
        }

        my_line = c(el_name, el_format)

        if(el_top){
            if(el_where == "coef"){
                res = rbind(my_line, res)
            } else if(el_where == "fixef"){
                before_fixef = rbind(before_fixef, my_line)
            } else {
                before_stat = rbind(before_stat, my_line)
            }
        } else {
            if(el_where == "coef"){
                res = rbind(res, my_line)
            } else if(el_where == "fixef"){
                after_fixef = rbind(after_fixef, my_line)
            } else {
                after_stat = rbind(after_stat, my_line)
            }
        }
    }

    #
    # Depvar
    #

    # The line with the dependent variable => defined here to get the width
    preamble = c()
    dep_width = 0
    if(depvar){
        preamble = rbind(c(style$depvar.title, depvar_list), preamble)
        dep_width = nchar(as.vector(preamble))
    }

    # the subtitles
    if(isSubtitles){
        # we need to provide unique names... sigh...

        n_sub = length(subtitles)
        sub_names = names(subtitles)
        if(is.null(sub_names)){
            sub_names = character(n_sub)
        }

        for(i in 1:n_sub){
            preamble = rbind(c(sfill(sub_names[i], i), subtitles[[i]]), preamble)
        }
    }

    # Used to draw lines
    longueur = apply(res, 2, function(x) max(nchar(as.character(x))))
    longueur = pmax(dep_width, longueur)

    #
    # fixed-effects
    #

    if(length(fe_names) > 0 && !"fixef" %in% drop.section){

        for(m in 1:n_models) {
            quoi = is_fe[[m]][fe_names]
            quoi[is.na(quoi)] = style$yesNo[2]
            is_fe[[m]] = quoi
        }

        fe_names = rename_fe(fe_names, dict)

        if(nchar(style$fixef.prefix) > 0){
            fe_names = paste0(style$fixef.prefix, fe_names)
        }

        if(nchar(style$fixef.suffix) > 0){
            fe_names = paste0(fe_names, style$fixef.suffix)
        }

        all_fe = matrix(c(is_fe, recursive=TRUE), nrow = length(fe_names))
        all_fe = cbind(fe_names, all_fe)

        if(nchar(style$fixef.title) > 0){
            myLine = paste(rep(style$fixef.line, 30), collapse = "")
            res = rbind(res, c(style$fixef.title, sprintf("%.*s", longueur[-1], myLine)))
        }

        if(length(before_fixef) > 0){
            res = rbind(res, before_fixef)
        }

        res = rbind(res, all_fe)
    } else {

        if(length(before_fixef) > 0 || length(after_fixef) > 0){
            # We create the fixed-effects title if needed
            if(!"fixef" %in% drop.section && nchar(style$fixef.title) > 0){
                myLine = paste(rep(style$fixef.line, 30), collapse = "")
                res = rbind(res, c(style$fixef.title, sprintf("%.*s", longueur[-1], myLine)))
            }
        }

        if(length(before_fixef) > 0){
            res = rbind(res, before_fixef)
        }
    }

    if(length(after_fixef) > 0){
        res = rbind(res, after_fixef)
    }

    #
    # The slopes
    #

    if(length(slope_names) > 0 && !"slopes" %in% drop.section){

        # reformatting the yes/no
        for(m in 1:n_models) {
            quoi = slope_flag_list[[m]][slope_names]
            quoi[is.na(quoi)] = style$yesNo[2]
            slope_flag_list[[m]] = quoi
        }

        all_slopes = matrix(c(slope_flag_list, recursive=TRUE), nrow = length(slope_names))
        all_slopes = cbind(slope_names, all_slopes)

        if(nchar(style$slopes.title) > 0){
            myLine = paste(rep(style$slopes.line, 30), collapse = "")
            res = rbind(res, c(style$slopes.title, sprintf("%.*s", longueur[-1], myLine)))
        }

        res = rbind(res, all_slopes)

    }

    # preamble created before because used to set the width
    if(length(preamble) > 0){
        preamble = rbind(preamble, rep(" ", length(longueur)))
        res = rbind(preamble, res)
    }

    if(nchar(style$stats.title) > 0){

        if(nchar(style$stats.title) == 1){
            n = max(nchar(res[, 1]), 12)
            if(all(grepl("(", unlist(se_type_list), fixed = TRUE))){
                n = max(n, 15)
            }

            fit_names = attr(fitstat_list, "format_names")
            if(length(fit_names) > 0){
                n = max(n, max(nchar(fit_names)))
            }

            myLine = paste(rep(style$stats.title, 40), collapse = "")
            stats_title = sprintf("%.*s", n, myLine)

        } else {
            stats_title = style$stats.title
        }

        myLine = paste(rep(style$stats.line, 30), collapse = "")
        res = rbind(res, c(stats_title, sprintf("%.*s", longueur[-1], myLine)))
    }

    # the line with the families
    if(family){
        res = rbind(res, c("Family", unlist(family_list)))
    }

    if(coefstat == "se"){
        coefstat_sentence = "S.E. type"
    } else {
        coefstat_sentence = "VCOV type"
    }

    se_type_format = c()
    for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]], by = TRUE)

    main_type = ""
    if(all(grepl("(", unlist(se_type_list), fixed = TRUE))){
        main_type = ": Clustered"
        coefstat_sentence = gsub(" type", "", coefstat_sentence)
    }

    res = rbind(res, c(paste0(coefstat_sentence, main_type), c(se_type_format, recursive = TRUE)))

    # convergence status
    if(convergence){
        res = rbind(res, c("Convergence", c(convergence_list, recursive = TRUE)))
    }

    #
    # Fit statistics
    #

    if(length(before_stat) > 0){
        res = rbind(res, before_stat)
    }

    if(!"stats" %in% drop.section && !all(sapply(fitstat_list, function(x) all(is.na(x))))){

        fit_names = attr(fitstat_list, "format_names")
        nb = length(fit_names)
        for(fit_id in 1:nb){
            fit = sapply(fitstat_list, function(x) x[[fit_id]])
            if(all(is.na(fit))) next
            fit[is.na(fit)] = "--"
            res = rbind(res, c(fit_names[fit_id], fit))
        }
    }

    if(length(after_stat) > 0){
        res = rbind(res, after_stat)
    }

    # if titles
    modelNames = model_names

    # we shorten the model names to fit the width
    for(m in 1:n_models) modelNames[m] = charShorten(modelNames[[m]], longueur[[1+m]])

    # I need to do that otherwise as.data.frame goes wild
    res_list = list()
    for(i in 1:ncol(res)) res_list[[i]] = unlist(res[, i])
    names(res_list) = paste0("x", 1:length(res_list))

    res = as.data.frame(res_list)
    names(res) = c("variables", modelNames)
    tvar = table(res$variables)
    if(any(tvar > 1)){
        qui = which(res$variables %in% names(tvar)[tvar > 1])
        add_space = c("", " ")
        if(length(qui) > 2) for(i in 3:length(qui)) add_space[i] = paste(rep(" ", i), collapse = "")
        res$variables = as.character(res$variables)
        res$variables[qui] = paste0(res$variables[qui], add_space)
    }

    row.names(res) = uniquify_names(unlist(res$variables))
    res$variables = NULL

    # We rename theta when NB is used
    quiTheta = which(row.names(res) == ".theta")
    row.names(res)[quiTheta] = "Over-dispersion"

    if(!is.null(file)){
        sink(file = file, append = !replace)
        on.exit(sink())

        print(res)

        return(invisible(res))
    } else {
        return(res)
    }

}


#' @rdname etable
setFixest_etable = function(digits = 4, digits.stats = 5, fitstat, coefstat = c("se", "tstat", "confint"), ci = 0.95, sdBelow = TRUE, keep, drop, order, dict, signifCode, float, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ ", depvar, style.tex = NULL, style.df = NULL, notes = NULL, group = NULL, extraline = NULL, fixef.group = NULL, placement = "htbp", drop.section = NULL, postprocess.tex = NULL, postprocess.df = NULL, fit_format = "__var__", reset = FALSE){

    # cat(names(formals(setFixest_etable)), sep = '", "')
    arg_list = c("digits", "digits.stats", "fitstat", "coefstat", "ci", "sdBelow", "keep", "drop", "order", "dict", "signifCode", "float", "fixef_sizes", "fixef_sizes.simplify", "family", "powerBelow", "interaction.combine", "depvar", "style.tex", "style.df", "notes", "group", "extraline", "placement", "drop.section", "postprocess.tex", "postprocess.df", "fit_format", "fixef.group", "reset")

    #
    # Argument checking => strong since these will become default values
    #

    check_set_digits(digits)
    check_set_digits(digits.stats)

    # fitstat (chiant) => controle reporte a fitstat_validate
    if(!missing(fitstat)){
        fitstat = fitstat_validate(fitstat)
    }


    check_arg_plus(coefstat, "match")
    check_arg(ci, "numeric scalar GT{0.5} LT{1}")

    check_arg("logical scalar", sdBelow, fixef_sizes, fixef_sizes.simplify, float, family, depvar, reset)

    check_arg(keep, drop, order, "character vector no na NULL", .message = "The arg. '__ARG__' must be a vector of regular expressions (see help(regex)).")

    check_arg_plus(signifCode, "NULL NA | match(letters) | named numeric vector no na GE{0} LE{1}")

    check_arg(interaction.combine, "character scalar")

    check_arg(notes, "character vector no na")

    check_arg(powerBelow, "integer scalar LE{-1}")

    check_arg(dict, "NULL logical scalar | named character vector no na")

    check_arg_plus(group, extraline, "NULL{list()} named list l0")

    check_arg_plus(fixef.group, "NULL{list()} logical scalar | named list l0")

    check_arg(placement, "character scalar")
    if(!missing(placement)){
        if(nchar(placement) == 0) stop("Argument 'placement' cannot be the empty string.")
        p_split = strsplit(placement, "")[[1]]
        check_value(p_split, "strict multi charin(h, t, b, p, H, !)", .message = "Argument 'placement' must be a character string containing only the following characters: 'h', 't', 'b', 'p', 'H', and '!'.")
    }

    check_arg_plus(drop.section, "NULL multi match(fixef, slopes, stats)")

    check_arg(style.tex, "NULL class(fixest_style_tex)")
    if(length(style.tex) > 0){
        # We ensure we always have ALL components provided
        basic_style = fixest::style.tex(main = "base")
        basic_style[names(style.tex)] = style.tex
        style.tex = basic_style
    }

    check_arg(style.df, "NULL class(fixest_style_df)")

    check_arg(postprocess.tex, postprocess.df, "NULL function arg(1,)")

    check_arg(fit_format, "character scalar")
    if(!grepl("__var__", fit_format, fixed = TRUE)){
        stop("The argument 'fit_format' should include the special name '__var__' that will be replaced by the variable name. So far it does not contain it.")
    }

    #
    # Setting the defaults
    #

    # Getting the existing defaults
    opts = getOption("fixest_etable")

    if(is.null(opts)){
        opts = list()
    } else if(!is.list(opts)){
        warning("Wrong formatting of option 'fixest_etable', all options are reset.")
        opts = list()
    } else if(reset){
        opts = list()
    }

    # Saving the default values
    mc = match.call()
    args_default = intersect(names(mc), arg_list)

    # NOTA: we don't allow delayed evaluation => all arguments must have hard values
    for(v in args_default){
        opts[[v]] = eval(as.name(v))
    }

    options(fixest_etable = opts)

}

#' @rdname etable
getFixest_etable = function(){
    opts = getOption("fixest_etable")
    if(!is.list(opts)){
        warning("Wrong formatting of option 'fplot_distr', all options are reset.")
        opts = list()
        options(fixest_etable = opts)
    }
    opts
}



tex_multicol = function(x){

    if(length(x) == 1) return(x)

    nb_multi = c(1)
    index = 1
    names_multi = x_current = x[1]

    for(val in x[-1]){
        if(val == x_current){
            nb_multi[index] = nb_multi[index] + 1
        } else {
            index = index + 1
            nb_multi[index] = 1
            names_multi[index] = x_current = val
        }
    }

    for(i in seq_along(nb_multi)){
        if(nb_multi[i] > 1){
            names_multi[i] = paste0("\\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
        }
    }

    res = paste(names_multi, collapse = " & ")

    return(res)
}

rename_fe = function(fe_all, dict){
    # Function used to rename the FEs

    res = c()

    for(i in seq_along(fe_all)){

        fe = fe_all[i]
        if(fe %in% names(dict)){
            res[i] = dict[fe]

        } else if(grepl("\\^", fe)){
            fe_split = strsplit(fe, "\\^")[[1]]
            who = fe_split %in% names(dict)
            fe_split[who] = dict[fe_split[who]]
            res[i] = paste(fe_split, collapse = "-")

        } else {
            res[i] = fe
        }
    }

    res
}



#' Style definitions for Latex tables
#'
#' This function describes the style of Latex tables to be exported with the function \code{\link[fixest]{etable}}.
#'
#' @param main Either "base", "aer" or "qje". Defines the basic style to start from. The styles "aer" and "qje" are almost identical and only differ on the top/bottom lines.
#' @param depvar.title A character scalar. The title of the line of the dependent variables (defaults to \code{"Dependent variable(s):"} if \code{main = "base"} (the 's' appears only if just one variable) and to \code{""} if \code{main = "aer"}).
#' @param model.title A character scalar. The title of the line of the models (defaults to \code{"Model:"} if \code{main = "base"} and to \code{""} if \code{main = "aer"}).
#' @param model.format A character scalar. The value to appear on top of each column. It defaults to \code{"(1)"}. Note that 1, i, I, a and A are special characters: if found, their values will be automatically incremented across columns.
#' @param line.top A character scalar. The line at the top of the table (defaults to \code{"\\tabularnewline\\toprule\\toprule"} if \code{main = "base"} and to \code{"\\toprule"} if \code{main = "aer"}).
#' @param line.bottom A character scalar. The line at the bottom of the table (defaults to \code{""} if \code{main = "base"} and to \code{"\\bottomrule"} if \code{main = "aer"}).
#' @param var.title A character scalar. The title line appearing before the variables (defaults to \code{"\\midrule \\emph{Variables}"} if \code{main = "base"} and to \code{"\\midrule"} if \code{main = "aer"}). Note that the behavior of \code{var.title = " "} (a space) is different from \code{var.title = ""} (the empty string): in the first case you will get an empty row, while in the second case you get no empty row. To get a line without an empty row, use \code{"\\midrule"} (and not \code{"\\midrule "}!--the space!).
#' @param fixef.title A character scalar. The title line appearing before the fixed-effects (defaults to \code{"\\midrule \\emph{Fixed-effects}"} if \code{main = "base"} and to \code{" "} if \code{main = "aer"}). Note that the behavior of \code{fixef.title = " "} (a space) is different from \code{fixef.title = ""} (the empty string): in the first case you will get an empty row, while in the second case you get no empty row. To get a line without an empty row, use \code{"\\midrule"} (and not \code{"\\midrule "}!--the space!).
#' @param fixef.prefix A prefix to add to the fixed-effects names. Defaults to \code{""} (i.e. no prefix).
#' @param fixef.suffix A suffix to add to the fixed-effects names. Defaults to \code{""} if \code{main = "base"}) and to \code{"fixed-effects"} if \code{main = "aer"}).
#' @param fixef.where Either "var" or "stats". Where to place the fixed-effects lines? Defaults to \code{"var"}, i.e. just after the variables, if \code{main = "base"}) and to \code{"stats"}, i.e. just after the statistics, if \code{main = "aer"}).
#' @param slopes.title A character scalar. The title line appearing before the variables with varying slopes (defaults to \code{"\\midrule \\emph{Varying Slopes}"} if \code{main = "base"} and to \code{""} if \code{main = "aer"}). Note that the behavior of \code{slopes.title = " "} (a space) is different from \code{slopes.title = ""} (the empty string): in the first case you will get an empty row, while in the second case you get no empty row. To get a line without an empty row, use \code{"\\midrule"} (and not \code{"\\midrule "}!--the space!).
#' @param slopes.format Character scalar representing the format of the slope variable name. There are two special characters: "__var__" and "__slope__", placeholers for the variable and slope names. Defaults to \code{"__var__ (__slope__)"} if \code{main = "base"}) and to \code{"__var__ $\\times $ __slope__"} if \code{main = "aer"}).
#' @param fixef_sizes.prefix A prefix to add to the fixed-effects names. Defaults to \code{"# "}.
#' @param fixef_sizes.suffix A suffix to add to the fixed-effects names. Defaults to \code{""} (i.e. no suffix).
#' @param stats.title A character scalar. The title line appearing before the statistics (defaults to \code{"\\midrule \\emph{Fit statistics}"} if \code{main = "base"} and to \code{" "} if \code{main = "aer"}). Note that the behavior of \code{stats.title = " "} (a space) is different from \code{stats.title = ""} (the empty string): in the first case you will get an empty row, while in the second case you get no empty row. To get a line without an empty row, use \code{"\\midrule"} (and not \code{"\\midrule "}!--the space!).
#' @param notes.title A character scalar. The title appearing just before the notes, defaults to \code{"\\medskip \\emph{Notes:} "}.
#' @param tablefoot A logical scalar. Whether or not to display a footer within the table. Defaults to \code{TRUE} if \code{main = "aer"}) and \code{FALSE} if \code{main = "aer"}).
#' @param tablefoot.title A character scalar. Only if \code{tablefoot = TRUE}, value to appear before the table footer. Defaults to \code{"\\bottomrule\\bottomrule"} if \code{main = "base"}.
#' @param tablefoot.value A character scalar. The notes to be displayed in the footer. Defaults to \code{"default"} if \code{main = "base"}, which leads to custom footers informing on the type of standard-error and significance codes, depending on the estimations.
#' @param yesNo A character vector of length 1 or 2. Defaults to \code{"Yes"} if \code{main = "base"} and to \code{"$\\checkmark$"} if \code{main = "aer"} (from package \code{amssymb}). This is the message displayed when a given fixed-effect is (or is not) included in a regression. If \code{yesNo} is of length 1, then the second element is the empty string.
#' @param tabular Character scalar equal to "normal" (default), "*" or "X". Represents the type of tabular to export.
#'
#' @details
#' The \code{\\checkmark} command, used in the "aer" style (in argument \code{yesNo}), is in the \code{amssymb} package.
#'
#' The commands \code{\\toprule}, \code{\\midrule} and \code{\\bottomrule} are in the \code{booktabs} package. You can set the width of the top/bottom rules with \\setlength\\heavyrulewidth\{wd\}, and of the midrule with \\setlength\\lightrulewidth\{wd\}.
#'
#' @return
#' Returns a list containing the style parameters.
#'
#' @seealso
#' \code{\link[fixest]{etable}}
#'
#' @examples
#'
#' # Multiple estimations => see details in feols
#' aq = airquality
#' est = feols(c(Ozone, Solar.R) ~
#'                 Wind + csw(Temp, Temp^2, Temp^3) | Month + Day,
#'             data = aq)
#'
#' # Playing a bit with the styles
#' etable(est, tex = TRUE)
#' etable(est, tex = TRUE, style.tex = style.tex("aer"))
#'
#' etable(est, tex = TRUE, style.tex = style.tex("aer",
#'                                       var.title = "\\emph{Expl. Vars.}",
#'                                       model.format = "[i]",
#'                                       yesNo = "x",
#'                                       tabular = "*"))
#'
style.tex = function(main = "base", depvar.title, model.title, model.format, line.top, line.bottom, var.title, fixef.title, fixef.prefix, fixef.suffix, fixef.where, slopes.title, slopes.format, fixef_sizes.prefix, fixef_sizes.suffix, stats.title, notes.title, tablefoot, tablefoot.title, tablefoot.value, yesNo, tabular = "normal"){

    # To implement later:
    # fixef_sizes.where = "obs"
    # se.par = "parentheses"
    # check_arg_plus(se.par, "match(parentheses, brackets, none)")
    # align
    # check_arg(align, "character vector no na len(,2)")

    # Checking
    check_arg_plus(main, "match(base, aer, qje)")

    check_arg("character scalar", depvar.title, model.title, line.top, line.bottom, var.title)
    check_arg("character scalar", fixef.title, fixef.prefix, fixef.suffix, slopes.title, slopes.format)
    check_arg("character scalar", fixef_sizes.prefix, fixef_sizes.suffix, stats.title)
    check_arg("character scalar", notes.title, tablefoot.title)

    check_arg(tablefoot.value, "character vector no na")
    check_arg(tablefoot, "logical scalar")
    check_arg_plus(fixef.where, "match(var, stats)")
    check_arg_plus(tabular, "match(normal, *, X)")

    check_arg(yesNo, "character vector len(,2) no na")
    if(!missing(yesNo) && length(yesNo) == 1){
        yesNo = c(yesNo, "")
    }

    mc = match.call()

    if("main" %in% names(mc)){
        if(main == "base"){
            res = list(depvar.title = "Dependent Variable(s):", model.title = "Model:", model.format = "(1)",
                       line.top = "\\tabularnewline\\midrule\\midrule", line.bottom = "",
                       var.title = "\\midrule \\emph{Variables}",
                       fixef.title = "\\midrule \\emph{Fixed-effects}", fixef.prefix = "", fixef.suffix = "", fixef.where = "var",
                       slopes.title = "\\midrule \\emph{Varying Slopes}", slopes.format = "__var__ (__slope__)",
                       fixef_sizes.prefix = "# ", fixef_sizes.suffix = "",
                       stats.title = "\\midrule \\emph{Fit statistics}", notes.title = "\\medskip \\emph{Notes:} ",
                       tablefoot = TRUE, tablefoot.title = "\\midrule\\midrule", tablefoot.value = "default", yesNo = c("Yes", ""))

        } else {
            res = list(depvar.title = "", model.title = "", model.format = "(1)",
                       line.top = "\\toprule", line.bottom = "\\bottomrule",
                       var.title = "\\midrule",
                       fixef.title = " ", fixef.prefix = "", fixef.suffix = " fixed effects", fixef.where = "stats",
                       slopes.title = "", slopes.format = "__var__ $\\times $ __slope__",
                       fixef_sizes.prefix = "# ", fixef_sizes.suffix = "",
                       stats.title = " ", notes.title = "\\medskip \\emph{Notes:} ",
                       tablefoot = FALSE, tablefoot.title = "", tablefoot.value = "", yesNo = c("$\\checkmark$", ""))

            if(main == "aer"){
                # just set

            } else if(main == "qje"){
                res$line.top = "\\tabularnewline\\toprule\\toprule"
                res$line.bottom = "\\bottomrule\\bottomrule & \\tabularnewline"
            }
        }

        res$tabular = tabular
    } else {
        res = list()
    }

    args2set = setdiff(names(mc)[-1], "main")

    for(var in args2set){
        res[[var]] = eval(as.name(var))
    }

    class(res) = "fixest_style_tex"

    return(res)
}



#' Style of data.frames created by etable
#'
#' This function describes the style of data.frames created with the function \code{\link[fixest]{etable}}.
#'
#' @param depvar.title Character scalar. Default is \code{"Dependent Var.:"}. The row name of the dependent variables.
#' @param fixef.title Character scalar. Default is \code{"Fixed-Effects:"}. The header preceding the fixed-effects. If equal to the empty string, then this line is removed.
#' @param fixef.line A single character. Default is \code{"-"}. A character that will be used to create a line of separation for the fixed-effects header. Used only if \code{fixef.title} is not the empty string.
#' @param fixef.prefix Character scalar. Default is \code{""}. A prefix to appear before each fixed-effect name.
#' @param fixef.suffix Character scalar. Default is \code{""}. A suffix to appear after each fixed-effect name.
#' @param slopes.title Character scalar. Default is \code{"Varying Slopes:"}. The header preceding the variables with varying slopes. If equal to the empty string, then this line is removed.
#' @param slopes.line Character scalar. Default is \code{"-"}. A character that will be used to create a line of separation for the variables with varying slopes header. Used only if \code{slopes.line} is not the empty string.
#' @param slopes.format Character scalar. Default is \code{"__var__ (__slope__)"}. The format of the name of the varying slopes. The values \code{__var__} and \code{__slope__} are special characters that will be replaced by the value of the variable name and slope name, respectively.
#' @param stats.title Character scalar. Default is \code{"_"}. The header preceding the statistics section. If equal to the empty string, then this line is removed. If equal to single character (like in the default), then this character will be expanded to take the full column width.
#' @param stats.line Character scalar. Default is \code{"_"}. A character that will be used to create a line of separation for the statistics header. Used only if \code{stats.title} is not the empty string.
#' @param yesNo Character vector of length 1 or 2. Default is \code{c("Yes", "No")}. Used to inform on the presence or absence of fixed-effects in the estimation. If of length 1, then automatically the second value is considered as the empty string.
#'
#' @details
#' The title elements (\code{depvar.title}, \code{fixef.title}, \code{slopes.title} and \code{stats.title}) will be the row names of the returned data.frame. Therefore keep in mind that any two of them should not be identical (since identical row names are forbidden in data.frames).
#'
#' @return
#' It returns an object of class \code{fixest_style_df}.
#'
#' @examples
#'
#' # Multiple estimations => see details in feols
#' aq = airquality
#' est = feols(c(Ozone, Solar.R) ~
#'                 Wind + csw(Temp, Temp^2, Temp^3) | Month + Day,
#'             data = aq)
#'
#'
#' # Default result
#' etable(est)
#'
#' # Playing a bit with the styles
#' etable(est, style_df = style.df(fixef.title = "", fixef.suffix = " FE",
#'                                  stats.line = " ", yesNo = "yes"))
#'
#'
style.df = function(depvar.title = "Dependent Var.:", fixef.title = "Fixed-Effects:", fixef.line = "-", fixef.prefix = "", fixef.suffix = "", slopes.title = "Varying Slopes:", slopes.line = "-", slopes.format = "__var__ (__slope__)", stats.title = "_", stats.line = "_", yesNo = c("Yes", "No")){

    # Checking

    check_arg("character scalar", depvar.title, fixef.title, fixef.line, fixef.prefix, fixef.suffix)
    check_arg("character scalar", slopes.title, slopes.line, slopes.format, stats.title, stats.line)

    check_arg(yesNo, "character vector len(,2) no na")
    if(length(yesNo) == 1){
        yesNo = c(yesNo, "")
    }

    if(nchar(fixef.line) != 1){
        stop("The argument 'fixef.line' must be a singe character! It's currently a string of length ", nchar(fixef.line), ".")
    }
    if(nchar(slopes.line) != 1){
        stop("The argument 'slopes.line' must be a singe character! It's currently a string of length ", nchar(slopes.line), ".")
    }
    if(nchar(stats.line) != 1){
        stop("The argument 'stats.line' must be a singe character! It's currently a string of length ", nchar(stats.line), ".")
    }

    res = list(depvar.title = depvar.title, fixef.title = fixef.title, fixef.line = fixef.line, fixef.prefix = fixef.prefix, fixef.suffix = fixef.suffix, slopes.title = slopes.title, slopes.line = slopes.line, slopes.format = slopes.format, stats.title = stats.title, stats.line = stats.line, yesNo = yesNo)

    class(res) = "fixest_style_df"

    return(res)
}

uniquify_names = function(x){
    # x: vector of names
    # we make each value of x unique by adding white spaces

    if(length(x) == 0) return(NULL)

    x_unik = unique(x)

    if(length(x_unik) == length(x)) return(x)

    x = gsub(" +$", " ", x)
    x_unik = unique(x)
    tab = rep(0, length(x_unik))
    names(tab) = x_unik

    x_new = x
    for(i in seq_along(x)){
        n = tab[x[i]]
        if(n > 0){
            x_new[i] = paste0(x_new[i], sprintf("% *s", n, ""))
        }
        tab[x[i]] = n + 1
    }

    x_new
}


extraline_extractor = function(x, name = NULL, tex = FALSE){
    # x must be a one sided formula
    # name: name of the listed element (empty = "")
    # => returns a named list

    is_name = !is.null(name) && nchar(name) > 0

    # extraline registered
    el_default = getOption("fixest_extraline")
    key_registered = names(el_default)

    # fitstat
    fitstat_fun_types = fitstat(give_types = TRUE)
    fitstat_type_allowed = fitstat_fun_types$types
    type_alias = fitstat_fun_types$type_alias

    # The variable(s) requested
    current_vars = attr(terms(x), "term.labels")

    if(length(current_vars) > 1 && is_name){
        stop_up("You cannot give list names in 'extraline' when several values are summoned via a formula. Simply remove the name associated to the formula to make it work (concerns '", name, "').", up = 2)
    }

    valid_keys = c(key_registered, fitstat_type_allowed)

    if(!all(current_vars %in% valid_keys)){
        pblm = setdiff(current_vars, valid_keys)
        stop_up("Argument 'extraline' can be a one-sided formula whose variables refer to the macros registered using the function 'extraline_register' or fit statistics (valid fitstat keywords). Problem: the values", enumerate_items(pblm, "s.were"), " not valid.", up = 2)
    }

    if(length(current_vars) == 0){
        extraline = list()
    } else {
        el_tmp = list()
        for(i in seq_along(current_vars)){

            key = current_vars[i]

            if(key %in% key_registered){
                element = el_default[[key]]

                if(is_name){
                    el_tmp[[name]] = element$fun
                } else {
                    el_tmp[[element$alias]] = element$fun
                }

            } else {
                # this is a fit statistic
                fun = eval(str2lang(paste0("function(est) fitstat(est, '", key, "', simplify = TRUE)[[1]]")))

                if(is_name){
                    el_tmp[[name]] = fun
                } else {
                    my_alias = if(tex) fitstat_fun_types$tex_alias else fitstat_fun_types$R_alias
                    el_tmp[[my_alias[[key]]]] = fun
                }
            }
        }
    }

    el_tmp
}


#' Register \code{extraline} macros to be used in \code{etable}
#'
#' This function is used to create \code{extraline} (which is an argument of \code{\link[fixest]{etable}}) macros that can be easily summoned in \code{\link[fixest]{etable}}.
#'
#' @param type A character scalar giving the type-name.
#' @param fun A function to be applied to a \code{fixest} estimation. It must return a scalar.
#' @param alias A character scalar. This is the alias to be used in lieu of the type name to form the row name.
#'
#' @details
#' You can register as many macros as you wish, the only constraint is that the type name should not conflict with a \code{\link[fixest]{fitstat}} type name.
#'
#'
#' @examples
#'
#'
#' # We register a function computing the standard-deviation of the dependent variable
#' my_fun = function(x) sd(model.matrix(x, type = "lhs"))
#' extraline_register("sdy", my_fun, "SD(y)")
#'
#' # An estimation
#' data(iris)
#' est = feols(Petal.Length ~ Sepal.Length | Species, iris)
#'
#' # Now we can easily create a row with the mean of y.
#' # We just "summon" it in a one-sided formula
#' etable(est, extraline = ~ sdy)
#'
#' # We can change the alias on the fly:
#' etable(est, extraline = list("_Standard deviation of the dep. var." = ~ sdy))
#'
#'
#'
#'
extraline_register = function(type, fun, alias){
    check_arg(type, "character scalar mbt")
    check_arg(fun, "function mbt")
    check_arg(alias, "character scalar mbt")

    # We check the type is not conflicting
    existing_types = fitstat(give_types = TRUE)$types

    opts = getOption("fixest_extraline")

    if(type %in% setdiff(existing_types, names(opts))){
        stop("The type name '", type, "' is the same as one of fitstat's built-in type. Please choose another one.")
    }

    if(missnull(alias)){
        alias = type
    }

    # We test the function on a simple estimation
    base = data.frame(y = rnorm(100), x = rnorm(100))
    est = feols(y ~ x, base)
    mc = match.call()
    fun_name = deparse_long(mc$fun)
    value = error_sender(fun(est), "The function '", fun_name, "' could not evaluated on a simple fixest object. Please try to improve it.")

    if(length(value) != 1){
        stop("The value returned by ", fun_name, " should be exactly of length 1. This is actually not the case (the result is of length ", length(value), ").")
    }

    res = list(fun = fun, alias = alias)

    opts[[type]] = res

    options(fixest_extraline = opts)

    invisible(NULL)
}

















































