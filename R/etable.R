#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Thu Jul 09 09:52:31 2020
# ~: etable
#----------------------------------------------#


#' Estimations table (export the results of multiples estimations to a DF or to Latex)
#'
#' Aggregates the results of multiple estimations and displays them in the form of either a Latex table or a \code{data.frame}. Note that you will need the \code{booktabs} package for the Latex table to render properly.
#'
#' @inheritParams summary.fixest
#' @inheritParams setFixest_nthreads
#'
#' @param ... Used to capture different \code{fixest} estimation objects (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). Note that any other type of element is discarded. Note that you can give a list of \code{fixest} objects.
#' @param digits Integer or character scalar. Default is 4 and represents the number of significant digits to be displayed for the coefficients and standard-errors. To apply rounding instead of significance use, e.g., \code{digits = "r3"} which will round at the first 3 decimals. If character, it must be of the form \code{"rd"} or \code{"sd"} with \code{d} a digit (\code{r} is for round and \code{s} is for significance). For the number of digits for the fit statistics, use \code{digits.stats}. Note that when significance is used it does not exactly display the number of significant digits: see details for its exact meaning.
#' @param digits.stats Integer or character scalar. Default is 5 and represents the number of significant digits to be displayed for the fit statistics. To apply rounding instead of significance use, e.g., \code{digits = "r3"} which will round at the first 3 decimals. If character, it must be of the form \code{"rd"} or \code{"sd"} with \code{d} a digit (\code{r} is for round and \code{s} is for significance). Note that when significance is used it does not exactly display the number of significant digits: see details for its exact meaning.
#' @param tex Logical: whether the results should be a data.frame or a Latex table. By default, this argument is \code{TRUE} if the argument \code{file} (used for exportation) is not missing; it is equal to \code{FALSE} otherwise.
#' @param fitstat A character vector or a one sided formula (both with only lowercase letters). A vector listing which fit statistics to display. The valid types are 'n', 'll', 'aic', 'bic' and r2 types like 'r2', 'pr2', 'war2', etc (see all valid types in \code{\link[fixest]{r2}}). Also accepts valid types from the function \code{\link[fixest]{fitstat}}. The default value depends on the models to display. Example of use: \code{fitstat=c('n', 'cor2', 'ar2', 'war2')}, or \code{fitstat=~n+cor2+ar2+war2} using a formula. You can use the dot to refer to default values:\code{ ~ . + ll} would add the log-likelihood to the default fit statistics.
#' @param title (Tex only.) Character scalar. The title of the Latex table.
#' @param float (Tex only.) Logical. By default, if the argument \code{title} or \code{label} is provided, it is set to \code{TRUE}. Otherwise, it is set to \code{FALSE}.
#' @param se.below Logical or \code{NULL} (default). Should the standard-errors be displayed below the coefficients? If \code{NULL}, then this is \code{TRUE} for Latex and \code{FALSE} otherwise.
#' @param se.row Logical scalar, default is \code{NULL}. Whether should be displayed the row with the type of standard-error for each model. When \code{tex = FALSE}, the default is \code{TRUE}. When \code{tex = FALSE}, the row is showed only when there is a table-footer and the types of standard-errors differ across models.
#' @param keep Character vector. This element is used to display only a subset of variables. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be kept. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{keep = "x[[:digit:]]$"}. If the first character is an exclamation mark, the effect is reversed (e.g. keep = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is kept). See details.
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be discarded. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{drop = "x[[:digit:]]{2}"}. If the first character is an exclamation mark, the effect is reversed (e.g. drop = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is dropped). See details.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions. This argument is applied post aliasing (see argument \code{dict}). Example: you have the following variables: \code{month1} to \code{month6}, then \code{x1} to \code{x5}, then \code{year1} to \code{year6}. If you want to display first the x's, then the years, then the months you could use: \code{order = c("x", "year")}. If the first character is an exclamation mark, the effect is reversed (e.g. order = "!Intercept" means: every variable that does not contain \dQuote{Intercept} goes first).  See details.
#' @param dict A named character vector or a logical scalar. It changes the original variable names to the ones contained in the \code{dict}ionary. E.g. to change the variables named \code{a} and \code{b3} to (resp.) \dQuote{$log(a)$} and to \dQuote{$bonus^3$}, use \code{dict=c(a="$log(a)$",b3="$bonus^3$")}. By default, it is equal to \code{getFixest_dict()}, a default dictionary which can be set with \code{\link[fixest]{setFixest_dict}}. You can use \code{dict = FALSE} to disable it. By default \code{dict} modifies the entries in the global dictionary, to disable this behavior, use "reset" as the first element (ex: \code{dict=c("reset", mpg="Miles per gallon")}).
#' @param file A character scalar. If provided, the Latex (or data frame) table will be saved in a file whose path is \code{file}. If you provide this argument, then a Latex table will be exported, to export a regular \code{data.frame}, use argument \code{tex = FALSE}.
#' @param replace Logical, default is \code{FALSE}. Only used if option \code{file} is used. Should the exported table be written in a new file that replaces any existing file?
#' @param convergence Logical, default is missing. Should the convergence state of the algorithm be displayed? By default, convergence information is displayed if at least one model did not converge.
#' @param signif.code Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)} for a Latex table and \code{c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)} for a data.frame (to conform with R's default). To suppress the significance codes, use \code{signif.code=NA} or \code{signif.code=NULL}. Can also be equal to \code{"letters"}, then the default becomes \code{c("a"=0.01, "b"=0.05, "c"=0.10)}.
#' @param label (Tex only.) Character scalar. The label of the Latex table.
#' @param headers Character vector or list. Adds one or more header lines in the table. A header line can be represented by a character vector or a named list of numbers where the names are the cell values and the numbers are the span. Example: \code{headers=list("M"=2, "F"=3)} will create a row with 2 times "M" and three time "F" (this is identical to \code{headers=rep(c("M", "F"), c(2, 3))}). You can stack header lines within a list, in that case the list names will be displayed in the leftmost cell. Example: \code{headers=list(Gender=list("M"=2, "F"=3), Country="US"} will create two header lines. When \code{tex = TRUE}, you can add a rule to separate groups by using \code{":_:"} somewhere in the row name (ex: \code{headers=list(":_:Gender"=list("M"=2, "F"=3))}. You can monitor the placement by inserting a special character in the row name: "^" means at the top, "-" means in the middle (default) and "_" means at the bottom. Example: \code{headers=list("_Country"="US")} will add the country row as the very last header row (after the model row). Finally, you can use the special value "auto" to include automatic headers when the data contains split sample estimations. By default it is equal to \code{list("auto")}. You can use \code{.()} instead of \code{list()}.
#' @param fixef_sizes (Tex only.) Logical, default is \code{FALSE}. If \code{TRUE} and fixed-effects were used in the models, then the number of "units" per fixed-effect dimension is also displayed.
#' @param fixef_sizes.simplify Logical, default is \code{TRUE}. Only used if \code{fixef_sizes = TRUE}. If \code{TRUE}, the fixed-effects sizes will be displayed in parentheses instead of in a separate line if there is no ambiguity (i.e. if the size is constant across models).
#' @param family Logical, default is missing. Whether to display the families of the models. By default this line is displayed when at least two models are from different families.
#' @param keepFactors Logical, default is \code{TRUE}. If \code{FALSE}, then factor variables are displayed as fixed-effects and no coefficient is shown.
#' @param powerBelow (Tex only.) Integer, default is -5. A coefficient whose value is below \code{10**(powerBelow+1)} is written with a power in Latex. For example \code{0.0000456} would be written \code{4.56$\\times 10^{-5}$} by default. Setting \code{powerBelow = -6} would lead to \code{0.00004} in Latex.
#' @param interaction.combine Character scalar, defaults to \code{" $\\times$ "} for Tex and to \code{" = "} otherwise. When the estimation contains interactions, then the variables names (after aliasing) are combined with this argument. For example: if \code{dict = c(x1="Wind", x2="Rain")} and you have the following interaction \code{x1:x2}, then it will be renamed (by default) \code{Wind $\\times$ Rain} -- using \code{interaction.combine = "*"} would lead to \code{Wind*Rain}.
#' @param interaction.order Character vector of regular expressions. Only affects variables that are interacted like x1 and x2 in \code{feols(y ~ x1*x2, data)}. You can change the order in which the interacted variables are displayed: e.g. \code{interaction.order = "x2"} would lead to "x1 x x2" instead of "x1 x x2". Please look at the argument 'order' and the dedicated section in the help page for more information.
#' @param i.equal Character scalar, defaults to \code{" $=$ "} when \code{tex = TRUE} and \code{" = "} otherwise. Only affects factor variables created with the function \code{\link[fixest]{i}}, tells how the variable should be linked to its value. For example if you have the \code{Species} factor from the \code{iris} data set, by default the display of the variable is \code{Species = Setosa}, etc. If \code{i.equal = ": "} the display becomes \code{Species: Setosa}.
#' @param depvar Logical, default is \code{TRUE}. Whether a first line containing the dependent variables should be shown.
#' @param coefstat One of \code{"se"} (default), \code{"tstat"} or \code{"confint"}. The statistic to report for each coefficient: the standard-error, the t-statistics or the confidence interval. You can adjust the confidence interval with the argument \code{ci}.
#' @param ci Level of the confidence interval, defaults to \code{0.95}. Only used if \code{coefstat = confint}.
#' @param style.tex An object created by the function \code{\link[fixest]{style.tex}}. It represents the style of the Latex table, see the documentation of \code{\link[fixest]{style.tex}}.
#' @param style.df An object created by the function \code{\link[fixest]{style.df}. }It represents the style of the data frame returned (if \code{tex = FALSE}), see the documentation of \code{\link[fixest]{style.df}}.
#' @param notes (Tex only.) Character vector. If provided, a \code{"notes"} section will be added at the end right after the end of the table, containing the text of this argument. If it is a vector, it will be collapsed with new lines. If \code{tpt = TRUE}, the behavior is different: each element of the vector is an item. If the first element of the vector starts with \code{"@"}, then it will be included verbatim, and in case of \code{tpt = TRUE}, right before the first item. If that element is provided, it will replace the value defined in \code{style.tex(notes.intro)} or \code{style.tex(notes.tpt.intro)}.
#' @param group A list. The list elements should be vectors of regular expressions. For each elements of this list: A new line in the table is created, all variables that are matched by the regular expressions are discarded (same effect as the argument \code{drop}) and \code{TRUE} or \code{FALSE} will appear in the model cell, depending on whether some of the previous variables were found in the model. Example: \code{group=list("Controls: personal traits"=c("gender", "height", "weight"))} will create an new line with \code{"Controls: personal traits"} in the leftmost cell, all three variables gender, height and weight are discarded, \code{TRUE} appearing in each model containing at least one of the three variables (the style of \code{TRUE}/\code{FALSE} is governed by the argument \code{yesNo}). You can control the placement of the new row by using 1 or 2 special characters at the start of the row name. The meaning of these special characters are: 1) \code{"^"}: coef., \code{"-"}: fixed-effect, \code{"_"}: stats, section; 2) \code{"^"}: 1st, \code{"_"}: last, row. For example: \code{group=list("_^Controls"=stuff)} will place the line at the top of the 'stats' section, and using \code{group=list("^_Controls"=stuff)} will make the row appear at the bottom of the coefficients section. For details, see the dedicated section.
#' @param extralines A vector, a list or a one sided formula. The list elements should be either a vector representing the value of each cell, a list of the form \code{list("item1" = #item1, "item2" = #item2, etc)}, or a function. This argument can be many things, please have a look at the dedicated help section; a simplified description follows. For each elements of this list: A new line in the table is created, the list name being the row name and the vector being the content of the cells. Example: \code{extralines=list("Sub-sample"=c("<20 yo", "all", ">50 yo"))} will create an new line with \code{"Sub-sample"} in the leftmost cell, the vector filling the content of the cells for the three models. You can control the placement of the new row by using 1 or 2 special characters at the start of the row name. The meaning of these special characters are: 1) \code{"^"}: coef., \code{"-"}: fixed-effect, \code{"_"}: stats, section; 2) \code{"^"}: 1st, \code{"_"}: last, row. For example: \code{extralines=list("__Controls"=stuff)} will place the line at the bottom of the stats section, and using \code{extralines=list("^^Controls"=stuff)} will make the row appear at the top of the 'coefficients' section. For details, see the dedicated section. You can use \code{.()} instead of \code{list()}.
#' @param fixef.group Logical scalar or list (default is \code{NULL}). If equal to \code{TRUE}, then all fixed-effects always appearing jointly in models will be grouped in one row. If a list, its elements must be character vectors of regular expressions and the list names will be the row names. For ex. \code{fixef.group=list("Dates fixed-effects"="Month|Day")} will remove the \code{"Month"} and \code{"Day"} fixed effects from the display and replace them with a single row named "Dates fixed-effects". You can monitor the placement of the new row with the special characters telling where to place the row within a section: \code{"^"} (first), or \code{"_"} (last); and in which section it should appear: \code{"^"} (coef.), \code{"-"} (fixed-effects), or \code{"_"} (stat.). These two special characters must appear first in the row names. Please see the dedicated section
#' @param placement (Tex only.) Character string giving the position of the float in Latex. Default is "htbp". It must consist of only the characters 'h', 't', 'b', 'p', 'H' and '!'. Reminder: h: here; t: top; b: bottom; p: float page; H: definitely here; !: prevents Latex to look for other positions. Note that it can be equal to the empty string (and you'll get the default placement).
#' @param drop.section Character vector which can be of length 0 (i.e. equal to \code{NULL}). Can contain the values "fixef", "slopes" or "stats". It would drop, respectively, the fixed-effects section, the variables with varying slopes section or the fit statistics section.
#' @param reset (\code{setFixest_etable} only.) Logical, default is \code{FALSE}. If \code{TRUE}, this will reset all the default values that were already set by the user in previous calls.
#' @param .vcov A function to be used to compute the standard-errors of each fixest object. You can pass extra arguments to this function using the argument \code{.vcov_args}. See the example.
#' @param .vcov_args A list containing arguments to be passed to the function \code{.vcov}.
#' @param poly_dict Character vector, default is \code{c("", " square", " cube")}. When raw polynomials (\code{x^2}, etc) are used, the variables are automatically renamed and \code{poly_dict} rules the display of the power. For powers greater than the number of elements of the vector, the value displayed is \code{$^{pow}$} in Latex and \code{^ pow} in the R console.
#' @param postprocess.tex A function that will postprocess the character vector defining the latex table. Only when \code{tex = TRUE}. By default it is equal to \code{NULL}, meaning that there is no postprocessing. When \code{tex = FALSE}, see the argument \code{postprocess.df}. See details.
#' @param postprocess.df A function that will postprocess.tex the resulting data.frame. Only when \code{tex = FALSE}. By default it is equal to \code{NULL}, meaning that there is no postprocessing. When \code{tex = TRUE}, see the argument \code{postprocess.tex}.
#' @param fit_format Character scalar, default is \code{"__var__"}. Only used in the presence of IVs. By default the endogenous regressors are named \code{fit_varname} in the second stage. The format of the endogenous regressor to appear in the table is governed by \code{fit_format}. For instance, by default, the prefix \code{"fit_"} is removed, leading to only \code{varname} to appear. If \code{fit_format = "$\\\\hat{__var__}$"}, then \code{"$\\hat{varname}$"} will appear in the table.
#' @param coef.just (DF only.) Either \code{"."}, \code{"("}, \code{"l"}, \code{"c"} or \code{"r"}, default is \code{NULL}. How the coefficients should be justified. If \code{NULL} then they are right aligned if \code{se.below = FALSE} and aligned to the dot if \code{se.below = TRUE}. The keywords stand respectively for dot-, parenthesis-, left-, center- and right-aligned.
#' @param meta (Tex only.) A one-sided formula that shall contain the following elements: date or time, sys, author, comment and call. Default is \code{NULL}. This argument is a shortcut to controlling the meta information that can be displayed in comments before the table. Typically if the element is in the formula, it means that the argument will be equal to \code{TRUE}. Example: \code{meta = ~time+call} is equivalent to \code{meta.time = TRUE} and \code{meta.call = TRUE}. The "author" and "comment" elements are a bit special. Using \code{meta = ~author("Mark")} is equivalent to \code{meta.author = "Mark"} while \code{meta=~author} is equiv. to \code{meta.author = TRUE}. The "comment" must be used with a character string inside: \code{meta = ~comment("this is a comment")}. The order in the formula controls the order of appearance of the meta elements. It also has precedence over the \code{meta.XX} arguments.
#' @param meta.time (Tex only.) Either a logical scalar (default is \code{FALSE}) or "time" or "date". Whether to include the time (if \code{TRUE} or "time") or the date (if "date") of creation of the table in a comment right before the table.
#' @param meta.sys (Tex only.) A logical scalar, default is \code{FALSE}. Whether to include system information (from \code{Sys.info()}) in a comment right before the table.
#' @param meta.author (Tex only.) A logical scalar (default is \code{FALSE}) or a character vector. If \code{TRUE} then the identity of the author (deduced from the system user in \code{Sys.info()}) is inserted in a comment right before the table. If a character vector, then it should contain author names that will be inserted as comments before the table, prefixed with \code{"Created by:"}. For free-form comments see the argument \code{meta.comment}.
#' @param meta.comment (Tex only.) A character vector containing free-form comments to be inserted right before the table.
#' @param meta.call (Tex only.) Logical scalar, default is \code{FALSE}. If \code{TRUE} then the call to the function is inserted right before the table in a comment.
#' @param view Logical, default is \code{FALSE}. If \code{TRUE}, then the table generated in Latex by \code{etable} and then is displayed in the viewer pane. Note that for this option to work you need i) pdflatex, ii) imagemagick and iii) ghostscript. All three software must be installed and on the path.
#' @param x An object returned by \code{etable}.
#' @param tpt (Tex only.) Logical scalar, default is FALSE. Whether to use the \code{threeparttable} environment. If so, the \code{notes} will be integrated into the \code{tablenotes} environment.
#' @param arraystretch (Tex only.) A numeric scalar, default is \code{NULL}. If provided, the command \code{\\renewcommand*{\\arraystretch}{x}} is inserted, replacing \code{x} by the value of \code{arraystretch}. The changes are specific to the current table and do not affect the rest of the document.
#' @param fontsize (Tex only.) A character scalar, default is \code{NULL}. Can be equal to \code{tiny}, \code{scriptsize}, \code{footnotesize}, \code{small}, \code{normalsize}, \code{large}, or \code{Large}. The change affect the table only (and not the rest of the document).
#' @param adjustbox (Tex only.) A logical, numeric or character scalar, default is \code{NULL}. If not \code{NULL}, the table is inserted within the \code{adjustbox} environment. By default the options are \code{width = 1\\textwidth, center} (if \code{TRUE}). A numeric value changes the value before \code{\\textwidth}. You can also add a character of the form \code{"x tw"} or \code{"x th"} with \code{x} a number and where tw (th) stands for text-width (text-height). Finally any other character value is passed verbatim as an \code{adjustbox} option.
#' @param tabular (Tex only.) Character scalar equal to "normal" (default), "*" or "X". Represents the type of tabular environment to use: either \code{tabular}, \code{tabular*} or \code{tabularx}.
#' @param export Character scalar giving the path to a PNG file to be created, default is \code{NULL}. If provided, the Latex table will be converted to PNG and copied to the \code{export} location. Note that for this option to work you need a working distribution of \code{pdflatex}, \code{imagemagick} and \code{ghostscript}.
#' @param markdown Character scalar giving the location of a directory, or a logical scalar. Default is \code{NULL}. This argument only works in Rmarkdown documents, when knitting the document. If provided: two behaviors depending on context. A) if the output document is Latex, the table is exported in Latex. B) if the output document is not Latex, the table will be exported to PNG at the desired location and inserted in the document via a markdown link. If equal to \code{TRUE}, the default location of the PNGs is a temporary folder for \code{R > 4.0.0}, or to \code{"images/etable/"} for earlier versions.
#' @param view.cache Logical, default is \code{FALSE}. Only used when \code{view = TRUE}. Whether the PNGs of the tables should be cached.
#' @param type Character scalar equal to 'pdflatex' (default), 'magick', 'dir' or 'tex'. Which log file to report; if 'tex', the full source code of the tex file is returned, if 'dir': the directory of the log files is returned.
#' @param highlight List containing coefficients to highlight. Highlighting is of the form \code{.("options1" = "coefs1", "options2" = "coefs2", etc)}.
#' The coefficients to be highlighted can be written in three forms: 1) row, eg \code{"x1"} will highlight the full row of the variable \code{x1}; 2) cells, use \code{'@'} after the coefficient name to give the column, it accepts ranges, eg \code{"x1@2, 4-6, 8"} will highlight only the columns 2, 4, 5, 6, and 8 of the variable \code{x1}; 3) range, by giving the top-left and bottom-right values separated with a semi-colon, eg \code{"x1@2 ; x3@5"} will highlight from the column 2 of \code{x1} to the 5th column of \code{x3}. Coefficient names are partially matched, use a \code{'\%'} first to refer to the original name (before dictionary) and use \code{'@'} first to use a regular expression. You can add a vector of row/cell/range.
#' The options are a comma-separated list of items. By default the highlighting is done with a frame (a thick box) around the coefficient, use \code{'rowcol'} to highlight with a row color instead. Here are the other options: \code{'se'} to highlight the standard-errors too; \code{'square'} to have a square box (instead of rounded); \code{'thick1'} to \code{'thick6'} to monitor the width of the box; \code{'sep0'} to \code{'sep9'} to monitor the inner spacing. Finally the remaining option is the color: simply add an R color (it must be a valid R color!). You can use \code{"color!alpha"} with "alpha" a number between 0 to 100 to change the alpha channel of the color.
#' @param coef.style Named list containing styles to be applied to the coefficients. It must be of the form \code{.("style1" = "coefs1", "style2" = "coefs2", etc)}. The style must contain the string \code{":coef:"} (or \code{":coef_se:"} to style both the coefficient and its standard-error). The string \code{:coef:} will be replaced verbatim by the coefficient value. For example use \code{"\\textbf{:coef:}"} to put the coefficient in bold. Note that markdown markup is enabled so \code{"**:coef:**"} would also put it in bold. The coefficients to be styled can be written in three forms: 1) row, eg \code{"x1"} will style the full row of the variable \code{x1}; 2) cells, use \code{'@'} after the coefficient name to give the column, it accepts ranges, eg \code{"x1@2, 4-6, 8"} will style only the columns 2, 4, 5, 6, and 8 of the variable \code{x1}; 3) range, by giving the top-left and bottom-right values separated with a semi-colon, eg \code{"x1@2 ; x3@5"} will style from the column 2 of \code{x1} to the 5th column of \code{x3}. Coefficient names are partially matched, use a \code{'\%'} first to refer to the original name (before dictionary) and use \code{'@'} first to use a regular expression. You can add a vector of row/cell/range.
#' @param page.width Character scalar equal to \code{'fit'} (default), \code{'a4'} or \code{'us'}; or a single Latex measure (like \code{'17cm'}) or a double one (like \code{"21, 2cm"}). Only used when the Latex table is to be viewed (\code{view = TRUE}), exported (\code{export != NULL}) or displayed in Rmarkdown (\code{markdown != NULL}). It represents the text width of the page in which the Latex table will be inserted. By default, \code{'fit'}, the page fits exactly the table (i.e. text width = table width). If \code{'a4'} or \code{'us'}, two times 2cm is removed from the page width to account for margins. Providing a page width and a margin width, like in \code{"17in, 1in"}, enables a correct display of the argument \code{adjustbox}. Note that the margin width represent the width of a single side margin (and hence will be doubled).
#' @param div.class Character scalar, default is \code{"etable"}. Only used in Rmarkdown documents when \code{markdown = TRUE}. The table in an image format is embedded in a \code{<div>} container, and that container is of class \code{div.class}.
#'
#'
#' @details
#' The function \code{esttex} is equivalent to the function \code{etable} with argument \code{tex = TRUE}.
#'
#' The function \code{esttable} is equivalent to the function \code{etable} with argument \code{tex = FALSE}.
#'
#' To display the table, you will need the Latex package \code{booktabs} which contains the \code{\\toprule}, \code{\\midrule} and \code{\\bottomrule} commands.
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
#' @section The argument \code{extralines}:
#'
#' The argument \code{extralines} adds well... extra lines to the table. It accepts either a list, or a one-sided formula.
#'
#' For each line, you can define the values taken by each cell using 4 different ways: a) a vector, b) a list, c) a function, and d) a formula.
#'
#' If a vector, it should represent the values taken by each cell. Note that if the length of the vector is smaller than the number of models, its values are recycled across models, but the length of the vector is required to be a divisor of the number of models.
#'
#' If a list, it should be of the form \code{list("item1" = #item1, "item2" = #item2, etc)}. For example \code{list("A"=2, "B"=3)} leads to \code{c("A", "A", "B", "B", "B")}. Note that if the number of items is 1, you don't need to add \code{= 1}. For example \code{list("A"=2, "B")} is valid and leads to \code{c("A", "A", "B"}. As for the vector the values are recycled if necessary.
#'
#' If a function, it will be applied to each model and should return a scalar (\code{NA} values returned are accepted).
#'
#' If a formula, it must be one-sided and the elements in the formula must represent either \code{extralines} macros, either fit statistics (i.e. valid types of the function \code{\link[fixest]{fitstat}}). One new line will be added for each element of the formula. To register \code{extralines} macros, you must first register them in \code{\link[fixest]{extralines_register}}.
#'
#' Finally, you can combine as many lines as wished by nesting them in a list. The names of the nesting list are the row titles (values in the leftmost cell). For example \code{extralines = list(~r2, Controls = TRUE, Group = list("A"=2, "B"))} will add three lines, the titles of which are "R2", "Controls" and "Group".
#'
#'
#' @section Controlling the placement of extra lines:
#'
#' The arguments \code{group}, \code{extralines} and \code{fixef.group} allow to add customized lines in the table. They can be defined via a list where the list name will be the row name. By default, the placement of the extra line is right after the coefficients (except for \code{fixef.group}, covered in the last paragraph). For instance, \code{group = list("Controls" = "x[[:digit:]]")} will create a line right after the coefficients telling which models contain the control variables.
#'
#' But the placement can be customized. The previous example (of the controls) will be used for illustration (the mechanism for \code{extralines} and \code{fixef.group} is identical).
#'
#' The row names accept 2 special characters at the very start. The first character tells in which section the line should appear: it can be equal to \code{"^"}, \code{"-"}, or \code{"_"}, meaning respectively the coefficients, the fixed-effects and the statistics section (which typically appear at the top, mid and bottom of the table). The second one governs the placement of the new line within the section: it can be equal to \code{"^"}, meaning first line, or \code{"_"}, meaning last line.
#'
#' Let's have some examples. Using the previous example, writing \code{"_^Controls"} would place the new line at the top of the statistics section. Writing \code{"-_Controls"} places it as the last row of the fixed-effects section; \code{"^^Controls"} at the top row of the coefficients section; etc...
#'
#' The second character is optional, the default placement being in the bottom. This means that \code{"_Controls"} would place it at the bottom of the statistics section.
#'
#' The placement in \code{fixef.group} is defined similarly, only the default placement is different. Its default placement is at the top of the fixed-effects section.
#'
#' @section Escaping special Latex characters:
#'
#' By default on all instances (with the notable exception of the elements of \code{\link[fixest]{style.tex}}) special Latex characters are escaped. This means that \code{title="Exports in million $."} will be exported as "Exports in million \\$.": the dollar sign will be escaped. This is true for the following characters: &, $, \%, _, ^ and #.
#'
#' Note, importantly, that equations are NOT escaped. This means that \code{title="Functional form $a_i \\times x^b$, variation in \%."} will be displayed as: \code{"Functional form $a_i \\times x^b$, variation in \\\%."}: only the last percentage will be escaped.
#'
#' If for some reason you don't want the escaping to take place, the arguments \code{headers} and \code{extralines} are the only ones allowing that. To disable escaping, add the special token ":tex:" in the row names. Example: in \code{headers=list(":tex:Row title"="weird & & \%\\n tex stuff\\\\")}, the elements will be displayed verbatim. Of course, since it can easily ruin your table, it is only recommended to super users.
#'
#' @section Markdown markup:
#'
#' Within anything that is Latex-escaped (see previous section), you can use a markdown-style markup to put the text in italic and/or bold. Use \code{*text*}, \code{**text**} or \code{***text***} to put some text in, respectively, italic (with \code{\\textit}), bold (with \code{\\textbf}) and italic-bold.
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
#'
#' est1 = feols(Ozone ~ i(Month) / Wind + Temp, data = airquality)
#' est2 = feols(Ozone ~ i(Month, Wind) + Temp | Month, data = airquality)
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
#' # signif.code
#' #
#'
#' etable(est1, est2, signif.code = c(" A"=0.01, " B"=0.05, " C"=0.1, " D"=0.15, " F"=1))
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
#' # Group and extralines
#' #
#'
#' # Sometimes it's useful to group control variables into a single line
#' # You can achieve that with the group argument
#'
#' setFixest_fml(..ctrl = ~ poly(Wind, 2) + poly(Temp, 2))
#' est_c0 = feols(Ozone ~ Solar.R, data = airquality)
#' est_c1 = feols(Ozone ~ Solar.R + ..ctrl, data = airquality)
#' est_c2 = feols(Ozone ~ Solar.R + Solar.R^2 + ..ctrl, data = airquality)
#'
#' etable(est_c0, est_c1, est_c2, group = list(Controls = "poly"))
#'
#' # 'group' here does the same as drop = "poly", but adds an extra line
#' # with TRUE/FALSE where the variables were found
#'
#' # 'extralines' adds an extra line, where you can add the value for each model
#' est_all  = feols(Ozone ~ Solar.R + Temp + Wind, data = airquality)
#' est_sub1 = feols(Ozone ~ Solar.R + Temp + Wind, data = airquality,
#'                  subset = ~ Month %in% 5:6)
#' est_sub2 = feols(Ozone ~ Solar.R + Temp + Wind, data = airquality,
#'                  subset = ~ Month %in% 7:8)
#' est_sub3 = feols(Ozone ~ Solar.R + Temp + Wind, data = airquality,
#'                  subset = ~ Month == 9)
#'
#' etable(est_all, est_sub1, est_sub2, est_sub3,
#'        extralines = list("Sub-sample" = c("All", "May-June", "Jul.-Aug.", "Sept.")))
#'
#' # You can monitor the placement of the new lines with two special characters
#' # at the beginning of the row name.
#' # 1) "^", "-" or "_" which mean the coefficients, the fixed-effects or the
#' # statistics section.
#' # 2) "^" or "_" which mean first or last line of the section
#' #
#' # Ex: starting with "_^" will place the line at the top of the stat. section
#' #     starting with "-_" will place the line at the bottom of the FEs section
#' #     etc.
#' #
#' # You can use a single character which will represent the section,
#' # the line would then appear at the bottom of the section.
#'
#' # Examples
#' etable(est_c0, est_c1, est_c2, group = list("_Controls" = "poly"))
#' etable(est_all, est_sub1, est_sub2, est_sub3,
#'        extralines = list("^^Sub-sample" = c("All", "May-June", "Jul.-Aug.", "Sept.")))
#'
#'
#' #
#' # headers
#' #
#'
#'
#' # You can add header lines with 'headers'
#' # These lines will appear at the top of the table
#'
#' # first, 3 estimations
#' est_header = feols(c(Ozone, Solar.R, Wind) ~  poly(Temp, 2), airquality)
#'
#' # header => vector: adds a line w/t title
#' etable(est_header, headers = c("A", "A", "B"))
#'
#' # header => list: identical way to do the previous header
#' # The form is: list(item1 = #item1, item2 = #item2,  etc)
#' etable(est_header, headers = list("A" = 2, "B" = 1))
#'
#' # Adding a title +
#' # when an element is to be repeated only once, you can avoid the "= 1":
#' etable(est_header, headers = list(Group = list("A" = 2, "B")))
#'
#' # To change the placement, add as first character:
#' # - "^" => top
#' # - "-" => mid (default)
#' # - "_" => bottom
#' # Note that "mid" and "top" are only distinguished when tex = TRUE
#'
#' # Placing the new header line at the bottom
#' etable(est_header, headers = list("_Group" = c("A", "A", "B"),
#'                                   "^Currency" = list("US $" = 2, "CA $" = 1)))
#'
#'
#' # In Latex, you can add "grouped underlines" (cmidrule from the booktabs package)
#' # by adding ":_:" in the title:
#' etable(est_header, tex = TRUE,
#'        headers = list("^:_:Group" = c("A", "A", "B")))
#'
#' #
#' # extralines and headers: .() for list()
#' #
#'
#' # In the two arguments extralines and headers, .() can be used for list()
#' # For example:
#' etable(est_header, headers = .("^Currency" = .("US $" = 2, "CA $" = 1)))
#'
#'
#'
#' #
#' # fixef.group
#' #
#'
#' # You can group the fixed-effects line with fixef.group
#'
#' est_0fe = feols(Ozone ~ Solar.R + Temp + Wind, airquality)
#' est_1fe = feols(Ozone ~ Solar.R + Temp + Wind | Month, airquality)
#' est_2fe = feols(Ozone ~ Solar.R + Temp + Wind | Month + Day, airquality)
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
#' # Using customized placement => as with 'group' and 'extralines',
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
#' # You can use external functions to compute the VCOVs
#' # by feeding functions in the 'vcov' argument.
#' # Let's use some covariances from the sandwich package
#'
#' etable(est_c0, est_c1, est_c2, vcov = sandwich::vcovHC)
#'
#' # To add extra arguments to vcovHC, you need to write your wrapper:
#' etable(est_c0, est_c1, est_c2, vcov = function(x) sandwich::vcovHC(x, type = "HC0"))
#'
#'
#' #
#' # Customize which fit statistic to display
#' #
#'
#' # You can change the fit statistics with the argument fitstat
#' # and you can rename them with the dictionary
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
#' est = feols(Ozone ~ Solar.R + Wind + Temp, data = airquality)
#'
#' #
#' # Method 1: use summary
#'
#' s1 = summary(est, "iid")
#' s2 = summary(est, cluster = ~ Month)
#' s3 = summary(est, cluster = ~ Day)
#' s4 = summary(est, cluster = ~ Day + Month)
#'
#' etable(list(s1, s2, s3, s4))
#'
#' #
#' # Method 2: using a list in the argument 'vcov'
#'
#' est_bis = feols(Ozone ~ Solar.R + Wind + Temp | Month, data = airquality)
#' etable(est, est_bis, vcov = list("hetero", ~ Month))
#'
#' # When you have only one model, this model is replicated
#' # along the elements of the vcov list.
#' etable(est, vcov = list("hetero", ~ Month))
#'
#' #
#' # Method 3: Using "each" or "times" in vcov
#'
#' # If the first element of the list in 'vcov' is "each" or "times",
#' # then all models will be replicated and all the VCOVs will be
#' # applied to each model. The order in which they are replicated
#' # are governed by the each/times keywords.
#'
#'
#' # each
#' etable(est, est_bis, vcov = list("each", "iid", ~ Month, ~ Day))
#'
#' # times
#' etable(est, est_bis, vcov = list("times", "iid", ~ Month, ~ Day))
#'
#' #
#' # Notes and markup
#' #
#'
#' # Notes can be also be set in a dictionary
#' # You can use markdown markup to put text into italic/bold
#'
#' dict = c("note 1" = "*Notes:* This data is not really random.",
#'          "source 1" = "**Source:** the internet?")
#'
#' est = feols(Ozone ~ csw(Solar.R, Wind, Temp), data = airquality)
#'
#' etable(est, dict = dict, tex = TRUE, notes = c("note 1", "source 1"))
#'
#'
#'
etable = function(..., vcov = NULL, stage = 2, agg = NULL,
                  se = NULL, ssc = NULL, cluster = NULL,
                  .vcov = NULL, .vcov_args = NULL, digits = 4, digits.stats = 5, tex,
                  fitstat = NULL, title = NULL, coefstat = "se", ci = 0.95,
                  se.row = NULL, se.below = NULL,
                  keep = NULL, drop = NULL, order = NULL,
                  dict = TRUE, file = NULL, replace = FALSE, convergence = NULL,
                  signif.code = NULL, label = NULL, float = NULL,
                  headers = list("auto"), fixef_sizes = FALSE,
                  fixef_sizes.simplify = TRUE, keepFactors = TRUE,
                  family = NULL, powerBelow = -5,
                  interaction.combine = NULL, interaction.order = NULL,
                  i.equal = NULL, depvar = TRUE, style.tex = NULL,
                  style.df = NULL, notes = NULL, group = NULL, extralines = NULL,
                  fixef.group = NULL, placement = "htbp", drop.section = NULL,
                  poly_dict = c("", " square", " cube"), postprocess.tex = NULL,
                  postprocess.df = NULL, tpt = FALSE, arraystretch = NULL, adjustbox = NULL,
                  fontsize = NULL, fit_format = "__var__", coef.just = NULL,
                  tabular = "normal", highlight = NULL, coef.style = NULL,
                  meta = NULL, meta.time = NULL, meta.author = NULL, meta.sys = NULL,
                  meta.call = NULL, meta.comment = NULL, view = FALSE,
                  export = NULL, markdown = NULL, page.width = "fit",
                  div.class = "etable"){

    #
    # Checking the arguments
    #

    # Need to check for the presence of the se
    useSummary = TRUE
    if(missnull(vcov) && missnull(se) && missnull(cluster) && missing(.vcov) && missing(stage) && missnull(agg)){
        useSummary = FALSE
    }

    # The depvar
    if(missing(depvar) && !missing(file)){
        depvar = TRUE
    }

    check_arg(tex, view, "logical scalar")
    if(missing(tex)){
        if(!missing(file)) {
            tex = TRUE
        } else {
            tex = FALSE
        }
    }

    # Float or not
    check_arg(float, "NULL logical scalar")
    if(missnull(float)){
        if(!missing(title) || !missing(label)){
            float = TRUE
        } else {
            float = FALSE
        }
    } else if(!float && (!missnull(title) || !missnull(label))) {
        what = c("title", "label")[c(!missing(title), !missing(label))]
        warning("Since float = FALSE, the argument", enumerate_items(what, "s.is"), " ignored.",
                immediate. = TRUE, call. = FALSE)
    }

    check_value(div.class, "character scalar")

    # NOTA: now that I allow the use of .(stuff) for headers and extralines
    # list(...) will raise an error if subtitle (now deprec) is used with .()
    # And there's no clear error message => I think it's a fair limitation
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

    #
    # Deprecated items
    #

    if("subtitles" %in% names(dots)){
        if(is.null(getOption("fixest_etable_arg_subtitles"))){
            warning("The argument 'subtitles' is deprecated. Please use 'headers' instead.")
            options(fixest_etable_arg_subtitles = TRUE)
        }
        headers = dots$subtitles
        dots$subtitles = NULL
    }

    if("extraline" %in% names(dots)){
        if(is.null(getOption("fixest_etable_arg_extraline"))){
            warning("The argument 'extraline' is deprecated. Please use 'extralines' instead (note the last 's'!).")
            options(fixest_etable_arg_extraline = TRUE)
        }
        extralines = dots$extraline
        dots$extraline = NULL
    }

    if("sdBelow" %in% names(dots)){
        if(is.null(getOption("fixest_etable_arg_sdBelow"))){
            warning("The argument 'sdBelow' is deprecated. Please use 'se.below' instead.")
            options(fixest_etable_arg_sdBelow = TRUE)
        }
        se.below = dots$sdBelow
        dots$sdBelow = NULL
    }

    if("signifCode" %in% names(dots)){
        if(is.null(getOption("fixest_etable_arg_signifCode"))){
            warning("The argument 'signifCode' is deprecated. Please use 'signif.code' instead.")
            options(fixest_etable_arg_signifCode = TRUE)
        }
        signif.code = dots$signifCode
        dots$signifCode = NULL
    }

    # Getting the model names
    if(.up == 2){
        # it's pain in the necky
        sysOrigin = sys.parent()
        mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
        dots_call = mc[["..."]]
    } else {
        mc = match.call(expand.dots = FALSE)
        dots_call = mc[["..."]]
    }

    # vcov
    vcov = oldargs_to_vcov(se, cluster, vcov, .vcov)

    if(is_function_in_it(vcov) && missnull(.vcov_args)){
        vcov_fun = if(is.function(vcov)) vcov else vcov[[1]]
        .vcov_args = catch_fun_args(vcov_fun, dots, exclude_args = "vcov", erase_args = TRUE)
        for(var in intersect(names(.vcov_args), names(dots))) dots[[var]] = NULL
    }


    # Arguments that can be set globally
    opts = getOption("fixest_etable")

    args_global = c("postprocess.tex", "postprocess.df", "view", "markdown", "page.width")
    for(arg in setdiff(args_global, names(mc))){
        if(arg %in% names(opts)){
            assign(arg, opts[[arg]])
        }
    }

    # argument only in setFixest_etable: cache
    cache = opts$view.cache

    check_arg(markdown, "NULL scalar(logical, character)")
    if(is.logical(markdown)){
        # R > 4.0.0: we use R_user_dir() to store the files across sessions
        # else images/etable/

        # The bookkeeping is handled in the dedicated function build_tex_png
        if(isFALSE(markdown)) markdown = NULL
    }

    tex_origin = tex
    is_md = !is.null(markdown)
    if(is_md){
        if(!is_Rmarkdown()){
            if("markdown" %in% names(mc)){
                warning("The argument 'markdown' only works when knitting Rmarkdown documents. It is currently ignored.")
            }
            markdown = NULL
            is_md = FALSE
        } else {
            tex = TRUE
            export = NULL
            view = FALSE
            if(knitr::is_latex_output()){
                is_md = FALSE
            }
        }
    }

    is_export = !is.null(export)

    # We can display the table in the R console AND at the same time showing
    # the latex table in the viewer
    do_df = isFALSE(tex)
    do_tex = tex || view || is_export

    #
    # postprocess.tex
    #

    # Note that I need to catch the arguments from the two pp functions
    # => so that the same call lead to valid evaluations

    check_arg(postprocess.tex, "NULL function arg(1,)")
    pp_tex = postprocess.tex
    is_pp_tex = !is.null(pp_tex)

    check_arg(postprocess.df, "NULL function arg(1,)")
    pp_df = postprocess.df
    is_pp_df = !is.null(pp_df)


    # Catching the arguments
    pp_tex_args = pp_df_args = list()
    if(!is.null(names(dots))){
        if(is_pp_tex){
            fm = formalArgs(pp_tex)
            qui = names(dots) %in% fm
            pp_tex_args = dots[qui]
            dots = dots[!qui]
        }

        if(is_pp_df){
            fm = formalArgs(pp_df)
            qui = names(dots) %in% fm
            pp_df_args = dots[qui]
            dots = dots[!qui]
        }
    }

    if(length(dots) == 0){
        stop("After cleaning the arguments to the postprocessing functions, there is no element left in '...'. Please provide at least one.")
    }

    for(i in seq_along(dots)){
        if(!is_fixest_model(dots[[i]]) && !(is.list(dots[[i]]) && is_fixest_model(dots[[i]][[1]]))){
            msg = ""
            if(!is.null(names(dots))){
                msg = paste0(" (named '", names(dots)[i], "')")
            }
            stop("The ", n_th(i), " element of '...'", msg, " is not valid: it should be a fixest object or a list of fixest objects, it is neither.")
        }
    }

    # internal argument
    caption.number = TRUE

    build_etable_list = function(TEX){
        results2formattedList(
            dots = dots, vcov = vcov, ssc = ssc, fitstat_all = fitstat,
            stage = stage, agg = agg,
            .vcov_args = .vcov_args, digits = digits, digits.stats = digits.stats,
            se.row = se.row, se.below = se.below,
            signif.code = signif.code, coefstat = coefstat,
            ci = ci, title = title, float = float, headers = headers,
            keepFactors = keepFactors, tex = TEX, useSummary = useSummary,
            dots_call = dots_call, powerBelow = powerBelow, dict = dict,
            interaction.combine = interaction.combine, interaction.order = interaction.order,
            i.equal = i.equal, convergence = convergence,
            family = family, keep = keep, drop = drop, file = file, order = order,
            label = label, fixef_sizes = fixef_sizes,
            fixef_sizes.simplify = fixef_sizes.simplify,
            depvar = depvar, style.tex = style.tex, style.df = style.df,
            replace = replace, notes = notes, group = group, extralines = extralines,
            fixef.group = fixef.group, placement = placement,
            drop.section = drop.section, poly_dict = poly_dict,
            tex_tag = TRUE, fit_format = fit_format,
            coef.just = coef.just, meta = meta, meta.time = meta.time,
            meta.author = meta.author, meta.sys = meta.sys,
            meta.call = meta.call, meta.comment = meta.comment,
            tpt = tpt, arraystretch = arraystretch, adjustbox = adjustbox,
            fontsize = fontsize, tabular = tabular, highlight = highlight,
            coef.style = coef.style, caption.number = caption.number,
            mc = mc, .up = .up + 1)
    }


    # Checking the requirements for exports to png
    if(view || is_md || is_export){
        build_ok = check_build_available()
        if(!isTRUE(build_ok)){
            args = c("view", "export", "markdown")
            what = args[args %in% names(mc)]
            if(length(what) > 0){
                warning("The argument", enumerate_items(what, "s.require.quote"), " ", build_ok, " to work. So they are currently ignored.")
            }

            view = is_md = is_export = FALSE
            markdown = export = NULL
            do_tex = tex = tex_origin
        }
    }

    if(do_tex){

        # if we're in markdown, we remove the table number
        if(is_md){
            caption.number = FALSE
        }

        info_tex = build_etable_list(TRUE)
        res_tex = etable_internal_latex(info_tex)
        n_models = attr(res_tex, "n_models")
        attr(res_tex, "n_models") = NULL
    }

    if(do_df){
        info_df = build_etable_list(FALSE)
        res_df = etable_internal_df(info_df)
    }

    make_png = function(x) NULL
    is_png = view || is_md || is_export
    if(is_png){
       make_png = function(x) build_tex_png(x, view = view, export = export,
                                            markdown = markdown, cache = cache,
                                            page.width = page.width)
    }

    # DF requested, but also png, OK user
    if(isFALSE(tex) && is_png){
        # Some exporting functions only print stuff!
        # I need to take care of that! PAIN IN THE NECK!

        ok = TRUE
        if(is_pp_tex){
            pp_tex_args_all = list(res_tex)
            if(length(pp_tex_args) > 0){
                pp_tex_args_all[names(pp_tex_args)] = pp_tex_args
            }

            tex_output = capture.output(res_tex <- do.call(pp_tex, pp_tex_args_all))
            if(length(tex_output) > 0){
                ok = FALSE
                args = c("view", "export")[c(view, is_export)]
                warning("The argument", enumerate_items(args, "s.don't.quote"),
                        " work with the current postprocessing function for tex. Try to remove it first.")
            }
        }

        if(ok) make_png(res_tex)
    }

    # Export to file
    is_file = !missnull(file)
    if(is_file){
        error_sender(sink(file = file, append = !replace),
                     "Argument 'file': error when creating the document in ", file)

        on.exit(sink())
    }


    # Output to the requested format
    if(tex){

        if(is_pp_tex){
            pp_tex_args_all = list(res_tex)
            if(length(pp_tex_args) > 0){
                pp_tex_args_all[names(pp_tex_args)] = pp_tex_args
            }

            res_tex = do.call(pp_tex, pp_tex_args_all)
        }

        if(is.null(res_tex)){
            if(is_png){
                args = c("view", "export", "markdown")[c(view, is_export, is_md)]
                warning("The argument", enumerate_items(args, "s.don't.quote"),
                        " work with the current postprocessing function for tex. Try to remove it first.")
            }

            return(invisible(NULL))
        }

        res_tex = res_tex[!res_tex %in% c("%start:tab\n", "%end:tab\n")]

        res_tex = tex.nice(res_tex, n_models) # we wait after PP to nicify

        path = make_png(res_tex)

        if(is_md){
            if(!knitr::is_latex_output()){
                path = path_to_relative(path)
                cat(.dsb0('<div class = ".[div.class]"><img src = ".[path]"></div>\n'))
                return(invisible(NULL))
            }
        }

        if(is_file){
            # We add extra whitespaces => otherwise the Latex file is a bit cluttered
            cat("\n")
            cat(res_tex, sep = "\n")
            cat("\n\n")
        }

        if(identical(class(res_tex), "character")){
            class(res_tex) = "etable_tex"
        }

        if(is_file){
            return(invisible(res_tex))
        } else {
            return(res_tex)
        }


    } else {
        # Some postptocessing functions return nothing (just print)
        # other return a character vector => I need to take care of that

        if(is_pp_df){
            pp_df_args_all = list(res_df)
            if(length(pp_df_args) > 0){
                pp_df_args_all[names(pp_df_args)] = pp_df_args
            }

            res_df = do.call(pp_df, pp_df_args_all)
        }

        if(is.null(res_df)){
            return(invisible(NULL))
        }

        if(is_file){
            if(is.data.frame(res_df)){
                print(res_df)
            } else {
                cat(res_df)
            }
            return(invisible(res_df))
        } else {
            if(is.data.frame(res_df)){
                return(res_df)
            } else {
                cat(res_df)
                return(invisible(res_df))
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

    qui_df = !arg_name %in% c("tex", "title", "label", "float", "style.tex",
                              "notes", "placement", "postprocess.tex",
                              "meta", "meta.time", "meta.author", "meta.sys",
                              "meta.call", "meta.comment", "tpt", "arraystretch",
                              "adjustbox", "fontsize", "view", "tabular", "markdown")

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
    intro = c("# Do not edit by hand\n# => aliases to the function etable\n\n\n")

    s = "\n\n\n\n"
    text = c(intro, s, esttable_rox, esttable_fun, s, esttex_rox, esttex_fun, s)

    update_file("R/etable_aliases.R", text)
}

results2formattedList = function(dots, vcov = NULL, ssc = getFixest_ssc(), stage = 2,
                                 agg = NULL, .vcov_args = NULL, digits = 4,
                                 digits.stats = 5, fitstat_all, se.row = NULL, se.below = NULL, dict,
                                 signif.code = c("***"=0.01, "**"=0.05, "*"=0.10),
                                 coefstat = "se", ci = 0.95, label, headers, title,
                                 float = FALSE, replace = FALSE, keepFactors = FALSE,
                                 tex = FALSE, useSummary, dots_call, powerBelow = -5,
                                 interaction.combine, interaction.order, i.equal,
                                 convergence, family, drop, order,
                                 keep, file, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE,
                                 depvar = FALSE, style.tex = NULL, style.df=NULL,
                                 notes = NULL, group = NULL, extralines=NULL,
                                 fixef.group = NULL, placement = "htbp", drop.section = NULL,
                                 poly_dict = c("", " square", " cube"), tex_tag = FALSE,
                                 fit_format = "__var__", coef.just = NULL,
                                 meta = NULL, meta.time = NULL, meta.author = NULL,
                                 meta.call = NULL, meta.sys = NULL, meta.comment = NULL,
                                 tpt = FALSE, arraystretch = NULL, adjustbox = NULL,
                                 fontsize = NULL, tabular = "normal", highlight = NULL,
                                 coef.style = NULL, caption.number = TRUE,
                                 mc = NULL, .up = 1){


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
    sysOrigin = sys.parent(.up)
    if(length(opts) > 0){
        args_usr = setdiff(names(mc), c("style.tex", "style.df"))

        # We modify only non-user provided arguments
        for(v in names(opts)){
            if(!v %in% args_usr){
                assign(v, opts[[v]])
            }
        }
    }

    # Getting the default style values
    if(tex){
        if(!"style.tex" %in% names(opts)){
            style = fixest::style.tex(main = "base")
        } else {
            style = style.tex
        }
    } else if(!tex){
        if(!"style.df" %in% names(opts)){
            style = fixest::style.df(default = TRUE)
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

    # Arguments both in style AND in etable
    args_dual = c("tpt", "arraystretch", "fontsize", "adjustbox",
                  "tabular", "signif.code")
    for(arg in setdiff(args_dual, names(mc))){
        # We set to the default in style (only if NOT user-provided)
        if(arg %in% names(style)){
            assign(arg, style[[arg]])
        }
    }

    #
    # Full control
    #

    check_arg(title, "NULL character scalar")
    check_arg_plus(coefstat, "match(se, tstat, confint)")

    check_arg_plus(notes, "NULL character vector no na")
    if(length(notes) > 0) notes = notes[nchar(notes) > 0]

    check_arg("logical scalar", replace, fixef_sizes, fixef_sizes.simplify,
              keepFactors, tex, depvar)

    check_arg("NULL logical scalar", convergence, family, se.below, se.row, tpt)
    if(is.null(tpt)) tpt = FALSE

    isTex = tex
    if(missnull(family)){
        show_family = NULL
    } else {
        show_family = family
    }

    if(is.null(se.below)){
        se.below = isTex
    }

    # start: coef.just
    check_arg(coef.just, "NULL charin", .choices = c(".", "(", "l", "c", "r"))

    if(!isTex){
        # this parameter is only used in DF
        if(is.null(coef.just) && se.below){
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

    check_arg(keep, drop, order, "character vector no na NULL",
              .message = "The arg. '__ARG__' must be a vector of regular expressions (see help(regex)).")

    check_arg(file, label, interaction.combine, i.equal, "NULL character scalar")

    check_arg_plus(tabular, "match(normal, *, X)")

    # interaction.combine: Default depends on type
    if(is.null(interaction.combine)){
        # interaction.combine = if(isTex) " $\\times $ " else " x "
        interaction.combine = style$interaction.combine
    }

    # i.equal
    if(is.null(i.equal)){
        # i.equal = if(isTex) " $=$ " else " = "
        i.equal = style$i.equal
    }

    check_arg_plus(signif.code, "NULL NA | match(letters) | named numeric vector no na GE{0} LE{1}")

    if(missnull(signif.code)){
        if(isTex){
            signif.code = c("***"=0.01, "**"=0.05, "*"=0.10)
        } else {
            signif.code = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)
        }
    }

    add_signif = TRUE
    if(identical(signif.code, "letters")){
        signif.code = c("a"=0.01, "b"=0.05, "c"=0.10)
    } else if(length(signif.code) == 0 || (length(signif.code) == 1 && is.na(signif.code))){
        add_signif = FALSE
    } else {
        signif.code = sort(signif.code)
    }

    # eval_dot: headers + extralines + highlight
    headers = error_sender(eval_dot(headers, .up + 1), arg_name = "headers", up = .up)
    extralines = error_sender(eval_dot(extralines, .up + 1), arg_name = "extralines", up = .up)
    highlight = error_sender(eval_dot(highlight, .up + 1), arg_name = "highlight", up = .up)
    coef.style = error_sender(eval_dot(coef.style, .up + 1), arg_name = "coef.style", up = .up)

    check_arg(headers, "NULL character vector no na | NA | list")
    if(is.null(headers)) headers = list()

    check_arg(highlight, "NULL character vector no na | list")
    if(is.character(highlight)){
        # We always want HL to be a list!
        # hl = "Petal.L>1-2"
        # hl = c("rowcol, lightblue" = "Petal.L>1-2")
        if(length(highlight) == 1){
            highlight = as.list(highlight)
        } else {
            # length 2 => can't have names
            # hl = c("Petal.L@1", "Sepal.L@3")
            highlight = list(highlight)
        }
    }

    check_arg(coef.style, "NULL named list")
    if(is.null(coef.style)){
        coef.style = list()
    }


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
    check_arg_plus(extralines, "NULL{list()} list l0 | os formula | vector")
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

    check_arg(arraystretch, "NULL numeric scalar GT{0}")

    check_arg_plus(fontsize, "NULL match(tiny, scriptsize, footnotesize, small, normalsize, large, Large)")

    # adjustbox + default
    adjustbox = check_set_adjustbox(adjustbox, .up)

    #
    # meta
    #

    meta_txt = ""
    if(isTex){
        check_arg(meta, "NULL os formula")
        # Note that meta CANNOT be set globally, so if not null it has **always** precedence!
        meta_order = c("time", "author", "sys", "call", "comment")
        if(!is.null(meta)){
            # meta = ~time + author("hahaha") + sys + comment("this is nuts")
            meta_vars = all.vars(meta, functions = TRUE)
            meta_vars = meta_vars[grepl("[[:alpha:]]{2,}", meta_vars)]

            check_value_plus(meta_vars, "multi match(time, date, sys, author, comment, call)",
                             .message = paste0("The argument 'meta' must be a one-sided formula composed of the following variables: \n",
                             "  time, date, author, sys, call or comment."))

            # date/time/call/sys cannot be directly edited
            if("date" %in% meta_vars) meta.time = "date"
            if("time" %in% meta_vars) meta.time = "time"
            if("call" %in% meta_vars) meta.call = TRUE
            if("sys" %in% meta_vars) meta.sys = TRUE

            meta_vars[meta_vars == "date"] = "time"
            meta_order = order_apply(meta_order, paste0("^", meta_vars, "$"))

            meta_terms = attr(terms(meta), "term.labels")
            if("author" %in% meta_vars) {
                if(!"author" %in% all.vars(meta, functions = FALSE)){
                    # means author is a function
                    my_term = meta_terms[grepl("^aut.*\\(", meta_terms)]
                    my_term = gsub("^[^\\(]+\\(", "I(", my_term)
                    meta.author = eval(str2lang(my_term))
                } else {
                    meta.author = TRUE
                }
            }

            if("comment" %in% meta_vars){
                if(!"comment" %in% all.vars(meta, functions = FALSE)){
                    my_term = meta_terms[grepl("^com.*\\(", meta_terms)]
                    my_term = gsub("^[^\\(]+\\(", "I(", my_term)
                    meta.comment = eval(str2lang(my_term))
                } else {
                    # Comments MUST be provided explicitly
                    # => we rely on the default comment then (if provided)
                }
            }
        }


        check_arg_plus(meta.time, "NULL match(date, time) | logical scalar")
        check_arg(meta.call, meta.sys, "NULL logical scalar")
        check_arg(meta.author, "NULL logical scalar | character vector no na")
        check_arg(meta.comment, "NULL character vector no na")

        for(meta_value in meta_order){

            if(meta_value == "time" && !is.null(meta.time)){
                if(isTRUE(meta.time)) meta.time = "time"
                my_date = switch(meta.time, "date" = Sys.Date(), "time" = Sys.time())
                meta_txt = paste0(meta_txt, "% Created: ", my_date, "\n")

            } else if(meta_value == "sys" && isTRUE(meta.sys)){
                my_sys = Sys.info()[1:5]
                meta_txt = paste0(meta_txt, "% System info: ",
                                  paste0(names(my_sys), ": ", my_sys, collapse = ", "), "\n")

            } else if(meta_value == "call" && isTRUE(meta.call)){
                sc = sys.call(sysOrigin)
                my_call = paste0("% ", deparse(sc), collapse = "\n")
                meta_txt = paste0(meta_txt, "% Call:\n", my_call, "\n")

            } else if(meta_value == "author" && !is.null(meta.author) && !isFALSE(meta.author)){

                if(isTRUE(meta.author)){

                    user_all = Sys.info()[c("login", "user", "effective_user")]

                    author_txt = "% Created by: "
                    if(length(unique(user_all)) == 1){
                        # Windows
                        author_txt = paste0(author_txt, user_all[1])
                    } else {
                        # Unix
                        user_login = user_all["login"]
                        user_no_login = user_all[-1]

                        if(length(unique(user_all)) == 2){
                            if(user_login %in% user_no_login){
                                i = which(user_no_login == user_login)
                                user_login = paste0(user_login, " (login, ", names(user_no_login)[i], ")")
                                user_rest = paste0(user_no_login[-i], " (", names(user_no_login)[-i], ")")
                            } else {
                                user_login = paste0(user_login, " (login)")
                                user_rest = paste0(user_all[2], "(user, effective_user)")
                            }

                            author_txt = paste0(author_txt, user_login, ", ", user_rest)
                        } else {
                            author_txt = paste0(author_txt, user_all[1], " (login), ", user_all[2], "(user), ", user_all[3], " (effective_user)")
                        }
                    }

                } else {
                    author_txt = meta.author

                    if(length(author_txt) > 1){
                        author_txt = paste0(author_txt, collapse = "\n%             ")
                    }

                    if(!grepl("^%", author_txt)){
                        author_txt = paste0("% Created by: ", author_txt)
                    }
                }

                meta_txt = paste0(meta_txt, author_txt, "\n")

            } else if(meta_value == "comment" && !is.null(meta.comment)){

                if(length(meta.comment) > 1){
                    meta.comment = paste0(meta.comment, collapse = "\n% ")
                }

                if(!grepl("^%", meta.comment)){
                    meta.comment = paste0("% ", meta.comment)
                }

                meta_txt = paste0(meta_txt, meta.comment, "\n")
            }
        }
    }


    # yesNo => used later when setting the fixed-effect line
    yesNo = style$yesNo

    # default values for dict
    dict_global = getFixest_dict()
    if(missing(dict) || isTRUE(dict)) {
        dict = dict_global
    } else if(isFALSE(dict)) {
        dict = NULL
    } else {
        # dict changes the values set in the global dict

        if(dict[1] == "reset"){
            dict_global = c()
            dict = dict[-1]
        }

        if(length(dict) > 0){
            dict_global[names(dict)] = dict
        }

        dict = dict_global
    }

    # headers => must be a list
    # We get the automatic headers, if split is used
    AUTO_HEADERS = FALSE
    i_auto_headers = 1
    if(is.list(headers)){
        if(length(headers) > 0){
            qui = sapply(headers, function(x) identical(x, "auto"))
            if(any(qui)){
                i_auto_headers = which(qui)[1]
                headers = headers[!qui]
                AUTO_HEADERS = TRUE
            }

            # We expand the headers if needed
            if(length(headers) > 0){
                # ex: headers = list(Gender = list("M"=2, "F"=2))

                if(length(headers[[1]]) == 1){
                    # ex: headers = list("M"=2, "F"=2)
                    # or list("M", "F" = 2)
                    if(is.numeric(headers[[1]]) || # case list("M" = 2, "F" = 3)
                       (!is.null(names(headers)) && nchar(names(headers)[1]) == 0 && is.character(headers[[1]]))){ # case list("M", "F" = 2)
                        headers = list(headers)
                    }
                }

                for(i in seq_along(headers)){
                    # expand_list_vector: list("A"=2, "B") => c("A", "A", "B")
                    headers[[i]] = expand_list_vector(headers[[i]])
                }

                # We ensure headers have names
                if(is.null(names(headers))){
                    names(headers) = character(length(headers))
                }

            }
        }
    } else if(anyNA(headers)){
        headers = list()
    } else {
        # It's a character vector
        if(identical(headers, "auto")){
            headers = list()
            AUTO_HEADERS = TRUE
        } else {
            # we need names
            headers = list(headers)
            names(headers) = ""
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

    n_dots = length(dots)

    if(n_dots == 0) stop_up("Not any estimation as argument.")


    all_models = list()
    model_names = list()
    auto_headers = list()
    model_id = NULL
    k = 1
    for(i in 1:n_dots){
        di = dots[[i]]

        if("fixest_multi" %in% class(di)){
            meta = attr(di, "meta")
            di = attr(di, "data")

            if(AUTO_HEADERS){
                if(!is.null(meta$all_names[["sample"]])){
                    my_headers = list()

                    if(tex == FALSE){
                        # We need to rename to avoid FE problem duplicate row names
                        my_headers$title = paste0("Sample (", dict_apply(meta$all_names$split.name, dict), ")")
                    } else {
                        my_headers$title = dict_apply(meta$all_names$split.name, dict)
                    }

                    n_mod = length(di)

                    if(!"sample" %in% names(meta$index)){
                        my_headers$value = rep(meta$all_names$sample, n_mod)
                    } else {
                        my_headers$value = meta$all_names$sample[meta$tree[, "sample"]]
                    }

                    my_headers$index = seq(k, length.out = n_mod)

                    auto_headers[[length(auto_headers) + 1]] = my_headers

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
                if(n_dots > 1){
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

    IS_MULTI_VCOV = FALSE
    IS_EACH = FALSE
    if(!missnull(vcov) && identical(class(vcov), "list") && length(vcov) > 1){
        IS_MULTI_VCOV = TRUE

        vcov_1 = vcov[[1]]

        is_rep = identical(vcov_1, "times") || identical(vcov_1, "each")
        if(is_rep || n_models == 1){

            if(is_rep){
                vcov[[1]] = NULL
            }

            n_vcov = length(vcov)
            if(n_vcov == 0){
                stop_up("The argument 'vcov' is not valid.")
            }

            n_total = n_models * n_vcov

            if(vcov_1 == "times"){
                id_mod = rep(1:n_models, n_vcov)
                id_vcov = rep(1:n_vcov, each = n_models)
            } else {
                IS_EACH = TRUE
                id_mod = rep(1:n_models, each = n_vcov)
                id_vcov = rep(1:n_vcov, n_models)
            }

            mega_model_names = vector("list", n_total)
            mega_models = vector("list", n_total)
            mega_vcov = vector("list", n_total)

            for(i in 1:n_total){
                m_name = model_names[[id_mod[i]]]
                if(length(m_name) == 0 || m_name == ""){
                    mega_model_names[[i]] = paste0("model ", id_mod[i], ".", id_vcov[i])
                } else {
                    mega_model_names[[i]] = m_name
                }

                mega_models[[i]] = all_models[[id_mod[i]]]
                mega_vcov[[i]] = vcov[[id_vcov[i]]]
            }

            if(length(auto_headers) > 0){
                # I need to find out the new indexes... pain in the neck....
                if(vcov_1 == "times"){
                    for(i in seq_along(auto_headers)){
                        index = auto_headers[[i]]
                        n_index = length(index)
                        auto_headers[[i]] = rep(index, n_vcov) + n_models * rep(0:(n_vcov-1), each = n_index)
                    }
                } else {
                    for(i in seq_along(auto_headers)){
                        index = auto_headers[[i]]
                        new_index = (index - 1) * n_vcov + 1
                        n_index = length(index)

                        auto_headers[[i]] = rep(new_index, each = n_vcov) + rep(1:n_vcov, n_index)
                    }
                }
            }

            model_names = mega_model_names
            all_models = mega_models
            n_models = length(all_models)
            vcov = mega_vcov

        } else {
            check_value(vcov, "list len(value)", .value = n_models, .message = "If 'vcov' is a list, it must be of the same length as the number of models, or you should add the 'each' or 'times' keyword as the first element of the list.")
        }
    }

    auto_headers_clean = list()
    if(length(auto_headers) > 0){
        # We reconstruct the headers properly

        for(i in seq_along(auto_headers)){
            my_sub = auto_headers[[i]]
            if(my_sub$title %in% names(auto_headers_clean)){
                value = auto_headers_clean[[my_sub$title]]
            } else {
                value = character(n_models)
            }

            value[my_sub$index] = my_sub$value
            auto_headers_clean[[my_sub$title]] = value
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

    if(is_function_in_it(vcov)){
        # finding the name

        sysOrigin = sys.parent(.up)
        mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
        # must be either in .vcov (backward comp) or in vcov
        if(".vcov" %in% names(mc)){
            vcov_name = deparse_long(mc$.vcov)
        } else {
            vcov_name = deparse_long(mc$vcov)
        }

        if(missnull(.vcov_args)) .vcov_args = list()
        check_arg_plus(.vcov_args, "list", .message = "The argument '.vcov_args' must be a list of arguments to be passed to the function in '.vcov'.")
    }

    # If vcov is provided, we use summary
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
    # we'll use iv headers only if length(stage) > 1
    iv_sub = c()
    for(m in 1:n_models){
        if(useSummary){
            if(is_function_in_it(vcov)){
                x = summary(all_models[[m]], vcov = vcov, stage = stage, .vcov_args = .vcov_args, vcov_name = vcov_name, agg = agg)
            } else {
                if(IS_MULTI_VCOV){
                    x = summary(all_models[[m]], vcov = vcov[[m]], ssc = ssc, stage = stage, agg = agg)
                } else {
                    x = summary(all_models[[m]], vcov = vcov, ssc = ssc, stage = stage, agg = agg)
                }
            }

        } else {
            # What do we do if vcov not provided?
            # we apply summary only to the ones that are not summaries
            x = all_models[[m]]
            if(!isTRUE(x$summary)){
                # not a summary => we apply summary to trigger default behavior
                x = summary(x, ssc = ssc)
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
        # We've got multiple estimations from summary => we expand the headers
        all_models = all_models_bis

        model_names = rep(model_names, iv_times)
        i_max = rep(iv_times, iv_times)
        suffix = paste0(".", unlist(lapply(iv_times, seq)))
        suffix[i_max == 1] = ""
        model_names = paste0(model_names, suffix)

        if(length(auto_headers) > 0){
            # We expand the headers
            for(i in seq_along(auto_headers_clean)){
                auto_headers_clean[[i]] = rep(auto_headers_clean[[i]], iv_times)
            }
        }

        if(AUTO_HEADERS){
            my_stages = sapply(all_models, function(x) ifelse(is.null(x$iv_stage), 3, x$iv_stage))
            if(length(unique(my_stages)) > 1){
                dict_stage = c("First", "Second", " ")
                auto_headers_clean[["IV stages"]] = dict_stage[my_stages]
            }
        }

        n_models = length(all_models)
    }

    if(length(auto_headers_clean) > 0){
        # We place the auto headers in the right location

        auto_headers_format = list()
        for(i in seq_along(auto_headers_clean)){
            my_title = names(auto_headers_clean)[i]
            auto_headers_format[[dict_apply(my_title, dict)]] = dict_apply(auto_headers_clean[[i]], dict)
        }

        headers = insert(headers, auto_headers_format, i_auto_headers)
    }

    # if there are headers
    if(length(headers) == 0){
        isHeaders = FALSE
    } else {

        n_all = lengths(headers)

        if(any(n_models %% n_all != 0)){
            i_pblm = which(n_models %% n_all != 0)[1]
            info = if(length(n_all) == 1) "" else paste0(" (", n_th(i_pblm), " header)")
            stop_up("If argument 'headers' is provided, its elements must be of the same length as, or a divisor of, the number of models. Current lengths: ", n_all[i_pblm], " vs ", n_models, " models", info, ".")
        }


        for(i in which(n_all != n_models)){
            if(IS_EACH){
                headers[[i]] = rep(headers[[i]], each = n_models/n_all[i])
            } else {
                headers[[i]] = rep(headers[[i]], n_models/n_all[i])
            }
        }

        isHeaders = TRUE

        if(isTex){
            h_names = names(headers)
            for(i in seq_along(headers)){
                h_i = h_names[i]
                if(grepl(":tex:", h_i)){
                    # no escaping
                    h_i = sub(":tex:", "", h_i, fixed = TRUE)
                } else {
                    h_i = escape_latex(h_i)
                    h_i = gsub("^\\\\(\\^|\\_)", "\\1", h_i)
                    h_i = gsub(":\\_:", ":_:", h_i, fixed = TRUE)
                    headers[[i]] = escape_latex(headers[[i]])
                }
                h_names[i] = h_i
            }
            names(headers) = h_names
        } else {
            # we clean the possible Latex markup
            h_names = gsub(":_:|:tex:", "", names(headers))
            names(headers) = h_names
        }

    }

    #
    # We check the group and extralines arguments
    #

    if(missing(drop)) drop = NULL
    for(i in seq_along(group)){
        check_value(group[[i]], "character vector", .message = "The elements of argument 'group' must be character vectors of regular expressions.")
        drop = unique(c(drop, group[[i]]))
    }

    #
    # ... extralines ####
    #

    if(inherits(extralines, "formula")){
        # If a formula => to summon registered stats
        extralines = extralines_extractor(extralines, tex = isTex)

    } else if(!is.list(extralines)){
        # => vector
        extralines = list(extralines)

    } else if(length(extralines) > 0 && all(lengths(extralines) == 1) &&
              !inherits(extralines[[1]], "formula") && !is.function(extralines[[1]])){
        # list("A" = 2, "B")
        extralines = list(extralines)
    }


    el_new = list() # I need it to cope with list(~f+ivf+macro, "my vars" = TRUE)
    # => the first command will create several lines
    el_names = names(extralines)
    if(is.null(el_names)) el_names = character(length(extralines))
    el_names = uniquify_names(el_names)
    for(i in seq_along(extralines)){
        check_value(extralines[[i]], "vector | function | os formula | list", .message = paste0("The elements of argument 'extralines' must be vectors of length ", n_models, ", logical scalars, functions, one-sided formulas, or lists."),
                    .value = n_models)

        el = extralines[[i]]
        if("formula" %in% class(el)){
            el_tmp = extralines_extractor(el, el_names[i], tex = isTex)
            for(k in seq_along(el_tmp)){
                el_new[[names(el_tmp)[k]]] = el_tmp[[k]]
            }
        } else {

            if(!is.function(el)){

                if(is.list(el)){
                    el = expand_list_vector(el)
                }

                if(length(el) < n_models){
                    # we extend

                    n_el = length(el)
                    if(n_models %% n_el != 0){
                        stop_up("In 'extralines', the number of elements in the ", n_th(i), " line is not a divisor of the number of models (", n_el, " vs ", n_models, " models).")
                    }

                    if(IS_EACH){
                        el = rep(el, each = n_models/n_el)
                    } else {
                        el = rep(el, n_models/n_el)
                    }
                }
            }

            el_new[[el_names[i]]] = el
        }
    }

    extralines = el_new

    # Now we catch the functions + normalization of the names
    el_fun_id = NULL
    if(length(extralines) > 0){
        el_fun_id = which(sapply(extralines, is.function)) # set of ID such that el[id] is a function
        if(length(el_fun_id) > 0){
            el_origin = extralines
        }

        if(length(unique(names(extralines))) != length(extralines)){
            new_names = uniquify_names(names(extralines))
            names(extralines) = new_names
        }
    }


    # end: extralines

    # we keep track of the SEs
    se_type_list = list()

    check_interaction_reorder = FALSE
    var_list = var_reorder_list = coef_list = coef_below = sd_below = list()
    depvar_list = obs_list = fitstat_list = list()
    r2_list = aic_list = bic_list = loglik_list = convergence_list = list()
    sqCor_list = family_list = list()

    # To take care of factors
    fe_names = c()
    is_fe = vector(mode = "list", n_models)
    nb_fe = vector(mode = "list", n_models) # the number of items per factor

    slope_names = c()
    slope_flag_list = vector(mode = "list", n_models)

    if(!is.null(dict) && isTex){
        dict = escape_latex(dict)
    }

    #
    # ... fitstat ####
    #

    if("fitstat" %in% names(opts)){
        fitstat_all = opts$fitstat
    }

    if(missnull(fitstat_all)){
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
        depvar = gsub(" ", "", as.character(x$fml)[[2]])

        a = x$coeftable
        if(!is.data.frame(a)){
            class(a) = NULL
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

        #
        # Now we rename the variables
        #

        # on enleve les espaces dans les noms de variables
        var = var_origin = c(gsub(" ", "", row.names(a)))
        # renaming
        if(TRUE){
            # Now I clean white spaces in dict_apply
            qui = gsub(" ", "", var, fixed = TRUE) %in% gsub(" ", "", names(dict), fixed = TRUE)
            var = dict_apply(var, dict)
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

                                res[i] = paste(value_split, collapse = i.equal)
                            }
                        }

                        # We put the factors on the right
                        res = res[base::order(qui_factor)]

                    } else {
                        res = x
                    }

                    res = dict_apply(res, dict)
                    res = order_apply(res, interaction.order)

                    res = paste0(res, collapse = interaction.combine)

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
            qui_pow = grepl("^I\\([^\\^]+\\^[[:digit:]]+\\)$", var_left)
            if(any(qui_pow)){
                # We clean only I(var^d)
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

            # poly(xx, d)p => poly(xx)p
            if(any(grepl("^poly\\(", var_left))){

                qui_poly = grepl("^poly\\([^,]+,[[:digit:]]\\)", var_left)
                if(any(qui_poly)){
                    poly_new = gsub("^(poly\\([^,]+),[[:digit:]]\\)", "\\1)", var_left[qui_poly])
                    var[!qui][qui_poly] = new_var[!qui][qui_poly] = poly_new
                }

                # poly(var) => poly(var)1
                qui_poly_clean = grepl("^poly\\([^,]+\\)$", var[!qui])
                if(any(qui_poly_clean)){
                    poly_new = gsub("^(poly\\([^,]+\\))$", "\\11", var[!qui][qui_poly_clean])
                    var[!qui][qui_poly_clean] = new_var[!qui][qui_poly_clean] = poly_new
                }
            }


            names(new_var) = names(var) = var_origin
            var_reorder_list[[m]] = new_var

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

            var_reorder_list[[m]] = new_var
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
                pval = cut(a[, 4], breaks = c(-1, signif.code, 100), labels = c(tex_star(names(signif.code)), ""))
            } else {
                pval = cut(a[, 4], breaks = c(-1, signif.code, 100), labels = c(names(signif.code), ""))
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
        var_list[[m]] = var
        names(structured_coef) = var
        coef_list[[m]] = structured_coef
        if(se.below){
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
        depvar_list[[m]] = depvar

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
        # extralines, when function
        #

        for(i in el_fun_id){
            f = el_origin[[i]]

            if(m == 1){
                extralines[[i]] = f(x)

            } else {
                extralines[[i]] = c(extralines[[i]], f(x))
            }
        }


    }

    if(check_interaction_reorder){
        if(length(unique(unlist(var_reorder_list))) < length(unique(unlist(var_list)))){
            var_list = var_reorder_list
            for(m in 1:length(var_list)){
                names(coef_list[[m]]) = var_list[[m]]
                if(se.below){
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
        if(missnull(title)){
            title = "no title"
        } else {
            title = escape_latex(title, makecell = FALSE)
        }
    } else {
        if(missnull(title)){
            title = NULL
        }
    }


    if((missnull(convergence) && any(convergence_list == FALSE)) || (!missnull(convergence) && convergence)){
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

                } else {
                    # we add the default placement
                    el_name = paste0("-^", el_name)
                }

                fixef.extralines[[el_name]] = is_there

            }

            if(length(fixef.extralines) > 0){
                # We add it to extra lines
                extralines[names(fixef.extralines)] = fixef.extralines
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
                warn_up("In 'fixef.group', ", msg, deparse_long(fixef.group[[i]]), ").\nTo create inconsistent rows: use drop.section = 'fixef' combined with the arghument 'extralines'.")
            }

        }

    }

    res = list(se_type_list = se_type_list, var_list = var_list, coef_list = coef_list,
               coef_below = coef_below, se.row = se.row, sd_below = sd_below, depvar_list = depvar_list,
               obs_list = obs_list, convergence_list = convergence_list, fe_names = fe_names,
               is_fe = is_fe, nb_fe = nb_fe, slope_flag_list = slope_flag_list,
               slope_names = slope_names, useSummary = useSummary, model_names = model_names,
               family_list = family_list, fitstat_list = fitstat_list, headers = headers,
               isHeaders = isHeaders, title = title, convergence = convergence, family = family,
               keep = keep, drop = drop, order = order, file = file, label = label, se.below = se.below,
               signif.code = signif.code, fixef_sizes = fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify,
               depvar = depvar, useSummary = useSummary, dict = dict, yesNo = yesNo, add_signif = add_signif,
               float = float, coefstat = coefstat, ci = ci, style = style, notes = notes, group = group,
               extralines = extralines, placement = placement, drop.section = drop.section,
               tex_tag = tex_tag, fun_format = fun_format, coef.just = coef.just, meta = meta_txt,
               tpt = tpt, arraystretch = arraystretch, adjustbox = adjustbox, fontsize = fontsize,
               tabular = tabular, highlight = highlight, coef.style = coef.style,
               caption.number = caption.number)

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
    headers = info$headers
    isHeaders = info$isHeaders
    title = info$title
    label = info$label
    keep = info$keep
    drop = info$drop
    order = info$order
    file = info$file
    family = info$family
    convergence = info$convergence
    se.below = info$se.below
    signif.code = info$signif.code
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
    extralines = info$extralines
    placement = info$placement
    drop.section = info$drop.section
    tex_tag = info$tex_tag
    fun_format = info$fun_format
    meta = info$meta
    se.row = info$se.row
    tpt = info$tpt
    arraystretch = info$arraystretch
    adjustbox = info$adjustbox
    fontsize = info$fontsize
    tabular = info$tabular
    highlight = info$highlight
    coef.style = info$coef.style
    caption.number = info$caption.number

    # top line
    style$line.top = switch(style$line.top,
                            "simple" = "\\toprule",
                            "double" = "\\tabularnewline \\midrule \\midrule",
                            style$line.top)

    if(nchar(style$line.top) > 1){
        style$line.top = paste0(style$line.top, "\n")
    }

    # bottom line
    if(style$tablefoot){
        style$line.bottom = switch(style$line.bottom,
                                   "simple" = "\\midrule",
                                   "double" = "\\midrule \\midrule",
                                   style$line.bottom)
    } else {
        style$line.bottom = switch(style$line.bottom,
                                   "simple" = "\\bottomrule",
                                   "double" = "\\midrule \\midrule & \\tabularnewline",
                                   style$line.bottom)
    }



    if(nchar(style$line.bottom) > 1){
        style$line.bottom = paste0(style$line.bottom, "\n")
    }

    #
    # prompting the infos gathered
    #

    # Starting the table
    myTitle = title
    if(!is.null(label)){
        myTitle = paste0("\\label{", label, "} ", myTitle)
    }

    caption = ""
    info_center = "\\centering\n"
    if(float){
        if(nchar(placement) > 0){
            placement = paste0("[", placement, "]")
        }

        table_begin = paste0("\\begin{table}", placement, "\n")

        caption = paste0("\\caption{",  myTitle, "}\n")
        if(nchar(style$caption.after) > 0){
            caption = paste0(caption, style$caption.after, "\n")
        }

        if(!caption.number){
            table_begin = paste0("\\renewcommand{\\thetable}{}\n", table_begin)
        }

        table_end = "\\end{table}"
    } else {
        table_begin = "\\begingroup\n"
        table_end = "\\par\\endgroup\n"
    }

    info_width = ""
    if(!is.null(style$rules_width) && !all(is.na(style$rules_width))){
        w = style$rules_width
        info_width = c()
        if(!is.na(w[1])){
            info_width = paste0("\\setlength\\heavyrulewidth{", w[1], "}\n")
        }

        if(!is.na(w[2])){
            info_width = c(info_width, paste0("\\setlength\\lightrulewidth{", w[2], "}\n"))
        }
    }


    space = if(style$no_border) "@{}" else ""
    if(tabular == "normal"){
        tabular_begin = paste0("\\begin{tabular}{", space, "l",
                               paste0(rep("c", n_models), collapse = ""), space,
                               "}\n", style$line.top)

        tabular_end = "\\end{tabular}\n"
    } else if(tabular == "*"){
        tabular_begin = paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}",
                               space, "l", paste0(rep("c", n_models), collapse = ""),
                               space, "}\n", style$line.top)

        tabular_end = "\\end{tabular*}\n"
    } else if(tabular == "X"){

        all_cols = .dsb("l *.[n_models]{>{\\centering\\arraybackslash}X}")

        tabular_begin = paste0("\\begin{tabularx}{\\textwidth}{",
                               space, all_cols, space, "}\n", style$line.top)

        tabular_end = "\\end{tabularx}\n"
    }

    # 1st lines => dep vars
    depvar_list = dict_apply(c(depvar_list, recursive = TRUE), dict)
    depvar_list = escape_latex(depvar_list)

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

    if(style$depvar.style != ""){
        if(style$depvar.style == "*"){
            names_multi = paste0("\\textit{", names_multi, "}")
        } else if(style$depvar.style == "**"){
            names_multi = paste0("\\textbf{", names_multi, "}")
        } else {
            names_multi = paste0("\\textbf{\\textit{", names_multi, "}}")
        }
    }

    # now the proper format
    first_line = style$depvar.title
    if(length(nb_multi) == 1){
        first_line = gsub("(s)", "", first_line, fixed = TRUE)
    } else {
        first_line = gsub("(s)", "s", first_line, fixed = TRUE)
    }

    for(i in 1:length(nb_multi)){
        if(nb_multi[i] == 1){
            # no multi column
            first_line = paste0(first_line, " & ", names_multi[i])
        } else {
            first_line = paste0(first_line, " & \\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
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
    model_line = paste0(style$model.title, " & ", paste0(model_format, collapse = " & "), "\\\\\n")

    # a simple line with only "variables" written in the first cell
    if(nchar(style$var.title) == 0){
        coef_title = ""
    } else if(style$var.title == "\\midrule"){
        coef_title = "\\midrule "
    } else {
        coef_title = paste0(style$var.title, "\\\\\n")
    }

    # Coefficients, the tricky part
    coef_lines = list()

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

    # The names are set in results2formattedList
    coef_names = escape_latex(all_vars)
    names(coef_names) = all_vars

    if(se.below){
        coef_mat = c()
        for(v in all_vars){

            myCoef = mySd = myLine = c()
            for(m in 1:n_models){
                myCoef = c(myCoef, coef_below[[m]][v])
                mySd = c(mySd, sd_below[[m]][v])
            }

            myCoef[is.na(myCoef)] = "  "
            mySd[is.na(mySd)] = "  "

            coef_mat = rbind(coef_mat, c(coef_names[v], myCoef), c(" ", mySd))
        }

    } else {

        coef_mat = all_vars
        for(m in 1:n_models){
            coef_mat = cbind(coef_mat, coef_list[[m]][all_vars])
        }
        coef_mat[is.na(coef_mat)] = "  "

        coef_mat[, 1] = coef_names
    }


    coef_mat = style_apply(coef.style, coef_mat, all_vars)
    info_mat = highlight_apply(highlight, coef_mat, all_vars)
    coef_mat = info_mat$coef_mat
    info_tikz = info_mat$preamble

    coef_lines = paste0(apply(coef_mat, 1, paste, collapse = " & "), "\\\\ \n")
    coef_lines = paste0(coef_lines, collapse = "")

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
            fixef_title = paste0(style$fixef.title, "\\\\\n")
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
                nb_FE_lines = paste0(paste0(apply(all_nb_FEs, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")
            }

        }

        all_fe = cbind(fe_names, all_fe)
        FE_lines = paste0(paste0(apply(all_fe, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")

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
            slope_intro = paste0(style$slopes.title, "\\\\\n")
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
        slope_lines = paste0(paste0(apply(all_slopes, 1, paste0, collapse = " & "), collapse="\\\\\n"), "\\\\\n")

    } else {
        slope_intro = NULL
        slope_lines = NULL
    }

    # Headers
    if(isHeaders){
        n_head = length(headers)
        h_names = names(headers)


        headers_top = c()
        headers_mid = c()
        headers_bottom = c()

        for(i in 1:n_head){
            h_i = h_names[i]

            add_rule = FALSE
            if(grepl(":_:", h_i)){
                add_rule = TRUE
                h_i = gsub(":_:", "", h_i, fixed = TRUE)
            }

            h_placement = substr(h_i, 1, 1)
            if(h_placement %in% c("^", "-", "_")){
                h_i = substr(h_i, 2, nchar(h_i))
            } else {
                h_placement = "-"
            }

            h_value = paste0(h_i, " & ", tex_multicol(headers[[i]], add_rule = add_rule), "\n")

            if(h_placement == "^") {
                headers_top[length(headers_top) + 1] = h_value
            } else if(h_placement == "-"){
                headers_mid[length(headers_mid) + 1] = h_value
            } else {
                headers_bottom[length(headers_bottom) + 1] = h_value
            }

        }

        headers_top = paste(headers_top, collapse = "")
        headers_mid = paste(headers_mid, collapse = "")
        headers_bottom = paste(headers_bottom, collapse = "")

    } else {
        headers_top = ""
        headers_mid = ""
        headers_bottom = ""
    }

    # Convergence information
    info_convergence = ""
    if(convergence){
        info_convergence = paste0("Convergence &", paste(convergence_list, collapse = " & "), "\\\\\n")
    }

    # information on family
    if(family){
        info_family = paste0(" &  ", paste(family_list, collapse = " & "), "\\\\\n")
    } else {
        info_family = ""
    }


    # The standard errors => if tablefoot = TRUE

    # We go through this in order to get the multi se if needed
    isUniqueSD = length(unique(unlist(se_type_list))) == 1
    nb_col = length(obs_list) + 1
    foot_intro = paste0("\\multicolumn{", nb_col, "}{l}{\\emph{")

    if(is.null(se.row)){
        # we show the SE row if:
        # - NO FOOTER & MULTIPLE SEs

        se.row = isFALSE(style$tablefoot) && !isUniqueSD
    }

    # Automatic footer information
    info_SE_footer = ""
    if(style$tablefoot){

        bottom_line = style$line.bottom
        style$line.bottom = ""

        if(identical(style$tablefoot.value, "default")){

            if(isUniqueSD){
                my_se = unique(unlist(se_type_list)) # it comes from summary
                # every model has the same type of SE

                # Now we modify the names of the clusters if needed
                my_se = format_se_type_latex(my_se, dict)

                if(coefstat == "se"){
                    coefstat_sentence = " standard-errors in parentheses"
                } else if(coefstat == "tstat"){
                    coefstat_sentence = " co-variance matrix, t-stats in parentheses"
                } else {
                    coefstat_sentence = paste0(" co-variance matrix, ", round(ci*100), "\\% confidence intervals in brackets")
                }
                info_SE_footer = paste0(bottom_line, foot_intro, my_se,
                                        coefstat_sentence, "}}\\\\\n")

                if(add_signif){
                    info_SE_footer = paste0(info_SE_footer, foot_intro, "Signif. Codes: ",
                                            paste(names(signif.code), signif.code, sep=": ", collapse = ", "), "}}\\\\\n")
                }

            } else {
                # Multiple types of SEs: we're more general
                all_se_type = sapply(se_type_list, format_se_type_latex, dict = dict, inline = TRUE)

                if(coefstat == "se"){
                    coefstat_sentence = "Standard-Errors"
                } else {
                    coefstat_sentence = "Co-variance"
                }

                info_muli_se = paste0(coefstat_sentence, " & ", paste(all_se_type, collapse = " & "), "\\\\\n")

                if(add_signif){
                    info_SE_footer = paste0(bottom_line, foot_intro, "Signif. Codes: ",
                                            paste(names(signif.code), signif.code,
                                                  sep = ": ", collapse = ", "),
                                            "}}\\\\\n")
                } else {
                    myAmpLine = paste0(paste0(rep(" ", length(depvar_list) + 1), collapse = " & "), "\\tabularnewline\n")
                    info_SE_footer = paste0(bottom_line, myAmpLine, "\\\\\n")
                }

            }

        } else {

            value = style$tablefoot.value
            if(isUniqueSD){
                my_se = unique(unlist(se_type_list))
                my_se = format_se_type_latex(my_se, dict)

                value = gsub("__se_type__", my_se, value)
            } else {
                # not super elegant, but can't to much more // else: enumerate the SEs?
                value = gsub("__se_type__", "", value)
            }

            info_SE_footer = paste0(bottom_line, paste(foot_intro, value, "}}\\\\\n", collapse = ""))

        }
    }

    info_muli_se = ""
    if(se.row){
        all_se_type = sapply(se_type_list, format_se_type_latex, dict = dict, inline = TRUE)

        if(coefstat == "se"){
            coefstat_sentence = "Standard-Errors"
        } else {
            coefstat_sentence = "Co-variance"
        }

        info_muli_se = paste0(coefstat_sentence, " & ", tex_multicol(all_se_type), "\n")
    }

    #
    # Fit statistics
    #

    if(!"stats" %in% drop.section){
        if(nchar(style$stats.title) == 0){
            stat_title = ""
        } else if(style$stats.title == "\\midrule"){
            stat_title = "\\midrule "
        } else {
            stat_title = paste0(style$stats.title, "\\\\\n")
        }

        stat_lines = paste0(nb_FE_lines, info_convergence, info_muli_se)

        if(!all(sapply(fitstat_list, function(x) all(is.na(x))))){

            fit_names = attr(fitstat_list, "format_names")
            nb = length(fit_names)
            for(fit_id in 1:nb){
                fit = sapply(fitstat_list, function(x) x[[fit_id]])
                if(all(is.na(fit))) next
                fit[is.na(fit)] = ""
                stat_lines = paste0(stat_lines, fit_names[fit_id], " & ", paste0(fit, collapse = " & "), "\\\\\n")
            }
        }
    } else {
        stat_title = stat_lines = ""
    }


    # Notes
    info_notes = ""
    if(length(notes) > 0){

        notes_intro = if(tpt) style$notes.tpt.intro else style$notes.intro

        if(grepl("^@", notes[1])){
            notes_intro = paste0(gsub("^@", "", notes[1]), " ")
            notes = notes[-1]
        }

        if(length(notes) > 0){
            notes = dict_apply(notes, dict)
            notes = escape_latex(notes, makecell = FALSE)

            if(tpt){
                if(nchar(trimws(notes_intro)) > 0){
                    notes_intro = paste0(trimws(notes_intro), "\n")
                }

                notes = paste0("\\item ", notes, collapse = "\n")

                # The note intro is placed right after the } so that you can pass options
                # like [flushleft]
                info_notes = paste0("\n\\begin{tablenotes}", notes_intro,
                                    notes,
                                    "\n\\end{tablenotes}\n")
            } else {
                info_notes = paste0("\n", notes_intro, paste0(notes, collapse = "\\\\\n"), "\n")
            }
        }
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
            # first character: section: ^|-|_
            # second character: the location: ^|_

            # sec: section

            sec = substr(gi_name, 1, 1)
            gi_name = substr(gi_name, 2, nchar(gi_name))
            if(grepl("^(\\^|_)", gi_name)){
                row = substr(gi_name, 1, 1)
                gi_name = substr(gi_name, 2, nchar(gi_name))
                gi_top = row == "^"
            }

            gi_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat")
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

    for(i in seq_along(extralines)){
        el = extralines[[i]]

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

        el_name = names(extralines)[i]

        el_full = ""
        el_where = "coef"

        el_top = FALSE
        if(grepl("^(\\^|_|-)", el_name)){
            # first character: section: ^|-|_
            # second character: the location: ^|_

            # sec: section

            sec = substr(el_name, 1, 1)
            el_name = substr(el_name, 2, nchar(el_name))
            if(grepl("^(\\^|_)", el_name)){
                row = substr(el_name, 1, 1)
                el_name = substr(el_name, 2, nchar(el_name))
                el_top = row == "^"
            }

            el_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat")
        }

        if(grepl(":tex:", el_name, fixed = TRUE)){
            # No escaping
            el_name = gsub(":tex:", "", el_name, fixed = TRUE)
        } else {
            # escaping
            el_name = escape_latex(el_name)
            el_format = escape_latex(el_format)
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
            fixef_title = paste0(style$fixef.title, "\\\\\n")
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
        tag_tabular_before = "%start:tab\n"
        tag_tabular_end = "%end:tab\n"
    } else {
        tag_tabular_before = tag_tabular_end = ""
    }

    #
    # Other info:
    #

    tpt_begin = tpt_end = tpt_caption = ""
    if(tpt){
        tpt_begin = "\\begin{threeparttable}[b]\n"
        tpt_end = "\\end{threeparttable}\n"
        tpt_caption = caption
        caption = ""
    }

    info_stretch = ""
    if(!is.null(arraystretch)){
        info_stretch = paste0("\\renewcommand*{\\arraystretch}{", arraystretch, "}")
    }

    adj_box_begin = adj_box_end = ""
    if(!is.null(adjustbox)){
        adj_box_begin = paste0("\\begin{adjustbox}{", adjustbox, "}\n")
        adj_box_end = "\\end{adjustbox}\n"
    }

    info_font = ""
    if(!is.null(fontsize)){
        info_font = paste0("\\", fontsize, "\n")
    }

    # meta information: has been computed in results2formattedList

    res = c(meta, table_begin, info_tikz, caption, info_center, info_width, info_font, adj_box_begin,
            tpt_begin, tpt_caption, info_stretch,
            tag_tabular_before, tabular_begin,
            headers_top, first_line, headers_mid, model_line, info_family, headers_bottom,
            coef_stack, stat_stack, info_SE_footer, style$line.bottom,
            tabular_end, tag_tabular_end,
            info_notes, tpt_end, adj_box_end, table_end)

    res = res[nchar(res) > 0]

    attr(res, "n_models") = n_models

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
    se.below = info$se.below
    signif.code = info$signif.code
    add_signif = info$add_signif
    fixef_sizes = info$fixef_sizes
    depvar = info$depvar
    useSummary = info$useSummary
    model_names = info$model_names
    coefstat = info$coefstat
    dict = info$dict
    group = info$group
    extralines = info$extralines
    style = info$style
    fun_format = info$fun_format
    drop.section = info$drop.section
    coef.just = info$coef.just
    se.row = info$se.row

    # naming differences
    headers = info$headers
    isHeaders = info$isHeaders

    # The coefficients

    depvar_list = dict_apply(unlist(depvar_list), dict)

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

    se.below = info$se.below
    if(se.below){
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
        coef_mat = all_vars
        for(m in 1:n_models) coef_mat = cbind(coef_mat, coef_list[[m]][all_vars])
        coef_mat[is.na(coef_mat)] = "  "
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
            # first character: section: ^|-|_
            # second character: the location: ^|_

            # sec: section

            sec = substr(gi_name, 1, 1)
            gi_name = substr(gi_name, 2, nchar(gi_name))
            if(grepl("^(\\^|_)", gi_name)){
                row = substr(gi_name, 1, 1)
                gi_name = substr(gi_name, 2, nchar(gi_name))
                gi_top = row == "^"
            }

            gi_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat")
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

    for(i in seq_along(extralines)){
        el = extralines[[i]]
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

        el_name = names(extralines)[i]
        el_top = FALSE

        if(grepl("^(\\^|_|-)", el_name)){
            # first character: section: ^|-|_
            # second character: the location: ^|_

            # sec: section

            sec = substr(el_name, 1, 1)
            el_name = substr(el_name, 2, nchar(el_name))
            if(grepl("^(\\^|_)", el_name)){
                row = substr(el_name, 1, 1)
                el_name = substr(el_name, 2, nchar(el_name))
                el_top = row == "^"
            }

            el_where = switch(sec, "^" = "coef", "-" = "fixef", "_" = "stat")
        }

        # we clean possible tex markup
        el_name = gsub(":tex:", "", el_name, fixed = TRUE)

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

    # the headers
    if(isHeaders){
        # we need to provide unique names... sigh...

        n_head = length(headers)
        h_names = names(headers)

        headers_top = c()
        headers_bottom = c()

        for(i in 1:n_head){
            h_i = h_names[i]

            h_placement = substr(h_i, 1, 1)
            if(h_placement %in% c("^", "-", "_")){
                h_i = substr(h_i, 2, nchar(h_i))
            } else {
                h_placement = "-"
            }

            h_value = c(sfill(h_i, i), headers[[i]])

            if(h_placement %in% c("^", "-")){
                headers_top = rbind(headers_top, h_value)
            } else {
                headers_bottom = rbind(headers_bottom, h_value)
            }

        }

        preamble = rbind(headers_top, preamble, headers_bottom)

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
        if(style$headers.sep){
            preamble = rbind(preamble, rep(" ", length(longueur)))
        }

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

    if(!isFALSE(se.row)){
        # default is to always show the SE row
        if(coefstat == "se"){
            coefstat_sentence = "S.E. type"
        } else {
            coefstat_sentence = "VCOV type"
        }

        se_type_format = c()
        for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]], by = TRUE)

        main_type = ""
        if(all(grepl("Clustered", unlist(se_type_list), fixed = TRUE))){
            main_type = ": Clustered"
            coefstat_sentence = gsub(" type", "", coefstat_sentence)
        }

        res = rbind(res, c(paste0(coefstat_sentence, main_type), c(se_type_format, recursive = TRUE)))
    }

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

    class(res) = c("etable_df", "data.frame")

    if(add_signif){
        attr(res, "signif") = signif.code
    }

    return(res)
}

####
#### User-level ####
####



#' @rdname etable
setFixest_etable = function(digits = 4, digits.stats = 5, fitstat,
                            coefstat = c("se", "tstat", "confint"),
                            ci = 0.95, se.below = TRUE, keep, drop, order, dict,
                            float,
                            fixef_sizes = FALSE, fixef_sizes.simplify = TRUE,
                            family, powerBelow = -5,
                            interaction.order = NULL, depvar, style.tex = NULL,
                            style.df = NULL, notes = NULL, group = NULL, extralines = NULL,
                            fixef.group = NULL, placement = "htbp", drop.section = NULL,
                            view = FALSE, markdown = NULL, view.cache = FALSE,
                            page.width = "fit",
                            postprocess.tex = NULL, postprocess.df = NULL,
                            fit_format = "__var__", meta.time = NULL,
                            meta.author = NULL, meta.sys = NULL,
                            meta.call = NULL, meta.comment = NULL,
                            reset = FALSE, save = FALSE){


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

    check_arg("logical scalar", se.below, fixef_sizes, fixef_sizes.simplify,
              float, family, depvar, reset, view)

    check_arg(keep, drop, order, "character vector no na NULL",
              .message = "The arg. '__ARG__' must be a vector of regular expressions (see help(regex)).")

    check_arg(interaction.order, "NULL character scalar")

    check_arg(notes, "character vector no na")

    check_arg(powerBelow, "integer scalar LE{-1}")

    check_arg(dict, "NULL logical scalar | named character vector no na")

    check_arg_plus(group, extralines, "NULL{list()} named list l0")

    check_arg_plus(fixef.group, "NULL{list()} logical scalar | named list l0")

    check_arg(placement, "character scalar")
    if(!missing(placement)){
        if(nchar(placement) == 0) stop("Argument 'placement' cannot be the empty string.")
        p_split = strsplit(placement, "")[[1]]
        check_value(p_split, "strict multi charin(h, t, b, p, H, !)", .message = "Argument 'placement' must be a character string containing only the following characters: 'h', 't', 'b', 'p', 'H', and '!'.")
    }

    check_arg_plus(drop.section, "NULL multi match(fixef, slopes, stats)")

    check_arg(style.tex, "NULL class(fixest_style_tex)")

    check_arg(postprocess.tex, postprocess.df, "NULL function arg(1,)")

    check_arg(fit_format, "character scalar")
    if(!grepl("__var__", fit_format, fixed = TRUE)){
        stop("The argument 'fit_format' should include the special name '__var__' that will be replaced by the variable name. So far it does not contain it.")
    }

    # meta
    check_arg_plus(meta.time, "NULL match(date, time) | logical scalar")
    check_arg(meta.call, meta.sys, "NULL logical scalar")
    check_arg(meta.author, "NULL logical scalar | character vector no na")
    check_arg(meta.comment, "NULL character vector no na")

    check_arg(view, view.cache, "logical scalar")
    check_arg(markdown, "NULL scalar(logical, character)")

    page.width = check_set_page_width(page.width)

    #
    # Setting the defaults
    #

    # Getting the existing defaults
    opts = getOption("fixest_etable")

    if(is.null(opts)){
        # We first look at the "root" default
        root_default = renvir_get("fixest_etable")
        if(is.null(root_default)){
            opts = list()
        } else {
            opts = root_default
        }

    } else if(!is.list(opts)){
        warning("Wrong formatting of option 'fixest_etable', all options are reset.")
        opts = list()

    } else if(reset){
        opts = list()
    }

    # Style setting
    if(length(style.tex) > 0){
        # We ensure we always have ALL components provided
        if(length(opts$style.tex) == 0){
            basic_style = fixest::style.tex(main = "base")
        } else {
            basic_style = opts$style.tex
        }

        basic_style[names(style.tex)] = style.tex
        style.tex = basic_style

    }

    check_arg(style.df, "NULL class(fixest_style_df)")
    if(length(style.df) > 0){
        # We ensure we always have ALL components provided
        if(length(opts$style.df) == 0){
            basic_style = fixest::style.df(default = TRUE)
        } else {
            basic_style = opts$style.df
        }

        basic_style[names(style.df)] = style.df
        style.df = basic_style
    }

    # Saving the default values
    mc = match.call()
    args_default = setdiff(names(mc)[-1], c("reset", "save"))

    # NOTA: we don't allow delayed evaluation => all arguments must have hard values
    for(v in args_default){
        opts[[v]] = eval(as.name(v))
    }

    options(fixest_etable = opts)

    # Saving at the project level if needed
    check_arg_plus(save, "logical scalar | match(reset)")
    if(isTRUE(save)){
        renvir_update("fixest_etable", opts)

    } else if(identical(save, "reset")){
        renvir_update("fixest_etable", NULL)
    }

}

#' @rdname etable
getFixest_etable = function(){
    opts = getOption("fixest_etable")
    if(!is.list(opts)){
        warning("Wrong formatting of option 'fixest_etable', all options are reset.")
        opts = list()
        options(fixest_etable = opts)
    }
    opts
}


#' Style definitions for Latex tables
#'
#' This function describes the style of Latex tables to be exported with the function \code{\link[fixest]{etable}}.
#'
#' @inheritParams etable
#'
#' @param main Either "base", "aer" or "qje". Defines the basic style to start from. The styles "aer" and "qje" are almost identical and only differ on the top/bottom lines.
#' @param depvar.title A character scalar. The title of the line of the dependent variables (defaults to \code{"Dependent variable(s):"} if \code{main = "base"} (the 's' appears only if just one variable) and to \code{""} if \code{main = "aer"}).
#' @param model.title A character scalar. The title of the line of the models (defaults to \code{"Model:"} if \code{main = "base"} and to \code{""} if \code{main = "aer"}).
#' @param model.format A character scalar. The value to appear on top of each column. It defaults to \code{"(1)"}. Note that 1, i, I, a and A are special characters: if found, their values will be automatically incremented across columns.
#' @param line.top A character scalar equal to \code{"simple"}, \code{"double"}, or anything else. The line at the top of the table (defaults to \code{"double"} if \code{main = "base"} and to \code{"simple"} if \code{main = "aer"}). \code{"simple"} is equivalent to \code{"\\toprule"}, and \code{"double"} to \code{"\\tabularnewline \\midrule \\midrule"}.
#' @param line.bottom A character scalar equal to \code{"simple"}, \code{"double"}, or anything else. The line at the bottom of the table (defaults to \code{"double"} if \code{main = "base"} and to \code{"simple"} if \code{main = "aer"}). \code{"simple"} is equivalent to \code{"\\bottomrule"}, and \code{"double"} to \code{"\\midrule \\midrule & \\tabularnewline"}.
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
#' @param notes.intro A character scalar. Some tex code appearing just before the notes, defaults to \code{"\\par \\raggedright \n"}.
#' @param tablefoot A logical scalar. Whether or not to display a footer within the table. Defaults to \code{TRUE} if \code{main = "aer"}) and \code{FALSE} if \code{main = "aer"}).
#' @param tablefoot.value A character scalar. The notes to be displayed in the footer. Defaults to \code{"default"} if \code{main = "base"}, which leads to custom footers informing on the type of standard-error and significance codes, depending on the estimations.
#' @param yesNo A character vector of length 1 or 2. Defaults to \code{"Yes"} if \code{main = "base"} and to \code{"$\\checkmark$"} if \code{main = "aer"} (from package \code{amssymb}). This is the message displayed when a given fixed-effect is (or is not) included in a regression. If \code{yesNo} is of length 1, then the second element is the empty string.
#' @param interaction.combine Character scalar, defaults to \code{" $\\times$ "}. When the estimation contains interactions, then the variables names (after aliasing) are combined with this argument. For example: if \code{dict = c(x1="Wind", x2="Rain")} and you have the following interaction \code{x1:x2}, then it will be renamed (by default) \code{Wind $\\times$ Rain} -- using \code{interaction.combine = "*"} would lead to \code{Wind*Rain}.
#' @param i.equal Character scalar, defaults to \code{" $=$ "}. Only affects factor variables created with the function \code{\link[fixest]{i}}, tells how the variable should be linked to its value. For example if you have the Species factor from the iris data set, by default the display of the variable is \code{Species $=$ Setosa}, etc. If \code{i.equal = ": "} the display becomes \code{Species: Setosa}.
#' @param notes.tpt.intro Character scalar. Only used if \code{tpt = TRUE}, it is some tex code that is passed before any \code{threeparttable} item (can be used for, typically, the font size). Default is the empty string.
#' @param depvar.style Character scalar equal to either \code{" "} (default), \code{"*"} (italic), \code{"**"} (bold), \code{"***"} (italic-bold). How the name of the dependent variable should be displayed.
#' @param no_border Logical, default is \code{FALSE}. Whether to remove any side border to the table (typically adds \code{@\{\}} to the sides of the tabular).
#' @param caption.after Character scalar. Tex code that will be placed right after the caption. Defaults to \code{""} for \code{main = "base"} and \code{"\\medskip"} for \code{main = "aer"}.
#' @param rules_width Character vector of length 1 or 2. This vector gives the width of the \code{booktabs} rules: the first element the heavy-width, the second element the light-width. NA values mean no modification. If of length 1, only the heavy rules are modified. The width are in Latex units (ex: \code{"0.1 em"}, etc).
#' @param signif.code Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)}. To suppress the significance codes, use \code{signif.code=NA} or \code{signif.code=NULL}. Can also be equal to \code{"letters"}, then the default becomes \code{c("a"=0.01, "b"=0.05, "c"=0.10)}.
#'
#' @details
#' The \code{\\checkmark} command, used in the "aer" style (in argument \code{yesNo}), is in the \code{amssymb} package.
#'
#' The commands \code{\\toprule}, \code{\\midrule} and \code{\\bottomrule} are in the \code{booktabs} package. You can set the width of the top/bottom rules with \\setlength\\heavyrulewidth\{wd\}, and of the midrule with \code{\\setlength\\lightrulewidth\{wd\}}.
#'
#' Note that all titles (\code{depvar.title}, \code{depvar.title}, etc) are not escaped, so they must be valid Latex expressions.
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
style.tex = function(main = "base", depvar.title, model.title, model.format, line.top,
                     line.bottom, var.title, fixef.title, fixef.prefix, fixef.suffix,
                     fixef.where, slopes.title, slopes.format, fixef_sizes.prefix,
                     fixef_sizes.suffix, stats.title, notes.intro,
                     notes.tpt.intro, tablefoot, tablefoot.value,
                     yesNo, tabular = "normal", depvar.style, no_border,
                     caption.after, rules_width, signif.code,
                     tpt, arraystretch, adjustbox = NULL, fontsize,
                     interaction.combine = " $\\times$ ", i.equal = " $=$ "){

    # To implement later:
    # fixef_sizes.where = "obs"
    # se.par = "parentheses"
    # check_arg_plus(se.par, "match(parentheses, brackets, none)")
    # align
    # check_arg(align, "character vector no na len(,2)")

    # Checking
    check_arg_plus(main, "match(base, aer, qje)")

    check_arg_plus(line.top, line.bottom, "match(simple, double) | character scalar")

    check_arg("character scalar", depvar.title, model.title, var.title)
    check_arg("character scalar", fixef.title, fixef.prefix, fixef.suffix, slopes.title, slopes.format)
    check_arg("character scalar", fixef_sizes.prefix, fixef_sizes.suffix, stats.title)
    check_arg("character scalar", notes.intro, interaction.combine, i.equal)
    check_arg("character scalar", notes.tpt.intro, depvar.style, caption.after)
    check_arg("character vector len(1,2)", rules_width)

    if(!missing(depvar.style)){
        depvar.style = trimws(depvar.style)
        if(!depvar.style %in% c("", "*", "**", "***")){
            stop("Argument 'depvar.style' must be one of ' ', '*', '**', '***', which means regular, italic, bold and italic-bold.")
        }
    }

    check_arg(tablefoot.value, "character vector no na")
    check_arg(tablefoot, tpt, no_border, "logical scalar")
    check_arg_plus(fixef.where, "match(var, stats)")
    check_arg_plus(tabular, "match(normal, *, X)")

    check_arg(yesNo, "character vector len(,2) no na")
    if(!missing(yesNo) && length(yesNo) == 1){
        yesNo = c(yesNo, "")
    }

    check_arg(arraystretch, "numeric scalar GT{0}")
    adjustbox = check_set_adjustbox(adjustbox)
    check_arg_plus(fontsize, "NULL match(tiny, scriptsize, footnotesize, small, normalsize, large, Large)")

    check_arg_plus(signif.code, "NULL NA | match(letters) | named numeric vector no na GE{0} LE{1}")

    mc = match.call()

    if("main" %in% names(mc)){
        if(main == "base"){
            res = list(depvar.title = "Dependent Variable(s):", model.title = "Model:",
                       model.format = "(1)",
                       line.top = "\\tabularnewline \\midrule \\midrule",
                       line.bottom = "double",
                       var.title = "\\midrule\n\\emph{Variables}",
                       fixef.title = "\\midrule\n\\emph{Fixed-effects}", fixef.prefix = "",
                       fixef.suffix = "", fixef.where = "var",
                       slopes.title = "\\midrule\n\\emph{Varying Slopes}",
                       slopes.format = "__var__ (__slope__)",
                       fixef_sizes.prefix = "\\# ", fixef_sizes.suffix = "",
                       stats.title = "\\midrule\n\\emph{Fit statistics}",
                       notes.intro = "\\par \\raggedright \n",
                       notes.tpt.intro = "",
                       tablefoot = TRUE,
                       caption.after = "",
                       tablefoot.value = "default", yesNo = c("Yes", ""),
                       depvar.style = "", no_border = FALSE)

            if(!missing(tablefoot) && isFALSE(tablefoot)){
                res$tablefoot = FALSE
                res$line.bottom = "\\midrule \\midrule & \\tabularnewline"
            }

        } else {
            res = list(depvar.title = "", model.title = "", model.format = "(1)",
                       line.top = "\\toprule", line.bottom = "\\bottomrule",
                       var.title = "\\midrule",
                       fixef.title = " ", fixef.prefix = "", fixef.suffix = " fixed effects",
                       fixef.where = "stats",
                       slopes.title = "", slopes.format = "__var__ $\\times $ __slope__",
                       fixef_sizes.prefix = "\\# ", fixef_sizes.suffix = "",
                       stats.title = " ", notes.intro = "\\par \\raggedright \n",
                       notes.tpt.intro = "", caption.after = "\\bigskip",
                       tablefoot = FALSE, tablefoot.value = "",
                       yesNo = c("$\\checkmark$", ""), depvar.style = "",
                       no_border = FALSE)

            if(main == "aer"){
                # just set

            } else if(main == "qje"){
                res$line.top = "\\tabularnewline\\midrule\\midrule"
                res$line.bottom = "\\midrule \\midrule & \\tabularnewline"
            }
        }

        res$tabular = tabular
        res$interaction.combine = interaction.combine
        res$i.equal = i.equal
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
#'  @inheritParams etable
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
#' @param headers.sep Logical, default is \code{TRUE}. Whether to add a line of separation between the headers and the coefficients.
#' @param interaction.combine Character scalar, defaults to \code{" x "}. When the estimation contains interactions, then the variables names (after aliasing) are combined with this argument. For example: if \code{dict = c(x1="Wind", x2="Rain")} and you have the following interaction \code{x1:x2}, then it will be renamed (by default) \code{Wind x Rain} -- using \code{interaction.combine = "*"} would lead to \code{Wind*Rain}.
#' @param i.equal Character scalar, defaults to \code{" = "}. Only affects factor variables created with the function \code{\link[fixest]{i}}, tells how the variable should be linked to its value. For example if you have the Species factor from the iris data set, by default the display of the variable is \code{Species = Setosa}, etc. If \code{i.equal = ": "} the display becomes \code{Species: Setosa}.
#' @param default Logical, default is \code{FALSE}. If \code{TRUE}, all the values not provided by the user are set to their default.
#' @param signif.code Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)}. To suppress the significance codes, use \code{signif.code=NA} or \code{signif.code=NULL}. Can also be equal to \code{"letters"}, then the default becomes \code{c("a"=0.01, "b"=0.05, "c"=0.10)}.
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
#' etable(est, style.df = style.df(fixef.title = "", fixef.suffix = " FE",
#'                                  stats.line = " ", yesNo = "yes"))
#'
#'
style.df = function(depvar.title = "Dependent Var.:", fixef.title = "Fixed-Effects:",
                    fixef.line = "-", fixef.prefix = "", fixef.suffix = "",
                    slopes.title = "Varying Slopes:", slopes.line = "-",
                    slopes.format = "__var__ (__slope__)", stats.title = "_",
                    stats.line = "_", yesNo = c("Yes", "No"), headers.sep = TRUE,
                    signif.code = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10),
                    interaction.combine = " x ", i.equal = " = ", default = FALSE){

    # Checking

    check_arg("character scalar", depvar.title, fixef.title, fixef.line, fixef.prefix, fixef.suffix)
    check_arg("character scalar", slopes.title, slopes.line, slopes.format, stats.title, stats.line)
    check_arg("character scalar", interaction.combine, i.equal)

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

    check_arg(headers.sep, default, "logical scalar")

    check_arg_plus(signif.code, "NULL NA | match(letters) | named numeric vector no na GE{0} LE{1}")

    res = list()

    mc = match.call()
    args2set = if(default) formalArgs(style.df) else names(mc)[-1]

    args2set = setdiff(args2set, "default")
    for(var in args2set){
        res[[var]] = eval(as.name(var))
    }

    class(res) = "fixest_style_df"

    return(res)
}


#' Register \code{extralines} macros to be used in \code{etable}
#'
#' This function is used to create \code{extralines} (which is an argument of \code{\link[fixest]{etable}}) macros that can be easily summoned in \code{\link[fixest]{etable}}.
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
#' extralines_register("sdy", my_fun, "SD(y)")
#'
#' # An estimation
#' data(iris)
#' est = feols(Petal.Length ~ Sepal.Length | Species, iris)
#'
#' # Now we can easily create a row with the mean of y.
#' # We just "summon" it in a one-sided formula
#' etable(est, extralines = ~ sdy)
#'
#' # We can change the alias on the fly:
#' etable(est, extralines = list("_Standard deviation of the dep. var." = ~ sdy))
#'
#'
#'
#'
extralines_register = function(type, fun, alias){
    check_arg(type, "character scalar mbt")
    check_arg(fun, "function mbt")
    check_arg(alias, "character scalar mbt")

    # We check the type is not conflicting
    existing_types = fitstat(give_types = TRUE)$types

    opts = getOption("fixest_extralines")

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

    options(fixest_extralines = opts)

    invisible(NULL)
}


#' @rdname etable
print.etable_tex = function(x, ...){
    cat(x, sep = "\n")
}

#' @rdname etable
print.etable_df = function(x, ...){
    # Almost equivalent to print.data.frame

    width = 0.99
    MAX_WIDTH = getOption("width") * width

    ROWNAMES = row.names(x)
    labels_width = max(nchar(ROWNAMES))
    MAX_WIDTH = MAX_WIDTH - labels_width - 1

    # names => first row
    xx = rbind(names(x), x)
    col_width = apply(xx, 2, function(x) max(nchar(x)))
    ROWNAMES = sprintf("%- *s", labels_width, c(" ", ROWNAMES))

    # We print on the fly
    my_col = 1
    nc = ncol(x)

    while(my_col <= nc){

        # We jump one line for readability
        if(my_col > 1) cat("\n")

        # Finding the nber of columns to print out
        cs = cumsum(col_width + 1) - 1
        if(cs[1] > MAX_WIDTH){
            n_max = 1
        } else {
            n_max = max(which(cs <= MAX_WIDTH))
        }

        m = cbind(ROWNAMES, xx[, seq(my_col, length.out = n_max)])

        for(j in 2:ncol(m)){
            m[, j] = sprintf("% *s", col_width[j - 1], m[, j])
        }

        table_print = apply(m, 1, paste, collapse = " ")
        cat(table_print, sep = "\n")

        # Next
        my_col = my_col + n_max
        col_width = col_width[-(1:n_max)]

    }

    signif.code = attr(x, "signif")
    if(!is.null(signif.code)){
        s = signif.code
        s_fmt = paste0(insert_in_between(.dsb("q?names(s)"), s), collapse = " ")
        cat("---\nSignif. codes: 0 ", s_fmt, " ' ' 1\n", sep = "")
    }

}

####
#### Viewer ####
####

check_build_available = function(){

    opt = getOption("fixest_build_available")

    if(!isTRUE(opt)){
        outcome = suppressWarnings(system2("pdflatex", "-help", FALSE, FALSE))
        if(outcome == 127 && !requireNamespace("tinytex", quietly = TRUE)){
            warn_up("The functionality you want to use requires the package 'tinytex' which is not installed or a working pdflatex installation which wasn't found.")

            options(fixest_build_available = "pdflatex")
            return("pdflatex")
        }

        outcome = suppressWarnings(system2("magick", "-help", FALSE, FALSE))
        if(outcome == 127 && !requireNamespace("pdftools", quietly = TRUE)){
            warn_up("The functionality you want to use requires the package 'pdftools' which is not installed or a working imagemagick + ghostscript installation which wasn't found.")

            options(fixest_build_available = "magick")
            return("magick")
        }

        options(fixest_build_available = TRUE)
    }

    return(TRUE)
}

build_tex_png = function(x, view = FALSE, export = NULL, markdown = NULL,
                         cache = FALSE, page.width = "fit", up = 0){

    up = up + 1
    set_up(up)
    check_value(view, "logical scalar")
    check_value(export, "NULL character scalar")
    check_value(markdown, "NULL scalar(logical, character)")
    page.width = check_set_page_width(page.width)

    # numbers used to uniquely identify the images
    NMAX = 10

    cache = isTRUE(cache)

    dir_cache = NULL
    if(cache || isTRUE(markdown)){
        is_global_dir = !is.null(find_config_path()) && !is.null(find_project_path())
        if(is_global_dir){
            dir_cache = config_get("cache_dir")

            if(!DIR_EXISTS(dir_cache)){
                # reset
                dir_cache = file.path(normalizePath(tempdir(), "/"), "etable")
                dir.create(dir_cache)
                config_update("cache_dir", dir_cache)
            }
        }
    }

    if(isFALSE(markdown)){
        markdown = NULL
    } else if(isTRUE(markdown)){
        # note that it's very important to have the default saving path to be relative
        # to the working directory and not (as I initially wanted) in a temp dir.
        # Why? To make it work properly in webpages, it should be self contained
        # and should not use things out of the current working dir.
        #
        # (In regular markdown, the temp dir is a bit more elegant since
        #  it does not clutter the workspace)

        markdown = "./images/etable/"
    }

    # we need to clean all the tmp tags otherwise the caching does not work
    x_clean = paste0(c(page.width, x), collapse = "\n")
    if(any(grepl("definecolor", x_clean, fixed = TRUE))){
        my_colors = dsb("'(?<=definecolor\\{)[[:alpha:]]+(?=\\})'X, '|'c ? x_clean")
        x_clean = gsub(my_colors, "", x_clean)
    }

    if(any(grepl("\\mark", x, fixed = TRUE))){
        x_clean = gsub("mark[[:alpha:]]+", "", x_clean)
    }

    do_build = TRUE
    export_markdown = id = NULL
    if(!is.null(markdown)){
        markdown_path = check_set_path(markdown, "w, dir", create = TRUE, up = up)

        all_files = list.files(markdown_path, "\\.png$", full.names = TRUE)
        id_all = gsub("^.+_|\\.png$", "", all_files)
        id = substr(cpp_hash_string(x_clean), 1, NMAX)

        if(!id %in% id_all){
            time = gsub(" .+", "", Sys.time())
            name = .dsb("etable_tex_.[time]_.[id].png")
            export_markdown = file.path(markdown_path, name)
        } else {
            do_build = FALSE
            export_markdown = png_name = normalizePath(all_files[id_all == id][1], "/")
        }
    }

    if(!is.null(export)){
        export_path = check_set_path(export, "w", up = up)
    }

    dir = NULL
    if(do_build){

        dir = if(cache && !is.null(dir_cache)) dir_cache else normalizePath(tempdir(), "/")

        if(cache){
            if(is.null(id)){
                id = substr(cpp_hash_string(x_clean), 1, NMAX)
            }

            # If the file already exists, we don't recreate it
            all_files = list.files(dir, "^etable.+\\.png$", full.names = TRUE)
            id_all = gsub("^.+_|\\.png$", "", all_files)
            if(id %in% id_all){
                do_build = FALSE
                png_name = normalizePath(all_files[id_all == id][1], "/")
            } else {
                time = gsub(" .+", "", Sys.time())
                png_name = .dsb("etable_tex_.[time]_.[id].png")
            }
        } else {
            png_name = "etable.png"
        }
    }

    if(view || do_build){
        # we set the working directory properly
        # used to:
        #  - compile .tex to .pdf to .png
        #  - view the HTML document (needs to be on the WD)

        if(is.null(dir)) dir = normalizePath(tempdir(), "/")

        current_dir = getwd()
        on.exit(setwd(current_dir))
        setwd(dir)
    }

    if(do_build){

        #
        # Latex document
        #

        # packages increase build time, so we load them sparingly
        # p: package ; pn: package name ; x: tex vector ; y: tex packages
        add_pkg = function(p, x, y, pn = p, opt = "", fixed = TRUE){
            if(any(grepl(p, x, fixed = fixed))){
                c(y, .dsb("\\usepackage.[opt]{.[pn]}"))
            } else {
                y
            }
        }

        minipage_start = minipage_end = ""

        w = "35cm"
        if(!identical(as.vector(page.width), "fit", )){
            # page.width is guaranteed to be of length 2
            w = page.width[1]
            minipage_start = .dsb("\\begin{minipage}{.[page.width[1]]} ",
                                  "\\centering ",
                                  "\\begin{minipage}{.[page.width[2]]}\n ", sep = "\n")

            minipage_end = "\n\\end{minipage}\n\\end{minipage}"
        }

        intro = .dsb0("\\documentclass[varwidth=.[w], border={ 10 5 10 5 }]{standalone}",
                      "\\usepackage[dvipsnames,table]{xcolor}",
                      "\\usepackage{.[/array, booktabs, multirow, helvet, amsmath, amssymb]}",
                      "\\renewcommand{\\familydefault}{\\sfdefault}", vectorize = TRUE)

        tex_pkg = c()
        tex_pkg = add_pkg("(row|cell)color", x, tex_pkg, "colortbl", fixed = FALSE)
        tex_pkg = add_pkg("threeparttable", x, tex_pkg, opt = "[flushleft]")
        tex_pkg = add_pkg("adjustbox", x, tex_pkg)
        tex_pkg = add_pkg("tabularx", x, tex_pkg)
        tex_pkg = add_pkg("makecell", x, tex_pkg)

        do_rerun = FALSE
        if(any(grepl("tikz", x, fixed = TRUE))){
            # If tikz is present, we need to run pdflatex twice to make it work properly
            do_rerun = TRUE
            tex_pkg = c(tex_pkg, "\\usepackage{tikz}",
                        "\\usetikzlibrary{matrix, shapes, arrows, fit, tikzmark}")
        }

        doc_full = c(intro, tex_pkg,
                     "\n\n\\begin{document}\n", minipage_start,
                     x,
                     minipage_end, "\n\\end{document}\n")

        #
        # Compiling + exporting
        #

        tex_file = file("etable.tex", "w", encoding = "UTF-8")
        writeLines(doc_full, tex_file)
        close(tex_file)

        options(fixest_log_dir = dir)

        # We compile the document
        draft = if(do_rerun) "-draftmode" else ""

        # we delete the aux that can generate problems
        if(file.exists("etable.aux")){
            unlink("etable.aux")
        }

        # I keep the CMD because it is faster
        # WITH: 2.5s
        # SANS: 3.5s

        DO_CMD = TRUE

        warn_msg = NULL
        ok_cmd_tex = TRUE && DO_CMD
        if(DO_CMD){
            outcome = suppressWarnings(system2("pdflatex", .dsb("-halt-on-error -interaction=nonstopmode .[draft] etable.tex"),
                              "etable_shell_pdf.log", "etable_shell_pdf.err"))

            if(outcome == 127){
                ok_cmd_tex = FALSE
            } else if(outcome == 1){

                # Sometimes recompiling works!!!
                outcome = suppressWarnings(system2("pdflatex", .dsb("-halt-on-error -interaction=nonstopmode .[draft] etable.tex"),
                                  "etable_shell_pdf.log", "etable_shell_pdf.err"))

                if(outcome == 1){
                    warn_msg = "pdflatex: error when compiling -- sorry! Check the log file with log_etable('pdflatex')."
                    ok_cmd_tex = FALSE
                }
            }
        }

        if(ok_cmd_tex && do_rerun){
            outcome = suppressWarnings(system2("pdflatex", "-halt-on-error -interaction=nonstopmode etable.tex",
                              "etable_shell_pdf.log", "etable_shell_pdf.err"))

            # No reason for the 2nd run to fail, but I add it anyway just to be safe
            if(outcome == 1){
                warn_msg = "pdflatex: error when compiling -- sorry! Check the log file with log_etable('pdflatex')."
                ok_cmd_tex = FALSE
            }
        }

        if(!ok_cmd_tex){
            # We use tinytex
            if(!requireNamespace("tinytex", quietly = TRUE)){

                if(!is.null(warn_msg)){
                    # we still give the warning to the user due to CMD pblm
                    warn_up(warn_msg)
                    return(NULL)
                }

                warn_up("The functionality you want to use requires the package 'tinytex' which is not installed or a working pdflatex installation which wasn't found.")
                return(NULL)
            }

            info_compile = try(suppressMessages(tinytex::pdflatex("etable.tex",
                                                                  clean = FALSE, min_times = 1 + do_rerun)))
            if(inherits(info_compile, "try-error")){
                warn_up("pdflatex: error when compiling -- sorry! Check the log file with log_etable('pdflatex').")
                return(NULL)
            }

        }

        #
        # Creating the PNG
        #

        ok_cmd_magick = TRUE && DO_CMD
        warn_msg = NULL
        if(DO_CMD){
            outcome = suppressWarnings(system2("magick", .dsb("-density 600 etable.pdf -colorspace RGB .[png_name]"),
                              "etable_shell_magick.log", "etable_shell_magick.err"))

            if(outcome == 127){
                ok_cmd_magick = FALSE

            } else if(outcome == 1){
                warn_msg = "magick: error when converting pdf to png. Check install of ghostscript? Check the log file with log_etable('magick')."
                ok_cmd_magick = FALSE
            }
        }


        if(!ok_cmd_magick){
            if(!requireNamespace("pdftools", quietly = TRUE)){

                if(!is.null(warn_msg)){
                    warn_up(warn_msg)
                    return(NULL)
                }

                warn_up("The functionality you want to use requires the package 'pdftools' which is not installed or a working imagemagick+ghostscript installation which wasn't found.")
                return(NULL)
            }

            # We use pdftools, but to avoid an ugly warning, we need to provide a name we don't want
            # and hence rename the file later...

            pdftools::pdf_convert("etable.pdf", dpi = 600, verbose = FALSE,
                                  filenames = .dsb(".[png_name]_%d.%s"))
            old_name = .dsb(".[png_name]_1.png")
            file.rename(old_name, png_name)
        }

    }

    # Now porting to the viewer

    if(view){
        my_viewer = getOption("viewer")
        if(is.null(my_viewer)){
            warning("To preview the table, we need RStudio's viewer -- which wasn't found.")
        } else {
            # setting up the html document

            # NOTE: all the viewer's data must be in tempdir
            # => when we cache, that can cause problems
            # hence we copy the file there if necessary

            tmp_dir = normalizePath(tempdir(), "/")

            if(normalizePath(getwd(), "/") != tmp_dir){
                old_name = png_name
                png_name = gsub(".+/", "", png_name)
                file.copy(old_name, file.path(tmp_dir, png_name))
                setwd(tmp_dir)
            }

            html_file = viewer_html_template(png_name)

            writeLines(html_file, "etable.html")

            my_viewer("etable.html")
        }
    }

    # And exporting
    if(!is.null(export)){
        file.copy(png_name, export_path, overwrite = TRUE)
    }

    if(!is.null(export_markdown) && export_markdown != png_name){
        file.copy(png_name, export_markdown, overwrite = TRUE)
    }

    return(export_markdown)
}


check_set_path = function(x, type = "", create = TRUE, up = 0){
    # type:
    # + r: read (file or dir must exists), w (file is to be created)
    # + dir: directory and not a document
    # create:
    # - if file: creates the parent dir if the grand parent exists
    # - if dir: creates the dir only if grand parent exists

    set_up(up + 1)

    flags = strsplit(type, ", *")[[1]]

    x_dp = deparse(substitute(x))

    path = try(normalizePath(x, "/", mustWork = FALSE))
    if("try-error" %in% class(path)){
        path = try(normalizePath(paste0("./", x), "/", mustWork = FALSE))
        if("try-error" %in% class(path)){
            stop_up("The path ", x, " is not valid, please revise.")
        }
    }

    # If path exists: fine!
    is_dir = "dir" %in% flags
    if(is_dir){
        if(dir.exists(path)){
            return(path)
        }
    } else if(file.exists(path)){
        return(path)
    }

    # Now the path does not exist already, so we don't need to check

    # here: it must be an error
    if("r" %in% flags){
        msg = if("dir" %in% flags) "directory" else "file"
        stop_up("Argument '", x_dp, "' should be a path to a ", msg,
                " that exists. \n  Problem: '", path, "' does not exist.")
    }

    # Here we're in write

    file_name = gsub(".+/", "", path)
    path_dir = str_trim(path, -nchar(file_name))
    if(nchar(path_dir) == 0) path_dir = "."

    path_parent = dirname(path)
    if(dir.exists(path_parent)){
        if(is_dir && create){
            dir.create(path)
        }

        return(path)
    }

    if(create){
        path_grand_parent = dirname(path_parent)
        if(dir.exists(path_grand_parent)){
            dir.create(path_parent)
            if(is_dir){
                dir.create(path)
            }

            return(path)
        }
    }

    msg = if("dir" %in% flags) "directory" else "file"
    stop_up("Argument '", x_dp, "' should be a path to a ", msg, ". \n  Problem: '",
            path_parent, "' does not exist.")

}

viewer_html_template = function(png_name){
    # I really wanted to see the full table all the time, so I had to add some JS.
    # There must be some straightforward way in CSS, but I don't know it...
    .dsb0('
<!DOCTYPE html>
<html> <head>

<style>

#container {
 width: 100%;
 display: block;
 height: 96vh;
 text-align: center;
}

img {
  max-width: 100%;
  max-height: 100%;
}

</style>

</head>

<body>

<div id="container" class = "etable">
	<img src = ".[png_name]" alt="etable preview">
</div>

</body> </html>
')
}

#' @rdname etable
log_etable = function(type = "pdflatex"){
    check_arg_plus(type, "match(pdflatex, magick, tex, dir)")

    dir = getOption("fixest_log_dir")

    if(length(dir) == 0){
        return("No log currently exists")
    }

    if(type == "pdflatex"){
        path = file.path(dir, "etable.log")
    } else if(type == "magick"){
        path = file.path(dir, "etable_shell_magick.log")
    } else if(type == "dir"){
        return(dir)
    } else {
        path = file.path(dir, "etable.tex")
    }

    if(!file.exists(path)){
        message(.dsb("No log currently exists for '.[type]'."))
        return(invisible(NULL))
    }

    res = readLines(path)

    class(res) = "etable_tex"
    res
}

DIR_EXISTS = function(x){
    if(length(x) != 1 || is.na(x) || !is.character(x)){
        return(FALSE)
    }

    dir.exists(x)
}


fix_pkgwdown_path = function(){
    # https://github.com/r-lib/pkgdown/issues/1218
    # just because I use google drive... it seems pkgdown cannot convert to relative path...

    # This is to ensure it only works for me
    if(!isTRUE(renvir_get("fixest_ROOT"))) return(NULL)
    # we check we're in the right directory (otherwise there can be prblms with Rmakdown)
    if(!isTRUE(file.exists("R/VCOV_aliases.R"))) return(NULL)

    all_files = list.files("docs/articles/", full.names = TRUE, pattern = "html$")

    for(f in all_files){
        my_file = file(f, "r", encoding = "UTF-8")
        text = readLines(f)
        close(my_file)
        if(any(grepl("../../../", text, fixed = TRUE))){
            # We embed the images directly: safer

            # A) we get the path
            # B) we transform to URI
            # C) we replace the line

            pat = "<img.+\\.\\./.+/fixest/.+/images/"
            qui = which(grepl(pat, text))
            for(i in qui){
                # ex: line = "<img src = \"../../../Google drive/fixest/fixest/vignettes/images/etable/etable_tex_2021-12-02_1.05477838.png\">"
                line = text[i]
                line_split = strsplit(line, "src *= *\"")[[1]]
                path = gsub("\".*", "", line_split[2])
                # ROOT is always fixest
                path = gsub(".+fixest/", "", path)
                path = gsub("^articles", "vignettes", path)

                URI = knitr::image_uri(path)

                rest = gsub("^[^\"]+\"", "", line_split[2])
                new_line = dsb('.[line_split[1]] src = ".[URI]".[rest]')

                text[i] = new_line
            }

            my_file = file(f, "w", encoding = "UTF-8")
            writeLines(text, f)
            close(my_file)
        }
    }

}

####
#### Utilities ####
####


style_apply = function(coef.style, coef_mat, coef_names){
    # applies a style to specified coefficients

    set_up(2)

    if(length(coef.style) == 0){
        return(coef_mat)
    }

    res = coef_mat
    n_models = ncol(res) - 1

    has_row_se = any(" " %in% res[, 1])
    get_row_id = function(i) if(has_row_se) 1 + 2 * (i - 1) else i

    for(i in seq_along(coef.style)){
        cs = coef.style[[i]]
        cs_name = paste0(trimws(names(coef.style)[i]), " ")

        if(!grepl(":coef(_se)?:", cs_name)){
            warn_up("In the argument 'coef.style', the names must contain the ':coef:' (or :coef_se:) string. ",
                    "Problem: this is not the case for the ", n_th(i), " element (equal to '", cs_name, "').")
        }

        is_se = grepl(":coef_se:", cs_name)
        if(is_se){
            cs_name = sub(":coef_se:", ":coef:", cs_name, fixed = TRUE)
        }
        add_se = function(x) if(is_se) rbind(x, x + matrix(c(1, 0), nrow(x), 2, byrow = TRUE)) else x

        style_split = strsplit(cs_name, ":coef:", fixed = TRUE)[[1]]

        if(cs %in% c(".", "all")){

            id = which(matrix(grepl("[^ ]", coef_mat), nrow(coef_mat)), arr.ind = TRUE)
            # we remove the row names
            id = id[id[, 2] != 1, ]
            if(!is_se && has_row_se){
                # we remove the se row!
                id = id[id[, 1] %% 2 == 1, ]
            }

            # We apply the style
            values = res[id]
            new_values = paste0(style_split[1])
            for(k in seq_along(style_split)[-1]){
                new_values = paste0(new_values, values, style_split[k])
            }
            new_values = trimws(new_values)
            res[id] = escape_latex(new_values)

        } else {
            coef_location_all = coef_location(cs, coef_names, n_models, "coef.style", pool = TRUE)
            cell = coef_location_all$cell

            # cell: list of n x 2 matrices
            for(j in seq_along(cell)){
                cell_mat = cell[[j]]

                id = cbind(get_row_id(cell_mat[, 1]), cell_mat[, 2] + 1)
                id = add_se(id)

                id = id[grepl("[^ ]", res[id]), , drop = FALSE]

                if(length(id) > 0){
                    values = res[id]
                    new_values = paste0(style_split[1])
                    for(k in seq_along(style_split)[-1]){
                        new_values = paste0(new_values, values, style_split[k])
                    }
                    new_values = trimws(new_values)
                    res[id] = escape_latex(new_values)
                }
            }
        }
    }

    return(res)
}

highlight_apply = function(highlight, coef_mat, coef_names){
    # applies coefficients highlighting

    set_up(2)

    if(length(highlight) == 0){
        return(list(preamble = "", coef_mat = coef_mat))
    }

    if(is.null(names(highlight))){
        names(highlight) = " "
    }

    # default color
    COLOR = "#c80815"

    res = coef_mat
    n_models = ncol(res) - 1

    # \definecolor{Mycolor2}{HTML}{00F9DE}
    # \definecolor{mypink1}{rgb}{0.858, 0.188, 0.478}
    # \definecolor{mypink2}{RGB}{219, 48, 122}

    has_row_se = any(" " %in% res[, 1])
    get_row_id = function(i) if(has_row_se) 1 + 2 * (i - 1) else i

    # frame id
    f_id = 0

    all_colors = c()
    tex_preamble = c()

    for(i in seq_along(highlight)){
        hl = highlight[[i]]
        hl_name = trimws(strsplit(names(highlight)[i], "[, ]+")[[1]])
        if(identical(hl_name, '')){
            hl_name = character(0)
        }

        #
        # Parsing the parameters
        #

        # TRUE/FALSE flags
        is_se = "se" %in% hl_name && has_row_se
        add_se = function(x) if(is_se) rbind(x, x + matrix(c(1, 0), nrow(x), 2, byrow = TRUE)) else x

        is_rowcol = "rowcol" %in% hl_name
        is_round = !"square" %in% hl_name

        hl_name = setdiff(hl_name, c("se", "rowcol", "square"))

        # Values
        thick_raw = grep("^thick\\d$", hl_name, value = TRUE)
        thick = if(length(thick_raw) == 0) 6 else as.numeric(sub("thick", "", thick_raw))
        if(thick > 6){
            stop_up("In the argument 'highlight', the name of the ", n_th(i),
                 " element (equal to '", hl_name, "') is ill formed. ",
                 "The item 'thick' can go only from 1 to 6.")
        }

        sep_raw = grep("^sep\\d+$", hl_name, value = TRUE)
        sep = if(length(sep_raw) == 0) 3 else sub("sep", "", sep_raw)

        hl_name = setdiff(hl_name, c(thick_raw, sep_raw))

        # Color: the remaining
        color = hl_name
        is_col = length(color) == 1

        if(length(color) > 1){
            stop_up("In the argument 'highlight', the name of the ", n_th(i),
                 " element (equal to '", hl_name, "') is ill formed. ",
                 "It should be a comma separated list of options, which include: 'rowcol', 'square', 'thickd' (with d from 0 to 6), 'sepd' (d: 0-9), 'se', and 'color!alpha' with 'color' a valid R color and, optionnaly, 'alpha' in 0-100.")
        }

        # conversion of the color
        is_alpha = FALSE
        if(is_col){
            if(grepl("!", color, fixed = TRUE)){
                color_split = strsplit(color, "!", fixed = TRUE)[[1]]
                is_alpha = TRUE
                alpha = color_split[2]
                color = color_split[1]
            }
        } else {
            color = COLOR
            if(is_rowcol){
                # Default row color is not solid
                alpha = "25"
                is_alpha = TRUE
            }
        }

        tex_color = error_sender(col2rgb(color),
                                 "In the argument 'highlight', the color, in the name of the ", n_th(i),
                                 " element (equal to '", color, "'), could not be converted to RGB. ")
        tex_color = paste0(as.vector(tex_color), collapse = ", ")

        if(tex_color %in% all_colors){
            tag_color = names(all_colors)[all_colors == tex_color]

        } else {
            tag_color = paste0("col", tag_gen())
            all_colors[tag_color] = tex_color
            tex_preamble = c(tex_preamble,
                             paste0("\\definecolor{", tag_color, "}{RGB}{", tex_color, "}"))
        }

        if(is_alpha){
            tag_color = paste0(tag_color, "!", alpha)
        }

        coef_location_all = coef_location(hl, coef_names, n_models, "highlight")
        cell = coef_location_all$cell
        row = coef_location_all$row
        range = coef_location_all$range

        if(is_rowcol){
            # cell: list of n x 2 matrices
            for(j in seq_along(cell)){
                cell_mat = cell[[j]]
                for(k in 1:nrow(cell_mat)){
                    my_cell = cell_mat[k, ]

                    id = cbind(get_row_id(my_cell[1]), my_cell[2] + 1)
                    id = add_se(id)

                    res[id] = paste0("\\cellcolor{", tag_color, "} ", res[id])
                }
            }

            # row: list of integers
            for(j in seq_along(row)){
                my_row = row[[j]]

                id = cbind(get_row_id(my_row), 1)
                id = add_se(id)

                res[id] = paste0("\\rowcolor{", tag_color, "} ", res[id])
            }

            # range: list of vectors of length 4 (i1, k1, i2, k2)
            for(j in seq_along(range)){
                my_range = range[[j]]

                id = expand.grid(my_range[1]:my_range[3], my_range[2]:my_range[4])
                id = add_se(id)

                res[id] = paste0("\\cellcolor{", tag_color, "} ", res[id])
            }

        } else {
            #
            # FRAME
            #

            tag_frame = tag_gen()

            tag_frame_NW = paste0("\\markNW", tag_frame)
            tag_frame_SE = paste0("\\markSE", tag_frame)

            # thickness
            # thin, very thin, ultra thin
            # thick, very thick, ultra thick
            thickness = .dsb(".[/thin, very thin, ultra thin, thick, very thick, ultra thick]")[thick]
            rounded = if(is_round) "rounded corners, " else ""

            # On the frame:
            # https://tex.stackexchange.com/questions/240542/adding-a-rectangular-box-with-tikz-to-table-beamer
            #

            tex_preamble = c(tex_preamble, .dsb(
                "
\\newcommand.[tag_frame_NW][1]{%
   \\tikz[overlay,remember picture]
      \\node (marker-#1-a) at (0, .3em) {};%
}

\\newcommand.[tag_frame_SE][1]{%
   \\tikz[overlay,remember picture]
      \\node (marker-#1-b) at (0, .3em) {};%
   \\tikz[overlay,remember picture,inner sep=.[sep]pt, .[thickness]]
      \\node[draw=.[tag_color],.[rounded]fit=(marker-#1-a.north west) (marker-#1-b.south east)] {};%
}\n"))


            # cell: list of n x 2 matrices
            for(j in seq_along(cell)){
                cell_mat = cell[[j]]
                n_cells = nrow(cell_mat)
                k_done = c()
                for(k in 1:n_cells){

                    if(k %in% k_done) next

                    my_cell = cell_mat[k, ]

                    # we find the rightmost index
                    cell_right = my_cell
                    while(k + 1 <= n_cells &&
                          cell_mat[k + 1, 1] == my_cell[1] &&
                          cell_mat[k + 1, 2] == my_cell[2] + 1){
                        k_right = k + 1
                        cell_right = cell_mat[k_right, ]
                        k_done[length(k_done) + 1] = k_right
                        k = k + 1
                    }

                    id_NW = cbind(get_row_id(my_cell[1]), my_cell[2] + 1)

                    if(is_se){
                        id_SE = cbind(get_row_id(my_cell[1]) + 1, cell_right[2] + 1)
                    } else {
                        id_SE = cbind(get_row_id(my_cell[1]), cell_right[2] + 1)
                    }

                    f_id = f_id + 1

                    res[id_NW] = .dsb(".[tag_frame_NW]{f.[f_id]} .[res[id_NW]]")
                    res[id_SE] = .dsb(".[res[id_SE]] .[tag_frame_SE]{f.[f_id]}")
                }
            }

            # row: list of integers
            for(j in seq_along(row)){
                my_row = row[[j]]

                id_NW = cbind(get_row_id(my_row), 2)

                if(is_se){
                    # We need to end at the line just below is is_se
                    id_SE = cbind(get_row_id(my_row) + 1, n_models + 1)
                } else {
                    id_SE = cbind(get_row_id(my_row), n_models + 1)
                }

                f_id = f_id + 1

                res[id_NW] = .dsb(".[tag_frame_NW]{f.[f_id]} .[res[id_NW]]")
                res[id_SE] = .dsb(".[res[id_SE]] .[tag_frame_SE]{f.[f_id]}")

            }

            # range: list of vectors of length 4 (i1, k1, i2, k2)
            for(j in seq_along(range)){
                my_range = range[[j]]

                id_NW = cbind(get_row_id(my_range[1]), my_range[2] + 1)
                id_SE = cbind(get_row_id(my_range[3]), my_range[4] + 1)

                if(is_se){
                    # We need to end at the line just below is is_se
                    id_SE[1] = id_SE[1] + 1
                }

                f_id = f_id + 1

                res[id_NW] = .dsb(".[tag_frame_NW]{f.[f_id]} .[res[id_NW]]")
                res[id_SE] = .dsb(".[res[id_SE]] .[tag_frame_SE]{f.[f_id]}")
            }
        }
    }

    return(list(coef_mat = res, preamble = tex_preamble))
}


coef_location = function(x, coef_names, n_models, arg_name, pool = FALSE){
    # x: vector of coef location
    # cell: "x1@2", "x1@1, 3, 5-8"
    # range: "x2@2; x3@.N", "x2; x3"
    # row: "x1"

    cell = range = row = list()

    for(i in seq_along(x)){
        xi = x[i]

        if(grepl(";", xi, fixed = TRUE)){
            # => range
            xi_split = trimws(strsplit(xi, ";", fixed = TRUE)[[1]])

            if(length(xi_split) != 2){
                stop_up(up = 3, "In the argument '", arg_name, "', the value of the ", n_th(i),
                        " element (equal to '", xi, "') is not valid. ",
                        "It should contain at most one semi-comma.")
            }

            xi_1 = xi_split[[1]]
            xi_2 = xi_split[[2]]

            xi_1_loc = coef_pos_parse(xi_1, coef_names, n_models, i, arg_name)
            xi_2_loc = coef_pos_parse(xi_2, coef_names, n_models, i, arg_name)

            # We always normalize the range, it must be two single positions
            xi_1_pos = if(xi_1_loc$is_row) 1 else min(xi_1_loc$pos)
            xi_2_pos = if(xi_2_loc$is_row) 1 else max(xi_2_loc$pos)

            if(pool){
                cell[[length(cell) + 1]] = expand.grid(xi_1_loc$coef:xi_2_loc$coef,
                                                       xi_1_pos:xi_2_pos)
            } else {
                range[[length(range) + 1]] = c(xi_1_loc$coef, xi_1_pos, xi_2_loc$coef, xi_2_pos)
            }

        } else {
            xi_loc = coef_pos_parse(xi, coef_names, n_models, i, arg_name)

            if(xi_loc$is_row){
                if(pool){
                    all_cells = cbind(xi_loc$coef, 1:n_models)
                    cell[[length(cell) + 1]] = all_cells
                } else {
                    row[[length(row) + 1]] = xi_loc$coef
                }
            } else {
                all_cells = cbind(xi_loc$coef, xi_loc$pos)
                cell[[length(cell) + 1]] = all_cells
            }
        }
    }


    res = list(cell = cell, row = row, range = range)
    return(res)
}

coef_pos_parse = function(x, coef_names, n_models, i, arg_name){
    # x: "x1@2", "x1", "x1@2, 5, 7-10, 13"

    set_up(4)

    x_split = strsplit(x, "@", fixed = TRUE)[[1]]
    if(length(x_split) > 2){
        stop_up("In the argument '", arg_name, "', the location of the ", n_th(i),
                " element (equal to '", x, "') is not valid. ",
                "It should contain at most one '@'.")
    }

    x_coef = x_split[1]
    x_pos = if(length(x_split) == 2) x_split[2] else NULL

    #
    # Position
    #

    x_pos = gsub(".N", n_models, x_pos, fixed = TRUE)

    if(length(x_pos) > 0){
        # x: 3, 5, 6-8, 10

        if(grepl("[^ ,-[[:digit:]]]", x_pos)){
            stop_up("In the argument '", arg_name, "', the location of the ", n_th(i),
                    " element (equal to '", x, "') contains non valid characters. Please have a look at the help/example.")
        }

        # new DSB power
        x_txt = paste0("c(", gsub("-", ":", x_pos), ")")
        x_call = error_sender(str2lang(x_txt), "In argument '", arg_name,
                              "' the position in '", x,
                              "' is not valid. Please have a look at the help/example.",
                              up = 4)

        x_pos = eval(x_call)
    }

    #
    # Coefficient
    #

    coef_id = pmatch_varname(x_coef, coef_names, arg_name)

    res = list(coef = coef_id, pos = x_pos, is_row = length(x_pos) == 0)

    return(res)
}

pmatch_varname = function(x, coef_names, arg_name){
    # We want a SINGLE match
    # coef_names: vector of names
    # names(coef_names): original names
    # coef_names = setNames(c("Petal length", "Petal width"), c("Petal.Length", "Petal.Width"))
    # x = "Petal.L"

    is_regex = is_original = FALSE
    if(grepl("^@", x)){
        is_regex = TRUE
        x = str_trim(x, 1)
        if(grepl("^%", x)){
            is_original = TRUE
            x = str_trim(x, 1)
        }
    } else if(grepl("^%", x)){
        is_original = TRUE
        x = str_trim(x, 1)
        if(grepl("^@", x)){
            is_regex = TRUE
            x = str_trim(x, 1)
        }
    }

    # Pattern of match finding
    if(is_regex && is_original){
        do_regex = TRUE
        do_original = TRUE
    } else if(is_regex){
        do_regex = c(TRUE, TRUE)
        do_original = c(FALSE, TRUE)
    } else if(is_original){
        do_regex = c(FALSE, TRUE)
        do_original = c(TRUE, TRUE)
    } else {
        do_regex = c(FALSE, FALSE, TRUE, TRUE)
        do_original = c(FALSE, TRUE, FALSE, TRUE)
    }

    ok = FALSE
    for(i in seq_along(do_regex)){

        if(do_regex[i]){
            if(do_original[i]){
                qui = grepl(x, names(coef_names))
            } else {
                qui = grepl(x, coef_names)
            }

            if(sum(qui) == 1){
                ok = TRUE
                qui = which(qui)
            }
        } else {
            if(do_original[i]){
                qui = charmatch(x, names(coef_names))
            } else {
                qui = charmatch(x, coef_names)
            }

            if(!is.na(qui) && qui != 0) ok = TRUE
        }

        if(ok) break
    }

    if(!ok){
        stop_up(up = 5, "In the argument '", arg_name, "', the value '", x,
                "' does not match any variable name.")
    }

    qui
}

# x = strsplit("*bonjour $x^2$*", "$", fixed = TRUE)[[1]]
# x = strsplit("*bonjour $x^2*3$*", "$", fixed = TRUE)[[1]]
# x = "et **bonsoir**!"
markup_apply = function(x){
    # x: if len > 1, means it contains equations that have been split
    # the function is slow, but we don't apply it to many things anyway, so that's OK

    if(!any(grepl("*", x, fixed = TRUE))){
        return(x)
    }

    is_eq = is_eq_star = FALSE
    n_x = length(x)
    if(n_x > 1){
        is_eq = TRUE
        id_eq = (1:n_x)[(1:n_x) %% 2 == 0]
        is_eq_star = any(grepl("*", x[id_eq], fixed = TRUE))
        if(is_eq_star){
            x[id_eq] = gsub("*", "_@_", x[id_eq], fixed = TRUE)
        }

        x = paste0(x, collapse = "|$|")
    }

    res = x

    # patterns
    dict_pat = c("***" = "\\textbf{\\textit{___}}",
                 "**" = "\\textbf{___}",
                 "*" = "\\textit{___}")
    failed = rep(FALSE, 3)

    for(i in 1:3){
        pat = names(dict_pat)[i]

        if(grepl(pat, x, fixed = TRUE)){
            markup = dict_pat[i]

            n_match = length(gregexpr(pat, x, fixed = TRUE)[[1]])

            if(!n_match %% 2 == 0){
                # failed markup
                res = gsub(pat, .dsb("__.[`i`*c!$]__"), res, fixed = TRUE)
                failed[i] = TRUE
            } else {
                x_split = strsplit(res, pat, fixed = TRUE)[[1]]
                n = length(x_split)
                for(j in (1:n)[(1:n) %% 2 == 0]){
                    x_split[j] = sub("___", x_split[j], markup, fixed = TRUE)
                }
                res = paste0(x_split, collapse = "")
            }
        }
    }

    for(i in which(failed)){
        pat = names(dict_pat)[i]
        res = gsub(.dsb("__.[`i`*c!$]__"), pat, res, fixed = TRUE)
    }

    if(is_eq){
        if(is_eq_star){
            res = gsub("_@_", "*", res, fixed = TRUE)
        }
        res = strsplit(res, "|$|", fixed = TRUE)[[1]]
    }

    res
}

escape_all = function(x){
    # we escape all
    res = gsub("((?<=[^\\\\])|(?<=^))(\\$|_|%|&|\\^|#)", "\\\\\\2", x, perl = TRUE)
    res
}

# escape_latex("Voici **une** *equation*: $5! = 5*4*3*2*1$. Est-ce que ***ca marche***?$^*$")
escape_latex = function(x_all, makecell = TRUE){
    # This is super tricky to escape properly!
    # We do NOT escape within equations

    x_name = deparse(substitute(x_all))

    res = c()

    for(index in seq_along(x_all)){
        x = x_all[index]

        # 1) finding out equations, ie non escaped dollar signs
        dollars = gregexpr("((?<=[^\\\\])|(?<=^))\\$", x, perl = TRUE)[[1]]

        is_eq = FALSE
        if(length(dollars) > 1){
            is_eq = TRUE
            if(length(dollars) %% 2 != 0){
                stop_up(up = 2, "There are ", length(dollars), " dollar signs in the following character string:\n", x, "\nIt will raise a Latex error (which '$' means equation? which means dollar-sign?): if you want to use a regular dollar sign, please escape it like that: \\\\$.")
            }
        }

        # 2) Escaping but conditionally on not being in an equation
        if(is_eq){
            # Finding out the equations
            all_items = strsplit(paste0(x, " "), "((?<=[^\\\\])|(?<=^))\\$", perl = TRUE)[[1]]
            for(i in seq_along(all_items)){
                if(i %% 2 == 1){
                    all_items[i] = escape_all(all_items[i])
                }
            }

            all_items = markup_apply(all_items)

            res[index] = gsub(" $", "", paste(all_items, collapse = "$"))
        } else {
            res[index] = markup_apply(escape_all(x))
        }
    }

    if(makecell){
        who = grepl("\n", res, fixed = TRUE)
        if(any(who)){
            res[who] = paste0("\\makecell{", gsub("\n", "\\\\", res[who], fixed = TRUE), "}")
        }
    }

    if(!is.null(names(x_all))){
        names(res) = names(x_all)
    }

    res
}

format_se_type = function(x, width, by = FALSE){
    # we make 'nice' se types
    # format_se_type("Clustered (species & fe2)", 10, by = TRUE)
    # format_se_type("vcovHC(x, type = \"HC0\")", 10, by = TRUE)
    # format_se_type("Newey-West (L=10)", 10, by = TRUE)
    # format_se_type("Driscoll-Kraay (L=10)", 10, by = TRUE)
    # format_se_type("Conley (100km)", 10, by = TRUE)

    if(!grepl("\\(", x) || !grepl("Clustered", x, fixed = TRUE)){
        # means not clustered
        if(nchar(x) <= width) return(x)

        # Special case: sandwich
        if(grepl("^vcov[^\\(]+\\(", x)){
            x = gsub("\\(x, ", "(", x)
            x = gsub("\\(\\)", "", x)

            if(nchar(x) > width){
                # we still reduce it
                x = gsub(" = ", "=", x)
            }

            if(nchar(x) > width){
                # ...even further
                x = paste0(substr(x, 1, width - 2), "..")
            }

            return(x)
        }

        # We reduce each word to 3 letters (if needed)
        x_split = c("$", strsplit(x, "")[[1]]) # we add a non-letter flag, marking the beginning
        x_split_new = x_split
        end_word = length(x_split)
        non_letter_flag = grepl("[^[:alpha:]]", x_split) * (1:end_word)
        letter_flag = grepl("[[:alpha:]]", x_split) * (1:end_word)
        while(TRUE){
            start_word = which.max(non_letter_flag[1:end_word]) + 1
            # we truncate
            word_length = end_word - start_word + 1
            slack = length(x_split_new) - (width + 1)
            letters_to_rm = min(word_length - 4, slack)
            if(letters_to_rm > 0){
                i_max = end_word - letters_to_rm
                x_split_new = x_split_new[-((i_max+1):end_word)]
                x_split_new[i_max] = "."
            }

            lf = letter_flag[1:(start_word - 1)]
            if(all(lf == 0)) break

            # new end_word
            end_word = which.max(lf)
        }

        return(paste(x_split_new[-1], collapse = ""))
    } else if(x == "NA (not-available)"){
        return("not available")
    }

    # Now the FEs
    all_fe = gsub(".+\\((.+)\\)", "\\1", x)

    all_fe_split = gsub(" ", "", strsplit(all_fe, "&")[[1]])
    n_fe = length(all_fe_split)
    n_char = nchar(all_fe_split)

    if(n_fe == 1 && !grepl("\\^", all_fe_split[1])){
        if(by){
            se_formatted = paste0("by: ", all_fe_split[1])
        } else {
            se_formatted = paste0("1-way: ", all_fe_split[1])
        }

        if(nchar(se_formatted) > width){
            se_formatted = paste0(substr(se_formatted, 1, width - 2), "..")
        }
        return(se_formatted)
    }

    nb = ifelse(by, 3, 6)

    if(width < nb + sum(n_char) + (n_fe-1) * 3){
        qui = n_char > 5
        for(i in which(qui)){
            if(grepl("\\^", all_fe_split[i])){
                single_split = strsplit(all_fe_split[i], "\\^")[[1]]
                qui_bis = nchar(single_split) > 4
                single_split[qui_bis] = paste0(substr(single_split[qui_bis], 1, 3), ".")
                all_fe_split[i] = paste(single_split, collapse = "^")
            } else {
                all_fe_split[i] = paste0(substr(all_fe_split[i], 1, 4), ".")
            }
        }
    }

    if(by){
        se_formatted = paste0("by: ", paste(all_fe_split, collapse = " & "))
    } else {
        se_formatted = paste0(n_fe, "-way: ", paste(all_fe_split, collapse = " & "))
    }


    # NOTA:
    # we do not trim if still too large because the SE-type IS informative!
    # A table without that information is useless, it's trimmed enough already

    # if(nchar(se_formatted) > width){
    #     # se_formatted = gsub("-way: ", "way: ", se_formatted)
    #     se_formatted = paste0(substr(se_formatted, 1, width - 2), "..")
    # }

    se_formatted
}

format_se_type_latex = function(x, dict = c(), inline = FALSE){
    # we make 'nice' se types

    if(!grepl("\\(", x)){
        # means not clustered
        # we escape all
        return(escape_all(x))
    }

    # Now the FEs
    main_type = gsub(" \\(.*", "", x)
    all_fe = gsub(".+\\((.+)\\)", "\\1", x)

    all_fe_split = gsub(" ", "", strsplit(all_fe, "&")[[1]])
    n_fe = length(all_fe_split)

    # Renaming the FEs

    all_fe_format = c()
    for(i in 1:length(all_fe_split)){
        fe = all_fe_split[i]

        if(fe %in% names(dict)){
            all_fe_format[i] = dict[fe]
        } else if(grepl("\\^", fe)){
            fe_split = strsplit(fe, "\\^")[[1]]
            who = fe_split %in% names(dict)
            fe_split[who] = dict[fe_split[who]]
            all_fe_format[i] = paste(fe_split, collapse = "-")
        } else {
            all_fe_format[i] = fe
        }
    }

    fe_format = paste(all_fe_format, collapse = " \\& ")

    # We add some flexibility: anticipation of more VCOV types
    main_type_dict = c("Clustered" = "Clustered", "Two-way" = "Clustered",
                       "Three-way" = "Clustered", "Four-way" = "Clustered")

    if(main_type %in% names(main_type_dict)){
        main_type = main_type_dict[main_type]
    }

    if(inline){
        # The fact that it is clustered is deduced
        se_formatted = fe_format
    } else {
        se_formatted = paste0(main_type, " (", fe_format, ")")
    }

    escape_latex(se_formatted)
}

tex_star = function(x){
    qui = nchar(x) > 0
    x[qui] = paste0("$^{", x[qui], "}$")
    x
}




tex_multicol = function(x, add_rule = FALSE){

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

    names_multi_raw = names_multi
    for(i in seq_along(nb_multi)){
        if(nb_multi[i] > 1){
            names_multi[i] = paste0("\\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
        }
    }

    res = paste0(paste(names_multi, collapse = " & "), " \\\\ ")

    if(add_rule){
        my_rule = c()
        start = 2 + c(0, cumsum(nb_multi))
        end = 1 + cumsum(nb_multi)
        for(i in seq_along(nb_multi)){
            if(grepl("[^ ]", names_multi_raw[i])){
                my_rule[i] = paste0("\\cmidrule(lr){", start[i], "-", end[i], "}")
            }
        }

        res = paste0(res, paste0(my_rule, collapse = " "))
    }

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



uniquify_names = function(x){
    # x: vector of names
    # we make each value of x unique by adding white spaces

    if(length(x) == 0) return(NULL)

    x_unik = unique(x)

    if(length(x_unik) == length(x)) return(x)

    x[nchar(x) == 0] = " "

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


extralines_extractor = function(x, name = NULL, tex = FALSE){
    # x must be a one sided formula
    # name: name of the listed element (empty = "")
    # => returns a named list

    is_name = !is.null(name) && nchar(name) > 0

    # extralines registered
    el_default = getOption("fixest_extralines")
    key_registered = names(el_default)

    # fitstat
    fitstat_fun_types = fitstat(give_types = TRUE)
    fitstat_type_allowed = fitstat_fun_types$types
    type_alias = fitstat_fun_types$type_alias

    # The variable(s) requested
    current_vars = attr(terms(x), "term.labels")

    if(length(current_vars) > 1 && is_name){
        stop_up("You cannot give list names in 'extralines' when several values are summoned via a formula. Simply remove the name associated to the formula to make it work (concerns '", name, "').", up = 2)
    }

    valid_keys = c(key_registered, fitstat_type_allowed)

    if(!all(current_vars %in% valid_keys)){
        pblm = setdiff(current_vars, valid_keys)
        stop_up("Argument 'extralines' can be a one-sided formula whose variables refer to the macros registered using the function 'extralines_register' or fit statistics (valid fitstat keywords). Problem: the values", enumerate_items(pblm, "s.were"), " not valid.", up = 2)
    }

    if(length(current_vars) == 0){
        extralines = list()
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



insert = function(x, y, i){
    # x = list(a = 1, b = 2, c = 3) ; y = list(u = 33, v = 55) ; i = 2 ; insert(x, y, i)
    # we insert y into x in location i
    # we don't lose the names!

    mode = mode(x)
    if(!mode %in% c("list", "numeric", "logical", "character", "integer")){
        stop("Internal error: the current mode (", mode, ") is node supported in insert().")
    }

    n_x = length(x)

    if(n_x == 0){
        return(y)
    }

    n_y = length(y)
    res = vector(mode, n_x + n_y)

    names_x = names(x)
    if(is.null(names_x)) names_x = character(n_x)
    names_y = names(y)
    if(is.null(names_y)) names_y = character(n_y)

    if(i > n_x){
        res[1:n_x] = x
        res[n_x + 1:n_y] = y
        names(res) = c(names_x, names_y)

    } else if(i == 1){
        res[1:n_y] = y
        res[n_y + 1:n_x] = x
        names(res) = c(names_y, names_x)

    } else {
        res[1:(i-1)] = x[1:(i-1)]
        res[(i-1) + 1:n_y] = y
        res[(i+n_y):(n_x + n_y)] = x[i:n_x]
        names(res) = c(names_x[1:(i-1)], names_y, names_x[i:n_x])

    }

    return(res)
}



is_fixest_model = function(x){
    any(c("fixest", "fixest_list", "fixest_multi") %in% class(x))
}


expand_list_vector = function(x){
    # we transform list("A" = 2, "B", "C" = 3) into c("A", "A", "B", "C", "C", "C")

    x_names = names(x)

    if(is.null(x_names)){
        # either list("m", "f") or c("a", "a", "c")
        x_new = unlist(x)
    } else {
        # ex: list("M", "F"=2)
        x_new = c()
        for(j in seq_along(x)){

            if(nchar(x_names[j]) == 0){
                x_new = c(x_new, x[j])
            } else {
                x_new = c(x_new, rep(x_names[j], x[j]))
            }
        }
    }

    return(x_new)
}


tex.nice = function(x, n_models){


    x = unlist(strsplit(x, "\n"))

    n = n_models + 1

    #
    # I) dealing with amps
    #

    x_split_amp = strsplit(x, "&", fixed = TRUE)

    qui_amp = lengths(x_split_amp) == n & !sapply(x_split_amp, function(x) any(grepl("midrule", x)))

    mat_amp = matrix(unlist(x_split_amp[qui_amp]), ncol = n, byrow = TRUE)
    mat_amp[, 1:(n-1)] = apply(mat_amp[, 1:(n-1), drop = FALSE], 2, format)

    if(any(grepl("\\", mat_amp, fixed = TRUE))){
        # we fix the \\ problem that will count for one character eventually
        who_slash = which(grepl("\\", mat_amp, fixed = TRUE))
        slash_count = lengths(gregexpr("\\", mat_amp[who_slash], fixed = TRUE))
        mat_amp[who_slash] = paste0(mat_amp[who_slash], sprintf("% *s", slash_count, " "))
    }

    amp_new = apply(mat_amp, 1, paste, collapse = "&")

    x[qui_amp] = amp_new

    #
    # II) dealing with tabs
    #

    # We assume everything is properly formatted
    x_begin = pmax(lengths(strsplit(x, "\\begin{", fixed = TRUE)) - 1, 0)
    x_end = pmax(lengths(strsplit(x, "\\end{", fixed = TRUE)) - 1, 0)

    tabs = pmax(cumsum(x_begin) - cumsum(x_end) - x_begin, 0)

    res = paste0(sprintf("% *s", tabs*3, " "), x)
    res = gsub("^ ([^ ])", "\\1", res)

    res
}

check_set_adjustbox = function(adjustbox, up = 0){
    set_up(up + 1)

    check_arg(adjustbox, "NULL scalar(character, strict logical, numeric) GT{0}")

    if(is.null(adjustbox)){
        # nothing
    } else if(isTRUE(adjustbox)){
        adjustbox = "width = \\textwidth, center"
    } else if(isFALSE(adjustbox)){
        adjustbox = NULL
    } else if(is.numeric(adjustbox)){
        if(adjustbox > 3){
            warn_up("When 'adjustbox' is a number, the unit is the text-width. Hence a value of ", fsignif(adjustbox), " may be too large.")
        }
        adjustbox = paste0("width = ", adjustbox, "\\textwidth, center")
    } else if(grepl("(?i)^[[:digit:]\\.]+ *t(w|h)$", trimws(adjustbox))){
        adj_nbr = gsub("[^[:digit:]\\.]", "", adjustbox)

        if(!is_numeric_in_char(adj_nbr)){
            stop_up("The number in the argument 'adjustbox' (equal to '", adjustbox, "') could not be parsed. Please revise.")
        }

        adj_unit = if(grepl("(?i)tw", adjustbox)) "\\textwidth" else "\\textheight"

        adjustbox = paste0("width = ", adj_nbr, adj_unit, ", center")
    }

    adjustbox
}



tag_gen = function(){
    # What's the fuss here?
    # We want to ensure that the tags are unique even if
    # the seed is reset, that's why we make all that.
    # Why? because otherwise, if two identical tags are used in a tex
    # document it will fail to compile.
    # => granted that's almost impossible, but that's still possible.
    # (ex: two identical code sections with set.seed in Rmarkdown)
    #
    # we only use letters because tex macro names only allow letters

    id = getOption("fixest_tag")
    if(is.null(id)) id = 1
    options(fixest_tag = id + 1)

    items = letters

    # Below: to ensure unicity even if the seed is the same
    n_shifts = id %% 26
    n_reshuffle = id %/% 26

    if(n_reshuffle > 0){
        for(i in 1:n_reshuffle){
            items = sample(items)
        }
    }

    if(n_shifts > 0){
        items = c(items[-(1:n_shifts)], items[1:n_shifts])
    }

    # 26 ** 5 = 11,881,376
    tag = paste0(sample(items, 6, replace = TRUE), collapse = "")

    return(tag)
}



check_set_page_width = function(page.width, up = 0){

    set_up(up + 1)
    check_value_plus(page.width, "match(fit, a4, us) | vector(character, numeric) GT{0} len(1, 2) no NA")

    if(isTRUE(attr(page.width, "no-margin"))){
        # Nothing => already set
    } else if(identical(page.width, "fit")){
        # nothing
    } else if(identical(page.width, "a4")){
        page.width = c("21cm", "17cm")

    } else if(identical(page.width, "us")){
        page.width = c("21.6cm", "16.5cm")

    } else if(is.numeric(page.width)){
        if(length(page.width) == 1){
            page.width = paste0(c(page.width, page.width - 4), "cm")
        } else {
            page.width = paste0(c(page.width[1], page.width[1] - 2*page.width[2]), "cm")
        }
    } else {
        # page.fit is a character scalar
        # format: 21:2cm
        #         21 is the total width
        #         2 is "optional". The margin width of ONE side
        #         cm: the unit

        if(length(page.width) > 1){
            page.width = paste0(page.width, collapse = ", ")
        }

        unit = .dsb("'[[:alpha:] ]+$'x, w, L?page.width")
        if(!unit %in% c("cm", "mm", "in", "pt", "em", "ex")){
            stop_up("The argument 'page.width' must represent a 'Latex' measure, ex: 12cm. ",
                    "It can be optionally followed by a side margin width, as in '21, 2cm'. ",
                    "Problem: a valid unit could not be found in '", unit,"'.")
        }

        # Now the numbers
        numbers = .dsb("'[[:alpha:]]'R, '[:,;]'S, w, s ? page.width")
        if(!length(numbers) %in% 1:2 || !all(sapply(numbers, is_numeric_in_char))){

            stop_up("The argument 'page.width' must represent a 'Latex' measure, ex: 12cm. ",
                    "It can be optionally followed by a side margin width, as in '21, 2cm'. ",
                    "Problem: the current format could not be parsed ('", page.width, "').")
        }

        if(length(numbers) == 1){
            numbers[2] = numbers[1]
        } else {
            numbers = as.numeric(numbers)
            numbers[2] = numbers[1] - 2*numbers[2]
        }

        page.width = paste0(c(numbers[1], numbers[2]), unit)
    }

    attr(page.width, "no-margin") = TRUE

    page.width
}


is_Rmarkdown = function(){
    "knitr" %in% loadedNamespaces() && !is.null(knitr::pandoc_to())
}

path_to_relative = function(x){
    # orig = "C:/Users/berge028/Google Drive/R_packages/fixest/fixest"
    # dest = "C:/Users/berge028/Google Drive/R_packages/automake/automake/NAMESPACE"

    # I'm not sure it works perfectly well on linux...

    dest = normalizePath(x, "/", mustWork = FALSE)
    orig = normalizePath(".", "/", mustWork = FALSE)

    if(dest == orig) return(".")

    dest_split = strsplit(dest, "/")[[1]]
    orig_split = strsplit(orig, "/")[[1]]

    n_d = length(dest_split)
    n_o = length(orig_split)
    n_od = min(n_d, n_o)

    i_common = which.max(dest_split[1:n_od] != orig_split[1:n_od])

    # I'm not sure of that on linux....
    if(i_common == 1){
        if(dest_split[1] == orig_split[1]){
            # means all are the same => inclusion
            i_common = n_od

        } else if(dest_split[1] != orig_split[1]){
        # different roots? => can't do much
            return(dest)
        }
    } else {
        i_common = i_common - 1
    }

    # common = paste0(dest_split[1:i_common], collapse = "/")

    # From orig, we need to go up to the common dir
    if(i_common == n_o){
        up = "./"
    } else {
        n_up = n_o - i_common
        up = dsb("`n_up`*c!../")
    }

    if(i_common == n_d){
        down = ""
    } else {
        down = paste(dest_split[(i_common + 1):n_d], collapse = "/")
    }

    relpath = paste0(up, down)

    relpath
}












