#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Thu Jul 09 09:52:31 2020
# ~: etable
#----------------------------------------------#


#' Estimations table (export the results of multiples estimations to a DF or to Latex)
#'
#' Aggregates the results of multiple estimations and displays them in the form of either a Latex table or a \code{data.frame}.
#'
#' @inheritParams summary.fixest
#'
#' @param ... Used to capture different \code{fixest} estimation objects (obtained with \code{\link[fixest]{femlm}}, \code{\link[fixest]{feols}} or \code{\link[fixest]{feglm}}). Note that any other type of element is discarded. Note that you can give a list of \code{fixest} objects.
#' @param digits Integer, default is 4. The number of digits to be displayed.
#' @param tex Logical: whether the results should be a data.frame or a Latex table. By default, this argument is \code{TRUE} if the argument \code{file} (used for exportation) is not missing; it is equal to \code{FALSE} otherwise.
#' @param fitstat A character vector or a one sided formula. A vector listing which fit statistics to display. The valid types are 'll', 'aic', 'bic' and r2 types like 'r2', 'pr2', 'war2', etc (see all valid types in \code{\link[fixest]{r2}}). The default value depends on the models to display. Example of use: \code{fitstat=c('sq.cor', 'ar2', 'war2')}, or \code{fitstat=~sq.cor+ar2+war2} using a formula.
#' @param title (Tex only.) Character scalar. The title of the Latex table.
#' @param float (Tex only.) Logical. By default, if the argument \code{title} or \code{label} is provided, it is set to \code{TRUE}. Otherwise, it is set to \code{FALSE}.
#' @param sdBelow (Tex only.) Logical, default is \code{TRUE}. Should the standard-errors be displayed below the coefficients?
#' @param keep Character vector. This element is used to display only a subset of variables. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be kept. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{keep = "x[[:digit:]]$"}. If the first character is an exclamation mark, the effect is reversed (e.g. keep = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is kept). See details.
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). Each variable satisfying any of the regular expressions will be discarded. This argument is applied post aliasing (see argument \code{dict}). Example: you have the variable \code{x1} to \code{x55} and want to display only \code{x1} to \code{x9}, then you could use \code{drop = "x[[:digit:]]{2}"}. If the first character is an exclamation mark, the effect is reversed (e.g. drop = "!Intercept" means: every variable that does not contain \dQuote{Intercept} is dropped). See details.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a vector of regular expressions (see \code{\link[base]{regex}} help for more info). The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions. This argument is applied post aliasing (see argument \code{dict}). Example: you have the following variables: \code{month1} to \code{month6}, then \code{x1} to \code{x5}, then \code{year1} to \code{year6}. If you want to display first the x's, then the years, then the months you could use: \code{order = c("x", "year")}. If the first character is an exclamation mark, the effect is reversed (e.g. order = "!Intercept" means: every variable that does not contain \dQuote{Intercept} goes first).  See details.
#' @param dict A named character vector or a logical scalar. It changes the original variable names to the ones contained in the \code{dict}ionary. E.g. to change the variables named \code{a} and \code{b3} to (resp.) \dQuote{$log(a)$} and to \dQuote{$bonus^3$}, use \code{dict=c(a="$log(a)$",b3="$bonus^3$")}. By default, if Tex output is requested or if argument \code{file} is not missing, it is equal to \code{getFixest_dict()}, a default dictionary which can be set with \code{\link[fixest]{setFixest_dict}}. The default is not to change names if a \code{data.frame} is requested (i.e. \code{tex = FALSE}); if so, you can use \code{dict = TRUE} to use the dictionary you've set globally with \code{setFixest_dict()}.
#' @param file A character scalar. If provided, the Latex (or data frame) table will be saved in a file whose path is \code{file}. If you provide this argument, then a Latex table will be exported, to export a regular \code{data.frame}, use argument \code{tex = FALSE}.
#' @param replace Logical, default is \code{FALSE}. Only used if option \code{file} is used. Should the exported table be written in a new file that replaces any existing file?
#' @param convergence Logical, default is missing. Should the convergence state of the algorithm be displayed? By default, convergence information is displayed if at least one model did not converge.
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)} for a Latex table and \code{c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10)} for a data.frame (to conform with R's default). To supress the significance codes, use \code{signifCode=NA} or \code{signifCode=NULL}. Can also be equal to \code{"letters"}, then the default becomes \code{c("a"=0.01, "b"=0.05, "c"=0.10)}.
#' @param label (Tex only.) Character scalar. The label of the Latex table.
#' @param subtitles Character vector of the same length as the number of models to be displayed. If provided, subtitles are added underneath the dependent variable name.
#' @param fixef_sizes (Tex only.) Logical, default is \code{FALSE}. If \code{TRUE} and fixed-effects were used in the models, then the number of "individuals" per fixed-effect dimension is also displayed.
#' @param fixef_sizes.simplify Logical, default is \code{TRUE}. Only used if \code{fixef_sizes = TRUE}. If \code{TRUE}, the fixed-effects sizes will be displayed in parentheses instead of in a separate line if there is no ambiguity (i.e. if the size is constant across models).
#' @param yesNo (Tex only.) A character vector of length 1 or 2. Default is \code{c("Yes", "No")}. This is the message displayed when a given fixed-effect is (or is not) included in a regression. If \code{yesNo} is of length 1, then the second element is the empty string.
#' @param family Logical, default is missing. Whether to display the families of the models. By default this line is displayed when at least two models are from different families.
#' @param keepFactors Logical, default is \code{TRUE}. If \code{FALSE}, then factor variables are displayed as fixed-effects and no coefficient is shown.
#' @param powerBelow (Tex only.) Integer, default is -5. A coefficient whose value is below \code{10**(powerBelow+1)} is written with a power in Latex. For example \code{0.0000456} would be written \code{4.56$\\times 10^{-5}$} by default. Setting \code{powerBelow = -6} would lead to \code{0.00004} in Latex.
#' @param interaction.combine (Tex only.) Character scalar, defaults to \code{" $\\times$ "}. When the estimation contains interactions, then the variables names (after aliasing) are combined with this argument. For example: if \code{dict = c(x1="Wind", x2="Rain")} and you have the following interaction \code{x1:x2}, then it will be renamed (by default) \code{Wind $\\times$ Rain} -- using \code{interaction.combine = "*"} would lead to \code{Wind*Rain}.
#' @param depvar (Data frame only.) Logical, default is missing. Whether a first line containing the dependent variables should be shown. By default, the dependent variables are shown only if they differ across models or if the argumen \code{file} is not missing.
#' @param coefstat One of \code{"se"} (default), \code{"tstat"} or \code{"confint"}. The statistic to report for each coefficient: the standard-error, the t-statistics or the confidence interval. You can adjust the confidence interval with the argument \code{ci}.
#' @param ci Level of the confidence interval, defaults to \code{0.95}. Only used if \code{coefstat = confint}.
#' @param style A list. You can change the general style of the table with this argument. It should be of the form \code{style = list(keyword="key1:value1;key2:value2")} etc. The available keywords are \code{lines} (to manage the type of lines appearing in the table), and \code{depvar}, \code{model}, \code{var}, \code{fixef}, \code{slopes}, \code{fixef.sizes}, \code{stats} and \code{notes}. Most of these keywords accept the key \code{title:} which affects the title appearing just before the section. Eg to drop the \emph{Variables} header, just use \code{style=list("title:")}. Note that if you use \code{style=list("title: ")} (note the space after ":"), then an empty line will still be there. The keywords fixef, slopes and fixef.sizes also accept the keys \code{prefix} and \code{suffix}. E.g. if \code{style=list(fixef="suffix: FE")}, then there will be no header showing but the text " FE" will be appended to the ficed-effects variable names. The keys accepted in \code{lines} are \code{top}, \code{bottom}, \code{foot} and \code{sep}.
#' @param notes Character vector. If provided, a \code{"notes"} section will be added at the end right after the end of the table, containing the text of this argument. Note that if it is a vector, it will be collapsed with new lines.
#' @param group A list. The list elements should be vectors of regular expressions. For each elements of this list: A new line in the table is created, all variables that are matched by the regular expressions are discarded (same effect as the argument \code{drop}) and \code{TRUE} or \code{FALSE} will appear in the model cell, depending on whether some of the previous variables were found in the model. Example: \code{group=list("Controls: personal traits"=c("gender", "height", "weight"))} will create an new line with \code{"Controls: personal traits"} in the leftmost cell, all three variables gender, height and weight are discared, TRUE appearing in each model containing at least one of the three variables (the style of TRUE/FALSE is governed by the argument \code{yesNo}). You can control the style with the \code{title} and \code{where} keywords in curly brackets. For example \code{group=list("{title:Controls; where:stats}Personal traits"=c("gender", "height", "weight"))} will add an extra line right before with "Control" written in it, and the group information will appear after the statistics. The keyword where can be equal to either \code{var} (default), \code{fixef} or \code{stats}.
#' @param extraline A list. The list elements should be either a single logical or a vector of the same length as the number of models. For each elements of this list: A new line in the table is created, the list name being the row name and the vector being the content of the cells. Example: \code{extraline=list("Sub-sample"=c("<20 yo", "all", ">50 yo"))} will create an new line with \code{"Sub-sample"} in the leftmost cell, the vector filling the content of the cells for the three models. You can control the style with the \code{title} and \code{where} keywords in curly brackets. For example \code{group=list("{title:Sub-sample; where:stats}By age"=c("<20 yo", "all", ">50 yo"))} will add an extra line right before with "Sub-sample" written in it, and the extraline information will appear after the statistics section. The keyword where can be equal to either \code{var} (default), \code{fixef} or \code{stats}.
#' @param tablefoot Logical, default is \code{TRUE}. Whether to display the table footer containing the information on the way the standard-errors where computed and the meaning of the significance codes.
#' @param reset (\code{setFixest_etable} only.) Logical, default is \code{FALSE}. If \code{TRUE}, this will reset all the default values that were already set by the user in previous calls.
#'
#' @details
#' The function \code{esttex} is equivalent to the function \code{etable} with argument \code{tex = TRUE}.
#'
#' The function \code{esttable} is equivalent to the function \code{etable} with argument \code{tex = FALSE}.
#'
#' You can permanently change the way your table looks in Latex by using \code{setFixest_etable}. The following vignette gives and example as well as illustrates how to use the argument \code{style}: \href{https://cran.r-project.org/package=fixest/vignettes/exporting_tables.html}{Exporting estimation tables}.
#'
#' @section Arguments keep, drop and order:
#' The arguments \code{keep}, \code{drop} and \code{order} use regular expressions. If you are not aware of regular expressions, I urge you to learn it, since it is an extremely powerful way to manipulate character strings (and it exists across most programming languages).
#'
#' For example drop = "Wind" would drop any variable whose name contains "Wind". Note that variables such as "Temp:Wind" or "StrongWind" do contain "Wind", so would be dropped. To drop only the variable named "Wind", you need to use \code{drop = "^Wind$"} (with "^" meaning beginning, resp. "$" meaning end, of the string => this is the language of regular expressions).
#'
#' Although you can combine several regular expressions in a single character string using pipes, \code{drop} also accepts a vector of regular expressions.
#'
#' You can use the special character "!" (exclamation mark) to reverse the effect of the regular expression (this feature is specific to this fonction). For example \code{drop = "!Wind"} would drop any variable that does not contain "Wind".
#'
#' You can use the special character "%" (percentage) to make reference to the original variable name instead of the aliased name. For example, you have a variable named \code{"Month6"}, and use a dictionary \code{dict = c(Month6="June")}. Thus the variable will be displayed as \code{"June"}. If you want to delete that variable, you can use either \code{drop="June"}, or \code{drop="%Month6"} (which makes reference to its original name).
#'
#' The argument \code{order} takes in a vector of regular expressions, the order will follow the elments of this vector. The vector gives a list of priorities, on the left the elements with highest priority. For example, order = c("Wind", "!Inter", "!Temp") would give highest priorities to the variables containing "Wind" (which would then appear first), second highest priority is the variables not containing "Inter", last, with lowest priority, the variables not containing "Temp". If you had the following variables: (Intercept), Temp:Wind, Wind, Temp you would end up with the following order: Wind, Temp:Wind, Temp, (Intercept).
#'
#'
#'
#' @return
#' If \code{tex = TRUE}, the lines composing the Latex table are returned invisibly while the table is directly prompted on the console.
#'
#' If \code{tex = FALSE}, the data.frame is directly returned. If the argument \code{file} is not missing, the \code{data.frame} is returned invisibly.
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
#' aq$Month = factor(aq$Month)
#'
#' est1 = feols(Ozone ~ Month / Wind + Temp, data = aq)
#' est2 = feols(Ozone ~ Wind + Temp | Month, data = aq)
#'
#' # Displaying the two results in a single table
#' etable(est1, est2)
#'
#' # keep/drop: keeping only interactions
#' etable(est1, est2, keep = ":")
#' # or using drop  (see regexp help):
#' etable(est1, est2, drop = "^[[:alnum:]]+$")
#'
#' # keep/drop: dropping interactions
#' etable(est1, est2, drop = ":")
#' # or using keep ("!" reverses the effect):
#' etable(est1, est2, keep = "!:")
#'
#' # order: Wind variable first, intercept last
#' etable(est1, est2, order = c("Wind", "Month"))
#' etable(est1, est2, order = c("^Wind", "Wind", "Month"))
#' # Interactions, then Intercept, last ("!" reverses the effect)
#' etable(est1, est2, order = c("!Int", "!:"))
#'
#' # dict + keep/drop/order: using "%" to match the original names
#' dict = c("Month5"="May", "Month6"="Jun", "Month7"="Jul",
#'          "Month8"="Aug", "Month9"="Sep")
#' etable(est1, est2, tex = TRUE, dict = dict)
#' # keeping only June and July
#' etable(est1, est2, tex = TRUE, dict = dict, keep = c("%Month6", "Jul"))
#' # All months variabes first
#' etable(est1, est2, tex = TRUE, dict = dict, order = c("%Month"))
#'
#' # signifCode
#' etable(est1, est2, signifCode = c(" A"=0.01, " B"=0.05, " C"=0.1,
#'                                   " D"=0.15, " F"=1))
#'
#' # fitstat
#' etable(est1, est2, fitstat = ~r2+ar2+apr2+war2)
#'
#' # Adding a dictionnary (Tex only)
#' dict = c(Month5="May", Month6="Jun", Month7="Jul", Month8="Aug", Month9="Sep")
#' etable(est1, est2, dict = dict, tex = TRUE)
#'
#' #
#' # Using the argument style to customize Latex exports
#' #
#'
#' # If you don't like the default layout of the table, no worries!
#' # You can modify many parameters with the argument style
#'
#' # To drop the headers before each section, use:
#' style_noHeaders = list(var="", fixef="suffix: FE", stats = "")
#' etable(est1, est2, dict = dict, tex = TRUE, style = style_noHeaders)
#'
#' # To change the lines of the table
#' style_lines = list(lines = "top:\\toprule;bottom:\\bottomrule;sep:\\midrule;foot:\\midrule")
#' etable(est1, est2, dict = dict, tex = TRUE, style = style_lines)
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
#' est_c2 = feols(Ozone ~ Solar.R + I(Solar.R**2) + ..ctrl, data = aq)
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
#' # When exporting to Latex, you can add meta arguments to 'group' and 'extraline'
#' # Two keywords are allowed: 'title' and 'where'
#' # 'title' adds a line just before with the content of 'title' in the leftmost cell
#' # 'where' governs the location of the line. It can be equal to 'var', 'stats' or 'fixef'.
#' # The syntax is: {title:Controls; where:stats}Group name.
#' # Note that starting with curly braces is mandatory.
#'
#' # Examples
#' etable(est_c0, est_c1, est_c2, tex = TRUE, group = list("{where:stats}Controls" = "poly"))
#' etable(est_all, est_sub1, est_sub2, est_sub3, tex = TRUE,
#'        extraline = list("{title:\\midrule}Sub-sample" =
#'                           c("All", "May-June", "Jul.-Aug.", "Sept.")))
#'
etable = function(..., se = c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = getFixest_dof(), cluster, digits=4, tex, fitstat, title, coefstat = c("se", "tstat", "confint"), ci = 0.95, sdBelow = TRUE, keep, drop, order, dict, file, replace=FALSE, convergence, signifCode, label, float, subtitles, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, yesNo = "Yes", keepFactors = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ ", depvar, style = list(), notes = NULL, group = NULL, extraline = NULL, tablefoot = TRUE){

    #
    # Checking the arguments
    #

    # Need to check for the presence of the se
    useSummary = TRUE
    if(missing(se) && missing(cluster)){
        useSummary = FALSE
    }

    if(!missing(se)){
        check_arg_plus(se, "match")
    } else {
        se = NULL
    }

    check_arg_plus(coefstat, "match")

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


    # to get the model names
    dots_call = match.call(expand.dots = FALSE)[["..."]]

    info = results2formattedList(..., se=se, dof=dof, fitstat=fitstat, cluster=cluster, digits=digits, sdBelow=sdBelow, signifCode=signifCode, coefstat = coefstat, ci = ci, title=title, float=float, subtitles=subtitles, yesNo=yesNo, keepFactors=keepFactors, tex = tex, useSummary=useSummary, dots_call=dots_call, powerBelow=powerBelow, dict=dict, interaction.combine=interaction.combine, convergence=convergence, family=family, keep=keep, drop=drop, file=file, order=order, label=label, fixef_sizes=fixef_sizes, fixef_sizes.simplify=fixef_sizes.simplify, depvar=depvar, style=style, replace=replace, notes = notes, group = group, tablefoot = tablefoot, extraline=extraline)

    if(tex){
        res = etable_internal_latex(info)
    } else {
        res = etable_internal_df(info)
    }

    if(!missnull(file)){
        sink(file = file, append = !replace)
        on.exit(sink())

        if(tex){
            cat(res, sep = "")
        } else {
            print(res)
        }

        return(invisible(res))
    } else {
        if(tex){
            cat(res, sep = "")
            return(invisible(res))
        } else {
            return(res)
        }
    }
}


#' @describeIn etable Exports the results of multiple \code{fixest} estimations in a Latex table.
esttex <- function(..., se = c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = getFixest_dof(), cluster, digits=4, fitstat, coefstat = c("se", "tstat", "confint"), ci = 0.95, title, float = float, sdBelow=TRUE, keep, drop, order, dict, file, replace=FALSE, convergence, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, subtitles, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, yesNo = "Yes", keepFactors = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ ", style = list(), notes = NULL, group = NULL, tablefoot = TRUE, extraline = NULL){
    # drop: a vector of regular expressions
    # order: a vector of regular expressions
    # dict: a 'named' vector
    # file: a character string

    useSummary = TRUE
    if(missing(se) && missing(cluster)){
        useSummary = FALSE
    }

    if(!missing(se)){
        check_arg_plus(se, "match")
    } else {
        se = NULL
    }

    check_arg_plus(coefstat, "match")

    # Float or not
    check_arg(float, "logical scalar")
    if(missing(float)){
        if(!missing(title) || !missing(label)){
            float = TRUE
        } else {
            float = FALSE
        }
    }

    # to get the model names
    dots_call = match.call(expand.dots = FALSE)[["..."]]

    info = results2formattedList(..., tex = TRUE, useSummary=useSummary, se = se, dof = dof, cluster = cluster, digits = digits, fitstat = fitstat, title = title, float = float, sdBelow = sdBelow, keep=keep, drop = drop, order = order, dict = dict, file = file, replace = replace, convergence = convergence, signifCode = signifCode, coefstat = coefstat, ci = ci, label = label, subtitles = subtitles, fixef_sizes = fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, yesNo = yesNo, keepFactors = keepFactors, family = family, powerBelow = powerBelow, interaction.combine = interaction.combine, dots_call = dots_call, depvar = TRUE, style=style, notes = notes, group = group, tablefoot = tablefoot, extraline=extraline)

    res = etable_internal_latex(info)

    if(!missnull(file)){
        sink(file = file, append = !replace)
        on.exit(sink())
    }

    cat(res, sep = "")
    return(invisible(res))

}

#' @describeIn etable Facility to display the results of multiple \code{fixest} estimations.
esttable <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), dof = getFixest_dof(), cluster, coefstat = c("se", "tstat", "confint"), ci = 0.95, depvar, keep, drop, dict, order, digits=4, fitstat, convergence, signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.10), subtitles, keepFactors = FALSE, family, group = NULL, extraline = NULL){

    # Need to check for the presence of the se
    useSummary = TRUE
    if(missing(se) && missing(cluster)){
        useSummary = FALSE
    }

    if(!missing(se)){
        se = match.arg(se)
    } else {
        se = NULL
    }

    coefstat = match.arg(coefstat)

    # to get the model names
    dots_call = match.call(expand.dots = FALSE)[["..."]]

    info = results2formattedList(..., se=se, dof = dof, cluster=cluster, digits=digits, signifCode=signifCode, coefstat = coefstat, ci = ci, subtitles=subtitles, keepFactors=keepFactors, useSummary=useSummary, dots_call=dots_call, fitstat=fitstat, yesNo = c("Yes", "No"), depvar = depvar, family = family, keep = keep, drop = drop, order = order, dict = dict, interaction.combine = ":", group = group, extraline = extraline)

    res = etable_internal_df(info)

    return(res)
}

results2formattedList = function(..., se, dof = getFixest_dof(), cluster, digits = 4, fitstat, sdBelow=TRUE, dict, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), coefstat = "se", ci = 0.95, label, subtitles, title, float = FALSE, replace = FALSE, yesNo = c("Yes", "No"), keepFactors = FALSE, tex = FALSE, useSummary, dots_call, powerBelow, interaction.combine, convergence, family, drop, order, keep, file, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, depvar = FALSE, style = list(), notes = NULL, group = NULL, tablefoot = TRUE, extraline=NULL){
    # This function is the core of the functions esttable and esttex


    # Setting the default values (we take extra care for "style")
    if(tex){
        check_arg_plus(style, "NULL{list()} named list l0")
        all_sections = strsplit("depvar, model, var, fixef, fixef.sizes, slopes, stats, lines, notes", ", ")[[1]]
        check_value(names(style), "NULL multi charin", .message = paste0("The names of the argument 'style' must be one of ", enumerate_items(all_sections, "or.quote", nmax = 100), "."), .choices = all_sections)

        # The variable style will be changed via the defaults
        style_user = style
    }

    #
    # Setting the default
    #

    opts = getOption("fixest_etable")
    if(length(opts) > 0){
        sysOrigin = sys.parent()
        mc = match.call(definition = sys.function(sysOrigin), call = sys.call(sysOrigin), expand.dots = FALSE)
        args_in = setdiff(names(mc), "style")

        # Arguments for which the defaults should not be changed in etable, tex = FALSE
        if(!tex) args_in = c(args_in, "signifCode")

        # We modify only non-user provided arguments
        for(v in names(opts)){
            if(!v %in% args_in){
                assign(v, opts[[v]])
            }
        }
    }


    #
    # Full control
    #

    set_up(1)
    # => we also allow lists (of fixest objects)
    check_arg(..., "class(fixest) | list mbt")

    check_arg(digits, "integer scalar GE{1}")
    check_arg(title, "character scalar")
    check_arg_plus(coefstat, "match(se, tstat, confint)")

    check_arg_plus(notes, "NULL{''} character vector no na")
    if(length(notes) > 1) notes = paste(notes, collapse = "\n")

    check_arg("logical scalar", sdBelow, replace, convergence, fixef_sizes, fixef_sizes.simplify, keepFactors, family, tex, depvar)
    check_arg("logical scalar", tablefoot)
    isTex = tex
    if(missing(family)){
        show_family = NULL
    } else {
        show_family = family
    }

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
    check_arg(subtitles, "character vector no na")
    check_arg(yesNo, "character vector no na len(,2)")

    if(isTex == FALSE) yesNo = c("Yes", "No")

    check_arg(powerBelow, "integer scalar LE{-1}")

    check_arg(ci, "numeric scalar GT{0.5} LT{1}")

    check_arg(dict, "NULL logical scalar | named character vector no na")

    check_arg_plus(group, extraline, "NULL{list()} named list l0")
    # we check it more in depth later


    # Setting the style defaults
    if(tex){
        # We check the argument style + set it // all checks are made here

        # We modify the default set with setFixest_etable()
        if(length(style_user) > 0){
            style[names(style_user)] = style_user
        }

        # setting the default + parsing
        check_value_plus(style$depvar, "NULL{'title:Dependent Variable(s):'} character scalar")
        style$depvar = parse_style(style$depvar, "title")

        check_value_plus(style$model, "NULL{'title:Model:;format:(1)'} character scalar")
        style$model = parse_style(style$model, c("title", "format"))

        default = 'top:\\tabularnewline\\toprule\\toprule;bottom:;foot:\\bottomrule\\bottomrule;sep:\\midrule'
        check_value_plus(style$lines, "NULL{default} character scalar")
        style$lines = parse_style(style$lines, c("top", "bottom", "sep", "foot"))

        # problem parsing { in dreamerr
        default = 'title:\\emph{Variables}'
        check_value_plus(style$var, "NULL{default} character scalar")
        style$var = parse_style(style$var, "title")

        default = 'title:\\emph{Fixed-effects}'
        check_value_plus(style$fixef, "NULL{default} character scalar")
        style$fixef = parse_style(style$fixef, c("title", "prefix", "suffix", "where"))

        default = 'title:\\emph{Varying Slopes}'
        check_value_plus(style$slopes, "NULL{default} character scalar")
        style$slopes = parse_style(style$slopes, c("title", "format"))

        check_value_plus(style$fixef.sizes, "NULL{'prefix:# ;where:obs'} character scalar")
        style$fixef.sizes = parse_style(style$fixef.sizes, c("title", "prefix", "suffix", "where"))

        default = "title:\\emph{Fit statistics}"
        check_value_plus(style$stats, "NULL{default} character scalar")
        style$stats = parse_style(style$stats, "title")

        default = "title:\\emph{\\medskip Notes:} "
        check_value_plus(style$notes, "NULL{default} character scalar")
        style$notes = parse_style(style$notes, "title")
    }



    # default values for dict
    if(missing(dict)){
        if(isTex || !missing(file)){
            dict = getFixest_dict()
        } else {
            dict = NULL
        }
    } else if(isTRUE(dict)) {
        dict = getFixest_dict()
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

    if(length(yesNo) == 1){
        yesNo = c(yesNo, "")
    }

    # at the moment: only fixest allowed
    allowed_types = "fixest"

    # We get all the models
    dots <- list(...)

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
    k = 1
    for(i in 1:n){
        di = dots[[i]]

        if(any(allowed_types %in% class(di))){
            all_models[[k]] = di
            if(any(class(dots_call[[i]]) %in% c("call", "name"))){
                model_names[[k]] = deparse_long(dots_call[[i]])
            } else {
                model_names[[k]] = as.character(dots_call[[i]])
            }

            k = k+1
        } else if(length(class(di))==1 && class(di)=="list"){
            # we get into this list to get the fixest objects
            types = sapply(di, class)
            qui = which(types %in% allowed_types)
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
                }

                k = k+1
            }
        }

    }

    if(length(all_models)==0) stop_up("Not any proper model (fixest) as argument!")

    n_models <- length(all_models)


    # We check the group and extraline arguments
    if(missing(drop)) drop = NULL
    for(i in seq_along(group)){
        check_value(group[[i]], "character vector", .message = "The elements of argument 'group' must be character vectors of regular expressions.")
        drop = unique(c(drop, group[[i]]))
    }

    for(i in seq_along(extraline)){
        check_value(extraline[[i]], "logical scalar | vector(character, numeric, logical) len(value)", .message = paste0("The elements of argument 'extraline' must be vectors of length ", n_models, " or logical scalars."), .value = n_models)
    }

    # formatting the names (if needed)
    alternative_names = paste0("model ", 1:n_models)
    who2replace = sapply(model_names, function(x) length(x) == 0 || x == "")
    model_names[who2replace] = alternative_names[who2replace]

    # we keep track of the SEs
    se_type_list = list()

    check_interaction_reorder = FALSE
    var_list <- var_reorder_list <- coef_list <- coef_below <- sd_below <- list()
    depvar_list <- obs_list <- fitstat_list <- list()
    r2_list <- aic_list <- bic_list <- loglik_list <- convergence_list <- list()
    sqCor_list = family_list = theta_list = list()

    # To take care of factors
    fe_names = c()
    is_fe = vector(mode = "list", n_models)
    nb_fe = vector(mode = "list", n_models) # the number of items per factor

    slope_names = c()
    slope_flag_list = vector(mode = "list", n_models)

    # if there are subtitles
    if(!missing(subtitles)){
        if(length(subtitles) != n_models){
            stop_up("If argument 'subtitles' is provided, it must be of the same length as the number of models. Current lengths: ", length(subtitles), " vs ", n_models, " models.")
        } else {
            isSubtitles = TRUE
        }

        if(isTex){
            subtitles = escape_latex(subtitles, up = 2)
        }

    } else {
        subtitles = NULL
        isSubtitles = FALSE
    }

    if(!is.null(dict) && isTex){
        dict = escape_latex(dict, up = 2)
    }

    #
    # fitstat: which R2 to display?
    #

    if(missing(fitstat)){
        # Default values:
        #   - if all OLS: typical R2
        #   - if any non-OLS: pseudo R2 + squared cor.
        is_ols = sapply(all_models, function(x) deparse(x$call[[1]]) == "feols")

        if(all(is_ols)){
            if(any(sapply(all_models, function(x) "fixef_vars" %in% names(x)))){
                # means any FE model
                fitstat = c("r2", "wr2")
            } else {
                fitstat = c("r2", "ar2")
            }
        } else {
            fitstat = c("sq.cor", "pr2", "bic")
        }


    } else if(isFALSE(fitstat) || (length(fitstat) == 1 && (is.na(fitstat) || fitstat == ""))){
        fitstat = NULL
    } else if("formula" %in% class(fitstat)){
        check_arg(fitstat, "os formula", .message = "Argument 'fitstat' must be a one sided formula (or a character vector) containing 'aic', 'bic', 'll', or valid r2 types names (see function r2). ")
        fitstat = attr(terms(fitstat), "term.labels")
    } else {
        check_arg(fitstat, "character vector no na", .message = "Argument 'fitstat' must be a character vector (or a one sided formula) containing 'aic', 'bic', 'll', or valid r2 types names (see function r2). ")
    }

    # checking the types
    fitstat_type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "par2", "wr2", "war2", "awr2", "wpr2", "pwr2", "wapr2", "wpar2", "awpr2", "apwr2", "pawr2", "pwar2", "ll", "aic", "bic")
    fitstat = unique(fitstat)

    pblm = setdiff(fitstat, fitstat_type_allowed)
    if(length(pblm) > 0){
        stop_up("Argument 'fitstat' must be a character vector (or a one sided formula) containing 'aic', 'bic', 'll', or valid r2 types names. ", enumerate_items(pblm, "is.quote"), " not valid (see function r2).")
    }

    fitstat_dict_tex = c("sq.cor"="Squared Correlation", r2="R$^2$", ar2="Adjusted R$^2$", pr2="Pseudo R$^2$", apr2="Adjusted Pseudo R$^2$", wr2="Within R$^2$", war2="Within Adjusted R$^2$", wpr2="Within Pseudo R$^2$", wapr2="Whithin Adjusted Pseudo R$^2$", aic = "AIC", bic = "BIC", ll = "Log-Likelihood")

    fitstat_dict_R = c("sq.cor"="Squared Corr.", r2="R2", ar2="Adjusted R2", pr2="Pseudo R2", apr2="Adj. Pseudo R2", wr2="Within R2", war2="Within Adj. R2", wpr2="Within Pseudo R2", wapr2="Whithin Adj. Pseudo R2", aic = "AIC", bic = "BIC", ll = "Log-Likelihood")

    fitstat_dict = fitstat_dict_R
    if(isTex) fitstat_dict = fitstat_dict_tex

    # end: fitstat

    for(m in 1:n_models){

        # If se or cluster is provided, we use summary
        if(useSummary){
            x = summary(all_models[[m]], se=se, cluster, dof = dof, nframes_up = 2)
        } else {
            # What do we do if se not provided?
            # we apply summary only to the ones that are not summaries
            x = all_models[[m]]
            if(!"cov.scaled" %in% names(x)){
                # not a summary => we apply summary to trigger default behavior
                x = summary(x, dof = dof)
            }

        }
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


        # Negbin parameter
        theta = all_models[[m]]$theta
        theta_list[[m]] = ifelse(is.null(theta), "", numberFormatNormal(theta))

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

        # now the normal FEs
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

            if(tex && nchar(style$slopes$format) > 0){
                slope_format = style$slopes$format
                slope_var_full = c()
                for(i in seq_along(slope_vars_name)){
                    slope_var_full[i] = gsub("\\_\\_slope\\_\\_", slope_fe_name[i], gsub("\\_\\_var\\_\\_", slope_vars_name[i], slope_format, fixed = TRUE), fixed = TRUE)
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
            if(length(var_left) > 0 && any(grepl(":", var_left))){
                check_interaction_reorder = TRUE


                qui_inter = grepl(":", var_left)
                inter = strsplit(var_left[qui_inter], "(?<=[^:]):(?=[^:])", perl = TRUE)

                fun_rename = function(x){
                    # We put the factors on the right
                    qui_factor = grepl("::", x)
                    if(any(qui_factor)){
                        res = x[base::order(qui_factor)]
                        res = gsub("::", " $=$ ", res)
                    } else {
                        res = x
                    }

                    who = res %in% names(dict)

                    res[who] = dict[res[who]]
                    paste0(res, collapse = interaction.combine)
                }

                inter_named = sapply(inter, fun_rename)
                new_inter = sapply(inter, function(x) fun_rename(sort(x)))

                var[!qui][qui_inter] = inter_named
                new_var[!qui][qui_inter] = new_inter

            }
            names(new_var) = names(var) = var_origin
            var_reorder_list[[m]] <- new_var
        } else {
            # We reorder the interaction terms alphabetically
            new_var = var
            qui = grepl(":", new_var)
            if(any(qui)){
                check_interaction_reorder = TRUE
                inter = strsplit(new_var[qui], ":")
                new_inter = sapply(inter, function(x) paste0(sort(x), collapse = ":"))
                new_var[qui] = new_inter
            }

            names(new_var) = names(var) = var_origin

            var_reorder_list[[m]] <- new_var
        }

        if(isTex){
            fun_format = function(x) coefFormatLatex(x, digits = digits, power = abs(powerBelow))
        } else {
            fun_format = function(x) as.character(mysignif(x, d = digits))
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

        # if(isTex){
        #     coef = coefFormatLatex(a[, 1], digits = digits, power = abs(powerBelow))
        #     se_value = coefFormatLatex(a[, 2], digits = digits, power = abs(powerBelow))
        # } else {
        #     coef = as.character(round(a[, 1], digits))
        #     se_value = as.character(myRound(a[, 2], digits))
        # }

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


        # If the coefficient is bounded, we supress the 'stars'
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

        # Pseudo-R2 // AIC // BIC // N
        n <- nobs(x)
        obs_list[[m]] <- n
        convergence_list[[m]] = ifelse(is.null(x$convStatus), TRUE, x$convStatus)

        K <- x$nparams

        if(length(fitstat) == 0){
            fitstat_list[[m]] = NA
        } else {
            fistat_format = list()

            fun_format = ifelse(isTex, numberFormatLatex, numberFormatNormal)

            if("aic" %in% fitstat) fistat_format[["aic"]] = fun_format(round(AIC(x), 3))
            if("bic" %in% fitstat) fistat_format[["bic"]] = fun_format(round(BIC(x), 3))
            if("ll" %in% fitstat) fistat_format[["ll"]] = fun_format(logLik(x))

            # regular r2s
            r2_type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "wr2", "war2", "wpr2", "wapr2")
            if(any(fitstat %in% r2_type_allowed)){
                all_r2 = r2(x, intersect(fitstat, r2_type_allowed))
                for(r2_val in names(all_r2)){
                    nb = 5 - (r2_val == "sq.cor") * 2
                    fistat_format[[r2_val]] = round(as.vector(all_r2[r2_val]), nb)
                }
            }

            fitstat_list[[m]] = fistat_format[fitstat]
        }

    }

    if(check_interaction_reorder){
        if(length(unique(unlist(var_reorder_list))) < length(unique(unlist(var_list)))){
            var_list = var_reorder_list
            for(m in 1:length(var_list)){
                names(coef_list[[m]]) <- var_list[[m]]
            }
        }
    }


    if(length(fitstat) > 0){
        attr(fitstat_list, "format_names") = fitstat_dict[fitstat]
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

    res = list(se_type_list=se_type_list, var_list=var_list, coef_list=coef_list, coef_below=coef_below, sd_below=sd_below, depvar_list=depvar_list, obs_list=obs_list, convergence_list=convergence_list, fe_names=fe_names, is_fe=is_fe, nb_fe=nb_fe, slope_flag_list = slope_flag_list, slope_names=slope_names, useSummary=useSummary, model_names=model_names, family_list=family_list, theta_list=theta_list, fitstat_list=fitstat_list, subtitles=subtitles, isSubtitles=isSubtitles, title=title, convergence=convergence, family=family, keep=keep, drop=drop, order=order, file=file, label=label, sdBelow=sdBelow, signifCode=signifCode, fixef_sizes=fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, depvar=depvar, useSummary=useSummary, dict=dict, yesNo=yesNo, add_signif=add_signif, float=float, coefstat=coefstat, ci=ci, style=style, notes=notes, group=group, tablefoot=tablefoot, extraline=extraline)

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
    obs_list = info$obs_list
    convergence_list = info$convergence_list
    fe_names = info$fe_names
    is_fe = info$is_fe
    nb_fe = info$nb_fe
    slope_names = info$slope_names
    slope_flag_list = info$slope_flag_list
    family_list = info$family_list
    theta_list = info$theta_list
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
    tablefoot = info$tablefoot
    extraline = info$extraline

    # Formatting the searating lines
    if(nchar(style$lines$top) > 1) style$lines$top = paste0(style$lines$top, "\n")
    if(nchar(style$lines$bottom) > 1) style$lines$bottom = paste0(style$lines$bottom, "\n")
    if(nchar(style$lines$sep) > 1) style$lines$sep = paste0(style$lines$sep, "\n")
    if(nchar(style$lines$foot) > 1) style$lines$foot = paste0(style$lines$foot, "\n")

    #
    # prompting the infos gathered
    #

    # Starting the table
    myTitle = title
    if(!is.null(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
    if(float){
        start_table = paste0("\\begin{table}[htbp]\\centering\n\\caption{",  myTitle, "}\n")
        end_table = "\\end{table}"
    } else {
        start_table = ""
        end_table = ""
    }


    # intro and outro Latex tabular
    intro_latex <- paste0("\\begin{tabular}{l", paste0(rep("c", n_models), collapse=""), "}\n", style$lines$top)

    outro_latex <- "\\end{tabular}\n"

    # 1st lines => dep vars
    depvar_list = dict_apply(c(depvar_list, recursive = TRUE), dict)
    depvar_list = escape_latex(depvar_list, up = 2)

    # We write the dependent variables properly, with multicolumn when necessary
    # to do that, we count the number of occurences of each variable (& we respect the order provided by the user)
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
    first_line <- style$depvar$title
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

    # Model line
    if(nchar(style$model$format) > 0){
        m = style$model$format
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
    model_line = paste0(style$model$title, "&", paste0(model_format, collapse = " & "), "\\\\\n")

    # a simple line with only "variables" written in the first cell
    if(nchar(style$var$title) == 0){
        variable_line = style$lines$sep
    } else if(style$var$title == "\\midrule"){
        variable_line = "\\midrule "
    } else {
        # variable_line = paste0(style$lines$sep, style$var$title, "\\tabularnewline\n")
        variable_line = paste0(style$lines$sep, style$var$title, "& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
    }

    # Coefficients,  the tricky part
    coef_lines <- list()

    # we need to loop not to lose names
    all_vars = c()
    for(vars in var_list){
        all_vars = c(all_vars, vars[!vars %in% all_vars])
    }
    # all_vars <- unique(unlist(var_list))

    for(i in seq_along(group)){
        gi = group[[i]]
        present = rep(FALSE, n_models)
        for(v in gi){
            if(grepl("^%", v)){
                # original names
                v = substr(v, 2, nchar(v))
                present = present | sapply(var_list, function(x) any(grepl(v, names(x))))
            } else {
                present = present | sapply(var_list, function(x) any(grepl(v, x)))
            }
        }

        group[[i]] = present
    }

    # keeping some coefs
    all_vars = keep_apply(all_vars, keep)

    # dropping some coefs
    all_vars = drop_apply(all_vars, drop)

    # ordering the coefs
    all_vars = order_apply(all_vars, order)

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
            myCoef = mySd= myLine = c()
            for(m in 1:n_models){
                myCoef = c(myCoef, coef_below[[m]][v])
                mySd = c(mySd, sd_below[[m]][v])
            }

            myCoef[is.na(myCoef)] = "  "
            mySd[is.na(mySd)] = "  "
            myCoef = paste0(aliasVars[v], "&", paste0(myCoef, collapse="&"))
            mySd = paste0("  &", paste0(mySd, collapse="&"))
            myLines = paste0(myCoef, "\\\\\n", mySd, "\\\\\n")
            coef_lines = c(coef_lines, myLines)
        }
        coef_lines = paste0(coef_lines, collapse="")
    } else {
        coef_lines = paste0(paste0(apply(coef_mat, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")
    }

    # Fixed-effects (if needed)
    nb_FE_lines = ""
    if(length(fe_names) > 0){

        if(nchar(style$fixef$title) == 0){
            dumIntro = style$lines$sep
        } else if(style$fixef$title == "\\midrule"){
            dumIntro = "\\midrule "
        } else {
            dumIntro = paste0(style$lines$sep, style$fixef$title, "& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
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
        for(i in seq_along(fe_names)){
            fe = fe_names[i]

            if(fe %in% names(dict)){
                fe_names[i] = dict[fe]
            } else if(grepl("\\^", fe)){
                fe_split = strsplit(fe, "\\^")[[1]]
                who = fe_split %in% names(dict)
                fe_split[who] = dict[fe_split[who]]
                fe_names[i] = paste(fe_split, collapse = "$\\times$")
            } else {
                fe_names[i] = fe
            }
        }
        fe_names_raw = escape_latex(fe_names)

        fe_names = paste0(style$fixef$prefix, fe_names_raw, style$fixef$suffix)

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
                fe_names_nbItems = paste0(style$fixef.sizes$prefix, fe_names_raw[is_complex], style$fixef.sizes$suffix)
                all_nb_FEs = cbind(fe_names_nbItems, all_nb_FEs[is_complex, , drop = FALSE])
                nb_FE_lines <- paste0(paste0(apply(all_nb_FEs, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")
            }

        }

        all_fe = cbind(fe_names, all_fe)
        FE_lines <- paste0(paste0(apply(all_fe, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

    } else {
        FE_lines = NULL
        dumIntro = NULL
    }

    # Slopes (if needed)
    if(length(slope_names) > 0){

        if(nchar(style$slopes$title) == 0){
            slope_intro = style$lines$sep
        } else if(style$slopes$title == "\\midrule"){
            slope_intro = "\\midrule "
        } else {
            slope_intro = paste0(style$lines$sep, style$slopes$title, "& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
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
        slope_lines <- paste0(paste0(apply(all_slopes, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

    } else {
        slope_intro = NULL
        slope_lines = NULL
    }

    # Subtitles
    if(isSubtitles){
        info_subtitles = paste0("  & ", paste(subtitles, collapse="&"), "\\\\\n")
    } else {
        info_subtitles = ""
    }

    # Convergence information
    info_convergence = ""
    if(convergence){
        info_convergence = paste0("Convergence &", paste(convergence_list, collapse="&"), "\\\\\n")
    }

    info_theta <- paste0("Overdispersion& ", paste(theta_list, collapse="&"), "\\\\\n")

    # information on family
    if(family){
        info_family <- paste0("Family &  ", paste(family_list, collapse=" & "), "\\\\\n")
    } else {
        info_family = ""
    }


    # The standard errors => if tablefoot = TRUE
    if(tablefoot){
        isUniqueSD = length(unique(unlist(se_type_list))) == 1
        nb_col = length(obs_list) + 1
        sd_intro = paste0("\\multicolumn{", nb_col, "}{l}{\\emph{")
        if(isUniqueSD){
            my_se = unique(unlist(se_type_list)) # it comes from summary
            # every model has the same type of SE
            if(my_se == "Standard") my_se = "Normal"
            if(my_se == "White") my_se = "White-corrected"

            # Now we modify the names of the clusters if needed
            my_se = format_se_type_latex(my_se, dict)

            if(coefstat == "se"){
                coefstat_sentence = " standard-errors in parentheses"
            } else if(coefstat == "tstat"){
                coefstat_sentence = " co-variance matrix, t-stats in parentheses"
            } else {
                coefstat_sentence = paste0(" co-variance matrix, ", round(ci*100), "\\% confidence intervals in brackets")
            }
            info_SD = paste0(style$lines$foot, sd_intro, my_se, coefstat_sentence, "}}\\\\\n")

            if(add_signif){
                info_SD = paste0(info_SD, sd_intro, "Signif. Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
            }

            info_muli_se = ""
        } else {
            all_se_type = sapply(se_type_list, format_se_type_latex, dict = dict, inline = TRUE)

            if(coefstat == "se"){
                coefstat_sentence = "Standard-Errors type"
            } else {
                coefstat_sentence = "Co-variance type"
            }

            info_muli_se = paste0(coefstat_sentence, "& ", paste(all_se_type, collapse = "&"), "\\\\\n")

            if(add_signif){
                info_SD = paste0(style$lines$foot, sd_intro, "Signif. Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
            } else {
                myAmpLine = paste0(paste0(rep(" ", length(depvar_list)+1), collapse="&"), "\\tabularnewline\n")
                info_SD = paste0(style$lines$foot, myAmpLine, "\\\\\n")
            }

        }

    } else {
        info_SD = ""
        info_muli_se = ""
        if(nchar(style$lines$bottom) == 0) style$lines$bottom = style$lines$foot
    }


    # Information on number of items
    supplemental_info = ""

    if(all(theta_list == "")) info_theta = ""

    #
    # Fit statistics
    #

    if(nchar(style$stats$title) == 0){
        fit_info = style$lines$sep
    } else if(style$stats$title == "\\midrule"){
        fit_info = "\\midrule "
    } else {
        fit_info = paste0(style$lines$sep, style$stats$title, "& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
    }

    fit_info = paste0(fit_info, "Observations& ", paste(addCommas(obs_list), collapse="&"), "\\\\\n")
    fit_info = paste0(fit_info, nb_FE_lines, info_convergence, info_muli_se)
    if(!all(sapply(fitstat_list, function(x) all(is.na(x))))){

        fit_names = attr(fitstat_list, "format_names")
        nb = length(fit_names)
        for(fit_id in 1:nb){
            fit = sapply(fitstat_list, function(x) x[[fit_id]])
            fit[is.na(fit)] = "--"
            fit_info = paste0(fit_info, fit_names[fit_id], " & ", paste0(fit, collapse = "&"), "\\\\\n")
        }
    }

    # Notes
    info_notes = ""
    if(nchar(notes) > 0){
        info_notes = paste0("\n", style$notes$title, notes, "\n")
    }

    # Stacking var and stat
    var_stack = c(variable_line, coef_lines, info_theta)
    stat_stack = c(fit_info)
    dum_stack = c(dumIntro, FE_lines, slope_intro, slope_lines)

    # Group
    for(i in seq_along(group)){
        gi = group[[i]]
        gi_format = yesNo[2 - gi]

        gi_name = names(group)[i]

        gi_full = ""
        gi_where = "var"

        if(grepl("^\\{", gi_name)){
            # style component
            gi_style = gsub("\\{|\\}.+", "", gi_name)
            gi_style_parsed = parse_style(gi_style, c("title", "where"))
            from = nchar(strsplit(gi_name, "}", fixed = TRUE)[[1]][1])
            gi_name = substr(gi_name, from + 2, nchar(gi_name))

            if(nchar(gi_style_parsed$where) == 0){
                gi_where = "var"
            } else {
                gi_where = check_value_plus(gi_style_parsed$where, "match(var, stats, fixef)", .message = "In 'group', the keyword 'where' accepts only three values: either 'var', 'fixef' or 'stats'.")
            }


            if(nchar(gi_style_parsed$title) > 0){
                if(gi_style_parsed$title == "\\midrule"){
                    gi_full = "\\midrule "
                } else {
                    gi_full = paste0(gi_style_parsed$title, "& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
                }

            }
        }

        gi_full = paste0(gi_full, gi_name, " & ", paste0(gi_format, collapse = " & "), "\\\\\n")

        if(gi_where == "var"){
            var_stack = c(var_stack, gi_full)
        } else if(gi_where == "fixef"){
            dum_stack = c(dum_stack, gi_full)
        } else {
            stat_stack = c(stat_stack, gi_full)
        }

    }

    # Extra lines
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
            el_format = numberFormatLatex(el)
        } else {
            el_format = el
        }

        el_name = names(extraline)[i]

        el_full = ""
        el_where = "var"

        if(grepl("^\\{", el_name)){
            # style component
            el_style = gsub("\\{|\\}.+", "", el_name)
            el_style_parsed = parse_style(el_style, c("title", "where"))
            from = nchar(strsplit(el_name, "}", fixed = TRUE)[[1]][1])
            el_name = substr(el_name, from + 2, nchar(el_name))

            if(nchar(el_style_parsed$where) == 0){
                el_where = "var"
            } else {
                el_where = check_value_plus(el_style_parsed$where, "match(var, stats, fixef)", .message = "In 'extraline', the keyword 'where' accepts only three values: either 'var', 'fixef' or 'stats'.")
            }

            if(nchar(el_style_parsed$title) > 0){
                if(el_style_parsed$title == "\\midrule"){
                    el_full = "\\midrule "
                } else {
                    el_full = paste0(el_style_parsed$title, "& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
                }
            }
        }

        el_full = paste0(el_full, el_name, " & ", paste0(el_format, collapse = " & "), "\\\\\n")

        if(el_where == "var"){
            var_stack = c(var_stack, el_full)
        } else if(el_where == "fixef"){
            dum_stack = c(dum_stack, el_full)
        } else {
            stat_stack = c(stat_stack, el_full)
        }

    }

    # Now we place the fixed-effects
    if(style$fixef$where %in% c("", "var")){
        var_stack = c(var_stack, dum_stack)
    } else {
        stat_stack = c(stat_stack, dum_stack)
    }

    res = c(supplemental_info, start_table, intro_latex, first_line, info_subtitles, model_line, info_family, var_stack, stat_stack, info_SD, style$lines$bottom, outro_latex, info_notes, end_table)

    # res = c(supplemental_info, start_table, intro_latex, first_line, info_subtitles, model_line, info_family, variable_line, coef_lines, info_theta, dumIntro, FE_lines, slope_intro, slope_lines, fit_info, info_SD, style$lines$bottom, outro_latex, info_notes, end_table)

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
    obs_list = info$obs_list
    convergence_list = info$convergence_list
    fe_names = info$fe_names
    is_fe = info$is_fe
    nb_fe = info$nb_fe
    slope_names = info$slope_names
    slope_flag_list = info$slope_flag_list
    family_list = info$family_list
    theta_list = info$theta_list
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

    # naming differences
    titles = info$subtitles
    isTitles = info$isSubtitles

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
        present = rep(FALSE, n_models)
        for(v in gi){
            if(grepl("^%", v)){
                # original names
                v = substr(v, 2, nchar(v))
                present = present | sapply(var_list, function(x) any(grepl(v, names(x))))
            } else {
                present = present | sapply(var_list, function(x) any(grepl(v, x)))
            }
        }

        group[[i]] = present
    }

    # keeping some coefs
    all_vars = keep_apply(all_vars, keep)

    # dropping some coefs
    all_vars = drop_apply(all_vars, drop)

    # ordering the coefs
    all_vars = order_apply(all_vars, order)

    coef_mat <- all_vars
    for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
    coef_mat[is.na(coef_mat)] <- "  "
    res = coef_mat

    # Group
    for(i in seq_along(group)){
        gi = group[[i]]
        gi_format = c("Yes", "No")[2 - gi]

        gi_name = names(group)[i]
        from = nchar(strsplit(gi_name, "}", fixed = TRUE)[[1]][1])
        gi_name = substr(gi_name, from + 2, nchar(gi_name))

        my_line = c(gi_name, gi_format)

        res = rbind(res, my_line)
    }

    # Extra lines
    for(i in seq_along(extraline)){
        el = extraline[[i]]
        yesNo = c("Yes", "No")

        # The format depends on the type
        if(is.logical(el)){
            if(length(el) == 1){
                el_format = rep(yesNo[2 - el], n_models)
            } else {
                el_format = yesNo[2 - el]
            }
        } else if(is.numeric(el)){
            el_format = numberFormatLatex(el)
        } else {
            el_format = el
        }

        el_name = names(extraline)[i]
        from = nchar(strsplit(el_name, "}", fixed = TRUE)[[1]][1])
        el_name = substr(el_name, from + 2, nchar(el_name))

        my_line = c(el_name, el_format)

        res = rbind(res, my_line)
    }

    if("Neg. Bin." %in% family_list){
        theta_line = c("Overdispersion:", unlist(theta_list))
        res = rbind(res, theta_line)
    }

    # The line with the dependent variable => defined here to get the width
    preamble = c()
    dep_width = 0
    if(depvar){
        preamble = rbind(c("Dependent Var.:", depvar_list), preamble)
        dep_width = nchar(as.vector(preamble))
    }

    # Used to draw a line
    myLine = "______________________________________"
    longueur = apply(res, 2, function(x) max(nchar(as.character(x))))
    longueur = pmax(dep_width, longueur)
    theLine = sapply(longueur, function(x) sprintf("%.*s", x, myLine))
    theLine[1] = sprintf("%.*s", max(nchar(theLine[1]), 19), myLine)

    # The FEs
    if(length(fe_names)>0){

        for(m in 1:n_models) {
            quoi = is_fe[[m]][fe_names]
            quoi[is.na(quoi)] = "No"
            is_fe[[m]] = quoi
        }
        all_fe = matrix(c(is_fe, recursive=TRUE), nrow = length(fe_names))
        all_fe = cbind(fe_names, all_fe)
        FE_lines <- paste0(paste0(apply(all_fe, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

        myLine = "-------------------------------"

        res = rbind(res, c("Fixed-Effects:", sprintf("%.*s", longueur[-1], myLine)))
        factmat = matrix(c(strsplit(strsplit(FE_lines, "\n")[[1]], "&"), recursive = TRUE), ncol=n_models+1, byrow=TRUE)
        factmat[, ncol(factmat)]=gsub("\\", "", factmat[, ncol(factmat)], fixed = TRUE)
        res = rbind(res, factmat)
    }

    # The slopes
    if(length(slope_names) > 0){
        # reformatting the yes/no
        for(m in 1:n_models) {
            quoi = slope_flag_list[[m]][slope_names]
            quoi[is.na(quoi)] = "No"
            slope_flag_list[[m]] = quoi
        }
        all_slopes = matrix(c(slope_flag_list, recursive=TRUE), nrow = length(slope_names))
        all_slopes = cbind(slope_names, all_slopes)
        slope_lines <- paste0(paste0(apply(all_slopes, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

        myLine = "-------------------------------"

        res = rbind(res, c("Varying Slopes:", sprintf("%.*s", longueur[-1], myLine)))
        slope_mat = matrix(c(strsplit(strsplit(slope_lines, "\n")[[1]], "&"), recursive = TRUE), ncol=n_models+1, byrow = TRUE)
        slope_mat[, ncol(slope_mat)] = gsub("\\", "", slope_mat[, ncol(slope_mat)], fixed = TRUE)
        res = rbind(res, slope_mat)

    }


    # preamble created before because used to set the width
    if(length(preamble) > 0){
        # preamble = rbind(preamble, c("  ", theLine[-1]))
        preamble = rbind(preamble, rep("   ", length(theLine)))
        res <- rbind(preamble, res)
    }

    res <- rbind(res, theLine)

    # the line with the families
    if(family){
        # preamble = rbind(c("Family:", unlist(family_list)), preamble)
        res = rbind(res, c("Family", unlist(family_list)))
    }

    if(coefstat == "se"){
        coefstat_sentence = "S.E. type"
    } else {
        coefstat_sentence = "VCOV type"
    }

    res <- rbind(res, c("Observations", addCommas(obs_list)))
    if(!useSummary || !any(grepl("\\(", unlist(se_type_list)))){
        se_type_format = c()

        if(length(unique(gsub(" \\(.+", "", se_type_list))) == 1){
            # All identical
            type = paste0(": ", gsub(" \\(.+", "", se_type_list[1]))
            for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]], by = TRUE)
        } else {
            type = ""
            for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]])
        }

        res <- rbind(res, c(paste0(coefstat_sentence, type), c(se_type_format, recursive = TRUE)))
    } else {
        main_type = gsub(" \\(.+", "", se_type_list[[1]])
        se_type_format = c()
        for(m in 1:n_models) se_type_format[m] = format_se_type(se_type_list[[m]], longueur[[1+m]], by = TRUE)
        res <- rbind(res, c(paste0(coefstat_sentence, ": ", main_type), c(se_type_format, recursive = TRUE)))
    }

    # convergence status
    if(convergence){
        res <- rbind(res, c("Convergence", c(convergence_list, recursive = TRUE)))
    }

    #
    # Fit statistics
    #

    if(!identical(fitstat_list, NA)){

        fit_names = attr(fitstat_list, "format_names")
        nb = length(fit_names)
        for(fit_id in 1:nb){
            fit = sapply(fitstat_list, function(x) x[[fit_id]])
            fit[is.na(fit)] = "--"
            res <- rbind(res, c(fit_names[fit_id], fit))
        }
    }

    # if titles
    if(isTitles){
        modelNames = titles
    } else {
        # modelNames = paste0("model ", 1:n_models)
        modelNames = model_names
    }

    # we shorten the model names to fit the width
    for(m in 1:n_models) modelNames[m] = charShorten(modelNames[[m]], longueur[[1+m]])

    res <- as.data.frame(res)
    names(res) <- c("variables", modelNames)
    tvar = table(res$variables)
    if(any(tvar > 1)){
        qui = which(res$variables %in% names(tvar)[tvar > 1])
        add_space = c("", " ")
        if(length(qui) > 2) for(i in 3:length(qui)) add_space[i] = paste(rep(" ", i), collapse = "")
        res$variables = as.character(res$variables)
        res$variables[qui] = paste0(res$variables[qui], add_space)
    }
    row.names(res) = res$variables
    res$variables = NULL

    # We rename theta when NB is used
    quiTheta = which(row.names(res) == ".theta")
    row.names(res)[quiTheta] = "Dispersion Parameter"

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
setFixest_etable = function(digits=4, fitstat, coefstat = c("se", "tstat", "confint"), ci = 0.95, sdBelow = TRUE, keep, drop, order, dict, signifCode, float, fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, yesNo = c("Yes", "No"), family, powerBelow = -5, interaction.combine = " $\\times $ ", depvar, style = list(), notes = NULL, group = NULL, extraline = NULL, tablefoot = TRUE, reset = FALSE){

    # cat(names(formals(setFixest_etable)), sep = '", "')
    arg_list = c("digits", "fitstat", "coefstat", "ci", "sdBelow", "keep", "drop", "order", "dict", "signifCode", "float", "fixef_sizes", "fixef_sizes.simplify", "yesNo", "family", "powerBelow", "interaction.combine", "depvar", "style", "notes", "group", "extraline", "tablefoot", "reset")

    #
    # Argument checking => strong since these will become default values
    #

    check_arg(digits, "integer scalar GE{1}")

    # fitstat (chiant)
    if(!missing(fitstat)){
        check_arg(fitstat, "NA | os formula | charin(FALSE) | character vector no na")
        if("formula" %in% class(fitstat)){
            fitstat = attr(terms(fitstat), "term.labels")
        } else if (length(fitstat) == 1 && (isFALSE(fitstat) || is.na(fitstat))){
            fitstat = c()
        }

        # checking the types
        fitstat_type_allowed = c("sq.cor", "r2", "ar2", "pr2", "apr2", "par2", "wr2", "war2", "awr2", "wpr2", "pwr2", "wapr2", "wpar2", "awpr2", "apwr2", "pawr2", "pwar2", "ll", "aic", "bic")
        fitstat = unique(fitstat)

        pblm = setdiff(fitstat, fitstat_type_allowed)
        if(length(pblm) > 0){
            stop_up("Argument 'fitstat' must be a character vector (or a one sided formula) containing 'aic', 'bic', 'll', or valid r2 types names. ", enumerate_items(pblm, "is.quote"), " not valid (see function r2).")
        }

        if(length(fitstat) == 0) fitstat = NA
    }


    check_arg_plus(coefstat, "match")
    check_arg(ci, "numeric scalar GT{0.5} LT{1}")

    check_arg("logical scalar", sdBelow, fixef_sizes, fixef_sizes.simplify, float, family, depvar, tablefoot, reset)

    check_arg(keep, drop, order, "character vector no na NULL", .message = "The arg. '__ARG__' must be a vector of regular expressions (see help(regex)).")

    check_arg_plus(signifCode, "NULL NA | match(letters) | named numeric vector no na GE{0} LE{1}")

    check_arg(interaction.combine, "character scalar")

    check_arg(notes, "character vector no na")

    check_arg(yesNo, "character vector no na len(,2)")

    check_arg(powerBelow, "integer scalar LE{-1}")

    check_arg(dict, "NULL logical scalar | named character vector no na")

    check_arg_plus(group, extraline, "NULL{list()} named list l0")

    if(!missnull(style)){
        # We check the argument style + set it // all checks are made here
        check_arg_plus(style, "NULL{list()} named list l0")

        all_sections = strsplit("depvar, model, var, fixef, fixef.sizes, slopes, stats, lines, notes", ", ")[[1]]
        check_value(names(style), "NULL multi charin", .message = paste0("The names of the argument 'style' must be one of ", enumerate_items(all_sections, "or.quote", nmax = 100), "."), .choices = all_sections)

        if(!is.null(style$depvar)){
            check_value(style$depvar, "character scalar")
            parse_style(style$depvar, "title")
        }

        if(!is.null(style$model)){
            check_value(style$model, "character scalar")
            parse_style(style$model, c("title", "format"))
        }

        if(!is.null(style$lines)){
            check_value(style$lines, "character scalar")
            parse_style(style$lines, c("top", "bottom", "sep", "foot"))
        }

        if(!is.null(style$var)){
            check_value(style$var, "character scalar")
            parse_style(style$var, "title")
        }

        if(!is.null(style$fixef)){
            check_value(style$fixef, "character scalar")
            parse_style(style$fixef, c("title", "prefix", "suffix", "where"))
        }

        if(!is.null(style$slopes)){
            check_value(style$slopes, "character scalar")
            parse_style(style$slopes, c("title", "format"))
        }

        if(!is.null(style$fixef.sizes)){
            check_value(style$fixef.sizes, "character scalar")
            parse_style(style$fixef.sizes, c("title", "prefix", "suffix", "where"))
        }

        if(!is.null(style$stats)){
            check_value(style$stats, "character scalar")
            parse_style(style$stats, "title")
        }

        if(!is.null(style$notes)){
            check_value(style$notes, "character scalar")
            parse_style(style$notes, "title")
        }
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




