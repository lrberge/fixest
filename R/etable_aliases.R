# Do not edit by hand
# => aliases to the function etable








#' @describeIn etable Exports the results of multiple \code{fixest} estimations in a Latex table.
esttable = function(..., se = NULL, dof = NULL, cluster = NULL, stage = 2, agg = NULL, .vcov, .vcov_args = NULL, digits = 4, digits.stats = 5, fitstat, coefstat = "se", ci = 0.95, sdBelow = NULL, keep, drop, order, dict, file, replace = FALSE, convergence, signifCode, subtitles = list("auto"), fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, keepFactors = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ ", depvar = TRUE, style.df = NULL, group = NULL, extraline = NULL, fixef.group = NULL, drop.section = NULL, poly_dict = c("", " square", " cube"), postprocess.df = NULL, fit_format = "__var__", coef.just = NULL){

	etable(..., se = se, dof = dof, cluster = cluster, stage = stage, agg = agg, .vcov = .vcov, .vcov_args = .vcov_args, digits = digits, digits.stats = digits.stats, fitstat = fitstat, coefstat = coefstat, ci = ci, sdBelow = sdBelow, keep = keep, drop = drop, order = order, dict = dict, file = file, replace = replace, convergence = convergence, signifCode = signifCode, subtitles = subtitles, fixef_sizes = fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, keepFactors = keepFactors, family = family, powerBelow = powerBelow, interaction.combine = interaction.combine, depvar = depvar, style.df = style.df, group = group, extraline = extraline, fixef.group = fixef.group, drop.section = drop.section, poly_dict = poly_dict, postprocess.df = postprocess.df, fit_format = fit_format, coef.just = coef.just, tex = FALSE, .up = 2)
}





#' @describeIn etable Exports the results of multiple \code{fixest} estimations in a Latex table.
esttex = function(..., se = NULL, dof = NULL, cluster = NULL, stage = 2, agg = NULL, .vcov, .vcov_args = NULL, digits = 4, digits.stats = 5, fitstat, title, coefstat = "se", ci = 0.95, sdBelow = NULL, keep, drop, order, dict, file, replace = FALSE, convergence, signifCode, label, float, subtitles = list("auto"), fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, keepFactors = TRUE, family, powerBelow = -5, interaction.combine = " $\\times $ ", depvar = TRUE, style.tex = NULL, notes = NULL, group = NULL, extraline = NULL, fixef.group = NULL, placement = "htbp", drop.section = NULL, poly_dict = c("", " square", " cube"), postprocess.tex = NULL, fit_format = "__var__"){

	etable(..., se = se, dof = dof, cluster = cluster, stage = stage, agg = agg, .vcov = .vcov, .vcov_args = .vcov_args, digits = digits, digits.stats = digits.stats, fitstat = fitstat, title = title, coefstat = coefstat, ci = ci, sdBelow = sdBelow, keep = keep, drop = drop, order = order, dict = dict, file = file, replace = replace, convergence = convergence, signifCode = signifCode, label = label, float = float, subtitles = subtitles, fixef_sizes = fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, keepFactors = keepFactors, family = family, powerBelow = powerBelow, interaction.combine = interaction.combine, depvar = depvar, style.tex = style.tex, notes = notes, group = group, extraline = extraline, fixef.group = fixef.group, placement = placement, drop.section = drop.section, poly_dict = poly_dict, postprocess.tex = postprocess.tex, fit_format = fit_format, tex = TRUE, .up = 2)
}





