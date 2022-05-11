# Do not edit by hand
# => aliases to the function etable








#' @describeIn etable Exports the results of multiple \code{fixest} estimations in a Latex table.
esttable = function(..., vcov = NULL, stage = 2, agg = NULL, se = NULL, ssc = NULL, cluster = NULL, .vcov = NULL, .vcov_args = NULL, digits = 4, digits.stats = 5, fitstat = NULL, coefstat = "se", ci = 0.95, se.row = NULL, se.below = NULL, keep = NULL, drop = NULL, order = NULL, dict = TRUE, file = NULL, replace = FALSE, convergence = NULL, signif.code = NULL, headers = list("auto"), fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, keepFactors = TRUE, family = NULL, powerBelow = -5, interaction.combine = NULL, interaction.order = NULL, i.equal = NULL, depvar = TRUE, style.df = NULL, group = NULL, extralines = NULL, fixef.group = NULL, drop.section = NULL, poly_dict = c("", " square", " cube"), postprocess.df = NULL, fit_format = "__var__", coef.just = NULL, highlight = NULL, coef.style = NULL, export = NULL, page.width = "fit", div.class = "etable"){

	etable(..., vcov = vcov, stage = stage, agg = agg, se = se, ssc = ssc, cluster = cluster, .vcov = .vcov, .vcov_args = .vcov_args, digits = digits, digits.stats = digits.stats, fitstat = fitstat, coefstat = coefstat, ci = ci, se.row = se.row, se.below = se.below, keep = keep, drop = drop, order = order, dict = dict, file = file, replace = replace, convergence = convergence, signif.code = signif.code, headers = headers, fixef_sizes = fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, keepFactors = keepFactors, family = family, powerBelow = powerBelow, interaction.combine = interaction.combine, interaction.order = interaction.order, i.equal = i.equal, depvar = depvar, style.df = style.df, group = group, extralines = extralines, fixef.group = fixef.group, drop.section = drop.section, poly_dict = poly_dict, postprocess.df = postprocess.df, fit_format = fit_format, coef.just = coef.just, highlight = highlight, coef.style = coef.style, export = export, page.width = page.width, div.class = div.class, tex = FALSE, .up = 2)
}





#' @describeIn etable Exports the results of multiple \code{fixest} estimations in a Latex table.
esttex = function(..., vcov = NULL, stage = 2, agg = NULL, se = NULL, ssc = NULL, cluster = NULL, .vcov = NULL, .vcov_args = NULL, digits = 4, digits.stats = 5, fitstat = NULL, title = NULL, coefstat = "se", ci = 0.95, se.row = NULL, se.below = NULL, keep = NULL, drop = NULL, order = NULL, dict = TRUE, file = NULL, replace = FALSE, convergence = NULL, signif.code = NULL, label = NULL, float = NULL, headers = list("auto"), fixef_sizes = FALSE, fixef_sizes.simplify = TRUE, keepFactors = TRUE, family = NULL, powerBelow = -5, interaction.combine = NULL, interaction.order = NULL, i.equal = NULL, depvar = TRUE, style.tex = NULL, notes = NULL, group = NULL, extralines = NULL, fixef.group = NULL, placement = "htbp", drop.section = NULL, poly_dict = c("", " square", " cube"), postprocess.tex = NULL, tpt = FALSE, arraystretch = NULL, adjustbox = NULL, fontsize = NULL, fit_format = "__var__", tabular = "normal", highlight = NULL, coef.style = NULL, meta = NULL, meta.time = NULL, meta.author = NULL, meta.sys = NULL, meta.call = NULL, meta.comment = NULL, view = FALSE, export = NULL, markdown = NULL, page.width = "fit", div.class = "etable"){

	etable(..., vcov = vcov, stage = stage, agg = agg, se = se, ssc = ssc, cluster = cluster, .vcov = .vcov, .vcov_args = .vcov_args, digits = digits, digits.stats = digits.stats, fitstat = fitstat, title = title, coefstat = coefstat, ci = ci, se.row = se.row, se.below = se.below, keep = keep, drop = drop, order = order, dict = dict, file = file, replace = replace, convergence = convergence, signif.code = signif.code, label = label, float = float, headers = headers, fixef_sizes = fixef_sizes, fixef_sizes.simplify = fixef_sizes.simplify, keepFactors = keepFactors, family = family, powerBelow = powerBelow, interaction.combine = interaction.combine, interaction.order = interaction.order, i.equal = i.equal, depvar = depvar, style.tex = style.tex, notes = notes, group = group, extralines = extralines, fixef.group = fixef.group, placement = placement, drop.section = drop.section, poly_dict = poly_dict, postprocess.tex = postprocess.tex, tpt = tpt, arraystretch = arraystretch, adjustbox = adjustbox, fontsize = fontsize, fit_format = fit_format, tabular = tabular, highlight = highlight, coef.style = coef.style, meta = meta, meta.time = meta.time, meta.author = meta.author, meta.sys = meta.sys, meta.call = meta.call, meta.comment = meta.comment, view = view, export = export, markdown = markdown, page.width = page.width, div.class = div.class, tex = TRUE, .up = 2)
}





