# Do not edit by hand
# => iplot calls coefplot internally








#' @describeIn coefplot Plots the coefficients generated with i()
iplot = function(object, ..., style = NULL, sd, ci_low, ci_high, x, x.shift = 0, horiz = FALSE, dict = getFixest_dict(), keep, drop, order, ci.width = "1%", ci_level = 0.95, add = FALSE, pt.pch = c(20, 17, 15, 21, 24, 22), pt.bg = NULL, cex = 1, pt.cex = cex, col = 1:8, pt.col = col, ci.col = col, lwd = 1, pt.lwd = lwd, ci.lwd = lwd, ci.lty = 1, grid = TRUE, grid.par = list(lty = 3, col = "gray"), zero = TRUE, zero.par = list(col = "black", lwd = 1), pt.join = FALSE, pt.join.par = list(col = pt.col, lwd = lwd), ci.join = FALSE, ci.join.par = list(lwd = lwd, col = col, lty = 2), ci.fill = FALSE, ci.fill.par = list(col = "lightgray", alpha = 0.5), ref = "auto", ref.line = "auto", ref.line.par = list(col = "black", lty = 2), lab.cex, lab.min.cex = 0.85, lab.max.mar = 0.25, lab.fit = "auto", xlim.add, ylim.add, only.params = FALSE, sep, as.multiple = FALSE, bg, group = "auto", group.par = list(lwd = 2, line = 3, tcl = 0.75), main = "Effect on __depvar__", value.lab = "Estimate and __ci__ Conf. Int.", ylab = NULL, xlab = NULL, sub = NULL){

	coefplot(object = object, ..., style = style, sd = sd, ci_low = ci_low, ci_high = ci_high, x = x, x.shift = x.shift, horiz = horiz, dict = dict, keep = keep, drop = drop, order = order, ci.width = ci.width, ci_level = ci_level, add = add, pt.pch = pt.pch, pt.bg = pt.bg, cex = cex, pt.cex = pt.cex, col = col, pt.col = pt.col, ci.col = ci.col, lwd = lwd, pt.lwd = pt.lwd, ci.lwd = ci.lwd, ci.lty = ci.lty, grid = grid, grid.par = grid.par, zero = zero, zero.par = zero.par, pt.join = pt.join, pt.join.par = pt.join.par, ci.join = ci.join, ci.join.par = ci.join.par, ci.fill = ci.fill, ci.fill.par = ci.fill.par, ref = ref, ref.line = ref.line, ref.line.par = ref.line.par, lab.cex = lab.cex, lab.min.cex = lab.min.cex, lab.max.mar = lab.max.mar, lab.fit = lab.fit, xlim.add = xlim.add, ylim.add = ylim.add, only.params = only.params, sep = sep, as.multiple = as.multiple, bg = bg, group = group, group.par = group.par, main = main, value.lab = value.lab, ylab = ylab, xlab = xlab, sub = sub, internal.only.i = TRUE)
}





