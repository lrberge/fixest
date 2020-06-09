#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sun Mar 22 09:04:17 2020
# ~: Only edit in fplot
#----------------------------------------------#


# Input of line.min / line.max => the default bottom margin goes from 0 to 4 (5 lines)
# Starts at 0!
# Beware: we return line => but in the axis sense! (input 0 is output -1)
xaxis_labels = function(at, labels, line.min = 0, line.max = 2, minCex = 0.8, add_ticks = FALSE, trunc = 20, trunc.method = "auto", only.params = FALSE, ...){
    # This function automates the placement of labels into the axis 1
    # It first put the cex to the appropritate level to insert the 1st
    # label into the frame, then
    # It uses the vertical placement to fit all the labels
    # if only.params => nothing is plotted

    # only 1 value => we do nothing
    if(length(at) == 1){
        if(only.params){
            res = list(cex = 1, line = line.min - 1)
            return(res)
        }
        axis(1, at=at, labels = labels, tick = add_ticks, line = line.min - 1)
        return(invisible(NULL))
    }

    # To send into a function: myat & lab
    myOrder = order(at)
    myLabels = labels[myOrder]
    myAt = at[myOrder]

    n = length(myAt)

    # # Truncation of the items
    # myLabels = substr(myLabels, 1, trunc)
    # # we replace the trucated element with dots to show that there is more
    # qui = nchar(myLabels) == trunc & nchar(labels[myOrder]) > trunc
    # myLabels[qui] = gsub("..$", "\\.\\.", myLabels[qui])
    if(!"call" %in% class(myLabels)){
        if(any(nchar(myLabels) > trunc)){
            # myLabels = truncate_string(myLabels, trunc, method = trunc.method)
        }

    } else {
        myLabels = gsub(" *phantom\\([\\)]*\\) *", "", deparse(myLabels))
    }


    # We compute the space that is left to display the label
    # 1st into the plot
    largeur = diff(par("usr")[1:2])
    half_plotSpace_in = (myAt[1] - par("usr")[1]) / largeur * par("pin")[1]
    # 2nd using the margin
    total_half_space = half_plotSpace_in + min(par("mai")[c(2,4)])

    # If it is too large, we reduce it:
    myCex = 1
    while(myCex > minCex && strwidth(myLabels[1], units = "in", cex = myCex)/2 > total_half_space){
        myCex = myCex * 0.95
    }

    if(myCex < minCex) myCex = minCex

    # Now we order the vertical space where the labels will appear
    # line_step = max(strheight(myLabels, units = "in", cex = myCex)) / max(strheight(myLabels, units = "in"))
    # we use the fact that:
    # 1) smaller characters can be more stacked vertically
    # 2) a line height is almost equivalent to character height

    # all_lines = -1:line.max
    ok = FALSE
    failed = FALSE
    while(!ok){
        ok = TRUE # if !ok there's a problem, we need to reduce cex

        all_width = strwidth(myLabels, units = "in", cex = myCex)

        # there can be more lines than expected depending on cex
        line_step = max(strheight(myLabels, units = "in", cex = myCex)) / max(strheight(myLabels, units = "in"))
        all_lines = seq(line.min, line.max, by = line_step)

        myLine = current_Line = line.min
        for(i in 2:n){
            # for each element, we find out the line where to write it
            for(line_index in all_lines){
                # we get the closest index with that Line level
                if(line_index %in% myLine){
                    index = max(which(myLine == line_index))
                    # we look at the distance between the two elements
                    at_first = myAt[index]
                    at_second = myAt[i]
                    # the distance in inches between the two 'at' + the space of one letter
                    dist_in = (at_second - at_first) * par("pin")[1] / largeur - strwidth("O", units = "in", cex = 1)
                    # the half sizes of the two strings
                    half_sums_in = (all_width[index] + all_width[i]) / 2
                    if(half_sums_in > dist_in){
                        # we go to the next line_index
                    } else {
                        myLine = c(myLine, line_index)
                        break
                    }
                } else {
                    # this line item has not been used already
                    myLine = c(myLine, line_index)
                    break
                }

                if(line_index == max(all_lines)) {
                    # Means it does not fit => we need to reduce cex
                    if(myCex <= minCex){
                        # already at the minimum, we keep going then
                        myLine = c(myLine, line_index)
                        failed = TRUE
                    } else {
                        # we get out totally and reduce the cex
                        ok = FALSE
                        myCex = myCex * 0.95
                    }
                }
            }

            if(!ok){
                # This means we've been into a non solvable situation
                break
            }
        }
    }

    if(only.params){
        # we substract 1 to make it "axis compatible"
        res = list(cex = myCex, line = myLine - 1 + 0.2, failed = failed)
        line_height = par("mai")[1] / par("mar")[1]
        res$height_in = (max(myLine) + 1.2) * line_height
        res$height_line = max(myLine) + 1.2
        return(res)
    }

    # We draw the ticks
    if(add_ticks){
        # 1) drawing them
        for(line in unique(myLine)){
            qui = which(myLine == line)
            axis(1, at = myAt[qui], labels = NA, tcl = -1-line, lwd = 0, lwd.ticks = 1)

        }

        # 2) "Cleaning" them
        n = length(myLabels)
        mid_width = strwidth(myLabels, units = "user") / 2

        # ceux qui 'debordent' a droite
        qui = which(mid_width[-n] > diff(myAt))
        if(length(qui) > 0){
            for(i in qui){
                # axis(1, at = myAt[i+1], labels = NA, lwd = 0, col = "white", line = myLine[i], lwd.ticks = 1.5)
                axis(1, at = myAt[i+1], labels = "|", lwd = 0, col.axis = "white", cex.axis = 1.5, line = myLine[i], lwd.ticks = 0)
            }
        }

        # ceux qui 'debordent' a gauche
        qui = which(mid_width[-1] > diff(myAt))
        if(length(qui) > 0){
            for(i in qui){
                # axis(1, at = myAt[i+1], labels = NA, lwd = 0, col = "white", line = myLine[i], lwd.ticks = 1.5)
                axis(1, at = myAt[i], labels = "|", lwd = 0, col.axis = "white", cex.axis = 1.5, line = myLine[i+1], lwd.ticks = 0)
            }
        }

    }

    # We draw the labels
    myLine = myLine + 0.2 - 1 # -1 to make it axis compatible
    for(line in unique(myLine)){
        qui = which(myLine == line)
        # the labels
        axis(1, at = myAt[qui], labels = myLabels[qui], line = line, cex.axis = myCex, lwd = 0, gap.axis = 0.01)

    }

    res = list(cex = myCex, line = myLine)
    line_height = par("mai")[1] / par("mar")[1]
    res$height_in = (max(myLine) + 2) * line_height
    res$height_line = max(myLine) + 2
    return(res)

}

# line: goes from 0 to 4 in a standard plot
xaxis_biased = function(at, labels, angle, cex, line.min = 0, line.max = 2, yadj = 0.5, trunc = 20, trunc.method = "auto", only.params = FALSE, ...){

    check_arg(angle, "null numeric vector no na")
    check_arg(cex, "null numeric vector no na")

    if(line.max < line.min){
        message("xaxis_biased: line.max < line.min (i.e. ", line.max, " < ", line.min, ") line.max set to ", line.min, ".")
        line.max = line.min
    }


    dots = list(...)
    dots$x = at

    if(!"call" %in% class(labels)){
        if(any(nchar(labels) > trunc)){
            # labels_trunc = truncate_string(labels, trunc = trunc, method = trunc.method)
            labels_trunc = labels
        } else {
            labels_trunc = labels
        }
    } else {
        labels_trunc = gsub(" *phantom\\([\\)]*\\) *", "", deparse(labels))
    }

    dots$labels = labels_trunc

    # setting automatically the cex and angle
    DO_ALGO = FALSE
    if(missnull(angle)){
        angle2check = c(45, 40, 35)
        DO_ALGO = TRUE
    } else {
        angle2check = angle
        DO_ALGO = length(angle) > 1
    }

    if(missnull(cex)){
        cex2check = c(1, 0.9, 0.8)
        DO_ALGO = TRUE
    } else {
        cex2check = cex
        DO_ALGO = length(cex) > 1
    }

    line_height = par("mai")[1] / par("mar")[1]

    if(DO_ALGO){
        lab_max = labels_trunc[which.max(strwidth(labels_trunc))]
        n_angle = length(angle2check)
        w_all = rep(sapply(cex2check, function(x) strwidth(lab_max, "in", cex = x)), n_angle)*1.05
        SH_all = rep(sapply(cex2check, function(x) strheight("W", "in", cex = x)), n_angle)
        angle_all = rep(angle2check, each = length(cex2check))
        h_all = SH_all / cos(angle_all / 360 * 2 * pi)
        longueur_cote = sin(angle_all / 360 * 2 * pi) * w_all + h_all

        total_height = line_height * (line.max - line.min + 1 - yadj)

        qui = longueur_cote <= total_height
        if(any(qui)){
            i = which.max(qui)
            angle = angle_all[i]
            cex = rep(cex2check, n_angle)[i]
            height_in = longueur_cote[i]
        } else {
            # message("xaxis_biased: Labels could not fit.")
            angle = tail(angle_all, 1)
            cex = tail(cex2check, 1)
            height_in = tail(longueur_cote, 1)
        }

    } else if(only.params) {
        lab_max = labels_trunc[which.max(strwidth(labels_trunc))]
        n_angle = length(angle2check)
        w_all = rep(sapply(cex2check, function(x) strwidth(lab_max, "in", cex = x)), n_angle)*1.05
        SH_all = rep(sapply(cex2check, function(x) strheight("W", "in", cex = x)), n_angle)
        angle_all = rep(angle2check, each = length(cex2check))
        h_all = SH_all / cos(angle_all / 360 * 2 * pi)
        longueur_cote = sin(angle_all / 360 * 2 * pi) * w_all + h_all

        height_in = longueur_cote
    }

    line_height_usr = line_height / par("pin")[2] * diff(par("usr")[3:4])

    if(only.params){
        res = list(cex = cex, angle = angle, height_in = height_in + strheight("W", units = "in", cex = cex))
        res$height_usr = height_in / par("pin")[2] * diff(par("usr")[3:4])
        res$height_line = height_in / line_height
        return(res)
    }

    dots$cex = cex
    SH = strheight("WWW", units = "user")
    dots$y = par("usr")[3] - yadj*SH - (line.min) * line_height_usr
    dots$srt = angle
    dots$xpd = TRUE
    dots$adj = 1

    # browser()

    do.call("text", dots)

    return(invisible(list(cex=cex, angle=angle)))
}



