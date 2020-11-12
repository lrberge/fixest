#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Nov 07 09:05:26 2020
# ~: fixest_multi
#----------------------------------------------#


setup_multi = function(index, all_names, data){
    # Basic setup function


    if("fixest_multi" %in% class(data[[1]])){
        # We need to "grow" the results

        new_data = list()
        for(i in seq_along(data)){
            my_res = data[[i]]

            if(i == 1){
                # setup in the first iteration
                meta = attr(my_res, "meta")
                old_index = meta$index
                old_all_names = meta$all_names
            }

            my_data = attr(my_res, "data")

            for(j in seq_along(my_data)){
                new_data[[length(new_data) + 1]] = my_data[[j]]
            }

        }

        #
        # we create the new information
        #

        # the length of all_names need not be the same as index
        new_all_names = all_names
        for(i in seq_along(old_all_names)){
            new_all_names[[names(old_all_names)[i]]] = old_all_names[[i]]
        }

        new_index = index
        for(i in seq_along(old_index)){
            new_index[[names(old_index)[i]]] = old_index[[i]]
        }

        # now we're ready to go
        index = new_index
        all_names = new_all_names
        data = new_data

    }


    # We drop non-used dimensions
    n_all = sapply(index, max)
    qui = n_all == 1
    if(any(qui)){
        index = index[!qui]
        # all_names = all_names[!qui]
        # We don't drop the names dimension => we always keep it as an imprint
    }

    meta = list(index = index, all_names = all_names)
    tree = as.data.frame(do.call("expand.grid", rev(lapply(index, seq))))
    meta$tree = tree[, rev(1:ncol(tree)), drop = FALSE]
    meta$tree$obs = 1:nrow(tree)

    res_multi = list()
    first_dim = all_names[[names(index)[1]]]
    for(i in seq_along(first_dim)){
        res_multi[[first_dim[i]]] = 1
    }

    attr(res_multi, "meta") = meta
    attr(res_multi, "data") = data

    class(res_multi) = "fixest_multi"

    return(res_multi)
}

#' Summary for fixest_multi objects
#'
#' Summary information for fixest_multi objects. In particular, this is used to specify the type of standard-errors to be computed.
#'
#' @method summary fixest_multi
#'
#' @inherit print.fixest_multi seealso
#'
#' @param object A \code{fixest_multi} object, obtained from a \code{fixest} estimation leading to multiple results.
#' @param type A character either equal to \code{"short"}, \code{"long"}, \code{"compact"}, or \code{"se_compact"}. If \code{short}, only the table of coefficients is displayed for each estimation. If \code{long}, then the full results are displayed for each estimation. If \code{compact}, a \code{data.frame} is returned with one line per model and the formatted coefficients + standard-errors in the columns. If \code{se_compact}, a \code{data.frame} is returned with one line per model, one numeric column for each coefficient and one numeric column for each standard-error.
#' @param ... Other arguments to be passed to \code{\link[fixest]{summary.fixest}}.
#'
#' @return
#' It returns either an object of class \code{fixest_multi} (if \code{type} equals \code{short} or \code{long}), either a \code{data.frame} (if type equals \code{compact} or \code{se_compact}).
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # Multiple estimation
#' res = feols(y ~ csw(x1, x2, x3), base, split = ~species)
#'
#' # By default, the type is "short"
#' # You can stil use the arguments from summary.fixest
#' summary(res, cluster = ~ species)
#'
#' summary(res, type = "long")
#'
#' summary(res, type = "compact")
#'
#' summary(res, type = "se_compact")
#'
#'
summary.fixest_multi = function(object, type = "short", ...){
    dots = list(...)
    data = attr(object, "data")

    check_arg_plus(type, "match(short, long, compact, se_compact)")

    if(!missing(type) || is.null(attr(object, "print_request"))){
        attr(object, "print_request") = type
    }

    est_1 = data[[1]]
    if(is.null(est_1$cov.scaled) || !isTRUE(dots$fromPrint)){

        for(i in 1:length(data)){
            data[[i]] = summary(data[[i]], nframes_up = 2, ...)
        }

        attr(object, "data") = data
    }

    if(type %in% c("compact", "se_compact")){
        meta = attr(object, "meta")
        data = attr(object, "data")

        index = meta$index
        all_names = meta$all_names
        tree = meta$tree

        res = data.frame(i = tree$obs)

        if(!"sample" %in% names(index) && !is.null(all_names$sample)){
            res$sample = all_names$sample
        }

        if(!"lhs" %in% names(index)){
            res$lhs = sapply(data, function(x) as.character(x$fml[[2]]))
        }
        for(my_dim in names(index)){
            res[[my_dim]] = sfill(as.character(factor(tree[[my_dim]], labels = all_names[[my_dim]])), right = TRUE)
        }
        res$i = NULL

        signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.1)

        ct_all = list()
        for(i in seq_along(data)){
            ct = data[[i]]$coeftable
            vname = row.names(ct)

            if(type == "compact"){
                stars = cut(ct[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
                stars[is.na(stars)] = ""

                value = paste0(signif_plus(ct[, 1], 3), stars, " (", signif_plus(ct[, 2], 3), ")")
                names(value) = vname

            } else if(type == "se_compact"){
                n = length(vname)
                vname_tmp = character(2 * n)
                qui_coef = seq(1, by = 2, length.out = n)
                qui_se   = seq(2, by = 2, length.out = n)
                vname_tmp[qui_coef] = vname
                vname_tmp[qui_se]   = paste0(vname, "__se")
                vname = vname_tmp

                value = numeric(2 * n)
                value[qui_coef] = ct[, 1]
                value[qui_se]   = ct[, 2]
                names(value) = vname
            }

            ct_all[[i]] = value
        }

        vname_all = unique(unlist(lapply(ct_all, names)))

        tmp = lapply(ct_all, function(x) x[vname_all])
        my_ct = do.call("rbind", tmp)
        colnames(my_ct) = vname_all
        if(type == "compact"){
            my_ct[is.na(my_ct)] = ""
        }

        for(i in seq_along(vname_all)){
            if(type == "compact"){
                res[[vname_all[i]]] = sfill(my_ct[, i], anchor = "(", right = TRUE)
            } else {
                res[[vname_all[i]]] = my_ct[, i]
            }
        }

        return(res)
    }

    return(object)
}


#' Print method for fixest_multi objects
#'
#' Displays summary information on fixest_multi objects in the R console.
#'
#' @method print fixest_multi
#'
#' @param x A \code{fixest_multi} object, obtained from a \code{fixest} estimation leading to multiple results.
#' @param ... Other arguments to be passed to \code{\link[fixest]{summary.fixest_multi}}.
#'
#' @seealso
#' The main fixest estimation functions: \code{\link[fixest]{feols}}, \code{\link[fixest:feglm]{fepois}}, \code{\link[fixest:femlm]{fenegbin}}, \code{\link[fixest]{feglm}}, \code{\link[fixest]{feNmlm}}. Tools for mutliple fixest estimations: \code{\link[fixest]{summary.fixest_multi}}, \code{\link[fixest]{print.fixest_multi}}, \code{\link[fixest]{as.list.fixest_multi}}, \code{\link[fixest]{sub-sub-.fixest_multi}}, \code{\link[fixest]{sub-.fixest_multi}}, \code{\link[fixest]{cash-.fixest_multi}}.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # Multiple estimation
#' res = feols(y ~ csw(x1, x2, x3), base, split = ~species)
#'
#' # Let's print all that
#' res
#'
print.fixest_multi = function(x, ...){

    x = summary(x, fromPrint = TRUE, ...)

    # Type = compact
    if(is.data.frame(x)){
        return(x)
    }

    is_short = identical(attr(x, "print_request"), "short")

    meta = attr(x, "meta")
    data = attr(x, "data")

    index = meta$index
    all_names = meta$all_names
    tree = meta$tree

    cat(signif_plus(do.call(prod, index)), " fixest estimations\n", sep = "")

    # Finding out the type of SEs
    if(is_short){
        est_1 = data[[1]]
        coeftable = est_1$coeftable
        se.type = attr(coeftable, "type")
        cat("Standard-errors:", se.type, "\n")
    }

    depth = length(index)

    dict_title = c("sample" = "Sample", "lhs" = "Dep. var.", "rhs" = "Expl. vars.")

    headers = list()
    headers[[1]] = function(d, i) cat(dict_title[d], ": ", all_names[[d]][i], "\n", sep = "")
    headers[[2]] = function(d, i) cat("\n### ", dict_title[d], ": ", all_names[[d]][i], "\n\n", sep = "")
    headers[[3]] = function(d, i) cat("\n\n# ", toupper(dict_title[d]), ": ", all_names[[d]][i], "\n\n", sep = "")
    headers[[4]] = function(d, i) cat("\n\n#\n# ", toupper(dict_title[d]), ": ", all_names[[d]][i], "\n#\n\n", sep = "")

    for(i in 1:nrow(tree)){
        for(j in 1:depth){
            d = names(index)[j]
            if(i == 1 || tree[i - 1, j] != tree[i, j]){
                headers[[depth - j + 1]](d, tree[i, j])
            }
        }

        if(is_short){
            myPrintCoefTable(coeftable = coeftable(data[[i]], ...), show_signif = FALSE)
            if(tree[i, depth] != index[[depth]]) cat("---\n")
        } else {
            print(data[[i]])
            if(tree[i, depth] != index[[depth]]) cat("\n")
        }

    }

}

#' Extracts one element from a \code{fixest_multi} object
#'
#' Extracts single elements from multiple \code{fixest} estimations.
#'
#' @method sub-sub- fixest_multi
#'
#' @inherit print.fixest_multi seealso
#' @inheritParams print.fixest_multi
#'
#' @param i An integer scalar. The identifier of the estimation to extract.
#'
#' @return
#' A \code{fixest} object is returned.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # Multiple estimation
#' res = feols(y ~ csw(x1, x2, x3), base, split = ~species)
#'
#' # The first estimation
#' res[[1]]
#'
#' # The second one, etc
#' res[[2]]
#'
"[[.fixest_multi" = function(x, i){

    data = attr(x, "data")

    check_arg_plus(i, "evalset integer scalar mbt", .data = list(.N = length(data)))

    return(data[[i]])
}


#' Subset a fixest_multi object
#'
#' Subset a fixest_multi object using different keys.
#'
#' @method sub- fixest_multi
#'
#' @inherit print.fixest_multi seealso
#' @inheritParams print.fixest_multi
#'
#' @param sample An integer vector of a character vector. It represents the \code{sample} identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation was a split sample. You can use \code{.N} to refer to the last element.
#' @param lhs An integer vector of a character vector. It represents the left-hand-sides identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation contained multiple left-hand-sides. You can use \code{.N} to refer to the last element.
#' @param rhs An integer vector. It represents the right-hand-sides identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation contained multiple right-hand-sides. You can use \code{.N} to refer to the last element.
#' @param i An integer vector. Represents the estimations to extract.
#' @param I An integer vector. Represents the root element to extract.
#'
#' @details
#' The order with we we use the keys matter. Every time a key sample, lhs or rhs is used, a reordering is performed to consider the leftmost-side key to be the new root.
#'
#' @return
#' It returns a \code{fixest_multi} object. If there is only one estimation left in the object, then the result is simplified into a \code{fixest} object.
#'
#' @examples
#'
#' # Estimation with multiple samples/LHS/RHS
#' aq = airquality[airquality$Month %in% 5:6, ]
#' est_split = feols(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
#'                   aq, split = ~ Month)
#'
#' # By default: sample is the root
#' etable(est_split)
#'
#' # Let's reorder, by considering lhs the root
#' etable(est_split[lhs = 1:.N])
#'
#' # Selecting only one LHS and RHS
#' etable(est_split[lhs = "Ozone", rhs = 1])
#'
#' # Taking the first root (here sample = 5)
#' etable(est_split[I = 1])
#'
#' # The first and last estimations
#' etable(est_split[i = c(1, .N)])
#'
"[.fixest_multi" = function(x, sample, lhs, rhs, i, I){

    mc = match.call()
    if(!any(c("sample", "lhs", "rhs", "i", "I") %in% names(mc))){
        return(x)
    }

    use_i = "i" %in% names(mc)
    if(use_i && any(c("sample", "lhs", "rhs", "I") %in% names(mc))){
        stop("The index 'i' cannot be used with any other index ('sample', 'lhs', 'rhs', or 'I').")
    }

    use_I = "I" %in% names(mc)
    if(use_I && any(c("sample", "lhs", "rhs") %in% names(mc))){
        stop("The index 'I' cannot be used with any other index ('i', 'sample', 'lhs', or 'rhs').")
    }

    # We get the meta information
    meta = attr(x, "meta")
    index = meta$index
    all_names = meta$all_names
    tree = meta$tree
    nc = ncol(tree)
    n = nrow(tree)

    data = attr(x, "data")

    if(!use_i && !use_I){
        pblm = setdiff(names(mc)[-(1:2)], names(index))
        if(length(pblm) > 0){
            stop("The ", ifsingle(pblm, "index", "indices"), " ", enumerate_items(pblm, "is"),
                 " not valid for this list of results (the valid one", plural_len(index, "s.is"), " ",  enumerate_items(names(index)), ").")
        }
    }

    if(use_i){
        check_arg_plus(i, "evalset integer vector l0 no na", .data = list(.N = n))

        if(length(i) == 0) return(list())

        if(any(abs(i) > n)){
            stop("The index 'i' cannot have values greater than ", n, ". Currently ", i[which.max(abs(i))], " is not valid.")
        }

        obs = (1:n)[i]

        return(data[obs])
    }

    if(use_I){
        I_max = index[[1]]
        check_arg_plus(I, "evalset integer vector no na", .data = list(.N = I_max))

        if(any(abs(I) > I_max)){
            stop("The index 'I' refers to root elements (here ", names(index)[1], "), and thus cannot be greater than ", I_max, ". Currently ", I[which.max(abs(I))], " is not valid.")
        }

        obs = (1:I_max)[I]

        # If only one dimension => result is simplified
        if(length(index) == 1 && length(obs) == 1){
            res = data[[obs]]

        } else {
            new_all_names = all_names
            new_all_names[[1]] = all_names[[1]][obs]
            obs_keep = tree[tree[, 1] %in% obs, "obs"]

            if(length(obs) == 1){
                # The first dimension is dissolved
                res = setup_multi(index[-1], new_all_names, data[obs_keep])
            } else {
                new_index = index
                new_index[[1]] = length(obs)
                res = setup_multi(new_index, new_all_names, data[obs_keep])
            }
        }

        return(res)
    }

    # Here:
    # We take care of reordering properly

    is_sample = !missing(sample)
    is_lhs    = !missing(lhs)
    is_rhs    = !missing(rhs)

    selection = list()

    s_max = index[["sample"]]
    if(is_sample){
        check_arg_plus(sample, "evalset multi charin | integer vector no na", .choices = all_names[["sample"]], .data = list(.N = s_max))


        if(is.character(sample)){
            dict = 1:s_max
            names(dict) = all_names[["sample"]]
            sample = as.integer(dict[sample])
        }

        if(any(abs(sample) > s_max)){
            stop("The index 'sample' cannot be greater than ", s_max, ". Currently ", I[which.max(abs(sample))], " is not valid.")
        }

        selection$sample = (1:s_max)[sample]
    } else if("sample" %in% names(index)){
        selection$sample = 1:s_max
    }

    lhs_max = index[["lhs"]]
    if(is_lhs){
        check_arg_plus(lhs, "evalset multi charin | integer vector no na", .choices = all_names[["lhs"]], .data = list(.N = lhs_max))

        if(is.character(lhs)){
            dict = 1:lhs_max
            names(dict) = all_names[["lhs"]]
            lhs = as.integer(dict[lhs])
        }

        if(any(abs(lhs) > lhs_max)){
            stop("The index 'lhs' cannot be greater than ", lhs_max, ". Currently ", I[which.max(abs(lhs))], " is not valid.")
        }

        selection$lhs = (1:lhs_max)[lhs]
    } else if("lhs" %in% names(index)){
        selection$lhs = 1:lhs_max
    }

    rhs_max = index[["rhs"]]
    if(is_rhs){
        check_arg_plus(rhs, "evalset multi charin | integer vector no na", .choices = all_names[["rhs"]], .data = list(.N = rhs_max))

        if(is.character(rhs)){
            dict = 1:rhs_max
            names(dict) = all_names[["rhs"]]
            rhs = as.integer(dict[rhs])
        }

        if(any(abs(rhs) > rhs_max)){
            stop("The index 'rhs' cannot be greater than ", rhs_max, ". Currently ", I[which.max(abs(rhs))], " is not valid.")
        }

        selection$rhs = (1:rhs_max)[rhs]
    } else if("rhs" %in% names(index)){
        selection$rhs = 1:rhs_max
    }

    # We keep the order of the user!!!!!
    sc = sys.call()
    user_order = names(sc)[-(1:2)]
    if(any(!user_order %in% names(index))){
        user_order = names(index)
    } else {
        user_order = c(user_order, setdiff(names(index), user_order))
    }

    new_index = list()
    new_all_names = list()
    new_data = list()

    new_index = list()
    new_all_names = list()
    for(d in user_order){
        new_index[[d]] = length(selection[[d]])
        new_all_names[[d]] = all_names[[d]][selection[[d]]]
    }
    new_all_names$split.name = all_names$split.name

    for(my_dim in rev(user_order)){
        # 1) we prune the tree
        obs_keep = tree[[my_dim]] %in% selection[[my_dim]]
        tree = tree[obs_keep, , drop = FALSE]
        # 2) we reorder it
        tree = tree[order(tree[[my_dim]]), ]
    }

    new_data = list()
    for(i in 1:nrow(tree)){
        new_data[[i]] = data[[tree$obs[i]]]
    }

    if(length(new_data) == 1){
        return(new_data[[1]])
    }

    res_multi = setup_multi(new_index, new_all_names, new_data)

    return(res_multi)
}


#' Extracts the root of a fixest_multi object
#'
#' Extracts an element at the root of a fixest_multi object.
#'
#' @method cash- fixest_multi
#'
#' @inherit print.fixest_multi seealso
#' @inheritParams print.fixest_multi
#'
#' @param name The name of the root element to select.
#'
#' @return
#' It either returns a \code{fixest_multi} object or a \code{fixest} object it there is only one estimation associated to the root element.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # Multiple estimation
#' res = feols(y ~ csw(x1, x2, x3), base, split = ~species)
#'
#' # Let's the results for the setosa species
#' res$setosa
#'
#' # now for versicolor
#' etable(res$versicolor)
#'
"$.fixest_multi" = function(x, name){

    # real variables
    if(!name %in% names(x)){
        mc = match.call()
        stop(name, " is not at the root of the multiple fixest estimations object.")
    }

    qui = which(name == names(x))

    return(x[I = qui])
}


#' Transforms a fixest_multi object into a list
#'
#' Extracts the results from a \code{fixest_multi} object and place them into a list.
#'
#' @inheritParams print.fixest_multi
#' @inherit print.fixest_multi seealso
#'
#' @method as.list fixest_multi
#'
#' @param ... Not currently used.
#'
#' @return
#' Returns a list containing all the results of the multiple estimations.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # Multiple estimation
#' res = feols(y ~ csw(x1, x2, x3), base, split = ~species)
#'
#' # All the results at once
#' as.list(res)
#'
#'
as.list.fixest_multi = function(x, ...){
    attr(x, "data")
}















