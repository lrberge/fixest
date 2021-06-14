#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Nov 07 09:05:26 2020
# ~: fixest_multi
#----------------------------------------------#


setup_multi = function(index, all_names, data, simplify = TRUE){
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
    if(simplify){
        n_all = sapply(index, max)
        qui = n_all == 1
        if(any(qui)){
            index = index[!qui]
            # all_names = all_names[!qui]
            # We don't drop the names dimension => we always keep it as an imprint
        }
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
#' @inheritParams summary.fixest
#'
#' @inherit print.fixest_multi seealso
#'
#' @param object A \code{fixest_multi} object, obtained from a \code{fixest} estimation leading to multiple results.
#' @param type A character either equal to \code{"short"}, \code{"long"}, \code{"compact"}, or \code{"se_compact"}. If \code{short}, only the table of coefficients is displayed for each estimation. If \code{long}, then the full results are displayed for each estimation. If \code{compact}, a \code{data.frame} is returned with one line per model and the formatted coefficients + standard-errors in the columns. If \code{se_compact}, a \code{data.frame} is returned with one line per model, one numeric column for each coefficient and one numeric column for each standard-error.
#' @param ... Not currently used.
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
#' # You can still use the arguments from summary.fixest
#' summary(res, se = "hetero")
#'
#' summary(res, type = "long")
#'
#' summary(res, type = "compact")
#'
#' summary(res, type = "se_compact")
#'
#'
summary.fixest_multi = function(object, type = "short", se = NULL, cluster = NULL, dof = NULL,
                                .vcov, stage = 2, lean = FALSE, n = 1000, ...){
    dots = list(...)
    data = attr(object, "data")

    check_arg_plus(type, "match(short, long, compact, se_compact)")

    if(!missing(type) || is.null(attr(object, "print_request"))){
        attr(object, "print_request") = type
    }

    est_1 = data[[1]]
    if(is.null(est_1$cov.scaled) || !isTRUE(dots$fromPrint)){

        for(i in 1:length(data)){
            data[[i]] = summary(data[[i]], se = se, cluster = cluster, dof = dof,
                                .vcov = .vcov, stage = stage, lean = lean, n = n, ...)
        }

        # In IV: multiple estimations can be returned

        if("fixest_multi" %in% class(data[[1]])){
            meta = attr(object, "meta")

            object = setup_multi(meta$index, meta$all_names, data)
        } else {
            attr(object, "data") = data
        }

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

    # Finding out the type of SEs
    if(is_short){

        all_se = unique(unlist(sapply(data, function(x) attr(x$cov.scaled, "type"))))

        if(length(all_se) > 1){
            cat("Standard-errors: mixed (use summary() with arg. 'se' or 'cluster' to harmonize them) \n")
        } else if(length(all_se) == 1){
            cat("Standard-errors:", all_se, "\n")
        }
    }

    dict_title = c("sample" = "Sample", "lhs" = "Dep. var.", "rhs" = "Expl. vars.", "iv" = "IV", "fixef" = "Fixed-effects")

    qui_drop = apply(tree, 2, max) == 1
    if(any(qui_drop)){
        var2drop = names(tree)[qui_drop]
        for(d in var2drop){
            cat(dict_title[d], ": ", all_names[[d]][1], "\n", sep = "")
        }

        tree = tree[, !qui_drop, drop = FALSE]
        index = index[!names(index) %in% var2drop]
    }

    depth = length(index)

    headers = list()
    headers[[1]] = function(d, i) cat(dict_title[d], ": ", all_names[[d]][i], "\n", sep = "")
    headers[[2]] = function(d, i) cat("\n### ", dict_title[d], ": ", all_names[[d]][i], "\n\n", sep = "")
    headers[[3]] = function(d, i) cat("\n\n# ", toupper(dict_title[d]), ": ", all_names[[d]][i], "\n\n", sep = "")
    headers[[4]] = function(d, i) cat("\n\n#\n# ", toupper(dict_title[d]), ": ", all_names[[d]][i], "\n#\n\n", sep = "")
    headers[[5]] = function(d, i) cat("\n\n#===\n# ", toupper(dict_title[d]), ": ", all_names[[d]][i], "\n#===\n\n", sep = "")

    for(i in 1:nrow(tree)){
        for(j in 1:depth){
            d = names(index)[j]
            if(i == 1 || tree[i - 1, j] != tree[i, j]){
                headers[[depth - j + 1]](d, tree[i, j])
            }
        }

        if(is_short){
            if(isTRUE(data[[i]]$onlyFixef)){
                cat("No variable (only the fixed-effects).\n")
            } else {
                myPrintCoefTable(coeftable = coeftable(data[[i]]), show_signif = FALSE)
            }
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
#'
#' @inherit print.fixest_multi seealso
#' @inheritParams print.fixest_multi
#'
#' @param sample An integer vector, a logical scalar, or a character vector. It represents the \code{sample} identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation was a split sample. You can use \code{.N} to refer to the last element. If logical, all elements are selected in both cases, but \code{FALSE} leads \code{sample} to become the rightmost key (just try it out).
#' @param lhs An integer vector, a logical scalar, or a character vector. It represents the left-hand-sides identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation contained multiple left-hand-sides. You can use \code{.N} to refer to the last element. If logical, all elements are selected in both cases, but \code{FALSE} leads \code{lhs} to become the rightmost key (just try it out).
#' @param rhs An integer vector or a logical scalar. It represents the right-hand-sides identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation contained multiple right-hand-sides. You can use \code{.N} to refer to the last element. If logical, all elements are selected in both cases, but \code{FALSE} leads \code{rhs} to become the rightmost key (just try it out).
#' @param fixef An integer vector or a logical scalar. It represents the fixed-effects identifiers for which the results should be extracted. Only valid when the \code{fixest} estimation contained fixed-effects in a stepwise fashion. You can use \code{.N} to refer to the last element. If logical, all elements are selected in both cases, but \code{FALSE} leads \code{fixef} to become the rightmost key (just try it out).
#' @param iv An integer vector or a logical scalar. It represent the stages of the IV. Note that the length can be greater than 2 when there are multiple endogenous regressors (the first stage corresponding to multiple estimations). Note that the order of the stages depends on the \code{stage} argument from \code{\link[fixest]{summary.fixest}}. If logical, all elements are selected in both cases, but \code{FALSE} leads \code{iv} to become the rightmost key (just try it out).
#' @param i An integer vector. Represents the estimations to extract.
#' @param I An integer vector. Represents the root element to extract.
#' @param reorder Logical, default is \code{TRUE}. Indicates whether reordering of the results should be performed depending on the user input.
#' @param drop Logical, default is \code{TRUE}. If the result contains only one estimation, then if \code{drop = TRUE} it will be transformed into a \code{fixest} object (instead of \code{fixest_multi}).
#'
#' @details
#' The order with we we use the keys matter. Every time a key \code{sample}, \code{lhs}, \code{rhs}, \code{fixef} or \code{iv} is used, a reordering is performed to consider the leftmost-side key to be the new root.
#'
#' Use logical keys to easily reorder. For example, say the object \code{res} contains a multiple estimation with multiple left-hand-sides, right-hand-sides and fixed-effects. By default the results are ordered as follows: \code{lhs}, \code{fixef}, \code{rhs}. If you use \code{res[lhs = FALSE]}, then the new order is: \code{fixef}, \code{rhs}, \code{lhs}. With \code{res[rhs = TRUE, lhs = FALSE]} it becomes: \code{rhs}, \code{fixef}, \code{lhs}. In both cases you keep all estimations.
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
"[.fixest_multi" = function(x, i, sample, lhs, rhs, fixef, iv, I, reorder = TRUE, drop = TRUE){

    core_args = c("sample", "lhs", "rhs", "fixef", "iv")
    check_arg(reorder, drop, "logical scalar")

    mc = match.call()
    if(!any(c(core_args, "i", "I") %in% names(mc))){
        return(x)
    }

    use_i = "i" %in% names(mc)
    if(use_i && any(c(core_args, "I") %in% names(mc))){
        stop("The index 'i' cannot be used with any other index (", enumerate_items(c(core_args, "I"), "quote.or"), ").")
    }

    use_I = "I" %in% names(mc)
    if(use_I && any(core_args %in% names(mc))){
        stop("The index 'I' cannot be used with any other index (", enumerate_items(c("i", core_args), "quote.or"), ").")
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
    is_fixef  = !missing(fixef)
    is_iv     = !missing(iv)

    selection = list()

    last = c()

    s_max = index[["sample"]]
    if(is_sample){
        check_arg_plus(sample, "evalset logical scalar | vector(character, integer) no na", .data = list(.N = s_max))
        sample = set_index_multi(sample, s_max, all_names)

        selection$sample = (1:s_max)[sample]

    } else if("sample" %in% names(index)){
        selection$sample = 1:s_max
    }

    lhs_max = index[["lhs"]]
    if(is_lhs){
        check_arg_plus(lhs, "evalset logical scalar | vector(character, integer) no na", .data = list(.N = lhs_max))
        lhs = set_index_multi(lhs, lhs_max, all_names)

        selection$lhs = (1:lhs_max)[lhs]
    } else if("lhs" %in% names(index)){
        selection$lhs = 1:lhs_max
    }

    rhs_max = index[["rhs"]]
    if(is_rhs){
        check_arg_plus(rhs, "evalset logical scalar | vector(character, integer) no na", .data = list(.N = rhs_max))
        rhs = set_index_multi(rhs, rhs_max, all_names)

        selection$rhs = (1:rhs_max)[rhs]
    } else if("rhs" %in% names(index)){
        selection$rhs = 1:rhs_max
    }

    fixef_max = index[["fixef"]]
    if(is_fixef){
        check_arg_plus(fixef, "evalset logical scalar | vector(character, integer) no na", .data = list(.N = fixef_max))
        fixef = set_index_multi(fixef, fixef_max, all_names)

        selection$fixef = (1:fixef_max)[fixef]
    } else if("fixef" %in% names(index)){
        selection$fixef = 1:fixef_max
    }

    iv_max = index[["iv"]]
    if(is_iv){
        check_arg_plus(iv, "evalset logical scalar | vector(character, integer) no na", .data = list(.N = iv_max))
        iv = set_index_multi(iv, iv_max, all_names)

        selection$iv = (1:iv_max)[iv]
    } else if("iv" %in% names(index)){
        selection$iv = 1:iv_max
    }

    # We keep the order of the user!!!!!
    sc = sys.call()
    user_order = names(sc)[-(1:2)]
    if(reorder == FALSE){
        user_order = names(index)
    } else {
        user_order = c(user_order, setdiff(names(index), user_order))
        if(length(last) > 0){
            user_order = c(setdiff(user_order, last), last)
        }
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

        # 2) we reorder it according to the order of the user
        new_tree = list()
        dim_order = selection[[my_dim]]
        for(i in seq_along(dim_order)){
            new_tree[[i]] = tree[tree[[my_dim]] == dim_order[i], ]
        }
        tree = do.call(base::rbind, new_tree)
    }

    new_data = list()
    for(i in 1:nrow(tree)){
        new_data[[i]] = data[[tree$obs[i]]]
    }

    if(length(new_data) == 1 && drop){
        return(new_data[[1]])
    }

    res_multi = setup_multi(new_index, new_all_names, new_data, simplify = FALSE)

    return(res_multi)
}


#' Extracts the root of a fixest_multi object
#'
#' Extracts an element at the root of a fixest_multi object.
#'
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


set_index_multi = function(x, vmax, all_names){
    # Function specific to [.fixest_multi => global assignments!!!
    arg = deparse(substitute(x))

    if(is.logical(x)){
        if(isFALSE(x)){
            last = get("last", parent.frame())
            last[length(last) + 1] = arg
            assign("last", last, parent.frame())
        }
        x = 1:vmax
    } else if(is.character(x)){
        dict = 1:vmax
        names(dict) = all_names[[arg]]
        vars = keep_apply(all_names[[arg]], x)
        vars = order_apply(vars, x)

        x = as.integer(dict[vars])

        if(length(x) == 0){
            stop_up("The set of regular expressions in '", arg, "' didn't match any choice.")
        }

    } else if(any(abs(x) > vmax)){
        stop_up("The index '", arg, "' cannot be greater than ", vmax, ". Currently ", x[which.max(abs(x))], " is not valid.")
    }

    x
}



#' Extracts the coefficients of fixest_multi objects
#'
#' Utility to extract the coefficients of multiple estimations and rearrange them into a matrix.
#'
#' @inheritParams etable
#'
#' @param object A \code{fixest_multi} object. Obtained from a multiple estimation.
#' @param ... Not currently used.
#'
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # A multiple estimation
#' est = feols(y ~ x1 + csw0(x2, x3), base)
#'
#' # Getting all the coefficients at once,
#' # each row is a model
#' coef(est)
#'
#' # Example of keep/drop/order
#' coef(est, keep = "Int|x1", order = "x1")
#'
#'
#' # To change the order of the model, use fixest_multi
#' # extraction tools:
#' coef(est[rhs = .N:1])
#'
#'
coef.fixest_multi = function(object, keep, drop, order, ...){
    # row: model
    # col: coefficient

    check_arg(keep, drop, order, "NULL character vector no na")

    data = attr(object, "data")

    res_list = list()
    for(i in seq_along(data)){
        res_list[[i]] = coef(data[[i]])
    }

    all_names = unique(unlist(lapply(res_list, names)))

    if(!missnull(keep)) all_names = keep_apply(all_names, keep)
    if(!missnull(drop)) all_names = drop_apply(all_names, drop)
    if(!missnull(order)) all_names = order_apply(all_names, order)

    if(length(all_names) == 0) return(NULL)

    nr = length(res_list)
    nc = length(all_names)

    res_list = lapply(res_list, function(x) x[all_names])
    res = do.call(rbind, res_list)
    colnames(res) = all_names

    res
}

#' @rdname coef.fixest_multi
coefficients.fixest_multi <- coef.fixest_multi


#' Extracts the residuals from a \code{fixest_multi} object
#'
#' Utility to extract the residuals from multiple \code{fixest} estimations. If possible, all the residuals are coerced into
#'
#' @inheritParams resid.fixest
#'
#' @param object A \code{fixes_multi} object.
#' @param na.rm Logical, default is \code{FALSE}. Should the NAs be kept? If \code{TRUE}, they are removed.
#'
#' @return
#' If all the models return residuals of the same length, a matrix is returned. Otherwise, a \code{data.frame} is returned.
#'
#' @examples
#'
#' base = iris
#' names(base) = c("y", "x1", "x2", "x3", "species")
#'
#' # A multiple estimation
#' est = feols(y ~ x1 + csw0(x2, x3), base)
#'
#' # We can get all the residuals at once,
#' # each column is a model
#' head(resid(est))
#'
#' # We can select/order the model using fixest_multi extraction
#' head(resid(est[rhs = .N:1]))
#'
resid.fixest_multi = function(object, type = c("response", "deviance", "pearson", "working"), na.rm = FALSE, ...){

    # Je fais un prototype pour le moment, je l'ameliorerai apres (07-04-2021)

    check_arg_plus(type, "match")
    check_arg_plus(na.rm, "logical scalar")

    data = attr(object, "data")

    res_list = list()
    for(i in seq_along(data)){
        res_list[[i]] = resid(data[[i]], type = type, na.rm = na.rm)
    }

    n_all = sapply(res_list, length)

    if(all(n_all == n_all[1])){
        res = do.call(cbind, res_list)
    } else {
        res = res_list
    }

    res
}


#' @rdname resid.fixest_multi
residuals.fixest_multi <- resid.fixest_multi


