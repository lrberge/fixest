#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Nov 07 09:05:26 2020
# ~: fixest_multi
#----------------------------------------------#


setup_multi = function(index, all_names, data){
    # Basic setup function

    # We drop non-used dimensions
    n_all = sapply(index, max)
    qui = n_all == 1
    if(any(qui)){
        index = index[!qui]
        all_names = all_names[!qui]
    }

    meta = list(index = index, all_names = all_names)
    tree = as.data.frame(do.call("expand.grid", rev(lapply(index, seq))))
    meta$tree = tree[, rev(1:ncol(tree)), drop = FALSE]
    meta$tree$obs = 1:nrow(tree)

    res_multi = list()
    first_dim = all_names[[1]]
    for(i in seq_along(first_dim)){
        res_multi[[first_dim[i]]] = 1
    }

    attr(res_multi, "meta") = meta
    attr(res_multi, "data") = data

    class(res_multi) = "fixest_multi"

    return(res_multi)
}

summary.fixest_multi = function(object, ...){
    dots = list(...)
    data = attr(object, "data")

    est_1 = data[[1]]
    if(!is.null(est_1$cov.scaled) && isTRUE(dots$fromPrint)) return(object)

    for(i in 1:length(data)){
        data[[i]] = summary(data[[i]], nframes_up = 2, ...)
    }

    attr(object, "data") = data

    return(object)
}


print.fixest_multi = function(x, n, type = getFixest_print.type(), ...){

    x = summary(x, fromPrint = TRUE, ...)

    meta = attr(x, "meta")
    data = attr(x, "data")

    index = meta$index
    all_names = meta$all_names
    tree = meta$tree

    cat(signif_plus(do.call(prod, index)), " fixest estimations\n", sep = "")

    # Finding out the type of SEs
    est_1 = data[[1]]
    coeftable = est_1$coeftable
    se.type = attr(coeftable, "type")
    cat("Standard-errors:", se.type, "\n")

    if(type == "compact"){
        res = data.frame(i = tree$obs)
        for(my_dim in names(index)){
            res[[my_dim]] = sfill(as.character(factor(tree[[my_dim]], labels = all_names[[my_dim]])), right = TRUE)
        }
        res$i = NULL

        signifCode = c(" ***"=0.001, " **"=0.01, " *"=0.05, " ."=0.1)

        ct_all = list()
        for(i in seq_along(data)){
            ct = data[[i]]$coeftable
            vname = row.names(ct)

            stars = cut(ct[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
            stars[is.na(stars)] = ""

            value = paste0(signif_plus(ct[, 1], 3), " (", signif_plus(ct[, 2], 3), ")", stars)
            names(value) = vname

            ct_all[[i]] = value
        }

        vname_all = unique(unlist(lapply(ct_all, names)))

        tmp = lapply(ct_all, function(x) x[vname_all])
        my_ct = do.call("rbind", tmp)
        colnames(my_ct) = vname_all
        my_ct[is.na(my_ct)] = ""

        for(i in seq_along(vname_all)){
            res[[vname_all[i]]] = sfill(my_ct[, i], anchor = "(", right = TRUE)
        }

        return(res)
    }

    depth = length(index)

    dict_title = c("sample" = "Sample", "lhs" = "Dep. var.", "rhs" = "Expl. vars.")

    if(depth == 3){

        for(i in 1:index[[1]]){
            cat("\n\n# ", toupper(dict_title[names(index)[1]]), ": ", all_names[[1]][i], "\n\n", sep = "")

            for(j in 1:index[[2]]){
                cat("\n### ", dict_title[names(index)[2]], ": ", all_names[[2]][j], "\n\n", sep = "")

                for(k in 1:index[[3]]){
                    cat(dict_title[names(index)[3]], ": ", all_names[[3]][k], "\n", sep = "")

                    obs = tree$obs[tree[, 1] == i & tree[, 2] == j & tree[, 3] == k]
                    myPrintCoefTable(coeftable = coeftable(data[[obs]], ...), show_signif = FALSE)
                    if(k != index[[3]]) cat("---\n")
                }
            }
        }

    } else if(depth == 2){
        for(i in 1:index[[1]]){
            cat("\n### ", dict_title[names(index)[1]], ": ", all_names[[1]][i], "\n\n", sep = "")

            for(j in 1:index[[2]]){
                cat(dict_title[names(index)[2]], ": ", all_names[[2]][j], "\n", sep = "")

                obs = tree$obs[tree[, 1] == i & tree[, 2] == j]
                myPrintCoefTable(coeftable = coeftable(data[[obs]], ...), show_signif = FALSE)
                if(j != index[[2]]) cat("---\n")
            }
        }
    } else {
        for(i in 1:index[[1]]){
            cat(dict_title[names(index)[1]], ": ", all_names[[1]][i], "\n", sep = "")

            obs = tree$obs[tree[, 1] == i]
            myPrintCoefTable(coeftable = coeftable(data[[obs]], ...), show_signif = FALSE)
            if(i != index[[1]]) cat("---\n")
        }
    }

}

"[[.fixest_multi" = function(x, i){

    data = attr(x, "data")

    check_arg_plus(i, "evalset integer scalar mbt", .data = list(.N = length(data)))

    return(data[[i]])
}


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
            stop("The index 'I' refers to root elements, and thus cannot be greater than ", I_max, ". Currently ", I[which.max(abs(I))], " is not valid.")
        }

        obs = (1:I_max)[I]

        # If only one dimension => result is simplified
        if(length(index) == 1 && length(obs) == 1){
            res = data[[obs]]

        } else if(length(obs) == 1){
            # The first dimension is dissolved
            obs_keep = tree[tree[, 1] == obs, "obs"]
            res = setup_multi(index[-1], all_names[-1], data[obs_keep])

        } else {
            new_index = index
            new_index[[1]] = length(obs)

            new_all_names = all_names
            new_all_names[[1]] = all_names[[1]][obs]

            obs_keep = tree[tree[, 1] %in% obs, "obs"]
            res = setup_multi(new_index, new_all_names, data[obs_keep])
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


"$.fixest_multi" = function(x, name){

    # real variables
    if(!name %in% names(x)){
        mc = match.call()
        stop(name, " is not at the root of the multiple fixest estimations object.")
    }

    qui = which(name == names(x))

    return(x[I = qui])
}


















