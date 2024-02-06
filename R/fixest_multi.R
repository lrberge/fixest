#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Sat Nov 07 09:05:26 2020
# ~: fixest_multi
#----------------------------------------------#

setup_multi = function(data, values, var = NULL, tree = NULL){
  # the incoming data is ALWAYS strongly structured
  # => they all have the same number of elements
  # data:
  # either a list of fixest objects
  # either a list of fixest_multi objects
  #
  # values: must be strongly and properly formatted
  # its length is the nber of objects (length(data)), with the appropriate names
  # var: to keep the information on the variable (sample, subset)

  # We also add the $model_info variable in each model

  # To remove after development
  check_arg(data, "list")
  check_arg(values, "named list")
  check_arg(var, "NULL character vector no na")
  check_arg(tree, "NULL data.frame")

  n_models = length(data)
  IS_VAR = !is.null(var)
  IS_TREE = !is.null(tree)

  if(!IS_TREE){
    stopifnot(identical(class(data), "list"))
  }

  var_label = NULL
  if(IS_VAR){
    stopifnot(length(values) == 1)
    var_label = names(values)[1]
    if(length(var) == 1){
      var = rep(var, n_models)
    }
  }

  IS_NESTED = inherits(data[[1]], "fixest_multi")

  if(IS_TREE){
    # This is an internal call from [.fixest_multi
    # data = the final data
    # values = the new tree

    res = data
    tree$id = NULL # we re-create it later

    if(IS_NESTED){

      # We allow non balanced data lists
      res = vector("list", sum(lengths(data)))
      tree_left = list()
      tree_right = list()
      index = 1
      for(i in 1:n_models){
        data_i = data[[i]]

        # updating the tree
        tree_nested = attr(data_i, "tree")[, -1, drop = FALSE]
        n_nested = nrow(tree_nested)
        tree_left[[i]] = rep_df(tree[i, , drop = FALSE], each = n_nested)
        tree_right[[i]] = tree_nested

        for(j in 1:n_nested){
          res[[index]] = data_i[[j]]
          index = index + 1
        }
      }

      tree_left = do.call(rbind, tree_left)
      tree_right = do.call(rbind, tree_right)

      tree = cbind(tree_left, tree_right)
    }

  } else {

    v_names = names(values)
    tree = as.data.frame(values)

    if(IS_NESTED){
      # bookkeeping needed: note that we ensure beforehand that each element is STRONGLY consistent

      tree_left = list()
      tree_right = list()
      for(i in 1:n_models){
        # updating the tree
        tree_nested = attr(data[[i]], "tree")[, -1, drop = FALSE]
        n_nested = nrow(tree_nested)
        tree_left[[i]] = rep_df(tree[i, , drop = FALSE], each = n_nested)
        tree_right[[i]] = tree_nested
      }

      tree_left = do.call(rbind, tree_left)
      tree_right = do.call(rbind, tree_right)

      tree = cbind(tree_left, tree_right)
    } else if(!inherits(data[[1]], "fixest")){
      stop("Internal error: the current data type is not supportded by setup_multi.")
    }

    # res: a plain list containing all the models
    res = vector("list", nrow(tree))
    index = 1
    for(i in 1:n_models){
      data_i = data[[i]]

      # new model information
      new_info = list()
      for(v in v_names){
        if(IS_VAR){
          new_info[[v]] = list(var = var[i], value = values[[v]][i])
        } else {
          new_info[[v]] = values[[v]][i]
        }
      }

      n_j = if(IS_NESTED) length(data_i) else 1
      for(j in 1:n_j){

        if(IS_NESTED){
          mod = data_i[[j]]
        } else {
          mod = data_i
        }

        # updating the model information
        model_info = mod$model_info
        for(v in names(new_info)){
          model_info[[v]] = new_info[[v]]
        }

        mod$model_info = model_info
        res[[index]] = mod

        index = index + 1
      }
    }
  }

  if(IS_VAR){
    tree = cbind(var, tree)
    names(tree)[1] = paste0(var_label, ".var")
  }

  tree_names = mapply(function(x, y) paste0(x, ": ", y), names(tree), tree)
  if(is.vector(tree_names)){
    model_names = paste(tree_names, collapse = "; ")
  } else {
    tree_names = as.data.frame(tree_names)
    model_names = apply(tree_names, 1, paste0, collapse = "; ")
  }

  # indexes
  info = index_from_tree(tree)
  index_names = info$index_names
  tree_index = info$tree_index

  tree = cbind(id = 1:nrow(tree), tree)

  # Shouldn't I remove tree_index and index_names since they can be built from the tree?
  # It seems it can be useful if they're directly computed... We'll see.
  names(res) = model_names
  class(res) = "fixest_multi"
  attr(res, "tree") = tree
  attr(res, "tree_index") = tree_index
  attr(res, "index_names") = index_names

  res
}

index_from_tree = function(tree){
  index_names = list()
  tree_index = list()
  names_keep = names(tree)[!grepl("\\.var$|^id$", names(tree))]
  for(v in names_keep){
    z = tree[[v]]
    fact = factor(z, levels = unique(z))
    index_names[[v]] = levels(fact)
    tree_index[[v]] = as.integer(unclass(fact))
  }
  tree_index = as.data.frame(tree_index)

  list(index_names = index_names, tree_index = tree_index)
}

reshape_multi = function(x, obs, colorder = NULL){
  # x: fixest_multi object
  # obs: indexes to keep

   tree = attr(x, "tree")
   new_tree = tree[obs, , drop = FALSE]

   if(!is.null(colorder)){
     new_tree_list = list()
     for(i in seq_along(colorder)){
       # I use grep to catch ".var" variables
       qui = grepl(colorder[i], names(tree))
       new_tree_list[[i]] = new_tree[, qui, drop = FALSE]
     }
     new_tree_list[[i + 1]] = new_tree["id"]

     new_tree = do.call(cbind, new_tree_list)

   }

   n_models = nrow(new_tree)
   new_data = vector("list", n_models)
   for(i in 1:n_models){
     new_data[[i]] = x[[new_tree$id[i]]]
   }

  setup_multi(new_data, tree = new_tree)
}


set_index_multi = function(x, index_names){
  # Function specific to [.fixest_multi => global assignments!!!
  arg = deparse(substitute(x))

  NAMES = index_names[[arg]]
  vmax = length(NAMES)

  if(is.logical(x)){
    if(isFALSE(x)){
      last = get("last", parent.frame())
      last[length(last) + 1] = arg
      assign("last", last, parent.frame())
    }
    res = 1:vmax
  } else if(is.character(x)){
    dict = 1:vmax
    names(dict) = NAMES
    vars = keep_apply(NAMES, x)
    vars = order_apply(vars, x)

    res = as.integer(dict[vars])

    if(length(res) == 0){
      stop_up("The set of regular expressions (equal to: ", x, ") in '", arg, 
              "' didn't match any choice.")
    }

  } else if(any(abs(x) > vmax)){
    stop_up("The index '", arg, "' cannot be greater than ", vmax, 
            ". Currently ", x[which.max(abs(x))], " is not valid.")
  } else {
    res = x
  }

  res
}




rep_df = function(x, times = 1, each = 1, ...){
  if(identical(times, 1) && identical(each, 1)){
    return(x)
  }

  as.data.frame(lapply(x, rep, times = times, each = each))
}

####
#### USER LEVEL ####
####

#' Extracts the models tree from a `fixest_multi` object
#'
#' Extracts the meta information on all the models contained in a `fixest_multi` estimation.
#'
#' @inheritParams print.fixest_multi
#' @param simplify Logical, default is `FALSE`. The default behavior is to display all the meta 
#' information, even if they are identical across models. By using `simplify = TRUE`, only the 
#' information with some variation is kept.
#'
#' @return
#' It returns a `data.frame` whose first column (named `id`) is the index of the models and 
#' the other columns contain the information specific to each model (e.g. which sample, 
#' which RHS,  which dependent variable, etc).
#'
#' @examples
#'
#' # a multiple estimation
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' est = feols(y ~ csw(x.[, 1:3]), base, fsplit = ~species)
#'
#' # All the meta information
#' models(est)
#'
#' # Illustration: Why use simplify
#' est_sub = est[sample = 2]
#' models(est_sub)
#' models(est_sub, simplify = TRUE)
#'
#'
#'
models = function(x, simplify = FALSE){
  check_arg(x, "class(fixest_multi)")

  res = attr(x, "tree")
  if(simplify){
    who_keep = sapply(res, function(x) length(unique(x)) != 1)

    if(!all(who_keep)){
      # we need to handle the behavior with the .var thing
      names_keep = names(res)[who_keep]
      pattern = sma("^({'|'c?names_keep})")

      res = res[, grepl(pattern, names(res)), drop = FALSE]
    }

  }

  res
}


####
#### METHODS ####
####



#' Summary for fixest_multi objects
#'
#' Summary information for fixest_multi objects. In particular, this is used to specify the 
#' type of standard-errors to be computed.
#'
#' @method summary fixest_multi
#'
#' @inheritParams summary.fixest
#'
#' @inherit print.fixest_multi seealso
#'
#' @param object A `fixest_multi` object, obtained from a `fixest` estimation leading to 
#' multiple results.
#' @param type A character either equal to `"short"`, `"long"`, `"compact"`, `"se_compact"` 
#' or `"se_long"`. If `short`, only the table of coefficients is displayed for each estimation. 
#' If `long`, then the full results are displayed for each estimation. If `compact`, 
#' a `data.frame` is returned with one line per model and the formatted 
#' coefficients + standard-errors in the columns. If `se_compact`, a `data.frame` is 
#' returned with one line per model, one numeric column for each coefficient and one numeric 
#' column for each standard-error. If `"se_long"`, same as `"se_compact"` but the data is in a 
#' long format instead of wide.
#' @param ... Not currently used.
#'
#' @return
#' It returns either an object of class `fixest_multi` (if `type` equals `short` or `long`), 
#' either a `data.frame` (if type equals `compact` or `se_compact`).
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
#' summary(res, type = "se_long")
#'
#'
summary.fixest_multi = function(object, type = "short", vcov = NULL, se = NULL, 
                                cluster = NULL, ssc = NULL,
                                .vcov = NULL, stage = 2, lean = FALSE, n = 1000, ...){
  dots = list(...)

  check_set_arg(type, "match(short, long, compact, se_compact, se_long)")

  if(!missing(type) || is.null(attr(object, "print_request"))){
    attr(object, "print_request") = type
  }

  if(is_user_level_call()){
    validate_dots(suggest_args = c("type", "vcov"),
            valid_args = c("agg", "forceCovariance", "keepBounded", "nthreads"))
  }

  est_1 = object[[1]]
  if(is.null(est_1$cov.scaled) || !isTRUE(dots$fromPrint)){

    for(i in 1:length(object)){
      object[[i]] = summary(object[[i]], vcov = vcov, se = se, cluster = cluster, ssc = ssc,
                            .vcov = .vcov, stage = stage, lean = lean, n = n, ...)
    }

    # In IV: multiple estimations can be returned

    if("fixest_multi" %in% class(object[[1]])){
      tree = attr(object, "tree")
      object = setup_multi(object, tree = tree)
    }

  }

  if(type %in% c("compact", "se_compact", "se_long")){
    tree = attr(object, "tree")
    tree_index = attr(object, "tree_index")

    res = data.frame(i = tree$id)

    if(!"lhs" %in% names(tree_index)){
      res$lhs = sapply(object, function(x) as.character(x$fml[[2]]))
    }

    for(my_dim in names(tree_index)){
      res[[my_dim]] = sfill(tree[[my_dim]], right = TRUE)
    }
    res$i = NULL

    if(type == "se_long"){
      res$type = "coef"
    }

    n_start = ncol(res)

    signifCode = c("***"=0.001, "**"=0.01, "*"=0.05, "."=0.1)

    ct_all = list()
    for(i in seq_along(object)){
      ct = object[[i]]$coeftable
      vname = row.names(ct)

      if(type == "compact"){
        stars = cut(ct[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
        stars[is.na(stars)] = ""

        value = paste0(format_number(ct[, 1], 3), stars, " (", format_number(ct[, 2], 3), ")")
        names(value) = vname

      } else if(type %in% c("se_compact", "se_long")){
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

    if(type == "se_long"){
      # clumsy... but works
      who_se = which(grepl("__se", names(res)))
      se_all = res[, c(1:n_start, who_se)]
      se_all$type = "se"
      names(se_all) = gsub("__se$", "", names(se_all))
      coef_all = res[, -who_se]

      quoi = rbind(coef_all, se_all)
      n = nrow(coef_all)
      res = quoi[rep(1:n, each = 2) + rep(c(0, n), n), ]
      row.names(res) = NULL
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
#' @param x A `fixest_multi` object, obtained from a `fixest` estimation leading to 
#' multiple results.
#' @param ... Other arguments to be passed to [`summary.fixest_multi`].
#'
#' @seealso
#' The main fixest estimation functions: [`feols`], [`fepois`][fixest::feglm], 
#' [`fenegbin`][fixest::femlm], [`feglm`], [`feNmlm`]. Tools for mutliple fixest 
#' estimations: [`summary.fixest_multi`], [`print.fixest_multi`], [`as.list.fixest_multi`], 
#' \code{\link[fixest]{sub-sub-.fixest_multi}}, \code{\link[fixest]{sub-.fixest_multi}}.
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

  if(is_user_level_call()){
    validate_dots(valid_args = dsb("/type, vcov, se, cluster, ssc, stage, lean, agg, forceCovariance, keepBounded, n, nthreads"))
  }

  x = summary(x, fromPrint = TRUE, ...)

  # Type = compact
  if(is.data.frame(x)){
    return(x)
  }

  is_short = identical(attr(x, "print_request"), "short")

  tree = attr(x, "tree")
  tree_index = attr(x, "tree_index")

  # Finding out the type of SEs
  if(is_short){

    all_se = unique(unlist(sapply(x, function(x) attr(x$cov.scaled, "type"))))

    if(length(all_se) > 1){
      cat("Standard-errors: mixed (use summary() with arg. 'vcov' to harmonize them) \n")
    } else if(length(all_se) == 1){
      cat("Standard-errors:", all_se, "\n")
    }
  }

  dict_title = c("sample" = "Sample", "lhs" = "Dep. var.", "rhs" = "Expl. vars.",
                 "iv" = "IV", "fixef" = "Fixed-effects", sample.var = "Sample var.")

  qui_drop = apply(tree_index, 2, max) == 1

  if(any(qui_drop) && !all(qui_drop)){
    var2drop = names(tree_index)[qui_drop]
    for(d in var2drop){
      cat(dict_title[d], ": ", tree[[d]][1], "\n", sep = "")
    }

    tree = tree[, !names(tree) %in% var2drop, drop = FALSE]
    tree_index = tree_index[, !qui_drop, drop = FALSE]
  }

  depth = ncol(tree_index)

  headers = list()
  headers[[1]] = function(d, i) cat(dict_title[d], ": ", tree[[d]][i], "\n", sep = "")
  headers[[2]] = function(d, i) cat("\n### ", dict_title[d], ": ", tree[[d]][i], "\n\n", sep = "")
  headers[[3]] = function(d, i) cat("\n\n# ", toupper(dict_title[d]), ": ", tree[[d]][i], "\n\n", sep = "")
  headers[[4]] = function(d, i) cat("\n\n#\n# ", toupper(dict_title[d]), ": ", tree[[d]][i], "\n#\n\n", sep = "")
  headers[[5]] = function(d, i) cat("\n\n#===\n# ", toupper(dict_title[d]), ": ", tree[[d]][i], "\n#===\n\n", sep = "")

  for(i in 1:nrow(tree)){
    for(j in 1:depth){
      d = names(tree_index)[j]
      if(i == 1 || tree_index[i - 1, j] != tree_index[i, j]){
        headers[[depth - j + 1]](d, i)
      }
    }

    if(is_short){
      if(isTRUE(x[[i]]$onlyFixef)){
        cat("No variable (only the fixed-effects).\n")
      } else {
        print_coeftable(coeftable = coeftable(x[[i]]), show_signif = FALSE)
      }
      if(tree_index[i, depth] != max(tree_index[, depth])) cat("---\n")
    } else {
      print(x[[i]])
      if(tree_index[i, depth] != max(tree_index[, depth])) cat("\n")
    }

  }

}

#' Extracts one element from a `fixest_multi` object
#'
#' Extracts single elements from multiple `fixest` estimations.
#'
#' @inherit print.fixest_multi seealso
#' @inheritParams print.fixest_multi
#'
#' @param i An integer scalar. The identifier of the estimation to extract.
#'
#' @return
#' A `fixest` object is returned.
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

  n = length(x)
  check_set_arg(i, "evalset integer scalar mbt", .data = list(.N = n))
  if(i < 0 || i > length(x)){
    stop("Index 'i' must lie within [1; ", n, "]. Problem: it is equal to ", i, ".")
  }

  `[[.data.frame`(x, i)
}


#' Subsets a fixest_multi object
#'
#' Subsets a fixest_multi object using different keys.
#'
#'
#' @inherit print.fixest_multi seealso
#' @inheritParams print.fixest_multi
#'
#' @param sample An integer vector, a logical scalar, or a character vector. It represents 
#' the `sample` identifiers for which the results should be extracted. Only valid when the 
#' `fixest` estimation was a split sample. You can use `.N` to refer to the last element. 
#' If logical, all elements are selected in both cases, but `FALSE` leads `sample` to become 
#' the rightmost key (just try it out).
#' @param lhs An integer vector, a logical scalar, or a character vector. It represents 
#' the left-hand-sides identifiers for which the results should be extracted. Only valid when 
#' the `fixest` estimation contained multiple left-hand-sides. You can use `.N` to refer to 
#' the last element. If logical, all elements are selected in both cases, but `FALSE` 
#' leads `lhs` to become the rightmost key (just try it out).
#' @param rhs An integer vector or a logical scalar. It represents the right-hand-sides 
#' identifiers for which the results should be extracted. Only valid when the `fixest` 
#' estimation contained multiple right-hand-sides. You can use `.N` to refer to the last 
#' element. If logical, all elements are selected in both cases, but `FALSE` leads `rhs` to 
#' become the rightmost key (just try it out).
#' @param fixef An integer vector or a logical scalar. It represents the fixed-effects 
#' identifiers for which the results should be extracted. Only valid when the `fixest` 
#' estimation contained fixed-effects in a stepwise fashion. You can use `.N` to refer to the 
#' last element. If logical, all elements are selected in both cases, but `FALSE` leads `fixef` 
#' to become the rightmost key (just try it out).
#' @param iv An integer vector or a logical scalar. It represent the stages of the IV. Note 
#' that the length can be greater than 2 when there are multiple endogenous regressors (the 
#' first stage corresponding to multiple estimations). Note that the order of the stages depends 
#' on the `stage` argument from [`summary.fixest`]. If logical, all elements are selected in 
#' both cases, but `FALSE` leads `iv` to become the rightmost key (just try it out).
#' @param i An integer vector. Represents the estimations to extract.
#' @param I An integer vector. Represents the root element to extract.
#' @param reorder Logical, default is `TRUE`. Indicates whether reordering of the results 
#' should be performed depending on the user input.
#' @param drop Logical, default is `FALSE`. If the result contains only one estimation, 
#' then if `drop = TRUE` it will be transformed into a `fixest` object (instead of `fixest_multi`).
#'
#' @details
#' The order with we we use the keys matter. Every time a key `sample`, `lhs`, `rhs`, 
#' `fixef` or `iv` is used, a reordering is performed to consider the leftmost-side key 
#' to be the new root.
#'
#' Use logical keys to easily reorder. For example, say the object `res` contains a 
#' multiple estimation with multiple left-hand-sides, right-hand-sides and fixed-effects. 
#' By default the results are ordered as follows: `lhs`, `fixef`, `rhs`. 
#' If you use `res[lhs = FALSE]`, then the new order is: `fixef`, `rhs`, `lhs`. 
#' With `res[rhs = TRUE, lhs = FALSE]` it becomes: `rhs`, `fixef`, `lhs`. In both cases 
#' you keep all estimations.
#'
#' @return
#' It returns a `fixest_multi` object. If there is only one estimation left in the object, then 
#' the result is simplified into a `fixest` object.
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
"[.fixest_multi" = function(x, i, sample, lhs, rhs, fixef, iv, I, reorder = TRUE, drop = FALSE){

  core_args = c("sample", "lhs", "rhs", "fixef", "iv")
  check_arg(reorder, drop, "logical scalar")
  extra_args = c("reorder", "drop")

  mc = match.call()
  if(!any(c(core_args, "i", "I") %in% names(mc))){
    return(x)
  }

  use_i = "i" %in% names(mc)
  if(use_i && any(c(core_args, "I") %in% names(mc))){
    stopi("The index 'i' cannot be used with any other index ({enum.q.or ? c(core_args, 'I')}).")
  }

  use_I = "I" %in% names(mc)
  if(use_I && any(core_args %in% names(mc))){
    stopi("The index 'I' cannot be used with any other index ({enum.q.or ? c('i', core_args)}).")
  }

  # We get the meta information
  tree = attr(x, "tree")
  tree_index = attr(x, "tree_index")
  index_names = attr(x, "index_names")
  index_n = lapply(index_names, length)

  # tree_index does not contain extra info like id or .var
  args = c(names(tree_index), extra_args)

  nc = ncol(tree)
  n = nrow(tree)

  if(!use_i && !use_I){
    pblm = setdiff(names(mc)[-(1:2)], args)
    if(length(pblm) > 0){
      stopi("The ind{$(ex;ices), enum.bq, is ? pblm} not valid for this list of results (the valid one{s, are, enum.bq ? index_n}).")
    }
  }

  if(use_i){
    check_set_arg(i, "evalset integer vector l0 no na", .data = list(.N = n))

    if(length(i) == 0) return(list())

    if(any(abs(i) > n)){
      stop("The index 'i' cannot have values greater than ", n, ". Currently ", i[which.max(abs(i))], " is not valid.")
    }

    obs = (1:n)[i]

    res = reshape_multi(x, obs)

    return(res)
  }

  if(use_I){
    I_max = index_n[[1]]
    check_set_arg(I, "evalset integer vector no na", .data = list(.N = I_max))

    if(any(abs(I) > I_max)){
      stop("The index 'I' refers to root elements (here ", names(index_n)[1], "), and thus cannot be greater than ", I_max, ". Currently ", I[which.max(abs(I))], " is not valid.")
    }

    obs = (1:I_max)[I]

    tree_index$obs = 1:nrow(tree_index)
    new_tree = list()
    for(i in seq_along(obs)){
      new_tree[[i]] = tree_index[tree_index[[1]] == obs[i], ]
    }
    tree_index = do.call(base::rbind, new_tree)

    res = reshape_multi(x, tree_index$obs)

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

  s_max = index_n[["sample"]]
  if(is_sample){
    check_set_arg(sample, "evalset logical scalar | vector(character, integer) no na", 
                  .data = list(.N = s_max))
    sample = set_index_multi(sample, index_names)

    selection$sample = (1:s_max)[sample]

  } else if("sample" %in% names(index_n)){
    selection$sample = 1:s_max
  }

  lhs_max = index_n[["lhs"]]
  if(is_lhs){
    check_set_arg(lhs, "evalset logical scalar | vector(character, integer) no na", 
                  .data = list(.N = lhs_max))
    lhs = set_index_multi(lhs, index_names)

    selection$lhs = (1:lhs_max)[lhs]
  } else if("lhs" %in% names(index_n)){
    selection$lhs = 1:lhs_max
  }

  rhs_max = index_n[["rhs"]]
  if(is_rhs){
    check_set_arg(rhs, "evalset logical scalar | vector(character, integer) no na", 
                  .data = list(.N = rhs_max))
    rhs = set_index_multi(rhs, index_names)

    selection$rhs = (1:rhs_max)[rhs]
  } else if("rhs" %in% names(index_n)){
    selection$rhs = 1:rhs_max
  }

  fixef_max = index_n[["fixef"]]
  if(is_fixef){
    check_set_arg(fixef, "evalset logical scalar | vector(character, integer) no na", 
                  .data = list(.N = fixef_max))
    fixef = set_index_multi(fixef, index_names)

    selection$fixef = (1:fixef_max)[fixef]
  } else if("fixef" %in% names(index_n)){
    selection$fixef = 1:fixef_max
  }

  iv_max = index_n[["iv"]]
  if(is_iv){
    check_set_arg(iv, "evalset logical scalar | vector(character, integer) no na", 
                  .data = list(.N = iv_max))
    iv = set_index_multi(iv, index_names)

    selection$iv = (1:iv_max)[iv]
  } else if("iv" %in% names(index_n)){
    selection$iv = 1:iv_max
  }

  # We keep the order of the user!!!!!
  sc = sys.call()
  user_order = setdiff(names(sc)[-(1:2)], extra_args)
  if(reorder == FALSE){
    user_order = names(index_n)
  } else {
    user_order = c(user_order, setdiff(names(index_n), user_order))
    if(length(last) > 0){
      user_order = c(setdiff(user_order, last), last)
    }
  }

  tree_index$obs = 1:nrow(tree_index)
  for(my_dim in rev(user_order)){
    # 1) we prune the tree
    obs_keep = tree_index[[my_dim]] %in% selection[[my_dim]]

    if(!any(obs_keep)){
      stop("No models ended up selected: revise selection?")
    }

    tree_index = tree_index[obs_keep, , drop = FALSE]

    # 2) we reorder it according to the order of the user
    new_tree = list()
    dim_order = selection[[my_dim]]
    for(i in seq_along(dim_order)){
      new_tree[[i]] = tree_index[tree_index[[my_dim]] == dim_order[i], , drop = FALSE]
    }
    tree_index = do.call(base::rbind, new_tree)
  }

  n_models = nrow(tree_index)
  if(n_models == 1 && drop){
    return(x[[tree_index$obs]])
  }

  # Reshaping a fixest_multi object properly
  res = reshape_multi(x, tree_index$obs, user_order)

  return(res)
}


#' Transforms a fixest_multi object into a list
#'
#' Extracts the results from a `fixest_multi` object and place them into a list.
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
  nm = names(x)
  attributes(x) = NULL
  names(x) = nm
  x
}

#' Extracts the coefficients of fixest_multi objects
#'
#' Utility to extract the coefficients of multiple estimations and rearrange them into a matrix.
#'
#' @inheritParams etable
#' @inheritParams coef.fixest
#'
#' @param object A `fixest_multi` object. Obtained from a multiple estimation.
#' @param long Logical, default is `FALSE`. Whether the results should be displayed 
#' in a long format.
#' @param na.rm Logical, default is `TRUE`. Only applies when `long = TRUE`: whether to remove 
#' the coefficients with `NA` values.
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
#' # collin + long + na.rm
#' base$x1_bis = base$x1 # => collinear
#' est = feols(y ~ x1_bis + csw0(x1, x2, x3), base, split = ~species)
#'
#' # does not display x1 since it is always collinear
#' coef(est)
#' # now it does
#' coef(est, collin = TRUE)
#'
#' # long
#' coef(est, long = TRUE)
#'
#' # long but balanced (with NAs then)
#' coef(est, long = TRUE, na.rm = FALSE)
#'
#'
coef.fixest_multi = function(object, keep, drop, order, collin = FALSE,
                             long = FALSE, na.rm = TRUE, ...){
  # row: model
  # col: coefficient

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(collin, "logical scalar")

  res_list = list()
  for(i in seq_along(object)){
    res_list[[i]] = coef(object[[i]], collin = collin)
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

  # model information
  mod = models(object)

  if(long){
    res_long = c(t(res), recursive = TRUE)
    tmp = data.frame(coefficient =  res_long)
    mod_long = rep_df(mod, each = ncol(res))
    res = cbind(mod_long, coefficient = rep(all_names, NROW(res)),
          estimate = res_long)

    if(na.rm && anyNA(res$estimate)){
      res = res[!is.na(res$estimate), , drop = FALSE]
    }

  } else {
    res = cbind(mod, res)
  }

  res
}

#' @rdname coef.fixest_multi
coefficients.fixest_multi <- coef.fixest_multi

#' Extracts the coefficients tables from `fixest_multi` estimations
#'
#' Series of methods to extract the coefficients table or its sub-components from a 
#' `fixest_multi` objects (i.e. the outcome of multiple estimations).
#'
#' @inheritParams etable
#'
#' @param object A `fixest_multi` object, coming from a `fixest` multiple estimation.
#' @param wide A logical scalar, default is `FALSE`. If `TRUE`, then a list is returned: 
#' the elements of the list are coef/se/tstat/pvalue. Each element of the list is a wide 
#' table with a column per coefficient.
#' @param long Logical scalar, default is `FALSE`. If `TRUE`, then all the information 
#' is stacked, with two columns containing the information: `"param"` and `"value"`. 
#' The column `param` contains the values `coef`/`se`/`tstat`/`pvalue`.
#' @param ... Other arguments to be passed to [`summary.fixest`].
#'
#' @return
#' It returns a `data.frame` containing the coefficients tables (or just the se/pvalue/tstat) 
#' along with the information on which model was estimated.
#'
#' If `wide = TRUE`, then a list is returned. The elements of the list are 
#' coef/se/tstat/pvalue. Each element of the list is a wide table with a column per coefficient.
#'
#' If `long = TRUE`, then all the information is stacked. This removes the 4 columns 
#' containing the coefficient estimates to the p-values, and replace them with two 
#' new columns: `"param"` and `"value"`. The column `param` contains the 
#' values `coef`/`se`/`tstat`/`pvalue`, and the column `values` the 
#' associated numerical information.
#'
#' @examples
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' est_multi = feols(y ~ csw(x.[,1:3]), base, split = ~species)
#'
#' # we get all the coefficient tables at once
#' coeftable(est_multi)
#'
#' # Now just the standard-errors
#' se(est_multi)
#'
#' # wide = TRUE => leads toa  list of wide tables
#' coeftable(est_multi, wide = TRUE)
#'
#' # long = TRUE, all the information is stacked
#' coeftable(est_multi, long = TRUE)
#'
#'
#'
coeftable.fixest_multi = function(object, vcov = NULL, keep = NULL, drop = NULL, 
                                  order = NULL, long = FALSE, wide = FALSE, ...){

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(wide, "logical scalar | charin(se, pvalue, tstat)")
  check_arg(long, "logical scalar")

  if(long && !isFALSE(wide)){
    stop("You cannot have 'wide = TRUE' with 'long = TRUE', please choose.")
  }

  mod = models(object)

  res_list = list()
  for(i in seq_along(object)){
    ct = coeftable(object[[i]], vcov = vcov, keep = keep, drop = drop, order = order, ...)
    if(is.null(ct)){
      next
    }

    ct = cbind(coefficient = row.names(ct), as.data.frame(ct))

    mod_current = rep_df(mod[i, , drop = FALSE], each = nrow(ct))
    res_list[[length(res_list) + 1]] = cbind(mod_current, ct)
  }

  if(length(res_list) == 0){
    stop("Not any variable was selected: revise you 'keep'/'drop' arguments?")
  }

  res = do.call(rbind, res_list)
  row.names(res) = NULL

  if(!isFALSE(wide)){
    # we return a list of wide tables
    res_list = list()

    roots = c("coef", "se", "tstat", "pvalue")
    if(isTRUE(wide)) wide = roots

    i_coef = which(names(res) == "coefficient")

    for(i_select in which(roots %in% wide)){

      all_coef = unique(res$coefficient)
      all_id = unique(res$id)

      key = paste0(res$id, "; ", res$coefficient)
      value = res[, i_coef + i_select]
      names(value) = key

      df_xpd = expand.grid(coef = all_coef, id = all_id)
      new_key = paste0(df_xpd$id, "; ", df_xpd$coef)


      res_wide = matrix(value[new_key], nrow = length(all_id), ncol = length(all_coef),
                dimnames = list(NULL, all_coef), byrow = TRUE)

      item = cbind(mod[all_id, , drop = FALSE], as.data.frame(res_wide))
      if(length(wide) == 1){
        # internal call
        return(item)
      }

      res_list[[roots[i_select]]] = item
    }

    res = res_list
  }

  if(long){
    i_coef = which(names(res) == "coefficient")
    values = res[, i_coef + 1:4]

    values_all = c(t(values), recursive = TRUE)

    params = data.frame(param = rep(c("coef", "se", "tstat", "pvalue"), nrow(res)))

    info = rep_df(res[, 1:i_coef, drop = FALSE], each = 4)

    res = cbind(info, params, value = values_all)

  }


  res
}


#' @describeIn coeftable.fixest_multi Extracts the standard-errors from `fixest_multi` estimations
se.fixest_multi = function(object, vcov = NULL, keep = NULL, drop = NULL, 
                           order = NULL, long = FALSE, ...){
  # Default is wide format => same as with coef

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(long, "logical scalar")

  mc = match.call()
  if("wide" %in% names(mc)){
    stop("The argument 'wide' is not a valid argument.")
  }

  wide = if(long) FALSE else "se"

  res = coeftable(object, vcov = vcov, keep = keep, drop = drop, 
                  order = order, wide = wide, ...)

  if(long){
    i_coef = which(names(res) == "coefficient")
    res = res[, c(1:i_coef, i_coef + 2)]
  }

  res
}

#' @describeIn coeftable.fixest_multi Extracts the t-stats from `fixest_multi` estimations
tstat.fixest_multi = function(object, vcov = NULL, keep = NULL, drop = NULL, 
                              order = NULL, long = FALSE, ...){
  # Default is wide format => same as with coef

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(long, "logical scalar")

  mc = match.call()
  if("wide" %in% names(mc)){
    stop("The argument 'wide' is not a valid argument.")
  }

  wide = if(long) FALSE else "tstat"

  res = coeftable(object, vcov = vcov, keep = keep, drop = drop, 
                  order = order, wide = wide, ...)

  if(long){
    i_coef = which(names(res) == "coefficient")
    res = res[, c(1:i_coef, i_coef + 3)]
  }

  res
}

#' @describeIn coeftable.fixest_multi Extracts the p-values from `fixest_multi` estimations
pvalue.fixest_multi = function(object, vcov = NULL, keep = NULL, drop = NULL, 
                               order = NULL, long = FALSE, ...){
  # Default is wide format => same as with coef

  check_arg(keep, drop, order, "NULL character vector no na")
  check_arg(long, "logical scalar")

  mc = match.call()
  if("wide" %in% names(mc)){
    stop("The argument 'wide' is not a valid argument.")
  }

  wide = if(long) FALSE else "pvalue"

  res = coeftable(object, vcov = vcov, keep = keep, drop = drop, 
                  order = order, wide = wide, ...)

  if(long){
    i_coef = which(names(res) == "coefficient")
    res = res[, c(1:i_coef, i_coef + 4)]
  }

  res
}


#' Extracts the residuals from a `fixest_multi` object
#'
#' Utility to extract the residuals from multiple `fixest` estimations. If possible, 
#' all the residuals are coerced into a matrix.
#'
#' @inheritParams resid.fixest
#'
#' @param object A `fixes_multi` object.
#' @param na.rm Logical, default is `FALSE`. Should the NAs be kept? If `TRUE`, they are removed.
#'
#' @return
#' If all the models return residuals of the same length, a matrix is returned. Otherwise, 
#' a `list` is returned.
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
resid.fixest_multi = function(object, type = c("response", "deviance", "pearson", "working"), 
                              na.rm = FALSE, ...){

  # Je fais un prototype pour le moment, je l'ameliorerai apres (07-04-2021)

  check_set_arg(type, "match")
  check_set_arg(na.rm, "logical scalar")

  res_list = list()
  for(i in seq_along(object)){
    res_list[[i]] = resid(object[[i]], type = type, na.rm = na.rm)
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
residuals.fixest_multi = resid.fixest_multi




#' Confidence intervals for `fixest_multi` objects
#'
#' Computes the confidence intervals of parameter estimates for `fixest`'s multiple 
#' estimation objects (aka `fixest_multi`).
#'
#' @inheritParams confint.fixest
#'
#' @param object A `fixest_multi` object obtained from a multiple estimation in `fixest`.
#'
#' @return
#' It returns a data frame whose first columns indicate which model has been estimated. 
#' The last three columns indicate the coefficient name, and the lower and upper 
#' confidence intervals.
#'
#' @examples
#'
#' base = setNames(iris, c("y", "x1", "x2", "x3", "species"))
#' est = feols(y ~ csw(x.[,1:3]) | sw0(species), base, vcov = "iid")
#'
#' confint(est)
#'
#' # focusing only on the coefficient 'x3'
#' confint(est, "x3")
#'
#' # the 'id' provides the index of the estimation
#' est[c(3, 6)]
#'
confint.fixest_multi = function(object, parm, level = 0.95, vcov = NULL, se = NULL,
                                cluster = NULL, ssc = NULL, ...){

  n = length(object)
  confint_all = vector("list", n)
  for(i in 1:n){
    confint_all[[i]] = confint(object[[i]], parm = parm, level = level,
                   vcov = vcov, se = se, cluster = cluster,
                   ssc = ssc, coef.col = TRUE, internal = TRUE, ...)
  }

  n_all = sapply(confint_all, NROW)

  if(max(n_all) == 0){
    stop("No coefficient could be selected. Revise the argument `parm`?")
  }

  mod = models(object)

  # Formatting
  mod_all = rep_df(mod, times = n_all)
  res = do.call(base::rbind, confint_all)
  res = cbind(mod_all, res)

  attr(res, "type") = attr(confint_all[n_all > 0][[1]], "type")

  res
}











