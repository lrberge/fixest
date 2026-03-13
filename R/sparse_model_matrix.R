# Contributed by @kylebutts; please tag him in github issues

#' Design matrix of a `fixest` object returned in sparse format
#'
#' This function creates the left-hand-side or the right-hand-side(s) of a [`femlm`], [`feols`] or [`feglm`] estimation.
#'
#' @inheritParams nobs.fixest
#' @inheritParams model.matrix.fixest
#'
#' @param data If missing (default) then the original data is obtained by evaluating the `call`.
#'  Otherwise, it should be a `data.frame`.
#' @param type Character vector or one sided formula, default is "rhs". 
#' Contains the type of matrix/data.frame to be returned. Possible values are: 
#' "lhs", "rhs", "fixef", "iv.rhs1" (1st stage RHS), "iv.rhs2" (2nd stage RHS), 
#' "iv.endo" (endogenous vars.), "iv.exo" (exogenous vars), "iv.inst" (instruments).
#' @param na.rm Default is `FALSE`. Should observations with NAs be removed from the matrix?
#' @param collin.rm Logical scalar. Whether to remove variables that were 
#' found to be collinear during the estimation. Beware: it does not perform a 
#' collinearity check and bases on the `coef(object)`. 
#' Default is TRUE if object is a `fixest` object, or `FALSE` if object is a formula.
#' @param combine Logical scalar, default is `TRUE`. Whether to combine each 
#' resulting sparse matrix.
#' @param ... Not currently used.
#'
#' @return
#' It returns either a single sparse matrix a list of matrices, 
#' depending whether `combine` is `TRUE` or `FALSE`. 
#' The sparse matrix is of class `dgCMatrix` from the `Matrix` package.
#'
#' @seealso
#' See also the main estimation functions [`femlm`], [`feols`] or [`feglm`]. [`formula.fixest`], [`update.fixest`], [`summary.fixest`], [`vcov.fixest`].
#'
#'
#' @author
#' Laurent Berge, Kyle Butts
#'
#' @examples
#'
#' est = feols(wt ~ i(vs) + hp | cyl, mtcars)
#' sparse_model_matrix(est)
#' sparse_model_matrix(wt ~ i(vs) + hp | cyl, mtcars)
#'
sparse_model_matrix = function(object, data, type = "rhs", sample = "estimation", 
                               na.rm = FALSE, collin.rm = NULL, 
                               combine = TRUE, ...) {
  # We evaluate the formula with the past call
  # type: lhs, rhs, fixef, iv.endo, iv.inst, iv.rhs1, iv.rhs2
  # if fixef => return a DF

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The `sparse_model_matrix` method requires the `Matrix` package.")
  }

  # Checking the arguments
  if (is_user_level_call()) {
    validate_dots(suggest_args = c("data", "type"))
  }
  
  check_set_arg(sample, "match(estimation, original)")

  # We allow type to be used in the location of data if data is missing
  if (!missing(data) && missing(type)) {
    sc = sys.call()
    if (!"data" %in% names(sc)) {
      if (!is.null(data) && (is.character(data) || "formula" %in% class(data))) {
        # data is in fact the type
        type = data
        data = NULL
      }
    }
  }

  type = check_set_types(type, c("lhs", "rhs", "fixef", "iv.endo", "iv.inst", 
                                 "iv.exo", "iv.rhs1", "iv.rhs2"))

  if (isTRUE(object$is_fit)) {
    stop("sparse_model_matrix method not available for fixest estimations obtained from fit methods.")
  }

  if (any(grepl("^iv", type)) && !isTRUE(object$iv)) {
    stop("The type", enumerate_items(grep("^iv", type, value = TRUE), "s.is"), " only valid for IV estimations.")
  }

  check_arg(subset, "logical scalar | character vector no na")

  if (missing(collin.rm)) {
    collin.rm = if (inherits(object, "formula")) FALSE else TRUE
  } else {
    check_arg_plus(collin.rm, "logical scalar")
  }

  # Evaluation with the data
  is_original_data = FALSE
  if (missnull(data)) {
    is_original_data = TRUE
    data = fetch_data(object, "To apply 'sparse_model_matrix', ")
  }

  # control of the data
  if (is.matrix(data)) {
    if (is.null(colnames(data))) {
      stop("If argument 'data' is to be a matrix, its columns must be named.")
    }
    data = as.data.frame(data)
  }

  if (!"data.frame" %in% class(data)) {
    stop("The argument 'data' must be a data.frame or a matrix.")
  }


  # Allows passage of formula to sparse_model_matrix. A bit inefficient, but it works.
  isFormula = FALSE
  if (inherits(object, "formula")) {
    fml_parts = get_fml_parts(object, parts_as_formula = TRUE)
    if (isTRUE(collin.rm)) {
      message("Formula passed to sparse_model_matrix with collin.rm == TRUE. Estimating feols with formula.")
      object = feols(object, data = data)
    } else if (!is.null(fml_parts$inst_fml)) {
      message("Formula passed to sparse_model_matrix with an iv formula. Estimating feols with formula.")
      object = feols(object, data = data)
    } else {
      isFormula = TRUE
    }
  } 

  # na.rm = FALSE doesn't work with type = "fixef" (which FE col gets NA?)
  if (("fixef" %in% type & !na.rm)) {
    # na.rm = TRUE
    message("na.rm = FALSE doesn't work with type = 'fixef'. It has been set to TRUE.")
  }

  # Panel setup
  if (!isFormula) {
    panel__meta__info = set_panel_meta_info(object, data)
  }

  # The formulas
  if (isFormula) {
    fml_parts = get_fml_parts(object, parts_as_formula = TRUE)
  } else {
    fml_parts = get_fml_parts(formula(object), parts_as_formula = TRUE)
  }  
  fml_0 = attr(stats::terms(fml_parts$exo_fml), "intercept") == 0
  # we kick out the intercept if there is presence of fixed-effects
  fake_intercept = !is.null(fml_parts$fe_fml) | fml_0
  res = list()

  if ("lhs" %in% type) {
    lhs = eval(fml_parts$y_fml[[2]], data) 
    lhs = Matrix::Matrix(lhs, sparse = TRUE, ncol = 1)
    colnames(lhs) = deparse_long(fml_parts$y_fml[[2]])

    res[["lhs"]] = lhs
  }

  if ("rhs" %in% type && !isTRUE(object$onlyFixef)) {
    fml = .xpd(~ ..endo + ..rhs, ..rhs = fml_parts$exo_fml, ..endo = fml_parts$endo_fml)
    vars = attr(stats::terms(fml), "term.labels")
    linear.mat = vars_to_sparse_mat(vars = vars, data = data, object = object, 
                                    collin.rm = collin.rm, type = "rhs", 
                                    add_intercept = !fake_intercept)

    res[["rhs"]] = linear.mat
  }

  if ("fixef" %in% type) {

    if (is.null(fml_parts$fe_fml)) {
      mat_FE = NULL
    } else { 
      fixef_fml = .xpd(~ ..fe, ..fe = fml_parts$fe_fml)
      fixef_terms_full = fixef_terms(fixef_fml)
      fe_vars = fixef_terms_full$fe_vars
      slope_var_list = fixef_terms_full$slope_vars_list

      fixef_df = prepare_df(fe_vars, data, fixef.keep_names = TRUE)

      # Check for slope vars
      if (any(fixef_terms_full$slope_flag > 0)) {
        slope_df = prepare_df(unlist(slope_var_list), data)
      }

      cols_lengths = c(0)
      total_cols = 0
      n_FE = 0
      nrows = nrow(data)
      for (i in seq_along(fe_vars)) {
        fe_var = fe_vars[i]
        slope_vars = slope_var_list[[i]]
        n_slope_vars = if (is.null(slope_vars)) 0 else length(slope_vars)

        fe = fixef_df[[fe_var]]
        unique_fe = unique(fe)
        n_cols = length(unique_fe[!is.na(unique_fe)])

        total_cols = total_cols + n_cols * (1 + n_slope_vars)
        cols_lengths = c(cols_lengths, rep(n_cols, 1 + n_slope_vars))
        n_FE = n_FE + 1 + n_slope_vars
      }
      running_cols = cumsum(cols_lengths)

      id_all = names_all = val_all = vector("list", n_FE)
      rowid = 1:nrows
      j = 1
      for (i in seq_along(fe_vars)) {

        fe_var = fe_vars[i]
        xi = fixef_df[[fe_var]]

        keep = which(!is.na(xi))
        if (length(keep) == 0) stop("All values of the fixed-effect variable '", fe_var, "' are NA.")

        xi = xi[keep]
        index_info = to_index_internal(xi, add_items = TRUE)

        col_id = index_info$x
        col_levels = as.character(index_info$items)

        slope_vars = slope_var_list[[i]]
        n_slope_vars = if (is.null(slope_vars)) 0 else length(slope_vars)

        # Be careful with NAs
        # First fixed-effect by itself
        val_all[[j]] = c(rep(1, length(col_id)), rep(NA, nrows - length(keep)))
        
        id_all[[j]] = cbind(
          c(
            rowid[keep], 
            rowid[-keep]
          ), 
          c( 
            running_cols[j] + col_id,  
            rep(running_cols[j] + 1, nrows - length(keep))
          )
        )
        
        names_all[[j]] = paste0(fe_var, "::", col_levels)
        j = j + 1

        for (k in seq_along(slope_vars)) {
          slope_var = slope_vars[k]
          slope = slope_df[[slope_var]]

          val_all[[j]] = c(as.numeric(slope[keep]), rep(NA, nrows - length(keep)))
          id_all[[j]] = cbind(
            c(
              rowid[keep],
              rowid[-keep]
            ), 
            c(
              running_cols[j] + col_id, 
              rep(running_cols[j] + 1, nrows - length(keep))
            )
          )
          
          names_all[[j]] = paste0(fe_var, "[[", slope_var, "]]", "::", col_levels)
          j = j + 1
        }
      }

      id_mat = do.call(rbind, id_all)
      val_vec = unlist(val_all)
      names_vec = unlist(names_all)

      mat_FE = Matrix::sparseMatrix(
        i = id_mat[, 1],
        j = id_mat[, 2],
        x = val_vec,
        dimnames = list(NULL, names_vec)
      )


      # Keep non-zero FEs
      if (collin.rm == TRUE) {
        fixefs = fixef(object, sorted = TRUE)

        select =	lapply(
          names(fixefs), 
          function(var) {
            names = names(fixefs[[var]])
            names = names[fixefs[[var]] != 0]

            paste0(var, "::", names)
          }
        )
        select = unlist(select)

        # When is_original_data isn't used, some FEs may not be in the new dataset, add them in
        missing_cols = setdiff(select, colnames(mat_FE))
        mat_FE = cbind(
          mat_FE, 
          Matrix::Matrix(0, ncol = length(missing_cols), nrow = nrow(mat_FE), 
                         sparse = TRUE, dimnames = list(NULL, missing_cols))
        )

        # Subset fixef and order columns to match fixef(object)
        idx = which(colnames(mat_FE) %in% select)
        mat_FE = mat_FE[, idx, drop = FALSE]

        # Reorder columns to match fixef(object)
        reorder = match(unlist(select), colnames(mat_FE))
        mat_FE = mat_FE[, reorder[!is.na(reorder)], drop = FALSE]
      }

    }

    res[["fixef"]] = mat_FE
  }

  #
  # IV
  #
  if ("iv.endo" %in% type) {
    fml = .xpd(~ ..endo, ..endo = fml_parts$endo_fml)
    vars = attr(stats::terms(fml), "term.labels")
    endo.mat = vars_to_sparse_mat(vars = vars, data = data, object, 
                                  collin.rm = collin.rm, type = "iv.endo")

    res[["iv.endo"]] = endo.mat
  }

  if ("iv.inst" %in% type) {
    fml = .xpd(~ ..inst, ..inst = fml_parts$inst_fml)
    vars = attr(stats::terms(fml), "term.labels")
    inst.mat = vars_to_sparse_mat(vars = vars, data = data, object, 
                                  collin.rm = collin.rm, type = "iv.inst")

    res[["iv.inst"]] = inst.mat
  }

  if ("iv.exo" %in% type) {
    fml = .xpd(~ ..exo, ..exo = fml_parts$exo_fml)
    vars = attr(stats::terms(fml), "term.labels")
    exo.mat = vars_to_sparse_mat(vars = vars, data = data, object = object, 
                                 collin.rm = collin.rm, type = "iv.exo", 
                                 add_intercept = !fake_intercept)

    res[["iv.exo"]] = exo.mat
  }

  if ("iv.rhs1" %in% type) {
    fml = .xpd(~ ..inst + ..rhs, ..rhs = fml_parts$exo_fml, ..inst = fml_parts$inst_fml)
    vars = attr(stats::terms(fml), "term.labels")
    iv_rhs1 = vars_to_sparse_mat(vars = vars, data = data, object = object, 
                                    collin.rm = collin.rm, type = "iv.rhs1", 
                                    add_intercept = !fake_intercept)

    res[["iv.rhs1"]] = iv_rhs1
  }

  if ("iv.rhs2" %in% type) {
    # Second stage
    if (!isFormula && !object$iv_stage == 2) {
      stop("In model.matrix, the type 'iv.rhs2' is only valid for second stage IV models. This estimation is the first stage.")
    }

    # I) we get the fitted values
    stage_1 = object$iv_first_stage

    fit_vars = c()
    for (i in seq_along(stage_1)) {
      fit_vars[i] = v = paste0("fit_", names(stage_1)[i])
      data[[v]] = predict(stage_1[[i]], newdata = data)
    }

    # II) we create the variables
    fml = object$fml
    fml = .xpd(~ ..fit + ..rhs, ..fit = fit_vars, ..rhs = fml_parts$exo_fml)
    vars = attr(stats::terms(fml), "term.labels")
    iv_rhs2 = vars_to_sparse_mat(vars = vars, data = data, object = object, 
                                 collin.rm = collin.rm, type = "iv.rhs2", 
                                 add_intercept = !fake_intercept)

    res[["iv.rhs2"]] = iv_rhs2
  }
  
  if(is_original_data){
    if(sample == "estimation"){
      
      if(length(object$obs_selection) > 0){
        for(var in names(res)){
          mat = res[[var]]
          for(obs in object$obs_selection){
            mat = mat[obs, , drop = FALSE]
          }
          res[[var]] = mat
        }
      }
      
      na.rm = FALSE
    }
  }

  if (na.rm) {
    na_rows = lapply(res, function(mat) {
      # Get rows with NA 
      1 + mat@i[is.na(mat@x)]
    })

    na_rows = unlist(na_rows, use.names = FALSE)
    na_rows = unique(na_rows)

    if (length(na_rows) > 0) {
      res = lapply(res, function(mat) { 
        mat[-na_rows, , drop = FALSE]
      })
    }
  } 

  # Formatting res
  if (length(res) == 0) {
    return(NULL)
    
  } else if (length(type) > 1) {
    if (combine) {
      res = res[type]
      res = do.call(cbind, unname(res))
    }
  } else {
    res = res[[1]]
  }

  res
}



# Internal: modifies the calls so that each variable/interaction is evaluated with mult_sparse
mult_wrap = function(x) {
  # x: character string of a variable to be evaluated
  # ex: "x1" => mult_sparse(x1)
  #     "x1:factor(x2):x3" => mult_sparse(x3, factor(x2), x1)
  #
  # We also add the argument sparse to i()
  #     "x1:i(species, TRUE)" => mult_sparse(x1, i(species, TRUE, sparse = TRUE))

  x_call = str2lang(x)

  res = (~ mult_sparse())[[2]]

  if (length(x_call) == 1 || x_call[[1]] != ":") {
    res[[2]] = x_call

  } else {
    res[[2]] = x_call[[3]]
    tmp = x_call[[2]]

    while (length(tmp) == 3 && tmp[[1]] == ":") {
      res[[length(res) + 1]] = tmp[[3]]
      tmp = tmp[[2]]
    }

    res[[length(res) + 1]] = tmp
  }

  # We also add sparse to i() if found
  for (i in 2:length(res)) {
    ri = res[[i]]
    if (length(ri) > 1 && ri[[1]] == "i") {
      ri[["sparse"]] = TRUE
      res[[i]] = ri
    }
  }

  if (length(res) > 2) {
    # we restore the original order
    res[-1] = rev(res[-1])
  }

  return(res)
}

# Internal function to evaluate the variables (and interactions) in a sparse way
mult_sparse = function(...) {
  # Only sparsifies factor variables
  # Takes care of interactions

  dots = list(...)
  n = length(dots)
  mc = match.call()
  var_call = sapply(seq_len(n), function(i) deparse_long(mc[[i+1]]))

  num_var = NULL
  factor_list = list()
  factor_idx = c()
  i_idx = c()
  info_i = NULL
  is_i = is_factor = FALSE
  # You can't have interactions between i and factors, it's either

  for (i in 1:n) {
    xi = dots[[i]]
    vari = mc[[i+1]]
    if (is.numeric(xi)) {
      # We stack the product
      num_var = if (is.null(num_var)) xi else xi * num_var
      
    } else if (inherits(xi, "i_sparse")) {
      is_i = TRUE
      info_i = xi
      i_idx = c(i_idx, i)
      
    } else {
      is_factor = TRUE
      factor_list[[length(factor_list) + 1]] = xi
      factor_idx = c(factor_idx, i)
    }
  }

  if (is_factor && is_i) {
    stop("Unfortunately, can not interact factor with `i()` in sparse_model_matrix")
  }

  # numeric
  if (!is_i && !is_factor) {
    return(num_var)
  }
  # factor()
  if (is_factor) {
    factor_list$add_items = TRUE
    factor_list$items.list = TRUE
    factor_list$multi.join = " ; "
    fact_as_int = do.call(to_integer, factor_list)

    values = if (is.null(num_var)) rep(1, length(fact_as_int$x)) else num_var
    rowid = seq_along(values)

    # messy, but need this to get things like `factor(am)1:hp:factor(cyl)6`
    if(is.character(fact_as_int$items)){
      # multiple elements are combined with " ; " ...
      items_mat = do.call(rbind, strsplit(fact_as_int$items, " ; "))
    } else {
      # ... single elements are not
      items_mat = cbind(as.character(fact_as_int$items))
    }
    
    # intersplice var_call with items 
    col_names = sapply(seq_len(nrow(items_mat)), function(i) {
      pasteable = var_call 
      pasteable[factor_idx] = paste0(var_call[factor_idx], items_mat[i, ])
      paste(pasteable, collapse = ":")
    }, USE.NAMES = FALSE)

    res = list(rowid = rowid, colid = fact_as_int$x, values = values,
               col_names = col_names, n_cols = length(fact_as_int$items))
  # i()
  } else {

    values = info_i$values
    if (!is.null(num_var)) {
      num_var = num_var[info_i$rowid]
      values = values * num_var
      col_names = sapply(info_i$col_names, function(item) {
        pasteable = var_call
        pasteable[i_idx] = item
        paste(pasteable, collapse = ":")
      }, USE.NAMES = FALSE)
    } else {
      col_names = info_i$col_names
    }

    res = list(rowid = info_i$rowid, colid = info_i$colid,
               values = values[info_i$rowid],
               col_names = col_names,
               n_cols = length(info_i$col_names))
  }

  class(res) = "sparse_var"

  res
}

# Takes a vector of strings denoting the variables (including terms like `poly()`, `i()`, `I()`, `lag()`, etc.) and returns a sparse matrix of the variables extracted from `data`.
vars_to_sparse_mat = function(vars, data, collin.rm = FALSE, object = NULL, 
                              type = NULL, add_intercept = FALSE) {
  
  if (collin.rm & is.null(object)) {
    stop("You need to provide the 'object' argument to use 'collin.rm = TRUE'.")
  }

  if (length(vars) == 0) {
    # Case only FEs
    mat = NULL
  } else {

    # Since we don't want to evaluate the factors,
    # the code is a bit intricate because we have to catch them before
    # any interaction takes place
    #
    # that's why I wrap interactions in a function (mult_sparse())
    #

    # Below, we evaluate all the variables in a "sparse" way

    vars_calls = lapply(vars, mult_wrap)

    n = length(vars)
    variables_list = vector("list", n)
    for (i in 1:n) {
      variables_list[[i]] = eval(vars_calls[[i]], data)
    }

    # To create the sparse matrix, we need the column indexes
    total_cols = 0
    running_cols = c(0)
    for (i in 1:n) {
      xi = variables_list[[i]]
      if (inherits(xi, "sparse_var")) {
        total_cols = total_cols + xi$n_cols
      } else {
        total_cols = total_cols + NCOL(xi)
      }
      running_cols[i + 1] = total_cols
    }

    # We just create a sparse matrix and fill it

    # 1) creating the indexes + names

    # NOTA: I use lists to avoid creating copies
    rowid = 1:nrow(data)
    id_all = values_all = names_all = vector("list", n)
    for (i in 1:n) {
      xi = variables_list[[i]]

      if (inherits(xi, "sparse_var")) {

        id_all[[i]] = cbind(xi$rowid, running_cols[i] + xi$colid)
        values_all[[i]] = xi$values
        names_all[[i]] = xi$col_names

      } else if (NCOL(xi) == 1) {

        id_all[[i]] = cbind(rowid, running_cols[i] + 1)
        values_all[[i]] = xi
        names_all[[i]] = vars[[i]]

      } else {

        colid = rep(1:NCOL(xi), each = nrow(data))
        id_all[[i]] = cbind(rep(rowid, NCOL(xi)), running_cols[i] + colid)
        values_all[[i]] = as.vector(xi)
        if (!is.null(colnames(xi))) {
          names_all[[i]] = paste0(vars[[i]], colnames(xi))
        } else {
          names_all[[i]] = paste0(vars[[i]], 1:NCOL(xi))
        }

      }
    }

    id_mat = do.call(rbind, id_all)
    values_vec = unlist(values_all)
    names_vec = unlist(names_all)

    # 2) filling the matrix: one shot, no copies

    mat = Matrix::Matrix(0, nrow(data), total_cols, dimnames = list(NULL, names_vec))
    mat[id_mat] = values_vec

  }

  if (collin.rm) {

    if (isTRUE(object$iv)) {
      fml_iv = object$fml_all$iv
      endo = fml_iv[[2]]        
      # Trick to get the rhs variables as a character vector
      endo = .xpd(~ ..endo, ..endo = endo)
      endo = attr(stats::terms(endo), "term.labels")

      exo = lapply(object$iv_first_stage, function(x) names(stats::coef(x)))
      exo = unique(unlist(exo, use.names = FALSE))
 
      inst = setdiff(exo, names(object$coefficients))
    }
    
    # Custom subsetting for collin.rm depending on `type`
    if (is.null(type)) { 
      coef_names = names(object$coefficients)
      add_intercept = "(Intercept)" %in% coef_names
    } else if (type == "rhs") {
      if (isTRUE(object$iv)) {
        coef_names = c(endo, names(object$coefficients))
        coef_names = coef_names[!grepl("^fit_", coef_names)]
        add_intercept = "(Intercept)" %in% coef_names
      } else {
        coef_names = names(object$coefficients)
        add_intercept = "(Intercept)" %in% coef_names
      }
    } else if (type == "iv.exo") { 
      coef_names = exo
      add_intercept = "(Intercept)" %in% coef_names
    } else if (type == "iv.inst") {
      coef_names = inst
      add_intercept = FALSE
    } else if (type == "iv.endo") {
      coef_names = endo
      add_intercept = FALSE
    } else if (type == "iv.rhs1") {
      coef_names = exo
      add_intercept = "(Intercept)" %in% coef_names
    } else if (type == "iv.rhs2") {
      coef_names = names(object$coefficients)
      add_intercept = "(Intercept)" %in% coef_names
    }
    keep = colnames(mat) %in% coef_names
    if (length(keep) == 0) {
      mat = NULL
    } else {
      mat =  mat[, keep, drop = FALSE]
    }
  }

  # add intercept if needed and not already present
  # when collin.rm = TRUE, we change add_intercept depending on if it is needed
  if (add_intercept && !("(Intercept)" %in% colnames(mat))) {
    mat = cbind(1, mat)
    colnames(mat)[1] = "(Intercept)"
  }

  return(mat)
}

get_fml_parts <- function(formula, parts_as_formula = FALSE) {
  has_lhs <- !is_rhs_only(formula)
  fml_split_tilde <- sm_fml_breaker(formula, "~")

  res <- list(
    y_fml = NULL,
    exo_fml = NULL,
    fe_fml = NULL,
    endo_fml = NULL,
    inst_fml = NULL
  )

  # LHS
  if (has_lhs) {
    res$y_fml <- fml_split_tilde[[length(fml_split_tilde)]]
    # Drop y and treat everything as RHS only formula
    fml_split_tilde <- fml_split_tilde[-length(fml_split_tilde)]
  }

  # OLS
  if (length(fml_split_tilde) == 1) {
    fml_split <- sm_fml_breaker(fml_split_tilde[[1]], "|")
    if (length(fml_split) == 2) {
      res$exo_fml <- fml_split[[2]]
      res$fe_fml <- fml_split[[1]]
    } else {
      res$exo_fml <- fml_split[[1]]
    }
  }

  # IV
  if (length(fml_split_tilde) == 2) {
    res$inst_fml <- fml_split_tilde[[1]]

    # This works b/c we know there is an IV
    fml_W_T_split <- sm_fml_breaker(fml_split_tilde[[2]], "|")
    if (length(fml_W_T_split) == 3) {
      res$endo_fml <- fml_W_T_split[[1]]
      res$fe_fml <- fml_W_T_split[[2]]
      res$exo_fml <- fml_W_T_split[[3]]
    } else {
      res$endo_fml <- fml_W_T_split[[1]]
      res$exo_fml <- fml_W_T_split[[2]]
    }
  }

  if (parts_as_formula) { 
    res = lapply(res, function(rhs) {
      if (is.null(rhs)) return(NULL)
      fml <- ~.
      fml[[2]] <- rhs
      return(fml)
    })
  }
  return(res)
}

is_rhs_only <- function(fml) {
  # e.g. fml = ~ x
  if (length(fml) == 2) {
    return(TRUE)
  }
  # e.g. fml = ~ x | t ~ z
  if (length(fml[[2]]) == 2) {
    return(TRUE)
  }
  return(FALSE)
}

# fml_breaker but works with RHS only
sm_fml_breaker <- function(fml, op) {
  res <- list()
  k <- 1
  while (is_operator(fml, op) & length(fml) == 3) {
    res[[k]] <- fml[[3]]
    k <- k + 1
    fml <- fml[[2]]
  }

  if (length(fml) == 2) {
    res[[k]] <- fml[[2]]
  } else {
    res[[k]] <- fml
  }
  res
}
