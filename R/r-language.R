#------------------------------------------------------------------------------#
# Author: Laurent R. Bergé
# Created: 2025-07-24
# ~: functions related to R-langage manipulation
#------------------------------------------------------------------------------#


is_formula = function(x){
  is.call(x) && length(x[[1]]) == 1 && x[[1]] == "~"
}


get_all_vars = function(expr, sort = FALSE, ignore_formula = FALSE, 
                        interpol = FALSE, i.prefix = FALSE){
  
  if(length(expr) == 1){
    if(is.name(expr)){
      return(as.character(expr))
      
    } else if(interpol && is.character(expr)){
      vars = stringmagic::get_interpolated_vars(expr)
      return(vars)
    }
    
    return(NULL)
  }
  
  res = character(0)
  i_start = if(is.expression(expr)) 1 else 2
  
  pos_prefix = 0
  if(i.prefix && i_start == 2 && length(expr) >= 3 && is_operator(expr, "i")){
    pos_prefix = 3
    if(!is.null(names(expr)) && "var" %in% names(expr)){
      pos_prefix = which("var" %in% names(expr))
    }
  }
  
  for(i in i_start:length(expr)){
    element = expr[[i]]
    if(length(element) == 1){
      if(is.name(element)){
        
        varname = as.character(element)
        if(i.prefix && i == pos_prefix && substr(varname, 1, 2) == "i."){
          varname = substr(varname, 3, 500)
        }
        
        res = c(res, varname)
      }
    } else if(ignore_formula && length(element[[1]]) == 1 && element[[1]] == "~"){
      # we ignore
      
    } else {
      vars = get_all_vars(element, ignore_formula = ignore_formula, 
                          interpol = interpol, i.prefix = i.prefix)
      
      if(length(vars) > 0){
        res = c(res, vars)
      }
    }
    
  }
  
  if(sort && length(res) > 1){
    return(unique(sort(res)))
  }
  
  res
}

get_all_vars_from_formula = function(fml){
  
  left = fml[[2]]
  right = fml[[3]]
  
  if(is_formula(left)){
    all_vars = c(
      get_all_vars(left[[2]], ignore_formula = TRUE),
      get_all_vars(left[[3]], ignore_formula = TRUE, i.prefix = TRUE),
      get_all_vars(right, ignore_formula = TRUE, i.prefix = TRUE)
    )
  } else {
    all_vars = c(
      get_all_vars(left, ignore_formula = TRUE),
      get_all_vars(right, ignore_formula = TRUE, i.prefix = TRUE)
    )
  }
  
  unique(all_vars)
}


####
#### formula ####
####

is_operator = function(x, op){
  if(length(x) <= 1) FALSE else x[[1]] == op
}

fml_breaker = function(fml, op){
  res = list()
  k = 1
  while(is_operator(fml, op)){
    res[[k]] = fml[[3]]
    k = k + 1
    fml = fml[[2]]
  }
  res[[k]] = fml

  res
}

fml_maker = function(lhs, rhs){

  while(is_operator(lhs, "(")){
    lhs = lhs[[2]]
  }

  if(missing(rhs)){
    if(is_operator(lhs, "~")){
      return(lhs)
    }
    res = ~ .
    res[[2]] = lhs
  } else {

    while(is_operator(rhs, "(")){
      rhs = rhs[[2]]
    }

    res = . ~ .
    res[[2]] = lhs
    res[[3]] = rhs
  }

  res
}

fml_split_internal = function(fml, split.lhs = FALSE){

  fml_split_tilde = fml_breaker(fml, "~")
  k = length(fml_split_tilde)

  # NOTA in fml_breaker: to avoid copies, the order of elements returned is reversed

  # currently res is the LHS
  res = list(fml_split_tilde[[k]])

  if(k == 2){
    rhs = fml_breaker(fml_split_tilde[[1]], "|")
    l = length(rhs)

  } else if(k == 3){
    rhs  = fml_breaker(fml_split_tilde[[2]], "|")
    l = length(rhs)
    rhs_right = fml_breaker(fml_split_tilde[[1]], "|")

    if(length(rhs_right) > 1){
      stop_up("Problem in the formula: the formula in the RHS (expressing the IVs) cannot be multipart.")
    }

    # The rightmost element of the RHS is in position 1!!!!
    iv_fml = fml_maker(rhs[[1]], rhs_right[[1]])

    rhs[[1]] = iv_fml

  } else {
    # This is an error
    stop_up("Problem in the formula: you cannot have more than one RHS part containing a formula.")
  }

  if(!split.lhs){
    new_fml = fml_maker(res[[1]], rhs[[l]])

    res[[1]] = new_fml
    if(l == 1) return(res)
    res[2:l] = rhs[(l - 1):1]

  } else {
    res[1 + (1:l)] = rhs[l:1]
  }

  res
}


fml_split = function(fml, i, split.lhs = FALSE, text = FALSE, raw = FALSE){
  # I had to create that function to cope with the following formula:
  #
  #                       y ~ x1 | fe1 | u ~ z
  #
  # For Formula to work well, one would need to write insstead: y | x1 | fe1 | (u ~ z)
  # that set of parentheses are superfluous

  my_split = fml_split_internal(fml, split.lhs)

  if(raw){
    return(my_split)
  } else if(text){

    if(!missing(i)){
      return(deparse_long(my_split[[i]]))
    } else {
      return(sapply(my_split, deparse_long))
    }

  } else if(!missing(i)) {
    return(fml_maker(my_split[[i]]))

  } else {
    res = lapply(my_split, fml_maker)

    return(res)
  }

}


# fml_char = "x + y + u^factor(v1, v2) + x5"
fml_combine = function(fml_char, fixef.keep_names, vars = FALSE){
  # function that transforms "hat" interactions into a proper function call:
  # Origin^Destination^Product + Year becomes ~combine_clusters(Origin, Destination, Product) + Year
  
  if(is.null(fixef.keep_names)){
    fun2combine = "combine_fixef"
  } else if(isTRUE(fixef.keep_names)){
    fun2combine = "combine_fixef_keep_names"
  } else {
    fun2combine = "combine_fixef_drop_names"
  }

  # we need to change ^ into %^% otherwise terms sends error
  labels = attr(terms(.xpd(rhs = gsub("\\^(?=[^0-9])", "%^%", fml_char, perl = TRUE))), 
                "term.labels")

  # now we work this out
  for(i in seq_along(labels)){
    lab = labels[i]
    if(grepl("^", lab, fixed = TRUE)){
      lab_split = trimws(strsplit(lab, "%^%", fixed = TRUE)[[1]])
      if(grepl("(", lab, fixed = TRUE)){
        # we add some error control -- imperfect, but... it's enough
        lab_collapsed = gsub("\\([^\\)]+\\)", "", lab)
        if(length(lab_split) != length(strsplit(lab_collapsed, "%^%", fixed = TRUE)[[1]])){
          msg = "Wrong formatting of the fixed-effects interactions. The `^` operator should not be within parentheses."
          stop(msg)
        }
      }
      labels[i] = paste0(fun2combine, "(", paste0(lab_split, collapse = ", "), ")")
    }
  }

  if(vars){
    return(labels)
  }

  fml = .xpd(rhs = labels)

  fml
}

# x = c('combine_clusters(bin(fe1, "!bin::2"), fe2)', 'fe3')
rename_hat = function(x){

  qui = grepl("combine_clusters", x, fixed = TRUE)
  if(!any(qui)) return(x)

  for(i in which(qui)){
    xi_new = gsub("combine_clusters(_fast)?", "sw", x[i])
    sw_only = extract_fun(xi_new, "sw")
    sw_eval = eval(str2lang(sw_only$fun))
    new_var = paste0(sw_only$before, paste0(sw_eval, collapse = "^"), sw_only$after)
    x[i] = new_var
  }

  x
}


fml2varnames = function(fml, combine_fun = FALSE){
  # This function transforms a one sided formula into a
  # character vector for each variable

  # In theory, I could just use terms.formula to extract the variable.
  # But I can't!!!! Because of this damn ^ argument.
  # I need to apply a trick

  # combine_fun: whether to add a call to combine_clusters_fast

  # Only the ^ users "pay the price"

  if("^" %in% all.vars(fml, functions = TRUE)){
    # new algo using fml_combine, more robust
    fml_char = as.character(fml)[2]
    all_var_names = fml_combine(fml_char, TRUE, vars = TRUE)
    if(!combine_fun){
      all_var_names = rename_hat(all_var_names)
    }

  } else {
    t = terms(fml)
    all_var_names = attr(t, "term.labels")
    # => maybe I should be more cautious below???? Sometimes `:` really means stuff like 1:5
    all_var_names = gsub(":", "*", all_var_names) # for very special cases
  }


  all_var_names
}


is_fml_inside = function(fml){
  # we remove parentheses first

  while(is_operator(fml, "(")){
    fml = fml[[2]]
  }

  is_operator(fml, "~")
}


merge_fml = function(fml_linear, fml_fixef = NULL, fml_iv = NULL){

  is_fe = length(fml_fixef) > 0
  is_iv = length(fml_iv) > 0

  if(!is_fe && !is_iv){
    res = fml_linear
  } else {
    fml_all = deparse_long(fml_linear)

    if(is_fe){
      # we add parentheses if necessary
      if(is_operator(fml_fixef[[2]], "|")){
        fml_all[[2]] = paste0("(", as.character(fml_fixef)[2], ")")
      } else {
        fml_all[[2]] = as.character(fml_fixef)[2]
      }
    }

    if(is_iv) fml_all[[length(fml_all) + 1]] = deparse_long(fml_iv)

     res = as.formula(paste(fml_all, collapse = "|"), .GlobalEnv)
  }

  res
}

lag_expand = function(x, k = 1, fill = NA){
  mc = match.call()
  
  if(length(k) > 1){
    args = lapply(as.numeric(k), function(i) {mc$k = i ; names(mc) = NULL ; deparse_long(mc)})
    args$sep = "+"
    res = do.call(paste, args)
    res = str2lang(res)
    return(res)
  }
  
  names(mc) = NULL
  
  mc
}

protect_powers_expand_lags = function(expr){
  # x1^2 + l(x2, 1:2) => I(x1^2) + l(x2, 1) + l(x2, 2)
  
  if(length(expr) == 1){
    return(expr)
  }
  
  if(length(expr[[1]]) == 1 && as.character(expr[[1]]) %in% c("l", "f", "d")){
    # lags
    
    envlist = list(lag_expand)
    names(envlist) = as.character(expr[[1]])
    new_expr = eval(expr, envlist)
    attr(new_expr, "lag_expand") = TRUE
    
    return(new_expr)
    
  } else if(length(expr[[1]]) == 1 && expr[[1]] == "^"){
    # power
    
    if(is.numeric(expr[[3]])){
      new_expr = quote(I(x))
      new_expr[[2]] = expr
      return(new_expr)
    }
    
  } else {
    is_plus = is_operator(expr, "+")
    fix_plus = FALSE
    for(i in 2:length(expr)){
      new_element = protect_powers_expand_lags(expr[[i]])
      
      if(!identical(new_element, expr[[i]])){
        if(is_plus && i == 3 && isTRUE(attr(new_element, "lag_expand"))){
          fix_plus = TRUE
        }
      }
      
      expr[[i]] = new_element
    }
    
    if(fix_plus){
      expr = str2lang(paste0(deparse_long(expr[[2]]), "+", deparse_long(expr[[3]])))
    }
    
  }
  
  expr
}

lag_expand_comma = function(x, k = 1, fill = NA){
  # returns a string
  mc = match.call()
  
  if(length(k) > 1){
    args = lapply(as.numeric(k), function(i) {mc$k = i ; names(mc) = NULL ; deparse_long(mc)})
    args$sep = ","
    res = do.call(paste, args)
    return(res)
  }
  
  names(mc) = NULL
  
  mc
}

expand_lags_lhs = function(expr, is_root = TRUE){
  # l(y, 1:2) => c(l(y, 1), l(y, 2))
  
  if(length(expr) == 1){
    return(expr)
  }
  
  if(length(expr[[1]]) == 1 && as.character(expr[[1]]) %in% c("l", "f", "d")){
    # lags
    
    envlist = list(lag_expand_comma)
    names(envlist) = as.character(expr[[1]])
    new_expr = eval(expr, envlist)
    attr(new_expr, "lag_expand") = TRUE
    
    if(is_root && !is.call(new_expr)){
      new_expr = str2lang(paste0("c(", new_expr, ")"))
      return(new_expr)
    }
    
    return(new_expr)
    
  } else {
    is_binary = length(expr[[1]]) == 1 && as.character(expr[[1]]) %in% c("+", "-", "*", "/")
    
    if(is_binary){
      # we avoid errors even if what the user wants to do is nonsensical
      
      is_plus = expr[[1]] == "+"
      fix_plus = FALSE
      for(i in 2:length(expr)){
        new_element = protect_powers_expand_lags(expr[[i]])
        
        if(!identical(new_element, expr[[i]])){
          if(is_plus && isTRUE(attr(new_element, "lag_expand"))){
            fix_plus = TRUE
          }
        }
        
        expr[[i]] = new_element
      }
      
      if(fix_plus){
        expr = str2lang(paste0(deparse_long(expr[[2]]), "+", deparse_long(expr[[3]])))
      }
      
    } else {
      # in the left hand side the lags are extended with commas
      # c(x^2, l(y, 1:2)) => c(x^2, l(y, 1), l(y, 2))
      
      if(is_root && length(expr[[1]]) == 1 && as.character(expr[[1]]) %in% c("sw", "sw0", "csw", "csw0")){
        expr[[1]] = as.name("c")
      }
      
      rebuild = FALSE
      for(i in 2:length(expr)){
        new_element = expand_lags_lhs(expr[[i]], is_root = FALSE)
        
        if(!identical(new_element, expr[[i]]) && isTRUE(attr(new_element, "lag_expand"))){
          rebuild = TRUE
          new_expr = as.list(expr)
          new_expr[[i]] = new_element
          break
        }
        
        expr[[i]] = new_element
      }
      
      if(rebuild){
        i_start = i + 1
        if(i_start <= length(expr)){
          for(i in i_start:length(expr)){
            new_expr[[i]] = expand_lags_lhs(expr[[i]], is_root = FALSE)
          }
        }
        
        # we rebuild the call
        args_list = lapply(new_expr[-1], 
                           function(x) ifelse(is.character(x), x, deparse_long(x)))
        args_list$sep = ", "
        core = do.call(paste, args_list)
        call_txt = paste0(deparse_long(expr[[1]]), "(", core, ")")
        new_expr = str2lang(call_txt)
        
        return(new_expr)
        
      }
    }
  }
  
  expr
}


fixest_fml_rewriter = function(fml){
  # Currently performs the following
  # - expands lags
  # - protects powers: x^3 => I(x^3)
  #
  # fml = sw(f(y, 1:2)) ~ x1 + l(x2, 1:2) + x2^2 | fe1 | y ~ z::e + g^3
  # fml = y ~ 1 | id + period | l(x_endo, -1:1) ~ l(x_exo, -1:1)
  
  
  lhs = fml[[2]]
  if(is_formula(lhs)){
    left_part = fml[[2]]
    
    lhs = expand_lags_lhs(left_part[[2]])
    rhs = protect_powers_expand_lags(left_part[[3]])
    iv = protect_powers_expand_lags(fml[[3]])
    
    left_part[[2]] = lhs
    left_part[[3]] = rhs
    
    fml[[2]] = left_part
    fml[[3]] = iv
    
  } else {
    lhs = expand_lags_lhs(lhs)
    rhs = protect_powers_expand_lags(fml[[3]])
    
    fml[[2]] = lhs
    fml[[3]] = rhs
  }
  
  all_funs = setdiff(all.vars(fml, functions = TRUE, unique = FALSE),
                     all.vars(fml, functions = FALSE, unique = FALSE))
  isPanel = any(c("l", "d", "f") %in% all_funs)

  res = list(fml = fml, isPanel = isPanel)

  return(res)
}

replace_dot_with_expr = function(expr, replacement){
  
  if(!'.' %in% all.vars(expr)){
    return(expr)
  }
  
  replace_target_with_expr(expr, quote(.), replacement)
}

replace_target_with_expr = function(expr, target, replacement){
  
  if(length(expr) == 1){
    if(is.name(expr) && expr == target){
      return(replacement)
    } else {
      return(expr)
    }
  }
  
  for(i in 2:length(expr)){
    expr[[i]] = replace_target_with_expr(expr[[i]], target, replacement)
  }
  
  return(expr)
}

fixest_upadte_formula = function(fml_new, fml_old){
  
  all_parts_old = fml_split(fml_old)
  all_parts_new = fml_split(fml_new)
  
  #
  # linear no FE
  #
  
  linear_old = all_parts_old[[1]]
  linear_new = all_parts_new[[1]]
  
  lhs_old = linear_old[[2]]
  lhs_new = linear_new[[2]]
  
  rhs_old = linear_old[[3]]
  rhs_new = linear_new[[3]]
  
  fml_main = xpd(lhs = replace_dot_with_expr(lhs_new, lhs_old),
                 rhs = replace_dot_with_expr(rhs_new, rhs_old))
  
  #
  # fe
  #
  
  fe_old = NULL
  if(length(all_parts_old) > 1 && length(all_parts_old[[2]]) == 2){
    fe_old = all_parts_old[[2]][[2]]
  }
  
  fe_new = NULL
  if(length(all_parts_new) > 1 && length(all_parts_new[[2]]) == 2){
    fe_new = all_parts_new[[2]][[2]]
  }
  
  fml_fe = fe_old
  if(!is.null(fe_new)){
    if(identical(fe_new, 0)){
      fml_fe = NULL
    } else {
      fml_fe = replace_dot_with_expr(fe_new, fe_old)
    }
  }
  
  #
  # iv
  #
  
  iv_old = NULL
  if(length(all_parts_old) > 1 && length(all_parts_old[[length(all_parts_old)]]) == 3){
    iv_old = all_parts_old[[length(all_parts_old)]]
  }
  
  iv_new = NULL
  if((length(all_parts_new) > 1 && length(all_parts_new[[length(all_parts_new)]]) == 3) || 
     length(all_parts_new) == 3){
    iv_new = all_parts_new[[length(all_parts_new)]]
    
  }
  
  fml_iv = iv_old
  if(!is.null(iv_new)){
    if(identical(iv_new[[2]], 0)){
      # y ~ . | . | 0
      fml_iv = NULL
    } else if(length(iv_new) == 2 && identical(iv_new[[2]], quote(.))){
      fml_iv = iv_old
    } else {
      if(is.null(iv_old)){
        fml_iv = iv_new
      } else {
        
        endo_old = iv_old[[2]]
        endo_new = iv_new[[2]]
        
        inst_old = iv_old[[3]]
        inst_new = iv_new[[3]]
        
        fml_iv = xpd(lhs = replace_dot_with_expr(endo_new, endo_old),
                     rhs = replace_dot_with_expr(inst_new, inst_old))
        
      }
    }
  }
  
  #
  # combining
  #
  
  res = fml_main
  
  if(!is.null(fml_fe)){
    res = xpd(res, add.after_pipe = fml_fe)
  }
  
  if(!is.null(fml_iv)){
    res = xpd(res, add.after_pipe = fml_iv)
  }
  
  res
  
}

