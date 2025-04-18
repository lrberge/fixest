


collinear_set = function(x, tol = 1e-6, nthreads = getFixest_nthreads(), 
                         out.full = FALSE, scale = TRUE){
  
  # TODO:
  # - handle IVs
  check_arg(x, "numeric matrix | class(fixest, fixest_multi)")
  
  if(inherits(x, "fixest_multi")){
    x = x[[1]]
  }
  
  if(inherits(x, "fixest")){
    if(!is.null(x$fixef_id)){
      # we demean and remove the dep var
      x = demean(x)[, -1]
    } else {
      x = model.matrix(x)
    }
  }
  
  vars = colnames(x)
  K = ncol(x)
  if(is.null(vars)){
    vars = paste0("x", 1:K)
  }
  
  # we need to scale x so that the results can be interpreted
  X = scale(x)
  is_constant = attr(X, "scaled:scale") == 0
  if(any(is_constant)){
    X[, is_constant] = 1
  }
  
  XtX = Xty = crossprod(X)
  
  if(out.full){
    res_full = vector("list", K)
  }
  
  is_collinear = rep(FALSE, K)
  for(i in 1:K){
    est = ols_fit(y = X[, i], X = X[, -i], w = 1, collin.tol = 1e-7, nthreads = nthreads, 
                  xwx = XtX[-i, -i], xwy = XtX[-i, i])
    
    r = est$residuals
    stat = max(abs(r))
    is_col = stat < tol
    
    if(out.full){
      df = data.frame(var = vars[i], collin = is_col, max_abs_resid = stat, 
                      collin.min_norm = est$collin.min_norm)
      df$est = list(est)
      res_full[[i]] = df
    } else {
      if(is_col){
        is_collinear[i] = TRUE
      }
    }
    
  }
  
  if(out.full){
    res = do.call(rbind, res_full)
    return(res)
  }
  
  vars[is_collinear]
}


