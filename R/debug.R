#------------------------------------------------------------------------------#
# Author: Laurent R. Bergé
# Created: 2025-07-21
# ~: debug tools
#------------------------------------------------------------------------------#

debug_path = function(filename){
  wd = getOption("fixest_onload_WD", ".")
  file.path(wd, filename)
}

debug_msg = function(...){
  # writes a message in debug.txt
  msg = sma(..., .envir = parent.frame())
  message(msg)
  f = file(debug_path("./../debug.txt"), "a")
  writeLines(msg, f)
  close(f)
}

debug_save = function(path = NULL, up = 0){
  # saves all the variables from the current function
  vars = ls(envir = parent.frame(up + 1))
  
  env = new.env(parent = emptyenv())
  for(v in vars){
    value = try(get(v, envir = parent.frame(up + 1)), silent = TRUE)
    if(!is_error(value)){
      assign(v, value, env)
    }
  }
  
  check_arg(path, "NULL path create")
  if(is.null(path)){
    path = debug_path("./../debug.RData")
  }
  
  save(list = names(env), envir = env, file = path)
}

debug_clear = function(){
  unlink(debug_path("./../debug.RData"))
}

debug_load = function(path = NULL, env = parent.frame(), return_env = FALSE){
  
  check_arg(path, "NULL path create")
  if(is.null(path)){
    path = debug_path("./../debug.RData")
  }
  
  if(return_env){
    env = new.env(parent = emptyenv())
  }
  
  load(path, env)
  
  return(invisible(env))
}

debug_any_variable_different_from_saved = function(path = NULL, ignore = character()){
  # checks if the variables from an environment (a closure) are the same
  # as the one which were previously saved
  # if there is no file => saving is done here
  
  check_arg(path, "NULL path create")
  if(is.null(path)){
    path = debug_path("./../debug.RData")
  }
  
  if(!file.exists(path)){
    debug_save(path, up = 1)
    return(FALSE)
  }
  
  env_old = new.env()
  load(path, envir = env_old)
  
  env_new = parent.frame()
  
  ignore = stringmagic::string_vec(ignore)
  
  all_vars = setdiff(names(env_old), ignore)
  
  is_different = FALSE
  for(v in all_vars){
    x = env_old[[v]]
    y = env_new[[v]]
    if(!is_similar(x, y, obj_name = v, msg = TRUE, do_all = TRUE)){
      is_different = TRUE
    }
  }
  
  return(is_different)
}


is_similar = function(x, y, env_x_done = list(), obj_name = "", msg = FALSE, do_all = FALSE){
  # we focus on the values and not the attributes
  
  if(identical(x, y)){
    return(TRUE)
  }
  
  if(length(x) != length(y)){
    mema("Object {bq ? obj_name} has different lengths: {len ? x} vs {len ? y}", .trigger = msg)
    return(FALSE)
  }
  
  if(!identical(class(x), class(y))){
    mema("Object {bq ? obj_name} has different classes: {enum.bq ? class(x)} vs {enum.bq ? class(y)}", .trigger = msg)
    return(FALSE)
  }
  
  if(is.function(x)){
    # functions are OK, but we check their environment
    if(!identical(formalArgs(x), formalArgs(y))){
      mema("Object {bq ? obj_name} has different arguments:\n",
           "| {enum.bq ? formalArgs(x)}\n", 
           "| vs \n",
           "| {enum.bq ? formalArgs(y)}", .trigger = msg)
      return(FALSE)
    }
    
    env_x = environment(x)
    env_y = environment(y)
    
    if(any(sapply(env_x_done, function(e) identical(e, env_x)))){
      return(TRUE)
    }
    
    obj_name = paste0(obj_name, "$env")
    return(is_similar(env_x, env_y, env_x_done, obj_name, msg, do_all))
  }
  
  x_clean = x
  attributes(x_clean) = NULL
  y_clean = y
  attributes(y_clean) = NULL
  
  if(identical(x_clean, y_clean)){
    return(TRUE)
  }
  
  if(!(is.list(x) || is.environment(x))){
    
    qui_pblm = which(x_clean != y_clean)
    id = qui_pblm[1]
    mema("Object {bq ? obj_name} has {len ? qui_pblm} different value{$s} in position: {enum ? qui_pblm}.\n", 
         "| x[{id}]: {x_clean[id]}\n", 
         "| y[{id}]: {y_clean[id]}", 
         .trigger = msg && length(qui_pblm) > 0)
    
    return(FALSE)
  }
  
  if(is.environment(x)){
    env_x_done[[length(env_x_done) + 1]] = x
  }
  
  for(v in names(x)){
    if(!is_similar(x[[v]], y[[v]], env_x_done, paste0(obj_name, "$", v), msg, do_all)){
      if(!do_all){
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}


compare = function(x, y){
  check_arg(x, y, "list")
  
  if(!identical(names(x), names(y))){
    x_solo = setdiff(names(x), names(y))
    y_solo = setdiff(names(y), names(x))
    
    if(length(x_solo) > 0){
      mema("Object{$s} in x not in y: {enum.bq ? x_solo}.")
    }
    
    if(length(y_solo) > 0){
      mema("Object{$s} in y not in x: {enum.bq ? y_solo}.")
    }
    
  }
  
  vars_ok = intersect(names(x), names(y))
  
  is_similar(x[vars_ok], y[vars_ok], msg = TRUE, do_all = TRUE)
  
  invisible(NULL)
}


