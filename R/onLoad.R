


#
# compatibility issues
#

# Nota:
# - gregexec, used only in format_help(), does not exist in R 3.5.0
#   => we can live without it
#


if(!exists("str2lang", asNamespace("base"), inherits = FALSE)){
  str2lang = function(x){
    parse(text = x, keep.source = FALSE)[[1]]
  }
}

if(!exists("str2expression", asNamespace("base"), inherits = FALSE)){
  str2expression = function(x){
    parse(text = x, keep.source = FALSE)
  }
}

change_defaults = function(fun_name, ...){
  # I don't rewrite the function by attaching the body of the fun_name
  # because I'm a bit wary of namespaces
  #
  # it would work for base R stuff but would be dangerous if extended to
  # functions from imported packages
  # it would also be a problem if some functions I use are conflicted with internal base R funs
  # (although I think they're all exposed)
  #
  # So here it's pretty innocuous, it's a simple rewrite of the call
  # but a bit less efficient
  # 
  
  defaults = list(...)
  for(i in seq_along(defaults)){
    arg_val = parse(text = deparse(defaults[[i]], width.cutoff = 500), keep.source = FALSE)
    if(is.expression(arg_val)){
      arg_val = arg_val[[1]]
    }
    defaults[[i]] = arg_val
  }
  
  fun = parse(text = fun_name, keep.source = FALSE)
  if(is.expression(fun)){
    fun = fun[[1]]
  }
  
  function(...){
    mc = match.call()
    mc_names = names(mc)
    for(arg_name in setdiff(names(defaults), mc_names)){
      mc[[arg_name]] = defaults[[i]]
    }
    
    mc[[1]] = fun
    eval(mc, parent.frame())
  }
  
}

if(is.factor(data.frame(x = "bonjour")$x)){
  data.frame = change_defaults("base::data.frame", stringsAsFactors = FALSE)
  
  as.data.frame = change_defaults("base::as.data.frame", stringsAsFactors = FALSE)
}


#
# startup
#


.onLoad = function(libname, pkgname){
  # setting some options
  
  options("fixest_notes" = TRUE)
  options("fixest_print" = list(type = "table"))
  options("fixest_fl_authorized" = FALSE)

  setFixest_coefplot("all", reset = TRUE)
  setFixest_etable()
  setFixest_vcov(all = "iid")
  setFixest_dict(c("(Intercept)" = "Constant"))

  # nthreads
  if(is_r_check()){
    # limiting the number of threads during R cmd check
    if(requireNamespace("data.table", quietly = TRUE)){
      data.table::setDTthreads(1)
    }
    setFixest_nthreads(1)
  } else {
    setFixest_nthreads()
  }	
  
  # Setup of builtin VCOVs
  vcov_setup()

  # Aliases must come after the VCOV
  create_aliases()

  # To circumvent a peculiar behavior from pkgdown
  fix_pkgwdown_path()

  # register emmeans methods
  if(requireNamespace("emmeans", quietly = TRUE)){
    emmeans::.emm_register(c("fixest", "fixest_multi"), pkgname)
  }

	invisible()
}

.onAttach = function(libname, pkgname) {

  # The startup message mechanism ends up being a bit complex because I try to avoid
  # annoyance as much as possible.
  # I also want to keep track of all the breaking changes so that someone that didn't update for
  # a while is fully informed on how to change his/her old code

  all_startup_msg = list(
    startup_msg(
      "0.13.0",
      "Beware, several BREAKING CHANGES:",
      "- By default the VCOV becomes `iid` for all estimations. Now you need to use the argument `vcov` explicitly to have clustered standard-errors.",
      "  To change the default to the way it was, add this to your .Rprofile:",
      "  `setFixest_vcov(all = \"cluster\", no_FE = \"iid\")`",
      "- Now by default fixed-effect singletons are removed: the coefficients of the estimations will not change but the standard-errors may differ (if you had singletons in your estimations).",
      "  Use the argument `fixef.rm = \"infinite_coef\"`` to get back to the previous results.",
      "  To change the default to the way it was, add this in your .Rprofile:",
      "  `setFixest_estimation(fixef.rm = \"infinite_coef\")`",
      "- `coefplot` and `iplot`: the argument `object` is removed, now all models need to be passed in `...`"
    ),
    startup_msg(
      "0.10.0",
      " - vcov: new argument 'vcov', it replaces the arguments 'se' and 'cluster' in all functions (retro compatibility is ensured)",
      "- BREAKING: arg. 'vcov' now comes after 'data'/'family', hence code using an offset w/t the arg. name may break (just use 'offset = stuff' now)",
      "- function 'dof()' has been renamed into 'ssc()' (i.e. small sample correction)"
    ),
    startup_msg(
      "0.9.0",
      "- in i():",
      "  + the first two arguments have been swapped! Now it's i(factor_var, continuous_var) for interactions",
      "  + argument 'drop' has been removed (put everything in 'ref' now)",
      "- in feglm(): the default family becomes 'gaussian' to be in line with glm(). Hence, for Poisson estimations, please use fepois() instead.",
      trigger = "\\bi\\(|feglm\\("
    )
  )

  fixest_startup_msg = initialize_startup_msg(all_startup_msg)

  if(length(fixest_startup_msg) > 0) {
    msg = c("(Permanently remove the following message with fixest_startup_msg(FALSE).)", 
            fixest_startup_msg)
    packageStartupMessage(fit_screen(paste(msg, collapse = "\n"), 1))
  }

}


is_r_check = function(){
  sum(grepl("_R_CHECK", names(Sys.getenv()), fixed = TRUE)) > 10
}


