

local({
  # 0-dep function to silently load (and install if needed) packages at startup
  load_pkg = function(...){
    mc = match.call(expand.dots = FALSE)
    mc_dots = mc[["..."]]
    
    for(i in seq_along(mc_dots)){
      pkg_i = mc_dots[[i]]
      pkg_name = if(is.character(pkg_i)) pkg_i else deparse(pkg_i)
      
      if(!requireNamespace(pkg_name, quietly = TRUE)){
        ok = try(install.packages(pkg_name, repos = "https://cloud.r-project.org/"))
        if(inherits(ok, "try-error")){
          stop("Could not install package `", pkg_name, "`. Please fix the problem manually.")
        }
      }
      
      suppressWarnings(suppressPackageStartupMessages(library(pkg_name, character.only = TRUE)))
    }
  }
  
  # add packages to load here: eg load_pkg(fixest, stringmagic)
  
  # Replacing View with datadive
  if(requireNamespace("datadive", quietly = TRUE)){
    assign("View", datadive::app_explore, .GlobalEnv)
  }
  
  # VSCode specific
  # adapted from: https://github.com/REditorSupport/vscode-R/wiki/Interacting-with-R-terminals
  # => no error if file is missing
  if (interactive() && Sys.getenv("RSTUDIO") == "") {
    Sys.setenv(TERM_PROGRAM = "vscode")
    HOME = if(.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
    my_file = file.path(Sys.getenv(HOME), ".vscode-R", "init.R")
    if(file.exists(my_file)){
      source(my_file)
      options(vsc.rstudioapi = TRUE)
    }
  }
})


