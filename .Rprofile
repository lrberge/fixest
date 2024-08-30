

local({
  
  # Replacing View with datadive
  if(requireNamespace("datadive", quietly = TRUE)){
    assign("View", datadive::app_explore, .GlobalEnv)
  }
  
  if(requireNamespace("microbenchmark", quietly = TRUE)){
    assign("mb", microbenchmark::microbenchmark, .GlobalEnv)
  }
  
  # VSCode specific
  # adapted from: https://github.com/REditorSupport/vscode-R/wiki/Interacting-with-R-terminals
  # => no error if file is missing + we hard check vscode
  if(interactive() && Sys.getenv("RSTUDIO") == "" && any(grepl("VSCODE", names(Sys.getenv())))){
    Sys.setenv(TERM_PROGRAM = "vscode")
    HOME = if(.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
    my_file = file.path(Sys.getenv(HOME), ".vscode-R", "init.R")
    if(file.exists(my_file)){
      source(my_file)
      options(vsc.rstudioapi = TRUE)
    }
  }
})


