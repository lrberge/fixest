

# adapted from: https://github.com/REditorSupport/vscode-R/wiki/Interacting-with-R-terminals
if (interactive() && Sys.getenv("RSTUDIO") == "") {
  Sys.setenv(TERM_PROGRAM = "vscode")
  HOME = if(.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
  my_file = file.path(Sys.getenv(HOME), ".vscode-R", "init.R")
  if(file.exists(my_file)){
    source(my_file)
    options(vsc.rstudioapi = TRUE)
  }
}


