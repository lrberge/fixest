#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Jul 10 08:53:10 2020
# ~: basic test function
#----------------------------------------------#


update_test_counter = function(reset = FALSE){
  n = getOption("TEST_COUNTER", default = 0)
  if(reset){
    options(TEST_COUNTER = 0)
  } else {
    options(TEST_COUNTER = n + 1)
  }  
}

get_test_counter = function(){
  getOption("TEST_COUNTER", default = 0)
}

test = function(x, y, type = "=", tol = 1e-6){
  
  update_test_counter()

  mc = match.call()
  IS_Y = TRUE
  if(missing(type) && length(y) == 1 && is.character(y) && y %in% c("err", "=", "~", "warn")){
    IS_Y = FALSE
    type = y
  }

  type = switch(type, "error" = "err", "warning" = "warn", type)

  plural_core = function(PLURAL, type, s, verb = FALSE, past = FALSE){

    if(!missing(type)){
      args = strsplit(type, "\\.")[[1]]
      s = "s" %in% args
      past = "past" %in% args
      verb = setdiff(args, c("s", "past"))
      if(length(verb) == 0){
        verb = FALSE
      } else {
        verb = verb[1]
      }
    }

    if(isFALSE(verb)){
      res = ifelse(PLURAL, "s", "")
    } else {

      if(verb %in% c("be", "are")){
        verb = "is"
      } else if(verb == "have"){
        verb = "has"
      } else if(verb == "does"){
        verb = "do"
      } else if(verb %in% c("do not", "does not", "don't", "doesn't")){
        verb = "do not"
      } else if(verb %in% c("is not", "are not", "isn't", "aren't")){
        verb = "is not"
      } else if(verb %in% c("was", "were")){
        verb = "is"
        past = TRUE
      }

      if(past){
        if(verb %in% c("is", "is not", "has", "do", "do not")){
          verb_format = switch(verb, is = ifelse(!PLURAL, "was", "were"), "is not" = ifelse(!PLURAL, "wasn't", "weren't"), has = "had", do = "did", "do not" = "didn't")
        } else {
          verb_format = paste0(verb, "ed")
        }
      } else {
        if(verb %in% c("is", "is not", "has", "do", "do not")){
          verb_format = switch(verb, is = ifelse(!PLURAL, "is", "are"), "is not" = ifelse(!PLURAL, "isn't", "aren't"), has = ifelse(!PLURAL, "has", "have"), do = ifelse(!PLURAL, "does", "do"), "do not" = ifelse(!PLURAL, "doesn't", "don't"))
        } else {
          verb_format = ifelse(PLURAL, verb, paste0(verb, "s"))
        }
      }

      if(!missing(s) && isTRUE(s)){
        res = paste0(ifelse(PLURAL, "s ", " "), verb_format)
      } else {
        res = verb_format
      }

    }

    res
  }


  plural = function(x, type, s, verb = FALSE, past = FALSE){
    # adds s if x > 1, can alternatively add a verb
    PLURAL = x[1] > 1
    plural_core(PLURAL, type, s, verb, past)
  }

  n_th = function(n){
    if(length(n) == 0) return("")
    dict = c("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "nineth", "tenth", "eleventh", "twelfth", "thirteenth")
    res = as.character(n)
    qui = n <= 13
    res[qui] = dict[n[qui]]
    if(any(!is.na(qui) & !qui)){
      other = n[!qui]
      rest = other %% 10
      rest[rest == 0 | rest >= 4] = 4
      postfix = c("st", "nd", "rd", "th")
      other = paste0(other, postfix[rest])
      res[!qui] = other
    }

    res
  }

  if(type == "=" && is.numeric(x)){
    type = "~"
    if(missing(tol)) tol = 1e-12
  }
  
  if(type  == "warn"){
    # we expect an warning
    m = tryCatch(x, warning = function(w) structure(conditionMessage(w), class = "try-warning"))
    if(!"try-warning" %in% class(m)){
      stop("Expected a warning that did not occur.")
    } else if(IS_Y && !grepl(tolower(y), tolower(m), fixed = TRUE)){
      stop("This is a warning, but the messages don't match:\nEXPECTED: \n", y, "\n..ACTUAL: \n", x)
    }
  } else if(type == "err"){
    # we expect an error
    m = tryCatch(x, error = function(e) structure(conditionMessage(e), class = "try-error"))
    if(!"try-error" %in% class(m)){
      stop("Expected an error that did not occur.")
    } else if(IS_Y && !grepl(tolower(y), tolower(m), fixed = TRUE)){
      stop("This is an error, but the messages don't match:\nEXPECTED: \n", y, "\n..ACTUAL: \n", x)
    }
  } else if(length(x) == 0){
    stop("Argument 'x' is of length 0. This is not allowed.")

  } else if(type %in% c("equal", "=")){
    if(length(x) != length(y)){
      stop("Lengths differ: ", 
            "\nEXPECTED: ", length(y), "\n",
            "\n..ACTUAL: ", length(x))

    } else if(anyNA(x)){
      if(!all(is.na(x) == is.na(y))){
        if(sum(is.na(x)) != sum(is.na(y))){
          stop("Number of NA values differ: ",
          "\nEXPECTED: ", sum(is.na(y)), "\n",
          "\n..ACTUAL: ", sum(is.na(x)))
        }

        i = which(is.na(x) != is.na(y))[1]

        na_y = 1 + is.na(y)[i]
        info = c("none", "one NA")

        stop("Position of the NA values differ. ", 
              "\nEXPECTED: ", info[na_y], " in position ", i, "\n",
              "\n..ACTUAL: ", info[3 - na_y], ".")

      } else if(!all(is_na <- is.na(x)) && any(qui_pblm <- x[!is_na] != y[!is_na])){

        if(all(qui_pblm)){
          if(length(x) == 1){
              stop("Non-NA values differ: ", 
                    "\nEXPECTED: ", y[!is_na], "\n",
                    "\n..ACTUAL: ", x[!is_na])
          } else {
              stop("All non-NA values differ: 1st elem.: ", 
                    "\nEXPECTED: ", y[!is_na][1], "\n",
                    "\n..ACTUAL: ", x[!is_na][1])
          }
        } else {
          n = sum(qui_pblm)
          i = which(qui_pblm)[1]
          stop(n, " non-NA value", plural(n, "s.differ"), ": ", n_th(i), " elem.: \nEXPECTED: ", y[!is_na][i], "\n..ACTUAL: ", x[!is_na][i])
        }
      }

    } else if(anyNA(y)){
      stop("Number of NA values differ: ", 
            "\nEXPECTED: ", sum(is.na(y)), "\n",
            "\n..ACTUAL: ", sum(is.na(x)))

    } else if(!all(x == y)){
      qui_pblm = x != y

      if(all(qui_pblm)){
        if(length(x) == 1){
          stop("Values differ: ", 
                "\nEXPECTED: ", y, "\n",
                "\n..ACTUAL: ", x)
        } else {
          stop("All values differ: 1st elem.: ", 
                "\nEXPECTED: ", y[1], "\n",
                "\n..ACTUAL: ", x[1])
        }
      } else {
        n = sum(qui_pblm)
        i = which(qui_pblm)[1]
        stop(n, " value", plural(n, "s.differ"), ": ", n_th(i), " elem.: \nEXPECTED: ", y[i], "\n..ACTUAL: ", x[i])
      }
    }
  } else if(type %in% c("~", "approx")){
    if(length(x) != length(y)){
      stop("Lengths differ: ", 
            "\nEXPECTED: ", length(y), "\n",
            "\n..ACTUAL: ", length(x))

    } else if(anyNA(x)){
      if(!all(is.na(x) == is.na(y))){
        if(sum(is.na(x)) != sum(is.na(y))){
          stop("Number of NA values differ: ", 
                "\nEXPECTED: ", sum(is.na(y)), "\n",
                "\n..ACTUAL: ", sum(is.na(x)))
        }

        i = which(is.na(x) != is.na(y))[1]

        na_y = 1 + is.na(y)[i]
        info = c("none", "one NA")

        stop("Position of the NA values differ. ", 
              "\nEXPECTED: ", info[na_y], " in position ", i, "\n",
              "\n..ACTUAL: ", info[3 - na_y], ".")

      } else if(max(abs(x - y), na.rm = TRUE) > tol){
        stop("Difference > tol: Max abs. diff: ", max(abs(x - y), na.rm = TRUE), " (in position ", which.max(abs(x - y)), ")")
      }
    } else if(anyNA(y)){
      stop("Number of NA values differ: ", 
            "\nEXPECTED: ", sum(is.na(y)), "\n",
            "\n..ACTUAL: ", sum(is.na(x)))

    } else if(max(abs(x - y)) > tol){
      stop("Difference > tol: Max abs. diff: ", max(abs(x - y)), " (in position ", which.max(abs(x - y)), ")")
    }
  }

}


chunk = function(x) message(toupper(x), "\n")

test_contains = function(x, pattern){
  update_test_counter()
  
  if(length(x) > 1){
    stop("Internal error: `test_contains` only works with vectors of length 1.")
  }
  
  if(!grepl(pattern, x)){
    stop("The pattern `", pattern, "` was not found in the following string:", 
         "\n\n", x)
  }
  
}

test_err_contains = function(x, pattern){
  update_test_counter()
  
  err = try(x, silent = TRUE)
  if(!inherits(err, "try-error")){
    x_dp = deparse(substitute(x))[1]
    stop("The expression should lead to an error.", 
            "\nPROBLEM: `", x_dp, "` is not an error.")
  }
  
  if(!grepl(pattern, err)){
    stop("The pattern `", pattern, "` was not found in the following error message:", 
         "\n\n", err)
  }
  
}

run_tests = function(chunk, from = 1, source = FALSE, use_devtools = TRUE){
  update_test_counter(reset = TRUE)
  
  wd = getwd()
  
  location = "tests"
  if(grepl("/tests$", wd)){
    location = "."
  }
  
  dir_tests = list.dirs(location, recursive = FALSE)
  if(length(dir_tests) > 0){
    # either all files in the dirs
    all_files = list.files(dir_tests, pattern = ".R$", full.names = TRUE)
  } else {
    # either the files in the main test dir
    all_files = list.files(location, pattern = ".R$", full.names = TRUE)
  }
  
  if(length(all_files) == 0){
    stop("No test file was found in the folder './tests/'.")
  }
  
  if(isTRUE(source)){
    assign("chunk", get("chunk", mode = "function"), parent.frame())
    assign("test", get("test", mode = "function"), parent.frame())
    assign("test_contains", get("test_contains", mode = "function"), parent.frame())
    assign("test_err_contains", test_err_contains, parent.frame())
    
    for(f in all_files){
      source(f)
    }
    return(NULL)
  }  

  test_code = c()
  lines_per_file = c()
  file_names = c()
  offset = 2
  for(f in all_files){
    file_names = c(file_names, f)
    f_open = file(f, encoding = "UTF-8")
    new_code = readLines(f_open)
    close(f_open)
    test_code = c(test_code, new_code)
    lines_per_file = c(lines_per_file, length(new_code))
  }

  pkg = gsub(".+/", "", normalizePath(".", "/"))

  # we don't load the package with library
  i_pkg = grepl(paste0("library(", pkg, ")"), test_code, fixed = TRUE)
  test_code[i_pkg] = ""

  # A) adding the line numbers

  lines = paste0("LINE_COUNTER = ", seq_along(test_code))

  code_LINE = rep(lines, each = 2)
  code_LINE[seq_along(test_code)  * 2] = test_code

  # We remove the counters for the lines with open parenthesis, like in
  # sum(x,
  #     y)
  # dt[, a :=
  #        32 + b]
  # since this leads to parsing errors

  i_open = which(grepl("[([,=][[:blank:]]*$", code_LINE))

  if(length(i_open)) i_open = i_open + 1

  # We remove the counters just before closing brackets => otherwise equivalent to a return statement!
  i_closing_bracket = which(grepl("^[[:blank:]]*[]})]", code_LINE))

  if(length(i_closing_bracket)) i_closing_bracket = i_closing_bracket - 1
  
  # we remove the counter for empty lines or commented lines
  i_empty = which(grepl("^\\s*$|^\\s*#", code_LINE))

  if(length(i_empty)) i_empty = i_empty - 1

  # removing
  i_rm = sort(unique(c(i_open, i_empty, i_closing_bracket)))
  if(length(i_rm) > 0){
      code_LINE = code_LINE[-(i_rm)]
  }

  # B) Adding the FOR loops

  bracket_close = which(grepl("^\\}", code_LINE))
  for_start = which(grepl("^for\\(", code_LINE))
  for_end = c()
  for(i in for_start){
      for_end = c(for_end, bracket_close[which.max(bracket_close > i)])
  }

  n = length(code_LINE)
  my_rep = rep(1, n)
  my_rep[c(for_start, for_end)] = 2

  code_LINE_FOR = rep(code_LINE, my_rep)
  n_loop = length(for_start)
  code_LINE_FOR[for_start + 2 * (0:(n_loop - 1))] = "INSIDE_LOOP = TRUE ; INDEX_LOOP = list()"
  code_LINE_FOR[for_end + 2 * 1:n_loop] = "INSIDE_LOOP = FALSE"

  qui_for = which(grepl("\\bfor\\(", code_LINE_FOR))
  index = gsub("^.*for\\(([^ ]+).+", "\\1", code_LINE_FOR[qui_for])

  n = length(code_LINE_FOR)
  my_rep = rep(1, n)
  my_rep[qui_for] = 2

  code_LINE_FOR = rep(code_LINE_FOR, my_rep)
  n_loop = length(qui_for)
  code_LINE_FOR[qui_for + 1 + (0:(n_loop - 1))] = paste0("INDEX_LOOP[[\"", index, "\"]] = ", index)

  test_code = code_LINE_FOR

  # C) Chunk selection
  
  if(!missing(chunk) || !missing(from)){

    qui = which(grepl("^chunk\\(", test_code))
    all_chunks = test_code[qui]
    chunk_names = tolower(gsub(".+\\(\"|\".*", "", all_chunks))
    n_chunks = length(qui)

    if(!missing(from)){
      if(is.character(from)){
        from = check_set_options(from, chunk_names, "test_fun")[1]
      } else {
        stopifnot(length(from) == 1, is.numeric(from), from %% 1 == 0)
      }

      if(is.numeric(from)){
        if(any(from > n_chunks)){
          stop("There are maximum ", n_chunks, " chunks.")
        }
        chunk_select = from:n_chunks
      } else {
        chunk_select = which(chunk_names %in% from):n_chunks
      }

    } else {
      if(is.character(chunk)){
        chunk = check_set_options(chunk, chunk_names, "test_fun")
      } else {
        stopifnot(length(chunk) == 1, is.numeric(chunk), chunk %% 1 == 0)
      }

      if(is.numeric(chunk)){
        if(any(chunk > n_chunks)){
          stop("There are maximum ", n_chunks, " chunks.")
        }
        chunk_select = sort(unique(chunk))
      } else {
        chunk_select = which(chunk_names %in% chunk)
      }
    }

    qui = c(qui, length(test_code))

    new_test_code = c()
    for(i in chunk_select){
      new_test_code = c(new_test_code, test_code[qui[i] : (qui[i + 1] - 1)])
    }

    # We add the preamble
    preamble = test_code[1:(qui[1] - 1)]

    test_code = c(preamble, new_test_code)
  }

  # D) Evaluation

  loader = if(use_devtools) "devtools::load_all(export_all = FALSE)" else ""
  test_code = c(loader, test_code)

  parsed_code = parse(text = test_code)

  env = new.env()
  assign("INSIDE_LOOP", FALSE, env)
  assign("chunk", get("chunk", mode = "function"), env)
  assign("test", test, env)
  assign("test_contains", get("test_contains", mode = "function"), parent.frame())
  assign("test_err_contains", test_err_contains, env)

  my_eval = try(eval(parsed_code, env), silent = TRUE)

  # E) Message

  if("try-error" %in% class(my_eval)){
    line_fail = get("LINE_COUNTER", env)
    inside_loop = get("INSIDE_LOOP", env)
    
    nfiles = length(all_files)
    line_start = c(0, cumsum(lines_per_file)[-nfiles])
    line_end = cumsum(lines_per_file)
    
    index = which(line_fail >= line_start & line_fail <= line_end)
    my_file_full = file_names[index]
    my_file_short = gsub(".+/", "", my_file_full)
    file_info = paste0("In file \"", my_file_short, "\"\n==> ")
    
    line_in_file = line_fail - line_start[index]
    
    message("\n---------------------------\n")

    message(file_info, "Test failed at line: ", line_in_file)
    
    i = which(test_code == paste0("LINE_COUNTER = ", line_fail))
    my_line = test_code[i + 1]
    message("\nHere is the line:\n", line_in_file, "| ", my_line, "\n")
    if(inside_loop){
      index_loop = get("INDEX_LOOP", env)
      index_names = names(index_loop)
      index_values = unlist(index_loop)
      msg = paste0(index_names, " = ", index_values, collapse = "; ")
      message("Loop values: ", msg, " (all currently assigned in GlobalEnv).\n")
    }

    # We assign the variables to the global env to facilitate debugging
    for(var in names(env)){
      assign(var, get(var, env), parent.frame())
    }
    
    error_clean = gsub("^[^:]+: *", "", my_eval)
    message("Here is the error:\n", error_clean)
    
    if(interactive() && use_devtools){
      command = paste0("rstudioapi::navigateToFile(\"", my_file_full, "\", ", line_in_file, ")")
      command = str2lang(command)
      
      eval(command)
    }    

  } else {
    n_tests = get_test_counter()
    message("tests performed successfully (", n_tests, ")")
  }
}


non_ascii = function(folder = "R"){
  # Finds non ascii characters lurking

  all_files = list.files(folder, full.names = TRUE)
  all_R_text = lapply(all_files, function(x) readLines(x, encoding = "UTF-8"))

  i_non_ascii = which(sapply(all_R_text, function(x) any(grepl("[^ -~\t]", x))))
  n = length(i_non_ascii)

  if(n == 0) return("No non-ASCII character found.")

  cat("Tip: Type, e.g., 1750G to go to the line in VIM\n\n")

  for(id in seq(n)){
    i = i_non_ascii[id]
    cat("File: ", gsub("R/", "", all_files[i]), "\n")

    text = all_R_text[[i]]
    all_lines = which(grepl("[^ -~\t]", text))
    for(line in all_lines){
      cat("-> line ", sfill(line, max(nchar(all_lines))), ":\n===|", text[line], "\n")
      cat("===|", gsub("[^ -~\t]", "__HERE__", text[line]), "\n")

      if(line != tail(all_lines, 1)) cat("\n")
    }

    if(id < n) cat("\n ---------- \n\n")
  }
}


check_set_options = function(x, options, op = NULL, free = FALSE, case = FALSE){
  # x: always a character vector
  # options: the options to match
  if(length(x) == 0) return(x)

  n = length(x)
  res = x
  for(i in 1:n){
    v = x[i]

    pm = pmatch(v, options)
    if(is.na(pm) && !case){
      pm = pmatch(tolower(v), tolower(options))
    }
    
    if(is.na(pm) && !free){
      # absence of a match
      stop("The option `", v, "` is not valid for the current operation.\n",
           "FYI the option available are ", paste(options, collapse = ", "), ".")
    }

    if(!is.na(pm)){
      res[i] = options[pm]
    }
  }

  res
}






























