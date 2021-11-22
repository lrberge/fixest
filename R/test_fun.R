#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Jul 10 08:53:10 2020
# ~: basic test function
#----------------------------------------------#





test = function(x, y, type = "=", tol = 1e-6){
    mc = match.call()
    IS_Y = TRUE
    if(missing(type) && length(y) == 1 && is.character(y) && y %in% c("err", "=", "~")){
        IS_Y = FALSE
        type = y
    }

    type = switch(type, "error" = "err", type)

    if(type == "=" && is.numeric(x)){
        type = "~"
        if(missing(tol)) tol = 1e-12
    }

    if(type == "err"){
        # we expect an error
        m = tryCatch(x, error = function(e) structure(conditionMessage(e), class = "try-error"))
        if(!"try-error" %in% class(m)){
            stop("Expected an error that did not occur.")
        } else if(IS_Y && !grepl(tolower(y), tolower(m), fixed = TRUE)){
            stop("This is an error, but the messages don't match:\nEXPECTED: \n", y, "\nACTUAL: \n", x)
        }
    } else if(length(x) == 0){
        stop("Argument 'x' is of length 0. This is not allowed.")

    } else if(type %in% c("equal", "=")){
        if(length(x) != length(y)){
            stop("Lengths differ: EXPECTED: ", length(y), "\n",
                 "                    ACTUAL: ", length(x))

        } else if(anyNA(x)){
            if(!all(is.na(x) == is.na(y))){
                if(sum(is.na(x)) != sum(is.na(y))){
                    stop("Number of NA values differ: EXPECTED: ", sum(is.na(y)), "\n",
                         "                                ACTUAL: ", sum(is.na(x)))
                }

                i = which(is.na(x) != is.na(y))[1]

                na_y = 1 + is.na(y)[i]
                info = c("none", "one NA")

                stop("Position of the NA values differ. EXPECTED: ", info[na_y], " in position ", i, "\n",
                     "                                      ACTUAL: ", info[3 - na_y], ".")

            } else if(!all(is_na <- is.na(x)) && any(qui_pblm <- x[!is_na] != y[!is_na])){

                if(all(qui_pblm)){
                    if(length(x) == 1){
                        stop("Non-NA values differ: EXPECTED: ", y[!is_na], "\n",
                             "                          ACTUAL: ", x[!is_na])
                    } else {
                        stop("All non-NA values differ: 1st elem.: EXPECTED: ", y[!is_na][1], "\n",
                             "                                         ACTUAL: ", x[!is_na][1])
                    }
                } else {
                    n = sum(qui_pblm)
                    i = which(qui_pblm)[1]
                    stop(n, " non-NA value", plural(n, "s.differ"), ": ", n_th(i), " elem.: \nEXPECTED: ", y[!is_na][i], "\nACTUAL: ", x[!is_na][i])
                }
            }

        } else if(anyNA(y)){
            stop("Number of NA values differ: EXPECTED: ", sum(is.na(y)), "\n",
                 "                                ACTUAL: ", sum(is.na(x)))

        } else if(!all(x == y)){
            qui_pblm = x != y

            if(all(qui_pblm)){
                if(length(x) == 1){
                    stop("Values differ: EXPECTED: ", y, "\n",
                         "                   ACTUAL: ", x)
                } else {
                    stop("All values differ: 1st elem.: EXPECTED: ", y[1], "\n",
                         "                                  ACTUAL: ", x[1])
                }
            } else {
                n = sum(qui_pblm)
                i = which(qui_pblm)[1]
                stop(n, " value", plural(n, "s.differ"), ": ", n_th(i), " elem.: \nEXPECTED: ", y[i], "\nACTUAL: ", x[i])
            }
        }
    } else if(type %in% c("~", "approx")){
        if(length(x) != length(y)){
            stop("Lengths differ: EXPECTED: ", length(y), "\n",
                 "                    ACTUAL: ", length(x))

        } else if(anyNA(x)){
            if(!all(is.na(x) == is.na(y))){
                if(sum(is.na(x)) != sum(is.na(y))){
                    stop("Number of NA values differ: EXPECTED: ", sum(is.na(y)), "\n",
                         "                                ACTUAL: ", sum(is.na(x)))
                }

                i = which(is.na(x) != is.na(y))[1]

                na_y = 1 + is.na(y)[i]
                info = c("none", "one NA")

                stop("Position of the NA values differ. EXPECTED: ", info[na_y], " in position ", i, "\n",
                     "                                      ACTUAL: ", info[3 - na_y], ".")

            } else if(max(abs(x - y), na.rm = TRUE) > tol){

                stop("Difference > tol: Max abs. diff: ", max(abs(x - y), na.rm = TRUE), " (in position ", which.max(abs(x - y)), ")")
            }
        } else if(anyNA(y)){
            stop("Number of NA values differ: EXPECTED: ", sum(is.na(y)), "\n",
                 "                                ACTUAL: ", sum(is.na(x)))

        } else if(max(abs(x - y)) > tol){
            stop("Difference > tol: Max abs. diff: ", max(abs(x - y)), " (in position ", which.max(abs(x - y)), ")")
        }
    }

}


chunk = function(x) cat(toupper(x), "\n\n")


run_test = function(chunk, from){

    test_code = readLines("tests/fixest_tests.R")

    i = grepl("library(fixest)", test_code, fixed = TRUE)
    test_code[i] = "devtools::load_all(export_all = FALSE)"

    # A) adding the line numbers

    lines = paste0("LINE_COUNTER = ", seq_along(test_code))

    code_LINE = rep(lines, each = 2)
    code_LINE[seq_along(test_code)  * 2] = test_code

    # We remove the counters for the lines with open parenthesis, like in
    # sum(x,
    #     y)
    # since this leads to parsing errors

    i_open = which(grepl(",[[:blank:]]*$", code_LINE))

    if(length(i_open)) i_open = i_open + 1

    # We remove the counters just before closing brackets => otherwise equivalent to a return statement!
    i_closing_bracket = which(grepl("^[[:blank:]]*\\}", code_LINE))

    if(length(i_closing_bracket)) i_closing_bracket = i_closing_bracket - 1

    # removing
    i_rm = c(i_open, i_closing_bracket)
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
            check_value_plus(from, "match | integer scalar no na", .choices = chunk_names)

            if(is.numeric(from)){
                if(any(from > n_chunks)){
                    stop("There are maximum ", n_chunks, " chunks.")
                }
                chunk_select = from:n_chunks
            } else {
                chunk_select = which(chunk_names %in% from):n_chunks
            }

        } else {
            check_value_plus(chunk, "multi match | integer vector no na", .choices = chunk_names)

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

    parsed_code = parse(text = test_code)

    env = new.env()
    assign("INSIDE_LOOP", FALSE, env)

    my_eval = try(eval(parsed_code, env))

    # E) Message

    if("try-error" %in% class(my_eval)){
        line_fail = get("LINE_COUNTER", env)
        inside_loop = get("INSIDE_LOOP", env)

        message("Test failed at line: ", line_fail)
        if(inside_loop){
            index_loop = get("INDEX_LOOP", env)
            index_names = names(index_loop)
            index_values = unlist(index_loop)
            msg = paste0(index_names, ":", index_values, collapse = ", ")
            message("Loop values: ", msg, ".")
        }
    } else {
        print("tests performed successfully")
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


# To circumvent an Rstudio/Google Drive bug preventing
# the source files to remain open when closing/opening the project
open_all = function(){

    # Opening R files

    R_not_open = c("R/alias_generator.R", "R/RcppExports.R", "R/Deprecated_funs.R", "R/index",
                   "R/etable_aliases.R", "R/ML_Families.R", "R/VCOV_aliases.R", "R/xaxis.R")

    R_files = list.files("R/", full.names = TRUE)

    for(f in R_files){
        if(f %in% R_not_open) next
        eval(str2lang(.dsb("rstudioapi::navigateToFile('.[f]', moveCursor = FALSE)")))
    }

    # Opening extra files

    extra_files = c("../PROBLEMS.R", "../Problems.Rmd", "../Development.R",
                    "../social.rmd", "NEWS.md", "DESCRIPTION", "NAMESPACE", "tests/fixest_tests.R",
                    list.files("vignettes/", pattern = "Rmd$", full.names = TRUE),
                    "../todo.txt")

    for(f in extra_files){
        eval(str2lang(.dsb("rstudioapi::navigateToFile('.[f]', moveCursor = FALSE)")))
    }


}












































