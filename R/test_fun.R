#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Jul 10 08:53:10 2020
# ~: basic test function
#----------------------------------------------#



error_catcher = function(expr) tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "try-error"))

test = function(x, y, type = "=", tol = 1e-6){
    mc = match.call()
    IS_Y = TRUE
    if(missing(type) && is.character(y) && y %in% c("err", "=", "~")){
        IS_Y = FALSE
        type = y
    }

    type = switch(type, "error" = "err", type)

    if(type == "=" && is.numeric(x)){
        type = "~"
        tol = 1e-12
    }

    if(type == "err"){
        # we expect an error
        m = error_catcher(x)
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


run_test = function(chunk){
    test_code = readLines("tests/fixest_tests.R")[-(1:17)]

    if(!missing(chunk)){
        qui = which(grepl("^chunk\\(", test_code))
        all_chunks = test_code[qui]
        chunk_names = tolower(gsub(".+\\(\"|\".*", "", all_chunks))
        check_value_plus(chunk, "multi match | integer vector no na", .choices = chunk_names)
        if(is.numeric(chunk)){
            if(any(chunk > length(qui))){
                stop("There are maximum ", length(qui), " chunks.")
            }
            chunk_select = sort(unique(chunk))
        } else {
            chunk_select = which(chunk_names %in% chunk)
        }

        qui = c(qui, length(test_code))

        new_test_code = c()
        for(i in chunk_select){
            new_test_code = c(new_test_code, test_code[qui[i] : (qui[i + 1] - 1)])
        }

        test_code = new_test_code
    }

    setFixest_notes(FALSE)

    parsed_code = parse(text = test_code)

    eval(parsed_code, .GlobalEnv)

    "tests performed successfully"
}

