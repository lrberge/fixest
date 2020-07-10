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
        tol = 1e-15
    }

    if(type == "err"){
        # we expect an error
        m = error_catcher(x)
        if(!"try-error" %in% class(m)){
            stop("Expected an error that did not occur.")
        } else if(IS_Y && !grepl(tolower(y), tolower(m), fixed = TRUE)){
            stop("This is an error, but the messages don't match:\nEXPECTED: \n", y, "\nACTUAL: \n", x)
        }
    } else if(type %in% c("equal", "=")){
        if(length(x) != length(y)){
            stop("Lengths differ: EXPECTED: ", length(y), "\nACTUAL: ", length(x))
        } else if(!all(x == y)){
            qui_pblm = x != y

            if(all(qui_pblm)){
                if(length(x) == 1) stop("Values differ: EXPECTED: ", y, "\nACTUAL: ", x)
                else stop("All values differ: 1st elem.: EXPECTED: ", y[1], "\nACTUAL: ", x[1])
            } else {
                n = sum(qui_pblm)
                i = which(qui_pblm)[1]
                stop(n, " value", plural(n, "s.differ"), ": ", n_th(i), " elem.: \nEXPECTED: ", y[i], "\nACTUAL: ", x[i])
            }
        }
    } else if(type %in% c("~", "approx")){
        if(length(x) != length(y)){
            stop("Lengths differ: EXPECTED: ", length(y), "\nACTUAL: ", length(x))
        } else if(max(abs(x - y)) > tol){
            stop("Difference > tol: Max abs. diff: ", max(abs(x - y)))
        }
    }

}

