


print.dsb = function(x, ...){
    cat(x, sep = "\n")
}



dsb = function(..., frame = parent.frame(), sep = "", concat = FALSE, nest = TRUE){
    check_arg(concat, nest, "logical scalar")
    check_arg(sep, "character scalar")
    check_arg(frame, "class(environment)")
    check_arg(..., "vector len(1)")

    if(...length() == 0) return("")

    sc = sys.call()
    if(identical(sc[[2]], "--help")){
        msg = c(
            "Welcome to dsb help\nUsage: dsb(s) with 's' a character string",
            " ",
            "BASIC usage ------------|",
            "    dsb evaluates anything in '.[]' and inserts it in 's'.",
            '    Ex: if x = "John", then dsb("Hi .[x]!") -> "Hi John!"',
            " ",
            "STRING OPERATIONS ------|",
            "    Each .[] instance supports one or more string operations. The syntax is .['arg'op?x], with:",
            "      - 'arg' a quoted string used as argument,",
            "      - op an operator code,",
            "      - ? (question mark) marks the end of the operation(s),",
            "      - x a variable to be evaluated.",
            '    Ex: dsb(".[\' + \'c?1:3] = 6") -> "1 + 2 + 3 = 6". 1:3 is collapsed (c) with \' + \'.',
            "",
            "    Using ! instead of ? applies the operation to the *verbatim* of the expression.",
            '    Ex: dsb(".[\': => 2\'r!1:3] = 6") -> "123 = 6".',
            "        In the string '1:3', ':' is replaced (r) with '2'.",
            "",
            "    Operations can be chained using comma-separation. The syntax is: .['s1'op1, 's2'op2?x]",
            "    Evaluations are from left to right.",
            '    Ex: dsb(".[\': => 2\'r, \'\'s, \' + \'c!1:3] = 6") -> "1 + 2 + 3 = 6',
            "        1) '1:3'            -> ':' is replaced (r) with '2'  -> '123',",
            "        2) '123'            -> is split (s) with ''          -> c('1', '2', '3')",
            "        3) c('1', '2', '3') -> is collapsed (c) with ' + '   -> '1 + 2 + 3'",
            "",
            "    Nesting works, but only in verbatim components.",
            "    Ex: x = c(\"Doe\", \"Smith\")",
            "        dsb(\"Hi .[' and 'c!John .[x]]\") -> \"Hi John Doe and John Smith\"",
            "",
            "    If nest = TRUE (default), the initial string is implicitly embedded in .[] ",
            "  and considered as verbatim. This means that operations can be applied from the start.",
            "    Ex: dsb(\"', 's!a1, b3, d8\") -> c(\"a1\", \"b3\", \"d8\")",
            "",
            "    Operators have default values, so the quoted argument is optional.",
            '    Ex: dsb("c?1:3") -> "123". 1:3 is collapsed (c) with \'\', its default.',
            "",
            "OPERATORS --------------|",
            "    Below is the list of operators and their default argument when relevant.",
            "    OPERATOR DEFAULT                   VALUE ",
            "",
            "    s        ', '                      split, fixed = TRUE",
            "    S        ', *'                     split, perl = TRUE",
            "    x        '[[:alnum:]]+'            extracts the first pattern, perl = TRUE",
            "    X        '[[:alnum:]]+'            extracts all patterns, perl = TRUE",
            "",
            "  Collapse:",
            "    - collapses a vector with paste(x, collapse = 's'),",
            "    - use a double pipe to apply a special collapse to the last values,",
            "    - the syntax is 's1||s2', s2 will be applied to the last 2 values,",
            "    - .[', || and 'c!1:3] -> \"1, 2 and 3\".",
            "    c        ''                        collapse",
            "    C        ', '                      collapse",
            "",
            "  Replication:",
            "    - use 5* to replicate 5 times,",
            "    - this is the only operator that does not require a quoted argument,",
            "    - but quoted numbers also work.",
            "    *        '1'                       replicates n times",
            "    *c       '1'                       replicates n times, then collapses with ''",
            "",
            "  Replacement: ",
            "    - the syntax is 'old => new', 'old=>new' or 'old_=>_new',",
            "    - if new is missing, it is considered as the empty string,",
            "    - the default for 'R' is: removing trailing spaces.",
            "    r        '\\n'                      replacement, fixed = TRUE ",
            "    R        '^[ \\t\\r\\n]+|[ \\t\\r\\n]+$' replacement, perl = TRUE ",
            "",
            "  Operators without arguments:",
            "    u                                  puts the first letter of the string to uppercase",
            "    U                                  puts the complete string to uppercase",
            "    l, L                               puts the complete string to lowercase",
            "    q                                  adds single quotes",
            "    Q                                  adds double quotes",
            "    f                                  applies format(x)",
            "    F                                  applies format(x, justify = 'right')",
            "",
            "  sprintf formatting:",
            "    - applies a formatting via sprintf,",
            "    - .['.3f'%?x] is equivalent to sprintf('%.3f', x),",
            "    - .['.2f'% ? pi] -> '3.14'.",
            "    %        No default                applies sprintf formatting",
            "",
            "SPECIALS ---------------|",
            "    Use '/' first to split the character with commas:",
            '    Ex: dsb("/x1, x2")             -> c("x1", "x2")',
            '        dsb("Hi .[/John, Linda]!") ->  c("Hi John!", "Hi Linda!")',
            "",
            "    In quoted arguments, use backticks to evaluate them from the frame.",
            '    Ex: n = 3 ; dsb("`n`*c!$") -> "$$$". The \'$\' is replicated n times, then collapsed.'
            )

        message(paste(msg, collapse = "\n"))
        return(invisible(NULL))
    }

    res = .dsb(..., frame = frame, sep = sep, concat = concat, nest = nest, check = TRUE)

    class(res) = c("dsb", "character")
    res
}

.dsb0 = function(..., frame = parent.frame(), sep = "", concat = FALSE){
    .dsb(..., frame = frame, nest = FALSE, sep = sep, concat = concat, check = FALSE)
}

.dsb = function(..., frame = parent.frame(), sep = "", concat = FALSE,
                nest = TRUE, check = FALSE){

    if(...length() == 0){
        return("")
    } else if(...length() == 1){
        x = as.character(..1)
    }else {
        # Note: using paste(..1, ..2, sep = sep) explicitly only saves 2us vav do.call
        # not worth it.

        dots = list(...)
        dots$sep = sep
        x = do.call(paste, dots)
    }

    if(length(x) == 1 && is.na(x)){
        return(NA_character_)
    }

    if(nest){
        x_parsed = list(cpp_dsb_full_string(x))

    } else if(!grepl(".[", x, fixed = TRUE)){
        return(x)

    } else {
        x_parsed = cpp_dsb(x)

    }

    n_x_all = lengths(x_parsed)

    n = length(x_parsed)
    res = character(0)
    i_done = FALSE
    for(i in 1:n){

        if(i_done){
            i_done = FALSE
            next
        }

        xi = x_parsed[[i]]

        if(length(xi) == 1){

            if(i == 1){
                res = xi
            } else {
                res = paste0(res, xi)
            }

        } else {
            operators = xi[[1]]
            xi = xi[[2]]

            if(length(operators) == 0){

                if(nest){
                    # means verbatim => no evaluation
                    if(grepl(".[", xi, fixed = TRUE)){
                        xi = .dsb(xi, frame = frame, nest = FALSE, check = check)
                    }

                } else {
                    # we need to evaluate xi
                    if(check){
                        xi_call = error_sender(str2lang(xi), "The value '", xi,
                                               "' could not be parsed.")
                    } else {
                        xi_call = str2lang(xi)
                    }

                    if(is.character(xi_call)){
                        # if a string literal => it's nesting
                        if(grepl(".[", xi_call, fixed = TRUE)){
                            xi = .dsb(xi_call, frame = frame, nest = FALSE, check = check)
                        }
                    } else {
                        if(check){
                            xi = error_sender(eval(xi_call, frame), "The value '", xi,
                                              "' could not be evaluated.")
                        } else {
                            xi = eval(xi_call, frame)
                        }
                    }
                }

            } else {

                n_op = length(operators)
                verbatim = operators[n_op] %in% c("!", "/")

                if(operators[n_op] == "/"){
                    operators = "', *'S"
                } else {
                    operators = operators[-n_op]
                    # The two separators ? and ! have no default operation
                }

                # If split operator, we concatenate
                concat_nested = length(operators) > 0 && grepl("(s|S)$", operators[[1]])

                if(verbatim && grepl(".[", xi, fixed = TRUE)){
                    xi = .dsb(xi, frame = frame, nest = FALSE, concat = concat_nested, check = check)

                } else if(!verbatim){
                    # evaluation
                    if(check){
                        xi_call = error_sender(str2lang(xi), "The value '", xi,
                                               "' could not be parsed.")
                        xi = error_sender(eval(xi_call, frame), "The value '", xi,
                                          "' could not be evaluated.")
                    } else {
                        xi = eval(str2lang(xi), frame)
                    }

                    if(is.function(xi)){
                        stop_up("dsb cannot coerce functions into strings. Problem: '",
                                trimws(x_parsed[[i]][[2]]), "' is a function.")
                    }

                }

                # Now we apply the operators
                for(j in seq_along(operators)){
                    opi = operators[[j]]
                    op_parsed = dsb_char2operator(opi, verbatim)

                    if(op_parsed$do_eval){
                        if(check){
                            quoted_call = error_sender(str2lang(op_parsed$quoted),
                                                       "In operation '", opi, "', the value '",
                                                       op_parsed$quoted, "' could not be parsed.")

                            quoted = error_sender(eval(quoted_call, frame),
                                                  "In operation '", opi, "', the value '",
                                                  op_parsed$quoted, "' could not be evaluated.")
                        } else {
                            quoted = eval(str2lang(op_parsed$quoted), frame)
                        }

                    } else {
                        quoted = op_parsed$quoted
                    }

                    if(check){
                        xi = error_sender(dsb_operators(xi, quoted, op_parsed$op),
                                          "The operation '", opi, "' failed. Please revise your call.")
                    } else {
                        xi = dsb_operators(xi, quoted, op_parsed$op)
                    }
                }
            }

            if(i == 1){
                res = xi
            } else{
                if(i < n && n_x_all[i + 1] == 1){
                    if(concat){
                        res = c(res, xi, x_parsed[[i + 1]])
                    } else {
                        res = paste0(res, xi, x_parsed[[i + 1]])
                    }

                    i_done = TRUE
                } else {
                    if(concat){
                        res = c(res, xi)
                    } else {
                        res = paste0(res, xi)
                    }
                }
            }
        }
    }

    return(res)
}


dsb_char2operator = function(x, verbatim){

    quote = substr(x, 1, 1)

    OPERATORS = c("s", "S", "x", "X", "c", "C", "r", "R", "*", "*c",
                  "u", "U", "l", "L", "q", "Q", "f", "F", "%",
                  "a", "A", "k", "K", "d", "D", "if")

    if(quote %in% c("'", "\"", "`")){
        in_quote = sub(paste0(quote, "[^", quote, "]*", "$"), "", x)
        in_quote = str_trim(in_quote, 1)
        op_abbrev = str_trim(x, nchar(in_quote) + 2)

        if(nchar(op_abbrev) == 0){
            op_abbrev = if(verbatim) "S" else "c"
        }

    } else if(x %in% OPERATORS){
        # default values
        in_quote = switch(x,
                          s = " ", S = ", *",
                          x = "[[:alnum:]]+", X = "[[:alnum:]]+",
                          c = "", C = ", ",
                          "*" = "1", "*c" = "1",
                          r = "\n", R = "^[ \t\r\n]+|[ \t\r\n]+$",
                          "")
        op_abbrev = x

        if(op_abbrev %in% c("%", "k", "K", "a", "A", "if")){
            ex = c("%" = ".['.3f'%?pi]", "k" = ".[4k ! longuest word ]", "K" = ".[2K ? 1:5]",
                   "a" = ".['|s'a ! cat]", "A" = ".['|two'A ! one]",
                   "if" = ".['len<3'if(d) ? c('a', 'abc')]")
            stop_up("The operator '", op_abbrev,
                    "' has no default value, you must provide values explicitly. Like in ",
                    ex[op_abbrev], " for instance.")
        }

    } else {
        last1 = substr(x, nchar(x), nchar(x))
        if(last1 %in% c("*", "k", "K")){
            in_quote = str_trim(x, -1)
            op_abbrev = last1
        } else {
            last2 = substr(x, nchar(x) - 1, nchar(x))
            if(last2 == "*c"){
                in_quote = str_trim(x, -2)
                op_abbrev = last2
            } else {
                op_abbrev = "problem"
            }
        }
    }

    if(!op_abbrev %in% OPERATORS){

        msg = c("The operation '", x, "' is not valid. It must be something quoted followed by a valid operator.",
                "\n  Valid operators: to split: s, S / to replace: r, R  / to collapse: c, C / to extract: x, X",
                "\n                   to replicate: * / to replicate and collapse with the empty string: *c",
                "\n                   to upper/lower case: u, U, L / to single/double quote: q, Q",
                "\n                   to format f, F / to apply sprintf format: %",
                "\n  NOTA: upper cases enable regular expressions.",
                "\n        type dsb('--help') for more help.",
                "\n  Example: .[', *'S, 'a => b'r? var] first splits the variable var by commas then replaces every 'a' with a 'b'.")

        message(msg)

        stop_up("In dsb, the operation is not valid, see upper message.")
    }

    res = list(quoted = in_quote, do_eval = quote == "`", op = op_abbrev)
    res
}


dsb_operators = function(x, quoted, op){

    if(op %in% c("s", "S")){
        # Split is always applied on verbatim stuff => length 1
        if(op == "s"){
            res = unlist(strsplit(x, quoted, fixed = TRUE))
        } else {
            res = unlist(strsplit(x, quoted, perl = TRUE))
        }

        res = res[nchar(res) > 0]

    } else if(op %in% c("c", "C")){
        # collapse

        n_x = length(x)
        if(n_x > 1 && grepl("||", quoted, fixed = TRUE)){
            # This is the "last" operator
            quoted_split = strsplit(quoted, "||", fixed = TRUE)[[1]]
            if(n_x == 2){
                res = paste(x, collapse = quoted_split[[2]])
            } else {
                res = paste(x[-n_x], collapse = quoted_split[[1]])
                res = paste0(res, quoted_split[[2]], x[n_x])
            }
        } else {
            res = paste(x, collapse = quoted)
        }

    } else if(op %in% c("r", "R")){
        new = ""
        if(grepl("=>", quoted, fixed = TRUE)){
            pat = "=>"

            if(grepl(" => ", quoted, fixed = TRUE)){
                pat = " => "
            } else if(grepl("_=>_", quoted, fixed = TRUE)){
                pat = "_=>_"
            }

            quoted_split = strsplit(quoted, pat, fixed = TRUE)[[1]]
            quoted = quoted_split[[1]]
            new = quoted_split[[2]]
        }

        if(op == "r"){
            res = gsub(quoted, new, x, fixed = TRUE)
        } else {
            res = gsub(quoted, new, x, perl = TRUE)
        }

    } else if(op %in% c("*", "*c", "*C")){
        if(!is_numeric_in_char(quoted)){
            stop("In dsb: the operator '", op, "' must have numeric arguments, '", quoted, "' is not numeric.")
        }

        res = rep(x, as.numeric(quoted))

        if(op == "*c"){
            res = paste(res, collapse = "")
        } else if(op == "*C"){
            res = paste(res, collapse = ", ")
        }

    } else if(op == "x"){
        # extract the first pattern

        x_pat = regexpr(quoted, x, perl = TRUE)

        res = substr(x, x_pat, attr(x_pat, "match.length"))

    } else if(op == "X"){
        # extract all patterns

        x_pat = gregexpr(quoted, x, perl = TRUE)

        my_extract = function(str, pat){
            if(length(pat) == 1 && pat == -1) return("")
            len = attr(pat, "match.length")
            n_pat = length(pat)
            res = character(n_pat)
            for(i in 1:n_pat){
                res[i] = substr(str, pat[i], pat[i] + len[i] - 1)
            }
            res
        }

        res = unlist(lapply(seq_along(x_pat), function(i) my_extract(x[i], x_pat[[i]])))

    } else if(op == "U"){
        res = toupper(x)
    } else if(op == "u"){
        # First letter only, if relevant
        res = x
        substr(res, 1, 1) = toupper(substr(x, 1, 1))
    } else if(op %in% c("l", "L")){
        res = tolower(x)
    } else if(op == "q"){
        res = paste0("'", x, "'")
    } else if(op == "Q"){
        res = paste0("\"", x, "\"")
    } else if(op == "f"){
        res = format(x)
    } else if(op == "F"){
        res = format(x, justify = "right")
    } else if(op == "%"){
        res = sprintf(paste0("%", quoted), x)
    } else if(op == "k"){

    } else if(op == "K"){

    } else if(op == "a"){

    } else if(op == "A"){

    } else if(op == "d"){

    } else if(op == "D"){

    } else if(op == "if"){
        stop("in progress")
    } else {
        stop("In dsb: the operator '", op, "' is not recognized. Internal error: this problem should have been spotted beforehand.")
    }

    return(res)
}


