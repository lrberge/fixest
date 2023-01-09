#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Wed Apr 13 10:05:46 2022
# ~: Data related tools
#----------------------------------------------#


#' Bins the values of a variable (typically a factor)
#'
#' Tool to easily group the values of a given variable.
#'
#' @param x A vector whose values have to be grouped. Can be of any type but must be atomic.
#' @param bin A list of values to be grouped, a vector, a formula, or the special values `"bin::digit"` or `"cut::values"`. To create a new value from old values, use `bin = list("new_value"=old_values)` with `old_values` a vector of existing values. You can use `.()` for `list()`.
#' It accepts regular expressions, but they must start with an `"@"`, like in `bin="@Aug|Dec"`. It accepts one-sided formulas which must contain the variable `x`, e.g. `bin=list("<2" = ~x < 2)`.
#' The names of the list are the new names. If the new name is missing, the first value matched becomes the new name. In the name, adding `"@d"`, with `d` a digit, will relocate the value in position `d`: useful to change the position of factors. Use `"@"` as first item to make subsequent items be located first in the factor.
#' Feeding in a vector is like using a list without name and only a single element. If the vector is numeric, you can use the special value `"bin::digit"` to group every `digit` element.
#' For example if `x` represents years, using `bin="bin::2"` creates bins of two years.
#' With any data, using `"!bin::digit"` groups every digit consecutive values starting from the first value.
#' Using `"!!bin::digit"` is the same but starting from the last value.
#' With numeric vectors you can: a) use `"cut::n"` to cut the vector into `n` equal parts, b) use `"cut::a]b["` to create the following bins: `[min, a]`, `]a, b[`, `[b, max]`.
#' The latter syntax is a sequence of number/quartile (q0 to q4)/percentile (p0 to p100) followed by an open or closed square bracket. You can add custom bin names by adding them in the character vector after `'cut::values'`. See details and examples. Dot square bracket expansion (see [`dsb`]) is enabled.
#'
#' @section "Cutting" a numeric vector:
#'
#' Numeric vectors can be cut easily into: a) equal parts, b) user-specified bins.
#'
#' Use `"cut::n"` to cut the vector into `n` (roughly) equal parts. Percentiles are used to partition the data, hence some data distributions can lead to create less than `n` parts (for example if P0 is the same as P50).
#'
#' The user can specify custom bins with the following syntax: `"cut::a]b]c]"`. Here the numbers `a`, `b`, `c`, etc, are a sequence of increasing numbers, each followed by an open or closed square bracket. The numbers can be specified as either plain numbers (e.g. \code{"cut::5]12[32["}), quartiles (e.g. \code{"cut::q1]q3["}), or percentiles (e.g. `"cut::p10]p15]p90]"`). Values of different types can be mixed: \code{"cut::5]q2[p80["} is valid provided the median (`q2`) is indeed greater than `5`, otherwise an error is thrown.
#'
#' The square bracket right of each number tells whether the numbers should be included or excluded from the current bin. For example, say `x` ranges from 0 to 100, then `"cut::5]"` will create two  bins: one from 0 to 5 and a second from 6 to 100. With `"cut::5["` the bins would have been 0-4 and 5-100.
#'
#' A factor is always returned. The labels always report the min and max values in each bin.
#'
#' To have user-specified bin labels, just add them in the character vector following `'cut::values'`. You don't need to provide all of them, and `NA` values fall back to the default label. For example, `bin = c("cut::4", "Q1", NA, "Q3")` will modify only the first and third label that will be displayed as `"Q1"` and `"Q3"`.
#'
#' @section `bin` vs `ref`:
#'
#' The functions [`bin`] and [`ref`] are able to do the same thing, then why use one instead of the other? Here are the differences:
#'
#' \itemize{
#' \item{}{`ref` always returns a factor. This is in contrast with `bin` which returns, when possible, a vector of the same type as the vector in input.}
#' \item{}{`ref` always places the values modified in the first place of the factor levels. On the other hand, `bin` tries to not modify the ordering of the levels. It is possible to make `bin` mimic the behavior of `ref` by adding an `"@"` as the first element of the list in the argument `bin`.}
#'  \item{}{when a vector (and not a list) is given in input, `ref` will place each element of the vector in the first place of the factor levels. The behavior of `bin` is totally different, `bin` will transform all the values in the vector into a single value in `x` (i.e. it's binning).}
#' }
#'
#' @return
#' It returns a vector of the same length as `x`.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' To re-factor variables: [`ref`].
#'
#' @examples
#'
#' data(airquality)
#' month_num = airquality$Month
#' table(month_num)
#'
#' # Grouping the first two values
#' table(bin(month_num, 5:6))
#'
#' # ... plus changing the name to '10'
#' table(bin(month_num, list("10" = 5:6)))
#'
#' # ... and grouping 7 to 9
#' table(bin(month_num, list("g1" = 5:6, "g2" = 7:9)))
#'
#' # Grouping every two months
#' table(bin(month_num, "bin::2"))
#'
#' # ... every 2 consecutive elements
#' table(bin(month_num, "!bin::2"))
#'
#' # ... idem starting from the last one
#' table(bin(month_num, "!!bin::2"))
#'
#' # Using .() for list():
#' table(bin(month_num, .("g1" = 5:6)))
#'
#'
#' #
#' # with non numeric data
#' #
#'
#' month_lab = c("may", "june", "july", "august", "september")
#' month_fact = factor(month_num, labels = month_lab)
#'
#' # Grouping the first two elements
#' table(bin(month_fact, c("may", "jun")))
#'
#' # ... using regex
#' table(bin(month_fact, "@may|jun"))
#'
#' # ...changing the name
#' table(bin(month_fact, list("spring" = "@may|jun")))
#'
#' # Grouping every 2 consecutive months
#' table(bin(month_fact, "!bin::2"))
#'
#' # ...idem but starting from the last
#' table(bin(month_fact, "!!bin::2"))
#'
#' # Relocating the months using "@d" in the name
#' table(bin(month_fact, .("@5" = "may", "@1 summer" = "@aug|jul")))
#'
#' # Putting "@" as first item means subsequent items will be placed first
#' table(bin(month_fact, .("@", "aug", "july")))
#'
#' #
#' # "Cutting" numeric data
#' #
#'
#' data(iris)
#' plen = iris$Petal.Length
#'
#' # 3 parts of (roughly) equal size
#' table(bin(plen, "cut::3"))
#'
#' # Three custom bins
#' table(bin(plen, "cut::2]5]"))
#'
#' # .. same, excluding 5 in the 2nd bin
#' table(bin(plen, "cut::2]5["))
#'
#' # Using quartiles
#' table(bin(plen, "cut::q1]q2]q3]"))
#'
#' # Using percentiles
#' table(bin(plen, "cut::p20]p50]p70]p90]"))
#'
#' # Mixing all
#' table(bin(plen, "cut::2[q2]p90]"))
#'
#' # NOTA:
#' # -> the labels always contain the min/max values in each bin
#'
#' # Custom labels can be provided, just give them in the char. vector
#' # NA values lead to the default label
#' table(bin(plen, c("cut::2[q2]p90]", "<2", "]2; Q2]", NA, ">90%")))
#'
#'
#'
#' #
#' # With a formula
#' #
#'
#' data(iris)
#' plen = iris$Petal.Length
#'
#' # We need to use "x"
#' table(bin(plen, list("< 2" = ~x < 2, ">= 2" = ~x >= 2)))
#'
#'
bin = function(x, bin){
    check_arg(x, "vector mbt")

    bin = error_sender(eval_dot(bin), arg_name = "bin")

    check_arg(bin, "list | vector mbt")

    varname = deparse(substitute(x))[1]
    bin_factor(bin, x, varname)
}

#' Refactors a variable
#'
#' Takes a variables of any types, transforms it into a factors, and modifies the values of the factors. Useful in estimations when you want to set some value of a vector as a reference.
#'
#' @inheritSection bin "Cutting" a numeric vector
#' @inheritSection bin `bin` vs `ref`
#'
#'
#' @param x A vector of any type (must be atomic though).
#' @param ref A vector or a list, or special binning values (explained later). If a vector, it must correspond to (partially matched) values of the vector `x`. The vector `x` which will be transformed into a factor and these values will be placed first in the levels. That's the main usage of this function. You can also bin on-the-fly the values of `x`, using the same syntax as the function [`bin`]. To create a new value from old values, use `ref = list("new_value"=old_values)` with `old_values` a vector of existing values. You can use `.()` for `list()`.
#' It accepts regular expressions, but they must start with an `"@"`, like in `ref="@Aug|Dec"`. It accepts one-sided formulas which must contain the variable `x`, e.g. `ref=list("<2" = ~x < 2)`.
#' The names of the list are the new names. If the new name is missing, the first value matched becomes the new name. In the name, adding `"@d"`, with `d` a digit, will relocate the value in position `d`: useful to change the position of factors.
#' If the vector `x` is numeric, you can use the special value `"bin::digit"` to group every `digit` element.
#' For example if `x` represents years, using `ref="bin::2"` creates bins of two years.
#' With any data, using `"!bin::digit"` groups every digit consecutive values starting from the first value.
#' Using `"!!bin::digit"` is the same but starting from the last value.
#' With numeric vectors you can: a) use `"cut::n"` to cut the vector into `n` equal parts, b) use `"cut::a]b["` to create the following bins: `[min, a]`, `]a, b[`, `[b, max]`.
#' The latter syntax is a sequence of number/quartile (q0 to q4)/percentile (p0 to p100) followed by an open or closed square bracket. You can add custom bin names by adding them in the character vector after `'cut::values'`. See details and examples. Dot square bracket expansion (see [`dsb`]) is enabled.
#'
#' @return
#' It returns a factor of the same length as `x`, where levels have been modified according to the argument `ref`.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' To bin the values of a vector: [`bin`].
#'
#' @examples
#'
#' data(airquality)
#'
#' # A vector of months
#' month_num = airquality$Month
#' month_lab = c("may", "june", "july", "august", "september")
#' month_fact = factor(month_num, labels = month_lab)
#' table(month_num)
#' table(month_fact)
#'
#' #
#' # Main use
#' #
#'
#' # Without argument: equivalent to as.factor
#' ref(month_num)
#'
#' # Main usage: to set a level first:
#' # (Note that partial matching is enabled.)
#' table(ref(month_fact, "aug"))
#'
#' # You can rename the level on-the-fly
#' # (Northern hemisphere specific!)
#' table(ref(month_fact, .("Hot month"="aug",
#'                         "Late summer" = "sept")))
#'
#'
#' # Main use is in estimations:
#' a = feols(Petal.Width ~ Petal.Length + Species, iris)
#'
#' # We change the reference
#' b = feols(Petal.Width ~ Petal.Length + ref(Species, "vers"), iris)
#'
#' etable(a, b)
#'
#'
#' #
#' # Binning
#' #
#'
#' # You can also bin factor values on the fly
#' # Using @ first means a regular expression will be used to match the values.
#' # Note that the value created is placed first.
#' # To avoid that behavior => use the function "bin"
#' table(ref(month_fact, .(summer = "@jul|aug|sep")))
#'
#' # Please refer to the example in the bin help page for more example.
#' # The syntax is the same.
#'
#'
#' #
#' # Precise relocation
#' #
#'
#' # You can place a factor at the location you want
#' #  by adding "@digit" in the name first:
#' table(ref(month_num, .("@5"=5)))
#'
#' # Same with renaming
#' table(ref(month_num, .("@5 five"=5)))
#'
#'
ref = function(x, ref){
    check_arg(x, "vector mbt")

    if(missing(ref)){
        return(as.factor(x))
    }

    ref = error_sender(eval_dot(ref), arg_name = "ref")
    check_arg(ref, "list | vector mbt")

    varname = deparse(substitute(x))[1]

    IS_SPECIAL = FALSE
    if(!is.list(ref)){
        if(is.character(ref[1]) && grepl("^(cut|bin)", ref[1])){
            IS_SPECIAL = TRUE
            if(!is.numeric(x)){
                stop(.dsb("To use the special binning `.[ref[1]]` the variable ",
                          "`.[varname]` must be numeric. Currently this is not the case ",
                          "(it is of class .[3KO, C?class(x)] instead)."))
            }
        } else {
            ref = as.list(ref)
        }
    }

    if(!IS_SPECIAL && !is.factor(x)){
        x = as.factor(x)
    }

    if(!IS_SPECIAL && ref[[1]] != "@"){
        ref_new = list("@")
        index = 1:length(ref) + 1
        ref_new[index] = ref
        names(ref_new)[index] = names(ref)
        ref = ref_new
    }

    bin_factor(ref, x, varname)
}



#' Prints the number of unique elements in a data set
#'
#' This utility tool displays the number of unique elements in one or multiple data.frames as well as their number of NA values.
#'
#' @param x A formula, with data set names on the LHS and variables on the RHS, like `data1 + data2 ~ var1 + var2`. The following special variables are admitted: `"."` to get default values, `".N"` for the number of observations, `".U"` for the number of unique rows, `".NA"` for the number of rows with at least one NA. Variables can be combined with `"^"`, e.g. `df~id^period`; use `id%^%period` to also include the terms on both sides. Note that using `:` and `*` is equivalent to `^` and `%^%`. Sub select with `id[cond]`, when doing so `id` is automatically included. Conditions can be chained, as in `id[cond1, cond2]`. Use `NA(x, y)` in conditions instead of `is.na(x) | is.na(y)`. Use the `!!` operator to have both a condition and its opposite. To compare the keys in two data sets, use `data1:data2`. If not a formula, `x` can be: a vector (displays the # of unique values); a `data.frame` (default values are displayed), or a "sum" of data sets like in `x = data1 + data2`, in that case it is equivalent to `data1 + data2 ~ .`.
#' @param ... Not currently used.
#'
#' @section Special values and functions:
#'
#' In the formula, you can use the following special values: `"."`, `".N"`, `".U"`, and `".NA"`.
#'
#' \itemize{
#'
#' \item{`"."`}{Accesses the default values. If there is only one data set and the data set is *not* a `data.table`, then the default is to display the number of observations and the number of unique rows. If the data is a `data.table`, the number of unique items in the key(s) is displayed instead of the number of unique rows (if the table has keys of course). If there are two or more data sets, then the default is to display the unique items for: a) the variables common across all data sets, if there's less than 4, and b) if no variable is shown in a), the number of variables common across at least two data sets, provided there are less than 5. If the data sets are data tables, the keys are also displayed on top of the common variables. In any case, the number of observations is always displayed.}
#'
#' \item{`".N"`}{Displays the number of observations.}
#'
#' \item{`".U"`}{Displays the number of unique rows.}
#'
#' \item{`".NA"`}{Displays the number of rows with at least one NA.}
#'
#' }
#'
#' @section The `NA` function:
#'
#' The special function `NA` is an equivalent to `is.na` but can handle several variables. For instance, `NA(x, y)` is equivalent to `is.na(x) | is.na(y)`. You can add as many variables as you want as arguments. If no argument is provided, as in `NA()`, it is identical to having all the variables of the data set as argument.
#'
#' @section Combining variables:
#'
#' Use the "hat", `"^"`, operator to combine several variables. For example `id^period` will display the number of unique values of id x period combinations.
#'
#' Use the "super hat", `"%^%"`, operator to also include the terms on both sides. For example, instead of writing `id + period + id^period`, you can simply write `id%^%period`.
#'
#' Alternatively, you can use `:` for `^` and `*` for `%^%`.
#'
#' @section Sub-selections:
#'
#' To show the number of unique values for sub samples, simply use `[]`. For example, `id[x > 10]` will display the number of unique `id` for which `x > 10`.
#'
#' Simple square brackets lead to the inclusion of both the variable and its subset. For example `id[x > 10]` is equivalent to `id + id[x > 10]`. To include only the sub selection, use double square brackets, as in `id[[x > 10]]`.
#'
#' You can add multiple sub selections at once, only separate them with a comma. For example `id[x > 10, NA(y)]` is equivalent to `id[x > 10] + id[NA(y)]`.
#'
#' Use the double negative operator, i.e. `!!`, to include both a condition and its opposite at once. For example `id[!!x > 10]` is equivalent to `id[x > 10, !x > 10]`. Double negative operators can be chained, like in `id[!!cond1 & !!cond2]`, then the cardinal product of all double negatived conditions is returned.
#'
#' @return
#' It returns a vector containing the number of unique values per element. If several data sets were provided, a list is returned, as long as the number of data sets, each element being a vector of unique values.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' data = base_did
#' data$x1.L1 = round(lag(x1~id+period, 1, data))
#'
#' # By default, just the formatted number of observations
#' n_unik(data)
#'
#' # Or the nber of unique elements of a vector
#' n_unik(data$id)
#'
#' # number of unique id values and id x period pairs
#' n_unik(data ~.N + id + id^period)
#'
#' # use the %^% operator to include the terms on the two sides at once
#' # => same as id*period
#' n_unik(data ~.N + id %^% period)
#'
#' # using sub selection with []
#' n_unik(data ~.N + period[!NA(x1.L1)])
#'
#' # to show only the sub selection: [[]]
#' n_unik(data ~.N + period[[!NA(x1.L1)]])
#'
#' # you can have multiple values in [],
#' # just separate them with a comma
#' n_unik(data ~.N + period[!NA(x1.L1), x1 > 7])
#'
#' # to have both a condition and its opposite,
#' # use the !! operator
#' n_unik(data ~.N[!!NA(x1.L1)])
#'
#' # the !! operator works within condition chains
#' n_unik(data ~.N[!!NA(x1.L1) & !!x1 > 7])
#'
#' # Conditions can be distributed
#' n_unik(data ~ (id + period)[x1 > 7])
#'
#' #
#' # Several data sets
#' #
#'
#' # Typical use case: merging
#' # Let's create two data sets and merge them
#'
#' data(base_did)
#' base_main = base_did
#' base_extra = sample_df(base_main[, c("id", "period")], 100)
#' base_extra$id[1:10] = 111:120
#' base_extra$period[11:20] = 11:20
#' base_extra$z = rnorm(100)
#'
#' # You can use db1:db2 to compare the common keys in two data sets
#'  n_unik(base_main:base_extra)
#'
#' tmp = merge(base_main, base_extra, all.x = TRUE, by = c("id", "period"))
#'
#' # You can show unique values for any variable, as before
#' n_unik(tmp + base_main + base_extra ~ id[!!NA(z)] + id^period)
#'
#'
#'
n_unik = function(x){
    # returns a vector with the nber of unique values
    # attr("na.info") => nber of NA values, vector

    if(missing(x)){
        stop("Argument `x` must be provided. Problem: it is missing.")
    }

    # Non standard-evaluation
    x_dp = deparse_long(substitute(x))
    if(!grepl("~", x_dp, fixed = TRUE)){
        if(grepl("[+:*]", x_dp)){
            # Non standard-evaluation
            x = as.formula(paste0(x_dp, "~ ."))
        } else {
            check_arg(x, "data.frame | vector l0 | ts formula")
        }
    } else {
        check_arg(x, "ts formula")
    }

    nthreads = getFixest_nthreads()

    # If vector
    comp_pairs = list()
    if(is.vector(x)){

        x_name = gsub("^[[:alpha:]\\.][[:alnum:]\\._]*\\$", "", x_dp)

        na.info = 0
        if(anyNA(x)){
            who_NA = is.na(x)
            na.info = sum(who_NA)
            x = x[!who_NA]
        }

        res = len_unique(x, nthreads)
        names(res) = x_name

        attr(res, "na.info") = na.info

        class(res) = "vec_n_unik"
        return(res)

    } else if(is.data.frame(x)){

        n_x = 1
        x_all = list(x)
        fml = ~.

    } else if(inherits(x, "formula")){

        x_terms = terms(x[1:2])

        fact_mat = attr(x_terms, "factors")

        x_all_names = rownames(fact_mat)

        # We get the comparison pairs
        vars = colnames(fact_mat)
        for(pair in grep(":", vars, fixed = TRUE, value = TRUE)){
            dict = setNames(1:length(x_all_names), x_all_names)
            pair_split = strsplit(pair, ":", fixed = TRUE)[[1]]
            comb = combn(pair_split, 2)
            for(i in 1:ncol(comb)){
                comp_pairs[[length(comp_pairs) + 1]] = unname(dict[comb[, i]])
            }
        }

        if(length(comp_pairs) > 0){
            comp_pairs = unique(comp_pairs)
        }

        # Construction of the list  + sanity check
        n_x = length(x_all_names)
        x_all = vector("list", n_x)
        for(i in 1:n_x){
            x_all[[i]] = error_sender(eval(str2lang(x_all_names[i]), parent.frame()),
                                      "The left-hand-side of the formula must contain valid data.frames. Problem in the evaluation of `", x_all_names[i], "`:")
            check_value(x_all[[i]], "data.frame",
                        .message = paste0("The value `", x_all_names[i], "` (in the LHS of the formula) must be a data.frame."))
        }

        fml = x[c(1, 3)]

    }

    # If DF + formula
    # ex: fml = ~ id^year + author_id[sw0(is.na(author_name), is.na(affil_name), year == min_year)]
    # fml = ~ id^sw0(fe)

    # check variable names
    naked_vars = origin_vars = all.vars(fml)
    valid_vars = c(".N", ".U", ".NA", ".", unique(unlist(lapply(x_all, names))))

    check_value_plus(naked_vars, "multi match", .choices = valid_vars,
                     .message = paste0("The formula must only use variables in the data set", plural_len(x_all), "."))

    if("." %in% naked_vars){
        # The default values

        IS_DT = requireNamespace("data.table", quietly = TRUE)

        dot_default = ".N"
        x_keys = c()
        x_common_vars = c()

        # keys
        if(IS_DT){
            for(I in 1:n_x){
                if(inherits(x_all[[I]], "data.table")){
                    x_keys = c(x_keys, data.table::key(x_all[[I]]))
                }
            }
        }

        # common variables
        if(n_x > 1){

            # Common to all
            x_names_current = names(x_all[[1]])
            for(i in 2:n_x){
                x_names_current = intersect(x_names_current, names(x_all[[i]]))
                if(length(x_names_current) == 0) break
            }

            if(length(x_names_current) %in% 1:4){
                # No more than 4 common variables by default
                x_common_vars = x_names_current

                if(length(comp_pairs) > 0 && length(x_names_current) <= 3){
                    x_common_vars = paste0(x_names_current, collapse = "*")
                }

            } else if(n_x > 2){
                # Common to at least 2 data sets
                for(i in 1:(n_x - 1)){
                    x_names_current = names(x_all[[i]])
                    for(j in (i + 1):n_x){
                        qui_common = x_names_current %in% names(x_all[[j]])
                        if(any(qui_common)){
                            x_common_vars = c(x_common_vars, x_names_current[qui_common])
                        }
                    }
                }

                x_common_vars = unique(x_common_vars)

                if(length(x_common_vars) > 5){
                    #  we keep max 5
                    x_common_vars = c()
                }
            }
        }

        if(length(x_keys) > 0 || length(x_common_vars) > 0){
            dot_default = unique(c(dot_default, x_keys, x_common_vars))
        }

        if(length(dot_default) == 1){
            dot_default = c(".N", ".U")
        }

        rhs_txt = as.character(fml)[[2]]
        rhs_txt = gsub("(^| )\\.( |$)", dsb(".['+'c? dot_default]"), rhs_txt)
        fml = .xpd(rhs = rhs_txt)
    }

    if(any(naked_vars != origin_vars)){
        # we complete the partial matching
        fml_txt = deparse_long(fml)

        var_diff = which(naked_vars != origin_vars)

        for(i in var_diff){
            re = paste0("(?!<[[:alnum:]._])", origin_vars[i], "(?!=[[:alnum:]._])")
            fml_txt = gsub(re, naked_vars[i], fml_txt, perl = TRUE)
        }

        fml = as.formula(fml_txt)
    }

    tm = terms_hat(fml, n_unik = TRUE)
    all_vars = attr(tm, "term.labels")
    var_final = c()

    # stepwise function
    for(var in all_vars){
        # var = all_vars[1]

        if(grepl("combine_clusters(_fast)?", var)){
            var = gsub("combine_clusters(_fast)?", "to_integer", var)
        }

        IS_HAT_SW = grepl("sw0?\\)\\(", var)
        if(IS_HAT_SW){
            var = paste0(gsub("(sw0?)\\)\\(", "\\1(", var), ")")
        }

        if(grepl("^[[:alpha:].][[:alnum:]._]*\\[", var, perl = TRUE) && !grepl("sw0?\\(", var)){
            # We use sw to make id[cond1, cond2] work
            # That's what's called path dependency!

            var_new = gsub("^([[:alpha:].][[:alnum:]._]*)\\[", "\\1[sw0(", var)
            if(grepl("sw0([", var_new, fixed = TRUE)){
                var_new = sub("sw0([", "sw(", str_trim(var_new, -1), fixed = TRUE)
            }
            var_new = sub("\\]$", ")]", var_new)

            var = var_new
        }

        if(is_fun_in_char(var, "sw0?")){

            info = extract_fun(var, "sw0?", err_msg = "The stepwise function can be used only once per variable.")

            sw_value = eval(str2lang(info$fun))

            qui_double = grepl("!!", trimws(sw_value))
            while(any(qui_double)){
                i = which(qui_double)[1]

                value = sw_value[i]
                value_split = strsplit(value, "!!")[[1]]

                n_split = length(value_split)
                n_s = 2 ** (n_split - 1)
                new_values = rep(value_split[1], n_s)

                prefixes = do.call(expand.grid, rep(list(c("", "!")), n_split - 1))

                for(j in 1:ncol(prefixes)){
                    new_values = paste0(new_values, prefixes[, j], value_split[j + 1])
                }

                sw_value = insert(sw_value[-i], new_values, i)
                qui_double = grepl("!!", sw_value)
            }

            var_char_new = paste0(info$before, sw_value, info$after)

            if(IS_HAT_SW){
                var_char_new[1] = gsub(", (,|\\))", "\\1", var_char_new[1])

                qui_nested = grepl("^(to_integer\\(.+){2,}", var_char_new)
                if(any(qui_nested)){
                    # later => add check problem if another function is used
                    var_char_new[qui_nested] = paste0("to_integer(", gsub("to_integer\\(|\\)", "", var_char_new[qui_nested]), ")")
                }
            }

            if(grepl("[]", var_char_new[1], fixed = TRUE)){
                var_char_new[1] = gsub("[]", "", var_char_new[1], fixed = TRUE)
            }

            var_final = c(var_final, var_char_new)
        } else {
            var_final = c(var_final, var)
        }
    }

    var_final = unique(var_final)

    var_final_names = var_final
    # we take care of id^year cases
    if(any(who_combine <- is_fun_in_char(var_final, "to_integer"))){

        for(i in which(who_combine)){
            info = extract_fun(var_final[i], "to_integer")

            # we use sw() to get the vars
            sw_char = gsub("^[^\\(]+\\(", "sw(", info$fun)
            sw_value = eval(str2lang(sw_char))

            hat_value = paste(sw_value, collapse = "^")

            var_final_names[i] = paste0(info$before, hat_value, info$after)
        }
    }

    # NA counting + unique counting

    res_all = list()
    for(I in 1:n_x){

        x = x_all[[I]]
        x_list = unclass(x)
        x_list[["NA_fun"]] = function(...) NA_fun(..., df = x)

        vars_legit = c(".N", ".U", ".NA", names(x))

        res = c()
        na.info = c()

        for(i in seq_along(var_final)){

            vf = var_final[i]
            vf_name = var_final_names[i]
            na_i = 0

            vf = gsub("((?<=[^[:alnum:]_\\.])|^)NA\\(", "NA_fun(", vf, perl = TRUE)

            vf_call = str2lang(vf)

            if(!all(all.vars(vf_call) %in% vars_legit)){
                res_i = NA_real_
                vf_name = ""

            } else if(grepl("^\\.(N|U)($|\\[)", vf)){


                if(vf %in% c(".N", ".N[]")){
                    res_i = nrow(x)
                    vf_name = "# Observations"

                } else if(vf %in% c(".U", ".U[]")){
                    res_i = nrow(unique(x))
                    vf_name = "# Unique rows"

                } else {
                    # Other methods are faster but are less general

                    vf_new = gsub("(^\\.(N|U)\\[)|(\\]$)", "", vf)
                    do_unik = grepl("^\\.U", vf)

                    val = eval(str2lang(vf_new), x_list)

                    # we want to drop the NAs for indices
                    if(anyNA(val)){
                        val = val[!is.na(val)]
                    }

                    if(length(val) == 0){
                        res_i = 0

                    } else if(do_unik){
                        res_i = NROW(unique(x[val, ]))

                    } else {
                        res_i = if(is.logical(val)) sum(val) else length(val)
                    }

                    msg = if(do_unik) "# Unique rows" else "# Obs."
                    vf_new = gsub("NA_fun", "NA", vf_new, fixed = TRUE)
                    vf_name = paste0(msg, " with ", vf_new)
                }

            } else if(grepl("^(\\.NA|NA_fun\\(\\))", vf)) {

                if(vf %in% c(".NA", ".NA[]", "NA_fun()", "NA_fun()[]")){
                    res_i = sum(!complete.cases(x))
                    vf_name = "# Rows with NAs"

                } else {

                    subselect = NULL
                    sub_text = ""
                    if(grepl("\\]$", vf)){
                        sub_text_split = strsplit(vf, "[", fixed = TRUE)[[1]]
                        sub_text = paste0(sub_text_split[-1], collapse = "[")
                        sub_text = gsub("\\]$", "", sub_text)

                        subselect = eval(str2lang(sub_text), x_list)
                    }

                    if(length(subselect) == 0){
                        res_i = 0
                    } else {
                        res_i = sum(!complete.cases(x[subselect, ]))
                    }

                    vf_name = dsb("# Rows with NAs, with .[sub_text]")
                }
            } else {
                val = eval(vf_call, x_list)

                if(anyNA(val)){
                    who_NA = is.na(val)
                    na_i = sum(who_NA)
                    val = val[!who_NA]
                }

                if(grepl("^NA_fun\\(", vf)){
                    res_i = sum(val)

                    # now the message
                    leftover = sub("^[^\\)]+\\)", "", vf)
                    begin = substr(vf, 1, nchar(vf) - nchar(leftover))

                    my_fun = list("NA_fun" = function(x) enumerate_items(sapply(sys.call()[-1], deparse_long)))
                    na_vars = eval(str2lang(begin), my_fun)
                    msg_start = paste0("# NAs in ", na_vars)

                    msg_end = ""
                    if(grepl("[", vf, fixed = TRUE)){
                        msg_end = paste0(" with ", gsub("^\\[|\\]$", "", leftover))
                    }

                    vf_name = paste0(msg_start, msg_end)

                } else {
                    res_i = len_unique(val, nthreads)
                }

            }

            res[vf_name] = res_i
            na.info[i] = na_i
        }

        attr(res, "na.info") = na.info
        class(res) = "vec_n_unik"

        if(n_x == 1){
            return(res)
        }

        res_all[[x_all_names[I]]] = res
    }


    # Specific data set comparisons
    info_pairs = list()
    for(pair in comp_pairs){

        i_x = pair[1]
        i_y = pair[2]

        x = x_all[[i_x]]
        x_list = unclass(x)
        x_list[["NA_fun"]] = function(...) NA_fun(..., df = x)

        y = x_all[[i_y]]
        y_list = unclass(y)
        y_list[["NA_fun"]] = function(...) NA_fun(..., df = y)

        vars_legit = intersect(names(x), names(y))

        # two last elements: id and common
        all_rows_id_common = c()
        row_temp = rep(NA_real_, n_x + 2)

        for(i in seq_along(var_final)){

            vf = var_final[i]
            vf = gsub("((?<=[^[:alnum:]_\\.])|^)NA\\(", "NA_fun(", vf, perl = TRUE)
            vf_call = str2lang(vf)

            if(!all(all.vars(vf_call) %in% vars_legit) || grepl("^NA_fun", vf)){
                next

            } else {

                if(grepl("to_integer", vf, fixed = TRUE)){
                    # specific scheme for combinations

                    if(grepl("NA_fun()", vf, fixed = TRUE)){
                        # we don't allow that
                        next
                    }

                    vars_keep = all.vars(vf_call)

                    xy_list = list()
                    for(v in vars_keep){
                        xy_list[[v]] = c(x_list[[v]], y_list[[v]])
                    }

                    xy_list[["NA_fun"]] = function(...) NA_fun(..., df = stop("Internal error."))


                    val_xy = eval(vf_call, xy_list)

                    nrow_x = length(x_list[[1]])
                    nrow_y = length(y_list[[1]])
                    val_x = val_xy[1:nrow_x]
                    val_y = val_xy[(nrow_x + 1):(nrow_x + nrow_y)]

                } else {
                    val_x = eval(vf_call, x_list)
                    val_y = eval(vf_call, y_list)
                }


                if(anyNA(val_x)){
                    val_x = val_x[!is.na(val_x)]
                }

                if(anyNA(val_y)){
                    val_y = val_y[!is.na(val_y)]
                }

                # first col is ID
                val_x_unik = unique(val_x)
                val_y_unik = unique(val_y)

                n_x_not_in_y = sum(!val_x_unik %in% val_y_unik)
                n_y_not_in_x = sum(!val_y_unik %in% val_x_unik)
                n_common = sum(val_x_unik %in% val_y_unik)

                row = row_temp
                row[i_x] = n_x_not_in_y
                row[i_y] = n_y_not_in_x

                row[n_x + 1] = i
                row[n_x + 2] = n_common

                all_rows_id_common = rbind(all_rows_id_common, row)

            }
        }

        if(length(all_rows_id_common) > 0){
            attr(all_rows_id_common, "data_id") = c(i_x, i_y)
            info_pairs[[length(info_pairs) + 1]] = all_rows_id_common
        }
    }

    if(length(info_pairs) > 0){
        attr(res_all, "info_pairs") = info_pairs
    }

    class(res_all) = "list_n_unik"

    return(res_all)
}

#' @rdname n_unik
print.vec_n_unik = function(x, ...){

    hash = "## "

    x_names = sfill(names(x))
    na.info = attr(x, "na.info")
    na.info_format = sfill(fsignif(na.info))
    x_value = sfill(fsignif(x))

    n = length(x)
    na_col = paste0("(# NAs: ", na.info_format, ")")
    na_col[na.info == 0] = ""

    res = paste0(hash, x_names, ": ", x_value, " ", na_col)
    cat(res, sep = "\n")
}

#' @rdname n_unik
print.list_n_unik = function(x, ...){


    # I can't use #> anymore!!! The auto complete by the new version
    # of Rstudio drives me nuts!!!
    hash = "## "

    x_all = x
    n_x = length(x_all)

    if(n_x == 1) return(x_all[[1]])

    info_pairs = attr(x, "info_pairs")
    IS_PAIR = !is.null(info_pairs)

    # First, the variable names

    all_names = names(x_all[[1]])
    qui_0 = nchar(all_names) == 0
    i = 2
    while(any(qui_0) && i <= n_x){
        names_i = names(x_all[[i]])
        all_names[qui_0] = names_i[qui_0]
        qui_0 = nchar(all_names) == 0
        i = i + 1
    }

    # id_vars: used in info_pairs
    id_vars = seq_along(all_names)

    if(all(qui_0)){
        stop("Not any valid information to display: please check that all your variables exist in all data sets.")

    } else if(any(qui_0)){
        warning("Some variables could not be evaluated in any data set: please check that all your variables exist in all data sets.")

        all_names = all_names[!qui_0]
        id_vars = id_vars[!qui_0]

        for(i in 1:n_x){
            x = x_all[[i]]
            na.info = attr(x, "na.info")

            x = x[!qui_0]
            na.info = na.info[!qui_0]
            attr(x, "na.info") = na.info

            x_all[[i]] = x
        }
    }


    x_names = sfill(all_names)

    # If ambiguous pairs: we add data set suffixes
    if(n_x > 2 && IS_PAIR){
        add_suffix = function(x, i) paste0(x, " (", letters[i], ")")
        var_width = max(nchar("Exclusive "), nchar(x_names[1]))
    } else {
        add_suffix = function(x, i) x
        var_width = nchar(x_names[1])
    }

    var_col = paste0(hash, x_names, ":")
    na_intro = paste0(hash, sfill(" ", var_width), "|NAs:")
    var_col = insert_in_between(var_col, na_intro)
    # we add the first row of the data sets names + format
    var_col = sfill(c(hash, var_col), right = TRUE)

    print_mat = var_col

    KEEP_NA_ROW = rep(FALSE, length(x_names))

    for(i in 1:n_x){

        data_name = add_suffix(names(x_all)[i], i)

        x = x_all[[i]]

        na.info = attr(x, "na.info")
        KEEP_NA_ROW = KEEP_NA_ROW | na.info > 0
        na.info_format = fsignif(na.info)
        x_value = fsignif(x)
        x_value[is.na(x)] = "--"
        na.info_format[is.na(x)] = "--"


        width = max(nchar(na.info_format), nchar(x_value))

        na.info_format = sfill(na.info_format, width)
        x_value = sfill(x_value, width)

        x_col = insert_in_between(x_value, na.info_format)
        x_col = sfill(c(data_name, x_col))

        print_mat = cbind(print_mat, x_col)
    }

    if(!any(KEEP_NA_ROW)){
        print_mat[, 1] = substr(print_mat[, 1], 1, nchar(print_mat[1, 1]) - 4)
    }

    keep = c(TRUE, insert_in_between(TRUE, KEEP_NA_ROW))

    print_mat = print_mat[keep, ]

    # PAIR information
    if(IS_PAIR){

        # identifiers used for insertion
        id_vars_all = c(0, insert_in_between(id_vars, 0))
        id_vars_all = id_vars_all[keep]

        # we add two columns: id and common
        print_mat = cbind(print_mat, id_vars_all, "")

        insert_row_after_id = function(mat, row, col_id){

            for(i in 1:nrow(row)){
                i_id_mat = max(which(mat[, col_id] == row[i, col_id]))
                tmp_before = mat[1:i_id_mat, ]

                if(i_id_mat < nrow(mat)){
                    tmp_after = mat[(i_id_mat + 1):nrow(mat), ]
                } else {
                    tmp_after = NULL
                }

                mat = rbind(tmp_before, row[i, ], tmp_after)
            }

            mat
        }

        for(pair in info_pairs){

            data_id = attr(pair, "data_id")
            pair = as.data.frame(pair)

            if(n_x == 2){
                # data_col = paste0(hash, sfill(" ", var_width, right = TRUE), "|Excl:")
                data_col = paste0(hash, sfill("# Exclusive ", var_width, right = TRUE), "|")
            } else {
                data_col = paste0(hash, sfill("# Exclusive ", var_width, right = TRUE),
                                  "|", paste0(letters[data_id], collapse = ":"))
            }

            # formatting the NAs
            for(i in 1:(ncol(pair) - 2)){
                pair[[i]] = fsignif(pair[[i]])
            }
            pair[is.na(pair)] = " "
            pair = as.matrix(pair)

            # adding the data col
            pair = cbind(data_col, pair)

            print_mat = insert_row_after_id(print_mat, pair, n_x + 2)

        }

        # Formatting
        print_mat = print_mat[, -(n_x + 2)]
        print_mat[, 1] = sfill(print_mat[, 1], right = TRUE)
        for(j in 2:(n_x + 1)){
            print_mat[, j] = sfill(print_mat[, j])
        }

        # Last column
        common = print_mat[, n_x + 2]
        is_num = grepl("\\d", common)
        common[is_num] = sfill(fsignif(as.numeric(common[is_num])))
        common[is_num] = paste0("# Common: ", common[is_num])
        print_mat[, n_x + 2] = common
    }

    vec2print = apply(print_mat, 1, paste, collapse = " ")

    cat(vec2print, sep = "\n")

}


#' Formatted object size
#'
#' Tools that returns a formatted object size, where the appropriate unit is automatically chosen.
#'
#' @param x Any R object.
#' @param ... Not currently used.
#'
#' @return
#' Returns a character scalar.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' osize(iris)
#'
#' data(trade)
#' osize(trade)
#'
#'
osize = function(x){

    size = as.numeric(utils::object.size(x))
    n = log10(size)

    if (n < 3) {
        res = paste0(size, " Octets.")
    } else if (n < 6) {
        res = paste0(fsignif(size/1000), " Ko.")
    } else {
        res = paste0(fsignif(size/1e+06), " Mo.")
    }

    class(res) = "osize"
    res
}

#' @rdname osize
print.osize = function(x, ...){
    cat(x, "\n")
}


#' Randomly draws observations from a data set
#'
#' This function is useful to check a data set. It gives a random number of rows of the input data set.
#'
#' @param x A data set: either a vector, a matrix or a data frame.
#' @param n The number of random rows/elements to sample randomly.
#' @param previous Logical scalar. Whether the results of the previous draw should be returned.
#'
#' @return
#' A data base (resp vector) with `n` rows (resp elements).
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' sample_df(iris)
#'
#' sample_df(iris, previous = TRUE)
#'
sample_df = function(x, n = 10, previous = FALSE){

    check_arg(n, "integer scalar")
    check_arg(previous, "logical scalar")

    if(MISSNULL(x)){
        if(missing(x)) stop("The argument `x` cannot be missing, please provide it.")
        if(is.null(x)) stop("The argument `x` must be a vector, matrix or data.frame. Currently it is NULL.")
    }

    all_draws = getOption("fixest_sample_df")
    x_dp = deparse_long(substitute(x))

    make_draw = TRUE
    if(previous){
        if(x_dp %in% names(all_draws)){
            draw__ = all_draws[[x_dp]]
            make_draw = FALSE
        } else {
            warning("No previous draw found for this data. Making a new draw.")
        }
    }

    is_unidim = is.null(dim(x)) || inherits(x, "table")
    n_data = if(is_unidim) length(x) else nrow(x)


    n = min(n, n_data)

    if(make_draw){
        # complicated var name to avoid issue with data.table variables
        draw__ = sample(n_data, n)
    }

    # saving
    all_draws[[x_dp]] = draw__
    options(fixest_sample_df = all_draws)

    # returning
    if(is_unidim){
        return(x[draw__])
    } else {
        return(x[draw__, ])
    }
}

len_unique = function(x, nthreads = getFixest_nthreads()){

    if(!is.vector(x)){
        return(length(unique(x)))
    }

    ANY_NA = FALSE
    if(is.numeric(x)){
        info_na = cpp_which_na_inf(x, nthreads)
        if(info_na$any_na_inf){
            ANY_NA = TRUE
            x = x[!info_na$is_na_inf]
        }
    } else {
        if(anyNA(x)){
            ANY_NA = TRUE
            x = x[!is.na(x)]
        }
    }

    if(length(x) == 0){
        return(0)
    }

    x_quf = quickUnclassFactor(x, addItem = TRUE, sorted = FALSE)

    length(x_quf$items) + ANY_NA
}


#' Formatted dimension
#'
#' Prints the dimension of a data set, in an user-readable way
#'
#' @param x An R object, usually a data.frame (but can also be a vector).
#'
#' @return
#' It does not return anything, the output is directly printed on the console.
#'
#' @author Laurent Berge
#'
#' @examples
#'
#' fdim(iris)
#'
#' fdim(iris$Species)
#'
#'
#'
fdim = function(x){
    if(!missing(x)){
        if(!is.null(dim(x))) {
            info = fsignif(dim(x))
            cat(info[1], " row", plural(info[1]), " and ", info[2], " column", plural(info[2]), "\n", sep = "")
        } else if (!is.null(length(x))) {
            cat(fsignif(length(x)), " obs\n")
        }
    }
}


#' Fast transform of any type of vector(s) into an integer vector
#'
#' Tool to transform any type of vector, or even combination of vectors, into an integer vector ranging from 1 to the number of unique values. This actually creates an unique identifier vector.
#'
#' @param ... Vectors of any type, to be transformed in integer.
#' @param sorted Logical, default is `FALSE`. Whether the integer vector should make reference to sorted values?
#' @param add_items Logical, default is `FALSE`. Whether to add the unique values of the original vector(s). If requested, an attribute `items` is created containing the values (alternatively, they can appear in a list if `items.list=TRUE`).
#' @param items.list Logical, default is `FALSE`. Only used if `add_items=TRUE`. If `TRUE`, then a list of length 2 is returned with `x` the integer vector and `items` the vector of items.
#' @param multi.df Logical, default is `FALSE`. If `TRUE` then a data.frame listing the unique elements is returned in the form of a data.frame. Ignored if `add_items = FALSE`.
#' @param multi.join Character scalar used to join the items of multiple vectors. The default is `"_"`. Ignored if `add_items = FALSE`.
#' @param internal Logical, default is `FALSE`. For programming only. If this function is used within another function, setting `internal = TRUE` is needed to make the evaluation of `...` valid. End users of `to_integer` should not care.
#'
#'
#' @return
#' Reruns a vector of the same length as the input vectors.
#' If `add_items=TRUE` and `items.list=TRUE`, a list of two elements is returned: `x` being the integer vector and `items` being the unique values to which the values in `x` make reference.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' x1 = iris$Species
#' x2 = as.integer(iris$Sepal.Length)
#'
#' # transforms the species vector into integers
#' to_integer(x1)
#'
#' # To obtain the "items":
#' to_integer(x1, add_items = TRUE)
#' # same but in list form
#' to_integer(x1, add_items = TRUE, items.list = TRUE)
#'
#' # transforms x2 into an integer vector from 1 to 4
#' to_integer(x2, add_items = TRUE)
#'
#' # To have the sorted items:
#' to_integer(x2, add_items = TRUE, sorted = TRUE)
#'
#' # The result can safely be used as an index
#' res = to_integer(x2, add_items = TRUE, sorted = TRUE, items.list = TRUE)
#' all(res$items[res$x] == x2)
#'
#'
#' #
#' # Multiple vectors
#' #
#'
#' to_integer(x1, x2, add_items = TRUE)
#'
#' # You can use multi.join to handle the join of the items:
#' to_integer(x1, x2, add_items = TRUE, multi.join = "; ")
#'
to_integer = function(..., sorted = FALSE, add_items = FALSE, items.list = FALSE,
                      multi.df = FALSE, multi.join = "_", internal = FALSE){

    if(!internal) check_arg(..., "vector mbt")
    check_arg(sorted, add_items, items.list, "logical scalar")
    check_arg(multi.join, "character scalar")

    dots = list(...)

    # Removing NAs
    Q = length(dots)
    n_all = lengths(dots)
    n = n_all[1]

    if(length(unique(n_all)) != 1) stop("All elements in `...` should be of the same length (current lenghts are ", enumerate_items(n_all), ").")

    is_na = is.na(dots[[1]])
    for(q in seq(from = 2, length.out = Q - 1)){
        is_na = is_na | is.na(dots[[q]])
    }

    ANY_NA = FALSE
    if(any(is_na)){
        ANY_NA = TRUE

        if(all(is_na)){
            message("NOTE: All values are NA.")
            res = rep(NA, n)
            if(add_items){
                if(items.list){
                    res = list(x = res, items = NA)
                } else {
                    attr(res, "items") = NA
                }
            }

            return(res)
        }

        for(q in 1:Q) dots[[q]] = dots[[q]][!is_na]
    }

    #
    # Creating the ID
    #

    if(Q == 1){
        if(sorted && !is.numeric(dots[[1]]) && !is.character(dots[[1]])){
            # general way => works for any type with a sort method
            f = dots[[1]]
            res_raw = quickUnclassFactor(f, addItem = TRUE, sorted = FALSE)
            obs_1st = cpp_get_first_item(res_raw$x, length(res_raw$items))
            f_unik = f[obs_1st]
            f_order = order(f_unik)
            x_new = order(f_order)[res_raw$x]
            if(add_items){
                items_new = f_unik[f_order]
                res = list(x = x_new, items = items_new)
            } else {
                res = x_new
            }

        } else {
            res = quickUnclassFactor(dots[[1]], addItem = add_items, sorted = sorted)
        }

    } else {

        QUF_raw = list()
        for(q in 1:Q){
            QUF_raw[[q]] = quickUnclassFactor(dots[[q]], sorted = FALSE, addItem = TRUE)
        }

        # Then we combine
        power = floor(1 + log10(sapply(QUF_raw, function(x) length(x$items))))

        is_large = sum(power) > 14
        if(is_large){
            # 15 Aug 2021, finally found a solution. It was so obvious with hindsight...
            QUF_raw_value = lapply(QUF_raw, `[[`, 1)
            order_index = do.call(order, QUF_raw_value)
            index = cpp_combine_clusters(QUF_raw_value, order_index)
        } else {
            # quicker, but limited by the precision of doubles
            index = QUF_raw[[1]]$x
            for(q in 2:Q){
                index = index + QUF_raw[[q]]$x*10**sum(power[1:(q-1)])
            }
        }

        res = quickUnclassFactor(index, addItem = add_items || sorted, sorted = sorted)

        if(add_items || sorted){
            # we re order appropriately
            # f prefix means factor

            obs_1st = cpp_get_first_item(res$x, length(res$items))

            f_all = list()
            for(q in 1:Q){
                f_all[[q]] = dots[[q]][obs_1st]
            }

            f_order = do.call("order", f_all)

            x_new = order(f_order)[res$x]

            if(multi.df){
                # Putting into a DF => we take care of names
                mc_dots = match.call(expand.dots = FALSE)[["..."]]
                n_dots = length(mc_dots)
                mc_dots_names = names(mc_dots)
                if(is.null(mc_dots_names)) mc_dots_names = character(n_dots)

                my_names = character(n_dots)
                for(q in 1:n_dots){
                    if(nchar(mc_dots_names[q]) > 0){
                        my_names[q] = mc_dots_names[q]
                    } else {
                        my_names[q] = deparse_long(mc_dots[[q]])
                    }
                }

                names(f_all) = my_names

                f_df = as.data.frame(f_all)
                items_new = f_df[f_order, , drop = FALSE]
                row.names(items_new) = 1:nrow(items_new)
            } else {
                # we "paste" them
                arg_list = f_all
                arg_list$sep = multi.join
                f_char = do.call("paste", arg_list)
                items_new = f_char[f_order]
            }

            if(add_items){
                res = list(x = x_new, items = items_new)
            } else {
                res = x_new
            }
        }
    }

    if(ANY_NA){
        if(is.list(res)){
            x = res$x
        } else {
            x = res
        }

        x_na = rep(NA, n)
        x_na[!is_na] = x

        if(is.list(res)){
            res$x = x_na
        } else {
            res = x_na
        }

    }

    if(add_items && isFALSE(items.list)){
        res_tmp = res$x
        attr(res_tmp, "items") = res$items
        res = res_tmp
    }

    res
}






















































