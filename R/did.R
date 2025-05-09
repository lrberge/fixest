#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Thu Apr 29 11:24:49 2021
# ~: DiD functions
#----------------------------------------------#




#' Sun and Abraham interactions
#'
#' User-level method to implement staggered difference-in-difference estimations a la Sun 
#' and Abraham (Journal of Econometrics, 2021).
#'
#'
#' @param cohort A vector representing the cohort. It should represent the period at 
#' which the treatment has been received (and thus be fixed for each unit).
#' @param period A vector representing the period. It can be either a relative time period 
#' (with negative values representing the before the treatment and positive values 
#' after the treatment), or a regular time period. In the latter case, the relative 
#' time period will be created from the cohort information (which represents the time at 
#' which the treatment has been received).
#' @param ref.c A vector of references for the cohort. By default the never treated 
#' cohorts are taken as reference and the always treated are excluded from the estimation. 
#' You can add more references with this argument, which means that dummies will not be 
#' created for them (but they will remain in the estimation).
#' @param ref.p A vector of references for the (relative!) period. By default the 
#' first relative period (RP) before the treatment, i.e. -1, is taken as reference. 
#' You can instead use your own references (i.e. RPs for which dummies will not be 
#' created -- but these observations remain in the sample). Please note that you will 
#' need at least two references. You can use the special variables `.F` and `.L` to 
#' access the first and the last relative periods.
#' @param att Logical, default is `FALSE`. If `TRUE`: then the total average treatment 
#' effect for the treated is computed (instead of the ATT for each relative period).
#' @param no_agg Logical, default is `FALSE`. If `TRUE`: then there is no aggregation, 
#' leading to the estimation of all `cohort x time to treatment` coefficients.
#' @param bin A list of values to be grouped, a vector, or the special value `"bin::digit"`. 
#' The binning will be applied to both the cohort and the period (to bin them separately, 
#' see `bin.c` and `bin.p`). To create a new value from old values, 
#' use `bin = list("new_value"=old_values)` with `old_values` a vector of 
#' existing values. It accepts regular expressions, but they must start with an `"@"`, 
#' like in `bin="@Aug|Dec"`. The names of the list are the new names. If the new 
#' name is missing, the first value matched becomes the new name. Feeding in a vector is 
#' like using a list without name and only a single element. If the vector is numeric, 
#' you can use the special value `"bin::digit"` to group every `digit` element. 
#' For example if `x` represent years, using `bin="bin::2"` create bins of two years. 
#' Using `"!bin::digit"` groups every digit consecutive values starting from the first value. 
#' Using `"!!bin::digit"` is the same bu starting from the last value. In both cases, 
#' `x` is not required to be numeric.
#' @param bin.rel A list or a vector defining which values to bin. Only applies to the 
#' relative periods and *not* the cohorts. Please refer to the help of the argument 
#' `bin` to understand the different ways to do the binning (or look at the help 
#' of [`bin`]).
#' @param bin.c A list or a vector defining which values to bin. Only applies to the cohort. 
#' Please refer to the help of the argument `bin` to understand the different ways to 
#' do the binning (or look at the help of [`bin`]).
#' @param bin.p A list or a vector defining which values to bin. Only applies to the period. 
#' Please refer to the help of the argument `bin` to understand the different ways to 
#' do the binning (or look at the help of [`bin`]).
#'
#' @details
#' This function creates a matrix of `cohort x relative_period` interactions, and if used within 
#' a `fixest` estimation, the coefficients will automatically be aggregated to obtain the ATT 
#' for each relative period. In practice, the coefficients are aggregated with the 
#' [`aggregate.fixest`] function whose argument `agg` is automatically set to the appropriate 
#' value.
#'
#' The SA method requires relative periods (negative/positive for before/after the treatment). 
#' Either the user can compute the RP (relative periods) by his/her own, either the RPs 
#' are computed on the fly from the periods and the cohorts (which then should represent 
#' the treatment period).
#'
#' The never treated, which are the cohorts displaying only negative RPs are used as references 
#' (i.e. no dummy will be constructed for them). On the other hand, the always treated are 
#' removed from the estimation, by means of adding NAs for each of their observations.
#'
#' If the RPs have to be constructed on the fly, any cohort that is not present in the 
#' period is considered as never treated. This means that if the period ranges from 
#' 1995 to 2005, `cohort = 1994` will be considered as never treated, although it 
#' should be considered as always treated: so be careful.
#'
#' If you construct your own relative periods, the controls cohorts should have only negative RPs.
#'
#' @section Binning:
#'
#' You can bin periods with the arguments `bin`, `bin.c`, `bin.p` and/or `bin.rel`.
#'
#' The argument `bin` applies both to the original periods and cohorts (the cohorts will also 
#' be binned!). This argument only works when the `period` represent "calendar" periods 
#' (not relative ones!).
#'
#' Alternatively you can bin the periods with `bin.p` (either "calendar" or relative); or 
#' the cohorts with `bin.c`.
#'
#' The argument `bin.rel` applies only to the relative periods (hence not to the cohorts) once 
#' they have been created.
#'
#' To understand how binning works, please have a look at the help and examples of the 
#' function [`bin`].
#'
#' Binning can be done in many different ways: just remember that it is not because it is 
#' possible that it does makes sense!
#'
#' @author
#' Laurent Berge
#'
#' @return
#' If not used within a `fixest` estimation, this function will return a matrix of 
#' interacted coefficients.
#'
#' @examples
#'
#' # Simple DiD example
#' data(base_stagg)
#' head(base_stagg)
#'
#' # Note that the year_treated is set to 1000 for the never treated
#' table(base_stagg$year_treated)
#' table(base_stagg$time_to_treatment)
#'
#' # The DiD estimation
#' res_sunab = feols(y ~ x1 + sunab(year_treated, year) | id + year, base_stagg)
#' etable(res_sunab)
#'
#' # By default the reference periods are the first year and the year before the treatment
#' # i.e. ref.p = c(-1, .F); where .F is a shortcut for the first period.
#' # Say you want to set as references the first three periods on top of -1
#'
#' res_sunab_3ref = feols(y ~ x1 + sunab(year_treated, year, ref.p = c(.F + 0:2, -1)) |
#'                          id + year, base_stagg)
#'
#' # Display the two results
#' iplot(list(res_sunab, res_sunab_3ref))
#'
#' # ... + show all refs
#' iplot(list(res_sunab, res_sunab_3ref), ref = "all")
#'
#'
#' #
#' # ATT
#' #
#'
#' # To get the total ATT, you can use summary with the agg argument:
#' summary(res_sunab, agg = "ATT")
#'
#' # You can also look at the total effect per cohort
#' summary(res_sunab, agg = "cohort")
#'
#'
#' #
#' # Binning
#' #
#'
#' # Binning can be done in many different ways
#'
#' # binning the cohort
#' est_bin.c   = feols(y ~ x1 + sunab(year_treated, year, bin.c = 3:2) | id + year, base_stagg)
#'
#' # binning the period
#' est_bin.p   = feols(y ~ x1 + sunab(year_treated, year, bin.p = 3:1) | id + year, base_stagg)
#'
#' # binning both the cohort and the period
#' est_bin     = feols(y ~ x1 + sunab(year_treated, year, bin = 3:1) | id + year, base_stagg)
#'
#' # binning the relative period, grouping every two years
#' est_bin.rel = feols(y ~ x1 + sunab(year_treated, year, bin.rel = "bin::2") | id + year, base_stagg)
#'
#' etable(est_bin.c, est_bin.p, est_bin, est_bin.rel, keep = "year")
#'
#'
sunab = function(cohort, period, ref.c = NULL, ref.p = -1, bin, bin.rel, 
                 bin.c, bin.p, att = FALSE, no_agg = FALSE){
  # LATER:
  # - add id or indiv argument, just to remove always treated
  # - add argument bin.p

  check_arg(cohort, "mbt vector")
  check_arg(period, "mbt vector len(data)", .data = cohort)
  check_arg(ref.c, "NULL vector no na")
  check_arg(att, no_agg, "logical scalar")
  check_arg(bin, bin.c, bin.p, bin.rel, "NULL list | vector")

  cohort_name = deparse_long(substitute(cohort))
  period_name = deparse_long(substitute(period))
  period_name = gsub("^[[:alpha:]][[:alpha:]_\\.]*\\$", "", period_name)

  # Finding out what kind of data that is

  is_bin = !missnull(bin)
  is_bin.c = !missnull(bin.c)
  is_bin.p = !missnull(bin.p)

  if(is_bin && (is_bin.c || is_bin.p)){
    stop("You cannot have the argument 'bin' with the arguments 'bin.p' or 'bin.c' at the same time. Use only the latter.")
  }

  # NAness (big perf hit)
  n_origin = length(cohort)
  IS_NA = which(is.na(cohort) | is.na(period))
  ANY_NA = length(IS_NA) > 0
  if(ANY_NA){
    cohort = cohort[-IS_NA]
    period = period[-IS_NA]
  }

  n = length(cohort)

  # Case 1
  # cohort: typically the year of treatment
  # period: relative period
  # we don't need to do anything


  # Case 2
  # cohort: year of treatment
  # period: the year
  # => we compute the relative period, we exclude the never/always treated

  # We find out the never/always treated with:
  #  absence of variance in the relative time

  period_unik = unique(period)
  cohort_unik = unique(cohort)

  if(is_bin.c){
    cohort = bin_factor(bin.c, cohort, cohort_name)
    cohort_unik = unique(cohort)
  }

  if(is_bin.p){
    period = bin_factor(bin.p, period, period_name)
    period_unik = unique(period)
  }

  # CASE 1
  is_CASE_1 = FALSE
  if(is.numeric(period) && 0 %in% period_unik && min(period_unik) < 0 && max(period_unik) > 0){
    # Case 1 => we don't need to do anything
    is_CASE_1 = TRUE

    if(is_bin){
      stop("You cannot use 'bin' when the argument 'period' contains relative periods. To use 'bin', 'period' should represent \"calendar\" periods.")
    }

  } else {

    if(is_bin){
      period = bin_factor(bin, period, period_name)
      cohort = bin_factor(bin, cohort, cohort_name, no_error = TRUE)

      period_unik = unique(period)
      cohort_unik = unique(cohort)
    }

    # CASE 2 => construction of the relative period
    refs = setdiff(cohort_unik, period_unik)

    if(length(refs) == length(cohort_unik)){
      stop("Problem in the creation of the relative time periods. We expected the cohort to be the treated period, yet not a single 'cohort' value was found in 'period'.")
    }

    qui_keep = which(!cohort %in% refs)
    cohort_valid = cohort[qui_keep]
    period_valid = period[qui_keep]

    if(is.numeric(period_valid) || is.numeric(cohort_valid)){

      # simplest case
      if(!is.numeric(cohort_valid)) cohort_valid = as.numeric(cohort_valid)
      if(!is.numeric(period_valid)) period_valid = as.numeric(period_valid)

      rel_period = period_valid - cohort_valid
    } else {
      # difficult case
      sunik_period = sort(unique(period_valid))
      dict_period = seq_along(sunik_period)
      names(dict_period) = sunik_period

      period_valid = dict_period[as.character(period_valid)]
      cohort_valid = dict_period[as.character(cohort_valid)]

      rel_period = period_valid - cohort_valid
    }

    # I put a negative number so that they are not considered as always treated
    new_period = rep(-1, n)
    new_period[qui_keep] = rel_period

    period = new_period
  }

  # Now the reference expressed in relative periods
  .F = period_min = min(period)
  .L = period_max = max(period)
  period_list = list(.F = period_min, .L = period_max)
  check_set_arg(ref.p, "evalset integer vector no na", .data = period_list)
  if(missing(ref.p)) ref.p = ref.p # One of the oddest line of code I ever wrote ;-)

  #
  #  we find out the never/always treated
  #

  cohort_int = quickUnclassFactor(cohort)
  c_order = order(cohort_int)
  info = cpp_find_never_always_treated(cohort_int[c_order], period[c_order])

  if(!is.null(ref.c)){
    qui_drop = which(cohort_int %in% info$ref | period %in% ref.p | cohort %in% ref.c)
  } else {
    qui_drop = which(cohort_int %in% info$ref | period %in% ref.p)
  }
  qui_NA = info$always_treated

  # All references have been removed => pure i() without ref
  cohort = cohort[-qui_drop]
  period = period[-qui_drop]

  if(!missing(bin.rel)){
    # we bin on the relative period
    period = bin_factor(bin.rel, period, "relative period")
  }

  res_raw = i(factor_var = period, f2 = cohort, f_name = period_name)

  # We extend the matrix

  if(ANY_NA){
    # We DON'T remove the references from the estimation
    res = matrix(NA_real_, nrow = n_origin, ncol = ncol(res_raw), 
                 dimnames = list(NULL, colnames(res_raw)))
    res[-IS_NA, ][qui_drop, ] = 0
    res[-IS_NA, ][-qui_drop, ] = res_raw
    if(length(qui_NA) > 0){
      res[-IS_NA, ][qui_NA, ] = NA_real_
    }
  } else {
    res = matrix(0, nrow = n_origin, ncol = ncol(res_raw), 
                 dimnames = list(NULL, colnames(res_raw)))
    res[-qui_drop, ] = res_raw
    if(length(qui_NA) > 0){
      res[qui_NA, ] = NA_real_
    }
  }


  # We add the agg argument to GLOBAL_fixest_mm_info
  if(!no_agg){
    is_GLOBAL = FALSE
    for(where in 1:min(8, sys.nframe())){
      if(exists("GLOBAL_fixest_mm_info", parent.frame(where))){
        GLOBAL_fixest_mm_info = get("GLOBAL_fixest_mm_info", parent.frame(where))
        is_GLOBAL = TRUE
        break
      }
    }

    if(is_GLOBAL){
      agg_att = c("ATT" = paste0("\\Q", period_name, "\\E::[[:digit:]]+:cohort"))
      agg_period = paste0("(\\Q", period_name, "\\E)::(-?[[:digit:]]+):cohort")

      if(att){
        agg = agg_att
      } else {
        agg = agg_period

        # We add the attribute containing the appropriate model_matrix_info
        info = list()
        period_unik = sort(unique(c(period, ref.p)))
        info$coef_names_full = paste0(period_name, "::", period_unik)
        info$items = period_unik

        if(length(ref.p) > 0){
          info$ref_id = c(which(info$items %in% ref.p[1]), which(info$items %in% ref.p[-1]))
          info$ref = info$items[info$ref_id]
        }

        info$f_name = period_name

        info$is_num = TRUE
        info$is_inter_num = info$is_inter_fact = FALSE

        attr(agg, "model_matrix_info") = info
      }

      GLOBAL_fixest_mm_info$sunab = list(agg = agg, agg_att = agg_att, 
                                         agg_period = agg_period, ref.p = ref.p)
      # re assignment
      assign("GLOBAL_fixest_mm_info", GLOBAL_fixest_mm_info, parent.frame(where))
    }
  }

  res
}

#' @rdname sunab
sunab_att = function(cohort, period, ref.c = NULL, ref.p = -1){
  sunab(cohort, period, ref.c, ref.p, att = TRUE)
}



#' Aggregates the values of DiD coefficients a la Sun and Abraham
#'
#' Simple tool that aggregates the value of CATT coefficients in staggered 
#' difference-in-difference setups (see details).
#'
#' @param x A `fixest` object.
#' @param agg A character scalar describing the variable names to be aggregated, 
#' it is pattern-based. For [`sunab`] estimations, the following keywords work: "att", 
#' "period", "cohort" and `FALSE` (to have full disaggregation). All variables that 
#' match the pattern will be aggregated. It must be of the form `"(root)"`, the parentheses 
#' must be there and the resulting variable name will be `"root"`. You can add another 
#' root with parentheses: `"(root1)regex(root2)"`, in which case the resulting 
#' name is `"root1::root2"`. To name the resulting variable differently you can pass 
#' a named vector: `c("name" = "pattern")` or `c("name" = "pattern(root2)")`. It's a 
#' bit intricate sorry, please see the examples.
#' @param full Logical scalar, defaults to `FALSE`. If `TRUE`, then all coefficients 
#' are returned, not only the aggregated coefficients.
#' @param use_weights Logical, default is `TRUE`. If the estimation was weighted, 
#' whether the aggregation should take into account the weights. Basically if the 
#' weights reflected frequency it should be `TRUE`.
#' @param ... Arguments to be passed to [`summary.fixest`].
#'
#' @details
#' This is a function helping to replicate the estimator from Sun and Abraham (2021). 
#' You first need to perform an estimation with cohort and relative periods dummies 
#' (typically using the function [`i`]), this leads to estimators of the cohort 
#' average treatment effect on the treated (CATT). Then you can use this function to 
#' retrieve the average treatment effect on each relative period, or for any other way 
#' you wish to aggregate the CATT.
#'
#' Note that contrary to the SA article, here the cohort share in the sample is 
#' considered to be a perfect measure for the cohort share in the population.
#'
#' @return
#' It returns a matrix representing a table of coefficients.
#'
#' @references
#' Liyang Sun and Sarah Abraham, 2021, "Estimating Dynamic Treatment Effects in 
#' Event Studies with Heterogeneous Treatment Effects". Journal of Econometrics.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' #
#' # DiD example
#' #
#'
#' data(base_stagg)
#'
#' # 2 kind of estimations:
#' # - regular TWFE model
#' # - estimation with cohort x time_to_treatment interactions, later aggregated
#'
#' # Note: the never treated have a time_to_treatment equal to -1000
#'
#' # Now we perform the estimation
#' res_twfe = feols(y ~ x1 + i(time_to_treatment, treated,
#'                             ref = c(-1, -1000)) | id + year, base_stagg)
#'
#' # we use the "i." prefix to force year_treated to be considered as a factor
#' res_cohort = feols(y ~ x1 + i(time_to_treatment, i.year_treated,
#'                               ref = c(-1, -1000)) | id + year, base_stagg)
#'
#' # Displaying the results
#' iplot(res_twfe, ylim = c(-6, 8))
#' att_true = tapply(base_stagg$treatment_effect_true,
#'                   base_stagg$time_to_treatment, mean)[-1]
#' points(-9:8 + 0.15, att_true, pch = 15, col = 2)
#'
#' # The aggregate effect for each period
#' agg_coef = aggregate(res_cohort, "(ti.*nt)::(-?[[:digit:]]+)")
#' x = c(-9:-2, 0:8) + .35
#' points(x, agg_coef[, 1], pch = 17, col = 4)
#' ci_low = agg_coef[, 1] - 1.96 * agg_coef[, 2]
#' ci_up = agg_coef[, 1] + 1.96 * agg_coef[, 2]
#' segments(x0 = x, y0 = ci_low, x1 = x, y1 = ci_up, col = 4)
#'
#' legend("topleft", col = c(1, 2, 4), pch = c(20, 15, 17),
#'        legend = c("TWFE", "True", "Sun & Abraham"))
#'
#'
#' # The ATT
#' aggregate(res_cohort, c("ATT" = "treatment::[^-]"))
#' with(base_stagg, mean(treatment_effect_true[time_to_treatment >= 0]))
#'
#' # The total effect for each cohort
#' aggregate(res_cohort, c("cohort" = "::[^-].*year_treated::([[:digit:]]+)"))
#'
#'
aggregate.fixest = function(x, agg, full = FALSE, use_weights = TRUE, ...){
  # Aggregates the value of coefficients

  check_arg(x, "class(fixest) mbt")
  if(isTRUE(x$is_sunab)){
    check_arg(agg, "scalar(character, logical)")
  } else {
    check_arg(agg, "character scalar")
  }

  check_arg(full, "logical scalar")
  # => later => extend it to more than one set of vars to agg

  dots = list(...)
  from_summary = isTRUE(dots$from_summary)

  no_agg = FALSE
  agg_rm = NULL
  check_set_value(agg, "match(att, period, cohort, TRUE) | scalar")
  if(agg %in% c("att", "period", "cohort", "TRUE")){
    if(isTRUE(x$is_sunab)){
      agg_name = names(agg)
      if(agg == "att"){
        agg = x$model_matrix_info$sunab$agg_att
        # we also remove the previous vars
        agg_rm = gsub("E::", "E::-?", agg, fixed = TRUE)
      } else if(agg == "cohort"){
        agg = c("cohort" = "::[^-].*:cohort::(.+)")
        agg_rm = gsub("E::", "E::-?", x$model_matrix_info$sunab$agg_att, fixed = TRUE)
      } else {
        agg = x$model_matrix_info$sunab$agg_period
      }
      if(!is.null(agg_name)) names(agg) = agg_name
    }
  } else if(isFALSE(agg)){
    agg = c("nothing to remove" = "we want all the coefficients")
  }

  is_name = !is.null(names(agg))

  if(!is_name && !grepl("(", agg, fixed = TRUE)){
    stop("Argument 'agg' must be a character in which the pattern to match must be in between parentheses. So far there are no parenthesis: please have a look at the examples.")
  }

  coef = x$coefficients
  cname = names(coef)

  qui = grepl(agg, cname, perl = TRUE)
  if(!any(qui)){
    if(from_summary){
      # We make it silent when aggregate is used in summary
      # => this way we can pool calls to agg even for models that don't have it
      # ==> useful in etable eg
      return(list(coeftable = x$coeftable, model_matrix_info = x$model_matrix_info))
    } else if(no_agg){
      x = summary(x, agg = FALSE, ...)
      return(x$coeftable)
    } else {
      stop("The argument 'agg' does not match any variable.")
    }
  }

  if(!isTRUE(x$summary)){
    x = summary(x, ...)
  }

  cname_select = cname[qui]
  if(is_name){
    root = rep(names(agg), length(cname_select))
    val = gsub(paste0(".*", agg, ".*"), "\\1", cname_select, perl = TRUE)
  } else {
    root = gsub(paste0(".*", agg, ".*"), "\\1", cname_select, perl = TRUE)
    val = gsub(paste0(".*", agg, ".*"), "\\2", cname_select, perl = TRUE)
  }

  V = x$cov.scaled

  mm = model.matrix(x)

  name_df = unique(data.frame(root, val, stringsAsFactors = FALSE))

  c_all = c()
  se_all = c()
  for(i in 1:nrow(name_df)){
    r = name_df[i, 1]
    v = name_df[i, 2]
    v_names = cname_select[root == r & val == v]

    if(use_weights && !is.null(x$weights)){
      shares = colSums(x$weights * sign(mm[, v_names, drop = FALSE]))
    } else {
      shares = colSums(sign(mm[, v_names, drop = FALSE]))
    }

    shares = shares / sum(shares)

    # The coef
    c_value = sum(shares * coef[v_names])

    # The variance
    n = length(v_names)
    s1 = matrix(shares, n, n)
    s2 = matrix(shares, n, n, byrow = TRUE)

    var_value = sum(s1 * s2 * V[v_names, v_names])
    se_value = sqrt(var_value)

    c_all[length(c_all) + 1] = c_value
    se_all[length(se_all) + 1] = se_value
  }

  # th z & p values
  zvalue = c_all/se_all
  pvalue = fixest_pvalue(x, zvalue, V)

  res = cbind(c_all, se_all, zvalue, pvalue)
  if(max(nchar(val)) == 0){
    rownames(res) = name_df[[1]]
  } else {
    rownames(res) = apply(name_df, 1, paste, collapse = "::")
  }

  colnames(res) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  if(full){
    if(!is.null(agg_rm)){
      qui = grepl(agg_rm, cname, perl = TRUE)
    }

    table_origin = x$coeftable
    i_min = min(which(qui)) - 1
    before = if(i_min > 0) table_origin[1:i_min, , drop = FALSE] else NULL

    i_after = (1:nrow(table_origin)) > i_min & !qui
    after = if(any(i_after)) table_origin[i_after, , drop = FALSE] else NULL

    res = rbind(before, res, after)

    attr(res, "type") = attr(table_origin, "type")
  }

  if(from_summary){
    # We add the model_matrix_info needed in iplot()
    mm_info = x$model_matrix_info
    mm_info_agg = attr(agg, "model_matrix_info")
    if(!is.null(mm_info_agg)){
      tmp = list(mm_info_agg)
      for(i in seq_along(mm_info)){
        my_name = names(mm_info)[i]
        if(my_name != ""){
          tmp[[my_name]] = mm_info[[i]]
        } else {
          tmp[[1 + i]] = mm_info[[i]]
        }

      }
      mm_info = tmp
    }

    res = list(coeftable = res, model_matrix_info = mm_info)
  }


  res
}


































