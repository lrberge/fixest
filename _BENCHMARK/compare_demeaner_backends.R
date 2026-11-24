#----------------------------------------------#
# Author: fixest contributors
# Purpose: Compare MAP vs within demeaner by
#          printing etable() and timings
#----------------------------------------------#

parse_args = function(args){
  out = list(
    seed = 1L,
    pow = 6L
  )

  for(a in args){
    if(grepl("^--seed=", a)){
      out$seed = as.integer(sub("^--seed=", "", a))
      next
    }
    if(grepl("^--pow=", a)){
      out$pow = as.integer(sub("^--pow=", "", a))
      next
    }
    stop("Unknown argument: ", a)
  }

  if(!is.finite(out$seed)){
    stop("--seed must be finite.")
  }
  if(!is.finite(out$pow) || !(out$pow %in% 4:7)){
    stop("--pow must be an integer in 4:7.")
  }

  out
}

load_fixest = function(){
  if(requireNamespace("pkgload", quietly = TRUE) &&
     file.exists("DESCRIPTION") && dir.exists("R")){
    pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
  } else {
    library(fixest)
  }
}

make_data_benchmark_difficult = function(seed, pow){
  # Adapted from _BENCHMARK/Data generation.R ("Difficult Data")
  set.seed(seed)
  n = 10^pow
  nb_indiv = n / 20
  nb_firm = round(n / 160)
  nb_year = round(n^0.3)

  id_indiv = sample.int(nb_indiv, n, replace = TRUE)
  id_firm = pmin(sample(0:20, n, TRUE) + pmax(1, id_indiv %/% 8 - 10), nb_firm)
  id_year = sample.int(nb_year, n, replace = TRUE)

  x1 = 5 * cos(id_indiv) + 5 * sin(id_firm) + 5 * sin(id_year) + runif(n)
  x2 = cos(id_indiv) + sin(id_firm) + sin(id_year) + rnorm(n)
  y = 3 * x1 + 5 * x2 + cos(id_indiv) + cos(id_firm)^2 + sin(id_year) + rnorm(n)
  w = 0.5 + runif(n)

  data.frame(
    y = y,
    x1 = x1,
    x2 = x2,
    id_indiv = id_indiv,
    id_firm = id_firm,
    id_year = id_year,
    w = w
  )
}

elapsed_seconds = function(expr){
  t0 = proc.time()
  val = force(expr)
  t1 = proc.time()
  list(value = val, elapsed = as.numeric((t1 - t0)[["elapsed"]]))
}

run_pair = function(label, fml, data, weights = NULL){
  args = list(fml = fml, data = data, demeaner = MapDemeaner())
  if(!is.null(weights)){
    args$weights = weights
  }
  map_run = elapsed_seconds(do.call(feols, args))

  args$demeaner = LsmrDemeaner(preconditioner = "additive")
  within_run = elapsed_seconds(do.call(feols, args))

  cat("\n==============================\n")
  cat("Model:", label, "\n")
  cat("==============================\n")
  print(etable(map_run$value, within_run$value,
               headers = c("MAP", "within"), digits = 6))

  data.frame(
    model = label,
    backend = c("MAP", "within"),
    elapsed_seconds = c(map_run$elapsed, within_run$elapsed),
    stringsAsFactors = FALSE
  )
}

main = function(){
  cfg = parse_args(commandArgs(trailingOnly = TRUE))

  if(!requireNamespace("withinr", quietly = TRUE) ||
     !exists("LsmrOptions", envir = asNamespace("withinr"), inherits = FALSE)){
    stop("A current 'withinr' with LsmrOptions() is required. Install it first.")
  }

  load_fixest()
  setFixest_notes(FALSE)

  dat = make_data_benchmark_difficult(cfg$seed, cfg$pow)

  cat("DGP: benchmark difficult data\n")
  cat("seed=", cfg$seed, ", n=10^", cfg$pow, "\n", sep = "")

  timings = list()
  i = 1L

  timings[[i]] = run_pair(
    label = "feols_1fe_weighted",
    fml = y ~ x1 + x2 | id_indiv,
    data = dat,
    weights = ~w
  ); i = i + 1L

  timings[[i]] = run_pair(
    label = "feols_2fe_weighted",
    fml = y ~ x1 + x2 | id_indiv + id_firm,
    data = dat,
    weights = ~w
  ); i = i + 1L

  timings[[i]] = run_pair(
    label = "feols_3fe_weighted",
    fml = y ~ x1 + x2 | id_indiv + id_firm + id_year,
    data = dat,
    weights = ~w
  ); i = i + 1L

  timings[[i]] = run_pair(
    label = "feols_3fe_unweighted",
    fml = y ~ x1 + x2 | id_indiv + id_firm + id_year,
    data = dat
  )

  timing_df = do.call(rbind, timings)
  rownames(timing_df) = NULL

  cat("\n==============================\n")
  cat("Timings (seconds)\n")
  cat("==============================\n")
  print(timing_df, row.names = FALSE, digits = 6)
}

main()
