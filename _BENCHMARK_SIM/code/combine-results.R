#------------------------------------------------------------------------------#
# Created: 2026-01-12
# ~: gathers all the individual results and saves them in a single table
#    after computing the average computing time
#------------------------------------------------------------------------------#


library(data.table)
all_results_list = list()
for(lang in c("R", "jl", "py")){
  for(estimator in c("ols", "poisson")){
    base_results = fread(paste0("./results/results_", estimator, "_", lang, ".csv"))
    if(lang %in% c("jl", "py")){
      # we need to drop the "burn in"
      base_results = base_results[-1]
    }
    base_results[, estimator := estimator]
    all_results_list[[length(all_results_list) + 1]] = base_results
  }
}

base_all_results = rbindlist(all_results_list)

base_avg_time = base_all_results[, .(time = mean(time, na.rm = TRUE)), 
                                 keyby = .(estimator, dgp, n_fe, n, method)]

fwrite(base_avg_time, "./results/all_results.csv")


