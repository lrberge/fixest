#------------------------------------------------------------------------------#
# Created: 2026-01-12
# ~: combine all results into a single table
#------------------------------------------------------------------------------#

library(data.table)

all_results_list = list()
for(lang in c("R", "jl", "py")){
  all_results_list[[lang]] = fread(paste0("./results/results_taxi_", lang, ".csv"))
}

base_results = rbindlist(all_results_list)

fwrite(base_results, "./results/all_results_taxi.csv")

