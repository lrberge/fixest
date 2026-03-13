#------------------------------------------------------------------------------#
# Created: 2026-01-12
# ~: Taxi benchmark for R
#------------------------------------------------------------------------------#

# needed packages
if(getOption("repos") == "@CRAN@") options(repos = "https://cran.rstudio.com")
if(!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
if(!requireNamespace("fixest", quietly = TRUE)) install.packages("fixest")

#
# Loading the data
#

base_taxi = arrow::read_parquet("./data/nyc_taxi.parquet")

#
# Estimations
#

library(fixest)

setFixest_estimation(data = base_taxi)

#
# Base (3 FEs)
time_start = Sys.time()
est_base = feols(tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type)
time_base = Sys.time() - time_start

#
# Varying slopes
time_start = Sys.time()
est_VS = feols(tip_amount ~ passenger_count | dofw + vendor_id + payment_type[trip_distance])
time_VS = Sys.time() - time_start

#
# Multiple y
base_taxi$tip_large = +(base_taxi$tip_amount >= 5)
time_start = Sys.time()
est_multy = feols(c(tip_amount, tip_large) ~ trip_distance + passenger_count | dofw + vendor_id + payment_type)
time_multy = Sys.time() - time_start

#
# Multiple VCOVs
time_start = Sys.time()
est_multvcov = feols(tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type)
est_multvcov_clust = summary(est_multvcov, vcov = vcov_cluster("vendor_id"))
time_multvcov = Sys.time() - time_start


#
# Save
#

base_results = data.frame(
  method = "fixest::feols",
  type = c("base", "varying-slopes", "multiple-y", "multiple-vcov"),
  time = c(time_base, time_VS, time_multy, time_multvcov)
)

write.csv(base_results, "./results/results_taxi_R.csv", row.names = FALSE)







