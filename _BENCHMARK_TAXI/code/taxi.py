

# NOTA: the module pyarrow is needed to read parquet files and should be installed
import time
import pandas
import pyfixest

#
# Loading the data
#

base_taxi = pandas.read_parquet("./data/nyc_taxi.parquet")

#
# Estimations
#

#
# Base (3 FEs)
time_start = time.perf_counter()
est_base = pyfixest.feols("tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type", base_taxi)
time_base = time.perf_counter() - time_start

#
# Varying slopes
# => no known equivalent in pyfixest

#
# Multiple y
base_taxi['tip_large'] = 1 * (base_taxi["tip_amount"] >= 5)
time_start = time.perf_counter()
est_multy = pyfixest.feols("tip_amount + tip_large ~ trip_distance + passenger_count | dofw + vendor_id + payment_type", base_taxi)
time_multy = time.perf_counter() - time_start

#
# Multiple VCOVs
time_start = time.perf_counter()
est_multvcov = pyfixest.feols("tip_amount ~ trip_distance + passenger_count | dofw + vendor_id + payment_type", base_taxi)
est_multvcov_clust = est_multvcov.vcov({"CRV1": "vendor_id"})
time_multvcov = time.perf_counter() - time_start


#
# Save
#

base_results = pandas.DataFrame({
  'method': ["pyfixest.feols" for i in range(3)],
  'type': ["base", "multiple-y", "multiple-vcov"],
  'time': [time_base, time_multy, time_multvcov]
})

base_results.to_csv("./results/results_taxi_py.csv", index = False)

