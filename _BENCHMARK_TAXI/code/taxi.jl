

#
# Loading the data 
#

# reading parquet.... there is a bug that prevents the direct read
# https://github.com/JuliaIO/Parquet.jl/issues/191

using DuckDB, Parquet, CSV, DataFrames

path = "./data/nyc_taxi.parquet"
con = DuckDB.DBInterface.connect(DuckDB.DB, ":memory:")
base_taxi = DataFrame(DuckDB.DBInterface.execute(con, "SELECT * FROM read_parquet('$path')"));

#
# Estimations
#

using FixedEffectModels, Vcov

# burn-in
est = reg(base_taxi, @formula(tip_amount ~ trip_distance))

#
# Base (3 FEs)
time_start = time()
est_base = reg(base_taxi, @formula(tip_amount ~ trip_distance + passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)));
time_base = time() - time_start

#
# Varying slopes
time_start = time()
est_VS = reg(base_taxi, @formula(tip_amount ~ passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)&trip_distance));
time_VS = time() - time_start

#
# Multiple y
base_taxi.tip_large = .+(base_taxi.tip_amount .> 5);
time_start = time()
est_multy_1 = reg(base_taxi, @formula(tip_amount ~ trip_distance + passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)));
est_multy_2 = reg(base_taxi, @formula(tip_large ~ trip_distance + passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)));
time_multy = time() - time_start

#
# Multiple VCOVs
time_start = time()
est_multvcov_1 = reg(base_taxi, @formula(tip_amount ~ trip_distance + passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)));
est_multvcov_2 = reg(base_taxi, @formula(tip_amount ~ trip_distance + passenger_count + fe(dofw) + fe(vendor_id) + fe(payment_type)),
                     Vcov.cluster(:vendor_id));
time_multvcov = time() - time_start


#
# Save
#

all_results = DataFrame(
  method = "FixedEffectModels::reg",
  type = ["base", "varying-slopes", "multiple-y", "multiple-vcov"],
  time = [time_base, time_VS, time_multy, time_multvcov]
)

CSV.write("results/results_taxi_jl.csv", all_results)
