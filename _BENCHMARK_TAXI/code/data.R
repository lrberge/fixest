#------------------------------------------------------------------------------#
# Created: 2026-01-14
# ~: Import the taxi data
#------------------------------------------------------------------------------#

# The origin of the data is this website:
# https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page
# 
# The direct download link of the files was access in 2026-01-14.
# If the download link becomes broken (unlikely but possible), please 
# manually download at the link above.


# needed packages
if(getOption("repos") == "@CRAN@") options(repos = "https://cran.rstudio.com")
if(!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")

#
# downloading the data source
#

# large timeout for downloading a large file
options(timeout = 200)

# We download the first three months of 2012
link = "https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2012-0{MONTH}.parquet"
destination_file = "./data/taxi/data_taxi_2012-0{MONTH}.parquet"

for(month in 1:3){
  
  dest_file = sub("{MONTH}", month, destination_file, fixed = TRUE)
  if(!file.exists(dest_file)){
    # we avoid downloading several times
    message("downloading month = ", month)
    download.file(sub("{MONTH}", month, link, fixed = TRUE),
                  dest_file, mode = "wb")
    
    Sys.sleep(1)
  }
  
}


#
# formatting the data
#

# required library
library(data.table)

# we format the three tables into a single file
all_tables = list()
for(month in 1:3){
  filename = sub("{MONTH}", month, destination_file, fixed = TRUE)
  message("reading ", filename)
  all_tables[[month]] = arrow::read_parquet(filename)
}

base_taxi = rbindlist(all_tables)
rm(all_tables)

# We only keep the variables we will use:
# - vendor_id = VendorID, either :
#   + 1 = Creative Mobile Technologies, LLC
#   + 2 = Curb Mobility, LLC
# - dofw = day of the week, extracted from the date tpep_pickup_datetime
# - trip_distance
# - passenger_count
# - payment_type, either:
#   + 1 = Credit card
#   + 2 = Cash
#   + 3 = No charge
#   + 4 = Dispute
#   + 5 = Unknown
# - tip_amount


base_taxi_small = base_taxi[, .(
  vendor_id = c("CMT", "CM")[VendorID],
  dofw = lubridate::wday(tpep_pickup_datetime),
  trip_distance,
  passenger_count,
  payment_type = c("Credit card", "Cash", "No charge", "Dispute", "Unknown")[payment_type],
  tip_amount
)]


# saving
arrow::write_parquet(base_taxi_small, "./data/nyc_taxi.parquet")

