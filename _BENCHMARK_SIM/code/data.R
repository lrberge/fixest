#------------------------------------------------------------------------------#
# Created: 2026-01-08
# ~: Data generation (1st step of the benchmarks)
#------------------------------------------------------------------------------#

# dependencies:
# - data.table, used to write compressed CSVs

####
#### Function to generate the data ####
####

base_dgp = function(n = 1000) {
  
  nb_year = 10
  nb_indiv_per_firm = 23
  
  nb_indiv = round(n / nb_year)
  nb_firm = round(nb_indiv / nb_indiv_per_firm)
  indiv_id = rep(1:nb_indiv, each = nb_year)
  year = rep(1:nb_year, times = nb_indiv)

  firm_id_simple = sample(1:nb_firm, n, TRUE)
  firm_id_difficult = rep(1:nb_firm, length.out = n)

  unit_fe = rnorm(nb_indiv)[indiv_id]
  year_fe = rnorm(nb_year)[year]
  firm_fe = rnorm(nb_firm)[firm_id_simple]
  
  x1 = rnorm(n)
  y = 1 * x1 + 0.05 * x1^2 + firm_fe + unit_fe + year_fe + rnorm(n)
  df = data.frame(
    indiv_id = indiv_id,
    year = year,
    firm_id = firm_id_simple,
    firm_id_difficult = firm_id_difficult,
    x1 = x1,
    x2 = x1^2,
    y = y
  )
  
  # we round everything at the 2nd digit to save disk space when compressing
  # => it divides the size of the data on disk by 3
  df = as.data.frame(lapply(df, round, digits = 2))
  
  return(df)
}

####
#### Loop to construct and save on disk all the data sets ####
####


# number of repetitions
N_REP = 3

# the seed (taken from as.numeric(Sys.time()) at the time of writing)
set.seed(1767871768)

for(n_obs in c(1e4, 1e5, 1e6, 1e7)){
  message(n_obs, appendLF = FALSE)
  # we replicate N_REP times
  for(id_rep in 1:N_REP){
    message(".", appendLF = FALSE)
    
    # we generate the data
    base = base_dgp(n = n_obs)
    
    # we save in a compressed csv file
    # (NOTE: by saving the variables n_obs and id_rep in the 
    #        file name and not as a variables of a larger data set, 
    #        we save some disk space [useful for the 1e7 case])
    
    filename = paste0("./data/data_n-1e", log10(n_obs), "_rep-", id_rep, ".gz")
    data.table::fwrite(base, filename, showProgress = FALSE)
    
  }
  message()
}





