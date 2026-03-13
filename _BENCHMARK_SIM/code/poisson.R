#------------------------------------------------------------------------------#
# Created: 2026-01-08
# ~: Poisson benchmarks for R functions
#------------------------------------------------------------------------------#


if(getOption("repos") == "@CRAN@") options(repos = "https://cran.rstudio.com")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("fixest", quietly = TRUE)) install.packages("fixest")
if(!requireNamespace("alpaca", quietly = TRUE)) install.packages("alpaca")

all_results_list = list()
for(n_obs in c(1e4, 1e5, 1e6)){
  message(n_obs, appendLF = FALSE)
  
  for(id_rep in 1:3){
    message("|", appendLF = FALSE)
    
    # loading the data 
    filename = paste0("./data/data_n-1e", log10(n_obs), "_rep-", id_rep, ".gz")
    base = data.table::fread(filename, showProgress = FALSE)
    
    # we exponentiate the dep. var
    base$y = exp(base$y)
    
    for(dgp in c("simple", "difficult")){
      message(substr(dgp, 1, 1), appendLF = FALSE)
      
      for(n_fe in 2:3){
        message(".", appendLF = FALSE)
        
        # NOTE: alpaca and fixest use the same formula
        if(dgp == "simple"){
          fml = switch(as.character(n_fe), 
                       "2" = y ~ x1 + x2 | indiv_id + firm_id,
                       "3" = y ~ x1 + x2 | indiv_id + firm_id + year)
        } else if(dgp == "difficult"){
          fml = switch(as.character(n_fe), 
                       "2" = y ~ x1 + x2 | indiv_id + firm_id_difficult,
                       "3" = y ~ x1 + x2 | indiv_id + firm_id_difficult + year)
        }
        
        #
        # fixest::fepois
        #
        
        time_start = Sys.time()
        est = fixest::fepois(fml, base)
        time_fepois = as.numeric(Sys.time() - time_start, units = "secs")
        
        #
        # alpaca::feglm
        #
        
        # we don't run for difficult large data sets because the run time is > 1h
        if(dgp == "difficult" && n_obs >= 1e6){
          time_feglm = NA
          
        } else {
          time_start = Sys.time()
          est = alpaca::feglm(fml, base, family = poisson())
          time_feglm = as.numeric(Sys.time() - time_start, units = "secs")

        } 
        
        
        #
        # saving
        #
        
        result = data.frame(
          n = n_obs, id_rep = id_rep, dgp = dgp, n_fe = n_fe, 
          method = c("fixest::fepois", "alpaca::feglm"),
          time = c(time_fepois, time_feglm)
        )
        
        all_results_list[[length(all_results_list) + 1]] = result
      }
      
      
    }
    
  }
  message("|")
}

all_results = do.call(rbind, all_results_list)

data.table::fwrite(all_results, "./results/results_poisson_R.csv")




