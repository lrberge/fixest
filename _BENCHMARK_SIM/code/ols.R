#------------------------------------------------------------------------------#
# Created: 2026-01-08
# ~: OLS benchmarks for R functions
#------------------------------------------------------------------------------#

# Three dependencies:

if(getOption("repos") == "@CRAN@") options(repos = "https://cran.rstudio.com")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("fixest", quietly = TRUE)) install.packages("fixest")
if(!requireNamespace("lfe", quietly = TRUE)) install.packages("lfe")

# All OLS simulations in a quadruple loop

all_results_list = list()
for(n_obs in c(1e4, 1e5, 1e6, 1e7)){
  message(n_obs, appendLF = FALSE)
  
  for(id_rep in 1:3){
    message("|", appendLF = FALSE)
    
    # loading the data 
    filename = paste0("./data/data_n-1e", log10(n_obs), "_rep-", id_rep, ".gz")
    base = data.table::fread(filename, showProgress = FALSE)
    
    for(dgp in c("simple", "difficult")){
      message(substr(dgp, 1, 1), appendLF = FALSE)
      
      for(n_fe in 2:3){
        message(".", appendLF = FALSE)
        
        #
        # fixest::feols
        #
        
        # NOTE: lfe and fixest use the same formula
        if(dgp == "simple"){
          fml = switch(as.character(n_fe), 
                       "2" = y ~ x1 + x2 | indiv_id + firm_id,
                       "3" = y ~ x1 + x2 | indiv_id + firm_id + year)
        } else if(dgp == "difficult"){
          fml = switch(as.character(n_fe), 
                       "2" = y ~ x1 + x2 | indiv_id + firm_id_difficult,
                       "3" = y ~ x1 + x2 | indiv_id + firm_id_difficult + year)
        }
        
        time_start = Sys.time()
        est = fixest::feols(fml, base)
        time_feols = as.numeric(Sys.time() - time_start, units = "secs")
        
        #
        # lfe::felm
        #
        
        # we don't run felm for difficult large data sets because the run time is > 1h
        if(dgp == "difficult" && n_obs >= 1e6){
          time_felm = NA
          
        } else {
          time_start = Sys.time()
          est = lfe::felm(fml, base)
          time_felm = as.numeric(Sys.time() - time_start, units = "secs")

        }
        
        #
        # saving
        #
        
        result = data.frame(
          n = n_obs, id_rep = id_rep, dgp = dgp, n_fe = n_fe, 
          method = c("fixest::feols", "lfe::felm"),
          time = c(time_feols, time_felm)
        )
        
        all_results_list[[length(all_results_list) + 1]] = result
        
      }
    }
    
  }
  message("|")
}

# saving
all_results = do.call(rbind, all_results_list)

data.table::fwrite(all_results, "./results/results_ols_R.csv")



