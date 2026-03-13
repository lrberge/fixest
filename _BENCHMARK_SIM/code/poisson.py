#------------------------------------------------------------------------------#
# Created: 2026-01-08
# ~: OLS benchmarks for python functions
#------------------------------------------------------------------------------#

# NOTA: on 2026-01-08 pyfixest could not be installed on python 3.14 because 
#       of an incompatible scipy version.
#       One may need to use python <= 3.13 for pyfixest to be installed properly.
#       

# Quadruple loop over:
# 1) number of observations
# 2) number of repetitions
# 3) type of Data Generating Process (dgp)
# 4) number of fixed effects
# 

import time
import pyfixest
import pandas
import math

all_results_list = list()
for n_obs in (1e4, 1e5, 1e6):
  print(f"n = 1e{int(math.log10(n_obs))}")
  
  for id_rep in (1, 2, 3):
    print(f" rep = {id_rep}")
    
    # loading the data 
    filename = f"./data/data_n-1e{int(math.log10(n_obs))}_rep-{id_rep}.gz"
    base = pandas.read_csv(filename, compression='gzip')
    
    # we exponentiate the dep.var
    base['y'] = [math.exp(x) for x in base['y']]
      
    for dgp in ("simple", "difficult"):
      print(f"  dgp = {dgp}", flush = True, end = "")
      
      for n_fe in (2, 3):
        print(".", flush = True, end = "")
        
        # setting the formula depending on the cases
        if dgp == "simple":
          fml = "y ~ x1 + x2 | indiv_id + firm_id" if n_fe == 2 else "y ~ x1 + x2 | indiv_id + firm_id + year"
        elif dgp == "difficult":
          # we don't run for large obs bc run time > 1h
          if n_obs >= 1e6: 
            continue
          fml = "y ~ x1 + x2 | indiv_id + firm_id_difficult" if n_fe == 2 else "y ~ x1 + x2 | indiv_id + firm_id_difficult + year"
        
        
        #
        # pyfixest.fepois
        #
        
        time_start = time.perf_counter()
        est = pyfixest.fepois(fml, base)
        time_fepois = time.perf_counter() - time_start
        
        #
        # saving
        #
        
        result = pandas.DataFrame({
          'n': [n_obs], 'id_rep': [id_rep], 'dgp': [dgp], 'n_fe': [n_fe],
          'method': ["pyfixest.fepois"], 'time': [time_fepois]
        })
        
        all_results_list.append(result)
      
      print("")
    

# saving
all_results = pandas.concat(all_results_list)

all_results.to_csv("./results/results_poisson_py.csv", index = False)







