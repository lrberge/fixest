# If you are running Julia for the first time, please read Readme.md
# 

# estimation
using CSV
using DataFrames
using GLFixedEffectModels, GLM, Distributions

all_results = DataFrame(n = Int[], id_rep = Int[], dgp = String[], 
                        n_fe = Int[], method = String[], time = Float64[])

for n_obs in [1e4, 1e5, 1e6]
  println("n = 1e$(Int(log10(n_obs)))")
  
  for id_rep in [1, 2, 3]
    println(" rep = $(id_rep)")
    
    # loading the data 
    filename = "./data/data_n-1e$(Int(log10(n_obs)))_rep-$id_rep.gz";
    base = CSV.read(filename, DataFrame);
    base.y = exp.(base.y);
    
    for dgp in ["simple", "difficult"]
      print("  dgp = $(dgp)")
    
      for n_fe in [2, 3]
        print(".")
        
        # setting up the formula
        if dgp == "simple"
          if n_fe == 2
            fml = @formula(y ~ x1 + x2 + fe(indiv_id) + fe(firm_id))
          else
            fml = @formula(y ~ x1 + x2 + fe(indiv_id) + fe(firm_id) + fe(year))
          end
        elseif dgp == "difficult"
          if n_fe == 2
            fml = @formula(y ~ x1 + x2 + fe(indiv_id) + fe(firm_id_difficult))
          else
            fml = @formula(y ~ x1 + x2 + fe(indiv_id) + fe(firm_id_difficult) + fe(year))
          end
        end
        
        #
        # GLFixedEffectModels.nlreg
        #
        
        time_start = time()
        est = nlreg(base, fml, Poisson(), LogLink());
        time_nlreg = time() - time_start
        
        
        #
        # saving
        #
        
        push!(all_results, (n_obs, id_rep, dgp, n_fe, "GLFixedEffectModels.nlreg", time_nlreg))
        
      end
      println("")
    end
  end
end


CSV.write("results/results_poisson_jl.csv", all_results)





