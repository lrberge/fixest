/***********************************************************************
 * ____________________________                                        *
 * || Fixed-effects Indexing ||                                        *
 * ----------------------------                                        *
 *                                                                     *
 * Author: Laurent R. Berge                                            *
 *                                                                     *
 * Turns the fixed-effects indicators into integers, and drops the FEs *
 * that need to be dropped because of singleton or perfect fit         *
 *                                                                     *
 **********************************************************************/

#include "Rcpp.h"
#include "to_index.h"
#include "util.h"

inline std::vector<int> seq(int from, int to){
  const int n = to - from + 1;
  std::vector<int> res(n);
  for(int i = 0 ; i < n ; ++i){
    res[i] = i + from;
  }
  
  return res;
}

void mark_obs_to_remove(std::vector<char> &removed_flag, bool &any_removed,
                        std::vector<int> &firstobs_rm,
                        const indexthis::IndexedVector &index_info, 
                        const bool rm_0, const bool rm_1, const bool rm_single){
  
  const std::vector<int> &firstobs = index_info.get_firstobs();
  const std::vector<int> &table = index_info.get_table();
  const std::vector<double> &sum_y = index_info.get_sum();
  
  //
  // step 1: we check if at least one is removed 
  //
  
  // G: number of groups
  const int G = table.size();
  bool any_to_remove = false;
  int i_start = 0;
  for(int g = 0 ; g < G ; ++g){
    if(rm_single && table[g] == 1){
      any_to_remove = true;
      i_start = firstobs[g] - 1;
      break;
    } else if(rm_0 && sum_y[g] == 0){
      any_to_remove = true;
      i_start = firstobs[g] - 1;
      break;
    } else if(rm_1 && sum_y[g] == table[g]){
      any_to_remove = true;
      i_start = firstobs[g] - 1;
      break;
    }
  }
  
  if(!any_to_remove){
    return;
  }
  
  //
  // step 2: we mark each observation 
  //
  
  // we will remove at least one
  any_removed = true;
  
  int *p_index = index_info.get_p_index();
  const int n = index_info.size();
  std::vector<bool> is_done(G, false);
  for(int i = i_start ; i < n ; ++i){
    const int g = p_index[i] - 1;
    if((rm_single && table[g] == 1) || (rm_0 && sum_y[g] == 0) || (rm_1 && sum_y[g] == table[g])){
      removed_flag[i] = 1;
      if(!is_done[g]){
        is_done[g] = true;
        firstobs_rm.push_back(firstobs[g]);
      }
    }
  }
  
}


// [[Rcpp::export]]
SEXP cpp_index_table_sum(SEXP fixef_list, SEXP y, const bool save_sum_y, 
                         const bool rm_0, const bool rm_1, const bool rm_single, 
                         Rcpp::IntegerVector only_slope, const int nthreads){
  
  int Q = Rf_length(fixef_list);
  SEXP fixef_vec = VECTOR_ELT(fixef_list, 0);
  const int n_obs = Rf_length(fixef_vec);
  
  // do_removal => whether we check for removal
  std::vector<bool> do_removal(Q, true);
  bool any_to_check_for_removal = rm_0 || rm_1 || rm_single;
  
  if(any_to_check_for_removal && only_slope.length() == Q){
    bool any_do_removal = false;
    for(int q=0 ; q<Q ; ++q){
      // we check for pblm only if NOT only slope
      do_removal[q] = only_slope[q] == false;
      any_do_removal = any_do_removal || do_removal[q];
    }
    
    any_to_check_for_removal = any_do_removal;
  }
  
  if(!any_to_check_for_removal){
    do_removal = std::vector<bool>(Q, false);
  }
  
  
  double *p_y = nullptr;
  if(TYPEOF(y) == REALSXP){
    p_y = REAL(y);
  } else {
    // => there will be no use of y, so nullptr is OK
    // but I must ensure that beforehand: save_sum_y = rm_0 = rm_1 = false
    if(save_sum_y || rm_0 || rm_1){
      Rcpp::stop("Internal error, cpp_index_table_sum: you need to give y to apply save_sum_y/rm_0/rm_1.");
    }
  }
  
  const bool do_sum_y = save_sum_y || rm_0 || rm_1;
  
  if(do_sum_y && Rf_length(y) != n_obs){
    Rcpp::stop("Internal error, cpp_index_table_sum: the length of y is different from the length of the fixed-effects.");
  }
  
  // we intialize the information on the indexes to be computed
  
  // we pre allocate the index output vector (the largest bit)
  Rcpp::List all_indexes_sexp(Q);
  for(int q = 0 ; q < Q ; ++q){
    all_indexes_sexp[q] = PROTECT(Rf_allocVector(INTSXP, n_obs));
  }
  
  std::vector<indexthis::IndexInputVector> all_input_vectors(Q);
  std::vector<indexthis::IndexedVector> all_index_info(Q);
  for(int q = 0 ; q < Q ; ++q){
    const SEXP &fixef_vec = VECTOR_ELT(fixef_list, q);
    all_input_vectors[q].initialize(fixef_vec);
    
    all_index_info[q].initialize(all_indexes_sexp[q]);
  }
  
  // below: only used when some observations were removed
  std::vector< std::vector<int> > all_raw_input_vectors(Q);
  std::vector<double> y_new;
  
  std::vector<char> removed_flag;
  std::vector<int> obs_keep;
  std::vector<int> obs_removed;
  std::vector< std::vector<int> > all_firstobs_rm(Q);
  
  bool do_sort_obs_removed = false;
  bool keep_running = true;
  bool first_iter = true;
  while(keep_running){
    keep_running = false;
    
    const int n_current = all_input_vectors[0].size();
    bool any_removed = false;
    std::vector< std::vector<int> > all_firstobs_rm_new(Q);
    
    if(any_to_check_for_removal){
      removed_flag = std::vector<char>(n_current, 0);
    }
    
    #pragma omp parallel for num_threads(nthreads)
    for(int q = 0 ; q < Q ; ++q){
      
      const indexthis::IndexInputVector &fixef_vec = all_input_vectors[q];
      indexthis::IndexedVector &index_info = all_index_info[q];
      
      indexthis::to_index_main(fixef_vec, index_info, do_sum_y, p_y);
      
      if(do_removal[q]){
        mark_obs_to_remove(removed_flag, any_removed, all_firstobs_rm_new[q],
                           index_info, rm_0, rm_1, rm_single);
      }
    }
    
    if(any_removed){
      keep_running = true;
      
      if(first_iter){
        first_iter = false;
        // This is the initialization: first iteration of the while loop
        for(int i = 0 ; i < n_current ; ++i){
          if(removed_flag[i] == 1){
            obs_removed.push_back(i + 1);
          } else {
            obs_keep.push_back(i + 1);
            if(do_sum_y){
              y_new.push_back(p_y[i]);
            }
          }
        }
        
        // firstobs removed
        for(int q = 0 ; q < Q ; ++q){
          all_firstobs_rm[q] = std::move(all_firstobs_rm_new[q]);
        }
        
      } else {
        // here we use the information on the obs_keep bc we nee to track 
        // which observation was removed
        
        // firstobs removed: we need to do it before obs_keep is modified
        for(int q = 0 ; q < Q ; ++q){
          auto &firstobs_rm = all_firstobs_rm[q];
          for(const auto &obs_rm : all_firstobs_rm_new[q]){
            firstobs_rm.push_back(obs_keep[obs_rm - 1]);
          }
        }
        
        // since we erase, we need to start from the end
        do_sort_obs_removed = true;
        for(int i = n_current - 1 ; i >= 0 ; --i){
          if(removed_flag[i] == 1){
            obs_removed.push_back(obs_keep[i]);
            obs_keep.erase(obs_keep.begin() + i);
            if(do_sum_y){
              y_new.erase(y_new.begin() + i);
            }
          }
        }
        
      }
      
      p_y = y_new.data();
      
      if(obs_keep.empty()){
        UNPROTECT(Q);
        Rcpp::List res = Rcpp::List::create(Rcpp::Named("all_removed") = true);
        return res;
      }
      
      // we take the index just computed (in index_info) as the new input
      const int n_new = obs_keep.size();
      #pragma omp parallel for num_threads(nthreads)
      for(int q = 0 ; q < Q ; ++q){
        
        // 1) we use the computed index as new input
        int *p_index = all_index_info[q].get_p_index();
        std::vector<int> new_input(n_new);
        int index = 0;
        for(int i = 0 ; i < n_current ; ++i){
          if(removed_flag[i] == 0){
            new_input[index++] = p_index[i];
          }
        }
        
        all_raw_input_vectors[q] = std::move(new_input);
        all_input_vectors[q].initialize(all_raw_input_vectors[q]);
        
      }
      
      // NOTA: we cannot call R API functions in a multi threaded loop
      for(int q = 0 ; q < Q ; ++q){
        // 2) we reset the containers of the future indexes
        UNPROTECT(1);
        all_indexes_sexp[q] = PROTECT(Rf_allocVector(INTSXP, n_new));
        all_index_info[q].initialize(all_indexes_sexp[q]);
      }
      
    }
    
  }
  
  //
  // building the return object 
  //
  
  /* RETURNED OBJECT
  * 
  * - index: 1-G values of length n_new (over the number of observations, 
  *   once we remove the ones to be removed)
  * - firstobs: first occurrence of the FEs that are kept
  * - table: the number of items for the FEs that are kept
  * - sum_y: the sum of y for the FEs that are kept
  * - obs_removed: the ID of the observations that were removed
  * - firstobs_removed: first occurence of the FEs that were removed
  * 
  * */
  
  Rcpp::List res;
  
  res["index"] = all_indexes_sexp;
  
  // table
  Rcpp::List all_tables(Q);
  Rcpp::List all_sum_y(Q);
  Rcpp::List all_firstobs(Q);
  Rcpp::List all_firstobs_removed(Q);
  const bool any_removed = !obs_removed.empty();
  for(int q = 0 ; q < Q ; ++q){
    
    all_tables[q] = all_index_info[q].get_table();
    
    if(save_sum_y){
      all_sum_y[q] = all_index_info[q].get_sum();
    } else {
      all_sum_y[q] = 0;
    }
    
    if(any_removed){
      auto firstobs_new = all_index_info[q].get_firstobs();
      for(auto &firstobs : firstobs_new){
        firstobs = obs_keep[firstobs - 1];
      }
      
      all_firstobs[q] = firstobs_new;
    } else {
      all_firstobs[q] = all_index_info[q].get_firstobs();
    }
    
    if(any_removed){
      all_firstobs_removed[q] = all_firstobs_rm[q];
    }
  }
  
  res["table"] = all_tables;
  res["sum_y"] = all_sum_y;
  res["firstobs"] = all_firstobs;
  if(any_removed){
    res["firstobs_removed"] = all_firstobs_removed;
  }
  
  if(any_removed){
    if(do_sort_obs_removed){
      std::sort(obs_removed.begin(), obs_removed.end());
    }
    res["obs_removed"] = obs_removed;
  }
  
  UNPROTECT(Q);
  
  return res;
  
}





