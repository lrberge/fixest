/**********************************************************************
 * ______________                                                     *
 * || Indexing ||                                                     *
 * --------------                                                     *
 *                                                                    *
 * Author: Laurent R. Berge                                           *
 *                                                                    *
 * Turn any vector of any type into an integer index                  *
 * Original code: package indexthis (lrberge)                         *
 *                                                                    *
 *********************************************************************/


#include "util.h"
#include "to_index.h"

namespace indexthis {
  
using std::vector;


void IndexInputVector::initialize(const SEXP &x){
  
  if(is_initialized){
    reset();
  }
  
  is_initialized = true;
  
  int n = Rf_length(x);
  this->n = n;
  
  bool IS_INT = false;
  
  if(TYPEOF(x) == STRSXP){
    // character
    this->type = T_STR;
    this->px_intptr = (intptr_t *) STRING_PTR_RO(x);
    
  } else if(Rf_isNumeric(x) || Rf_isFactor(x) || TYPEOF(x) == LGLSXP){
      
    if(TYPEOF(x) == REALSXP){
      // we check if the underlying structure is int
      this->px_dbl = REAL(x);
      IS_INT = true;
      double *px = REAL(x);
      double x_min = 0, x_max = 0, x_tmp;
      
      // taking care of NA corner cases
      int i_start = 0;
      while(i_start < n && std::isnan(px[i_start])){
        ++i_start;
      }
      
      bool any_na = i_start > 0;
      if(i_start < n){
        x_min = px[i_start];
        x_max = px[i_start];
        
        for(int i=i_start ; i<n ; ++i){
          x_tmp = px[i];
          
          if(std::isnan(x_tmp)){
            any_na = true;
          } else if(!(x_tmp == (int) x_tmp)){
            IS_INT = false;
            break;
          } else if(x_tmp > x_max){
            x_max = x_tmp;
          } else if(x_tmp < x_min){
            x_min = x_tmp;
          }
        }
      }      
      
      this->any_na = any_na;

      this->x_min = static_cast<int>(x_min);
      // +1 for the NAs
      this->x_range = x_max - x_min + 2;
      
      this->type = IS_INT ? T_DBL_INT : T_DBL;
    } else {
      // logical, factor and integer are all integers
      IS_INT = true;
      this->px_int = INTEGER(x);
      this->type = T_INT;
      
      if(TYPEOF(x) == INTSXP){
        int *px = INTEGER(x);
        int x_min = 0, x_max = 0, x_tmp;
        
        // taking care of NA corner cases
        int i_start = 0;
        while(i_start < n && px[i_start] == NA_INTEGER){
          ++i_start;
        }
        bool any_na = i_start > 0;
        if(i_start < n){
          x_min = px[i_start];
          x_max = px[i_start];
          
          for(int i=i_start ; i<n ; ++i){
            x_tmp = px[i];
            
            if(x_tmp > x_max){
              x_max = x_tmp;
            } else if(x_tmp < x_min){
              // NA integer is the smallest int, defined as -2147483648
              if(x_tmp == NA_INTEGER){
                any_na = true;
              } else {
                x_min = x_tmp;
              }            
            }
          }
        }
        
        this->any_na = any_na;
        this->x_min = x_min;
        // +1 for the NAs
        this->x_range = x_max - x_min + 2;
      } else if(TYPEOF(x) == LGLSXP){
        this->x_min = 0;
        // 0, 1, NA
        this->x_range = 3;
      } else {
        // factor
        SEXP labels = Rf_getAttrib(x, R_LevelsSymbol);
        // factors always start at 1
        this->x_min = 1;
        // we add 1 for the NAs
        this->x_range = Rf_length(labels) + 1;
      }
    }
    
    if(IS_INT){
      // finding out if we're in the easy case
      this->x_range_bin = power_of_two(this->x_range);    
      this->is_fast_int = this->x_range < 100000 || this->x_range <= 2*n;
      this->NA_value = this->x_range - 1;
    }
    
  } else {
    // error
    
    Rcpp::stop("Internal error: Only vectors of type character, numeric, logical or factor can be converted to an index.");
    
  }
}


void IndexInputVector::initialize(const vector<int> &x){
  
  if(is_initialized){
    reset();
  }
  
  is_initialized = true;
  
  this->n = x.size();
  
  this->px_int = x.data();
  this->type = T_INT;
  int x_max = *std::max_element(x.begin(), x.end());
  int x_min = *std::min_element(x.begin(), x.end());
  
  this->x_min = x_min;
  this->x_range = x_max - x_min + 1;
  this->x_range_bin = power_of_two(this->x_range);    
  this->is_fast_int = this->x_range < 100000 || this->x_range <= 2*n;
  
  this->any_na = false;
  this->NA_value = this->x_range - 1;
  
}


SEXP std_string_to_r_string(std::vector<std::string> x){
  
  int n = x.size();
  
  SEXP res = PROTECT(Rf_allocVector(STRSXP, n));
  
  for(int i=0 ; i<n ; ++i){
    SET_STRING_ELT(res, i, Rf_mkCharCE(x[i].c_str(), CE_UTF8));
  }
  
  UNPROTECT(1);
  
  return res;
}

void general_type_to_index_single(const IndexInputVector *x, int *__restrict p_index, int &n_groups,
                                  const bool is_final,
                                  vector<int> &vec_firstobs, vector<int> &vec_table, 
                                  vector<double> &vec_sum,
                                  const bool do_sum, const double *p_vec_to_sum){
  
  const size_t n = x->n;
  
  // we hash the vectors successively to turn them into "sparse" int32
  // we combine values from different vectors with xor
  // we cut the hashed value to fit into 2**(log2(n + 1))
  // we fill the group vector and check for collision all the time (this is costly)
  
  
  // we find out the number of bits (see shifter)
  // we find the first multiple of 2 greater than n
  int shifter = power_of_two(2.0 * n + 1.0);
  if(shifter < 8) shifter = 8;
  size_t larger_n = std::pow(2, shifter);
  
  // hashed_obs_vec:
  // - assume hash(value) leads to ID
  // - then hashed_obs_vec[ID] is the observation id of the first observation with that hash
  // - note that using an array makes the algo twice faster
  int *hashed_obs_vec = new int[larger_n + 1];
  std::fill_n(hashed_obs_vec, larger_n + 1, 0);
  
  const int *px_int = (int *) x->px_int;
  const double *px_dbl = (double *) x->px_dbl;
  const intptr_t *px_intptr = (intptr_t *) x->px_intptr;
  
  const int x_type = x->type;
  
  int g = 0;
  uint32_t id = 0;
  int obs = 0;
  if(x_type == T_STR){
    for(size_t i=0 ; i<n ; ++i){
      
      id = hash_single(px_intptr[i] & 0xffffffff, shifter);
      
      bool does_exist = false;
      while(hashed_obs_vec[id] != 0){
        obs = hashed_obs_vec[id] - 1;
        if(px_intptr[obs] == px_intptr[i]){
          p_index[i] = p_index[obs];
          does_exist = true;
          break;
        } else {
          ++id;
          if(id > larger_n){
            id %= larger_n;
          }
        }
      }

      if(!does_exist){
        // hash never seen => ok
        hashed_obs_vec[id] = i + 1;
        p_index[i] = ++g;
        if(is_final){
          vec_firstobs.push_back(i + 1);
          vec_table.push_back(1);
          if(do_sum){
            vec_sum.push_back(p_vec_to_sum[i]);
          }
        }
      } else if(is_final){
        // NOTA: the index is 1-indexed, we need to 0-index it
        const int group_id = p_index[i] - 1;
        vec_table[group_id]++;
        if(do_sum){
          vec_sum[group_id] += p_vec_to_sum[i];
        }
      }
      
    }
  } else if(x_type == T_INT){
    for(size_t i=0 ; i<n ; ++i){
      id = hash_single(px_int[i], shifter);
      
      bool does_exist = false;
      while(hashed_obs_vec[id] != 0){
        obs = hashed_obs_vec[id] - 1;
        if(px_int[obs] == px_int[i]){
          p_index[i] = p_index[obs];
          does_exist = true;
          break;
        } else {
          ++id;
          if(id > larger_n){
            id %= larger_n;
          }
        }
      }

      if(!does_exist){
        hashed_obs_vec[id] = i + 1;
        p_index[i] = ++g;
        if(is_final){
          vec_firstobs.push_back(i + 1);
          vec_table.push_back(1);
          if(do_sum){
            vec_sum.push_back(p_vec_to_sum[i]);
          }
        }
      } else if(is_final){
        // NOTA: the index is 1-indexed, we need to 0-index it
        const int group_id = p_index[i] - 1;
        vec_table[group_id]++;
        if(do_sum){
          vec_sum[group_id] += p_vec_to_sum[i];
        }
      }
    }
  } else {
    // NOTA: the compiler will take the if() out of the loop
    // not possible if we also have an if() inside the while() loop
    // => explains why I needed to repeat all the for loops
    const bool any_na = x->any_na;
    const int NA_value = x->NA_value;
    for(size_t i=0 ; i<n ; ++i){
      
      if(x_type == T_DBL_INT){
        if(any_na){
          if(std::isnan(px_dbl[i])){
            id = hash_single(NA_value, shifter);
          } else {
            id = hash_single((int) px_dbl[i], shifter);
          }
        } else {
          id = hash_single((int) px_dbl[i], shifter);
        }
      } else {
        id = hash_single(double_to_uint32(px_dbl[i]), shifter);
      }      
      
      bool does_exist = false;
      while(hashed_obs_vec[id] != 0){
        obs = hashed_obs_vec[id] - 1;
        if(is_equal_dbl(px_dbl[obs], px_dbl[i])){
          p_index[i] = p_index[obs];
          does_exist = true;
          break;
        } else {
          ++id;
          if(id > larger_n){
            id %= larger_n;
          }
        }
      }

      if(!does_exist){
        hashed_obs_vec[id] = i + 1;
        p_index[i] = ++g;
        if(is_final){
          vec_firstobs.push_back(i + 1);
          vec_table.push_back(1);
          if(do_sum){
            vec_sum.push_back(p_vec_to_sum[i]);
          }
        }
      } else if(is_final){
        // NOTA: the index is 1-indexed, we need to 0-index it
        const int group_id = p_index[i] - 1;
        vec_table[group_id]++;
        if(do_sum){
          vec_sum[group_id] += p_vec_to_sum[i];
        }
      }
    }
  }
  
  n_groups = g;
  delete[] hashed_obs_vec;
  
}

void general_type_to_index_double(const IndexInputVector *x, int *__restrict p_index_in, 
                                  int *__restrict p_index_out, int &n_groups,
                                  const bool is_final,
                                  vector<int> &vec_firstobs, vector<int> &vec_table,
                                  vector<double> &vec_sum,
                                  const bool do_sum, const double *p_vec_to_sum){
  // Two differences with the *_single version:
  // - when hashing and checking for collision => we use the extra index
  // - we include the possibility of fast ints
  //
  
  // general information
  const size_t n = x->n;
  
  const int *px_int = (int *) x->px_int;
  const double *px_dbl = (double *) x->px_dbl;
  const intptr_t *px_intptr = (intptr_t *) x->px_intptr;
  
  const int x_type = x->type;
  int g = 0;
  
  bool do_fast_int = false;
  if(x->is_fast_int){
    int sum_range_bin = x->x_range_bin + power_of_two(n_groups);
    do_fast_int = sum_range_bin < 17 || sum_range_bin <= power_of_two(5.0 * n);
  }
  
  if(do_fast_int){
    
    int n_groups_bin = power_of_two(n_groups);
    size_t lookup_size = std::pow(2, x->x_range_bin + n_groups_bin + 1);
    int *int_array = new int[lookup_size];
    std::fill_n(int_array, lookup_size, 0);
    
    const bool is_x_int = x_type == T_INT;
    const int x_min = x->x_min;
    const int offset = n_groups_bin;
    int id = 0;
    const bool any_na = x->any_na;
    const int NA_value = x->NA_value;
    for(size_t i=0 ; i<n ; ++i){
      if(is_x_int){
        if(any_na){
          if(px_int[i] == NA_INTEGER){
            id = p_index_in[i] + ((NA_value) << offset);
          } else {
            id = p_index_in[i] + ((px_int[i] - x_min) << offset);
          }
        } else {
          id = p_index_in[i] + ((px_int[i] - x_min) << offset);
        }
      } else {
        if(any_na){
          if(std::isnan(px_dbl[i])){
            id = p_index_in[i] + ((NA_value) << offset);
          } else {
            id = p_index_in[i] + ((static_cast<int>(px_dbl[i]) - x_min) << offset);
          }
        } else {
          id = p_index_in[i] + ((static_cast<int>(px_dbl[i]) - x_min) << offset);
        }
      }      
      
      if(int_array[id] == 0){
        ++g;
        int_array[id] = g;
        p_index_out[i] = g;
        if(is_final){
          vec_firstobs.push_back(i + 1);
          vec_table.push_back(1);
          if(do_sum){
            vec_sum.push_back(p_vec_to_sum[i]);
          }
        }
      } else {
        p_index_out[i] = int_array[id];
        
        if(is_final){
          // NOTA: the index is 1-indexed, we need to 0-index it
          const int group_id = p_index_out[i] - 1;
          vec_table[group_id]++;
          if(do_sum){
            vec_sum[group_id] += p_vec_to_sum[i];
          }
        }
      }
    }
    
    delete[] int_array;
    
  } else {  
    // we hash the vectors successively to turn them into "sparse" int32
    // we combine values from different vectors with xor
    // we cut the hashed value to fit into 2**(log2(n + 1))
    // we fill the group vector and check for collision all the time (this is costly)
    
    // we find out the number of bits (see shifter)
    // we find the first multiple of 2 greater than n
    int shifter = power_of_two(2.0 * n + 1.0);
    if(shifter < 8) shifter = 8;
    size_t larger_n = std::pow(2, shifter);
    
    // hashed_obs_vec:
    // - assume hash(value) leads to ID
    // - then hashed_obs_vec[ID] is the observation id of the first observation with that hash
    // - note that using an array makes the algo twice faster
    int *hashed_obs_vec = new int[larger_n + 1];
    std::fill_n(hashed_obs_vec, larger_n + 1, 0);
    
    uint32_t id = 0;
    int obs = 0;
    if(x_type == T_STR){
      for(size_t i=0 ; i<n ; ++i){
        
        id = hash_double(px_intptr[i] & 0xffffffff, p_index_in[i], shifter);
        
        bool does_exist = false;
        while(hashed_obs_vec[id] != 0){
          obs = hashed_obs_vec[id] - 1;
          if(px_intptr[obs] == px_intptr[i] && p_index_in[obs] == p_index_in[i]){
            p_index_out[i] = p_index_out[obs];
            does_exist = true;
            break;
          } else {
            ++id;
            if(id > larger_n){
              id %= larger_n;
            }
          }
        }

        if(!does_exist){
          // hash never seen => ok
          hashed_obs_vec[id] = i + 1;
          p_index_out[i] = ++g;
          if(is_final){
            vec_firstobs.push_back(i + 1);
            vec_table.push_back(1);
            if(do_sum){
              vec_sum.push_back(p_vec_to_sum[i]);
            }
          }
        } else if(is_final){
          // NOTA: the index is 1-indexed, we need to 0-index it
          const int group_id = p_index_out[i] - 1;
          vec_table[group_id]++;
          if(do_sum){
            vec_sum[group_id] += p_vec_to_sum[i];
          }
        }
        
      }
    } else if(x_type == T_INT){
      for(size_t i=0 ; i<n ; ++i){
        id = hash_double(px_int[i], p_index_in[i], shifter);
        
        bool does_exist = false;
        while(hashed_obs_vec[id] != 0){
          obs = hashed_obs_vec[id] - 1;
          if(px_int[obs] == px_int[i] && p_index_in[obs] == p_index_in[i]){
            p_index_out[i] = p_index_out[obs];
            does_exist = true;
            break;
          } else {
            ++id;
            if(id > larger_n){
              id %= larger_n;
            }
          }
        }

        if(!does_exist){
          hashed_obs_vec[id] = i + 1;
          p_index_out[i] = ++g;
          if(is_final){
            vec_firstobs.push_back(i + 1);
            vec_table.push_back(1);
            if(do_sum){
              vec_sum.push_back(p_vec_to_sum[i]);
            }
          }
        } else if(is_final){
          // NOTA: the index is 1-indexed, we need to 0-index it
          const int group_id = p_index_out[i] - 1;
          vec_table[group_id]++;
          if(do_sum){
            vec_sum[group_id] += p_vec_to_sum[i];
          }
        }
      }
    } else {
      // NOTA: the compiler will take the if() out of the loop
      // not possible if we also have an if() inside the while() loop
      // => explains why I needed to repeat all the for loops
      const bool any_na = x->any_na;
      const int NA_value = x->NA_value;
      for(size_t i=0 ; i<n ; ++i){
        
        if(x_type == T_DBL_INT){
          if(any_na){
            if(std::isnan(px_dbl[i])){
              id = hash_double(NA_value, p_index_in[i], shifter);
            } else {
              id = hash_double((int) px_dbl[i], p_index_in[i], shifter);
            }
          } else {
            id = hash_double((int) px_dbl[i], p_index_in[i], shifter);
          }
        } else {
          id = hash_double(double_to_uint32(px_dbl[i]), p_index_in[i], shifter);
        }
        
        bool does_exist = false;
        while(hashed_obs_vec[id] != 0){
          obs = hashed_obs_vec[id] - 1;
          if(is_equal_dbl(px_dbl[obs], px_dbl[i]) && p_index_in[obs] == p_index_in[i]){
            p_index_out[i] = p_index_out[obs];
            does_exist = true;
            break;
          } else {
            ++id;
            if(id > larger_n){
              id %= larger_n;
            }
          }
        }

        if(!does_exist){
          hashed_obs_vec[id] = i + 1;
          p_index_out[i] = ++g;
          if(is_final){
            vec_firstobs.push_back(i + 1);
            vec_table.push_back(1);
            if(do_sum){
              vec_sum.push_back(p_vec_to_sum[i]);
            }
          }
        } else if(is_final){
          // NOTA: the index is 1-indexed, we need to 0-index it
          const int group_id = p_index_out[i] - 1;
          vec_table[group_id]++;
          if(do_sum){
            vec_sum[group_id] += p_vec_to_sum[i];
          }
        }
      }
    }
    
    delete[] hashed_obs_vec;
  }
  
  n_groups = g;
}


inline void update_index_intarray_g_obs(int id, size_t i, int &g, int * &int_array, 
                                        int *__restrict &p_index, const bool &is_final, 
                                        vector<int> &vec_firstobs, vector<int> &vec_table,
                                        vector<double> &vec_sum, const bool do_sum, 
                                        const double *p_vec_to_sum){
  
  if(int_array[id] == 0){
    ++g;
    int_array[id] = g;
    p_index[i] = g;
    if(is_final){
      vec_firstobs.push_back(i + 1);
      vec_table.push_back(1);
      if(do_sum){
        vec_sum.push_back(p_vec_to_sum[i]);
      }
    }
  } else {
    p_index[i] = int_array[id];
    if(is_final){
      // NOTA: the index is 1-indexed, we need to 0-index it
      const int group_id = p_index[i] - 1;
      vec_table[group_id]++;
      if(do_sum){
        vec_sum[group_id] += p_vec_to_sum[i];
      }
    }
  }  
}

void multiple_ints_to_index(const vector<IndexInputVector> &all_vecs, vector<int> &all_k, 
                            int *__restrict p_index, int &n_groups, const bool is_final,
                            vector<int> &vec_firstobs, vector<int> &vec_table,
                            vector<double> &vec_sum,
                            const bool do_sum, const double *p_vec_to_sum){
  
  int sum_bin_ranges = 0;
  int K = all_k.size();
    
  for(auto &&k : all_k){
    sum_bin_ranges += all_vecs[k].x_range_bin;
  }  
  
  int k0 = all_k[0];
  const IndexInputVector *x0 = &all_vecs[k0];
  const size_t n = x0->n;
  const int * px0_int = (int *) x0->px_int;
  const double * px0_dbl = (double *) x0->px_dbl;
  
  const int x0_type = x0->type;
  const bool is_x0_int = x0_type == T_INT;
  const int x0_min = x0->x_min;
  
  size_t lookup_size = K == 1 ? x0->x_range + 1 : std::pow(2, sum_bin_ranges + K - 1);
  int *int_array = new int[lookup_size];
  std::fill_n(int_array, lookup_size, 0);
  
  int g = 0;                    
  
  if(K == 1){
    int id = 0;
    
    if(is_x0_int){
      if(x0->any_na){
        int NA_value = x0->NA_value;
        for(size_t i=0 ; i<n ; ++i){
          if(px0_int[i] == NA_INTEGER){
            id = NA_value;
          } else {
            id = px0_int[i] - x0_min;
          }
          
          update_index_intarray_g_obs(id, i, g, int_array, p_index, is_final, 
                                      vec_firstobs, vec_table, vec_sum, do_sum, p_vec_to_sum);
        }
      } else {
        for(size_t i=0 ; i<n ; ++i){
          
          id = px0_int[i] - x0_min;
          update_index_intarray_g_obs(id, i, g, int_array, p_index, is_final, 
                                      vec_firstobs, vec_table, vec_sum, do_sum, p_vec_to_sum);
        }
      }
    } else {
      // Double as int
      if(x0->any_na){
        int NA_value = x0->NA_value;
        for(size_t i=0 ; i<n ; ++i){
          if(std::isnan(px0_dbl[i])){
            id = NA_value;
          } else {
            id = static_cast<int>(px0_dbl[i]) - x0_min;
          }
          
          update_index_intarray_g_obs(id, i, g, int_array, p_index, is_final, 
                                      vec_firstobs, vec_table, vec_sum, do_sum, p_vec_to_sum);
        }
      } else {
        for(size_t i=0 ; i<n ; ++i){
          id = static_cast<int>(px0_dbl[i]) - x0_min;
          update_index_intarray_g_obs(id, i, g, int_array, p_index, is_final, 
                                      vec_firstobs, vec_table, vec_sum, do_sum, p_vec_to_sum);
        }
      }
    }
  } else {
    
    int k1 = all_k[1];
    const IndexInputVector *x1 = &all_vecs[k1];
    const int *px1_int = (int *) x1->px_int;
    const double *px1_dbl = (double *) x1->px_dbl;
    
    const bool x1_type = x1->type;
    const bool is_x1_int = x1_type == T_INT;   
    const int x1_min = x1->x_min;
    
    // ad hoc macro to update v depending on the vector type
    #define UPDATE_V(is_xk_int, any_na_xk, pxk_int, pxk_dbl, vk, NA_value_xk, xk_min)  \
      if(is_xk_int){                                                                   \
        if(any_na_xk){                                                                 \
          if(pxk_int[i] == NA_INTEGER){                                                \
            vk = NA_value_xk;                                                          \
          } else {                                                                     \
            vk = pxk_int[i] - xk_min;                                                  \
          }                                                                            \
        } else {                                                                       \
          vk = pxk_int[i] - xk_min;                                                    \
        }                                                                              \
      } else {                                                                         \
        if(any_na_xk){                                                                 \
          if(std::isnan(pxk_dbl[i])){                                                  \
            vk = NA_value_xk;                                                          \
          } else {                                                                     \
            vk = static_cast<int>(pxk_dbl[i]) - xk_min;                                \
          }                                                                            \
        } else {                                                                       \
          vk = static_cast<int>(pxk_dbl[i]) - xk_min;                                  \
        }                                                                              \
      }
    
    if(K == 2){
      int id = 0;
      int offset = x0->x_range_bin;
      int v0 = 0, v1 = 0;
      const bool any_na_x0 = x0->any_na, any_na_x1 = x1->any_na;
      const int NA_value_x0 = x0->NA_value, NA_value_x1 = x1->NA_value;
      for(size_t i=0 ; i<n ; ++i){
        
        UPDATE_V(is_x0_int, any_na_x0, px0_int, px0_dbl, v0, NA_value_x0, x0_min)
        UPDATE_V(is_x1_int, any_na_x1, px1_int, px1_dbl, v1, NA_value_x1, x1_min)
        
        id = v0 + (v1 << offset);
        
        if(int_array[id] == 0){
          ++g;
          int_array[id] = g;
          p_index[i] = g;
          if(is_final){
            vec_firstobs.push_back(i + 1);
            vec_table.push_back(1);
            if(do_sum){
              vec_sum.push_back(p_vec_to_sum[i]);
            }
          }
        } else {
          p_index[i] = int_array[id];
          if(is_final){
            // NOTA: the index is 1-indexed, we need to 0-index it
            const int group_id = p_index[i] - 1;
            vec_table[group_id]++;
            if(do_sum){
              vec_sum[group_id] += p_vec_to_sum[i];
            }
          }
        }
      }
    } else {
      int *sum_vec = new int[n];
      
      // we first initialize sum_vec to the sum of the two first elements
      
      int offset = x0->x_range_bin;
      int v0 = 0, v1 = 0;
      const bool any_na_x0 = x0->any_na, any_na_x1 = x1->any_na;
      const int NA_value_x0 = x0->NA_value, NA_value_x1 = x1->NA_value;
      // the `if`s below in the macro are *not* unrolled by the compiler, sigh
      // unrolling them by hand is too costly and will create unreadable code (16 differents cases!!!!, it would be over 300 lines of code!). 
      // But, on the benchmarks, we could gain 15% perf. So be it.
      for(size_t i=0 ; i<n ; ++i){
        UPDATE_V(is_x0_int, any_na_x0, px0_int, px0_dbl, v0, NA_value_x0, x0_min)
        UPDATE_V(is_x1_int, any_na_x1, px1_int, px1_dbl, v1, NA_value_x1, x1_min)
        
        sum_vec[i] = v0 + (v1 << offset);
      }
      
      // we sum
      offset += x1->x_range_bin;
      for(int ind=2 ; ind<K-1 ; ++ind){
        int k = all_k[ind];
        const IndexInputVector *xk = &all_vecs[k];
        const int *pxk_int = (int *) xk->px_int;
        const double *pxk_dbl = (double *) xk->px_dbl;
        
        const int x_type = xk->type;
        const bool is_xk_int = x_type == T_INT;
        const int xk_min = xk->x_min;
        int vk = 0;
        const bool any_na_xk = xk->any_na;
        const int NA_value_xk = xk->NA_value;
        for(size_t i=0 ; i<n ; ++i){
          UPDATE_V(is_xk_int, any_na_xk, pxk_int, pxk_dbl, vk, NA_value_xk, xk_min)
          
          sum_vec[i] += (vk << offset);
        }
        offset += xk->x_range_bin;
      }
            
      // last element + group creation
      int k = all_k[K - 1];
      const IndexInputVector *xk = &all_vecs[k];
      const int *pxk_int = (int *) xk->px_int;
      const double *pxk_dbl = (double *) xk->px_dbl;
      
      const int x_type = xk->type;
      const bool is_xk_int = x_type == T_INT;
      const int xk_min = xk->x_min;
      
      int id = 0;
      int vk = 0;
      const bool any_na_xk = xk->any_na;
      const int NA_value_xk = xk->NA_value;
      for(size_t i=0 ; i<n ; ++i){
        
        UPDATE_V(is_xk_int, any_na_xk, pxk_int, pxk_dbl, vk, NA_value_xk, xk_min)
        
        id = sum_vec[i] + (vk << offset);
        
        if(int_array[id] == 0){
          // new item: we save the group id
          ++g;
          int_array[id] = g;
          p_index[i] = g;
          if(is_final){
            vec_firstobs.push_back(i + 1);
            vec_table.push_back(1);
            if(do_sum){
              vec_sum.push_back(p_vec_to_sum[i]);
            }
          }
        } else {
          p_index[i] = int_array[id];
          if(is_final){
            // NOTA: the index is 1-indexed, we need to 0-index it
            const int group_id = p_index[i] - 1;
            vec_table[group_id]++;
            if(do_sum){
              vec_sum[group_id] += p_vec_to_sum[i];
            }
          }
        }
      }
      
      delete[] sum_vec;
    }
  }
  
  n_groups = g;
  
  delete[] int_array;
}

std::vector<IndexInputVector> SEXP_to_vec_r_vector(const SEXP &x){
  
  size_t n = 0;
  int K = 0;
  std::vector<IndexInputVector> all_vecs;
  
  // we set up the info with the rvec class. It makes it easy to pass across functions
  if(TYPEOF(x) == VECSXP){
    K = Rf_length(x);
    for(int k=0; k<K; ++k){
      IndexInputVector rvec{VECTOR_ELT(x, k)};
      all_vecs.push_back(rvec);
      
      if(k == 0){
        n = Rf_length(VECTOR_ELT(x, 0));
      } else if((size_t) Rf_length(VECTOR_ELT(x, k)) != n){
        Rcpp::stop("All the vectors to turn into an index must be of the same length. This is currently not the case.");
      }
    }
    
  } else {
    IndexInputVector rvec{x};
    all_vecs.push_back(rvec);
  }
  
  return all_vecs;
}

void to_index_main(const SEXP &x, IndexedVector &output, 
                   const bool do_sum, const double *p_vec_to_sum){
  
  std::vector<IndexInputVector> all_vecs = SEXP_to_vec_r_vector(x);
  
  to_index_main(all_vecs, output, do_sum, p_vec_to_sum);
}

void to_index_main(const SEXP &x, IndexedVector &output){
  const bool do_sum = false;
  const double *p_vec_to_sum = nullptr;
  to_index_main(x, output, do_sum, p_vec_to_sum);
}

void to_index_main(const IndexInputVector &x, IndexedVector &output, 
                   const bool do_sum, const double *p_vec_to_sum){
  std::vector<IndexInputVector> all_vecs;
  all_vecs.push_back(x);
  to_index_main(all_vecs, output, do_sum, p_vec_to_sum);
}

void to_index_main(const IndexInputVector &x, IndexedVector &output){
  const bool do_sum = false;
  const double *p_vec_to_sum = nullptr;
  to_index_main(x, output, do_sum, p_vec_to_sum);
}

void to_index_main(const std::vector<IndexInputVector> &x, IndexedVector &output){
  const bool do_sum = false;
  const double *p_vec_to_sum = nullptr;
  to_index_main(x, output, do_sum, p_vec_to_sum);
}

void to_index_main(const std::vector<IndexInputVector> &all_vecs, IndexedVector &output, 
                   const bool do_sum, const double *p_vec_to_sum){
  
  // x: vector or list of vectors of the same length (n)
  // returns:
  // - index: vector of length n, from 1 to the numberof unique values of x (g)
  // - first_obs: vector of length g of the first observation belonging to each group
  
  int K = all_vecs.size();
  int n = all_vecs.at(0).size();
  
  if(n != output.size()){
    Rcpp::stop("Internal error `to_index_main`: The index size allocated in output is different from the input!");
  }
  
  
  // the result to be returned
  int *p_index = output.get_p_index();
  
  // vector of the first observation of the group
  std::vector<int> &vec_firstobs = output.get_firstobs();
  std::vector<int> &vec_table = output.get_table();
  std::vector<double> &vec_sum = output.get_sum();
  
  // finding out the fast cases
  // Note that partial fast ordering is enabled and 
  // we stop at the first feasible possibility
  int sum_bin_ranges = 0;
  vector<int> id_fast_int;
  for(int k=0 ; k<K ; ++k){
    const IndexInputVector *x = &all_vecs[k];
    if(x->is_fast_int){
      int new_bin_range = sum_bin_ranges + x->x_range_bin;
      if(new_bin_range < 17 || (K >= 2 && new_bin_range <= power_of_two(5 * n))){
        id_fast_int.push_back(k);
        sum_bin_ranges = new_bin_range;
      } else {
        break;
      }      
    }
  }
  
  int n_groups;
  
  //
  // STEP 1: taking care of fast indexing of ints
  //
  
  bool is_final = false;
  bool init_done = false;
  if(!id_fast_int.empty()){
    init_done = true;
    
    is_final = (size_t) K == id_fast_int.size();
    multiple_ints_to_index(all_vecs, id_fast_int, p_index, n_groups, is_final, 
                           vec_firstobs, vec_table, vec_sum, do_sum, p_vec_to_sum);
  }
  
  if(!is_final){
    
    // 
    // STEP 2: general algorithm
    //
    
    // note: here only if not all vectors are "fast"
    //  
    
    // first we find out who is left
    vector<int> all_k_left;
    for(int k=0 ; k<K ; ++k){
      if(std::find(id_fast_int.begin(), id_fast_int.end(), k) == id_fast_int.end()){
        // i.e. if k not in id_fast_int (which has been done)
        all_k_left.push_back(k);
      }
    }
    
    if(!init_done){
      int k0 = all_k_left[0];
      // we remove that element
      all_k_left.erase(all_k_left.begin());
      
      is_final = all_k_left.empty();
      general_type_to_index_single(&all_vecs[k0], p_index, n_groups, is_final, 
                                   vec_firstobs, vec_table, vec_sum, do_sum, p_vec_to_sum);
    }
    
    if(!is_final){
      // here: p_index is an index
      // we will loop over all remainining items and create the index sequentially
      
      int *p_extra_index = new int[n];
      
      bool is_res_updated_index = true;
      for(size_t ind=0 ; ind<all_k_left.size() ; ++ind){
        int k = all_k_left[ind];
        is_final = ind == all_k_left.size() - 1;
        if(is_res_updated_index){
          general_type_to_index_double(&all_vecs[k], p_index, p_extra_index, n_groups, 
                                       is_final, vec_firstobs, vec_table, vec_sum, 
                                       do_sum, p_vec_to_sum);
          is_res_updated_index = false;
        } else {
          general_type_to_index_double(&all_vecs[k], p_extra_index, p_index, n_groups, 
                                       is_final, vec_firstobs, vec_table, vec_sum, 
                                       do_sum, p_vec_to_sum);
          is_res_updated_index = true;
        }
      }
      
      if(!is_res_updated_index){
        std::memcpy(p_index, p_extra_index, sizeof(int) * n);
      }
      delete[] p_extra_index;
    }
  }
}
 
} // end namespace indexthis
 
 
 // [[Rcpp::export]]
SEXP cpp_to_index(SEXP &x){
  
  // x can be a vector or a list of vectors of the same length
  std::vector<indexthis::IndexInputVector> all_vecs = indexthis::SEXP_to_vec_r_vector(x);
  
  const int n = all_vecs.at(0).size();
  SEXP index = PROTECT(Rf_allocVector(INTSXP, n));
  
  indexthis::IndexedVector index_info(index);
  
  //
  // computing the index 
  //
  
  to_index_main(all_vecs, index_info);
  
  //
  // saving the results 
  //
  
  // we copy the first observations into an R vector
  SEXP first_obs = PROTECT(indexthis::to_r_vector(index_info.get_firstobs()));
  SEXP table = PROTECT(indexthis::to_r_vector(index_info.get_table()));
  
  // we save the results into a list
  SEXP res = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, index);
  SET_VECTOR_ELT(res, 1, first_obs);
  SET_VECTOR_ELT(res, 2, table);
  
  // names
  Rf_setAttrib(res, R_NamesSymbol, indexthis::std_string_to_r_string({"index", "first_obs", "table"}));
    
  UNPROTECT(4);
  
  return res;
  
}



