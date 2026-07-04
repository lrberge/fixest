
#include <stdint.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>


namespace indexthis {
  
using std::vector;

enum {T_INT, T_DBL_INT, T_DBL, T_STR};

inline int double_to_uint32(const double &x){
  uint32_t y[2];
  std::memcpy(y, &x, sizeof(y));
  return y[0] + y[1];
}

inline int power_of_two(double x){
  return std::ceil(std::log2(x + 1));
}

inline uint32_t hash_single(uint32_t value, int shifter){
  return (3141592653U * value >> (32 - shifter));
}

inline uint32_t hash_double(uint32_t v1, uint32_t v2, int shifter){
  return (((3141592653U * v1) ^ (3141592653U * v2)) >> (32 - shifter));
}

inline bool is_equal_dbl(double x, double y){
  return std::isnan(x) ? std::isnan(y) : x == y;
}

inline SEXP to_r_vector(const vector<int> &x){
  const int n = x.size();
  SEXP r_vec = PROTECT(Rf_allocVector(INTSXP, n));
  int *p_int = INTEGER(r_vec);
  std::memcpy(p_int, x.data(), sizeof(int) * n);
  UNPROTECT(1);
  
  return r_vec;
}

//
// IndexInputVector --------------------------------------------------------------------
//


// Class very useful to pass around the data on R vectors 
class IndexInputVector {
  
  bool is_initialized = false;
  
  void reset(){
    // we reset to the default values
    
    is_fast_int = false;
    x_range = 0;
    x_range_bin = 0;
    x_min = 0;
    type = 0;
    any_na = true;
    NA_value = -1;
    px_int = nullptr;
    px_dbl = nullptr;
    px_intptr = nullptr;
  }
  
public:
  IndexInputVector() = default;
  IndexInputVector(const SEXP &x) { initialize(x); }
  IndexInputVector(const vector<int> &x) { initialize(x); }
  
  void initialize(const SEXP &x);
  void initialize(const vector<int> &x);
  
  int size() const { return n; }
  
  // public properties
  int n = 0;
  bool is_fast_int = false;
  int x_range = 0;
  int x_range_bin = 0;
  int x_min = 0;
  int type = 0;
  
  // this is only used in the quick ints algorithm
  // for factors and bool we assume there are NAs since we don't traverse the data
  //  to find the range, contrary to ints or dbl_ints
  bool any_na = true;
  int NA_value = -1;
  
  // pointers: only the valid one ends up non-null
  const int *px_int = (int *) nullptr;
  const double *px_dbl = (double *) nullptr;
  const intptr_t *px_intptr = (intptr_t *) nullptr;
  
};

//
//  IndexedVector --------------------------------------------------------------
//


class IndexedVector {
  
  vector<int> firstobs;
  vector<int> table;
  vector<double> sum;
  
  int *p_index = nullptr;
  int n = 0;
  
  bool is_initialized = false;
  
  void check_init() const {
    if(!is_initialized){
      Rf_error("IndexedVector: trying to use IndexedVector while not initialized!!! Please fix.");
    }
  }
  
  void reset(){
    firstobs.clear();
    table.clear();
    sum.clear();
    n = 0;
    p_index = nullptr;
  }
  
public:
  
  IndexedVector() = default;
  IndexedVector(SEXP x){ initialize(x); }
  IndexedVector(vector<int> &x){ initialize(x); }
  
  void initialize(SEXP x){
    if(TYPEOF(x) != INTSXP){
      Rf_error("The class IndexedVector can only be initialized with INTEGER R objects.");
    }
    
    if(is_initialized){
      reset();
    }
    
    p_index = INTEGER(x);
    n = Rf_length(x);
    is_initialized = true;
  }
  
  void initialize(vector<int> &x){
    
    if(is_initialized){
      reset();
    }
    
    p_index = x.data();
    n = x.size();
    is_initialized = true;
  }
  
  vector<int>& get_firstobs(){
    check_init();
    return firstobs;
  }
  
  const vector<int>& get_firstobs() const {
    check_init();
    return firstobs;
  }
  
  vector<int>& get_table(){
    check_init();
    return table;
  }
  
  const vector<int>& get_table() const {
    check_init();
    return table;
  }
  
  vector<double>& get_sum(){
    check_init();
    return sum;
  }
  
  const vector<double>& get_sum() const {
    check_init();
    return sum;
  }
  
  int* get_p_index() const {
    check_init();
    return p_index;
  }
  
  int size() const { 
    check_init();
    return n;
  }
  
};


//
// function declarations -------------------------------------------------------
//

void to_index_main(const SEXP &x, IndexedVector &output);

void to_index_main(const SEXP &x, IndexedVector &output,
                   const bool do_sum, const double *p_vec_to_sum);

void to_index_main(const IndexInputVector &x, IndexedVector &output);

void to_index_main(const IndexInputVector &x, IndexedVector &output, 
                   const bool do_sum, const double *p_vec_to_sum);

void to_index_main(const std::vector<IndexInputVector> &x, IndexedVector &output);

void to_index_main(const std::vector<IndexInputVector> &x, IndexedVector &output, 
                   const bool do_sum, const double *p_vec_to_sum);

} // namespace indexthis




