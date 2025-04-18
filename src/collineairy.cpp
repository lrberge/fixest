
#include "fixest_main.h"
using namespace Rcpp;

inline int get_sexp_ncol(SEXP &x){
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if(!Rf_isNull(dim)){
    return INTEGER(dim)[1];
  } else if(Rf_inherits(x, "data.frame")){
    return Rf_length(x);
  } else {
    return 1;
  }
}


/*
// use 128bit for scaling? Would avoid issues. 
SEXP cpp_scale(SEXP x){
  // x: numeric matrix
  //    garanteed without NA
  
  int N = Rf_length(x);
  int K = get_sexp_ncol(x);
  
  bool is_int = TYPEOF(x) == INTSXP;
  if(!is_int && TYPEOF(x) == INTSXP){
    Rf_error("Internal error: cpp_scale only accepts numeric matrices.");
  }
  
  // step 1: we check if we can center
  IntegerVector is_constant(K, 1);
  
  for(int k=0 ; k<K ; ++k){
    if(is_int){
      int *px = INTEGER(x);
      int val = px[k * N];
      for(int i=1 ; i<N ; ++i){
        if(px[i] != val){
          is_constant[k] = 0;
          break;
        }
      }
    } else {
      double *px = REAL(x);
      double val = px[k * N];
      for(int i=1 ; i<N ; ++i){
        if(px[i] != val){
          is_constant[k] = 0;
          break;
        }
      }
    }
  }
  
}

*/







