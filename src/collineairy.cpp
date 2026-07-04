
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






