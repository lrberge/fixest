
#include "fixest_main.h"

/**
 * @brief Set the parallel scheme for a 1 dimentional matrix
 * so that parallel computation can be leveraged for 1D matrices
 * 
 * @param N number of obs
 * @param nthreads 
 * @return vector length = nber of threads + 1 giving the start/stop of each thread
 */
std::vector<int> set_parallel_scheme(int N, int nthreads){

  std::vector<int> res(nthreads + 1, 0);
  double N_rest = N;

  for(int i=0 ; i<nthreads ; ++i){
    res[i + 1] = ceil(N_rest / (nthreads - i));
    N_rest -= res[i + 1];
    res[i + 1] += res[i];
  }

  return res;
}


//
// RealVec 
//

RealVec::RealVec(SEXP x){
  if(TYPEOF(x) == REALSXP){
    is_int = false;
    p_dble = REAL(x);
  } else if(TYPEOF(x) == INTSXP){
    is_int = true;
    p_int = INTEGER(x);
  } else {
    Rf_error("The current SEXP type is not supported by the RealVec class.");
  }
}

//
// RealMat 
//

RealMat::RealMat(SEXP x, bool single_obs){
  // NOTA: you only default values to declarations, not implementations

  if(TYPEOF(x) == VECSXP){
    // x can be a list of either vectors or matrices

    int L = Rf_length(x);

    for(int l=0 ; l<L ; ++l){
      SEXP xx = VECTOR_ELT(x, l);
      SEXP dim = Rf_getAttrib(xx, R_DimSymbol);

      int n_tmp = 0, K_tmp = 0;

      if(Rf_length(dim) == 0){
        // vector
        n_tmp = Rf_length(xx);
        K_tmp = 1;
      } else {
        int *pdim = INTEGER(dim);
        n_tmp = pdim[0];
        K_tmp = pdim[1];
      }

      // we set the number of rows at the first iteration
      if(l == 0){
        n = n_tmp;
      } else {
        if(n != n_tmp) Rf_error("When setting up the class RealMat: The number of observations in the list is not coherent across columns.");
      }

      K += K_tmp;

      if(TYPEOF(xx) == REALSXP){
        double *p_x = REAL(xx);
        for(int k=0 ; k<K_tmp ; ++k){
          p_RealVec.push_back(RealVec(p_x));
          if(k + 1 < K_tmp) p_x += n;
        }

      } else if(TYPEOF(xx) == INTSXP){
        int *p_x = INTEGER(xx);
        for(int k=0 ; k<K_tmp ; ++k){
          p_RealVec.push_back(RealVec(p_x));
          if(k + 1 < K_tmp) p_x += n;
        }
      } else {
        Rf_error("The current SEXP type is not supported by the RealMat class.");
      }
    }


  } else {
    // Matrix or vector

    SEXP dim = Rf_getAttrib(x, R_DimSymbol);

    if(Rf_length(dim) == 0){
      // vector
      n = Rf_length(x);
      K = 1;
    } else {
      const int *pdim = INTEGER(dim);
      n = pdim[0];
      K = pdim[1];
    }

    if(!single_obs && (n == 1 && K == 1)){
      // => absence of data
      n = 0;
      K = 0;

    } else if(TYPEOF(x) == REALSXP){
      double *p_x = REAL(x);
      for(int k=0 ; k<K ; ++k){
        p_RealVec.push_back(RealVec(p_x));
        if(k + 1 < K) p_x += n;
      }

    } else if(TYPEOF(x) == INTSXP){
      int *p_x = INTEGER(x);
      for(int k=0 ; k<K ; ++k){
        p_RealVec.push_back(RealVec(p_x));
        if(k + 1 < K) p_x += n;
      }
    } else {
      Rf_error("The current SEXP type is not supported by the RealMat class.");
    }
  }
}
