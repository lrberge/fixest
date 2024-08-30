
#include "fixest_main.hpp"

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
