
#pragma once

#include <Rcpp.h>
#include <vector>
#include <cstdint>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
  #define omp_get_max_threads() 0
#endif


std::vector<int> set_parallel_scheme(int N, int nthreads);

//
// inline 
//


/**
 * @brief checks whether the algorithm should continue
 * 
 * @param a numeric scalar
 * @param b numeric scalar
 * @param diffMax max absolute difference
 * @return true if the absolute difference between a and b is greater than diffMax
 */
inline bool continue_crit(double a, double b, double diffMax){
  // continuing criterion of the algorithm
  double diff = fabs(a - b);
  return ( (diff > diffMax) && (diff/(0.1 + fabs(a)) > diffMax) );
}

/**
 * @brief checks whether the algorithm should stop
 * 
 * @param a numeric scalar
 * @param b numeric scalar
 * @param diffMax max absolute difference
 * @return true if the absolute difference between a and b is lower than diffMax
 */
inline bool stopping_crit(double a, double b, double diffMax){
  // stopping criterion of the algorithm
  double diff = fabs(a - b);
  return ( (diff < diffMax) || (diff/(0.1 + fabs(a)) < diffMax) );
}

//
// classes 
//

//
// We introduce a class that handles varying types of SEXP and behaves as a regular matrix
//


/**
 * @brief treats an int/double vector as a double vector without copy
 * 
 */
class RealVec{
  double *p_dble = nullptr;
  int *p_int = nullptr;

public:
  // several constructors

  // is_int public member
  bool is_int = false;

  RealVec(){};
  RealVec(SEXP);
  RealVec(double *p_x): p_dble(p_x), is_int(false){};
  RealVec(int *p_x): p_int(p_x), is_int(true){};
  RealVec(std::nullptr_t){};


  inline double operator[](int i){
    if(is_int) return static_cast<double>(p_int[i]);
    return p_dble[i];
  }

};

/**
 * @brief treats a matrix/R DF as a matrix of doubles without copy
 * 
 */
class RealMat{

  std::vector<RealVec> p_RealVec;
  int n = 0;
  int K = 0;

  RealMat() = delete;

public:
  RealMat(SEXP x, bool single_obs = false);

  inline int nrow(){return n;};
  inline int ncol(){return K;};

  inline RealVec operator[](int k){
    return p_RealVec[k];
  };
  inline double operator()(int i, int k){
    return p_RealVec[k][i];
  };
};



//
// fork detection
//

// Regarding fork detection => I don't know enough yet
// Is it worth the pain finding out?
// 

// #include <pthread.h>
// static bool fixest_in_fork = false;
//
// // Trick taken from data.table to detect forking
// // Actually I don't know if it works...
// 
// void when_fork() {
//   fixest_in_fork = true;
// }
//
// void after_fork() {
//   fixest_in_fork = false;
// }
//
// // [[Rcpp::export]]
// void cpp_setup_fork_presence() {
//  // Called only once at startup
//  #ifdef _OPENMP
//     pthread_atfork(&when_fork, &after_fork, NULL);
//  #endif
// }
//
// // [[Rcpp::export]]
// bool cpp_is_in_fork(){
//     return fixest_in_fork;
// }


