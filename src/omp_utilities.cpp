#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 0
#endif

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
int get_nb_threads(){
	int res = omp_get_max_threads();
	return(res);
}

