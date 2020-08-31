/***************************************************************************
 * _________________                                                       *
 * || Convergence ||                                                       *
 * -----------------                                                       *
 *                                                                         *
 * Author: Laurent R. Berge                                                *
 *                                                                         *
 * Compute the optimal set of cluster coefficients based on Berge (2018),  *
 * in a ML framework.                                                      *
 *                                                                         *
 * The code concerns both the cluster coefficients (i.e. fixed-effects)    *
 * and the coefficients of the derivatives of the cluster coefficients.    *
 *                                                                         *
 * 2-FEs trick:                                                            *
 * For both the Poisson and the Gaussian methods I use a trick to hasten   *
 * computation. Instead of looping on the number of obsertvations, the     *
 * loop is on the number of unique cases. For balanced panels, this is     *
 * not useful, but for strongly unbalanced panels, this makes a big        *
 * difference.                                                             *
 *                                                                         *
 * Logit/Negbin:                                                           *
 * To obtain the cluster coefficients for these two likelihoods, I use     *
 * a Newton-Raphson + dichotomy algorithm. This ensures convergence        *
 * even when the NR algo would go astray (which may be the case in some    *
 * situations).                                                            *
 *                                                                         *
 * Drawback:                                                               *
 * The big drawback of the ML methods is that, as opposed to GLM metohods, *
 * parallel computing cannot be leveraged.                                 *
 *                                                                         *
 *                                                                         *
 **************************************************************************/

#include <Rcpp.h>
#include <math.h>
#include <vector>
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using std::vector;

// Stopping / continuing criteria
// Functions used inside all loops
inline bool continue_criterion(double a, double b, double diffMax){
    // continuing criterion of the algorithm
    double diff = fabs(a - b);
    return ( (diff > diffMax) && (diff/(0.1 + fabs(a)) > diffMax) );
}

inline bool stopping_criterion(double a, double b, double diffMax){
    // stopping criterion of the algorithm
    double diff = fabs(a - b);
    return ( (diff < diffMax) || (diff/(0.1 + fabs(a)) < diffMax) );
}

// List of objects that will be used to
// lighten the writting of the functions
// Contains only parameters that are fixed throguhout ALL the process
struct PARAM_CCC{
	int family;
	int n_obs;
	int K;
	double theta;
	double diffMax_NR;
	int nthreads;

	// vectors from R
	double *mu_init;
	int *pcluster;
	double *lhs;

	// vectors of pointers
	vector<int*> pdum;
	vector<int*> ptable;
	vector<double*> psum_y;
	vector<int*> pobsCluster;
	vector<int*> pcumtable;

	// value that will vary
	double *mu_with_coef;

};

// IT update + returns numerical convergence indicator
bool update_X_IronsTuck(int nb_coef_no_K, vector<double> &X,
                        const vector<double> &GX, const vector<double> &GGX,
                        vector<double> &delta_GX, vector<double> &delta2_X){

	for(int i=0 ; i<nb_coef_no_K ; ++i){
	    double GX_tmp = GX[i];
	    delta_GX[i] = GGX[i] - GX_tmp;
	    delta2_X[i] = delta_GX[i] - GX_tmp + X[i];
		// delta_GX[i] = GGX[i] - GX[i];
		// delta2_X[i] = delta_GX[i] - GX[i] + X[i];
	}

	// delta_GX %*% delta2_X and crossprod(delta2_X)
	double vprod = 0, ssq = 0;
	for(int i=0 ; i<nb_coef_no_K ; ++i){
	    double delta2_X_tmp = delta2_X[i];
	    vprod += delta_GX[i] * delta2_X_tmp;
	    ssq += delta2_X_tmp * delta2_X_tmp;
		// vprod += delta_GX[i] * delta2_X[i];
		// ssq += delta2_X[i] * delta2_X[i];
	}

	bool res = false;

	if(ssq == 0){
		res = true;
	} else {
		double coef = vprod/ssq;

		// update of X:
		for(int i=0 ; i<nb_coef_no_K ; ++i){
			X[i] = GGX[i] - coef * delta_GX[i];
		}
	}

	return(res);
}

void CCC_poisson(int n_obs, int nb_cluster,
                 double *cluster_coef, double *exp_mu,
                 double *sum_y, int *dum){
	// compute cluster coef, poisson
	// Rprintf("in gaussian\n");

	// initialize cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = 0;
	}

	// looping sequentially over exp_mu
	for(int i=0 ; i<n_obs ; ++i){
		cluster_coef[dum[i]] += exp_mu[i];
	}

	// calculating cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = sum_y[m] / cluster_coef[m];
	}

	// "output" is the update of my_cluster_coef
}

void CCC_poisson_log(int n_obs, int nb_cluster,
                     double *cluster_coef, double *mu,
                     double *sum_y, int *dum){
	// compute cluster coef, poisson
	// This is log poisson => we are there because classic poisson did not work
	// Thus high chance there are very high values of the cluster coefs (in abs value)
	// we need to take extra care in computing it => we apply trick of substracting the max in the exp

	vector<double> mu_max(nb_cluster);
	vector<bool> doInit(nb_cluster);

	// initialize cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = 0;
		doInit[m] = true;
	}

	// finding the max mu for each cluster
	int d;
	for(int i=0 ; i<n_obs ; ++i){
		d = dum[i];
		if(doInit[d]){
			mu_max[d] = mu[i];
			doInit[d] = false;
		} else if(mu[i] > mu_max[d]){
			mu_max[d] = mu[i];
		}
	}

	// looping sequentially over exp_mu
	for(int i=0 ; i<n_obs ; ++i){
		d = dum[i];
		cluster_coef[d] += exp(mu[i] - mu_max[d]);
	}

	// calculating cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = log(sum_y[m]) - log(cluster_coef[m]) - mu_max[m];
	}

	// "output" is the update of my_cluster_coef
}


void CCC_gaussian(int n_obs, int nb_cluster,
                  double *cluster_coef, double *mu,
                  double *sum_y, int *dum, int *table){
	// compute cluster coef, gaussian

	// initialize cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = 0;
	}

	// looping sequentially over mu
	for(int i=0 ; i<n_obs ; ++i){
		cluster_coef[dum[i]] += mu[i];
	}

	// calculating cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = (sum_y[m] - cluster_coef[m]) / table[m];
	}

	// "output" is the update of my_cluster_coef
}

void CCC_negbin(int nthreads, int nb_cluster, double theta, double diffMax_NR,
                double *cluster_coef, double *mu,
                double *lhs, double *sum_y, int *obsCluster, int *table, int *cumtable){
	// compute cluster coefficients negbin
	// This is not direct: needs to apply dichotomy+NR algorithm

	// first we find the min max for each cluster to get the bounds
	int iterMax = 100, iterFullDicho = 10;

	// finding the max/min values of mu for each cluster
	vector<double> borne_inf(nb_cluster);
	vector<double> borne_sup(nb_cluster);
	// attention borne_inf => quand mu est maximal
	int u0;
	double value, mu_min, mu_max;

	for(int m=0 ; m<nb_cluster ; ++m){
		// the min/max of mu
		u0 = (m == 0 ? 0 : cumtable[m - 1]);
		mu_min = mu[obsCluster[u0]];
		mu_max = mu[obsCluster[u0]];
		for(int u = 1+u0 ; u<cumtable[m] ; ++u){
			value = mu[obsCluster[u]];
			if(value < mu_min){
				mu_min = value;
			} else if(value > mu_max){
				mu_max = value;
			}
		}

		// computing the bounds
		borne_inf[m] = log(sum_y[m]) - log(static_cast<double>(table[m])) - mu_max;
		borne_sup[m] = borne_inf[m] + (mu_max - mu_min);
	}

	// Rprintf("inf: %f -- sup: %f -- middle: %f\n", borne_inf[0], borne_sup[0], (borne_inf[0] + borne_sup[0])/2);

    #pragma omp parallel for num_threads(nthreads)
	for(int m=0 ; m<nb_cluster ; ++m){
		// we loop over each cluster

		// we initialise the cluster coefficient at 0 (it should converge to 0 at some point)
		double x1 = 0;
		bool keepGoing = true;
		int iter = 0, i;
		int u0 = (m == 0 ? 0 : cumtable[m - 1]);

		double value, x0, derivee = 0, exp_mu;

		// the bounds
		double lower_bound = borne_inf[m];
		double upper_bound = borne_sup[m];

		// Update of the value if it goes out of the boundaries
		// because we dont know ex ante if 0 is within the bounds
		if(x1 >= upper_bound || x1 <= lower_bound){
			x1 = (lower_bound + upper_bound)/2;
		}

		while(keepGoing){
			++iter;

			// 1st step: initialisation des bornes

			// computing the value of f(x)
			value = sum_y[m];
			for(int u = u0 ; u<cumtable[m] ; ++u){
				i = obsCluster[u];
				value -= (theta + lhs[i]) / (1 + theta*exp(-x1 - mu[i]));
			}

			// update of the bounds.
			if(value > 0){
				lower_bound = x1;
			} else {
				upper_bound = x1;
			}

			// 2nd step: NR iteration or Dichotomy
			x0 = x1;
			if(value == 0){
				keepGoing = false;
			} else if(iter <= iterFullDicho){
				// computing the derivative
				derivee = 0;
				for(int u = u0 ; u<cumtable[m] ; ++u){
					i = obsCluster[u];
					exp_mu = exp(x1 + mu[i]);
					derivee -= theta * (theta + lhs[i]) / ( (theta/exp_mu + 1) * (theta + exp_mu) );
				}

				x1 = x0 - value / derivee;
				// Rprintf("x1: %5.2f\n", x1);

				// 3rd step: dichotomy (if necessary)
				// Update of the value if it goes out of the boundaries
				if(x1 >= upper_bound || x1 <= lower_bound){
					x1 = (lower_bound + upper_bound)/2;
				}
			} else {
				x1 = (lower_bound + upper_bound)/2;
			}

			// the stopping criterion
			if(iter == iterMax){
				keepGoing = false;
			    if(omp_get_thread_num() == 0){
			        Rprintf("[Getting cluster coefficients nber %i] max iterations reached (%i).\n", m, iterMax);
			        Rprintf("Value Sum Deriv (NR) = %f. Difference = %f.\n", value, fabs(x0-x1));
			    }
			}

			// if(fabs(x0-x1) / (0.1 + fabs(x1)) < diffMax_NR){
			if(stopping_criterion(x0, x1, diffMax_NR)){
				keepGoing = false;
			}

		}

		// if(m == 0) Rprintf("diffMax = %f, value = %f, test = %f\n", fabs(x0-x1), value, fabs(5.1 - 5));
		// if(m == 0) Rprintf("total iter = %i, final x1 = %f\n", iter, x1);

		// after convegence: only update of cluster coef
		cluster_coef[m] = x1;
	}

}

void CCC_logit(int nthreads, int nb_cluster, double diffMax_NR,
               double *cluster_coef, double *mu,
               double *sum_y, int *obsCluster, int *table, int *cumtable){
	// compute cluster coefficients negbin
	// This is not direct: needs to apply dichotomy+NR algorithm

	// first we find the min max for each cluster to get the bounds
	int iterMax = 100, iterFullDicho = 10;

	// finding the max/min values of mu for each cluster
	vector<double> borne_inf(nb_cluster);
	vector<double> borne_sup(nb_cluster);
	// attention borne_inf => quand mu est maximal

	int u0;
	double value, mu_min, mu_max;
	for(int m=0 ; m<nb_cluster ; ++m){
		// the min/max of mu
		u0 = (m == 0 ? 0 : cumtable[m - 1]);
		mu_min = mu[obsCluster[u0]];
		mu_max = mu[obsCluster[u0]];
		for(int u = 1+u0 ; u<cumtable[m] ; ++u){
			value = mu[obsCluster[u]];
			if(value < mu_min){
				mu_min = value;
			} else if(value > mu_max){
				mu_max = value;
			}
		}

		// computing the "bornes"
		borne_inf[m] = log(sum_y[m]) - log(table[m]-sum_y[m]) - mu_max;
		borne_sup[m] = borne_inf[m] + (mu_max - mu_min);
	}


	//
	// Parallel loop
	//

    #pragma omp parallel for num_threads(nthreads)
	for(int m=0 ; m<nb_cluster ; ++m){
		// we loop over each cluster

		// we initialise the cluster coefficient at 0 (it should converge to 0 at some point)
		double x1 = 0;
		bool keepGoing = true;
		int iter = 0;
		int u0 = (m == 0 ? 0 : cumtable[m - 1]);

		double value, x0, derivee = 0, exp_mu;

		// the bounds
		double lower_bound = borne_inf[m];
		double upper_bound = borne_sup[m];

		// Update of the value if it goes out of the boundaries
		// because we dont know ex ante if 0 is within the bounds
		if(x1 >= upper_bound || x1 <= lower_bound){
			x1 = (lower_bound + upper_bound)/2;
		}

		while(keepGoing){
			++iter;

			// 1st step: initialisation des bornes

			// computing the value of f(x)
			value = sum_y[m];
			for(int u = u0 ; u<cumtable[m] ; ++u){
				value -= 1 / (1 + exp(-x1 - mu[obsCluster[u]]));
			}

			// update of the bounds.
			if(value > 0){
				lower_bound = x1;
			} else {
				upper_bound = x1;
			}

			// 2nd step: NR iteration or Dichotomy
			x0 = x1;
			if(value == 0){
				keepGoing = false;
			} else if(iter <= iterFullDicho){
				// computing the derivative
				derivee = 0;
				for(int u = u0 ; u<cumtable[m] ; ++u){
					exp_mu = exp(x1 + mu[obsCluster[u]]);
					derivee -= 1 / ( (1/exp_mu + 1) * (1 + exp_mu) );
				}

				x1 = x0 - value / derivee;
				// Rprintf("x1: %5.2f\n", x1);

				// 3rd step: dichotomy (if necessary)
				// Update of the value if it goes out of the boundaries
				if(x1 >= upper_bound || x1 <= lower_bound){
					x1 = (lower_bound + upper_bound)/2;
				}
			} else {
				x1 = (lower_bound + upper_bound)/2;
			}

			// the stopping criteria
			if(iter == iterMax){
				keepGoing = false;
				Rprintf("[Getting cluster coefficients nber %i] max iterations reached (%i).\n", m, iterMax);
				Rprintf("Value Sum Deriv (NR) = %f. Difference = %f.\n", value, fabs(x0-x1));
			}

			// if(fabs(x0-x1) / (0.1 + fabs(x1)) < diffMax_NR){
		    if(stopping_criterion(x0, x1, diffMax_NR)){
				keepGoing = false;
			}

		}
		// Rprintf("iter=%i.\n", iter);
		// res[k] = iter;
		// res[k] = value;

		// after convegence: only update of cluster coef
		cluster_coef[m] = x1;
	}

}

void computeClusterCoef_single(int family, int n_obs, int nb_cluster, double theta, double diffMax_NR,
                               double *cluster_coef, double *mu,
                               double *lhs, double *sum_y,
                               int *dum, int *obsCluster, int *table, int *cumtable, int nthreads){

	// Leads to the appropriate function
	// we update the cluster "in place" (ie using pointers)
	// switch is peculiar... we need break too

	// comment on mutlithreading:
	// I know I could use only the CCC_par_negbin/CCC_par_logit functions but I'm scared that
	// in some instance openmp creates some bug, thus when single core, there is no reference
	// to openmp at all


	switch(family){
	case 1:
		CCC_poisson(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum);
		break;
	case 2: // Negbin
		CCC_negbin(nthreads, nb_cluster, theta, diffMax_NR, cluster_coef, mu, lhs, sum_y, obsCluster, table, cumtable);
		break;
	case 3: // logit
		CCC_logit(nthreads, nb_cluster, diffMax_NR, cluster_coef, mu, sum_y, obsCluster, table, cumtable);
		break;
	case 4: // Gaussian
		CCC_gaussian(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum, table);
		break;
	case 5: // log poisson
		CCC_poisson_log(n_obs, nb_cluster, cluster_coef, mu, sum_y, dum);
		break;
	}

}

// Function to delete => only for drbugging
// [[Rcpp::export]]
SEXP compute_cluster_coef_r(int family, int nb_coef, double theta, double diffMax_NR,
                            SEXP r_mu, SEXP r_lhs, SEXP r_sum_y, SEXP r_dum,
                            SEXP r_obsCluster, SEXP r_table, SEXP r_cumtable, int nthreads){

	int n_obs = Rf_length(r_mu);

	// pointers to R values
	double *mu = REAL(r_mu);
	double *lhs = REAL(r_lhs);
	double *sum_y = REAL(r_sum_y);
	int *dum = INTEGER(r_dum);
	int *obsCluster = INTEGER(r_obsCluster);
	int *table = INTEGER(r_table);
	int *cumtable = INTEGER(r_cumtable);

	SEXP res = PROTECT(Rf_allocVector(REALSXP, nb_coef));
	double *pcoef = REAL(res);

	computeClusterCoef_single(family, n_obs, nb_coef, theta, diffMax_NR,
                           pcoef, mu, lhs, sum_y, dum, obsCluster, table, cumtable, nthreads);

	UNPROTECT(1);

	return(res);
}

// [[Rcpp::export]]
SEXP update_mu_single_cluster(int family, int nb_cluster, double theta, double diffMax_NR, SEXP mu_in,
                              SEXP lhs, SEXP sum_y, SEXP dum, SEXP obsCluster, SEXP table,
                              SEXP cumtable, int nthreads){
	// Function used to compute the cluster coefficient for ONE fixed-effect only
	// and then add it to the existing mu

	int n_obs = Rf_length(mu_in);

	// Getting the pointers
	int *pdum = INTEGER(dum);
	int *pobsCluster = INTEGER(obsCluster);
	int *ptable = INTEGER(table);
	int *pcumtable = INTEGER(cumtable);
	double *plhs = REAL(lhs);
	double *psum_y = REAL(sum_y);
	double *pmu_in = REAL(mu_in);

	// The cluster coeffficients
	vector<double> cluster_coef(nb_cluster);

	computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR, cluster_coef.data(), pmu_in,
                           plhs, psum_y, pdum, pobsCluster, ptable, pcumtable, nthreads);

	// the result
	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);

	if(family == 1){
		// Poisson, exp form
		for(int i=0 ; i<n_obs ; ++i){
			pmu[i] = pmu_in[i] * cluster_coef[pdum[i]];
		}
	} else {
		for(int i=0 ; i<n_obs ; ++i){
			pmu[i] = pmu_in[i] + cluster_coef[pdum[i]];
		}
	}

	UNPROTECT(1);

	return(mu);
}

void computeClusterCoef(vector<double*> &pcluster_origin, vector<double*> &pcluster_destination,
                        PARAM_CCC *args){
	// update of the cluster coefficients
	// first we update mu, then we update the cluster coefficicents

	//
	// Loading the variables
	//

	int family = args->family;
	int n_obs = args->n_obs;
	int K = args->K;
	int nthreads = args->nthreads;
	double theta = args->theta;
	double diffMax_NR = args->diffMax_NR;

	int *pcluster = args->pcluster;
	double *lhs = args->lhs;
	double *mu_init = args->mu_init;

	vector<int*> &pdum = args->pdum;
	vector<int*> &ptable = args->ptable;
	vector<double*> &psum_y = args->psum_y;
	vector<int*> &pobsCluster = args->pobsCluster;
	vector<int*> &pcumtable = args->pcumtable;

	// value that will be modified
	double *mu_with_coef = args->mu_with_coef;

	// We update each cluster coefficient, starting from K

	// we first set the value of mu_with_coef
	for(int i=0 ; i<n_obs ; ++i){
		mu_with_coef[i] = mu_init[i];
	}

	for(int k=0 ; k<(K-1) ; ++k){
		int *my_dum = pdum[k];
		double *my_cluster_coef = pcluster_origin[k];

		if(family == 1){ // Poisson
			for(int i=0 ; i<n_obs ; ++i){
				mu_with_coef[i] *= my_cluster_coef[my_dum[i]];
			}
		} else {
			for(int i=0 ; i<n_obs ; ++i){
				mu_with_coef[i] += my_cluster_coef[my_dum[i]];
			}
		}
	}


	for(int k=K-1 ; k>=0 ; k--){
		R_CheckUserInterrupt();

		// computing the optimal cluster coef -- given mu_with_coef
		double *my_cluster_coef = pcluster_destination[k];
		int *my_table = ptable[k];
		double *my_sum_y = psum_y[k];
		int *my_dum = pdum[k];
		int *my_cumtable = pcumtable[k];
		int *my_obsCluster = pobsCluster[k];
		int nb_cluster = pcluster[k];

		// update of the cluster coefficients
		computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR,
                            my_cluster_coef, mu_with_coef, lhs, my_sum_y,
                            my_dum, my_obsCluster, my_table, my_cumtable, nthreads);


		// updating the value of mu_with_coef (only if necessary)
		if(k != 0){

			//
			// Test 1 => add then withdraw
			//

			// // we add the computed value
			// if(family == 1){
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] *= my_cluster_coef[my_dum[i]];
			// 	}
			// } else {
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] += my_cluster_coef[my_dum[i]];
			// 	}
			// }
			//
			// // and withdraw the next value (awaiting to be computed)
			// int *my_dum = pdum[k-1];
			// double *my_cluster_coef = pcluster_origin[k-1];
			//
			// if(family == 1){
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] /= my_cluster_coef[my_dum[i]];
			// 	}
			// } else {
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] -= my_cluster_coef[my_dum[i]];
			// 	}
			// }



			//
			// Test 2 => withdraw then add
			//

			// // we take out the next
			// int *my_dum = pdum[k-1];
			// double *my_cluster_coef = pcluster_origin[k-1];
			//
			// if(family == 1){
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] /= my_cluster_coef[my_dum[i]];
			// 	}
			// } else {
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] -= my_cluster_coef[my_dum[i]];
			// 	}
			// }
			//
			// // and add the current
			// my_dum = pdum[k];
			// my_cluster_coef = pcluster_destination[k];
			//
			// if(family == 1){
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] *= my_cluster_coef[my_dum[i]];
			// 	}
			// } else {
			// 	for(i=0 ; i<n_obs ; ++i){
			// 		mu_with_coef[i] += my_cluster_coef[my_dum[i]];
			// 	}
			// }

			//
			// Test 3: recompute from scratch => we keep this one
			//


			for(int i=0 ; i<n_obs ; ++i){
				mu_with_coef[i] = mu_init[i];
			}

			int *my_dum;
			double *my_cluster_coef;
			for(int h=0 ; h<K ; h++){
				if(h == k-1) continue;

				my_dum = pdum[h];

				if(h < k-1){
					my_cluster_coef = pcluster_origin[h];
				} else {
					my_cluster_coef = pcluster_destination[h];
				}

				if(family == 1){ // Poisson
					for(int i=0 ; i<n_obs ; ++i){
						mu_with_coef[i] *= my_cluster_coef[my_dum[i]];
					}
				} else {
					for(int i=0 ; i<n_obs ; ++i){
						mu_with_coef[i] += my_cluster_coef[my_dum[i]];
					}
				}
			}

		}
	}

	// In the end, the array pcluster_coef is fully updated, starting from K to 1

}

// [[Rcpp::export]]
List cpp_conv_acc_gnl(int family, int iterMax, double diffMax, double diffMax_NR, double theta, SEXP nb_cluster_all,
                 SEXP lhs, SEXP mu_init, SEXP dum_vector, SEXP tableCluster_vector,
                 SEXP sum_y_vector, SEXP cumtable_vector, SEXP obsCluster_vector, int nthreads){

	//initial variables
	int K = Rf_length(nb_cluster_all);
	int *pcluster = INTEGER(nb_cluster_all);
	int n_obs = Rf_length(mu_init);
	double *pmu_init = REAL(mu_init);

	int nb_coef=0;
	for(int k=0 ; k<K ; ++k){
		nb_coef += pcluster[k];
	}

	// setting the pointers
	// table, coef, sum_y: length nb_coef
	// dum_all: length K*n_obs
	vector<int*> ptable(K);
	vector<double*> psum_y(K);
	ptable[0] = INTEGER(tableCluster_vector);
	psum_y[0] = REAL(sum_y_vector);

	for(int k=1 ; k<K ; ++k){
		ptable[k] = ptable[k - 1] + pcluster[k - 1];
		psum_y[k] = psum_y[k - 1] + pcluster[k - 1];
	}

	// cumtable (points to nothing for non negbin/logit families)
	vector<int*> pcumtable(K);
	if(family == 2 || family == 3){
		pcumtable[0] = INTEGER(cumtable_vector);
		for(int k=1 ; k<K ; ++k){
			pcumtable[k] = pcumtable[k - 1] + pcluster[k - 1];
		}
	}

	// obsCluster (points to nothing for non negbin/logit families)
	vector<int*> pobsCluster(K);
	if(family == 2 || family == 3){
		pobsCluster[0] = INTEGER(obsCluster_vector);
		for(int k=1 ; k<K ; ++k){
			pobsCluster[k] = pobsCluster[k - 1] + n_obs;
		}
	}

	// cluster of each observation
	vector<int*> pdum(K);
	pdum[0] = INTEGER(dum_vector);
	for(int k=1 ; k<K ; ++k){
		pdum[k] = pdum[k - 1] + n_obs;
	}

	// lhs (only negbin will use it)
	double *plhs = REAL(lhs);

	//
	// Sending variables to envir
	//

	PARAM_CCC args;

	args.family = family;
	args.n_obs = n_obs;
	args.K = K;
	args.nthreads = nthreads;
	args.theta = (family == 2 ? theta : 1); // theta won't be used if family not negbin
	args.diffMax_NR = diffMax_NR;
	args.pdum = pdum;
	args.mu_init = pmu_init;
	args.ptable = ptable;
	args.psum_y = psum_y;
	args.pcluster = pcluster;
	args.pcumtable = pcumtable;
	args.pobsCluster = pobsCluster;
	args.lhs = plhs;

	// value that will be modified
	vector<double> mu_with_coef(n_obs);
	args.mu_with_coef = mu_with_coef.data();

	//
	// IT iteration (preparation)
	//

	// variables on 1:K
	vector<double> X(nb_coef);
	vector<double> GX(nb_coef);
	vector<double> GGX(nb_coef);
	// pointers:
	vector<double*> pX(K);
	vector<double*> pGX(K);
	vector<double*> pGGX(K);
	pX[0] = X.data();
	pGX[0] = GX.data();
	pGGX[0] = GGX.data();

	for(int k=1 ; k<K ; ++k){
		pX[k] = pX[k - 1] + pcluster[k - 1];
		pGX[k] = pGX[k - 1] + pcluster[k - 1];
		pGGX[k] = pGGX[k - 1] + pcluster[k - 1];
	}

	// variables on 1:(K-1)
	int nb_coef_no_K = 0;
	for(int k = 0 ; k<(K-1) ; ++k){
		nb_coef_no_K += pcluster[k];
	}
	vector<double> delta_GX(nb_coef_no_K);
	vector<double> delta2_X(nb_coef_no_K);

	//
	// the main loop
	//

	// initialisation of X and then pGX
	if(family == 1){
		for(int i=0 ; i<nb_coef ; ++i){
			X[i] = 1;
		}
	} else {
		for(int i=0 ; i<nb_coef ; ++i){
			X[i] = 0;
		}
	}


	// first iteration
	computeClusterCoef(pX, pGX, &args);

	// Rprintf("init: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", GX[i]);
	// }
	// Rprintf("\n");

	// stop("kkk");

	// flag for problem with poisson
	bool any_negative_poisson = false;

	// check whether we should go into the loop
	bool keepGoing = false;
	for(int i=0 ; i<nb_coef ; ++i){
		// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
		if(continue_criterion(X[i], GX[i], diffMax)){
			keepGoing = true;
			break;
		}
	}

	int iter = 0;
	bool numconv = false;
	while(keepGoing && iter<iterMax){
		++iter;

		// GGX -- origin: GX, destination: GGX
		computeClusterCoef(pGX, pGGX, &args);

		// Rprintf("iter: %i ; ", iter);
		//
		// if(iter >= iterMax - 3){
		// 	Rprintf("GGX: ");
		// 	for(int i=0 ; i<8 ; ++i){
		// 		Rprintf("%2.7f ", GGX[i]);
		// 	}
		// 	Rprintf("\n");
		// }


		// X ; update of the cluster coefficient
		numconv = update_X_IronsTuck(nb_coef_no_K, X, GX, GGX, delta_GX, delta2_X);
		if(numconv) break;

		// if(iter >= iterMax - 3){
		// 	Rprintf(" IT: ");
		// 	for(int i=0 ; i<8 ; ++i){
		// 		Rprintf("%2.7f ", X[i]);
		// 	}
		// 	Rprintf("\n");
		// }


		if(family == 1){
			// We control for possible problems with poisson
			for(int i=0 ; i<nb_coef_no_K ; ++i){
				if(X[i] <= 0){
					any_negative_poisson = true;
					break;
				}
			}

			if(any_negative_poisson){
				break; // we quit the loop
				// update of mu is OK, it's only IT iteration that leads to negative values
			}
		}

		// GX -- origin: X, destination: GX
		computeClusterCoef(pX, pGX, &args);

		keepGoing = false;
		for(int i=0 ; i<nb_coef_no_K ; ++i){
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if(continue_criterion(X[i], GX[i], diffMax)){
				keepGoing = true;
				break;
			}
		}

	}

	//
	// We update mu => result
	//

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);
	for(int i=0 ; i<n_obs ; ++i){
		pmu[i] = pmu_init[i];
	}

	// // regular code (using last iteration)
	// for(int k=0 ; k<K ; ++k){
	// 	int *my_dum = pdum[k];
	// 	double *my_cluster_coef = pGX[k]; // GX is last updated
	// 	if(family == 1){ // "quick" Poisson case
	// 		for(int i=0 ; i<n_obs ; ++i){
	// 			pmu[i] *= my_cluster_coef[my_dum[i]];
	// 		}
	// 	} else {
	// 		for(int i=0 ; i<n_obs ; ++i){
	// 			pmu[i] += my_cluster_coef[my_dum[i]];
	// 		}
	// 	}
	// }

	// To obtain stg identical to the R code:
	computeClusterCoef(pGX, pGGX, &args);
	for(int k=0 ; k<K ; ++k){
		int *my_dum = pdum[k];
		double *my_cluster_coef = pGGX[k];

		if(family == 1){ // "quick" Poisson case
			for(int i=0 ; i<n_obs ; ++i){
				pmu[i] *= my_cluster_coef[my_dum[i]];
			}
		} else {
			for(int i=0 ; i<n_obs ; ++i){
				pmu[i] += my_cluster_coef[my_dum[i]];
			}
		}
	}

	UNPROTECT(1);

	List res;
	res["mu_new"] = mu;
	res["iter"] = iter;
	res["any_negative_poisson"] = any_negative_poisson;

	return(res);
}


// [[Rcpp::export]]
List cpp_conv_seq_gnl(int family, int iterMax, double diffMax, double diffMax_NR, double theta, SEXP nb_cluster_all,
                 SEXP lhs, SEXP mu_init, SEXP dum_vector, SEXP tableCluster_vector,
                 SEXP sum_y_vector, SEXP cumtable_vector, SEXP obsCluster_vector, int nthreads){

	//initial variables
	int K = Rf_length(nb_cluster_all);
	int *pcluster = INTEGER(nb_cluster_all);
	int n_obs = Rf_length(mu_init);
	double *pmu_init = REAL(mu_init);

	int nb_coef=0;
	for(int k=0 ; k<K ; ++k){
		nb_coef += pcluster[k];
	}

	// setting the pointers
	// table, coef, sum_y: length nb_coef
	// dum_all: length K*n_obs
	vector<int*> ptable(K);
	vector<double*> psum_y(K);
	ptable[0] = INTEGER(tableCluster_vector);
	psum_y[0] = REAL(sum_y_vector);

	for(int k=1 ; k<K ; ++k){
		ptable[k] = ptable[k - 1] + pcluster[k - 1];
		psum_y[k] = psum_y[k - 1] + pcluster[k - 1];
	}

	// cumtable (points to nothing for non negbin/logit families)
	vector<int*> pcumtable(K);
	if(family == 2 || family == 3){
		pcumtable[0] = INTEGER(cumtable_vector);
		for(int k=1 ; k<K ; ++k){
			pcumtable[k] = pcumtable[k - 1] + pcluster[k - 1];
		}
	}

	// obsCluster (points to nothing for non negbin/logit families)
	vector<int*> pobsCluster(K);
	if(family == 2 || family == 3){
		pobsCluster[0] = INTEGER(obsCluster_vector);
		for(int k=1 ; k<K ; ++k){
			pobsCluster[k] = pobsCluster[k - 1] + n_obs;
		}
	}

	// cluster of each observation
	vector<int*> pdum(K);
	pdum[0] = INTEGER(dum_vector);
	for(int k=1 ; k<K ; ++k){
		pdum[k] = pdum[k - 1] + n_obs;
	}

	// lhs (only negbin will use it)
	double *plhs = REAL(lhs);

	// the variable with the value of mu
	vector<double> mu_with_coef(n_obs);
	// initialization of mu_with_coef
	for(int i=0 ; i<n_obs ; ++i){
		mu_with_coef[i] = pmu_init[i];
	}

	// The cluster coefficients
	// variables on 1:K
	vector<double> cluster_coef(nb_coef);
	// pointers:
	vector<double*> pcluster_coef(K);
	pcluster_coef[0] = cluster_coef.data();
	for(int k=1 ; k<K ; ++k){
		pcluster_coef[k] = pcluster_coef[k - 1] + pcluster[k - 1];
	}

	//
	// the main loop
	//

	// initialisation of the cluster coefficients
	if(family == 1){ // Poisson
		for(int i=0 ; i<nb_coef ; ++i){
			cluster_coef[i] = 1;
		}
	} else {
		for(int i=0 ; i<nb_coef ; ++i){
			cluster_coef[i] = 0;
		}
	}

	// the main loop
	bool keepGoing = true;
	int iter = 1;
	while(keepGoing && iter <= iterMax){
		// Rprintf("\n%i.", iter);
		++iter;
		keepGoing = false;

		/// we loop over all clusters => from K to 1
		for(int k=(K-1) ; k>=0 ; k--){
			// Rprintf("k=%i ", k);
			R_CheckUserInterrupt();

			//
			// 1) computing the cluster coefficient
			//

			// computing the optimal cluster coef -- given mu_with_coef
			double *my_cluster_coef = pcluster_coef[k];
			int *my_table = ptable[k];
			double *my_sum_y = psum_y[k];
			int *my_dum = pdum[k];
			int *my_cumtable = pcumtable[k];
			int *my_obsCluster = pobsCluster[k];
			int nb_cluster = pcluster[k];

			// update of the cluster coefficients
			computeClusterCoef_single(family, n_obs, nb_cluster, theta, diffMax_NR,
                             my_cluster_coef, mu_with_coef.data(), plhs, my_sum_y,
                             my_dum, my_obsCluster, my_table, my_cumtable, nthreads);

			//
			// 2) Updating the value of mu
			//

			// we add the computed value
			if(family == 1){
				for(int i=0 ; i<n_obs ; ++i){
					mu_with_coef[i] *= my_cluster_coef[my_dum[i]];
				}
			} else {
				for(int i=0 ; i<n_obs ; ++i){
					mu_with_coef[i] += my_cluster_coef[my_dum[i]];
				}
			}

			// Stopping criterion
			if(keepGoing == false){

				// Rprintf("nb: %i ", nb_cluster);
				// Rprintf("value: %f", my_cluster_coef[nb_cluster]);
				// stop("hhhh");

				if(family == 1){
					for(int m=0 ; m<nb_cluster ; ++m){
						if(fabs(my_cluster_coef[m]-1) > diffMax){
							keepGoing = true;
							break;
						}
					}
				} else {
					for(int m=0 ; m<nb_cluster ; ++m){
						if(fabs(my_cluster_coef[m]) > diffMax){
							keepGoing = true;
							break;
						}
					}
				}

			}

		}
	}

	//
	// We update mu => result
	//

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);
	for(int i=0 ; i<n_obs ; ++i){
		pmu[i] = mu_with_coef[i];
	}

	UNPROTECT(1);

	List res;
	res["mu_new"] = mu;
	res["iter"] = iter;

	return(res);
}

// [[Rcpp::export]]
int get_n_cells(IntegerVector index_i, IntegerVector index_j){

	int n = index_i.length();

	// we count the nber of different elements
	int index_current = 0;

	for(int i=1 ; i<n ; ++i){
		if(index_j[i] != index_j[i-1] || index_i[i] != index_i[i-1]){
			// new row
			index_current++;
		}
	}

	index_current++;

	return(index_current);
}



void CCC_poisson_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
                   int n_i, int n_j, int n_cells,
                   const vector<int> &mat_row, vector<int> &mat_col, vector<double> &mat_value,
                   const vector<double> &ca, const vector<double> &cb,
                   vector<double> &alpha){

	// alpha = ca / (Ab %m% (cb / (Ab %tm% alpha)))

	double *beta = pcluster_destination.data() + n_i;

	for(int i=0 ; i<n_i ; ++i){
		alpha[i] = 0;
	}

	for(int j=0 ; j<n_j ; ++j){
		beta[j] = 0;
	}

	for(int obs=0 ; obs<n_cells ; ++obs){
		beta[mat_col[obs]] += mat_value[obs]*pcluster_origin[mat_row[obs]];
	}

	for(int j=0 ; j<n_j ; ++j){
		beta[j] = cb[j] / beta[j];
	}

	for(int obs=0 ; obs<n_cells ; ++obs){
		alpha[mat_row[obs]] += mat_value[obs]*beta[mat_col[obs]];
	}

	for(int i=0 ; i<n_i ; ++i){
		pcluster_destination[i] = ca[i] / alpha[i];
	}

}

// [[Rcpp::export]]
List cpp_conv_acc_poi_2(int n_i, int n_j, int n_cells, SEXP index_i, SEXP index_j,
                        SEXP dum_vector, SEXP sum_y_vector,
                        int iterMax, double diffMax, SEXP exp_mu_in, SEXP order){



	// values that will be used later
	vector<double> alpha(n_i);

	// We compute the matrix Ab
	int index_current = 0;
	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);

	int n_obs = Rf_length(exp_mu_in);

	int *new_order = INTEGER(order);
	double *pexp_mu_in = REAL(exp_mu_in);
	// double value = pexp_mu_in[0];
	double value = pexp_mu_in[new_order[0]];

	for(int i=1 ; i<n_obs ; ++i){
		if(pindex_j[i] != pindex_j[i-1] || pindex_i[i] != pindex_i[i-1]){
			// save the value if we change index
			mat_row[index_current] = pindex_i[i-1];
			mat_col[index_current] = pindex_j[i-1];
			mat_value[index_current] = value;

			// new row
			index_current++;

			// the new value
			value = pexp_mu_in[new_order[i]];
		} else {
			value += pexp_mu_in[new_order[i]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value[index_current] = value;

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i+n_j);
	vector<double> GX(n_i+n_j);
	vector<double> GGX(n_i+n_j);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);

	//
	// the main loop
	//

	// initialisation of X and then GX
	for(int i=0 ; i<n_i ; ++i){
		X[i] = 1;
	}

	// ca and cb
	double *psum_y = REAL(sum_y_vector);

	vector<double> ca(n_i);
	vector<double> cb(n_j);
	for(int i=0 ; i<n_i ; ++i){
		ca[i] = psum_y[i];
	}

	double *psum_y_j = psum_y + n_i;
	for(int j=0 ; j<n_j ; ++j){
		cb[j] = psum_y_j[j];
	}

	// first iteration
	CCC_poisson_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);

	// Rprintf("init: ");
	// for(i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", GX[i]);
	// }
	// Rprintf("\n");

	// stop("kkk");

	// flag for problem with poisson
	bool any_negative_poisson = false;

	// check whether we should go into the loop => always 1 iter
	bool keepGoing = true;
	for(int i=0 ; i<n_i ; ++i){
		// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
		if(continue_criterion(X[i], GX[i], diffMax)){
			keepGoing = true;
			break;
		}
	}

	bool numconv = false;
	int iter = 0;
	while(keepGoing && iter<iterMax){
		// Rprintf("%i.", iter);
		++iter;

		// GGX -- origin: GX, destination: GGX
		CCC_poisson_2(GX, GGX, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);

		// Rprintf("ggx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f -- ", GGX[i]);
		// }
		// Rprintf("\n");

		// X ; update of the cluster coefficient
		numconv = update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
		if(numconv) break;

		// Control for negative values
		for(int i=0 ; i<n_i ; ++i){
			if(X[i] <= 0){
				any_negative_poisson = true;
				break;
			}
		}

		if(any_negative_poisson){
			break; // we quit the loop
			// update of mu is OK, it's only IT iteration that leads to negative values
		}

		// Rprintf("  x: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f -- ", X[i]);
		// }
		// Rprintf("\n");

		// GX -- origin: X, destination: GX
		CCC_poisson_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);

		// Rprintf("   gx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f\t", GX[i]);
		// }
		// Rprintf("\n");

		keepGoing = false;
		for(int i=0 ; i<n_i ; ++i){
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if(continue_criterion(X[i], GX[i], diffMax)){
				keepGoing = true;
				break;
			}
		}

	}

	SEXP exp_mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(exp_mu);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	// double *beta = GX + n_i;
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	pmu[obs] = pexp_mu_in[obs] * GX[dum_i[obs]] * beta[dum_j[obs]];
	// }

	// pour avoir identique a acc_pois
	CCC_poisson_2(GX, X, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);
	double *beta = X.data() + n_i;
	for(int obs=0 ; obs<n_obs ; ++obs){
		pmu[obs] = pexp_mu_in[obs] * X[dum_i[obs]] * beta[dum_j[obs]];
	}

	// Rprintf("alpha: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", X[i]);
	// }
	// Rprintf("\n");

	UNPROTECT(1);

	List res;
	res["mu_new"] = exp_mu;
	res["iter"] = iter;
	res["any_negative_poisson"] = any_negative_poisson;

	return(res);
}


// [[Rcpp::export]]
List cpp_conv_seq_poi_2(int n_i, int n_j, int n_cells, SEXP index_i, SEXP index_j,
                        SEXP dum_vector, SEXP sum_y_vector,
                        int iterMax, double diffMax, SEXP exp_mu_in, SEXP order){



	// values that will be used later
	vector<double> alpha(n_i);

	// We compute the matrix Ab
	int index_current = 0;
	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);

	int n_obs = Rf_length(exp_mu_in);

	int *new_order = INTEGER(order);
	double *pexp_mu_in = REAL(exp_mu_in);
	// double value = pexp_mu_in[0];
	double value = pexp_mu_in[new_order[0]];

	for(int i=1 ; i<n_obs ; ++i){
		if(pindex_j[i] != pindex_j[i-1] || pindex_i[i] != pindex_i[i-1]){
			// save the value if we change index
			mat_row[index_current] = pindex_i[i-1];
			mat_col[index_current] = pindex_j[i-1];
			mat_value[index_current] = value;

			// new row
			index_current++;

			// the new value
			value = pexp_mu_in[new_order[i]];
		} else {
			value += pexp_mu_in[new_order[i]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value[index_current] = value;

	// X, X_new => the vector of coefficients
	vector<double> X_new(n_i+n_j);
	vector<double> X(n_i+n_j);

	//
	// the main loop
	//

	// initialisation of X and then GX
	for(int i=0 ; i<n_i ; ++i){
		X[i] = 1;
	}

	// ca and cb
	double *psum_y = REAL(sum_y_vector);

	vector<double> ca(n_i);
	vector<double> cb(n_j);
	for(int i=0 ; i<n_i ; ++i){
		ca[i] = psum_y[i];
	}

	double *psum_y_j = psum_y + n_i;
	for(int j=0 ; j<n_j ; ++j){
		cb[j] = psum_y_j[j];
	}

	bool keepGoing = true;
	int iter = 0;
	while(keepGoing && iter<iterMax){
		// Rprintf("%i.", iter);
		++iter;

		// This way I don't need to update the values of X with a loop
		if(iter % 2 == 1){
			CCC_poisson_2(X, X_new, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);
		} else {
			CCC_poisson_2(X_new, X, n_i, n_j, n_cells, mat_row, mat_col, mat_value, ca, cb, alpha);
		}

		// double *X_current = (iter % 2 == 1 ? X_new : X);
		// Rprintf("alpha: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f -- ", X_current[i]);
		// }
		// Rprintf("\n");


		keepGoing = false;
		for(int i=0 ; i<n_i ; ++i){
			// if(fabs(X[i] - X_new[i]) / (0.1 + fabs(X_new[i])) > diffMax){
			if(continue_criterion(X[i], X_new[i], diffMax)){
				keepGoing = true;
				break;
			}
		}

	}

	double *X_final = (iter % 2 == 1 ? X_new.data() : X.data());

	SEXP exp_mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(exp_mu);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	double *beta = X_final + n_i;
	for(int obs=0 ; obs<n_obs ; ++obs){
		pmu[obs] = pexp_mu_in[obs] * X_final[dum_i[obs]] * beta[dum_j[obs]];
	}

	// Rprintf(" beta: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", beta[i]);
	// }
	// Rprintf("\n");
	//
	// Rprintf("alpha: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", X_final[i]);
	// }
	// Rprintf("\n");

	UNPROTECT(1);

	List res;
	res["mu_new"] = exp_mu;
	res["iter"] = iter;

	return(res);
}

// [[Rcpp::export]]
List cpp_fixed_cost_gaussian(int n_i, int n_cells, SEXP index_i, SEXP index_j, SEXP order,
                             SEXP invTableCluster_vector, SEXP dum_vector){

	// conversion of R objects
	int n_obs = Rf_length(index_i);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;
	double *invTable_i = REAL(invTableCluster_vector);
	double *invTable_j = invTable_i + n_i;

	// We compute the matrix Ab and Ba
	int index_current = 0;

	// the objects returned
	SEXP r_mat_row = PROTECT(Rf_allocVector(INTSXP, n_cells));
	SEXP r_mat_col = PROTECT(Rf_allocVector(INTSXP, n_cells));
	SEXP r_mat_value_Ab = PROTECT(Rf_allocVector(REALSXP, n_cells));
	SEXP r_mat_value_Ba = PROTECT(Rf_allocVector(REALSXP, n_cells));
	int *mat_row = INTEGER(r_mat_row);
	int *mat_col = INTEGER(r_mat_col);
	double *mat_value_Ab = REAL(r_mat_value_Ab);
	double *mat_value_Ba = REAL(r_mat_value_Ba);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);
	int *new_order = INTEGER(order);

	// double value_Ab = invTable_i[dum_i[0]];
	// double value_Ba = invTable_j[dum_j[0]];
	double value_Ab = invTable_i[dum_i[new_order[0]]];
	double value_Ba = invTable_j[dum_j[new_order[0]]];

	for(int obs=1 ; obs<n_obs ; ++obs){
		if(pindex_j[obs] != pindex_j[obs-1] || pindex_i[obs] != pindex_i[obs-1]){
			// save the value if we change index
			mat_row[index_current] = pindex_i[obs-1];
			mat_col[index_current] = pindex_j[obs-1];
			mat_value_Ab[index_current] = value_Ab;
			mat_value_Ba[index_current] = value_Ba;

			// new row
			index_current++;

			// the new value
			value_Ab = invTable_i[dum_i[new_order[obs]]];
			value_Ba = invTable_j[dum_j[new_order[obs]]];
		} else {
			value_Ab += invTable_i[dum_i[new_order[obs]]];
			value_Ba += invTable_j[dum_j[new_order[obs]]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs - 1];
	mat_col[index_current] = pindex_j[n_obs - 1];
	mat_value_Ab[index_current] = value_Ab;
	mat_value_Ba[index_current] = value_Ba;

	// the object returned
	List res;

	res["mat_row"] = r_mat_row;
	res["mat_col"] = r_mat_col;
	res["mat_value_Ab"] = r_mat_value_Ab;
	res["mat_value_Ba"] = r_mat_value_Ba;

	UNPROTECT(4);

	return(res);
}

void CCC_gaussian_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
                    int n_i, int n_j, int n_cells,
                    int *mat_row, int *mat_col,
                    double *mat_value_Ab, double *mat_value_Ba,
                    const vector<double> &a_tilde, vector<double> &beta){

	// alpha = a_tilde + (Ab %m% (Ba %m% alpha))

	// for(int i=0 ; i<n_i ; ++i){
	// 	alpha[i] = 0;
	// }
	//
	// for(int j=0 ; j<n_j ; ++j){
	// 	beta[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_cells ; ++obs){
	// 	beta[mat_col[obs]] += mat_value_Ba[obs]*pcluster_origin[mat_row[obs]];
	// }
	//
	// for(int obs=0 ; obs<n_cells ; ++obs){
	// 	alpha[mat_row[obs]] += mat_value_Ab[obs]*beta[mat_col[obs]];
	// }
	//
	// for(int i=0 ; i<n_i ; ++i){
	// 	pcluster_destination[i] = a_tilde[i] + alpha[i];
	// }

	for(int i=0 ; i<n_i ; ++i){
		pcluster_destination[i] = a_tilde[i];
	}

	for(int j=0 ; j<n_j ; ++j){
		beta[j] = 0;
	}

	for(int obs=0 ; obs<n_cells ; ++obs){
		beta[mat_col[obs]] += mat_value_Ba[obs]*pcluster_origin[mat_row[obs]];
	}

	for(int obs=0 ; obs<n_cells ; ++obs){
		pcluster_destination[mat_row[obs]] += mat_value_Ab[obs]*beta[mat_col[obs]];
	}

}

// [[Rcpp::export]]
List cpp_conv_acc_gau_2(int n_i, int n_j, int n_cells,
                        SEXP r_mat_row, SEXP r_mat_col, SEXP r_mat_value_Ab, SEXP r_mat_value_Ba,
                        SEXP dum_vector, SEXP lhs, SEXP invTableCluster_vector,
                        int iterMax, double diffMax, SEXP mu_in){

	//
	// Setting up
	//

	int n_obs = Rf_length(mu_in);

	int *mat_row = INTEGER(r_mat_row);
	int *mat_col = INTEGER(r_mat_col);
	double *mat_value_Ab = REAL(r_mat_value_Ab);
	double *mat_value_Ba = REAL(r_mat_value_Ba);

	vector<double> resid(n_obs);
	double *plhs = REAL(lhs), *pmu_in = REAL(mu_in);
	for(int obs=0 ; obs<n_obs ; ++obs){
		resid[obs] = plhs[obs] - pmu_in[obs];
	}

	//
	// const_a and const_b
	vector<double> const_a(n_i, 0);
	vector<double> const_b(n_j, 0);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;
	double *invTable_i = REAL(invTableCluster_vector);
	double *invTable_j = invTable_i + n_i;

	// for(int i=0 ; i<n_i ; ++i){
	// 	// init
	// 	const_a[i] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_a[dum_i[obs]] += resid[obs] * invTable_i[dum_i[obs]];
	// }
	//
	// for(int j=0 ; j<n_j ; ++j){
	// 	// init
	// 	const_b[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_b[dum_j[obs]] += resid[obs] * invTable_j[dum_j[obs]];
	// }

	for(int obs=0 ; obs<n_obs ; ++obs){
	    double resid_tmp = resid[obs];
	    int d_i = dum_i[obs];
	    int d_j = dum_j[obs];

	    const_a[d_i] += resid_tmp * invTable_i[d_i];

	    const_b[d_j] += resid_tmp * invTable_j[d_j];
	}


	// Rprintf("const_a: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", const_a[i]);
	// }
	// Rprintf("\n");
	//
	//
	// Rprintf("const_b: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f -- ", const_b[i]);
	// }
	// Rprintf("\n");


	// values that will be used later
	vector<double> beta(n_j);


	//
	// some tests
	//

	// double test_Ab[n_i], test_Ba[n_j];
	// for(int i=0 ; i<n_i ; ++i){
	// 	test_Ab[i] = 0;
	// }
	// for(int j=0 ; j<n_j ; ++j){
	// 	test_Ba[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_cells ; ++obs){
	// 	test_Ab[mat_row[obs]] += mat_value_Ab[obs];
	// 	test_Ba[mat_col[obs]] += mat_value_Ba[obs];
	// }

	// Rprintf("Ab x 1: \n");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f \t", test_Ab[i]);
	// }
	// Rprintf("\n");
	//
	// Rprintf("Ba x 1: \n");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f \t", test_Ba[i]);
	// }
	// Rprintf("\n");

	//
	// alpha_tilde:
	// a_tilde = const_a - (Ab %m% const_b)
	// vector<double> a_tilde(n_i);
	//
	// for(int i=0 ; i<n_i ; ++i){
	// 	a_tilde[i] = const_a[i];
	// }

	vector<double> a_tilde(const_a); // init at const_a

	for(int obs=0 ; obs<n_cells ; ++obs){
		a_tilde[mat_row[obs]] -= mat_value_Ab[obs]*const_b[mat_col[obs]];
	}

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i);
	vector<double> GX(n_i);
	vector<double> GGX(n_i);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);

	//
	// the main loop
	//

	// initialisation of X and then GX
	for(int i=0 ; i<n_i ; ++i){
		// X[i] = a_tilde[i];
		X[i] = 0;
	}

	// first iteration
	CCC_gaussian_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);

	// Rprintf("  X: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", X[i]);
	// }
	// Rprintf("\n");
	//
	//
	// Rprintf(" GX: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", GX[i]);
	// }
	// Rprintf("\n");

	// stop("kkk");

	bool numconv = false;
	bool keepGoing = true;
	int iter = 0;
	while(keepGoing && iter<iterMax){
		// Rprintf("%i.", iter);
		++iter;

		// GGX -- origin: GX, destination: GGX
		CCC_gaussian_2(GX, GGX, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);

		// Rprintf("ggx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f \t", GGX[i]);
		// }
		// Rprintf("\n");

		// X ; update of the cluster coefficient
		numconv = update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
		if(numconv) break;

		// Rprintf("  x: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f \t", X[i]);
		// }
		// Rprintf("\n");

		// GX -- origin: X, destination: GX
		CCC_gaussian_2(X, GX, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);

		// Rprintf("   gx: ");
		// for(int i=0 ; i<5 ; ++i){
		// 	Rprintf("%f\t", GX[i]);
		// }
		// Rprintf("\n");

		keepGoing = false;
		for(int i=0 ; i<n_i ; ++i){
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if(continue_criterion(X[i], GX[i], diffMax)){
			    keepGoing = true;
				break;
			}
		}

	}

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);

	// we need to compute beta, and then alpha

	//	beta = const_b - (Ba %m% alpha)
	// vector<double> beta_final(n_j);
	// for(int j=0 ; j<n_j ; ++j){
	// 	beta_final[j] = const_b[j];
	// }
	vector<double> beta_final(const_b);

	for(int obs=0 ; obs<n_cells ; ++obs){
		beta_final[mat_col[obs]] -= mat_value_Ba[obs]*GX[mat_row[obs]];
	}

	// alpha = const_a - (Ab %m% beta)
	// vector<double> alpha_final(n_i);
	// for(int i=0 ; i<n_i ; ++i){
	// 	alpha_final[i] = const_a[i];
	// }
	vector<double> alpha_final(const_a);

	for(int obs=0 ; obs<n_cells ; ++obs){
		alpha_final[mat_row[obs]] -= mat_value_Ab[obs]*beta_final[mat_col[obs]];
	}

	// mu final
	for(int obs=0 ; obs<n_obs ; ++obs){
		pmu[obs] = pmu_in[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
	}

	// Rprintf("alpha: ");
	// for(int i=0 ; i<5 ; ++i){
	// 	Rprintf("%f\t", alpha_final[i]);
	// }
	// Rprintf("\n");

	UNPROTECT(1);

	List res;
	res["mu_new"] = mu;
	res["iter"] = iter;

	return(res);
}


// [[Rcpp::export]]
List cpp_conv_seq_gau_2(int n_i, int n_j, int n_cells,
                        SEXP r_mat_row, SEXP r_mat_col, SEXP r_mat_value_Ab, SEXP r_mat_value_Ba,
                        SEXP dum_vector, SEXP lhs, SEXP invTableCluster_vector,
                        int iterMax, double diffMax, SEXP mu_in){

	//
	// Setting up
	//

	int n_obs = Rf_length(mu_in);

	int *mat_row = INTEGER(r_mat_row);
	int *mat_col = INTEGER(r_mat_col);
	double *mat_value_Ab = REAL(r_mat_value_Ab);
	double *mat_value_Ba = REAL(r_mat_value_Ba);

	vector<double> resid(n_obs);
	double *plhs = REAL(lhs), *pmu_in = REAL(mu_in);
	for(int obs=0 ; obs<n_obs ; ++obs){
		resid[obs] = plhs[obs] - pmu_in[obs];
	}

	//
	// const_a and const_b
	vector<double> const_a(n_i, 0);
	vector<double> const_b(n_j, 0);
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;
	double *invTable_i = REAL(invTableCluster_vector);
	double *invTable_j = invTable_i + n_i;

	// for(int i=0 ; i<n_i ; ++i){
	// 	// init
	// 	const_a[i] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_a[dum_i[obs]] += resid[obs] * invTable_i[dum_i[obs]];
	// }
	//
	// for(int j=0 ; j<n_j ; ++j){
	// 	// init
	// 	const_b[j] = 0;
	// }
	//
	// for(int obs=0 ; obs<n_obs ; ++obs){
	// 	const_b[dum_j[obs]] += resid[obs] * invTable_j[dum_j[obs]];
	// }

	for(int obs=0 ; obs<n_obs ; ++obs){
	    double resid_tmp = resid[obs];
	    int d_i = dum_i[obs];
	    int d_j = dum_j[obs];

	    const_a[d_i] += resid_tmp * invTable_i[d_i];

	    const_b[d_j] += resid_tmp * invTable_j[d_j];
	}


	// values that will be used later
	vector<double> beta(n_j);

	//
	// alpha_tilde:
	// a_tilde = const_a - (Ab %m% const_b)
	// vector<double> a_tilde(n_i);
	//
	// for(int i=0 ; i<n_i ; ++i){
	// 	a_tilde[i] = const_a[i];
	// }

	vector<double> a_tilde(const_a);

	for(int obs=0 ; obs<n_cells ; ++obs){
		a_tilde[mat_row[obs]] -= mat_value_Ab[obs]*const_b[mat_col[obs]];
	}

	//
	// the main loop
	//


	// X, X_new => the vector of coefficients
	vector<double> X(n_i);
	vector<double> X_new(n_i);

	// initialisation of X
	for(int i=0 ; i<n_i ; ++i){
		X[i] = a_tilde[i];
	}

	bool keepGoing = true;
	int iter = 0;
	while(keepGoing && iter<iterMax){
		++iter;

		// This way I don't need to update the values of X with a loop
		if(iter % 2 == 1){
			CCC_gaussian_2(X, X_new, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);
		} else {
			CCC_gaussian_2(X_new, X, n_i, n_j, n_cells, mat_row, mat_col, mat_value_Ab, mat_value_Ba, a_tilde, beta);
		}

		keepGoing = false;
		for(int i=0 ; i<n_i ; ++i){
			// if(fabs(X[i] - X_new[i]) / (0.1 + fabs(X_new[i])) > diffMax){
			if(continue_criterion(X[i], X_new[i], diffMax)){
			    keepGoing = true;
				break;
			}
		}

	}

	double *X_final = (iter % 2 == 1 ? X_new.data() : X.data());

	SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_obs));
	double *pmu = REAL(mu);

	// we need to compute beta, and then alpha

	//	beta = const_b - (Ba %m% alpha)
	// vector<double> beta_final(n_j);
	// for(int j=0 ; j<n_j ; ++j){
	// 	beta_final[j] = const_b[j];
	// }
	vector<double> beta_final(const_b);

	for(int obs=0 ; obs<n_cells ; ++obs){
		beta_final[mat_col[obs]] -= mat_value_Ba[obs]*X_final[mat_row[obs]];
	}

	// alpha = const_a - (Ab %m% beta)
	// vector<double> alpha_final(n_i);
	// for(int i=0 ; i<n_i ; ++i){
	// 	alpha_final[i] = const_a[i];
	// }
	vector<double> alpha_final(const_a);

	for(int obs=0 ; obs<n_cells ; ++obs){
		alpha_final[mat_row[obs]] -= mat_value_Ab[obs]*beta_final[mat_col[obs]];
	}

	// mu final
	for(int obs=0 ; obs<n_obs ; ++obs){
		pmu[obs] = pmu_in[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
	}

	UNPROTECT(1);

	List res;
	res["mu_new"] = mu;
	res["iter"] = iter;

	return(res);
}


//
// Maintenant la convergence des derivees
//



// [[Rcpp::export]]
List cpp_derivconv_seq_gnl(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all, SEXP ll_d2,
                                    SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector){

	int n_obs = Rf_length(ll_d2);
	int K = Rf_length(nb_cluster_all);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);

	int nb_coef = 0;
	for(int k=0 ; k<K ; ++k){
		nb_coef += pcluster[k];
	}

	// Setting up the vectors on variables
	vector<double*> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pjac[v] = pjac[v-1] + n_obs;
	}

	// setting up the vectors on clusters
	vector<int*> pdum(K);
	pdum[0] = INTEGER(dum_vector);

	for(int k=1 ; k<K ; ++k){
		pdum[k] = pdum[k-1] + n_obs;
	}

	// Setting up the target => deriv
	vector<double> deriv(n_obs*n_vars);
	double *my_init = REAL(deriv_init_vector);
	for(int i=0 ; i<n_obs*n_vars ; ++i){
		deriv[i] = my_init[i];
	}
	// pointers to deriv
	vector<double*> pderiv(n_vars);
	pderiv[0] = deriv.data();
	for(int v=1 ; v<n_vars ; ++v){
		pderiv[v] = pderiv[v-1] + n_obs;
	}

	// the deriv coefficients & sum_ll_d2
	vector<double> deriv_coef(nb_coef);
	vector<double> sum_ll_d2(nb_coef);

	vector<double*> pderiv_coef(K);
	vector<double*> psum_ll_d2(K);

	pderiv_coef[0] = deriv_coef.data();
	psum_ll_d2[0] = sum_ll_d2.data();
	for(int k=1 ; k<K ; ++k){
		pderiv_coef[k] = pderiv_coef[k-1] + pcluster[k-1];
		psum_ll_d2[k] = psum_ll_d2[k-1] + pcluster[k-1];
	}

	// computing the sum_ll_d2
	for(int m=0 ; m<nb_coef ; ++m){
		// init
		sum_ll_d2[m] = 0;
	}

	for(int k=0 ; k<K ; ++k){
		double *my_sum_ll_d2 = psum_ll_d2[k];
		int *my_dum = pdum[k];
		for(int i=0 ; i<n_obs ; ++i){
			my_sum_ll_d2[my_dum[i]] += pll_d2[i];
		}
	}

	//
	// Loop
	//

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for(int v=0 ; v<n_vars ; ++v){

		// loading the required data
		double *my_deriv = pderiv[v];
		double *my_jac = pjac[v];

		keepGoing = true;
		iter = 0;

		while(keepGoing && iter < iterMax){
			++iter;
			keepGoing = false;

			// we update the clusters sequentially
			// we loop over all clusters => from K to 1
			for(int k=(K-1) ; k>=0 ; k--){
				// Rprintf("k=%i ", k);
				R_CheckUserInterrupt();

				// loading the required info
				double *my_deriv_coef = pderiv_coef[k];
				int *my_dum = pdum[k];
				double *my_sum_ll_d2 = psum_ll_d2[k];
				int nb_cluster = pcluster[k];

				// init the deriv coef
				for(int m=0 ; m<nb_cluster ; ++m){
					my_deriv_coef[m] = 0;
				}

				// sum the jac and deriv
				for(int i=0 ; i<n_obs ; ++i){
					my_deriv_coef[my_dum[i]] += (my_jac[i] + my_deriv[i])*pll_d2[i];
				}

				// divide by the LL sum
				for(int m=0 ; m<nb_cluster ; ++m){
					my_deriv_coef[m] /= -my_sum_ll_d2[m];
				}

				// we update deriv
				for(int i=0 ; i<n_obs ; ++i){
					my_deriv[i] += my_deriv_coef[my_dum[i]];
				}


				// the stopping criterion
				if(keepGoing == false){
					for(int m=0 ; m<nb_cluster ; ++m){
						if(fabs(my_deriv_coef[m]) > diffMax){
							keepGoing = true;
							break;
						}
					}
				}
			}
		}

		if(iter > iter_all_max){
			iter_all_max = iter;
		}

	}

	//
	// The result
	//

	NumericMatrix dxi_dbeta(n_obs, n_vars);

	for(int v=0 ; v<n_vars ; ++v){
		double *my_deriv = pderiv[v];
		for(int i=0 ; i<n_obs ; ++i){
			dxi_dbeta(i, v) = my_deriv[i];
		}
	}

	List res;
	res["dxi_dbeta"] = dxi_dbeta;
	res["iter"] = iter_all_max;

	return(res);
}

// easier to handle in the acceleration
struct PARAM_DERIV_COEF{
	int n_obs;
	int K;

	// vectors of pointers (length K)
	vector<int*> pdum;
	vector<double*> psum_ll_d2;
	vector<double*> psum_jac_lld2;

	// simple vectors using pointers
	int *pcluster;
	double *ll_d2;

	// only value that will vary
	double *deriv_with_coef; // => vector length n_obs

};


void computeDerivCoef(vector<double*> &pcoef_origin, vector<double*> &pcoef_destination,
                      double *my_deriv_init, PARAM_DERIV_COEF *args){

	//
	// Loading the arguments
	//

	int n_obs = args->n_obs;
	int K = args->K;
	vector<int*> &pdum = args->pdum;
	vector<double*> &psum_ll_d2 = args->psum_ll_d2;
	vector<double*> &psum_jac_lld2 = args->psum_jac_lld2;
	int *pcluster = args->pcluster;
	double *ll_d2 = args->ll_d2;
	double *deriv_with_coef = args->deriv_with_coef;

	// compute all the first two coefficients
	// thus you need to start by the last item

	for(int i=0 ; i<n_obs ; ++i){
		deriv_with_coef[i] = my_deriv_init[i];
	}

	for(int k=0 ; k<(K-1) ; ++k){
		int *my_dum = pdum[k];
		double *my_deriv_coef = pcoef_origin[k];

		for(int i=0 ; i<n_obs ; ++i){
			deriv_with_coef[i] += my_deriv_coef[my_dum[i]];
		}
	}


	for(int k=K-1 ; k>=0 ; k--){
		R_CheckUserInterrupt();

		// computing the optimal cluster coef -- given mu_with_coef
		double *my_deriv_coef = pcoef_destination[k];
		double *my_sum_ll_d2 = psum_ll_d2[k];
		double *my_sum_jac_lld2 = psum_jac_lld2[k];
		int *my_dum = pdum[k];
		int nb_cluster = pcluster[k];

		//
		// update of the deriv coefficients
		//

		// init the deriv coef
		for(int m=0 ; m<nb_cluster ; ++m){
			my_deriv_coef[m] = my_sum_jac_lld2[m];
		}

		// sum the jac and deriv
		for(int i=0 ; i<n_obs ; ++i){
			my_deriv_coef[my_dum[i]] += deriv_with_coef[i]*ll_d2[i];
		}

		// divide by the LL sum
		for(int m=0 ; m<nb_cluster ; ++m){
			my_deriv_coef[m] /= -my_sum_ll_d2[m];
		}

		// updating the value of deriv_with_coef (only if necessary)
		if(k != 0){

			for(int i=0 ; i<n_obs ; ++i){
				deriv_with_coef[i] = my_deriv_init[i];
			}

			int *my_dum;
			double *my_deriv_coef;
			for(int h=0 ; h<K ; h++){
				if(h == k-1) continue;

				my_dum = pdum[h];

				if(h < k-1){
					my_deriv_coef = pcoef_origin[h];
				} else {
					my_deriv_coef = pcoef_destination[h];
				}

				for(int i=0 ; i<n_obs ; ++i){
					deriv_with_coef[i] += my_deriv_coef[my_dum[i]];
				}
			}

		}
	}

}




// [[Rcpp::export]]
List cpp_derivconv_acc_gnl(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all, SEXP ll_d2,
                                    SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector){


	int n_obs = Rf_length(ll_d2);
	int K = Rf_length(nb_cluster_all);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);

	int nb_coef = 0;
	for(int k=0 ; k<K ; ++k){
		nb_coef += pcluster[k];
	}

	// Setting up the vectors on variables
	vector<double*> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pjac[v] = pjac[v-1] + n_obs;
	}

	// setting up the vectors on clusters
	vector<int*> pdum(K);
	pdum[0] = INTEGER(dum_vector);
	for(int k=1 ; k<K ; ++k){
		pdum[k] = pdum[k-1] + n_obs;
	}

	// pointers to deriv_init_vector (length n_vars*n_obs)
	vector<double*> pderiv_init(n_vars);
	pderiv_init[0] = REAL(deriv_init_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pderiv_init[v] = pderiv_init[v-1] + n_obs;
	}

	// the deriv coefficients & sum_ll_d2 & sum_jac_lld2
	vector<double> deriv_coef(nb_coef);
	vector<double> sum_ll_d2(nb_coef);
	vector<double> sum_jac_lld2(nb_coef);
	vector<double*> pderiv_coef(K);
	vector<double*> psum_ll_d2(K);
	vector<double*> psum_jac_lld2(K);
	pderiv_coef[0] = deriv_coef.data();
	psum_ll_d2[0] = sum_ll_d2.data();
	psum_jac_lld2[0] = sum_jac_lld2.data();
	for(int k=1 ; k<K ; ++k){
		pderiv_coef[k] = pderiv_coef[k-1] + pcluster[k-1];
		psum_ll_d2[k] = psum_ll_d2[k-1] + pcluster[k-1];
		psum_jac_lld2[k] = psum_jac_lld2[k-1] + pcluster[k-1];
	}

	// computing the sum_ll_d2
	for(int m=0 ; m<nb_coef ; ++m){
		// init
		sum_ll_d2[m] = 0;
	}

	for(int k=0 ; k<K ; ++k){
		double *my_sum_ll_d2 = psum_ll_d2[k];
		int *my_dum = pdum[k];
		for(int i=0 ; i<n_obs ; ++i){
			my_sum_ll_d2[my_dum[i]] += pll_d2[i];
		}
	}

	vector<double> deriv_with_coef(n_obs);

	//
	// Sending the information
	//

	PARAM_DERIV_COEF args;
	args.n_obs = n_obs;
	args.K = K;
	args.pdum = pdum;
	args.psum_ll_d2 = psum_ll_d2;
	args.psum_jac_lld2 = psum_jac_lld2;
	args.pcluster = pcluster;
	args.ll_d2 = pll_d2;
	args.deriv_with_coef = deriv_with_coef.data();

	//
	// IT iteration (preparation)
	//

	// variables on 1:K
	vector<double> X(nb_coef);
	vector<double> GX(nb_coef);
	vector<double> GGX(nb_coef);
	// pointers:
	vector<double*> pX(K);
	vector<double*> pGX(K);
	vector<double*> pGGX(K);

	pX[0] = X.data();
	pGX[0] = GX.data();
	pGGX[0] = GGX.data();

	for(int k=1 ; k<K ; ++k){
		pX[k] = pX[k - 1] + pcluster[k - 1];
		pGX[k] = pGX[k - 1] + pcluster[k - 1];
		pGGX[k] = pGGX[k - 1] + pcluster[k - 1];
	}

	// variables on 1:(K-1)
	int nb_coef_no_K = 0;
	for(int k = 0 ; k<(K-1) ; ++k){
		nb_coef_no_K += pcluster[k];
	}
	vector<double> delta_GX(nb_coef_no_K);
	vector<double> delta2_X(nb_coef_no_K);

	//
	// Loop
	//

	NumericMatrix dxi_dbeta(n_obs, n_vars);

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for(int v=0 ; v<n_vars ; ++v){

		// loading the required data
		double *my_deriv_init = pderiv_init[v];
		double *my_jac = pjac[v];

		//
		// we update the values of sum_jac_lld2
		//

		for(int k=0 ; k<K ; ++k){
			int *my_dum = pdum[k];
			double *my_sum_jac_lld2 = psum_jac_lld2[k];

			for(int m=0 ; m<pcluster[k] ; ++m){
				my_sum_jac_lld2[m] = 0;
			}

			// sum the jac and deriv
			for(int i=0 ; i<n_obs ; ++i){
				my_sum_jac_lld2[my_dum[i]] += my_jac[i]*pll_d2[i];
			}
		}

		//
		// The IT loop
		//

		// we initialize the deriv coefficients
		for(int m=0 ; m<nb_coef ; ++m){
			X[m] = 0;
		}

		computeDerivCoef(pX, pGX, my_deriv_init, &args);

		keepGoing = true;
		iter = 0;
		bool numconv = false;

		while(keepGoing && iter < iterMax){
			++iter;

			// origin: GX, destination: GGX
			computeDerivCoef(pGX, pGGX, my_deriv_init, &args);

			// X ; update of the cluster coefficient
			numconv = update_X_IronsTuck(nb_coef_no_K, X, GX, GGX, delta_GX, delta2_X);
			if(numconv) break;

			// origin: X, destination: GX
			computeDerivCoef(pX, pGX, my_deriv_init, &args);

			keepGoing = false;
			// the stopping criterion
			for(int m=0 ; m<nb_coef_no_K ; ++m){
				// if(fabs(X[m] - GX[m]) / (0.1 + fabs(GX[m])) > diffMax){
				if(continue_criterion(X[m], GX[m], diffMax)){
				    // Rprintf("diffMax = %f\n", fabs(X[m] - GX[m]));
					keepGoing = true;
					break;
				}
			}
		}

		// save iteration info
		if(iter > iter_all_max){
			iter_all_max = iter;
		}

		//
		// Save the deriv into res
		//

		// we compute the deriv based on the cluster coefs
		for(int i=0 ; i<n_obs ; ++i){
			deriv_with_coef[i] = my_deriv_init[i];
		}

		for(int k=0 ; k<K ; ++k){
			int *my_dum = pdum[k];
			double *my_deriv_coef = pGX[k];
			for(int i=0 ; i<n_obs ; ++i){
				deriv_with_coef[i] += my_deriv_coef[my_dum[i]];
			}
		}

		// save
		for(int i=0 ; i<n_obs ; ++i){
			dxi_dbeta(i, v) = deriv_with_coef[i];
		}

	}


	List res;
	res["dxi_dbeta"] = dxi_dbeta;
	res["iter"] = iter_all_max;

	return(res);
}

void computeDerivCoef_2(vector<double> &alpha_origin, vector<double> &alpha_destination,
                        int n_i, int n_j, int n_cells,
                        const vector<double> &a_tilde,
                        const vector<int> &mat_row, const vector<int> &mat_col,
                        const vector<double> &mat_value_Ab, const vector<double> &mat_value_Ba,
                        vector<double> &beta){

	// a_tile + Ab * Ba * alpha

	for(int m=0 ; m<n_i ; ++m){
		alpha_destination[m] = a_tilde[m];
	}

	for(int m=0 ; m<n_j ; ++m){
		beta[m] = 0;
	}

	for(int obs=0 ; obs<n_cells ; ++obs){
		beta[mat_col[obs]] += mat_value_Ba[obs] * alpha_origin[mat_row[obs]];
	}

	for(int obs=0 ; obs<n_cells ; ++obs){
		alpha_destination[mat_row[obs]] += mat_value_Ab[obs] * beta[mat_col[obs]];
	}

}


// [[Rcpp::export]]
List cpp_derivconv_acc_2(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all,
                                  int n_cells, SEXP index_i, SEXP index_j, SEXP ll_d2, SEXP order,
                                  SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector){

	int n_obs = Rf_length(ll_d2);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);
	int n_i = pcluster[0], n_j = pcluster[1];

	// Setting up the vectors on variables
	vector<double*> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pjac[v] = pjac[v-1] + n_obs;
	}

	// setting up the vectors on clusters
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	// pointers to deriv_init_vector (length n_vars*n_obs)
	vector<double*> pderiv_init(n_vars);
	pderiv_init[0] = REAL(deriv_init_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pderiv_init[v] = pderiv_init[v-1] + n_obs;
	}

	// the sum_ll_d2
	vector<double> sum_ll_d2_i(n_i);
	vector<double> sum_ll_d2_j(n_j);

	// computing the sum_ll_d2
	for(int m=0 ; m<n_i ; ++m){
		sum_ll_d2_i[m] = 0;
	}

	for(int m=0 ; m<n_j ; ++m){
		sum_ll_d2_j[m] = 0;
	}

	for(int obs=0 ; obs<n_obs ; ++obs){
		sum_ll_d2_i[dum_i[obs]] += pll_d2[obs];
		sum_ll_d2_j[dum_j[obs]] += pll_d2[obs];
	}

	//
	// Setting up the matrices A and B, and vector a_tilde
	//

	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value_Ab(n_cells);
	vector<double> mat_value_Ba(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);
	int *new_order = INTEGER(order);

	int index_current = 0;

	double value_Ab = pll_d2[new_order[0]];
	double value_Ba = pll_d2[new_order[0]];

	for(int obs=1 ; obs<n_obs ; ++obs){
		if(pindex_j[obs] != pindex_j[obs-1] || pindex_i[obs] != pindex_i[obs-1]){
			// save the value if we change index
			mat_row[index_current] = pindex_i[obs-1];
			mat_col[index_current] = pindex_j[obs-1];
			mat_value_Ab[index_current] = value_Ab / -sum_ll_d2_i[pindex_i[obs-1]];
			mat_value_Ba[index_current] = value_Ba / -sum_ll_d2_j[pindex_j[obs-1]];

			// new row
			index_current++;

			// the new value
			value_Ab = pll_d2[new_order[obs]];
			value_Ba = pll_d2[new_order[obs]];
		} else {
			value_Ab += pll_d2[new_order[obs]];
			value_Ba += pll_d2[new_order[obs]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs-1];
	mat_col[index_current] = pindex_j[n_obs-1];
	mat_value_Ab[index_current] = value_Ab / -sum_ll_d2_i[pindex_i[n_obs-1]];;
	mat_value_Ba[index_current] = value_Ba / -sum_ll_d2_j[pindex_j[n_obs-1]];;

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i);
	vector<double> GX(n_i);
	vector<double> GGX(n_i);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);
	vector<double> beta(n_j);
	vector<double> alpha_final(n_i);
	vector<double> beta_final(n_j);

	//
	// Loop
	//

	NumericMatrix dxi_dbeta(n_obs, n_vars);

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for(int v=0 ; v<n_vars ; ++v){

		// loading the required data
		double *my_deriv_init = pderiv_init[v];
		double *my_jac = pjac[v];

		//
		// we update the constants and then alpha tilde
		//

		// we compute the constants
		vector<double> a(n_i);
		vector<double> b(n_j);

		for(int m=0 ; m<n_i ; ++m){
			a[m] = 0;
		}
		for(int m=0 ; m<n_j ; ++m){
			b[m] = 0;
		}
		for(int obs=0 ;obs<n_obs ; ++obs){
			a[dum_i[obs]] += (my_jac[obs]+my_deriv_init[obs]) * pll_d2[obs];
			b[dum_j[obs]] += (my_jac[obs]+my_deriv_init[obs]) * pll_d2[obs];
		}
		for(int m=0 ; m<n_i ; ++m){
			a[m] /= -sum_ll_d2_i[m];
		}
		for(int m=0 ; m<n_j ; ++m){
			b[m] /= -sum_ll_d2_j[m];
		}


		// a_tilde: a_tilde = a + (Ab %m% b)
		vector<double> a_tilde(n_i);
		for(int m=0 ; m<n_i ; ++m){
			a_tilde[m] = a[m];
		}
		for(int i=0 ; i<n_cells ; ++i){
			a_tilde[mat_row[i]] += mat_value_Ab[i] * b[mat_col[i]];
		}

		//
		// The IT loop
		//

		// we initialize the deriv coefficients
		for(int m=0 ; m<n_i ; ++m){
			X[m] = 0;
		}

		computeDerivCoef_2(X, GX, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);

		keepGoing = true;
		iter = 0;
		bool numconv = false;

		while(keepGoing && iter < iterMax){
			++iter;

			// origin: GX, destination: GGX
			computeDerivCoef_2(GX, GGX, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);

			// X ; update of the cluster coefficient
			numconv = update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
			if(numconv) break;

			// origin: X, destination: GX
			computeDerivCoef_2(X, GX, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);

			keepGoing = false;
			// the stopping criterion
			for(int m=0 ; m<n_i ; ++m){
				// if(fabs(X[m] - GX[m]) / (0.1 + fabs(GX[m])) > diffMax){
				if(continue_criterion(X[m], GX[m], diffMax)){
					keepGoing = true;
					break;
				}
			}
		}

		// save iteration info
		if(iter > iter_all_max){
			iter_all_max = iter;
		}

		//
		// Save the deriv into res
		//

		// we compute the last alpha and beta
		for(int m=0 ; m<n_i ; ++m){
			alpha_final[m] = a[m];
		}

		for(int m=0 ; m<n_j ; ++m){
			beta_final[m] = b[m];
		}

		for(int obs=0 ; obs<n_cells ; ++obs){
			beta_final[mat_col[obs]] += mat_value_Ba[obs] * GX[mat_row[obs]];
		}

		for(int obs=0 ; obs<n_cells ; ++obs){
			alpha_final[mat_row[obs]] += mat_value_Ab[obs] * beta_final[mat_col[obs]];
		}

		// save
		for(int obs=0 ; obs<n_obs ; ++obs){
			dxi_dbeta(obs, v) = my_deriv_init[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
		}

	}


	List res;
	res["dxi_dbeta"] = dxi_dbeta;
	res["iter"] = iter_all_max;

	return(res);
}


// [[Rcpp::export]]
List cpp_derivconv_seq_2(int iterMax, double diffMax, int n_vars, SEXP nb_cluster_all,
                                  int n_cells, SEXP index_i, SEXP index_j, SEXP order, SEXP ll_d2,
                                  SEXP jacob_vector, SEXP deriv_init_vector, SEXP dum_vector){

	int n_obs = Rf_length(ll_d2);

	int *pcluster = INTEGER(nb_cluster_all);
	double *pll_d2 = REAL(ll_d2);
	int n_i = pcluster[0], n_j = pcluster[1];

	// Rprintf("n_i = %i ; n_j = %i ; n_obs = %i ; n_cells = %i\n", n_i, n_j, n_obs, n_cells);

	// Setting up the vectors on variables
	vector<double*> pjac(n_vars);
	pjac[0] = REAL(jacob_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pjac[v] = pjac[v-1] + n_obs;
	}

	// setting up the vectors on clusters
	int *dum_i = INTEGER(dum_vector);
	int *dum_j = dum_i + n_obs;

	// pointers to deriv_init_vector (length n_vars*n_obs)
	vector<double*> pderiv_init(n_vars);
	pderiv_init[0] = REAL(deriv_init_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pderiv_init[v] = pderiv_init[v-1] + n_obs;
	}

	// the sum_ll_d2
	vector<double> sum_ll_d2_i(n_i);
	vector<double> sum_ll_d2_j(n_j);


	// computing the sum_ll_d2
	for(int m=0 ; m<n_i ; ++m){
		sum_ll_d2_i[m] = 0;
	}

	for(int m=0 ; m<n_j ; ++m){
		sum_ll_d2_j[m] = 0;
	}

	for(int obs=0 ; obs<n_obs ; ++obs){
		sum_ll_d2_i[dum_i[obs]] += pll_d2[obs];
		sum_ll_d2_j[dum_j[obs]] += pll_d2[obs];
	}

	//
	// Setting up the matrices A and B, and vector a_tilde
	//

	vector<int> mat_row(n_cells);
	vector<int> mat_col(n_cells);
	vector<double> mat_value_Ab(n_cells);
	vector<double> mat_value_Ba(n_cells);

	int *pindex_i = INTEGER(index_i);
	int *pindex_j = INTEGER(index_j);
	int *new_order = INTEGER(order);

	int index_current = 0;

	double value_current = pll_d2[new_order[0]];

	for(int obs=1 ; obs<n_obs ; ++obs){
		if(pindex_j[obs] != pindex_j[obs-1] || pindex_i[obs] != pindex_i[obs-1]){
			// save the value if we change index
			mat_row[index_current] = pindex_i[obs-1];
			mat_col[index_current] = pindex_j[obs-1];
			mat_value_Ab[index_current] = value_current / -sum_ll_d2_i[pindex_i[obs-1]];
			mat_value_Ba[index_current] = value_current / -sum_ll_d2_j[pindex_j[obs-1]];

			// new row
			index_current++;

			// the new value
			value_current = pll_d2[new_order[obs]];
		} else {
			value_current += pll_d2[new_order[obs]];
		}
	}

	// last save
	mat_row[index_current] = pindex_i[n_obs-1];
	mat_col[index_current] = pindex_j[n_obs-1];
	mat_value_Ab[index_current] = value_current / -sum_ll_d2_i[pindex_i[n_obs-1]];
	mat_value_Ba[index_current] = value_current / -sum_ll_d2_j[pindex_j[n_obs-1]];

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i);
	vector<double> X_new(n_i);
	double *X_final;
	vector<double> beta(n_j);
	vector<double> alpha_final(n_i);
	vector<double> beta_final(n_j);

	//
	// Loop
	//

	NumericMatrix dxi_dbeta(n_obs, n_vars);

	bool keepGoing = true;
	int iter = 0, iter_all_max = 0;

	// We first loop on each variable
	for(int v=0 ; v<n_vars ; ++v){

		// loading the required data
		double *my_deriv_init = pderiv_init[v];
		double *my_jac = pjac[v];

		//
		// we update the constants and then alpha tilde
		//

		// we compute the constants
		vector<double> a(n_i);
		vector<double> b(n_j);
		for(int m=0 ; m<n_i ; ++m){
			a[m] = 0;
		}
		for(int m=0 ; m<n_j ; ++m){
			b[m] = 0;
		}
		for(int obs=0 ;obs<n_obs ; ++obs){
			a[dum_i[obs]] += (my_jac[obs]+my_deriv_init[obs]) * pll_d2[obs];
			b[dum_j[obs]] += (my_jac[obs]+my_deriv_init[obs]) * pll_d2[obs];
		}
		for(int m=0 ; m<n_i ; ++m){
			a[m] /= -sum_ll_d2_i[m];
		}
		for(int m=0 ; m<n_j ; ++m){
			b[m] /= -sum_ll_d2_j[m];
		}


		// a_tilde: a_tilde = a + (Ab %m% b)
		vector<double> a_tilde(n_i);
		for(int m=0 ; m<n_i ; ++m){
			a_tilde[m] = a[m];
		}
		for(int i=0 ; i<n_cells ; ++i){
			a_tilde[mat_row[i]] += mat_value_Ab[i] * b[mat_col[i]];
		}

		// init of X (alpha)
		for(int m=0 ; m<n_i ; ++m){
			X[m] = 0;
		}

		keepGoing = true;
		iter = 0;

		while(keepGoing && iter < iterMax){
			++iter;

			if(iter % 2 == 1){
				computeDerivCoef_2(X, X_new, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);
			} else {
				computeDerivCoef_2(X_new, X, n_i, n_j, n_cells, a_tilde, mat_row, mat_col, mat_value_Ab, mat_value_Ba, beta);
			}


			keepGoing = false;
			// the stopping criterion
			for(int m=0 ; m<n_i ; ++m){
				// if(fabs(X[m] - X_new[m]) / (0.1 + fabs(X_new[m])) > diffMax){
				if(continue_criterion(X[m], X_new[m], diffMax)){
				    keepGoing = true;
					break;
				}
			}
		}

		// save iteration info
		if(iter > iter_all_max){
			iter_all_max = iter;
		}

		//
		// Save the deriv into res
		//


		X_final = (iter % 2 == 1) ? X_new.data() : X.data();

		// we compute the last alpha and beta
		for(int m=0 ; m<n_i ; ++m){
			alpha_final[m] = a[m];
		}

		for(int m=0 ; m<n_j ; ++m){
			beta_final[m] = b[m];
		}

		for(int obs=0 ; obs<n_cells ; ++obs){
			beta_final[mat_col[obs]] += mat_value_Ba[obs] * X_final[mat_row[obs]];
		}

		for(int obs=0 ; obs<n_cells ; ++obs){
			alpha_final[mat_row[obs]] += mat_value_Ab[obs] * beta_final[mat_col[obs]];
		}

		// save
		for(int obs=0 ; obs<n_obs ; ++obs){
			dxi_dbeta(obs, v) = my_deriv_init[obs] + alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
		}

	}


	List res;
	res["dxi_dbeta"] = dxi_dbeta;
	res["iter"] = iter_all_max;

	return(res);
}

// [[Rcpp::export]]
NumericMatrix update_deriv_single(int n_vars, int nb_coef,
                         SEXP r_ll_d2, SEXP r_jacob_vector, SEXP r_dum_vector){

	int n_obs = Rf_length(r_ll_d2);

	// loading variables
	double *ll_d2 = REAL(r_ll_d2);
	int *dum = INTEGER(r_dum_vector);

	vector<double*> pjac(n_vars);
	pjac[0] = REAL(r_jacob_vector);
	for(int v=1 ; v<n_vars ; ++v){
		pjac[v] = pjac[v-1] + n_obs;
	}

	// vector sum_ll_d2
	vector<double> sum_ll_d2(nb_coef, 0);
	for(int obs=0 ; obs<n_obs ; ++obs){
		sum_ll_d2[dum[obs]] += ll_d2[obs];
	}

	// the vector of coefficients
	vector<double> coef_deriv(nb_coef);

	// the result
	NumericMatrix res(n_obs, n_vars); // init at 0

	for(int v=0 ; v<n_vars ; ++v){
		double *my_jac = pjac[v];

		//
		// 1) we compute the coef
		//

		for(int m=0 ; m<nb_coef ; ++m){
			coef_deriv[m] = 0;
		}

		for(int obs=0 ; obs<n_obs ; ++obs){
			coef_deriv[dum[obs]] += my_jac[obs]*ll_d2[obs];
		}

		// divide is costly, we do it only nb_coef times
		for(int m=0 ; m<nb_coef ; ++m){
			coef_deriv[m] /= -sum_ll_d2[m];
		}

		//
		// 2) We save the value in res
		//

		for(int obs=0 ; obs<n_obs ; ++obs){
			res(obs, v) = coef_deriv[dum[obs]];
		}

	}


	return(res);
}




