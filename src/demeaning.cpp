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


// Stopping / continuing criteria:
// Functions used inside all loops
inline bool continue_crit(double a, double b, double diffMax){
    // continuing criterion of the algorithm
    double diff = fabs(a - b);
    return ( (diff > diffMax) && (diff/(0.1 + fabs(a)) > diffMax) );
}

inline bool stopping_crit(double a, double b, double diffMax){
    // continuing criterion of the algorithm
    double diff = fabs(a - b);
    return ( (diff < diffMax) || (diff/(0.1 + fabs(a)) < diffMax) );
}

//
// Demeans each variable in input
// The method is based on obtaining the optimal cluster coefficients
//

// for interruption
int pending_interrupt() {
	return !(R_ToplevelExec(Rcpp::checkInterruptFn, NULL));
}
// Works but only the master thread can call that function!
// what happens if the master thread has finished its job but the lower thread is in an "infinite" loop?
// this is tricky, as such we cannot stop it
// solution: create a function keeping the threads idle waiting for the complete job to be done
// BUT I need to add static allocation of threads => performance cost


// List of objects, used to
// lighten the writting of the functions
struct PARAM_DEMEAN{
	int n_obs;
	int Q;
	int nb_coef;
	int iterMax;
	double diffMax;

	// number of clusters & iterations
	int *pcluster;
	int *piterations_all;

	// weights + slopes:
	bool isWeight;
	// double *obs_weights;
	vector<double*> psum_weights;
	vector<double*> all_obs_weights;

	bool isSlope;
	int *slope_flag;
	vector<double*> all_slope_vars;

	// vectors of pointers
	vector<int*> pdum;
	vector<int*> ptable; // to remove
	vector<double*> pinput;
	vector<double*> poutput;

	// saving the fixed effects
	bool save_fixef;
	double *fixef_values;

	// soptflag
	bool *stopnow;
	int *jobdone;
};

bool dm_update_X_IronsTuck(int nb_coef_no_Q, vector<double> &X,
                        const vector<double> &GX, const vector<double> &GGX,
                        vector<double> &delta_GX, vector<double> &delta2_X){

	for(int i=0 ; i<nb_coef_no_Q ; ++i){
	    double GX_tmp = GX[i];
		delta_GX[i] = GGX[i] - GX_tmp;
		delta2_X[i] = delta_GX[i] - GX_tmp + X[i];
	}

	double vprod = 0, ssq = 0;
	for(int i=0 ; i<nb_coef_no_Q ; ++i){
	    double delta2_X_tmp = delta2_X[i];
		vprod += delta_GX[i] * delta2_X_tmp;
		ssq += delta2_X_tmp * delta2_X_tmp;
	}

	bool res = false;

	if(ssq == 0){
		res = true;
	} else {
		double coef = vprod/ssq;

		// update of X:
		for(int i=0 ; i<nb_coef_no_Q ; ++i){
			X[i] = GGX[i] - coef * delta_GX[i];
		}
	}

	return(res);
}



void demean_single_1(int v, PARAM_DEMEAN* args){
	// v: variable identifier to demean

	// Q == 1: nothing to say, just compute the closed form

	// loading the data
	int n_obs = args->n_obs;
	int nb_coef = args->nb_coef;
	vector<int*> &pdum = args->pdum;
	vector<double*> &psum_weights = args->psum_weights;
	vector<double*> &pinput = args->pinput;
	vector<double*> &poutput = args->poutput;

	// vector of cluster coefficients initialized at 0
	vector<double> cluster_coef(nb_coef, 0);

	// interruption handling
	bool isMaster = omp_get_thread_num() == 0;
	bool *pStopNow = args->stopnow;
	if(isMaster){
		if(pending_interrupt()){
			*pStopNow = true;
		}
	}

	// dum and table
	int *dum = pdum[0];
	double *sum_weights = psum_weights[0];

	// weights
	bool isWeight = args->isWeight;
	// double *obs_weights = args->obs_weights;
	double *obs_weights = args->all_obs_weights[0];

	// slopes
	bool isSlope = args->isSlope;

	// the input & output
	double *input = pinput[v];
	double *output = poutput[v];

	// The sum
	if(isWeight){
		for(int obs=0 ; obs<n_obs ; ++obs){
			cluster_coef[dum[obs]] += obs_weights[obs] * input[obs];
		}
	} else {
		for(int obs=0 ; obs<n_obs ; ++obs){
			cluster_coef[dum[obs]] += input[obs];
		}
	}

	// calculating cluster coef
	for(int m=0 ; m<nb_coef ; ++m){
		cluster_coef[m] /= sum_weights[m];
	}

	// Output:
	if(isSlope){
	    double *my_slope_var = args->all_slope_vars[0];
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        output[obs] = my_slope_var[obs] * cluster_coef[dum[obs]];
	    }
	} else {
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        output[obs] = cluster_coef[dum[obs]];
	    }
	}

	// saving the fixef coefs
	double *fixef_values = args->fixef_values;
	if(args->save_fixef){
	    for(int m=0 ; m<nb_coef ; ++m){
	        fixef_values[m] = cluster_coef[m];
	    }
	}


}

void CCC_gaussian_2(const vector<double> &pcluster_origin, vector<double> &pcluster_destination,
                    int n_i, int n_j,
                    int n_obs, int *dum_i, int *dum_j,
                    bool isWeight, double *obs_weights_i, double *obs_weights_j,
                    bool isSlope_i, double *slope_var_i, bool isSlope_j, double *slope_var_j,
                    double *sum_weights_i, double *sum_weights_j,
                    const vector<double> &a_tilde, vector<double> &beta){

	// alpha = a_tilde + (Ab %m% (Ba %m% alpha))
	for(int i=0 ; i<n_i ; ++i){
		pcluster_destination[i] = 0;
	}

	for(int j=0 ; j<n_j ; ++j){
		beta[j] = 0;
	}

	if(isSlope_i){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        beta[dum_j[obs]] += obs_weights_j[obs] * slope_var_i[obs] * pcluster_origin[dum_i[obs]];
	    }
	} else if(isWeight){
		for(int obs=0 ; obs<n_obs ; ++obs){
			beta[dum_j[obs]] += obs_weights_j[obs] * pcluster_origin[dum_i[obs]];
		}
	} else {
		for(int obs=0 ; obs<n_obs ; ++obs){
			beta[dum_j[obs]] += pcluster_origin[dum_i[obs]];
		}
	}

	for(int j=0 ; j<n_j ; ++j){
		beta[j] /= sum_weights_j[j];
	}


	if(isSlope_j){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        pcluster_destination[dum_i[obs]] += obs_weights_i[obs] * slope_var_j[obs] * beta[dum_j[obs]];
	    }
	} else if(isWeight){
		for(int obs=0 ; obs<n_obs ; ++obs){
			pcluster_destination[dum_i[obs]] += obs_weights_i[obs] * beta[dum_j[obs]];
		}
	} else {
		for(int obs=0 ; obs<n_obs ; ++obs){
			pcluster_destination[dum_i[obs]] += beta[dum_j[obs]];
		}
	}

	for(int i=0 ; i<n_i ; ++i){
		pcluster_destination[i] /= sum_weights_i[i];
		pcluster_destination[i] += a_tilde[i];
	}

	// Rprintf("pcluster_destination: %.3f, %.3f, %.3f, %.3f\n", pcluster_destination[0], pcluster_destination[1], pcluster_destination[2], pcluster_destination[3]);

}

void demean_acc_2(int v, int iterMax, PARAM_DEMEAN *args){

	//
	// Setting up
	//

	// loading data
	int n_obs = args->n_obs;
	int *nb_cluster = args->pcluster;
	int n_i = nb_cluster[0], n_j = nb_cluster[1];
	double diffMax = args->diffMax;

	// input output
	vector<double*> &pinput = args->pinput;
	vector<double*> &poutput = args->poutput;
	double *input = pinput[v];
	double *output = poutput[v];

	vector<double*> &psum_weights = args->psum_weights;
	double *sum_weights_i = psum_weights[0];
	double *sum_weights_j = psum_weights[1];

	bool isWeight = args->isWeight;
	double *obs_weights_i = args->all_obs_weights[0];
	double *obs_weights_j = args->all_obs_weights[1];

	bool isSlope = args->isSlope;
	bool isSlope_i = args->slope_flag[0];
	bool isSlope_j = args->slope_flag[1];

	double *slope_var_i = args->all_slope_vars[0];
	double *slope_var_j = args->all_slope_vars[1];

	int *dum_i = args->pdum[0];
	int *dum_j = args->pdum[1];

	int *iterations_all = args->piterations_all;

	//
	// const_a and const_b
	vector<double> const_a(n_i, 0);
	vector<double> const_b(n_j, 0);

	if(isSlope){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        double resid_tmp = input[obs] - output[obs];
	        const_a[dum_i[obs]] += obs_weights_i[obs]*resid_tmp;
	        const_b[dum_j[obs]] += obs_weights_j[obs]*resid_tmp;
	    }
	} else if(isWeight){
	    // weights are identical if there is no slope
        for(int obs=0 ; obs<n_obs ; ++obs){
            double resid_tmp = obs_weights_i[obs]*(input[obs] - output[obs]);
            const_a[dum_i[obs]] += resid_tmp;
            const_b[dum_j[obs]] += resid_tmp;
        }
    } else {
        for(int obs=0 ; obs<n_obs ; ++obs){
            double resid_tmp = input[obs] - output[obs];
	        const_a[dum_i[obs]] += resid_tmp;
	        const_b[dum_j[obs]] += resid_tmp;
        }
    }


	for(int i=0 ; i<n_i ; ++i){
		const_a[i] /= sum_weights_i[i];
	}

	for(int j=0 ; j<n_j ; ++j){
		const_b[j] /= sum_weights_j[j];
	}

	// Rprintf("const_a: %.3f, %.3f, %.3f, %.3f\n", const_a[0], const_a[1], const_a[2], const_a[3]);
	// Rprintf("const_b: %.3f, %.3f, %.3f, %.3f\n", const_b[0], const_b[1], const_b[2], const_b[3]);


	// values that will be used later
	vector<double> beta(n_j);

	// alpha_tilde
	vector<double> a_tilde(n_i, 0);

	if(isSlope_j){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        a_tilde[dum_i[obs]] -= obs_weights_i[obs] * slope_var_j[obs] * const_b[dum_j[obs]];
	    }
	} else if(isWeight){
		for(int obs=0 ; obs<n_obs ; ++obs){
			a_tilde[dum_i[obs]] -= obs_weights_i[obs] * const_b[dum_j[obs]];
		}
	} else {
		for(int obs=0 ; obs<n_obs ; ++obs){
			a_tilde[dum_i[obs]] -= const_b[dum_j[obs]];
		}
	}


	for(int i=0 ; i<n_i ; ++i){
		a_tilde[i] /= sum_weights_i[i];
		a_tilde[i] += const_a[i];
	}

	// Rprintf("a_tilde: %.3f, %.3f, %.3f, %.3f\n", a_tilde[0], a_tilde[1], a_tilde[2], a_tilde[3]);

	// interruption handling
	bool isMaster = omp_get_thread_num() == 0;
	bool *pStopNow = args->stopnow;
	double flop = 20 * n_obs; // rough estimate nber operation per iter
	int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second

	//
	// IT iteration (preparation)
	//

	vector<double> X(n_i, 0);
	vector<double> GX(n_i);
	vector<double> GGX(n_i);
	vector<double> delta_GX(n_i);
	vector<double> delta2_X(n_i);

	//
	// the main loop
	//

	// first iteration => update GX
	CCC_gaussian_2(X, GX, n_i, n_j, n_obs, dum_i, dum_j, isWeight, obs_weights_i, obs_weights_j,
                isSlope_i, slope_var_i, isSlope_j, slope_var_j,
                sum_weights_i, sum_weights_j, a_tilde, beta);

	// For the stopping criterion on total addition
	vector<double> mu_last(n_obs, 0);
	double input_mean = 0;

	bool numconv = false;
	bool keepGoing = true;
	int iter = 1;
	while(!*pStopNow && keepGoing && iter<=iterMax){

		if(isMaster && iter % iterSecond == 0){
			if(pending_interrupt()){
				*pStopNow = true;
				break;
			}
		}

		++iter;

		// GGX -- origin: GX, destination: GGX
		CCC_gaussian_2(GX, GGX, n_i, n_j, n_obs, dum_i, dum_j, isWeight, obs_weights_i, obs_weights_j,
                 isSlope_i, slope_var_i, isSlope_j, slope_var_j,
                 sum_weights_i, sum_weights_j, a_tilde, beta);

		// X ; update of the cluster coefficient
		numconv = dm_update_X_IronsTuck(n_i, X, GX, GGX, delta_GX, delta2_X);
		if(numconv) break;

		// GX -- origin: X, destination: GX
		CCC_gaussian_2(X, GX, n_i, n_j, n_obs, dum_i, dum_j, isWeight, obs_weights_i, obs_weights_j,
                 isSlope_i, slope_var_i, isSlope_j, slope_var_j,
                 sum_weights_i, sum_weights_j, a_tilde, beta);

		keepGoing = false;
		for(int i=0 ; i<n_i ; ++i){
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if(continue_crit(X[i], GX[i], diffMax)){
				keepGoing = true;
				break;
			}
		}

		// Other stopping criterion: change to output very small
		if(iter % 50 == 0){

			// we need to compute beta
			for(int j=0 ; j<n_j ; ++j){
				beta[j] = 0;
			}

			if(isSlope_i){
			    for(int obs=0 ; obs<n_obs ; ++obs){
			        beta[dum_j[obs]] += obs_weights_j[obs] * slope_var_i[obs] * GX[dum_i[obs]];
			    }
			} else if(isWeight){
				for(int obs=0 ; obs<n_obs ; ++obs){
					beta[dum_j[obs]] += obs_weights_j[obs] * GX[dum_i[obs]];
				}
			} else {
				for(int obs=0 ; obs<n_obs ; ++obs){
					beta[dum_j[obs]] += GX[dum_i[obs]];
				}
			}

			for(int j=0 ; j<n_j ; ++j){
				beta[j] /= sum_weights_j[j];
				beta[j] += const_b[j];
			}

			vector<double> mu_current(n_obs);
			if(isSlope){
			    for(int obs=0 ; obs<n_obs ; ++obs){
			        mu_current[obs] = slope_var_i[obs] * GX[dum_i[obs]] + slope_var_j[obs] * beta[dum_j[obs]];
			    }
			} else {
			    for(int obs=0 ; obs<n_obs ; ++obs){
			        mu_current[obs] = GX[dum_i[obs]] + beta[dum_j[obs]];
			    }
			}


			// init mu_last if iter == 50 / otherwise, comparison
			if(iter == 50){
				for(int i=0 ; i<n_obs ; ++i){
					mu_last[i] = mu_current[i];
				}

				for(int i=0 ; i<n_obs ; ++i){
					input_mean += input[i];
				}
				input_mean /= n_obs;
				input_mean = fabs(input_mean);

			} else {
				double mu_diff = 0;
				for(int i=0 ; i<n_obs ; ++i){
					mu_diff += fabs(mu_last[i] - mu_current[i]);
				}
				mu_diff = mu_diff / n_obs; // we scale it to have a meaningful value
				// Rprintf("iter %i -- mu_diff = %f (mean = %f)\n", iter, mu_diff, input_mean);

				// if(mu_diff < diffMax || mu_diff / (0.1 + input_mean) < diffMax){
				if(stopping_crit(input_mean, mu_diff + input_mean, diffMax)){
					// Rprintf("mu_diff BREAK\n");
					break;
				}

				// we update mu_last
				for(int i=0 ; i<n_obs ; ++i){
					mu_last[i] = mu_current[i];
				}

			}
		}

	}

	//
	// we update the result (output)
	//

	// we need to compute beta, and then alpha
	vector<double> beta_final(n_j, 0);

	if(isSlope_i){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        beta_final[dum_j[obs]] -= obs_weights_j[obs] * slope_var_i[obs] * GX[dum_i[obs]];
	    }
	} else if(isWeight){
		for(int obs=0 ; obs<n_obs ; ++obs){
			beta_final[dum_j[obs]] -= obs_weights_j[obs]*GX[dum_i[obs]];
		}
	} else {
		for(int obs=0 ; obs<n_obs ; ++obs){
			beta_final[dum_j[obs]] -= GX[dum_i[obs]];
		}
	}

	for(int j=0 ; j<n_j ; ++j){
		beta_final[j] /= sum_weights_j[j];
		beta_final[j] += const_b[j];
	}

	// alpha = const_a - (Ab %m% beta)
	vector<double> alpha_final(n_i, 0);

	if(isSlope_j){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        alpha_final[dum_i[obs]] -= obs_weights_i[obs] * slope_var_j[obs] * beta_final[dum_j[obs]];
	    }
	} else if(isWeight){
		for(int obs=0 ; obs<n_obs ; ++obs){
			alpha_final[dum_i[obs]] -= obs_weights_i[obs]*beta_final[dum_j[obs]];
		}
	} else {
		for(int obs=0 ; obs<n_obs ; ++obs){
			alpha_final[dum_i[obs]] -= beta_final[dum_j[obs]];
		}
	}

	for(int i=0 ; i<n_i ; ++i){
		alpha_final[i] /= sum_weights_i[i];
		alpha_final[i] += const_a[i];
	}

	// mu final
	if(isSlope){
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        output[obs] += slope_var_i[obs] * alpha_final[dum_i[obs]] + slope_var_j[obs] * beta_final[dum_j[obs]];
	    }
	} else {
	    for(int obs=0 ; obs<n_obs ; ++obs){
	        output[obs] += alpha_final[dum_i[obs]] + beta_final[dum_j[obs]];
	    }
	}

	// keeping track of iterations
	iterations_all[v] += iter;

	// saving the fixef coefs
	double *fixef_values = args->fixef_values;
	if(args->save_fixef){

	    for(int i=0 ; i<n_i ; ++i){
	        fixef_values[i] += alpha_final[i];
	    }

	    for(int j=0 ; j<+n_j ; ++j){
	        fixef_values[n_i + j] += beta_final[j];
	    }
	}

}

void compute_mean(int n_obs, int nb_cluster, double *cluster_coef, const vector<double> &sum_other_means,
                  double *sum_in_out, int *dum, double *sum_weights){

    // NOTA: sum_in_out and sum_other_means are already "weighted" (see in computeMeans and in demean_acc_gnl)

	// initialize cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = 0;
	}

	// looping sequentially over the sum of other coefficients
	for(int i=0 ; i<n_obs ; ++i){
		cluster_coef[dum[i]] += sum_other_means[i];
	}

	// calculating cluster coef
	for(int m=0 ; m<nb_cluster ; ++m){
		cluster_coef[m] = (sum_in_out[m] - cluster_coef[m]) / sum_weights[m];
	}

	// "output" is the update of cluster_coef
}

void computeMeans(vector<double*> &pcluster_origin, vector<double*> &pcluster_destination,
                  vector<double> &sum_other_means, vector<double*> &psum_input_output, PARAM_DEMEAN *args){
	// update of the cluster coefficients
	// first we update mu, then we update the cluster coefficicents

	//
	// Loading the variables
	//

	int n_obs = args->n_obs;
	int Q = args->Q;

	int *pcluster = args->pcluster;

	vector<int*> &pdum = args->pdum;

	// weights:
	vector<double*> &psum_weights = args->psum_weights;
	bool isWeight = args->isWeight;
	vector<double*> &all_obs_weights = args->all_obs_weights;

	// slopes:
	int *slope_flag = args->slope_flag;
	vector<double*> &all_slope_vars = args->all_slope_vars;

	// We update each cluster coefficient, starting from Q (the smallest one)

	for(int i=0 ; i<n_obs ; ++i){
		sum_other_means[i] = 0;
	}

    // we start with Q-1
	double *obs_weights_current = all_obs_weights[Q-1];
	for(int q=0 ; q<(Q-1) ; ++q){
		int *my_dum = pdum[q];
		double *my_cluster_coef = pcluster_origin[q];

		if(slope_flag[q]){
		    double *my_slope_var = all_slope_vars[q];
		    for(int obs=0 ; obs<n_obs ; ++obs){
		        sum_other_means[obs] += obs_weights_current[obs] * my_slope_var[obs] * my_cluster_coef[my_dum[obs]];
		    }
		} else if(isWeight){
			for(int obs=0 ; obs<n_obs ; ++obs){
				sum_other_means[obs] += obs_weights_current[obs] * my_cluster_coef[my_dum[obs]];
			}
		} else {
			for(int obs=0 ; obs<n_obs ; ++obs){
				sum_other_means[obs] += my_cluster_coef[my_dum[obs]];
			}
		}
	}


	for(int q=Q-1 ; q>=0 ; q--){

		// computing the optimal cluster coef -- given mu_with_coef
		double *my_cluster_coef = pcluster_destination[q];
		double *my_sum_weights = psum_weights[q];
		double *my_sum_in_out = psum_input_output[q];
		int *my_dum = pdum[q];
		int nb_cluster = pcluster[q];

		// update of the cluster coefficients
		compute_mean(n_obs, nb_cluster, my_cluster_coef, sum_other_means, my_sum_in_out, my_dum, my_sum_weights);

		// if(int q == 0){
		// 	Rprintf("pcluster_destination: %.3f, %.3f, %.3f, %.3f\n", my_cluster_coef[0], my_cluster_coef[1], my_cluster_coef[2], my_cluster_coef[3]);
		// }


		// updating the value of sum_other_means (only if necessary)
		if(q != 0){

			// We recompute it from scratch (only way -- otherwise precision problems arise)

			for(int i=0 ; i<n_obs ; ++i){
				sum_other_means[i] = 0;
			}

			int *my_dum;
			double *my_cluster_coef;
			double *obs_weights_current = all_obs_weights[q-1];
			for(int h=0 ; h<Q ; h++){
				if(h == q-1) continue;

				my_dum = pdum[h];

				if(h < q-1){
					my_cluster_coef = pcluster_origin[h];
				} else {
					my_cluster_coef = pcluster_destination[h];
				}

				if(slope_flag[h]){
				    double *my_slope_var = all_slope_vars[h];
				    for(int obs=0 ; obs<n_obs ; ++obs){
				        sum_other_means[obs] += obs_weights_current[obs] * my_slope_var[obs] * my_cluster_coef[my_dum[obs]];
				    }
				} else if(isWeight){
					for(int obs=0 ; obs<n_obs ; ++obs){
						sum_other_means[obs] += obs_weights_current[obs] * my_cluster_coef[my_dum[obs]];
					}
				} else {
					for(int obs=0 ; obs<n_obs ; ++obs){
						sum_other_means[obs] += my_cluster_coef[my_dum[obs]];
					}
				}
			}

		}
	}

	// In the end, the array pcluster_destination is fully updated, starting from Q to 1

}

bool demean_acc_gnl(int v, int iterMax, PARAM_DEMEAN *args){

	//
	// data
	//

	int n_obs = args->n_obs;
	int nb_coef = args->nb_coef;
	int Q = args->Q;
	double diffMax = args->diffMax;

	int *pcluster = args->pcluster;

	vector<int*> pdum = args->pdum;

	// weights
	bool isWeight = args->isWeight;
	// double *obs_weights = args->obs_weights;
	vector<double*> &all_obs_weights = args->all_obs_weights;

	// slopes
	int *slope_flag = args->slope_flag;
	vector<double*> &all_slope_vars = args->all_slope_vars;

	// input output
	vector<double*> &pinput = args->pinput;
	vector<double*> &poutput = args->poutput;
	double *input = pinput[v];
	double *output = poutput[v];

	// temp var:
	vector<double> sum_other_means(n_obs);

	// conditional sum of input minus output
	vector<double> sum_input_output(nb_coef, 0);
	vector<double*> psum_input_output(Q);
	psum_input_output[0] = sum_input_output.data();
	for(int q=1 ; q<Q ; ++q){
		psum_input_output[q] = psum_input_output[q - 1] + pcluster[q - 1];
	}


	for(int q=0 ; q<Q ; ++q){
		double *my_sum_input_output = psum_input_output[q];
		int *my_dum = pdum[q];

		if(isWeight){
		    double *obs_weights_current = all_obs_weights[q];
			for(int obs=0 ; obs<n_obs ; ++obs){
				my_sum_input_output[my_dum[obs]] += obs_weights_current[obs] * (input[obs] - output[obs]);
			}
		} else {
			for(int obs=0 ; obs<n_obs ; ++obs){
				my_sum_input_output[my_dum[obs]] += input[obs] - output[obs];
			}
		}
	}

	// interruption handling
	bool isMaster = omp_get_thread_num() == 0;
	bool *pStopNow = args->stopnow;
	double flop = 4*(5 + 12*(Q-1) + 4*(Q-1)*(Q-1))*n_obs; // rough estimate nber operation per iter
	int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second

	//
	// IT iteration (preparation)
	//

	// variables on 1:K
	vector<double> X(nb_coef, 0);
	vector<double> GX(nb_coef);
	vector<double> GGX(nb_coef);
	// pointers:
	vector<double*> pX(Q);
	vector<double*> pGX(Q);
	vector<double*> pGGX(Q);
	pX[0] = X.data();
	pGX[0] = GX.data();
	pGGX[0] = GGX.data();

	for(int q=1 ; q<Q ; ++q){
		pX[q] = pX[q - 1] + pcluster[q - 1];
		pGX[q] = pGX[q - 1] + pcluster[q - 1];
		pGGX[q] = pGGX[q - 1] + pcluster[q - 1];
	}

	// variables on 1:(Q-1)
	int nb_coef_no_Q = 0;
	for(int q = 0 ; q<(Q-1) ; ++q){
		nb_coef_no_Q += pcluster[q];
	}
	vector<double> delta_GX(nb_coef_no_Q);
	vector<double> delta2_X(nb_coef_no_Q);

	//
	// the main loop
	//

	// first iteration
	computeMeans(pX, pGX, sum_other_means, psum_input_output, args);

	// check whether we should go into the loop
	bool keepGoing = false;
	for(int i=0 ; i<nb_coef ; ++i){
		// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
		if(continue_crit(X[i], GX[i], diffMax)){
			keepGoing = true;
			break;
		}
	}

	// For the stopping criterion on total addition
	vector<double> mu_last(n_obs, 0);
	double input_mean = 0;

	int iter = 0;
	bool numconv = false;
	while(!*pStopNow && keepGoing && iter<iterMax){

		if(isMaster && iter % iterSecond == 0){
			if(pending_interrupt()){
				*pStopNow = true;
				break;
			}
		}

		iter++;

		// GGX -- origin: GX, destination: GGX
		computeMeans(pGX, pGGX, sum_other_means, psum_input_output, args);

		// X ; update of the cluster coefficient
		numconv = dm_update_X_IronsTuck(nb_coef_no_Q, X, GX, GGX, delta_GX, delta2_X);
		if(numconv) break;

		// GX -- origin: X, destination: GX
		computeMeans(pX, pGX, sum_other_means, psum_input_output, args);

		keepGoing = false;
		for(int i=0 ; i<nb_coef_no_Q ; ++i){
			// diffmax: nber of significant digits
			// if(fabs(X[i] - GX[i]) / (0.1 + fabs(GX[i])) > diffMax){
			if(continue_crit(X[i], GX[i], diffMax)){
				keepGoing = true;
				break;
			}
		}

		// Other stopping criterion: change to output very small
		if(iter % 50 == 0){

			vector<double> mu_current(n_obs, 0);
			for(int q=0 ; q<Q ; ++q){
				int *my_dum = pdum[q];
				double *my_cluster_coef = pGX[q];

				if(slope_flag[q]){
				    double *my_slope_var = all_slope_vars[q];
				    for(int obs=0 ; obs<n_obs ; ++obs){
				        mu_current[obs] += my_slope_var[obs] * my_cluster_coef[my_dum[obs]];
				    }
				} else {
				    for(int obs=0 ; obs<n_obs ; ++obs){
				        mu_current[obs] += my_cluster_coef[my_dum[obs]];
				    }
				}

			}

			// init mu_last if iter == 50 / otherwise, comparison
			if(iter == 50){
				for(int i=0 ; i<n_obs ; ++i){
					mu_last[i] = mu_current[i];
				}

				for(int i=0 ; i<n_obs ; ++i){
					input_mean += input[i];
				}
				input_mean /= n_obs;
				input_mean = fabs(input_mean);

			} else {
				double mu_diff = 0;
				for(int i=0 ; i<n_obs ; ++i){
					mu_diff += fabs(mu_last[i] - mu_current[i]);
				}
				mu_diff = mu_diff / n_obs; // we scale it to have a meaningful value
				// Rprintf("iter %i -- mu_diff = %f (mean = %f)\n", iter, mu_diff, input_mean);

				// if(mu_diff / (0.1 + input_mean) < diffMax){
				if(stopping_crit(input_mean, mu_diff + input_mean, diffMax)){
					// Rprintf("mu_diff BREAK\n");
					break;
				}

				// we update mu_last
				for(int i=0 ; i<n_obs ; ++i){
					mu_last[i] = mu_current[i];
				}

			}
		}


	}

	//
	// Updating the output
	//

	for(int q=0 ; q<Q ; ++q){
		int *my_dum = pdum[q];
		double *my_cluster_coef = pGX[q];

		if(slope_flag[q]){
		    double *my_slope_var = all_slope_vars[q];
		    for(int obs=0 ; obs<n_obs ; ++obs){
		        output[obs] += my_slope_var[obs] * my_cluster_coef[my_dum[obs]];
		    }
		} else {
		    for(int obs=0 ; obs<n_obs ; ++obs){
		        output[obs] += my_cluster_coef[my_dum[obs]];
		    }
		}
	}

	// keeping track of iterations
	int *iterations_all = args->piterations_all;
	iterations_all[v] += iter;

	// saving the fixef coefs
	double *fixef_values = args->fixef_values;
	if(args->save_fixef){
	    for(int m=0 ; m<nb_coef ; ++m){
	        fixef_values[m] += GX[m];
	    }
	}

	bool conv = iter == iterMax ? false : true;

	return(conv);
}

void demean_single_gnl(int v, PARAM_DEMEAN* args){
	// v: variable identifier to demean

	// Algorithm to quickly get the means
	// Q >= 3 => acceleration for 15 iter
	// if no convergence: conv 2 clusters
	// then acceleration again

	// data
	int iterMax = args->iterMax;
	int Q = args->Q;

	if(Q == 2){
		demean_acc_2(v, iterMax, args);
		// demean_acc_gnl(v, iterMax, args);
	} else {
		// 15 iterations
		bool conv = demean_acc_gnl(v, 15, args);

		if(conv == false){
			// 2 convergence
			demean_acc_2(v, iterMax / 2 - 15, args);

			if(Q > 2){
				// re-acceleration
				demean_acc_gnl(v, iterMax / 2, args);
			}
		}
	}

	int *jobdone = args->jobdone;
	jobdone[v] = 1;

}

void stayIdleCheckingInterrupt(bool *stopnow, vector<int> &jobdone, int n_vars, int *counterInside){
	// function that keeps the master thread busy until everything is done

	int nbDone = 0, iter = 1;
	bool isMaster = omp_get_thread_num() == 0;


	while(isMaster && nbDone < n_vars && !(*stopnow)){
		++iter;

		if(iter % 500000000 == 0){
			if(pending_interrupt()){
				(*counterInside)++;
				*stopnow = true;
				break;
			} else {
				// to avoid int overflow:
				iter = 0;
			}
		}

		if(iter % 1000000 == 0){
			nbDone = 0;
			for(int v=0 ; v<n_vars ; v++){
				nbDone += jobdone[v];
			}
		}

	}

}

// Loop over demean_single
// [[Rcpp::export]]
List cpp_demean(SEXP y, SEXP X_raw, SEXP r_weights, int iterMax, double diffMax, SEXP nb_cluster_all,
                SEXP dum_vector, SEXP tableCluster_vector, SEXP slope_flag, SEXP slope_vars,
                SEXP r_init, int checkWeight, int nthreads, bool save_fixef = false){
	// main fun that calls demean_single
	// preformat all the information needed on the clusters
	// y: the dependent variable
	// X_raw: the matrix of the explanatory variables -- can be "empty"

	// when including weights: recreate table values
	// export weights and isWeight bool

	// slope_flag: whether a FE is a varying slope
	// slope_var: the associated variables with varying slopes

	//initial variables
	int Q = Rf_length(nb_cluster_all);
	int *pcluster = INTEGER(nb_cluster_all);
	int n_obs = Rf_length(y);
	bool isWeight = Rf_length(r_weights) != 1;

	// whether we use X_raw
	int n_X = Rf_length(X_raw);
	int n_vars;
	bool useX;
	if(n_X == 1){
		// means X_raw not needed
		n_vars = 1; // only y
		useX = false;
	} else {
		n_vars = n_X / n_obs + 1;
		useX = true;
	}

	// initialisation if needed
	bool isInit = Rf_length(r_init) != 1;
	double *init = REAL(r_init);
	bool saveInit = isInit || init[0] != 0;

	int nb_coef = 0;
	for(int q=0 ; q<Q ; ++q){
		nb_coef += pcluster[q];
	}

	// cluster id for each observation
	vector<int*> pdum(Q);
	pdum[0] = INTEGER(dum_vector);
	for(int q=1 ; q<Q ; ++q){
		pdum[q] = pdum[q - 1] + n_obs;
	}

	// Handling slopes
	int nb_slopes = 0;
	int *pslope_flag = INTEGER(slope_flag);
	for(int q=0 ; q<Q ; ++q){
        nb_slopes += pslope_flag[q];
	}
	bool isSlope = nb_slopes > 0;
	vector<double*> all_slope_vars(Q);
	vector<double> neutral_var(isSlope ? n_obs : 1, 1);

    // we initialize all_slope_vars to the values of slope_vars
    // to neutral_var if not slope
    int index = 0;
    for(int q=0 ; q<Q ; ++q){
        if(pslope_flag[q]){
            all_slope_vars[q] = REAL(slope_vars) + index;
            index += n_obs;
        } else {
            all_slope_vars[q] = neutral_var.data();
        }
    }

	// Handling weights
	// if there are weights: sum_weights
	vector<double> sum_weights(nb_coef, 0);
	vector<double*> psum_weights(Q);
	psum_weights[0] = sum_weights.data();
	for(int q=1 ; q<Q ; ++q){
		psum_weights[q] = psum_weights[q - 1] + pcluster[q - 1];
	}

	double *obs_weights = REAL(r_weights);
	// all_obs_weights: weights for each FE (I created it to take care of slopes)
	vector<double*> all_obs_weights(Q);
	int *table_vector = INTEGER(tableCluster_vector);

	// slope_weights_vector will contain obs_weight[obs]*vars[obs] for slopes, and obs_weight[obs]
	//    for non slopes [thus we initialize at 1 -- default if no weights no slope]
	vector<double> slope_weights_vector(isSlope ? Q * n_obs : 1, 1);

	if(isSlope){

	    // all_obs_weights refer to the values in slope_weights_vector
	    all_obs_weights[0] = slope_weights_vector.data();
	    for(int q=1 ; q<Q ; ++q){
	        all_obs_weights[q] = all_obs_weights[q - 1] + n_obs;
	    }

	    // we compute the values of slope_weights_vector and sum of weights

	    for(int q=0 ; q<Q ; ++q){
	        int *my_dum = pdum[q];
	        double *my_SW = psum_weights[q];
	        double *my_slope_weights = all_obs_weights[q];

	        if(pslope_flag[q]){
	            double *my_slope_var = all_slope_vars[q];
	            if(isWeight){
	                for(int obs=0 ; obs<n_obs ; ++obs){
	                    double var = my_slope_var[obs];
	                    double weight_var = obs_weights[obs] * var;
	                    my_slope_weights[obs] = weight_var;
	                    my_SW[my_dum[obs]] += weight_var * var;
	                }
	            } else {
	                for(int obs=0 ; obs<n_obs ; ++obs){
	                    double var = my_slope_var[obs];
	                    my_slope_weights[obs] = var;
	                    my_SW[my_dum[obs]] += var * var;
	                }
	            }
	        } else {
	            // if NOT a slope var: my_slope_weights eq to weights or 1
	            if(isWeight){
	                for(int obs=0 ; obs<n_obs ; ++obs){
	                    double weight = obs_weights[obs];
	                    my_slope_weights[obs] = weight;
	                    my_SW[my_dum[obs]] += weight;
	                }
	            } else {
	                for(int obs=0 ; obs<n_obs ; ++obs){
	                    my_SW[my_dum[obs]]++;
	                }
	            }
	        }
	    }

	} else {
	    // NO SLOPE

	    // this is always the same weights
	    for(int q=0 ; q<Q ; ++q){
	        all_obs_weights[q] = REAL(r_weights);
	    }

	    if(isWeight){
	        for(int q=0 ; q<Q ; ++q){
	            int *my_dum = pdum[q];
	            double *my_SW = psum_weights[q];
	            for(int obs=0 ; obs<n_obs ; ++obs){
	                my_SW[my_dum[obs]] += obs_weights[obs];
	            }
	        }
	    } else {
	        // we pass the value of table_vector to sum_weights
	        for(int i=0 ; i<nb_coef ; ++i){
	            sum_weights[i] = table_vector[i];
	        }
	    }
	}

	// We update the weight information => slope is (almost) like using weights
	isWeight = isWeight || isSlope;

	// We avoid 0 weight clusters => (otherwise division by 0 leads to NA)
	if(checkWeight || isSlope){
		for(int coef=0 ; coef<nb_coef ; ++coef){
			if(sum_weights[coef] == 0){
				sum_weights[coef] = 1;
			}
		}
	}

	// we put all input variables into a single vector
	// the dep var is the last one
	vector<double> input_values(n_obs*n_vars);

	if(useX){
		double* pX_raw = REAL(X_raw);
		for(int i = 0 ; i < (n_obs*(n_vars - 1)) ; ++i){
			input_values[i] = pX_raw[i];
		}
	}

	double* py = REAL(y);
	int y_start = n_obs*(n_vars - 1);
	for(int i = 0 ; i < n_obs ; ++i){
		input_values[y_start + i] = py[i];
	}

	// output vector:
	vector<double> output_values(n_obs*n_vars, 0);

	if(isInit){
		for(int i=0 ; i<(n_obs*n_vars) ; ++i){
			output_values[i] = init[i];
		}
	}

	// vector of pointers: input/output
	vector<double*> pinput(n_vars);
	vector<double*> poutput(n_vars);
	pinput[0] = input_values.data();
	poutput[0] = output_values.data();
	for(int v=1 ; v<n_vars ; v++){
		pinput[v] = pinput[v - 1] + n_obs;
		poutput[v] = poutput[v - 1] + n_obs;
	}

	// keeping track of iterations
	vector<int> iterations_all(n_vars, 0);
	int *piterations_all = iterations_all.data();

	// save fixef option
	if(useX && save_fixef){
	    stop("save_fixef can be used only when there is no Xs.");
	}

	vector<double> fixef_values(save_fixef ? nb_coef : 1, 0);
	double *pfixef_values = fixef_values.data();

	//
	// Sending variables to envir
	//

	PARAM_DEMEAN args;

	args.n_obs = n_obs;
	args.iterMax = iterMax;
	args.diffMax = diffMax;
	args.Q = Q;
	args.nb_coef = nb_coef;
	args.pdum = pdum;
	args.pcluster = pcluster;
	args.pinput = pinput;
	args.poutput = poutput;
	args.piterations_all = piterations_all;

	// weights + slope:
	args.isWeight = isWeight;
	args.psum_weights = psum_weights;
	args.all_obs_weights = all_obs_weights;
	args.all_slope_vars = all_slope_vars;
	args.slope_flag = pslope_flag;
	args.isSlope = isSlope;

	// save fixef:
	args.save_fixef = save_fixef;
	args.fixef_values = pfixef_values;

	// stopping flag + indicator that job is finished
	bool stopnow = false;
	args.stopnow = &stopnow;
	vector<int> jobdone(n_vars, 0);
	int *pjobdone = jobdone.data();
	args.jobdone = pjobdone;

	int counter = 0;
	int *pcounter = &counter;


	//
	// the main loop
	//

	// enlever les rprintf dans les nthreads jobs
#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
	for(int v = 0 ; v<(n_vars+nthreads) ; ++v){
		// demean_single is the workhorse
		// you get the "mean"

		if(!*(args.stopnow)){
			if(v < n_vars){
				if(Q == 1){
					demean_single_1(v, &args);
				} else {
					demean_single_gnl(v, &args);
				}
			} else if(true && Q != 1){
				stayIdleCheckingInterrupt(&stopnow, jobdone, n_vars, pcounter);
			}
		}

	}


	if(*(args.stopnow)){
		stop("cpp_demean: User interrupt.");
	}

	// Rprintf("Master checking: %i\n", *pcounter);

	//
	// save
	//

	int nrow = useX ? n_obs : 1;
	int ncol = useX ? n_vars - 1 : 1;
	NumericMatrix X_demean(nrow, ncol);
	for(int k=0 ; k<(n_vars - 1) ; ++k){
		int start = k*n_obs;
		for(int i=0 ; i < n_obs ; ++i){
			X_demean(i, k) = input_values[start + i] - output_values[start + i];
		}
	}

	NumericVector y_demean(n_obs);
	for(int i=0 ; i < n_obs ; ++i){
		y_demean[i] = input_values[y_start + i] - output_values[y_start + i];
	}

	// iterations
	IntegerVector iter_final(n_vars);
	for(int v=0 ; v<n_vars ; ++v){
		iter_final[v] = piterations_all[v];
	}

	// if save is requested
	int n = saveInit ? n_obs*n_vars : 1;
	NumericVector saved_output(n);
	for(int i=0 ; i < n ; ++i){
		saved_output[i] = output_values[i];
	}

	// save fixef coef
	n = save_fixef ? nb_coef : 1;
	NumericVector saved_fixef_coef(n);
	for(int i=0 ; i < n ; ++i){
	    saved_fixef_coef[i] = fixef_values[i];
	}


	List res; // a vector and a matrix
	res["X_demean"] = X_demean;
	res["y_demean"] = y_demean;
	res["iterations"] = iter_final;
	res["means"] = saved_output;
	res["fixef_coef"] = saved_fixef_coef;

	return(res);
}




