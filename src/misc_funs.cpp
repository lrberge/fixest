#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <Rcpp/Benchmark/Timer.h>
#include <stdio.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]


// Small utility to time -- only for development
double get_time(Timer &timer, bool us = false){
    timer.step("stop");

    int div = us ? 1000 : 1000000;

    NumericVector res_v(timer);
    int n = res_v.length();
    if(n < 2) stop("Problem timer length");
    double res = (res_v[n - 1] - res_v[n - 2]) / div;
    return(res);
}

// [[Rcpp::export]]
bool is_little_endian(){
    int num = 1;
    bool res = *(char *)&num == 1;
    return res;
}

// [[Rcpp::export]]
NumericVector cpp_lgamma(NumericVector x){
	// simple function to compute lgamma of a vector

	int n = x.length();
	NumericVector res(n);

	for(int i=0 ; i<n ; i++){
		res[i] = lgamma(x[i]);
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpp_log_a_exp(double a, NumericVector mu, NumericVector exp_mu){
	// faster this way

	int n = mu.length();
	NumericVector res(n);

	for(int i=0 ; i<n ; i++){
		if(mu[i] < 200){
			res[i] = log(a + exp_mu[i]);
		} else {
			res[i] = mu[i];
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpp_partialDerivative_other(int iterMax, int Q, int N, double epsDeriv, NumericVector ll_d2,	NumericVector dx_dother, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster // must be in cpp index!!!
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// dx_dother: the vector of dx_dother

	int iter;

	int i, q, c;
	int index;
	int sum_cases=0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		// the total number of clusters (eg if man/woman and 10 countries: total of 12 cases)
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericVector clusterDeriv(sum_cases); // the derivatives per cluster, stacked in one vector
	NumericVector sum_lld2(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			index = start[q] + dumMat(i, q);
			sum_lld2[index] += ll_d2(i);
		}
	}

	// the result to return
	NumericVector S(N);
	for(i=0 ; i<N ; i++){
		S[i] = init(i);
	}

	ok = true;
	iter = 0;
	while( ok & (iter<iterMax) ){
		iter++;
		ok = false;

		for(q=0 ; q<Q ; q++){
			R_CheckUserInterrupt();

			// init of the vector
			for(c=start[q] ; c<end[q] ; c++){
				clusterDeriv(c) = 0;
			}

			for(i=0 ; i<N ; i++){
				index = start[q] + dumMat(i, q);
				clusterDeriv(index) += dx_dother(i) + S(i)*ll_d2(i);
			}

			// on finit de calculer clusterDeriv + controle
			for(c=start[q] ; c<end[q] ; c++){
				new_value = -clusterDeriv(c)/sum_lld2[c];
				clusterDeriv(c) = new_value;
				if(fabs(new_value) > epsDeriv){
					ok = true;
				}
			}

			// on ajoute la derivee a S:
			for(i=0 ; i<N ; i++){
				index = start[q] + dumMat(i, q);
				S[i] += clusterDeriv(index);
			}

		}
	}

	// Rprintf("other, nb iter=%i\n", iter);
	if(iter == iterMax){
		Rprintf("[Getting cluster deriv. other] Max iterations reached (%i)\n", iterMax);
	}

	return(S);
}

// Function to get the conditional sum of a matrix
// [[Rcpp::export]]
NumericMatrix cpp_tapply_sum(int Q, NumericMatrix x, IntegerVector dum){
	// Q: nber of classes
	// N: nber of observations
	// x: a matrix
	// dum: the N vector of clusters

	int N = x.nrow();
	int K = x.ncol();

	NumericMatrix res(Q, K);
	int i, q, k;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C

		for(k=0 ; k<K ; k++){
			res(q, k) += x(i, k);
		}
	}

	return(res);
}

// Function to get the conditional sum of a vector
// [[Rcpp::export]]
NumericVector cpp_tapply_vsum(int Q, NumericVector x, IntegerVector dum){
	// Q: nber of classes
	// x: a matrix
	// dum: the N vector of clusters

	int N = x.size();

	NumericVector res(Q);
	int i, q;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C
		res(q) += x(i);
	}

	return(res);
}

// similar a table but faster
// [[Rcpp::export]]
NumericVector cpp_table(int Q, IntegerVector dum){
	// Q: nber of classes
	// dum: the N vector of clusters

	int N = dum.size();

	NumericVector res(Q);
	int i, q;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C
		res(q)++;
	}

	return(res);
}


// [[Rcpp::export]]
IntegerVector cpp_unclassFactor(NumericVector x){
	// x: a sorted integer vector

	int N = x.size();

	IntegerVector res(N);
	int k=1;
	res[0] = 1;

	for(int i=1 ; i<N ; i++){
		if(x(i-1)!=x(i)) k++;
		res[i] = k;
	}

	return(res);
}

// [[Rcpp::export]]
IntegerVector cpp_unik(NumericVector x_sorted, int k_max){
	// x_sorted: a sorted numeric vector

	int n = x_sorted.size();

	IntegerVector res(k_max);
	int k = 1;
	res[0] = x_sorted[0];

	for(int i=1 ; i<n ; i++){
		if(x_sorted(i - 1) != x_sorted(i)){
			res[k] = x_sorted(i);
			k++;
			if(k == k_max) break;
		}
	}

	return(res);
}



//
// Getting the cluster coefficients
//

// [[Rcpp::export]]
NumericMatrix cpp_get_fe_2(SEXP clusterSize,
                           SEXP i_sorted_index_j, SEXP i_sorted_sumFE,
                           SEXP j_sorted_index_i, SEXP j_sorted_sumFE,
                           SEXP r_cumtable_i, SEXP r_cumtable_j){

	int *psize = INTEGER(clusterSize);
	int n_i = psize[0], n_j = psize[1];
	int nb_coef = n_i + n_j;


	int *is_ind_j = INTEGER(i_sorted_index_j);
	double *is_sumFE = REAL(i_sorted_sumFE);
	int *cumtable_i = INTEGER(r_cumtable_i);

	int *js_ind_i = INTEGER(j_sorted_index_i);
	double *js_sumFE = REAL(j_sorted_sumFE);
	int *cumtable_j = INTEGER(r_cumtable_j);

	// vectors
	vector<bool> isRef_i(n_i, false);
	vector<bool> isRef_j(n_j, false);
	vector<double> cluster_coef_i(n_i);
	vector<double> cluster_coef_j(n_j);
	vector<bool> to_visit_i(n_i, true);
	vector<bool> to_visit_j(n_j, true);
	vector<int> cluster_pending_i(n_i);
	vector<int> cluster_pending_j(n_j);

	int n_done = 0, n_pending_i = 0, n_pending_j = 0, i, j_start = 0, j, m, u;
	while(n_done < nb_coef){

		if(n_pending_j == 0){
			// we set the first guy encountered to 0
			for(j=j_start ; to_visit_j[j] == false ; j++){}
			j_start = j + 1;
			isRef_j[j] = true;
			to_visit_j[j] = false;
			n_done++;
			cluster_coef_j[j] = 0;
			n_pending_j = 1;
			cluster_pending_j[0] = j;
		}

		// if(n_pending_i == 0){
		// 	// we set the first guy encountered to 0
		// 	for(i=i_start ; to_visit_i[i] == false ; i++){}
		// 	i_start = i + 1;
		// 	isRef_i[i] = true;
		// 	to_visit_i[i] = false;
		// 	n_done++;
		// 	cluster_coef_i[i] = 0;
		// 	n_pending_i = 1;
		// 	cluster_pending_i[0] = i;
		// }

		n_pending_i = 0;
		for(m=0 ; m<n_pending_j ; m++){
			j = cluster_pending_j[m];
			for(u = j == 0 ? 0 : cumtable_j[j-1] ; u<cumtable_j[j] ; u++){
				i = js_ind_i[u];
				if(to_visit_i[i]){
					cluster_coef_i[i] = js_sumFE[u] - cluster_coef_j[j];

					to_visit_i[i] = false;
					n_done++;
					cluster_pending_i[n_pending_i] = i;
					n_pending_i++;
				}
			}
		}

		n_pending_j = 0;
		for(m=0 ; m<n_pending_i ; m++){
			i = cluster_pending_i[m];
			for(u = i == 0 ? 0 : cumtable_i[i-1] ; u<cumtable_i[i] ; u++){
				j = is_ind_j[u];
				if(to_visit_j[j]){
					cluster_coef_j[j] = is_sumFE[u] - cluster_coef_i[i];

					to_visit_j[j] = false;
					n_done++;
					cluster_pending_j[n_pending_j] = j;
					n_pending_j++;
				}
			}
		}

	}

	NumericMatrix res(nb_coef, 4);

	for(i=0 ; i<n_i ; i++){
		res(i, 0) = 1;
		res(i, 1) = i + 1;
		res(i, 2) = cluster_coef_i[i];
		res(i, 3) = isRef_i[i];
	}

	for(j=0 ; j<n_j ; j++){
		res(j+n_i, 0) = 2;
		res(j+n_i, 1) = j + 1;
		res(j+n_i, 2) = cluster_coef_j[j];
		res(j+n_i, 3) = isRef_j[j];
	}

	return(res);
}



// [[Rcpp::export]]
List cpp_get_fe_gnl(int Q, int N, NumericVector sumFE, IntegerMatrix dumMat, IntegerVector cluster_sizes, IntegerVector obsCluster){
	// This function returns a list of the cluster coefficients for each cluster
	// dumMat: the matrix of cluster ID for each observation, with cpp index style
	// Q, N: nber of clusters / obs
	// cluster_sizes: vector of the number of cases per cluster
	// obsCluster: the integer vector that is equal to order(dum[[g]])
	// RETURN:
	// a list for each cluster of the cluster coefficient value
	// + the last element is the number of clusters that have been set as references (nb_ref)
	// => much optimized version May 2019


	int iter = 0, iterMax = 10000;
	int iter_loop = 0, iterMax_loop = 10000;


	// Creation of the indices to put all the cluster values into a single vector
	int nb_coef = 0;
	IntegerVector nb_ref(Q); // nb_ref takes the nb of elements set as ref

	for(int q=0 ; q<Q ; q++){
		// the total number of clusters (eg if c1: man/woman and c2: 10 countries: total of 12 cases)
		nb_coef += cluster_sizes(q);
	}

	NumericVector cluster_values(nb_coef);
	IntegerVector cluster_visited(nb_coef); // whether a value has been already assigned

	// index of the cluster
	std::vector<int*> pindex_cluster(Q);
	std::vector<int> index_cluster(nb_coef);

	for(int i=0 ; i<nb_coef ; i++){
		index_cluster[i] = i;
	}

	pindex_cluster[0] = index_cluster.data();
	for(int q=1 ; q<Q ; q++){
		pindex_cluster[q] = pindex_cluster[q - 1] + cluster_sizes(q - 1);
	}


	// Now we create the vector of observations for each cluster
	// we need a strating and an end vector as well
	IntegerVector start_cluster(nb_coef), end_cluster(nb_coef);

	int index;
	int k;

	for(int q=0 ; q<Q ; q++){

		// table cluster: nber of elements for each cluster class
		IntegerVector tableCluster(cluster_sizes(q));
		for(int i=0 ; i<N ; i++){
			k = dumMat(i, q);
			tableCluster(k) += 1; // the number of obs per case
		}

		// now creation of the start/end vectors
		for(int k=0 ; k<cluster_sizes(q) ; k++){
			int *pindex = pindex_cluster[q];
			index = pindex[k];
			// index = start(q) + k;
			if(k == 0){
				start_cluster[index] = 0;
				end_cluster[index] = tableCluster[k];
			} else {
				start_cluster[index] = end_cluster[index - 1];
				end_cluster[index] = end_cluster[index - 1] + tableCluster[k];
			}
		}
	}

	// matrix of the clusters that have been computed
	IntegerMatrix mat_done(N, Q);
	IntegerVector rowsums(N, 0);

	// vector of elements to loop over
	IntegerVector id2do(N);
	IntegerVector id2do_next(N);
	int nb2do = N, nb2do_next = N;
	for(int i=0 ; i<nb2do ; i++){
		id2do(i) = i;
		id2do_next(i) = i;
	}

	// Other indices and variables
	int qui_max, obs;
	int rs, rs_max; // rs: row sum
	int id_cluster;
	double other_value;
	bool first;

	//
	// THE MAIN LOOP
	//

	while(iter < iterMax){
		iter++;

		//
		// Finding the row where to put the 0s
		//


		if(iter == 1){
			// 1st iter, we select the first element
			qui_max = 0;
		} else {
			// we find the row that has the maximum of items done

			qui_max = 0;
			rs_max = 0;
			for(int i=0 ; i<nb2do ; i++){
				obs = id2do[i];

				rs = rowsums[obs];

				if(rs == Q-2){
					// if rs == Q-2 => its the maximum possible, no need to go further
					qui_max = obs;
					break;
				} else if(rs<Q && rs>rs_max){
					// this is to handle complicated cases with more than 3+ clusters
					qui_max = obs;
					rs_max = rs;
				} else if(qui_max == 0 && rs == 0){
					// used to initialize qui_max
					qui_max = obs;
				}
			}
		}

		// Rprintf("Obs selected: %i\n", qui_max);

		//
		// Putting the 0s, ie setting the references
		//

		// the first element is spared
		first = true;
		for(int q=0 ; q<Q ; q++){
			if(mat_done(qui_max, q) == 0){
				if(first){
					// we spare the first element
					first = false;
				} else {
					// we set the cluster to 0
					// 1) we find the cluster
					id_cluster = dumMat(qui_max, q);
					// Rprintf("Cluster: %i\n", id_cluster + 1);
					// 2) we get the index of the cluster vector
					int *pindex = pindex_cluster[q];
					index = pindex[id_cluster];
					// 3) we set the cluster value to 0
					cluster_values(index) = 0;
					// 4) we update the mat_done matrix for the elements of this cluster
					for(int i=start_cluster[index] ; i<end_cluster[index] ; i++){
						obs = obsCluster(i, q);
						mat_done(obs, q) = 1;
						rowsums[obs]++;
					}
					// 5) => we save the information on which cluster was set as a reference
					nb_ref(q)++;
				}
			}
		}

		//
		// LOOP OF ALL OTHER UPDATES (CRITICAL)
		//

		iter_loop = 0;
		while(iter_loop < iterMax_loop){
			iter_loop++;

			// Rprintf("nb2do_next: %i -- nb2do: %i\n", nb2do_next, nb2do);

			R_CheckUserInterrupt();

			//
			// Selection of indexes (new way) to be updated
			//

			// initialisation of the observations to cover (before first loop the two are identical)
			if(iter_loop != 1){
				nb2do = nb2do_next;
				for(int i=0 ; i<nb2do ; i++){
					id2do[i] = id2do_next[i];
				}
			}

			nb2do_next = 0;

			for(int i=0 ; i<nb2do ; i++){
				// we compute the rowsum of the obs that still needs to be done
				obs = id2do[i];

				rs = rowsums[obs];

				if(rs < Q-1){
					// you need to check it next time!
					id2do_next[nb2do_next] = obs;
					nb2do_next++;
				} else if(rs == Q-1){
					// means: needs to be updated
					int q;
					for(q=0 ; mat_done(obs, q)!=0 ; q++){
						// selection of the q that is equal to 0
					}

					int *pindex = pindex_cluster[q];
					int index_select = pindex[dumMat(obs, q)];

					// Finding the value of the cluster coefficient
					other_value = 0;
					// Computing the sum of the other cluster values
					// and finding the cluster to be updated (q_select)
					for(int l=0 ; l<Q ; l++){
						// we can loop over all q because cluster_values is initialized to 0
						index = pindex_cluster[l][dumMat(obs, l)];
						other_value += cluster_values(index);
					}

					// the index to update
					cluster_values(index_select) = sumFE(obs) - other_value;

					// Update of the mat_done
					for(int i=start_cluster[index_select] ; i<end_cluster[index_select] ; i++){
						obs = obsCluster(i, q);
						mat_done(obs, q) = 1;
						rowsums[obs]++;
					}
				}
			}

			if(nb2do_next == nb2do) break;

		}

		// out of this loop:
		//	- nb2do == nb2do_next
		// - id2do == id2do_next

		// Check that everything is all right
		if(iter_loop == iterMax_loop){
			Rprintf("Problem getting FE, maximum iterations reached (2nd order loop).");
		}

		// if no more obs to be covered
		if(nb2do_next == 0) break;

	}

	if(iter == iterMax){
		Rprintf("Problem getting FE, maximum iterations reached (1st order loop).");
	}

	// final formatting and save
	List res(Q + 1);

	int *pindex;
	for(int q=0 ; q<Q ; q++){
		NumericVector quoi(cluster_sizes(q));
		pindex = pindex_cluster[q];
		for(k=0 ; k<cluster_sizes(q) ; k++){
			index = pindex[k];
			// index = start(q) + k;
			quoi(k) = cluster_values(index);
		}
		res(q) = quoi;
	}

	res(Q) = nb_ref;

	return(res);
}



// [[Rcpp::export]]
double cpp_ssr_null(NumericVector y){
	// simple fun to compute the ssr of the null ols model
	// 2/3 times faster than pure r

	int n = y.length();

	// the mean
	double y_mean = 0;
	for(int i=0 ; i<n ; i++){
		y_mean += y[i];
	}

	y_mean = y_mean/n;

	double res = 0, value;
	for(int i=0 ; i<n ; i++){
		value = y[i] - y_mean;
		res += value * value;
	}

	return(res);
}

// [[Rcpp::export]]
double cpp_ssq(NumericVector x){
	// simple fun to compute the sum of the square of the elt of a vector
	// 30% fatser than pure r

	int n = x.length();

	// the mean
	double res = 0;
	for(int i=0 ; i<n ; i++){
		res += x[i] * x[i];
	}

	return(res);
}

// [[Rcpp::export]]
List cpp_update_dum(IntegerVector dum, int k_max){
	// much faster than redoing unclas(factor(dum))

	int n = dum.length();
	IntegerVector update_value(k_max);
	IntegerVector dum_new(n);
	int total = 0;
	int k;

	// we find out who's still there
	for(int i=0 ; i<n ; i++){
		k = dum[i] - 1;
		if(update_value[k] == 0){
			update_value[k] = 1;
			total++;
			if(total == k_max){
				break;
			}
		}
	}

	if(total == k_max){
		for(int i=0 ; i<n ; i++){
			dum_new[i] = dum[i];
		}
	} else {
		IntegerVector adjust(k_max, 0);
		int a = 0;

		for(int k=0 ; k<k_max ; k++){
			if(update_value[k] == 0){
				a++;
			}
			adjust[k] = a;
		}

		int k;
		for(int i=0 ; i<n ; i++){
			k = dum[i] - 1;
			dum_new[i] = dum[i] - adjust[k];
		}
	}

	List res;
	res["dum_new"] = dum_new;
	res["keep"] = update_value;

	return(res);
}


// [[Rcpp::export]]
bool cpp_isConstant(NumericVector x){
	// simple fun to see whether a variable is constant
	// it is unexpensive -- not the most useful function however:
	//		for 1e7 obs, you gain 50ms over var(x)==0, but it's still 50ms!

	int n = x.length();
	bool res = true;
	double value = x[0];
	for(int i=1 ; i<n ; i++){
		if(x[i] != value){
			res = false;
			break;
		}
	}

	return(res);
}

// [[Rcpp::export]]
bool cpp_any_na_null(SEXP x){
    // > twice faster than testing the two separately

    int n = Rf_length(x);
    double *px = REAL(x);

    for(int i=0 ; i<n ; ++i){
        double x_tmp = px[i];
        if(std::isnan(x_tmp) || x_tmp == 0){
            return true;
        }
    }

    return false;
}

// [[Rcpp::export]]
int cpp_constant_dum(int k, NumericVector x, IntegerVector dum, bool only_0 = false){
    // number of values of dum for which x is constant

    int n_obs = dum.length();

    double ref = x[0];
    int dum_current = dum[0];

    bool found_different = only_0 ? ref != 0 : false;
    int nb_constant = 0;

    for(int i=1 ; i<n_obs ; ++i){
        if(dum[i] != dum_current){
            // new guy
            dum_current = dum[i];
            if(found_different == false){
                ++nb_constant;
            }

            ref = x[i];
            found_different = only_0 ? ref != 0 : false;

        } else if(!found_different){
            if(x[i] != ref){
                found_different = true;
            }
        }

    }

    if(found_different == false){
        ++nb_constant;
    }

    return nb_constant;
}


//
// Lag related functions //
//

// [[Rcpp::export]]
List cpp_find_duplicates(IntegerVector id, IntegerVector time){
    // we check whether there are duplicated rows
    // if so, we provide information

    int n = id.length();
    int n_dup = 0;
    int obs_dup = 0;
    bool any_dup = false;

    /// finding the first duplicate value
    int i = 0;
    for(i=1 ; i<n ; ++i){
        if(time[i - 1] == time[i]){
            if(id[i - 1] == id[i]){
                any_dup = true;
                break;
            }
        }
    }

    // if dup: we find out the nber
    if(any_dup){
        obs_dup = i; // the 1 is implicitely added to make it r style
        int id_dup = id[i];
        int time_dup = time[i];
        n_dup = 2;
        while(++i<n && id_dup == id[i] && time_dup == time[i]) n_dup++;
    }

    List res;
    res["n_dup"] = n_dup;
    res["obs_dup"] = obs_dup;

    return(res);
}

// [[Rcpp::export]]
int cpp_pgcd(IntegerVector x){
    // quick and dirty, but does not matter

    int n = x.length();

    if(n == 1){
        return(x[0]);
    }

    bool ok = false;
    int pgcd = x[0];

    // the min
    for(int i=1 ; i<n ; ++i){
        if(pgcd > x[i]){
            pgcd = x[i];
        }
    }

    // the denom
    while(!ok && pgcd > 1){
        ok = true;
        for(int i=0 ; i<n ; ++i){
            if(x[i] % pgcd != 0){
                pgcd--;
                ok = false;
                break;
            }
        }
    }


    return pgcd;
}





inline unsigned long long float_to_ull(void *u, int i) {
    unsigned long long *pu_ll = reinterpret_cast<unsigned long long*>(u);
    unsigned long long u_ull = pu_ll[i];
    unsigned long long mask = -(u_ull >> 63) | 0x8000000000000000;
    return (u_ull ^ mask);
}

inline double ull_to_float(unsigned long long *u, int i) {
    unsigned long long u_ull = u[i];
    unsigned long long mask = ((u_ull >> 63) - 1) | 0x8000000000000000;
    unsigned long long res_ull = u_ull ^ mask;
    unsigned long long *pres_ull = &res_ull;
    double *pres = reinterpret_cast<double*>(pres_ull);
    double res = *pres;
    return res;
}

void quf_double(vector<int> &x_uf, double *px, vector<double> &x_unik){

    // x_uf: x unclassed factor
    // px: pointer to x vector (R vector) -- READ ONLY!!!
    // x_unik: empty vector

    int n = x_uf.size();

    if(n == 0) stop("x_unclass: size problem");

    bool debug = false;

    Timer timer;
    if(debug) timer.step("start");

    // variables
    unsigned long long x_ull_current = 0;
    // char: 2 octets: 0-255
    // one double is made of 8 char
    unsigned char x_char = 0;
    vector<unsigned long long> x_ulong(n), x_tmp(n);

    // radix array
    int radix_table[8][256] = { {0} };

    // 1) Counting
    for(int i=0 ; i<n ; ++i){
        x_ull_current = float_to_ull(px, i);
        x_ulong[i] = x_ull_current;

        for(int b=0 ; b<8 ; ++b){
            x_char = ((unsigned char *)&x_ull_current)[b];
            ++radix_table[b][x_char];
        }
    }

    // 1') skipping
    vector<bool> skip_flag(8);
    for(int b = 0 ; b < 8 ; ++b){
        x_char = ((unsigned char *)&x_ull_current)[b];
        skip_flag[b] = radix_table[b][x_char] == n;
    }

    if(debug) Rprintf("Counting and skipping in %.3fms\n", get_time(timer));

    // 2) Cumulating
    for(int b = 0 ; b < 8 ; ++b){
        for(int d = 1 ; d < 256 ; ++d){
            radix_table[b][d] += radix_table[b][d - 1];
        }
    }

    // 3) Sorting
    unsigned long long *x_read = x_ulong.data();
    unsigned long long *x_write = x_tmp.data();
    unsigned long long *p_tmp; // used for swapping

    vector<int> x_order(n), x_unclass(n);
    int *o_read = x_unclass.data();
    int *o_write = x_order.data();
    int *o_tmp;

    int k = 0;
    for(auto &p : x_unclass) p = k++;

    for(int b=0 ; b < 8 ; ++b){
        if(skip_flag[b] == false){
            int index;
            for(int i = n - 1 ; i >= 0 ; --i){
                x_ull_current = x_read[i];
                x_char = ((unsigned char *)&x_ull_current)[b];
                index = --radix_table[b][x_char];
                x_write[index] = x_read[i];
                o_write[index] = o_read[i];
            }

            p_tmp = x_read;
            x_read = x_write;
            x_write = p_tmp;

            o_tmp = o_read;
            o_read = o_write;
            o_write = o_tmp;
        }
    }

    if(debug) Rprintf("Cumulating and sorting in %.3fms\n", get_time(timer));


    if(o_read != x_order.data()){
        // we copy the result
        memcpy(x_order.data(), o_read, sizeof(int) * n);
    }

    // We unclass, starting at 1
    k=1;
    x_unclass[0] = k;

    double xi, xim1 = ull_to_float(x_read, 0);
    x_unik.push_back(xim1);

    for(int i=1 ; i<n ; ++i){
        xi = ull_to_float(x_read, i);
        if(xi != xim1){
            ++k;
            x_unik.push_back(xi);
        }
        x_unclass[i] = k;
        xim1 = xi;
    }

    if(debug) Rprintf("Unclassing in %.3fms\n", get_time(timer));

    // We put into the right order
    for(int i=0 ; i<n ; ++i){
        x_uf[x_order[i]] = x_unclass[i];
    }

    if(debug) Rprintf("Reordering in %.3fms\n", get_time(timer));


}

void quf_int_gnl(vector<int> &x_uf, int *px, vector<double> &x_unik, int x_min){
    // we can sort a range up to 2**31 => ie not the full int range
    // for ranges > 2**31 => as double
    // px: pointer to the values of x

    int n = x_uf.size();

    if(static_cast<int>(x_uf.size()) != n) stop("x_unclass: size problem");

    bool debug = false;

    Timer timer;
    if(debug) timer.step("start");

    // variables
    int x_uint_current = 0;
    vector<int> x_uint(n), x_tmp(n);

    // radix array
    int radix_table[4][256] = { {0} };

    // 1) Counting
    for(int i=0 ; i<n ; ++i){
        x_uint_current = px[i] - x_min;
        x_uint[i] = x_uint_current;

        for(int b=0 ; b<4 ; ++b){
            ++radix_table[b][(x_uint_current >> 8*b) & 0xFF];
        }
    }

    // 1') skipping
    vector<bool> skip_flag(4);
    for(int b = 0 ; b < 4 ; ++b){
        skip_flag[b] = radix_table[b][(x_uint_current >> 8*b) & 0xFF] == n;
    }

    if(debug) Rprintf("Counting and skipping in %.3fms\n", get_time(timer));

    // 2) Cumulating
    for(int b = 0 ; b < 4 ; ++b){
        for(int d = 1 ; d < 256 ; ++d){
            radix_table[b][d] += radix_table[b][d - 1];
        }
    }

    // 3) Sorting
    int *x_read = x_uint.data();
    int *x_write = x_tmp.data();
    int *p_tmp; // used for swapping

    vector<int> x_order(n), x_unclass(n);
    int *o_read = x_unclass.data();
    int *o_write = x_order.data();
    int *o_tmp;

    int k = 0;
    for(auto &p : x_unclass) p = k++;

    for(int b=0 ; b < 4 ; ++b){
        if(skip_flag[b] == false){
            // Rprintf("bin: %i\n", b);
            int index;
            for(int i = n - 1 ; i >= 0 ; --i){
                index = --radix_table[b][(x_read[i] >> 8*b) & 0xFF];
                x_write[index] = x_read[i];
                o_write[index] = o_read[i];
            }

            p_tmp = x_read;
            x_read = x_write;
            x_write = p_tmp;

            o_tmp = o_read;
            o_read = o_write;
            o_write = o_tmp;
        }
    }

    if(debug) Rprintf("Cumulating and sorting in %.3fms\n", get_time(timer));

    if(o_read != x_order.data()){
        // we copy the result, if needed
        memcpy(x_order.data(), o_read, sizeof(int) * n);
    }

    // We unclass
    k=1;
    x_unclass[0] = k;

    double xi, xim1 = x_read[0] + x_min;
    x_unik.push_back(xim1);

    for(int i=1 ; i<n ; ++i){
        xi = x_read[i] + x_min;
        if(xi != xim1){
            ++k;
            x_unik.push_back(static_cast<double>(xi));
        }
        x_unclass[i] = k;
        xim1 = xi;
    }

    if(debug) Rprintf("Unclassing in %.3fms\n", get_time(timer));

    // We put into the right order
    for(int i=0 ; i<n ; ++i){
        x_uf[x_order[i]] = x_unclass[i];
    }

    if(debug) Rprintf("Reordering in %.3fms\n", get_time(timer));


}

void quf_int(vector<int> &x_uf, int *px, vector<double> &x_unik, int x_min, int max_value){
    // Principle:
    // we go through x only once
    // we keep a table of all x values
    // whenever we find a new value of x, we save it

    if(max_value > 1000000) stop("max_value must be lower than 100000.");

    int n = x_uf.size();
    int n_unik = 0; // nber of uniques minus one

    // radix array
    vector<int> x_lookup(max_value + 1, 0);

    int x_tmp, x_pos;
    for(int i=0 ; i<n ; ++i){
        x_tmp = px[i] - x_min;

        if(x_lookup[x_tmp] == 0){
            ++n_unik;
            x_uf[i] = n_unik;
            x_unik.push_back(static_cast<double>(px[i]));
            x_lookup[x_tmp] = n_unik;
        } else {
            x_pos = x_lookup[x_tmp];
            x_uf[i] = x_pos;
        }
    }

}


// [[Rcpp::export]]
List cpp_quf_gnl(SEXP x){

    int n = Rf_length(x);

    vector<int> x_uf(n);
    vector<double> x_unik;

    if(TYPEOF(x) == INTSXP){
        // integer

        int *px = INTEGER(x);
        int x_min = px[0], x_max = px[0], x_tmp;
        for(int i=1 ; i<n ; ++i){
            x_tmp = px[i];
            if(x_tmp > x_max) x_max = x_tmp;
            if(x_tmp < x_min) x_min = x_tmp;
        }

        int max_value = x_max - x_min;

        // Rprintf("max_value = %i, n = %i\n", max_value, n);

        if(max_value < n && max_value < 1000000){
            quf_int(x_uf, px, x_unik, x_min, max_value);
        } else if(max_value < 0x10000000){
            // we don't cover ranges > 2**31
            quf_int_gnl(x_uf, px, x_unik, x_min);
        } else {
            // ranges > 2**31 => as double
            // we need to create a vector of double, otherwise: pointer issue
            vector<double> x_dble(n);
            for(int i=0 ; i<n ; ++i) x_dble[i] = static_cast<double>(px[i]);
            quf_double(x_uf, x_dble.data(), x_unik);
        }

    } else {
        // double
        double *px = REAL(x);
        quf_double(x_uf, px, x_unik);
    }

    List res;
    res["x_uf"] = x_uf;
    res["x_unik"] = x_unik;

    return res;
}







// [[Rcpp::export]]
IntegerVector cpp_lag_obs(IntegerVector id, IntegerVector time, int nlag){
    // in case of ties, we sum
    // must be two consecutive years
    // returns an observation nber of where to find the lagged obs
    // note that we forbid duplicate rows!!!!! => in case of duplicate we use a rule of thumb

    int nobs = id.length();
    int obs;
    int id_current;
    int time_current;
    IntegerVector res(nobs, NA_INTEGER);
    int i, j, diff_time;

    if(nlag > 0){
        i = 0;
        while(i < nobs){
            // R_CheckUserInterrupt(); // this is (too) costly
            id_current = id[i];
            time_current = time[i];
            obs = i + 1; // observation, R style
            j = i + 1;
            while(j < nobs){
                diff_time = time[j] - time_current;
                if(id[j] != id_current){
                    // we start someone else
                    i = j - 1; // minus 1 because after we indent i
                    break;
                } else if(diff_time > nlag){
                    // we are too far => stop
                    break;
                } else if(diff_time == 0){
                    // it's the same time/id => we skip
                    ++i;
                } else if(diff_time < nlag){
                    // not far enough
                } else if(diff_time == nlag){
                    // match!
                    res[j] = obs;
                }

                ++j;
            }
            ++i;
        }
    } else if(nlag < 0){
        /* NOTA:I could have tweaked the previous if() to get rid of the condition
         //      but the code would have lost in clarity.
         // For the lead: opposite to what is done before */
         int nlead = -nlag;
        i = nobs;
        while(i >= 0){
            // R_CheckUserInterrupt(); // this is (too) costly
            id_current = id[i];
            time_current = time[i];
            obs = i + 1; // observation, R style
            j = i - 1;
            while(j >= 0){
                diff_time = time_current - time[j];
                if(id[j] != id_current){
                    // we start someone else
                    i = j + 1; // plus 1 because after we dedent i
                    break;
                } else if(diff_time > nlead){
                    // we are too far => stop
                    break;
                } else if(diff_time == 0){
                    // it's the same time/id => we skip
                    --i;
                } else if(diff_time < nlead){
                    // not far enough
                } else if(diff_time == nlead){
                    // match!
                    res[j] = obs;
                }

                --j;
            }
            --i;
        }
    } else {
        for(int i=0 ; i<nobs ; ++i){
            res[i] = i + 1;
        }
    }

    return(res);
}



