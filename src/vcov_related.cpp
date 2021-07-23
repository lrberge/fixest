#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using std::vector;


// [[Rcpp::export]]
NumericMatrix cpp_newey_west_meat(NumericMatrix S, NumericVector w, int nthreads){
    // Basic NeweyWest for time series
    // note that the data MUST be sorted by period beforehand
    // S: scores
    // w: weights
    // Note that the first weight needs to be halved

    int N = S.nrow();
    int K = S.ncol();

    int L = w.size();
    if(w[L - 1] == 0) L -= 1;
    if(L > N - 1) L = N - 1;

    // We set the parallel scheme depending on the data
    bool par_on_col = K >= L || nthreads == 1;

    int K_sq = K * K;
    std::vector<int> all_k1, all_k2;
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=0 ; k2<K ; ++k2){
            all_k1.push_back(k1);
            all_k2.push_back(k2);
        }
    }

    NumericMatrix meat(K, K);

    if(par_on_col){

        NumericMatrix mat_prod(K, K);

        for(int l=0 ; l<L ; ++l){

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
            for(int index=0 ; index<K_sq ; ++index){
                int k1 = all_k1[index];
                int k2 = all_k2[index];

                double cp_sum = 0;
                for(int i=0 ; i<(N - l) ; ++i){
                    cp_sum += S(i, k1) * S(i + l, k2);
                }

                meat(k1, k2) += w[l] * cp_sum;
            }
        }
    } else {
        // parallelize on each lag
        int n_steps = ceil(1.0 * L / nthreads);
        int step_size = L / n_steps;

        // we avoid race conditions that way

        // I still don't know how to conveniently do the same with Rcpp::NumericMatrix...
        int K_sq = K * K;
        std::vector<double> mat_all_stacked(K_sq * step_size);
        std::vector<double*> p_mat_prods(step_size);
        p_mat_prods[0] = mat_all_stacked.data();
        for(int s=1 ; s<step_size ; ++s){
            p_mat_prods[s] = p_mat_prods[s - 1] + K_sq;
        }

        int L_start = 0;
        int L_end = step_size;
        for(int s=0 ; s<n_steps ; ++s){

#pragma omp parallel for num_threads(nthreads)
            for(int l=L_start ; l<L_end ; ++l){

                // pointer matrix (not a copy)
                double * mat_prod = p_mat_prods[l - L_start];

                for(int k1=0 ; k1<K ; ++k1){
                    for(int k2=0 ; k2<K ; ++k2){
                        double cp_sum = 0;
                        for(int i=0 ; i<(N - l) ; ++i){
                            cp_sum += S(i, k1) * S(i + l, k2);
                        }

                        mat_prod[k1 + k2*K] = cp_sum;
                    }
                }
            }


            // Adding all the matrices
            for(int l=L_start ; l<L_end ; ++l){
                double * mat_prod = p_mat_prods[l - L_start];
#pragma omp parallel for num_threads(nthreads)
                for(int k1=0 ; k1<K ; ++k1){
                    for(int k2=0 ; k2<K ; ++k2){
                        meat(k1, k2) += w[l] * mat_prod[k1 + k2*K];
                    }
                }
            }

            // updating counters
            L_start += step_size;
            L_end += step_size;
            if(L_end > L) L_end = L;

        }


    }

    // Finishing
    // we add the transpose
    NumericMatrix res = clone(meat);
#pragma omp parallel for num_threads(nthreads)
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=0 ; k2<K ; ++k2){
            res(k1, k2) += meat(k2, k1);
        }
    }

    return res;
}

// [[Rcpp::export]]
NumericMatrix cpp_newey_west_panel_meat(NumericMatrix S, NumericVector w, IntegerVector unit,
                                        int G, IntegerVector time, int T, int nthreads){
    // Newey West,  but for panels
    // S: scores
    // w: weights
    // unit: IDs
    // time: must be int
    // the data MUST be sorted by unit and time (in that order)
    // Note that the first weight needs to be halved

    int N = S.nrow();
    int K = S.ncol();

    int L = w.size();
    if(w[L - 1] == 0) L -= 1;
    if(L > T - 1) L = T - 1;

    NumericMatrix meat(K, K);

    // utilities
    NumericVector time_table(T);
    for(int i=0 ; i<N ; ++i){
        ++time_table[time[i] - 1];
    }

    NumericVector time_start(T);
    NumericVector time_end(T);
    time_end[0] = time_table[0];
    for(int t=1 ; t<T ; ++t){
        time_start[t] = time_start[t - 1] + time_table[t - 1];
        time_end[t] = time_end[t - 1] + time_table[t];
    }

    // checking the balance
    bool balanced = true;
    if(unit[0] != 1){
        balanced = false;
    } else {
        int n_unit = 1;
        int t_current = time[0];

        for(int i=1 ; i<N ; ++i){
            if(t_current != time[i]){
                // we change time period
                // we check all is fine
                if(n_unit != G){
                    balanced = false;
                    break;
                }

                // OK
                n_unit = 0;
                t_current = time[i];
            } else if(unit[i] - unit[i - 1] != 1){
                balanced = false;
                break;
            }

            ++n_unit;
        }
    }

    // Rcout << "balanced: " << balanced << "\n";

    int K_sq = K * K;
    std::vector<int> all_k1, all_k2;
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=0 ; k2<K ; ++k2){
            all_k1.push_back(k1);
            all_k2.push_back(k2);
        }
    }

    // l == 0 => easy
#pragma omp parallel for num_threads(nthreads)
    for(int index=0 ; index<K_sq ; ++index){
        int k1 = all_k1[index];
        int k2 = all_k2[index];

        if(k1 > k2) continue;

        double tmp = 0;
        for(int i=0 ; i<N ; ++i){
            tmp += S(i, k1) * S(i, k2);
        }

        meat(k1, k2) = w[0] * tmp;
        if(k1 != k2) meat(k2, k1) = w[0] * tmp;
    }

    if(balanced){
        // computationally easy

        for(int l=1 ; l<L ; ++l){

            int s1 = time_start[l];
            int s2 = time_start[0];

            int nmax = 0;
            for(int t=l ; t<T ; ++t){
                nmax += time_table[t];
            }

#pragma omp parallel for num_threads(nthreads)
            for(int index=0 ; index<K_sq ; ++index){
                int k1 = all_k1[index];
                int k2 = all_k2[index];

                double tmp = 0;
                for(int id=0 ; id<nmax ; ++id){
                    tmp += S(s1 + id, k1) * S(s2 + id, k2);
                }

                meat(k1, k2) += w[l] * tmp;
            }
        }

    } else {
        // we need to do some extra legwork
        // we only make take the product of matching units

        // the rest
        for(int l=1 ; l<L ; ++l){

#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
            for(int index=0 ; index<K_sq ; ++index){
                int k1 = all_k1[index];
                int k2 = all_k2[index];

                double tmp = 0;

                for(int t=l ; t<T ; ++t){

                    // we only make the product of matching units

                    int obs_left = time_start[t];
                    int obs_right = time_start[t - l];

                    int obs_left_max = obs_left + time_table[t];
                    int obs_right_max = obs_right + time_table[t - l];

                    while(obs_left < obs_left_max && obs_right < obs_right_max){
                        if(unit[obs_left] == unit[obs_right]){
                            // Match!
                            tmp += S(obs_left, k1) * S(obs_right, k2);

                            ++obs_left;
                            ++obs_right;
                        } else if(unit[obs_left] < unit[obs_right]){
                            ++obs_left;
                        } else {
                            ++obs_right;
                        }
                    }
                }

                meat(k1, k2) += w[l] * tmp;
            }
        }
    }


    // Finishing
    // we add the transpose
    NumericMatrix res = clone(meat);
#pragma omp parallel for num_threads(nthreads)
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=0 ; k2<K ; ++k2){
            res(k1, k2) += meat(k2, k1);
        }
    }

    return res;

}


// [[Rcpp::export]]
NumericMatrix cpp_driscoll_kraay_meat(NumericMatrix S, NumericVector w,
                                      IntegerVector time, int T, int nthreads){
    // Driscoll and Kraay
    // S: scores
    // w: weights
    // time: must be int
    // the data MUST be sorted by time
    // Note that the first weight needs to be halved

    int N = S.nrow();
    int K = S.ncol();

    int L = w.size();
    if(w[L - 1] == 0) L -= 1;
    if(L > T - 1) L = T - 1;

    NumericMatrix meat(K, K);

    // Scores
    NumericMatrix time_scores(T, K);

    // we sum the scores by period
#pragma omp parallel for num_threads(nthreads)
    for(int k=0 ; k<K ; ++k){
        for(int i=0 ; i<N ; ++i){
            time_scores(time[i] - 1, k) += S(i, k);
        }
    }

    int K_sq = K * K;
    std::vector<int> all_k1, all_k2;
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=0 ; k2<K ; ++k2){
            all_k1.push_back(k1);
            all_k2.push_back(k2);
        }
    }

    for(int l=0 ; l<L ; ++l){
        // X_t' %*% X_t+l
#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for(int index=0 ; index<K_sq ; ++index){
            int k1 = all_k1[index];
            int k2 = all_k2[index];

            if(l == 0 && k1 > k2) continue;

            double tmp = 0;
            for(int t=0 ; t<T-l ; ++t){
                tmp += time_scores(t, k1) * time_scores(t + l, k2);
            }

            meat(k1, k2) += w[l] * tmp;
            if(l == 0 && k1 != k2) meat(k2, k1) += w[l] * tmp;
        }
    }


    // Finishing
    // we add the transpose
    NumericMatrix res = clone(meat);
#pragma omp parallel for num_threads(nthreads)
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=0 ; k2<K ; ++k2){
            res(k1, k2) += meat(k2, k1);
        }
    }

    return res;

}
