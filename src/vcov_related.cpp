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


// [[Rcpp::export]]
NumericMatrix cpp_newey_west(NumericMatrix S, NumericVector w, int nthreads){
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
NumericMatrix cpp_newey_west_panel(NumericMatrix S, NumericVector w, IntegerVector unit,
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
NumericMatrix cpp_driscoll_kraay(NumericMatrix S, NumericVector w,
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


/******************************************
 *                                        *
 *                                        *
 *                 Conley                 *
 *                                        *
 *                                        *
 *****************************************/

// we define a class to access the data since access will be pretty intensive
class mat_row_scheme{

    int64_t K = 0;
    int64_t N = 0;
    int64_t n_total = 0;

public:
    std::vector<double> mat;

    mat_row_scheme() = delete;
    mat_row_scheme(mat_row_scheme&);
    mat_row_scheme(Rcpp::NumericMatrix&);
    mat_row_scheme(int, int);

    double& operator()(int64_t i, int64_t k){
        return mat[i * K + k];
    }

    double *data() {return mat.data();}

    void scale(double s){
        for(int64_t i=0 ; i<n_total ; ++i){
            mat[i] *= s;
        }
    }

    void check(){
        Rcout << "CHECK!\n";

        int cpt = 0;
        for(int64_t i=0 ; i<N ; ++i){
            for(int64_t k=0 ; k<K ; ++k){
                if(cpt++ <= 10) Rcout << "i = " << i << ", k = " << k << ", x = " << mat[i * K + k] << "\n";
            }
        }
        Rcout << "...end\n";
    }

    int nrow() {return N;}
    int ncol() {return K;}

};

mat_row_scheme::mat_row_scheme(Rcpp::NumericMatrix &x){

    this->N = x.nrow();
    this->K = x.ncol();
    n_total = N * K;

    // N and K are int64, so product is OK
    mat.resize(n_total);

    // filling scheme is slow, memcpy cannot be leveraged
    for(int64_t i=0 ; i<N ; ++i){
        for(int64_t k=0 ; k<K ; ++k){
            mat[i * K + k] = x(i, k);
        }
    }
}

mat_row_scheme::mat_row_scheme(mat_row_scheme &x){
    // copy

    this->N = x.nrow();
    this->K = x.ncol();

    // N and K are int64, so product is OK
    n_total = N * K;
    mat.resize(n_total);

    // std::copy can be used => fast
    std::copy(x.mat.begin(), x.mat.end(), mat.begin());
}

mat_row_scheme::mat_row_scheme(int N_in, int K_in){

    this->N = N_in;
    this->K = K_in;
    n_total = N * K;

    // N and K are int64, so product is OK
    mat.resize(n_total);

    std::fill(mat.begin(), mat.end(), 0);
}


inline double to_sq(double x){
    return x * x;
}

inline double dist_km(double lon_1, double lat_1, double cos_lat_1,
                      double lon_2, double lat_2, double cos_lat_2){

    double delta_lon = (lon_2 - lon_1) / 2;
    double delta_lat = (lat_2 - lat_1) / 2;

    double a = to_sq(sin(delta_lat)) + cos_lat_1 * cos_lat_2 * to_sq(sin(delta_lon));
    double res = 12752 * asin(fmin(1, sqrt(a)));

    return res;
}

inline double degree_to_radian(double x){
    return x * 3.14159 / 180;
}

inline double fabs_lon(double x, double y){
    // there is a border problem that we take care of

    // this is in radians
    double diff = fabs(x - y);
    // in degrees it would be: diff < 180 ? diff : 360 - diff;
    return diff < 3.14159 ? diff : 6.28318 - diff;
}

inline double fabs_lat(double x, double y){
    // There is no border problem wrt latitude
    return fabs(x - y);
}


// [[Rcpp::export]]
NumericMatrix cpp_vcov_conley(NumericMatrix S, NumericVector lon_rad, NumericVector lat_rad,
                                   const int distance, const double cutoff, int nthreads){
    // S: scores
    // lon_rad/lat_rad: longitude/latitude
    //  => the lat_rad and lon_rad **must** be in radians!!!!
    //
    // IMPORTANT: S, lon_rad, lat_rad must be sorted by lat
    //
    // cutoff, in km
    // kernel type: uniform only, it does not matter
    // distance:
    // - 1: spherical
    // - 2: triangular
    //

    // Notes to self:
    // [border problem]
    // - when we use latitude as main dimension, there are no border problems
    // - border problems (that -170 is close to 170) happen only for longitude
    // - using only the latitude is the way to go, even though the data is usually more scattered across longitude
    // and it would save time to use longitude as the main cutoff, that would imply just too many computations
    //
    // [longitude cutoff]
    // - the longitude cutoff depends on the latitude of the points
    // - I initially used only the latitude of i
    // - but there are two points with different latitudes => thus that breaks symmetry!
    // - in the end I took the avg latitudes
    // - I think that in some rare instances there can still be cases where:
    //   + dist_ij < cutoff, but..
    //   + dist_lon_rad_ij > cutoff_lon_mean_ij (and possibly: dist_lon_rad_ij > cutoff_lon_i OR dist_lon_rad_ij > cutoff_lon_j)
    // - but it should be super edge cases, at the very limit... so not really important
    //
    // [using longitude]
    // - given the problems mentioned above, I should NEVER use longitude as the main dimension
    // - Laurent, remember you first started with longitude in main and although you can't mention here
    // all the problems that you encountered, you did encounter them. So don't.
    //

    if(distance <= 0 || distance > 2){
        Rcpp::stop("'distance' is not valid (internal error).");
    }

    int K = S.ncol();
    int N = S.nrow();

    mat_row_scheme scores(S);

    // utilities
    NumericVector cos_lat(N);
    for(int i=0 ; i<N ; ++i){
        cos_lat[i] = cos(lat_rad[i]);
    }

    mat_row_scheme cum_scores(scores);

    cum_scores.scale(0.5);

    // 1 lat degree is approx 111km
    const double lat_cutoff_rad = degree_to_radian(cutoff / 111);
    const double lon_cutoff_rad_factor = degree_to_radian(cutoff / 111);

    // cutoff_rad_sq used when distance == 2
    const double cutoff_rad_sq = to_sq(degree_to_radian(cutoff) / 111);

#pragma omp parallel for num_threads(nthreads)
    for(int i=0 ; i<N ; ++i){

        const double lon_rad_i = lon_rad[i];
        const double lat_rad_i = lat_rad[i];
        const double cos_lat_i = cos_lat[i];

        // 1 lon degree is approx cos(lat)*111km
        double lon_cutoff_rad = 0;

        bool ok = false;
        double dist_lon_rad = 0;
        double dist_lat_rad = 0;
        double dist = 0;
        double cos_lat_mean = 0;

        for(int i2=i + 1; i2<N ; ++i2){
            dist_lat_rad = fabs_lat(lat_rad[i2], lat_rad_i);
            if(dist_lat_rad > lat_cutoff_rad){
                break;
            }

            dist_lon_rad = fabs_lon(lon_rad[i2], lon_rad_i);

            // we take the avg lat for symmetry
            cos_lat_mean = cos((lat_rad_i + lat_rad[i2]) / 2);
            lon_cutoff_rad = lon_cutoff_rad_factor / cos_lat_mean;

            if(dist_lon_rad > lon_cutoff_rad){
                continue;
            }

            if(distance == 1){
                dist = dist_km(lon_rad_i, lat_rad_i, cos_lat_i,
                               lon_rad[i2], lat_rad[i2], cos_lat[i2]);
                ok = dist <= cutoff;

            } else if(distance == 2){
                // dist: in radian, and put to square
                dist = to_sq(dist_lat_rad) + to_sq(cos_lat_mean * dist_lon_rad);
                ok = dist <= cutoff_rad_sq;
            }

            // if(cpt++ < 10) Rcout << "i = " << i << ", j = " << i2 << ", dij_lon = " << dist_lon_rad << ", dij_lat = " << dist_lat_rad << ", d_ij = " << dist << ", w = " << weight << "\n";

            if(ok){
                // ++n_done;
                for(int k=0 ; k<K ; ++k){
                    cum_scores(i, k) += scores(i2, k);
                }
            }
        }

    }

    // Rcout << "scores:\n";
    // scores.check();
    // Rcout << "cum scores:\n";
    // cum_scores.check();

    // Now the matrix multiplications

    // Rcout << "total done: " << ++n_done << "\n";

    NumericMatrix res(K, K);

    for(int i=0 ; i<N ; ++i){
        for(int k1=0 ; k1<K ; ++k1){
            for(int k2=0 ; k2<K ; ++k2){
                res(k1, k2) += scores(i, k1) * cum_scores(i, k2);
            }
        }
    }

    // we add the transpose
    double tmp = 0;
    for(int k1=0 ; k1<K ; ++k1){
        for(int k2=k1 ; k2<K ; ++k2){
            if(k1 == k2){
                res(k1, k2) *= 2;
            } else {
                tmp = res(k1, k2);
                res(k1, k2) += res(k2, k1);
                res(k2, k1) += tmp;
            }
        }
    }


    return res;
}

































































































