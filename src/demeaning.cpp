/*********************************************************************
 * _______________                                                   *
 * || Demeaning ||                                                   *
 * ---------------                                                   *
 *                                                                   *
 * Author: Laurent R. Berge                                          *
 *                                                                   *
 * Workhorse for feols and feglm.                                    *
 *                                                                   *
 * It demeans any variable given in input, the algortihm is not      *
 * a real demeaning algorithm.                                       *
 * It is in fact identical as obtaining the optimal set of           *
 * cluster coefficients (i.e. fixed-effects) in a ML context with    *
 * a Gaussian likelihood (as described in Berge, 2018).              *
 *                                                                   *
 * This way we can leverage the powerful Irons and Tuck acceleration *
 * algorithm. For simple cases, this doesn't matter so much.         *
 * But for messy data (employee-company for instance), it matters    *
 * quite a bit.                                                      *
 *                                                                   *
 * In terms of functionnality it accomodate weights and coefficients *
 * with varying slopes.                                              *
 *                                                                   *
 * Of course any input is **strongly** checked before getting into   *
 * this function.                                                    *
 *                                                                   *
 * I had to apply a trick to accomodate user interrupt in a parallel *
 * setup. It costs a bit, but it's clearly worth it.                 *
 *                                                                   *
 ********************************************************************/

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


/* CHANGELOG:
 *  - 12/09/2020: the vector of fixed-effects id becomes a list, the same that is returned in R (it avoids a deep copy).
 *  So I had to change the R-style indexes into C-style indexes by substratcting 1.
 *  - October 2020:
 *  Complete rewriting of the function. Now there is a class, FEClass, that takes care of computing the FE related stuff.
 *  I included the handling of varying slopes with closed-form.
 *  => this leads to much clearer code.
 */

/* CONVENTION:
 * - suffix _Q means the vector is of length Q
 * - suffix _I (for Identifier) means the vector is of length the total number of FE identifiers
 * - suffix _T means a scalar representing the Total of something, usually the sum of a _Q object
 * - suffix _C means the vector is of length the number of coefficients
 * - suffix VS_C means the vector is of length the number of coefficients for the varying slopes
 * - suffix noVS_C means the vector is of length the number of coefficients for regular fixed-effects (no VS!)
 * - suffix _N means the vector is of length n_obs
 * - prefix p_ means a pointer
 * - nb_id: refers to the fixed-effects IDs
 * - nb_coef: refers to the fixed-effects coefficients. Coef and id are identical in **absence** of varying slopes.
 * - vs: means varying slopes
 */





// Stopping / continuing criteria:
// Functions used inside all loops
inline bool continue_crit(double a, double b, double diffMax){
    // continuing criterion of the algorithm
    double diff = fabs(a - b);
    return ( (diff > diffMax) && (diff/(0.1 + fabs(a)) > diffMax) );
}

inline bool stopping_crit(double a, double b, double diffMax){
    // stopping criterion of the algorithm
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
// What happens if the master thread has finished its job but the lower thread is in an "infinite" loop?
// this is tricky, as such we cannot stop it.
// Solution: create a function keeping the threads idle waiting for the complete job to be done
// BUT I need to add static allocation of threads => performance cost


//
// We introduce a class that handles varying types of SEXP and behaves as a regular matrix
//


// Later I may use the number of observations
// I don't at the moment because everything is strongly checked beforehand
class sVec{
    double *p_dble = nullptr;
    int *p_int = nullptr;

public:
    // several constructors

    // is_int public member
    bool is_int = false;

    sVec(){};
    sVec(SEXP);
    sVec(double *p_x): p_dble(p_x), is_int(false){};
    sVec(int *p_x): p_int(p_x), is_int(true){};
    sVec(std::nullptr_t){};


    double operator[](int i){
        if(is_int) return static_cast<double>(p_int[i]);
        return p_dble[i];
    }

};

sVec::sVec(SEXP x){
    if(TYPEOF(x) == REALSXP){
        is_int = false;
        p_dble = REAL(x);
    } else if(TYPEOF(x) == INTSXP){
        is_int = true;
        p_int = INTEGER(x);
    } else {
        stop("The current SEXP type is not supported by the sVec class.");
    }
}


class sMat{

    std::vector<sVec> p_sVec;
    int n = 0;
    int K = 0;

    sMat() = delete;

public:
    sMat(SEXP);

    int nrow(){return n;};
    int ncol(){return K;};

    sVec operator[](int);
    double operator()(int, int);
};

sMat::sMat(SEXP x){

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
                if(n != n_tmp) stop("When setting up the class sMat: The number of observations in the list is not coherent across columns.");
            }

            K += K_tmp;

            if(TYPEOF(xx) == REALSXP){
                double *p_x = REAL(xx);
                for(int k=0 ; k<K_tmp ; ++k){
                    p_sVec.push_back(sVec(p_x));
                    if(k + 1 < K_tmp) p_x += n;
                }

            } else if(TYPEOF(xx) == INTSXP){
                int *p_x = INTEGER(xx);
                for(int k=0 ; k<K_tmp ; ++k){
                    p_sVec.push_back(sVec(p_x));
                    if(k + 1 < K_tmp) p_x += n;
                }
            } else {
                stop("The current SEXP type is not supported by the sMat class.");
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

        if(n == 1 && K == 1){
            // => absence of data
            n = 0;
            K = 0;

        } else if(TYPEOF(x) == REALSXP){
            double *p_x = REAL(x);
            for(int k=0 ; k<K ; ++k){
                p_sVec.push_back(sVec(p_x));
                if(k + 1 < K) p_x += n;
            }

        } else if(TYPEOF(x) == INTSXP){
            int *p_x = INTEGER(x);
            for(int k=0 ; k<K ; ++k){
                p_sVec.push_back(sVec(p_x));
                if(k + 1 < K) p_x += n;
            }
        } else {
            stop("The current SEXP type is not supported by the sMat class.");
        }
    }
}

sVec sMat::operator[](int k){
    return p_sVec[k];
}

double sMat::operator()(int i, int k){
    return p_sVec[k][i];
}



//
// Now we start a big chunk => computing the varying slopes coefficients
// That's a big job. To simplify it, I created the class FEClass that takes care of it.
//


class simple_mat_with_id{
    // => Access to one of the n_coef matrices of size n_vs x n_vs; all stacked in a single vector
    //
    // px0: origin of the vector (which is of length n_coef * n_vs * n_vs)
    // px_current: working location of the n_vs x n_vs matrix
    // n: n_vs
    // n2: explicit
    // id_current: current identifier. The identifiers range from 0 to (n_coef - 1)

    simple_mat_with_id() = delete;

    double *px0;
    double *px_current;
    int nrow, ncol, n_total, id_current = 0;

public:
    simple_mat_with_id(double* px_in, int nrow_in):
        px0(px_in), px_current(px_in), nrow(nrow_in), ncol(nrow_in), n_total(nrow * ncol) {};
    simple_mat_with_id(double* px_in, int nrow_in, int ncol_in):
        px0(px_in), px_current(px_in), nrow(nrow_in), ncol(ncol_in), n_total(nrow * ncol) {};
    double& operator()(int id, int i, int j);
    double& operator()(int id, int i);
};

inline double& simple_mat_with_id::operator()(int id, int i, int j){
    if(id != id_current){
        id_current = id;
        px_current = px0 + n_total * id;
    }

    return px_current[i + nrow * j];
}

inline double& simple_mat_with_id::operator()(int id, int i){
    if(id != id_current){
        id_current = id;
        px_current = px0 + n_total * id;
    }

    return px_current[i];
}

class FEClass{

    int Q;
    int n_obs;
    int slope_index;
    bool is_weight;
    bool is_slope;

    // Dense vectors that we populate in the class, and their associated pointers
    vector<double> eq_systems_VS_C;
    vector<double*> p_eq_systems_VS_C;

    vector<double> sum_weights_noVS_C;
    vector<double*> p_sum_weights_noVS_C;

    // p_fe_id: pointers to the fe_id vectors
    // p_vs_vars: pointers to the VS variables
    // p_weights: pointer to the weight vector
    // eq_systems_VS_C: vector stacking all the systems of equations (each system is of size n_coef * n_vs * n_vs)
    // p_eq_systems_VS_C: pointer to the right equation system. Of length Q.
    vector<int*> p_fe_id;
    vector<sVec> p_vs_vars;
    double *p_weights = nullptr;

    vector<bool> is_slope_Q;
    vector<bool> is_slope_fe_Q;

    vector<int> nb_vs_Q;
    vector<int> nb_vs_noFE_Q;

    int *nb_id_Q;

    vector<int> coef_start_Q;

    // internal functions
    void compute_fe_coef_internal(int, double *, bool, sVec, double *, double *);
    void compute_fe_coef_2_internal(double *, double *, double *, bool);
    void add_wfe_coef_to_mu_internal(int, double *, double *, bool);

public:

    // Utility class: Facilitates the access to the VS variables
    class simple_mat_of_vs_vars{
        int K_fe;
        vector<sVec> pvars;

    public:
        simple_mat_of_vs_vars(const FEClass*, int);
        double operator()(int, int);
    };

    int nb_coef_T;
    vector<int> nb_coef_Q;

    // constructor:
    FEClass(int n_obs, int Q, SEXP r_weights, SEXP fe_id_list, SEXP r_nb_id_Q, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list);

    // functions
    void compute_fe_coef(double *fe_coef, sVec &mu_in_N);
    void compute_fe_coef(int q, double *fe_coef, double *sum_other_coef_N, double *in_out_C);

    void add_wfe_coef_to_mu(int q, double *fe_coef_C, double *out_N);
    void add_fe_coef_to_mu(int q, double *fe_coef_C, double *out_N);

    void compute_fe_coef_2(double *fe_coef_in_C, double *fe_coef_out_C, double *fe_coef_tmp, double *in_out_C);

    void add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b, double *in_out_C, double *out_N, bool update_beta);

    void compute_in_out(int q, double *in_out_C, sVec &in_N, double *out_N);
};

FEClass::FEClass(int n_obs, int Q, SEXP r_weights, SEXP fe_id_list, SEXP r_nb_id_Q, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list){
    // The constructor does the job of creating the pre-solved system of equations

    // Information on p_slope_flag_Q:
    // - vector of length Q
    // - say sf = p_slope_flag_Q[q]
    //  * abs(sf) number of VARIABLES with varying slope
    //  * sf > 0: we add the fixed-effect
    //  * sf < 0: the FE should NOT be included

    this->n_obs = n_obs;
    this->Q = Q;

    //
    // Step 0: General information
    //

    // fixed-effect id for each observation
    p_fe_id.resize(Q);
    // New version => dum_vector (a vector) is replaced by fe_id_list (a list)
    for(int q=0 ; q<Q ; ++q){
        p_fe_id[q] = INTEGER(VECTOR_ELT(fe_id_list, q));
    }

    nb_id_Q = INTEGER(r_nb_id_Q);


    //
    // Step 1: we check if slopes are needed
    //

    // nb_slopes: number of variables with varying slopes (the FE does not count!)

    int nb_slopes = 0;
    int sf = 0; // slope flag
    int *p_slope_flag_Q = INTEGER(slope_flag_Q);
    vector<bool> is_slope_Q(Q, false);
    vector<bool> is_slope_fe_Q(Q, false);
    vector<int> nb_vs_Q(Q, 0);
    vector<int> nb_vs_noFE_Q(Q, 0);
    vector<int> nb_coef_Q(Q);
    int nb_coef_T = 0;

    for(int q=0 ; q<Q ; ++q){
        //   0: no slope
        // < 0: slope but no fixed-effect
        // > 0: slope WITH fixed-effect
        // here we count the number of slopes only, we exclude the FEs (that's why there's the substraction)

        sf = p_slope_flag_Q[q];
        if(sf != 0){
            nb_slopes += abs(sf);
            is_slope_Q[q] = true;
            nb_vs_Q[q] = abs(sf);
            nb_vs_noFE_Q[q] = abs(sf);

            if(sf > 0){
                ++nb_vs_Q[q];
                is_slope_fe_Q[q] = true;
            }

            // There is n_vs coefficients (and not n_vs squared => this is only in the systems of eq)
            nb_coef_Q[q] = nb_vs_Q[q] * nb_id_Q[q];
        } else {
            nb_coef_Q[q] = nb_id_Q[q];
        }

        nb_coef_T += nb_coef_Q[q];
    }

    // where to start the coefficients
    vector<int> coef_start_Q(Q, 0);
    for(int q=1 ; q<Q ; ++q) coef_start_Q[q] = coef_start_Q[q - 1] + nb_coef_Q[q - 1];

    // Copying (tiny objects)
    this->is_slope_Q = is_slope_Q;
    this->is_slope_fe_Q = is_slope_fe_Q;
    this->nb_vs_Q = nb_vs_Q;
    this->nb_vs_noFE_Q = nb_vs_noFE_Q;
    this->nb_coef_Q = nb_coef_Q;
    this->nb_coef_T = nb_coef_T;
    this->coef_start_Q = coef_start_Q;

    //
    // Step 2: precomputed stuff for non slopes and slopes
    //

    // Weights
    bool is_weight = Rf_length(r_weights) != 1;
    this->is_weight = is_weight;
    p_weights = REAL(r_weights);

    is_slope = nb_slopes > 0;

    // First the non slope coefficients
    // We create the sum of weights

    int nb_coef_noVS_T = 0;
    for(int q=0; q<Q ; ++q){
        if(is_slope_Q[q] == false){
            nb_coef_noVS_T += nb_id_Q[q];
        }
    }

    sum_weights_noVS_C.resize(nb_coef_noVS_T > 0 ? nb_coef_noVS_T : 1);
    std::fill(sum_weights_noVS_C.begin(), sum_weights_noVS_C.end(), 0);

    p_sum_weights_noVS_C.resize(Q);
    p_sum_weights_noVS_C[0] = sum_weights_noVS_C.data();
    for(int q=1 ; q<Q ; ++q){
        p_sum_weights_noVS_C[q] = p_sum_weights_noVS_C[q - 1] + (is_slope_Q[q - 1] == true ? 0 : nb_id_Q[q - 1]);
    }

    // table is already computed
    vector<int*> p_table_id_I(Q);
    p_table_id_I[0] = INTEGER(table_id_I);
    for(int q=1 ; q<Q ; ++q){
        p_table_id_I[q] = p_table_id_I[q - 1] + nb_id_Q[q - 1];
    }


    for(int q=0 ; q<Q ; ++q){
        if(is_slope_Q[q] == true){
            continue;
        }

        double *my_SW = p_sum_weights_noVS_C[q];

        if(is_weight){
            int *my_fe = p_fe_id[q];
            for(int obs=0 ; obs<n_obs ; ++obs){
                my_SW[my_fe[obs] - 1] += p_weights[obs];
            }
        } else {
            int nb_coef = nb_id_Q[q];
            int *my_table = p_table_id_I[q];
            for(int i=0 ; i<nb_coef ; ++i){
                my_SW[i] = my_table[i];
            }
        }

    }


    if(is_weight && nb_coef_noVS_T > 0){
        // Checking the 0-weights => we set them to 1 to wavoid division by 0
        for(int c=0 ; c<nb_coef_noVS_T ; ++c){
            if(sum_weights_noVS_C[c] == 0){
                sum_weights_noVS_C[c] = 1;
            }
        }

    }

    // Then the slopes
    if(nb_slopes > 0){

        // A) Meta variables => the ones containing the main information

        // slope_vars_list: R list
        // p_vs_vars.resize(nb_slopes);
        // for(int v=0 ; v<nb_slopes ; ++v){
        //     p_vs_vars[v] = REAL(VECTOR_ELT(slope_vars_list, v));
        // }
        sMat m_slopes(slope_vars_list);
        p_vs_vars.resize(nb_slopes);
        for(int v=0 ; v<nb_slopes ; ++v){
            p_vs_vars[v] = m_slopes[v];
        }

        // B) Computing the coefficients of the systems of equations

        int nb_vs_coef_T = 0;
        for(int q=0 ; q<Q ; ++q){
            nb_vs_coef_T += nb_vs_Q[q] * nb_vs_Q[q] * nb_id_Q[q];
        }

        // The sys of eqs, all coefs to 0
        eq_systems_VS_C.resize(nb_vs_coef_T);
        std::fill(eq_systems_VS_C.begin(), eq_systems_VS_C.end(), 0);

        p_eq_systems_VS_C.resize(Q);
        p_eq_systems_VS_C[0] = eq_systems_VS_C.data();
        for(int q=1 ; q<Q ; ++q){
            p_eq_systems_VS_C[q] = p_eq_systems_VS_C[q - 1] + nb_vs_Q[q - 1] * nb_vs_Q[q - 1] * nb_id_Q[q - 1];
        }

        for(int q=0 ; q<Q ; ++q){
            if(is_slope_Q[q] == false) continue;

            simple_mat_of_vs_vars VS_mat(this, q);
            simple_mat_with_id my_system(p_eq_systems_VS_C[q], nb_vs_Q[q]);

            int V = nb_vs_Q[q];
            int *my_fe = p_fe_id[q];
            int nb_coef = nb_id_Q[q];

            for(int i=0 ; i<n_obs ; ++i){
                for(int v1=0 ; v1<V ; ++v1){
                    for(int v2=0 ; v2<=v1 ; ++v2){
                        if(is_weight){
                            my_system(my_fe[i] - 1, v1, v2) += VS_mat(i, v1) * VS_mat(i, v2) * p_weights[i];
                        } else {
                            my_system(my_fe[i] - 1, v1, v2) += VS_mat(i, v1) * VS_mat(i, v2);
                        }
                    }
                }
            }

            // Finishing the computation of the system (symmetry)
            for(int c=0 ; c<nb_coef ; ++c){
                for(int v1=0 ; v1<V ; ++v1){
                    for(int v2=0 ; v2<v1 ; ++v2){
                        my_system(c, v2, v1) = my_system(c, v1, v2);
                    }
                }
            }


            // Precomputing the solver coefficients
            double my_row_coef = 0;
            for(int c=0 ; c<nb_coef ; ++c){
                for(int v=0 ; v<V ; ++v){
                    if(my_system(c, v, v) == 0){
                        // The pivot is equal to 0 => overidentified system
                        for(int i=v + 1 ; i<V ; ++i){
                            my_system(c, i, v) = 0;
                        }
                    } else {
                        for(int i=v + 1 ; i<V ; ++i){
                            my_row_coef = my_system(c, i, v) / my_system(c, v, v);
                            my_system(c, i, v) = my_row_coef;
                            for(int j=v + 1 ; j<V ; ++j){
                                my_system(c, i, j) -= my_row_coef * my_system(c, v, j);
                            }
                        }
                    }
                }
            }

            // We end up with all the (pre-solved) systems of equations
        }

    }
}


// Overloaded versions
void FEClass::compute_fe_coef(double *fe_coef_C, sVec &mu_in_N){
    // mu: length n_obs, vector giving sum_in_out
    // fe_coef: vector receiving the cluster coefficients

    // single FE version
    compute_fe_coef_internal(0, fe_coef_C, true, mu_in_N, nullptr, nullptr);
}

void FEClass::compute_fe_coef(int q, double *fe_coef_C, double *sum_other_coef_N, double *in_out_C){
    // mu: length n_obs, vector giving sum_in_out
    // fe_coef: vector receiving the cluster coefficients

    // multiple FE version
    compute_fe_coef_internal(q, fe_coef_C, false, nullptr, sum_other_coef_N, in_out_C);

}

void FEClass::compute_fe_coef_internal(int q, double *fe_coef_C, bool is_single, sVec mu_in_N, double *sum_other_coef_N, double *in_out_C){
    // mu: length n_obs, vector giving sum_in_out
    // fe_coef: vector receiving the cluster coefficients

    int V = nb_vs_Q[q];
    int *my_fe = p_fe_id[q];
    int nb_coef = nb_coef_Q[q];
    int nb_id = nb_id_Q[q];

    double *my_fe_coef = fe_coef_C + coef_start_Q[q];

    if(is_slope_Q[q] == false){

        double *my_SW = p_sum_weights_noVS_C[q];

        if(is_single){

            for(int obs=0 ; obs<n_obs ; ++obs){
                if(is_weight){
                    my_fe_coef[my_fe[obs] - 1] += p_weights[obs] * mu_in_N[obs];
                } else {
                    my_fe_coef[my_fe[obs] - 1] += mu_in_N[obs];
                }
            }

        } else {

            double *sum_in_out = in_out_C + coef_start_Q[q];

            // initialize cluster coef
            for(int m=0 ; m<nb_coef ; ++m){
                my_fe_coef[m] = sum_in_out[m];
            }

            // looping sequentially over the sum of other coefficients
            for(int i=0 ; i<n_obs ; ++i){
                my_fe_coef[my_fe[i] - 1] -= sum_other_coef_N[i];
            }

        }

        // Rcout << "is_weight:" << is_weight << ", is_single: " << is_single << "\n";

        // Rcout << "Coefs:\n ";
        for(int m=0 ; m<nb_coef ; ++m){
            // Rcout << my_fe_coef[m] << " ("  << my_SW[m] << "), ";

            my_fe_coef[m] /= my_SW[m];
        }
        // Rcout << "\n\n";

    } else {
        // Setting up => computing the raw coefficient of the last column
        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_system(p_eq_systems_VS_C[q], nb_vs_Q[q]);

        // We use the vector of FE coefficients
        simple_mat_with_id rhs(my_fe_coef, nb_vs_Q[q], 1);

        if(is_single){

            for(int i=0 ; i<n_obs ; ++i){
                for(int v1=0 ; v1<V ; ++v1){
                    if(is_weight){
                        rhs(my_fe[i] - 1, v1) += VS_mat(i, v1) * mu_in_N[i] * p_weights[i];
                    } else {
                        rhs(my_fe[i] - 1, v1) += VS_mat(i, v1) * mu_in_N[i];
                    }
                }
            }

        } else {

            // initialization
            double *sum_in_out = in_out_C + coef_start_Q[q];
            for(int m=0 ; m<nb_coef ; ++m){
                my_fe_coef[m] = sum_in_out[m];
            }

            for(int i=0 ; i<n_obs ; ++i){
                for(int v=0 ; v<V ; ++v){
                    rhs(my_fe[i] - 1, v) -= VS_mat(i, v) * sum_other_coef_N[i];
                }
            }

        }


        // "Correcting" the last column
        for(int c=0 ; c<nb_id ; ++c){
            for(int v=0 ; v<V ; ++v){
                for(int v1=v + 1 ; v1<V ; ++v1){
                    // -= row_coef * value from the v row
                    rhs(c, v1) -= my_system(c, v1, v) * rhs(c, v);
                }
            }
        }

        // Solving => we store the solution in rhs
        double val = 0;
        for(int c=0 ; c<nb_id ; ++c){
            // We backward solve
            for(int v=V - 1 ; v>=0 ; --v){
                if(my_system(c, v, v) == 0){
                    rhs(c, v) = 0;
                } else {
                    val = rhs(c, v);
                    for(int v_done=v + 1 ; v_done<V ; ++v_done){
                        val -= rhs(c, v_done) * my_system(c, v, v_done);
                    }
                    rhs(c, v) = val / my_system(c, v, v);
                }
            }
        }

    }

}


void FEClass::compute_fe_coef_2_internal(double *fe_coef_in_out_C, double *fe_coef_tmp, double *in_out_C, bool step_2 = false){
    // The aim here is to update the value of b given the value of a
    // fe_coef_tmp is always of size the 2nd FE
    // fe_in/out are always the size the 1st FE
    // a: origin
    // b: destination

    int index_a = 0;
    int index_b = 1;

    double *my_fe_coef_a;
    double *my_fe_coef_b;

    if(step_2){
        index_a = 1;
        index_b = 0;

        my_fe_coef_a = fe_coef_tmp;
        my_fe_coef_b = fe_coef_in_out_C;
    } else {

        my_fe_coef_a = fe_coef_in_out_C;
        my_fe_coef_b = fe_coef_tmp;
    }

    int V_a = nb_vs_Q[index_a];
    int V_b = nb_vs_Q[index_b];

    int *my_fe_a = p_fe_id[index_a];
    int *my_fe_b = p_fe_id[index_b];

    // int nb_coef_a = nb_coef_Q[index_a];
    int nb_coef_b = nb_coef_Q[index_b];

    // int nb_id_a = nb_id_Q[index_a];
    int nb_id_b = nb_id_Q[index_b];

    // double *my_in_out_a = in_out_C + coef_start_Q[index_a];
    double *my_in_out_b = in_out_C + coef_start_Q[index_b];

    bool is_slope_a = is_slope_Q[index_a];
    bool is_slope_b = is_slope_Q[index_b];

    // VSlope utilities
    simple_mat_of_vs_vars VS_mat_a(this, index_a);
    simple_mat_with_id my_vs_coef_a(my_fe_coef_a, V_a, 1);

    simple_mat_of_vs_vars VS_mat_b(this, index_b);
    simple_mat_with_id my_vs_coef_b(my_fe_coef_b, V_b, 1);


    // initialize cluster coef
    for(int m=0 ; m<nb_coef_b ; ++m){
        my_fe_coef_b[m] = my_in_out_b[m];
    }

    if(is_slope_b == false){

        // Adding the values of the first FE
        for(int i=0 ; i<n_obs ; ++i){
            if(is_slope_a){
                for(int v=0 ; v<V_a ; ++v){
                    if(is_weight){
                        my_fe_coef_b[my_fe_b[i] - 1] -= my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v) * p_weights[i];
                    } else {
                        my_fe_coef_b[my_fe_b[i] - 1] -= my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v);
                    }
                }
            } else if(is_weight){
                my_fe_coef_b[my_fe_b[i] - 1] -= my_fe_coef_a[my_fe_a[i] - 1] * p_weights[i];
            } else {
                my_fe_coef_b[my_fe_b[i] - 1] -= my_fe_coef_a[my_fe_a[i] - 1];
            }
        }

        double *my_SW = p_sum_weights_noVS_C[index_b];
        for(int m=0 ; m<nb_coef_b ; ++m){
            my_fe_coef_b[m] /= my_SW[m];
        }

    } else {
        simple_mat_with_id my_system_b(p_eq_systems_VS_C[index_b], V_b);
        simple_mat_with_id rhs_b(my_fe_coef_b, V_b, 1);

        double tmp = 0;

        for(int i=0 ; i<n_obs ; ++i){

            // if(is_slope_a){
            //     tmp = 0;
            //     for(int v=0 ; v<V_a ; ++v){
            //         if(is_weight){
            //             tmp += my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v) * p_weights[i];
            //         } else {
            //             tmp += my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v);
            //         }
            //     }
            //
            // } else if(is_weight){
            //     tmp = my_fe_coef_a[my_fe_a[i] - 1] * p_weights[i];
            // } else {
            //     tmp = my_fe_coef_a[my_fe_a[i] - 1];
            // }
            if(is_slope_a){
                tmp = 0;
                for(int v=0 ; v<V_a ; ++v){
                    tmp += my_vs_coef_a(my_fe_a[i] - 1, v) * VS_mat_a(i, v);
                }

            } else {
                tmp = my_fe_coef_a[my_fe_a[i] - 1];
            }

            for(int v=0 ; v<V_b ; ++v){
                if(is_weight){
                    rhs_b(my_fe_b[i] - 1, v) -= VS_mat_b(i, v) * tmp * p_weights[i];
                } else {
                    rhs_b(my_fe_b[i] - 1, v) -= VS_mat_b(i, v) * tmp;
                }
            }
        }

        // "Correcting" the last column
        for(int c=0 ; c<nb_id_b ; ++c){
            for(int v=0 ; v<V_b ; ++v){
                for(int v1=v + 1 ; v1<V_b ; ++v1){
                    // -= row_coef * value from the v row
                    rhs_b(c, v1) -= my_system_b(c, v1, v) * rhs_b(c, v);
                }
            }
        }

        // Solving => we store the solution in rhs
        double val = 0;
        for(int c=0 ; c<nb_id_b ; ++c){
            // We backward solve
            for(int v=V_b - 1 ; v>=0 ; --v){
                if(my_system_b(c, v, v) == 0){
                    rhs_b(c, v) = 0;
                } else {
                    val = rhs_b(c, v);
                    for(int v_done=v + 1 ; v_done<V_b ; ++v_done){
                        val -= rhs_b(c, v_done) * my_system_b(c, v, v_done);
                    }
                    rhs_b(c, v) = val / my_system_b(c, v, v);
                }
            }
        }
    }

}

void FEClass::compute_fe_coef_2(double *fe_coef_in_C, double *fe_coef_out_C, double *fe_coef_tmp, double *in_out_C){
    // Specific to the 2-FEs case
    // This way we avoid creating and using a temp object of length N



    //
    // Step 1: Updating b
    //

    compute_fe_coef_2_internal(fe_coef_in_C, fe_coef_tmp, in_out_C);

    // Rcout << "Coefs IN:\nFirst dim: ";
    // for(int i=0 ; i<nb_coef_Q[0] ; ++i){
    //     Rcout << fe_coef_in_C[i] << ", ";
    // }
    // Rcout << "\n";

    //
    // Step 2: Updating a
    //

    compute_fe_coef_2_internal(fe_coef_out_C, fe_coef_tmp, in_out_C, true);

    // Rcout << "Coefs OUT:\nFirst dim: ";
    // for(int i=0 ; i<nb_coef_Q[0] ; ++i){
    //     Rcout << fe_coef_out_C[i] << ", ";
    // }
    // Rcout << "\n";

}

void FEClass::add_wfe_coef_to_mu(int q, double *fe_coef_C, double *out_N){
    // We add the weighted FE coefficient to each observation

    add_wfe_coef_to_mu_internal(q, fe_coef_C, out_N, true);

}

void FEClass::add_fe_coef_to_mu(int q, double *fe_coef_C, double *out_N){
    // We add the FE coefficient to each observation -- NO WEIGHTS!!!

    // Single FE version
    add_wfe_coef_to_mu_internal(q, fe_coef_C, out_N, false);

}

void FEClass::add_wfe_coef_to_mu_internal(int q, double *fe_coef_C, double *out_N, bool add_weights){
    // We just add the weighted (or not) FE coefficient to each observation
    // We need to be careful for VS

    int V = nb_vs_Q[q];
    int *my_fe = p_fe_id[q];

    double *my_fe_coef = fe_coef_C + coef_start_Q[q];

    bool use_weights = add_weights && is_weight;

    if(is_slope_Q[q] == false){

        for(int i=0 ; i<n_obs ; ++i){
            if(use_weights){
                out_N[i] += my_fe_coef[my_fe[i] - 1] * p_weights[i];
            } else {
                out_N[i] += my_fe_coef[my_fe[i] - 1];
            }
        }

    } else {
        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_vs_coef(my_fe_coef, nb_vs_Q[q], 1);

        for(int i=0 ; i<n_obs ; ++i){
            for(int v=0 ; v<V ; ++v){
                if(use_weights){
                    out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v) * p_weights[i];
                } else {
                    out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v);
                }
            }
        }

    }

}

void FEClass::add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b, double *in_out_C, double *out_N, bool update_beta = true){
    // We add the value of the FE coefficients to each observation


    //
    // Step 1: we update the coefficients of b
    //

    if(update_beta){
        compute_fe_coef_2_internal(fe_coef_a, fe_coef_b, in_out_C, out_N);
    }

    //
    // Step 2: we add the value of each coef
    //

    for(int q=0 ; q<2 ; ++q){

        double *my_fe_coef = q == 0 ? fe_coef_a : fe_coef_b;
        int *my_fe = p_fe_id[q];
        bool is_slope = is_slope_Q[q];
        int V = nb_vs_Q[q];

        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_vs_coef(my_fe_coef, nb_vs_Q[q], 1);

        for(int i=0 ; i<n_obs ; ++i){
            if(is_slope){
                for(int v=0 ; v<V ; ++v){
                    out_N[i] += my_vs_coef(my_fe[i] - 1, v) * VS_mat(i, v);
                }
            } else {
                out_N[i] += my_fe_coef[my_fe[i] - 1];
            }
        }

    }

}


void FEClass::compute_in_out(int q, double *in_out_C, sVec &in_N, double *out_N){
    // output: vector of length the number of coefficients

    int V = nb_vs_Q[q];
    int *my_fe = p_fe_id[q];

    double *sum_in_out = in_out_C + coef_start_Q[q];

    if(is_slope_Q[q] == false){

        for(int i=0 ; i<n_obs ; ++i){
            if(is_weight){
                sum_in_out[my_fe[i] - 1] += (in_N[i] - out_N[i]) * p_weights[i];
            } else {
                sum_in_out[my_fe[i] - 1] += (in_N[i] - out_N[i]);
            }
        }

    } else {
        simple_mat_of_vs_vars VS_mat(this, q);
        simple_mat_with_id my_vs_sum_in_out(sum_in_out, nb_vs_Q[q], 1);

        for(int i=0 ; i<n_obs ; ++i){
            for(int v=0 ; v<V ; ++v){
                if(is_weight){
                    my_vs_sum_in_out(my_fe[i] - 1, v) += (in_N[i] - out_N[i]) * VS_mat(i, v) * p_weights[i];
                } else {
                    my_vs_sum_in_out(my_fe[i] - 1, v) += (in_N[i] - out_N[i]) * VS_mat(i, v);
                }
            }
        }

    }

}

FEClass::simple_mat_of_vs_vars::simple_mat_of_vs_vars(const FEClass *FE_info, int q){
    // We set up the matrix
    int start = 0;
    for(int l=0 ; l<q ; ++l) start += FE_info->nb_vs_noFE_Q[l];

    int K = FE_info->nb_vs_noFE_Q[q];
    pvars.resize(K);
    for(int k=0 ; k<K ; ++k){
        pvars[k] = FE_info->p_vs_vars[start + k];
    }

    if(FE_info->is_slope_fe_Q[q]){
        K_fe = K;
    } else {
        K_fe = -1;
    }
}

inline double FEClass::simple_mat_of_vs_vars::operator()(int i, int k){
    if(k == K_fe){
        return 1;
    }

    return pvars[k][i];
}



//
// END OF FE Class
//


// List of objects, used to
// lighten the writting of the functions
struct PARAM_DEMEAN{
    int n_obs;
    int Q;
    int nb_coef_T;
    int iterMax;
    double diffMax;

    // iterations
    int *p_iterations_all;

    // vectors of pointers
    vector<sVec> p_input;
    vector<double*> p_output;

    // saving the fixed effects
    bool save_fixef;
    double *fixef_values;

    // FE information
    FEClass *p_FE_info;

    // stopflag
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
    int nb_coef_T = args->nb_coef_T;

    vector<sVec> &p_input = args->p_input;
    vector<double*> &p_output = args->p_output;

    // fe_info
    FEClass &FE_info = *(args->p_FE_info);

    // vector of fixed-effects coefficients initialized at 0
    vector<double> fe_coef(nb_coef_T, 0);
    double *p_fe_coef = fe_coef.data();

    // interruption handling
    bool isMaster = omp_get_thread_num() == 0;
    bool *pStopNow = args->stopnow;
    if(isMaster){
        if(pending_interrupt()){
            *pStopNow = true;
        }
    }

    // the input & output
    sVec &input = p_input[v];
    double *output = p_output[v];

    // We compute the FEs
    FE_info.compute_fe_coef(p_fe_coef, input);

    // Output:
    FE_info.add_fe_coef_to_mu(0, p_fe_coef, output);

    // saving the fixef coefs
    double *fixef_values = args->fixef_values;
    if(args->save_fixef){
        for(int m=0 ; m<nb_coef_T ; ++m){
            fixef_values[m] = fe_coef[m];
        }
    }


}

void demean_acc_2(int v, int iterMax, PARAM_DEMEAN *args){

    //
    // Setting up
    //

    // loading data
    int n_obs = args->n_obs;
    double diffMax = args->diffMax;

    // input output
    vector<sVec> &p_input = args->p_input;
    vector<double*> &p_output = args->p_output;
    sVec &input = p_input[v];
    double *output = p_output[v];

    int *iterations_all = args->p_iterations_all;

    // FE info
    FEClass &FE_info = *(args->p_FE_info);

    // vector of FE coefs => needs not be equal to nb_coef_T (because Q can be higher than 2)
    int n_a = FE_info.nb_coef_Q[0];
    int n_b = FE_info.nb_coef_Q[1];
    int nb_coef_2 = n_a + n_b;

    // conditional sum of input minus output
    vector<double> sum_input_output(nb_coef_2, 0);
    double *p_sum_in_out = sum_input_output.data();

    for(int q=0 ; q<2 ; ++q){
        FE_info.compute_in_out(q, p_sum_in_out, input, output);
    }

    // Rcout << "Value of in_out:\n";
    // for(int c=0 ; c<nb_coef_2 ; ++c){
    //     if(c == n_a) Rcout << "\n";
    //     Rcout << p_sum_in_out[c] << ", ";
    //
    // }
    // Rcout << "\n";

    // Rprintf("a_tilde: %.3f, %.3f, %.3f, %.3f\n", a_tilde[0], a_tilde[1], a_tilde[2], a_tilde[3]);

    // interruption handling
    bool isMaster = omp_get_thread_num() == 0;
    bool *pStopNow = args->stopnow;
    double flop = 20.0 * static_cast<double>(n_obs); // rough estimate nber operation per iter
    int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second

    //
    // IT iteration (preparation)
    //

    vector<double> X(n_a, 0);
    vector<double> GX(n_a);
    vector<double> GGX(n_a);

    double *p_X   = X.data();
    double *p_GX  = GX.data();
    double *p_GGX = GGX.data();

    vector<double> delta_GX(n_a);
    vector<double> delta2_X(n_a);

    vector<double> coef_beta(n_b);
    double *p_coef_beta = coef_beta.data();

    //
    // the main loop
    //

    // first iteration => update GX
    FE_info.compute_fe_coef_2(p_X, p_GX, p_coef_beta, p_sum_in_out);

    // Rcout << "Coefs at first iteration:\na: ";
    // for(int i=0 ; i<n_a ; ++i){
    //     Rcout << p_GX[i] << ", ";
    // }
    // Rcout << "\nb: ";
    // for(int i=0 ; i<n_b ; ++i){
    //     Rcout << p_coef_beta[i] << ", ";
    // }
    // Rcout << "\n";

    // For the stopping criterion on total addition
    // vector<double> mu_last(n_obs, 0);
    // double input_mean = 0;
    double ssr = 0;

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
        FE_info.compute_fe_coef_2(p_GX, p_GGX, p_coef_beta, p_sum_in_out);

        // X ; update of the fixed-effects coefficients
        numconv = dm_update_X_IronsTuck(n_a, X, GX, GGX, delta_GX, delta2_X);
        if(numconv) break;

        // GX -- origin: X, destination: GX
        FE_info.compute_fe_coef_2(p_X, p_GX, p_coef_beta, p_sum_in_out);

        keepGoing = false;
        for(int i=0 ; i<n_a ; ++i){
            if(continue_crit(X[i], GX[i], diffMax)){
                keepGoing = true;
                break;
            }
        }

        // Other stopping criterion: change to SSR very small
        if(iter % 50 == 0){

            vector<double> mu_current(n_obs, 0);
            double *p_mu = mu_current.data();

            FE_info.add_2_fe_coef_to_mu(p_GX, p_coef_beta, p_sum_in_out, p_mu);

            // init ssr if iter == 50 / otherwise, comparison
            if(iter == 50){
                ssr = 0;
                double resid;
                for(int i=0 ; i<n_obs ; ++i){
                    resid = input[i] - mu_current[i];
                    ssr += resid*resid;
                }

            } else {
                double ssr_old = ssr;

                // we compute the new SSR
                ssr = 0;
                double resid;
                for(int i=0 ; i<n_obs ; ++i){
                    resid = input[i] - mu_current[i];
                    ssr += resid*resid;
                }

                // if(isMaster) Rprintf("iter %i -- SSR = %.0f (diff = %.0f)\n", iter, ssr, ssr_old - ssr);

                if(stopping_crit(ssr_old, ssr, diffMax)){
                    break;
                }

            }
        }

    }

    //
    // we update the result (output)
    //

    // we end with a last iteration
    double *p_alpha_final = p_GX;
    double *p_beta_final = p_coef_beta;

    FE_info.compute_fe_coef_2(p_alpha_final, p_alpha_final, p_beta_final, p_sum_in_out);

    // Rcout << "Output before:\n";
    // for(int i=0 ; i<6 ; ++i){
    //     Rcout << output[i] << ",";
    // }

    // we update the output
    FE_info.add_2_fe_coef_to_mu(p_alpha_final, p_beta_final, p_sum_in_out, output, false);

    // Rcout << "\nOutput after:\n";
    // for(int i=0 ; i<6 ; ++i){
    //     Rcout << output[i] << ",";
    // }
    // Rcout << "\n";
    //
    // Rcout << "Final coefs:\na: ";
    // for(int i=0 ; i<n_a ; ++i){
    //     Rcout << p_alpha_final[i] << ", ";
    // }
    // Rcout << "\nb: ";
    // for(int i=0 ; i<n_b ; ++i){
    //     Rcout << p_beta_final[i] << ", ";
    // }
    // Rcout << "\n\n";

    // keeping track of iterations
    iterations_all[v] += iter;

    // saving the fixef coefs
    double *fixef_values = args->fixef_values;
    if(args->save_fixef){

        for(int i=0 ; i<n_a ; ++i){
            fixef_values[i] += p_alpha_final[i];
        }

        for(int j=0 ; j<n_b ; ++j){
            fixef_values[n_a + j] += p_beta_final[j];
        }
    }

}

void compute_fe_gnl(double *p_fe_coef_origin, double *p_fe_coef_destination,
                    double *p_sum_other_means, double *p_sum_in_out, PARAM_DEMEAN *args){
    // update of the cluster coefficients
    // first we update mu, then we update the cluster coefficicents

    //
    // Loading the variables
    //

    int n_obs = args->n_obs;
    int Q = args->Q;

    // fe_info
    FEClass &FE_info = *(args->p_FE_info);

    // We update each cluster coefficient, starting from Q (the smallest one)

    std::fill_n(p_sum_other_means, n_obs, 0);

    // we start with Q-1
    for(int q=0 ; q<(Q-1) ; ++q){
        FE_info.add_wfe_coef_to_mu(q, p_fe_coef_origin, p_sum_other_means);
    }

    // Rcout << "Head of sum_other_means: ";
    // for(int i=0 ; i<10 ; ++i) Rcout << p_sum_other_means[i] << ", ";
    // Rcout << "\n";

    for(int q=Q-1 ; q>=0 ; q--){


        FE_info.compute_fe_coef(q, p_fe_coef_destination, p_sum_other_means, p_sum_in_out);

        // if(int q == 0){
        // 	Rprintf("p_fe_coef_destination: %.3f, %.3f, %.3f, %.3f\n", my_fe_coef[0], my_fe_coef[1], my_fe_coef[2], my_fe_coef[3]);
        // }


        // updating the value of p_sum_other_means (only if necessary)
        if(q != 0){

            // We recompute it from scratch (only way -- otherwise precision problems arise)

            std::fill_n(p_sum_other_means, n_obs, 0);

            double *my_fe_coef;
            for(int h=0 ; h<Q ; h++){
                if(h == q-1) continue;

                if(h < q-1){
                    my_fe_coef = p_fe_coef_origin;
                } else {
                    my_fe_coef = p_fe_coef_destination;
                }

                FE_info.add_wfe_coef_to_mu(h, my_fe_coef, p_sum_other_means);
            }

        }
    }

    // In the end, the array p_fe_coef_destination is fully updated, starting from Q to 1

}



bool demean_acc_gnl(int v, int iterMax, PARAM_DEMEAN *args){

    //
    // data
    //

    int n_obs = args->n_obs;
    int nb_coef_T = args->nb_coef_T;
    int Q = args->Q;
    double diffMax = args->diffMax;

    // fe_info
    FEClass &FE_info = *(args->p_FE_info);

    // input output
    vector<sVec> &p_input = args->p_input;
    vector<double*> &p_output = args->p_output;
    sVec &input = p_input[v];
    double *output = p_output[v];

    // temp var:
    vector<double> sum_other_means(n_obs);
    double *p_sum_other_means = sum_other_means.data();

    // conditional sum of input minus output
    vector<double> sum_input_output(nb_coef_T, 0);
    double *p_sum_in_out = sum_input_output.data();

    for(int q=0 ; q<Q ; ++q){
        FE_info.compute_in_out(q, p_sum_in_out, input, output);
    }

    // Rcout << "Sum in out:\n";
    // int c_next = 0;
    // int q_current = 0;
    // for(int c=0 ; c<nb_coef_T ; ++c){
    //     if(c == c_next){
    //         if(c != 0) Rcout << "\n";
    //
    //         Rcout << "q = " << q_current << ": ";
    //         c_next += FE_info.nb_coef_Q[q_current];
    //         q_current++;
    //     }
    //
    //     Rcout << p_sum_in_out[c] << ", ";
    // }
    // Rcout << "\n";



    // interruption handling
    bool isMaster = omp_get_thread_num() == 0;
    bool *pStopNow = args->stopnow;
    // I overcast to remember the lesson
    double flop = 4.0*(5 + 12*(Q-1) + 4*(Q-1)*(Q-1))*static_cast<double>(n_obs); // rough estimate nber operation per iter
    int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second

    //
    // IT iteration (preparation)
    //

    // variables on 1:K
    vector<double> X(nb_coef_T, 0);
    vector<double> GX(nb_coef_T);
    vector<double> GGX(nb_coef_T);
    // pointers:
    double *p_X   = X.data();
    double *p_GX  = GX.data();
    double *p_GGX = GGX.data();

    // variables on 1:(Q-1)
    int nb_coef_no_Q = 0;
    for(int q = 0 ; q<(Q-1) ; ++q){
        nb_coef_no_Q += FE_info.nb_coef_Q[q];
    }
    vector<double> delta_GX(nb_coef_no_Q);
    vector<double> delta2_X(nb_coef_no_Q);

    //
    // the main loop
    //

    // first iteration
    compute_fe_gnl(p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

    // Rcout << "Coef first iter:\n";
    // c_next = 0;
    // q_current = 0;
    // for(int c=0 ; c<nb_coef_T ; ++c){
    //     if(c == c_next){
    //         if(c != 0) Rcout << "\n";
    //
    //         Rcout << "q = " << q_current << ": ";
    //         c_next += FE_info.nb_coef_Q[q_current];
    //         q_current++;
    //     }
    //
    //     Rcout << p_GX[c] << ", ";
    // }
    // Rcout << "\n\n\n";

    // check whether we should go into the loop
    bool keepGoing = false;
    for(int i=0 ; i<nb_coef_T ; ++i){
        if(continue_crit(X[i], GX[i], diffMax)){
            keepGoing = true;
            break;
        }
    }

    // For the stopping criterion on total addition
    // vector<double> mu_last(n_obs, 0);
    // double input_mean = 0;
    double ssr = 0;

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
        compute_fe_gnl(p_GX, p_GGX, p_sum_other_means, p_sum_in_out, args);

        // X ; update of the cluster coefficient
        numconv = dm_update_X_IronsTuck(nb_coef_no_Q, X, GX, GGX, delta_GX, delta2_X);
        if(numconv) break;

        // GX -- origin: X, destination: GX
        compute_fe_gnl(p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

        keepGoing = false;
        for(int i=0 ; i<nb_coef_no_Q ; ++i){
            if(continue_crit(X[i], GX[i], diffMax)){
                keepGoing = true;
                break;
            }
        }

        // Other stopping criterion: change to SSR very small
        if(iter % 50 == 0){

            // mu_current is the vector of means
            vector<double> mu_current(n_obs, 0);
            double *p_mu = mu_current.data();
            for(int q=0 ; q<Q ; ++q){

                FE_info.add_fe_coef_to_mu(q, p_GX, p_mu);

            }

            // init ssr if iter == 50 / otherwise, comparison
            if(iter == 50){

                ssr = 0;
                double resid;
                for(int i=0 ; i<n_obs ; ++i){
                    resid = input[i] - mu_current[i];
                    ssr += resid*resid;
                }

            } else {
                double ssr_old = ssr;

                // we compute the new SSR
                ssr = 0;
                double resid;
                for(int i=0 ; i<n_obs ; ++i){
                    resid = input[i] - mu_current[i];
                    ssr += resid*resid;
                }

                // if(isMaster) Rprintf("iter %i -- SSR = %.0f (diff = %.0f)\n", iter, ssr, ssr_old - ssr);

                if(stopping_crit(ssr_old, ssr, diffMax)){
                    break;
                }

            }
        }


    }

    //
    // Updating the output
    //

    for(int q=0 ; q<Q ; ++q){
        FE_info.add_fe_coef_to_mu(q, p_GX, output);
    }

    // keeping track of iterations
    int *iterations_all = args->p_iterations_all;
    iterations_all[v] += iter;

    // saving the fixef coefs
    double *fixef_values = args->fixef_values;
    if(args->save_fixef){
        for(int m=0 ; m<nb_coef_T ; ++m){
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

        if(conv == false && iterMax > 15){
            // 2 convergence
            demean_acc_2(v, iterMax / 2 - 15, args);

            if(Q > 2){
                // re-acceleration
                demean_acc_gnl(v, iterMax / 2, args);
            }
        }
    }

    // if(args->p_iterations_all[v] > 50) balance_FEs(v, args);

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
List cpp_demean(SEXP y, SEXP X_raw, SEXP r_weights, int iterMax, double diffMax, SEXP r_nb_id_Q,
                   SEXP fe_id_list, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list,
                   SEXP r_init, int nthreads, bool save_fixef = false){
    // main fun that calls demean_single
    // preformats all the information needed on the fixed-effects
    // y: the dependent variable
    // X_raw: the matrix of the explanatory variables -- can be "empty"

    // when including weights: recreate table values
    // export weights and is_weight bool

    // slope_flag: whether a FE is a varying slope
    // slope_var: the associated variables with varying slopes

    //initial variables
    int Q = Rf_length(r_nb_id_Q);

    // info on y
    sMat m_y(y);
    int n_vars_y = m_y.ncol();
    bool useY = n_vars_y > 0;
    bool is_y_list = n_vars_y > 1 || TYPEOF(y) == VECSXP;
    int n_obs = m_y.nrow();

    // info on X
    sMat m_X(X_raw);
    int n_vars_X = m_X.ncol();
    if(useY == false){
        n_obs = m_X.nrow();
    }
    bool useX = n_vars_X > 0;

    int n_vars = n_vars_y + n_vars_X;

    // initialisation if needed (we never initialize when only one FE, except if asked explicitly)
    bool isInit = Rf_xlength(r_init) != 1 && Q > 1;
    double *init = REAL(r_init);
    bool saveInit = ((isInit || init[0] != 0) && Q > 1) || init[0] == 666;

    // Creating the object containing all information on the FEs
    FEClass FE_info(n_obs, Q, r_weights, fe_id_list, r_nb_id_Q, table_id_I, slope_flag_Q, slope_vars_list);
    int nb_coef_T = FE_info.nb_coef_T;

    // output vector: (Note that if the means are provided, we use that vector and will modify it in place)
    int64_t n_total = static_cast<int64_t>(n_obs) * n_vars;
    SEXP output_values = PROTECT(Rf_allocVector(REALSXP, isInit ? 1 : n_total));

    double *p_output_origin = isInit ? init : REAL(output_values);
    if(isInit == false){
        std::fill_n(p_output_origin, n_total, 0);
    }

    //
    // vector of pointers: input/output
    //

    vector<double*> p_output(n_vars);
    p_output[0] = p_output_origin;
    for(int v=1 ; v<n_vars ; v++){
        p_output[v] = p_output[v - 1] + n_obs;
    }

    vector<sVec> p_input(n_vars);

    for(int k=0 ; k<n_vars_X ; ++k){
        p_input[k] = m_X[k];
    }

    for(int k=0 ; k<n_vars_y ; ++k){
        p_input[n_vars_X + k] = m_y[k];
    }

    // keeping track of iterations
    vector<int> iterations_all(n_vars, 0);
    int *p_iterations_all = iterations_all.data();

    // save fixef option
    if(useX && save_fixef){
        stop("save_fixef can be used only when there is no Xs.");
    }

    vector<double> fixef_values(save_fixef ? nb_coef_T : 1, 0);
    double *p_fixef_values = fixef_values.data();

    //
    // Sending variables to envir
    //

    PARAM_DEMEAN args;

    args.n_obs = n_obs;
    args.iterMax = iterMax;
    args.diffMax = diffMax;
    args.Q = Q;
    args.nb_coef_T = nb_coef_T;
    args.p_input = p_input;
    args.p_output = p_output;
    args.p_iterations_all = p_iterations_all;

    // save FE_info
    args.p_FE_info = &FE_info;

    // save fixef:
    args.save_fixef = save_fixef;
    args.fixef_values = p_fixef_values;

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

    int nthreads_current = nthreads > n_vars ? n_vars : nthreads;

    // enlever les rprintf dans les nthreads jobs
    #pragma omp parallel for num_threads(nthreads_current) schedule(static, 1)
    for(int v = 0 ; v<(n_vars+nthreads_current) ; ++v){
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
        UNPROTECT(1);
        stop("cpp_demean: User interrupt.");
    }

    // Rprintf("Master checking: %i\n", *pcounter);

    //
    // save
    //

    List res; // a vector and a matrix

    int nrow = useX ? n_obs : 1;
    int ncol = useX ? n_vars_X : 1;
    NumericMatrix X_demean(nrow, ncol);

    sVec p_input_tmp;
    double *p_output_tmp;
    for(int k=0 ; k<n_vars_X ; ++k){
        p_input_tmp = p_input[k];
        p_output_tmp = p_output[k];

        for(int i=0 ; i < n_obs ; ++i){
            X_demean(i, k) = p_input_tmp[i] - p_output_tmp[i];
        }
    }

    res["X_demean"] = X_demean;



    if(is_y_list && useY){
        List y_demean(n_vars_y);

        for(int v=0 ; v<n_vars_y ; ++v){
            p_input_tmp  = p_input[n_vars_X + v];
            p_output_tmp = p_output[n_vars_X + v];

            NumericVector y_demean_tmp(n_obs);
            for(int i=0 ; i < n_obs ; ++i){
                y_demean_tmp[i] = p_input_tmp[i] - p_output_tmp[i];
            }

            y_demean[v] = y_demean_tmp;
        }

        res["y_demean"] = y_demean;

    } else {
        NumericVector y_demean(useY ? n_obs : 1);
        if(useY){
            // y is always the last variable
            p_input_tmp = p_input[n_vars - 1];
            p_output_tmp = p_output[n_vars - 1];
            for(int i=0 ; i < n_obs ; ++i){
                y_demean[i] = p_input_tmp[i] - p_output_tmp[i];
            }
        }
        res["y_demean"] = y_demean;
    }



    // iterations
    IntegerVector iter_final(n_vars);
    for(int v=0 ; v<n_vars ; ++v){
        iter_final[v] = p_iterations_all[v];
    }

    res["iterations"] = iter_final;

    // if save is requested
    if(saveInit){
        if(isInit){
            res["means"] = r_init;
        } else {
            res["means"] = output_values;
        }
    } else {
        res["means"] = 0.0;
    }

    // save fixef coef
    if(save_fixef){
        res["fixef_coef"] = fixef_values;
    }

    UNPROTECT(1);

    return(res);
}

//
// Next version => clean c++ code, use only sMat, create file with common functions
//


std::vector<int> set_parallel_scheme_ter(int N, int nthreads){
    // => this concerns only the parallel application on a 1-Dimensional matrix
    // takes in the nber of observations of the vector and the nber of threads
    // gives back a vector of the length the nber of threads + 1 giving the start/stop of each threads

    std::vector<int> res(nthreads + 1, 0);
    double N_rest = N;

    for(int i=0 ; i<nthreads ; ++i){
        res[i + 1] = ceil(N_rest / (nthreads - i));
        N_rest -= res[i + 1];
        res[i + 1] += res[i];
    }

    return res;
}


// [[Rcpp::export]]
List cpp_which_na_inf(SEXP x, int nthreads){
    // x: vector, matrix, data.frame // double or integer

    /*
     This function takes a matrix and looks at whether it contains NA or infinite values
     return: flag for na/inf + logical vector of obs that are Na/inf
     std::isnan, std::isinf are OK since cpp11 required
     do_any_na_inf: if high suspicion of NA present: we go directly constructing the vector is_na_inf
     in the "best" case (default expected), we need not construct is_na_inf
     */

    sMat mat(x);

    int nobs = mat.nrow();
    int K = mat.ncol();
    bool anyNAInf = false;
    bool any_na = false;    // return value
    bool any_inf = false;   // return value

    /*
     we make parallel the anyNAInf loop
     why? because we want that when there's no NA (default) it works as fast as possible
     if there are NAs, single threaded mode is faster, but then we circumvent with the do_any_na_inf flag
     */

    // no need to care about the race condition
    // "trick" to make a break in a multi-threaded section

    std::vector<int> bounds = set_parallel_scheme_ter(nobs, nthreads);

    #pragma omp parallel for num_threads(nthreads)
    for(int t=0 ; t<nthreads ; ++t){
        for(int k=0 ; k<K ; ++k){
            for(int i=bounds[t]; i<bounds[t + 1] && !anyNAInf ; ++i){

                if(mat[k].is_int){
                    if(mat(i, k) == -2147483648.0){
                        anyNAInf = true;
                    }
                } else if(std::isnan(mat(i, k)) || std::isinf(mat(i, k))){
                    anyNAInf = true;
                }
            }
        }
    }

    // object to return: is_na_inf
    LogicalVector is_na_inf(anyNAInf ? nobs : 1);

    if(anyNAInf){
        #pragma omp parallel for num_threads(nthreads)
        for(int i=0 ; i<nobs ; ++i){
            double x_tmp = 0;
            for(int k=0 ; k<K ; ++k){
                x_tmp = mat(i, k);
                if(mat[k].is_int){
                    if(mat(i, k) == -2147483648.0){
                        is_na_inf[i] = true;
                        any_na = true;
                        break;
                    }
                } else if(std::isnan(x_tmp)){
                    is_na_inf[i] = true;
                    any_na = true;
                    break;
                } else if(std::isinf(x_tmp)){
                    is_na_inf[i] = true;
                    any_inf = true;
                    break;
                }
            }
        }
    }

    // Return
    List res;
    res["any_na"] = any_na;
    res["any_inf"] = any_inf;
    res["any_na_inf"] = any_na || any_inf;
    res["is_na_inf"] = is_na_inf;

    return res;
}

