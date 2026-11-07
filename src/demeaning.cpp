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

#include "fixest_main.h"
#include <math.h>

using namespace Rcpp;
using std::vector;
using std::int64_t;


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


//
// Demeans each variable in input
// The method is based on obtaining the optimal cluster coefficients
//

// for interruption
int pending_interrupt() {
  return !(R_ToplevelExec(Rcpp::checkInterruptFn, NULL));
}

inline int omp_num_tasks(int64_t n, int64_t min_grain = 32768){
#ifdef _OPENMP
  int nt = 4 * omp_get_num_threads();
  int64_t max_t = n / min_grain;
  if(max_t < 1) max_t = 1;
  return max_t < nt ? static_cast<int>(max_t) : nt;
#else
  return 1;
#endif
}
// Works but only the master thread can call that function!
// What happens if the master thread has finished its job but the lower thread is in an "infinite" loop?
// this is tricky, as such we cannot stop it.
// Solution: create a function keeping the threads idle waiting for the complete job to be done
// BUT I need to add static allocation of threads => performance cost





//
// Now we start a big chunk => computing the varying slopes coefficients
// That's a big job. To simplify it, I created the class FEClass that takes care of it.
//


class simple_mat_with_id{
  // => Access to one of the n_coef matrices of size n_vs x n_vs; all stacked in a single vector
  //
  // px0: origin of the vector (which is of length n_coef * n_vs * n_vs)
  // n: n_vs
  // n2: explicit

  simple_mat_with_id() = delete;

  double *px0;
  int nrow, ncol, n_total;

public:
  simple_mat_with_id(double* px_in, int nrow_in):
    px0(px_in), nrow(nrow_in), ncol(nrow_in), n_total(nrow * ncol) {};
  simple_mat_with_id(double* px_in, int nrow_in, int ncol_in):
    px0(px_in), nrow(nrow_in), ncol(ncol_in), n_total(nrow * ncol) {};
  double& operator()(int id, int i, int j);
  double& operator()(int id, int i);
};

inline double& simple_mat_with_id::operator()(int id, int i, int j){
  return px0[n_total * id + i + nrow * j];
}

inline double& simple_mat_with_id::operator()(int id, int i){
  return px0[n_total * id + i];
}

class FEClass{

  int Q;
  int n_obs;
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
  vector<RealVec> p_vs_vars;
  double *p_weights = nullptr;

  vector<bool> is_slope_Q;
  vector<bool> is_slope_fe_Q;

  vector<int> nb_vs_Q;
  vector<int> nb_vs_noFE_Q;

  int *nb_id_Q;

  vector<int> coef_start_Q;

  // internal functions
  void compute_fe_coef_internal(int, double *, bool, RealVec, double *, double *);
  void compute_fe_coef_2_internal(double *, double *, double *, bool);
  void add_wfe_coef_to_mu_internal(int, double *, double *, bool);

public:

  // Utility class: Facilitates the access to the VS variables
  class simple_mat_of_vs_vars{
    int K_fe;
    vector<RealVec> pvars;

  public:
    simple_mat_of_vs_vars(const FEClass*, int);
    double operator()(int, int);
  };

  int nb_coef_T;
  vector<int> nb_coef_Q;

  // constructor:
  FEClass(int n_obs, int Q, SEXP r_weights, SEXP fe_id_list, SEXP r_nb_id_Q, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list);

  // functions
  void compute_fe_coef(double *fe_coef, RealVec &mu_in_N);
  void compute_fe_coef(int q, double *fe_coef, double *sum_other_coef_N, double *in_out_C);

  void add_wfe_coef_to_mu(int q, double *fe_coef_C, double *out_N);
  void add_fe_coef_to_mu(int q, double *fe_coef_C, double *out_N);
  void add_all_fe_coef_to_mu(double *fe_coef_C, double *out_N);
  void compute_sum_other_means(int q, double *fe_coef_origin, double *fe_coef_destination, double *out_N);
  double compute_ssr(double *fe_coef_C, RealVec &in_N);
  double compute_2_fe_ssr(double *fe_coef_a, double *fe_coef_b, double *in_out_C, RealVec &in_N);

  void compute_fe_coef_2(double *fe_coef_in_C, double *fe_coef_out_C, double *fe_coef_tmp, double *in_out_C);

  void add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b, double *in_out_C, double *out_N, bool update_beta);

  void compute_in_out(int q, double *in_out_C, RealVec &in_N, double *out_N);
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

    RealMat m_slopes(slope_vars_list);
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
void FEClass::compute_fe_coef(double *fe_coef_C, RealVec &mu_in_N){
  // mu: length n_obs, vector giving sum_in_out
  // fe_coef: vector receiving the fixed-effects coefficients

  // single FE version
  compute_fe_coef_internal(0, fe_coef_C, true, mu_in_N, nullptr, nullptr);
}

void FEClass::compute_fe_coef(int q, double *fe_coef_C, double *sum_other_coef_N, double *in_out_C){
  // mu: length n_obs, vector giving sum_in_out
  // fe_coef: vector receiving the fixed-effects coefficients

  // multiple FE version
  compute_fe_coef_internal(q, fe_coef_C, false, nullptr, sum_other_coef_N, in_out_C);

}

void FEClass::compute_fe_coef_internal(int q, double *fe_coef_C, bool is_single, RealVec mu_in_N, 
                                       double *sum_other_coef_N, double *in_out_C){
  // mu: length n_obs, vector giving sum_in_out
  // fe_coef: vector receiving the fixed-effects coefficients

  int V = nb_vs_Q[q];
  int *my_fe = p_fe_id[q];
  int nb_coef = nb_coef_Q[q];
  int nb_id = nb_id_Q[q];
  const bool is_weight = this->is_weight;
  const double *const p_weights = this->p_weights;

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

      const double *sum_in_out = in_out_C + coef_start_Q[q];

      // initialize cluster coef
      #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef)) default(shared) if(nb_coef >= 32768)
      for(int m=0 ; m<nb_coef ; ++m){
        my_fe_coef[m] = sum_in_out[m];
      }

      // looping sequentially over the sum of other coefficients
      for(int i=0 ; i<n_obs ; ++i){
        my_fe_coef[my_fe[i] - 1] -= sum_other_coef_N[i];
      }

    }
    
    #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef)) default(shared) if(nb_coef >= 32768)
    for(int m=0 ; m<nb_coef ; ++m){
      my_fe_coef[m] /= my_SW[m];
    }

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
      const double *sum_in_out = in_out_C + coef_start_Q[q];
      #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef)) default(shared) if(nb_coef >= 32768)
      for(int m=0 ; m<nb_coef ; ++m){
        my_fe_coef[m] = sum_in_out[m];
      }

      for(int i=0 ; i<n_obs ; ++i){
        for(int v=0 ; v<V ; ++v){
          rhs(my_fe[i] - 1, v) -= VS_mat(i, v) * sum_other_coef_N[i];
        }
      }

    }


    // Forward elimination and backward solve per group c
    #pragma omp taskloop num_tasks(omp_num_tasks(nb_id, 16384)) default(shared) if(nb_id >= 16384)
    for(int c=0 ; c<nb_id ; ++c){
      for(int v=0 ; v<V ; ++v){
        for(int v1=v + 1 ; v1<V ; ++v1){
          // -= row_coef * value from the v row
          rhs(c, v1) -= my_system(c, v1, v) * rhs(c, v);
        }
      }

      // We backward solve
      for(int v=V - 1 ; v>=0 ; --v){
        if(my_system(c, v, v) == 0){
          rhs(c, v) = 0;
        } else {
          double val = rhs(c, v);
          for(int v_done=v + 1 ; v_done<V ; ++v_done){
            val -= rhs(c, v_done) * my_system(c, v, v_done);
          }
          rhs(c, v) = val / my_system(c, v, v);
        }
      }
    }

  }

}


void FEClass::compute_fe_coef_2_internal(double *fe_coef_in_out_C, double *fe_coef_tmp, 
                                         double *in_out_C, bool step_2 = false){
  // The aim here is to update the value of b given the value of a
  // fe_coef_tmp is always of size the 2nd FE
  // fe_in/out are always the size the 1st FE
  // a: origin
  // b: destination

  int index_a = 0;
  int index_b = 1;

  double *my_fe_coef_a;
  double *my_fe_coef_b;
  
  const bool is_weight = this->is_weight;

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

  int nb_coef_b = nb_coef_Q[index_b];

  int nb_id_b = nb_id_Q[index_b];

  double *my_in_out_b = in_out_C + coef_start_Q[index_b];

  const bool is_slope_a = is_slope_Q[index_a];
  const bool is_slope_b = is_slope_Q[index_b];

  // VSlope utilities
  simple_mat_of_vs_vars VS_mat_a(this, index_a);
  simple_mat_with_id my_vs_coef_a(my_fe_coef_a, V_a, 1);

  simple_mat_of_vs_vars VS_mat_b(this, index_b);
  simple_mat_with_id my_vs_coef_b(my_fe_coef_b, V_b, 1);


  // initialize cluster coef
  #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef_b)) default(shared) if(nb_coef_b >= 32768)
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
    #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef_b)) default(shared) if(nb_coef_b >= 32768)
    for(int m=0 ; m<nb_coef_b ; ++m){
      my_fe_coef_b[m] /= my_SW[m];
    }

  } else {
    simple_mat_with_id my_system_b(p_eq_systems_VS_C[index_b], V_b);
    simple_mat_with_id rhs_b(my_fe_coef_b, V_b, 1);

    double tmp = 0;

    for(int i=0 ; i<n_obs ; ++i){
      
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

    // Forward elimination and backward solve per group c
    #pragma omp taskloop num_tasks(omp_num_tasks(nb_id_b, 16384)) default(shared) if(nb_id_b >= 16384)
    for(int c=0 ; c<nb_id_b ; ++c){
      for(int v=0 ; v<V_b ; ++v){
        for(int v1=v + 1 ; v1<V_b ; ++v1){
          // -= row_coef * value from the v row
          rhs_b(c, v1) -= my_system_b(c, v1, v) * rhs_b(c, v);
        }
      }

      // We backward solve
      for(int v=V_b - 1 ; v>=0 ; --v){
        if(my_system_b(c, v, v) == 0){
          rhs_b(c, v) = 0;
        } else {
          double val = rhs_b(c, v);
          for(int v_done=v + 1 ; v_done<V_b ; ++v_done){
            val -= rhs_b(c, v_done) * my_system_b(c, v, v_done);
          }
          rhs_b(c, v) = val / my_system_b(c, v, v);
        }
      }
    }
  }

}

void FEClass::compute_fe_coef_2(double *fe_coef_in_C, double *fe_coef_out_C, 
                                double *fe_coef_tmp, double *in_out_C){
  // Specific to the 2-FEs case
  // This way we avoid creating and using a temp object of length N



  //
  // Step 1: Updating b
  //

  compute_fe_coef_2_internal(fe_coef_in_C, fe_coef_tmp, in_out_C);

  //
  // Step 2: Updating a
  //

  compute_fe_coef_2_internal(fe_coef_out_C, fe_coef_tmp, in_out_C, true);

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

  const bool use_weights = add_weights && is_weight;

  if(is_slope_Q[q] == false){

    #pragma omp taskloop num_tasks(omp_num_tasks(n_obs)) default(shared) if(n_obs >= 20000)
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

    #pragma omp taskloop num_tasks(omp_num_tasks(n_obs)) default(shared) if(n_obs >= 20000)
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

void FEClass::compute_sum_other_means(int q, double *fe_coef_origin, double *fe_coef_destination, double *out_N){
  // Fused Step 1 of compute_fe_gnl: sums all h != q into a register per observation i
  // and writes out_N[i] in a single pass instead of Q taskloop passes over DRAM.
  const bool use_weights = is_weight;
  vector<double*> my_fe_coef_H(Q);
  vector<simple_mat_of_vs_vars> vs_mat_H;
  vs_mat_H.reserve(Q);
  vector<simple_mat_with_id> vs_coef_H;
  vs_coef_H.reserve(Q);

  for(int h=0 ; h<Q ; ++h){
    double *base_coef = (h < q) ? fe_coef_origin : fe_coef_destination;
    my_fe_coef_H[h] = base_coef + coef_start_Q[h];
    vs_mat_H.emplace_back(this, h);
    vs_coef_H.emplace_back(my_fe_coef_H[h], is_slope_Q[h] ? nb_vs_Q[h] : 1, 1);
  }

  #pragma omp taskloop num_tasks(omp_num_tasks(n_obs)) default(shared) if(n_obs >= 20000)
  for(int i=0 ; i<n_obs ; ++i){
    double sum_i = 0.0;
    const double w_i = use_weights ? p_weights[i] : 1.0;
    for(int h=0 ; h<Q ; ++h){
      if(h == q) continue;
      int id = p_fe_id[h][i] - 1;
      if(is_slope_Q[h] == false){
        if(use_weights){
          sum_i += my_fe_coef_H[h][id] * w_i;
        } else {
          sum_i += my_fe_coef_H[h][id];
        }
      } else {
        int V = nb_vs_Q[h];
        simple_mat_of_vs_vars &VS_mat = vs_mat_H[h];
        simple_mat_with_id &my_vs_coef = vs_coef_H[h];
        for(int v=0 ; v<V ; ++v){
          if(use_weights){
            sum_i += my_vs_coef(id, v) * VS_mat(i, v) * w_i;
          } else {
            sum_i += my_vs_coef(id, v) * VS_mat(i, v);
          }
        }
      }
    }
    out_N[i] = sum_i;
  }
}

void FEClass::add_all_fe_coef_to_mu(double *fe_coef_C, double *out_N){
  // Fused pass adding all Q FE contributions to out_N[i] in a single taskloop
  vector<double*> my_fe_coef_Q(Q);
  vector<simple_mat_of_vs_vars> vs_mat_Q;
  vs_mat_Q.reserve(Q);
  vector<simple_mat_with_id> vs_coef_Q;
  vs_coef_Q.reserve(Q);

  for(int q=0 ; q<Q ; ++q){
    my_fe_coef_Q[q] = fe_coef_C + coef_start_Q[q];
    vs_mat_Q.emplace_back(this, q);
    vs_coef_Q.emplace_back(my_fe_coef_Q[q], is_slope_Q[q] ? nb_vs_Q[q] : 1, 1);
  }

  #pragma omp taskloop num_tasks(omp_num_tasks(n_obs)) default(shared) if(n_obs >= 20000)
  for(int i=0 ; i<n_obs ; ++i){
    double sum_i = 0.0;
    for(int q=0 ; q<Q ; ++q){
      int id = p_fe_id[q][i] - 1;
      if(is_slope_Q[q] == false){
        sum_i += my_fe_coef_Q[q][id];
      } else {
        int V = nb_vs_Q[q];
        simple_mat_of_vs_vars &VS_mat = vs_mat_Q[q];
        simple_mat_with_id &my_vs_coef = vs_coef_Q[q];
        for(int v=0 ; v<V ; ++v){
          sum_i += my_vs_coef(id, v) * VS_mat(i, v);
        }
      }
    }
    out_N[i] += sum_i;
  }
}

double FEClass::compute_ssr(double *fe_coef_C, RealVec &in_N){
  // Computes sum_{i=0}^{n_obs-1} (in_N[i] - sum_q FE_q(i))^2 in a single fused taskloop
  // without allocating a temporary n_obs vector.
  vector<double*> my_fe_coef_Q(Q);
  vector<simple_mat_of_vs_vars> vs_mat_Q;
  vs_mat_Q.reserve(Q);
  vector<simple_mat_with_id> vs_coef_Q;
  vs_coef_Q.reserve(Q);

  for(int q=0 ; q<Q ; ++q){
    my_fe_coef_Q[q] = fe_coef_C + coef_start_Q[q];
    vs_mat_Q.emplace_back(this, q);
    vs_coef_Q.emplace_back(my_fe_coef_Q[q], is_slope_Q[q] ? nb_vs_Q[q] : 1, 1);
  }

  int n_tasks = omp_num_tasks(n_obs);
  vector<double> chunk_ssr(n_tasks, 0.0);

  #pragma omp taskloop num_tasks(n_tasks) default(shared) if(n_obs >= 20000)
  for(int t=0 ; t<n_tasks ; ++t){
    int i_start = static_cast<int>((static_cast<int64_t>(t) * n_obs) / n_tasks);
    int i_end   = static_cast<int>((static_cast<int64_t>(t + 1) * n_obs) / n_tasks);
    double local_ssr = 0.0;
    for(int i=i_start ; i<i_end ; ++i){
      double mu_i = 0.0;
      for(int q=0 ; q<Q ; ++q){
        int id = p_fe_id[q][i] - 1;
        if(is_slope_Q[q] == false){
          mu_i += my_fe_coef_Q[q][id];
        } else {
          int V = nb_vs_Q[q];
          simple_mat_of_vs_vars &VS_mat = vs_mat_Q[q];
          simple_mat_with_id &my_vs_coef = vs_coef_Q[q];
          for(int v=0 ; v<V ; ++v){
            mu_i += my_vs_coef(id, v) * VS_mat(i, v);
          }
        }
      }
      double resid = in_N[i] - mu_i;
      local_ssr += resid * resid;
    }
    chunk_ssr[t] = local_ssr;
  }

  double ssr = 0.0;
  for(int t=0 ; t<n_tasks ; ++t){
    ssr += chunk_ssr[t];
  }
  return ssr;
}

double FEClass::compute_2_fe_ssr(double *fe_coef_a, double *fe_coef_b, double *in_out_C, RealVec &in_N){
  compute_fe_coef_2_internal(fe_coef_a, fe_coef_b, in_out_C, false);

  simple_mat_of_vs_vars VS_mat_0(this, 0);
  simple_mat_of_vs_vars VS_mat_1(this, 1);
  simple_mat_with_id my_vs_coef_0(fe_coef_a, is_slope_Q[0] ? nb_vs_Q[0] : 1, 1);
  simple_mat_with_id my_vs_coef_1(fe_coef_b, is_slope_Q[1] ? nb_vs_Q[1] : 1, 1);
  int *my_fe_0 = p_fe_id[0];
  int *my_fe_1 = p_fe_id[1];
  bool is_slope_0 = is_slope_Q[0];
  bool is_slope_1 = is_slope_Q[1];
  int V_0 = nb_vs_Q[0];
  int V_1 = nb_vs_Q[1];

  int n_tasks = omp_num_tasks(n_obs);
  vector<double> chunk_ssr(n_tasks, 0.0);

  #pragma omp taskloop num_tasks(n_tasks) default(shared) if(n_obs >= 20000)
  for(int t=0 ; t<n_tasks ; ++t){
    int i_start = static_cast<int>((static_cast<int64_t>(t) * n_obs) / n_tasks);
    int i_end   = static_cast<int>((static_cast<int64_t>(t + 1) * n_obs) / n_tasks);
    double local_ssr = 0.0;
    for(int i=i_start ; i<i_end ; ++i){
      double mu_i = 0.0;
      if(is_slope_0){
        for(int v=0 ; v<V_0 ; ++v) mu_i += my_vs_coef_0(my_fe_0[i] - 1, v) * VS_mat_0(i, v);
      } else {
        mu_i += fe_coef_a[my_fe_0[i] - 1];
      }
      if(is_slope_1){
        for(int v=0 ; v<V_1 ; ++v) mu_i += my_vs_coef_1(my_fe_1[i] - 1, v) * VS_mat_1(i, v);
      } else {
        mu_i += fe_coef_b[my_fe_1[i] - 1];
      }
      double resid = in_N[i] - mu_i;
      local_ssr += resid * resid;
    }
    chunk_ssr[t] = local_ssr;
  }

  double ssr = 0.0;
  for(int t=0 ; t<n_tasks ; ++t){
    ssr += chunk_ssr[t];
  }
  return ssr;
}

void FEClass::add_2_fe_coef_to_mu(double *fe_coef_a, double *fe_coef_b, double *in_out_C, 
                                  double *out_N, bool update_beta = true){
  // We add the value of the FE coefficients to each observation


  //
  // Step 1: we update the coefficients of b
  //

  if(update_beta){
    compute_fe_coef_2_internal(fe_coef_a, fe_coef_b, in_out_C, false);
  }

  //
  // Step 2: we add the value of each coef in a single fused taskloop
  //

  simple_mat_of_vs_vars VS_mat_0(this, 0);
  simple_mat_of_vs_vars VS_mat_1(this, 1);
  simple_mat_with_id my_vs_coef_0(fe_coef_a, is_slope_Q[0] ? nb_vs_Q[0] : 1, 1);
  simple_mat_with_id my_vs_coef_1(fe_coef_b, is_slope_Q[1] ? nb_vs_Q[1] : 1, 1);
  int *my_fe_0 = p_fe_id[0];
  int *my_fe_1 = p_fe_id[1];
  bool is_slope_0 = is_slope_Q[0];
  bool is_slope_1 = is_slope_Q[1];
  int V_0 = nb_vs_Q[0];
  int V_1 = nb_vs_Q[1];

  #pragma omp taskloop num_tasks(omp_num_tasks(n_obs)) default(shared) if(n_obs >= 20000)
  for(int i=0 ; i<n_obs ; ++i){
    double sum_i = 0.0;
    if(is_slope_0){
      for(int v=0 ; v<V_0 ; ++v){
        sum_i += my_vs_coef_0(my_fe_0[i] - 1, v) * VS_mat_0(i, v);
      }
    } else {
      sum_i += fe_coef_a[my_fe_0[i] - 1];
    }
    if(is_slope_1){
      for(int v=0 ; v<V_1 ; ++v){
        sum_i += my_vs_coef_1(my_fe_1[i] - 1, v) * VS_mat_1(i, v);
      }
    } else {
      sum_i += fe_coef_b[my_fe_1[i] - 1];
    }
    out_N[i] += sum_i;
  }

}


void FEClass::compute_in_out(int q, double *in_out_C, RealVec &in_N, double *out_N){
  // q: the index of the FE
  // *in_out_C: pointer to the vector of cumulated (input-output), of length C (nber of FE coefs)
  //            => this is a precomputation of the target quantity
  // in_N: vector of input, of length N
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
  
  int algo_extraProj;
  int algo_iter_warmup;
  int algo_iter_projAfterAcc;
  int algo_iter_grandAcc;

  // iterations
  int *p_iterations_all;

  // vectors of pointers
  vector<RealVec> p_input;
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

  #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef_no_Q)) default(shared) if(nb_coef_no_Q >= 20000)
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
    #pragma omp taskloop num_tasks(omp_num_tasks(nb_coef_no_Q)) default(shared) if(nb_coef_no_Q >= 20000)
    for(int i=0 ; i<nb_coef_no_Q ; ++i){
      X[i] = GGX[i] - coef * delta_GX[i];
    }
  }

  return(res);
}


void compute_fe_gnl(double *p_fe_coef_origin, double *p_fe_coef_destination,
                    double *p_sum_other_means, double *p_sum_in_out, PARAM_DEMEAN *args){
  // update of the cluster coefficients
  // p_fe_coef_origin: the FEs vector of origin
  // p_fe_coef_destination: the FEs of destination (which will receive the updated values)
  // p_sum_other_means: vector of length N, which will be equal to the sum of the FEs but on
  //   one dimension
  // p_sum_in_out: vector of length C (# of FE coefs) which contains the precomputed sum of 
  //   the target variable per FE
  //
  // first we update mu, then we update the cluster coefficicents

  //
  // Loading the variables
  //

  int Q = args->Q;

  // fe_info
  FEClass &FE_info = *(args->p_FE_info);

  // We update each cluster coefficient, starting from Q (the smallest one)

  for(int q=Q-1 ; q>=0 ; q--){
    
    // STEP 1: we compute the sum of the other coef values in a single fused taskloop
    FE_info.compute_sum_other_means(q, p_fe_coef_origin, p_fe_coef_destination, p_sum_other_means);

    // STEP 2: computing the FE coef
    FE_info.compute_fe_coef(q, p_fe_coef_destination, p_sum_other_means, p_sum_in_out);
    
  }

  // In the end, the array p_fe_coef_destination is fully updated, starting from Q to 1

}


void compute_fe(int Q, double *p_fe_coef_origin, double *p_fe_coef_destination,
                double *p_sum_other_means, double *p_sum_in_out, 
                PARAM_DEMEAN *args){
  if(Q == 2){
    FEClass &FE_info = *(args->p_FE_info);
    FE_info.compute_fe_coef_2(p_fe_coef_origin, p_fe_coef_destination, 
                              p_sum_other_means, p_sum_in_out);
  } else {
    compute_fe_gnl(p_fe_coef_origin, p_fe_coef_destination, p_sum_other_means, 
                   p_sum_in_out, args);
  }
}


void demean_single_1(int v, PARAM_DEMEAN* args){
  // v: variable identifier to demean

  // Q == 1: nothing to say, just compute the closed form

  // loading the data
  int nb_coef_T = args->nb_coef_T;

  vector<RealVec> &p_input = args->p_input;
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
  RealVec &input = p_input[v];
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


bool demean_acc_gnl(int v, int iterMax, PARAM_DEMEAN *args, bool two_fe = false){

  //
  // data
  //
  
  // fe_info
  FEClass &FE_info = *(args->p_FE_info);
  
  // algo info
  const int n_extraProj = args->algo_extraProj;
  const int iter_projAfterAcc = args->algo_iter_projAfterAcc;
  const int iter_grandAcc = args->algo_iter_grandAcc;

  int n_obs = args->n_obs;
  int nb_coef_T = args->nb_coef_T;
  int Q = args->Q;
  double diffMax = args->diffMax;
  
  int nb_coef_all = nb_coef_T;
  const bool two_fe_algo = two_fe || Q == 2;
  if(two_fe_algo){
    Q = 2;
    // special case: 2 FEs, below will contain only the size of the first FE
    //   this is because 2-FE is a special case which allows to avoid creating 
    //   a N-size vector
    nb_coef_T = FE_info.nb_coef_Q[0];
    // nb_coef_all: only used in in-out
    nb_coef_all = FE_info.nb_coef_Q[0] + FE_info.nb_coef_Q[1];
  }  

  // input output
  vector<RealVec> &p_input = args->p_input;
  vector<double*> &p_output = args->p_output;
  RealVec &input = p_input[v];
  double *output = p_output[v];

  // temp var:
  int size_other_means = two_fe_algo ? FE_info.nb_coef_Q[1] : n_obs;
  vector<double> sum_other_means_or_second_coef(size_other_means);
  double *p_sum_other_means = sum_other_means_or_second_coef.data();

  // conditional sum of input minus output
  vector<double> sum_input_output(nb_coef_all, 0);
  double *p_sum_in_out = sum_input_output.data();

  for(int q=0 ; q<Q ; ++q){
    FE_info.compute_in_out(q, p_sum_in_out, input, output);
  }

  // interruption handling
  bool isMaster = omp_get_thread_num() == 0;
  bool *pStopNow = args->stopnow;
  // I overcast to remember the lesson
  // rough estimate nber operation per iter
  double flop = 4.0*(5 + 12*(Q-1) + 4*(Q-1)*(Q-1))*static_cast<double>(n_obs); 
  if(two_fe_algo){
    flop = 20.0 * static_cast<double>(n_obs);
  }
  int iterSecond = ceil(2000000000 / flop / 5); // nber iter per 1/5 second at 2GHz

  //
  // IT iteration (preparation)
  //

  // variables on 1:nb_coef
  // note that in the case of two_fe_algo, the length is equal to the 1st FE coefs
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
  
  // additional vectors used to save the coefs
  vector<double> Y(nb_coef_T);
  vector<double> GY(nb_coef_T);
  vector<double> GGY(nb_coef_T);
  // pointers:
  double *p_Y   = Y.data();
  double *p_GY  = GY.data();
  double *p_GGY = GGY.data();
  
  int grand_acc = 0;

  //
  // the main loop
  //

  // first iteration
  compute_fe(Q, p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

  // check whether we should go into the loop
  bool keepGoing = false;
  for(int i=0 ; i<nb_coef_T ; ++i){
    if(continue_crit(X[i], GX[i], diffMax)){
      keepGoing = true;
      break;
    }
  }

  // For the stopping criterion on total addition
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
    
    for(int rep=0 ; rep<n_extraProj ; ++rep){
      // simple projections, at the request of the user
      // may be useful to hasten convergence on special cases
      // default is 0
      compute_fe(Q, p_GX, p_GGX, p_sum_other_means, p_sum_in_out, args);
      compute_fe(Q, p_GGX, p_X, p_sum_other_means, p_sum_in_out, args);
      compute_fe(Q, p_X, p_GX, p_sum_other_means, p_sum_in_out, args);
    }

    // origin: GX, destination: GGX
    compute_fe(Q, p_GX, p_GGX, p_sum_other_means, p_sum_in_out, args);

    // X: outcome of the acceleration
    numconv = dm_update_X_IronsTuck(nb_coef_no_Q, X, GX, GGX, delta_GX, delta2_X);
    if(numconv) break;
    
    if(iter >= iter_projAfterAcc){
      std::memcpy(p_Y, p_X, nb_coef_T * sizeof(double));
      compute_fe(Q, p_Y, p_X, p_sum_other_means, p_sum_in_out, args);
    }    

    // origin: X, destination: GX
    compute_fe(Q, p_X, p_GX, p_sum_other_means, p_sum_in_out, args);

    keepGoing = false;
    for(int i=0 ; i<nb_coef_no_Q ; ++i){
      if(continue_crit(X[i], GX[i], diffMax)){
        keepGoing = true;
        break;
      }
    }
    
    if(iter % iter_grandAcc == 0){
      ++grand_acc;
      if(grand_acc == 1){
        std::memcpy(p_Y, p_GX, nb_coef_T * sizeof(double));
      } else if(grand_acc == 2){
        std::memcpy(p_GY, p_GX, nb_coef_T * sizeof(double));
      } else {
        std::memcpy(p_GGY, p_GX, nb_coef_T * sizeof(double));
        numconv = dm_update_X_IronsTuck(nb_coef_no_Q, Y, GY, GGY, delta_GX, delta2_X);
        if(numconv) break;
        compute_fe(Q, p_Y, p_GX, p_sum_other_means, p_sum_in_out, args);
        grand_acc = 0;
      }
    }

    // Other stopping criterion: change to SSR very small
    if(iter % 40 == 0){
      double ssr_old = ssr;

      // we compute the new SSR without allocating a temporary n_obs vector
      if(two_fe_algo){
        ssr = FE_info.compute_2_fe_ssr(p_GX, p_sum_other_means, p_sum_in_out, input);
      } else {
        ssr = FE_info.compute_ssr(p_GX, input);
      }

      if(stopping_crit(ssr_old, ssr, diffMax)){
        break;
      }
      
    }

  }

  //
  // Updating the output
  //

  double *p_beta_final = nullptr;
  if(two_fe_algo){
    // we end with a last iteration
    double *p_alpha_final = p_GX;
    p_beta_final = p_sum_other_means;

    FE_info.compute_fe_coef_2(p_alpha_final, p_alpha_final, p_beta_final, p_sum_in_out);
    FE_info.add_2_fe_coef_to_mu(p_alpha_final, p_beta_final, p_sum_in_out, output, false);
  } else {
    FE_info.add_all_fe_coef_to_mu(p_GX, output);
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
    
    if(two_fe_algo){
      // we add the other coefficients
      int n_coefs_FE1 = nb_coef_T;
      int n_coefs_FE2 = size_other_means;
      for(int m=0 ; m<n_coefs_FE2; ++m){
        fixef_values[n_coefs_FE1 + m] += p_beta_final[m];
      }
    }
  }

  bool conv = iter == iterMax ? false : true;

  return(conv);
}

void demean_single_gnl(int v, PARAM_DEMEAN* args){
  // v: integer identifying the variable to demean

  // Algorithm to quickly get the means
  // Q >= 3 => acceleration for 15 iter
  // if no convergence: conv 2 FEs
  // then acceleration again

  // data
  int iterMax = args->iterMax;
  
  // Note that the historical default was 15 iterations (still the case on 2024-02-14)
  int iter_warmup = args->algo_iter_warmup;
  int Q = args->Q;

  if(Q == 2){
    demean_acc_gnl(v, iterMax, args);
  } else {
    bool conv = false;
    
    // iter_warmup <= 0 means no warmup and we start directly with 2FE algo
    if(iter_warmup > 0){
      conv = demean_acc_gnl(v, iter_warmup, args);
      iter_warmup = 0;
    }

    if(conv == false && iterMax > iter_warmup){
      // convergence for the first 2 FEs
      int iter_max_2FE = iterMax / 2 - iter_warmup;
      if(iter_max_2FE > 0){
        demean_acc_gnl(v, iter_max_2FE, args, true);
      }
      
      // re-acceleration
      int iter_previous = args->p_iterations_all[v];
      demean_acc_gnl(v, iterMax - iter_previous, args);
    }
  }

  int *jobdone = args->jobdone;
  jobdone[v] = 1;

}

// Loop over demean_single
// [[Rcpp::export]]
List cpp_demean(SEXP y, SEXP X_raw, SEXP r_weights, int iterMax, double diffMax, SEXP r_nb_id_Q,
                SEXP fe_id_list, SEXP table_id_I, SEXP slope_flag_Q, SEXP slope_vars_list,
                SEXP r_init, int nthreads, 
                int algo_extraProj, int algo_iter_warmup, int algo_iter_projAfterAcc, 
                int algo_iter_grandAcc, bool save_fixef = false){
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
  RealMat m_y(y);
  int n_vars_y = m_y.ncol();
  bool useY = n_vars_y > 0;
  bool is_y_list = n_vars_y > 1 || TYPEOF(y) == VECSXP;
  int n_obs = m_y.nrow();

  // info on X
  RealMat m_X(X_raw);
  int n_vars_X = m_X.ncol();
  if(useY == false){
    n_obs = m_X.nrow();
  }
  bool useX = n_vars_X > 0;
  
  if(n_obs == 0 || n_obs == 1){
    // The data set is of length 1!!!!
    n_obs = 1;
    
    m_y = RealMat(y, true);
    n_vars_y = m_y.ncol();
    useY = true;
    
    m_X = RealMat(X_raw, true);
    n_vars_X = m_X.ncol();
    useX = n_vars_X > 0;
    
  }

  int n_vars = n_vars_y + n_vars_X;

  // initialisation if needed (we never initialize when only one FE, except if asked explicitly)
  bool isInit = Rf_xlength(r_init) != 1 && Q > 1;
  double *init = REAL(r_init);
  bool saveInit = ((isInit || init[0] != 0) && Q > 1) || init[0] == 666;

  // Creating the object containing all information on the FEs
  FEClass FE_info(n_obs, Q, r_weights, fe_id_list, r_nb_id_Q, table_id_I, slope_flag_Q, 
                  slope_vars_list);
  int nb_coef_T = FE_info.nb_coef_T;

  // output vector: (Note that if the means are provided, we use that vector and will modify it in place)
  int64_t n_total = static_cast<int64_t>(n_obs) * n_vars;
  SEXP output_values = PROTECT(Rf_allocVector(REALSXP, isInit ? 1 : n_total));

  double *p_output_origin = isInit ? init : REAL(output_values);

  //
  // vector of pointers: input/output
  //

  vector<double*> p_output(n_vars);
  p_output[0] = p_output_origin;
  for(int v=1 ; v<n_vars ; v++){
    p_output[v] = p_output[v - 1] + n_obs;
  }

  if(isInit == false){
    #pragma omp parallel for num_threads(nthreads) schedule(static, 1) if(nthreads > 1 && n_total >= 20000)
    for(int v=0 ; v<n_vars ; ++v){
      std::fill_n(p_output[v], n_obs, 0.0);
    }
  }

  vector<RealVec> p_input(n_vars);

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
    stop("save_fixef can only be used when there is no Xs.");
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
  args.algo_extraProj = algo_extraProj;
  args.algo_iter_warmup = algo_iter_warmup;
  args.algo_iter_projAfterAcc = algo_iter_projAfterAcc;
  // negative number or 0 means never
  args.algo_iter_grandAcc = algo_iter_grandAcc <= 0 ? 1000000 : algo_iter_grandAcc;

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


  //
  // the main loop
  //

  if(nthreads == 1){
    for(int v = 0 ; v < n_vars ; ++v){
      if(!*(args.stopnow)){
        if(Q == 1){
          demean_single_1(v, &args);
        } else {
          demean_single_gnl(v, &args);
        }
      }
    }
  } else {
    #pragma omp parallel num_threads(nthreads) default(shared)
    {
      #pragma omp master
      {
        for(int v = 1 ; v < n_vars ; ++v){
          #pragma omp task firstprivate(v) default(shared)
          {
            if(!*(args.stopnow)){
              if(Q == 1){
                demean_single_1(v, &args);
              } else {
                demean_single_gnl(v, &args);
              }
            }
          }
        }

        if(!*(args.stopnow)){
          if(Q == 1){
            demean_single_1(0, &args);
          } else {
            demean_single_gnl(0, &args);
          }
        }
      }
    }
  }


  if(*(args.stopnow)){
    UNPROTECT(1);
    stop("cpp_demean: User interrupt.");
  }

  //
  // save
  //

  List res; // a vector and a matrix

  int nrow = useX ? n_obs : 1;
  int ncol = useX ? n_vars_X : 1;
  NumericMatrix X_demean(nrow, ncol);

  double *p_X_demean = REAL(X_demean);
  #pragma omp parallel for num_threads(nthreads) if(nthreads > 1 && n_vars_X > 1)
  for(int k=0 ; k<n_vars_X ; ++k){
    RealVec p_input_tmp = p_input[k];
    double *p_output_tmp = p_output[k];
    double *p_X_col = p_X_demean + static_cast<size_t>(k) * n_obs;

    for(int i=0 ; i < n_obs ; ++i){
      p_X_col[i] = p_input_tmp[i] - p_output_tmp[i];
    }
  }

  res["X_demean"] = X_demean;


  if(is_y_list && useY){
    List y_demean(n_vars_y);
    vector<double*> p_y_demean_vec(n_vars_y);

    for(int v=0 ; v<n_vars_y ; ++v){
      NumericVector y_demean_tmp(n_obs);
      p_y_demean_vec[v] = REAL(y_demean_tmp);
      y_demean[v] = y_demean_tmp;
    }

    #pragma omp parallel for num_threads(nthreads) if(nthreads > 1 && n_vars_y > 1)
    for(int v=0 ; v<n_vars_y ; ++v){
      RealVec p_input_tmp  = p_input[n_vars_X + v];
      double *p_output_tmp = p_output[n_vars_X + v];
      double *p_y_col = p_y_demean_vec[v];

      for(int i=0 ; i < n_obs ; ++i){
        p_y_col[i] = p_input_tmp[i] - p_output_tmp[i];
      }
    }

    res["y_demean"] = y_demean;

  } else {
    NumericVector y_demean(useY ? n_obs : 1);
    if(useY){
      // y is always the last variable
      RealVec p_input_tmp = p_input[n_vars - 1];
      double *p_output_tmp = p_output[n_vars - 1];
      double *p_y_col = REAL(y_demean);
      #pragma omp parallel for num_threads(nthreads) if(nthreads > 1 && n_obs >= 65536)
      for(int i=0 ; i < n_obs ; ++i){
        p_y_col[i] = p_input_tmp[i] - p_output_tmp[i];
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




