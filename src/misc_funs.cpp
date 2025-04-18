/***********************************************************************************
 * ________________                                                                *
 * || Misc. funs ||                                                                *
 * ----------------                                                                *
 *                                                                                 *
 * Author: Laurent R. Berge                                                        *
 *                                                                                 *
 * Three groups of non-parallel functions:                                         *
 *   1) simple functions that are faster than base R (because more specific)       *
 *   2) function to obtain the FE coefficients after the estimation is done        *
 *   3) functions to lag variables                                                 *
 *                                                                                 *
 * Nothing much to say, everything is pretty explicit.                             *
 *                                                                                 *
 **********************************************************************************/


#include "fixest_main.h"
#include <math.h>
#include <stdio.h>
#include <cmath>
#include <functional>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector cpp_partialDerivative_other(int iterMax, int Q, int N, double epsDeriv, 
                                          NumericVector ll_d2,	NumericVector dx_dother, 
                                          NumericVector init, IntegerMatrix dumMat, 
                                          IntegerVector nbCluster){
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



//
// Getting the cluster coefficients
//

// [[Rcpp::export]]
List cpp_get_fe_gnl(int Q, int N, NumericVector sumFE, IntegerMatrix dumMat, 
                    IntegerVector cluster_sizes, IntegerVector obsCluster){
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
          for(int j=start_cluster[index_select] ; j<end_cluster[index_select] ; j++){
            obs = obsCluster(j, q);
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
double cpp_ssr_null(NumericVector y, NumericVector w = NumericVector(0)){
  // simple fun to compute the ssr of the null ols model
  // 2/3 times faster than pure r

  bool is_weight = w.length() > 1;

  int n = y.length();
  double denom = 0;

  // the "mean"
  double y_mean = 0;
  for(int i=0 ; i<n ; ++i){
    if(is_weight){
      y_mean += y[i] * w[i];
      denom += w[i];
    } else {
      y_mean += y[i];
    }
  }

  if(is_weight){
    y_mean = y_mean/denom;
  } else {
    y_mean = y_mean/n;
  }

  double res = 0, value = 0;
  for(int i=0 ; i<n ; ++i){
    value = y[i] - y_mean;
    if(is_weight){
      res += value * value * w[i];
    } else {
      res += value * value;
    }
  }

  return(res);
}

// [[Rcpp::export]]
double cpp_ssq(NumericVector x, NumericVector w = NumericVector(0)){
  // simple fun to compute the sum of the square of the elt of a vector
  // 30% faster than pure r (twice faster with weights)

  bool is_weight = w.length() > 1;

  int n = x.length();

  // the mean
  double res = 0;
  for(int i=0 ; i<n ; i++){
    if(is_weight){
      res += x[i] * x[i] * w[i];
    } else {
      res += x[i] * x[i];
    }
  }

  return res;
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
  // x is a vector

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
    /**************************************************************************
    * NOTA: I could have tweaked the previous if() to get rid of the condition
    *       but the code would have lost in clarity.
    *       For the lead: opposite to what is done before
    ***************************************************************************/
    int nlead = -nlag;
    i = nobs - 1;
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


// [[Rcpp::export]]
IntegerVector cpp_check_nested(SEXP fe_list, SEXP cluster_list, IntegerVector fe_sizes, int n){
  // Returns boolean vector of whether each FE is nested in the clusters

  int Q = Rf_length(fe_list);
  int G = Rf_length(cluster_list);

  // SEXP x0 = VECTOR_ELT(x, 0);

  IntegerVector res(Q);

  for(int q=0 ; q<Q ; ++q){

    int *pfe = INTEGER(VECTOR_ELT(fe_list, q));

    for(int g=0 ; g<G ; ++g){
      vector<int> fe_clust(fe_sizes[q]);

      int *pclust =INTEGER(VECTOR_ELT(cluster_list, g));

      bool nested = true;
      int fe_value = 0;
      int clust_value = 0;
      for(int i=0 ; i<n ; ++i){
        fe_value = pfe[i] - 1;
        clust_value = fe_clust[fe_value];
        if(clust_value == 0){
          fe_clust[fe_value] = pclust[i];
        } else if(clust_value != pclust[i]){
          nested = false;
          break;
        }
      }

      if(nested){
        res[q] = 1;
        break;
      }
    }
  }

  return res;
}


// [[Rcpp::export]]
NumericVector cpp_diag_XUtX(NumericMatrix X, NumericMatrix U){
  // computes the diagonal of X %*% U %*% t(X)

  int n = X.nrow();
  int K = X.ncol();

  NumericVector res(n);

  for(int i=0 ; i<n ; ++i){

    double res_i = 0;
    for(int k=0 ; k<K ; ++k){

      double xk = 0;
      for(int k2=0 ; k2<K ; ++k2){
        xk += X(i, k2) * U(k, k2);
      }

      res_i += xk * X(i,k);
    }

    res[i] = res_i;
  }

  return res;
}


class simple_vec_double{
  simple_vec_double() = delete;
  double *px_double = nullptr;
  int *px_int = nullptr;
  int n;
  bool is_real;
public:
  simple_vec_double(SEXP x);
  double operator[](int);
};

simple_vec_double::simple_vec_double(SEXP x){

  n = Rf_length(x);

  if(TYPEOF(x) == REALSXP){
    px_double = REAL(x);
    is_real = true;

  } else if(TYPEOF(x) == INTSXP){
    px_int = INTEGER(x);
    is_real = false;

  } else {
    stop("Error: Wrong argument type in cpp_factor_matrix.");
  }
}

double simple_vec_double::operator[](int i){
  if(i >= n){
    return 1;
  } else if(is_real){
    return px_double[i];
  } else {
    return static_cast<double>(px_int[i]);
  }
}


// [[Rcpp::export]]
NumericVector cpp_factor_matrix(IntegerVector fact, LogicalVector is_na_all, 
                                IntegerVector who_is_dropped, SEXP var, 
                                CharacterVector col_names){
  // fact: integer vector from 1 (!) to K, can contain NAs
  // Checking Na is cheap as opposed to populating the matrix, but having an argument avoids creating a new object

  int n = fact.length();
  int K = 0;

  // Finding out the TOTAL nber of cols (before removal)
  for(int i=0 ; i<n ; ++i){
    if(!is_na_all[i] && K < fact[i]){
      K = fact[i];
    }
  }

  // Are there values to remove? If so, we create a mapping
  int n_drop = who_is_dropped[0] == -1 ? 0 : Rf_length(who_is_dropped);
  bool IS_REMOVAL = n_drop > 0;
  std::vector<int> mapping;

  if(IS_REMOVAL){
    mapping.resize(K);

    for(int k=0 ; k<K ; ++k){
      mapping[k] = k;
    }

    // In mapping: values equal to -1 mean dropped value

    int k_drop = 0;
    for(int k=0 ; k<K ; ++k){
      if(k_drop < n_drop && k + 1 == who_is_dropped[k_drop]){
        ++k_drop;
        mapping[k] = -1;
      } else {
        mapping[k] -= k_drop;
      }
    }

    K -= k_drop;
  }

  NumericMatrix res(n, K);

  // The interacted var
  simple_vec_double my_var(var);

  // Filling the matrix
  for(int i=0 ; i<n ; ++i){
    if(is_na_all[i]){
      // we fill the row
      for(int k=0 ; k<K ; ++k){
        res(i, k) = NA_REAL;
      }

    } else if(IS_REMOVAL) {
      if(mapping[fact[i] - 1] != -1){
        res(i, mapping[fact[i] - 1]) = my_var[i];
      }

    } else {
      res(i, fact[i] - 1) = my_var[i];
    }
  }

  colnames(res) = col_names;

  return res;
}

// [[Rcpp::export]]
std::string cpp_add_commas(double x, int r = 1, bool whole = true){
  // a bit like (but not exactly equal to) format(x, nsmall = 1, big.mark = ",") but about 40-100 times faster
  // for whole numbers => no trailing digits
  // does not accept vectors, although super easy to expand to vectors

  double xr = round(x * pow(10, r)) / pow(10, r);

  std::string x_str = std::to_string(static_cast<int>(abs(xr)));
  std::string res;

  if(xr < 0){
    res.push_back('-');
    xr = -xr;
  }

  if(xr < 1000){
    res.insert(res.size(), x_str);

  } else {
    int n = x_str.size();
    int e = n; // e: exponent

    while(e > 0){
      res.push_back(x_str[n - e]);
      --e;
      if(e > 1 && e % 3 == 0) {
        res.push_back(',');
      }
    }
  }

  double rest = xr - floor(xr);
  bool is_whole = whole && (x - floor(x)) == 0;
  if(!is_whole && r > 0){
    // not a whole number

    res.push_back('.');

    std::string rest_str = std::to_string(rest);
    int nr = rest_str.size();

    for(int i=2 ; i < nr && i-1 <= r ; ++i){
      res.push_back(rest_str[i]);
    }
  }

  return res;
}


// [[Rcpp::export]]
List cpp_find_never_always_treated(IntegerVector cohort, NumericVector period){
  // Note that both cohort and period are sorted according to cohort

  IntegerVector always_treated;
  // cohort_ref => the guys that will be out of the interaction
  IntegerVector cohort_ref;

  int n = cohort.size();

  bool is_pos = false;
  bool is_neg = false;
  bool is_ok = false;

  int j = cohort[0];
  if(period[0] < 0){
    is_neg = true;
  } else {
    is_pos = true;
  }

  for(int i=1 ; i<n ; ++i){

    if(j == cohort[i]){
      if(!is_ok){
        // condition avoids checking period
        // only worth for very (very) long vectors

        if(period[i] < 0){
          is_neg = true;
          is_ok = is_pos;
        } else {
          is_pos = true;
          is_ok = is_neg;
        }
      }
    } else {
      // we change IDs
      if(!is_ok){
        if(is_pos) always_treated.push_back(j);
        cohort_ref.push_back(j);
      }

      // re-init
      j = cohort[i];
      is_neg = false;
      is_pos = false;
      is_ok = false;

    }
  }

  // Last element
  if(!is_ok){
    if(is_pos) always_treated.push_back(j);
    cohort_ref.push_back(j);
  }

  List res;
  res["always_treated"] = always_treated;
  res["ref"] = cohort_ref;

  return res;
}


// [[Rcpp::export]]
IntegerVector cpp_get_first_item(IntegerVector x, int n_items){
  // observation id of the first occurrence
  // x ranges from 1 to n_items
  // we return indexes R style

  int n = x.size();
  IntegerVector res(n_items);

  for(int i=0 ; i<n ; ++i){
    if(res[x[i] - 1] == 0){
      res[x[i] - 1] = i + 1;
    }
  }

  return res;
}


// [[Rcpp::export]]
IntegerVector cpp_combine_clusters(SEXP cluster_list, IntegerVector index){
  // cluster: list of integer vectors, each ranging from 1 to the number of cases
  // index: result of order() on the clusters

  if(TYPEOF(cluster_list) != VECSXP){
    stop("Internal error: Only lists are accepted!");
  }

  int Q = Rf_length(cluster_list);

  int n = index.size();
  IntegerVector res(n);

  // Loading the data
  vector<int*> pcluster(Q);
  for(int q=0 ; q<Q ; ++q){
    SEXP cluster_q = VECTOR_ELT(cluster_list, q);

    pcluster[q] = INTEGER(cluster_q);
  }

  // the observation ID
  int obs = index[0] - 1;

  // vector holding the current value
  vector<int> current_value(Q);

  // initialization
  int counter = 1;
  res[obs] = counter;
  for(int q=0 ; q<Q ; ++q){
    current_value[q] = pcluster[q][obs];
  }

  // we loop on the vector and flag values that are different
  int q = 0;
  for(int i=1 ; i<n ; ++i){
    obs = index[i] - 1;

    for(q=0 ; q<Q ; ++q){
      if(pcluster[q][obs] != current_value[q]){
        break;
      }
    }

    // if the condition holds => means the values are different
    if(q < Q){
      ++counter;
      // we save the new values
      for(; q<Q ; ++q){
        current_value[q] = pcluster[q][obs];
      }
    }

    res[obs] = counter;
  }

  return res;
}




// [[Rcpp::export]]
List cpp_cut(NumericVector x_sorted, NumericVector cut_points, IntegerVector is_included){
  // x_sorted: no NA, sorted
  // cut_points: bounds
  // is_included: for each bound, if it is included or not
  //

  int N = x_sorted.length();
  int n_cuts = cut_points.length();

  bool is_int = true;
  for(int i=0 ; i<N ; ++i){
    if(fabs(x_sorted[i] - round(x_sorted[i])) > 0.00000000001){
      is_int = false;
      break;
    }
  }

  IntegerVector x_int(N, n_cuts + 1);
  IntegerVector isnt_empty(n_cuts + 1);
  NumericVector value_min(n_cuts + 1);
  NumericVector value_max(n_cuts + 1);

  int index = 0;
  bool first = true;
  bool include = is_included[0];
  double cutoff = cut_points[0];
  int i = 0;
  while(i < N){

    if(include ? x_sorted[i] <= cutoff : x_sorted[i] < cutoff){

      if(first){
        isnt_empty[index] = true;
        value_min[index] = x_sorted[i];
        first = false;
      }

      x_int[i] = index + 1;
      ++i;

    } else {
      // we increment index
      // bins can be empty

      if(isnt_empty[index] && i > 0){
        value_max[index] = x_sorted[i - 1];
      }

      ++index;

      if(index == n_cuts){
        // last bin: we've done it at initialization
        isnt_empty[index] = true;
        value_min[index] = x_sorted[i];
        value_max[index] = x_sorted[N - 1];
        break;
      }

      include = is_included[index];
      cutoff = cut_points[index];
      first = true;
    }
  }

  if(index != n_cuts){
    value_max[index] = x_sorted[N - 1];
  }

  List res;
  res["x_int"] = x_int;
  res["isnt_empty"] = isnt_empty;
  res["value_min"] = value_min;
  res["value_max"] = value_max;
  res["is_int"] = is_int;

  return res;
}

// [[Rcpp::export]]
bool cpp_is_int(SEXP x){

  if(TYPEOF(x) == INTSXP){
    return true;
  }

  if(TYPEOF(x) != REALSXP){
    return false;
  }

  int N = Rf_length(x);
  double *px = REAL(x);

  bool is_int = true;
  for(int i=0 ; i<N ; ++i){
    if(fabs(px[i] - round(px[i])) > 0.00000000001){
      is_int = false;
      break;
    }
  }

  return is_int;
}

// [[Rcpp::export]]
double cpp_hash_string(std::string x){
  // simple function hashing a string
  // used to identify tables in etable

  double res = std::hash<std::string>{}(x);

  return res;
}

inline bool is_md_markup(const char * x, int i, int n){
  if(x[i] != '*') return false;

  if(i + 1 >= n || x[i + 1] != '*') return true;
  if(i + 2 >= n ||  x[i + 2] != '*') return true;
  if(i + 3 < n && x[i + 3] == '*') return false;
  return true;
}

inline bool is_special_char(const char x){
  return x == '&' || x == '%' || x == '_' || x == '^' || x == '#';
}

inline int find_id_markup(const char * x, int i, int n){
  if(i + 1 >= n || x[i + 1] != '*') return 1;
  if(i + 2 >= n ||  x[i + 2] != '*') return 2;
  return 3;
}

std::string apply_escape_markup(const char * x){
  //
  //
  //
  // element 0 is the main character vector: it is never closed

  int n = std::strlen(x);
  std::vector<std::string> tmp_all(4, "");
  std::string tmp = "", res = "";

  std::vector<bool> is_open(4, false);
  std::vector<int> id_open;
  id_open.push_back(0);

  std::vector<std::string> open_tags{"", "\\textit{", "\\textbf{", "\\textbf{\\textit{"};
  std::vector<std::string> closing_tags{"", "}", "}", "}}"};
  std::vector<std::string> stars{"", "*", "**", "***"};

  int id_mkp = 0;

  int i = 0, i_save = 0;
  while(i < n){
    if(x[i] == '$'){
      if(i > 0 && x[i - 1] == '\\'){
        tmp_all[id_mkp] += '$';
      } else {
        // This is an equation
        i_save = i;
        ++i;
        tmp = "$";
        while(i < n && (x[i] != '$' || x[i - 1] == '\\')){
          tmp += x[i];
          ++i;
        }

        if(i == n){
          // we went all the way without a closing $
          // => we come back, escape it and continue
          tmp_all[id_mkp] += "\\$";
          i = i_save;

        } else {
          tmp_all[id_mkp] += tmp + "$";
        }
      }
    } else if(is_special_char(x[i])){
      if(i > 0 && x[i - 1] == '\\'){
        // we do nothing
        tmp_all[id_mkp] += x[i];
      } else {
        // we escape
        tmp_all[id_mkp] += "\\";
        tmp_all[id_mkp] += x[i];
      }

    } else if(x[i] == '\\'){
      if(i + 1 < n && x[i + 1] == '*'){
        // escape of MD markup
        ++i;
        while(i < n && x[i] == '*'){
          tmp_all[id_mkp] += '*';
          ++i;
        }
        --i; // compensated at the end of the main while
      } else {
        tmp_all[id_mkp] += '\\';
      }
    } else if(is_md_markup(x, i, n)){
      id_mkp = find_id_markup(x, i, n);
      i += id_mkp - 1;

      if(is_open[id_mkp]){

        for(int j=id_open.size() - 1 ; j >= 1 ; --j){

          int id_j = id_open[j];
          int id_prev = id_open[j - 1];

          if(id_j == id_mkp){
            // we apply the markup
            tmp = open_tags[id_mkp] + tmp_all[id_mkp] + closing_tags[id_mkp];

            tmp_all[id_prev] += tmp;

            is_open[id_mkp] = false;
            id_open.pop_back();
            tmp_all[id_mkp] = "";
            id_mkp = id_prev;
            break;

          } else {
            // we flush: the markup is invalid
            // example:  "** bonjour * les gens **"

            tmp = stars[id_j] + tmp_all[id_j];
            tmp_all[id_j] = "";
            tmp_all[id_prev] += tmp;

            is_open[id_j] = false;
            id_open.pop_back();
          }
        }

      } else {
        is_open[id_mkp] = true;
        id_open.push_back(id_mkp);
      }

    } else {
      tmp_all[id_mkp] += x[i];
    }
    ++i;
  }

  // we flush the rest if needed
  for(int j=id_open.size() - 1 ; j >= 1 ; --j){
    int id_j = id_open[j];
    int id_prev = id_open[j - 1];

    tmp = stars[id_j] + tmp_all[id_j];
    tmp_all[id_j] = "";
    tmp_all[id_prev] += tmp;
  }

  return tmp_all[0];
}


// [[Rcpp::export]]
StringVector cpp_escape_markup(SEXP Rstr){
  // cpp_escape_markup("**bonjour** *les* ***gens * \\***heureux*** ")
  // cpp_escape_markup("stars: 10%: *, 5%: **, 1%: ***")

  int n = LENGTH(Rstr);
  StringVector res(n);

  if(n == 0){
    return res;
  }

  for(int i=0 ; i<n ; ++i){
    res[i] = apply_escape_markup(CHAR(STRING_ELT(Rstr, i)));
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

  RealMat mat(x);

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

  std::vector<int> bounds = set_parallel_scheme(nobs, nthreads);

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



































