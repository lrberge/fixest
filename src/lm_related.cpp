/**************************************************************
* ___________________                                         *
* || OLS Functions ||                                         *
* -------------------                                         *
*                                                             *
* Set of functions to perform OLS estimations.                *
*                                                             *
* In general, the functions are slower than BLAS ones         *
* but they have the merit of doing exacly what I want.        *
*                                                             *
**************************************************************/

#include <Rcpp.h>
#include <cmath>
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#endif
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
NumericMatrix cpp_invert_tri(const NumericMatrix &R_mat, int nthreads = 1){
    // Parallel is not super efficient but it's better than nothing: 2 threads: 40%; 4 threads: 60%; speedup
    //

    // Startegy: we invert by bands (b)
    int K = R_mat.ncol();
    NumericMatrix res(K, K);
    NumericMatrix R_prime(K, K);

    // initialization of R prime
    for(int i=0 ; i<K ; ++i){
        for(int j=i+1 ; j<K ; ++j){
            R_prime(j, i) = R_mat(i, j);
        }
    }

    for(int i=0 ; i<K ; ++i){
        res(i, i) = 1/R_mat(i, i);
    }


    for(int b=1 ; b<K ; ++b){

        #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for(int i=0 ; i<K-b ; ++i){
            double numerator = 0;
            int col = i+b;
            for(int k=i+1 ; k<=i+b ; ++k){
                numerator -= R_prime(k, i) * res(k, col);
            }
            res(i, col) = numerator * res(i, i);
        }
    }

    return res;
}

// [[Rcpp::export]]
NumericMatrix cpp_tri_tprod(const NumericMatrix &R_mat, int nthreads = 1){

    int K = R_mat.ncol();

    NumericMatrix R_prime(K, K);
    // initialization of R prime
    for(int i=0 ; i<K ; ++i){
        for(int j=i ; j<K ; ++j){
            R_prime(j, i) = R_mat(i, j);
        }
    }

    NumericMatrix RRt(K, K);

    #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for(int i=0 ; i<K ; ++i){
        for(int j=i ; j<K ; ++j){

            double value = 0;

            for(int k=i ; k<K ; ++k){
                value += R_prime(k, j) * R_prime(k, i);
            }
            RRt(i, j) = value;
            RRt(j, i) = value;
        }
    }

    return RRt;
}

// [[Rcpp::export]]
List cpp_cholesky(NumericMatrix X, int nthreads = 1){
    // X est symetrique, semi definie positive
    // rank-revealing on-the-fly

    int K = X.ncol();

    NumericMatrix L(K, K);
    LogicalVector id_excl(K);
    int n_excl = 0;
    double tol = 1.0/100000.0/100000.0;

    for(int j=0 ; j<K ; ++j){

        // implicit pivoting (it's like 0-rank variables are stacked in the end)

        double L_jj = X(j, j);
        for(int k=0 ; k<j ; ++k){

            if(id_excl[k]) continue;

            L_jj -= L(j, k) * L(j, k);
        }

        if(L_jj < tol){
            n_excl++;
            id_excl[j] = true;
            // Rcout << "excluded: " << j << ", L_jj: " << L_jj << "\n";

            // Corner case, may happen:
            if(n_excl == K){
                List res;
                res["all_removed"] = true;
                return res;
            }

            continue;
        }

        L_jj = sqrt(L_jj);

        L(j, j) = L_jj;

        #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for(int i=j+1 ; i<K ; ++i){

            double value = X(i, j);
            for(int k=0 ; k<j ; ++k){
                if(id_excl[k]) continue;

                value -= L(i, k) * L(j, k);
            }
            L(i, j) = value/L_jj;
        }
    }

    // We recontruct a R matrix and compute (RtR)-1

    // we reconstruct the R matrix, otherwise it's just a mess to handle excluded variables on the fly
    int K_excl = K - n_excl;
    NumericMatrix R_mat(K_excl, K_excl);
    // column fill scheme
    int n_j_excl = 0;
    for(int j=0 ; j<K ; ++j){

        if(id_excl[j]){
            ++n_j_excl;
            continue;
        }

        int n_i_excl = n_j_excl;
        for(int i=j ; i<K ; ++i){
            if(id_excl[i]){
                ++n_i_excl;
                continue;
            }
            R_mat(j - n_j_excl, i - n_i_excl) = L(i, j);
        }
    }


    NumericMatrix R_inv = cpp_invert_tri(R_mat, nthreads);
    NumericMatrix XtX_inv = cpp_tri_tprod(R_inv, nthreads);

    List res;
    // res["R"] = R_mat;
    // res["L"] = L;
    res["XtX_inv"] = XtX_inv;
    res["id_excl"] = id_excl;


    return res;
}


static bool sparse_check(const NumericMatrix &X){
    // Super cheap sparsity test
    // works even for small K large N

    int N = X.nrow();
    int K = X.ncol();

    if(K < 5) return false;
    if(N < 1000 && K < 100) return false;
    if(N < 100) return false; // avoids corner cases where you have more vars than obs

    // Let's look at some rows
    // If there are factors, there should be 0s
    //
    // Let's put it differently: if there are no 0s, then it's not factors for sure

    int n0_first_row = 0;
    int n0_middle_row = 0;
    int n0_last_row = 0;

    for(int k=0 ; k<K ; ++k){
        if(X(0, k) == 0) ++n0_first_row;
        if(X(N/2, k) == 0) ++n0_middle_row;
        if(X(N, k) == 0) ++n0_last_row;
    }

    // Rcout << "1st: " << n0_first_row << ", mid: " << n0_middle_row << ", last: " << n0_last_row << "\n";

    if(n0_first_row > K/2 && n0_middle_row > K/2 && n0_last_row > K/2) return true;

    return false;
}

// [[Rcpp::export]]
List cpp_sparse_products(NumericMatrix X, NumericVector w, NumericVector y, bool correct_0w = false, int nthreads = 1){

    int N = X.nrow();
    int K = X.ncol();

    bool isWeight = false;
    if(w.length() > 1){
        isWeight = true;
    }

    NumericMatrix XtX(K, K);
    NumericVector Xty(K);

    if(sparse_check(X) == false){
        // NOT SPARSE

        if(isWeight){
            NumericMatrix wX(Rcpp::clone(X));
            for(int k=0 ; k<K ; ++k){
                for(int i=0 ; i<N ; ++i){
                    wX(i, k) *= sqrt(w[i]);
                }
            }

            NumericVector wy(Rcpp::clone(y));
            for(int i=0 ; i<N ; ++i){
                wy(i) *= sqrt(w[i]);
            }

            //
            // XtX
            //

            // specific scheme for large N and K == 1
            if(K == 1){
                std::vector<double> all_values(nthreads, 0);
                #pragma omp parallel num_threads(nthreads)
                {
                    int i = omp_get_thread_num()*N/omp_get_num_threads();
                    int stop = (omp_get_thread_num()+1)*N/omp_get_num_threads();
                    double val = 0;
                    for(; i<stop ; ++i){
                        val += wX(i, 0) * wX(i, 0);
                    }
                    all_values[omp_get_thread_num()] = val;
                }

                double value = 0;
                for(int m=0 ; m<nthreads ; ++m){
                    value += all_values[m];
                }

                XtX(0, 0) = value;

            } else {
                // We use this trick to even out the load on the threads
                int nValues = K * K;
                NumericVector all_values(nValues);

                #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        double val = 0;
                        for(int i=0 ; i<N ; ++i){
                            val += wX(i, k_row) * wX(i, k_col);
                        }

                        all_values[index] = val;
                    }

                }

                // save
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        XtX(k_row, k_col) = all_values[index];
                        XtX(k_col, k_row) = all_values[index];
                    }

                }
            }


            //
            // Xty
            //

            if(K == 1){
                std::vector<double> all_values(nthreads, 0);
                #pragma omp parallel num_threads(nthreads)
                {
                    int i = omp_get_thread_num()*N/omp_get_num_threads();
                    int stop = (omp_get_thread_num()+1)*N/omp_get_num_threads();
                    double val = 0;
                    for(; i<stop ; ++i){
                        val += wX(i, 0) * wy[i];
                    }
                    all_values[omp_get_thread_num()] = val;
                }

                double value = 0;
                for(int m=0 ; m<nthreads ; ++m){
                    value += all_values[m];
                }

                Xty[0] = value;
            } else  {
                #pragma omp parallel for num_threads(nthreads)
                for(int j=0 ; j<K ; ++j){
                    double val = 0;
                    for(int i=0 ; i<N ; ++i){
                        val += wX(i, j) * wy[i];
                    }

                    Xty[j] = val;
                }
            }

        } else {
            // Identique, mais pas besoin de faire une copie de X ni de y qui peuvent etre couteuses
            //
            // XtX
            //

            if(K == 1){
                std::vector<double> all_values(nthreads, 0);
                #pragma omp parallel num_threads(nthreads)
                {
                    int i = omp_get_thread_num()*N/omp_get_num_threads();
                    int stop = (omp_get_thread_num()+1)*N/omp_get_num_threads();
                    double val = 0;
                    for(; i<stop ; ++i){
                        val += X(i, 0) * X(i, 0);
                    }
                    all_values[omp_get_thread_num()] = val;
                }

                double value = 0;
                for(int m=0 ; m<nthreads ; ++m){
                    value += all_values[m];
                }

                XtX(0, 0) = value;

            } else {
                // We use this trick to even out the load on the threads
                int nValues = K * K;
                NumericVector all_values(nValues);

                #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        double val = 0;
                        for(int i=0 ; i<N ; ++i){
                            val += X(i, k_row) * X(i, k_col);
                        }

                        all_values[index] = val;
                    }

                }

                // save
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        XtX(k_row, k_col) = all_values[index];
                        XtX(k_col, k_row) = all_values[index];
                    }

                }
            }


            //
            // Xty
            //

            if(K == 1){
                std::vector<double> all_values(nthreads, 0);
                #pragma omp parallel num_threads(nthreads)
                {
                    int i = omp_get_thread_num()*N/omp_get_num_threads();
                    int stop = (omp_get_thread_num()+1)*N/omp_get_num_threads();
                    double val = 0;
                    for(; i<stop ; ++i){
                        val += X(i, 0) * y[i];
                    }
                    all_values[omp_get_thread_num()] = val;
                }

                double value = 0;
                for(int m=0 ; m<nthreads ; ++m){
                    value += all_values[m];
                }

                Xty[0] = value;

            } else  {
                #pragma omp parallel for num_threads(nthreads)
                for(int j=0 ; j<K ; ++j){
                    double val = 0;
                    for(int i=0 ; i<N ; ++i){
                        val += X(i, j) * y[i];
                    }

                    Xty[j] = val;
                }
            }
        }

        List res;
        res["XtX"] = XtX;
        res["Xty"] = Xty;

        return res;
    }

    //
    // SPARSE case
    //


    std::vector<int> n_j(K, 0);
    std::vector<int> all_i;
    std::vector<double> x;

    // CAN'T BE PARALLELIZED
    for(int k=0 ; k<K ; ++k){
        for(int i=0 ; i<N ; ++i){

            if(X(i, k) != 0){
                if(isWeight){
                    if(correct_0w && w[i] == 0){
                        continue;
                    }
                    ++n_j[k];
                    all_i.push_back(i);

                    // BEWARE: here there is no square root for w!
                    // (as opposed to the non-sparse case)
                    x.push_back(X(i, k) * w[i]);
                } else {
                    // No weights
                    ++n_j[k];
                    all_i.push_back(i);
                    x.push_back(X(i, k));
                }
            }
        }
    }


    //
    // The cross-product
    //


    // starting values of j's
    std::vector<int> start_j(K, 0);
    for(int k=1 ; k<K ; ++k){
        start_j[k] = start_j[k-1] + n_j[k-1];
    }
    start_j.push_back(start_j[K-1] + n_j[K-1]);



    #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for(int j1=0 ; j1<K ; ++j1){

        int k = 0, start = 0, end = 0;
        double tmp = 0;

        // Upper triangle
        for(int j2=j1 ; j2<K ; ++j2){

            if(n_j[j1] < n_j[j2]){
                start = start_j[j1];
                end = start_j[j1 + 1];
                k = j2;
            } else {
                start = start_j[j2];
                end = start_j[j2 + 1];
                k = j1;
            }

            double value = 0;
            for(int index=start ; index<end ; ++index){
                tmp = X(all_i[index], k);
                if(tmp != 0) value += tmp * x[index];
            }

            if(value == 0) continue;

            XtX(j1, j2) = value;
            XtX(j2, j1) = value;

        }
    }


    //
    // now Xty
    //

    #pragma omp parallel for num_threads(nthreads)
    for(int j=0 ; j<K ; ++j){

        int start = start_j[j];
        int end = start_j[j + 1];

        double value = 0;
        for(int index=start ; index<end ; ++index){
            value += y[all_i[index]] * x[index];
        }

        if(value == 0) continue;

        Xty[j] = value;
    }

    List res;
    res["XtX"] = XtX;
    res["Xty"] = Xty;

    return res;
}



// [[Rcpp::export]]
NumericMatrix cpppar_crossprod(NumericMatrix X, NumericVector w, int nthreads){

    int N = X.nrow();
    int K = X.ncol();

    bool isWeight = false;
    if(w.length() > 1){
        isWeight = true;
    }

    NumericMatrix XtX(K, K);

    if(sparse_check(X) == false){
        // NOT SPARSE

        if(isWeight){
            // I don't use the sqrt, because I use the function when weights are negative too (ll_d2 used as 'weight')
            NumericMatrix wX(Rcpp::clone(X));
            for(int k=0 ; k<K ; ++k){
                for(int i=0 ; i<N ; ++i){
                    wX(i, k) *= w[i];
                }
            }

            //
            // XtX
            //

            // specific scheme for large N and K == 1
            if(K == 1){
                std::vector<double> all_values(nthreads, 0);
                #pragma omp parallel num_threads(nthreads)
                {
                    int i = omp_get_thread_num()*N/omp_get_num_threads();
                    int stop = (omp_get_thread_num()+1)*N/omp_get_num_threads();
                    double val = 0;
                    for(; i<stop ; ++i){
                        val += X(i, 0) * wX(i, 0);
                    }
                    all_values[omp_get_thread_num()] = val;
                }

                double value = 0;
                for(int m=0 ; m<nthreads ; ++m){
                    value += all_values[m];
                }

                XtX(0, 0) = value;

            } else {
                // We use this trick to even out the load on the threads
                int nValues = K * K;
                NumericVector all_values(nValues);

                #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        double val = 0;
                        for(int i=0 ; i<N ; ++i){
                            val += X(i, k_row) * wX(i, k_col);
                        }

                        all_values[index] = val;
                    }

                }

                // save
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        XtX(k_row, k_col) = all_values[index];
                        XtX(k_col, k_row) = all_values[index];
                    }

                }
            }

        } else {
            // Identique, mais pas besoin de faire une copie de X ni de y qui peuvent etre couteuses
            //
            // XtX
            //

            if(K == 1){
                std::vector<double> all_values(nthreads, 0);
                #pragma omp parallel num_threads(nthreads)
                {
                    int i = omp_get_thread_num()*N/omp_get_num_threads();
                    int stop = (omp_get_thread_num()+1)*N/omp_get_num_threads();
                    double val = 0;
                    for(; i<stop ; ++i){
                        val += X(i, 0) * X(i, 0);
                    }

                    all_values[omp_get_thread_num()] = val;
                }

                double value = 0;
                for(int m=0 ; m<nthreads ; ++m){
                    value += all_values[m];
                }

                XtX(0, 0) = value;

            } else {
                // We use this trick to even out the load on the threads
                int nValues = K * K;
                NumericVector all_values(nValues);

                #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        double val = 0;
                        for(int i=0 ; i<N ; ++i){
                            val += X(i, k_row) * X(i, k_col);
                        }

                        all_values[index] = val;
                    }

                }

                // save
                for(int index=0 ; index<nValues ; index++){
                    int k_row = index % K;
                    int k_col = index / K;

                    if(k_row <= k_col){
                        XtX(k_row, k_col) = all_values[index];
                        XtX(k_col, k_row) = all_values[index];
                    }

                }
            }
        }

        return XtX;
    }

    //
    // SPARSE case
    //


    std::vector<int> n_j(K, 0);
    std::vector<int> all_i;
    std::vector<double> x;

    // CAN'T BE PARALLELIZED
    for(int k=0 ; k<K ; ++k){
        for(int i=0 ; i<N ; ++i){

            if(X(i, k) != 0){
                if(isWeight){
                    ++n_j[k];
                    all_i.push_back(i);

                    // BEWARE: here there is no square root for w!
                    // (as opposed to the non-sparse case)
                    x.push_back(X(i, k) * w[i]);
                } else {
                    // No weights
                    ++n_j[k];
                    all_i.push_back(i);
                    x.push_back(X(i, k));
                }
            }
        }
    }


    //
    // The cross-product
    //


    // starting values of j's
    std::vector<int> start_j(K, 0);
    for(int k=1 ; k<K ; ++k){
        start_j[k] = start_j[k-1] + n_j[k-1];
    }
    start_j.push_back(start_j[K-1] + n_j[K-1]);



    #pragma omp parallel for num_threads(nthreads) schedule(static, 1)
    for(int j1=0 ; j1<K ; ++j1){

        int k = 0, start = 0, end = 0;
        double tmp = 0;

        // Upper triangle
        for(int j2=j1 ; j2<K ; ++j2){

            if(n_j[j1] < n_j[j2]){
                start = start_j[j1];
                end = start_j[j1 + 1];
                k = j2;
            } else {
                start = start_j[j2];
                end = start_j[j2 + 1];
                k = j1;
            }

            double value = 0;
            for(int index=start ; index<end ; ++index){
                tmp = X(all_i[index], k);
                if(tmp != 0) value += tmp * x[index];
            }

            if(value == 0) continue;

            XtX(j1, j2) = value;
            XtX(j2, j1) = value;

        }
    }

    return XtX;
}




// [[Rcpp::export]]
NumericMatrix cpp_mat_reconstruct(NumericMatrix X, Rcpp::LogicalVector id_excl){

    int K = id_excl.length();
    int K_small = X.ncol();

    NumericMatrix res(K, K);

    int n_col_excl = 0;
    for(int j=0 ; j<K_small ; ++j){

        while(id_excl[j + n_col_excl]) ++n_col_excl;

        int col = j + n_col_excl;
        int n_row_excl = 0;
        for(int i=0 ; i<K_small ; ++i){
            while(id_excl[i + n_row_excl]) ++n_row_excl;
            res(i+n_row_excl, col) = X(i, j);
        }
    }

    return res;
}
