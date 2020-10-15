/**********************************************************************
 * _____________                                                      *
 * || QUFing ||                                                       *
 * -------------                                                      *
 *                                                                    *
 * Author: Laurent R. Berge                                           *
 *                                                                    *
 * QUF stands for "quick unclass factor", an operation consisting     *
 * in transforming a vector of arbitrary values into an integer       *
 * vector ranging from 1 to the number of unique values (i.e:         *
 *   a)  (50, 55, 32, 12, 12) => (3, 4, 2, 1, 1)                      *
 *   b)  ("a", "d", "a", "b") => (1, 3, 1, 2)                         *
 *  )                                                                 *
 *                                                                    *
 *  The code here is used to quf vectors of integers, floats          *
 *  or strings. I convert any other type of identifier to             *
 *  character if they're not numeric before getting into this         *
 *  function.                                                         *
 *                                                                    *
 *  All we care in qufing is to get a vector from 1 to n the number   *
 *  of unique values. We don't care about the order (e.g., in         *
 *  example a) the result (1,2,3,4,4) would have been fine).          *
 *                                                                    *
 *  Therefore, qufing integer vectors of small range is both          *
 *  very simple and very fast: we just have to set a lookup table.    *
 *  For integers of large range, setting up a large                   *
 *  lookup table can be too costly, therefore I sort the vector       *
 *  first and then create the new values.                             *
 *  The process is similar for doubles: first sort, then get the      *
 *  new values.                                                       *
 *  Finally for string vectors, I use R's feature of stocking         *
 *  character strings in a unique location. They are therefore        *
 *  uniquely identified by their pointer. I then tranform the         *
 *  pointer to int, and sort accordingly.                             *
 *                                                                    *
 *  The sorting, when sorting is required, is radix-based. My code    *
 *  has greatly benefited from two sources which clarified a lot      *
 *  of implicit things:                                               *
 *  - http://codercorner.com/RadixSortRevisited.htm (Pierre Terdiman) *
 *  - http://stereopsis.com/radix.html (Michael Herf)                 *
 *                                                                    *
 *  As said before, all that matters for qufing is that identical     *
 *  values are consecutive after the sort. Thus, the endianness of    *
 *  the system is of no importantce: whatever the order of bytes      *
 *  on which we sort, we obtain what we want.                         *
 *                                                                    *
 *                                                                    *
 *********************************************************************/

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using std::vector;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]




inline unsigned long long float_to_ull(void *u, int i) {
    unsigned long long *pu_ll = reinterpret_cast<unsigned long long*>(u);
    unsigned long long u_ull = pu_ll[i];
    unsigned long long mask = -(u_ull >> 63) | 0x8000000000000000;
    return (u_ull ^ mask);
}

inline double ull_to_float(unsigned long long u_ull) {
    unsigned long long mask = ((u_ull >> 63) - 1) | 0x8000000000000000;
    unsigned long long res_ull = u_ull ^ mask;
    unsigned long long *pres_ull = &res_ull;
    double *pres = reinterpret_cast<double*>(pres_ull);
    double res = *pres;
    return res;
}

void quf_double(vector<int> &x_uf, void *px, vector<double> &x_unik, bool is_string = false){

    // x_uf: x unclassed factor
    // px: pointer to x vector (R vector) -- READ ONLY!!!
    // x_unik: empty vector
    // px: either double or ULL (in case of strings)

    int n = x_uf.size();

    double *px_dble = (double *)px;
    unsigned long long *px_ull = (unsigned long long *)px;

    // variables
    unsigned long long x_ull_current = 0;
    // char: 2 octets: 0-255
    // one double is made of 8 char
    unsigned char x_char = 0;
    vector<unsigned long long> x_ulong(is_string ? 1 : n), x_tmp(n);

    // in case the data is string, we use px as the x_ulong
    unsigned long long *px_ulong = is_string ? px_ull : x_ulong.data();


    // radix array
    int radix_table[8][256] = { {0} };

    // 1) Counting
    for(int i=0 ; i<n ; ++i){

        if(is_string){
            x_ull_current = px_ull[i];
        } else {
            // we change x double into ulong after flipping to keep order
            x_ull_current = float_to_ull(px_dble, i);
            x_ulong[i] = x_ull_current;
        }

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

    // 2) Cumulating
    for(int b = 0 ; b < 8 ; ++b){
        for(int d = 1 ; d < 256 ; ++d){
            radix_table[b][d] += radix_table[b][d - 1];
        }
    }

    // 3) Sorting
    unsigned long long *x_read = px_ulong;
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


    if(o_read != x_order.data()){
        // we copy the result
        memcpy(x_order.data(), o_read, sizeof(int) * n);
    }

    // We unclass, starting at 1
    k=1;
    x_unclass[0] = k;

    unsigned long long xi, xim1 = x_read[0];
    x_unik.push_back(is_string ? static_cast<double>(x_order[0]  + 1) : ull_to_float(xim1));

    for(int i=1 ; i<n ; ++i){
        xi = x_read[i];
        if(xi != xim1){
            ++k;
            x_unik.push_back(is_string ? static_cast<double>(x_order[i]  + 1) : ull_to_float(xi));
        }
        x_unclass[i] = k;
        xim1 = xi;
    }

    // We put into the right order
    for(int i=0 ; i<n ; ++i){
        x_uf[x_order[i]] = x_unclass[i];
    }

}

void quf_int_gnl(vector<int> &x_uf, void *px, vector<double> &x_unik, int x_min, bool is_double){
    // we can sort a range up to 2**31 => ie not the full int range
    // for ranges > 2**31 => as double
    // px: pointer to the values of x -- R vector READ ONLY!!!
    // px: either a R vector of type integer // either a R vector of type double
    // px: we know it thanks to the value of 'is_double'

    int *px_int = (int *)px;
    double *px_dble = (double *)px;

    int n = x_uf.size();

    // variables
    int x_uint_current = 0;
    vector<int> x_uint(n), x_tmp(n);

    // radix array
    int radix_table[4][256] = { {0} };

    // 1) Counting
    for(int i=0 ; i<n ; ++i){
        // x_uint_current = px[i] - x_min;
        x_uint_current = is_double ? static_cast<int>(px_dble[i] - x_min) : px_int[i] - x_min;
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

    if(o_read != x_order.data()){
        // we copy the result, if needed
        memcpy(x_order.data(), o_read, sizeof(int) * n);
    }

    // We unclass, starting at 1
    k=1;
    x_unclass[0] = k;

    double xi, xim1 = x_read[0];
    x_unik.push_back(xim1 + x_min);

    for(int i=1 ; i<n ; ++i){
        xi = x_read[i];
        if(xi != xim1){
            ++k;
            x_unik.push_back(static_cast<double>(xi + x_min));
        }
        x_unclass[i] = k;
        xim1 = xi;
    }

    // We put into the right order
    for(int i=0 ; i<n ; ++i){
        x_uf[x_order[i]] = x_unclass[i];
    }


}

void quf_int(vector<int> &x_uf, void *px, vector<double> &x_unik, int x_min, int max_value, bool is_double = false){
    // Principle:
    // we go through x only once
    // we keep a table of all x values
    // whenever we find a new value of x, we save it
    // px: pointer to the values of x -- R vector READ ONLY!!!
    // px: either a R vector of type integer // either a R vector of type double
    // px: we know it thanks to the value of 'is_double'

    int *px_int = (int *)px;
    double *px_dble = (double *)px;

    int n = x_uf.size();
    int n_unik = 0; // nber of uniques minus one

    // radix array
    vector<int> x_lookup(max_value + 1, 0);

    int x_tmp, x_pos;
    for(int i=0 ; i<n ; ++i){
        x_tmp = is_double ? static_cast<int>(px_dble[i]) - x_min : px_int[i] - x_min;
        // x_tmp = x[i] - x_min;

        if(x_lookup[x_tmp] == 0){
            ++n_unik;
            x_uf[i] = n_unik;
            // x_unik.push_back(static_cast<double>(px_int[i]));
            x_unik.push_back(is_double ? px_dble[i] : static_cast<double>(px_int[i]));
            x_lookup[x_tmp] = n_unik;
        } else {
            x_pos = x_lookup[x_tmp];
            x_uf[i] = x_pos;
        }
    }

}

// [[Rcpp::export]]
List cpp_quf_str(SEXP x){

    int n = Rf_length(x);

    vector<int> x_uf(n);
    vector<double> x_unik;

    // we change the string vector into a uint64 vector
    // since we don't care about the sorting order, we can just
    // recast into double (an not uint64)

    vector<double> x_ull(n);
    std::uintptr_t xi_uintptr;
    unsigned long long xi_ull;


    for(int i=0 ; i<n ; ++i){
        const char *pxi = CHAR(STRING_ELT(x, i));
        xi_uintptr = reinterpret_cast<std::uintptr_t>(pxi);
        xi_ull = static_cast<unsigned long long>(xi_uintptr);

        x_ull[i] = xi_ull;

        // Rcout << xi_uintptr << "  ----  " << xi_ull << "\n";
    }

    quf_double(x_uf, x_ull.data(), x_unik, true);

    List res;
    res["x_uf"] = x_uf;
    res["x_unik"] = x_unik;

    return res;
}




// [[Rcpp::export]]
List cpp_quf_gnl(SEXP x){

    // INT: we try as possible to send the data to quf_int, the most efficient function
    // for data of large range, we have a separate algorithms that avoids the creation
    // of a possibly too large lookup table.
    // for extreme ranges (> 2**31) we use the algorithm for double

    // DOUBLE: quf_int is the most efficient function.
    // even if we have a vector of double, we check whether the underlying structure
    // is in fact int, so that we can send it to quf_int
    // if the vector is really double, then the checking cost is negligible
    // if the vector is in fact int, the cost of checking is largely compensated by the
    //   efficiency of the algorithm

    // STRING: for string vectors, we first transform them into ULL using their pointer
    // before applying them the algorithm for doubles
    // note that the x_unik returned is NOT of type string but is equal to
    // the location of the unique elements

    int n = Rf_length(x);

    vector<int> x_uf(n);
    vector<double> x_unik;

    // preparation for strings
    bool IS_STR = false;
    vector<unsigned long long> x_ull;

    bool IS_INT = false;
    bool is_int_in_double = false;
    if(TYPEOF(x) == REALSXP){
        // we check if underlying structure is int
        IS_INT = true;
        double *px = REAL(x);
        for(int i=0 ; i<n ; ++i){
            if(!(px[i] == (int) px[i])){
                IS_INT = false;
                break;
            }
        }

        is_int_in_double = IS_INT; // true: only if x is REAL + INT test OK
    } else if(TYPEOF(x) == STRSXP){

        IS_STR = true;
        std::uintptr_t xi_uintptr;

        for(int i=0 ; i<n ; ++i){
            const char *pxi = CHAR(STRING_ELT(x, i));
            xi_uintptr = reinterpret_cast<std::uintptr_t>(pxi);

            x_ull.push_back(static_cast<unsigned long long>(xi_uintptr));

            // Rcout << xi_uintptr << "  ----  " << xi_ull << "\n";
        }
    } else {
        IS_INT = true;
    }

    if(IS_INT){
        // integer

        double max_value;
        void *px_generic;
        int X_MIN;
        if(is_int_in_double){
            double *px = REAL(x);
            px_generic = REAL(x);
            double x_min = px[0], x_max = px[0], x_tmp;
            for(int i=1 ; i<n ; ++i){
                x_tmp = px[i];
                if(x_tmp > x_max) x_max = x_tmp;
                if(x_tmp < x_min) x_min = x_tmp;
            }
            X_MIN = static_cast<int>(x_min);
            max_value = x_max - x_min;
        } else {
            int *px = INTEGER(x);
            px_generic = INTEGER(x);
            int x_min = px[0], x_max = px[0], x_tmp;
            for(int i=1 ; i<n ; ++i){
                x_tmp = px[i];
                if(x_tmp > x_max) x_max = x_tmp;
                if(x_tmp < x_min) x_min = x_tmp;
            }
            X_MIN = x_min;
            max_value = x_max - x_min;
        }

        // creating + copying a 10**5 vector takes about 0.5ms which is OK even if n << 10**5
        // quf_int is really extremely efficient, so we use it whenever appropriate
        // In quf_int_gnl we create two n-size vectors + the method is less efficient by construction
        // so we go into quf_int whenever max_value <= 2.5*n

        if(max_value < 100000 || max_value <= 2.5*n){
            quf_int(x_uf, px_generic, x_unik, X_MIN, static_cast<int>(max_value), is_int_in_double);
        } else if(max_value < 0x10000000){
            // we don't cover ranges > 2**31 (uints are pain in the neck)
            quf_int_gnl(x_uf, px_generic, x_unik, X_MIN, is_int_in_double);
        } else {
            // ranges > 2**31 => as double

            if(is_int_in_double){
                quf_double(x_uf, (double *)px_generic, x_unik);
            } else {
                // we need to create a vector of double, otherwise: pointer issue
                vector<double> x_dble(n);
                int *px = INTEGER(x);
                for(int i=0 ; i<n ; ++i) x_dble[i] = static_cast<double>(px[i]);
                quf_double(x_uf, x_dble.data(), x_unik);
            }
        }

    } else if(IS_STR){
        // string -- beforehand transformed as ULL
        quf_double(x_uf, x_ull.data(), x_unik, true);
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


void quf_single(void *px_in, std::string &x_type, int n, vector<int> &x_uf, vector<double> &x_unik){

    // preparation for strings
    bool IS_STR = x_type == "string";

    bool IS_INT = false;
    bool is_int_in_double = false;
    if(x_type == "double"){
        // we check if underlying structure is int
        IS_INT = true;
        double *px = (double *) px_in;
        for(int i=0 ; i<n ; ++i){
            if(!(px[i] == (int) px[i])){
                IS_INT = false;
                break;
            }
        }

        is_int_in_double = IS_INT; // true: only if x is REAL + INT test OK

    } else if(x_type == "int"){
        IS_INT = true;

    }

    if(IS_INT){
        // integer

        double max_value;
        void *px_generic;
        int X_MIN;
        if(is_int_in_double){
            double *px = (double *) px_in;
            px_generic = (double *) px_in;
            double x_min = px[0], x_max = px[0], x_tmp;
            for(int i=1 ; i<n ; ++i){
                x_tmp = px[i];
                if(x_tmp > x_max) x_max = x_tmp;
                if(x_tmp < x_min) x_min = x_tmp;
            }
            X_MIN = static_cast<int>(x_min);
            max_value = x_max - x_min;
        } else {
            int *px = (int *) px_in;
            px_generic = (int *) px_in;
            int x_min = px[0], x_max = px[0], x_tmp;
            for(int i=1 ; i<n ; ++i){
                x_tmp = px[i];
                if(x_tmp > x_max) x_max = x_tmp;
                if(x_tmp < x_min) x_min = x_tmp;
            }
            X_MIN = x_min;
            max_value = x_max - x_min;
        }

        // creating + copying a 10**5 vector takes about 0.5ms which is OK even if n << 10**5
        // quf_int is really extremely efficient, so we use it whenever appropriate
        // In quf_int_gnl we create two n-size vectors + the method is less efficient by construction
        // so we go into quf_int whenever max_value <= 2.5*n

        if(max_value < 100000 || max_value <= 2.5*n){
            quf_int(x_uf, px_generic, x_unik, X_MIN, static_cast<int>(max_value), is_int_in_double);
        } else if(max_value < 0x10000000){
            // we don't cover ranges > 2**31 (uints are pain in the neck)
            quf_int_gnl(x_uf, px_generic, x_unik, X_MIN, is_int_in_double);
        } else {
            // ranges > 2**31 => as double

            if(is_int_in_double){
                quf_double(x_uf, (double *)px_generic, x_unik);
            } else {
                // we need to create a vector of double, otherwise: pointer issue
                vector<double> x_dble(n);
                int *px = (int *) px_in;
                for(int i=0 ; i<n ; ++i) x_dble[i] = static_cast<double>(px[i]);
                quf_double(x_uf, x_dble.data(), x_unik);
            }
        }

    } else if(IS_STR){
        // string -- beforehand transformed as ULL
        quf_double(x_uf, px_in, x_unik, true);
    } else {
        // double
        quf_double(x_uf, px_in, x_unik);
    }

}

void quf_table_sum_single(void *px_in, std::string &x_type, int n, int q, vector<int> &x_quf,
                          vector<double> &x_unik, vector<int> &x_table,
                          double *py, vector<double> &sum_y, bool do_sum_y,
                          int type, vector<bool> &any_pblm, vector<bool> &id_pblm,
                          bool check_pblm){

    //

    // Problems depending on the types:
    // 0: no problem
    // 1: pblm if sum y = 0
    // 2: pblm if sum y = 0 or sum_y == table

    x_quf.resize(n);

    // Rcout << "q = " << q << ", n = " << n;

    // UFing
    quf_single(px_in, x_type, n, x_quf, x_unik);

    // table + sum_y

    // You need sum_y to remove clusters (even if you don't need sum_y in the end)
    bool compute_sum_y = do_sum_y || type > 0;

    int D = x_unik.size();

    x_table.resize(D);

    sum_y.resize(compute_sum_y > 0 ? D : 1);
    std::fill(sum_y.begin(), sum_y.end(), 0);

    int obs;
    for(int i=0 ; i<n ; ++i){
        obs = x_quf[i] - 1;
        ++x_table[obs];
        if(compute_sum_y) sum_y[obs] += py[i];
    }

    if(type > 0 && check_pblm){

        // 1) We check whether there is any problem
        bool is_logit = type == 2;
        int d_end = 0;
        for(int d=0 ; d<D ; ++d){
            if(sum_y[d] == 0 || (is_logit && sum_y[d] == x_table[d])){
                any_pblm[q] = true;
                d_end = d;
                break;
            }
        }

        // 2) If so: we find them
        if(any_pblm[q]){
            id_pblm.resize(D);
            std::fill(id_pblm.begin(), id_pblm.end(), false);

            for(int d=d_end ; d<D ; ++d){
                if(sum_y[d] == 0 || (is_logit && sum_y[d] == x_table[d])){
                    id_pblm[d] = true;
                }
            }

            // int sum_problems = std::accumulate(id_pblm.begin(), id_pblm.end(), 0);
            // Rcout << ", sum_problems = " << sum_problems;
        }

    }

    // Rcout << "\n";

}

void quf_refactor_table_sum_single(vector<int> &quf_old, vector<int> &quf_new, vector<bool> &obs_removed,
                                   vector<double> &x_unik, vector<double> &x_unik_new, vector<double> &x_removed,
                                   vector<int> &x_table, double *py, vector<double> &sum_y, bool do_sum_y,
                                   int type, vector<bool> &id_pblm, bool check_pblm, bool *pstop_now){
    // takes in the old quf, the observations removed,
    // then recreates the vectors of:
    // quf, table, sum_y, x_unik_new

    // IDs range from 1 to n_values

    // Type: 0: no restriction -- 1: only 0s -- 2: 0s and 1s

    int D = x_unik.size();
    int n = quf_old.size();

    // Only if type == 2 (0s and 1s removed) that we recreate the pblmatic IDs
    // because if only 0s, only the ones in the original id_pblm are removed
    // we find out which ones are only 0s or only 1s
    // NOTA: I don't reiterate the observation removal => just one round
    //  (otherwise it could go several rounds => too time consuming, no big value added)

    if((type == 2 || !check_pblm) && !*pstop_now){
        // If we haven't excluded observations because of 0 only values (ie check_pblm == false)
        // then we here need to update id_pblm, because some IDs could have been taken
        // out due to other clusters being removed

        // the new ID problem => what are the IDs that still exists?
        vector<bool> id_still_exists(D, false);

        // We also check that some observations are still left!

        // we recompute the table and sum_y
        std::fill(x_table.begin(), x_table.end(), 0);
        std::fill(sum_y.begin(), sum_y.end(), 0);

        int id;
        for(int i=0 ; i<n ; ++i){
            if(!obs_removed[i]){
                id = quf_old[i] - 1;
                if(!id_still_exists[id]) id_still_exists[id] = true;
                ++x_table[id];
                sum_y[id] += py[i];
            }
        }

        // recreating id_pblm
        id_pblm.resize(D); // we resize because we could have had length = 0
        for(int d=0 ; d<D ; ++d) id_pblm[d] = !id_still_exists[d];

        // Loop finding out the problems
        vector<bool> id_pblm_check(D, false);
        for(int d=0 ; d<D ; ++d){
            if(sum_y[d] == 0 || sum_y[d] == x_table[d]){
                id_pblm_check[d] = true;
            }
        }

        int sum_problems = std::accumulate(id_pblm_check.begin(), id_pblm_check.end(), 0);
        // Rcout << "sum_problems = " << sum_problems << "\n";
        *pstop_now = sum_problems == D;

    }

    if(*pstop_now == false){

        // is there a problem?
        bool any_pblm = false;
        // we check only if needed (if type<2 && length(id_pblm) == 0 => means no problem)
        if(id_pblm.size() > 0){
            for(int d=0 ; d<D ; ++d){
                if(id_pblm[d]){
                    any_pblm = true;
                    break;
                }
            }
        }

        // Creating the new IDs
        vector<int> id_new;
        int nb_pblm = 0;
        if(any_pblm){
            // we create these new IDs only if there is a pblm
            // if no problem then the elements of quf_new are the same as in quf_old

            id_new.resize(D);
            std::fill(id_new.begin(), id_new.end(), 0);

            for(int d=0 ; d<D ; ++d){
                if(id_pblm[d]){
                    ++nb_pblm;
                } else {
                    id_new[d] = 1 + d - nb_pblm;
                }
            }
        }

        // New quf + table + sum_y

        int D_new = D - nb_pblm;

        x_table.resize(D_new);
        std::fill(x_table.begin(), x_table.end(), 0);
        if(do_sum_y){
            sum_y.resize(D_new);
            std::fill(sum_y.begin(), sum_y.end(), 0);
        }

        int id;
        for(int i=0 ; i<n ; ++i){
            if(!obs_removed[i]){
                id = any_pblm ? id_new[quf_old[i] - 1] : quf_old[i];
                ++x_table[id - 1];
                if(do_sum_y) sum_y[id - 1] += py[i];
                quf_new.push_back(id);
            }
        }

        // x_unik + x_removed

        if(any_pblm){
            // we create x_unik_new and x_removed

            for(int d=0 ; d<D ; ++d){
                if(id_pblm[d]){
                    x_removed.push_back(x_unik[d]);
                } else {
                    x_unik_new.push_back(x_unik[d]);
                }
            }
        } else {
            // x_unik_new is identical to x_unik
            x_unik_new = x_unik;
        }
    }
}


// [[Rcpp::export]]
List cpppar_quf_table_sum(SEXP x, SEXP y, bool do_sum_y, int type, IntegerVector only_slope, int nthreads){
    // x: List of vectors of IDs (type int/num or char only)
    // y: dependent variable
    // Type:
    //  0: no removal of obs
    //  1: removal of 0 only clusters
    //  2: removal of 0 only and 1 only clusters
    // do_sum_y: should we compute the sum_y?

    if(type < 0 || type > 2) stop("Argument type must be equal to 0, 1, or 2.");

    int Q = Rf_length(x);
    SEXP xq = VECTOR_ELT(x, 0);
    int n = Rf_length(xq);

    vector<bool> check_pblm(Q, true);
    if(only_slope.length() == Q){
        for(int q=0 ; q<Q ; ++q){
            // we check for pblm only if NOT only slope
            check_pblm[q] = only_slope[q] == false;
        }
    }

    // Rcout << "Q = " << Q << "\n";

    double *py = REAL(y);

    // the vectors of qufed
    vector< vector<int> > x_quf_all(Q);
    vector< vector<int> > x_table_all(Q);
    vector< vector<double> > x_unik_all(Q);
    vector<bool> any_pblm(Q, false);
    vector< vector<bool> > id_pblm_all(Q);
    // The following may not be needed:
    vector< vector<double> > sum_y_all(Q);
    vector<bool> obs_removed;
    vector< vector<double> > x_removed_all(Q);

    // New code => we avoid passing SEXP within the parallel loop
    // We get the types of the R vectors
    // For strings => we convert into numeric without parallel (sigh)

    vector<void *> px_all(Q);
    vector<std::string> x_type_all(Q);
    // vector to store modified strings
    vector< vector<unsigned long long> > x_ull_all(Q);

    for(int q=0 ; q<Q ; ++q){

        xq = VECTOR_ELT(x, q);

        if(TYPEOF(xq) == INTSXP){
            x_type_all[q] = "int";
            px_all[q] = INTEGER(xq);

        } else if(TYPEOF(xq) == REALSXP){
            x_type_all[q] = "double";
            px_all[q] = REAL(xq);

        } else if(TYPEOF(xq) == STRSXP){
            // We make the conversion to unsigned long long
            x_type_all[q] = "string";

            std::uintptr_t xi_uintptr;

            for(int i=0 ; i<n ; ++i){
                const char *pxi = CHAR(STRING_ELT(xq, i));
                xi_uintptr = reinterpret_cast<std::uintptr_t>(pxi);

                x_ull_all[q].push_back(static_cast<unsigned long long>(xi_uintptr));

                // Rcout << xi_uintptr << "  ----  " << xi_ull << "\n";
            }

            px_all[q] = x_ull_all[q].data();
        } else {
            // We NEVER end here
            stop("Error: wrong type.");
        }
    }



#pragma omp parallel for num_threads(nthreads)
    for(int q=0 ; q<Q ; ++q){
        quf_table_sum_single(px_all[q], x_type_all[q], n, q, x_quf_all[q], x_unik_all[q], x_table_all[q],
                             py, sum_y_all[q], do_sum_y, type, any_pblm, id_pblm_all[q], check_pblm[q]);
    }


    bool is_pblm = false;
    for(int q=0 ; q<Q ; ++q){
        if(any_pblm[q]){
            is_pblm = true;
            break;
        }
    }

    if(is_pblm){

        // creating the obs2remove vector
        obs_removed.resize(n);
        std::fill(obs_removed.begin(), obs_removed.end(), false);

        // No need to take care of race conditions
#pragma omp parallel for num_threads(nthreads)
        for(int q=0 ; q<Q ; ++q){
            // skipping loop if no problem
            if(any_pblm[q]){
                vector<bool> id_pblm = id_pblm_all[q];
                int *pquf = x_quf_all[q].data();
                for(int i=0 ; i<n ; ++i){
                    if(id_pblm[pquf[i] - 1]){
                        obs_removed[i] = true;
                    }
                }
            }
        }

        // int nb_rm = 0;
        // for(int i=0 ; i<n ; ++i){
        //     nb_rm += obs_removed[i];
        // }
        //
        // Rcout << "#obs removed = " << nb_rm << "\n";
        // stop("end after removing obs.");


        // refactoring, recomputing all the stuff
        // ie x_quf, x_unik, x_table, sum_y (if needed)
        // but also: x_removed

        // The values of x_table and sum_y are changed in place
        // I create the addition new_quf and new_unik because we need the
        // them old items for the construction of the new.

        vector< vector<int> > x_new_quf_all(Q);
        vector< vector<double> > x_new_unik_all(Q);
        bool stop_now = false;
        bool *pstop_now = &stop_now;

#pragma omp parallel for num_threads(nthreads)
        for(int q=0 ; q<Q ; ++q){
            quf_refactor_table_sum_single(x_quf_all[q], x_new_quf_all[q], obs_removed,
                                          x_unik_all[q], x_new_unik_all[q], x_removed_all[q],
                                                                                         x_table_all[q], py, sum_y_all[q], do_sum_y,
                                                                                         type, id_pblm_all[q], check_pblm[q], pstop_now);
        }

        if(*pstop_now){
            stop("The dependent variable is fully explained by the fixed-effects.");
        }

        // Saving the values in the appropriate locations
        x_quf_all = x_new_quf_all;
        x_unik_all = x_new_unik_all;

    }

    // The object to be returned
    List res;

    // Rcout << "start res -- ";

    // x: UF
    List res_tmp;
    for(int q=0 ; q<Q ; ++q){
        res_tmp.push_back(x_quf_all[q]);
    }
    res["quf"] = clone(res_tmp);

    // Rcout << "quf: ok";

    // x: Unik
    for(int q=0 ; q<Q ; ++q){
        res_tmp[q] = x_unik_all[q];
    }
    res["items"] = clone(res_tmp);

    // Rcout << ", items: ok";

    // table
    for(int q=0 ; q<Q ; ++q){
        res_tmp[q] = x_table_all[q];
    }
    res["table"] = clone(res_tmp);

    // Rcout << ", table: ok";

    // sum y
    for(int q=0 ; q<Q ; ++q){
        if(do_sum_y){
            res_tmp[q] = sum_y_all[q];
        } else {
            res_tmp[q] = 0;
        }
    }
    res["sum_y"] = clone(res_tmp);

    // Rcout << ", sum_y: ok";

    //
    // IF PROBLEM ONLY
    //

    if(is_pblm){
        // The removed observations
        res["obs_removed"] = obs_removed;

        // Rcout << ", obs_removed: ok";

        // The removed clusters
        for(int q=0 ; q<Q ; ++q){
            res_tmp[q] = x_removed_all[q];
        }
        res["fe_removed"] = clone(res_tmp);

        // Rcout << ", fe_removed: ok";

    }

    // Rcout << "\n";

    return res;
}




