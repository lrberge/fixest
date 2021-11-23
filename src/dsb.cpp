#include <Rcpp.h>
#include <string>
#include <cstring>
#include <vector>
using namespace Rcpp;


inline bool is_dsb_open(const char * str, int i, int n){
    if(i + 2 < n){
        return str[i] == '.' && str[i + 1] == '[';
    } else {
        return false;
    }
}

inline bool is_dsb_bound(const char * str, int i, int n){
    if(str[i] == ']'){
        return true;
    } else if(i + 2 < n){
        return str[i] == '.' && str[i + 1] == '[';
    } else {
        return false;
    }
}

inline bool is_special_char(const char * str, int i){
    return str[i] == '\'' || str[i] == '"' || str[i] == '`' || str[i] == ':' ||
        str[i] == ';' || str[i] == '/';
}

inline bool is_quote(const char * str, int i){
    return str[i] == '\'' || str[i] == '"' || str[i] == '`';
}

inline bool is_basic_quote(const char * str, int i){
    return str[i] == '\'' || str[i] == '"';
}

inline bool is_separator(const char * str, int i){
    return str[i] == '!' || str[i] == '?';
}


void extract_quote(const char * str, int &i, int n,
                   std::string &operator_tmp){

    char quote = str[i++];
    operator_tmp += quote;

    while(i < n && str[i] != quote){
        operator_tmp += str[i++];
    }

    if(i < n){
        operator_tmp += quote;
        ++i;
    }
}

void extract_operator(const char * str, int &i, int n,
                      std::vector<std::string> &operator_vec,
                      bool &is_eval, bool no_whitespace = false){
    // modifies the operator and gives the updated i
    // the i should start where to evaluate post separator (if present)

    // we first get the operator, if there is an operator
    std::string operator_tmp = "";

    int i_start = i;
    bool any_operator = true;

    if(str[i] == '/'){
        // special case!!!!
        operator_tmp = "/";
        operator_vec.push_back(operator_tmp);
        ++i;

    } else {
        // normal case

        bool is_comma = false;

        if(str[i] == ' '){
            if(no_whitespace){
                // not an operator.
                // i = n to avoid the loop
                i = n;
                any_operator = false;
            } else {
                // skip the first whitespaces
                while(i < n && str[i] == ' ') ++i;
            }
        }

        while(i < n && !is_separator(str, i)){
            if(str[i] == '('){
                // we deal with quotes
                // Used (so far) only in if statements
                if(i >= 2 && str[i - 1] == 'f' && str[i - 2] == 'i'){

                    // There is no such thing as nested if(),
                    // so if we find another '(' it's a problem

                    operator_tmp += '(';
                    ++i;

                    while(i < n && str[i] != ')' && !is_separator(str, i)){
                        // we don't really care about parsing here
                        // we'll do that later

                        if(is_quote(str, i)){
                            extract_quote(str, i, n, operator_tmp);
                        } else {
                            // any other stuff gets in
                            operator_tmp += str[i];
                            ++i;
                        }
                    }

                    if(i < n && str[i] != ')'){
                        any_operator = false;
                        break;
                    }

                } else {
                    // Otherwise => problem
                    any_operator = false;
                    break;
                }
            }

            if(is_quote(str, i)){
                is_comma = false;
                if(operator_tmp.length() > 0){
                    // we save the existing command if needed
                    operator_vec.push_back(operator_tmp);
                    operator_tmp = "";
                }

                // we get the full quoted value
                extract_quote(str, i, n, operator_tmp);

            } else if(is_dsb_bound(str, i, n)){
                // there should be no dsb bound in the operator
                // if so, this is an error

                any_operator = false;
                break;

            } else {

                // comma: separation between operations
                if(str[i] == ','){
                    if(operator_tmp.length() > 0){
                        operator_vec.push_back(operator_tmp);
                        operator_tmp = "";
                    }
                    is_comma = true;
                } else if(str[i] == ' '){
                    // nothing, but if NOT after a comma => error
                    if(!is_comma){
                        // if the spaces are only trailing, OK
                        while(i < n && str[i] == ' ') ++i;
                        if(is_separator(str, i)){
                            // OK
                            break;
                        } else {
                            any_operator = false;
                            break;
                        }
                    }

                } else {
                    is_comma = false;
                    operator_tmp += str[i];
                }

                ++i;
            }
        }

        if(any_operator){
            if(operator_tmp.length() > 0){
                operator_vec.push_back(operator_tmp);
                operator_tmp = "";
            }

            if(i < n){
                // we end with a separator
                is_eval = str[i] == '?';
                operator_tmp = str[i];
                operator_vec.push_back(operator_tmp);
                ++i;
            } else {
                // end of the string with no separator => problem
                any_operator = false;
            }
        }

        if(!any_operator){
            std::vector<std::string> empty_vec;
            operator_vec = empty_vec;
            i = i_start;
        }
    }
}

// [[Rcpp::export]]
List cpp_dsb(SEXP Rstr){
    // Rstr: string from R of length 1

    List res;
    const char *str = CHAR(STRING_ELT(Rstr, 0));

    // DSB open flag
    int n_open = 0;

    std::string string_value = "";
    std::string dsb_value = "";

    int n = std::strlen(str);

    int i = 0;
    while(i < n){

        // if not currently open => we check until open
        if(n_open == 0){
            while(i < n && !is_dsb_open(str, i, n)){
                string_value += str[i];
                ++i;
            }

            res.push_back(string_value);

            if(i < n){
                // there was one open
                i += 2; // we increment i bc the opening string is 2 char long
                ++n_open;
                string_value = "";
            }

        } else {

            List dsb_element;
            std::vector<std::string> operator_vec;

            // modifies i and operator_vec "in place"
            bool is_eval = true;
            extract_operator(str, i, n, operator_vec, is_eval);

            dsb_element.push_back(operator_vec);

            // init
            dsb_value = "";

            // we now get the value to be evaluated or treated as verbatim

            if(is_eval){
                // we strip all the white spaces
                while(i < n && str[i] == ' ') ++i;

                if(is_basic_quote(str, i)){
                    // we take verbatim the full quote
                    // this means it can contain .[] without issue

                    char quote = str[i++];
                    dsb_value += quote;

                    while(i < n && str[i] != quote) dsb_value += str[i++];

                    if(i < n){
                        dsb_value += quote;
                        ++i;
                        // we then let it go, if there is a parsing problem, we'll
                        // spot that directly in R, no point to
                        // doing thorough error handling here
                    }
                }
            }

            while(i < n && n_open > 0){

                if(str[i] == '['){
                    ++n_open;
                } else if(str[i] == ']'){
                    --n_open;
                }

                if(n_open == 0){
                    ++i;
                    break;
                }

                dsb_value += str[i];
                ++i;
            }

            dsb_element.push_back(dsb_value);

            res.push_back(dsb_element);
        }
    }

    return res;
}

// [[Rcpp::export]]
List cpp_dsb_full_string(SEXP Rstr){
    // When we consider the full string a verbatim within dsb

    const char *str = CHAR(STRING_ELT(Rstr, 0));

    int n = std::strlen(str);

    List dsb_element;
    std::vector<std::string> operator_vec;

    // is eval is not used here but is required in extract_operator
    bool is_eval = false;

    // modifies i and operator_vec "in place"
    int i = 0;
    extract_operator(str, i, n, operator_vec, is_eval, true);

    dsb_element.push_back(operator_vec);

    // init
    std::string dsb_value = "";

    // Remember that the full string is verbatim
    for(; i < n ; ++i){
        dsb_value += str[i];
    }

    dsb_element.push_back(dsb_value);

    return dsb_element;
}

inline bool is_if_separator(const char * str, int i, int n, bool semicolon = false){
    if(semicolon){
        return i < n && str[i] == ':';
    } else {
        return i >= n;
    }
}

// [[Rcpp::export]]
List cpp_dsb_if_extract(SEXP Rstr){

    const char *str = CHAR(STRING_ELT(Rstr, 0));
    int n = std::strlen(str);

    List if_elements;
    std::vector<std::string> operator_vec;
    std::vector<std::string> empty_vec;
    std::string operator_tmp = "";

    // the code is close to extract_operator, but not identical....
    // it's a bit code duplication and I really don't like that
    // but otherwise, the code would get more ugly in extract_operator...

    bool is_comma = false;
    bool any_problem = false;

    // true
    int i = 0;
    int n_loop = 0;

    while(n_loop++ < 2){

        // 1st loop: semicolon is separator
        bool semicolon = n_loop == 1;
        // 2nd loop: EOL is separator

        while(i < n && !is_if_separator(str, i, n, semicolon)){

            if(is_quote(str, i)){
                is_comma = false;
                if(operator_tmp.length() > 0){
                    // we save the existing command if needed
                    operator_vec.push_back(operator_tmp);
                    operator_tmp = "";
                }

                // we get the full quoted value
                extract_quote(str, i, n, operator_tmp);

            } else {

                // comma: separation between operations
                if(str[i] == ','){
                    if(operator_tmp.length() > 0){
                        operator_vec.push_back(operator_tmp);
                        operator_tmp = "";
                    }
                    is_comma = true;
                } else if(str[i] == ' '){
                    // nothing, but if NOT after a comma => error
                    if(!is_comma){
                        // if the spaces are only trailing, OK
                        while(i < n && str[i] == ' ') ++i;
                        if(is_if_separator(str, i, n, semicolon)){
                            // OK
                            break;
                        } else if(i < n){
                            any_problem = true;
                            break;
                        }
                    }

                } else {
                    is_comma = false;
                    operator_tmp += str[i];
                }

                ++i;
            }
        }

        if(any_problem){
            List error;
            error.push_back(false);
            return(error);
        }

        while(i < n && (str[i] == ' ' || str[i] == ':')) ++i;

        if(operator_tmp.length() > 0){
            operator_vec.push_back(operator_tmp);
            operator_tmp = "";
        }

        if_elements.push_back(operator_vec);
        operator_vec = empty_vec;
    }


    return if_elements;
}



// [[Rcpp::export]]
StringVector cpp_paste_conditional(StringVector x, IntegerVector id, int n){

    StringVector res(n);

    int n_x = x.length();

    if(n_x == 0){
        return res;
    }

    std::string tmp = "";
    int id_current = id[0];

    for(int i=0 ; i<n_x ; ++i){
        if(id[i] == id_current){
            tmp += x[i];
        } else {
            res[id_current - 1] = tmp;
            tmp = "";
            id_current = id[i];
        }
    }

    // don't forget the last item
    res[id[n_x - 1] - 1] = tmp;

    return res;
}

