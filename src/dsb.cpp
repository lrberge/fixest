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

inline bool is_separator(const char * str, int i){
    return str[i] == '!' || str[i] == '?';
}


void extract_operator(const char * str, int &i, int n,
                      std::vector<std::string> &operator_list, bool no_whitespace = false){
    // modifies the operator and gives the updated i
    // the i should start where to evaluate post separator (if present)

    char quote = '\'';

    // we first get the operator, if there is an operator
    std::string operator_tmp = "";

    int i_start = i;
    bool any_operator = true;

    if(str[i] == '/'){
        // special case!!!!
        operator_tmp = "/";
        operator_list.push_back(operator_tmp);
        ++i;

    } else {
        // normal case

        bool is_comma = false;

        if(str[i] == ' '){
            if(no_whitespace){
                // not an operator. i = n to avoid the loop
                i = n;
                any_operator = false;
            } else {
                // skip the first whitespaces
                while(i < n && str[i] == ' ') ++i;
            }
        }

        while(i < n && !is_separator(str, i)){
            if(is_quote(str, i)){
                is_comma = false;
                if(operator_tmp.length() > 0){
                    // we save the existing command if needed
                    operator_list.push_back(operator_tmp);
                    operator_tmp = "";
                }

                // we get the full quoted value
                quote = str[i];
                operator_tmp += quote;
                ++i;

                while(i < n && str[i] != quote){
                    operator_tmp += str[i];
                    ++i;
                }

                if(i < n){
                    operator_tmp += quote;
                    ++i;
                }
            } else if(is_dsb_bound(str, i, n)){
                // there should be no dsb bound in the operator
                // if so, this is an error

                any_operator = false;
                break;

            } else {

                // comma: separation between operations
                if(str[i] == ','){
                    if(operator_tmp.length() > 0){
                        operator_list.push_back(operator_tmp);
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
                operator_list.push_back(operator_tmp);
                operator_tmp = "";
            }

            if(i < n){
                operator_tmp = str[i];
                operator_list.push_back(operator_tmp);
                ++i;
            } else {
                // end of the string with no separator => problem
                any_operator = false;
            }
        }

        if(!any_operator){
            std::vector<std::string> empty_list;
            operator_list = empty_list;
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
            std::vector<std::string> operator_list;

            // modifies i and operator_list "in place"
            extract_operator(str, i, n, operator_list);

            dsb_element.push_back(operator_list);

            // init
            dsb_value = "";

            // we now get the value to be evaluated or treated as verbatim
            while(i < n && n_open > 0){
                // if(is_dsb_open(str, i, n)){
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
    std::vector<std::string> operator_list;

    // modifies i and operator_list "in place"
    int i = 0;
    extract_operator(str, i, n, operator_list, true);

    dsb_element.push_back(operator_list);

    // init
    std::string dsb_value = "";

    // Remember that the full string is verbatim
    for(; i < n ; ++i){
        dsb_value += str[i];
    }

    dsb_element.push_back(dsb_value);

    return dsb_element;
}




