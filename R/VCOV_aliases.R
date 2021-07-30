# Do not edit by hand
# => aliases some VCOV functions




    #' @rdname vcov_hac
    NW = function(lag = NULL){
        extra_args = list(lag = lag)
        vcov_request = list(vcov = "NW", extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
        vcov_request
    }

    #' @rdname vcov_hac
    newey_west = function(lag = NULL){
        extra_args = list(lag = lag)
        vcov_request = list(vcov = "NW", extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
        vcov_request
    }

    #' @rdname vcov_hac
    DK = function(lag = NULL){
        extra_args = list(lag = lag)
        vcov_request = list(vcov = "DK", extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
        vcov_request
    }

    #' @rdname vcov_hac
    driscoll_kraay = function(lag = NULL){
        extra_args = list(lag = lag)
        vcov_request = list(vcov = "DK", extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
        vcov_request
    }

    #' @rdname vcov_conley
    conley = function(cutoff = NULL, pixel = NULL, distance = NULL){
        extra_args = list(cutoff = cutoff, pixel = pixel, distance = distance)
        vcov_request = list(vcov = "conley", extra_args = extra_args)
        class(vcov_request) = "fixest_vcov_request"
        vcov_request
    }
