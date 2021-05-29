


.onLoad <- function(libname, pkgname){
	# setting some options

	options("fixest_dict" = c())
	options("fixest_notes" = TRUE)
	options("fixest_print" = list(type = "table"))
	options("fixest_fl_authorized" = FALSE)

	setFixest_coefplot("all", reset = TRUE)
	setFixest_dof()

	# # To include later
	# cpp_setup_fork_presence()

	# nthreads
	nthreads = Sys.getenv("fixest_nthreads")
	if(identical(nthreads, "")) nthreads = NULL
	if(!is.null(nthreads)){
	    # we have extra strong checking in setFixest_nthreads
	    # so here we create a warning instead if nthreads is not valid
	    if(length(nthreads) != 1 || !isScalar(nthreads) || nthreads < 0 ||
	       !(nthreads %% 1 == 0 || nthreads < 1)){
            warning("Environment variable 'fixest_nthreads' not valid, using the default value of setFixest_nthreads() instead.", call. = FALSE)
	        nthreads = NULL
	    }
	}

	setFixest_nthreads(nthreads)



	invisible()
}



.onAttach = function(libname, pkgname) {

    # breaking message: don't know how long I'll keep it
    is_msg = !isFALSE(Sys.getenv("fixest_startup_msg"))
    if(is_msg) packageStartupMessage("fixest 0.9.0, BREAKING changes! \nin i():\n  a) the two first arguments have been swapped! Now it's i(factor_var, continuous_var) for interactions. \n  b) argument 'drop' has been removed (put everything in 'ref' now).")

}


