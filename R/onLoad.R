


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
	setFixest_nthreads()

	invisible()
}



.onAttach = function(libname, pkgname) {

    # breaking message: don't know how long I'll keep it

    do_msg = initialize_startup_msg()

    is_msg = do_msg || !isFALSE(renvir_get("fixest_startup_msg"))

    if(is_msg) packageStartupMessage("fixest 0.9.0, BREAKING changes! (Permanently remove this message with fixest_startup_msg(FALSE).) \n- In i():\n    + the first two arguments have been swapped! Now it's i(factor_var, continuous_var) for interactions. \n    + argument 'drop' has been removed (put everything in 'ref' now).\n- In feglm(): \n    + the default family becomes 'gaussian' to be in line with glm(). Hence, for Poisson estimations, please use fepois() instead.")

}


