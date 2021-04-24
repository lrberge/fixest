


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
	setFixest_nthreads()

	# Setup of builtin VCOVs
	vcov_setup()

	invisible()
}

