


.onLoad <- function(libname, pkgname){
	# setting some options

	options("fixest_dict" = c())
	options("fixest_notes" = TRUE)
	options("fixest_print.type" = "table")
	options("fixest_fl_authorized" = FALSE)

	setFixest_coefplot("all", reset = TRUE)
	setFixest_dof()
	setFixest_nthreads()

	invisible()
}

