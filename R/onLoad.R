


.onLoad <- function(libname, pkgname){
	# setting the two options

	options("fixest_dict" = c())
	options("fixest_notes" = TRUE)
	options("fixest_na_inf.rm" = TRUE)
	options("fixest_print.type" = "table")
	options("fixest_fl_authorized" = FALSE)
	setFixest_nthreads()

	invisible()
}

