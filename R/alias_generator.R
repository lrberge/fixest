#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Jul 02 11:58:49 2021
# ~: Generate aliases
#----------------------------------------------#


# The problem with aliased functions is that when the main-function changes its
# arguments, then we also need to make the changes in the aliased ones, we usually need to
# make those changes by hand, incurring substantial maintenance costs and highly susceptible to errors
#
# It is much safer to generate aliased functions automatically
# This cuts the maintenance costs to almost 0 and avoids any mistake
#



create_aliases = function(){
    # This function is triggered only at loading time
    # BUT ONLY FOR ME, THE DEVELOPPER!!!!!
    # The function knows it's "me" with the system environment variable fixest_ROOT

    if(!isTRUE(renvir_get("fixest_ROOT"))) return(NULL)
    # we check we're in the right directory (otherwise there can be prblms with Rmakdown)
    if(!isTRUE(file.exists("R/VCOV_aliases.R"))) return(NULL)

    gen_etable_aliases()
    gen_iplot()
    gen_vcov_aliases()

}





































































































