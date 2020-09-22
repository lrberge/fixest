#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Tue Sep 22 15:30:09 2020
# ~: Code snippets to share
#----------------------------------------------#




####
#### Imai & Song ####
####

# Equivalence between the "unit" wfe estimation and fixest
# A) function to compute the unit weights
# B) comparaison of the estimations
#

#
# A) weights
#

# Function to compute the weights (with full checks)

require(dreamerr) ; require(data.table)
IS_weights_unit = function(index, base){

    # ARGUMENT CHECK #

    check_arg(index, "mbt data.frame ncol(2) | os formula | character vector len(2)",
              .message = "The argument 'index' must contain the a) unit and b) the treatment variables. It can be either: i) a data.frame of two columns, ii) a one sided formula, or iii) a character vector of length 2.")

    if(!is.data.frame(index)){
        check_arg(base, "data.frame mbt", .message = "If the argument 'index' is not a data.frame, you must provide a data.frame in the argument 'base'.")

        if("formula" %in% class(index)){
            check_arg(index, "formula var(data)", .data = base)

            index = model.frame(index, base)

            check_arg(index, "data.frame ncol(2)", .message = "If a formula, then 'index' must lead to a data.frame with two variables.")

        } else {
            # flexible matching of names
            check_arg_plus(index, "multi match", .choices = names(base), .message = "If the argument 'index' is a character vector, it must match (at least partially) the names of the data.frame in argument 'base'.")

            index = as.data.frame(base)[, index]
        }

    }

    # Here 'index' is a data.frame with two variables
    res = as.data.table(index)
    names(res) = c("unit", "treatment")

    res[, n_treated := sum(treatment), by = unit]
    res[, n_not_treated := sum(1 - treatment), by = unit]
    res[, m_size := treatment * n_not_treated + (1 - treatment) * n_treated]
    res[, w := 1 + treatment * (m_size / n_treated) + (1 - treatment) * (m_size / n_not_treated)]
    res[is.na(w), w := 0]

    res$w
}


#
# B) Estimations
#

library(fixest) ; library(wfe)

# Data preparation

data(base_did)
set.seed(0)
base = data.table(base_did)
# we trim first/last periods of some guys to add variation in the weights (otherwise: all equal to 2 or 0)
id_first = sample(108, 20) ; id_last = sample(108, 20)
base = base[!(id %in% id_first & period <= 3) & !(id %in% id_last & period >= 7)]
base[, treat_post := treat*post]

# Estimations

system.time(res_wfe <- wfe(y ~ treat_post + x1, data = base, treat = "treat_post", unit.index = "id", time.index = "period", method = "unit",
                           qoi = "ate", hetero.se=TRUE, auto.se=TRUE))


system.time(res_feols <- feols(y ~ treat_post + x1 | id, base, weights = w <- IS_weights_unit(~ id + treat_post, base)))

# Comparing the weights
all.equal(w, summary(res_wfe)$W$W.it)

# Timings on my laptop:
#    wfe: 320ms
# fixest:  10ms


















































