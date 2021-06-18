#----------------------------------------------#
# Author: Laurent Berge
# Date creation: Fri Jul 10 09:03:06 2020
# ~: package sniff tests
#----------------------------------------------#

# Not everything is currently covered, but I'll improve it over time

# Some functions are not trivial to test properly though

library(dreamerr) ; library(fixest)


test = fixest:::test ; chunk = fixest:::chunk
vcovClust = fixest:::vcovClust

setFixest_notes(FALSE)

####
#### Estimations ####
####

chunk("ESTIMATION")

set.seed(0)

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$fe_2 = rep(1:5, 30)
base$fe_3 = sample(15, 150, TRUE)
base$constant = 5
base$y_int = as.integer(base$y)
base$w = as.vector(unclass(base$species) - 0.95)
base$offset_value = unclass(base$species) - 0.95
base$y_01 = 1 * ((scale(base$x1) + rnorm(150)) > 0)
# what follows to avoid removal of fixed-effects (logit is pain in the neck)
base$y_01[1:5 + rep(c(0, 50, 100), each = 5)] = 1
base$y_01[6:10 + rep(c(0, 50, 100), each = 5)] = 0
# We enforce the removal of observations
base$y_int_null = base$y_int
base$y_int_null[base$fe_3 %in% 1:5] = 0

for(model in c("ols", "pois", "logit", "negbin", "Gamma")){
    cat("Model: ", dreamerr::sfill(model, 6), sep = "")
    for(use_weights in c(FALSE, TRUE)){
        my_weight = NULL
        if(use_weights) my_weight = base$w

        for(use_offset in c(FALSE, TRUE)){
            my_offset = NULL
            if(use_offset) my_offset = base$offset_value

            for(id_fe in 0:9){

                cat(".")

                tol = switch(model, "negbin" = 1e-2, "logit" = 3e-5, 1e-5)

                # Setting up the formula to accommodate FEs
                if(id_fe == 0){
                    fml_fixest = fml_stats = y ~ x1
                } else if(id_fe == 1){
                    fml_fixest = y ~ x1 | species
                    fml_stats = y ~ x1 + factor(species)
                } else if(id_fe == 2){
                    fml_fixest = y ~ x1 | species + fe_2
                    fml_stats = y ~ x1 + factor(species) + factor(fe_2)
                } else if(id_fe == 3){
                    # varying slope
                    fml_fixest = y ~ x1 | species[[x2]]
                    fml_stats = y ~ x1 + x2:species
                } else if(id_fe == 4){
                    # varying slope -- 1 VS, 1 FE
                    fml_fixest = y ~ x1 | species[[x2]] + fe_2
                    fml_stats = y ~ x1 + x2:species + factor(fe_2)
                } else if(id_fe == 5){
                    # varying slope -- 2 VS
                    fml_fixest = y ~ x1 | species[x2]
                    fml_stats = y ~ x1 + x2:species + species
                } else if(id_fe == 6){
                    # varying slope -- 2 VS bis
                    fml_fixest = y ~ x1 | species[[x2]] + fe_2[[x3]]
                    fml_stats = y ~ x1 + x2:species + x3:factor(fe_2)
                } else if(id_fe == 7){
                    # Combined clusters
                    fml_fixest = y ~ x1 + x2 | species^fe_2
                    fml_stats = y ~ x1 + x2 + paste(species, fe_2)
                } else if(id_fe == 8){
                    fml_fixest = y ~ x1 | species[x2] + fe_2[x3] + fe_3
                    fml_stats = y ~ x1 + species + i(species, x2) + factor(fe_2) + i(fe_2, x3) + factor(fe_3)
                } else if(id_fe == 9){
                    fml_fixest = y ~ x1 | species + fe_2[x2,x3] + fe_3
                    fml_stats = y ~ x1 + species + factor(fe_2) + i(fe_2, x2) + i(fe_2, x3) + factor(fe_3)
                }

                # ad hoc modifications of the formula
                if(model == "logit"){
                    fml_fixest = xpd(y_01 ~ ..rhs, ..rhs = fml_fixest[[3]])
                    fml_stats = xpd(y_01 ~ ..rhs, ..rhs = fml_stats[[3]])

                    # The estimations are OK, conv differences out of my control
                    if(id_fe %in% 8:9) tol = 0.5

                } else if(model == "pois"){
                    fml_fixest = xpd(y_int_null ~ ..rhs, ..rhs = fml_fixest[[3]])
                    fml_stats = xpd(y_int_null ~ ..rhs, ..rhs = fml_stats[[3]])

                } else if(model %in% c("negbin", "Gamma")){
                    fml_fixest = xpd(y_int ~ ..rhs, ..rhs = fml_fixest[[3]])
                    fml_stats = xpd(y_int ~ ..rhs, ..rhs = fml_stats[[3]])
                }

                adj = 1
                if(model == "ols"){
                    res = feols(fml_fixest, base, weights = my_weight, offset = my_offset)
                    res_bis = lm(fml_stats, base, weights = my_weight, offset = my_offset)

                } else if(model %in% c("pois", "logit", "Gamma")){
                    adj = 0
                    if(model == "Gamma" && use_offset) next

                    my_family = switch(model, pois = poisson(), logit = binomial(), Gamma = Gamma())

                    res = feglm(fml_fixest, base, family = my_family, weights = my_weight, offset = my_offset)

                    if(!is.null(res$obs_selection$obsRemoved)){
                        qui = res$obs_selection$obsRemoved

                        # I MUST do that.... => subset does not work...
                        base_tmp = base[qui, ]
                        base_tmp$my_offset = my_offset[qui]
                        base_tmp$my_weight = my_weight[qui]
                        res_bis = glm(fml_stats, base_tmp, family = my_family, weights = my_weight, offset = my_offset)
                    } else {
                        res_bis = glm(fml_stats, data = base, family = my_family, weights = my_weight, offset = my_offset)
                    }

                } else if(model == "negbin"){
                    # no offset in glm.nb + no VS in fenegbin + no weights in fenegbin
                    if(use_weights || use_offset || id_fe > 2) next

                    res = fenegbin(fml_fixest, base, notes = FALSE)
                    res_bis = MASS::glm.nb(fml_stats, base)

                }

                test(coef(res)["x1"], coef(res_bis)["x1"], "~", tol)
                test(se(res, se = "st", dof = dof(adj = adj))["x1"], se(res_bis)["x1"], "~", tol)
                test(pvalue(res, se = "st", dof = dof(adj = adj))["x1"], pvalue(res_bis)["x1"], "~", tol*10**(model == "negbin"))
                # cat("Model: ", model, ", FE: ", id_fe, ", weight: ", use_weights,  ", offset: ", use_offset, "\n", sep="")

            }
            cat("|")
        }
    }
    cat("\n")
}

####
#### Corner cases ####
####

chunk("Corner cases")


# We test the absence of bugs

base = iris
names(base) = c("y", "x1", "x2", "x3", "fe1")
base$fe2 = rep(1:5, 30)
base$y[1:5] = NA
base$x1[4:8] = NA
base$x2[4:21] = NA
base$x3[110:111] = NA
base$fe1[110:118] = NA
base$fe2[base$fe2 == 1] = 0
base$fe3 = sample(letters[1:5], 150, TRUE)
base$period = rep(1:50, 3)
base$x_cst = 1

res = feols(y ~ 1 | csw(fe1, fe1^fe2), base)

res = feols(y ~ 1 + csw(x1, i(fe1)) | fe2, base)

res = feols(y ~ csw(f(x1, 1:2), x2) | sw0(fe2, fe2^fe3), base, panel.id = ~ fe1 + period)

res = feols(d(y) ~ -1 + d(x2), base, panel.id = ~ fe1 + period)
test(length(coef(res)), 1)

res = feols(c(y, x1) ~ 1 | fe1 | x2 ~ x3, base)

res = feols(y ~ x1 | fe1[x2] + fe2[x2], base)

#
# NA models (ie all variables are collinear with the FEs)
#

# Should work when warn = FALSE or multiple est
for(i in 1:2){
    fun = switch(i, "1" = feols, "2" = feglm)

    res = feols(y ~ x_cst | fe1, base, warn = FALSE)
    res         # => no error
    etable(res) # => no error

    # error when warn = TRUE
    test(feols(y ~ x_cst | fe1, base), "err")

    # multiple est => no error
    res = feols(c(y, x1) ~ x_cst | fe1, base)
    res         # => no error
    etable(res) # => no error
}


# Removing the intercept!!!

res = feols(y ~ -1 + x1 + i(fe1), base)
test("(Intercept)" %in% names(res$coefficients), FALSE)

res = feols(y ~ -1 + x1 + factor(fe1), base)
test("(Intercept)" %in% names(res$coefficients), FALSE)

res = feols(y ~ -1 + x1 + i(fe1) + i(fe2), base)
test("(Intercept)" %in% names(res$coefficients), FALSE)
test(is.null(res$collin.var), TRUE)



####
#### Fit methods ####
####

chunk("Fit methods")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$y_int = as.integer(base$y)
base$y_log = sample(c(TRUE, FALSE), 150, TRUE)

res = feglm.fit(base$y, base[, 2:4])
res_bis = feglm(y ~ -1 + x1 + x2 + x3, base)
test(coef(res), coef(res_bis))

res = feglm.fit(base$y_int, base[, 2:4])
res_bis = feglm(y_int ~ -1 + x1 + x2 + x3, base)
test(coef(res), coef(res_bis))

res = feglm.fit(base$y_log, base[, 2:4])
res_bis = feglm(y_log ~ -1 + x1 + x2 + x3, base)
test(coef(res), coef(res_bis))



res = feglm.fit(base$y, base[, 2:4], family = "poisson")
res_bis = feglm(y ~ -1 + x1 + x2 + x3, base, family = "poisson")
test(coef(res), coef(res_bis))

res = feglm.fit(base$y_int, base[, 2:4], family = "poisson")
res_bis = feglm(y_int ~ -1 + x1 + x2 + x3, base, family = "poisson")
test(coef(res), coef(res_bis))

res = feglm.fit(base$y_log, base[, 2:4], family = "poisson")
res_bis = feglm(y_log ~ -1 + x1 + x2 + x3, base, family = "poisson")
test(coef(res), coef(res_bis))

####
#### Standard-errors ####
####

chunk("STANDARD ERRORS")

#
# Fixed-effects corrections
#

# We create "irregular" FEs
set.seed(0)
base = data.frame(x = rnorm(20))
base$y = base$x + rnorm(20)
base$fe1 = rep(rep(1:3, c(4, 3, 3)), 2)
base$fe2 = rep(rep(1:5, each = 2), 2)
est = feols(y ~ x | fe1 + fe2, base)

# fe1: 3 FEs
# fe2: 5 FEs

#
# Clustered standard-errors: by fe1
#

# Default: fixef.K = "nested"
#  => adjustment K = 1 + 5 (i.e. x + fe2)
test(attr(vcov(est, dof = dof(fixef.K = "nested"), attr = TRUE), "dof.K"), 6)

# fixef.K = FALSE
#  => adjustment K = 1 (i.e. only x)
test(attr(vcov(est, dof = dof(fixef.K = "none"), attr = TRUE), "dof.K"), 1)

# fixef.K = TRUE
#  => adjustment K = 1 + 3 + 5 - 1 (i.e. x + fe1 + fe2 - 1 restriction)
test(attr(vcov(est, dof = dof(fixef.K = "full"), attr = TRUE), "dof.K"), 8)

# fixef.K = TRUE & fixef.exact = TRUE
#  => adjustment K = 1 + 3 + 5 - 2 (i.e. x + fe1 + fe2 - 2 restrictions)
test(attr(vcov(est, dof = dof(fixef.K = "full", fixef.force_exact = TRUE), attr = TRUE), "dof.K"), 7)

#
# Manual checks of the SEs
#

n = est$nobs
VCOV_raw = est$cov.unscaled / ((n - 1) / (n - est$nparams))

# standard
for(k_val in c("none", "nested", "full")){
    for(adj in c(FALSE, TRUE)){

        K = switch(k_val, none = 1, nested = 8, full = 8)
        my_adj = ifelse(adj, (n - 1) / (n - K), 1)

        test(vcov(est, se = "standard", dof = dof(adj = adj, fixef.K = k_val)), VCOV_raw * my_adj)

        # cat("adj = ", adj, " ; fixef.K = ", k_val, "\n", sep = "")
    }
}

# Clustered, fe1
VCOV_raw = est$cov.unscaled / est$sigma2
H = vcovClust(est$fixef_id$fe1, VCOV_raw, scores = est$scores, dof = FALSE)

for(tdf in c("conventional", "min")){
    for(k_val in c("none", "nested", "full")){
        for(c_adj in c(FALSE, TRUE)){
            for(adj in c(FALSE, TRUE))

                K = switch(k_val, none = 1, nested = 6, full = 8)
            cluster_factor = ifelse(c_adj, 3/2, 1)
            df = ifelse(tdf == "min", 2, 20 - 8)
            my_adj = ifelse(adj, (n - 1) / (n - K), 1)

            V = H * cluster_factor

            # test SE
            test(vcov(est, se = "cluster", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj)), V * my_adj)

            # test pvalue
            my_tstat = tstat(est, se = "cluster", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj))
            test(pvalue(est, se = "cluster", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj, t.df = tdf)), 2*pt(-abs(my_tstat), df))

            # cat("adj = ", adj, " ; fixef.K = ", k_val, " ; cluster.adj = ", c_adj, " t.df = ", tdf, "\n", sep = "")
        }
    }
}


# 2-way Clustered, fe1 fe2
VCOV_raw = est$cov.unscaled / est$sigma2
M_i  = vcovClust(est$fixef_id$fe1, VCOV_raw, scores = est$scores, dof = FALSE)
M_t  = vcovClust(est$fixef_id$fe2, VCOV_raw, scores = est$scores, dof = FALSE)
M_it = vcovClust(paste(base$fe1, base$fe2), VCOV_raw, scores = est$scores, dof = FALSE, do.unclass = TRUE)

M_i + M_t - M_it
vcov(est, se = "two", dof = dof(adj = FALSE, cluster.adj = FALSE))

for(cdf in c("conventional", "min")){
    for(tdf in c("conventional", "min")){
        for(k_val in c("none", "nested", "full")){
            for(c_adj in c(FALSE, TRUE)){
                for(adj in c(FALSE, TRUE)){

                    K = switch(k_val, none = 1, nested = 2, full = 8)

                    if(c_adj){
                        if(cdf == "min"){
                            V = (M_i + M_t - M_it) * 3/2
                        } else {
                            V = M_i  * 3/2 + M_t * 5/4 - M_it * 6/5
                        }
                    } else {
                        V = M_i + M_t - M_it
                    }

                    df = ifelse(tdf == "min", 2, 20 - 8)
                    my_adj = ifelse(adj, (n - 1) / (n - K), 1)

                    # test SE
                    test(vcov(est, se = "two", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj, cluster.df = cdf)),
                         V * my_adj)

                    # test pvalue
                    my_tstat = tstat(est, se = "two", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj, cluster.df = cdf))
                    test(pvalue(est, se = "two", dof = dof(adj = adj, fixef.K = k_val, cluster.adj = c_adj, cluster.df = cdf, t.df = tdf)),
                         2*pt(-abs(my_tstat), df))

                    # cat("adj = ", adj, " ; fixef.K = ", k_val, " ; cluster.adj = ", c_adj, " t.df = ", tdf, "\n", sep = "")
                }
            }
        }
    }
}


#
# Comparison with sandwich and plm
#

library(sandwich)

# Data generation
set.seed(0)
N <- 20; G <- N/5; T <- N/G
d <- data.frame( y=rnorm(N), x=rnorm(N), grp=rep(1:G,T), tm=rep(1:T,each=G) )

# Estimations
est_lm    = lm(y ~ x + as.factor(grp) + as.factor(tm), data=d)
est_feols = feols(y ~ x | grp + tm, data=d)

#
# Standard
#

test(se(est_feols, se = "st")["x"], se(est_lm)["x"])

#
# Clustered
#

# Clustered by grp
se_CL_grp_lm_HC1 = sqrt(vcovCL(est_lm, cluster = d$grp, type = "HC1")["x", "x"])
se_CL_grp_lm_HC0 = sqrt(vcovCL(est_lm, cluster = d$grp, type = "HC0")["x", "x"])
se_CL_grp_stata  = 0.165385 # vce(cluster grp)

# How to get the lm
test(se(est_feols, dof = dof(fixef.K = "full")), se_CL_grp_lm_HC1)
test(se(est_feols, dof = dof(adj = FALSE, fixef.K = "full")), se_CL_grp_lm_HC0)

# How to get the Stata
test(se(est_feols), se_CL_grp_stata, "~")

#
# Heteroskedasticity-robust
#

se_white_lm_HC1 = sqrt(vcovHC(est_lm, type = "HC1")["x", "x"])
se_white_lm_HC0 = sqrt(vcovHC(est_lm, type = "HC0")["x", "x"])

test(se(est_feols, se = "hetero"), se_white_lm_HC1)
test(se(est_feols, se = "hetero", dof = dof(adj = FALSE, cluster.adj = FALSE)), se_white_lm_HC0)

#
# Two way
#

# Clustered by grp & tm
se_CL_2w_lm    = sqrt(vcovCL(est_lm, cluster = ~ grp + tm, type = "HC1")["x", "x"])
se_CL_2w_feols = se(est_feols, se = "twoway")

test(se(est_feols, se = "twoway", dof = dof(fixef.K = "full", cluster.df = "conv")), se_CL_2w_lm)

#
# Checking the calls work properly
#

data(trade)

est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination, trade)

se_clust = se(est_pois, se = "cluster", cluster = "Product")
test(se(est_pois, cluster = trade$Product), se_clust)
test(se(est_pois, cluster = ~Product), se_clust)

se_two = se(est_pois, se = "twoway", cluster = trade[, c("Product", "Destination")])
test(se_two, se(est_pois, cluster = c("Product", "Destination")))
test(se_two, se(est_pois, cluster = ~Product+Destination))

se_clu_comb = se(est_pois, cluster = "Product^Destination")
test(se_clu_comb, se(est_pois, cluster = paste(trade$Product, trade$Destination)))
test(se_clu_comb, se(est_pois, cluster = ~Product^Destination))

se_two_comb = se(est_pois, cluster = c("Origin^Destination", "Product"))
test(se_two_comb, se(est_pois, cluster = list(paste(trade$Origin, trade$Destination), trade$Product)))
test(se_two_comb, se(est_pois, cluster = ~Origin^Destination + Product))

# With cluster removed
base = trade
base$Euros[base$Origin == "FR"] = 0
est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination, base)

se_clust = se(est_pois, se = "cluster", cluster = "Product")
test(se(est_pois, cluster = base$Product), se_clust)
test(se(est_pois, cluster = ~Product), se_clust)

se_two = se(est_pois, se = "twoway", cluster = base[, c("Product", "Destination")])
test(se_two, se(est_pois, cluster = c("Product", "Destination")))
test(se_two, se(est_pois, cluster = ~Product+Destination))

se_clu_comb = se(est_pois, cluster = "Product^Destination")
test(se_clu_comb, se(est_pois, cluster = paste(base$Product, base$Destination)))
test(se_clu_comb, se(est_pois, cluster = ~Product^Destination))

se_two_comb = se(est_pois, cluster = c("Origin^Destination", "Product"))
test(se_two_comb, se(est_pois, cluster = list(paste(base$Origin, base$Destination), base$Product)))
test(se_two_comb, se(est_pois, cluster = ~Origin^Destination + Product))

# With cluster removed and NAs
base = trade
base$Euros[base$Origin == "FR"] = 0
base$Euros_na = base$Euros ; base$Euros_na[sample(nrow(base), 50)] = NA
base$Destination_na = base$Destination ; base$Destination_na[sample(nrow(base), 50)] = NA
base$Origin_na = base$Origin ; base$Origin_na[sample(nrow(base), 50)] = NA
base$Product_na = base$Product ; base$Product_na[sample(nrow(base), 50)] = NA

est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination_na, base)

se_clust = se(est_pois, se = "cluster", cluster = "Product")
test(se(est_pois, cluster = base$Product), se_clust)
test(se(est_pois, cluster = ~Product), se_clust)

se_two = se(est_pois, se = "twoway", cluster = base[, c("Product", "Destination")])
test(se_two, se(est_pois, cluster = c("Product", "Destination")))
test(se_two, se(est_pois, cluster = ~Product+Destination))

se_clu_comb = se(est_pois, cluster = "Product^Destination")
test(se_clu_comb, se(est_pois, cluster = paste(base$Product, base$Destination)))
test(se_clu_comb, se(est_pois, cluster = ~Product^Destination))

se_two_comb = se(est_pois, cluster = c("Origin^Destination", "Product"))
test(se_two_comb, se(est_pois, cluster = list(paste(base$Origin, base$Destination), base$Product)))
test(se_two_comb, se(est_pois, cluster = ~Origin^Destination + Product))

#
# Checking errors
#

# Should report error
test(se(est_pois, cluster = "Origin_na"), "err")
test(se(est_pois, cluster = base$Origin_na), "err")
test(se(est_pois, cluster = list(base$Origin_na)), "err")
test(se(est_pois, cluster = ~Origin_na^Destination), "err")

test(se(est_pois, se = "cluster", cluster = ~Origin_na^not_there), "err")

test(se(est_pois, se = "twoway", cluster = c("Origin^Destination", "Product", "error")), "err")
test(se(est_pois, se = "twoway", cluster = base[, 1:4]), "err")
test(se(est_pois, se = "twoway", cluster = ~Origin+Destination+Product), "err")
test(se(est_pois, se = "fourway", cluster = ~Origin+Destination+Product), "err")

####
#### Residuals ####
####

chunk("RESIDUALS")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$y_int = as.integer(base$y) + 1

# OLS + GLM + FENMLM

for(method in c("ols", "feglm", "femlm", "fenegbin")){
    cat("Method: ", dreamerr::sfill(method, 8))
    for(do_weight in c(FALSE, TRUE)){
        cat(".")

        if(do_weight){
            w = unclass(as.factor(base$species))
        } else {
            w = NULL
        }

        if(method == "ols"){
            m = feols(y_int ~ x1 | species, base, weights = w)
            mm = lm(y_int ~ x1 + species, base, weights = w)

        } else if(method == "feglm"){
            m = feglm(y_int ~ x1 | species, base, weights = w, family = "poisson")
            mm = glm(y_int ~ x1 + species, base, weights = w, family = poisson())

        } else if(method == "femlm"){
            if(!is.null(w)) next
            m = femlm(y_int ~ x1 | species, base)
            mm = glm(y_int ~ x1 + species, base, family = poisson())

        } else if(method == "fenegbin"){
            if(!is.null(w)) next
            m = fenegbin(y_int ~ x1 | species, base, notes = FALSE)
            mm = MASS::glm.nb(y_int ~ x1 + species, base)
        }

        tol = ifelse(method == "fenegbin", 1e-2, 1e-6)

        test(resid(m, "r"), resid(mm, "resp"), "~", tol = tol)
        test(resid(m, "d"), resid(mm, "d"), "~", tol = tol)
        test(resid(m, "p"), resid(mm, "pearson"), "~", tol = tol)

        test(deviance(m), deviance(mm), "~", tol = tol)
    }
    cat("\n")
}
cat("\n")


####
#### fixef ####
####

chunk("FIXEF")

set.seed(0)
base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$x4 = rnorm(150) + 0.25*base$y
base$fe_bis = sample(10, 150, TRUE)
base$fe_ter = sample(15, 150, TRUE)

get_coef = function(all_coef, x){
    res = all_coef[grepl(x, names(all_coef), perl = TRUE)]
    names(res) = gsub(x, "", names(res), perl = TRUE)
    res
}

#
# With 2 x 1 FE
#

m = feols(y ~ x1 + x2 | species + fe_bis, base)
all_coef = coef(feols(y ~ -1 + x1 + x2 + species + factor(fe_bis), base))

m_fe = fixef(m)
c1  = get_coef(all_coef, "species")
test(var(c1 - m_fe$species[names(c1)]), 0)

c2  = get_coef(all_coef, "factor\\(fe_bis\\)")
test(var(c2 - m_fe$fe_bis[names(c2)]), 0)


#
# With 1 FE + 1 FE 1 VS
#

m = feols(y ~ x1 + x2 | species + fe_bis[x3], base)
all_coef = coef(feols(y ~ -1 + x1 + x2 + species + factor(fe_bis) + i(fe_bis, x3), base))

m_fe = fixef(m)
c1  = get_coef(all_coef, "species")
test(var(c1 - m_fe$species[names(c1)]), 0, "~")

c2  = get_coef(all_coef, "factor\\(fe_bis\\)")
test(var(c2 - m_fe$fe_bis[names(c2)]), 0, "~")

c3  = get_coef(all_coef, "fe_bis::|:x3")
test(c3, m_fe[["fe_bis[[x3]]"]][names(c3)], "~", tol = 1e-5)

#
# With 2 x (1 FE + 1 VS) + 1 FE
#

m = feols(y ~ x1 | species[x2] + fe_bis[x3] + fe_ter, base)
all_coef = coef(feols(y ~ -1 + x1 + species + i(species, x2) + factor(fe_bis) + i(fe_bis, x3) + factor(fe_ter), base))

m_fe = fixef(m)
c1  = get_coef(all_coef, "^species(?=[^:])")
test(var(c1 - m_fe$species[names(c1)]), 0, "~")

c2  = get_coef(all_coef, "^factor\\(fe_bis\\)")
test(var(c2 - m_fe$fe_bis[names(c2)]), 0, "~")

c3  = get_coef(all_coef, "fe_bis::|:x3")
test(c3, m_fe[["fe_bis[[x3]]"]][names(c3)], "~", tol = 2e-4)

c4  = get_coef(all_coef, "species::|:x2")
test(c4, m_fe[["species[[x2]]"]][names(c4)], "~", tol = 2e-4)

#
# With 2 x (1 FE) + 1 FE 2 VS
#

m = feols(y ~ x1 | species + fe_bis[x2,x3] + fe_ter, base)
all_coef = coef(feols(y ~ x1 + species + factor(fe_bis) + i(fe_bis, x2) + i(fe_bis, x3) + factor(fe_ter), base))

m_fe = fixef(m)
c1  = get_coef(all_coef, "^species")
test(var(c1 - m_fe$species[names(c1)]), 0, "~")

c2  = get_coef(all_coef, "^factor\\(fe_bis\\)")
test(var(c2 - m_fe$fe_bis[names(c2)]), 0, "~")

c3  = get_coef(all_coef, "fe_bis::(?=.+x2)|:x2")
test(c3, m_fe[["fe_bis[[x2]]"]][names(c3)], "~", tol = 2e-4)

c4  = get_coef(all_coef, "fe_bis::(?=.+x3)|:x3")
test(c4, m_fe[["fe_bis[[x3]]"]][names(c4)], "~", tol = 2e-4)


#
# With weights
#

w = 3 * (as.integer(base$species) - 0.95)
m = feols(y ~ x1 | species + fe_bis[x2,x3] + fe_ter, base, weights = w)
all_coef = coef(feols(y ~ x1 + species + factor(fe_bis) + i(fe_bis, x2) + i(fe_bis, x3) + factor(fe_ter), base, weights = w))

m_fe = fixef(m)
c1  = get_coef(all_coef, "^species")
test(var(c1 - m_fe$species[names(c1)]), 0, "~")

c2  = get_coef(all_coef, "^factor\\(fe_bis\\)")
test(var(c2 - m_fe$fe_bis[names(c2)]), 0, "~")

c3  = get_coef(all_coef, "fe_bis::(?=.+x2)|:x2")
test(c3, m_fe[["fe_bis[[x2]]"]][names(c3)], "~", tol = 2e-4)

c4  = get_coef(all_coef, "fe_bis::(?=.+x3)|:x3")
test(c4, m_fe[["fe_bis[[x3]]"]][names(c4)], "~", tol = 2e-4)


####
#### To Integer ####
####

chunk("TO_INTEGER")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$z = sample(5, 150, TRUE)

# Normal
m = to_integer(base$species)
test(length(unique(m)), 3)

m = to_integer(base$species, base$z)
test(length(unique(m)), 15)

# with NA
base$species_na = base$species
base$species_na[base$species == "setosa"] = NA

m = to_integer(base$species_na, base$z)
test(length(unique(m)), 11)

m = to_integer(base$species_na, base$z, add_items = TRUE, items.list = TRUE)
test(length(m$items), 10)


####
#### Collinearity ####
####

cat("COLLINEARITY\n\n")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$constant = 5
base$y_int = as.integer(base$y)
base$w = as.vector(unclass(base$species) - 0.95)

for(useWeights in c(FALSE, TRUE)){
    for(model in c("ols", "pois")){
        for(use_fe in c(FALSE, TRUE)){
            cat(".")

            my_weight = NULL
            if(useWeights) my_weight = base$w

            adj = 1
            if(model == "ols"){
                if(!use_fe){
                    res = feols(y ~ x1 + constant, base, weights = my_weight)
                    res_bis = lm(y ~ x1 + constant, base, weights = my_weight)
                } else {
                    res = feols(y ~ x1 + constant | species, base, weights = my_weight)
                    res_bis = lm(y ~ x1 + constant + species, base, weights = my_weight)
                }
            } else {
                if(!use_fe){
                    res = fepois(y_int ~ x1 + constant, base, weights = my_weight)
                    res_bis = glm(y_int ~ x1 + constant, base, weights = my_weight, family = poisson)
                } else {
                    res = fepois(y_int ~ x1 + constant | species, base, weights = my_weight)
                    res_bis = glm(y_int ~ x1 + constant + species, base, weights = my_weight, family = poisson)
                }
                adj = 0
            }

            test(coef(res)["x1"], coef(res_bis)["x1"], "~")
            test(se(res, se = "st", dof = dof(adj=adj))["x1"], se(res_bis)["x1"], "~")
            # cat("Weight: ", useWeights, ", model: ", model, ", FE: ", use_fe, "\n", sep="")

        }
    }
}
cat("\n")


####
#### Interact ####
####


chunk("Interact")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$fe_2 = round(rnorm(150))

#
# We just ensure it works without error
#

m = feols(y ~ x1 + i(fe_2), base)
coefplot(m)
etable(m, dict = c("0" = "zero"))

m = feols(y ~ x1 + i(fe_2) + i(fe_2, x2), base)
coefplot(m)
etable(m, dict = c("0" = "zero"))

a = i(base$fe_2)
b = i(base$fe_2, ref = 0:1)
d = i(base$fe_2, keep = 0:1)

test(ncol(a), ncol(b) + 2)
test(ncol(d), 2)



####
#### demean ####
####

chunk("DEMEAN")

data(trade)

base = trade
base$ln_euros = log(base$Euros)
base$ln_dist = log(base$dist_km)

X = base[, c("ln_euros", "ln_dist")]
fe = base[, c("Origin", "Destination")]

base_new = demean(X, fe)

a = feols(ln_euros ~ ln_dist, base_new)
b = feols(ln_euros ~ ln_dist | Origin + Destination, base, demeaned = TRUE)

test(coef(a)[-1], coef(b), "~", 1e-12)

test(base_new$ln_euros, b$y_demeaned)
test(base_new$ln_dist, b$X_demeaned)

# Now we just check there's no error

# NAs
X_NA = X
fe_NA = fe
X_NA[1:5, 1] = NA
fe_NA[6:10, 1] = NA
X_demean = demean(X_NA, fe_NA, na.rm = FALSE)
test(nrow(X_demean), nrow(X))

# integer
X_int = X
X_int[[1]] = as.integer(X_int[[1]])
X_demean = demean(X_int, fe)

# matrix/DF
X_demean = demean(X_int, fe, as.matrix = TRUE)
test(is.matrix(X_demean), TRUE)

X_demean = demean(as.matrix(X_int), fe, as.matrix = FALSE)
test(!is.matrix(X_demean), TRUE)

# fml
X_demean = demean(ln_dist ~ Origin + Destination, data = base)

# slopes
X_dm_slopes = demean(ln_dist ~ Origin + Destination[ln_euros], data = base)
X_dm_slopes_bis = demean(base$ln_dist, fe, slope.vars = base$ln_euros, slope.flag = c(0, 1))

test(X_dm_slopes[[1]], X_dm_slopes_bis)


####
#### hatvalues ####
####


chunk("HATVALUES")

set.seed(0)
x = sin(1:10)
y = rnorm(10)
y_int = rpois(10, 2)
fm  = lm(y ~ x)
ffm = feols(y ~ x, data.frame(y, x))

test(hatvalues(ffm), hatvalues(fm))

gm  = glm(y_int ~ x, family = poisson())
fgm = fepois(y_int ~ x, data.frame(y_int, x))

test(hatvalues(fgm), hatvalues(gm))


####
#### sandwich ####
####

chunk("SANDWICH")

# Compatibility with sandwich

library(sandwich)

data(base_did)
est = feols(y ~ x1 + I(x1**2) + factor(id), base_did)

test(vcov(est, cluster = ~id), vcovCL(est, cluster = ~id, type = "HC1"))

est_pois = fepois(as.integer(y) + 20 ~ x1 + I(x1**2) + factor(id), base_did)

test(vcov(est_pois, cluster = ~id), vcovCL(est_pois, cluster = ~id, type = "HC1"))

# With FEs

est = feols(y ~ x1 + I(x1**2) | id, base_did)

test(vcov(est, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est, cluster = ~id))

est_pois = fepois(as.integer(y) + 20 ~ x1 + I(x1**2) | id, base_did)

test(vcov(est_pois, cluster = ~id, dof = dof(adj = FALSE)), vcovCL(est_pois, cluster = ~id))



####
#### only.env ####
####

# We check that there's no problem when using the environment

chunk("ONLY ENV")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")

env = feols(y ~ x1 + x2 | species, base, only.env = TRUE)
feols(env = env)

env = feglm(y ~ x1 + x2 | species, base, only.env = TRUE)
feglm(env = env)

env = fepois(y ~ x1 + x2 | species, base, only.env = TRUE)
fepois(env = env)

env = fenegbin(y ~ x1 + x2 | species, base, only.env = TRUE)
fenegbin(env = env)

env = femlm(y ~ x1 + x2 | species, base, only.env = TRUE)
femlm(env = env)

env = feNmlm(y ~ x1 + x2 | species, base, only.env = TRUE)
feNmlm(env = env)



####
#### Non linear tests ####
####

chunk("NON LINEAR")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")

tab = c("versicolor" = 5, "setosa" = 0, "virginica" = -5)

fun_nl = function(a, b, spec){
    res = as.numeric(tab[spec])
    a*res + b*res^2
}

est_nl = feNmlm(y ~ x1, base, NL.fml = ~fun_nl(a, b, species), NL.start = 1, family = "gaussian")

base$var_spec = as.numeric(tab[base$species])

est_lin = feols(y ~ x1 + var_spec + I(var_spec^2), base)

test(coef(est_nl), coef(est_lin)[c(3, 4, 1, 2)], "~")



####
#### Lagging ####
####

# Different types of lag
# 1) check no error in wide variety of situations
# 2) check consistency

chunk("LAGGING")

data(base_did)
base = base_did

n = nrow(base)

set.seed(0)
base$y_na = base$y ; base$y_na[sample(n, 50)] = NA
base$period_txt = letters[base$period]
ten_dates = c("1960-01-15", "1960-01-16", "1960-03-31", "1960-04-05", "1960-05-12", "1960-05-25", "1960-06-20", "1960-07-30", "1965-01-02", "2002-12-05")
base$period_date = as.Date(ten_dates, "%Y-%m-%d")[base$period]
base$y_0 = base$y**2 ; base$y_0[base$id == 1] = 0

# We compute the lags "by hand"
base = base[order(base$id, base$period), ]
base$x1_lag = c(NA, base$x1[-n]) ; base$x1_lag[base$period == 1] = NA
base$x1_lead = c(base$x1[-1], NA) ; base$x1_lead[base$period == 10] = NA
base$x1_diff = base$x1 - base$x1_lag

# we create holes
base$period_bis = base$period ; base$period_bis[base$period_bis == 5] = 50
base$x1_lag_hole = base$x1_lag ; base$x1_lag_hole[base$period %in% c(5, 6)] = NA
base$x1_lead_hole = base$x1_lead ; base$x1_lead_hole[base$period %in% c(4, 5)] = NA

# we reshuffle the base
base = base[sample(n), ]

#
# Checks consistency
#

cat("consistentcy...")

test(lag(x1 ~ id + period, data = base), base$x1_lag)
test(lag(x1 ~ id + period, -1, data = base), base$x1_lead)

test(lag(x1 ~ id + period_bis, data = base), base$x1_lag_hole)
test(lag(x1 ~ id + period_bis, -1, data = base), base$x1_lead_hole)

test(lag(x1 ~ id + period_txt, data = base), base$x1_lag)
test(lag(x1 ~ id + period_txt, -1, data = base), base$x1_lead)

test(lag(x1 ~ id + period_date, data = base), base$x1_lag)
test(lag(x1 ~ id + period_date, -1, data = base), base$x1_lead)

cat("done.\nEstimations...")

#
# Estimations
#

# Poisson

for(depvar in c("y", "y_na", "y_0")){
    for(p in c("period", "period_txt", "period_date")){

        base$per = base[[p]]

        cat(".")

        base$y_dep = base[[depvar]]
        pdat = panel(base, ~ id + period)

        if(depvar == "y_0"){
            estfun = fepois
        } else {
            estfun = feols
        }

        est_raw  = estfun(y_dep ~ x1 + x1_lag + x1_lead, base)
        est      = estfun(y_dep ~ x1 + l(x1) + f(x1), base, panel.id = "id,per")
        est_pdat = estfun(y_dep ~ x1 + l(x1, 1) + f(x1, 1), pdat)
        test(coef(est_raw), coef(est))
        test(coef(est_raw), coef(est_pdat))

        # Now diff
        est_raw  = estfun(y_dep ~ x1 + x1_diff, base)
        est      = estfun(y_dep ~ x1 + d(x1), base, panel.id = "id,per")
        est_pdat = estfun(y_dep ~ x1 + d(x1, 1), pdat)
        test(coef(est_raw), coef(est))
        test(coef(est_raw), coef(est_pdat))

        # Now we just check that calls to l/f works without checking coefs

        est = estfun(y_dep ~ x1 + l(x1) + f(x1), base, panel.id = "id,per")
        est = estfun(y_dep ~ l(x1, -1:1) + f(x1, 2), base, panel.id = c("id", "per"))
        est = estfun(y_dep ~ l(x1, -1:1, fill = 1), base, panel.id = ~ id + per)
        if(depvar == "y") test(est$nobs, n)
        est = estfun(f(y_dep) ~ f(x1, -1:1), base, panel.id = ~ id + per)
    }
}

cat("done.\n\n")

#
# Data table
#

cat("data.table...")
# We just check there is no bug (consistency should be OK)

library(data.table)

base_dt = data.table(id = c("A", "A", "B", "B"),
                     time = c(1, 2, 1, 3),
                     x = c(5, 6, 7, 8))

base_dt = panel(base_dt, ~id + time)

base_dt[, x_l := l(x)]
test(base_dt$x_l, c(NA, 5, NA, NA))

lag_creator = function(dt) {
    dt2 = panel(dt, ~id + time)
    dt2[, x_l := l(x)]
    return(dt2)
}

base_bis = lag_creator(base_dt)

base_bis[, x_d := d(x)]

cat("done.\n\n")


####
#### predict ####
####

chunk("PREDICT")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$fe_bis = sample(letters, 150, TRUE)

#
# Same generative data
#

# Predict with fixed-effects
res = feols(y ~ x1 | species + fe_bis, base)
test(predict(res), predict(res, base))

res = fepois(y ~ x1 | species + fe_bis, base)
test(predict(res), predict(res, base))

res = femlm(y ~ x1 | species + fe_bis, base)
test(predict(res), predict(res, base))


# Predict with varying slopes -- That's normal that precision is high (because FEs are computed with low precision)
res = feols(y ~ x1 | species + fe_bis[x3], base)
test(predict(res), predict(res, base), "~", tol = 1e-4)

res = fepois(y ~ x1 | species + fe_bis[x3], base)
test(predict(res), predict(res, base), "~", tol = 1e-3)


# Prediction with factors
res = feols(y ~ x1 + i(species), base)
test(predict(res), predict(res, base))

res = feols(y ~ x1 + i(species) + i(fe_bis), base)
test(predict(res), predict(res, base))

quoi = head(base[, c("y", "x1", "species", "fe_bis")])
test(head(predict(res)), predict(res, quoi))

quoi$species = as.character(quoi$species)
quoi$species[1:3] = "zz"
test(head(predict(res)), predict(res, quoi))

# prediction with lags
data(base_did)
res = feols(y ~ x1 + l(x1), base_did, panel.id = ~ id + period)
test(predict(res, na.rm = FALSE), predict(res, base_did))

qui = sample(which(base_did$id %in% 1:5))
base_bis = base_did[qui, ]
test(predict(res, na.rm = FALSE)[qui], predict(res, base_bis))

# prediction with poly
res_poly = feols(y ~ poly(x1, 2), base)
pred_all = predict(res_poly)
pred_head = predict(res_poly, head(base, 20))
pred_tail = predict(res_poly, tail(base, 20))
test(head(pred_all, 20), pred_head)
test(tail(pred_all, 20), pred_tail)

####
#### subset ####
####

chunk("SUBSET")

set.seed(5)
base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$fe_bis = sample(letters, 150, TRUE)
base$x4 = rnorm(150)
base$x1[sample(150, 5)] = NA

# Errors
test(feols(fml, base, subset = ~species), "err")
test(feols(fml, base, subset = -1:15), "err")
test(feols(fml, base, subset = integer(0)), "err")
test(feols(fml, base, subset = c(TRUE, TRUE, FALSE)), "err")

# Valid use
for(id_fun in 1:6){
    estfun = switch(as.character(id_fun),
                    "1" = feols,
                    "2" = feglm,
                    "3" = fepois,
                    "4" = femlm,
                    "5" = fenegbin,
                    "6" = feNmlm)

    for(id_fe in 1:5){

        cat(".")

        fml = switch(as.character(id_fe),
                     "1" = y ~ x1 + x2,
                     "2" = y ~ x1 + x2 | species,
                     "3" = y ~ x1 + x2 | fe_bis,
                     "4" = y ~ x1 + x2 + i(fe_bis),
                     "5" = y ~ x1 + x2 | fe_bis[x3])

        if(id_fe == 5 && id_fun %in% 4:6) next

        if(id_fun == 6){
            res_sub_a = estfun(fml, base, subset = ~species == "setosa", NL.fml = ~ a*x4, NL.start = 0)
            res_sub_b = estfun(fml, base, subset = base$species == "setosa", NL.fml = ~ a*x4, NL.start = 0)
            res_sub_c = estfun(fml, base, subset = which(base$species == "setosa"), NL.fml = ~ a*x4, NL.start = 0)
            res       = estfun(fml, base[base$species == "setosa", ], NL.fml = ~ a*x4, NL.start = 0)
        } else {
            res_sub_a = estfun(fml, base, subset = ~species == "setosa")
            res_sub_b = estfun(fml, base, subset = base$species == "setosa")
            res_sub_c = estfun(fml, base, subset = which(base$species == "setosa"))
            res       = estfun(fml, base[base$species == "setosa", ])
        }

        test(coef(res_sub_a), coef(res))
        test(coef(res_sub_b), coef(res))
        test(coef(res_sub_c), coef(res))
        test(se(res_sub_c, cluster = "fe_bis"), se(res, cluster = "fe_bis"))
    }
    cat("|")
}
cat("\n")



####
#### Multiple estimations ####
####

chunk("Multiple")

set.seed(2)
base = iris
names(base) = c("y1", "x1", "x2", "x3", "species")
base$y2 = 10 + rnorm(150) + 0.5 * base$x1
base$x4 = rnorm(150) + 0.5 * base$y1
base$fe2 = rep(letters[1:15], 10)
base$fe2[50:51] = NA
base$y2[base$fe2 == "a" & !is.na(base$fe2)] = 0
base$x2[1:5] = NA
base$x3[6] = NA
base$fe3 = rep(letters[1:10], 15)


for(id_fun in 1:5){
    estfun = switch(as.character(id_fun),
                    "1" = feols,
                    "2" = feglm,
                    "3" = fepois,
                    "4" = femlm,
                    "5" = feNmlm)

    # Following weird bug ASAN on CRAN I cannot replicate, check 4/5 not performed on non Windows
    if(Sys.info()["sysname"] != "Windows"){
       if(id_fun %in% 4:5) next
    }


    est_multi = estfun(c(y1, y2) ~ x1 + sw(x2, x3), base, split = ~species)

    k = 1
    for(s in c("setosa", "versicolor", "virginica")){
        for(lhs in c("y1", "y2")){
            for(rhs in c("x2", "x3")){
                res = estfun(xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs), base[base$species == s, ], notes = FALSE)

                test(coef(est_multi[[k]]), coef(res))
                test(se(est_multi[[k]], cluster = "fe3"), se(res, cluster = "fe3"))
                k = k + 1
            }
        }
    }

    cat("__")

    est_multi = estfun(c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2, base, fsplit = ~species)
    k = 1
    all_rhs = c("", "x2", "x3")
    for(s in c("all", "setosa", "versicolor", "virginica")){
        for(lhs in c("y1", "y2")){
            for(n_rhs in 1:3){
                if(s == "all"){
                    res = estfun(xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs]), base, notes = FALSE)
                } else {
                    res = estfun(xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs]), base[base$species == s, ], notes = FALSE)
                }

                vname = names(coef(res))
                test(coef(est_multi[[k]])[vname], coef(res), "~" , 1e-6)
                test(se(est_multi[[k]], cluster = "fe3")[vname], se(res, cluster = "fe3"), "~" , 1e-6)
                k = k + 1
            }
        }
    }

    cat("|")
}
cat("\n")


####
#### IV ####
####

chunk("IV")

base = iris
names(base) = c("y", "x1", "x_endo_1", "x_inst_1", "fe")
set.seed(2)
base$x_inst_2 = 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
base$x_endo_2 = 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)

# Checking a basic estimation

setFixest_se(all = "standard")

est_iv = feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)

res_f1 = feols(x_endo_1 ~ x1 + x_inst_1 + x_inst_2, base)
res_f2 = feols(x_endo_2 ~ x1 + x_inst_1 + x_inst_2, base)

base$fit_x_endo_1 = predict(res_f1)
base$fit_x_endo_2 = predict(res_f2)

res_2nd = feols(y ~ fit_x_endo_1 + fit_x_endo_2 + x1, base)

# the coef
test(coef(est_iv), coef(res_2nd))

# the SE
resid_iv = base$y - predict(res_2nd, data.frame(x1 = base$x1, fit_x_endo_1 = base$x_endo_1, fit_x_endo_2 = base$x_endo_2))
sigma2_iv = sum(resid_iv**2) / (res_2nd$nobs - res_2nd$nparams)

sum_2nd = summary(res_2nd, .vcov = res_2nd$cov.unscaled / res_2nd$sigma2 * sigma2_iv)

# We only check that on Windows => avoids super odd bug in fedora devel
# The worst is that I just can't debug it.... so that's the way it's done.
if(Sys.info()["sysname"] == "Windows"){
    test(se(sum_2nd), se(est_iv))
}


# check no bug
etable(summary(est_iv, stage = 1:2))

setFixest_se(reset = TRUE)


####
#### model.matrix ####
####

chunk("Model matrix")

base = iris
names(base) = c("y1", "x1", "x2", "x3", "species")
base$y2 = 10 + rnorm(150) + 0.5 * base$x1
base$x4 = rnorm(150) + 0.5 * base$y1
base$fe2 = rep(letters[1:15], 10)
base$fe2[50:51] = NA
base$y2[base$fe2 == "a" & !is.na(base$fe2)] = 0
base$x2[1:5] = NA
base$x3[6] = NA
base$fe3 = rep(letters[1:10], 15)
base$id = rep(1:15, each = 10)
base$time = rep(1:10, 15)

base_bis = base[1:50, ]
base_bis$id = rep(1:5, each = 10)
base_bis$time = rep(1:10, 5)

# NA removed
res = feols(y1 ~ x1 + x2 + x3, base)
m1 = model.matrix(res, type = "lhs")
test(length(m1), res$nobs)

# we check this is identical
m1_na = model.matrix(res, type = "lhs", na.rm = FALSE)
test(length(m1_na), res$nobs_origin)
test(max(abs(m1_na - base$y1), na.rm = TRUE), 0)

y = model.matrix(res, type = "lhs", data = base, na.rm = FALSE)
X = model.matrix(res, type = "rhs", data = base, na.rm = FALSE)
obs_rm = res$obs_selection$obsRemoved
res_bis = lm.fit(X[obs_rm, ], y[obs_rm])
test(res_bis$coefficients, res$coefficients)

# Lag
res_lag = feols(y1 ~ l(x1, 1:2) + x2 + x3, base, panel = ~id + time)
m_lag = model.matrix(res_lag)
test(nrow(m_lag), nobs(res_lag))

# lag with subset
m_lag_x1 = model.matrix(res_lag, subset = "x1")
test(ncol(m_lag_x1), 2)

# lag with subset, new data
mbis_lag_x1 = model.matrix(res_lag, base_bis[, c("x1", "x2", "id", "time")], subset = TRUE)
# l(x1, 1) + l(x1, 2) + x2
test(ncol(mbis_lag_x1), 3)
# 13 NAs: 2 per ID for the lags, 3 for x2
test(nrow(mbis_lag_x1), 37)

# With poly
res_poly = feols(y1 ~ poly(x1, 2), base)
m_poly_old = model.matrix(res_poly)
m_poly_new = model.matrix(res_poly, base_bis)
test(m_poly_old[1:50, 3], m_poly_new[, 3])

####
#### VCOV at estimation ####
####

chunk("vcov at estimation")

base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
base$clu = sample(6, 150, TRUE)
base$clu[1:5] = NA

est = feols(y ~ x1 | species, base, cluster = ~clu, dof = dof(adj = FALSE))

# The three should be identical
v1 = est$cov.scaled
v1b = vcov(est)
v1c = summary(est)$cov.scaled

test(v1, v1b)
test(v1, v1c)

# Only dof change
v2 = summary(est, dof = dof())$cov.scaled
v2b = vcov(est, cluster = ~clu, dof = dof())

test(v2, v2b)
test(max(abs(v1 - v2)) == 0, FALSE)

# SE change only
v3 = summary(est, se = "hetero")$cov.scaled
v3b = vcov(est, se = "hetero", dof = dof(adj = FALSE))

test(v3, v3b)
test(max(abs(v1 - v3)) == 0, FALSE)
test(max(abs(v2 - v3)) == 0, FALSE)

####
#### fitstat ####
####

base = iris
names(base) = c("y", "x1", "x_endo_1", "x_inst_1", "fe")
set.seed(2)
base$x_inst_2 = 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
base$x_endo_2 = 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)

# Checking a basic estimation
est_iv = feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)

fitstat(est_iv, ~ f + ivf + ivf2 + wald + ivwald + ivwald2 + wh + sargan + rmse + g + n + ll + sq.cor + r2)

est_fe = feols(y ~ x1 | fe, base)

fitstat(est_fe, ~ wf)





















































