
base <- datab14()
Est <- list()
mods <- c("ols", "glm", "femlm", "feNmlm")
for (k in 1:length(mods)) {
  Est[[k]] <- fixest_mod_select(model = mods[k],
                                fmla = c(y1, y2) ~ x1 + sw(x2, x3),
                                base = base,
                                split = ~species,
                                famly = "poisson")
}
species_indx <- list(
  "setosa" = which(base$species == "setosa"),
  "versicolor" = which(base$species == "versicolor"),
  "virginica" = which(base$species == "virginica")
)

### Test that multiple estimation works properly for feols, feglm, femlm and feNmlm.
for (k in 1:length(mods)) { # For each estimation method (feols, feglm, fmlm, feNmlm)
  est_multi <- Est[[k]]

  patrick::with_parameters_test_that("multiple estimations works properly for feols, feglm, femlm and feNmlm",
                                     {
                                       est <- est_multi[[num_fmla]]

                                       base_aux <- base[species_indx[[s]], ]
                                       fmlas = xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs)
                                       res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)

                                       expect_equal(coef(est), coef(res))
                                       expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
                                     },
                                     .cases = multiple_cases(mods[k])
  )
}


# k = 1
# est_multi <- Est[[k]]
#
# K = 1
# casos = multiple_cases(mods[k])[K,]
#
# num_fmla =casos$num_fmla
# s = casos$s
# lhs = casos$lhs
# rhs = casos$rhs
# method = casos$method
# fmly = casos$fmly
#
# est <- est_multi[[ num_fmla ]]
#
# base_aux <- base[species_indx[[s]], ]
# fmlas = xpd(..lhs ~ x1 + ..rhs, ..lhs = lhs, ..rhs = rhs)
# res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
#
# expect_equal(coef(est), coef(res))
# expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))



### Test that multiple estimation with fixed effects works properly for feols, feglm, femlm and feNmlm

## Est multiple
Est <- list()
mods <- c("ols", "glm", "femlm", "feNmlm")
for (k in 1:length(mods)) {
  Est[[k]] <- fixest_mod_select(model = mods[k],
                                fmla = c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2,
                                base = base,
                                fsplit = ~species,
                                famly = "poisson")
}
all_rhs <- c("", "x2", "x3")
for (k in 1:length(mods)) {
  est_multi <- Est[[k]]
  patrick::with_parameters_test_that("multiple estimations with fixed effects works properly for feols, feglm, femlm and feNmlm",
            {
              fmlas = xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs])
              if(s == "all"){
                res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base, famly = fmly, notes = FALSE)
              } else {
                base_aux <- base[species_indx[[s]], ]
                res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)
              }

              est <- est_multi[[num_fmla]]

              vname <- names(coef(res))
              expect_equal2(coef(est_multi[[ num_fmla ]])[vname], coef(res), tolerance =  1e-6)
              expect_equal2(se(est_multi[[ num_fmla ]], cluster = "fe3")[vname], se(res, cluster = "fe3"), tolerance = 1e-6)
            },
            .cases = multiple_cases2(mods[k])
  )
}

# k = 1
# est_multi <- Est[[k]]
#
# K = 3
# casos = multiple_cases2(mods[k])[K,]
#
# method = casos$method
# s = casos$s
# famly = casos$fmly
# num_fmla = casos$num_fmla
#
# lhs = casos$lhs
# n_rhs = casos$n_rhs
# s = casos$s
# fmlas = xpd(..lhs ~ x1 + ..rhs + x4 | species + fe2, ..lhs = lhs, ..rhs = all_rhs[1:n_rhs])
#
#
# if(s == "all"){
#   res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base, famly = famly, notes = FALSE)
# } else {
#   base_aux <- base[species_indx[[s]], ]
#   res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base_aux, famly = famly, notes = FALSE)
# }
#
#
# est <- est_multi[[num_fmla]]
#
# vname <- names(coef(res))
# expect_equal2(coef(est_multi[[K]])[vname], coef(res), tolerance =  1e-6)
# expect_equal2(se(est_multi[[K]], cluster = "fe3")[vname], se(res, cluster = "fe3"), tolerance = 1e-6)
