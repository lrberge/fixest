
base <- datab14()
Est <- list()
mods <- c("ols", "glm", "femlm", "feNmlm")
for (k in 1:length(mods)) {
  Est[[k]] <- fixest_mod_select(model = mods[k], fmla = c(y1, y2) ~ x1 + sw(x2, x3), base = base, split = ~species, famly = "poisson")
}
names(Est) <- mods

species_indx <- list(
  "setosa" = which(base$species == "setosa"),
  "versicolor" = which(base$species == "versicolor"),
  "virginica" = which(base$species == "virginica")
)

for (k in 1:length(mods)) {
  est_multi <- Est[[k]]

  patrick::with_parameters_test_that("aux",
    {
      est <- est_multi[[num_fmla]]

      base_aux <- base[species_indx[[s]], ]
      res <- fixest_mod_select(model = mods[k], fmla = fmlas, base = base_aux, famly = fmly, notes = FALSE)

      expect_equal(coef(est), coef(res))
      expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
    },
    .cases = multiple_cases(mods[k])
  )
}

# K = 20
# casos = multiple_cases()[K,]
#
# method = casos$method
# fmla = casos$fmlas
# s = casos$s
# famly = casos$fmly
# num_fmla = casos$num_fmla
#
# est_multi = Est[[method]]
# est <- est_multi[[num_fmla]]
#
# base_aux <- base[species_indx[[s]],]
# res <- fixest_mod_select(model = method, fmla = fmla, base = base_aux, famly = fmly, notes = FALSE)
#
# expect_equal(coef(est), coef(res))
# expect_equal(se(est, cluster = "fe3"), se(res, cluster = "fe3"))
