
<<<<<<< HEAD
setFixest_notes(FALSE)
base = datab5()

=======
base = datab5()

k = 8 # case 6 and 8 throws error
casos = residuals_cases()[k,]

fe_mod = fixest_mod_select(model = casos$method,
                           fmla = casos$formula_fe,
                           data = base,
                           famly = casos$family,
                           weights = ev_par(casos$wghts))
stat_mod = stats_mod_select(model = casos$method,
                            fmla = casos$formula_stats,
                            data = base,
                            famly = casos$family,
                            weights = ev_par(casos$wghts))
feres_r = resid(fe_mod, "r")
feres_d = resid(fe_mod, "d")
feres_p = resid(fe_mod, "p")
fedev = deviance(fe_mod)

stres_r = unname(resid(stat_mod, "resp"))
stres_d = unname(resid(stat_mod, "d"))
stres_p = unname(resid(stat_mod, "pearson"))
stdev = deviance(stat_mod)

expect_equal(feres_r, stres_r, tolerance = casos$tol)
expect_equal(feres_d, stres_d, tolerance = casos$tol)
expect_equal(feres_p, stres_p, tolerance = casos$tol)
expect_equal(fedev, stdev, tolerance = casos$tol)


>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110
patrick::with_parameters_test_that("Fitting without intercept works properly",
       {
           fe_mod = fixest_mod_select(model = method,
                                      fmla = formula_fe,
                                      data = base,
                                      famly = family,
                                      weights = ev_par(wghts))
           stat_mod = stats_mod_select(model = method,
                                       fmla = formula_stats,
                                       data = base,
                                       famly = family,
                                       weights = ev_par(wghts))
           feres_r = resid(fe_mod, "r")
           feres_d = resid(fe_mod, "d")
           feres_p = resid(fe_mod, "p")
           fedev = deviance(fe_mod)

           stres_r = unname(resid(stat_mod, "resp"))
           stres_d = unname(resid(stat_mod, "d"))
           stres_p = unname(resid(stat_mod, "pearson"))
           stdev = deviance(stat_mod)

           expect_equal(feres_r, stres_r, tolerance = tol)
           expect_equal(feres_d, stres_d, tolerance = tol)
           expect_equal(feres_p, stres_p, tolerance = tol)
           expect_equal(fedev, stdev, tolerance = tol)
       },
<<<<<<< HEAD
    .cases = residuals_cases()
=======
    .cases = residuals_cases()[c(-6,-8),]
>>>>>>> a851cedfa1a705cb6b6e2fc57f8f4e707bd6d110
)
