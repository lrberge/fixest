
base <- datab11()
setFixest_notes(FALSE)
patrick::with_parameters_test_that("estimation with subset works properly",
  {
    res_sub_a <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base, subset = ~ species == "setosa")
    res_sub_b <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base, subset = base$species == "setosa")
    res_sub_c <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base, subset = which(base$species == "setosa"))
    res <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base[base$species == "setosa", ])

    test_that("coeficients are equal", {
      expect_equal(coef(res_sub_a), coef(res))
      expect_equal(coef(res_sub_b), coef(res))
      expect_equal(coef(res_sub_c), coef(res))
    })
    test_that("standard errors are equal", {
      expect_equal(se(res_sub_c, cluster = "fe_bis"), se(res, cluster = "fe_bis"))
    })
    expect_equal(1, 1)
  },
  .cases = subset_cases()[19:22, ]
)

# K = 20
# casos = subset_cases()[K,]
#
# method = casos$method
# fmlas = casos$fmlas
# famly = casos$famly
#
# res_sub_a <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base, subset = ~ species == "setosa")
# res_sub_b <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base, subset = base$species == "setosa")
# res_sub_c <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base, subset = which(base$species == "setosa"))
# res <- fixest_mod_select(model = method, fmla = fmlas, famly = famly, base = base[base$species == "setosa", ])
#
# test_that("coeficients are equal", {
#     expect_equal(coef(res_sub_a), coef(res))
#     expect_equal(coef(res_sub_b), coef(res))
#     expect_equal(coef(res_sub_c), coef(res))
# })
# test_that("standard errors are equal", {
#     expect_equal(se(res_sub_c, cluster = "fe_bis"), se(res, cluster = "fe_bis"))
# })
