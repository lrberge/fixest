#
#
# ##### Bug in model.matrix for fixest??
#
# base <- iris
# names(base) <- c("y", "x1", "x2", "x3", "species")
# mod <- feols(y ~ x1 + x2 | x3 + species, data = base)
# model.matrix(mod, type = "fixef") # Error -> Estimation does not contain fixed-effects
# mod$fixef_vars # There are fixed effects
#
# # Trying different models:
# mod <- feglm(y ~ x1 + x2 | x3 + species, data = base)
# model.matrix(mod, type = "fixef")
# mod <- femlm(y ~ x1 + x2 | x3 + species, data = base)
# model.matrix(mod, type = "fixef")
#
# mod <- feols(y ~ x1 + x2 | species, data = base)
# model.matrix(mod, type = "fixef")
# mod <- feglm(y ~ x1 + x2 | species, data = base)
# model.matrix(mod, type = "fixef")
# mod <- femlm(y ~ x1 + x2 | species, data = base)
# model.matrix(mod, type = "fixef")
#
#
# base$fe <- base$species
# mod <- feols(y ~ x1 | fe, data = base)
# model.matrix(mod, type = "fixef")
# mod <- feglm(y ~ x1 | fe, data = base)
# model.matrix(mod, type = "fixef")
# mod <- femlm(y ~ x1 | fe, data = base)
# model.matrix(mod, type = "fixef")
# ## Got the same error over and over
#
# ### lm fit
# base <- datab()
# k <- 1
# casos <- ols_cases()[k, ]
#
# fml_fixest <- casos$fml_fixest
# fml_stats <- casos$fml_stats
# my_weight <- casos$my_weight
# my_offset <- casos$my_offset
#
# res <- fixest::feols(as.formula(fml_fixest), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
# res_bis <- lm(as.formula(fml_stats), base, weights = ev_par(my_weight), offset = ev_par(my_offset))
# expect_ols_equal(res, res_bis)
#
# object <- res
# reference <- res_bis
#
# ### glm fit
# base <- datab()
# k <- 35 # 25,35,90
# casos <- feglm_cases()[k, ]
#
# fml_fixest <- casos$fml_fixest
# fml_stats <- casos$fml_stats
# my_family <- casos$my_family
# my_weight <- casos$my_weight
# my_offset <- casos$my_offset
#
#
# res <- feglm(as.formula(fml_fixest), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
# if (!is.null(res$obs_selection$obsRemoved)) {
#   qui <- res$obs_selection$obsRemoved
#   base <- base[qui, ]
#   res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight)[qui], offset = ev_par(my_offset)[qui])
# } else {
#   res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
# }
# local_edition(2)
# expect_glm_equal(res, res_bis)
#
#
# #### negbin fit
# K <- 1
# fml_fixest <- fenegbin_cases()[K, ]$fml_fixest
# fml_stats <- fenegbin_cases()[K, ]$fml_stats
#
# res <- fenegbin(as.formula(fml_fixest), base, notes = FALSE)
# res_bis <- MASS::glm.nb(as.formula(fml_stats), base)
# expect_negbin_equal(res, res_bis)
#
# object <- res
# reference <- res_bis
#
# ##########
#
#
# adj <- 0
# tol <- 3e-5
#
# abs((pvalue(object, se = "st", dof = dof(adj = adj))["x1"] - pvalue(reference)["x1"]) / pvalue(reference)["x1"])
# abs((pvalue(object, se = "st", dof = dof(adj = adj))["x1"] - pvalue(reference)["x1"]))
#
#
# local_edition(2)
# expect_equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
#   pvalue(reference)["x1"],
#   tolerance = tol,
#   scale = 1
# )
#
# test(
#   pvalue(res, se = "st", dof = dof(adj = adj))["x1"],
#   pvalue(res_bis)["x1"], "~",
#   tol
# )
#
#
# all.equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
#   pvalue(reference)["x1"],
#   tolerance = tol
# )
# all.equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
#   pvalue(reference)["x1"],
#   tolerance = tol,
#   scale = 1
# )
#
# compare(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
#   pvalue(reference)["x1"],
#   tolerance = tol,
#   scale = pvalue(reference)["x1"]
# )
