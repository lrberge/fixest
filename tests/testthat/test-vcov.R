# vcovClust <- fixest:::vcovClust
#
# ### Standard SEs
# base <- datab3()
#
# method <- c("ols", "glm", "femlm", "feNmlm")
# y_dep <- c("y", "y_int", "y_int", "y_int")
# fmly <- c("NULL", rep("poisson", 3))
#
# ## Fitting feols, feglm, femlm and feNmlm
# Est <- list()
# VCOVs_raw <- list()
# ns <- list()
# for (k in 1:length(method)) {
#   fmla <- xpd(fml = y_dep ~ x | fe1 + fe2, lhs = y_dep[k])
#   Est[[k]] <- fixest_mod_select(model = method[k], fmla = fmla, base = base, famly = fmly[k])
#   ns[[k]] <- Est[[k]]$nobs
#   VCOVs_raw[[k]] <- vcov(Est[[k]], se = "standard") / ((ns[[k]] - 1) / (ns[[k]] - Est[[k]]$nparams))
# }
# names(Est) <- method
# names(VCOVs_raw) <- method
#
# with_parameters_test_that("Standard VarCov estimations are correct",
#                           {
#                             est <- Est[[method]]
#                             VCOV_raw <- VCOVs_raw[[method]]
#                             vcov_est <- vcov(est, se = "standard", ssc = ssc(adj = adj, fixef.K = k_val))
#                             vcov_raw <- VCOV_raw * my_adj
#                             expect_equal(vcov_est, vcov_raw)
#                           },
#                           .cases = vcov_cases1()
# )
#
#
# ########## Clustered SEs
# VCOVs_raw <- list()
# Hs <- list()
# for (k in 1:length(method)) {
#   if (k == 1) {
#     VCOVs_raw[[k]] <- vcov(Est[[k]], se = "standard") / Est[[k]]$sigma2
#   } else {
#     VCOVs_raw[[k]] <- vcov(Est[[k]], se = "standard") # / est$dispersion
#   }
#   Hs[[k]] <- vcovClust(Est[[k]]$fixef_id$fe1, VCOVs_raw[[k]], scores = Est[[k]]$scores, ssc = FALSE)
# }
# names(VCOVs_raw) <- method
# names(Hs) <- method
#
# with_parameters_test_that("Clustered VarCov estimations are correct",
#                           {
#                             vcov_est <- vcov(Est[[method]], se = "cluster", ssc = ssc(adj = adj, fixef.K = k_val, cluster.adj = c_adj))
#                             vcov_raw <- Hs[[method]] * cluster_factor * my_adj
#                             expect_equal(as.numeric(vcov_est), as.numeric(vcov_raw))
#                           },
#                           .cases = vcov_cases2()[1:24, ] # works only for ols fitting
# )
#
# ############# Two-way clustered se
# VCOVs_raw <- list()
# Ms_i <- list()
# Ms_t <- list()
# Ms_it <- list()
# for (k in 1:length(method)) {
#   est <- Est[[k]]
#   if (k == 1) {
#     VCOVs_raw[[k]] <- vcov(est, se = "standard") / est$sigma2
#   } else {
#     VCOVs_raw[[k]] <- vcov(est, se = "standard") # / est$dispersion
#   }
#   Ms_i[[k]] <- vcovClust(est$fixef_id$fe1, VCOVs_raw[[k]], scores = est$scores, ssc = FALSE)
#   Ms_t[[k]] <- vcovClust(est$fixef_id$fe2, VCOVs_raw[[k]], scores = est$scores, ssc = FALSE)
#   Ms_it[[k]] <- vcovClust(paste(base$fe1, base$fe2), VCOVs_raw[[k]], scores = est$scores, ssc = FALSE, do.unclass = TRUE)
# }
# names(VCOVs_raw) <- method
# names(Ms_i) <- method
# names(Ms_t) <- method
# names(Ms_it) <- method
#
# with_parameters_test_that("Two-Way Clustered VarCov estimations are correct", {
#   est <- Est[[method]]
#   vcov_est <- vcov(est, se = "two", ssc = ssc(adj = adj, fixef.K = k_val, cluster.adj = c_adj, cluster.df = cdf))
#   vcov_raw <- V_matrix(Ms_i[[method]], Ms_t[[method]], Ms_it[[method]], c_adj, cdf) * my_adj
#   expect_equal(as.numeric(vcov_est), as.numeric(vcov_raw))
# },
# .cases = vcov_cases3()[1:48, ] # Only works for ols
# )
#
# base <- iris
# names(base) <- c("y", "x1", "x2", "x3", "species")
# base$clu <- sample(6, 150, TRUE)
# base$clu[1:5] <- NA
# base$y_int <- round(base$y)
# base$y_01 <- c(base$y > mean(base$y)) + 0
#
# with_parameters_test_that("vcov estimation from different sources are equal", {
#   fmla <- xpd(lhs ~ x1 | species, lhs = y_dep)
#   est <- fixest_mod_select(model = method,
#                            fmla = fmla,
#                            base = base,
#                            famly = fmly,
#                            cluster = ~clu, ssc = ssc(adj = FALSE))
#   v1 <- est$cov.scaled
#   v1b <- vcov(est)
#   v1c <- summary(est)$cov.scaled
#
#   expect_equal(as.numeric(v1), as.numeric(v1b))
#   expect_equal(as.numeric(v1), as.numeric(v1c))
#
#   # Only ssc change
#   v2 <- summary(est, ssc = ssc())$cov.scaled
#   v2b <- vcov(est, cluster = ~clu, ssc = ssc())
#
#   expect_equal(as.numeric(v2), as.numeric(v2b))
#   expect_true(max(abs(v1 - v2)) != 0)
#
#   # SE change only
#   v3 <- summary(est, se = "hetero")$cov.scaled
#   v3b <- vcov(est, se = "hetero", ssc = ssc(adj = FALSE))
#
#   expect_equal(as.numeric(v3), as.numeric(v3b))
#   expect_true(max(abs(v1 - v3)) != 0)
#   expect_true(max(abs(v2 - v3)) != 0)
# },
# .cases = vcov_cases4()
# )
