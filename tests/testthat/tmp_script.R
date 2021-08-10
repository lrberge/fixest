
feglm_cases()[19:20, ]
K <- 18
fml_fixest <- feglm_cases()[K, ]$fml_fixest
fml_stats <- feglm_cases()[K, ]$fml_stats
my_family <- feglm_cases()[K, ]$my_family
my_offset <- feglm_cases()[K, ]$my_offset
my_weight <- feglm_cases()[K, ]$my_weight

rm(fml_fixest, fml_stats, my_family, my_offset, my_weight)





k = 35 # 25,35,90
casos = feglm_cases()[k,]

fml_fixest = casos$fml_fixest
fml_stats = casos$fml_stats
my_family = casos$my_family
my_weight = casos$my_weight
my_offset = casos$my_offset


res <- feglm(as.formula(fml_fixest), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
if (!is.null(res$obs_selection$obsRemoved)) {
    qui <- res$obs_selection$obsRemoved
    base <- base[qui, ]
    res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight)[qui], offset = ev_par(my_offset)[qui])
} else {
    res_bis <- glm(as.formula(fml_stats), base, family = my_family, weights = ev_par(my_weight), offset = ev_par(my_offset))
}
local_edition(2)
expect_glm_equal(res, res_bis)



adj = 0
tol = 3e-5

abs((pvalue(object, se = "st", dof = dof(adj = adj))["x1"] - pvalue(reference)["x1"])/pvalue(reference)["x1"])
abs((pvalue(object, se = "st", dof = dof(adj = adj))["x1"] - pvalue(reference)["x1"]))


local_edition(2)
expect_equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
             pvalue(reference)["x1"],
             tolerance = tol,
             scale = 1
)

test(pvalue(res, se = "st", dof = dof(adj = adj))["x1"],
     pvalue(res_bis)["x1"], "~",
     tol)


all.equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
          pvalue(reference)["x1"],
          tolerance = tol)
all.equal(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
          pvalue(reference)["x1"],
          tolerance = tol,
          scale = 1)

compare(pvalue(object, se = "st", dof = dof(adj = adj))["x1"],
        pvalue(reference)["x1"],
        tolerance = tol,
        scale = pvalue(reference)["x1"])
