## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE-----------------------------------------------------------
library(fixest)
data(trade)


## ---- echo=FALSE, results='asis'-----------------------------------------
tab = head(trade)
knitr::kable(tab)

## ------------------------------------------------------------------------
gravity_results <- feglm(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade)

## ------------------------------------------------------------------------
print(gravity_results)

## ------------------------------------------------------------------------
summary(gravity_results, se = "twoway")

## ------------------------------------------------------------------------
# Equivalent ways of clustering the SEs:
# One-way clustering is deduced from the arguent 'cluster'
# - using the vector:
summary(gravity_results, cluster = trade$Product)
# - by reference:
summary(gravity_results, cluster = "Product")
# - with a formula:
summary(gravity_results, cluster = ~Product)

## ------------------------------------------------------------------------
gravity_simple = feglm(Euros ~ log(dist_km), trade)
# Two way clustering is deduced from the argument 'cluster'
# Using data:
summary(gravity_simple, cluster = trade[, c("Origin", "Destination")])
# Using a formula (note that the values of the variables are 
#  fetched directly in the original database):
summary(gravity_simple, cluster = ~Origin+Destination)

## ------------------------------------------------------------------------
gravity_results_ols <- feols(log(Euros) ~ log(dist_km)|Origin+Destination+Product+Year, trade)

## ------------------------------------------------------------------------
gravity_results_negbin <- fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade)


## ---- eval=FALSE---------------------------------------------------------
#  esttable(gravity_results, gravity_results_negbin, gravity_results_ols, se = "twoway", titles = c("Poisson", "Negative Binomial", "Gaussian"))

## ---- echo=FALSE, results='asis'-----------------------------------------
tab = esttable(gravity_results, gravity_results_negbin, gravity_results_ols, se = "twoway", titles = c("Poisson", "Negative Binomial", "Gaussian"))
# problem to display the second empty line in markdown
knitr::kable(tab[-2, ])

## ------------------------------------------------------------------------
gravity_subcluster = list()
all_clusters = c("Year", "Destination", "Origin", "Product")
for(i in 1:4){
	gravity_subcluster[[i]] = feglm(Euros ~ log(dist_km), trade, fixef = all_clusters[1:i])
}

## ---- eval=FALSE---------------------------------------------------------
#  esttable(gravity_subcluster, se = "twoway", cluster = trade[, c("Origin", "Destination")])

## ---- echo=FALSE, results='asis'-----------------------------------------
tab = esttable(gravity_subcluster, se = "twoway", cluster = trade[, c("Origin", "Destination")])
knitr::kable(tab)

## ------------------------------------------------------------------------
esttex(gravity_subcluster, se = "twoway", cluster = trade[, c("Origin", "Destination")])

## ---- eval=FALSE---------------------------------------------------------
#  # we set the dictionary once and for all
#  myDict = c("log(dist_km)" = "$\\ln (Distance)$", "(Intercept)" = "Constant")
#  # 1st export: we change the signif code and drop the intercept
#  esttex(gravity_subcluster, signifCode = c("a" = 0.01, "b" = 0.05), drop = "Int", dict = myDict, file = "Estimation Table.tex", replace = TRUE, title = "First export -- normal Standard-errors")
#  # 2nd export: clustered S-E + distance as the first coefficient
#  esttex(gravity_subcluster, se = "cluster", cluster = trade$Product, order = "dist", dict = myDict, file = "Estimation Table.tex", title = "Second export -- clustered standard-errors (on Product variable)")
#  

## ------------------------------------------------------------------------
fixedEffects <- fixef(gravity_results)
summary(fixedEffects)

## ------------------------------------------------------------------------
fixedEffects$Year

## ---- fig.width=7--------------------------------------------------------
plot(fixedEffects)

## ------------------------------------------------------------------------
base_vs = iris
names(base_vs) = c(paste0("x", 1:4), "species")

## ------------------------------------------------------------------------
est_vs = feols(x1 ~ x2 | species[x3], base_vs)
est_vs

## ------------------------------------------------------------------------
summary(fixef(est_vs))

## ------------------------------------------------------------------------
# we create another "fixed-effect"
base_vs$fe = rep(1:5, 30)
head(base_vs)

## ------------------------------------------------------------------------
est_comb = feols(x1 ~ x2 | species^fe, base_vs)
est_comb

## ------------------------------------------------------------------------
fixef(est_comb)[[1]]

## ------------------------------------------------------------------------
# Sample data illustrating the DiD
data(base_did)
head(base_did)
# Estimation of yearly effect (they are automatically added)
# We also add individual/time fixed-effects:
est_did = did_estimate_yearly_effects(y ~ x1 | id + period, base_did,
                                      treat_time = ~treat+period, reference = 5)
est_did

## ---- fig.width=7--------------------------------------------------------
did_plot_yearly_effects(est_did)

## ------------------------------------------------------------------------
base_lag = base_did
# we create a lagged value of the variable x1
base_lag$x1.l1 = lag(x1~id+period, 1, base_lag)
head(base_lag)

## ------------------------------------------------------------------------
library(data.table)
base_lag_dt = as.data.table(base_did)
# we create a lagged value of the variable x1
base_lag_dt[, x1.l1 := lag(x1~id+period, 1)]

## ------------------------------------------------------------------------
# Generating data:
n = 1000
# x and y: two positive random variables
x = rnorm(n, 1, 5)**2
y = rnorm(n, -1, 5)**2
# E(z) = 2*x + 3*y and some noise
z = rpois(n, 2*x + 3*y) + rpois(n, 1)
base = data.frame(x, y, z)

## ------------------------------------------------------------------------
result_NL = feNmlm(z~0, base, NL.fml = ~ log(a*x + b*y), NL.start = list(a=1, b=1), lower = list(a=0, b=0))

## ------------------------------------------------------------------------
print(result_NL)

## ------------------------------------------------------------------------
# the class of each observation
id = sample(20, n, replace = TRUE)
base$id = id
# the vector of fixed-effects
gamma = rnorm(20)**2
# the new vector z_bis
z_bis = rpois(n, gamma[id] * (2*x + 3*y)) + rpois(n, 1)
base$z_bis = z_bis

## ------------------------------------------------------------------------
# we add the fixed-effect in the formula
result_NL_fe = feNmlm(z_bis~0|id, base, NL.fml = ~ log(2*x + b*y), NL.start = list(b=1), lower = list(b=0))
# The coef should be around 3
coef(result_NL_fe)
# the gamma and the exponential of the fixed-effects should be similar
rbind(gamma, exp(fixef(result_NL_fe)$id))


## ---- eval = FALSE-------------------------------------------------------
#  # Sample of results:
#  # 1 nthreads: 3.13s
#  system.time(fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade, nthreads = 1))
#  # 2 nthreads: 1.82s
#  system.time(fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade, nthreads = 2))
#  # 4 nthreads: 1.17s
#  system.time(fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade, nthreads = 4))

## ------------------------------------------------------------------------
base_coll = trade
base_coll$constant_variable = 1
res <- femlm(Euros ~ log(dist_km) + constant_variable|Origin+Destination+Product+Year, base_coll)
collinearity(res)


