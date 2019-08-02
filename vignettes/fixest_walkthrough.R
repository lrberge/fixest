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

