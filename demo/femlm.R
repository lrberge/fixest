# Example of femlm using trade data.

# loading data
data(trade)

# Estimation with three sets of fixed-effects (Poisson)
est_pois = femlm(Euros ~ log(dist_km) | Origin + Destination + Product, trade)

# the result with two way clustered standard errors:
summary(est_pois, se = "twoway")

# estimation with Gaussian likelihood (in log-log)
est_gaus = femlm(log(Euros) ~ log(dist_km) | Origin + Destination + Product, trade, family = "gaussian")

# displaying the two sets of results
esttable(est_pois, est_gaus, se = "twoway")

# updating the first result (adding one variable and deleting one cluster)
update(est_pois, . ~ . + log(Year) | . - Product)

# obtaining the set of fixed-effects of the first estimation
fe_pois = fixef(est_pois)
# displaying the 6 (2*3) most notable fixed-effects:
plot(fe_pois, n = 3)


