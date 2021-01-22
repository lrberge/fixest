## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(fixest)
setFixest_nthreads(1)

## -----------------------------------------------------------------------------
base = iris
names(base) = c("y1", "y2", "x1", "x2", "species")

res_multi = feols(c(y1, y2) ~ x1 + csw(x2, x2^2) | sw0(species), base, fsplit = ~species)


## -----------------------------------------------------------------------------
summary(res_multi, "compact", se = "hetero")

## -----------------------------------------------------------------------------
etable(feols(c(y1, y2) ~ x1 + x2, base))

## -----------------------------------------------------------------------------
etable(feols(y1 ~ csw(x1, x2) | sw0(species), base, cluster = ~species))

## -----------------------------------------------------------------------------
etable(feols(y1 ~ x1 + x2, base, fsplit = ~species))

## -----------------------------------------------------------------------------
etable(res_multi[lhs = 1, fixef = 1, rhs = TRUE, sample = -1])

