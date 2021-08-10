
library(fixest)
base = datab2()

setFixest_notes(FALSE)
# ols
res = feols(y ~ -1 + x1 + i(fe1), base)
expect_model_nointercept(res)
res = feols(y ~ -1 + x1 + factor(fe1), base)
expect_model_nointercept(res)
res = feols(y ~ -1 + x1 + i(fe1) + i(fe2), base)
expect_model_nointercept(res)

## feglm
res = feglm(y ~ -1 + x1 + i(fe1), base)
expect_model_nointercept(res)
res = feglm(y ~ -1 + x1 + factor(fe1), base)
expect_model_nointercept(res)
res = feglm(y ~ -1 + x1 + i(fe1) + i(fe2), base)
expect_model_nointercept(res)

res = feglm(y_r ~ -1 + x1 + i(fe1), base, family = "poisson")
expect_model_nointercept(res)
res = feglm(y_r ~ -1 + x1 + factor(fe1), base, family = "poisson")
expect_model_nointercept(res)
res = feglm(y_r ~ -1 + x1 + i(fe1) + i(fe2), base, family = "poisson")
expect_model_nointercept(res)

res = feglm(y ~ -1 + x1 + i(fe1), base, family = "Gamma")
expect_model_nointercept(res)
res = feglm(y_r ~ -1 + x1 + factor(fe1), base, family = "Gamma")
expect_model_nointercept(res)
res = feglm(y_r ~ -1 + x1 + i(fe1) + i(fe2), base, family = "Gamma")
expect_model_nointercept(res)

## fenegbin
res = fenegbin(y ~ -1 + x1 + i(fe1), base)
expect_model_nointercept(res)
res = fenegbin(y ~ -1 + x1 + factor(fe1), base)
expect_model_nointercept(res)
res = fenegbin(y ~ -1 + x1 + i(fe1) + i(fe2), base)
expect_model_nointercept(res)


## Parametrize the previous tests?

setFixest_notes(TRUE)
