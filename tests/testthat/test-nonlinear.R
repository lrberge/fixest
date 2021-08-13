
base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$tab <- c("versicolor" = 5, "setosa" = 0, "virginica" = -5)
base$var_spec <- as.numeric(base$tab[base$species])

# a and b are te parameters to be estimated
fun_nl <- function(a, b, spec) {
    res <- as.numeric(base$tab[spec])
    a * res + b * res^2
}

test_that("estimation of a linear model between feNmlm and feols are equal",
          {
              est_nl <- feNmlm(y ~ x1, base, NL.fml = ~ fun_nl(a, b, species), NL.start = 1, family = "gaussian")
              est_lin <- feols(y ~ x1 + var_spec + I(var_spec^2), base)

              coef_nl = coef(est_nl)
              coef_lin = coef(est_lin)[c(3, 4, 1, 2)]
              expect_equal2(unname(coef_nl), unname(coef_lin))
          })

## standard errors are different!
