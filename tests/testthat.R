library(testthat)
library(fixest)
library(sandwich)
library(patrick)

test_check("fixest")
devtools::test()
