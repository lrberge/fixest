base <- datab8()

## All in one test
# test_that("to_integer creates the correct number of elements",
#           {
#               m1 <- to_integer(base$species)
#               m2 <- to_integer(base$species, base$z)
#               m3 <- to_integer(base$species_na, base$z)
#               m4 <- to_integer(base$species_na, base$z, add_items = TRUE, items.list = TRUE)
#               M = c(length(unique(m1)),
#                     length(unique(m2)),
#                     length(unique(m3)),
#                     length(unique(m4$items)))
#               expect_equal(M,c(3,15,11,10))
#
#           })


## Separate in different tests?

test_that("to_integer creates the correct number of elements", {
  m <- to_integer(base$species)
  expect_equal(length(unique(m)), 3)
})
test_that("to_integer creates the correct number of elements with combined factors", {
  m <- to_integer(base$species, base$z)
  expect_equal(length(unique(m)), 15)
})
test_that("to_integer creates the correct number of elements with NA presence", {
  m <- to_integer(base$species_na, base$z)
  expect_equal(length(unique(m)), 11)
})
test_that("to_integer creates the correct number of elements with NA presence and added items", {
  m <- to_integer(base$species_na, base$z, add_items = TRUE, items.list = TRUE)
  expect_equal(length(unique(m$items)), 10)
})
