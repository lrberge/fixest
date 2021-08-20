base <- datab8()

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
