test_that("basic_traits() works", {
  traits <- basic_traits()
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF") %in% names(traits)))
})

test_that("basic_traits() with sia_link = TRUE works", {
  expect_snapshot(traits <- basic_traits(sia_link = TRUE))
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF", "GE", "GL") %in% names(traits)))
})

test_that("all_traits() works", {
  traits <- all_traits()
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF", "A2Fc") %in% names(traits)))
})

test_that("all_traits() with sia_link = TRUE works", {
  expect_snapshot(traits <- all_traits(sia_link = TRUE))
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF", "A2Fc", "GE", "GL", "A2E", "A2L") %in% names(traits)))
})