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

test_that("traits_clerc_2018() works", {
  traits <- traits_clerc_2018()
  expect_true(is.list(traits))
  expect_true(all(c("A2F0GS", "A3F0GS", "A4F0GS") %in% names(traits)))
})

test_that("traits_clerc_2018() with sia_link = TRUE works", {
  expect_snapshot(traits <- traits_clerc_2018(sia_link = TRUE))
  expect_true(is.list(traits))
  expect_true(all(c("A2F0GS", "A3F0GS", "A4F0GS", "A2L0F", "A3L0F") %in% names(traits)))
})

test_that("traits_fu_2026() works", {
  traits <- traits_fu_2026()
  expect_true(is.list(traits))
  expect_true(all(c("TM", "THy", "TC") %in% names(traits)))
})

test_that("traits_li_2025() works", {
  traits <- traits_li_2025()
  expect_true(is.list(traits))
  expect_true(all(c("S1", "S2", "S3", "S4") %in% names(traits)))
})