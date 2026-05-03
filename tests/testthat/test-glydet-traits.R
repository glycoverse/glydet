test_that("traits_basic() works", {
  traits <- traits_basic()
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF") %in% names(traits)))
})

test_that("traits_basic() with sia_link = TRUE works", {
  expect_snapshot(traits <- traits_basic(sia_link = TRUE))
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF", "GE", "GL") %in% names(traits)))
})

test_that("traits_detailed() works", {
  traits <- traits_detailed()
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF", "A2Fc") %in% names(traits)))
})

test_that("traits_detailed() with sia_link = TRUE works", {
  expect_snapshot(traits <- traits_detailed(sia_link = TRUE))
  expect_true(is.list(traits))
  expect_true(all(c("TM", "CA2", "TF", "A2Fc", "GE", "GL", "A2E", "A2L") %in% names(traits)))
})

test_that("old trait set names remain aliases", {
  expect_identical(basic_traits(), traits_basic())
  expect_identical(suppressMessages(basic_traits(sia_link = TRUE)), suppressMessages(traits_basic(sia_link = TRUE)))
  expect_identical(all_traits(), traits_detailed())
  expect_identical(suppressMessages(all_traits(sia_link = TRUE)), suppressMessages(traits_detailed(sia_link = TRUE)))
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
