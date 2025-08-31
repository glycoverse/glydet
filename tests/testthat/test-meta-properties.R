test_that("get_meta_properties works with default meta-property functions", {
  glycans <- c(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  )
  result <- get_meta_properties(glycans)
  expected <- tibble::tibble(
    T = factor(
      c("paucimannose", "complex", "complex"),
      levels = c("paucimannose", "hybrid", "highmannose", "complex")
    ),
    B = c(FALSE, FALSE, FALSE),
    nA = c(0L, 2L, 1L),
    nFc = c(0L, 0L, 0L),
    nFa = c(0L, 1L, 0L),
    nG = c(0L, 0L, 0L),
    nGt = c(0L, 0L, 0L),
    nS = c(0L, 0L, 0L),
    nM = c(3L, 3L, 3L)
  )
  expect_equal(result, expected)
})

test_that("get_meta_properties works with custom meta-property functions", {
  glycans <- c(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  )
  result <- get_meta_properties(glycans, list(
    nN = ~ glyrepr::count_mono(.x, "HexNAc"),
    nH = ~ glyrepr::count_mono(.x, "Hex")
  ))
  expected <- tibble::tibble(
    nN = c(2L, 4L, 3L),
    nH = c(3L, 3L, 3L)
  )
  expect_equal(result, expected)
})