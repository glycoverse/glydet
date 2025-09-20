create_expr_mat <- function() {
  expr_mat <- matrix(
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  expr_mat
}

test_that("prop(nFc > 0) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(nFc = c(1, 0, 1))
  trait_fn <- prop(nFc > 0)
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(2/3, 1/2, 1/2))
})

test_that("prop(nFc > 0 & nS > 0) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    nFc = c(1, 0, 1),
    nS = c(1, 1, 0)
  )

  trait_fn <- prop(nFc > 0 & nS > 0)
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(1/3, 0, 1/2))
})

test_that("prop(nFc > 0, within = (T == 'complex')) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    nFc = c(1, 0, 1),
    T = c("complex", "complex", "hybrid")
  )

  trait_fn <- prop(nFc > 0, within = (T == "complex"))
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(1/2, 0, 1/2))
})

test_that("prop() handles Inf", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    nFc = c(1, 0, 1),
    T = c("hybrid", "hybrid", "hybrid")  # no complex here
  )

  trait_fn <- prop(nFc > 0, within = (T == "complex"))
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, rep(NA_real_, 3))
})

test_that("prop() handles NA through na_action = 'zero'", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    nFc = c(1, 0, 1),
    T = c("hybrid", "hybrid", "hybrid")  # no complex here
  )

  trait_fn <- prop(nFc > 0, within = (T == "complex"), na_action = "zero")
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, rep(0, 3))
})

test_that("prop() handles NA in expr_mat", {
  expr_mat <- matrix(
    c(NA, 0, 1,
      1, 1, 1,
      1, 1, 0),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")

  mp_tbl <- tibble::tibble(
    nFc = c(1, 0, 1),
  )

  trait_fn <- prop(nFc > 0)
  result <- trait_fn(expr_mat, mp_tbl)
  # The first sample is 1/2, because NA is ignored
  expect_equal(result, c(1/2, 1/2, 1/2))
})

test_that("prop() handles NA in cond and within", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(nFc = c(NA, 0, 1))
  trait_fn <- prop(nFc > 0)
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(1/3, 1/2, 0))
})

test_that("ratio(T == 'complex, T == 'hybrid') works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(T = c("complex", "complex", "hybrid"))
  trait_fn <- ratio(T == "complex", T == "hybrid")
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(2, 1, NA))
})

test_that("ratio(B, !B, within = (T == 'complex')) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    B = c(1, 0, 1),
    T = c("complex", "complex", "hybrid")
  )
  trait_fn <- ratio(B, !B, within = (T == "complex"))
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(1, 0, 1))
})

test_that("wmean(nA) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(nA = c(1, 2, 3))
  trait_fn <- wmean(nA)
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(2, 2.5, 1.5))
})

test_that("wmean(nA, within = (T == 'complex')) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    nA = c(1, 2, 3),
    T = c("complex", "complex", "hybrid")
  )
  trait_fn <- wmean(nA, within = (T == "complex"))
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(1.5, 2, 1.5))
})

test_that("wmean(nS / nA) works", {
  expr_mat <- create_expr_mat()
  mp_tbl <- tibble::tibble(
    nA = c(1, 2, 3),
    nS = c(0, 1, 3)
  )
  trait_fn <- wmean(nS / nA)
  result <- trait_fn(expr_mat, mp_tbl)
  expect_equal(result, c(0.5, 0.75, 0.25))
})

test_that("prop print", {
  expect_snapshot(print(prop(nFc > 0)))
  expect_snapshot(print(prop(nFc > 0, within = T == "complex")))
  expect_snapshot(print(prop(nFc > 0, within = (T == "complex"))))
  expect_snapshot(print(prop(nFc > 0, within = NULL)))
  expect_snapshot(print(prop(nFc > 0, na_action = "zero")))
})

test_that("ratio print", {
  expect_snapshot(print(ratio(T == "complex", T == "hybrid")))
  expect_snapshot(print(ratio(T == "complex", T == "hybrid", within = T == "complex")))
  expect_snapshot(print(ratio(T == "complex", T == "hybrid", within = (T == "complex"))))
  expect_snapshot(print(ratio(T == "complex", T == "hybrid", within = NULL)))
  expect_snapshot(print(ratio(T == "complex", T == "hybrid", na_action = "zero")))
})

test_that("wmean print", {
  expect_snapshot(print(wmean(nA)))
  expect_snapshot(print(wmean(nA, within = T == "complex")))
  expect_snapshot(print(wmean(nA, within = (T == "complex"))))
  expect_snapshot(print(wmean(nA, within = NULL)))
  expect_snapshot(print(wmean(nA, na_action = "zero")))
})