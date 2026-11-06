test_that("trait suffixes do not collide with existing or generated identifiers", {
  expect_identical(
    .ensure_unique_variables(c("a", "a", "a-1")),
    c("a-2", "a-3", "a-1")
  )
  expect_identical(
    .ensure_unique_variables(c("a-1", "a", "a", "a-2")),
    c("a-1", "a-3", "a-4", "a-2")
  )
  values <- c("a", "a", "a-1", "a-1", "a-1-1")
  result <- .ensure_unique_variables(values)
  expect_identical(anyDuplicated(result), 0L)
  expect_identical(result[length(result)], "a-1-1")
  expect_identical(.ensure_unique_variables(c("a", "a")), c("a-1", "a-2"))
  expect_identical(.ensure_unique_variables(c("a", "b")), c("a", "b"))
})

test_that("derive_traits keeps colliding protein-site-trait identifiers unique", {
  structures <- rep(glyparse::parse_iupac_condensed(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ), 2)
  mat <- matrix(c(1, 3, 2, 4), nrow = 2, dimnames = list(c("V1", "V2"), c("S1", "S2")))
  exp <- glyexp::GlycoproteomicSE(
    mat,
    rowData = S4Vectors::DataFrame(
      protein = c("P-1", "P"),
      protein_site = c(2L, 1L),
      glycan_structure = structures,
      glycan_composition = glyrepr::as_glycan_composition(structures)
    ),
    colData = S4Vectors::DataFrame(group = c("A", "B")),
    metadata = list(glycan_type = "N")
  )
  traits <- stats::setNames(rep(list(total(TRUE)), 3), c("T", "2-T", "T-1"))
  result <- derive_traits(exp, trait_fns = traits)
  info <- SummarizedExperiment::rowData(result)

  expect_equal(nrow(result), 6)
  expect_identical(anyDuplicated(rownames(result)), 0L)
  singleton <- which(info$protein == "P-1" & info$trait == "T-1")
  expect_identical(rownames(result)[singleton], "P-1-2-T-1")
  expected <- mat[match(info$protein, c("P-1", "P")), , drop = FALSE]
  expect_equal(unname(SummarizedExperiment::assay(result)), unname(expected))
  expect_identical(SummarizedExperiment::assay(exp), mat)
})
