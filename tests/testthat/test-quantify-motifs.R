glycomics_exp <- function() {
  expr_mat <- matrix(
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3"),
    glycan_structure = glyparse::parse_iupac_condensed(c(
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # No Sia
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # One Sia
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # Two Sia
    )),
    glycan_composition = glyrepr::as_glycan_composition(glycan_structure)
  )
  glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N")
}

glycoproteomics_exp <- function() {
  expr_mat <- matrix(
    c(1, 0, 1,  # Site 1, Glycan 1
      1, 1, 1,  # Site 1, Glycan 2
      1, 1, 0,  # Site 1, Glycan 3
      1, 1, 1,  # Site 2, Glycan 1
      1, 1, 1), # Site 2, Glycan 2
    nrow = 5, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3", "V4", "V5")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3", "V4", "V5"),
    protein = c("P1", "P1", "P1", "P1", "P1"),
    protein_site = c(1L, 1L, 1L, 2L, 2L),
    glycan_structure = glyparse::parse_iupac_condensed(c(
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Glycan 1: No Sia
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Glycan 2: One Sia
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Glycan 3: Two Sia
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Glycan 1: No Sia
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # Glycan 2: One Sia
    )),
    glycan_composition = glyrepr::as_glycan_composition(glycan_structure)
  )
  glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycoproteomics", glycan_type = "N")
}

test_that("quantify_motifs works for glycomics data absolutely", {
  exp <- glycomics_exp()
  motifs <- c(
    sia = "Neu5Ac(a2-3)Gal(b1-",  # 0, 1, 2
    core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # 1, 1, 1
  )
  result <- quantify_motifs(exp, motifs, method = "absolute")
  expected_expr_mat <- matrix(
    c(3, 3, 1,
      3, 2, 2),
    nrow = 2, byrow = TRUE
  )
  colnames(expected_expr_mat) <- colnames(exp$expr_mat)
  rownames(expected_expr_mat) <- c("V1", "V2")
  expect_equal(result$expr_mat, expected_expr_mat)
  expect_equal(result$var_info$motif, c("sia", "core"))
})

test_that("quantify_motifs works for glycomics data relatively", {
  exp <- glycomics_exp()
  motifs <- c(
    sia = "Neu5Ac(a2-3)Gal(b1-",  # 0, 1, 2
    core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # 1, 1, 1
  )
  result <- quantify_motifs(exp, motifs, method = "relative")
  expected_expr_mat <- matrix(
    c(3/3, 3/2, 1/2,
      3/3, 2/2, 2/2),
    nrow = 2, byrow = TRUE
  )
  colnames(expected_expr_mat) <- colnames(exp$expr_mat)
  rownames(expected_expr_mat) <- c("V1", "V2")
  expect_equal(result$expr_mat, expected_expr_mat)
  expect_equal(result$var_info$motif, c("sia", "core"))
})

test_that("quantify_motifs works for glycoproteomics data absolutely", {
  # Set up
  exp <- glycoproteomics_exp()
  motifs <- c(
    sia = "Neu5Ac(a2-3)Gal(b1-",  # 0, 1, 2
    core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # 1, 1, 1
  )

  # Perform
  result <- quantify_motifs(exp, motifs, method = "absolute")

  # Test
  expected_expr_mat <- matrix(
    c(3, 3, 1,
      3, 2, 2,
      1, 1, 1,
      2, 2, 2),
    nrow = 4, byrow = TRUE
  )
  colnames(expected_expr_mat) <- colnames(exp$expr_mat)
  rownames(expected_expr_mat) <- c("V1", "V2", "V3", "V4")
  expect_equal(result$expr_mat, expected_expr_mat)

  # Check non-structure columns match
  expect_equal(result$var_info$variable, c("V1", "V2", "V3", "V4"))
  expect_equal(result$var_info$protein, c("P1", "P1", "P1", "P1"))
  expect_equal(result$var_info$protein_site, c(1, 1, 2, 2))
  expect_equal(result$var_info$motif, c("sia", "core", "sia", "core"))

  # Check motif_structure column exists and has correct type
  expect_true("motif_structure" %in% colnames(result$var_info))
  expect_true(glyrepr::is_glycan_structure(result$var_info$motif_structure))
  expect_equal(as.character(result$var_info$motif_structure),
               as.character(glyparse::parse_iupac_condensed(c(
                 "Neu5Ac(a2-3)Gal(b1-",
                 "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
                 "Neu5Ac(a2-3)Gal(b1-",
                 "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
               ))))
})

test_that("quantify_motifs works for glycoproteomics data relatively", {
  # Set up
  exp <- glycoproteomics_exp()
  motifs <- c(
    sia = "Neu5Ac(a2-3)Gal(b1-",  # 0, 1, 2
    core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # 1, 1, 1
  )

  # Perform
  result <- quantify_motifs(exp, motifs, method = "relative")

  # Test
  expected_expr_mat <- matrix(
    c(3/3, 3/2, 1/2,
      3/3, 2/2, 2/2,
      1/2, 1/2, 1/2,
      2/2, 2/2, 2/2),
    nrow = 4, byrow = TRUE
  )
  colnames(expected_expr_mat) <- colnames(exp$expr_mat)
  rownames(expected_expr_mat) <- c("V1", "V2", "V3", "V4")
  expect_equal(result$expr_mat, expected_expr_mat)

  # Check non-structure columns match
  expect_equal(result$var_info$variable, c("V1", "V2", "V3", "V4"))
  expect_equal(result$var_info$protein, c("P1", "P1", "P1", "P1"))
  expect_equal(result$var_info$protein_site, c(1, 1, 2, 2))
  expect_equal(result$var_info$motif, c("sia", "core", "sia", "core"))

  # Check motif_structure column exists and has correct type
  expect_true("motif_structure" %in% colnames(result$var_info))
  expect_true(glyrepr::is_glycan_structure(result$var_info$motif_structure))
  expect_equal(as.character(result$var_info$motif_structure),
               as.character(glyparse::parse_iupac_condensed(c(
                 "Neu5Ac(a2-3)Gal(b1-",
                 "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
                 "Neu5Ac(a2-3)Gal(b1-",
                 "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
               ))))
})

test_that("quantify_motifs keeps additional columns for glycoproteomics data", {
  # Set up
  exp <- glycoproteomics_exp() |>
    glyexp::mutate_var(gene = c("G1", "G1", "G1", "G1", "G1"))
  motifs <- c(sia = "Neu5Ac(a2-3)Gal(b1-")

  # Perform
  result <- quantify_motifs(exp, motifs)

  # Test
  expect_setequal(colnames(result$var_info), c("variable", "protein", "protein_site", "motif", "motif_structure", "gene"))
})

test_that("quantify_motifs adds motif_structure column for glycomics data", {
  exp <- glycomics_exp()
  motifs <- c(
    sia = "Neu5Ac(a2-3)Gal(b1-",  # 0, 1, 2
    core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # 1, 1, 1
  )
  result <- quantify_motifs(exp, motifs, method = "absolute")

  # Test motif_structure column exists and has correct type
  expect_true("motif_structure" %in% colnames(result$var_info))
  expect_true(glyrepr::is_glycan_structure(result$var_info$motif_structure))

  # Test motif_structure values match the input motifs (compare string representations)
  expect_equal(as.character(result$var_info$motif_structure[1]), "Neu5Ac(a2-3)Gal(b1-")
  expect_equal(as.character(result$var_info$motif_structure[2]), "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
})

test_that("quantify_motifs adds motif_structure column for glycoproteomics data", {
  exp <- glycoproteomics_exp()
  motifs <- c(
    sia = "Neu5Ac(a2-3)Gal(b1-",  # 0, 1, 2
    core = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # 1, 1, 1
  )
  result <- quantify_motifs(exp, motifs, method = "absolute")

  # Test motif_structure column exists and has correct type
  expect_true("motif_structure" %in% colnames(result$var_info))
  expect_true(glyrepr::is_glycan_structure(result$var_info$motif_structure))

  # Check structures are in correct order (sia, core, sia, core)
  expect_equal(as.character(result$var_info$motif_structure[1]), "Neu5Ac(a2-3)Gal(b1-")
  expect_equal(as.character(result$var_info$motif_structure[2]), "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  expect_equal(as.character(result$var_info$motif_structure[3]), "Neu5Ac(a2-3)Gal(b1-")
  expect_equal(as.character(result$var_info$motif_structure[4]), "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
})

test_that("quantify_motifs works with known motif names", {
  exp <- glycomics_exp()
  motifs <- c("N-Glycan core basic", "N-Glycan high mannose")

  result <- quantify_motifs(exp, motifs, method = "absolute")

  # Test motif_structure column exists
  expect_true("motif_structure" %in% colnames(result$var_info))
  expect_true(glyrepr::is_glycan_structure(result$var_info$motif_structure))

  # Test motifs are correctly named
  expect_equal(result$var_info$motif, c("N-Glycan core basic", "N-Glycan high mannose"))
})

test_that("quantify_motifs works with glycan_structure vector", {
  exp <- glycomics_exp()
  motifs <- glyparse::parse_iupac_condensed(c(
    "Neu5Ac(a2-3)Gal(b1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))

  result <- quantify_motifs(exp, motifs, method = "absolute")

  # Test motif_structure column exists
  expect_true("motif_structure" %in% colnames(result$var_info))
  expect_true(glyrepr::is_glycan_structure(result$var_info$motif_structure))

  # Test motifs are correctly named (using structure strings as names)
  expect_equal(result$var_info$motif, c(
    "Neu5Ac(a2-3)Gal(b1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
})