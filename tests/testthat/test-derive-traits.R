test_that("derive_traits works for glycomics experiments", {
  # Construct a test experiment
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3"),
    glycan_structure = glyparse::parse_iupac_condensed(c(
      "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
    ))
  )
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  expr_mat <- matrix(
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N")

  # Calculate derived traits
  trait_fns <- list(
    TFc = prop(nFc > 0),
    TC = prop(T == "complex")
  )
  trait_exp <- derive_traits(exp, trait_fns)

  # Test expr_mat
  expected_expr_mat <- matrix(
    c(1/3, 0.5, 0.5,
      2/3, 1, 0.5),
    nrow = 2, ncol = 3, byrow = TRUE
  )
  rownames(expected_expr_mat) <- c("V1", "V2")
  colnames(expected_expr_mat) <- c("S1", "S2", "S3")
  expect_equal(trait_exp$expr_mat, expected_expr_mat)

  # Test var_info
  expected_var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    trait = c("TFc", "TC")
  )
  expect_equal(trait_exp$var_info, expected_var_info)

  # Test sample_info
  expect_equal(trait_exp$sample_info, exp$sample_info)

  # Test meta_data
  expect_equal(trait_exp$meta_data$exp_type, "traitomics")
})

test_that("derive_traits works for glycoproteomics experiments", {
  # Construct a test experiment
  var_info <- tibble::tibble(
    variable = paste0("V", 1:6),
    protein = rep(c("P1", "P2"), each = 3),
    protein_site = rep(c(10L, 20L), each = 3),
    glycan_structure = rep(glyparse::parse_iupac_condensed(c(
      "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
    )), 2)
  )
  sample_info <- tibble::tibble(sample = paste0("S", 1:3))
  expr_mat <- matrix(
    # Site 1
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0,
    # Site 2
      1, 1, 1,
      1, 1, 1,
      1, 1, 1),
    nrow = 6, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- paste0("V", 1:6)
  colnames(expr_mat) <- paste0("S", 1:3)
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycoproteomics", glycan_type = "N")

  # Calculate derived traits
  trait_fns <- list(
    TFc = prop(nFc > 0),
    TC = prop(T == "complex")
  )
  trait_exp <- derive_traits(exp, trait_fns)

  # Test expr_mat
  expected_expr_mat <- matrix(
    c(1/3, 0.5, 0.5,
      2/3, 1, 0.5,
      1/3, 1/3, 1/3,
      2/3, 2/3, 2/3
    ),
   nrow = 4, ncol = 3, byrow = TRUE
  )
  rownames(expected_expr_mat) <- paste0("V", 1:4)
  colnames(expected_expr_mat) <- paste0("S", 1:3)
  expect_equal(trait_exp$expr_mat, expected_expr_mat)

  # Test var_info
  expected_var_info <- tibble::tibble(
    variable = paste0("V", 1:4),
    protein = c("P1", "P1", "P2", "P2"),
    protein_site = c(10L, 10L, 20L, 20L),
    trait = c("TFc", "TC", "TFc", "TC")
  )
  expect_equal(trait_exp$var_info, expected_var_info)

  # Test sample_info
  expect_equal(trait_exp$sample_info, exp$sample_info)

  # Test meta_data
  expect_equal(trait_exp$meta_data$exp_type, "traitproteomics")
})

test_that("derive_traits works with custom meta-properties", {
  # Construct a test experiment
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3"),
    glycan_structure = glyparse::parse_iupac_condensed(c(
      "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
    ))
  )
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  expr_mat <- matrix(
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N")

  # Create custom meta-property functions
  mp_fns <- list(
    nN = ~ glyrepr::count_mono(.x, "HexNAc"),
    nH = ~ glyrepr::count_mono(.x, "Hex")
  )

  # Calculate derived traits
  trait_fns <- list(
    many_N = prop(nN > 2),
    many_H = prop(nH > 5)
  )
  trait_exp <- derive_traits(exp, trait_fns, mp_fns)

  # Test expr_mat
  expected_expr_mat <- matrix(
    c(2/3, 1, 1/2,
      1/3, 0, 1/2
    ),
    nrow = 2, ncol = 3, byrow = TRUE
  )
  rownames(expected_expr_mat) <- c("V1", "V2")
  colnames(expected_expr_mat) <- c("S1", "S2", "S3")
  expect_equal(trait_exp$expr_mat, expected_expr_mat)

  # Test var_info
  expected_var_info <- tibble::tibble(
    variable = c("V1", "V2"),
    trait = c("many_N", "many_H")
  )
  expect_equal(trait_exp$var_info, expected_var_info)

  # Test sample_info
  expect_equal(trait_exp$sample_info, exp$sample_info)

  # Test meta_data
  expect_equal(trait_exp$meta_data$exp_type, "traitomics")
})

test_that("derive_traits works with default traits", {
  # Construct a test experiment
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3"),
    glycan_structure = glyparse::parse_iupac_condensed(c(
      "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
    ))
  )
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  expr_mat <- matrix(
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N")

  # Calculate derived traits
  trait_exp <- derive_traits(exp)

  # Test var_info
  expect_setequal(trait_exp$var_info$trait, names(all_traits()))
})

test_that("derive_traits keeps glycosite descriptive columns in var_info", {
  strucs <- glyparse::parse_iupac_condensed(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
  ))
  comps <- glyrepr::as_glycan_composition(strucs)
  var_info <- tibble::tibble(
    variable = paste0("V", 1:6),  # index column
    protein = rep(c("P1", "P2"), each = 3),  # keep
    protein_site = rep(c(10L, 20L), each = 3),  # keep
    gene = rep(c("G1", "G2"), each = 3),  # keep
    glycan_structure = rep(strucs, 2),  # remove
    glycan_composition = rep(comps, 2)  # remove
  )
  sample_info <- tibble::tibble(sample = paste0("S", 1:3))
  expr_mat <- matrix(
    # Site 1
    c(1, 0, 1,
      1, 1, 1,
      1, 1, 0,
    # Site 2
      1, 1, 1,
      1, 1, 1,
      1, 1, 1),
    nrow = 6, ncol = 3, byrow = TRUE
  )
  rownames(expr_mat) <- paste0("V", 1:6)
  colnames(expr_mat) <- paste0("S", 1:3)
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycoproteomics", glycan_type = "N")

  # Calculate derived traits
  trait_fns <- list(
    TFc = prop(nFc > 0),
    TC = prop(T == "complex")
  )
  trait_exp <- derive_traits(exp, trait_fns)

  # Test var_info
  expect_setequal(
    colnames(trait_exp$var_info),
    c("variable", "protein", "protein_site", "gene", "trait")
  )
})

test_that("derive_traits raises error for missing columns", {
  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3")
  )
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  expr_mat <- matrix(c(1, 0, 1, 1, 1, 1, 1, 1, 0), nrow = 3, ncol = 3, byrow = TRUE)
  rownames(expr_mat) <- c("V1", "V2", "V3")
  colnames(expr_mat) <- c("S1", "S2", "S3")
  exp <- glyexp::experiment(expr_mat, sample_info, var_info, exp_type = "glycomics", glycan_type = "N")
  expect_error(derive_traits(exp), "Variable information must contain the following columns: glycan_structure.")
})

test_that("derive_traits_() works for glycomics experiments", {
  glycan_structure <- glyparse::parse_iupac_condensed(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
  ))
  tbl <- tibble::tibble(
    sample = rep(c("S1", "S2", "S3"), 3),
    glycan_structure = rep(glycan_structure, each = 3),
    value = c(1, 0, 1, 1, 1, 1, 1, 1, 0)
  )

  trait_fns <- list(
    TFc = prop(nFc > 0),
    TC = prop(T == "complex")
  )
  trait_tbl <- derive_traits_(tbl, "glycomics", trait_fns)

  # Test trait_tbl
  expected <- tibble::tibble(
    trait = rep(c("TFc", "TC"), each = 3),
    sample = rep(c("S1", "S2", "S3"), 2),
    value = c(1/3, 0.5, 0.5, 2/3, 1, 0.5)
  )
  expect_equal(trait_tbl, expected)
})

test_that("derive_traits_() works for glycoproteomics experiments", {
  glycan_structure <- glyparse::parse_iupac_condensed(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
  ))
  tbl <- tibble::tibble(
    sample = rep(c("S1", "S2", "S3"), 6),
    protein = rep(c("P1", "P2"), each = 9),
    protein_site = rep(c(10L, 20L), each = 9),
    glycan_structure = rep(rep(glycan_structure, each = 3), 2),
    value = c(1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  )

  trait_fns <- list(
    TFc = prop(nFc > 0),
    TC = prop(T == "complex")
  )
  trait_tbl <- derive_traits_(tbl, "glycoproteomics", trait_fns)

  # Test trait_tbl
  expected <- tibble::tibble(
    protein = rep(c("P1", "P2"), each = 6),
    protein_site = rep(c(10L, 20L), each = 6),
    trait = rep(rep(c("TFc", "TC"), each = 3), 2),
    sample = rep(c("S1", "S2", "S3"), 4),
    value = c(1/3, 0.5, 0.5, 2/3, 1, 0.5, 1/3, 1/3, 1/3, 2/3, 2/3, 2/3)
  )
  expect_equal(trait_tbl, expected)
})

test_that("derive_traits_() raises error for invalid data type", {
  tbl <- tibble::tibble(
    sample = c("S1", "S2", "S3"),
    variable = c("V1", "V2", "V3"),
    value = c(1, 0, 1)
  )
  expect_error(derive_traits_(tbl, "invalid"), 'must be "glycomics" or "glycoproteomics"')
})

test_that("derive_traits_ ignores other columns for glycomics experiments", {
  glycan_structure <- glyparse::parse_iupac_condensed(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
  ))
  tbl <- tibble::tibble(
    sample = c("S1", "S2", "S3"),
    glycan_structure = glycan_structure,
    value = c(1, 0, 1),
    other = c("A", "B", "C")
  )
  trait_tbl <- derive_traits_(tbl, "glycomics", trait_fns = list(TFc = prop(nFc > 0)))
  expect_equal(colnames(trait_tbl), c("trait", "sample", "value"))
})

test_that("derive_traits_ ignores other columns for glycoproteomics experiments", {
  glycan_structure <- glyparse::parse_iupac_condensed(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
  ))
  tbl <- tibble::tibble(
    sample = c("S1", "S2", "S3"),
    protein = c("P1", "P2", "P3"),
    protein_site = c(10L, 20L, 30L),
    glycan_structure = glycan_structure,
    value = c(1, 0, 1),
    other = c("A", "B", "C")
  )
  trait_tbl <- derive_traits_(tbl, "glycoproteomics", trait_fns = list(TFc = prop(nFc > 0)))
  expect_equal(colnames(trait_tbl), c("protein", "protein_site", "trait", "sample", "value"))
})

test_that("derive_traits_ works with default traits", {
  glycan_structure <- glyparse::parse_iupac_condensed(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # Man9
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",  # core-fuc, bi-antennary
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # bi-antennary
  ))
  tbl <- tibble::tibble(
    sample = rep(c("S1", "S2", "S3"), 3),
    glycan_structure = rep(glycan_structure, each = 3),
    value = c(1, 0, 1, 1, 1, 1, 1, 1, 0)
  )
  trait_tbl <- derive_traits_(tbl, "glycomics")
  expect_setequal(trait_tbl$trait, names(all_traits()))
})

test_that("derive_traits raises error for empty trait_fns", {
  exp <- glyexp::toy_experiment
  expect_error(derive_traits(exp, trait_fns = list()), "must be a non-empty named list or NULL.")
})

test_that("derive_traits_ raises error for empty trait_fns", {
  tbl <- tibble::as_tibble(glyexp::toy_experiment)
  expect_error(derive_traits_(tbl, "glycomics", trait_fns = list()), "must be a non-empty named list or NULL.")
})

test_that("derive_traits raises error when trait_fns has no names", {
  exp <- glyexp::toy_experiment
  expect_error(derive_traits(exp, trait_fns = list(prop(nFc > 0))), "must be a non-empty named list or NULL.")
})

test_that("derive_traits_ raises error when trait_fns has no names", {
  tbl <- tibble::as_tibble(glyexp::toy_experiment)
  expect_error(derive_traits_(tbl, "glycomics", trait_fns = list(prop(nFc > 0))), "must be a non-empty named list or NULL.")
})