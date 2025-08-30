# Test for N-glycan properties functions in glydet-mp.R

# Test glycan structures
paucimannose_H3N2 <- function() {
  glyparse::parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

paucimannose_H4N2a3 <- function() {
  glyparse::parse_iupac_condensed("Man(a1-3)Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

highmannose_H5N2 <- function() {
  glyparse::parse_iupac_condensed("Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

highmannose_H6N2 <- function() {
  glyparse::parse_iupac_condensed("Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

hybrid_H5N3 <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

hybrid_H4N3a3 <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H3N3 <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H3N4 <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H4N4 <- function() {
  glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H4N4S1 <- function() {
  glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H3N5_bisect <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H6N5a3 <- function() {
  glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H3N4F1_coreF <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-")
}

complex_H3N4F1_armF <- function() {
  glyparse::parse_iupac_condensed("Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

complex_H3N4F2_2coreF <- function() {
  glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)][Fuc(a1-8)]GlcNAc(?1-")
}

complex_H3N4F2_2armF <- function() {
  glyparse::parse_iupac_condensed("Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
}

o_glycan_core_1 <- function() {
  glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
}

# ========== N-glycan types ==========
test_that("n_glycan_type identifies paucimannose correctly", {
  expect_equal(n_glycan_type(paucimannose_H3N2()), factor("paucimannose", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
  expect_equal(n_glycan_type(paucimannose_H4N2a3()), factor("paucimannose", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
})

test_that("n_glycan_type identifies highmannose correctly", {
  expect_equal(n_glycan_type(highmannose_H5N2()), factor("highmannose", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
  expect_equal(n_glycan_type(highmannose_H6N2()), factor("highmannose", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
})

test_that("n_glycan_type identifies hybrid correctly", {
  expect_equal(n_glycan_type(hybrid_H5N3()), factor("hybrid", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
  expect_equal(n_glycan_type(hybrid_H4N3a3()), factor("hybrid", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
})

test_that("n_glycan_type identifies complex correctly", {
  expect_equal(n_glycan_type(complex_H3N3()), factor("complex", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
  expect_equal(n_glycan_type(complex_H3N4()), factor("complex", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
  expect_equal(n_glycan_type(complex_H4N4()), factor("complex", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
})

test_that("n_glycan_type works vectorizedly", {
  glycans <- c(
    paucimannose_H3N2(),
    complex_H4N4()
  )
  expected <- factor(c("paucimannose", "complex"), levels = c("paucimannose", "hybrid", "highmannose", "complex"))
  expect_equal(n_glycan_type(glycans), expected)
})

# ========== Bisecting N-glycan ==========
test_that("has_bisecting identifies bisecting correctly", {
  expect_true(has_bisecting(complex_H3N5_bisect()))
  expect_false(has_bisecting(complex_H4N4()))
})

test_that("has_bisecting works vectorizedly", {
  glycans <- c(
    complex_H3N5_bisect(),
    complex_H4N4()
  )
  expect_equal(has_bisecting(glycans), c(TRUE, FALSE))
})

# ========== Number of antennae ==========
test_that("n_antennae counts antennae correctly", {
  expect_identical(n_antennae(complex_H3N3()), 1L)
  expect_identical(n_antennae(complex_H4N4()), 2L)
  expect_identical(n_antennae(complex_H3N5_bisect()), 2L) # bisecting doesn't count as antenna
  expect_identical(n_antennae(complex_H6N5a3()), 3L)
})

test_that("n_antennae counts antenna motifs in non-complex glycans", {
  # In the simplified version, we directly count antenna motifs regardless of glycan type
  expect_identical(n_antennae(hybrid_H5N3()), 1L)  # Has one antenna branch
  expect_identical(n_antennae(highmannose_H5N2()), 0L)  # No antenna branches
  expect_identical(n_antennae(paucimannose_H3N2()), 0L)  # No antenna branches
})

test_that("n_antennae works vectorizedly", {
  glycans <- c(
    complex_H3N3(),
    complex_H4N4()
  )
  expect_equal(n_antennae(glycans), c(1L, 2L))
})

# ========== Number of core fucoses ==========
test_that("n_core_fuc counts core fucoses correctly", {
  expect_identical(n_core_fuc(complex_H3N3()), 0L)
  expect_identical(n_core_fuc(complex_H3N4F1_coreF()), 1L)
  expect_identical(n_core_fuc(complex_H3N4F2_2coreF()), 2L)
  expect_identical(n_core_fuc(complex_H3N4F1_armF()), 0L) # arm fucose, not core
})

test_that("n_core_fuc works vectorizedly", {
  glycans <- c(
    complex_H3N3(),
    complex_H3N4F1_coreF()
  )
  expect_equal(n_core_fuc(glycans), c(0L, 1L))
})

# ========== Number of arm fucoses ==========
test_that("n_arm_fuc counts arm fucoses correctly", {
  expect_identical(n_arm_fuc(complex_H3N3()), 0L)
  expect_identical(n_arm_fuc(complex_H3N4F1_armF()), 1L)
  expect_identical(n_arm_fuc(complex_H3N4F2_2armF()), 2L)
  expect_identical(n_arm_fuc(complex_H3N4F1_coreF()), 0L) # core fucose, not arm
})

test_that("n_arm_fuc works vectorizedly", {
  glycans <- c(
    complex_H3N3(),
    complex_H3N4F1_armF()
  )
  expect_equal(n_arm_fuc(glycans), c(0L, 1L))
})

# ========== Number of galactoses ==========
test_that("n_gal counts galactoses correctly", {
  expect_identical(n_gal(complex_H4N4S1()), 1L)
  expect_identical(n_gal(highmannose_H5N2()), 0L)
})

test_that("n_gal works vectorizedly", {
  glycans <- c(
    complex_H4N4S1(),
    highmannose_H5N2()
  )
  expect_equal(n_gal(glycans), c(1L, 0L))
})

# ========== Number of terminal galactoses ==========
test_that("n_terminal_gal counts terminal galactoses correctly", {
  expect_identical(n_terminal_gal(complex_H4N4()), 1L)
  expect_identical(n_terminal_gal(complex_H4N4S1()), 0L) # capped with sialic acid
  expect_identical(n_terminal_gal(highmannose_H5N2()), 0L)
})

test_that("n_terminal_gal works vectorizedly", {
  glycans <- c(
    complex_H4N4(),
    complex_H4N4S1()
  )
  expect_equal(n_terminal_gal(glycans), c(1L, 0L))
})

# ========== Non-N-glycan handling ==========
test_that("functions handle non-N-glycans properly", {
  # The new version should throw an error for non-N-glycans
  expect_error(n_glycan_type(o_glycan_core_1()), "not N-glycans")
  expect_error(has_bisecting(o_glycan_core_1()), "not N-glycans")
  expect_error(n_antennae(o_glycan_core_1()), "not N-glycans")
  expect_error(n_core_fuc(o_glycan_core_1()), "not N-glycans")
  expect_error(n_arm_fuc(o_glycan_core_1()), "not N-glycans")
  expect_error(n_gal(o_glycan_core_1()), "not N-glycans")
  expect_error(n_terminal_gal(o_glycan_core_1()), "not N-glycans")
})

# ========== Multi-format Support ==========
test_that("N-glycan functions support multiple structure formats", {
  # Test IUPAC-condensed format
  glycan_iupac <- "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  expect_equal(n_glycan_type(glycan_iupac), factor("paucimannose", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
})

# ========== Processing and Simplification ==========
test_that("processing handles linkages, generics, and substituents correctly", {
  # Test with a complex glycan that has linkages and potential substituents
  glycan_complex <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  
  # Should work fine and be processed as generic
  expect_equal(n_glycan_type(glycan_complex), factor("complex", levels = c("paucimannose", "hybrid", "highmannose", "complex")))
  expect_identical(n_antennae(glycan_complex), 2L)
  expect_identical(n_gal(glycan_complex), 1L)
})
