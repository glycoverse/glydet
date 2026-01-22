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
  glycans <- c(
    paucimannose_H3N2(),
    paucimannose_H4N2a3()
  )
  expected <- factor(
    c("paucimannose", "paucimannose"),
    levels = c("paucimannose", "hybrid", "highmannose", "complex")
  )
  expect_equal(n_glycan_type(glycans), expected)
})

test_that("n_glycan_type identifies highmannose correctly", {
  glycans <- c(
    highmannose_H5N2(),
    highmannose_H6N2()
  )
  expected <- factor(
    c("highmannose", "highmannose"),
    levels = c("paucimannose", "hybrid", "highmannose", "complex")
  )
  expect_equal(n_glycan_type(glycans), expected)
})

test_that("n_glycan_type identifies hybrid correctly", {
  glycans <- c(
    hybrid_H5N3(),
    hybrid_H4N3a3()
  )
  expected <- factor(
    c("hybrid", "hybrid"),
    levels = c("paucimannose", "hybrid", "highmannose", "complex")
  )
  expect_equal(n_glycan_type(glycans), expected)
})

test_that("n_glycan_type identifies complex correctly", {
  glycans <- c(
    complex_H3N3(),
    complex_H3N4(),
    complex_H4N4()
  )
  expected <- factor(
    c("complex", "complex", "complex"),
    levels = c("paucimannose", "hybrid", "highmannose", "complex")
  )
  expect_equal(n_glycan_type(glycans), expected)
})

# ========== Bisecting N-glycan ==========
test_that("has_bisecting identifies bisecting correctly", {
  glycans <- c(
    complex_H3N5_bisect(),
    complex_H4N4()
  )
  expect_equal(has_bisecting(glycans), c(TRUE, FALSE))
})

# ========== Number of antennae ==========
test_that("n_antennae counts antennae correctly", {
  glycans <- c(
    complex_H3N3(),
    complex_H4N4(),
    complex_H3N5_bisect(),
    complex_H6N5a3()
  )
  expect_equal(n_antennae(glycans), c(1L, 2L, 2L, 3L))
})

test_that("n_antennae counts antenna motifs in non-complex glycans", {
  glycans <- c(
    hybrid_H5N3(),
    highmannose_H5N2(),
    paucimannose_H3N2()
  )
  expect_equal(n_antennae(glycans), c(1L, 0L, 0L))
})

# ========== Number of core fucoses ==========
test_that("n_core_fuc counts core fucoses correctly", {
  glycans <- c(
    complex_H3N3(),
    complex_H3N4F1_coreF(),
    complex_H3N4F2_2coreF(),
    complex_H3N4F1_armF()
  )
  expect_equal(n_core_fuc(glycans), c(0L, 1L, 2L, 0L))
})

# ========== Number of arm fucoses ==========
test_that("n_arm_fuc counts arm fucoses correctly", {
  glycans <- c(
    complex_H3N3(),
    complex_H3N4F1_armF(),
    complex_H3N4F2_2armF(),
    complex_H3N4F1_coreF()
  )
  expect_equal(n_arm_fuc(glycans), c(0L, 1L, 2L, 0L))
})

# ========== Number of galactoses ==========
test_that("n_gal counts galactoses correctly", {
  glycans <- c(
    complex_H4N4S1(),
    highmannose_H5N2()
  )
  expect_equal(n_gal(glycans), c(1L, 0L))
})

# ========== Number of terminal galactoses ==========
test_that("n_terminal_gal counts terminal galactoses correctly", {
  glycans <- c(
    complex_H4N4(),
    complex_H4N4S1(),
    highmannose_H5N2()
  )
  expect_equal(n_terminal_gal(glycans), c(1L, 0L, 0L))
})

# ========== Number of sialic acids ==========
test_that("n_sia counts sialic acids correctly", {
  glycans <- c(
    complex_H4N4S1(),
    complex_H4N4()
  )
  expect_equal(n_sia(glycans), c(1L, 0L))
})

# ========== Number of mannoses ==========
test_that("n_man counts mannoses correctly", {
  glycans <- c(
    paucimannose_H3N2(),
    highmannose_H6N2(),
    hybrid_H5N3(),
    complex_H6N5a3()
  )
  expect_equal(n_man(glycans), c(3L, 6L, 5L, 3L))
})
