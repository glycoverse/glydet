scale_test_exp <- function(proteomics = FALSE, legacy = FALSE) {
  structures <- glyparse::parse_iupac_condensed(c("Man(b1-4)GlcNAc", "GlcNAc"))
  info <- tibble::tibble(
    variable = c("V1", "V2"), glycan_structure = structures,
    glycan_composition = glyrepr::as_glycan_composition(structures),
    selected = c(TRUE, FALSE), count = c(2, 4)
  )
  if (proteomics) {
    info$protein <- "P1"
    info$protein_site <- 10L
  }
  exp <- glyexp::experiment(
    matrix(c(1, 3, 2, 6), 2, dimnames = list(c("V1", "V2"), c("S1", "S2"))),
    tibble::tibble(sample = c("S1", "S2"), group = factor(c("A", "B"))),
    info, exp_type = if (proteomics) "glycoproteomics" else "glycomics",
    glycan_type = "N"
  )
  exp$meta_data$scale_test_note <- "preserved"
  if (legacy) exp else as_test_se(exp)
}

set_test_scale <- function(exp, scale) {
  if (inherits(exp, "glyexp_experiment")) {
    exp$meta_data["scale"] <- list(scale)
  } else {
    S4Vectors::metadata(exp)["scale"] <- list(scale)
  }
  exp
}

test_that("nonlinear and invalid scales are rejected before any computation", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  for (proteomics in c(FALSE, TRUE)) {
    for (legacy in c(FALSE, TRUE)) {
      exp <- scale_test_exp(proteomics, legacy)
      for (scale in list("log", "log2", "log10", "logx", "unknown", NA_character_,
                         NULL, character(), c("raw", "log"), 1)) {
        input <- set_test_scale(exp, scale)
        expect_error(derive_traits(input, trait_fns = list(p = prop(selected))),
                     "require linear abundances")
        # Invalid motif ensures scale validation happens before motif work.
        expect_error(quantify_motifs(input, motifs = NULL), "require linear abundances")
      }
    }
  }
})

test_that("raw mixed traits have explicit unlogged output and preserve sample metadata", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  traits <- list(p = prop(selected), r = ratio(selected, !selected),
                 m = wmean(count), s = wsum(count), t = total(count > 0))
  for (proteomics in c(FALSE, TRUE)) {
    for (legacy in c(FALSE, TRUE)) {
      input <- set_test_scale(scale_test_exp(proteomics, legacy), "raw")
      result <- derive_traits(input, trait_fns = traits)
      se <- if (legacy) glyexp::as_se(result) else result
      expect_equal(unname(SummarizedExperiment::assay(se)),
                   matrix(c(.25, 1 / 3, 3.5, 14, 4, .25, 1 / 3, 3.5, 28, 8), 5))
      expect_identical(S4Vectors::metadata(se)$scale, "raw")
      expect_identical(S4Vectors::metadata(se)$scale_test_note, "preserved")
      before <- if (legacy) glyexp::as_se(input) else input
      expect_identical(SummarizedExperiment::colData(se), SummarizedExperiment::colData(before))
      expect_identical(colnames(se), colnames(before))
      expect_identical(inherits(result, "glyexp_experiment"), legacy)
      if (proteomics) {
        expect_equal(SummarizedExperiment::rowData(se)$protein, rep("P1", 5))
        expect_equal(SummarizedExperiment::rowData(se)$protein_site, rep(10L, 5))
      }
    }
  }
})

test_that("absent scale is backward compatible but output scale is explicit", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  for (proteomics in c(FALSE, TRUE)) {
    for (legacy in c(FALSE, TRUE)) {
      exp <- scale_test_exp(proteomics, legacy)
      if (legacy) exp$meta_data$scale <- NULL else S4Vectors::metadata(exp)$scale <- NULL
      results <- list(
        derive_traits(exp, list(p = prop(selected))),
        quantify_motifs(exp, c(mannose = "Man"), method = "relative"),
        quantify_motifs(exp, c(mannose = "Man"), method = "absolute")
      )
      for (i in seq_along(results)) {
        se <- if (legacy) glyexp::as_se(results[[i]]) else results[[i]]
        expect_identical(S4Vectors::metadata(se)$scale, "raw")
        expect_identical(S4Vectors::metadata(se)$scale_test_note, "preserved")
        expected <- if (i == 3L) c(1, 2) else c(.25, .25)
        expect_equal(as.numeric(SummarizedExperiment::assay(se)), expected)
      }
    }
  }
})

test_that("the reported log2 proportion bug errors instead of returning zero", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  exp <- set_test_scale(scale_test_exp(), "log2")
  SummarizedExperiment::assay(exp) <- log2(SummarizedExperiment::assay(exp))
  expect_error(derive_traits(exp, list(p = prop(selected))), "require linear abundances")
  # The factory itself has no scale metadata: passing logged data directly is invalid.
  expect_equal(unname(prop(selected)(matrix(log2(c(1, 3)), 2),
                                     tibble::tibble(selected = c(TRUE, FALSE)))), 0)
})

test_that("relative and absolute motifs return unlogged values for both experiment types", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  for (proteomics in c(FALSE, TRUE)) {
    for (legacy in c(FALSE, TRUE)) {
      exp <- set_test_scale(scale_test_exp(proteomics, legacy), "raw")
      for (method in c("relative", "absolute")) {
        result <- quantify_motifs(exp, c(mannose = "Man"), method = method)
        se <- if (legacy) glyexp::as_se(result) else result
        expect_identical(S4Vectors::metadata(se)$scale, "raw")
        expected <- if (method == "relative") c(.25, .25) else c(1, 2)
        expect_equal(as.numeric(SummarizedExperiment::assay(se)), expected)
      }
    }
  }
})
