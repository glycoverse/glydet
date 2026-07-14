test_that("container functions accept GlycomicSE natively", {
  se <- glyexp::real_experiment2 |>
    glyexp::slice_head_var(n = 10) |>
    glyexp::slice_head_obs(n = 5) |>
    glyexp::as_glycomic_se()
  motif <- c(
    observed = as.character(
      SummarizedExperiment::rowData(se)$glycan_structure[1]
    )
  )

  with_mps <- add_meta_properties(se)
  traits <- derive_traits(
    se,
    trait_fns = list(TS = prop(nS > 0))
  )
  motifs <- quantify_motifs(
    se,
    motifs = motif
  )

  expect_s4_class(with_mps, "GlycomicSE")
  expect_true("nS" %in% colnames(SummarizedExperiment::rowData(with_mps)))
  expect_s4_class(traits, "SummarizedExperiment")
  expect_identical(S4Vectors::metadata(traits)$exp_type, "traitomics")
  expect_identical(SummarizedExperiment::rowData(traits)$trait, "TS")
  expect_s4_class(motifs, "SummarizedExperiment")
  expect_identical(S4Vectors::metadata(motifs)$exp_type, "traitomics")
  expect_identical(
    SummarizedExperiment::rowData(motifs)$trait,
    "observed"
  )
})

test_that("container functions accept GlycoproteomicSE natively", {
  se <- glyexp::real_experiment |>
    glyexp::slice_head_var(n = 10) |>
    glyexp::slice_head_obs(n = 5) |>
    glyexp::as_glycoproteomic_se()
  motif <- c(
    observed = as.character(
      SummarizedExperiment::rowData(se)$glycan_structure[1]
    )
  )

  with_mps <- add_meta_properties(se)
  traits <- derive_traits(
    se,
    trait_fns = list(TS = prop(nS > 0))
  )
  motifs <- quantify_motifs(
    se,
    motifs = motif
  )

  expect_s4_class(with_mps, "GlycoproteomicSE")
  expect_true("nS" %in% colnames(SummarizedExperiment::rowData(with_mps)))
  expect_s4_class(traits, "SummarizedExperiment")
  expect_identical(S4Vectors::metadata(traits)$exp_type, "traitproteomics")
  expect_true(all(
    c("protein", "protein_site", "trait") %in%
      colnames(SummarizedExperiment::rowData(traits))
  ))
  expect_s4_class(motifs, "SummarizedExperiment")
  expect_identical(S4Vectors::metadata(motifs)$exp_type, "traitproteomics")
})

test_that("legacy experiment support preserves legacy return types", {
  se <- as_test_se(glyexp::real_experiment2)[seq_len(10), seq_len(5)]
  motif <- c(
    observed = as.character(SummarizedExperiment::rowData(se)$glycan_structure[
      1
    ])
  )
  exp <- suppressWarnings(glyexp::from_se(se))

  expect_s3_class(add_meta_properties(exp), "glyexp_experiment")
  expect_s3_class(
    derive_traits(exp, trait_fns = list(TS = prop(nS > 0))),
    "glyexp_experiment"
  )
  expect_s3_class(
    quantify_motifs(exp, motifs = motif),
    "glyexp_experiment"
  )
})

test_that("derive_traits_ is removed from the public API", {
  expect_false("derive_traits_" %in% getNamespaceExports("glydet"))
})
