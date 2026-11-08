test_that("SE annotations never supply or replace axis identifiers", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  for (proteomics in c(FALSE, TRUE)) {
    for (values in list(c("V1", "V2"), c("annotation", "annotation"))) {
      exp <- scale_test_exp(proteomics)
      SummarizedExperiment::rowData(exp)$variable <- values
      SummarizedExperiment::colData(exp)$sample <- c("group", "group")
      before <- SummarizedExperiment::assay(exp, withDimnames = FALSE)
      result <- add_meta_properties(exp, list(n = ~ rep(1L, length(.x))))
      expect_identical(SummarizedExperiment::rowData(result)$variable, values)
      expect_identical(
        SummarizedExperiment::colData(result),
        SummarizedExperiment::colData(exp)
      )
      expect_identical(
        SummarizedExperiment::assay(result, withDimnames = FALSE),
        before
      )
      expect_identical(rownames(result), rownames(exp))
      expect_true(methods::validObject(result))
      if (proteomics && length(unique(values)) == 1L) {
        derived <- derive_traits(exp, list(p = prop(selected)))
        expect_identical(
          SummarizedExperiment::rowData(derived)$variable,
          "annotation"
        )
        expect_identical(rownames(derived), "P1-10-p")
      }
    }
  }
})

test_that("metadata enrichment repairs only absent scale through public conversion", {
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)
  for (proteomics in c(FALSE, TRUE)) {
    exp <- scale_test_exp(proteomics)
    S4Vectors::metadata(exp)$scale <- NULL
    result <- add_meta_properties(exp, list(n = ~ rep(1L, length(.x))))
    expect_identical(glyexp::get_scale(result), "raw")
    expect_null(glyexp::get_scale(exp))
    expect_identical(
      SummarizedExperiment::assay(result),
      SummarizedExperiment::assay(exp)
    )
    expect_true(methods::validObject(result))
    for (scale in c("log2", "log10", "logx")) {
      logged <- set_test_scale(exp, scale)
      result <- add_meta_properties(logged, list(n = ~ rep(1L, length(.x))))
      expect_identical(glyexp::get_scale(result), scale)
    }
  }
})
