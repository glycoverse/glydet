test_that("package code uses glyexp accessors instead of experiment internals", {
  r_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"))
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  source_lines <- purrr::map(r_files, readLines, warn = FALSE)
  names(source_lines) <- r_files

  direct_internal_access <- purrr::imap(
    source_lines,
    ~ grep(
      "[$](expr_mat|sample_info|var_info|meta_data)",
      .x,
      value = TRUE
    )
  )
  direct_internal_access <- purrr::discard(
    direct_internal_access,
    ~ length(.x) == 0
  )

  expect_length(direct_internal_access, 0)
})
