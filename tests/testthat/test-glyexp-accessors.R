test_that("package code handles data through SummarizedExperiment natively", {
  r_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"))
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  source_lines <- purrr::map(r_files, readLines, warn = FALSE)
  source_lines <- purrr::map(
    source_lines,
    ~ .x[!stringr::str_detect(.x, "^\\s*#")]
  )
  names(source_lines) <- r_files

  legacy_container_access <- purrr::imap(
    source_lines,
    ~ grep(
      paste0(
        "glyexp::(get_expr_mat|get_var_info|get_sample_info|get_meta_data|",
        "mutate_var|select_var|left_join_var|standardize_variable)"
      ),
      .x,
      value = TRUE
    )
  )
  legacy_container_access <- purrr::discard(
    legacy_container_access,
    ~ length(.x) == 0
  )

  expect_length(legacy_container_access, 0)
})
