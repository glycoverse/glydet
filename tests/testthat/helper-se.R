as_test_se <- function(exp) {
  if (methods::is(exp, "SummarizedExperiment")) {
    return(exp)
  }

  switch(
    glyexp::get_exp_type(exp),
    glycomics = glyexp::as_glycomic_se(exp),
    glycoproteomics = glyexp::as_glycoproteomic_se(exp)
  )
}

test_expr_mat <- function(exp) {
  SummarizedExperiment::assay(exp, 1)
}

test_var_info <- function(exp) {
  row_data <- SummarizedExperiment::rowData(exp)
  tibble::as_tibble(as.list(row_data), .name_repair = "minimal") |>
    tibble::add_column(variable = rownames(exp), .before = 1)
}

test_sample_info <- function(exp) {
  tibble::as_tibble(SummarizedExperiment::colData(exp), rownames = "sample")
}

test_meta_data <- function(exp) {
  S4Vectors::metadata(exp)
}
