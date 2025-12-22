.process_glycans <- function(glycans) {
  if (is.character(glycans)) {
    glycans <- glyparse::auto_parse(glycans)
  } else if (!glyrepr::is_glycan_structure(glycans)) {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} object or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(glycans)}}."
    ))
  }
  glycans
}

.check_var_info_cols <- function(exp, cols) {
  has_cols <- cols %in% colnames(exp$var_info)
  if (!all(has_cols)) {
    missing_cols <- cols[!has_cols]
    cli::cli_abort(c(
      "Variable information must contain the following columns: {.field {cols}}.",
      "x" = "Cannot find {.field {missing_cols}} in {.field var_info}."
    ))
  }
}

.get_api_key <- function() {
  api_key <- Sys.getenv("DEEPSEEK_API_KEY")
  if (api_key == "") {
    cli::cli_abort(c(
      "API key for DeepSeek chat model is not set.",
      "i" = "Please set the environment variable `DEEPSEEK_API_KEY` to your API key.",
      "i" = "You can obtain an API key from https://platform.deepseek.com."
    ))
  }
  api_key
}