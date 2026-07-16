#' Get Meta-Properties of Glycans
#'
#' This function calculates the meta-properties of the given glycans.
#' Meta-properties are properties describing certain structural characteristics of glycans.
#' For example, the number of antennae, the number of core fucoses, etc.
#'
#' @param glycans A `glyrepr::glycan_structure()` vector, or a character vector of glycan structure strings
#'   supported by `glyparse::auto_parse()`.
#' @param mp_fns A named list of meta-property functions.
#'   Names of the list are the names of the meta-properties. Default is [all_mp_fns()].
#'   A meta-property function should takes a `glyrepr::glycan_structure()` vector,
#'   and returns a vector of the meta-property values.
#'   purrr-style lambda functions are supported.
#'
#' @returns A tibble with the meta-properties. Column names are the names of the meta-properties.
#' @seealso [all_mp_fns()]
#'
#' @examples
#' glycans <- c(
#'   "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
#'   "Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",
#'   "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
#' )
#'
#' # Use default meta-property functions
#' get_meta_properties(glycans)
#'
#' # Use custom meta-property functions
#' fns <- list(
#'   nN = ~ glyrepr::count_mono(.x, "HexNAc"),  # purrr-style lambda function
#'   nH = ~ glyrepr::count_mono(.x, "Hex")
#' )
#' get_meta_properties(glycans, fns)
#'
#' @export
get_meta_properties <- function(glycans, mp_fns = NULL) {
  glycans <- .process_glycans(glycans)
  if (is.null(mp_fns)) {
    mp_fns <- all_mp_fns()
  } else {
    mp_fns <- purrr::map(mp_fns, rlang::as_function)
    checkmate::assert_named(mp_fns)
  }
  tibble::as_tibble(purrr::map(mp_fns, ~ .x(glycans)))
}

#' Add Meta-Properties to a Glyco SummarizedExperiment
#'
#' This function adds meta-properties to the variable information of a
#' [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()].
#' Under the hood, it uses [get_meta_properties()] to calculate the meta-properties
#' on the "glycan_structure" column (or column specified by `struc_col`) of the variable information tibble,
#' and then adds the result back as new columns.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object.
#' @param mp_fns A named list of meta-property functions.
#'   Names of the list are the names of the meta-properties. Default is [all_mp_fns()].
#'   A meta-property function should takes a `glyrepr::glycan_structure()` vector,
#'   and returns a vector of the meta-property values.
#'   purrr-style lambda functions are supported.
#' @param struc_col The column name of the glycan structures in the variable information tibble.
#'   Default is "glycan_structure".
#' @param overwrite Whether to overwrite the existing meta-property columns.
#'   Default is FALSE, raising an error if the existing columns are found.
#'
#' @return The input data container with meta-properties added to its variable
#'   information. The input container type is preserved.
#'
#' @examples
#' library(glyexp)
#' library(SummarizedExperiment)
#'
#' # Compare rowData columns before and after adding meta-properties
#' gp_se <- real_experiment |>
#'   slice_sample_row(n = 10)
#' colnames(rowData(gp_se))
#'
#' gp_se2 <- add_meta_properties(gp_se)
#' colnames(rowData(gp_se2))
#'
#' @seealso [get_meta_properties()], [glyexp::GlycomicSE()],
#'   [glyexp::GlycoproteomicSE()]
#'
#' @export
add_meta_properties <- function(
  exp,
  mp_fns = NULL,
  struc_col = "glycan_structure",
  overwrite = FALSE
) {
  .assert_data_container(exp)
  legacy <- inherits(exp, "glyexp_experiment")
  exp <- .as_glyco_se(exp)
  checkmate::assert_string(struc_col)
  checkmate::assert_flag(overwrite)
  .check_var_info_cols(exp, struc_col)

  var_info <- .get_var_info(exp)
  if (is.null(mp_fns)) {
    mp_names <- names(all_mp_fns())
  } else {
    mp_names <- names(mp_fns)
  }
  if (overwrite) {
    # Remove the existing meta-property columns if any
    var_info <- dplyr::select(var_info, -dplyr::any_of(mp_names))
  } else {
    # Check if any existing columns are the same as the new meta-property names
    exist_mp_names <- intersect(mp_names, colnames(var_info))
    if (length(exist_mp_names) > 0) {
      cli::cli_abort(c(
        "Variable information tibble must not contain columns with the same names as the meta-properties.",
        "x" = "The following columns already exist: {.field {exist_mp_names}}."
      ))
    }
  }

  meta_properties <- get_meta_properties(var_info[[struc_col]], mp_fns)
  exp <- .set_var_info(
    exp,
    dplyr::bind_cols(var_info, tibble::as_tibble(meta_properties))
  )
  .restore_data_container(exp, legacy)
}
