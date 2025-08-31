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

#' Add Meta-Properties to Experiment
#'
#' This function adds meta-properties to the variable information of a [glyexp::experiment()].
#' Under the hood, it uses [get_meta_properties()] to calculate the meta-properties
#' on the "glycan_structure" column (or column specified by `struc_col`) of the variable information tibble,
#' and then adds the result back as new columns.
#'
#' @param exp An [glyexp::experiment()] object.
#' @param mp_fns A named list of meta-property functions.
#'   Names of the list are the names of the meta-properties. Default is [all_mp_fns()].
#'   A meta-property function should takes a `glyrepr::glycan_structure()` vector,
#'   and returns a vector of the meta-property values.
#'   purrr-style lambda functions are supported.
#' @param struc_col The column name of the glycan structures in the variable information tibble.
#'   Default is "glycan_structure".
#'
#' @return An [glyexp::experiment()] object with meta-properties added to the variable information.
#' @seealso [get_meta_properties()], [glyexp::experiment()]
#'
#' @export
add_meta_properties <- function(exp, mp_fns = NULL, struc_col = "glycan_structure") {
  checkmate::assert_class(exp, "glyexp_experiment")
  checkmate::assert_string(struc_col)

  if (!struc_col %in% colnames(exp$var_info)) {
    cli::cli_abort(c(
      "Column {.field {struc_col}} not found in {.field var_info}.",
      "i" = "Got {.col {colnames(exp$var_info)}}."
    ))
  }

  meta_properties <- get_meta_properties(exp$var_info[[struc_col]], mp_fns)
  exp$var_info <- dplyr::bind_cols(exp$var_info, meta_properties)
  exp
}