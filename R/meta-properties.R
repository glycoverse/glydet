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
get_meta_properties <- function(glycans, mp_fns = all_mp_fns()) {
  glycans <- .process_glycans(glycans)
  mp_fns <- purrr::map(mp_fns, rlang::as_function)
  tibble::as_tibble(purrr::map(mp_fns, ~ .x(glycans)))
}