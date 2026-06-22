#' Get GlycoMotif database metadata
#'
#' @returns A tibble returned by [glymotif::db_motif_info()].
#' @noRd
.glymotif_db_motif_info <- memoise::memoise(function() {
  glymotif::db_motif_info()
})

#' Check GlycoMotif database motif names
#'
#' @param motif_names A character vector of motif names.
#'
#' @returns A scalar logical indicating whether all names are known database motifs.
#' @noRd
.are_db_motif_names <- function(motif_names) {
  checkmate::assert_character(motif_names, any.missing = FALSE)

  if ("db_motif_info" %in% getNamespaceExports("glymotif")) {
    all(motif_names %in% .glymotif_db_motif_info()$name)
  } else {
    all(glymotif::is_known_motif(motif_names))
  }
}

#' Get GlycoMotif database motif structures by name
#'
#' @param motif_names A character vector of motif names.
#'
#' @returns A `glyrepr::glycan_structure()` vector.
#' @noRd
.get_db_motif_structure <- function(motif_names) {
  checkmate::assert_character(motif_names, any.missing = FALSE)

  if ("db_motif_info" %in% getNamespaceExports("glymotif")) {
    motif_info <- .glymotif_db_motif_info()
    motif_idx <- match(motif_names, motif_info$name)
    missing_motifs <- is.na(motif_idx)

    if (any(missing_motifs)) {
      cli::cli_abort(c(
        "Unknown GlycoMotif database motif name{?s}: {.val {motif_names[missing_motifs]}}.",
        "i" = "Use `glymotif::db_motif_info()` to inspect supported motif names."
      ))
    }

    motif_info$glycan_structure[motif_idx]
  } else {
    glymotif::get_motif_structure(motif_names)
  }
}
