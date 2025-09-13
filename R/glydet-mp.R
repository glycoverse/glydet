# This file contains built-in meta-property functions.
# These functions all work with generic glycans. (e.g. "Hex", "HexNAc")


#' Get All Meta-Property Functions
#'
#' This function returns a named list of all meta-property functions:
#' - "T": `n_glycan_type()`: type of the glycan
#' - "B": `has_bisecting()`: whether the glycan has a bisecting GlcNAc
#' - "nA": `n_antennae()`: number of antennae
#' - "nFc": `n_core_fuc()`: number of core fucoses
#' - "nFa": `n_arm_fuc()`: number of arm fucoses
#' - "nG": `n_gal()`: number of galactoses
#' - "nGt": `n_terminal_gal()`: number of terminal galactoses
#' - "nS": `n_sia()`: number of sialic acids
#' - "nM": `n_man()`: number of mannoses
#'
#' @returns
#' A named list of meta-property functions.
#'
#' @export
all_mp_fns <- function() {
  list(
    "T" = n_glycan_type,
    "B" = has_bisecting,
    "nA" = n_antennae,
    "nFc" = n_core_fuc,
    "nFa" = n_arm_fuc,
    "nG" = n_gal,
    "nGt" = n_terminal_gal,
    "nS" = n_sia,
    "nM" = n_man
  )
}

#' Determine N-Glycan Key Properties
#'
#' @description
#' These functions check key properties of an N-glycan:
#' - `n_glycan_type()`: Determine the N-glycan type.
#' - `has_bisecting()`: Check if the glycan has a bisecting GlcNAc.
#' - `n_antennae()`: Count the number of antennae.
#' - `n_core_fuc()`: Count the number of core fucoses.
#' - `n_arm_fuc()`: Count the number of arm fucoses.
#' - `n_gal()`: Count the number of galactoses.
#' - `n_terminal_gal()`: Count the number of terminal galactoses.
#' - `n_sia()`: Count the number of sialic acids.
#' - `n_man()`: Count the number of mannoses.
#'
#' All functions assume the glycans are N-glycans without validation,
#' thus may return meaningless values for non-N-glycans.
#' Therefore, please make sure to pass in N-glycans only.
#'
#' All functions put minimum requirement on the glycans,
#' i.e. they work with glycans with generic monosaccharides (e.g. "Hex", "HexNAc")
#' and no linkage information.
#' This type of structures are common in glycoproteomics and glycomics studies.
#'
#' @details
#' # `n_glycan_type()`: N-Glycan Types
#'
#' Four types of N-glycans are recognized: high mannose, hybrid, complex, and paucimannose.
#' For more information about N-glycan types,
#' see [Essentials of Glycobiology](https://www.ncbi.nlm.nih.gov/books/NBK579964/#_s9_2_).
#'
#' # `has_bisecting()`: Bisecting GlcNAc
#'
#' Bisecting GlcNAc is a GlcNAc residue attached to the core mannose of N-glycans.
#' ```
#'      Man
#'         \
#' GlcNAc - Man - GlcNAc - GlcNAc -
#' ~~~~~~  /
#'      Man
#' ```
#'
#' # `n_antennae()`: Number of Antennae
#'
#' The number of antennae is the number of branching GlcNAc to the core mannoses.
#'
#' # `n_core_fuc()`: Number of Core Fucoses
#'
#' Core fucoses are those fucose residues attached to the core GlcNAc of an N-glycan.
#' ```
#' Man             Fuc  <- core fucose
#'    \             |
#'     Man - GlcNAc - GlcNAc -
#'    /
#' Man
#' ```
#'
#' # `n_arm_fuc()`: Number of Arm Fucoses
#'
#' Arm focuses are those focuse residues attached to the branching GlcNAc
#' of an N-glycan.
#' ```
#'  Fuc  <- arm fucose
#'   |
#' GlcNAc - Man
#'             \
#'              Man - GlcNAc - GlcNAc -
#'             /
#' GlcNAc - Man
#' ```
#'
#' # `n_gal()`: Number of Galactoses
#'
#' This function seems useless and silly.
#' It is, if you have a well-structured glycan with concrete monosaccharides.
#' However, if you only have "Hex" or "H" at hand,
#' it is tricky to know how many of them are "Gal" and how many are "Man".
#' This function makes a simply assumption that all the rightmost "H" in a
#' "H-H-N-H" unit is a galactose.
#' The two "H" on the left are mannoses of the N-glycan core.
#' The "N" is a GlcNAc attached to one core mannose.
#'
#' # `n_terminal_gal()`: Number of Terminal Galactoses
#'
#' Terminal galactoses are those galactose residues on the non-reducing end
#' without sialic acid capping.
#' ```
#'          Gal - GlcNAc - Man
#'          ~~~               \
#'      terminal Gal           Man - GlcNAc - GlcNAc -
#'                            /
#' Neu5Ac - Gal - GlcNAc - Man
#'          ~~~
#'    not terminal Gal
#' ```
#'
#' # `n_sia()`: Number of Sialic Acids
#'
#' Number of sialic acids (Neu5Ac). Neu5Gc is not counted.
#'
#' # `n_man()`: Number of Mannoses
#'
#' Number of mannoses. This function assumes the Hex of the N-glycan core is mannoses.
#' Also, for high-mannose and paucimannose glycans, all Hex are mannoses.
#' Finally, for hybrid glycans, all the rightmost (the side without branching HexNAc) are mannoses.
#'
#' @param glycans A `glyrepr::glycan_structure()` vector.
#'
#' @returns
#' - `n_glycan_type()`: A factor vector indicating the N-glycan type,
#'   either "highmannose", "hybrid", "complex", or "paucimannose".
#' - `has_bisecting()`: A logical vector indicating if the glycan has a bisecting GlcNAc.
#' - `n_antennae()`: An integer vector indicating the number of antennae.
#' - `n_core_fuc()`: An integer vector indicating the number of core fucoses.
#' - `n_arm_fuc()`: An integer vector indicating the number of arm fucoses.
#' - `n_gal()`: An integer vector indicating the number of galactoses.
#' - `n_terminal_gal()`: An integer vector indicating the number of terminal galactoses.
#' - `n_sia()`: An integer vector indicating the number of sialic acids.
#' - `n_man()`: An integer vector indicating the number of mannoses.
#'
#' @export
n_glycan_type <- function(glycans) {
  motif_names <- c("core", "pauciman", "hybrid", "highman")
  motifs <- purrr::map(motif_names, .get_n_glycan_motif)
  motifs <- do.call(c, motifs)
  alignments <- c("whole", "whole", "core", "core")
  have_motif_mat <- glymotif::have_motifs(glycans, motifs, alignments, ignore_linkages = TRUE)
  colnames(have_motif_mat) <- motif_names
  res <- dplyr::case_when(
    have_motif_mat[, "core"] ~ "paucimannose",
    have_motif_mat[, "pauciman"] ~ "paucimannose",
    have_motif_mat[, "hybrid"] ~ "hybrid",
    have_motif_mat[, "highman"] ~ "highmannose",
    TRUE ~ "complex"
  )
  res <- factor(res, levels = c("paucimannose", "hybrid", "highmannose", "complex"))
  res
}

#' @rdname n_glycan_type
#' @export
has_bisecting <- function(glycans) {
  bisect_motif <- .get_n_glycan_motif("bisect")
  glymotif::have_motif(glycans, bisect_motif, alignment = "core", ignore_linkages = TRUE)
}

#' @rdname n_glycan_type
#' @export
n_antennae <- function(glycans) {
  antenna_motif <- .get_n_glycan_motif("antenna")
  glymotif::count_motif(glycans, antenna_motif, alignment = "core", ignore_linkages = TRUE)
}

#' @rdname n_glycan_type
#' @export
n_core_fuc <- function(glycans) {
  core_fuc_motif <- .get_n_glycan_motif("core_fuc")
  glymotif::count_motif(glycans, core_fuc_motif, alignment = "core", ignore_linkages = TRUE)
}

#' @rdname n_glycan_type
#' @export
n_arm_fuc <- function(glycans) {
  arm_fuc_motif <- .get_n_glycan_motif("arm_fuc")
  glymotif::count_motif(glycans, arm_fuc_motif, alignment = "core", ignore_linkages = TRUE)
}

#' @rdname n_glycan_type
#' @export
n_gal <- function(glycans) {
  gal_motif <- .get_n_glycan_motif("gal")
  glymotif::count_motif(glycans, gal_motif, alignment = "substructure", ignore_linkages = TRUE)
}

#' @rdname n_glycan_type
#' @export
n_terminal_gal <- function(glycans) {
  gal_motif <- .get_n_glycan_motif("gal")
  glymotif::count_motif(glycans, gal_motif, alignment = "terminal", ignore_linkages = TRUE)
}

#' @rdname n_glycan_type
#' @export
n_sia <- function(glycans) {
  glyrepr::count_mono(glycans, "NeuAc")
}

#' @rdname n_glycan_type
#' @export
n_man <- function(glycans) {
  gal_motif <- .get_n_glycan_motif("gal")
  n_hex <- glyrepr::count_mono(glycans, "Hex")
  n_gal <- glymotif::count_motif(glycans, gal_motif, alignment = "substructure", ignore_linkages = TRUE)
  n_hex - n_gal
}

.get_n_glycan_motif <- function(name) {
  motif <- switch(
    name,
    core     = glymotif::get_motif_structure("N-Glycan core basic"),
    pauciman = glyparse::parse_iupac_condensed("Hex(??-?)[Hex(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-"),
    hybrid   = glymotif::get_motif_structure("N-Glycan hybrid"),
    highman  = glymotif::get_motif_structure("N-Glycan high mannose"),
    bisect   = glymotif::get_motif_structure("N-glycan core, bisected"),
    antenna  = glyparse::parse_iupac_condensed("HexNAc(??-?)Hex(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"),
    core_fuc = glymotif::get_motif_structure("N-Glycan core, core-fucosylated"),
    arm_fuc  = glymotif::get_motif_structure("N-Glycan core, arm-fucosylated"),
    gal      = glyparse::parse_iupac_condensed("Hex(??-?)HexNAc(??-?)Hex(??-?)Hex(??-")
  )
  motif <- glyrepr::convert_to_generic(motif)
  motif <- glyrepr::remove_linkages(motif)
  motif
}
.get_n_glycan_motif <- memoise::memoise(.get_n_glycan_motif)