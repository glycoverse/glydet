#' Quantify Motifs in an Experiment
#'
#' @description
#' This function quantifies motifs from glycomic or glycoproteomic profiles.
#' For glycomics data, it calculates motif quantifications directly.
#' For glycoproteomics data, each glycosite is treated as a separate glycome,
#' and motif quantifications are calculated in a site-specific manner.
#'
#' The function takes a `glyexp::GlycomicSE` or `glyexp::GlycoproteomicSE`
#' object and returns a plain `SummarizedExperiment` with motif quantifications.
#' Instead of containing quantifications of individual glycans on each glycosite
#' in each sample, the output assay contains quantifications of each motif on
#' each glycosite in each sample (for glycoproteomics data) or motif
#' quantifications in each sample (for glycomics data).
#'
#' @details
#' # Relative and Absolute Motif Quantification
#'
#' Motif quantification can be performed in two ways: absolute and relative.
#' Do not confuse this with the absolute and relative quantification of glycans or glycopeptides.
#'
#' In an omics context, absolute quantification means we can determine the actual concentration
#' (e.g., in mg/L) or counts (e.g., mRNA copy numbers),
#' while relative quantification means we can only compare the abundance of the same
#' molecule between different samples,
#' but the absolute values themselves are not meaningful.
#' Label-free quantification is a relative quantification method.
#'
#' In the context of motif quantification,
#' absolute and relative refer to whether we normalize the motif quantifications
#' by the total abundance of the glycome (or the mini-glycomes of individual glycosites).
#'
#' For example, let's say we have a glycomics dataset with two samples, A and B.
#' In each sample, three glycans (G1, G2, G3) are present,
#' with the following abundances:
#'
#' - Sample A: G1 = 10, G2 = 20, G3 = 30
#' - Sample B: G1 = 20, G2 = 40, G3 = 60
#'
#' For simplicity, let's say the motif we want to quantify happens to appear once in all glycans.
#'
#' For absolute motif quantification, we simply sum the abundances of the motif across all glycans:
#'
#' - Sample A: 10 × 1 + 20 × 1 + 30 × 1 = 60  (× 1 because the motif appears once in each glycan)
#' - Sample B: 20 × 1 + 40 × 1 + 60 × 1 = 120
#'
#' The results are clearly different.
#'
#' However, if we quantify the motif in a relative way:
#'
#' - Sample A: (10 × 1 + 20 × 1 + 30 × 1) / (10 + 20 + 30) = 60 / 60 = 1
#' - Sample B: (20 × 1 + 40 × 1 + 60 × 1) / (20 + 40 + 60) = 120 / 120 = 1
#'
#' The results are identical!
#'
#' The absolute motif quantification answers the question "how many motifs are there in the sample?",
#' while the relative motif quantification answers the question
#' "if I take out one glycan molecule, how many motifs are on it in average?"
#'
#' So which method should you use?
#' It depends on both your data type and research objectives.
#' Here are some general guidelines:
#'
#' - For glycomics data, you should typically use relative motif quantification,
#'   because glycomics data is inherently compositional (search "compositional data" for more details).
#' - For glycoproteomics data, the choice depends on your research question.
#'   If you want to compare observed motif abundance across samples, use absolute motif quantification.
#'   Note that many factors can cause one sample to have higher values than another,
#'   such as upregulation of enzymes responsible for the motif,
#'   or higher overall glycosylation site occupancy.
#'   If you want to understand the underlying regulatory mechanisms,
#'   use relative motif quantification to correct for differences in site occupancy.
#'
#' # Relationship with Derived Traits
#'
#' Motif quantification is a special type of derived trait.
#' It is simply a weighted sum ([wsum()]) or weighted mean ([wmean()]) of the motif counts,
#' for absolute and relative motif quantification, respectively.
#' You can perform motif quantification manually using [derive_traits()].
#'
#' Let's perform absolute motif quantification manually using [derive_traits()]
#' to better understand the process (we will use custom column meta-properties here):
#'
#' ```r
#' # Add the meta-properties to the variable information tibble
#' motifs <- c(
#'   nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
#'   nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
#' )
#' exp_with_mps <- exp |>
#'   glyexp::mutate_row(
#'     tibble::as_tibble(glymotif::count_motifs(glycan_structure, motifs))
#'   )
#'
#' # Define the traits
#' trait_fns <- list(Lx = wsum(nLx), SLx = wsum(nSLx))
#'
#' # Calculate the traits
#' derive_traits(exp_with_mps, trait_fns = trait_fns)
#' ```
#'
#' The code snippet above is equivalent to:
#'
#' ```r
#' # Quantify the motifs
#' quantify_motifs(exp, motifs, method = "absolute")
#' ```
#'
#' In fact, this is essentially the implementation of the `quantify_motifs()` function,
#' except the actual implementation is more robust and user-friendly.
#'
#' For relative motif quantification, simply replace `wsum()` with `wmean()`,
#' and everything else remains the same.
#'
#' @inheritParams derive_traits
#' @param motifs A character vector of motif names, glycan structure strings,
#'   a 'glyrepr_structure' object, or a motif specification from [glymotif::dynamic_motifs()]
#'   or [glymotif::branch_motifs()].
#'   For glycan structure strings, all formats supported by [glyparse::auto_parse()] are accepted,
#'   including IUPAC-condensed, WURCS, GlycoCT, and others.
#'   If the vector is named, the names will be used as motif names.
#'   Otherwise, IUPAC-condensed structure strings will be used as motif names.
#'   For motif specifications, motifs are extracted automatically from the glycan structures
#'   in the experiment, and their IUPAC-condensed strings are used as motif names.
#' @param method A character string specifying the quantification method.
#'   Must be either "absolute" or "relative". Default is "relative".
#'   See "Relative and Absolute Motif Quantification" section for details.
#' @param alignments A character vector specifying the alignment method for each motif.
#'   Can be "terminal", "substructure", "core", "whole" or "exact". Default is "substructure".
#'   See [glymotif::have_motifs()] for details.
#' @param ignore_linkages A logical value. If `TRUE`, linkages will be ignored in the comparison.
#'   See [glymotif::have_motifs()] for details.
#'
#' @returns
#' New `GlycomicSE` and `GlycoproteomicSE` inputs return a plain
#' `SummarizedExperiment`; compatible legacy glyexp inputs preserve their
#' legacy container type. The output contains motif quantifications.
#' Instead of containing quantifications of individual glycans on each glycosite in each sample,
#' the output assay contains quantifications of each motif on each glycosite in each sample
#' (for glycoproteomics data) or motif quantifications in each sample (for glycomics data).
#'
#' `rowData()` includes a `trait` column with motif names and a
#' `motif_structure` column containing the parsed glycan structure for each motif,
#' allowing traceability of motif definitions.
#'
#' For glycoproteomics data, with additional columns:
#' - `protein`: protein ID
#' - `protein_site`: the glycosite position on the protein
#'
#' Other columns in `rowData()` (e.g., `gene`) are retained if they have a "many-to-one"
#' relationship with glycosites (unique combinations of `protein`, `protein_site`).
#' That is, each glycosite cannot have multiple values for these columns.
#' `gene` is a common example, as a glycosite can only be associated with one gene.
#' Glycan descriptions are not such columns, as a glycosite can have multiple glycans,
#' thus having multiple descriptions.
#' Columns that do not have this relationship with glycosites will be dropped.
#' Don't worry if you cannot understand this logic—
#' just know that this function will do its best to preserve useful information.
#'
#' `colData()` and `metadata()` are not modified, except that the `exp_type`
#' field in `metadata()` is set to "traitomics" for glycomics data and
#' "traitproteomics" for glycoproteomics data.
#'
#' @examples
#' library(glyexp)
#' library(SummarizedExperiment)
#' library(glyclean)
#'
#' gp_se <- real_experiment |>
#'   auto_clean() |>
#'   slice_head_row(n = 10)
#'
#' motifs <- c(
#'   nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
#'   nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
#' )
#'
#' motif_se <- quantify_motifs(gp_se, motifs)
#' rowData(motif_se)
#' assay(motif_se)
#' colData(motif_se)
#'
#' # Using dynamic motifs (auto-extracted from data)
#' quantify_motifs(gp_se, glymotif::dynamic_motifs(max_size = 3))
#'
#' # Using branch motifs (auto-extracted from data)
#' quantify_motifs(gp_se, glymotif::branch_motifs())
#'
#' @seealso [derive_traits()], [glymotif::have_motifs()]
#' @export
quantify_motifs <- function(
  exp,
  motifs,
  method = "relative",
  alignments = NULL,
  ignore_linkages = FALSE
) {
  .assert_data_container(exp)
  legacy <- inherits(exp, "glyexp_experiment")
  exp <- .as_glyco_se(exp)
  checkmate::assert_choice(method, c("absolute", "relative"))

  # ----- Create motif lookup tibble -----
  # Handle three cases: glycan_structure vector, known motif names, or IUPAC strings
  if (glyrepr::is_glycan_structure(motifs)) {
    # Case: glycan_structure vector (no names possible)
    motif_structures <- motifs
  } else if (is.character(motifs)) {
    if (.are_db_motif_names(motifs)) {
      # Case: Known motif names from database
      motif_structures <- .get_db_motif_structure(motifs)
    } else {
      # Case: IUPAC structure strings
      motif_structures <- glyparse::auto_parse(motifs)
    }
  } else if (
    inherits(motifs, "dynamic_motifs_spec") ||
      inherits(motifs, "branch_motifs_spec")
  ) {
    # Case: Motif spec objects - resolved while counting motif meta-properties
    motif_structures <- NULL # Will be set from the resulting motif columns
  } else {
    rlang::abort(
      "`motifs` must be a character vector, a 'glyrepr_structure' object, or a motif specification from `dynamic_motifs()` or `branch_motifs()`."
    )
  }

  exp2 <- .add_motif_count_mps(
    exp,
    motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages
  )
  mp_cols <- setdiff(
    colnames(.get_var_info(exp2)),
    colnames(.get_var_info(exp))
  )

  # For motif specs, extract structures from column names (IUPAC strings)
  is_motif_spec <- inherits(motifs, "dynamic_motifs_spec") ||
    inherits(motifs, "branch_motifs_spec")
  if (is_motif_spec) {
    motif_structures <- glyparse::parse_iupac_condensed(mp_cols)
  }

  # Create lookup tibble: motif names -> motif structures
  motif_lookup <- tibble::tibble(
    trait = mp_cols,
    motif_structure = unname(motif_structures)
  )

  # Define the trait functions
  factory <- switch(method, absolute = wsum, relative = wmean)
  trait_fns <- vector("list", length(mp_cols))
  names(trait_fns) <- mp_cols
  for (i in seq_along(mp_cols)) {
    trait_fns[[i]] <- factory(!!rlang::sym(mp_cols[[i]]))
  }

  # Calculate the traits
  trait_exp <- derive_traits(exp2, trait_fns = trait_fns, mp_cols = mp_cols)
  # Remove meta-property columns (if any). `derive_traits()` has a special
  # logic of rescuing columns that have a "many-to-one" relationship with
  # glycosites, which can leave some meta-property columns in the result.
  result_var_info <- .get_var_info(trait_exp) |>
    dplyr::select(-dplyr::any_of(mp_cols)) |>
    dplyr::select(-dplyr::all_of("explanation")) |>
    dplyr::left_join(motif_lookup, by = "trait")
  result <- .set_var_info(trait_exp, result_var_info)
  .restore_data_container(result, legacy)
}

#' Add motif count meta-property columns
#'
#' @param exp A legacy glyexp data container.
#' @inheritParams quantify_motifs
#'
#' @returns A legacy glyexp data container with motif count columns added to its
#'   variable information.
#' @noRd
.add_motif_count_mps <- function(
  exp,
  motifs,
  alignments = NULL,
  ignore_linkages = FALSE
) {
  .check_var_info_cols(exp, "glycan_structure")

  motif_counts <- glymotif::count_motifs(
    .get_var_info(exp)[["glycan_structure"]],
    motifs,
    alignments = alignments,
    ignore_linkages = ignore_linkages
  )
  colnames(motif_counts) <- .motif_count_colnames(motifs, motif_counts)

  .set_var_info(
    exp,
    dplyr::bind_cols(.get_var_info(exp), tibble::as_tibble(motif_counts))
  )
}

#' Get motif count column names
#'
#' @param motifs Motif definitions passed to [quantify_motifs()].
#' @param motif_counts A count matrix returned by [glymotif::count_motifs()].
#'
#' @returns A character vector of column names for motif count meta-properties.
#' @noRd
.motif_count_colnames <- function(motifs, motif_counts) {
  count_names <- colnames(motif_counts)
  if (!is.null(count_names)) {
    return(count_names)
  }

  fallback_names <- NULL
  if (glyrepr::is_glycan_structure(motifs)) {
    fallback_names <- as.character(motifs)
  } else if (is.character(motifs)) {
    fallback_names <- motifs
  }

  if (is.null(fallback_names)) {
    cli::cli_abort(
      "Could not determine motif count column names from {.arg motifs}."
    )
  }

  motif_names <- names(motifs)
  if (is.null(motif_names)) {
    motif_names <- rep("", length(fallback_names))
  }
  unnamed <- !nzchar(motif_names)
  motif_names[unnamed] <- fallback_names[unnamed]
  motif_names
}
