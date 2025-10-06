#' Quantify Motifs in an Experiment
#'
#' @description
#' This function quantifies motifs from glycomic or glycoproteomic profiles.
#' For glycomics data, it calculates motif quantifications directly.
#' For glycoproteomics data, each glycosite is treated as a separate glycome,
#' and motif quantifications are calculated in a site-specific manner.
#'
#' The function takes a `glyexp::experiment()` object and returns a new `glyexp::experiment()`
#' object with motif quantifications. Instead of containing quantifications of individual glycans
#' on each glycosite in each sample, the new experiment contains quantifications
#' of each motif on each glycosite in each sample (for glycoproteomics data) or
#' motif quantifications in each sample (for glycomics data).
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
#' exp_with_mps <- glymotif::add_motifs_int(exp, motifs)
#'
#' # Define the traits
#' trait_fns <- list(Lx = wsum(nLx), SLx = wsum(nSLx))
#'
#' # Calculate the traits
#' derive_traits(exp_with_mps, trait_fns = trait_fns, mp_cols = c("nLx", "nSLx"))
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
#'   or a 'glyrepr_structure' object.
#'   For glycan structure strings, all formats supported by [glyparse::auto_parse()] are accepted,
#'   including IUPAC-condensed, WURCS, GlycoCT, and others.
#'   If the vector is named, the names will be used as motif names.
#'   Otherwise, IUPAC-condensed structure strings will be used as motif names.
#' @param method A character string specifying the quantification method.
#'   Must be either "absolute" or "relative". Default is "relative".
#'   See "Relative and Absolute Motif Quantification" section for details.
#'
#' @returns
#' A new [glyexp::experiment()] object containing motif quantifications.
#' Instead of containing quantifications of individual glycans on each glycosite in each sample,
#' the new experiment contains quantifications of each motif on each glycosite in each sample
#' (for glycoproteomics data) or motif quantifications in each sample (for glycomics data).
#'
#' For glycoproteomics data, with additional columns:
#' - `protein`: protein ID
#' - `protein_site`: the glycosite position on the protein
#'
#' Other columns in the `var_info` table (e.g., `gene`) are retained if they have a "many-to-one"
#' relationship with glycosites (unique combinations of `protein`, `protein_site`).
#' That is, each glycosite cannot have multiple values for these columns.
#' `gene` is a common example, as a glycosite can only be associated with one gene.
#' Glycan descriptions are not such columns, as a glycosite can have multiple glycans,
#' thus having multiple descriptions.
#' Columns that do not have this relationship with glycosites will be dropped.
#' Don't worry if you cannot understand this logic—
#' just know that this function will do its best to preserve useful information.
#'
#' The `sample_info` and `meta_data` tables are not modified,
#' except that the `exp_type` field in `meta_data` is set to "traitomics" for glycomics data
#' and "traitproteomics" for glycoproteomics data.
#'
#' @examples
#' library(glyexp)
#' library(glyclean)
#'
#' exp <- real_experiment |>
#'   auto_clean() |>
#'   slice_head_var(n = 10)
#'
#' motifs <- c(
#'   nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
#'   nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
#' )
#'
#' quantify_motifs(exp, motifs)
#'
#' @seealso [derive_traits()]
#' @export
quantify_motifs <- function(exp, motifs, method = "relative") {
  checkmate::assert_class(exp, "glyexp_experiment")
  # `motifs` is validated in `glymotif::add_motifs_int()`
  checkmate::assert_choice(method, c("absolute", "relative"))

  # Add meta-properties columns to the variable information tibble
  exp2 <- glymotif::add_motifs_int(exp, motifs)
  # `add_motifs_int()` has a complex logic of determining the column names,
  # so we use a simpler approach to get the column names.
  mp_cols <- setdiff(colnames(exp2$var_info), colnames(exp$var_info))

  # Define the trait functions
  factory <- switch(method, absolute = wsum, relative = wmean)
  trait_fns <- vector("list", length(mp_cols))
  names(trait_fns) <- mp_cols
  for (i in seq_along(mp_cols)) {
    trait_fns[[i]] <- factory(!!rlang::sym(mp_cols[[i]]))
  }

  # Calculate the traits
  trait_exp <- derive_traits(exp2, trait_fns = trait_fns, mp_cols = mp_cols)
  trait_exp |>
    glyexp::rename_var(dplyr::all_of(c("motif" = "trait"))) |>
    # Remove meta-property columns (if any)
    # `derive_traits()` has a special logic of rescuing columns
    # that have "many-to-one" relationship with glycosites.
    # This behavior might result in some meta-property columns being left over.
    # We remove them here.
    glyexp::select_var(-dplyr::any_of(mp_cols))
}