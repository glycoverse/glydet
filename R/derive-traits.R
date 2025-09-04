#' Calculate Derived Traits
#'
#' @description
#' This function calculates derived traits from glycomic or glycoproteomic profiles.
#' For glycomics data, it calculates the derived traits directly.
#' For glycoproteomics data, each glycosite is treated as a separate glycome,
#' and derived traits are calculated in a site-specific manner.
#'
#' - `derive_traits()`: Calculate derived traits from a [glyexp::experiment()] object.
#' - `derive_traits_()`: Calculate derived traits from a tibble in tidy format.
#'   Use this function if you are not using the `glyexp` package.
#'
#' @param exp A [glyexp::experiment()] object. Before using this function,
#'   you should preprocess the data using the `glyclean` package.
#'   For glycoproteomics data, the data should be aggregated to the
#'   "gfs" (glycoforms with structures) level using `glyclean::aggregate()`.
#'   Also, please make sure that the `glycan_structure` column is present in the `var_info` table,
#'   as not all glycoproteomics identification softwares provide this information.
#' @param tbl
#'   A tibble in tidy format, with the following columns:
#'   - `sample`: sample ID
#'   - `glycan_structure`: glycan structures, either a `glyrepr::glycan_structure()` vector
#'     or a character vector of glycan structure strings supported by `glyparse::auto_parse()`.
#'   - `value`: the quantification of the glycan in the sample.
#'
#'   For glycoproteomics data, additional columns are needed:
#'   - `protein`: protein ID
#'   - `protein_site`: the glycosite position on the protein
#'   The unique combination of `protein` and `protein_site` determines a glycosite.
#'
#'   Other columns are ignored.
#'
#'   Please make sure that the data has been properly preprocessed,
#'   including normalization, missing value handling, etc.
#'   Specifically, for glycoproteomics data, please make sure that the data has been aggregated to the
#'   "glycoforms with structures" level.
#'   That is the quantification of each glycan structure on each glycosite in each sample.
#' @param data_type Either "glycomics" or "glycoproteomics". Only needed for `derive_traits_()`.
#' @param trait_fns A named list of derived trait functions created by trait factories.
#'   Names of the list are the names of the derived traits.
#'   Default is `NULL`, which means all derived traits in the package are calculated.
#' @param mp_fns A named list of meta-property functions.
#'   This parameter is useful if your trait functions use custom meta-properties
#'   other than those in [all_mp_fns()].
#'   Default is `NULL`, which means all meta-properties in [all_mp_fns()] are used.
#'
#' @returns
#' - For `derive_traits()`: a new [glyexp::experiment()] object for derived traits with "traitomics" type.
#'   Instead of "quantification of each glycan on each glycosite in each sample",
#'   the new `experiment()` contains "the value of each derived trait on each glycosite in each sample",
#'   with the following columns in the `var_info` table:
#'   - `variable`: variable ID
#'   - `trait`: derived trait name
#'
#'   For glycoproteomics data, with additional columns:
#'   - `protein`: protein ID
#'   - `protein_site`: the glycosite position on the protein
#'
#'   Other columns in the `var_info` table (e.g. `gene`) are retained if they have "many-to-one"
#'   relationship with glycosites (unique combinations of `protein`, `protein_site`).
#'   That is, each glycosite cannot have multiple values for these columns.
#'   `gene` is a common example, as a glycosite can only be associate with one gene.
#'   Descriptions about glycans are not such a column, as a glycosite can have multiple glycans,
#'   thus having multiple descriptions.
#'   Columns not having this relationship with glycosites will be dropped.
#'   Don't worry if you cannot understand this logic,
#'   as long as you know that this function will try its best to preserve useful information.
#'
#'   `sample_info` and `meta_data` are not modified,
#'    except that the `exp_type` field of `meta_data` is set to "traitomics".
#'
#' - For `derive_traits_()`: a tidy tibble containing the following columns:
#'   - `sample`: sample ID
#'   - `trait`: derived trait name
#'   - `value`: the value of the derived trait
#'
#'   For glycoproteomics data, with additional columns:
#'   - `protein`: protein ID
#'   - `protein_site`: the glycosite position on the protein
#'
#'   Other columns in the original tibble are not included.
#'
#' @export
derive_traits <- function(exp, trait_fns = NULL, mp_fns = NULL) {
  checkmate::assert_class(exp, "glyexp_experiment")
  if (is.null(trait_fns)) {
    rlang::abort("Not implemented.")
  }

  switch(
    exp$meta_data$exp_type,
    glycomics = .derive_traits_glycomics(exp, trait_fns, mp_fns),
    glycoproteomics = .derive_traits_glycoproteomics(exp, trait_fns, mp_fns),
    cli::cli_abort(c(
      "{.arg exp} must be of type {.val glycomics} or {.val glycoproteomics}.",
      "x" = "Got {.val {exp$meta_data$exp_type}}."
    ))
  )
}

.derive_traits_glycomics <- function(exp, trait_fns, mp_fns) {
  expr_mat <- exp$expr_mat
  glycans <- exp$var_info[["glycan_structure"]]
  mp_tbl <- get_meta_properties(glycans, mp_fns)
  res_mat <- .derive_traits_mat(expr_mat, trait_fns, mp_tbl)
  var_info <- tibble::tibble(
    variable = paste0("V", seq_len(nrow(res_mat))),
    trait = names(trait_fns)
  )
  exp$expr_mat <- res_mat
  exp$var_info <- var_info
  exp$meta_data$exp_type <- "traitomics"
  exp
}

.derive_traits_glycoproteomics <- function(exp, trait_fns, mp_fns) {
  mp_tbl <- get_meta_properties(exp$var_info[["glycan_structure"]], mp_fns)
  glycosites <- stringr::str_c(exp$var_info[["protein"]], exp$var_info[["protein_site"]], sep = "@")
  splits <- split(seq_along(glycosites), glycosites)

  # For now, all other columns in the var_info table are dropped.
  # TODO: Handle other columns.

  derive_one_site <- function(site_idx) {
    site_expr_mat <- exp$expr_mat[site_idx, , drop = FALSE]
    site_mp_tbl <- mp_tbl[site_idx, ]

    res_mat <- .derive_traits_mat(site_expr_mat, trait_fns, site_mp_tbl)
    res_var_info <- tibble::tibble(
      variable = NA_character_,  # Placeholder, will be replaced later
      protein = unique(exp$var_info[["protein"]][site_idx]),
      protein_site = unique(exp$var_info[["protein_site"]][site_idx]),
      trait = names(trait_fns)
    )

    list(res_mat = res_mat, res_var_info = res_var_info)
  }

  res_list <- purrr::map(splits, derive_one_site)
  res_mat <- do.call(rbind, purrr::map(res_list, ~ .x$res_mat))
  res_var_info <- do.call(rbind, purrr::map(res_list, ~ .x$res_var_info))
  res_var_info$variable <- paste0("V", seq_len(nrow(res_var_info)))
  rownames(res_mat) <- res_var_info$variable

  exp$expr_mat <- res_mat
  exp$var_info <- res_var_info
  exp$meta_data$exp_type <- "traitomics"
  exp
}

#' Calculate Derived Traits from a Matrix
#'
#' @param expr_mat A matrix of expression values with samples as columns and glycans as rows.
#' @param trait_fns A named list of derived trait functions.
#' @param mp_tbl A tibble with the meta-properties of the glycans.
#'
#' @returns A matrix of derived traits with samples as columns and traits as rows.
#'   Row names are the names of the derived traits.
#'   Column names are the names of the samples.
#'
#' @noRd
.derive_traits_mat <- function(expr_mat, trait_fns, mp_tbl) {
  res_list <- purrr::map(trait_fns, ~ .x(expr_mat, mp_tbl))
  res_mat <- do.call(rbind, res_list)
  rownames(res_mat) <- names(trait_fns)
  colnames(res_mat) <- colnames(expr_mat)
  res_mat
}

#' @rdname derive_traits
#' @export
derive_traits_ <- function(tbl, data_type, trait_fns, mp_fns) {
  # TODO: implement
}