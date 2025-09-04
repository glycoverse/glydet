# This file contains helper functions for creating derived trait functions.
#
# A derived trait function has the following signature:
# `function(expr_mat, mp_tbl)`
# - `expr_mat`: a expression matrix (samples as columns and glycans as rows)
# - `mp_tbl`: a tibble with the meta-properties of the glycans


#' Create a Propotion Trait
#'
#' A proportion trait is the proportion of certain group of glycans within a larger group of glycans.
#' For example, the proportion of sialylated glycans within all glycans,
#' or the proportion of tetra-antennary glycans within all complex glycans.
#' This type of traits is the most common type of glycan derived traits.
#'
#' @details
#' You can use `prop()` to create proportion trait easily.
#'
#' For example:
#'
#' ```r
#' # Proportion of core-fucosylated glycans within all glycans
#' prop(nFc > 0)
#'
#' # Proportion of complex glycans within all glycans
#' prop(T == "complex")
#'
#' # Proportion of sialylated and fucosylated glycans within all glycans
#' prop(nS > 0 & nFc > 0)
#' ```
#'
#' Note that the last example uses `&` for logical AND.
#' Actually, you can use any logical operator in the expression in R (e.g., `|`, `!`, etc.).
#'
#' If you want to perform a pre-filtering before calculating the proportion,
#' for example, you want to calculate the proportion of core-fucosylated glycans within only complex glycans,
#' you can use `within` to define the denominator.
#'
#' ```r
#' # Proportion of core-fucosylated glycans within complex glycans
#' prop(nFc > 0, within = (T == "complex"))
#'
#' # Proportion of core-fucosylated glycans with tetra-antenary complex glycans
#' prop(nFc > 0, within = (T == "complex" & nA == 4))
#' ```
#'
#' The parentheses around the condition in `within` are optional,
#' but it is recommended to use them for clarity.
#'
#' @param cond Condition to use for defining the smaller group.
#'   An expression that evaluates to a logical vector.
#'   The names of all built-in meta-properties (see [all_mp_fns()]) and custom meta-properties
#'   can be used in the expression.
#' @param within Condition to use for defining the larger group,
#'   with the same format as `cond`.
#'   If `NULL` (default), all glycans are used as the larger group.
#' @param na_action How to handle missing values.
#'   - "keep" (default): keep the missing values as NA.
#'   - "zero": set the missing values to 0.
#'
#' @returns A derived trait function.
#' @export
prop <- function(cond, within = NULL, na_action = "keep") {
  cond <- rlang::enquo(cond)
  within <- rlang::enquo(within)
  checkmate::assert_choice(na_action, c("keep", "zero"))

  function(expr_mat, mp_tbl) {
    cond <- rlang::eval_tidy(cond, data = mp_tbl)
    within <- rlang::eval_tidy(within, data = mp_tbl)
    if (is.null(within)) {
      within <- rep(TRUE, nrow(expr_mat))
    }
    cond[is.na(cond)] <- FALSE
    within[is.na(within)] <- FALSE
    num <- colSums(expr_mat[cond & within, , drop = FALSE], na.rm = TRUE)
    denom <- colSums(expr_mat[within, , drop = FALSE], na.rm = TRUE)
    res <- num / denom
    res[!is.finite(res)] <- NA
    if (na_action == "zero") {
      res[is.na(res)] <- 0
    }
    res
  }
}