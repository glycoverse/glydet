# This file contains helper functions for creating derived trait functions.
#
# A derived trait function has the following signature:
# `function(expr_mat, mp_tbl)`
# - `expr_mat`: a expression matrix (samples as columns and glycans as rows)
# - `mp_tbl`: a tibble with the meta-properties of the glycans

#' Create a Proportion Trait
#'
#' A proportion trait is the proportion of certain group of glycans within a larger group of glycans.
#' For example, the proportion of sialylated glycans within all glycans,
#' or the proportion of tetra-antennary glycans within all complex glycans.
#' This type of traits is the most common type of glycan derived traits.
#' It can be regarded as a special case of the ratio trait (see [ratio()]).
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

  # Create calculation function for proportion
  calc_fn <- function(expr_mat, mp_tbl, within_cond) {
    cond_eval <- rlang::eval_tidy(cond, data = mp_tbl)
    cond_eval[is.na(cond_eval)] <- FALSE
    num <- colSums(expr_mat[cond_eval & within_cond, , drop = FALSE], na.rm = TRUE)
    denom <- colSums(expr_mat[within_cond, , drop = FALSE], na.rm = TRUE)
    list(num = num, denom = denom)
  }

  .create_trait_calculator(calc_fn, within_quo = within, na_action = na_action)
}

#' Create a Ratio Trait
#'
#' A ratio trait is the ratio of total abundance of two groups of glycans.
#' For example, the ratio of complex glycans and hybrid glycans,
#' or the ratio of bisecting and unbisecting glycans.
#'
#' @details
#' You can use `ratio()` to create ratio trait easily.
#'
#' For example:
#'
#' ```r
#' # Ratio of complex glycans and hybrid glycans
#' ratio(T == "complex", T == "hybrid")
#'
#' # Ratio of bisecting and unbisecting glycans
#' ratio(B, !B)
#'
#' # Ratio of core-fucosylated and non-core-fucosylated glycans within complex glycans
#' ratio(nFc > 0 & T == "complex", nFc == 0 & T == "complex")
#'
#' # The above example can be simplified as:
#' ratio(nFc > 0, nFc == 0, within = (T == "complex"))  # more readable
#' ```
#'
#' Note that the last example uses `&` for logical AND.
#' Actually, you can use any logical operator in the expression in R (e.g., `|`, `!`, etc.).
#'
#' `prop()` is a special case of `ratio()`,
#' i.e., `prop(cond, within)` is equivalent to `ratio(cond & within, within)`.
#' We recommend using `prop()` instead of `ratio()` for clarity if possible.
#'
#' @param num_cond Condition to use for defining the numerator.
#'   An expression that evaluates to a logical vector.
#'   The names of all built-in meta-properties (see [all_mp_fns()]) and custom meta-properties
#'   can be used in the expression.
#' @param denom_cond Condition to use for defining the denominator. Same format as `num_cond`.
#' @param within Condition to set a restriction for the glycans. Same format as `num_cond`.
#' @param na_action How to handle missing values.
#'   - "keep" (default): keep the missing values as NA.
#'   - "zero": set the missing values to 0.
#'
#' @returns A derived trait function.
#' @export
ratio <- function(num_cond, denom_cond, within = NULL, na_action = "keep") {
  num_cond <- rlang::enquo(num_cond)
  denom_cond <- rlang::enquo(denom_cond)
  within <- rlang::enquo(within)

  # Create calculation function for ratio
  calc_fn <- function(expr_mat, mp_tbl, within_cond) {
    num_cond_eval <- rlang::eval_tidy(num_cond, data = mp_tbl)
    denom_cond_eval <- rlang::eval_tidy(denom_cond, data = mp_tbl)
    num_cond_eval[is.na(num_cond_eval)] <- FALSE
    denom_cond_eval[is.na(denom_cond_eval)] <- FALSE
    num <- colSums(expr_mat[num_cond_eval & within_cond, , drop = FALSE], na.rm = TRUE)
    denom <- colSums(expr_mat[denom_cond_eval & within_cond, , drop = FALSE], na.rm = TRUE)
    list(num = num, denom = denom)
  }

  .create_trait_calculator(calc_fn, within_quo = within, na_action = na_action)
}

#' Create a Weighted-Mean Trait
#'
#' A weighted-mean trait is the average value of some quantitative property within a group of glycans,
#' weighted by the abundance of the glycans.
#' For example, the average number of antennae within all complex glycans,
#' or the average number of sialic acids within all glycans.
#'
#' @details
#' You can use `wmean()` to create weighted-mean trait easily.
#'
#' For example:
#'
#' ```r
#' # Weighted mean of the number of sialic acids within all glycans
#' wmean(nS)
#'
#' # Average degree of sialylation per antenna within all glycans
#' wmean(nS / nA)
#' ```
#'
#' Note that the last example uses `/` for division.
#' Actually, you can use any arithmetic operator in the expression in R (e.g., `*`, `+`, `-`, etc.).
#'
#' If you want to perform a pre-filtering before calculating the weighted-mean,
#' for example, you want to calculate the average degree of sialylation per antenna within only complex glycans,
#' you can use `within` to define the restriction.
#'
#' ```r
#' # Average number of antennae within complex glycans
#' wmean(nA, within = (T == "complex"))
#' ```
#'
#' @param val_cond Condition to use for defining the value.
#'   An expression that evaluates to a numeric vector.
#'   The names of all built-in meta-properties (see [all_mp_fns()]) and custom meta-properties
#'   can be used in the expression.
#' @param within Condition to set a restriction for the glycans. Same format as `val_cond`.
#' @param na_action How to handle missing values.
#'   - "keep" (default): keep the missing values as NA.
#'   - "zero": set the missing values to 0.
#'
#' @returns A derived trait function.
#' @export
wmean <- function(val_cond, within = NULL, na_action = "keep") {
  val_cond <- rlang::enquo(val_cond)
  within <- rlang::enquo(within)

  # Create calculation function for weighted mean
  calc_fn <- function(expr_mat, mp_tbl, within_cond) {
    val <- rlang::eval_tidy(val_cond, data = mp_tbl)
    val_mat <- expr_mat * val
    num <- colSums(val_mat[within_cond, , drop = FALSE], na.rm = TRUE)
    denom <- colSums(expr_mat[within_cond, , drop = FALSE], na.rm = TRUE)
    list(num = num, denom = denom)
  }

  .create_trait_calculator(calc_fn, within_quo = within, na_action = na_action)
}

# Internal helper function to create trait calculators with common logic
.create_trait_calculator <- function(calc_fn, within_quo, na_action = "keep") {
  checkmate::assert_choice(na_action, c("keep", "zero"))

  function(expr_mat, mp_tbl) {
    # Evaluate within condition
    within_cond <- rlang::eval_tidy(within_quo, data = mp_tbl)
    if (is.null(within_cond)) {
      within_cond <- rep(TRUE, nrow(expr_mat))
    }
    within_cond[is.na(within_cond)] <- FALSE

    # Calculate numerator and denominator using the provided calculation function
    result <- calc_fn(expr_mat, mp_tbl, within_cond)

    # Apply common post-processing
    res <- result$num / result$denom
    res[!is.finite(res)] <- NA
    if (na_action == "zero") {
      res[is.na(res)] <- 0
    }
    unname(res)
  }
}