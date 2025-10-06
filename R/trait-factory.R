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
#' @section How to use:
#'
#' You can use `prop()` to create proportion trait easily.
#'
#' For example:
#'
#' ```r
#' # Proportion of core-fucosylated glycans within all glycans
#' prop(nFc > 0)
#'
#' # Proportion of complex glycans within all glycans
#' prop(Tp == "complex")
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
#' prop(nFc > 0, within = (Tp == "complex"))
#'
#' # Proportion of core-fucosylated glycans with tetra-antenary complex glycans
#' prop(nFc > 0, within = (Tp == "complex" & nA == 4))
#' ```
#'
#' The parentheses around the condition in `within` are optional,
#' but it is recommended to use them for clarity.
#'
#' @section Note about NA:
#'
#' All the internal summation operations ignore NAs by default.
#' Therefore, NAs in the expression matrix and meta-property values will not result in NAs in the derived traits.
#' However, as all derived traits calculate a ratio of two values,
#' NAs will be introduced when:
#'
#' 1. The denominator is 0. This can happen when the `within` condition selects no glycans.
#' 2. Both the numerator and denominator are 0.
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
#'
#' @examples
#' # Proportion of core-fucosylated glycans within all glycans
#' prop(nFc > 0)
#'
#' # Proportion of bisecting glycans within all glycans
#' prop(B)
#'
#' # Proportion of sialylated and arm-fucosylated glycans within all glycans
#' prop(nS > 0 & nFa > 0)
#'
#' # Proportion of bi-antennary glycans within complex glycans
#' prop(nA == 2, within = (Tp == "complex"))
#'
#' # Proportion of sialylated glycans within core-fucosylated tetra-antennary glycans
#' prop(nS > 0, within = (nFc > 0 & nA == 4))
#'
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

  structure(
    .create_trait_calculator(calc_fn, within_quo = within, na_action = na_action),
    cond = rlang::quo_get_expr(cond),
    within = rlang::quo_get_expr(within),
    na_action = na_action,
    class = c("glydet_prop", "glydet_trait")
  )
}

#' @export
print.glydet_prop <- function(x, ...) {
  cond_expr <- rlang::expr_text(attr(x, "cond"))
  within_expr <- rlang::expr_text(attr(x, "within"))
  na_action <- attr(x, "na_action")
  if (within_expr == "NULL") {
    cli::cli_text("prop({.strong {cond_expr}}, na_action = \"{na_action}\")")
  } else {
    within_expr <- .format_within_expr(within_expr)
    cli::cli_text("prop({.strong {cond_expr}}, within = {.strong ({within_expr})}, na_action = \"{na_action}\")")
  }
  invisible(x)
}

#' Create a Ratio Trait
#'
#' A ratio trait is the ratio of total abundance of two groups of glycans.
#' For example, the ratio of complex glycans and hybrid glycans,
#' or the ratio of bisecting and unbisecting glycans.
#'
#' @section How to use:
#'
#' You can use `ratio()` to create ratio trait easily.
#'
#' For example:
#'
#' ```r
#' # Ratio of complex glycans and hybrid glycans
#' ratio(Tp == "complex", Tp == "hybrid")
#'
#' # Ratio of bisecting and unbisecting glycans
#' ratio(B, !B)
#'
#' # Ratio of core-fucosylated and non-core-fucosylated glycans within complex glycans
#' ratio(nFc > 0 & Tp == "complex", nFc == 0 & Tp == "complex")
#'
#' # The above example can be simplified as:
#' ratio(nFc > 0, nFc == 0, within = (Tp == "complex"))  # more readable
#' ```
#'
#' Note that the last example uses `&` for logical AND.
#' Actually, you can use any logical operator in the expression in R (e.g., `|`, `!`, etc.).
#'
#' `prop()` is a special case of `ratio()`,
#' i.e., `prop(cond, within)` is equivalent to `ratio(cond & within, within)`.
#' We recommend using `prop()` instead of `ratio()` for clarity if possible.
#'
#' @inheritSection prop Note about NA
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
#'
#' @examples
#' # Ratio of complex glycans and hybrid glycans
#' ratio(Tp == "complex", Tp == "hybrid")
#'
#' # Ratio of bisecting and unbisecting glycans within bi-antennary glycans
#' ratio(B, !B, within = (nA == 2))
#'
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

  structure(
    .create_trait_calculator(calc_fn, within_quo = within, na_action = na_action),
    num_cond = rlang::quo_get_expr(num_cond),
    denom_cond = rlang::quo_get_expr(denom_cond),
    within = rlang::quo_get_expr(within),
    na_action = na_action,
    class = c("glydet_ratio", "glydet_trait")
  )
}

#' @export
print.glydet_ratio <- function(x, ...) {
  num_cond_expr <- rlang::expr_text(attr(x, "num_cond"))
  denom_cond_expr <- rlang::expr_text(attr(x, "denom_cond"))
  within_expr <- rlang::expr_text(attr(x, "within"))
  na_action <- attr(x, "na_action")
  if (within_expr == "NULL") {
    cli::cli_text("ratio({.strong {num_cond_expr}}, {.strong {denom_cond_expr}}, na_action = \"{na_action}\")")
  } else {
    within_expr <- .format_within_expr(within_expr)
    cli::cli_text("ratio({.strong {num_cond_expr}}, {.strong {denom_cond_expr}}, within = {.strong ({within_expr})}, na_action = \"{na_action}\")")
  }
  invisible(x)
}

#' Create a Weighted-Mean Trait
#'
#' A weighted-mean trait is the average value of some quantitative property within a group of glycans,
#' weighted by the abundance of the glycans.
#' For example, the average number of antennae within all complex glycans,
#' or the average number of sialic acids within all glycans.
#'
#' @section How to use:
#'
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
#' wmean(nA, within = (Tp == "complex"))
#' ```
#'
#' @inheritSection prop Note about NA
#'
#' @param val Expression to use for defining the value.
#'   An expression that evaluates to a numeric vector.
#'   The names of all built-in meta-properties (see [all_mp_fns()]) and custom meta-properties
#'   can be used in the expression.
#' @param within Condition to set a restriction for the glycans. Same format as `val`.
#' @param na_action How to handle missing values.
#'   - "keep" (default): keep the missing values as NA.
#'   - "zero": set the missing values to 0.
#'
#' @returns A derived trait function.
#'
#' @examples
#' # Weighted mean of the number of sialic acids within all glycans
#' wmean(nS)
#'
#' # Average degree of sialylation per antenna within all glycans
#' wmean(nS / nA)
#'
#' # Average number of antennae within complex glycans
#' wmean(nA, within = (Tp == "complex"))
#'
#' @export
wmean <- function(val, within = NULL, na_action = "keep") {
  val <- rlang::enquo(val)
  within <- rlang::enquo(within)

  # Create calculation function for weighted mean
  calc_fn <- function(expr_mat, mp_tbl, within_cond) {
    val <- rlang::eval_tidy(val, data = mp_tbl)
    val_mat <- expr_mat * val
    num <- colSums(val_mat[within_cond, , drop = FALSE], na.rm = TRUE)
    denom <- colSums(expr_mat[within_cond, , drop = FALSE], na.rm = TRUE)
    list(num = num, denom = denom)
  }

  structure(
    .create_trait_calculator(calc_fn, within_quo = within, na_action = na_action),
    val = rlang::quo_get_expr(val),
    within = rlang::quo_get_expr(within),
    na_action = na_action,
    class = c("glydet_wmean", "glydet_trait")
  )
}

#' @export
print.glydet_wmean <- function(x, ...) {
  val_expr <- rlang::expr_text(attr(x, "val"))
  within_expr <- rlang::expr_text(attr(x, "within"))
  na_action <- attr(x, "na_action")
  if (within_expr == "NULL") {
    cli::cli_text("wmean({.strong {val_expr}}, na_action = \"{na_action}\")")
  } else {
    within_expr <- .format_within_expr(within_expr)
    cli::cli_text("wmean({.strong {val_expr}}, within = {.strong ({within_expr})}, na_action = \"{na_action}\")")
  }
  invisible(x)
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

#' Create a Total Abundance Trait
#'
#' A total abundance trait is the total abundance of a group of glycans.
#' For example, the total abundance of all complex glycans,
#' or the total abundance of all tetra-antennary glycans.
#'
#' @section How to use:
#'
#' You can use `total()` to create total abundance trait easily.
#'
#' For example:
#'
#' ```r
#' # Total abundance of all complex glycans
#' total(Tp == "complex")
#'
#' # Total abundance of all tetra-antennary glycans
#' total(nA == 4)
#' ```
#'
#' @param cond Condition to use for defining the group of glycans.
#'   An expression that evaluates to a logical vector.
#'   The names of all built-in meta-properties (see [all_mp_fns()]) and custom meta-properties
#'   can be used in the expression.
#'
#' @returns A derived trait function.
#'
#' @examples
#' # Total abundance of all complex glycans
#' total(Tp == "complex")
#'
#' # Total abundance of all tetra-antennary glycans
#' total(nA == 4)
#'
#' @export
total <- function(cond) {
  cond <- rlang::enquo(cond)

  # Create calculation function for total abundance
  f <- function(expr_mat, mp_tbl) {
    cond_eval <- rlang::eval_tidy(cond, data = mp_tbl)
    cond_eval[is.na(cond_eval)] <- FALSE
    res <- colSums(expr_mat[cond_eval, , drop = FALSE], na.rm = TRUE)
    unname(res)
  }

  structure(f, cond = rlang::quo_get_expr(cond), class = c("glydet_total", "glydet_trait"))
}

#' @export
print.glydet_total <- function(x, ...) {
  cond_expr <- rlang::expr_text(attr(x, "cond"))
  cli::cli_text("total({.strong {cond_expr}})")
  invisible(x)
}

#' Create a Weighted Sum Trait
#'
#' A weighted sum trait is the sum of a quantitative property within a group of glycans,
#' weighted by the abundance of the glycans.
#' For example, the sum of the number of sialic acids within all glycans,
#' or the sum of the number of Lewis x antigens within all glycans.
#'
#' @section How to use:
#'
#' You can use `wsum()` to create weighted sum trait easily.
#'
#' For example:
#'
#' ```r
#' # Weighted sum of the number of sialic acids within all glycans
#' wsum(nS)
#' ```
#'
#' This can be regarded as the quantification of sialic acids.
#' If some glycan has only one sialic acid,
#' its abundance is added to the results.
#' If another glycan has two sialic acids,
#' its abundance is doubled before being added to the results.
#'
#' You can also use `within` to restrict the weighted sum calculation to specific glycan subsets.
#' For example, you can calculate the weighted sum of the number of sialic acids within complex glycans:
#'
#' ```r
#' wsum(nS, within = (Tp == "complex"))
#' ```
#'
#' @param val Expression to use for defining the value.
#'   An expression that evaluates to a logical vector.
#'   The names of all built-in meta-properties (see [all_mp_fns()]) and custom meta-properties
#'   can be used in the expression.
#' @param within Condition to set a restriction for the glycans. Same format as `val`.
#'
#' @returns A derived trait function.
#'
#' @examples
#' # Weighted sum of the number of sialic acids within all glycans
#' wsum(nS)
#'
#' # Weighted sum of the number of sialic acids within complex glycans
#' wsum(nS, within = (Tp == "complex"))
#'
#' @export
wsum <- function(val, within = NULL) {
  val <- rlang::enquo(val)
  within <- rlang::enquo(within)

  # Create calculation function for weighted sum
  f <- function(expr_mat, mp_tbl) {
    val_eval <- rlang::eval_tidy(val, data = mp_tbl)
    val_mat <- expr_mat * val_eval
    within_eval <- rlang::eval_tidy(within, data = mp_tbl)
    if (is.null(within_eval)) {
      within_eval <- rep(TRUE, nrow(expr_mat))
    }
    within_eval[is.na(within_eval)] <- FALSE
    res <- colSums(val_mat[within_eval, , drop = FALSE], na.rm = TRUE)
    unname(res)
  }

  structure(
    f,
    val = rlang::quo_get_expr(val),
    within = rlang::quo_get_expr(within),
    class = c("glydet_wsum", "glydet_trait")
  )
}

#' @export
print.glydet_wsum <- function(x, ...) {
  val_expr <- rlang::expr_text(attr(x, "val"))
  within_expr <- rlang::expr_text(attr(x, "within"))
  if (within_expr == "NULL") {
    cli::cli_text("wsum({.strong {val_expr}})")
  } else {
    within_expr <- .format_within_expr(within_expr)
    cli::cli_text("wsum({.strong {val_expr}}, within = {.strong ({within_expr})})")
  }
  invisible(x)
}

#' Remove outer parentheses from within expression
#'
#' @param within The `within` character string.
#' @returns The formatted within character string.
#' @noRd
.format_within_expr <- function(within) {
  if (stringr::str_starts(within, stringr::fixed("("))) {
    within <- stringr::str_sub(within, 2, -2)
  }
  within
}