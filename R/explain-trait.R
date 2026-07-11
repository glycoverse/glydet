#' Explain Derived Traits
#'
#' @description
#' These functions provide human-readable English explanations of derived traits.
#' `explain_trait()` explains one trait, while `explain_traits()` explains a list
#' of traits. They work with trait functions created by [prop()], [ratio()],
#' [wmean()], [total()], and [wsum()], and work best with traits defined by
#' built-in meta-properties. When `use_ai = TRUE`, `explain_traits()` sends all
#' valid traits in one request so the shared prompt is only sent once. Traits
#' that cannot be explained are returned as `NA` with a warning.
#'
#' @param trait_fn A derived trait function created by one of the trait factories.
#' @param trait_fns A list of derived trait functions created by the trait
#'   factories.
#' @param use_ai `r lifecycle::badge("experimental")`
#'   Whether to use a Large Language Model (LLM) to explain the trait. Default is FALSE.
#'   To use this feature, you need to install the `ellmer` package.
#'   DeepSeek is used by default for backward compatibility. Other `ellmer`
#'   providers can be selected with `provider`, `model`, and provider-specific
#'   API key configuration.
#' @param custom_mp A named character vector of custom meta-properties.
#'   The names are the meta-property names, and the values are in the format
#'   "(type) description". Only used when `use_ai = TRUE`.
#' @param provider AI provider passed to `ellmer` when `use_ai = TRUE`.
#'   One of "deepseek", "openai", "anthropic", "gemini", "openrouter", or
#'   "openai_compatible". "google_gemini" is accepted as an alias for
#'   "gemini". Defaults to
#'   `getOption("glydet.ai_provider", "deepseek")`.
#' @param model Model to use when `use_ai = TRUE`. Defaults to
#'   `getOption("glydet.ai_model")`, or "deepseek-chat" for DeepSeek and the
#'   provider default for other providers.
#' @param api_key API key for the selected provider. If `NULL`, the provider
#'   specific environment variable is used. Defaults to
#'   `getOption("glydet.ai_api_key")`.
#' @param base_url Optional base URL for custom or OpenAI-compatible endpoints.
#'   Defaults to `getOption("glydet.ai_base_url")`.
#'
#' @returns `explain_trait()` returns a character string containing a concise
#'   English explanation. `explain_traits()` returns a character vector with
#'   input names preserved; entries that cannot be explained are `NA`.
#'
#' @examples
#' # Explain built-in traits
#' explain_trait(traits_basic()$TM)
#' explain_trait(traits_basic()$GS)
#'
#' # Explain custom traits
#' explain_trait(prop(nFc > 0))
#' explain_trait(prop(nFc > 0, within = (Tp == "complex")))
#' explain_trait(ratio(Tp == "complex", Tp == "hybrid"))
#' explain_trait(wmean(nA, within = (Tp == "complex")))
#' explain_trait(wmean(nS / nG, within = nA == 4 & nFc > 0))
#'
#' # Explain total and wsum traits
#' explain_trait(total(Tp == "complex"))
#' explain_trait(wsum(nS))
#' explain_trait(wsum(nS, within = (Tp == "complex")))
#'
#' # Explain multiple traits
#' explain_traits(traits_basic()[c("TM", "GS")])
#'
#' @export
explain_trait <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  UseMethod("explain_trait")
}

#' @rdname explain_trait
#' @export
explain_traits <- function(
  trait_fns,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  checkmate::assert_list(trait_fns)

  explanations <- rep(NA_character_, length(trait_fns))
  names(explanations) <- names(trait_fns)
  valid <- purrr::map_lgl(
    trait_fns,
    ~ is.function(.x) && inherits(.x, "glydet_trait")
  )

  if (any(valid) && !use_ai) {
    explanations[valid] <- purrr::map_chr(
      trait_fns[valid],
      ~ tryCatch(
        explain_trait(.x),
        error = function(error) NA_character_
      )
    )
  }

  if (any(valid) && use_ai) {
    formulas <- purrr::map_chr(
      trait_fns[valid],
      ~ tryCatch(
        .trait_to_formula(.x),
        error = function(error) NA_character_
      )
    )
    positions <- which(valid)[!is.na(formulas)]

    if (length(positions) > 0) {
      system_prompt <- paste0(
        .explain_sys_prompt("all", custom_mp),
        "\nThe inputs may use prop(), ratio(), wmean(), total(), or wsum().\n",
        "Explain every input independently. Return exactly one line per input ",
        "in the form `index<TAB>explanation`. Return `index<TAB><INVALID>` ",
        "when an input cannot be explained."
      )
      user_prompt <- paste0(
        "INPUTS:\n",
        paste0(positions, "\t", formulas[!is.na(formulas)], collapse = "\n"),
        "\nOUTPUT:"
      )
      response <- .ask_ai(
        system_prompt,
        user_prompt,
        api_key = api_key,
        model = model,
        provider = provider,
        base_url = base_url
      )
      explanations[positions] <- .parse_batch_response(
        response,
        length(trait_fns)
      )[positions]
    }
  }

  invalid_positions <- which(is.na(explanations))
  if (length(invalid_positions) > 0) {
    cli::cli_warn(
      "Could not explain trait(s) at position(s): {paste(invalid_positions, collapse = ', ')}."
    )
  }

  explanations
}

#' @export
explain_trait.default <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  cli::cli_abort(c(
    "Object must be a glydet trait function.",
    "x" = "Got an object of class {.cls {class(trait_fn)}}.",
    "i" = "Use trait factories like {.fn prop}, {.fn ratio}, or {.fn wmean} to create trait functions."
  ))
}

#' @export
explain_trait.glydet_prop <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  if (use_ai) {
    .explain_with_ai(
      trait_fn,
      custom_mp,
      provider = provider,
      model = model,
      api_key = api_key,
      base_url = base_url
    )
  } else {
    .explain_prop(trait_fn)
  }
}

#' @export
explain_trait.glydet_ratio <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  if (use_ai) {
    .explain_with_ai(
      trait_fn,
      custom_mp,
      provider = provider,
      model = model,
      api_key = api_key,
      base_url = base_url
    )
  } else {
    .explain_ratio(trait_fn)
  }
}

#' @export
explain_trait.glydet_wmean <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  if (use_ai) {
    .explain_with_ai(
      trait_fn,
      custom_mp,
      provider = provider,
      model = model,
      api_key = api_key,
      base_url = base_url
    )
  } else {
    .explain_wmean(trait_fn)
  }
}

#' @export
explain_trait.glydet_total <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  if (use_ai) {
    .explain_with_ai(
      trait_fn,
      custom_mp,
      provider = provider,
      model = model,
      api_key = api_key,
      base_url = base_url
    )
  } else {
    .explain_total(trait_fn)
  }
}

#' @export
explain_trait.glydet_wsum <- function(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  if (use_ai) {
    .explain_with_ai(
      trait_fn,
      custom_mp,
      provider = provider,
      model = model,
      api_key = api_key,
      base_url = base_url
    )
  } else {
    .explain_wsum(trait_fn)
  }
}

# ----- Hard-coded explanations for built-in traits --------------------

.explain_prop <- function(trait_fn) {
  cond_expr <- attr(trait_fn, "cond")
  within_expr <- attr(trait_fn, "within")

  cond_desc <- .expr_to_description(cond_expr)
  scope_desc <- .within_to_scope(within_expr)

  paste0("Proportion of ", cond_desc, " ", scope_desc, ".")
}

.explain_ratio <- function(trait_fn) {
  num_cond_expr <- attr(trait_fn, "num_cond")
  denom_cond_expr <- attr(trait_fn, "denom_cond")
  within_expr <- attr(trait_fn, "within")

  num_desc <- .expr_to_description(num_cond_expr)
  denom_desc <- .expr_to_description(denom_cond_expr)
  scope_desc <- .within_to_scope(within_expr)

  paste0("Ratio of ", num_desc, " to ", denom_desc, " ", scope_desc, ".")
}

.explain_wmean <- function(trait_fn) {
  val_expr <- attr(trait_fn, "val")
  within_expr <- attr(trait_fn, "within")

  val_desc <- .expr_to_value_description(val_expr)
  scope_desc <- .within_to_scope(within_expr)

  paste0("Abundance-weighted mean of ", val_desc, " ", scope_desc, ".")
}

.explain_total <- function(trait_fn) {
  cond_expr <- attr(trait_fn, "cond")
  cond_desc <- .expr_to_description(cond_expr)
  paste0("Total abundance of ", cond_desc, ".")
}

.explain_wsum <- function(trait_fn) {
  val_expr <- attr(trait_fn, "val")
  within_expr <- attr(trait_fn, "within")

  val_desc <- .expr_to_value_description(val_expr)
  scope_desc <- .within_to_scope(within_expr)

  paste0("Abundance-weighted sum of ", val_desc, " ", scope_desc, ".")
}

# ----- Internal helper functions --------------------------------------

.within_to_scope <- function(within_expr) {
  if (is.null(within_expr)) {
    return("among all glycans")
  }

  # Remove outer parentheses if present (consistent with print methods)
  within_text <- rlang::expr_text(within_expr)
  if (
    stringr::str_starts(within_text, stringr::fixed("(")) &&
      stringr::str_ends(within_text, stringr::fixed(")"))
  ) {
    within_expr <- rlang::parse_expr(stringr::str_sub(within_text, 2, -2))
  }

  within_desc <- .expr_to_description(within_expr)
  paste0("within ", within_desc)
}

.expr_to_description <- function(expr) {
  if (is.null(expr)) {
    return("all glycans")
  }

  # Handle special cases for common patterns
  desc <- .try_special_patterns(expr)
  if (!is.null(desc)) {
    return(desc)
  }

  # Handle logical operators recursively
  if (rlang::is_call(expr)) {
    op <- rlang::call_name(expr)
    args <- rlang::call_args(expr)

    if (op == "&" && length(args) == 2) {
      left_desc <- .expr_to_description(args[[1]])
      right_desc <- .expr_to_description(args[[2]])
      # Try to merge similar descriptions for better readability
      merged <- .try_merge_descriptions(left_desc, right_desc, "and")
      if (!is.null(merged)) {
        return(merged)
      }
      return(paste(left_desc, "and", right_desc))
    }

    if (op == "|" && length(args) == 2) {
      left_desc <- .expr_to_description(args[[1]])
      right_desc <- .expr_to_description(args[[2]])
      # Try to merge similar descriptions for better readability
      merged <- .try_merge_descriptions(left_desc, right_desc, "or")
      if (!is.null(merged)) {
        return(merged)
      }
      return(paste(left_desc, "or", right_desc))
    }

    if (op == "!" && length(args) == 1) {
      arg_desc <- .expr_to_description(args[[1]])
      # Handle special negation cases
      if (stringr::str_starts(arg_desc, "glycans with")) {
        return(stringr::str_replace(
          arg_desc,
          "glycans with",
          "glycans without"
        ))
      }
      if (stringr::str_ends(arg_desc, " glycans")) {
        return(paste("non-", arg_desc, sep = ""))
      }
      return(paste("not", arg_desc))
    }
  }

  # Fallback to generic description for complex expressions
  expr_text <- rlang::expr_text(expr)
  paste0("glycans satisfying '", expr_text, "'")
}

.expr_to_value_description <- function(expr) {
  if (is.null(expr)) {
    return("value")
  }

  # Handle special arithmetic patterns
  desc <- .try_special_value_patterns(expr)
  if (!is.null(desc)) {
    return(desc)
  }

  # Handle simple variable names
  if (rlang::is_symbol(expr)) {
    var_name <- rlang::as_string(expr)
    return(.translate_variable(var_name))
  }

  # Handle arithmetic expressions
  if (rlang::is_call(expr)) {
    op <- rlang::call_name(expr)
    args <- rlang::call_args(expr)

    if (op == "/" && length(args) == 2) {
      num_var <- rlang::as_string(args[[1]])
      denom_var <- rlang::as_string(args[[2]])

      # Special cases for ratios
      if (num_var == "nS" && denom_var == "nG") {
        return("degree of sialylation per galactose")
      }
      if (num_var == "nG" && denom_var == "nA") {
        return("degree of galactosylation per antenna")
      }
      if (num_var == "nS" && denom_var == "nA") {
        return("degree of sialylation per antenna")
      }

      num_desc <- .translate_variable(num_var)
      denom_desc <- .translate_variable(denom_var)
      return(paste(num_desc, "per", denom_desc))
    }
  }

  # Fallback to generic description
  expr_text <- rlang::expr_text(expr)
  paste0("value of '", expr_text, "'")
}

.try_special_patterns <- function(expr) {
  if (!rlang::is_call(expr)) {
    if (rlang::is_symbol(expr)) {
      var_name <- rlang::as_string(expr)
      if (var_name == "B") {
        return("glycans with bisecting GlcNAc")
      }
    }
    return(NULL)
  }

  op <- rlang::call_name(expr)
  args <- rlang::call_args(expr)

  # Handle comparisons
  if (op == "==" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]

    if (rlang::is_symbol(var) && rlang::as_string(var) == "Tp") {
      if (rlang::is_string(val)) {
        type_val <- rlang::eval_bare(val)
        if (type_val == "highmannose") {
          return("high-mannose glycans")
        }
        if (type_val == "hybrid") {
          return("hybrid glycans")
        }
        if (type_val == "complex") return("complex glycans")
      }
    }

    if (rlang::is_symbol(var) && rlang::as_string(var) == "nA") {
      if (rlang::is_syntactic_literal(val)) {
        antenna_count <- rlang::eval_bare(val)
        if (antenna_count == 1) {
          return("mono-antennary glycans")
        }
        if (antenna_count == 2) {
          return("bi-antennary glycans")
        }
        if (antenna_count == 3) {
          return("tri-antennary glycans")
        }
        if (antenna_count == 4) return("tetra-antennary glycans")
      }
    }

    if (rlang::is_symbol(var) && rlang::is_syntactic_literal(val)) {
      var_name <- rlang::as_string(var)
      val_num <- rlang::eval_bare(val)

      if (var_name == "nS" && val_num %in% 1:4) {
        return(paste0(.count_prefix(val_num), "-sialylated glycans"))
      }
      if (var_name == "nF" && val_num %in% 1:4) {
        return(paste0(.count_prefix(val_num), "-fucosylated glycans"))
      }
    }
  }

  if (op == ">" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]

    if (rlang::is_symbol(var) && rlang::is_syntactic_literal(val)) {
      val_num <- rlang::eval_bare(val)
      var_name <- rlang::as_string(var)

      # Handle nX > 0 patterns (presence/absence)
      if (val_num == 0) {
        if (var_name == "nS") {
          return("sialylated glycans")
        }
        if (var_name == "nF") {
          return("fucosylated glycans")
        }
        if (var_name == "nG") {
          return("galactosylated glycans")
        }
        if (var_name == "nFc") {
          return("core-fucosylated glycans")
        }
        if (var_name == "nFa") return("arm-fucosylated glycans")
      }

      # Handle nA > n patterns (antenna count comparison)
      if (var_name == "nA") {
        return(paste0("glycans with more than ", val_num, " antennae"))
      }
    }
  }

  if (op == "<=" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]

    if (
      rlang::is_symbol(var) &&
        rlang::as_string(var) == "nA" &&
        rlang::is_syntactic_literal(val)
    ) {
      antenna_count <- rlang::eval_bare(val)
      return(paste0("glycans with at most ", antenna_count, " antennae"))
    }
  }

  if (op == ">=" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]

    if (
      rlang::is_symbol(var) &&
        rlang::as_string(var) == "nA" &&
        rlang::is_syntactic_literal(val)
    ) {
      antenna_count <- rlang::eval_bare(val)
      return(paste0("glycans with at least ", antenna_count, " antennae"))
    }
  }

  if (op == "==" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]

    if (
      rlang::is_symbol(var) &&
        rlang::is_syntactic_literal(val) &&
        rlang::eval_bare(val) == 0
    ) {
      var_name <- rlang::as_string(var)
      if (var_name == "nF") {
        return("non-fucosylated glycans")
      }
      if (var_name == "nFc") {
        return("non-core-fucosylated glycans")
      }
      if (var_name == "nFa") {
        return("non-arm-fucosylated glycans")
      }
      if (var_name == "nS") return("non-sialylated glycans")
    }
  }

  if (op == "!=" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]

    if (
      rlang::is_symbol(var) &&
        rlang::is_syntactic_literal(val) &&
        rlang::eval_bare(val) == 0
    ) {
      var_name <- rlang::as_string(var)
      if (var_name == "nF") {
        return("fucosylated glycans")
      }
      if (var_name == "nFc") {
        return("core-fucosylated glycans")
      }
      if (var_name == "nFa") {
        return("arm-fucosylated glycans")
      }
      if (var_name == "nS") return("sialylated glycans")
    }
  }

  return(NULL)
}

#' Translate Small Counts to Glycan Description Prefixes
#'
#' @param count A count from 1 to 4.
#'
#' @returns A character prefix such as "mono", "bi", "tri", or "tetra".
#' @noRd
.count_prefix <- function(count) {
  prefixes <- c(
    "1" = "mono",
    "2" = "bi",
    "3" = "tri",
    "4" = "tetra"
  )

  prefixes[[as.character(count)]]
}

.try_special_value_patterns <- function(expr) {
  # Already handled in .expr_to_value_description
  return(NULL)
}

.translate_variable <- function(var_name) {
  translations <- list(
    "nM" = "mannose count",
    "nA" = "antenna count",
    "nF" = "fucose count",
    "nFc" = "core fucose count",
    "nFa" = "arm fucose count",
    "nG" = "galactose count",
    "nS" = "sialic acid count",
    "nGt" = "terminal galactose count",
    "B" = "bisecting GlcNAc",
    "Tp" = "type"
  )

  if (var_name %in% names(translations)) {
    return(translations[[var_name]])
  }

  # Fallback to the original variable name
  return(var_name)
}

.try_merge_descriptions <- function(left_desc, right_desc, connector) {
  # Handle special patterns first
  if (connector == "and") {
    # Special handling for "glycans with X" patterns
    if (
      stringr::str_starts(right_desc, "glycans with ") &&
        stringr::str_ends(left_desc, " glycans")
    ) {
      left_adj <- stringr::str_replace(left_desc, " glycans$", "")
      right_feature <- stringr::str_replace(right_desc, "^glycans with ", "")
      return(paste0(left_adj, " glycans with ", right_feature))
    } else if (
      stringr::str_starts(left_desc, "glycans with ") &&
        stringr::str_ends(right_desc, " glycans")
    ) {
      right_adj <- stringr::str_replace(right_desc, " glycans$", "")
      left_feature <- stringr::str_replace(left_desc, "^glycans with ", "")
      return(paste0(right_adj, " glycans with ", left_feature))
    }
  }

  # Extract adjectives from glycan descriptions for merging
  left_pattern <- "^(.*?) glycans$"
  right_pattern <- "^(.*?) glycans$"

  if (
    stringr::str_detect(left_desc, left_pattern) &&
      stringr::str_detect(right_desc, right_pattern)
  ) {
    left_adj <- stringr::str_extract(left_desc, left_pattern) |>
      stringr::str_replace(left_pattern, "\\1")
    right_adj <- stringr::str_extract(right_desc, right_pattern) |>
      stringr::str_replace(right_pattern, "\\1")

    # Handle merging differently for "and" vs "or"
    if (
      !stringr::str_detect(left_adj, " and | or ") &&
        !stringr::str_detect(right_adj, " and | or ")
    ) {
      if (connector == "and") {
        # For "and": compound characteristics of the same glycans
        # Convert adjective endings to noun forms when using "with"
        right_adj_noun <- .adjective_to_noun(right_adj)
        if (!is.null(right_adj_noun)) {
          return(paste0(left_adj, " glycans with ", right_adj_noun))
        } else {
          # Use parallel adjective structure if noun conversion not available
          return(paste0(left_adj, " and ", right_adj, " glycans"))
        }
      } else if (connector == "or") {
        # For "or": different glycan groups
        return(paste0(left_adj, " or ", right_adj, " glycans"))
      }
    }
  }

  return(NULL)
}

.adjective_to_noun <- function(adj) {
  # Convert common adjective forms to noun forms for "with X" constructions
  conversions <- list(
    "fucosylated" = "fucosylation",
    "core-fucosylated" = "core-fucosylation",
    "arm-fucosylated" = "arm-fucosylation",
    "sialylated" = "sialylation",
    "galactosylated" = "galactosylation",
    "fucosylated" = "fucosylation",
    "mannosylated" = "mannosylation"
  )

  if (adj %in% names(conversions)) {
    return(conversions[[adj]])
  }

  # For other adjectives ending in -ylated, try removing -ylated and adding -ylation
  if (stringr::str_ends(adj, "ylated")) {
    base <- stringr::str_replace(adj, "ylated$", "ylation")
    return(base)
  }

  return(NULL)
}

# ----- AI-assisted explanations --------------------------------------

.explain_sys_prompt <- function(trait_type, custom_mp = NULL) {
  checkmate::assert_choice(
    trait_type,
    c("prop", "ratio", "wmean", "total", "wsum", "all")
  )

  # Build custom meta-properties section
  custom_mp_lines <- ""
  if (!is.null(custom_mp) && length(custom_mp) > 0) {
    custom_mp_lines <- paste0(
      "\nHere are the definitions of user-defined custom meta-properties:\n",
      paste0("- ", names(custom_mp), ": ", custom_mp, collapse = "\n"),
      "\n"
    )
  }

  prompt <- paste0(
    "You are a professional glycobiologist.\n",
    "Your task is to explain derived traits expressions in one sentence.\n",
    "Derived traits are defined by expressions using meta-properties.\n",
    "Here are the definitions of all built-in meta-properties:\n",
    "- Tp: glycan type (complex, hybrid, highmannose, pausimannose)\n",
    "- B: glycans with bisecting GlcNAc\n",
    "- nA: number of antennae\n",
    "- nF: number of fucoses\n",
    "- nFc: number of core fucoses\n",
    "- nFa: number of arm fucoses\n",
    "- nG: number of galactoses\n",
    "- nGt: number of terminal galactoses\n",
    "- nS: number of sialic acids\n",
    "- nM: number of mannoses\n",
    custom_mp_lines,
    "When encountering a meta-property that is not listed here (and not in the custom list), use it as is.\n",
    "A special case is custom meta-properties with prefix 'n' should be interpreted as the number of the corresponding unit, e.g. 'nE' -> 'number of E'\n",
    "For `MP > 0` patterns, use words like 'fucosylated', 'sialylated', etc.\n",
    "For nA related patterns, use words like 'bi-antennary', 'tri-antennary', 'tetra-antennary', etc.\n",
    "For multiple adjectives, the order should be: fucosylation/sialylation/galactosylation/mannosylation/bisection -> antennae count -> glycan type."
  )

  prop_examples <- paste(
    "Here are some examples:",
    'INPUT: prop(Tp == "highmannose")',
    "OUTPUT: Proportion of highmannose glycans among all glycans.",
    'INPUT: CA2 = prop(nA == 2, within = (Tp == "complex"))',
    "OUTPUT: Proportion of bi-antennary glycans within complex glycans.",
    "INPUT: TF = prop(nF > 0)",
    "OUTPUT: Proportion of fucosylated glycans among all glycans.",
    "INPUT: TFc = prop(nFc > 0)",
    "OUTPUT: Proportion of core-fucosylated glycans among all glycans.",
    "INPUT: TB = prop(B)",
    "OUTPUT: Proportion of glycans with bisecting GlcNAc among all glycans.",
    "INPUT: prop(nS > 0, within = (nFc > 0 & nA == 4))",
    "OUTPUT: Proportion of sialylated glycans within core-fucosylated tetra-antennary glycans.",
    "INPUT: TE = prop(E > 0)",
    "OUTPUT: Proportion of glycans with at least one E among all glycans.",
    sep = "\n"
  )

  ratio_examples <- paste(
    "Here are some examples:",
    'INPUT: ratio(Tp == "complex", Tp == "hybrid")',
    "OUTPUT: Ratio of complex glycans to hybrid glycans.",
    "INPUT: ratio(B, !B)",
    "OUTPUT: Ratio of glycans with bisecting GlcNAc to glycans without bisecting GlcNAc.",
    'INPUT: ratio(nFc > 0, nFc == 0, within = (Tp == "complex"))',
    "OUTPUT: Ratio of core-fucosylated glycans to non-core-fucosylated glycans within complex glycans.",
    "INPUT: EL = ratio(nE / nL)",
    "OUTPUT: Ratio of number of E to number of L.",
    sep = "\n"
  )

  wmean_examples <- paste(
    "Here are some examples:",
    "INPUT: wmean(nA)",
    "OUTPUT: Average number of antennae among all glycans.",
    "INPUT: wmean(nS / nA)",
    "OUTPUT: Average degree of sialylation per antenna among all glycans.",
    'INPUT: wmean(nS / nA, within = (Tp == "complex"))',
    "OUTPUT: Average degree of sialylation per antenna within complex glycans.",
    "INPUT: wmean(nS / nG, within = (nA == 4))",
    "OUTPUT: Average degree of sialylation per galactose within tetra-antennary glycans.",
    "INPUT: wmean(nE / nG)",
    "OUTPUT: Average number of E per galactose among all glycans.",
    sep = "\n"
  )

  total_examples <- paste(
    "Here are some examples:",
    'INPUT: total(Tp == "complex")',
    "OUTPUT: Total abundance of complex glycans.",
    'INPUT: total(Tp == "complex" & nA == 2)',
    "OUTPUT: Total abundance of bi-antennary complex glycans.",
    "INPUT: total(nE > 0)",
    "OUTPUT: Total abundance of glycans with at least one E.",
    sep = "\n"
  )

  wsum_examples <- paste(
    "Here are some examples:",
    "INPUT: wsum(nS)",
    "OUTPUT: Abundance-weighted sum of number of sialic acids among all glycans.",
    'INPUT: wsum(nA, within = (Tp == "complex"))',
    "OUTPUT: Abundance-weighted sum of number of antennae within complex glycans.",
    'INPUT: wsum(nS, within = (Tp == "complex"))',
    "OUTPUT: Abundance-weighted sum of number of sialic acids within complex glycans.",
    "INPUT: wsum(nE / nG)",
    "OUTPUT: Abundance-weighted sum of number of E per galactose among all glycans.",
    sep = "\n"
  )

  example_prompt <- switch(
    trait_type,
    "prop" = prop_examples,
    "ratio" = ratio_examples,
    "wmean" = wmean_examples,
    "total" = total_examples,
    "wsum" = wsum_examples,
    "all" = paste(
      prop_examples,
      ratio_examples,
      wmean_examples,
      total_examples,
      wsum_examples,
      sep = "\n"
    )
  )
  paste(prompt, example_prompt, sep = "\n")
}

.explain_with_ai <- function(
  trait_fn,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  trait_str <- rlang::expr_text(trait_fn)
  str_type <- stringr::str_extract(trait_str, "prop|ratio|wmean|total|wsum")
  system_prompt <- .explain_sys_prompt(str_type, custom_mp)
  user_prompt <- paste0("INPUT: ", trait_str, "\nOUTPUT: ")
  .ask_ai(
    system_prompt,
    user_prompt,
    api_key = api_key,
    model = model,
    provider = provider,
    base_url = base_url
  )
}

#' Convert a Derived Trait Function to Its Factory Formula
#'
#' @param trait_fn A derived trait function.
#' @returns A one-line trait factory formula.
#' @noRd
.trait_to_formula <- function(trait_fn) {
  trait_class <- class(trait_fn)
  trait_type <- sub(
    "^glydet_",
    "",
    trait_class[grepl("^glydet_", trait_class)]
  )[[1]]
  within <- attr(trait_fn, "within")
  within_arg <- if (is.null(within)) {
    ""
  } else {
    paste0(", within = (", rlang::expr_text(within), ")")
  }

  switch(
    trait_type,
    prop = paste0(
      "prop(",
      rlang::expr_text(attr(trait_fn, "cond")),
      within_arg,
      ")"
    ),
    ratio = paste0(
      "ratio(",
      rlang::expr_text(attr(trait_fn, "num_cond")),
      ", ",
      rlang::expr_text(attr(trait_fn, "denom_cond")),
      within_arg,
      ")"
    ),
    wmean = paste0(
      "wmean(",
      rlang::expr_text(attr(trait_fn, "val")),
      within_arg,
      ")"
    ),
    total = paste0("total(", rlang::expr_text(attr(trait_fn, "cond")), ")"),
    wsum = paste0(
      "wsum(",
      rlang::expr_text(attr(trait_fn, "val")),
      within_arg,
      ")"
    ),
    cli::cli_abort("Object must be a supported glydet trait function.")
  )
}
