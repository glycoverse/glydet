#' Explain a Derived Trait
#'
#' @description
#' This function provides a human-readable English explanation of what a derived trait represents.
#' It works with trait functions created by the trait factories [prop()], [ratio()], and [wmean()],
#' and works best with traits defined by built-in meta-properties.
#'
#' @param trait_fn A derived trait function created by one of the trait factories.
#'
#' @returns A character string containing a concise English explanation of the trait.
#'
#' @examples
#' # Explain built-in traits
#' explain_trait(basic_traits()$TM)
#' explain_trait(basic_traits()$SG)
#'
#' # Explain custom traits
#' explain_trait(prop(nFc > 0))
#' explain_trait(prop(nFc > 0, within = (T == "complex")))
#' explain_trait(ratio(T == "complex", T == "hybrid"))
#' explain_trait(wmean(nA, within = (T == "complex")))
#' explain_trait(wmean(nS / nG, within = nA == 4 & nFc > 0))
#'
#' @export
explain_trait <- function(trait_fn) {
  UseMethod("explain_trait")
}

#' @export
explain_trait.default <- function(trait_fn) {
  cli::cli_abort(c(
    "Object must be a glydet trait function.",
    "x" = "Got an object of class {.cls {class(trait_fn)}}.",
    "i" = "Use trait factories like {.fn prop}, {.fn ratio}, or {.fn wmean} to create trait functions."
  ))
}

#' @export
explain_trait.glydet_prop <- function(trait_fn) {
  cond_expr <- attr(trait_fn, "cond")
  within_expr <- attr(trait_fn, "within")
  
  cond_desc <- .expr_to_description(cond_expr)
  scope_desc <- .within_to_scope(within_expr)
  
  paste0("Proportion of ", cond_desc, " ", scope_desc, ".")
}

#' @export
explain_trait.glydet_ratio <- function(trait_fn) {
  num_cond_expr <- attr(trait_fn, "num_cond")
  denom_cond_expr <- attr(trait_fn, "denom_cond")
  within_expr <- attr(trait_fn, "within")
  
  num_desc <- .expr_to_description(num_cond_expr)
  denom_desc <- .expr_to_description(denom_cond_expr)
  scope_desc <- .within_to_scope(within_expr)
  
  paste0("Ratio of ", num_desc, " to ", denom_desc, " ", scope_desc, ".")
}

#' @export
explain_trait.glydet_wmean <- function(trait_fn) {
  val_expr <- attr(trait_fn, "val")
  within_expr <- attr(trait_fn, "within")
  
  val_desc <- .expr_to_value_description(val_expr)
  scope_desc <- .within_to_scope(within_expr)
  
  paste0("Abundance-weighted mean of ", val_desc, " ", scope_desc, ".")
}

# Internal helper functions

.within_to_scope <- function(within_expr) {
  if (is.null(within_expr)) {
    return("among all glycans")
  }
  
  # Remove outer parentheses if present (consistent with print methods)
  within_text <- rlang::expr_text(within_expr)
  if (stringr::str_starts(within_text, stringr::fixed("(")) && 
      stringr::str_ends(within_text, stringr::fixed(")"))) {
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
        return(stringr::str_replace(arg_desc, "glycans with", "glycans without"))
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
        if (type_val == "highmannose") return("high-mannose glycans")
        if (type_val == "hybrid") return("hybrid glycans")
        if (type_val == "complex") return("complex glycans")
      }
    }
    
    if (rlang::is_symbol(var) && rlang::as_string(var) == "nA") {
      if (rlang::is_syntactic_literal(val)) {
        antenna_count <- rlang::eval_bare(val)
        if (antenna_count == 2) return("bi-antennary glycans")
        if (antenna_count == 3) return("tri-antennary glycans")
        if (antenna_count == 4) return("tetra-antennary glycans")
      }
    }
  }
  
  if (op == ">" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]
    
    if (rlang::is_symbol(var) && rlang::is_syntactic_literal(val) && rlang::eval_bare(val) == 0) {
      var_name <- rlang::as_string(var)
      if (var_name == "nS") return("sialylated glycans")
      if (var_name == "nF") return("fucosylated glycans")
      if (var_name == "nFc") return("core-fucosylated glycans")
      if (var_name == "nFa") return("arm-fucosylated glycans")
    }
  }
  
  if (op == "==" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]
    
    if (rlang::is_symbol(var) && rlang::is_syntactic_literal(val) && rlang::eval_bare(val) == 0) {
      var_name <- rlang::as_string(var)
      if (var_name == "nF") return("non-fucosylated glycans")
      if (var_name == "nFc") return("non-core-fucosylated glycans")
      if (var_name == "nFa") return("non-arm-fucosylated glycans")
      if (var_name == "nS") return("non-sialylated glycans")
    }
  }
  
  if (op == "!=" && length(args) == 2) {
    var <- args[[1]]
    val <- args[[2]]
    
    if (rlang::is_symbol(var) && rlang::is_syntactic_literal(val) && rlang::eval_bare(val) == 0) {
      var_name <- rlang::as_string(var)
      if (var_name == "nF") return("fucosylated glycans")
      if (var_name == "nFc") return("core-fucosylated glycans")
      if (var_name == "nFa") return("arm-fucosylated glycans")
      if (var_name == "nS") return("sialylated glycans")
    }
  }
  
  return(NULL)
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
    "nGt" = "galactose count",
    "nT" = "terminal count",
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
    if (stringr::str_starts(right_desc, "glycans with ") && 
        stringr::str_ends(left_desc, " glycans")) {
      left_adj <- stringr::str_replace(left_desc, " glycans$", "")
      right_feature <- stringr::str_replace(right_desc, "^glycans with ", "")
      return(paste0(left_adj, " glycans with ", right_feature))
    } else if (stringr::str_starts(left_desc, "glycans with ") && 
               stringr::str_ends(right_desc, " glycans")) {
      right_adj <- stringr::str_replace(right_desc, " glycans$", "")
      left_feature <- stringr::str_replace(left_desc, "^glycans with ", "")
      return(paste0(right_adj, " glycans with ", left_feature))
    }
  }
  
  # Extract adjectives from glycan descriptions for merging
  left_pattern <- "^(.*?) glycans$"
  right_pattern <- "^(.*?) glycans$"
  
  if (stringr::str_detect(left_desc, left_pattern) && 
      stringr::str_detect(right_desc, right_pattern)) {
    
    left_adj <- stringr::str_extract(left_desc, left_pattern) |>
      stringr::str_replace(left_pattern, "\\1")
    right_adj <- stringr::str_extract(right_desc, right_pattern) |>
      stringr::str_replace(right_pattern, "\\1")
    
    # Handle merging differently for "and" vs "or"
    if (!stringr::str_detect(left_adj, " and | or ") && 
        !stringr::str_detect(right_adj, " and | or ")) {
      
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
