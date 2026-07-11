#' Use a Large Language Model (LLM) to Create Derived Trait Functions
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' These functions create derived trait functions from natural-language
#' descriptions. `make_trait()` creates one trait and checks its consistency
#' with the requested description. `make_traits()` creates and validates a list
#' of traits in batched requests, reducing repeated prompt tokens. For batch
#' creation, descriptions that cannot be understood, produce an invalid formula,
#' or do not match the generated trait are returned as `NA` with a warning.
#' LLM-generated traits should always be verified manually.
#' Try to read the descriptions of built-in traits to get ideas.
#' Currently, only `prop()`, `ratio()`, and `wmean()` are supported.
#' To use this feature, you need to install the `ellmer` package.
#' DeepSeek is used by default for backward compatibility. Other `ellmer`
#' providers can be selected with `provider`, `model`, and provider-specific
#' API key configuration.
#'
#' @section Batch multi-agent workflow:
#' `make_traits()` uses one batch writer to generate formulas for all active
#' descriptions, one batch explainer to describe the generated formulas, and one
#' batch evaluator to compare those explanations with the original descriptions.
#' Successful positions are retained. Invalid or mismatched positions are sent
#' back to the writer with their validation error or generated explanation, and
#' only those positions are regenerated. This continues until all positions pass
#' or `max_retries` is reached.
#'
#' @param description A description of the trait in natural language.
#' @param descriptions A character vector of trait descriptions.
#' @param custom_mp A named character vector of custom meta-properties.
#'   The names are the meta-property names, and the values are in the format
#'   "(type) description". For example:
#'   `c(nE = "(integer) number of a2,6-linked sialic acids")`.
#'   These custom meta-properties will be available for the LLM to use.
#'   Note that defining the meta-properties here is not enough for you to use them.
#'   You need to define corresponding meta-property functions or specifying meta-property columns.
#'   For more information about custom meta-properties, see the vignette
#'   [Custom Meta-Properties](https://glycoverse.github.io/glydet/articles/custom-traits.html#using-make_trait).
#' @param max_retries Maximum number of retries after an invalid formula or an
#'   explanation that doesn't match the original description. In batch mode,
#'   only unresolved descriptions are retried. Default is 2.
#' @param verbose Whether to print verbose output. Default is FALSE.
#'   This is useful for inspecting how LLMs generate trait functions.
#' @param provider AI provider passed to `ellmer`. One of "deepseek",
#'   "openai", "anthropic", "gemini", "openrouter", or "openai_compatible".
#'   "google_gemini" is accepted as an alias for "gemini". Defaults to
#'   `getOption("glydet.ai_provider", "deepseek")`.
#' @param model Model to use. Defaults to `getOption("glydet.ai_model")`,
#'   or "deepseek-chat" for DeepSeek and the provider default for other providers.
#' @param api_key API key for the selected provider. If `NULL`, the provider
#'   specific environment variable is used. Defaults to
#'   `getOption("glydet.ai_api_key")`.
#' @param base_url Optional base URL for custom or OpenAI-compatible endpoints.
#'   Defaults to `getOption("glydet.ai_base_url")`.
#'
#' @returns `make_trait()` returns a derived trait function. `make_traits()`
#'   returns a list of derived trait functions with input names preserved; entries
#'   that cannot be created are `NA`.
#'
#' @examples
#' # Sys.setenv(DEEPSEEK_API_KEY = "your_api_key")
#' # my_traits <- list(
#' #   nS = make_trait("the average number of sialic acids"),
#' #   nG = make_trait("the average number of galactoses")
#' # )
#'
#' # The trait function can then be used in `derive_traits()`:
#' # derive_traits(exp, trait_fns = my_traits)
#'
#' \dontrun{
#' make_traits(c(
#'   sialylated = "proportion of sialylated glycans",
#'   galactose = "average number of galactoses"
#' ))
#' }
#'
#' @export
make_trait <- function(
  description,
  custom_mp = NULL,
  max_retries = 2,
  verbose = FALSE,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  checkmate::assert_string(description)
  checkmate::assert_character(custom_mp, names = "named", null.ok = TRUE)
  checkmate::assert_count(max_retries)
  provider <- .normalize_ai_provider(provider)
  model <- .normalize_optional_ai_string(model)
  api_key <- .normalize_optional_ai_string(api_key)
  base_url <- .normalize_optional_ai_string(base_url)
  checkmate::assert_string(model, null.ok = TRUE)
  checkmate::assert_string(api_key, null.ok = TRUE)
  checkmate::assert_string(base_url, null.ok = TRUE)
  description <- stringr::str_trim(description)
  rlang::check_installed("ellmer")

  model <- .resolve_ai_model(provider, model)
  api_key <- .get_api_key(provider = provider, api_key = api_key)
  system_prompt <- .make_trait_sys_prompt(custom_mp = custom_mp)

  chat <- .create_ai_chat(
    system_prompt = system_prompt,
    api_key = api_key,
    provider = provider,
    model = model,
    base_url = base_url
  )

  # Initial prompt
  current_prompt <- paste0("INPUT: ", description, "\nOUTPUT: ")

  for (i in 0:max_retries) {
    if (i > 0 && verbose) {
      cli::cli_alert_info(
        "Attempt {i}/{max_retries}: Retrying with feedback..."
      )
    }

    # Call AI to generate trait formula
    output <- as.character(chat$chat(current_prompt))
    if (verbose) {
      cli::cli_alert_info("AI generated formula: {cli::col_grey(output)}")
    }
    result <- .process_trait_response(output)

    if (!result$valid) {
      # Format error, retry
      if (i < max_retries) {
        current_prompt <- paste0(
          "The previous formula was invalid:\n",
          result$error,
          "\n",
          "Please fix the formula and return only the corrected expression."
        )
        next
      } else {
        cli::cli_abort(c(
          "Failed to create a derived trait function after {max_retries} retries.",
          "x" = "Last error: {result$error}",
          "i" = "Please try again with a different description."
        ))
      }
    }

    # Valid formula, now do reflection check
    trait_fn <- result$trait_fn
    reflection_result <- .check_trait_consistency(
      description,
      trait_fn,
      custom_mp,
      api_key = api_key,
      model = model,
      provider = provider,
      base_url = base_url
    )

    if (reflection_result$consistent) {
      return(trait_fn)
    }

    # Explanation doesn't match, retry with feedback
    if (i < max_retries) {
      current_prompt <- paste0(
        "The formula you generated was: ",
        result$formula,
        "\n",
        "The AI explanation of this formula is: ",
        reflection_result$explanation,
        "\n",
        "This does not match the original intent: ",
        description,
        "\n",
        "Please generate a corrected formula. Return only the expression."
      )
    } else {
      cli::cli_abort(c(
        "Failed to create a consistent derived trait function after {max_retries} retries.",
        "x" = "Generated formula: {result$formula}",
        "x" = "AI explanation: {reflection_result$explanation}",
        "x" = "Original description: {description}",
        "i" = "Please try again with a different description."
      ))
    }
  }
}

#' @rdname make_trait
#' @export
make_traits <- function(
  descriptions,
  custom_mp = NULL,
  max_retries = 2,
  verbose = FALSE,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  checkmate::assert_character(descriptions, any.missing = TRUE)
  checkmate::assert_character(custom_mp, names = "named", null.ok = TRUE)
  checkmate::assert_count(max_retries)
  checkmate::assert_flag(verbose)
  provider <- .normalize_ai_provider(provider)
  model <- .normalize_optional_ai_string(model)
  api_key <- .normalize_optional_ai_string(api_key)
  base_url <- .normalize_optional_ai_string(base_url)
  checkmate::assert_string(model, null.ok = TRUE)
  checkmate::assert_string(api_key, null.ok = TRUE)
  checkmate::assert_string(base_url, null.ok = TRUE)

  trait_fns <- rep(list(NA), length(descriptions))
  names(trait_fns) <- names(descriptions)
  descriptions <- stringr::str_squish(descriptions)
  valid <- !is.na(descriptions) & nzchar(descriptions)
  positions <- which(valid)

  if (length(positions) > 0) {
    api_key <- .get_api_key(provider = provider, api_key = api_key)
    system_prompt <- paste0(
      .make_trait_sys_prompt(custom_mp = custom_mp),
      "\nCreate every requested trait independently. Return exactly one line per ",
      "description in the form `index<TAB>formula`. Return `index<TAB><INVALID>` ",
      "when a description cannot be understood."
    )
    chat <- .create_ai_chat(
      system_prompt = system_prompt,
      api_key = api_key,
      provider = provider,
      model = model,
      base_url = base_url
    )
    current_prompt <- paste0(
      "INPUTS:\n",
      paste0(positions, "\t", descriptions[positions], collapse = "\n"),
      "\nOUTPUT:"
    )
    active_positions <- positions
    feedback <- rep(NA_character_, length(descriptions))

    for (attempt in 0:max_retries) {
      if (attempt > 0 && verbose) {
        cli::cli_alert_info(
          "Batch retry {attempt}/{max_retries}: regenerating {length(active_positions)} trait(s)."
        )
      }

      output <- as.character(chat$chat(current_prompt))
      formulas <- .parse_batch_response(output, length(descriptions))
      candidate_positions <- integer()
      candidate_formulas <- character()
      candidate_fns <- list()

      for (position in active_positions) {
        formula <- formulas[[position]]
        if (is.na(formula)) {
          feedback[[position]] <- paste(
            "The previous response did not contain a valid formula.",
            "Generate a formula for this description."
          )
          next
        }

        result <- .process_trait_response(formula)
        if (!result$valid) {
          feedback[[position]] <- paste0(
            "The previous formula was invalid: ",
            result$error,
            " Generate a corrected formula."
          )
          next
        }

        candidate_positions <- c(candidate_positions, position)
        candidate_formulas <- c(candidate_formulas, result$formula)
        candidate_fns <- c(candidate_fns, list(result$trait_fn))
      }

      if (length(candidate_positions) > 0) {
        consistency <- .check_traits_consistency(
          descriptions[candidate_positions],
          candidate_fns,
          custom_mp = custom_mp,
          api_key = api_key,
          model = model,
          provider = provider,
          base_url = base_url
        )

        for (index in seq_along(candidate_positions)) {
          position <- candidate_positions[[index]]
          if (isTRUE(consistency$consistent[[index]])) {
            trait_fns[[position]] <- candidate_fns[[index]]
            next
          }

          explanation <- consistency$explanations[[index]]
          if (is.na(explanation)) {
            explanation <- "The generated formula could not be explained."
          }
          feedback[[position]] <- paste0(
            "The previous formula was: ",
            candidate_formulas[[index]],
            "\nThe generated explanation was: ",
            explanation,
            "\nThis does not match the original description. Generate a corrected formula."
          )
        }
      }

      active_positions <- active_positions[
        !purrr::map_lgl(trait_fns[active_positions], is.function)
      ]
      if (length(active_positions) == 0 || attempt == max_retries) {
        break
      }

      current_prompt <- paste0(
        "RETRY INPUTS:\n",
        paste0(
          active_positions,
          "\tOriginal description: ",
          descriptions[active_positions],
          "\n\tFeedback: ",
          feedback[active_positions],
          collapse = "\n"
        ),
        "\nReturn exactly one corrected formula per input using the original ",
        "index in the form `index<TAB>formula`.\nOUTPUT:"
      )
    }
  }

  invalid_positions <- which(purrr::map_lgl(
    trait_fns,
    ~ is.atomic(.x) && length(.x) == 1 && is.na(.x)
  ))
  if (length(invalid_positions) > 0) {
    cli::cli_warn(
      "Could not make trait(s) at position(s): {paste(invalid_positions, collapse = ', ')}."
    )
  }

  trait_fns
}

.process_trait_response <- function(output) {
  # Clean up output
  output <- stringr::str_trim(output)
  output <- stringr::str_remove_all(output, "`")

  if (stringr::str_detect(output, "<INVALID>")) {
    return(list(valid = FALSE, error = "AI returned <INVALID> tag."))
  }

  if (!stringr::str_detect(output, "^(prop|ratio|wmean)\\((.*)\\)$")) {
    return(list(valid = FALSE, error = paste0("Invalid format: ", output)))
  }

  tryCatch(
    {
      expr <- rlang::parse_expr(output)

      # Evaluate in a clean environment to avoid capturing temporary objects
      clean_env <- rlang::new_environment(parent = asNamespace("glydet"))
      trait_fn <- eval(expr, envir = clean_env)

      # Basic validation: try explain without AI
      explain <- explain_trait(trait_fn)

      list(valid = TRUE, trait_fn = trait_fn, formula = output)
    },
    error = function(e) {
      list(valid = FALSE, error = paste0("Parse/eval error: ", e$message))
    }
  )
}

.check_trait_consistency <- function(
  description,
  trait_fn,
  custom_mp = NULL,
  api_key = getOption("glydet.ai_api_key", NULL),
  model = getOption("glydet.ai_model", NULL),
  provider = getOption("glydet.ai_provider", "deepseek"),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  # Get AI explanation of the generated trait
  explanation <- tryCatch(
    explain_trait(
      trait_fn,
      use_ai = TRUE,
      custom_mp = custom_mp,
      api_key = api_key,
      model = model,
      provider = provider,
      base_url = base_url
    ),
    error = function(e) NULL
  )

  if (is.null(explanation)) {
    # If explanation fails, assume inconsistent
    return(list(consistent = FALSE, explanation = "Failed to get explanation."))
  }

  # Ask AI to judge if the explanation matches the description
  is_consistent <- .ask_ai_consistency(
    description,
    explanation,
    api_key = api_key,
    model = model,
    provider = provider,
    base_url = base_url
  )

  list(consistent = is_consistent, explanation = explanation)
}

#' Check Multiple Generated Traits Against Their Descriptions
#'
#' @param descriptions Character vector of original trait descriptions.
#' @param trait_fns List of generated trait functions.
#' @param custom_mp A named character vector of custom meta-properties.
#' @param api_key API key for the selected provider.
#' @param model Model to use.
#' @param provider AI provider passed to `ellmer`.
#' @param base_url Optional base URL for custom or OpenAI-compatible endpoints.
#' @returns A list containing a logical consistency vector and the generated
#'   explanations. Unparseable responses are `NA`.
#' @noRd
.check_traits_consistency <- function(
  descriptions,
  trait_fns,
  custom_mp = NULL,
  api_key = getOption("glydet.ai_api_key", NULL),
  model = getOption("glydet.ai_model", NULL),
  provider = getOption("glydet.ai_provider", "deepseek"),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  explanations <- suppressWarnings(
    explain_traits(
      trait_fns,
      use_ai = TRUE,
      custom_mp = custom_mp,
      api_key = api_key,
      model = model,
      provider = provider,
      base_url = base_url
    )
  )
  positions <- which(!is.na(explanations))
  consistent <- rep(NA, length(trait_fns))

  if (length(positions) == 0) {
    return(list(consistent = consistent, explanations = explanations))
  }

  custom_mp_lines <- ""
  if (!is.null(custom_mp) && length(custom_mp) > 0) {
    custom_mp_lines <- paste0(
      "Custom meta-properties:\n",
      paste0("- ", names(custom_mp), ": ", custom_mp, collapse = "\n"),
      "\n"
    )
  }
  system_prompt <- paste0(
    "You are a professional glycobiologist.\n",
    "Your task is to judge whether each generated explanation is semantically ",
    "equivalent to its original trait description.\n",
    custom_mp_lines,
    "Return exactly one line per input in the form `index<TAB>YES` or ",
    "`index<TAB>NO`. Use YES only when the meanings are the same; minor wording ",
    "differences are acceptable."
  )
  user_prompt <- paste0(
    "INPUTS:\n",
    paste0(
      positions,
      "\tDescription: ",
      descriptions[positions],
      "\n\tGenerated explanation: ",
      explanations[positions],
      collapse = "\n"
    ),
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
  answers <- .parse_batch_response(response, length(trait_fns))
  consistent[positions] <- toupper(answers[positions]) == "YES"
  list(consistent = consistent, explanations = explanations)
}

.ask_ai_consistency <- function(
  description,
  explanation,
  api_key = getOption("glydet.ai_api_key", NULL),
  model = getOption("glydet.ai_model", NULL),
  provider = getOption("glydet.ai_provider", "deepseek"),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  system_prompt <- paste(
    "You are a professional glycobiologist.",
    "Your task is to judge if two statements about glycan traits are semantically equivalent.",
    "Answer only 'YES' if they mean the same thing, or 'NO' if they are different.",
    "Be lenient: minor wording differences are acceptable if the meaning is the same.",
    sep = "\n"
  )
  user_prompt <- paste0(
    "Original description: ",
    description,
    "\n",
    "Generated explanation: ",
    explanation,
    "\n",
    "Are these two statements semantically equivalent? Answer YES or NO only."
  )

  response <- .ask_ai(
    system_prompt,
    user_prompt,
    api_key = api_key,
    model = model,
    provider = provider,
    base_url = base_url
  )
  stringr::str_detect(toupper(response), "YES")
}

#' Build the System Prompt for Making Derived Traits
#'
#' @param description Deprecated and unused. Kept for internal compatibility.
#' @param custom_mp A named character vector of custom meta-properties.
#' @returns A system prompt for an AI chat.
#' @noRd
.make_trait_sys_prompt <- function(description = NULL, custom_mp = NULL) {
  # Build custom meta-properties section
  custom_mp_lines <- ""
  if (!is.null(custom_mp) && length(custom_mp) > 0) {
    custom_mp_lines <- paste0(
      "Here are the definitions of user-defined custom meta-properties:\n",
      paste0("- ", names(custom_mp), ": ", custom_mp, collapse = "\n"),
      "\n"
    )
  }

  paste0(
    "You are a professional glycobiologist.\n",
    "Your task is to create a derived trait function using natural language.\n",
    "Derived traits are defined by expressions using meta-properties.\n",
    "Here are the definitions of all built-in meta-properties:\n",
    "- Tp: (string) glycan type (complex, hybrid, highmannose, pausimannose)\n",
    "- B: (logical) glycans with bisecting GlcNAc\n",
    "- nA: (numeric) number of antennae\n",
    "- nF: (numeric) number of fucoses\n",
    "- nFc: (numeric) number of core fucoses\n",
    "- nFa: (numeric) number of arm fucoses\n",
    "- nG: (numeric) number of galactoses\n",
    "- nGt: (numeric) number of terminal galactoses (galactoses without further modifications)",
    "- nS: (numeric) number of sialic acids\n",
    "- nM: (numeric) number of mannoses\n",
    custom_mp_lines,
    "When encountering a meta-property that is not listed here (and not in the custom list), output an <INVALID> tag.\n",
    "To create a derived trait function, you need to choose one of the following factories:\n",
    "- prop(): for calculating the abundance proportion of a specific glycan subset\n",
    "- ratio(): for calculating the abundance ratio between two glycan subsets\n",
    "- wmean(): for calculating the weighted mean of a quantitative property, weighted by glycan abundances\n",
    "prop() accepts a logical expression as the first argument.\n",
    "ratio() accepts two logical expressions as the first two arguments.\n",
    "wmean() accepts a numeric expression as the first argument.\n",
    "Each function has a `within` parameter to restrict the calculation to a specific glycan subset, which is a logical expression.\n",
    "Here are some examples:\n",
    "INPUT: proportion of fucosylated glycans within mono-antennary glycans\n",
    "OUTPUT: prop(nF > 0, within = (nA == 1))\n",
    "INPUT: proportion of fucosylated glycans within mono-antennary glycans\n",
    "OUTPUT: prop(nF > 0, within = (nA == 1))\n",
    "INPUT: Proportion of core-fucosylated glycans within bi-antennary glycans\n",
    "OUTPUT: prop(nFc > 0, within = (nA == 2))\n",
    "INPUT: Proportion of arm-fucosylated glycans within asialylated tri-antennary glycans\n",
    "OUTPUT: prop(nFa > 0, within = (nA == 3 & nS == 0))\n",
    "INPUT: Propotion of bisecting glycans within a-core-fucosylated tetra-antennary glycans\n",
    "OUTPUT: prop(B, within = (nA == 4 & nFc == 0))\n",
    "INPUT: Ratio of complex glycans to hybrid glycans\n",
    'OUTPUT: ratio(Tp == "complex", Tp == "hybrid")\n',
    "INPUT: Ratio of bisecting glycans to unbisecting glycans within bi-antennary glycans\n",
    "OUTPUT: ratio(B, !B, within = (nA == 2))\n",
    "INPUT: Ratio of core-fucosylated glycans to non-core-fucosylated glycans within complex glycans\n",
    'OUTPUT: ratio(nFc > 0, nFc == 0, within = (Tp == "complex"))\n',
    "INPUT: Average degree of sialylation per antenna within bi-antennary glycans\n",
    "OUTPUT: wmean(nS / nA, within = (nA == 2))\n",
    "INPUT: Average degree of galactosylation per antenna within tri-antennary glycans\n",
    "OUTPUT: wmean(nG / nA, within = (nA == 3))\n",
    "INPUT: Average numbers of antennae within complex glycans\n",
    'OUTPUT: wmean(nA, within = (Tp == "complex"))\n',
    "INPUT: Average number of mannoses within highmannose glycans\n",
    'OUTPUT: wmean(nM, within = (Tp == "highmannose"))\n',
    "You need to decide:\n",
    "- Which factory to use\n",
    "- The logical expression to use as the first (and second) argument(s) of the factory\n",
    "- The logical expression to use as the `within` parameter of the factory\n",
    "- The complete derived trait function\n",
    "If you cannot make a trait function, output an <INVALID> tag."
  )
}
