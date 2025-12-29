#' Use a Large Language Model (LLM) to create a derived trait function
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' This function allows you to create a derived trait function using natural language.
#' Note that LLMs can be unreliable, so the result should be verified manually.
#' If the description is not clear, an error will be raised.
#' Try to read the descriptions of built-in traits to get ideas.
#' Currently, only `prop()`, `ratio()`, and `wmean()` are supported.
#' To use this feature, you need to install the `ellmer` package.
#' You also need to provide an API key for the DeepSeek chat model.
#' Please set the environment variable `DEEPSEEK_API_KEY` to your API key using `Sys.setenv()`.
#' You can obtain an API key from https://platform.deepseek.com.
#'
#' @param description A description of the trait in natural language.
#' @param custom_mp A named character vector of custom meta-properties.
#'   The names are the meta-property names, and the values are in the format
#'   "(type) description". For example:
#'   `c(nE = "(integer) number of a2,6-linked sialic acids")`.
#'   These custom meta-properties will be available for the LLM to use.
#' @param max_retries Maximum number of reflection retries when the AI-generated
#'   formula's explanation doesn't match the original description. Default is 2.
#' @param verbose Whether to print verbose output. Default is FALSE.
#'   This is useful for inspecting how LLMs generate trait functions.
#'
#' @returns A derived trait function.
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
#' @export
make_trait <- function(description, custom_mp = NULL, max_retries = 2, verbose = FALSE) {
  checkmate::assert_string(description)
  checkmate::assert_character(custom_mp, names = "named", null.ok = TRUE)
  checkmate::assert_count(max_retries)
  description <- stringr::str_trim(description)
  rlang::check_installed("ellmer")

  api_key <- .get_api_key()
  system_prompt <- .make_trait_sys_prompt(description, custom_mp)

  chat <- ellmer::chat_deepseek(
    system_prompt = system_prompt,
    model = "deepseek-chat",
    echo = "none",
    credentials = function() api_key
  )

  # Initial prompt
  current_prompt <- paste0("INPUT: ", description, "\nOUTPUT: ")

  for (i in 0:max_retries) {
    if (i > 0 && verbose) {
      cli::cli_alert_info("Attempt {i}/{max_retries}: Retrying with feedback...")
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
          result$error, "\n",
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
    reflection_result <- .check_trait_consistency(description, trait_fn, custom_mp)

    if (reflection_result$consistent) {
      return(trait_fn)
    }

    # Explanation doesn't match, retry with feedback
    if (i < max_retries) {
      current_prompt <- paste0(
        "The formula you generated was: ", result$formula, "\n",
        "The AI explanation of this formula is: ", reflection_result$explanation, "\n",
        "This does not match the original intent: ", description, "\n",
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

.check_trait_consistency <- function(description, trait_fn, custom_mp = NULL) {
  # Get AI explanation of the generated trait
  explanation <- tryCatch(
    explain_trait(trait_fn, use_ai = TRUE, custom_mp = custom_mp),
    error = function(e) NULL
  )

  if (is.null(explanation)) {
    # If explanation fails, assume inconsistent
    return(list(consistent = FALSE, explanation = "Failed to get explanation."))
  }

  # Ask AI to judge if the explanation matches the description
  is_consistent <- .ask_ai_consistency(description, explanation)

  list(consistent = is_consistent, explanation = explanation)
}

.ask_ai_consistency <- function(description, explanation) {
  system_prompt <- paste(
    "You are a professional glycobiologist.",
    "Your task is to judge if two statements about glycan traits are semantically equivalent.",
    "Answer only 'YES' if they mean the same thing, or 'NO' if they are different.",
    "Be lenient: minor wording differences are acceptable if the meaning is the same.",
    sep = "\n"
  )
  user_prompt <- paste0(
    "Original description: ", description, "\n",
    "Generated explanation: ", explanation, "\n",
    "Are these two statements semantically equivalent? Answer YES or NO only."
  )

  response <- .ask_ai(system_prompt, user_prompt)
  stringr::str_detect(toupper(response), "YES")
}

.make_trait_sys_prompt <- function(description, custom_mp = NULL) {
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