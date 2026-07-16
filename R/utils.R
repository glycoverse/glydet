.process_glycans <- function(glycans) {
  if (glyrepr::is_glycan_structure(glycans)) {
    glycans
  } else if (is.character(glycans)) {
    glyparse::auto_parse(glycans)
  } else {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} object or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(glycans)}}."
    ))
  }
}

#' Check a glydet data container
#'
#' @param exp A `GlycomicSE`, `GlycoproteomicSE`, or legacy glyexp data
#'   container.
#'
#' @returns `NULL` invisibly, or an error if `exp` is unsupported.
#' @noRd
.assert_data_container <- function(exp) {
  if (
    inherits(exp, "glyexp_experiment") ||
      glyexp::is_glycomic_se(exp) ||
      glyexp::is_glycoproteomic_se(exp)
  ) {
    return(invisible(NULL))
  }

  cli::cli_abort(
    paste0(
      "{.arg exp} must be a {.cls glyexp_experiment}, ",
      "{.cls GlycomicSE}, or {.cls GlycoproteomicSE} object."
    )
  )
}

#' Convert a supported data container to a glyco SummarizedExperiment
#'
#' @inheritParams .assert_data_container
#'
#' @returns A `GlycomicSE` or `GlycoproteomicSE` object.
#' @noRd
.as_glyco_se <- function(exp) {
  if (glyexp::is_glycomic_se(exp) || glyexp::is_glycoproteomic_se(exp)) {
    return(exp)
  }

  exp_type <- S4Vectors::metadata(glyexp::as_se(exp))[["exp_type"]]
  switch(
    exp_type,
    glycomics = glyexp::as_glycomic_se(exp),
    glycoproteomics = glyexp::as_glycoproteomic_se(exp),
    cli::cli_abort(
      c(
        "{.arg exp} must be of type {.val glycomics} or {.val glycoproteomics}.",
        "x" = "Got {.val {exp_type}}."
      ),
      call = NULL
    )
  )
}

#' Restore the legacy container type
#'
#' @param exp A `SummarizedExperiment` object.
#' @param legacy Whether the original input used the legacy glyexp container.
#'
#' @returns `exp`, converted back to the legacy glyexp container when `legacy`
#'   is `TRUE`.
#' @noRd
.restore_data_container <- function(exp, legacy) {
  if (legacy) {
    return(glyexp::from_se(exp))
  }
  exp
}

#' Extract the abundance matrix
#'
#' @param exp A `SummarizedExperiment` object.
#'
#' @returns A numeric matrix with variables in rows and samples in columns.
#' @noRd
.get_expr_mat <- function(exp) {
  as.matrix(SummarizedExperiment::assay(exp, 1))
}

#' Extract variable information
#'
#' @param exp A `SummarizedExperiment` object.
#'
#' @returns A tibble containing `rowData` columns.
#' @noRd
.get_var_info <- function(exp) {
  row_data <- SummarizedExperiment::rowData(exp)
  tibble::as_tibble(as.list(row_data), .name_repair = "minimal")
}

#' Replace variable information
#'
#' @param exp A `SummarizedExperiment` object.
#' @param var_info A data frame containing replacement `rowData`.
#'
#' @returns `exp` with updated `rowData`.
#' @noRd
.set_var_info <- function(exp, var_info) {
  variable <- NULL
  if ("variable" %in% colnames(var_info)) {
    variable <- var_info[["variable"]]
    var_info <- dplyr::select(var_info, -dplyr::all_of("variable"))
  }
  if (is.null(variable)) {
    variable <- rownames(exp)
  }

  SummarizedExperiment::rowData(exp) <- S4Vectors::DataFrame(
    var_info,
    row.names = variable,
    check.names = FALSE
  )
  exp
}

#' Determine the glyco experiment type
#'
#' @param exp A `GlycomicSE` or `GlycoproteomicSE` object.
#'
#' @returns Either `"glycomics"` or `"glycoproteomics"`.
#' @noRd
.get_exp_type <- function(exp) {
  if (glyexp::is_glycomic_se(exp)) {
    return("glycomics")
  }
  "glycoproteomics"
}

#' Check variable-information columns
#'
#' @param exp A `SummarizedExperiment` object.
#' @param cols Required column names.
#'
#' @returns `NULL` invisibly, or an error if columns are missing.
#' @noRd
.check_var_info_cols <- function(exp, cols) {
  var_info <- .get_var_info(exp)
  has_cols <- cols %in% colnames(var_info)
  if (!all(has_cols)) {
    missing_cols <- cols[!has_cols]
    cli::cli_abort(c(
      "Variable information must contain the following columns: {.field {cols}}.",
      "x" = "Cannot find {.field {missing_cols}} in {.field var_info}."
    ))
  }
}

#' Rebuild a SummarizedExperiment with Updated Data
#'
#' @param exp A `SummarizedExperiment` object.
#' @param expr_mat A matrix of expression values.
#' @param var_info A tibble of variable information.
#' @param exp_type Experiment type for the rebuilt experiment.
#'
#' @returns A `SummarizedExperiment` object with the updated data.
#' @noRd
.rebuild_experiment <- function(exp, expr_mat, var_info, exp_type) {
  meta_data <- S4Vectors::metadata(exp)
  meta_data[["exp_type"]] <- exp_type

  variable <- var_info[["variable"]]
  rownames(expr_mat) <- variable
  var_info <- dplyr::select(var_info, -dplyr::all_of("variable"))
  assay_name <- SummarizedExperiment::assayNames(exp)[[1]]

  SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(list(expr_mat), assay_name),
    rowData = S4Vectors::DataFrame(
      var_info,
      row.names = variable,
      check.names = FALSE
    ),
    colData = SummarizedExperiment::colData(exp),
    metadata = meta_data
  )
}

#' Standardize trait variable identifiers
#'
#' @param exp A derived-trait `SummarizedExperiment` object.
#'
#' @returns `exp` with standardized row names.
#' @noRd
.standardize_trait_variables <- function(exp) {
  var_info <- .get_var_info(exp)
  exp_type <- S4Vectors::metadata(exp)[["exp_type"]]

  variables <- switch(
    exp_type,
    traitomics = var_info[["trait"]],
    traitproteomics = stringr::str_c(
      var_info[["protein"]],
      dplyr::if_else(
        is.na(var_info[["protein_site"]]),
        "X",
        as.character(var_info[["protein_site"]])
      ),
      var_info[["trait"]],
      sep = "-"
    )
  )
  rownames(exp) <- .ensure_unique_variables(variables)
  exp
}

#' Make variable identifiers unique
#'
#' @param variables A character vector of proposed identifiers.
#'
#' @returns Unique identifiers, with `-N` suffixes added to every member of a
#'   duplicated group.
#' @noRd
.ensure_unique_variables <- function(variables) {
  if (length(unique(variables)) == length(variables)) {
    return(variables)
  }

  counts <- table(variables)
  seen <- stats::setNames(rep(0L, length(counts)), names(counts))
  purrr::map_chr(variables, function(variable) {
    if (counts[[variable]] == 1L) {
      return(variable)
    }
    seen[[variable]] <<- seen[[variable]] + 1L
    stringr::str_c(variable, seen[[variable]], sep = "-")
  })
}

#' Supported AI provider names
#'
#' @returns Character vector of provider choices accepted by glydet.
#' @noRd
.ai_provider_choices <- function() {
  c(
    "deepseek",
    "openai",
    "anthropic",
    "gemini",
    "google_gemini",
    "openrouter",
    "openai_compatible"
  )
}

#' Normalize AI provider aliases
#'
#' @param provider Provider name or supported alias.
#' @returns Canonical provider name.
#' @noRd
.normalize_ai_provider <- function(
  provider = getOption("glydet.ai_provider", "deepseek")
) {
  provider <- rlang::arg_match(provider, .ai_provider_choices())
  if (identical(provider, "google_gemini")) {
    return("gemini")
  }
  provider
}

#' Normalize optional AI string arguments
#'
#' @param value Optional string value.
#' @returns `NULL` for missing or empty values, otherwise the original value.
#' @noRd
.normalize_optional_ai_string <- function(value) {
  if (is.null(value) || identical(value, "")) {
    return(NULL)
  }
  value
}

#' Human-readable AI provider label
#'
#' @param provider Provider name.
#' @returns Provider label for messages.
#' @noRd
.ai_provider_label <- function(provider) {
  switch(
    .normalize_ai_provider(provider),
    deepseek = "DeepSeek",
    openai = "OpenAI",
    anthropic = "Anthropic",
    gemini = "Google Gemini",
    openrouter = "OpenRouter",
    openai_compatible = "OpenAI-compatible"
  )
}

#' Environment variable for an AI provider API key
#'
#' @param provider Provider name.
#' @returns Environment variable name.
#' @noRd
.ai_provider_envvar <- function(provider) {
  switch(
    .normalize_ai_provider(provider),
    deepseek = "DEEPSEEK_API_KEY",
    openai = "OPENAI_API_KEY",
    anthropic = "ANTHROPIC_API_KEY",
    gemini = "GEMINI_API_KEY",
    openrouter = "OPENROUTER_API_KEY",
    openai_compatible = "OPENAI_API_KEY"
  )
}

#' Resolve the default AI model for a provider
#'
#' @param provider Provider name.
#' @param model Optional model name supplied by the caller.
#' @returns Model name, or `NULL` to use the provider default.
#' @noRd
.resolve_ai_model <- function(
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL)
) {
  provider <- .normalize_ai_provider(provider)
  model <- .normalize_optional_ai_string(model)
  if (!is.null(model)) {
    return(model)
  }
  if (identical(provider, "deepseek")) {
    return("deepseek-chat")
  }
  NULL
}

#' Resolve an AI API key
#'
#' @param provider Provider name.
#' @param api_key Optional explicit API key.
#' @returns API key string.
#' @noRd
.get_api_key <- function(
  provider = getOption("glydet.ai_provider", "deepseek"),
  api_key = getOption("glydet.ai_api_key", NULL)
) {
  api_key <- .normalize_optional_ai_string(api_key)
  if (!is.null(api_key) && nzchar(api_key)) {
    return(api_key)
  }
  provider <- .normalize_ai_provider(provider)
  envvar <- .ai_provider_envvar(provider)
  api_key <- Sys.getenv(envvar)
  if (api_key == "") {
    label <- .ai_provider_label(provider)
    bullets <- c(
      "API key for {label} chat model is not set.",
      "i" = "Please set the environment variable `{envvar}` to your API key, or pass `api_key` directly."
    )
    if (identical(provider, "openai_compatible")) {
      bullets <- c(
        bullets,
        "i" = "For OpenAI-compatible endpoints, also pass `base_url`."
      )
    }
    cli::cli_abort(bullets)
  }
  api_key
}

#' Create an ellmer chat object for a configured provider
#'
#' @param system_prompt System prompt for the chat object.
#' @param api_key API key for the selected provider.
#' @param provider Provider name.
#' @param model Optional model name.
#' @param base_url Optional provider base URL.
#' @returns An `ellmer` chat object.
#' @noRd
.create_ai_chat <- function(
  system_prompt,
  api_key,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  rlang::check_installed("ellmer")
  provider <- .normalize_ai_provider(provider)
  model <- .resolve_ai_model(provider, model)
  base_url <- .normalize_optional_ai_string(base_url)

  args <- list(
    system_prompt = system_prompt,
    model = model,
    echo = "none",
    credentials = function() api_key
  )
  if (is.null(model)) {
    args$model <- NULL
  }

  chat_fun <- switch(
    provider,
    deepseek = ellmer::chat_deepseek,
    openai = ellmer::chat_openai,
    anthropic = ellmer::chat_anthropic,
    gemini = ellmer::chat_google_gemini,
    openrouter = ellmer::chat_openrouter,
    openai_compatible = ellmer::chat_openai_compatible
  )

  if (identical(provider, "openai_compatible")) {
    if (is.null(base_url)) {
      cli::cli_abort(
        "`base_url` is required when `provider = \"openai_compatible\"`."
      )
    }
    args <- c(list(base_url = base_url, name = "OpenAI-compatible"), args)
  } else if (!is.null(base_url)) {
    args$base_url <- base_url
  }

  do.call(chat_fun, args)
}

#' Send a text-only AI request
#'
#' @param system_prompt System prompt for the AI request.
#' @param user_prompt User prompt for the AI request.
#' @param api_key API key for the selected provider.
#' @param model Optional model name.
#' @param provider Provider name.
#' @param base_url Optional provider base URL.
#' @returns Character AI response.
#' @noRd
.ask_ai <- function(
  system_prompt,
  user_prompt,
  api_key = getOption("glydet.ai_api_key", NULL),
  model = getOption("glydet.ai_model", NULL),
  provider = getOption("glydet.ai_provider", "deepseek"),
  base_url = getOption("glydet.ai_base_url", NULL)
) {
  api_key <- .get_api_key(provider = provider, api_key = api_key)
  chat <- .create_ai_chat(
    system_prompt = system_prompt,
    api_key = api_key,
    provider = provider,
    model = model,
    base_url = base_url
  )
  as.character(chat$chat(user_prompt))
}

#' Parse a Line-Based Batch AI Response
#'
#' @param output Text returned by the AI.
#' @param n Number of requested items.
#' @returns A character vector with one entry per requested item. Invalid or
#'   missing responses are `NA`.
#' @noRd
.parse_batch_response <- function(output, n) {
  values <- rep(NA_character_, n)
  lines <- stringr::str_split(paste(output, collapse = "\n"), "\n")[[1]]

  purrr::walk(lines, function(line) {
    match <- stringr::str_match(stringr::str_trim(line), "^([0-9]+)\\t(.+)$")
    if (is.na(match[1, 1])) {
      return(invisible(NULL))
    }

    position <- as.integer(match[1, 2])
    value <- stringr::str_trim(match[1, 3])
    if (
      !is.na(position) &&
        position >= 1 &&
        position <= n &&
        is.na(values[[position]]) &&
        !identical(toupper(value), "<INVALID>")
    ) {
      values[[position]] <<- value
    }
    invisible(NULL)
  })

  values
}
