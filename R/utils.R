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

.check_var_info_cols <- function(exp, cols) {
  var_info <- glyexp::get_var_info(exp)
  has_cols <- cols %in% colnames(var_info)
  if (!all(has_cols)) {
    missing_cols <- cols[!has_cols]
    cli::cli_abort(c(
      "Variable information must contain the following columns: {.field {cols}}.",
      "x" = "Cannot find {.field {missing_cols}} in {.field var_info}."
    ))
  }
}

#' Rebuild an Experiment with Updated Data
#'
#' @param exp A [glyexp::experiment()] object.
#' @param expr_mat A matrix of expression values.
#' @param var_info A tibble of variable information.
#' @param exp_type Experiment type for the rebuilt experiment.
#'
#' @returns A [glyexp::experiment()] object with the updated data.
#' @noRd
.rebuild_experiment <- function(exp, expr_mat, var_info, exp_type) {
  meta_data <- glyexp::get_meta_data(exp)
  meta_data[["exp_type"]] <- exp_type
  args <- c(
    list(
      expr_mat = expr_mat,
      sample_info = glyexp::get_sample_info(exp),
      var_info = var_info
    ),
    meta_data
  )
  do.call(glyexp::experiment, args)
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
