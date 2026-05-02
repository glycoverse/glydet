test_that("optional AI strings treat empty values as missing", {
  old_key <- Sys.getenv("OPENAI_API_KEY", unset = NA_character_)
  Sys.unsetenv("OPENAI_API_KEY")
  on.exit(
    if (is.na(old_key)) {
      Sys.unsetenv("OPENAI_API_KEY")
    } else {
      Sys.setenv(OPENAI_API_KEY = old_key)
    },
    add = TRUE
  )

  expect_error(
    make_trait(
      "proportion of sialylated glycans",
      provider = "openai",
      model = "",
      api_key = "",
      base_url = ""
    ),
    "API key for OpenAI chat model"
  )
})

test_that("empty optional AI model and base URL are omitted from chat args", {
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    chat_openai = function(system_prompt, echo, credentials, ...) {
      captured$args <- list(...)
      list(chat = function(prompt) "AI response")
    },
    .package = "ellmer"
  )

  .create_ai_chat(
    system_prompt = "system",
    api_key = "openai-key",
    provider = "openai",
    model = "",
    base_url = ""
  )

  expect_false("model" %in% names(captured$args))
  expect_false("base_url" %in% names(captured$args))
})

test_that("missing API key guidance is provider-specific", {
  old_deepseek_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = NA_character_)
  old_openai_key <- Sys.getenv("OPENAI_API_KEY", unset = NA_character_)
  Sys.unsetenv("DEEPSEEK_API_KEY")
  Sys.unsetenv("OPENAI_API_KEY")
  on.exit(
    {
      if (is.na(old_deepseek_key)) {
        Sys.unsetenv("DEEPSEEK_API_KEY")
      } else {
        Sys.setenv(DEEPSEEK_API_KEY = old_deepseek_key)
      }
      if (is.na(old_openai_key)) {
        Sys.unsetenv("OPENAI_API_KEY")
      } else {
        Sys.setenv(OPENAI_API_KEY = old_openai_key)
      }
    },
    add = TRUE
  )

  deepseek_error <- tryCatch(
    .get_api_key(provider = "deepseek"),
    error = function(e) conditionMessage(e)
  )
  compatible_error <- tryCatch(
    .get_api_key(provider = "openai_compatible"),
    error = function(e) conditionMessage(e)
  )

  expect_match(deepseek_error, "DEEPSEEK_API_KEY")
  expect_no_match(deepseek_error, "OpenAI-compatible")
  expect_match(compatible_error, "OpenAI-compatible endpoints")
})

test_that("google_gemini is accepted as a Gemini provider alias", {
  expect_equal(.normalize_ai_provider("google_gemini"), "gemini")
})
