test_that("make_trait works when AI response is valid and consistent", {
  # Mock .get_api_key to avoid needing actual API key
  local_mocked_bindings(.get_api_key = function(...) "mock_key")

  # Mock the chat object
  mock_chat <- list(
    chat = function(prompt) "prop(nS > 0)"
  )

  # Mock ellmer::chat_deepseek to return our mock chat
  local_mocked_bindings(
    chat_deepseek = function(...) mock_chat,
    .package = "ellmer"
  )

  # Mock explain_trait with use_ai = TRUE to return consistent explanation
  local_mocked_bindings(
    explain_trait = function(trait_fn, use_ai = FALSE, custom_mp = NULL, ...) {
      if (use_ai) {
        return("Proportion of sialylated glycans among all glycans.")
      }
      "Proportion of sialylated glycans among all glycans."
    }
  )

  # Mock the consistency check to return YES
  local_mocked_bindings(
    .ask_ai = function(...) {
      "YES"
    }
  )

  res <- make_trait("proportion of sialylated glycans")
  expect_s3_class(res, "glydet_prop")
})

test_that("make_traits preserves names and marks invalid descriptions", {
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .get_api_key = function(...) "mock_key",
    .check_traits_consistency = function(...) {
      list(
        consistent = TRUE,
        explanations = "Proportion of sialylated glycans among all glycans."
      )
    },
    .create_ai_chat = function(system_prompt, ...) {
      captured$system_prompt <- system_prompt
      list(
        chat = function(user_prompt) {
          captured$user_prompt <- user_prompt
          "1\tprop(nS > 0)\n2\t<INVALID>"
        }
      )
    }
  )

  descriptions <- c(
    sialylated = "proportion of sialylated glycans",
    invalid = "the colour of glycans"
  )
  expect_warning(
    trait_fns <- make_traits(descriptions, max_retries = 0),
    "Could not make trait"
  )

  expect_named(trait_fns, names(descriptions))
  expect_s3_class(trait_fns[["sialylated"]], "glydet_prop")
  expect_true(is.na(trait_fns[["invalid"]]))
  expect_match(captured$user_prompt, "1\\tproportion of sialylated glycans")
  expect_match(captured$user_prompt, "2\\tthe colour of glycans")
})

test_that("make_traits retries inconsistent formulas with batch feedback", {
  writer_calls <- 0
  evaluator_calls <- 0
  prompts <- character()

  local_mocked_bindings(
    .get_api_key = function(...) "mock_key",
    .create_ai_chat = function(...) {
      list(chat = function(prompt) {
        writer_calls <<- writer_calls + 1
        prompts[[writer_calls]] <<- prompt
        if (writer_calls == 1) {
          "1\tprop(nS > 0)\n2\tprop(nS > 0)"
        } else {
          "2\tprop(nF > 0)"
        }
      })
    },
    .check_traits_consistency = function(...) {
      evaluator_calls <<- evaluator_calls + 1
      if (evaluator_calls == 1) {
        list(
          consistent = c(TRUE, FALSE),
          explanations = c(
            "Proportion of sialylated glycans among all glycans.",
            "Proportion of sialylated glycans among all glycans."
          )
        )
      } else {
        list(
          consistent = TRUE,
          explanations = "Proportion of fucosylated glycans among all glycans."
        )
      }
    }
  )

  trait_fns <- make_traits(
    c(
      "proportion of sialylated glycans",
      "proportion of fucosylated glycans"
    ),
    max_retries = 2
  )

  expect_s3_class(trait_fns[[1]], "glydet_prop")
  expect_s3_class(trait_fns[[2]], "glydet_prop")
  expect_equal(writer_calls, 2)
  expect_equal(evaluator_calls, 2)
  expect_no_match(prompts[[2]], "1\\tOriginal description")
  expect_match(prompts[[2]], "2\\tOriginal description")
  expect_match(prompts[[2]], "prop\\(nS > 0\\)")
  expect_match(prompts[[2]], "Proportion of sialylated glycans")
})

test_that("make_traits stops retrying after max_retries", {
  writer_calls <- 0
  local_mocked_bindings(
    .get_api_key = function(...) "mock_key",
    .create_ai_chat = function(...) {
      list(chat = function(...) {
        writer_calls <<- writer_calls + 1
        "1\tprop(nF > 0)"
      })
    },
    .check_traits_consistency = function(...) {
      list(
        consistent = FALSE,
        explanations = "Proportion of fucosylated glycans among all glycans."
      )
    }
  )

  expect_warning(
    trait_fns <- make_traits(
      "proportion of sialylated glycans",
      max_retries = 1
    ),
    "Could not make trait"
  )

  expect_true(is.na(trait_fns[[1]]))
  expect_equal(writer_calls, 2)
})

test_that("batch consistency uses an AI explainer and evaluator", {
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    explain_traits = function(trait_fns, use_ai = FALSE, ...) {
      captured$use_ai <- use_ai
      c(
        "Proportion of sialylated glycans among all glycans.",
        "Proportion of fucosylated glycans among all glycans."
      )
    },
    .ask_ai = function(system_prompt, user_prompt, ...) {
      captured$system_prompt <- system_prompt
      captured$user_prompt <- user_prompt
      "1\tYES\n2\tNO"
    }
  )

  result <- .check_traits_consistency(
    c(
      "proportion of sialylated glycans",
      "proportion of sialylated glycans"
    ),
    list(prop(nS > 0), prop(nF > 0))
  )

  expect_true(captured$use_ai)
  expect_equal(result$consistent, c(TRUE, FALSE))
  expect_length(result$explanations, 2)
  expect_match(captured$user_prompt, "Generated explanation")
})

test_that("make_trait routes explicit provider settings to ellmer", {
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    chat_openai = function(system_prompt, model, echo, credentials) {
      captured$system_prompt <- system_prompt
      captured$model <- model
      captured$echo <- echo
      captured$api_key <- credentials()
      list(chat = function(prompt) "prop(nS > 0)")
    },
    .package = "ellmer"
  )

  local_mocked_bindings(
    explain_trait = function(
      trait_fn,
      use_ai = FALSE,
      custom_mp = NULL,
      api_key = NULL,
      model = NULL,
      provider = NULL,
      base_url = NULL
    ) {
      if (use_ai) {
        captured$explain_api_key <- api_key
        captured$explain_model <- model
        captured$explain_provider <- provider
        captured$explain_base_url <- base_url
      }
      "Proportion of sialylated glycans among all glycans."
    }
  )

  local_mocked_bindings(
    .ask_ai = function(
      system_prompt,
      user_prompt,
      api_key = NULL,
      model = NULL,
      provider = NULL,
      base_url = NULL
    ) {
      captured$consistency_api_key <- api_key
      captured$consistency_model <- model
      captured$consistency_provider <- provider
      captured$consistency_base_url <- base_url
      "YES"
    }
  )

  res <- make_trait(
    "proportion of sialylated glycans",
    provider = "openai",
    model = "gpt-4o-mini",
    api_key = "openai-key"
  )

  expect_s3_class(res, "glydet_prop")
  expect_match(captured$system_prompt, "professional glycobiologist")
  expect_equal(captured$model, "gpt-4o-mini")
  expect_equal(captured$echo, "none")
  expect_equal(captured$api_key, "openai-key")
  expect_equal(captured$explain_provider, "openai")
  expect_equal(captured$explain_model, "gpt-4o-mini")
  expect_equal(captured$explain_api_key, "openai-key")
  expect_null(captured$explain_base_url)
  expect_equal(captured$consistency_provider, "openai")
  expect_equal(captured$consistency_model, "gpt-4o-mini")
  expect_equal(captured$consistency_api_key, "openai-key")
  expect_null(captured$consistency_base_url)
})

test_that("make_trait uses package-level AI provider options", {
  captured <- new.env(parent = emptyenv())
  old_options <- options(
    glydet.ai_provider = "openai",
    glydet.ai_model = "gpt-4o-mini",
    glydet.ai_api_key = "option-key"
  )
  on.exit(options(old_options), add = TRUE)

  local_mocked_bindings(
    chat_openai = function(system_prompt, model, echo, credentials) {
      captured$model <- model
      captured$api_key <- credentials()
      list(chat = function(prompt) "prop(nS > 0)")
    },
    .package = "ellmer"
  )

  local_mocked_bindings(
    explain_trait = function(trait_fn, use_ai = FALSE, custom_mp = NULL, ...) {
      "Proportion of sialylated glycans among all glycans."
    },
    .ask_ai = function(...) "YES"
  )

  res <- make_trait("proportion of sialylated glycans")

  expect_s3_class(res, "glydet_prop")
  expect_equal(captured$model, "gpt-4o-mini")
  expect_equal(captured$api_key, "option-key")
})

test_that("make_trait supports OpenAI-compatible endpoints", {
  captured <- new.env(parent = emptyenv())

  local_mocked_bindings(
    chat_openai_compatible = function(
      base_url,
      name,
      system_prompt,
      model = NULL,
      echo,
      credentials
    ) {
      captured$base_url <- base_url
      captured$name <- name
      captured$model <- model
      captured$api_key <- credentials()
      list(chat = function(prompt) "prop(nS > 0)")
    },
    .package = "ellmer"
  )

  local_mocked_bindings(
    explain_trait = function(trait_fn, use_ai = FALSE, custom_mp = NULL, ...) {
      "Proportion of sialylated glycans among all glycans."
    },
    .ask_ai = function(...) "YES"
  )

  res <- make_trait(
    "proportion of sialylated glycans",
    provider = "openai_compatible",
    api_key = "compatible-key",
    base_url = "https://example.test/v1"
  )

  expect_s3_class(res, "glydet_prop")
  expect_equal(captured$base_url, "https://example.test/v1")
  expect_equal(captured$name, "OpenAI-compatible")
  expect_null(captured$model)
  expect_equal(captured$api_key, "compatible-key")
})

test_that("make_trait retries when explanation is inconsistent", {
  # Mock .get_api_key
  local_mocked_bindings(.get_api_key = function(...) "mock_key")

  call_count <- 0

  # Mock the chat object that returns correct formula on second attempt
  mock_chat <- list(
    chat = function(prompt) {
      call_count <<- call_count + 1
      if (call_count == 1) {
        # First attempt: wrong formula (but syntactically valid)
        "prop(nF > 0)"
      } else {
        # Second attempt: correct formula
        "prop(nS > 0)"
      }
    }
  )

  local_mocked_bindings(
    chat_deepseek = function(...) mock_chat,
    .package = "ellmer"
  )

  # Mock explain_trait
  local_mocked_bindings(
    explain_trait = function(trait_fn, use_ai = FALSE, custom_mp = NULL, ...) {
      if (use_ai) {
        attr_val <- attr(trait_fn, "cond")
        if (!is.null(attr_val) && rlang::expr_text(attr_val) == "nF > 0") {
          return("Proportion of fucosylated glycans among all glycans.")
        }
        return("Proportion of sialylated glycans among all glycans.")
      }
      "some explanation"
    }
  )

  # Mock consistency check: NO for first, YES for second
  consistency_call_count <- 0
  local_mocked_bindings(
    .ask_ai = function(...) {
      consistency_call_count <<- consistency_call_count + 1
      if (consistency_call_count == 1) "NO" else "YES"
    }
  )

  res <- make_trait("proportion of sialylated glycans")
  expect_s3_class(res, "glydet_prop")
  expect_equal(call_count, 2)
})

test_that("make_trait raises an error when AI response is invalid", {
  # Mock .get_api_key
  local_mocked_bindings(.get_api_key = function(...) "mock_key")

  # Mock the chat object
  mock_chat <- list(
    chat = function(prompt) "<INVALID>"
  )

  local_mocked_bindings(
    chat_deepseek = function(...) mock_chat,
    .package = "ellmer"
  )

  expect_error(make_trait("proportion of sialylated glycans", max_retries = 0))
})

test_that("make_trait raises an error after max retries with inconsistent explanations", {
  # Mock .get_api_key
  local_mocked_bindings(.get_api_key = function(...) "mock_key")

  # Mock the chat object
  mock_chat <- list(
    chat = function(prompt) "prop(nF > 0)"
  )

  local_mocked_bindings(
    chat_deepseek = function(...) mock_chat,
    .package = "ellmer"
  )

  # Mock explain_trait
  local_mocked_bindings(
    explain_trait = function(trait_fn, use_ai = FALSE, custom_mp = NULL, ...) {
      if (use_ai) {
        return("Proportion of fucosylated glycans among all glycans.")
      }
      "some explanation"
    }
  )

  # Always return NO for consistency
  local_mocked_bindings(
    .ask_ai = function(...) "NO"
  )

  expect_error(
    make_trait("proportion of sialylated glycans", max_retries = 2),
    "consistent"
  )
})

test_that("custom_mp parameter includes custom meta-properties in system prompt", {
  custom_mp <- c(
    nE = "(integer) number of a2,6-linked sialic acids",
    nL = "(integer) number of a2,3-linked sialic acids"
  )

  prompt <- glydet:::.make_trait_sys_prompt("test", custom_mp)

  expect_match(prompt, "user-defined custom meta-properties")
  expect_match(prompt, "nE: \\(integer\\) number of a2,6-linked sialic acids")
  expect_match(prompt, "nL: \\(integer\\) number of a2,3-linked sialic acids")
})

test_that("trait prompt omits within for all-glycan traits", {
  prompt <- glydet:::.make_trait_sys_prompt()

  expect_match(
    prompt,
    "Only include `within` when the description names a subset; otherwise omit it."
  )
  expect_match(prompt, "Never use `within = TRUE`.")
})

test_that("custom_mp works with make_trait", {
  # Mock .get_api_key
  local_mocked_bindings(.get_api_key = function(...) "mock_key")

  captured_system_prompt <- NULL

  # Mock ellmer::chat_deepseek to capture the system prompt
  local_mocked_bindings(
    chat_deepseek = function(system_prompt, ...) {
      captured_system_prompt <<- system_prompt
      list(chat = function(prompt) "prop(nE > 0)")
    },
    .package = "ellmer"
  )

  # Mock explain_trait
  local_mocked_bindings(
    explain_trait = function(trait_fn, use_ai = FALSE, custom_mp = NULL, ...) {
      if (use_ai) {
        return("Proportion of glycans with a2,6-linked sialic acids.")
      }
      "some explanation"
    }
  )

  # Mock the consistency check to return YES
  local_mocked_bindings(
    .ask_ai = function(...) "YES"
  )

  custom_mp <- c(nE = "(integer) number of a2,6-linked sialic acids")

  res <- make_trait(
    "proportion of glycans with a2,6-linked sialic acids",
    custom_mp = custom_mp
  )

  # Verify system prompt contains custom meta-property
  expect_match(
    captured_system_prompt,
    "nE: \\(integer\\) number of a2,6-linked sialic acids"
  )
  expect_s3_class(res, "glydet_prop")
})
