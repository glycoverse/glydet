test_that("make_trait works when AI response is valid and consistent", {
  # Mock .get_api_key to avoid needing actual API key
  local_mocked_bindings(.get_api_key = function() "mock_key")

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
    explain_trait = function(trait_fn, use_ai = FALSE) {
      if (use_ai) {
        return("Proportion of sialylated glycans among all glycans.")
      }
      "Proportion of sialylated glycans among all glycans."
    }
  )

  # Mock the consistency check to return YES
  local_mocked_bindings(
    .ask_ai = function(system_prompt, user_prompt, model = "deepseek-chat") {
      "YES"
    }
  )

  res <- make_trait("proportion of sialylated glycans")
  expect_s3_class(res, "glydet_prop")
})

test_that("make_trait retries when explanation is inconsistent", {
  # Mock .get_api_key
  local_mocked_bindings(.get_api_key = function() "mock_key")

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
    explain_trait = function(trait_fn, use_ai = FALSE) {
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
    .ask_ai = function(system_prompt, user_prompt, model = "deepseek-chat") {
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
  local_mocked_bindings(.get_api_key = function() "mock_key")

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
  local_mocked_bindings(.get_api_key = function() "mock_key")

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
    explain_trait = function(trait_fn, use_ai = FALSE) {
      if (use_ai) {
        return("Proportion of fucosylated glycans among all glycans.")
      }
      "some explanation"
    }
  )

  # Always return NO for consistency
  local_mocked_bindings(
    .ask_ai = function(system_prompt, user_prompt, model = "deepseek-chat") "NO"
  )

  expect_error(
    make_trait("proportion of sialylated glycans", max_retries = 2),
    "consistent"
  )
})