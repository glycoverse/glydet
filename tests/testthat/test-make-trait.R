test_that("make_trait works when AI response is valid", {
  local_mocked_bindings(.ask_ai = function(system_prompt, user_prompt) "prop(nS > 0)")
  res <- make_trait("proportion of sialylated glycans")
  expect_snapshot(res)
})

test_that("make_trait raises an error when AI response is invalid", {
  local_mocked_bindings(.ask_ai = function(system_prompt, user_prompt) "<INVALID>")
  expect_error(make_trait("proportion of sialylated glycans"))

  local_mocked_bindings(.ask_ai = function(system_prompt, user_prompt) "not a valid expression")
  expect_error(make_trait("proportion of sialylated glycans"))

  local_mocked_bindings(.ask_ai = function(system_prompt, user_prompt) "prop(not good)")
  expect_error(make_trait("proportion of sialylated glycans"))
})