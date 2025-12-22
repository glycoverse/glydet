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
make_trait <- function(description) {
  checkmate::assert_string(description)
  description <- stringr::str_trim(description)
  system_prompt <- .make_trait_sys_prompt(description)
  user_prompt <- paste0("INPUT: ", description, "\nOUTPUT: ")
  output <- .ask_ai(system_prompt, user_prompt)
  if (stringr::str_detect(output, "<INVALID>")) {
    cli::cli_abort(c(
      "Failed to create a derived trait function using AI.",
      "x" = "The output from AI is: {.val {output}}",
      "i" = "Please try again with a different description."
    ))
  }
  tryCatch(
    expr <- rlang::parse_expr(output),
    error = function(e) {
      cli::cli_abort(c(
        "Failed to create a derived trait function using AI.",
        "x" = "The output from AI is: {.val {output}}",
        "i" = "Error: {conditionMessage(e)}",
        "i" = "Please try again with a different description."
      ))
    }
  )
  eval(expr)
}

.make_trait_sys_prompt <- function(description) {
  paste(
    "You are a professional glycobiologist.",
    "Your task is to create a derived trait function using natural language.",
    "Derived traits are defined by expressions using meta-properties.",
    "Here are the definitions of all built-in meta-properties:",
    "- Tp: (string) glycan type (complex, hybrid, highmannose, pausimannose)",
    "- B: (logical) glycans with bisecting GlcNAc",
    "- nA: (numeric) number of antennae",
    "- nF: (numeric) number of fucoses",
    "- nFc: (numeric) number of core fucoses",
    "- nFa: (numeric) number of arm fucoses",
    "- nG: (numeric) number of galactoses",
    "- nS: (numeric) number of sialic acids",
    "- nM: (numeric) number of mannoses",
    "When encountering a meta-property that is not listed here, output an <INVALID> tag.",
    "To create a derived trait function, you need to choose one of the following factories:",
    "- prop(): for calculating the abundance proportion of a specific glycan subset",
    "- ratio(): for calculating the abundance ratio between two glycan subsets",
    "- wmean(): for calculating the weighted mean of a quantitative property, weighted by glycan abundances",
    "prop() accepts a logical expression as the first argument.",
    "ratio() accepts two logical expressions as the first two arguments.",
    "wmean() accepts a numeric expression as the first argument.",
    "Each function has a `within` parameter to restrict the calculation to a specific glycan subset, which is a logical expression.",
    "Here are some examples:",
    "INPUT: proportion of fucosylated glycans within mono-antennary glycans",
    "OUTPUT: prop(nF > 0, within = (nA == 1))",
    "INPUT: proportion of fucosylated glycans within mono-antennary glycans",
    "OUTPUT: prop(nF > 0, within = (nA == 1))",
    "INPUT: Proportion of core-fucosylated glycans within bi-antennary glycans",
    "OUTPUT: prop(nFc > 0, within = (nA == 2))",
    "INPUT: Proportion of arm-fucosylated glycans within asialylated tri-antennary glycans",
    "OUTPUT: prop(nFa > 0, within = (nA == 3 & nS == 0))",
    "INPUT: Propotion of bisecting glycans within a-core-fucosylated tetra-antennary glycans",
    "OUTPUT: prop(B, within = (nA == 4 & nFc == 0))",
    "INPUT: Ratio of complex glycans to hybrid glycans",
    'OUTPUT: ratio(Tp == "complex", Tp == "hybrid")',
    "INPUT: Ratio of bisecting glycans to unbisecting glycans within bi-antennary glycans",
    "OUTPUT: ratio(B, !B, within = (nA == 2))",
    "INPUT: Ratio of core-fucosylated glycans to non-core-fucosylated glycans within complex glycans",
    'OUTPUT: ratio(nFc > 0, nFc == 0, within = (Tp == "complex"))',
    "INPUT: Average degree of sialylation per antenna within bi-antennary glycans",
    "OUTPUT: wmean(nS / nA, within = (nA == 2))",
    "INPUT: Average degree of galactosylation per antenna within tri-antennary glycans",
    "OUTPUT: wmean(nG / nA, within = (nA == 3))",
    "INPUT: Average numbers of antennae within complex glycans",
    'OUTPUT: wmean(nA, within = (Tp == "complex"))',
    "INPUT: Average number of mannoses within highmannose glycans",
    'OUTPUT: wmean(nM, within = (Tp == "highmannose"))',
    "You need to decide:",
    "- Which factory to use",
    "- The logical expression to use as the first (and second) argument(s) of the factory",
    "- The logical expression to use as the `within` parameter of the factory",
    "- The complete derived trait function",
    "If you cannot make a trait function, output an <INVALID> tag.",
    sep = "\n"
  )
}