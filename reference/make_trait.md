# Use a Large Language Model (LLM) to Create Derived Trait Functions

**\[experimental\]** These functions create derived trait functions from
natural-language descriptions. `make_trait()` creates one trait and
checks its consistency with the requested description. `make_traits()`
creates and validates a list of traits in batched requests, reducing
repeated prompt tokens. For batch creation, descriptions that cannot be
understood, produce an invalid formula, or do not match the generated
trait are returned as `NA` with a warning. LLM-generated traits should
always be verified manually. Try to read the descriptions of built-in
traits to get ideas. Currently, only
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md), and
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md) are
supported. To use this feature, you need to install the `ellmer`
package. DeepSeek is used by default for backward compatibility. Other
`ellmer` providers can be selected with `provider`, `model`, and
provider-specific API key configuration.

## Usage

``` r
make_trait(
  description,
  custom_mp = NULL,
  max_retries = 2,
  verbose = FALSE,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
)

make_traits(
  descriptions,
  custom_mp = NULL,
  max_retries = 2,
  verbose = FALSE,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
)
```

## Arguments

- description:

  A description of the trait in natural language.

- custom_mp:

  A named character vector of custom meta-properties. The names are the
  meta-property names, and the values are in the format "(type)
  description". For example:
  `c(nE = "(integer) number of a2,6-linked sialic acids")`. These custom
  meta-properties will be available for the LLM to use. Note that
  defining the meta-properties here is not enough for you to use them.
  You need to define corresponding meta-property functions or specifying
  meta-property columns. For more information about custom
  meta-properties, see the vignette [Custom
  Meta-Properties](https://glycoverse.github.io/glydet/articles/custom-traits.html#using-make_trait).

- max_retries:

  Maximum number of retries after an invalid formula or an explanation
  that doesn't match the original description. In batch mode, only
  unresolved descriptions are retried. Default is 2.

- verbose:

  Whether to print verbose output. Default is FALSE. This is useful for
  inspecting how LLMs generate trait functions.

- provider:

  AI provider passed to `ellmer`. One of "deepseek", "openai",
  "anthropic", "gemini", "openrouter", or "openai_compatible".
  "google_gemini" is accepted as an alias for "gemini". Defaults to
  `getOption("glydet.ai_provider", "deepseek")`.

- model:

  Model to use. Defaults to `getOption("glydet.ai_model")`, or
  "deepseek-chat" for DeepSeek and the provider default for other
  providers.

- api_key:

  API key for the selected provider. If `NULL`, the provider specific
  environment variable is used. Defaults to
  `getOption("glydet.ai_api_key")`.

- base_url:

  Optional base URL for custom or OpenAI-compatible endpoints. Defaults
  to `getOption("glydet.ai_base_url")`.

- descriptions:

  A character vector of trait descriptions.

## Value

`make_trait()` returns a derived trait function. `make_traits()` returns
a list of derived trait functions with input names preserved; entries
that cannot be created are `NA`.

## Batch multi-agent workflow

`make_traits()` uses one batch writer to generate formulas for all
active descriptions, one batch explainer to describe the generated
formulas, and one batch evaluator to compare those explanations with the
original descriptions. Successful positions are retained. Invalid or
mismatched positions are sent back to the writer with their validation
error or generated explanation, and only those positions are
regenerated. This continues until all positions pass or `max_retries` is
reached.

## Examples

``` r
# Sys.setenv(DEEPSEEK_API_KEY = "your_api_key")
# my_traits <- list(
#   nS = make_trait("the average number of sialic acids"),
#   nG = make_trait("the average number of galactoses")
# )

# The trait function can then be used in `derive_traits()`:
# derive_traits(exp, trait_fns = my_traits)

if (FALSE) { # \dontrun{
make_traits(c(
  sialylated = "proportion of sialylated glycans",
  galactose = "average number of galactoses"
))
} # }
```
