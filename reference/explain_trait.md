# Explain Derived Traits

These functions provide human-readable English explanations of derived
traits. `explain_trait()` explains one trait, while `explain_traits()`
explains a list of traits. They work with trait functions created by
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md),
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md),
[`total()`](https://glycoverse.github.io/glydet/reference/total.md), and
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md), and
work best with traits defined by built-in meta-properties. When
`use_ai = TRUE`, `explain_traits()` sends all valid traits in one
request so the shared prompt is only sent once. Traits that cannot be
explained are returned as `NA` with a warning.

## Usage

``` r
explain_trait(
  trait_fn,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
)

explain_traits(
  trait_fns,
  use_ai = FALSE,
  custom_mp = NULL,
  provider = getOption("glydet.ai_provider", "deepseek"),
  model = getOption("glydet.ai_model", NULL),
  api_key = getOption("glydet.ai_api_key", NULL),
  base_url = getOption("glydet.ai_base_url", NULL)
)
```

## Arguments

- trait_fn:

  A derived trait function created by one of the trait factories.

- use_ai:

  **\[experimental\]** Whether to use a Large Language Model (LLM) to
  explain the trait. Default is FALSE. To use this feature, you need to
  install the `ellmer` package. DeepSeek is used by default for backward
  compatibility. Other `ellmer` providers can be selected with
  `provider`, `model`, and provider-specific API key configuration.

- custom_mp:

  A named character vector of custom meta-properties. The names are the
  meta-property names, and the values are in the format "(type)
  description". Only used when `use_ai = TRUE`.

- provider:

  AI provider passed to `ellmer` when `use_ai = TRUE`. One of
  "deepseek", "openai", "anthropic", "gemini", "openrouter", or
  "openai_compatible". "google_gemini" is accepted as an alias for
  "gemini". Defaults to `getOption("glydet.ai_provider", "deepseek")`.

- model:

  Model to use when `use_ai = TRUE`. Defaults to
  `getOption("glydet.ai_model")`, or "deepseek-chat" for DeepSeek and
  the provider default for other providers.

- api_key:

  API key for the selected provider. If `NULL`, the provider specific
  environment variable is used. Defaults to
  `getOption("glydet.ai_api_key")`.

- base_url:

  Optional base URL for custom or OpenAI-compatible endpoints. Defaults
  to `getOption("glydet.ai_base_url")`.

- trait_fns:

  A list of derived trait functions created by the trait factories.

## Value

`explain_trait()` returns a character string containing a concise
English explanation. `explain_traits()` returns a character vector with
input names preserved; entries that cannot be explained are `NA`.

## Examples

``` r
# Explain built-in traits
explain_trait(traits_basic()$TM)
#> [1] "Proportion of high-mannose glycans among all glycans."
explain_trait(traits_basic()$GS)
#> [1] "Abundance-weighted mean of degree of sialylation per galactose among all glycans."

# Explain custom traits
explain_trait(prop(nFc > 0))
#> [1] "Proportion of core-fucosylated glycans among all glycans."
explain_trait(prop(nFc > 0, within = (Tp == "complex")))
#> [1] "Proportion of core-fucosylated glycans within complex glycans."
explain_trait(ratio(Tp == "complex", Tp == "hybrid"))
#> [1] "Ratio of complex glycans to hybrid glycans among all glycans."
explain_trait(wmean(nA, within = (Tp == "complex")))
#> [1] "Abundance-weighted mean of antenna count within complex glycans."
explain_trait(wmean(nS / nG, within = nA == 4 & nFc > 0))
#> [1] "Abundance-weighted mean of degree of sialylation per galactose within tetra-antennary glycans with core-fucosylation."

# Explain total and wsum traits
explain_trait(total(Tp == "complex"))
#> [1] "Total abundance of complex glycans."
explain_trait(wsum(nS))
#> [1] "Abundance-weighted sum of sialic acid count among all glycans."
explain_trait(wsum(nS, within = (Tp == "complex")))
#> [1] "Abundance-weighted sum of sialic acid count within complex glycans."

# Explain multiple traits
explain_traits(traits_basic()[c("TM", "GS")])
#>                                                                                  TM 
#>                             "Proportion of high-mannose glycans among all glycans." 
#>                                                                                  GS 
#> "Abundance-weighted mean of degree of sialylation per galactose among all glycans." 
```
