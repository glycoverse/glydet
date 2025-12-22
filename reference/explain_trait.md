# Explain a Derived Trait

This function provides a human-readable English explanation of what a
derived trait represents. It works with trait functions created by the
trait factories
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md),
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md),
[`total()`](https://glycoverse.github.io/glydet/reference/total.md), and
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md), and
works best with traits defined by built-in meta-properties.

## Usage

``` r
explain_trait(trait_fn, use_ai = FALSE)
```

## Arguments

- trait_fn:

  A derived trait function created by one of the trait factories.

- use_ai:

  **\[experimental\]** Whether to use a Large Language Model (LLM) to
  explain the trait. Default is FALSE. To use this feature, you need to
  install the `ellmer` package. You also need to provide an API key for
  the DeepSeek chat model. Please set the environment variable
  `DEEPSEEK_API_KEY` to your API key. You can obtain an API key from
  https://platform.deepseek.com.

## Value

A character string containing a concise English explanation of the
trait.

## Examples

``` r
# Explain built-in traits
explain_trait(basic_traits()$TM)
#> [1] "Proportion of high-mannose glycans among all glycans."
explain_trait(basic_traits()$GS)
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
```
