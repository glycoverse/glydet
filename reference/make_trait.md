# Use a Large Language Model (LLM) to create a derived trait function

**\[experimental\]** This function allows you to create a derived trait
function using natural language. Note that LLMs can be unreliable, so
the result should be verified manually. If the description is not clear,
an error will be raised. Try to read the descriptions of built-in traits
to get ideas. Currently, only
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md), and
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md) are
supported. To use this feature, you need to install the `ellmer`
package. You also need to provide an API key for the DeepSeek chat
model. Please set the environment variable `DEEPSEEK_API_KEY` to your
API key using [`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html).
You can obtain an API key from https://platform.deepseek.com.

## Usage

``` r
make_trait(description, custom_mp = NULL, max_retries = 2, verbose = FALSE)
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

  Maximum number of reflection retries when the AI-generated formula's
  explanation doesn't match the original description. Default is 2.

- verbose:

  Whether to print verbose output. Default is FALSE. This is useful for
  inspecting how LLMs generate trait functions.

## Value

A derived trait function.

## Examples

``` r
# Sys.setenv(DEEPSEEK_API_KEY = "your_api_key")
# my_traits <- list(
#   nS = make_trait("the average number of sialic acids"),
#   nG = make_trait("the average number of galactoses")
# )

# The trait function can then be used in `derive_traits()`:
# derive_traits(exp, trait_fns = my_traits)
```
