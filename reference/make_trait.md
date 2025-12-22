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
make_trait(description)
```

## Arguments

- description:

  A description of the trait in natural language.

## Value

A derived trait function.

## Examples

``` r
# Sys.setenv(DEEPSEEK_API_KEY = "your_api_key")
# make_trait("the average number of sialic acids")
```
