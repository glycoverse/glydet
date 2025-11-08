# Explain a Derived Trait

This function provides a human-readable English explanation of what a
derived trait represents. It works with trait functions created by the
trait factories
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md), and
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md), and
works best with traits defined by built-in meta-properties.

## Usage

``` r
explain_trait(trait_fn)
```

## Arguments

- trait_fn:

  A derived trait function created by one of the trait factories.

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
explain_trait(prop(nFc > 0, within = (T == "complex")))
#> [1] "Proportion of core-fucosylated glycans within glycans satisfying 'T == \"complex\"'."
explain_trait(ratio(T == "complex", T == "hybrid"))
#> [1] "Ratio of glycans satisfying 'T == \"complex\"' to glycans satisfying 'T == \"hybrid\"' among all glycans."
explain_trait(wmean(nA, within = (T == "complex")))
#> [1] "Abundance-weighted mean of antenna count within glycans satisfying 'T == \"complex\"'."
explain_trait(wmean(nS / nG, within = nA == 4 & nFc > 0))
#> [1] "Abundance-weighted mean of degree of sialylation per galactose within tetra-antennary glycans with core-fucosylation."
```
