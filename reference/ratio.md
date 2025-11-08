# Create a Ratio Trait

A ratio trait is the ratio of total abundance of two groups of glycans.
For example, the ratio of complex glycans and hybrid glycans, or the
ratio of bisecting and unbisecting glycans.

## Usage

``` r
ratio(num_cond, denom_cond, within = NULL, na_action = "keep")
```

## Arguments

- num_cond:

  Condition to use for defining the numerator. An expression that
  evaluates to a logical vector. The names of all built-in
  meta-properties (see
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md))
  and custom meta-properties can be used in the expression.

- denom_cond:

  Condition to use for defining the denominator. Same format as
  `num_cond`.

- within:

  Condition to set a restriction for the glycans. Same format as
  `num_cond`.

- na_action:

  How to handle missing values.

  - "keep" (default): keep the missing values as NA.

  - "zero": set the missing values to 0.

## Value

A derived trait function.

## How to use

You can use `ratio()` to create ratio trait easily.

For example:

    # Ratio of complex glycans and hybrid glycans
    ratio(Tp == "complex", Tp == "hybrid")

    # Ratio of bisecting and unbisecting glycans
    ratio(B, !B)

    # Ratio of core-fucosylated and non-core-fucosylated glycans within complex glycans
    ratio(nFc > 0 & Tp == "complex", nFc == 0 & Tp == "complex")

    # The above example can be simplified as:
    ratio(nFc > 0, nFc == 0, within = (Tp == "complex"))  # more readable

Note that the last example uses `&` for logical AND. Actually, you can
use any logical operator in the expression in R (e.g., `|`, `!`, etc.).

[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md) is a
special case of `ratio()`, i.e., `prop(cond, within)` is equivalent to
`ratio(cond & within, within)`. We recommend using
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md)
instead of `ratio()` for clarity if possible.

## Note about NA

All the internal summation operations ignore NAs by default. Therefore,
NAs in the expression matrix and meta-property values will not result in
NAs in the derived traits. However, as all derived traits calculate a
ratio of two values, NAs will be introduced when:

1.  The denominator is 0. This can happen when the `within` condition
    selects no glycans.

2.  Both the numerator and denominator are 0.

## Examples

``` r
# Ratio of complex glycans and hybrid glycans
ratio(Tp == "complex", Tp == "hybrid")
#> ratio(Tp == "complex", Tp == "hybrid", na_action = "keep")

# Ratio of bisecting and unbisecting glycans within bi-antennary glycans
ratio(B, !B, within = (nA == 2))
#> ratio(B, !B, within = (nA == 2), na_action = "keep")
```
