# Create a Weighted-Mean Trait

A weighted-mean trait is the average value of some quantitative property
within a group of glycans, weighted by the abundance of the glycans. For
example, the average number of antennae within all complex glycans, or
the average number of sialic acids within all glycans.

## Usage

``` r
wmean(val, within = NULL, na_action = "keep")
```

## Arguments

- val:

  Expression to use for defining the value. An expression that evaluates
  to a numeric vector. The names of all built-in meta-properties (see
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md))
  and custom meta-properties can be used in the expression.

- within:

  Condition to set a restriction for the glycans. Same format as `val`.

- na_action:

  How to handle missing values.

  - "keep" (default): keep the missing values as NA.

  - "zero": set the missing values to 0.

## Value

A derived trait function.

## How to use

You can use `wmean()` to create weighted-mean trait easily.

For example:

    # Weighted mean of the number of sialic acids within all glycans
    wmean(nS)

    # Average degree of sialylation per antenna within all glycans
    wmean(nS / nA)

Note that the last example uses `/` for division. Actually, you can use
any arithmetic operator in the expression in R (e.g., `*`, `+`, `-`,
etc.).

If you want to perform a pre-filtering before calculating the
weighted-mean, for example, you want to calculate the average degree of
sialylation per antenna within only complex glycans, you can use
`within` to define the restriction.

    # Average number of antennae within complex glycans
    wmean(nA, within = (Tp == "complex"))

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
# Weighted mean of the number of sialic acids within all glycans
wmean(nS)
#> wmean(nS, na_action = "keep")

# Average degree of sialylation per antenna within all glycans
wmean(nS / nA)
#> wmean(nS/nA, na_action = "keep")

# Average number of antennae within complex glycans
wmean(nA, within = (Tp == "complex"))
#> wmean(nA, within = (Tp == "complex"), na_action = "keep")
```
