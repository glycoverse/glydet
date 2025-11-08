# Create a Weighted Sum Trait

A weighted sum trait is the sum of a quantitative property within a
group of glycans, weighted by the abundance of the glycans. For example,
the sum of the number of sialic acids within all glycans, or the sum of
the number of Lewis x antigens within all glycans.

## Usage

``` r
wsum(val, within = NULL)
```

## Arguments

- val:

  Expression to use for defining the value. An expression that evaluates
  to a logical vector. The names of all built-in meta-properties (see
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md))
  and custom meta-properties can be used in the expression.

- within:

  Condition to set a restriction for the glycans. Same format as `val`.

## Value

A derived trait function.

## How to use

You can use `wsum()` to create weighted sum trait easily.

For example:

    # Weighted sum of the number of sialic acids within all glycans
    wsum(nS)

This can be regarded as the quantification of sialic acids. If some
glycan has only one sialic acid, its abundance is added to the results.
If another glycan has two sialic acids, its abundance is doubled before
being added to the results.

You can also use `within` to restrict the weighted sum calculation to
specific glycan subsets. For example, you can calculate the weighted sum
of the number of sialic acids within complex glycans:

    wsum(nS, within = (Tp == "complex"))

## Examples

``` r
# Weighted sum of the number of sialic acids within all glycans
wsum(nS)
#> wsum(nS)

# Weighted sum of the number of sialic acids within complex glycans
wsum(nS, within = (Tp == "complex"))
#> wsum(nS, within = (Tp == "complex"))
```
