# Create a Total Abundance Trait

A total abundance trait is the total abundance of a group of glycans.
For example, the total abundance of all complex glycans, or the total
abundance of all tetra-antennary glycans.

## Usage

``` r
total(cond)
```

## Arguments

- cond:

  Condition to use for defining the group of glycans. An expression that
  evaluates to a logical vector. The names of all built-in
  meta-properties (see
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md))
  and custom meta-properties can be used in the expression.

## Value

A derived trait function.

## How to use

You can use `total()` to create total abundance trait easily.

For example:

    # Total abundance of all complex glycans
    total(Tp == "complex")

    # Total abundance of all tetra-antennary glycans
    total(nA == 4)

## Examples

``` r
# Total abundance of all complex glycans
total(Tp == "complex")
#> total(Tp == "complex")

# Total abundance of all tetra-antennary glycans
total(nA == 4)
#> total(nA == 4)
```
