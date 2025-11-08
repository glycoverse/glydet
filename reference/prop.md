# Create a Proportion Trait

A proportion trait is the proportion of certain group of glycans within
a larger group of glycans. For example, the proportion of sialylated
glycans within all glycans, or the proportion of tetra-antennary glycans
within all complex glycans. This type of traits is the most common type
of glycan derived traits. It can be regarded as a special case of the
ratio trait (see
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)).

## Usage

``` r
prop(cond, within = NULL, na_action = "keep")
```

## Arguments

- cond:

  Condition to use for defining the smaller group. An expression that
  evaluates to a logical vector. The names of all built-in
  meta-properties (see
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md))
  and custom meta-properties can be used in the expression.

- within:

  Condition to use for defining the larger group, with the same format
  as `cond`. If `NULL` (default), all glycans are used as the larger
  group.

- na_action:

  How to handle missing values.

  - "keep" (default): keep the missing values as NA.

  - "zero": set the missing values to 0.

## Value

A derived trait function.

## How to use

You can use `prop()` to create proportion trait easily.

For example:

    # Proportion of core-fucosylated glycans within all glycans
    prop(nFc > 0)

    # Proportion of complex glycans within all glycans
    prop(Tp == "complex")

    # Proportion of sialylated and fucosylated glycans within all glycans
    prop(nS > 0 & nFc > 0)

Note that the last example uses `&` for logical AND. Actually, you can
use any logical operator in the expression in R (e.g., `|`, `!`, etc.).

If you want to perform a pre-filtering before calculating the
proportion, for example, you want to calculate the proportion of
core-fucosylated glycans within only complex glycans, you can use
`within` to define the denominator.

    # Proportion of core-fucosylated glycans within complex glycans
    prop(nFc > 0, within = (Tp == "complex"))

    # Proportion of core-fucosylated glycans with tetra-antenary complex glycans
    prop(nFc > 0, within = (Tp == "complex" & nA == 4))

The parentheses around the condition in `within` are optional, but it is
recommended to use them for clarity.

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
# Proportion of core-fucosylated glycans within all glycans
prop(nFc > 0)
#> prop(nFc > 0, na_action = "keep")

# Proportion of bisecting glycans within all glycans
prop(B)
#> prop(B, na_action = "keep")

# Proportion of sialylated and arm-fucosylated glycans within all glycans
prop(nS > 0 & nFa > 0)
#> prop(nS > 0 & nFa > 0, na_action = "keep")

# Proportion of bi-antennary glycans within complex glycans
prop(nA == 2, within = (Tp == "complex"))
#> prop(nA == 2, within = (Tp == "complex"), na_action = "keep")

# Proportion of sialylated glycans within core-fucosylated tetra-antennary glycans
prop(nS > 0, within = (nFc > 0 & nA == 4))
#> prop(nS > 0, within = (nFc > 0 & nA == 4), na_action = "keep")
```
