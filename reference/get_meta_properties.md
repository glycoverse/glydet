# Get Meta-Properties of Glycans

This function calculates the meta-properties of the given glycans.
Meta-properties are properties describing certain structural
characteristics of glycans. For example, the number of antennae, the
number of core fucoses, etc.

## Usage

``` r
get_meta_properties(glycans, mp_fns = NULL)
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- mp_fns:

  A named list of meta-property functions. Names of the list are the
  names of the meta-properties. Default is
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md).
  A meta-property function should takes a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, and returns a vector of the meta-property values. purrr-style
  lambda functions are supported.

## Value

A tibble with the meta-properties. Column names are the names of the
meta-properties.

## See also

[`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md)

## Examples

``` r
glycans <- c(
  "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  "Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-",
  "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
)

# Use default meta-property functions
get_meta_properties(glycans)
#> # A tibble: 3 × 10
#>   Tp           B        nA    nF   nFc   nFa    nG   nGt    nS    nM
#>   <fct>        <lgl> <int> <int> <int> <int> <int> <int> <int> <int>
#> 1 paucimannose FALSE     0     0     0     0     0     0     0     3
#> 2 complex      FALSE     2     1     0     1     0     0     0     3
#> 3 complex      FALSE     1     0     0     0     0     0     0     3

# Use custom meta-property functions
fns <- list(
  nN = ~ glyrepr::count_mono(.x, "HexNAc"),  # purrr-style lambda function
  nH = ~ glyrepr::count_mono(.x, "Hex")
)
get_meta_properties(glycans, fns)
#> # A tibble: 3 × 2
#>      nN    nH
#>   <int> <int>
#> 1     2     3
#> 2     4     3
#> 3     3     3
```
