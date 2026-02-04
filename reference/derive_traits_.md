# Calculate Derived Traits from Tidy Data

This function calculates derived traits from a tibble in tidy format.
Use this function if you are not using the `glyexp` package. For
glycomics data, it calculates the derived traits directly. For
glycoproteomics data, each glycosite is treated as a separate glycome,
and derived traits are calculated in a site-specific manner.

## Usage

``` r
derive_traits_(tbl, data_type, trait_fns = NULL, mp_fns = NULL)
```

## Arguments

- tbl:

  A tibble in tidy format, with the following columns:

  - `sample`: sample ID

  - `glycan_structure`: glycan structures, either a
    [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
    vector or a character vector of glycan structure strings supported
    by
    [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

  - `value`: the quantification of the glycan in the sample.

  For glycoproteomics data, additional columns are needed:

  - `protein`: protein ID

  - `protein_site`: the glycosite position on the protein The unique
    combination of `protein` and `protein_site` determines a glycosite.

  Other columns are ignored.

  Please make sure that the data has been properly preprocessed,
  including normalization, missing value handling, etc. Specifically,
  for glycoproteomics data, please make sure that the data has been
  aggregated to the "glycoforms with structures" level. That is the
  quantification of each glycan structure on each glycosite in each
  sample.

- data_type:

  Either "glycomics" or "glycoproteomics".

- trait_fns:

  A named list of derived trait functions created by trait factories.
  Names of the list are the names of the derived traits. Default is
  `NULL`, which means all derived traits in
  [`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md)
  are calculated.

- mp_fns:

  A named list of meta-property functions. This parameter is useful if
  your trait functions use custom meta-properties other than those in
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md).
  Default is `NULL`, which means all meta-properties in
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md)
  are used.

## Value

A tidy tibble containing the following columns:

- `sample`: sample ID

- `trait`: derived trait name

- `value`: the value of the derived trait

- `explanation`: a concise English explanation of the trait

For glycoproteomics data, with additional columns:

- `protein`: protein ID

- `protein_site`: the glycosite position on the protein

Other columns in the original tibble are not included.

## See also

[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md),
[`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md),
[`all_traits()`](https://glycoverse.github.io/glydet/reference/all_traits.md)

## Examples

``` r
# Create example tidy data
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(glyexp)
library(tibble)

tidy_data <- as_tibble(real_experiment2)

# Calculate traits
traits <- derive_traits_(tidy_data, data_type = "glycomics")
traits
#> # A tibble: 2,016 × 4
#>    trait explanation                                           sample  value
#>    <chr> <chr>                                                 <chr>   <dbl>
#>  1 TM    Proportion of high-mannose glycans among all glycans. S1     0.0322
#>  2 TM    Proportion of high-mannose glycans among all glycans. S2     0.0274
#>  3 TM    Proportion of high-mannose glycans among all glycans. S3     0.0215
#>  4 TM    Proportion of high-mannose glycans among all glycans. S4     0.0178
#>  5 TM    Proportion of high-mannose glycans among all glycans. S5     0.0238
#>  6 TM    Proportion of high-mannose glycans among all glycans. S6     0.0254
#>  7 TM    Proportion of high-mannose glycans among all glycans. S7     0.0234
#>  8 TM    Proportion of high-mannose glycans among all glycans. S8     0.0200
#>  9 TM    Proportion of high-mannose glycans among all glycans. S9     0.0170
#> 10 TM    Proportion of high-mannose glycans among all glycans. S10    0.0207
#> # ℹ 2,006 more rows
```
