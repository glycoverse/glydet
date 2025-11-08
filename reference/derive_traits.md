# Calculate Derived Traits

This function calculates derived traits from a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. For glycomics data, it calculates the derived traits directly.
For glycoproteomics data, each glycosite is treated as a separate
glycome, and derived traits are calculated in a site-specific manner.

## Usage

``` r
derive_traits(exp, trait_fns = NULL, mp_fns = NULL, mp_cols = NULL)
```

## Arguments

- exp:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object. Before using this function, you should preprocess the data
  using the `glyclean` package. For glycoproteomics data, the data
  should be aggregated to the "gfs" (glycoforms with structures) level
  using
  [`glyclean::aggregate()`](https://glycoverse.github.io/glyclean/reference/aggregate.html).
  Also, please make sure that the `glycan_structure` column is present
  in the `var_info` table, as not all glycoproteomics identification
  softwares provide this information. "glycan_structure" can be a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

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

- mp_cols:

  A character vector of column names in the `var_info` tibble to use as
  meta-properties. If names are provided, they will be used as names of
  the meta-properties, otherwise the column names will be used.
  Meta-properties specified in `mp_cols` will overwrite those introduced
  by `mp_fns` with the same names, including the built-in
  meta-properties. Default is `NULL`, which means no columns are used as
  meta-properties.

## Value

A new
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object for derived traits. Instead of "quantification of each glycan on
each glycosite in each sample", the new
[`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
contains "the value of each derived trait on each glycosite in each
sample", with the following columns in the `var_info` table:

- `variable`: variable ID

- `trait`: derived trait name

For glycoproteomics data, with additional columns:

- `protein`: protein ID

- `protein_site`: the glycosite position on the protein

Other columns in the `var_info` table (e.g. `gene`) are retained if they
have "many-to-one" relationship with glycosites (unique combinations of
`protein`, `protein_site`). That is, each glycosite cannot have multiple
values for these columns. `gene` is a common example, as a glycosite can
only be associate with one gene. Descriptions about glycans are not such
a column, as a glycosite can have multiple glycans, thus having multiple
descriptions. Columns not having this relationship with glycosites will
be dropped. Don't worry if you cannot understand this logic, as long as
you know that this function will try its best to preserve useful
information.

`sample_info` and `meta_data` are not modified, except that the
`exp_type` field of `meta_data` is set to "traitomics" for glycomics
data, and "traitproteomics" for glycoproteomics data.

## See also

[`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md),
[`all_traits()`](https://glycoverse.github.io/glydet/reference/all_traits.md)

## Examples

``` r
library(glyexp)
library(glyclean)
#> 
#> Attaching package: ‘glyclean’
#> The following object is masked from ‘package:stats’:
#> 
#>     aggregate

exp <- real_experiment |>
  auto_clean() |>
  slice_sample_var(n = 100)
#> ℹ Normalizing data (Median)
#> ✔ Normalizing data (Median) [191ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [19ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size <= 30, using sample minimum imputation
#> ℹ Imputing missing values
#> ✔ Imputing missing values [68ms]
#> 
#> ℹ Aggregating data
#> ✔ Aggregating data [1.6s]
#> 
#> ℹ Normalizing data again
#> ✔ Normalizing data again [16ms]
#> 
trait_exp <- derive_traits(exp)
trait_exp
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 896 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>

# By default, only basic traits are calculated
names(basic_traits())
#>  [1] "TM"  "TH"  "TC"  "MM"  "CA2" "CA3" "CA4" "TF"  "TFc" "TFa" "TB"  "GS" 
#> [13] "AG"  "TS" 

# You can calculate all traits in `all_traits()`
more_trait_exp <- derive_traits(exp, trait_fns = all_traits())
more_trait_exp
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3968 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>
```
