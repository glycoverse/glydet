# Add Meta-Properties to Experiment

This function adds meta-properties to the variable information of a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html).
Under the hood, it uses
[`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md)
to calculate the meta-properties on the "glycan_structure" column (or
column specified by `struc_col`) of the variable information tibble, and
then adds the result back as new columns.

## Usage

``` r
add_meta_properties(
  exp,
  mp_fns = NULL,
  struc_col = "glycan_structure",
  overwrite = FALSE
)
```

## Arguments

- exp:

  An
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object.

- mp_fns:

  A named list of meta-property functions. Names of the list are the
  names of the meta-properties. Default is
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md).
  A meta-property function should takes a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, and returns a vector of the meta-property values. purrr-style
  lambda functions are supported.

- struc_col:

  The column name of the glycan structures in the variable information
  tibble. Default is "glycan_structure".

- overwrite:

  Whether to overwrite the existing meta-property columns. Default is
  FALSE, raising an error if the existing columns are found.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with meta-properties added to the variable information.

## See also

[`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md),
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)

## Examples

``` r
library(glyexp)

# Compare the columns in the variable information before and after adding meta-properties
exp <- real_experiment |>
  slice_sample_var(n = 10)
colnames(get_var_info(exp))
#> [1] "variable"           "peptide"            "peptide_site"      
#> [4] "protein"            "protein_site"       "gene"              
#> [7] "glycan_composition" "glycan_structure"  

exp2 <- add_meta_properties(exp)
colnames(get_var_info(exp2))
#>  [1] "variable"           "peptide"            "peptide_site"      
#>  [4] "protein"            "protein_site"       "gene"              
#>  [7] "glycan_composition" "glycan_structure"   "Tp"                
#> [10] "B"                  "nA"                 "nF"                
#> [13] "nFc"                "nFa"                "nG"                
#> [16] "nGt"                "nS"                 "nM"                
```
