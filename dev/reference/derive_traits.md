# Calculate Derived Traits

This function calculates derived traits from a
[`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
or
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
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
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object. Before using this function, you should preprocess the data
  using the `glyclean` package. For glycoproteomics data, the data
  should be aggregated to the "gfs" (glycoforms with structures) level
  using
  [`glyclean::aggregate()`](https://glycoverse.github.io/glyclean/reference/aggregate.html).
  Also, please make sure that the `glycan_structure` column is present
  in
  [`SummarizedExperiment::rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
  as not all glycoproteomics identification softwares provide this
  information. "glycan_structure" can be a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector, or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).

- trait_fns:

  A named list of derived trait functions created by trait factories.
  Names of the list are the names of the derived traits. Default is
  `NULL`, which means all derived traits in
  [`traits_basic()`](https://glycoverse.github.io/glydet/dev/reference/traits_basic.md)
  are calculated.

- mp_fns:

  A named list of meta-property functions. This parameter is useful if
  your trait functions use custom meta-properties other than those in
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/dev/reference/all_mp_fns.md).
  Default is `NULL`, which means all meta-properties in
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/dev/reference/all_mp_fns.md)
  are used.

- mp_cols:

  A character vector of column names in
  [`SummarizedExperiment::rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  to use as meta-properties. If names are provided, they will be used as
  names of the meta-properties, otherwise the column names will be used.
  When `mp_cols` is specified, the selected columns overwrite
  meta-properties introduced by `mp_fns` with the same names, including
  built-in meta-properties. Default is `NULL`, which means all columns
  in
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  are available as meta-properties by their existing names. In this
  default mode, meta-properties introduced by `mp_fns` take precedence
  over
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  columns with the same names.

## Value

New `GlycomicSE` and `GlycoproteomicSE` inputs return a plain
`SummarizedExperiment`; compatible legacy glyexp inputs preserve their
legacy container type. Instead of "quantification of each glycan on each
glycosite in each sample", its assay contains "the value of each derived
trait on each glycosite in each sample", with the following columns in
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):

- `trait`: derived trait name

- `explanation`: a concise English explanation of the trait

For glycoproteomics data, with additional columns:

- `protein`: protein ID

- `protein_site`: the glycosite position on the protein

Other columns in
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
(e.g. `gene`) are retained if they have "many-to-one" relationship with
glycosites (unique combinations of `protein`, `protein_site`). That is,
each glycosite cannot have multiple values for these columns. `gene` is
a common example, as a glycosite can only be associate with one gene.
Descriptions about glycans are not such a column, as a glycosite can
have multiple glycans, thus having multiple descriptions. Columns not
having this relationship with glycosites will be dropped. Don't worry if
you cannot understand this logic, as long as you know that this function
will try its best to preserve useful information.

[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
and `metadata()` are not modified, except that the `exp_type` field of
`metadata()` is set to "traitomics" for glycomics data and
"traitproteomics" for glycoproteomics data.

## See also

[`traits_basic()`](https://glycoverse.github.io/glydet/dev/reference/traits_basic.md),
[`traits_detailed()`](https://glycoverse.github.io/glydet/dev/reference/traits_detailed.md)

## Examples

``` r
library(glyexp)
library(SummarizedExperiment)
library(glyclean)
#> 
#> Attaching package: ‘glyclean’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     aggregate
#> The following object is masked from ‘package:stats’:
#> 
#>     aggregate

gp_se <- real_experiment |>
  auto_clean() |>
  slice_sample_row(n = 100)
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 24 (0.56%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Normalizing data ──
#> 
#> ℹ Normalization method: `normalize_median()`
#> ℹ Reason: default for "glycoproteomics".
#> ✔ Normalization completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ Imputation method: `impute_min_prob()`
#> ℹ Reason: default for "glycoproteomics" with n_samples < 30.
#> ✔ Imputation completed.
#> 
#> ── Aggregating data ──
#> 
#> ℹ Aggregating to "gfs" level
#> ✔ Aggregation completed.
#> 
#> ── Normalizing data again ──
#> 
#> ℹ Normalization method: `normalize_median()`
#> ℹ Reason: default for "glycoproteomics".
#> ✔ Normalization completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column batch not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
trait_se <- derive_traits(gp_se)
rowData(trait_se)
#> DataFrame with 798 rows and 5 columns
#>                     protein protein_site       trait        gene
#>                 <character>    <integer> <character> <character>
#> O75882-1073-TM       O75882         1073          TM        ATRN
#> O75882-1073-TH       O75882         1073          TH        ATRN
#> O75882-1073-TC       O75882         1073          TC        ATRN
#> O75882-1073-MM       O75882         1073          MM        ATRN
#> O75882-1073-CA2      O75882         1073         CA2        ATRN
#> ...                     ...          ...         ...         ...
#> Q14624-517-TFa       Q14624          517         TFa       ITIH4
#> Q14624-517-TB        Q14624          517          TB       ITIH4
#> Q14624-517-GS        Q14624          517          GS       ITIH4
#> Q14624-517-AG        Q14624          517          AG       ITIH4
#> Q14624-517-TS        Q14624          517          TS       ITIH4
#>                            explanation
#>                            <character>
#> O75882-1073-TM  Proportion of high-m..
#> O75882-1073-TH  Proportion of hybrid..
#> O75882-1073-TC  Proportion of comple..
#> O75882-1073-MM  Abundance-weighted m..
#> O75882-1073-CA2 Proportion of bi-ant..
#> ...                                ...
#> Q14624-517-TFa  Proportion of arm-fu..
#> Q14624-517-TB   Proportion of glycan..
#> Q14624-517-GS   Abundance-weighted m..
#> Q14624-517-AG   Abundance-weighted m..
#> Q14624-517-TS   Proportion of sialyl..
assay(trait_se)[1:5, 1:5]
#>                 C1 C2 C3 H1 H2
#> O75882-1073-TM   0  0  0  0  0
#> O75882-1073-TH   1  1  1  1  1
#> O75882-1073-TC   0  0  0  0  0
#> O75882-1073-MM  NA NA NA NA NA
#> O75882-1073-CA2 NA NA NA NA NA
colData(trait_se)
#> DataFrame with 12 rows and 1 column
#>        group
#>     <factor>
#> C1         C
#> C2         C
#> C3         C
#> H1         H
#> H2         H
#> ...      ...
#> M2         M
#> M3         M
#> Y1         Y
#> Y2         Y
#> Y3         Y

# By default, only basic traits are calculated
names(traits_basic())
#>  [1] "TM"  "TH"  "TC"  "MM"  "CA2" "CA3" "CA4" "TF"  "TFc" "TFa" "TB"  "GS" 
#> [13] "AG"  "TS" 

# You can calculate detailed traits in `traits_detailed()`
more_trait_se <- derive_traits(gp_se, trait_fns = traits_detailed())
more_trait_se
#> class: SummarizedExperiment 
#> dim: 3534 12 
#> metadata(3): exp_type glycan_type quant_method
#> assays(1): abundance
#> rownames(3534): O75882-1073-TM O75882-1073-TH ... Q14624-517-A3GS
#>   Q14624-517-A4GS
#> rowData names(5): protein protein_site trait gene explanation
#> colnames(12): C1 C2 ... Y2 Y3
#> colData names(1): group
```
