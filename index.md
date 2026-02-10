# glydet

The goal of glydet is to describe glycosylation structural properties in
a site specific manner. In the field of glycomics, this analytical
approach is known as derived traits. glydet provides functions to
calculate derived traits well-defined in literature, and implements a
domain-specific language to define custom derived traits.

## Installation

You can install the latest release of glydet from
[r-universe](https://glycoverse.r-universe.dev/glydet)
(**recommended**):

``` r
# install.packages("pak")
pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::pkg_install("glydet")
```

Or from [GitHub](https://github.com/glycoverse/glydet):

``` r
pak::pkg_install("glycoverse/glydet@*release")
```

Or install the development version (NOT recommended):

``` r
pak::pkg_install("glycoverse/glydet")
```

## Documentation

- 🚀 Get started:
  [Here](https://glycoverse.github.io/glydet/articles/glydet.html)
- 🔧 Custom derived traits:
  [Here](https://glycoverse.github.io/glydet/articles/custom-traits.html)
- 📚 Reference:
  [Here](https://glycoverse.github.io/glydet/reference/index.html)

## Role in `glycoverse`

glydet is a high-level package in the `glycoverse` ecosystem. It is
designed to be used by glycomics or glycoproteomics researchers directly
to calculate derived traits. It is built on top of many other packages
in the `glycoverse` ecosystem, including `glyexp`, `glyrepr`,
`glyparse`, and `glymotif`.

## Example

First, let’s load necessary packages and get the data ready.

``` r
library(glyexp)
#> Warning: 程序包'glyexp'是用R版本4.5.2 来建造的
library(glyclean)
library(glydet)

exp <- auto_clean(real_experiment)
```

``` r
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3979 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>
```

Now, let’s calculate some derived traits!

``` r
trait_exp <- derive_traits(exp)
trait_exp
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3864 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>, explanation <chr>
```

Voilà! What you see is a brand new
[`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with “traitomics” type. Think of it as your original dataset’s
sophisticated cousin 🎭 — instead of tracking “quantification of each
glycan on each glycosite in each sample,” it now contains “the value of
each derived trait on each glycosite in each sample.”

``` r
get_var_info(trait_exp)
#> # A tibble: 3,864 × 6
#>    variable protein protein_site trait gene  explanation                        
#>    <chr>    <chr>          <int> <chr> <chr> <chr>                              
#>  1 V1       A6NJW9            49 TM    CD8B2 Proportion of high-mannose glycans…
#>  2 V2       A6NJW9            49 TH    CD8B2 Proportion of hybrid glycans among…
#>  3 V3       A6NJW9            49 TC    CD8B2 Proportion of complex glycans amon…
#>  4 V4       A6NJW9            49 MM    CD8B2 Abundance-weighted mean of mannose…
#>  5 V5       A6NJW9            49 CA2   CD8B2 Proportion of bi-antennary glycans…
#>  6 V6       A6NJW9            49 CA3   CD8B2 Proportion of tri-antennary glycan…
#>  7 V7       A6NJW9            49 CA4   CD8B2 Proportion of tetra-antennary glyc…
#>  8 V8       A6NJW9            49 TF    CD8B2 Proportion of fucosylated glycans …
#>  9 V9       A6NJW9            49 TFc   CD8B2 Proportion of core-fucosylated gly…
#> 10 V10      A6NJW9            49 TFa   CD8B2 Proportion of arm-fucosylated glyc…
#> # ℹ 3,854 more rows
```

``` r
get_expr_mat(trait_exp)[1:5, 1:5]
#>    C1 C2 C3 H1 H2
#> V1  0  0  0  0  0
#> V2  0  0  0  0  0
#> V3  1  1  1  1  1
#> V4 NA NA NA NA NA
#> V5  1  1  1  1  1
```
