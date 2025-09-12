
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glydet <a href="https://glycoverse.github.io/glydet/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/glydet)](https://CRAN.R-project.org/package=glydet)
[![R-CMD-check](https://github.com/glycoverse/glydet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glydet/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glydet/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glydet)
<!-- badges: end -->

The goal of glydet is to describe glycosylation structural properties in
a site specific manner. In the field of glycomics, this analytical
approach is known as derived traits. glydet provides functions to
calculate derived traits well-defined in literature, and implements a
domain-specific language to define custom derived traits.

## Installation

You can install the latest release of glydet from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("glycoverse/glydet@*release")
```

Or install the development version:

``` r
remotes::install_github("glycoverse/glydet")
```

## Documentation

-   🚀 Get started:
    [Here](https://glycoverse.github.io/glydet/articles/glydet.html)
-   🔧 Custom derived traits:
    [Here](https://glycoverse.github.io/glydet/articles/custom-traits.html)
-   📚 Reference:
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
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:stats':
#> 
#>     aggregate
library(glydet)

exp <- auto_clean(real_experiment)
#> ℹ Normalizing data (Median)
#> ✔ Normalizing data (Median) [71ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [38ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size <= 30, using sample minimum imputation
#> ℹ Imputing missing values✔ Imputing missing values [10ms]
#> 
#> ℹ Aggregating data
#> ✔ Aggregating data [349ms]
#> 
#> ℹ Normalizing data again
#> ✔ Normalizing data again [7ms]
exp
#> 
#> ── Experiment ──────────────────────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3880 variables
#> ℹ Sample information fields: group <chr>
#> ℹ Variable information fields: protein <chr>, gene <chr>, glycan_composition <glyrpr_c>, glycan_structure <glyrpr_s>, protein_site <int>
```

Now, let’s calculate some derived traits!

``` r
trait_exp <- derive_traits(exp)
trait_exp
#> 
#> ── Experiment ──────────────────────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3836 variables
#> ℹ Sample information fields: group <chr>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>
```

Voilà! What you see is a brand new `experiment()` object with
“traitomics” type. Think of it as your original dataset’s sophisticated
cousin 🎭 — instead of tracking “quantification of each glycan on each
glycosite in each sample,” it now contains “the value of each derived
trait on each glycosite in each sample.”

``` r
get_var_info(trait_exp)
#> # A tibble: 3,836 × 5
#>    variable protein protein_site trait gene 
#>    <chr>    <chr>          <int> <chr> <chr>
#>  1 V1       A6NJW9            49 TM    CD8B2
#>  2 V2       A6NJW9            49 TH    CD8B2
#>  3 V3       A6NJW9            49 TC    CD8B2
#>  4 V4       A6NJW9            49 MM    CD8B2
#>  5 V5       A6NJW9            49 CA2   CD8B2
#>  6 V6       A6NJW9            49 CA3   CD8B2
#>  7 V7       A6NJW9            49 CA4   CD8B2
#>  8 V8       A6NJW9            49 TF    CD8B2
#>  9 V9       A6NJW9            49 TFc   CD8B2
#> 10 V10      A6NJW9            49 TFa   CD8B2
#> # ℹ 3,826 more rows
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
