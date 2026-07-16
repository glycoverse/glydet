# Add Meta-Properties to a Glyco SummarizedExperiment

This function adds meta-properties to the variable information of a
[`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
or
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html).
Under the hood, it uses
[`get_meta_properties()`](https://glycoverse.github.io/glydet/dev/reference/get_meta_properties.md)
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

  A
  [`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
  or
  [`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
  object.

- mp_fns:

  A named list of meta-property functions. Names of the list are the
  names of the meta-properties. Default is
  [`all_mp_fns()`](https://glycoverse.github.io/glydet/dev/reference/all_mp_fns.md).
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

The input data container with meta-properties added to its variable
information. The input container type is preserved.

## See also

[`get_meta_properties()`](https://glycoverse.github.io/glydet/dev/reference/get_meta_properties.md),
[`glyexp::GlycomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html),
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)

## Examples

``` r
library(glyexp)
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
#> The following object is masked from ‘package:glyexp’:
#> 
#>     samples

# Compare rowData columns before and after adding meta-properties
gp_se <- real_experiment |>
  slice_sample_row(n = 10)
colnames(rowData(gp_se))
#> [1] "peptide"            "peptide_site"       "protein"           
#> [4] "protein_site"       "gene"               "glycan_composition"
#> [7] "glycan_structure"  

gp_se2 <- add_meta_properties(gp_se)
colnames(rowData(gp_se2))
#>  [1] "peptide"            "peptide_site"       "protein"           
#>  [4] "protein_site"       "gene"               "glycan_composition"
#>  [7] "glycan_structure"   "Tp"                 "B"                 
#> [10] "nA"                 "nF"                 "nFc"               
#> [13] "nFa"                "nG"                 "nGt"               
#> [16] "nS"                 "nM"                
```
