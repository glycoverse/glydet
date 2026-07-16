# Get Started with glydet

Derived traits are fabricated features that capture specific structural
characteristics of a glycome. For example, the proportion of core
fucosylated glycans within all glycans is a derived trait that reflects
the overall level of core fucosylation. The glycomics community has long
recognized the value of derived traits for interpreting complex
glycomics data, but calculating them has often been a manual,
error-prone process.

`glydet` makes this workflow more direct. The package provides tools for
automatically calculating a wide range of derived traits. It also
provides a flexible interface for defining custom traits. For
glycoproteomics data, `glydet` can calculate derived traits in a
site-specific manner.

## Important Notes Before You Start

### Prerequisites

This package is built on the
[glyrepr](https://github.com/glycoverse/glyrepr) package, and heavily
relies on the [glyexp](https://github.com/glycoverse/glyexp) package. If
you are not familiar with these two packages, we highly recommend
checking out their introductions first.

### Data Types

`glydet` is designed to work with untargeted glycomics and
glycoproteomics data. Label-free quantification data is readily
supported by `glyread` and `glyclean`.

For labeling quantification like TMT, ratios between the target and
reference channels (TMT ratios) must be converted to abundance matrix
following this procedure:

1.  Calculate the median of all reference channel MS2 summed intensities
    from the unnormalized intensity matrix for each glycopeptide
2.  Multiply the TMT ratios by the median MS2 summed intensities for
    each glycopeptide

This enables the quantification of different glycopeptides in the same
sample to be comparable.

## Your First Analysis

``` r

library(glydet)
library(glyrepr)
library(glyexp)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
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
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> The following object is masked from 'package:glyexp':
#> 
#>     samples
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:S4Vectors':
#> 
#>     aggregate
#> The following object is masked from 'package:stats':
#> 
#>     aggregate
library(glystats)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:Biobase':
#> 
#>     combine
#> The following objects are masked from 'package:GenomicRanges':
#> 
#>     intersect, setdiff, union
#> The following object is masked from 'package:Seqinfo':
#> 
#>     intersect
#> The following objects are masked from 'package:IRanges':
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from 'package:S4Vectors':
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from 'package:BiocGenerics':
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from 'package:generics':
#> 
#>     explain
#> The following object is masked from 'package:matrixStats':
#> 
#>     count
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

We’ll work with
[`glyexp::real_experiment2`](https://glycoverse.github.io/glyexp/reference/real_experiment2.html),
a real `GlycomicSE` dataset containing 144 samples and 67 glycans.

Before calculating derived traits, preprocess your data with `glyclean`.
This makes the trait analysis use the cleaned abundance data.

``` r

glyco_se <- auto_clean(real_experiment2)  # Preprocess the data
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 10 (14.93%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ Imputation method: `impute_miss_forest()`
#> ℹ Reason: default for "glycomics" with n_samples > 100.
#> ✔ Imputation completed.
#> 
#> ── Normalizing data ──
#> 
#> ℹ Normalization method: `normalize_total_area()`
#> ℹ Reason: default for "glycomics".
#> ✔ Normalization completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column batch not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
glyco_se
#> 
#> ── GlycomicSE ──────────────────────────────────────────────────────────────────
#> ℹ Abundance assay: 144 samples, 57 variables
#> ℹ Glycan type: N
#> ℹ Row data fields: glycan_composition <comp>, glycan_structure <struct>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>
```

Let’s inspect the dataset before calculating traits:

``` r

tibble::as_tibble(rowData(glyco_se), rownames = "variable")
#> # A tibble: 57 × 3
#>    variable                             glycan_composition      glycan_structure
#>    <chr>                                <comp>                  <struct>        
#>  1 Man(3)GlcNAc(3)                      Man(3)GlcNAc(3)         GlcNAc(?1-?)Man…
#>  2 Man(3)GlcNAc(7)                      Man(3)GlcNAc(7)         GlcNAc(?1-?)[Gl…
#>  3 Man(5)GlcNAc(2)                      Man(5)GlcNAc(2)         Man(?1-?)[Man(?…
#>  4 Man(4)Gal(2)GlcNAc(4)Neu5Ac(2)       Man(4)Gal(2)GlcNAc(4)N… Neu5Ac(?2-?)Gal…
#>  5 Man(3)Gal(1)GlcNAc(3)                Man(3)Gal(1)GlcNAc(3)   Gal(?1-?)GlcNAc…
#>  6 Man(3)Gal(2)GlcNAc(4)Fuc(2)          Man(3)Gal(2)GlcNAc(4)F… Gal(?1-?)GlcNAc…
#>  7 Man(3)GlcNAc(3)Fuc(1)                Man(3)GlcNAc(3)Fuc(1)   GlcNAc(?1-?)Man…
#>  8 Man(3)GlcNAc(4)                      Man(3)GlcNAc(4)         GlcNAc(?1-?)Man…
#>  9 Man(3)Gal(2)GlcNAc(5)Neu5Ac(1)       Man(3)Gal(2)GlcNAc(5)N… Neu5Ac(?2-?)Gal…
#> 10 Man(3)Gal(1)GlcNAc(5)Fuc(1)Neu5Ac(1) Man(3)Gal(1)GlcNAc(5)F… Neu5Ac(?2-?)Gal…
#> # ℹ 47 more rows
```

``` r

tibble::as_tibble(colData(glyco_se), rownames = "sample")
#> # A tibble: 144 × 2
#>    sample group
#>    <chr>  <fct>
#>  1 S1     H    
#>  2 S2     H    
#>  3 S3     Y    
#>  4 S4     C    
#>  5 S5     H    
#>  6 S6     C    
#>  7 S7     M    
#>  8 S8     C    
#>  9 S9     M    
#> 10 S10    M    
#> # ℹ 134 more rows
```

``` r

assay(glyco_se)[1:5, 1:5]
#>                                          S1           S2          S3
#> Man(3)GlcNAc(3)                0.0008620851 0.0010173209 0.001771775
#> Man(3)GlcNAc(7)                0.0021105912 0.0013498370 0.001590182
#> Man(5)GlcNAc(2)                0.0044181138 0.0031911698 0.002150341
#> Man(4)Gal(2)GlcNAc(4)Neu5Ac(2) 0.0028248014 0.0040683859 0.002618020
#> Man(3)Gal(1)GlcNAc(3)          0.0008569431 0.0008980034 0.001330593
#>                                          S4           S5
#> Man(3)GlcNAc(3)                0.0010535720 0.0008144561
#> Man(3)GlcNAc(7)                0.0016065736 0.0016060105
#> Man(5)GlcNAc(2)                0.0021080114 0.0024060222
#> Man(4)Gal(2)GlcNAc(4)Neu5Ac(2) 0.0024716069 0.0017383226
#> Man(3)Gal(1)GlcNAc(3)          0.0009142087 0.0006549279
```

Now let’s calculate derived traits:

``` r

trait_se <- derive_traits(glyco_se)
trait_se
#> class: SummarizedExperiment 
#> dim: 14 144 
#> metadata(2): exp_type glycan_type
#> assays(1): abundance
#> rownames(14): TM TH ... AG TS
#> rowData names(2): trait explanation
#> colnames(144): S1 S2 ... S143 S144
#> colData names(1): group
```

The result is a plain `SummarizedExperiment` with “traitomics” type.
Instead of storing the abundance of each glycan in each sample, it
stores the value of each derived trait in each sample.

``` r

tibble::as_tibble(rowData(trait_se), rownames = "variable")
#> # A tibble: 14 × 3
#>    variable trait explanation                                                   
#>    <chr>    <chr> <chr>                                                         
#>  1 TM       TM    Proportion of high-mannose glycans among all glycans.         
#>  2 TH       TH    Proportion of hybrid glycans among all glycans.               
#>  3 TC       TC    Proportion of complex glycans among all glycans.              
#>  4 MM       MM    Abundance-weighted mean of mannose count within high-mannose …
#>  5 CA2      CA2   Proportion of bi-antennary glycans within complex glycans.    
#>  6 CA3      CA3   Proportion of tri-antennary glycans within complex glycans.   
#>  7 CA4      CA4   Proportion of tetra-antennary glycans within complex glycans. 
#>  8 TF       TF    Proportion of fucosylated glycans among all glycans.          
#>  9 TFc      TFc   Proportion of core-fucosylated glycans among all glycans.     
#> 10 TFa      TFa   Proportion of arm-fucosylated glycans among all glycans.      
#> 11 TB       TB    Proportion of glycans with bisecting GlcNAc among all glycans.
#> 12 GS       GS    Abundance-weighted mean of degree of sialylation per galactos…
#> 13 AG       AG    Abundance-weighted mean of degree of galactosylation per ante…
#> 14 TS       TS    Proportion of sialylated glycans among all glycans.
```

``` r

# These are the trait values
assay(trait_se)[1:5, 1:5]
#>             S1         S2         S3         S4         S5
#> TM  0.03173244 0.02709974 0.02069125 0.01752109 0.02344732
#> TH  0.02049304 0.01848629 0.02646816 0.01943236 0.01528559
#> TC  0.94777452 0.95441397 0.95284059 0.96304655 0.96126709
#> MM  7.12257722 7.37659157 7.43377129 7.41549126 7.43498970
#> CA2 0.85120989 0.87515334 0.84044661 0.82167357 0.83514565
```

[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
calculated the following built-in derived traits:

- **`TM`**: Proportion of high-mannose glycans
- **`TH`**: Proportion of hybrid glycans  
- **`TC`**: Proportion of complex glycans
- **`MM`**: Average number of mannoses within high-mannose glycans
- **`CA2`**: Proportion of bi-antennary glycans within complex glycans
- **`CA3`**: Proportion of tri-antennary glycans within complex glycans
- **`CA4`**: Proportion of tetra-antennary glycans within complex
  glycans
- **`TF`**: Proportion of fucosylated glycans
- **`TFc`**: Proportion of core-fucosylated glycans
- **`TFa`**: Proportion of arm-fucosylated glycans
- **`TB`**: Proportion of glycans with bisecting GlcNAc
- **`GS`**: Average degree of sialylation per galactose
- **`AG`**: Average degree of galactosylation per antenna
- **`TS`**: Proportion of sialylated glycans

**Important note:** All built-in derived traits only work for N-glycans.
For other types of glycans, you need to define your own derived traits.
Check out the [Defining Custom
Traits](https://glycoverse.github.io/glydet/articles/custom-traits.html)
vignette to learn how to define your own derived traits. For glycan
classes with simpler structural patterns than N-glycans, such as
O-glycans, we recommend using [motif
quantification](https://glycoverse.github.io/glydet/articles/quantify-motifs.html)
instead.

Because
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
returns a `SummarizedExperiment`, functions in `glystats` can be applied
to the derived traits to perform statistical analyses.

``` r

anova_res <- gly_anova(trait_se)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
anova_res |>
  get_tidy_result("main_test") |>
  filter(p_adj < 0.05)
#> # A tibble: 8 × 12
#>   variable trait explanation        term     df  sumsq  meansq statistic   p_val
#>   <chr>    <chr> <chr>              <chr> <dbl>  <dbl>   <dbl>     <dbl>   <dbl>
#> 1 CA2      CA2   Proportion of bi-… group     3 0.0404 0.0135       5.71 1.02e-3
#> 2 CA3      CA3   Proportion of tri… group     3 1.75   0.584        5.38 1.56e-3
#> 3 CA4      CA4   Proportion of tet… group     3 1.12   0.372        4.82 3.19e-3
#> 4 TF       TF    Proportion of fuc… group     3 2.15   0.716        8.45 3.36e-5
#> 5 TFc      TFc   Proportion of cor… group     3 2.15   0.716        8.45 3.36e-5
#> 6 TFa      TFa   Proportion of arm… group     3 2.72   0.906        7.46 1.14e-4
#> 7 TB       TB    Proportion of gly… group     3 1.64   0.545        3.68 1.37e-2
#> 8 AG       AG    Abundance-weighte… group     3 0.0222 0.00739      3.61 1.51e-2
#> # ℹ 3 more variables: p_adj <dbl>, effect_size <dbl>, post_hoc <chr>
```

## Site-Specific Derived Traits in Glycoproteomics

Calculating derived traits for glycoproteomics data requires grouping by
glycosite.

Let’s first review the structure of glycoproteomics data. Usually, in
glycoproteomics, we analyze the quantification of glycoforms in samples.
A glycoform is a unique combination of a glycosite and a glycan
structure. Note that one glycosite can have many different glycan
structures, resulting in multiple glycoforms for the same site.

To calculate derived traits for glycoproteomics data, we calculate
derived traits separately within each glycosite. This captures
site-specific glycosylation patterns and their variations across
samples.

Operationally, there is no difference between calculating derived traits
for glycomics and glycoproteomics data.

``` r

# A GlycoproteomicSE dataset
gp_se <- auto_clean(glyexp::real_experiment)
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
gp_se
#> 
#> ── GlycoproteomicSE ────────────────────────────────────────────────────────────
#> ℹ Abundance assay: 12 samples, 3979 variables
#> ℹ Glycan type: N
#> ℹ Row data fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>, quant_method <chr>
```

``` r

gp_trait_se <- derive_traits(gp_se)
gp_trait_se
#> class: SummarizedExperiment 
#> dim: 3864 12 
#> metadata(3): exp_type glycan_type quant_method
#> assays(1): abundance
#> rownames(3864): A6NJW9-49-TM A6NJW9-49-TH ... Q9Y6W6-184-AG
#>   Q9Y6W6-184-TS
#> rowData names(5): protein protein_site trait gene explanation
#> colnames(12): C1 C2 ... Y2 Y3
#> colData names(1): group
```

The result is also a plain `SummarizedExperiment`, now with
“traitproteomics” type. The only difference is that
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
now contains a `protein_site` column, which indicates the glycosite for
each trait variable. The result therefore contains the value of each
derived trait for each glycosite in each sample.

You can again apply `glystats` functions to perform statistical
analyses. For example, we can identify glycosites with dysregulated core
fucosylation:

``` r

gly_anova(gp_trait_se) |>
  get_tidy_result("main_test") |>
  filter(trait == "TFc", p_adj < 0.05) |>
  select(protein, protein_site)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
#> Warning: 281 variables failed to fit the model
#> # A tibble: 12 × 2
#>    protein protein_site
#>    <chr>          <int>
#>  1 P00748           249
#>  2 P01591            71
#>  3 P02679            78
#>  4 P02765           176
#>  5 P02790           240
#>  6 P04004            86
#>  7 P05090            98
#>  8 P06681           621
#>  9 P0C0L4          1328
#> 10 P0C0L4          1391
#> 11 P19652           103
#> 12 P20851            64
```

## Understanding Meta-Properties

To understand how `glydet` calculates derived traits, it helps to look
at the meta-properties used by the package.

The key concept is **“meta-properties”**: properties that describe
individual glycans.

**What’s the difference?**

- **Derived traits** describe entire glycomes, or all glycans on a
  glycosite, and their values vary across samples
- **Meta-properties** describe individual glycans regardless of their
  abundance, such as the number of antennae, core fucoses, or sialic
  acids on a single structure

**The connection:** Derived traits are calculated from meta-properties
and abundance values. When you call
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md),
`glydet` automatically calculates meta-properties for all glycans first,
then uses this information to compute the derived traits you see.

You can also work with meta-properties directly through two functions:

- **[`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md)**:
  Calculate meta-properties for any set of glycans
- **[`add_meta_properties()`](https://glycoverse.github.io/glydet/reference/add_meta_properties.md)**:
  Add meta-properties to the
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  of a `GlycomicSE` or `GlycoproteomicSE` object

### get_meta_properties()

Let’s try
[`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md)
with a few glycan structures from the dataset:

``` r

glycans <- unique(rowData(glyco_se)$glycan_structure)[1:5]
glycans
#> <glycan_structure[5]>
#> [1] GlcNAc(?1-?)Man(?1-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-
#> [2] GlcNAc(?1-?)[GlcNAc(?1-?)]Man(?1-?)[GlcNAc(?1-?)[GlcNAc(?1-?)]Man(?1-?)][GlcNAc(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-
#> [3] Man(?1-?)[Man(?1-?)]Man(?1-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-
#> [4] Neu5Ac(?2-?)Gal(?1-?)GlcNAc(?1-?)[Neu5Ac(?2-?)Gal(?1-?)GlcNAc(?1-?)]Man(?1-?)[Man(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-
#> [5] Gal(?1-?)GlcNAc(?1-?)Man(?1-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-
#> # Unique structures: 5
```

**Note:** `glycans` is a
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector. These are standardized representations of glycan structures.

Now calculate their meta-properties:

``` r

get_meta_properties(glycans)
#> # A tibble: 5 × 10
#>   Tp          B        nA    nF   nFc   nFa    nG   nGt    nS    nM
#>   <fct>       <lgl> <int> <int> <int> <int> <int> <int> <int> <int>
#> 1 hybrid      FALSE     1     0     0     0     0     0     0     3
#> 2 complex     TRUE      4     0     0     0     0     0     0     3
#> 3 highmannose FALSE     0     0     0     0     0     0     0     5
#> 4 complex     FALSE     2     0     0     0     2     0     2     4
#> 5 hybrid      FALSE     1     0     0     0     1     1     0     3
```

### add_meta_properties()

When working with `GlycomicSE` or `GlycoproteomicSE` objects, you can
add meta-properties directly to
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):

``` r

se_with_mp <- add_meta_properties(glyco_se)
tibble::as_tibble(rowData(se_with_mp), rownames = "variable")
#> # A tibble: 57 × 13
#>    variable    glycan_composition glycan_structure Tp    B        nA    nF   nFc
#>    <chr>       <comp>             <struct>         <fct> <lgl> <int> <int> <int>
#>  1 Man(3)GlcN… Man(3)GlcNAc(3)    GlcNAc(?1-?)Man… hybr… FALSE     1     0     0
#>  2 Man(3)GlcN… Man(3)GlcNAc(7)    GlcNAc(?1-?)[Gl… comp… TRUE      4     0     0
#>  3 Man(5)GlcN… Man(5)GlcNAc(2)    Man(?1-?)[Man(?… high… FALSE     0     0     0
#>  4 Man(4)Gal(… Man(4)Gal(2)GlcNA… Neu5Ac(?2-?)Gal… comp… FALSE     2     0     0
#>  5 Man(3)Gal(… Man(3)Gal(1)GlcNA… Gal(?1-?)GlcNAc… hybr… FALSE     1     0     0
#>  6 Man(3)Gal(… Man(3)Gal(2)GlcNA… Gal(?1-?)GlcNAc… comp… FALSE     2     2     1
#>  7 Man(3)GlcN… Man(3)GlcNAc(3)Fu… GlcNAc(?1-?)Man… hybr… FALSE     1     1     1
#>  8 Man(3)GlcN… Man(3)GlcNAc(4)    GlcNAc(?1-?)Man… comp… FALSE     2     0     0
#>  9 Man(3)Gal(… Man(3)Gal(2)GlcNA… Neu5Ac(?2-?)Gal… comp… TRUE      2     0     0
#> 10 Man(3)Gal(… Man(3)Gal(1)GlcNA… Neu5Ac(?2-?)Gal… comp… TRUE      2     1     1
#> # ℹ 47 more rows
#> # ℹ 5 more variables: nFa <int>, nG <int>, nGt <int>, nS <int>, nM <int>
```

The variable information now contains multiple meta-property columns.
This supports filtering based on structural features.

For instance, filter for all glycoforms containing high-mannose glycans:

``` r

filter_row(se_with_mp, Tp == "highmannose")
#> 
#> ── GlycomicSE ──────────────────────────────────────────────────────────────────
#> ℹ Abundance assay: 144 samples, 5 variables
#> ℹ Glycan type: N
#> ℹ Row data fields: glycan_composition <comp>, glycan_structure <struct>, Tp <fct>, B <lgl>, nA <int>, nF <int>, nFc <int>, nFa <int>, nG <int>, nGt <int>, nS <int>, nM <int>
#> ℹ Column data fields: group <fct>
#> ℹ Metadata fields: exp_type <chr>, glycan_type <chr>
```

### Meta-Property Functions

Meta-properties are implemented as functions that take
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vectors and return corresponding property values. `glydet` includes a
set of built-in meta-property functions:

``` r

names(all_mp_fns())
#>  [1] "Tp"  "B"   "nA"  "nF"  "nFc" "nFa" "nG"  "nGt" "nS"  "nM"
```

Here is the full list of built-in meta-property functions:

| Name | Function | Description |
|----|----|----|
| `Tp` | [`n_glycan_type()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Type of the glycan, either “complex”, “hybrid”, “highmannose”, or “pausimannose” |
| `B` | [`has_bisecting()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Whether the glycan has a bisecting GlcNAc |
| `nA` | [`n_antennae()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of antennae |
| `nF` | [`n_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of fucoses |
| `nFc` | [`n_core_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of core fucoses |
| `nFa` | [`n_arm_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of arm fucoses |
| `nG` | [`n_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of galactoses |
| `nGt` | [`n_terminal_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of terminal galactoses |
| `nS` | [`n_sia()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of sialic acids |
| `nM` | [`n_man()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of mannoses |

Each function can be called directly:

``` r

n_glycan_type(glycans)
#> [1] hybrid      complex     highmannose complex     hybrid     
#> Levels: paucimannose hybrid highmannose complex
```

## Working with Structural Ambiguity

An important design principle of `glydet` is its ability to handle
glycan structures with varying levels of detail. All built-in
meta-properties and derived traits are designed to work with the
**minimum information typically available** for N-glycans in most
experimental scenarios.

### Generic vs. Specific Monosaccharides

`glydet` works with generic monosaccharide names, such as “Hex”,
“HexNAc”, and “dHex”, as well as structures lacking linkage information.
This level of structural resolution reflects what is commonly achievable
in glycoproteomics workflows, where complete structural determination is
often challenging.

For example, this ambiguous structure is supported:

``` r

# Generic monosaccharides with unknown linkages
ambiguous_glycan <- "HexNAc(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-"
```

### Handling Detailed Structures

This design philosophy doesn’t limit `glydet`’s applicability to
well-characterized structures. The package also handles glycans with
complete structural information:

``` r

# Fully specified structure with specific monosaccharides and linkages
detailed_glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-3)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
```

### Extending Functionality

When working with highly detailed structural information, you may want
to create specialized meta-property functions that leverage specific
monosaccharide identities or linkage patterns. This allows you to define
custom derived traits that use structural features beyond those covered
by the built-in functions.

## What’s Next?

Now you have the main concepts needed to use `glydet`. You can try
[`traits_detailed()`](https://glycoverse.github.io/glydet/reference/traits_detailed.md)
to calculate more detailed derived traits. You can also start to define
your own meta-property functions and derived traits. Check out the
[Custom
Traits](https://glycoverse.github.io/glydet/articles/custom-traits.html)
vignette to learn how to define your own derived traits. Or you can
check out a special type of derived traits: [motif
quantification](https://glycoverse.github.io/glydet/articles/quantify-motifs.html).
