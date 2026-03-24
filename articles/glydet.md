# Get Started with glydet

Glycan derived traits are summary features calculated from individual
glycan abundances. Instead of analyzing every single glycan structure,
we combine related glycans into biologically meaningful groups. For
example, traits describing the overall level of galactosylation,
sialylation, fucosylation, or branching. These derived traits capture
broader patterns in glycosylation and reduce noise from individual
measurements.

Compared with analyzing raw glycan abundances (the “direct traits”),
using derived traits has several advantages. It simplifies the data
while keeping key biological information, making it easier to interpret
and compare across samples. Derived traits also tend to be more robust
and less affected by technical variation, and they can better highlight
biological trends or associations with phenotypes. In short, derived
traits help us see the forest rather than just the trees 🌲🌲🌲.

Enter `glydet` 🚀—your toolkit for calculating derived traits with
unprecedented precision. But here’s where it gets exciting: `glydet`
brings the power of derived traits to the glycoproteomics community for
the first time, enabling **site-specific trait analysis**. Plus, it
features an intuitive domain-specific language that lets you define
custom traits tailored to your research needs.

## Important Notes Before You Start

### Prerequisites

This package is built on the
[glyrepr](https://github.com/glycoverse/glyrepr) package, and heavily
relies on the [glyexp](https://github.com/glycoverse/glyexp) package. If
you are not familiar with these two packages, we highly recommend
checking out their introductions first. Also, to fully understand the
concepts and functions in this package, it is recommended to have a
basic understanding of the
[glymotif](https://github.com/glycoverse/glymotif) package.

### Data Types

`glydet` is designed to work with untargeted glycomics and
glycoproteomics data. Label-free quantification data is readily
supported by `glyread` and `glyclean`. For labeling quantification like
TMT, ratios between the target and reference channels (TMT ratios) must
be converted to abundance matrix following this procedure:

1.  Calculate the median of all reference channel MS2 summed intensities
    from the unnormalized intensity matrix for each glycopeptide
2.  Multiply the TMT ratios by the median MS2 summed intensities for
    each glycopeptide

This enables the quantification of different glycopeptides in the same
sample to be comparable.

## 🎯 Dive Right In: Your First Analysis

``` r
library(glydet)
library(glyrepr)
library(glyexp)
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:stats':
#> 
#>     aggregate
library(glystats)
```

Ready to see `glydet` in action? Let’s jump straight into a real-world
example that demonstrates its power! We’ll work with
[`glyexp::real_experiment`](https://glycoverse.github.io/glyexp/reference/real_experiment.html)
— an authentic N-glycoproteomics dataset from 12 patients with varying
liver conditions.

⚠️ **Pro tip:** Always preprocess your data with `glyclean` before
diving into trait analysis. This ensures your results are as clean and
reliable as your data!

``` r
exp <- auto_clean(real_experiment)  # Preprocess the data
#> 
#> ── Normalizing data ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycoproteomics". Using `normalize_median()`.
#> ✔ Normalization completed.
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ No QC samples found. Using all samples.
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 24 (0.56%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ No QC samples found. Using default imputation method based on sample size.
#> ℹ Sample size <= 30, using `impute_sample_min()`.
#> ✔ Imputation completed.
#> 
#> ── Aggregating data ──
#> 
#> ℹ Aggregating to "gfs" level
#> ✔ Aggregation completed.
#> 
#> ── Normalizing data again ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycoproteomics". Using `normalize_median()`.
#> ✔ Normalization completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column  not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3979 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>
```

Let’s take a quick peek at our dataset to understand what we’re working
with:

``` r
get_var_info(exp)
#> # A tibble: 3,979 × 6
#>    variable       protein glycan_composition glycan_structure protein_site gene 
#>    <chr>          <chr>   <comp>             <struct>                <int> <chr>
#>  1 P08185-176-He… P08185  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          176 SERP…
#>  2 P04196-344-He… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#>  3 P04196-344-He… P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          344 HRG  
#>  4 P04196-344-He… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#>  5 P10909-291-He… P10909  Hex(6)HexNAc(5)    Hex(??-?)HexNAc…          291 CLU  
#>  6 P04196-344-He… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#>  7 P04196-345-He… P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          345 HRG  
#>  8 P04196-344-He… P04196  Hex(5)HexNAc(4)dH… dHex(??-?)Hex(?…          344 HRG  
#>  9 P04196-344-He… P04196  Hex(4)HexNAc(3)    Hex(??-?)HexNAc…          344 HRG  
#> 10 P04196-344-He… P04196  Hex(4)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG  
#> # ℹ 3,969 more rows
```

``` r
get_sample_info(exp)
#> # A tibble: 12 × 2
#>    sample group
#>    <chr>  <fct>
#>  1 C1     C    
#>  2 C2     C    
#>  3 C3     C    
#>  4 H1     H    
#>  5 H2     H    
#>  6 H3     H    
#>  7 M1     M    
#>  8 M2     M    
#>  9 M3     M    
#> 10 Y1     Y    
#> 11 Y2     Y    
#> 12 Y3     Y
```

``` r
get_expr_mat(exp)[1:5, 1:5]
#>                                                C1           C2         C3
#> P08185-176-Hex(5)HexNAc(4)NeuAc(2)   6.676837e+03 2.007225e+04      13368
#> P04196-344-Hex(5)HexNAc(4)NeuAc(1)-1 3.772892e+08 5.658012e+08   99052192
#> P04196-344-Hex(5)HexNAc(4)           5.300372e+08 5.611186e+08  210626085
#> P04196-344-Hex(5)HexNAc(4)NeuAc(1)-2 3.006477e+09 2.649997e+09 1201420056
#> P10909-291-Hex(6)HexNAc(5)-1         2.772362e+07 3.181527e+07    8016730
#>                                                H1           H2
#> P08185-176-Hex(5)HexNAc(4)NeuAc(2)   4.105520e+04 1.754469e+04
#> P04196-344-Hex(5)HexNAc(4)NeuAc(1)-1 2.391413e+04 1.408332e+07
#> P04196-344-Hex(5)HexNAc(4)           9.224067e+08 8.450856e+08
#> P04196-344-Hex(5)HexNAc(4)NeuAc(1)-2 3.438030e+09 3.879662e+09
#> P10909-291-Hex(6)HexNAc(5)-1         6.820649e+07 4.501783e+07
```

Now for the magic moment ✨—let’s calculate some derived traits!

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
#>    variable      protein protein_site trait gene  explanation                   
#>    <glue>        <chr>          <int> <chr> <chr> <chr>                         
#>  1 A6NJW9-49-TM  A6NJW9            49 TM    CD8B2 Proportion of high-mannose gl…
#>  2 A6NJW9-49-TH  A6NJW9            49 TH    CD8B2 Proportion of hybrid glycans …
#>  3 A6NJW9-49-TC  A6NJW9            49 TC    CD8B2 Proportion of complex glycans…
#>  4 A6NJW9-49-MM  A6NJW9            49 MM    CD8B2 Abundance-weighted mean of ma…
#>  5 A6NJW9-49-CA2 A6NJW9            49 CA2   CD8B2 Proportion of bi-antennary gl…
#>  6 A6NJW9-49-CA3 A6NJW9            49 CA3   CD8B2 Proportion of tri-antennary g…
#>  7 A6NJW9-49-CA4 A6NJW9            49 CA4   CD8B2 Proportion of tetra-antennary…
#>  8 A6NJW9-49-TF  A6NJW9            49 TF    CD8B2 Proportion of fucosylated gly…
#>  9 A6NJW9-49-TFc A6NJW9            49 TFc   CD8B2 Proportion of core-fucosylate…
#> 10 A6NJW9-49-TFa A6NJW9            49 TFa   CD8B2 Proportion of arm-fucosylated…
#> # ℹ 3,854 more rows
```

``` r
# These are the trait values!
get_expr_mat(trait_exp)[1:5, 1:5]
#>               C1 C2 C3 H1 H2
#> A6NJW9-49-TM   0  0  0  0  0
#> A6NJW9-49-TH   0  0  0  0  0
#> A6NJW9-49-TC   1  1  1  1  1
#> A6NJW9-49-MM  NA NA NA NA NA
#> A6NJW9-49-CA2  1  1  1  1  1
```

🎉 **Congratulations!** You’ve just calculated a comprehensive suite of
derived traits in a site-specific manner:

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

💡 **The key insight:** We treat the glycans on each glycosite as a
separate mini-glycome, then calculate derived traits for each one across
all samples. For instance, if a particular glycosite hosts 10 different
glycans, the `TFc` value represents the proportion of core-fucosylated
glycans within those 10 structures in each sample.

**Important Note:** All built-in derived traits only work for N-glycans.
For other types of glycans, you need to define your own derived traits.
Check out the [Defining Custom
Traits](https://glycoverse.github.io/glydet/articles/custom-traits.html)
vignette to learn how to define your own derived traits. In fact, for
other types of glycans not so complex as N-glycans, e.g. O-glycans, we
recommend using [motif
quantification](https://glycoverse.github.io/glydet/articles/quantify-motifs.html)
instead.

Now comes the fun part! 📊 You can leverage all the powerful functions
in the `glystats` package to analyze your derived traits. Let’s
demonstrate with an ANOVA analysis to identify glycosites with
significantly different levels of core-fucosylation across conditions:

``` r
anova_res <- gly_anova(trait_exp)
#> ℹ Number of groups: 4
#> ℹ Groups: "H", "M", "Y", and "C"
#> ℹ Pairwise comparisons will be performed, with levels coming first as reference groups.
#> Warning: 281 variables failed to fit the model
anova_res$tidy_result$main_test |>
  dplyr::filter(trait == "TFc", p_adj < 0.05)
#> # A tibble: 12 × 14
#>    variable     protein protein_site trait gene  explanation term     df   sumsq
#>    <glue>       <chr>          <int> <chr> <chr> <chr>       <chr> <dbl>   <dbl>
#>  1 P00748-249-… P00748           249 TFc   F12   Proportion… group     3 5.48e-4
#>  2 P01591-71-T… P01591            71 TFc   JCHA… Proportion… group     3 7.71e-2
#>  3 P02679-78-T… P02679            78 TFc   FGG   Proportion… group     3 3.65e-3
#>  4 P02765-176-… P02765           176 TFc   AHSG  Proportion… group     3 9.41e-5
#>  5 P02790-240-… P02790           240 TFc   HPX   Proportion… group     3 6.29e-2
#>  6 P03952-494-… P03952           494 TFc   KLKB1 Proportion… group     3 2.31e-3
#>  7 P04004-86-T… P04004            86 TFc   VTN   Proportion… group     3 6.50e-3
#>  8 P04278-396-… P04278           396 TFc   SHBG  Proportion… group     3 2.99e-2
#>  9 P05090-98-T… P05090            98 TFc   APOD  Proportion… group     3 1.70e-2
#> 10 P0C0L4-1328… P0C0L4          1328 TFc   C4A   Proportion… group     3 1.74e-2
#> 11 P19652-103-… P19652           103 TFc   ORM2  Proportion… group     3 6.44e-2
#> 12 P43652-33-T… P43652            33 TFc   AFM   Proportion… group     3 5.47e-3
#> # ℹ 5 more variables: meansq <dbl>, statistic <dbl>, p_val <dbl>, p_adj <dbl>,
#> #   post_hoc <chr>
```

🔍 **Discovery time!** We’ve identified several glycosites with
statistically significant differences in core-fucosylation levels across
our patient groups — exactly the kind of biological insights that make
derived traits so powerful!

## 🔧 Under the Hood: Understanding Meta-Properties

Curious about how the magic happens? Let’s lift the hood and explore
`glydet`’s inner workings—but don’t worry, we’ll keep things accessible!

The key concept you need to understand is **“meta-properties”** - think
of them as the molecular fingerprints of individual glycans.

🆚 **What’s the difference?**

- **Derived traits** describe entire glycomes (or all glycans on a
  glycosite) and their values fluctuate across samples
- **Meta-properties** describe individual glycans regardless of their
  abundance — like counting antennae,core fucoses, or sialic acids on a
  single structure

🧠 **The connection:** Meta-properties are the building blocks for
derived traits. When you call
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md),
`glydet` automatically calculates meta-properties for all glycans first,
then uses this information to compute the derived traits you see.

Want to work with meta-properties directly? 🛠️ You’re in luck! `glydet`
provides two handy functions:

- **[`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md)**:
  Calculate meta-properties for any set of glycans
- **[`add_meta_properties()`](https://glycoverse.github.io/glydet/reference/add_meta_properties.md)**:
  Enrich your
  [`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object by adding meta-properties to variable information

### 🔬 get_meta_properties()

Let’s see
[`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md)
in action! We’ll extract a few glycan structures from our dataset:

``` r
glycans <- unique(get_var_info(exp)$glycan_structure)[1:5]
glycans
#> <glycan_structure[5]>
#> [1] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [2] NeuAc(??-?)Hex(??-?)HexNAc(??-?)[HexNAc(??-?)]Hex(??-?)[Hex(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [3] Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [4] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [5] Hex(??-?)HexNAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> # Unique structures: 5
```

📝 **Note:** `glycans` is a
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector—these are standardized representations of glycan structures.

Now watch the magic happen as we calculate their meta-properties:

``` r
get_meta_properties(glycans)
#> # A tibble: 5 × 10
#>   Tp      B        nA    nF   nFc   nFa    nG   nGt    nS    nM
#>   <fct>   <lgl> <int> <int> <int> <int> <int> <int> <int> <int>
#> 1 complex FALSE     2     0     0     0     2     0     2     3
#> 2 complex FALSE     2     0     0     0     1     0     1     4
#> 3 complex FALSE     2     0     0     0     2     2     0     3
#> 4 complex FALSE     2     0     0     0     2     1     1     3
#> 5 complex FALSE     2     0     0     0     2     1     0     4
```

### 📈 add_meta_properties()

Working with
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
objects? Perfect! You can supercharge your variable information by
adding meta-properties directly:

``` r
exp_with_mp <- add_meta_properties(exp)
get_var_info(exp_with_mp)
#> # A tibble: 3,979 × 16
#>    variable protein glycan_composition glycan_structure protein_site gene  Tp   
#>    <chr>    <chr>   <comp>             <struct>                <int> <chr> <fct>
#>  1 P08185-… P08185  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          176 SERP… comp…
#>  2 P04196-… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#>  3 P04196-… P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          344 HRG   comp…
#>  4 P04196-… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#>  5 P10909-… P10909  Hex(6)HexNAc(5)    Hex(??-?)HexNAc…          291 CLU   comp…
#>  6 P04196-… P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#>  7 P04196-… P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          345 HRG   comp…
#>  8 P04196-… P04196  Hex(5)HexNAc(4)dH… dHex(??-?)Hex(?…          344 HRG   comp…
#>  9 P04196-… P04196  Hex(4)HexNAc(3)    Hex(??-?)HexNAc…          344 HRG   hybr…
#> 10 P04196-… P04196  Hex(4)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#> # ℹ 3,969 more rows
#> # ℹ 9 more variables: B <lgl>, nA <int>, nF <int>, nFc <int>, nFa <int>,
#> #   nG <int>, nGt <int>, nS <int>, nM <int>
```

✨ **Look at that transformation!** Your variable information is now
enriched with multiple meta-property columns. This opens up powerful
filtering possibilities based on structural features.

For instance, let’s filter for all glycoforms containing high-mannose
glycans:

``` r
exp_with_mp |>
  filter_var(Tp == "highmannose")
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 332 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>, Tp <fct>, B <lgl>, nA <int>, nF <int>, nFc <int>, nFa <int>, nG <int>, nGt <int>, nS <int>, nM <int>
```

### 🧰 Meta-Property Functions: Your Structural Toolkit

Behind the scenes, meta-properties are actually functions that take
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vectors and return corresponding property values. `glydet` comes packed
with a comprehensive library of built-in meta-property functions:

``` r
names(all_mp_fns())
#>  [1] "Tp"  "B"   "nA"  "nF"  "nFc" "nFa" "nG"  "nGt" "nS"  "nM"
```

📚 **Your complete toolkit:** Here’s the full roster of built-in
meta-property functions:

| Name  | Function                                                                             | Description                                                                      |
|-------|--------------------------------------------------------------------------------------|----------------------------------------------------------------------------------|
| `Tp`  | [`n_glycan_type()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)  | Type of the glycan, either “complex”, “hybrid”, “highmannose”, or “pausimannose” |
| `B`   | [`has_bisecting()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)  | Whether the glycan has a bisecting GlcNAc                                        |
| `nA`  | [`n_antennae()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)     | Number of antennae                                                               |
| `nF`  | [`n_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | Number of fucoses                                                                |
| `nFc` | [`n_core_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)     | Number of core fucoses                                                           |
| `nFa` | [`n_arm_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)      | Number of arm fucoses                                                            |
| `nG`  | [`n_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | Number of galactoses                                                             |
| `nGt` | [`n_terminal_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | Number of terminal galactoses                                                    |
| `nS`  | [`n_sia()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | Number of sialic acids                                                           |
| `nM`  | [`n_man()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | Number of mannoses                                                               |

Each function can be called directly for quick structural analysis:

``` r
n_glycan_type(glycans)
#> [1] complex complex complex complex complex
#> Levels: paucimannose hybrid highmannose complex
```

## 🧩 Working with Structural Ambiguity

An important design principle of `glydet` is its ability to handle
glycan structures with varying levels of detail. All built-in
meta-properties and derived traits are designed to work with the
**minimum information typically available** for N-glycans in most
experimental scenarios.

### 🔧 Generic vs. Specific Monosaccharides

`glydet` works seamlessly with generic monosaccharide names (e.g.,
“Hex”, “HexNAc”, “dHex”) and structures lacking linkage information.
This level of structural resolution reflects what is commonly achievable
in glycoproteomics workflows, where complete structural determination is
often challenging.

For example, this ambiguous structure works perfectly:

``` r
# Generic monosaccharides with unknown linkages ❓
ambiguous_glycan <- "HexNAc(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-"
```

### ✨ Handling Detailed Structures

This design philosophy doesn’t limit `glydet`’s applicability to
well-characterized structures. The package equally handles glycans with
complete structural information:

``` r
# Fully specified structure with specific monosaccharides and linkages ✅
detailed_glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-3)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
```

### 🚀 Extending Functionality

When working with highly detailed structural information, you may want
to create specialized meta-property functions that leverage specific
monosaccharide identities or linkage patterns. This allows you to define
custom derived traits that capture structural features beyond the
generic framework provided by the built-in functions.

## Working with Glycomics Data

Working with glycomics data has no difference from working with
glycoproteomics data, even more straightforward as the resulting
[`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
has a simpler structure. Here we briefly demonstrate how to work with
glycomics data using
[`glyexp::real_experiment2`](https://glycoverse.github.io/glyexp/reference/real_experiment2.html).

``` r
exp <- auto_clean(real_experiment2)
#> 
#> ── Removing variables with too many missing values ──
#> 
#> ℹ No QC samples found. Using all samples.
#> ℹ Applying preset "discovery"...
#> ℹ Total removed: 10 (14.93%) variables.
#> ✔ Variable removal completed.
#> 
#> ── Normalizing data ──
#> 
#> ℹ No QC samples found. Using default normalization method based on experiment type.
#> ℹ Experiment type is "glycomics". Using `normalize_median_quotient()` + `normalize_total_area()`.
#> ✔ Normalization completed.
#> 
#> ── Normalizing data (Total Area) ──
#> 
#> ✔ Total area normalization completed.
#> 
#> ── Imputing missing values ──
#> 
#> ℹ No QC samples found. Using default imputation method based on sample size.
#> ℹ Sample size > 100, using `impute_miss_forest()`.
#> ✔ Imputation completed.
#> 
#> ── Correcting batch effects ──
#> 
#> ℹ Batch column  not found in sample_info. Skipping batch correction.
#> ✔ Batch correction completed.
trait_exp <- derive_traits(exp)
trait_exp
#> 
#> ── Traitomics Experiment ───────────────────────────────────────────────────────
#> ℹ Expression matrix: 144 samples, 14 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: trait <chr>, explanation <chr>
```

## What’s Next?

Now you have a good understanding of `glydet` and how to use it. You can
try
[`all_traits()`](https://glycoverse.github.io/glydet/reference/all_traits.md)
to calculate more advanced and detailed derived traits. You can also
start to define your own meta-property functions and derived traits.
Check out the [Custom
Traits](https://glycoverse.github.io/glydet/articles/custom-traits.html)
vignette to learn how to define your own derived traits. Or you can
check out a special type of derived traits: [motif
quantification](https://glycoverse.github.io/glydet/articles/quantify-motifs.html).
