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
#> ℹ Normalizing data (Median)
#> ✔ Normalizing data (Median) [144ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [82ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size <= 30, using sample minimum imputation
#> ℹ Imputing missing values✔ Imputing missing values [27ms]
#> 
#> ℹ Aggregating data
#> ✔ Aggregating data [1.1s]
#> 
#> ℹ Normalizing data again
#> ✔ Normalizing data again [21ms]
exp
#> 
#> ── Glycoproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3880 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, glycan_composition <comp>, glycan_structure <struct>, protein_site <int>, gene <chr>
```

Let’s take a quick peek at our dataset to understand what we’re working
with:

``` r
get_var_info(exp)
#> # A tibble: 3,880 × 6
#>    variable protein glycan_composition      glycan_structure  protein_site gene 
#>    <chr>    <chr>   <comp>                  <struct>                 <int> <chr>
#>  1 V1       P08185  Hex(5)HexNAc(4)NeuAc(2) NeuAc(??-?)Hex(?…          176 SERP…
#>  2 V2       P04196  Hex(5)HexNAc(4)NeuAc(1) NeuAc(??-?)Hex(?…          344 HRG  
#>  3 V3       P04196  Hex(5)HexNAc(4)         Hex(??-?)HexNAc(…          344 HRG  
#>  4 V4       P04196  Hex(5)HexNAc(4)NeuAc(1) NeuAc(??-?)Hex(?…          344 HRG  
#>  5 V5       P10909  Hex(6)HexNAc(5)         Hex(??-?)HexNAc(…          291 CLU  
#>  6 V6       P04196  Hex(5)HexNAc(4)NeuAc(2) NeuAc(??-?)Hex(?…          344 HRG  
#>  7 V7       P04196  Hex(5)HexNAc(4)         Hex(??-?)HexNAc(…          345 HRG  
#>  8 V8       P04196  Hex(5)HexNAc(4)dHex(2)  dHex(??-?)Hex(??…          344 HRG  
#>  9 V9       P04196  Hex(4)HexNAc(3)         Hex(??-?)HexNAc(…          344 HRG  
#> 10 V10      P04196  Hex(4)HexNAc(4)NeuAc(1) NeuAc(??-?)Hex(?…          344 HRG  
#> # ℹ 3,870 more rows
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
#>              C1           C2           C3           H1           H2
#> V1 6.626760e+03 2.019159e+04      13432.7 4.072473e+04 1.771879e+04
#> V2 3.744595e+08 5.691652e+08   99531624.5 2.372164e+04 1.422307e+07
#> V3 5.260619e+08 5.644547e+08  211645556.7 9.149818e+08 8.534716e+08
#> V4 2.983928e+09 2.665752e+09 1207235166.5 3.410355e+09 3.918161e+09
#> V5 2.751569e+07 3.200443e+07    8055532.6 6.765746e+07 4.546455e+07
```

Now for the magic moment ✨—let’s calculate some derived traits!

``` r
trait_exp <- derive_traits(exp)
trait_exp
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 3836 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>
```

Voilà! What you see is a brand new
[`experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with “traitomics” type. Think of it as your original dataset’s
sophisticated cousin 🎭 — instead of tracking “quantification of each
glycan on each glycosite in each sample,” it now contains “the value of
each derived trait on each glycosite in each sample.”

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
# These are the trait values!
get_expr_mat(trait_exp)[1:5, 1:5]
#>    C1 C2 C3 H1 H2
#> V1  0  0  0  0  0
#> V2  0  0  0  0  0
#> V3  1  1  1  1  1
#> V4 NA NA NA NA NA
#> V5  1  1  1  1  1
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
#> Warning: 267 variables failed to fit the model
anova_res$tidy_result$main_test |>
  dplyr::filter(trait == "TFc", p_adj < 0.05)
#> # A tibble: 12 × 13
#>    variable term     df     sumsq    meansq statistic     p_val   p_adj post_hoc
#>    <chr>    <chr> <dbl>     <dbl>     <dbl>     <dbl>     <dbl>   <dbl> <chr>   
#>  1 V1115    group     3 0.0000941 0.0000314      26.2   1.72e-4 7.80e-3 H_vs_M;…
#>  2 V1227    group     3 0.0629    0.0210         14.0   1.50e-3 3.06e-2 H_vs_Y;…
#>  3 V1353    group     3 0.00231   0.000770       19.3   5.05e-4 1.61e-2 H_vs_C;…
#>  4 V1381    group     3 0.00640   0.00213        14.9   1.23e-3 2.69e-2 H_vs_Y;…
#>  5 V1661    group     3 0.0299    0.00998        14.3   1.40e-3 2.92e-2 H_vs_M;…
#>  6 V1675    group     3 0.0174    0.00581        43.1   2.78e-5 2.90e-3 H_vs_M;…
#>  7 V2165    group     3 0.0174    0.00581       172.    1.34e-7 1.01e-4 H_vs_M;…
#>  8 V2487    group     3 0.0644    0.0215         74.1   3.53e-6 8.91e-4 H_vs_M;…
#>  9 V2837    group     3 0.00547   0.00182        27.5   1.45e-4 7.00e-3 H_vs_M;…
#> 10 V457     group     3 0.000548  0.000183       52.2   1.34e-5 2.55e-3 H_vs_C;…
#> 11 V709     group     3 0.0771    0.0257         22.4   3.00e-4 1.22e-2 H_vs_M;…
#> 12 V919     group     3 0.00365   0.00122        31.9   8.46e-5 4.97e-3 H_vs_C;…
#> # ℹ 4 more variables: protein <chr>, protein_site <int>, trait <chr>,
#> #   gene <chr>
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
#> 2 hybrid  FALSE     2     0     0     0     1     0     1     4
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
#> # A tibble: 3,880 × 16
#>    variable protein glycan_composition glycan_structure protein_site gene  Tp   
#>    <chr>    <chr>   <comp>             <struct>                <int> <chr> <fct>
#>  1 V1       P08185  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          176 SERP… comp…
#>  2 V2       P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   hybr…
#>  3 V3       P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          344 HRG   comp…
#>  4 V4       P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#>  5 V5       P10909  Hex(6)HexNAc(5)    Hex(??-?)HexNAc…          291 CLU   comp…
#>  6 V6       P04196  Hex(5)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#>  7 V7       P04196  Hex(5)HexNAc(4)    Hex(??-?)HexNAc…          345 HRG   comp…
#>  8 V8       P04196  Hex(5)HexNAc(4)dH… dHex(??-?)Hex(?…          344 HRG   comp…
#>  9 V9       P04196  Hex(4)HexNAc(3)    Hex(??-?)HexNAc…          344 HRG   comp…
#> 10 V10      P04196  Hex(4)HexNAc(4)Ne… NeuAc(??-?)Hex(…          344 HRG   comp…
#> # ℹ 3,870 more rows
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
#> ℹ Expression matrix: 12 samples, 207 variables
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
#> [1] complex hybrid  complex complex complex
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
#> ℹ Normalizing data (Median Quotient)
#> ✔ Normalizing data (Median Quotient) [223ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [16ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size > 100, using MissForest imputation
#> ℹ Imputing missing values✔ Imputing missing values [6.9s]
#> 
#> ℹ Normalizing data (Total Area)
#> ✔ Normalizing data (Total Area) [15ms]
trait_exp <- derive_traits(exp)
trait_exp
#> 
#> ── Traitomics Experiment ───────────────────────────────────────────────────────
#> ℹ Expression matrix: 144 samples, 14 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: trait <chr>
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
