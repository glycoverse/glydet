# Quantifying Glycan Motifs

In this vignette, we’ll explore a special type of derived trait: motif
quantification. This powerful feature allows us to measure the abundance
of biologically meaningful glycan substructures across your samples.

``` r
library(glydet)
library(glyexp)
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:stats':
#> 
#>     aggregate

exp <- auto_clean(real_experiment)
#> ℹ Normalizing data (Median)
#> ✔ Normalizing data (Median) [149ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [79ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size <= 30, using sample minimum imputation
#> ℹ Imputing missing values✔ Imputing missing values [30ms]
#> 
#> ℹ Aggregating data
#> ✔ Aggregating data [1.2s]
#> 
#> ℹ Normalizing data again
#> ✔ Normalizing data again [19ms]
```

## What is motif quantification?

Glycan motifs are substructures with special biological or structural
significance. Think of them as functional building blocks that carry
specific biological messages. For example, the Lewis x antigen is a
fucosylated carbohydrate epitope commonly found on glycoproteins and
glycolipids. It plays crucial roles in cell–cell recognition, adhesion,
and immune responses. Understanding the abundance of Lewis x antigens
can provide valuable insights into immune system dynamics.

But here’s the challenge: how do we quantify Lewis x antigens across
samples? Different glycans contain varying numbers of Lewis x motifs,
and these glycans themselves have different abundances across samples.
Our solution elegantly combines both pieces of information to estimate
the true abundance of Lewis x antigens in your dataset.

## The `quantify_motifs()` function

Glydet provides the
[`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
function to handle this complex task seamlessly. This function takes a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object along with your motifs of interest and returns a new experiment
object enriched with motif quantifications. Just like
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md),
[`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
works beautifully with both glycomics and glycoproteomics data.

Let’s dive into a practical example:

``` r
# Define our motifs of interest using IUPAC-condensed format
motifs <- c(
  Lx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
  SLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
)

# Quantify the motifs in our dataset
motif_exp <- quantify_motifs(exp, motifs)
motif_exp
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 548 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, motif <chr>, gene <chr>
```

``` r
get_var_info(motif_exp)
#> # A tibble: 548 × 5
#>    variable protein protein_site motif gene  
#>    <chr>    <chr>          <int> <chr> <chr> 
#>  1 V1       A6NJW9            49 Lx    CD8B2 
#>  2 V2       A6NJW9            49 SLx   CD8B2 
#>  3 V3       O14786           150 Lx    NRP1  
#>  4 V4       O14786           150 SLx   NRP1  
#>  5 V5       O43866           226 Lx    CD5L  
#>  6 V6       O43866           226 SLx   CD5L  
#>  7 V7       O75437           244 Lx    ZNF254
#>  8 V8       O75437           244 SLx   ZNF254
#>  9 V9       O75581           281 Lx    LRP6  
#> 10 V10      O75581           281 SLx   LRP6  
#> # ℹ 538 more rows
```

Notice how the variable information now includes a `motif` column
instead of the usual `trait` column. Don’t worry—this is just a cosmetic
difference. Under the hood, the functionality works exactly the same way
as
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md).

## Absolute vs. relative motif quantification

Here’s where things get a little bit complex: motif quantification can
be performed in two distinct ways—absolute and relative. By default,
[`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
uses relative quantification, but understanding both approaches is
crucial for choosing the right method for your research.

**Important note:** Don’t confuse this with the absolute and relative
quantification of glycans or glycopeptides themselves— that’s a
different concept entirely.

In traditional omics contexts, “absolute quantification” means
determining actual concentrations (like mg/L) or counts (like mRNA copy
numbers), while “relative quantification” means comparing abundance
between samples without meaningful absolute values. Label-free
quantification, for instance, is a relative method.

However, in motif quantification, “absolute” and “relative” refer to
whether we normalize motif abundances by the total abundance of the
glycome (or individual glycosite mini-glycomes).

Let’s illustrate this with a concrete example. Imagine a glycomics
dataset with two samples, A and B, each containing three glycans (G1,
G2, G3) with these abundances:

- **Sample A:** G1 = 10, G2 = 20, G3 = 30
- **Sample B:** G1 = 20, G2 = 40, G3 = 60

For simplicity, let’s assume our target motif appears exactly once in
each glycan.

### Absolute motif quantification

We simply sum the motif abundances across all glycans:

- **Sample A:** 10 × 1 + 20 × 1 + 30 × 1 = 60
- **Sample B:** 20 × 1 + 40 × 1 + 60 × 1 = 120

The results clearly differ between samples.

### Relative motif quantification

We normalize by the total glycan abundance:

- **Sample A:** (10 × 1 + 20 × 1 + 30 × 1) ÷ (10 + 20 + 30) = 60 ÷ 60 =
  1
- **Sample B:** (20 × 1 + 40 × 1 + 60 × 1) ÷ (20 + 40 + 60) = 120 ÷ 120
  = 1

The results are identical in two samples!

This distinction matters because these methods answer different
biological questions:

- **Absolute quantification** asks: “How many motifs are present in this
  sample?”
- **Relative quantification** asks: “If I randomly select one glycan
  molecule, how many motifs would I expect to find on it in average?”

### Choosing the right approach

So which method should you use? The answer depends on both your data
type and research objectives. Here are some practical guidelines:

**For glycomics data:** Use relative motif quantification in most cases,
since glycomics data is inherently compositional (look up “compositional
data analysis” if you’re curious about the statistical theory behind
this).

**For glycoproteomics data:**

- **Choose absolute quantification** when you want to compare observed
  motif abundance across samples. Keep in mind that many factors can
  cause differences between samples, including upregulation of enzymes
  responsible for the motif or changes in overall glycosylation site
  occupancy.
- **Choose relative quantification** when you want to understand
  underlying regulatory mechanisms by correcting for differences in site
  occupancy. This approach helps isolate the biological signal from
  technical variation.

## Connection to derived traits

You might wonder why
[`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
lives in the `glydet` package rather than `glymotif`. There’s actually
an interesting history here! Before `glymotif` v0.9.0, there was indeed
a
[`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
function in `glymotif`. However, we eventually realized that motif
quantification is fundamentally a special case of derived traits, so we
reimplemented it in `glydet` with a more consistent and powerful
interface.

To demonstrate this connection, let’s manually perform motif
quantification using
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md).
This will help you understand that motif quantification is essentially
just a specialized application of trait derivation:
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md) traits
for absolute quantification and
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
traits for relative quantification.

``` r
# First, add the meta-properties to the variable information
motifs <- c(
  nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
  nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
)
exp_with_mps <- glymotif::add_motifs_int(exp, motifs)

# Define the traits using wsum() for absolute quantification
trait_fns <- list(Lx = wsum(nLx), SLx = wsum(nSLx))

# Calculate the traits
derive_traits(exp_with_mps, trait_fns = trait_fns, mp_cols = c("nLx", "nSLx"))
```

This code snippet is functionally equivalent to:

``` r
# The much simpler approach using quantify_motifs()
quantify_motifs(exp, motifs, method = "absolute")
```

Pretty neat, right? The
[`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
function is essentially doing all the heavy lifting shown in the manual
approach, but with additional robustness and user-friendly features.

For relative motif quantification, you’d simply replace
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md) with
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md) in
the manual approach, while everything else stays the same. This
flexibility showcases the power of the underlying trait derivation
system.
