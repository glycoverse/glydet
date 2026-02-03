# Quantify Motifs in an Experiment

This function quantifies motifs from glycomic or glycoproteomic
profiles. For glycomics data, it calculates motif quantifications
directly. For glycoproteomics data, each glycosite is treated as a
separate glycome, and motif quantifications are calculated in a
site-specific manner.

The function takes a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object and returns a new
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with motif quantifications. Instead of containing quantifications
of individual glycans on each glycosite in each sample, the new
experiment contains quantifications of each motif on each glycosite in
each sample (for glycoproteomics data) or motif quantifications in each
sample (for glycomics data).

## Usage

``` r
quantify_motifs(
  exp,
  motifs,
  method = "relative",
  alignments = NULL,
  ignore_linkages = FALSE
)
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

- motifs:

  A character vector of motif names, glycan structure strings, a
  'glyrepr_structure' object, or a motif specification from
  [`glymotif::dynamic_motifs()`](https://glycoverse.github.io/glymotif/reference/dynamic_motifs.html)
  or
  [`glymotif::branch_motifs()`](https://glycoverse.github.io/glymotif/reference/branch_motifs.html).
  For glycan structure strings, all formats supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html)
  are accepted, including IUPAC-condensed, WURCS, GlycoCT, and others.
  If the vector is named, the names will be used as motif names.
  Otherwise, IUPAC-condensed structure strings will be used as motif
  names. For motif specifications, motifs are extracted automatically
  from the glycan structures in the experiment, and their
  IUPAC-condensed strings are used as motif names.

- method:

  A character string specifying the quantification method. Must be
  either "absolute" or "relative". Default is "relative". See "Relative
  and Absolute Motif Quantification" section for details.

- alignments:

  A character vector specifying the alignment method for each motif. Can
  be "terminal", "substructure", "core", "whole" or "exact". Default is
  "substructure". See
  [`glymotif::have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.html)
  for details.

- ignore_linkages:

  A logical value. If `TRUE`, linkages will be ignored in the
  comparison. See
  [`glymotif::have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.html)
  for details.

## Value

A new
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object containing motif quantifications. Instead of containing
quantifications of individual glycans on each glycosite in each sample,
the new experiment contains quantifications of each motif on each
glycosite in each sample (for glycoproteomics data) or motif
quantifications in each sample (for glycomics data).

The `var_info` table includes a `motif_structure` column containing the
parsed glycan structure for each motif, allowing traceability of motif
definitions.

For glycoproteomics data, with additional columns:

- `protein`: protein ID

- `protein_site`: the glycosite position on the protein

Other columns in the `var_info` table (e.g., `gene`) are retained if
they have a "many-to-one" relationship with glycosites (unique
combinations of `protein`, `protein_site`). That is, each glycosite
cannot have multiple values for these columns. `gene` is a common
example, as a glycosite can only be associated with one gene. Glycan
descriptions are not such columns, as a glycosite can have multiple
glycans, thus having multiple descriptions. Columns that do not have
this relationship with glycosites will be dropped. Don't worry if you
cannot understand this logic— just know that this function will do its
best to preserve useful information.

The `sample_info` and `meta_data` tables are not modified, except that
the `exp_type` field in `meta_data` is set to "traitomics" for glycomics
data and "traitproteomics" for glycoproteomics data.

## Relative and Absolute Motif Quantification

Motif quantification can be performed in two ways: absolute and
relative. Do not confuse this with the absolute and relative
quantification of glycans or glycopeptides.

In an omics context, absolute quantification means we can determine the
actual concentration (e.g., in mg/L) or counts (e.g., mRNA copy
numbers), while relative quantification means we can only compare the
abundance of the same molecule between different samples, but the
absolute values themselves are not meaningful. Label-free quantification
is a relative quantification method.

In the context of motif quantification, absolute and relative refer to
whether we normalize the motif quantifications by the total abundance of
the glycome (or the mini-glycomes of individual glycosites).

For example, let's say we have a glycomics dataset with two samples, A
and B. In each sample, three glycans (G1, G2, G3) are present, with the
following abundances:

- Sample A: G1 = 10, G2 = 20, G3 = 30

- Sample B: G1 = 20, G2 = 40, G3 = 60

For simplicity, let's say the motif we want to quantify happens to
appear once in all glycans.

For absolute motif quantification, we simply sum the abundances of the
motif across all glycans:

- Sample A: 10 × 1 + 20 × 1 + 30 × 1 = 60 (× 1 because the motif appears
  once in each glycan)

- Sample B: 20 × 1 + 40 × 1 + 60 × 1 = 120

The results are clearly different.

However, if we quantify the motif in a relative way:

- Sample A: (10 × 1 + 20 × 1 + 30 × 1) / (10 + 20 + 30) = 60 / 60 = 1

- Sample B: (20 × 1 + 40 × 1 + 60 × 1) / (20 + 40 + 60) = 120 / 120 = 1

The results are identical!

The absolute motif quantification answers the question "how many motifs
are there in the sample?", while the relative motif quantification
answers the question "if I take out one glycan molecule, how many motifs
are on it in average?"

So which method should you use? It depends on both your data type and
research objectives. Here are some general guidelines:

- For glycomics data, you should typically use relative motif
  quantification, because glycomics data is inherently compositional
  (search "compositional data" for more details).

- For glycoproteomics data, the choice depends on your research
  question. If you want to compare observed motif abundance across
  samples, use absolute motif quantification. Note that many factors can
  cause one sample to have higher values than another, such as
  upregulation of enzymes responsible for the motif, or higher overall
  glycosylation site occupancy. If you want to understand the underlying
  regulatory mechanisms, use relative motif quantification to correct
  for differences in site occupancy.

## Relationship with Derived Traits

Motif quantification is a special type of derived trait. It is simply a
weighted sum
([`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md)) or
weighted mean
([`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)) of
the motif counts, for absolute and relative motif quantification,
respectively. You can perform motif quantification manually using
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md).

Let's perform absolute motif quantification manually using
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
to better understand the process (we will use custom column
meta-properties here):

    # Add the meta-properties to the variable information tibble
    motifs <- c(
      nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
      nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
    )
    exp_with_mps <- glymotif::add_motifs_int(exp, motifs)

    # Define the traits
    trait_fns <- list(Lx = wsum(nLx), SLx = wsum(nSLx))

    # Calculate the traits
    derive_traits(exp_with_mps, trait_fns = trait_fns, mp_cols = c("nLx", "nSLx"))

The code snippet above is equivalent to:

    # Quantify the motifs
    quantify_motifs(exp, motifs, method = "absolute")

In fact, this is essentially the implementation of the
`quantify_motifs()` function, except the actual implementation is more
robust and user-friendly.

For relative motif quantification, simply replace
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md) with
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md), and
everything else remains the same.

## See also

[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md),
[`glymotif::have_motifs()`](https://glycoverse.github.io/glymotif/reference/have_motif.html)

## Examples

``` r
library(glyexp)
library(glyclean)

exp <- real_experiment |>
  auto_clean() |>
  slice_head_var(n = 10)
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

motifs <- c(
  nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
  nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
)

quantify_motifs(exp, motifs)
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 8 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, motif <chr>, gene <chr>, motif_structure <struct>

# Using dynamic motifs (auto-extracted from data)
quantify_motifs(exp, glymotif::dynamic_motifs(max_size = 3))
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 80 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, motif <chr>, gene <chr>, motif_structure <struct>

# Using branch motifs (auto-extracted from data)
quantify_motifs(exp, glymotif::branch_motifs())
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 20 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, motif <chr>, gene <chr>, motif_structure <struct>
```
