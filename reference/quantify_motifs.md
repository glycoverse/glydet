# Quantify Motifs in an Experiment

This function quantifies motifs from glycomic or glycoproteomic
profiles. For glycomics data, it calculates motif quantifications
directly. For glycoproteomics data, each glycosite is treated as a
separate glycome, and motif quantifications are calculated in a
site-specific manner.

The function takes a
[`glyexp::GlycomicSE`](https://glycoverse.github.io/glyexp/reference/GlycomicSE.html)
or
[`glyexp::GlycoproteomicSE`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
object and returns a plain `SummarizedExperiment` with motif
quantifications. Instead of containing quantifications of individual
glycans on each glycosite in each sample, the output assay contains
quantifications of each motif on each glycosite in each sample (for
glycoproteomics data) or motif quantifications in each sample (for
glycomics data).

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

New `GlycomicSE` and `GlycoproteomicSE` inputs return a plain
`SummarizedExperiment`; compatible legacy glyexp inputs preserve their
legacy container type. The output contains motif quantifications.
Instead of containing quantifications of individual glycans on each
glycosite in each sample, the output assay contains quantifications of
each motif on each glycosite in each sample (for glycoproteomics data)
or motif quantifications in each sample (for glycomics data).

[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
includes a `trait` column with motif names and a `motif_structure`
column containing the parsed glycan structure for each motif, allowing
traceability of motif definitions.

For glycoproteomics data, with additional columns:

- `protein`: protein ID

- `protein_site`: the glycosite position on the protein

Other columns in
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
(e.g., `gene`) are retained if they have a "many-to-one" relationship
with glycosites (unique combinations of `protein`, `protein_site`). That
is, each glycosite cannot have multiple values for these columns. `gene`
is a common example, as a glycosite can only be associated with one
gene. Glycan descriptions are not such columns, as a glycosite can have
multiple glycans, thus having multiple descriptions. Columns that do not
have this relationship with glycosites will be dropped. Don't worry if
you cannot understand this logic— just know that this function will do
its best to preserve useful information.

[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
and `metadata()` are not modified, except that the `exp_type` field in
`metadata()` is set to "traitomics" for glycomics data and
"traitproteomics" for glycoproteomics data.

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
    exp_with_mps <- exp |>
      glyexp::mutate_row(
        tibble::as_tibble(glymotif::count_motifs(glycan_structure, motifs))
      )

    # Define the traits
    trait_fns <- list(Lx = wsum(nLx), SLx = wsum(nSLx))

    # Calculate the traits
    derive_traits(exp_with_mps, trait_fns = trait_fns)

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
library(SummarizedExperiment)
library(glyclean)

gp_se <- real_experiment |>
  auto_clean() |>
  slice_head_row(n = 10)
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

motifs <- c(
  nLx = "Hex(??-?)[dHex(??-?)]HexNAc(??-",  # Lewis x antigen
  nSLx = "NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-"  # Sialyl Lewis x antigen
)

motif_se <- quantify_motifs(gp_se, motifs)
rowData(motif_se)
#> DataFrame with 8 rows and 5 columns
#>                     protein protein_site       trait        gene
#>                 <character>    <integer> <character> <character>
#> P04196-344-nLx       P04196          344         nLx         HRG
#> P04196-344-nSLx      P04196          344        nSLx         HRG
#> P04196-345-nLx       P04196          345         nLx         HRG
#> P04196-345-nSLx      P04196          345        nSLx         HRG
#> P08185-176-nLx       P08185          176         nLx    SERPINA6
#> P08185-176-nSLx      P08185          176        nSLx    SERPINA6
#> P10909-291-nLx       P10909          291         nLx         CLU
#> P10909-291-nSLx      P10909          291        nSLx         CLU
#>                        motif_structure
#>                    <glyrepr_structure>
#> P04196-344-nLx  dHex(??-?)[Hex(??-?)..
#> P04196-344-nSLx NeuAc(??-?)Hex(??-?)..
#> P04196-345-nLx  dHex(??-?)[Hex(??-?)..
#> P04196-345-nSLx NeuAc(??-?)Hex(??-?)..
#> P08185-176-nLx  dHex(??-?)[Hex(??-?)..
#> P08185-176-nSLx NeuAc(??-?)Hex(??-?)..
#> P10909-291-nLx  dHex(??-?)[Hex(??-?)..
#> P10909-291-nSLx NeuAc(??-?)Hex(??-?)..
assay(motif_se)
#>                 C1 C2 C3 H1 H2 H3 M1 M2 M3 Y1 Y2 Y3
#> P04196-344-nLx   0  0  0  0  0  0  0  0  0  0  0  0
#> P04196-344-nSLx  0  0  0  0  0  0  0  0  0  0  0  0
#> P04196-345-nLx   0  0  0  0  0  0  0  0  0  0  0  0
#> P04196-345-nSLx  0  0  0  0  0  0  0  0  0  0  0  0
#> P08185-176-nLx   0  0  0  0  0  0  0  0  0  0  0  0
#> P08185-176-nSLx  0  0  0  0  0  0  0  0  0  0  0  0
#> P10909-291-nLx   0  0  0  0  0  0  0  0  0  0  0  0
#> P10909-291-nSLx  0  0  0  0  0  0  0  0  0  0  0  0
colData(motif_se)
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

# Using dynamic motifs (auto-extracted from data)
quantify_motifs(gp_se, glymotif::dynamic_motifs(max_size = 3))
#> class: SummarizedExperiment 
#> dim: 80 12 
#> metadata(3): exp_type glycan_type quant_method
#> assays(1): abundance
#> rownames(80): P04196-344-NeuAc(??- P04196-344-Hex(??- ...
#>   P10909-291-dHex(??-?)Hex(??- P10909-291-dHex(??-?)Hex(??-?)HexNAc(??-
#> rowData names(5): protein protein_site trait gene motif_structure
#> colnames(12): C1 C2 ... Y2 Y3
#> colData names(1): group

# Using branch motifs (auto-extracted from data)
quantify_motifs(gp_se, glymotif::branch_motifs())
#> class: SummarizedExperiment 
#> dim: 20 12 
#> metadata(3): exp_type glycan_type quant_method
#> assays(1): abundance
#> rownames(20): P04196-344-NeuAc(??-?)Hex(??-?)HexNAc(??-
#>   P04196-344-HexNAc(??- ...
#>   P10909-291-Hex(??-?)HexNAc(??-?)Hex(??-?)HexNAc(??-
#>   P10909-291-dHex(??-?)Hex(??-?)HexNAc(??-
#> rowData names(5): protein protein_site trait gene motif_structure
#> colnames(12): C1 C2 ... Y2 Y3
#> colData names(1): group
```
