# Changelog

## glydet (development version)

## glydet 0.8.0

### New features

- [`make_trait()`](https://glycoverse.github.io/glydet/reference/make_trait.md)
  and
  [`explain_trait()`](https://glycoverse.github.io/glydet/reference/explain_trait.md)
  now have a `custom_mp` argument to allow users to define custom
  meta-properties.

### Minor improvements and bug fixes

- Update documentation of
  [`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
  to add a new type of `alignment`: “exact”, to be consistent with
  glymotif 0.12.0.
- Add reflection mechanism to
  [`make_trait()`](https://glycoverse.github.io/glydet/reference/make_trait.md)
  to make it more accurate.
- Add a `verbose` argument to
  [`make_trait()`](https://glycoverse.github.io/glydet/reference/make_trait.md)
  to control the amount of information printed to the console.

## glydet 0.7.0

We introduce exciting new features to `glydet`, including explaining and
creating derived traits using a Large Language Model (LLM). We hope
these features will help you work with custom derived traits more
efficiently with less effort.

### New features

- [`explain_trait()`](https://glycoverse.github.io/glydet/reference/explain_trait.md)
  now supports
  [`total()`](https://glycoverse.github.io/glydet/reference/total.md)
  and [`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md)
  traits.
- [`explain_trait()`](https://glycoverse.github.io/glydet/reference/explain_trait.md)
  now supports a new argument `use_ai` to use a Large Language Model
  (LLM) to explain the trait.
- Add
  [`make_trait()`](https://glycoverse.github.io/glydet/reference/make_trait.md)
  to create a derived trait function using natural language with a Large
  Language Model (LLM).

### Minor improvements and bug fixes

- The result of
  [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  and
  [`derive_traits_()`](https://glycoverse.github.io/glydet/reference/derive_traits_.md)
  now have an `explanation` column in the output tibble, explaining the
  derived traits.

## glydet 0.6.5

### Minor improvements and fixes

- glydet now depends on the CRAN version of glyparse.

## glydet 0.6.4

### Minor improvements and fixes

- glydet now depends on the CRAN version of glyrepr.

## glydet 0.6.3

### Minor improvements and bug fixes

- Updated dependency on `glymotif (>= 0.10.0)` to ensure compatibility
  with recent changes in package `igraph v2.2.0`.

## glydet 0.6.2

### Minor improvements and bug fixes

- Add `alignments` and `ignore_linkages` parameters to
  [`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md).
  These two parameters were left behind when we migrated
  [`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
  from `glymotif` to `glydet` in glydet 0.6.0.

## glydet 0.6.1

### Minor improvements and bug fixes

- Fix bugs introduced by the breaking changes in `glyexp` 0.10.0.

## glydet 0.6.0

### Breaking changes

- The first parameter of
  [`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md) is
  renamed from `val_cond` to `val`.
- The built-in trait `SG` is renamed to `GS`, and `GA` to `AG`, to match
  the convention that the last element of a trait name should be the
  focused property.

### New features

- Add two new trait factories:
  [`total()`](https://glycoverse.github.io/glydet/reference/total.md)
  and [`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md).
- Add
  [`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md).
  This function was once in the `glymotif` package, but we realized that
  motif quantification is essentially a special case of derived traits,
  so we reimplemented it in `glydet` with a more consistent and powerful
  interface.

### Minor improvements and bug fixes

- Fix some improperly named traits in the “Defining Custom Traits”
  vignette.
- Fix an error in a code snippet in the “Defining Custom Traits”
  vignette: `prop(nFc + nFa > 0)` -\> `prop((nFc + nFa) > 0)`.
- Add a vignette about quantifying glycan motifs.
- Update the introduction part in the “Get Started with glydet”
  vignette.
- All trait factories now have a `glydet_trait` super class.

## glydet 0.5.0

### Breaking changes

- The meta-property `T` is renamed to `Tp`, to avoid confusion with the
  R alias `T` for `TRUE`.

### New features

- Add `mp_cols` parameter to
  [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  as a new way to define custom meta-properties.

### Minor improvements and bug fixes

- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  and
  [`derive_traits_()`](https://glycoverse.github.io/glydet/reference/derive_traits_.md)
  now have separate documentations.
- Update the “Defining Custom Traits” vignette to include the newly
  added `mp_cols` parameter.
- Better error message when custom derived traits use undefined
  meta-properties.
- Add a section about how NAs can appear in the documentations of
  [`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
  [`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md),
  and
  [`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md).

## glydet 0.4.1

### Minor improvements and bug fixes

- Add another example to the “Defining Custom Traits” vignette.
- Use smaller datasets in examples to reduce build time.

## glydet 0.4.0

### Breaking changes

- The old
  [`all_traits()`](https://glycoverse.github.io/glydet/reference/all_traits.md)
  is renamed to
  [`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md).
- A new
  [`all_traits()`](https://glycoverse.github.io/glydet/reference/all_traits.md)
  is added, which includes more advanced derived traits with more
  detailed `within` conditions.
- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  now uses
  [`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md)
  by default.
- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  now returns a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with the “traitomics” type for glycomics data, and
  “traitproteomics” type for glycoproteomics data, instead of
  “traitomics” for all input.

### Minor improvements and bug fixes

- Introduce the new
  [`all_traits()`](https://glycoverse.github.io/glydet/reference/all_traits.md)
  in the Get Started vignette.
- Add examples to some functions.

## glydet 0.3.1

### Minor improvements and bug fixes

- Add `missForest` to the Suggests field in DESCRIPTION. This package is
  used in vignettes to impute missing values for glycomics data.

## glydet 0.3.0

### New features

- Trait functions created by
  [`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
  [`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md),
  and
  [`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
  can be printed nicely into the console.
- Add
  [`explain_trait()`](https://glycoverse.github.io/glydet/reference/explain_trait.md)
  to provide a human-readable explanation of trait definitions.
- Add a `nF` built-in meta-property to count the number of fucoses. The
  `TF` trait is redefined as `prop(nF > 0)`.

### Minor improvements and bug fixes

- Emphasize in the documentation of
  [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  and
  [`derive_traits_()`](https://glycoverse.github.io/glydet/reference/derive_traits_.md)
  that the “glycan_structure” column can be either a
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector or a character vector of glycan structure strings supported by
  [`glyparse::auto_parse()`](https://glycoverse.github.io/glyparse/reference/auto_parse.html).
- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  and
  [`derive_traits_()`](https://glycoverse.github.io/glydet/reference/derive_traits_.md)
  now raise an error if the `trait_fns` parameter is not a named list.
- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  and
  [`derive_traits_()`](https://glycoverse.github.io/glydet/reference/derive_traits_.md)
  now raise an error if the `trait_fns` parameter is an empty list.
- A “Working with Glycomics Data” section is added to the “Get Started”
  vignette.
- A “Validating Trait Definitions” section is added to the “Defining
  Custom Traits” vignette.

## glydet 0.2.0

### New features

- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  and
  [`add_meta_properties()`](https://glycoverse.github.io/glydet/reference/add_meta_properties.md)
  now validates if the `struc_col` column exists.
- Add an `overwrite` parameter to
  [`add_meta_properties()`](https://glycoverse.github.io/glydet/reference/add_meta_properties.md)
  to deal with existing trait columns.

### Minor improvements and bug fixes

- Fix a serious bug in
  [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md),
  where the row names of the expression matrix were not set correctly
  for glycomics data.
- Fix typos in documentations.

## glydet 0.1.1

### Minor improvements and bug fixes

- Update dependencies to depend on release versions of glycoverse
  packages.

## glydet 0.1.0

- First release.
