# glydet 0.3.0

## New features

* Trait functions created by `prop()`, `ratio()`, and `wmean()` can be printed nicely into the console.
* Add `explain_trait()` to provide a human-readable explanation of trait definitions.
* Add a `nF` built-in meta-property to count the number of fucoses. The `TF` trait is redefined as `prop(nF > 0)`.

## Minor improvements and bug fixes

* Emphasize in the documentation of `derive_traits()` and `derive_traits_()` that the "glycan_structure" column can be either a `glyrepr::glycan_structure()` vector or a character vector of glycan structure strings supported by `glyparse::auto_parse()`.
* `derive_traits()` and `derive_traits_()` now raise an error if the `trait_fns` parameter is not a named list.
* `derive_traits()` and `derive_traits_()` now raise an error if the `trait_fns` parameter is an empty list.
* A "Working with Glycomics Data" section is added to the "Get Started" vignette.
* A "Validating Trait Definitions" section is added to the "Defining Custom Traits" vignette.

# glydet 0.2.0

## New features

* `derive_traits()` and `add_meta_properties()` now validates if the `struc_col` column exists.
* Add an `overwrite` parameter to `add_meta_properties()` to deal with existing trait columns.

## Minor improvements and bug fixes

* Fix a serious bug in `derive_traits()`, where the row names of the expression matrix were not set correctly for glycomics data.
* Fix typos in documentations.

# glydet 0.1.1

## Minor improvements and bug fixes

* Update dependencies to depend on release versions of glycoverse packages.

# glydet 0.1.0

* First release.
