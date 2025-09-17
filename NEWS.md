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
