# Changelog

## MOSuite 0.3.0

### New features

- New function
  [`write_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/dev/reference/write_multiOmicDataSet_properties.md).
  ([\#173](https://github.com/CCBR/MOSuite/issues/173),
  [@kelly-sovacool](https://github.com/kelly-sovacool))
  - Extracts all the properties from a multiOmicDataSet and writes any
    data frames as csv files, other objects are written as rds files.
- New utility functions to reduce code duplication across Code Ocean
  capsules: ([\#185](https://github.com/CCBR/MOSuite/issues/185),
  [@kelly-sovacool](https://github.com/kelly-sovacool))
  [`setup_capsule_environment()`](https://ccbr.github.io/MOSuite/dev/reference/setup_capsule_environment.md),
  [`load_moo_from_data_dir()`](https://ccbr.github.io/MOSuite/dev/reference/load_moo_from_data_dir.md),
  [`parse_optional_vector()`](https://ccbr.github.io/MOSuite/dev/reference/parse_optional_vector.md),
  [`parse_vector_with_default()`](https://ccbr.github.io/MOSuite/dev/reference/parse_vector_with_default.md),
  and
  [`parse_samples_to_rename()`](https://ccbr.github.io/MOSuite/dev/reference/parse_samples_to_rename.md).

### Bug fixes

- Fixed bug in
  [`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md)
  where duplicate gene rows were not being aggregated correctly.
  ([\#162](https://github.com/CCBR/MOSuite/issues/162),
  [@TJoshMeyer](https://github.com/TJoshMeyer))
- Fixed
  [`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md)
  to gracefully skip batch correction when only a single batch level
  exists; the function now warns without erroring.
  ([\#158](https://github.com/CCBR/MOSuite/issues/158),
  [@TJoshMeyer](https://github.com/TJoshMeyer))
- Fixed `multiOmicDataSet` validator to return character vector instead
  of using [`stop()`](https://rdrr.io/r/base/stop.html) per S7
  documentation. ([\#177](https://github.com/CCBR/MOSuite/issues/177),
  [@copilot](https://github.com/copilot),
  [@kelly-sovacool](https://github.com/kelly-sovacool))
- Fixed duplicate PCA figure files in
  [`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
  [`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
  and
  [`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md).
  ([\#180](https://github.com/CCBR/MOSuite/issues/180),
  [@kelly-sovacool](https://github.com/kelly-sovacool))
- Bug fixes: ([\#174](https://github.com/CCBR/MOSuite/issues/174),
  [@kelly-sovacool](https://github.com/kelly-sovacool))
  - Fixed bugs in
    [`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_summary.md),
    [`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_enhanced.md),
    and
    [`plot_pca_3d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_3d.md)
    when used with multiOmicDataSet objects.
  - Fixed bug in
    [`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md)
    when `filtering_mode = "all"` that was causing plot rendering
    errors.
  - Fixed
    [`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md)
    to save plots to disk correctly.
- Replaced deprecated
  [`arrange_()`](https://dplyr.tidyverse.org/reference/defunct-lazyeval.html)
  with [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)
  in
  [`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md).
  ([\#182](https://github.com/CCBR/MOSuite/issues/182),
  [@kelly-sovacool](https://github.com/kelly-sovacool))
- Improvements for use with Galaxy.
  ([\#168](https://github.com/CCBR/MOSuite/issues/168),
  [\#170](https://github.com/CCBR/MOSuite/issues/170),
  [\#171](https://github.com/CCBR/MOSuite/issues/171),
  [\#174](https://github.com/CCBR/MOSuite/issues/174),
  [@kelly-sovacool](https://github.com/kelly-sovacool))

## MOSuite 0.2.1

- A docker image is now available.
  ([\#134](https://github.com/CCBR/MOSuite/issues/134))
  - <https://hub.docker.com/r/nciccbr/mosuite>
- Minor documentation improvements.
  ([\#135](https://github.com/CCBR/MOSuite/issues/135))
- Improvements for use with Galaxy.
  ([\#149](https://github.com/CCBR/MOSuite/issues/149))
- Fixed bug where 3D PCA plots were not being saved.
  ([\#149](https://github.com/CCBR/MOSuite/issues/149))

## MOSuite 0.2.0

- Any user-facing function can now be called from the unix command line
  to support Galaxy.
  ([\#126](https://github.com/CCBR/MOSuite/issues/126),
  [\#127](https://github.com/CCBR/MOSuite/issues/127)) Usage:
  `mosuite [function] --json=path/to/args`
  - It is not recommended for most users to run MOSuite via the CLI;
    this is only intended for the Galaxy workflow.
- MOSuite is now archived in Zenodo with a DOI:
  [10.5281/zenodo.16371580](http://doi.org/10.5281/zenodo.16371580)

## MOSuite 0.1.0

This is the first release of MOSuite 🎉

- Note: at the start of development, this package was called reneeTools.
  Later it was renamed to MOSuite.
  ([\#76](https://github.com/CCBR/MOSuite/issues/76))

### Main functions & classes

- `multiOmicDataSet` ([\#16](https://github.com/CCBR/MOSuite/issues/16),
  [\#28](https://github.com/CCBR/MOSuite/issues/28))
  - [`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_files.md)
  - [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_dataframes.md)
- [`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md)
- [`calc_cpm()`](https://ccbr.github.io/MOSuite/dev/reference/calc_cpm.md)
  ([\#38](https://github.com/CCBR/MOSuite/issues/38))
- [`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md)
  ([\#38](https://github.com/CCBR/MOSuite/issues/38))
- [`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md)
  ([\#79](https://github.com/CCBR/MOSuite/issues/79))
- [`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md)
  ([\#82](https://github.com/CCBR/MOSuite/issues/82))
- [`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md)
  ([\#87](https://github.com/CCBR/MOSuite/issues/87))
- [`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md)
  ([\#102](https://github.com/CCBR/MOSuite/issues/102))
- [`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md)
  ([\#110](https://github.com/CCBR/MOSuite/issues/110))

#### visualization

- [`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md)
  ([\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md)
  ([\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md)
  ([\#90](https://github.com/CCBR/MOSuite/issues/90),
  [\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md)
  ([\#88](https://github.com/CCBR/MOSuite/issues/88),
  [\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_enhanced.md)
  ([\#112](https://github.com/CCBR/MOSuite/issues/112))
- [`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_summary.md)
  ([\#112](https://github.com/CCBR/MOSuite/issues/112))
- [`plot_venn_diagram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_venn_diagram.md)
  ([\#111](https://github.com/CCBR/MOSuite/issues/111))

### vignettes

- `intro`
- `visualization`
- `memory`
- `renee`
