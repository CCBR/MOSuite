# Changelog

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
  - [`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_files.md)
  - [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md)
- [`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md)
- [`calc_cpm()`](https://ccbr.github.io/MOSuite/reference/calc_cpm.md)
  ([\#38](https://github.com/CCBR/MOSuite/issues/38))
- [`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md)
  ([\#38](https://github.com/CCBR/MOSuite/issues/38))
- [`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md)
  ([\#79](https://github.com/CCBR/MOSuite/issues/79))
- [`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md)
  ([\#82](https://github.com/CCBR/MOSuite/issues/82))
- [`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md)
  ([\#87](https://github.com/CCBR/MOSuite/issues/87))
- [`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md)
  ([\#102](https://github.com/CCBR/MOSuite/issues/102))
- [`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md)
  ([\#110](https://github.com/CCBR/MOSuite/issues/110))

#### visualization

- [`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md)
  ([\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md)
  ([\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md)
  ([\#90](https://github.com/CCBR/MOSuite/issues/90),
  [\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md)
  ([\#88](https://github.com/CCBR/MOSuite/issues/88),
  [\#96](https://github.com/CCBR/MOSuite/issues/96))
- [`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/reference/plot_volcano_enhanced.md)
  ([\#112](https://github.com/CCBR/MOSuite/issues/112))
- [`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/reference/plot_volcano_summary.md)
  ([\#112](https://github.com/CCBR/MOSuite/issues/112))
- [`plot_venn_diagram()`](https://ccbr.github.io/MOSuite/reference/plot_venn_diagram.md)
  ([\#111](https://github.com/CCBR/MOSuite/issues/111))

### vignettes

- `intro`
- `visualization`
- `memory`
- `renee`
