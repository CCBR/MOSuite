# MOSuite development version

- Fixed duplicate PCA figure files in `filter_counts()`, `normalize_counts()`, and `batch_correct_counts()`.
- Fixed bug in `clean_raw_counts()` where duplicate gene rows were not being aggregated correctly. (#162, @TJoshMeyer)
- New function `write_multiOmicDataSet_properties()`. (#173, @kelly-sovacool)
  - Extracts all the properties from a multiOmicDataSet and writes any data frames as csv files, other objects are written as rds files.
- Bug fixes: (#174, @kelly-sovacool)
  - Fixed bugs in `plot_volcano_summary()`, `plot_volcano_enhanced()`, and `plot_pca_3d()` when used with multiOmicDataSet objects.
  - Fixed bug in `filter_diff()` when `filtering_mode = "all"` that was causing plot rendering errors.
  - Fixed `plot_pca_2d()` to save plots to disk correctly.
- Improvements for use with Galaxy. (#168, #170, #171, #174, @kelly-sovacool)

# MOSuite 0.2.1

- A docker image is now available. (#134)
  - <https://hub.docker.com/r/nciccbr/mosuite>
- Minor documentation improvements. (#135)
- Improvements for use with Galaxy. (#149)
- Fixed bug where 3D PCA plots were not being saved. (#149)

# MOSuite 0.2.0

- Any user-facing function can now be called from the unix command line to support Galaxy.  (#126, #127)
  Usage: `mosuite [function] --json=path/to/args`
  - It is not recommended for most users to run MOSuite via the CLI; this is only intended for the Galaxy workflow.
- MOSuite is now archived in Zenodo with a DOI: [10.5281/zenodo.16371580](http://doi.org/10.5281/zenodo.16371580)

# MOSuite 0.1.0

This is the first release of MOSuite 🎉

- Note: at the start of development, this package was called reneeTools.
  Later it was renamed to MOSuite. (#76)

## Main functions & classes

- `multiOmicDataSet` (#16, #28)
  - `create_multiOmicDataSet_from_files()`
  - `create_multiOmicDataSet_from_dataframes()`
- `run_deseq2()`
- `calc_cpm()` (#38)
- `filter_counts()` (#38)
- `clean_raw_counts()` (#79)
- `normalize_counts()` (#82)
- `batch_correct_counts()` (#87)
- `diff_counts()` (#102)
- `filter_diff()` (#110)

### visualization

- `plot_histogram()` (#96)
- `plot_corr_heatmap()` (#96)
- `plot_expr_heatmap()` (#90, #96)
- `plot_pca()` (#88, #96)
- `plot_volcano_enhanced()` (#112)
- `plot_volcano_summary()` (#112)
- `plot_venn_diagram()` (#111)

## vignettes

- `intro`
- `visualization`
- `memory`
- `renee`
