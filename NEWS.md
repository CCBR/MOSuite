# MOSuite (development version)

# MOSuite 0.2.0

- Any user-facing function can now be called from the unix command line.  (#126, #127)
  Usage: `mosuite [function] --json=path/to/args`
  See the new `cli` vignette for detailed instructions.
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
