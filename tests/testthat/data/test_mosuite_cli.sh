#!/usr/bin/env bash
set -euo pipefail

# set MOSuite options for plots
export MOO_SAVE_PLOTS=TRUE
export MOO_PLOTS_DIR=tests/testthat/data/figures
mkdir -p $MOO_PLOTS_DIR

# add mosuite executable to the path
mosuite=$(R -s -e "cat(system.file('exec','mosuite', package='MOSuite'))")
export PATH="$PATH:$(dirname $mosuite)"

mosuite create_multiOmicDataSet_from_files --json=tests/testthat/data/create_multiOmicDataSet_from_files.json
mosuite clean_raw_counts --json=tests/testthat/data/clean_raw_counts.json
mosuite filter_counts --json=tests/testthat/data/filter_counts.json
mosuite normalize_counts --json=tests/testthat/data/normalize_counts.json
mosuite batch_correct_counts --json=tests/testthat/data/batch_correct_counts.json
mosuite diff_counts --json=tests/testthat/data/diff_counts.json
mosuite filter_diff --json=tests/testthat/data/filter_diff.json
