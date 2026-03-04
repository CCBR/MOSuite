#!/usr/bin/env bash
set -euo pipefail

# set MOSuite options for plots
export MOO_SAVE_PLOTS=TRUE
export MOO_PLOTS_DIR=./figures
mkdir -p $MOO_PLOTS_DIR

# add mosuite executable to the path
mosuite=$(R -s -e "cat(system.file('exec','mosuite', package='MOSuite'))")
export PATH="$PATH:$(dirname $mosuite)"

mosuite create_multiOmicDataSet_from_files --json=json_args/common/create_multiOmicDataSet_from_files.json
mosuite clean_raw_counts --json=json_args/common/clean_raw_counts.json
mosuite filter_counts --json=json_args/common/filter_counts.json
mosuite normalize_counts --json=json_args/common/normalize_counts.json
mosuite batch_correct_counts --json=json_args/common/batch_correct_counts.json
mosuite diff_counts --json=json_args/common/diff_counts.json
mosuite filter_diff --json=json_args/common/filter_diff.json
