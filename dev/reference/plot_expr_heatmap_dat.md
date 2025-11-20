# Plot expression heatmap for counts dataframe

Plot expression heatmap for counts dataframe

## Arguments

- moo_counts:

  counts dataframe (**Required**)

- sample_metadata:

  sample metadata as a data frame or tibble (**Required**)

- sample_id_colname:

  The column from the sample metadata containing the sample names. The
  names in this column must exactly match the names used as the sample
  column names of your input Counts Matrix. (Default: `NULL` - first
  column in the sample metadata will be used.)

- feature_id_colname:

  The column from the counts dataa containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- group_colname:

  The column from the sample metadata containing the sample group
  information. This is usually a column showing to which experimental
  treatments each sample belongs (e.g. WildType, Knockout, Tumor,
  Normal, Before, After, etc.).

- label_colname:

  The column from the sample metadata containing the sample labels as
  you wish them to appear in the plots produced by this template. This
  can be the same Sample Names Column. However, you may desire different
  labels to display on your figure (e.g. shorter labels are sometimes
  preferred on plots). In that case, select the column with your
  preferred Labels here. The selected column should contain unique names
  for each sample. (Default: `NULL` – `sample_id_colname` will be used.)

- color_values:

  vector of colors as hex values or names recognized by R

- samples_to_include:

  Which samples would you like to include? Usually, you will choose all
  sample columns, or you could choose to remove certain samples. Samples
  excluded here will be removed in this step and from further analysis
  downstream of this step. (Default: `NULL` - all sample IDs in
  `moo@sample_meta` will be used.)

- include_all_genes:

  Set to TRUE if all genes are to be included. Set to FALSE if you want
  to filter genes by variance and/or provide a list of specific genes
  that will appear in the heatmap.

- filter_top_genes_by_variance:

  Set to TRUE if you want to only include the top genes by variance. Set
  to FALSE if you do not want to filter genes by variance.

- top_genes_by_variance_to_include:

  The number of genes to include if filtering genes by variance. This
  parameter is ignored if "Filter top genes by variance" is set to
  FALSE.

- specific_genes_to_include_in_heatmap:

  Enter the gene symbols to be included in the heatmap, with each gene
  symbol separated with a space from the others. Alternatively, paste in
  a column of gene names from any spreadsheet application. This
  parameter is ignored if "Include all genes" is set to TRUE.

- cluster_genes:

  Choose whether to cluster the rows (genes). If TRUE, rows will have
  clustering applied. If FALSE, clustering will not be applied to rows.

- gene_distance_metric:

  Distance metric to be used in clustering genes. (TODO document
  options)

- gene_clustering_method:

  Clustering method metric to be used in clustering samples. (TODO
  document options)

- display_gene_dendrograms:

  Set to TRUE to show gene dendrograms. Set to FALSE to hide
  dendrograms.

- display_gene_names:

  Set to TRUE to display gene names on the right side of the heatmap.
  Set to FALSE to hide gene names.

- center_and_rescale_expression:

  Center and rescale expression for each gene across all included
  samples.

- cluster_samples:

  Choose whether to cluster the columns (samples). If TRUE, columns will
  have clustering applied. If FALSE, clustering will not be applied to
  columns.

- arrange_sample_columns:

  If TRUE, arranges columns by annotation groups. If FALSE, and "Cluster
  Samples" is FALSE, samples will appear in the order of input (samples
  to include)

- order_by_gene_expression:

  If TRUE, set gene name below and direction for ordering

- gene_to_order_columns:

  Gene to order columns by expression levels

- gene_expression_order:

  Choose direction for gene order

- smpl_distance_metric:

  Distance metric to be used in clustering samples. (TODO document
  options)

- smpl_clustering_method:

  Clustering method to be used in clustering samples. (TODO document
  options)

- display_smpl_dendrograms:

  Set to TRUE to show sample dendrograms. Set to FALSE to hide
  dendrogram.

- reorder_dendrogram:

  If TRUE, set the order of the dendrogram (below)

- reorder_dendrogram_order:

  Reorder the samples (columns) of the dendrogram by name, e.g.
  “sample2”,“sample3",“sample1".

- display_sample_names:

  Set to TRUE if you want sample names to be displayed on the plot. Set
  to FALSE to hide sample names.

- group_columns:

  Columns containing the sample groups for annotation tracks

- assign_group_colors:

  If TRUE, set the groups assigned colors (below)

- assign_color_to_sample_groups:

  Enter each sample to color in the format: group_name: color This
  parameter is ignored if "Assign Colors" is set to FALSE.

- group_colors:

  Set group annotation colors.

- heatmap_color_scheme:

  color scheme (TODO document options)

- autoscale_heatmap_color:

  Set to TRUE to autoscale the heatmap colors between the maximum and
  minimum heatmap color parameters. If FALSE, set the heatmap colors
  between "Set max heatmap color" and "Set min heatmap color" (below).

- set_min_heatmap_color:

  If Autoscale heatmap color is set to FALSE, set the minimum heatmap
  z-score value

- set_max_heatmap_color:

  If Autoscale heatmap color is set to FALSE, set the maximum heatmap
  z-score value.

- aspect_ratio:

  Set figure Aspect Ratio. Ratio refers to entire figure including
  legend. If set to Auto figure size is based on number of rows and
  columns form counts matrix. default - Auto

- legend_font_size:

  Set Font size for figure legend. Default is 10.

- gene_name_font_size:

  Font size for gene names. If you don't want gene labels to show,
  toggle "Display Gene Names" below to FALSE

- sample_name_font_size:

  Font size for sample names. If you don't want to display samples
  names, toggle "Display sample names" (below) to FALSE

- display_numbers:

  Setting to FALSE (default) will not display numerical value of heat on
  heatmap. Set to TRUE if you want to see these numbers on the plot.

## See also

[plot_expr_heatmap](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md)
generic

Other heatmaps:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md)

Other plotters for counts dataframes:
[`plot_corr_heatmap_dat`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap_dat.md),
[`plot_histogram_dat`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram_dat.md),
[`plot_pca_dat`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_dat.md),
[`plot_read_depth_dat`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth_dat.md)
