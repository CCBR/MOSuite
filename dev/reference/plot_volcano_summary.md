# Volcano Plot - Summary

Produces one volcano plot for each tested contrast in the input DEG
table. It can be sorted by either fold change, t-statistic, or p-value.
The returned dataset includes one row for each significant gene in each
contrast, and contains columns from the DEG analysis of that contrast as
well as columns useful to the Venn diagram template downstream.

## Usage

``` r
plot_volcano_summary(
  moo_diff,
  feature_id_colname = NULL,
  signif_colname = "pval",
  signif_threshold = 0.05,
  change_threshold = 1,
  value_to_sort_the_output_dataset = "t-statistic",
  num_features_to_label = 30,
  add_features = FALSE,
  label_features = FALSE,
  custom_gene_list = "",
  default_label_color = "black",
  custom_label_color = "green3",
  label_x_adj = 0.2,
  label_y_adj = 0.2,
  line_thickness = 0.5,
  label_font_size = 4,
  label_font_type = 1,
  displace_feature_labels = FALSE,
  custom_gene_list_special_label_displacement = "",
  special_label_displacement_x_axis = 2,
  special_label_displacement_y_axis = 2,
  color_of_signif_threshold_line = "blue",
  color_of_non_significant_features = "black",
  color_of_logfold_change_threshold_line = "red",
  color_of_features_meeting_only_signif_threshold = "lightgoldenrod2",
  color_for_features_meeting_pvalue_and_foldchange_thresholds = "red",
  flip_vplot = FALSE,
  use_default_x_axis_limit = TRUE,
  x_axis_limit = 5,
  use_default_y_axis_limit = TRUE,
  y_axis_limit = 10,
  point_size = 2,
  add_deg_columns = c("FC", "logFC", "tstat", "pval", "adjpval"),
  graphics_device = grDevices::png,
  image_width = 15,
  image_height = 15,
  dpi = 300,
  use_default_grid_layout = TRUE,
  number_of_rows_in_grid_layout = 1,
  aspect_ratio = 0,
  plot_filename = "volcano_summary.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff"
)
```

## Arguments

- moo_diff:

  Differential expression analysis result from one or more contrasts.
  This must be a dataframe.

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- signif_colname:

  column name of significance values (e.g., adjusted p-values or FDR).
  This column will be used to determine which points are considered
  significant in the volcano plot.

- signif_threshold:

  Numeric value specifying the significance cutoff for p-values (i.e.
  filters on `signif_colname`)

- change_threshold:

  Numeric value specifying the fold change cutoff for significance (i.e.
  filters on `change_colname`)

- value_to_sort_the_output_dataset:

  How to sort the output dataset. Options are "fold-change" or
  "p-value".

- num_features_to_label:

  Number of top features/genes to label in the volcano plot. Default is
  30.

- add_features:

  Add custom_gene_list To Labels. Set TRUE when you want to label a
  specific set of features (features) in the "custom_gene_list"
  parameter" IN ADDITION to the number of features you set in the
  "Number of Features to Label" parameter.

- label_features:

  Select TRUE when you want to label ONLY a specific list of
  features(features) given in the "custom_gene_list" parameter.

- custom_gene_list:

  Provide a list of features (comma separated) to be labeled on the
  volcano plot. You must toggle one of the following ON to see these
  labels: "Add features" or "Label Only My Feature List".

- default_label_color:

  Set the color for the text used to add feature (gene) name labels to
  points.

- custom_label_color:

  Set the color for the specific list of features (features) provided in
  the "Feature List" parameter.

- label_x_adj:

  adjust position of the labels on the x-axis. Default: 0.2

- label_y_adj:

  adjust position of the labels on the y-axis. Default: 0.2

- line_thickness:

  Set the thickness of the lines in the plot. Default: 0.5

- label_font_size:

  Set the font size of the labels. Default: 4

- label_font_type:

  Set the font type of the labels. Default: 1

- displace_feature_labels:

  Set to TRUE to displace gene labels. Default: FALSE. Set TRUE if you
  want to displace the feature (gene) label for a specific set of
  features. Make sure to use custom x- and y- limits and give sufficient
  space for displacement; otherwise other labels than the desired ones
  will appear displaced.

- custom_gene_list_special_label_displacement:

  Provide a list of features (comma separated) for which you want
  special displacement of the feature label.

- special_label_displacement_x_axis:

  Displacement of the feature label on the x-axis. Default: 2

- special_label_displacement_y_axis:

  Displacement of the feature label on the y-axis. Default: 2

- color_of_signif_threshold_line:

  Color of the significance threshold line. Default: "blue"

- color_of_non_significant_features:

  Color of the non-significant features. Default: "black"

- color_of_logfold_change_threshold_line:

  Color of the log fold change threshold line. Default: "red"

- color_of_features_meeting_only_signif_threshold:

  Color of the features that meet only the significance threshold.
  Default: "lightgoldenrod2"

- color_for_features_meeting_pvalue_and_foldchange_thresholds:

  Color of the features that meet both the p-value and fold change
  thresholds. Default: "red"

- flip_vplot:

  Set to TRUE to flip the fold change values so that the volcano plot
  looks like a comparison was B-A. Default: FALSE

- use_default_x_axis_limit:

  Set to TRUE to use the default x-axis limit. Default: TRUE

- x_axis_limit:

  Custom x-axis limit. Default: c(-5, 5)

- use_default_y_axis_limit:

  Set to TRUE to use the default y-axis limit. Default: TRUE

- y_axis_limit:

  Custom y-axis limit. Default: c(0, 10)

- point_size:

  Size of the points in the plot. Default: 1

- add_deg_columns:

  Add additional columns from the DEG analysis to the output dataset.
  Default: FALSE

- graphics_device:

  passed to `ggsave(device)`. Default:
  [`grDevices::png`](https://rdrr.io/r/grDevices/png.html)

- image_width:

  output image width in pixels - only used if save_plots is TRUE

- image_height:

  output image height in pixels - only used if save_plots is TRUE

- dpi:

  dots-per-inch of the output image (see `ggsave()`) - only used if
  save_plots is TRUE

- use_default_grid_layout:

  Set to TRUE to use the default grid layout. Default: TRUE

- number_of_rows_in_grid_layout:

  Number of rows in the grid layout. Default: 1

- aspect_ratio:

  Aspect ratio of the output image. Default: 4/3

- plot_filename:

  Filename for the output plot. Default: "volcano_plot.png"

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `TRUE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- plots_subdir:

  subdirectory in `figures/` where plots will be saved if `save_plots`
  is `TRUE`

## Examples

``` r
plot_volcano_summary(nidap_deg_analysis, print_plots = TRUE)
#> Preparing table for contrast: B-A
#> Fold change column: B-A_logFC
#> pval column: B-A_pval
#> Total number of features included in volcano plot: 7943
#> Preparing table for contrast: C-A
#> Fold change column: C-A_logFC
#> pval column: C-A_pval
#> Total number of features included in volcano plot: 7943
#> Preparing table for contrast: B-C
#> Fold change column: B-C_logFC
#> pval column: B-C_pval
#> Total number of features included in volcano plot: 7943
#> Warning: ggrepel: 18 unlabeled data points (too many overlaps). Consider increasing max.overlaps
#> Warning: ggrepel: 27 unlabeled data points (too many overlaps). Consider increasing max.overlaps
#> Warning: ggrepel: 27 unlabeled data points (too many overlaps). Consider increasing max.overlaps

#> Saving 6.67 x 6.67 in image
#> Warning: ggrepel: 18 unlabeled data points (too many overlaps). Consider increasing max.overlaps
#> Warning: ggrepel: 27 unlabeled data points (too many overlaps). Consider increasing max.overlaps
#> Warning: ggrepel: 27 unlabeled data points (too many overlaps). Consider increasing max.overlaps
#>                   Gene Contrast          FC     logFC      tstat         pval
#> B-A.1             Dntt      B-A  -42.746586 -5.417737 -15.687975 3.159343e-09
#> B-A.2           Tmsb4x      B-A    3.850020  1.944866  12.910261 2.760555e-08
#> B-A.3             Flt3      B-A   -7.714394 -2.947553 -11.380840 1.093405e-07
#> B-A.4          Tspan13      B-A   -7.038498 -2.815268 -11.031274 1.531110e-07
#> B-A.5            Tapt1      B-A   -5.291816 -2.403763 -10.658467 2.214593e-07
#> B-A.6            Itgb7      B-A    8.873823  3.149556  10.561474 2.442070e-07
#> B-A.7            Nfkb1      B-A    6.153137  2.621322   9.778964 5.527550e-07
#> B-A.8             Bcl2      B-A    5.036402  2.332394   9.189193 1.059895e-06
#> B-A.9            Itga4      B-A    4.306962  2.106671   9.170168 1.082972e-06
#> B-A.10            Cnn3      B-A   -8.514926 -3.089994  -8.891098 1.491360e-06
#> B-A.11             Id2      B-A   71.935840  6.168639   8.385761 2.715658e-06
#> B-A.12          Il17ra      B-A   -3.537347 -1.822668  -8.163613 3.564467e-06
#> B-A.13           Eltd1      B-A  -12.867730 -3.685686  -7.697566 6.419952e-06
#> B-A.14           Xrcc6      B-A   -3.868517 -1.951781  -7.666373 6.683762e-06
#> B-A.15           Skp1a      B-A    2.899286  1.535698   7.654455 6.787602e-06
#> B-A.16           Csrp1      B-A   -4.209207 -2.073549  -7.641318 6.904062e-06
#> B-A.17            Ighm      B-A   -2.752411 -1.460696  -7.559218 7.682077e-06
#> B-A.18          Samsn1      B-A   -6.104075 -2.609773  -7.376627 9.769386e-06
#> B-A.19          Tagln2      B-A    2.150776  1.104857   7.355488 1.004767e-05
#> B-A.20            Cnn2      B-A    3.060121  1.613589   7.334952 1.032615e-05
#> B-A.21       Serpinb1a      B-A   -4.189829 -2.066891  -7.227236 1.192838e-05
#> B-A.22            Add3      B-A    3.503996  1.809001   6.989236 1.648969e-05
#> B-A.23          Fcer1g      B-A    4.829421  2.271850   6.834513 2.043187e-05
#> B-A.24          Notch1      B-A   -7.561324 -2.918639  -6.757727 2.275132e-05
#> B-A.25          Lgals9      B-A   -3.543556 -1.825198  -6.698467 2.473284e-05
#> B-A.26             Emb      B-A   -2.266181 -1.180263  -6.654169 2.633392e-05
#> B-A.27          Tmem51      B-A  245.180031  7.937698   6.600994 2.840297e-05
#> B-A.28            Grb2      B-A    2.887777  1.529959   6.457725 3.488815e-05
#> B-A.29          Ptp4a3      B-A   -3.061090 -1.614045  -6.452968 3.512883e-05
#> B-A.30        Ankrd13a      B-A   -2.044362 -1.031651  -6.440875 3.574869e-05
#> B-A.31           Runx3      B-A    4.723949  2.239993   6.339389 4.143323e-05
#> B-A.32             Tox      B-A   68.039168  6.088294   6.321678 4.252017e-05
#> B-A.33            Bin2      B-A    3.136470  1.649142   6.162287 5.378258e-05
#> B-A.34           Dusp6      B-A -192.941350 -7.592019  -6.080287 6.077557e-05
#> B-A.35           Coro7      B-A   -2.518113 -1.332343  -6.067403 6.195938e-05
#> B-A.36        Tmem176b      B-A   -2.343345 -1.228569  -5.879357 8.232286e-05
#> B-A.37           Ramp1      B-A   -2.456709 -1.296727  -5.874224 8.296966e-05
#> B-A.38             Myc      B-A   -2.399224 -1.262568  -5.849352 8.618073e-05
#> B-A.39             Cfp      B-A    5.060412  2.339255   5.826732 8.921559e-05
#> B-A.40           Prr13      B-A    3.504592  1.809246   5.819213 9.024934e-05
#> B-A.41           Anxa2      B-A    5.852724  2.549108   5.810872 9.141088e-05
#> B-A.42            Jak1      B-A    2.122836  1.085993   5.696077 1.091110e-04
#> B-A.43           Mef2c      B-A   -3.375360 -1.755041  -5.618635 1.230779e-04
#> B-A.44            Mycn      B-A    6.724322  2.749389   5.468122 1.559196e-04
#> B-A.45           Padi2      B-A    3.140490  1.650990   5.453858 1.594805e-04
#> B-A.46           Cpne2      B-A    6.875618  2.781489   5.412824 1.702116e-04
#> B-A.47            Dbnl      B-A    2.127868  1.089409   5.383098 1.784601e-04
#> B-A.48            Cd81      B-A   -3.810977 -1.930161  -5.323362 1.963385e-04
#> B-A.49          Atp2a3      B-A    2.204316  1.140331   5.234769 2.264137e-04
#> B-A.50           Rogdi      B-A    4.725617  2.240503   5.232174 2.273648e-04
#> B-A.51           Nol4l      B-A  204.580987  7.676528   5.229554 2.283291e-04
#> B-A.52            Cd52      B-A    2.808133  1.489611   5.216068 2.333619e-04
#> B-A.53            Cd34      B-A   -2.639522 -1.400277  -5.200682 2.392467e-04
#> B-A.54            Tcf4      B-A   -3.262787 -1.706105  -5.170808 2.511244e-04
#> B-A.55         Fam189b      B-A    8.540812  3.094373   5.162343 2.546014e-04
#> B-A.56           Esyt1      B-A    2.246173  1.167469   5.063560 2.991430e-04
#> B-A.57          Smim14      B-A   -2.550455 -1.350755  -5.049824 3.059580e-04
#> B-A.58           Egfl7      B-A  -18.634380 -4.219895  -5.023885 3.192769e-04
#> B-A.59           Mgat1      B-A   -2.140223 -1.097761  -4.952621 3.591015e-04
#> B-A.60           Ap3s1      B-A    3.447118  1.785391   4.910781 3.848858e-04
#> B-A.61            Cd93      B-A   -3.555863 -1.830200  -4.867770 4.134270e-04
#> B-A.62           Emid1      B-A    5.006098  2.323687   4.834254 4.372034e-04
#> B-A.63            Myh9      B-A    2.566767  1.359952   4.813335 4.527673e-04
#> B-A.64            Chn2      B-A    9.901175  3.307600   4.810549 4.548837e-04
#> B-A.65          Ifitm2      B-A   -4.640046 -2.214139  -4.803053 4.606293e-04
#> B-A.66          Tyrobp      B-A    3.749437  1.906674   4.802010 4.614347e-04
#> B-A.67             Xpc      B-A   -2.230507 -1.157372  -4.783767 4.757637e-04
#> B-A.68            Il16      B-A   -2.149369 -1.103913  -4.606592 6.417896e-04
#> B-A.69            Fgl2      B-A   10.274616  3.361013   4.598688 6.504805e-04
#> B-A.70          Fam73b      B-A  108.743767  6.764789   4.589438 6.608086e-04
#> B-A.71            Cdc6      B-A   -2.287589 -1.193828  -4.574200 6.781975e-04
#> B-A.72             Gyg      B-A    2.910230  1.541133   4.568512 6.848111e-04
#> B-A.73           Gfod1      B-A   53.047575  5.729215   4.561962 6.925094e-04
#> B-A.74            Ugcg      B-A    5.274843  2.399128   4.559181 6.958056e-04
#> B-A.75         Clec12a      B-A   -7.646484 -2.934797  -4.556716 6.987410e-04
#> B-A.76            Bin1      B-A   -2.426882 -1.279104  -4.544213 7.138285e-04
#> B-A.77         Tnfaip1      B-A  -42.411505 -5.406384  -4.542924 7.154036e-04
#> B-A.78            Cd96      B-A    2.544419  1.347336   4.513896 7.518448e-04
#> B-A.79            Ctsc      B-A    2.019511  1.014006   4.512494 7.536536e-04
#> B-A.80         Tcrg-C1      B-A  400.378672  8.645221   4.506327 7.616631e-04
#> B-A.81             Cd7      B-A  109.175408  6.770504   4.494462 7.773215e-04
#> B-A.82          Malat1      B-A    2.452171  1.294060   4.460130 8.245548e-04
#> B-A.83           Gata3      B-A  211.409117  7.723894   4.458686 8.266047e-04
#> B-A.84         Fam134b      B-A   -5.845783 -2.547396  -4.442935 8.493274e-04
#> B-A.85            Rora      B-A  239.421600  7.903410   4.437781 8.569026e-04
#> B-A.86           Oasl2      B-A  -63.269725 -5.983443  -4.424749 8.763734e-04
#> B-A.87            Fut7      B-A    2.796718  1.483735   4.412888 8.944949e-04
#> B-A.88            Lsm3      B-A   -2.333481 -1.222484  -4.387701 9.342848e-04
#> B-A.89         Plekhf1      B-A   72.941043  6.188659   4.379241 9.480610e-04
#> B-A.90           Rab44      B-A   12.807668  3.678936   4.374293 9.562148e-04
#> B-A.91            Rfc4      B-A   -2.886617 -1.529380  -4.358022 9.835507e-04
#> B-A.92           Cdk19      B-A   -3.422914 -1.775225  -4.308755 1.071390e-03
#> B-A.93         Tmem109      B-A   -3.864945 -1.950448  -4.268835 1.148537e-03
#> B-A.94           Itgb2      B-A    3.170260  1.664601   4.259454 1.167491e-03
#> B-A.95            Rbpj      B-A    2.803813  1.487390   4.249652 1.187643e-03
#> B-A.96           Arl11      B-A  -22.842209 -4.513630  -4.247925 1.191231e-03
#> B-A.97            Dtx4      B-A   -3.539979 -1.823741  -4.224743 1.240496e-03
#> B-A.98         Arhgef3      B-A  -46.655855 -5.543986  -4.204707 1.284783e-03
#> B-A.99            Nrgn      B-A    6.620665  2.726976   4.197008 1.302233e-03
#> B-A.100          Runx2      B-A   -6.249573 -2.643758  -4.182444 1.335921e-03
#> B-A.101          Prrc1      B-A    3.389171  1.760933   4.178637 1.344874e-03
#> B-A.102           Tab2      B-A   27.545118  4.783725   4.173710 1.356555e-03
#> B-A.103           Lmo4      B-A    5.959620  2.575220   4.148608 1.417725e-03
#> B-A.104          Thtpa      B-A   36.910559  5.205962   4.124381 1.479476e-03
#> B-A.105           Cd74      B-A    5.155806  2.366198   4.123526 1.481707e-03
#> B-A.106          Fnbp1      B-A    2.680391  1.422444   4.113368 1.508464e-03
#> B-A.107         Acot11      B-A   71.208895  6.153986   4.110528 1.516033e-03
#> B-A.108       Ccdc102a      B-A   63.919679  5.998188   4.108067 1.522624e-03
#> B-A.109          P2rx7      B-A   10.594439  3.405235   4.076516 1.609814e-03
#> B-A.110        Pip4k2a      B-A   -2.027260 -1.019531  -4.070607 1.626711e-03
#> B-A.111            Npl      B-A  -15.866297 -3.987894  -4.068490 1.632809e-03
#> B-A.112          Smad3      B-A    4.899903  2.292753   4.060342 1.656500e-03
#> B-A.113          Tgtp2      B-A   -3.400936 -1.765932  -4.058037 1.663267e-03
#> B-A.114         Ifitm3      B-A  -32.172476 -5.007755  -4.054787 1.672857e-03
#> B-A.115          Nfil3      B-A    7.200048  2.848007   4.032045 1.741582e-03
#> B-A.116        Zc3hav1      B-A   -2.483945 -1.312633  -4.017993 1.785498e-03
#> B-A.117          Mier3      B-A   -3.715714 -1.893640  -4.017395 1.787391e-03
#> B-A.118           Nkg7      B-A    3.023418  1.596181   4.007159 1.820142e-03
#> B-A.119          Dirc2      B-A  -10.328614 -3.368575  -4.001744 1.837716e-03
#> B-A.120             Hp      B-A   71.955312  6.169029   3.980435 1.908596e-03
#> B-A.121         Sh3bp5      B-A  -27.949023 -4.804726  -3.971700 1.938460e-03
#> B-A.122         Man2a1      B-A   11.026183  3.462862   3.970527 1.942507e-03
#> B-A.123        Ppp2r5a      B-A   -2.037206 -1.026592  -3.962646 1.969929e-03
#> B-A.124          Ero1l      B-A    2.711130  1.438894   3.946863 2.026049e-03
#> B-A.125          Ttc19      B-A   46.185827  5.529378   3.942012 2.043630e-03
#> B-A.126            Gsn      B-A    3.463276  1.792137   3.936011 2.065595e-03
#> B-A.127           Ctso      B-A   -3.810401 -1.929943  -3.934285 2.071959e-03
#> B-A.128           Chka      B-A   -9.035705 -3.175637  -3.934083 2.072703e-03
#> B-A.129          Il2rb      B-A 1209.478551 10.240169   3.921991 2.117868e-03
#> B-A.130         Mif4gd      B-A    5.010323  2.324904   3.916393 2.139119e-03
#> B-A.131           Gse1      B-A   22.218071  4.473662   3.903776 2.187825e-03
#> B-A.132           Snx9      B-A   -4.340030 -2.117705  -3.903573 2.188616e-03
#> B-A.133          Lims1      B-A   -2.160844 -1.111595  -3.902919 2.191174e-03
#> B-A.134           Gm2a      B-A   -2.445361 -1.290048  -3.886258 2.257358e-03
#> B-A.135          Acap2      B-A   -6.967154 -2.800569  -3.878450 2.289078e-03
#> B-A.136         Ndufa4      B-A   -2.214824 -1.147192  -3.853797 2.392286e-03
#> B-A.137          Hvcn1      B-A    3.803689  1.927399   3.836284 2.468512e-03
#> B-A.138           Tle3      B-A    2.638008  1.399449   3.835543 2.471790e-03
#> B-A.139          Aldh2      B-A   -2.246148 -1.167453  -3.816881 2.555888e-03
#> B-A.140           Vav3      B-A   -3.961958 -1.986214  -3.816336 2.558386e-03
#> B-A.141         Shcbp1      B-A   -3.023834 -1.596379  -3.812845 2.574461e-03
#> B-A.142          Egln1      B-A   -2.610038 -1.384071  -3.812177 2.577545e-03
#> B-A.143          Dcaf7      B-A    2.291884  1.196534   3.786309 2.700077e-03
#> B-A.144          Hif1a      B-A    2.020754  1.014894   3.786300 2.700118e-03
#> B-A.145        Aldh4a1      B-A   43.187410  5.432539   3.785835 2.702375e-03
#> B-A.146          Ckap5      B-A   -2.950729 -1.561071  -3.774641 2.757292e-03
#> B-A.147          Prdx4      B-A  -14.734362 -3.881113  -3.767972 2.790551e-03
#> B-A.148           Klf3      B-A  -12.372451 -3.629059  -3.743241 2.917535e-03
#> B-A.149           Apoe      B-A   61.806697  5.949691   3.734907 2.961655e-03
#> B-A.150          Susd1      B-A   -2.404836 -1.265939  -3.723314 3.024167e-03
#> B-A.151          Mpeg1      B-A   -3.673017 -1.876965  -3.723101 3.025330e-03
#> B-A.152           Sdc1      B-A -255.605993 -7.997778  -3.707830 3.109782e-03
#> B-A.153           Cttn      B-A    5.220442  2.384172   3.702531 3.139652e-03
#> B-A.154         Il10rb      B-A   -6.948561 -2.796714  -3.693182 3.193067e-03
#> B-A.155          Lmnb1      B-A   -2.125653 -1.087906  -3.682351 3.256117e-03
#> B-A.156           Gapt      B-A    2.588379  1.372049   3.676965 3.287944e-03
#> B-A.157          Phtf2      B-A    3.073799  1.620023   3.675274 3.298007e-03
#> B-A.158          Rnf14      B-A    3.213478  1.684135   3.674749 3.301134e-03
#> B-A.159          Sesn2      B-A   25.461823  4.670264   3.668941 3.335959e-03
#> B-A.160           Ubl4      B-A   -2.722295 -1.444823  -3.658724 3.398131e-03
#> B-A.161           St14      B-A  -13.594249 -3.764925  -3.657770 3.403998e-03
#> B-A.162         Cyb5r3      B-A    2.070435  1.049934   3.655274 3.419392e-03
#> B-A.163          Rbpms      B-A    3.136838  1.649311   3.644689 3.485485e-03
#> B-A.164           Exo1      B-A   -2.771384 -1.470606  -3.641018 3.508709e-03
#> B-A.165       Trp53i11      B-A   -5.829045 -2.543259  -3.638528 3.524553e-03
#> B-A.166          Ikzf2      B-A   53.533609  5.742373   3.636408 3.538100e-03
#> B-A.167            Cpq      B-A  -18.039568 -4.173093  -3.636202 3.539415e-03
#> B-A.168         Camk1d      B-A   -2.853254 -1.512608  -3.618624 3.653863e-03
#> B-A.169           Rin3      B-A    2.271398  1.183581   3.618373 3.655523e-03
#> B-A.170          Itgb3      B-A   45.994493  5.523389   3.610649 3.707029e-03
#> B-A.171          Plcg1      B-A   -2.862754 -1.517404  -3.600996 3.772444e-03
#> B-A.172           Cdk1      B-A   -2.482017 -1.311513  -3.596702 3.801921e-03
#> B-A.173        Tsc22d1      B-A   -7.528249 -2.912314  -3.593946 3.820964e-03
#> B-A.174            Ivd      B-A   -2.351548 -1.233611  -3.585524 3.879769e-03
#> B-A.175         Nsmce1      B-A   -2.324371 -1.216840  -3.580745 3.913545e-03
#> B-A.176          Cdyl2      B-A   38.815992  5.278579   3.579227 3.924333e-03
#> B-A.177          Kcnk5      B-A  -14.207329 -3.828563  -3.567373 4.009662e-03
#> B-A.178          Iigp1      B-A  -56.392525 -5.817432  -3.564816 4.028317e-03
#> B-A.179           Lyz2      B-A   -8.932998 -3.159145  -3.558808 4.072496e-03
#> B-A.180        Slc22a3      B-A    7.979566  2.996310   3.558434 4.075262e-03
#> B-A.181          Dapp1      B-A    2.160688  1.111491   3.555394 4.097818e-03
#> B-A.182           Wee1      B-A   -2.885280 -1.528711  -3.554705 4.102948e-03
#> B-A.183        Tsc22d3      B-A   -2.503284 -1.323822  -3.549042 4.145371e-03
#> B-A.184          Tpst2      B-A    2.104833  1.073706   3.545921 4.168939e-03
#> B-A.185         Camk2d      B-A   -3.329911 -1.735484  -3.527461 4.311162e-03
#> B-A.186          Stat3      B-A    2.322369  1.215597   3.519902 4.370826e-03
#> B-A.187         Zfp868      B-A    6.470819  2.693948   3.519659 4.372760e-03
#> B-A.188          Paqr8      B-A   38.409959  5.263409   3.506188 4.481232e-03
#> B-A.189         Ppfia1      B-A    2.442647  1.288445   3.504940 4.491415e-03
#> B-A.190           Cpa3      B-A  581.285503  9.183103   3.495718 4.567426e-03
#> B-A.191         Ms4a4c      B-A   32.754586  5.033625   3.487837 4.633427e-03
#> B-A.192          Cdc20      B-A   -2.793176 -1.481907  -3.486416 4.645436e-03
#> B-A.193        Laptm4a      B-A   -2.072610 -1.051449  -3.484759 4.659468e-03
#> B-A.194           Zeb2      B-A   -3.151028 -1.655823  -3.476810 4.727411e-03
#> B-A.195         Tmem43      B-A   -2.619861 -1.389490  -3.441946 5.037555e-03
#> B-A.196          Ccnb2      B-A   -2.276418 -1.186766  -3.436448 5.088318e-03
#> B-A.197         Cldn25      B-A    2.051266  1.036515   3.435542 5.096739e-03
#> B-A.198          Sepp1      B-A    2.396237  1.260770   3.435356 5.098460e-03
#> B-A.199          Abca7      B-A    2.294962  1.198470   3.435036 5.101445e-03
#> B-A.200           Gphn      B-A   -9.600645 -3.263131  -3.432452 5.125546e-03
#> B-A.201        Fam212a      B-A  -13.616238 -3.767256  -3.423441 5.210504e-03
#> B-A.202        Herpud1      B-A   -2.729730 -1.448758  -3.419904 5.244246e-03
#> B-A.203          Ddit4      B-A   14.779726  3.885548   3.418940 5.253484e-03
#> B-A.204        Il12rb2      B-A   16.068919  4.006201   3.416562 5.276333e-03
#> B-A.205         Fam69b      B-A   -5.810566 -2.538679  -3.415590 5.285705e-03
#> B-A.206          Snx14      B-A   60.684304  5.923252   3.413403 5.306845e-03
#> B-A.207          H2-Aa      B-A   34.261610  5.098521   3.412497 5.315633e-03
#> B-A.208          Slfn2      B-A   35.775375  5.160895   3.410746 5.332647e-03
#> B-A.209         Trim59      B-A   -3.001289 -1.585582  -3.407114 5.368127e-03
#> B-A.210           Rxra      B-A   27.350115  4.773475   3.404006 5.398680e-03
#> B-A.211           Cish      B-A  -25.482295 -4.671423  -3.391160 5.526832e-03
#> B-A.212         Map4k5      B-A  -10.974561 -3.456091  -3.389872 5.539852e-03
#> B-A.213        Ppp1r18      B-A   -2.174998 -1.121014  -3.389493 5.543689e-03
#> B-A.214         Sema4a      B-A    6.462744  2.692147   3.388481 5.553941e-03
#> B-A.215           Accs      B-A  -22.196295 -4.472247  -3.356660 5.886542e-03
#> B-A.216         Jmjd1c      B-A    2.808512  1.489806   3.355819 5.895600e-03
#> B-A.217         Ms4a4b      B-A   31.559832  4.980018   3.350416 5.954144e-03
#> B-A.218         Ankib1      B-A   23.407110  4.548875   3.349375 5.965493e-03
#> B-A.219           Ctsw      B-A   11.104279  3.473044   3.348368 5.976491e-03
#> B-A.220          Rab37      B-A    2.402005  1.264239   3.338221 6.088458e-03
#> B-A.221           Rgs3      B-A   28.183894  4.816799   3.333487 6.141426e-03
#> B-A.222           Pim2      B-A    7.076854  2.823108   3.333319 6.143310e-03
#> B-A.223       Rab3gap1      B-A   -3.334710 -1.737561  -3.332156 6.156399e-03
#> B-A.224           Ctc1      B-A   -2.367995 -1.243666  -3.321728 6.275020e-03
#> B-A.225          Hmha1      B-A    2.004662  1.003359   3.321464 6.278048e-03
#> B-A.226            Pbk      B-A   -2.706516 -1.436437  -3.317963 6.318411e-03
#> B-A.227           Lgmn      B-A   -7.914001 -2.984407  -3.313731 6.367548e-03
#> B-A.228          Clip1      B-A   -2.470676 -1.304906  -3.303349 6.489753e-03
#> B-A.229  1110059E24Rik      B-A   -2.867505 -1.519796  -3.302591 6.498763e-03
#> B-A.230          Abcb9      B-A   21.572481  4.431120   3.298821 6.543796e-03
#> B-A.231           Uap1      B-A   -4.001821 -2.000657  -3.295582 6.582726e-03
#> B-A.232           Cd63      B-A    4.984620  2.317484   3.283330 6.732136e-03
#> B-A.233        Ppapdc2      B-A    9.461879  3.242127   3.276787 6.813333e-03
#> B-A.234         Gnptab      B-A   -2.174773 -1.120865  -3.271243 6.882912e-03
#> B-A.235           Rfc3      B-A   -3.296641 -1.720997  -3.270029 6.898243e-03
#> B-A.236           Mafk      B-A   -2.049511 -1.035280  -3.266934 6.937479e-03
#> B-A.237          Ttc13      B-A   -2.049463 -1.035246  -3.260504 7.019729e-03
#> B-A.238          Abhd4      B-A   -2.649640 -1.405796  -3.254930 7.091837e-03
#> B-A.239          Pdcd1      B-A  157.624217  7.300345   3.242380 7.256918e-03
#> B-A.240          Nucb2      B-A  -10.997313 -3.459079  -3.241231 7.272222e-03
#> B-A.241         Nusap1      B-A   -2.736207 -1.452178  -3.236583 7.334486e-03
#> B-A.242         Ccp110      B-A  -10.550172 -3.399195  -3.222555 7.525658e-03
#> B-A.243           Abi2      B-A   -4.169837 -2.059991  -3.219765 7.564276e-03
#> B-A.244         B3gnt5      B-A   -3.721114 -1.895735  -3.219269 7.571176e-03
#> B-A.245          Tram1      B-A   -2.072318 -1.051245  -3.212210 7.669868e-03
#> B-A.246         Tspan2      B-A   -2.819381 -1.495379  -3.211469 7.680313e-03
#> B-A.247          Stag1      B-A   -2.075458 -1.053430  -3.208280 7.725387e-03
#> B-A.248  1110034G24Rik      B-A   21.667617  4.437469   3.201815 7.817593e-03
#> B-A.249        Ankrd44      B-A    2.406178  1.266744   3.200294 7.839436e-03
#> B-A.250         Hivep1      B-A   20.686853  4.370642   3.192819 7.947748e-03
#> B-A.251         Rassf2      B-A   -2.606168 -1.381930  -3.186089 8.046542e-03
#> B-A.252         Cdkn1a      B-A    2.394776  1.259891   3.180279 8.132836e-03
#> B-A.253         Pik3cb      B-A    9.849067  3.299987   3.177524 8.174077e-03
#> B-A.254          Snx30      B-A   -2.887846 -1.529994  -3.176265 8.192990e-03
#> B-A.255           Eya2      B-A   53.341462  5.737185   3.173883 8.228912e-03
#> B-A.256         Adssl1      B-A   -2.144670 -1.100756  -3.168008 8.318165e-03
#> B-A.257          Ssbp3      B-A    3.958112  1.984812   3.165093 8.362805e-03
#> B-A.258          Asap1      B-A   -5.366597 -2.424007  -3.164575 8.370777e-03
#> B-A.259           Rinl      B-A    3.084081  1.624841   3.158727 8.461160e-03
#> B-A.260          H2-Ob      B-A   -2.000929 -1.000670  -3.157021 8.487710e-03
#> B-A.261         Kif20b      B-A  -10.913545 -3.448048  -3.147963 8.630105e-03
#> B-A.262           Btg2      B-A    3.277826  1.712739   3.134252 8.850261e-03
#> B-A.263           Rhof      B-A    7.356213  2.878963   3.127501 8.960713e-03
#> B-A.264         Rassf4      B-A    2.365253  1.241995   3.125388 8.995571e-03
#> B-A.265          Il6st      B-A   -6.253893 -2.644754  -3.118901 9.103453e-03
#> B-A.266          Kif23      B-A   -2.711187 -1.438925  -3.114054 9.184893e-03
#> B-A.267           Prc1      B-A   -2.684647 -1.424733  -3.108971 9.271107e-03
#> B-A.268           Capg      B-A   -2.458168 -1.297584  -3.108032 9.287124e-03
#> B-A.269         Map3k5      B-A   -4.554623 -2.187332  -3.106388 9.315218e-03
#> B-A.270          Csf1r      B-A   -3.971865 -1.989817  -3.099652 9.431263e-03
#> B-A.271        Aldh7a1      B-A   -6.101495 -2.609163  -3.091912 9.566410e-03
#> B-A.272         Hspbp1      B-A   -2.096642 -1.068080  -3.091181 9.579272e-03
#> B-A.273           Ccr2      B-A   33.171551  5.051875   3.088621 9.624452e-03
#> B-A.274         Nudt19      B-A   -2.732073 -1.449996  -3.087787 9.639220e-03
#> B-A.275        Aldh1b1      B-A   -3.427682 -1.777233  -3.087432 9.645524e-03
#> B-A.276           Amd1      B-A   -2.123806 -1.086652  -3.085928 9.672227e-03
#> B-A.277        Wbscr27      B-A   22.768376  4.508960   3.085383 9.681925e-03
#> B-A.278         Lmbrd1      B-A  -20.287104 -4.342491  -3.084719 9.693740e-03
#> B-A.279          Ypel5      B-A    2.256966  1.174385   3.082840 9.727281e-03
#> B-A.280        Creb3l2      B-A   44.650409  5.480601   3.082660 9.730501e-03
#> B-A.281          Rgs14      B-A   24.609769  4.621159   3.074905 9.870224e-03
#> B-A.282         Osbpl3      B-A    3.234768  1.693662   3.069622 9.966562e-03
#> B-A.283          Klhl5      B-A   -2.099241 -1.069868  -3.058994 1.016323e-02
#> B-A.284         Sept11      B-A    2.850468  1.511199   3.057908 1.018356e-02
#> B-A.285        Slc29a1      B-A   -2.092017 -1.064895  -3.050277 1.032746e-02
#> B-A.286         Pkmyt1      B-A   -3.671524 -1.876379  -3.048998 1.035177e-02
#> B-A.287            Btk      B-A   -2.342020 -1.227753  -3.048667 1.035807e-02
#> B-A.288          Gna12      B-A    2.023146  1.016601   3.040240 1.051984e-02
#> B-A.289         Klrb1f      B-A    3.422111  1.774887   3.038345 1.055657e-02
#> B-A.290          Stap1      B-A   -5.226817 -2.385933  -3.036860 1.058544e-02
#> B-A.291          Gpsm3      B-A    5.978963  2.579895   3.036398 1.059443e-02
#> B-A.292           Acp2      B-A   -5.961757 -2.575738  -3.028039 1.075856e-02
#> B-A.293          Terf1      B-A   -4.428062 -2.146675  -3.027216 1.077486e-02
#> B-A.294         B3gnt1      B-A  -15.469373 -3.951343  -3.017847 1.096213e-02
#> B-A.295          Ggta1      B-A   -3.251296 -1.701015  -3.015935 1.100074e-02
#> B-A.296           Elf4      B-A    2.087564  1.061820   3.013539 1.104931e-02
#> B-A.297        Slc29a3      B-A   -2.121542 -1.085114  -3.011709 1.108656e-02
#> B-A.298           Lpxn      B-A    5.088668  2.347288   2.989348 1.155205e-02
#> B-A.299          Frmd6      B-A   -6.443002 -2.687733  -2.989287 1.155333e-02
#> B-A.300           Pdpr      B-A   39.417272  5.300756   2.987320 1.159522e-02
#> B-A.301           Glrx      B-A   -2.846656 -1.509268  -2.982528 1.169787e-02
#> B-A.302          Whamm      B-A   15.510292  3.955154   2.981336 1.172354e-02
#> B-A.303           Fut8      B-A    2.605270  1.381433   2.976882 1.181999e-02
#> B-A.304           Aspm      B-A   -8.923643 -3.157633  -2.976308 1.183246e-02
#> B-A.305          Acss2      B-A    7.298838  2.867667   2.970963 1.194937e-02
#> B-A.306          Uvrag      B-A   -2.333120 -1.222260  -2.967975 1.201523e-02
#> B-A.307         Lpcat4      B-A    5.472063  2.452085   2.964412 1.209424e-02
#> B-A.308           Cd48      B-A   -2.322931 -1.215946  -2.948220 1.245985e-02
#> B-A.309          Dgat1      B-A    3.447476  1.785541   2.947747 1.247068e-02
#> B-A.310          Gcnt1      B-A    5.243668  2.390576   2.946691 1.249493e-02
#> B-A.311          Anks1      B-A   -2.458536 -1.297799  -2.946108 1.250833e-02
#> B-A.312         Gnpda1      B-A    3.150336  1.655506   2.941495 1.261492e-02
#> B-A.313         Zfp945      B-A   10.454887  3.386106   2.940539 1.263711e-02
#> B-A.314           Dgkz      B-A   -2.207429 -1.142367  -2.940171 1.264567e-02
#> B-A.315         Atp8b2      B-A   -5.247670 -2.391677  -2.939500 1.266129e-02
#> B-A.316        Tmem194      B-A  -19.681353 -4.298757  -2.938830 1.267690e-02
#> B-A.317        Tmem55a      B-A  -12.004379 -3.585489  -2.936238 1.273747e-02
#> B-A.318        Nckipsd      B-A   21.074980  4.397459   2.931264 1.285455e-02
#> B-A.319        Eif2ak2      B-A   -5.154857 -2.365933  -2.925612 1.298887e-02
#> B-A.320         Osbpl5      B-A    4.997855  2.321309   2.923556 1.303807e-02
#> B-A.321          Trim8      B-A    2.077155  1.054609   2.914557 1.325566e-02
#> B-A.322        Abhd17b      B-A   -2.094179 -1.066385  -2.912277 1.331134e-02
#> B-A.323           Brd8      B-A   -2.065480 -1.046477  -2.908917 1.339385e-02
#> B-A.324         Hmbox1      B-A   12.678956  3.664364   2.903769 1.352126e-02
#> B-A.325           Vasp      B-A    2.771228  1.470525   2.903643 1.352440e-02
#> B-A.326           Mzt1      B-A    2.703497  1.434827   2.900292 1.360800e-02
#> B-A.327         Zfp763      B-A    4.762157  2.251615   2.896385 1.370612e-02
#> B-A.328         Gm1966      B-A   -3.006111 -1.587899  -2.894890 1.374384e-02
#> B-A.329         Chst15      B-A   -2.639746 -1.400399  -2.888181 1.391445e-02
#> B-A.330          N4bp3      B-A   -4.057466 -2.020579  -2.880183 1.412060e-02
#> B-A.331          Gtse1      B-A   -6.245817 -2.642890  -2.878695 1.415929e-02
#> B-A.332            Kmo      B-A    2.058710  1.041741   2.877984 1.417780e-02
#> B-A.333         Scamp3      B-A    2.055899  1.039770   2.876622 1.421334e-02
#> B-A.334         Pik3cg      B-A   -2.267710 -1.181236  -2.866273 1.448638e-02
#> B-A.335          Rnf13      B-A   -2.284322 -1.191766  -2.863336 1.456481e-02
#> B-A.336         Mboat1      B-A   -9.043712 -3.176915  -2.861665 1.460962e-02
#> B-A.337         Unc119      B-A   10.045409  3.328464   2.861561 1.461240e-02
#> B-A.338          Birc5      B-A   -2.176736 -1.122167  -2.859990 1.465466e-02
#> B-A.339           Rrs1      B-A   -2.059324 -1.042171  -2.857267 1.472822e-02
#> B-A.340           Faah      B-A    2.473889  1.306781   2.855384 1.477929e-02
#> B-A.341       Fam160b2      B-A   17.157551  4.100772   2.845047 1.506278e-02
#> B-A.342          Kdm4a      B-A    3.717693  1.894407   2.838655 1.524078e-02
#> B-A.343           Aff4      B-A    3.071668  1.619022   2.831615 1.543926e-02
#> B-A.344          Akip1      B-A   -3.946743 -1.980663  -2.828639 1.552392e-02
#> B-A.345          Atg9a      B-A    5.675415  2.504726   2.821983 1.571495e-02
#> B-A.346          Myo1c      B-A   -2.314516 -1.210711  -2.815652 1.589884e-02
#> B-A.347         Ddx19b      B-A   -5.619448 -2.490428  -2.813453 1.596322e-02
#> B-A.348          Arrb2      B-A    2.138804  1.096804   2.808580 1.610676e-02
#> B-A.349          Dock5      B-A    9.056004  3.178875   2.795719 1.649181e-02
#> B-A.350         Erlin2      B-A   -2.109614 -1.076979  -2.793972 1.654481e-02
#> B-A.351           Mafg      B-A    4.912622  2.296493   2.791997 1.660495e-02
#> B-A.352         Zfp664      B-A    3.795543  1.924306   2.786789 1.676452e-02
#> B-A.353         Fam64a      B-A   -4.858847 -2.280614  -2.786476 1.677416e-02
#> B-A.354         Ift172      B-A  -17.193293 -4.103774  -2.784232 1.684344e-02
#> B-A.355          Pycr1      B-A    4.109284  2.038887   2.783353 1.687066e-02
#> B-A.356          Hipk2      B-A    3.226803  1.690106   2.778944 1.700779e-02
#> B-A.357          Fbxo4      B-A   -8.567616 -3.098894  -2.776711 1.707767e-02
#> B-A.358          Kif15      B-A   -2.581608 -1.368270  -2.775563 1.711369e-02
#> B-A.359         Golga3      B-A    3.126279  1.644446   2.769815 1.729528e-02
#> B-A.360         Gimap5      B-A    5.150082  2.364595   2.767241 1.737718e-02
#> B-A.361        Golph3l      B-A   -2.367160 -1.243157  -2.756212 1.773257e-02
#> B-A.362  4632428N05Rik      B-A    2.427670  1.279572   2.755858 1.774412e-02
#> B-A.363         Incenp      B-A   -2.159690 -1.110824  -2.753966 1.780582e-02
#> B-A.364          Itm2c      B-A    2.256096  1.173828   2.753082 1.783473e-02
#> B-A.365            Tk1      B-A   -2.062246 -1.044216  -2.752834 1.784288e-02
#> B-A.366         Med12l      B-A    9.370915  3.228190   2.750991 1.790333e-02
#> B-A.367           Bcor      B-A    2.857943  1.514977   2.747686 1.801222e-02
#> B-A.368          Cdca2      B-A   -2.288857 -1.194627  -2.741844 1.820637e-02
#> B-A.369            Mn1      B-A   -3.641285 -1.864448  -2.739597 1.828158e-02
#> B-A.370          Cela1      B-A   15.706193  3.973262   2.735318 1.842567e-02
#> B-A.371            Pml      B-A   -5.162682 -2.368121  -2.734616 1.844942e-02
#> B-A.372          Chtf8      B-A   -2.056274 -1.040033  -2.730794 1.857924e-02
#> B-A.373          Sgpl1      B-A   -2.808808 -1.489958  -2.728659 1.865213e-02
#> B-A.374         Trim65      B-A    5.743630  2.521963   2.727398 1.869532e-02
#> B-A.375          Zfpm1      B-A   37.183496  5.216591   2.723703 1.882248e-02
#> B-A.376           Ect2      B-A   -2.280046 -1.189063  -2.720985 1.891653e-02
#> B-A.377          Neil3      B-A   -4.402805 -2.138423  -2.720644 1.892837e-02
#> B-A.378          Mki67      B-A   -2.644293 -1.402882  -2.719163 1.897986e-02
#> B-A.379        Smarca2      B-A   -2.723511 -1.445467  -2.718392 1.900670e-02
#> B-A.380         Sigirr      B-A   -2.389209 -1.256533  -2.717605 1.903416e-02
#> B-A.381          Cers5      B-A   -2.376995 -1.249139  -2.710970 1.926714e-02
#> B-A.382          Furin      B-A    2.097805  1.068881   2.709980 1.930215e-02
#> B-A.383          Sort1      B-A   -2.620611 -1.389903  -2.701080 1.961968e-02
#> B-A.384          Stim1      B-A    2.363805  1.241111   2.697307 1.975585e-02
#> B-A.385         Vps37a      B-A    3.975381  1.991093   2.696041 1.980171e-02
#> B-A.386         Tbxa2r      B-A   -2.454037 -1.295157  -2.692020 1.994817e-02
#> B-A.387          Gins3      B-A   -6.980337 -2.803297  -2.690564 2.000150e-02
#> B-A.388           Als2      B-A   39.175294  5.291872   2.690087 2.001897e-02
#> B-A.389           Tpx2      B-A   -2.206950 -1.142054  -2.680959 2.035663e-02
#> B-A.390           Gpd2      B-A    3.010032  1.589779   2.677372 2.049081e-02
#> B-A.391  D230025D16Rik      B-A    2.078956  1.055859   2.671393 2.071648e-02
#> B-A.392          Pcgf5      B-A   -3.613577 -1.853428  -2.667633 2.085961e-02
#> B-A.393          Ly6c2      B-A    3.465706  1.793149   2.664482 2.098035e-02
#> B-A.394         Knstrn      B-A   -2.397714 -1.261659  -2.663200 2.102967e-02
#> B-A.395        Pla2g15      B-A   -2.711712 -1.439204  -2.661823 2.108274e-02
#> B-A.396        Dennd5a      B-A   -7.057091 -2.819074  -2.661759 2.108522e-02
#> B-A.397           Cptp      B-A   13.409900  3.745227   2.659450 2.117454e-02
#> B-A.398        Rhobtb3      B-A  -30.631357 -4.936937  -2.657234 2.126063e-02
#> B-A.399          Wdr90      B-A   -2.427357 -1.279386  -2.652237 2.145601e-02
#> B-A.400           Hagh      B-A    4.789201  2.259785   2.650723 2.151557e-02
#> B-A.401          Uhrf1      B-A   -2.204915 -1.140723  -2.646688 2.167502e-02
#> B-A.402         Glt8d1      B-A  -11.162238 -3.480554  -2.645626 2.171722e-02
#> B-A.403         Cd2bp2      B-A   -2.846493 -1.509186  -2.644120 2.177715e-02
#> B-A.404          Prkch      B-A    2.242036  1.164809   2.644004 2.178176e-02
#> B-A.405           Lrmp      B-A   -2.234848 -1.160177  -2.638675 2.199521e-02
#> B-A.406       Tmem176a      B-A   -2.183637 -1.126733  -2.635221 2.213465e-02
#> B-A.407        Tcp11l2      B-A    3.262775  1.706100   2.635116 2.213889e-02
#> B-A.408           E4f1      B-A    4.193093  2.068015   2.634290 2.217237e-02
#> B-A.409           B9d2      B-A   -2.453140 -1.294630  -2.628578 2.240526e-02
#> B-A.410            Hk3      B-A    3.442565  1.783484   2.626284 2.249945e-02
#> B-A.411      Gabarapl1      B-A   -2.290881 -1.195903  -2.625290 2.254042e-02
#> B-A.412           Cln8      B-A    4.376930  2.129919   2.621891 2.268094e-02
#> B-A.413         Sapcd2      B-A  -11.124824 -3.475711  -2.621175 2.271065e-02
#> B-A.414       Mettl7a1      B-A   -2.432662 -1.282536  -2.620967 2.271932e-02
#> B-A.415          Rufy1      B-A   -2.360303 -1.238972  -2.620023 2.275854e-02
#> B-A.416           Rtp4      B-A   -7.563542 -2.919062  -2.619534 2.277890e-02
#> B-A.417          Fads2      B-A   -2.408432 -1.268094  -2.618604 2.281768e-02
#> B-A.418          Capn1      B-A    2.380345  1.251171   2.616604 2.290128e-02
#> B-A.419           Pfkp      B-A    2.219211  1.150047   2.611946 2.309713e-02
#> B-A.420           Pck2      B-A   -2.177932 -1.122959  -2.611893 2.309935e-02
#> B-A.421            Fgr      B-A   15.306101  3.936035   2.610225 2.316988e-02
#> B-A.422         Pou2f1      B-A   -5.392403 -2.430928  -2.606504 2.332803e-02
#> B-A.423        Tubgcp4      B-A   10.839601  3.438240   2.603879 2.344017e-02
#> B-A.424        Ehbp1l1      B-A    2.175719  1.121493   2.603456 2.345831e-02
#> B-A.425           Jak2      B-A    2.050501  1.035977   2.603177 2.347027e-02
#> B-A.426         Nmral1      B-A   -2.152882 -1.106269  -2.603131 2.347225e-02
#> B-A.427         Mrps10      B-A   -2.567407 -1.360312  -2.600326 2.359286e-02
#> B-A.428         Frmd4b      B-A  -14.686877 -3.876456  -2.587046 2.417221e-02
#> B-A.429           Nptn      B-A    3.595179  1.846064   2.586099 2.421408e-02
#> B-A.430       Tnfrsf22      B-A   10.618451  3.408501   2.583428 2.433244e-02
#> B-A.431          Zfp64      B-A    3.690855  1.883955   2.583395 2.433394e-02
#> B-A.432          Scrn2      B-A   -4.995239 -2.320554  -2.582335 2.438109e-02
#> B-A.433          Abhd5      B-A   11.417884  3.513223   2.581436 2.442111e-02
#> B-A.434        Gramd1b      B-A    7.264603  2.860884   2.581299 2.442722e-02
#> B-A.435           Ccnf      B-A   -2.171849 -1.118924  -2.579445 2.451006e-02
#> B-A.436         Ckap2l      B-A   -2.177221 -1.122488  -2.576067 2.466169e-02
#> B-A.437         Gm8995      B-A   -2.114830 -1.080542  -2.572435 2.482574e-02
#> B-A.438           Rrad      B-A   14.984391  3.905389   2.571749 2.485681e-02
#> B-A.439         Dcaf10      B-A   -4.213876 -2.075148  -2.571026 2.488966e-02
#> B-A.440  2610008E11Rik      B-A    6.583987  2.718961   2.568665 2.499714e-02
#> B-A.441       Slc25a24      B-A    4.177114  2.062507   2.567149 2.506638e-02
#> B-A.442         Polr3f      B-A    3.644759  1.865823   2.564753 2.517621e-02
#> B-A.443           Net1      B-A   -2.559028 -1.355596  -2.561261 2.533714e-02
#> B-A.444        Zdhhc21      B-A   10.567585  3.401574   2.559018 2.544103e-02
#> B-A.445          Kif2c      B-A   -2.912936 -1.542474  -2.556235 2.557049e-02
#> B-A.446         Pik3ca      B-A   -2.529420 -1.338807  -2.553937 2.567790e-02
#> B-A.447         Fam45a      B-A   -5.772683 -2.529242  -2.550635 2.583299e-02
#> B-A.448          Herc3      B-A    9.557757  3.256672   2.545328 2.608420e-02
#> B-A.449         Dopey1      B-A  -20.800052 -4.378515  -2.542633 2.621262e-02
#> B-A.450          Bahd1      B-A    4.719971  2.238778   2.535290 2.656582e-02
#> B-A.451         Man2a2      B-A   -2.137078 -1.095640  -2.535176 2.657133e-02
#> B-A.452         Atp1a1      B-A   -2.042449 -1.030300  -2.528760 2.688376e-02
#> B-A.453          Itpr1      B-A    2.312926  1.209719   2.521543 2.723951e-02
#> B-A.454          Ube2t      B-A   -4.882336 -2.287571  -2.516627 2.748444e-02
#> B-A.455          Slfn8      B-A   -8.426315 -3.074902  -2.515069 2.756255e-02
#> B-A.456         Pxylp1      B-A   -4.783554 -2.258083  -2.500616 2.829717e-02
#> B-A.457           Cd69      B-A   -4.284528 -2.099136  -2.500079 2.832479e-02
#> B-A.458          Trit1      B-A   -2.056814 -1.040411  -2.497788 2.844311e-02
#> B-A.459         Zdhhc8      B-A    4.325902  2.113001   2.497235 2.847174e-02
#> B-A.460          Bub1b      B-A   -2.072562 -1.051416  -2.495615 2.855577e-02
#> B-A.461         Setdb2      B-A    8.608696  3.105795   2.495229 2.857582e-02
#> B-A.462         Lpcat1      B-A   -4.394495 -2.135698  -2.486941 2.900972e-02
#> B-A.463        Zscan29      B-A    4.747780  2.247253   2.483487 2.919248e-02
#> B-A.464          Irak3      B-A   14.549212  3.862869   2.483161 2.920978e-02
#> B-A.465           Etv3      B-A    3.218620  1.686442   2.480722 2.933956e-02
#> B-A.466          Kdm7a      B-A    4.642765  2.214984   2.478607 2.945257e-02
#> B-A.467         Ahctf1      B-A    2.025441  1.018236   2.478232 2.947261e-02
#> B-A.468            Ttk      B-A   -4.042282 -2.015170  -2.476593 2.956056e-02
#> B-A.469            Smo      B-A   -2.781851 -1.476045  -2.475777 2.960439e-02
#> B-A.470          Anxa4      B-A   -3.211609 -1.683296  -2.471740 2.982231e-02
#> B-A.471         Zfp110      B-A    3.127126  1.644837   2.470764 2.987525e-02
#> B-A.472           Cnr2      B-A   26.947312  4.752069   2.468995 2.997141e-02
#> B-A.473          P2rx4      B-A    2.731495  1.449691   2.468807 2.998166e-02
#> B-A.474         Dusp11      B-A    3.093569  1.629272   2.467714 3.004124e-02
#> B-A.475         Maged1      B-A   -8.650551 -3.112792  -2.462637 3.031945e-02
#> B-A.476          Pskh1      B-A    2.030532  1.021858   2.461592 3.037705e-02
#> B-A.477         Rhbdf2      B-A    9.189789  3.200032   2.460207 3.045355e-02
#> B-A.478        Jakmip1      B-A   -4.285148 -2.099345  -2.458813 3.053070e-02
#> B-A.479       Arhgap26      B-A    2.962240  1.566689   2.458363 3.055565e-02
#> B-A.480          Pdcd7      B-A   -2.544972 -1.347650  -2.456283 3.067123e-02
#> B-A.481  5430427O19Rik      B-A   -2.223277 -1.152688  -2.455968 3.068877e-02
#> B-A.482         Kif13a      B-A  -12.183541 -3.606862  -2.453143 3.084652e-02
#> B-A.483          Tmem9      B-A   -2.632144 -1.396238  -2.450973 3.096824e-02
#> B-A.484         Zdhhc7      B-A    2.156039  1.108383   2.450072 3.101892e-02
#> B-A.485          Fnip1      B-A   15.713597  3.973942   2.446462 3.122276e-02
#> B-A.486         Supt16      B-A    2.800952  1.485917   2.437948 3.170859e-02
#> B-A.487           Tle6      B-A   14.501060  3.858086   2.437568 3.173049e-02
#> B-A.488           Irf5      B-A    2.701717  1.433877   2.434482 3.190850e-02
#> B-A.489         Zfp597      B-A   15.783303  3.980327   2.433957 3.193890e-02
#> B-A.490        Dennd2c      B-A   16.316431  4.028254   2.432279 3.203620e-02
#> B-A.491         Spata6      B-A   -7.761557 -2.956346  -2.424661 3.248153e-02
#> B-A.492        Arfgap3      B-A   -5.390538 -2.430429  -2.417764 3.288987e-02
#> B-A.493           G2e3      B-A   -4.711753 -2.236264  -2.415623 3.301767e-02
#> B-A.494         H2-Eb1      B-A   34.920172  5.125989   2.414570 3.308069e-02
#> B-A.495         Parp10      B-A   -3.143561 -1.652400  -2.414459 3.308731e-02
#> B-A.496          Atad5      B-A   -2.236697 -1.161370  -2.414071 3.311058e-02
#> B-A.497         Zfp808      B-A    8.718223  3.124034   2.412955 3.317756e-02
#> B-A.498        Fam210a      B-A    2.105309  1.074032   2.409478 3.338708e-02
#> B-A.499        Spata13      B-A   -3.516824 -1.814273  -2.409205 3.340355e-02
#> B-A.500          Nphp1      B-A   -3.844309 -1.942724  -2.408550 3.344318e-02
#> B-A.501         Glipr2      B-A   -4.470185 -2.160335  -2.405112 3.365194e-02
#> B-A.502  A830080D01Rik      B-A   -4.988523 -2.318613  -2.402394 3.381786e-02
#> B-A.503           Ern1      B-A   -3.942238 -1.979015  -2.399875 3.397229e-02
#> B-A.504         Ahcyl2      B-A    3.775768  1.916770   2.396395 3.418686e-02
#> B-A.505         Zfp608      B-A   -4.362595 -2.125187  -2.387918 3.471491e-02
#> B-A.506          Egln3      B-A    3.383397  1.758473   2.384229 3.494717e-02
#> B-A.507         Qtrtd1      B-A  -11.435481 -3.515445  -2.383558 3.498955e-02
#> B-A.508         Apol7e      B-A   17.357080  4.117452   2.383442 3.499692e-02
#> B-A.509          Gstm1      B-A   -2.817776 -1.494557  -2.376851 3.541612e-02
#> B-A.510          Brip1      B-A   -4.374652 -2.129168  -2.376292 3.545191e-02
#> B-A.511          Ppil4      B-A    2.167837  1.116256   2.372756 3.567908e-02
#> B-A.512           Mavs      B-A    2.226461  1.154753   2.368236 3.597147e-02
#> B-A.513          Vars2      B-A    2.408390  1.268069   2.366840 3.606220e-02
#> B-A.514          Prkd2      B-A   -6.263510 -2.646971  -2.365998 3.611709e-02
#> B-A.515         Sorbs3      B-A   -4.093843 -2.033456  -2.364090 3.624169e-02
#> B-A.516        Tsc22d2      B-A    3.490739  1.803532   2.362293 3.635940e-02
#> B-A.517          Snx11      B-A   -4.123106 -2.043732  -2.361165 3.643351e-02
#> B-A.518         Ptpn18      B-A   -2.516984 -1.331696  -2.361049 3.644110e-02
#> B-A.519          Rcor2      B-A    3.812438  1.930714   2.360868 3.645303e-02
#> B-A.520           Sos2      B-A    2.683699  1.424223   2.360743 3.646125e-02
#> B-A.521          Ndor1      B-A    2.365236  1.241984   2.359770 3.652532e-02
#> B-A.522        Slc27a1      B-A   -3.996565 -1.998760  -2.359133 3.656734e-02
#> B-A.523          Birc3      B-A  -11.620751 -3.538631  -2.355553 3.680423e-02
#> B-A.524           Ice2      B-A   -3.275053 -1.711518  -2.354399 3.688097e-02
#> B-A.525        Mettl14      B-A   -3.396491 -1.764045  -2.350766 3.712337e-02
#> B-A.526          Limk1      B-A    2.135015  1.094246   2.346129 3.743505e-02
#> B-A.527        Gm16039      B-A    4.496636  2.168846   2.341110 3.777523e-02
#> B-A.528           Pcnt      B-A   -2.269700 -1.182501  -2.337775 3.800293e-02
#> B-A.529          Brca2      B-A   -3.588345 -1.843319  -2.336798 3.806982e-02
#> B-A.530          Mdfic      B-A    8.781579  3.134480   2.334428 3.823273e-02
#> B-A.531        Tmem260      B-A   -2.050533 -1.035999  -2.333901 3.826904e-02
#> B-A.532          Kif11      B-A   -2.101100 -1.071145  -2.332412 3.837181e-02
#> B-A.533         Zfp109      B-A    6.629889  2.728985   2.328807 3.862170e-02
#> B-A.534         Rnf141      B-A   -2.535999 -1.342554  -2.328322 3.865544e-02
#> B-A.535            Nin      B-A    3.090534  1.627856   2.328026 3.867601e-02
#> B-A.536          Cyth3      B-A   15.780102  3.980035   2.326396 3.878971e-02
#> B-A.537           Sc5d      B-A    2.577061  1.365727   2.322791 3.904221e-02
#> B-A.538         Shkbp1      B-A    2.482620  1.311864   2.322290 3.907741e-02
#> B-A.539         Pik3r5      B-A    5.191969  2.376282   2.322128 3.908880e-02
#> B-A.540          Dtwd1      B-A   -4.342438 -2.118505  -2.320596 3.919674e-02
#> B-A.541           Rdm1      B-A   -6.418870 -2.682319  -2.319319 3.928691e-02
#> B-A.542        Ankrd32      B-A   -6.345679 -2.665775  -2.317625 3.940683e-02
#> B-A.543          Ckap2      B-A   -5.345770 -2.418398  -2.316849 3.946183e-02
#> B-A.544         Sptlc1      B-A    2.303694  1.203949   2.316416 3.949260e-02
#> B-A.545         Zmynd8      B-A   -2.048360 -1.034469  -2.314509 3.962832e-02
#> B-A.546         Ifnar1      B-A   -2.171501 -1.118693  -2.313426 3.970554e-02
#> B-A.547          Actr6      B-A   -3.417725 -1.773036  -2.306100 4.023202e-02
#> B-A.548          Fbxw4      B-A    4.025034  2.009001   2.305878 4.024801e-02
#> B-A.549          Mfsd4      B-A    2.265514  1.179838   2.304966 4.031405e-02
#> B-A.550          Prtn3      B-A    2.073758  1.052248   2.298839 4.076036e-02
#> B-A.551          Tiam1      B-A   12.041270  3.589916   2.294491 4.107992e-02
#> B-A.552           Sik1      B-A    4.094364  2.033640   2.293156 4.117854e-02
#> B-A.553         Klhl36      B-A    3.281136  1.714196   2.292000 4.126407e-02
#> B-A.554            Evl      B-A   -5.678417 -2.505489  -2.290805 4.135269e-02
#> B-A.555          Dhx29      B-A    2.751238  1.460081   2.290112 4.140420e-02
#> B-A.556        Racgap1      B-A   -2.219420 -1.150183  -2.289972 4.141459e-02
#> B-A.557        Elmsan1      B-A    6.329226  2.662029   2.289286 4.146562e-02
#> B-A.558          Dscr3      B-A   -2.706768 -1.436571  -2.288384 4.153278e-02
#> B-A.559          Bpnt1      B-A    3.331966  1.736374   2.284641 4.181274e-02
#> B-A.560         Mfsd11      B-A   -5.319529 -2.411299  -2.278772 4.225532e-02
#> B-A.561          Tdrd7      B-A  -14.105988 -3.818236  -2.276017 4.246465e-02
#> B-A.562        Rtn4rl1      B-A    2.096582  1.068039   2.275427 4.250955e-02
#> B-A.563           Pogk      B-A    6.713672  2.747102   2.271656 4.279790e-02
#> B-A.564     St6galnac6      B-A    2.915356  1.543672   2.265229 4.329370e-02
#> B-A.565          Pde6d      B-A   -7.443555 -2.895992  -2.263989 4.339003e-02
#> B-A.566           Mcm8      B-A  -10.420727 -3.381384  -2.262899 4.347482e-02
#> B-A.567            Ada      B-A   -3.596557 -1.846617  -2.260517 4.366070e-02
#> B-A.568           Pcnx      B-A   -3.918945 -1.970465  -2.258626 4.380884e-02
#> B-A.569            Pcx      B-A   15.130499  3.919388   2.258545 4.381516e-02
#> B-A.570           Ass1      B-A   -4.055432 -2.019856  -2.256588 4.396893e-02
#> B-A.571       Dync1li2      B-A   -2.195876 -1.134796  -2.256196 4.399982e-02
#> B-A.572         Snap29      B-A    2.669516  1.416578   2.255491 4.405541e-02
#> B-A.573           Dctd      B-A   -2.467847 -1.303253  -2.255302 4.407029e-02
#> B-A.574          Gnptg      B-A    3.153229  1.656830   2.254634 4.412300e-02
#> B-A.575         Acadvl      B-A    2.401319  1.263827   2.253199 4.423655e-02
#> B-A.576         Kif21b      B-A   10.190176  3.349107   2.253031 4.424982e-02
#> B-A.577          Stk10      B-A   -2.148610 -1.103404  -2.250175 4.447660e-02
#> B-A.578        Zkscan3      B-A    2.372289  1.246280   2.249713 4.451335e-02
#> B-A.579           Gbp2      B-A   -8.755519 -3.130193  -2.249517 4.452892e-02
#> B-A.580          Gins2      B-A   -2.314340 -1.210601  -2.245182 4.487563e-02
#> B-A.581          Abtb1      B-A    2.259918  1.176271   2.241077 4.520631e-02
#> B-A.582         B3gat3      B-A    3.845457  1.943155   2.238881 4.538412e-02
#> B-A.583         Map3k4      B-A   -2.280066 -1.189076  -2.236366 4.558858e-02
#> B-A.584           Plk4      B-A   -2.270757 -1.183174  -2.234294 4.575770e-02
#> B-A.585           Tet1      B-A   -2.743363 -1.455946  -2.232024 4.594371e-02
#> B-A.586       Mapkapk3      B-A    2.089427  1.063107   2.230445 4.607348e-02
#> B-A.587  D030056L22Rik      B-A   -3.285809 -1.716249  -2.229291 4.616854e-02
#> B-A.588         Sap30l      B-A    5.087269  2.346891   2.221685 4.679985e-02
#> B-A.589           Lnx2      B-A    4.430157  2.147358   2.221042 4.685365e-02
#> B-A.590         Dhrs11      B-A    3.435622  1.780571   2.218664 4.705290e-02
#> B-A.591           Pkig      B-A   -2.404827 -1.265933  -2.215805 4.729355e-02
#> B-A.592        Rtn4ip1      B-A    2.882524  1.527333   2.214062 4.744079e-02
#> B-A.593        Proser1      B-A    4.036045  2.012942   2.212825 4.754561e-02
#> B-A.594          Pde4b      B-A   -2.370156 -1.244982  -2.208155 4.794318e-02
#> B-A.595         Nufip2      B-A    2.464402  1.301238   2.207612 4.798969e-02
#> B-A.596          Pear1      B-A   -2.243324 -1.165638  -2.205971 4.813024e-02
#> B-A.597        Il18rap      B-A   -3.387196 -1.760091  -2.204037 4.829641e-02
#> B-A.598          Cbll1      B-A    3.847999  1.944108   2.203677 4.832744e-02
#> B-A.599          Pus7l      B-A    3.596367  1.846540   2.202055 4.846734e-02
#> B-A.600          Mrps6      B-A   -4.223550 -2.078456  -2.200031 4.864246e-02
#> B-A.601         Zfp830      B-A    2.417943  1.273780   2.199370 4.869973e-02
#> B-A.602         Marcks      B-A   -3.045408 -1.606635  -2.199177 4.871647e-02
#> B-A.603          Nsun3      B-A   12.522365  3.646435   2.198568 4.876940e-02
#> B-A.604          Zzef1      B-A   14.207578  3.828589   2.198434 4.878102e-02
#> B-A.605        Fastkd3      B-A   -2.493383 -1.318105  -2.197557 4.885730e-02
#> B-A.606        St3gal5      B-A    7.045211  2.816643   2.194663 4.910971e-02
#> B-A.607         Klhl42      B-A   -2.642711 -1.402019  -2.193438 4.921700e-02
#> B-A.608         Ddx26b      B-A    3.231596  1.692247   2.191509 4.938629e-02
#> B-A.609          Rad51      B-A   -2.081532 -1.057645  -2.187600 4.973107e-02
#> B-A.610  1110057K04Rik      B-A   -2.661698 -1.412347  -2.186139 4.986048e-02
#> B-A.611         Zfand1      B-A   -3.634121 -1.861607  -2.185397 4.992641e-02
#> B-A.612          Siah2      B-A   -5.169170 -2.369933  -2.184969 4.996441e-02
#> C-A.1            H2afy      C-A   -5.347700 -2.418918 -16.225491 2.163237e-09
#> C-A.2           Tmsb4x      C-A    4.079607  2.028430  14.396435 8.258691e-09
#> C-A.3            Itgb7      C-A    8.758115  3.130620  11.236885 1.254682e-07
#> C-A.4            Tapt1      C-A   -5.281643 -2.400987 -11.155166 1.357504e-07
#> C-A.5             Sell      C-A   -6.133608 -2.616736 -11.149835 1.364521e-07
#> C-A.6             Lmo4      C-A   87.188286  6.446062  10.980139 1.609614e-07
#> C-A.7            Il2rg      C-A    7.946997  2.990410  10.735937 2.049305e-07
#> C-A.8             Irf8      C-A  -15.186704 -3.924737 -10.680906 2.165288e-07
#> C-A.9             Ighm      C-A   -4.331836 -2.114979 -10.574966 2.408973e-07
#> C-A.10           Esyt1      C-A    4.576928  2.194380  10.550526 2.469287e-07
#> C-A.11           Gpr97      C-A   -8.565627 -3.098559 -10.532731 2.514223e-07
#> C-A.12            Cd82      C-A   12.290158  3.619432  10.426215 2.802303e-07
#> C-A.13             Id2      C-A  168.950218  7.400454  10.339187 3.064086e-07
#> C-A.14           Plac8      C-A  -15.583847 -3.961980 -10.288433 3.228823e-07
#> C-A.15           Prr13      C-A    7.306977  2.869275  10.175043 3.632347e-07
#> C-A.16           Stmn1      C-A   -2.970409 -1.570662 -10.131261 3.802408e-07
#> C-A.17             Myc      C-A   -5.071972 -2.342547 -10.002701 4.353246e-07
#> C-A.18            Cnn3      C-A  -12.796789 -3.677710  -9.938727 4.658845e-07
#> C-A.19          Ifngr1      C-A    3.827268  1.936315   9.774838 5.552173e-07
#> C-A.20            Ctr9      C-A   -2.884560 -1.528351  -9.746870 5.722241e-07
#> C-A.21           Mgat1      C-A   -5.590455 -2.482966  -9.696181 6.044921e-07
#> C-A.22           Rftn1      C-A    3.263421  1.706385   9.688406 6.096122e-07
#> C-A.23           Ero1l      C-A    8.571943  3.099622   9.565295 6.972198e-07
#> C-A.24            Jak1      C-A    3.134028  1.648018   9.521010 7.319680e-07
#> C-A.25            Map4      C-A    6.391533  2.676162   9.406590 8.306681e-07
#> C-A.26            Mcm6      C-A   -2.156713 -1.108834  -9.141664 1.118563e-06
#> C-A.27         Hnrnpa1      C-A   -2.939855 -1.555745  -9.126098 1.138528e-06
#> C-A.28        Slc25a12      C-A   -6.121088 -2.613788  -9.113054 1.155553e-06
#> C-A.29           Xrcc6      C-A   -5.288881 -2.402962  -9.050019 1.241773e-06
#> C-A.30           Actr3      C-A    2.307154  1.206114   9.046711 1.246484e-06
#> C-A.31          Zfp706      C-A   -2.892236 -1.532185  -8.903816 1.469526e-06
#> C-A.32             Msn      C-A    2.565530  1.359257   8.818238 1.623352e-06
#> C-A.33            Faah      C-A   12.694778  3.666163   8.811300 1.636559e-06
#> C-A.34          Ptp4a3      C-A   -5.478977 -2.453907  -8.599785 2.100085e-06
#> C-A.35          Dnajc3      C-A    3.464381  1.792598   8.525848 2.293864e-06
#> C-A.36        BC035044      C-A  -15.304945 -3.935926  -8.505792 2.349675e-06
#> C-A.37          Ndfip1      C-A    4.422077  2.144724   8.503177 2.357059e-06
#> C-A.38            Cnn2      C-A    3.313874  1.728519   8.431650 2.568986e-06
#> C-A.39            Bcl2      C-A    4.016658  2.005996   8.401388 2.664726e-06
#> C-A.40          Galnt1      C-A    2.822366  1.496905   8.331222 2.901747e-06
#> C-A.41          Shisa5      C-A    3.987397  1.995447   8.195301 3.427702e-06
#> C-A.42          Lgals9      C-A   -4.603995 -2.202886  -8.154617 3.604347e-06
#> C-A.43           Eltd1      C-A   -5.648828 -2.497952  -8.109046 3.813872e-06
#> C-A.44           Fnbp1      C-A    5.586884  2.482044   8.004267 4.346773e-06
#> C-A.45            Cd34      C-A  -48.030892 -5.585891  -7.975112 4.508863e-06
#> C-A.46            Ets1      C-A   18.542484  4.212763   7.953965 4.630473e-06
#> C-A.47           Itgal      C-A    3.799982  1.925993   7.949445 4.656923e-06
#> C-A.48          Tmem51      C-A  506.267203  8.983755   7.906886 4.914019e-06
#> C-A.49           Srsf7      C-A   -3.265755 -1.707416  -7.870247 5.147548e-06
#> C-A.50         Unc93b1      C-A   -9.413065 -3.234665  -7.828335 5.429335e-06
#> C-A.51           Prkcq      C-A   11.242011  3.490828   7.766664 5.874426e-06
#> C-A.52            Gpx1      C-A   -2.885628 -1.528886  -7.670398 6.649075e-06
#> C-A.53          Rnf145      C-A    4.722436  2.239531   7.643000 6.889035e-06
#> C-A.54            Dntt      C-A -356.676314 -8.478472  -7.622449 7.075097e-06
#> C-A.55            Idh2      C-A   -3.327400 -1.734395  -7.610765 7.183265e-06
#> C-A.56          Ahcyl2      C-A   48.434492  5.597963   7.536872 7.909717e-06
#> C-A.57           Mef2c      C-A  -19.703817 -4.300403  -7.512527 8.165960e-06
#> C-A.58           Il6st      C-A  -25.873648 -4.693412  -7.510631 8.186281e-06
#> C-A.59           Gpr56      C-A  -29.945649 -4.904275  -7.482604 8.493133e-06
#> C-A.60           Itm2c      C-A    6.951168  2.797255   7.435315 9.039101e-06
#> C-A.61           Paics      C-A   -2.659383 -1.411092  -7.416937 9.261300e-06
#> C-A.62            Pgls      C-A   -3.488409 -1.802569  -7.374891 9.791921e-06
#> C-A.63           Nrros      C-A   -3.205899 -1.680729  -7.277014 1.115721e-05
#> C-A.64            Ly6e      C-A    2.499572  1.321681   7.228573 1.190693e-05
#> C-A.65          Gpr171      C-A   -3.577725 -1.839042  -7.191166 1.252267e-05
#> C-A.66           Srsf3      C-A   -2.203580 -1.139849  -7.151501 1.321290e-05
#> C-A.67           Satb1      C-A   -7.830846 -2.969168  -7.147230 1.328961e-05
#> C-A.68          Il17rb      C-A   21.530353  4.428300   7.075642 1.464871e-05
#> C-A.69             Dtl      C-A   -4.351911 -2.121649  -7.051616 1.513750e-05
#> C-A.70           Stim1      C-A    7.653413  2.936103   7.022522 1.575281e-05
#> C-A.71          Notch1      C-A   -8.582268 -3.101359  -7.001101 1.622291e-05
#> C-A.72           Furin      C-A    5.399093  2.432717   6.953997 1.731007e-05
#> C-A.73           Itgb2      C-A    5.527340  2.466585   6.931669 1.785232e-05
#> C-A.74          Slain2      C-A    3.846747  1.943639   6.925600 1.800283e-05
#> C-A.75            Lsm3      C-A   -4.151033 -2.053471  -6.876373 1.927477e-05
#> C-A.76            Ctsc      C-A    2.664978  1.414124   6.849167 2.001860e-05
#> C-A.77             Tox      C-A   81.031298  6.340407   6.823742 2.074141e-05
#> C-A.78            Cd47      C-A   -2.078250 -1.055369  -6.800916 2.141407e-05
#> C-A.79           Hspa9      C-A   -2.138664 -1.096710  -6.763738 2.255999e-05
#> C-A.80            Fli1      C-A    3.316483  1.729654   6.750735 2.297604e-05
#> C-A.81           Aldoa      C-A    2.634664  1.397619   6.744622 2.317448e-05
#> C-A.82            Mcm7      C-A   -2.578116 -1.366317  -6.722387 2.391174e-05
#> C-A.83         Hnrnpab      C-A   -2.165614 -1.114776  -6.722270 2.391568e-05
#> C-A.84           Ramp1      C-A   -2.690344 -1.427791  -6.711373 2.428615e-05
#> C-A.85            Cd74      C-A   11.532442  3.527626   6.674495 2.558597e-05
#> C-A.86           Nup62      C-A   -2.671572 -1.417689  -6.597722 2.853582e-05
#> C-A.87            Btg1      C-A    6.035236  2.593410   6.590655 2.882500e-05
#> C-A.88            Aff4      C-A   10.807301  3.433934   6.551061 3.050404e-05
#> C-A.89           Dcaf7      C-A    3.629819  1.859897   6.493242 3.314557e-05
#> C-A.90           Hif1a      C-A    3.006195  1.587938   6.489582 3.332078e-05
#> C-A.91            Cd93      C-A   -5.749183 -2.523357  -6.466185 3.446440e-05
#> C-A.92           Esyt2      C-A    4.166303  2.058768   6.452029 3.517654e-05
#> C-A.93            Sox4      C-A   -2.961356 -1.566258  -6.437011 3.594920e-05
#> C-A.94          H2-DMa      C-A   -4.719131 -2.238521  -6.410622 3.735098e-05
#> C-A.95            Flna      C-A    2.156537  1.108716   6.370557 3.959149e-05
#> C-A.96           Ap3s1      C-A    4.401861  2.138114   6.340323 4.137672e-05
#> C-A.97            Tcf4      C-A   -4.405061 -2.139162  -6.335951 4.164190e-05
#> C-A.98            Rora      C-A 2234.021189 11.125427   6.325913 4.225756e-05
#> C-A.99           Hmgn2      C-A   -2.080367 -1.056838  -6.304178 4.362401e-05
#> C-A.100         Rasal3      C-A    2.113568  1.079681   6.246873 4.745694e-05
#> C-A.101           Flt3      C-A -347.254654 -8.439850  -6.237546 4.811399e-05
#> C-A.102          Cdca7      C-A   -3.313278 -1.728259  -6.233299 4.841632e-05
#> C-A.103          Csrp1      C-A   -2.785973 -1.478181  -6.206283 5.038762e-05
#> C-A.104          Gpr65      C-A    5.826992  2.542751   6.205062 5.047865e-05
#> C-A.105          Clic4      C-A   -3.658637 -1.871306  -6.189150 5.168183e-05
#> C-A.106         Fam65a      C-A   -3.893755 -1.961162  -6.170534 5.312815e-05
#> C-A.107          Cmtm6      C-A    2.647365  1.404557   6.165049 5.356248e-05
#> C-A.108           Lcp2      C-A    4.600274  2.201720   6.156116 5.427786e-05
#> C-A.109           Ncf1      C-A  -12.262691 -3.616204  -6.120042 5.727185e-05
#> C-A.110           Gclc      C-A    8.390543  3.068764   6.109517 5.817811e-05
#> C-A.111           Lpxn      C-A   21.881897  4.451666   6.106573 5.843430e-05
#> C-A.112          Dnmt1      C-A   -2.219127 -1.149992  -6.093313 5.960314e-05
#> C-A.113           Hes6      C-A   -3.063274 -1.615074  -6.062346 6.243069e-05
#> C-A.114           Rbl2      C-A    4.386948  2.133218   6.062223 6.244218e-05
#> C-A.115        Plekhb2      C-A    2.689334  1.427249   5.992168 6.937888e-05
#> C-A.116         Slc9a9      C-A    4.082867  2.029583   5.979628 7.070457e-05
#> C-A.117       Marcksl1      C-A   -3.089540 -1.627392  -5.961939 7.262059e-05
#> C-A.118          Itgb3      C-A  504.901400  8.979858   5.959900 7.284490e-05
#> C-A.119           Hn1l      C-A   -2.737934 -1.453088  -5.939660 7.511233e-05
#> C-A.120           Add3      C-A    2.739239  1.453775   5.937474 7.536172e-05
#> C-A.121        Fam102a      C-A    5.740591  2.521199   5.933714 7.579265e-05
#> C-A.122           Akna      C-A    2.497776  1.320644   5.919230 7.747716e-05
#> C-A.123        Smpdl3a      C-A    4.469794  2.160208   5.908410 7.876144e-05
#> C-A.124          Lpar6      C-A    7.536813  2.913955   5.895639 8.030639e-05
#> C-A.125        Ehbp1l1      C-A    4.902356  2.293475   5.891324 8.083570e-05
#> C-A.126         Il27ra      C-A    8.769486  3.132492   5.838548 8.761640e-05
#> C-A.127         Anapc5      C-A   -2.071155 -1.050435  -5.826920 8.918996e-05
#> C-A.128         Ifitm2      C-A  -21.870454 -4.450911  -5.819363 9.022858e-05
#> C-A.129          Gata3      C-A  892.769845  9.802144   5.803666 9.242727e-05
#> C-A.130        Plekhf1      C-A  188.952950  7.561883   5.774185 9.671140e-05
#> C-A.131         Il17ra      C-A   -2.187170 -1.129066  -5.762785 9.842402e-05
#> C-A.132          Nup85      C-A   -2.766646 -1.468138  -5.739844 1.019684e-04
#> C-A.133         Gnptab      C-A   -4.855378 -2.279584  -5.723658 1.045504e-04
#> C-A.134           Cdc6      C-A   -2.818517 -1.494936  -5.704151 1.077546e-04
#> C-A.135      Epb4.1l4b      C-A   -4.093668 -2.033394  -5.701108 1.082636e-04
#> C-A.136          Cela1      C-A  113.313563  6.824177   5.698634 1.086794e-04
#> C-A.137        Tcrg-C1      C-A 1565.558834 10.612462   5.690734 1.100184e-04
#> C-A.138          Nol4l      C-A  219.703260  7.779412   5.686738 1.107025e-04
#> C-A.139          Ints1      C-A   -3.249303 -1.700130  -5.685902 1.108462e-04
#> C-A.140            Sla      C-A    2.741255  1.454837   5.675427 1.126630e-04
#> C-A.141          Myo1g      C-A    2.390063  1.257049   5.665967 1.143309e-04
#> C-A.142           Cerk      C-A    4.055530  2.019890   5.634134 1.201381e-04
#> C-A.143         Stk17b      C-A    3.127046  1.644800   5.626692 1.215403e-04
#> C-A.144       Rap1gds1      C-A    4.216251  2.075961   5.610401 1.246704e-04
#> C-A.145           Ybx3      C-A   -2.716646 -1.441827  -5.602431 1.262328e-04
#> C-A.146          Napsa      C-A   -2.811666 -1.491425  -5.595710 1.275663e-04
#> C-A.147        Fam189b      C-A    8.797173  3.137040   5.584186 1.298877e-04
#> C-A.148        Slc29a3      C-A   -7.765538 -2.957086  -5.577292 1.312977e-04
#> C-A.149          Naa50      C-A   -2.703898 -1.435041  -5.570484 1.327061e-04
#> C-A.150         Pdgfrb      C-A   -6.552977 -2.712150  -5.562356 1.344084e-04
#> C-A.151          Phgdh      C-A   -2.225055 -1.153841  -5.534617 1.403941e-04
#> C-A.152          Idh3g      C-A   -2.240829 -1.164032  -5.522288 1.431444e-04
#> C-A.153         Malat1      C-A    2.794529  1.482605   5.512088 1.454627e-04
#> C-A.154          Ptpn6      C-A   -2.718546 -1.442835  -5.460388 1.578400e-04
#> C-A.155           Ugcg      C-A    6.285771  2.652090   5.457947 1.584511e-04
#> C-A.156         Nfatc3      C-A    2.828407  1.499990   5.453588 1.595488e-04
#> C-A.157          Cops4      C-A   -2.981113 -1.575851  -5.447982 1.609724e-04
#> C-A.158          Egfl7      C-A  -13.153774 -3.717405  -5.428833 1.659369e-04
#> C-A.159          Ikzf2      C-A  287.251928  8.166173   5.421547 1.678683e-04
#> C-A.160          Socs1      C-A   90.043090  6.492544   5.408240 1.714568e-04
#> C-A.161          Stag1      C-A   -4.845591 -2.276673  -5.408128 1.714874e-04
#> C-A.162         Ltb4r1      C-A   61.681479  5.946765   5.397161 1.745065e-04
#> C-A.163           Nkg7      C-A    3.930334  1.974652   5.396633 1.746532e-04
#> C-A.164          Tmed3      C-A   -4.075043 -2.026815  -5.380272 1.792659e-04
#> C-A.165         Ptpn22      C-A    4.519760  2.176246   5.357120 1.860143e-04
#> C-A.166           Tifa      C-A   -4.916578 -2.297654  -5.349809 1.882007e-04
#> C-A.167          H2-Ob      C-A   -3.364184 -1.750256  -5.348048 1.887312e-04
#> C-A.168          Hmga1      C-A   -3.688656 -1.883095  -5.344837 1.897028e-04
#> C-A.169          Sirt3      C-A   -4.627840 -2.210339  -5.322142 1.967230e-04
#> C-A.170           Lat2      C-A   -3.258050 -1.704009  -5.311547 2.000933e-04
#> C-A.171       Tmem229b      C-A   -2.959597 -1.565401  -5.279794 2.105647e-04
#> C-A.172            Srm      C-A   -2.744293 -1.456434  -5.277430 2.113669e-04
#> C-A.173         Lpcat4      C-A   16.363470  4.032407   5.251107 2.205220e-04
#> C-A.174           Sdhd      C-A   -2.181228 -1.125141  -5.247105 2.219502e-04
#> C-A.175  4632428N05Rik      C-A    4.661733  2.220866   5.246935 2.220110e-04
#> C-A.176         Atp2a3      C-A    2.105457  1.074134   5.246215 2.222690e-04
#> C-A.177          Ssbp3      C-A    8.171827  3.030659   5.235380 2.261904e-04
#> C-A.178          Tpd52      C-A   -2.690110 -1.427665  -5.211148 2.352269e-04
#> C-A.179         Adssl1      C-A   -4.092075 -2.032833  -5.205932 2.372213e-04
#> C-A.180          Hmha1      C-A    2.708330  1.437403   5.197278 2.405697e-04
#> C-A.181           Tab2      C-A   48.855244  5.610442   5.194747 2.415584e-04
#> C-A.182          Skap2      C-A   -4.170632 -2.060266  -5.184582 2.455725e-04
#> C-A.183            Dek      C-A   -2.393838 -1.259325  -5.180075 2.473749e-04
#> C-A.184       AU040320      C-A    3.022402  1.595695   5.159511 2.557763e-04
#> C-A.185        S100a10      C-A    3.232653  1.692718   5.157986 2.564110e-04
#> C-A.186          Cnpy3      C-A   -2.936742 -1.554217  -5.141542 2.633633e-04
#> C-A.187           Myh9      C-A    2.564449  1.358649   5.140070 2.639953e-04
#> C-A.188           Cdk4      C-A   -2.630627 -1.395407  -5.138171 2.648131e-04
#> C-A.189           Dok2      C-A    4.925918  2.300393   5.131710 2.676150e-04
#> C-A.190          Runx2      C-A  -10.102123 -3.336587  -5.129959 2.683798e-04
#> C-A.191           Rrm1      C-A   -2.934771 -1.553248  -5.128064 2.692098e-04
#> C-A.192         Smim14      C-A   -2.449332 -1.292388  -5.127309 2.695413e-04
#> C-A.193            Myb      C-A   -4.347423 -2.120160  -5.106272 2.789532e-04
#> C-A.194            Cd7      C-A  172.714766  7.432248   5.103873 2.800486e-04
#> C-A.195        Themis2      C-A  -42.969539 -5.425242  -5.103635 2.801574e-04
#> C-A.196         Fam73b      C-A  124.353937  6.958308   5.094030 2.845889e-04
#> C-A.197            Btk      C-A  -16.042121 -4.003793  -5.088362 2.872387e-04
#> C-A.198          Gfod1      C-A   57.682632  5.850065   5.087205 2.877831e-04
#> C-A.199          Wdfy4      C-A   -6.008352 -2.586969  -5.086571 2.880815e-04
#> C-A.200          Mdfic      C-A   66.360121  6.052245   5.077615 2.923347e-04
#> C-A.201            Bid      C-A   -4.995443 -2.320612  -5.069154 2.964133e-04
#> C-A.202        Plekha2      C-A    2.870247  1.521175   5.058576 3.015970e-04
#> C-A.203         Trip12      C-A    3.053603  1.610512   5.052886 3.044248e-04
#> C-A.204         Fcer1g      C-A    3.006495  1.588083   5.035467 3.132560e-04
#> C-A.205          Abcb9      C-A   74.513454  6.219429   5.028617 3.168023e-04
#> C-A.206           Cd52      C-A    2.550146  1.350580   5.022302 3.201091e-04
#> C-A.207           Chn2      C-A    9.494933  3.247158   5.013348 3.248601e-04
#> C-A.208         Il10ra      C-A   -7.935113 -2.988251  -5.003911 3.299474e-04
#> C-A.209          Spns3      C-A   -7.197287 -2.847453  -5.003263 3.302998e-04
#> C-A.210          Lonp2      C-A    3.185253  1.671408   4.984786 3.405165e-04
#> C-A.211      Trp53inp1      C-A    6.011104  2.587630   4.970688 3.485350e-04
#> C-A.212           Eya2      C-A  434.786676  8.764164   4.963302 3.528148e-04
#> C-A.213           Prr5      C-A  -14.960052 -3.903043  -4.960596 3.543963e-04
#> C-A.214        Clec12a      C-A  -25.056151 -4.647093  -4.959689 3.549285e-04
#> C-A.215           Elp2      C-A   -2.685111 -1.424982  -4.951895 3.595336e-04
#> C-A.216        Fam111a      C-A   -2.512128 -1.328910  -4.935708 3.692992e-04
#> C-A.217            Wls      C-A    3.593386  1.845344   4.926543 3.749513e-04
#> C-A.218          Gcnt1      C-A   13.411502  3.745399   4.926270 3.751209e-04
#> C-A.219          Sept1      C-A    2.248910  1.169226   4.921687 3.779828e-04
#> C-A.220         Ormdl3      C-A    9.795573  3.292130   4.911258 3.845812e-04
#> C-A.221           Cd84      C-A    3.703697  1.888966   4.902364 3.903040e-04
#> C-A.222           Snx9      C-A  -54.908258 -5.778951  -4.886542 4.007065e-04
#> C-A.223         Acot11      C-A  116.750253  6.867282   4.873679 4.093778e-04
#> C-A.224           Bin1      C-A   -2.506305 -1.325562  -4.873531 4.094785e-04
#> C-A.225            Syk      C-A  -15.252166 -3.930942  -4.872541 4.101541e-04
#> C-A.226          Plcg2      C-A   -2.258044 -1.175073  -4.872193 4.103921e-04
#> C-A.227          Anks1      C-A  -18.263378 -4.190882  -4.864818 4.154653e-04
#> C-A.228          Prmt1      C-A   -2.020278 -1.014554  -4.859873 4.189042e-04
#> C-A.229          Sfxn1      C-A   -3.589970 -1.843972  -4.847759 4.274547e-04
#> C-A.230        Erbb2ip      C-A    2.691585  1.428456   4.839421 4.334465e-04
#> C-A.231          Ap3d1      C-A    3.109489  1.636677   4.837320 4.349703e-04
#> C-A.232           Snx2      C-A    2.157389  1.109286   4.832464 4.385128e-04
#> C-A.233            Fh1      C-A   -2.736525 -1.452345  -4.827617 4.420798e-04
#> C-A.234          Ppm1l      C-A    9.204163  3.202287   4.817274 4.497930e-04
#> C-A.235           Pfkp      C-A    3.803500  1.927328   4.814348 4.520002e-04
#> C-A.236          Susd1      C-A   -3.494979 -1.805284  -4.783272 4.761591e-04
#> C-A.237         Apol7e      C-A  160.936669  7.330349   4.782347 4.768988e-04
#> C-A.238        Gramd1b      C-A   27.186977  4.764844   4.780310 4.785315e-04
#> C-A.239         Clint1      C-A    2.134993  1.094231   4.737774 5.139991e-04
#> C-A.240          Ttc13      C-A   -2.858155 -1.515084  -4.729707 5.210307e-04
#> C-A.241        Rasgrp2      C-A   -4.484432 -2.164925  -4.727167 5.232661e-04
#> C-A.242           Mcm2      C-A   -2.130381 -1.091112  -4.727080 5.233427e-04
#> C-A.243        Hnrnpdl      C-A   -4.943657 -2.305579  -4.714579 5.344904e-04
#> C-A.244          Trim8      C-A    2.953051  1.562206   4.710958 5.377660e-04
#> C-A.245          Rpl32      C-A   -2.011536 -1.008297  -4.700277 5.475506e-04
#> C-A.246         Bcl11a      C-A  -13.925537 -3.799661  -4.686371 5.605691e-04
#> C-A.247         Nsmce1      C-A   -3.213690 -1.684231  -4.684017 5.628050e-04
#> C-A.248           Exo1      C-A   -4.484188 -2.164847  -4.670907 5.754278e-04
#> C-A.249          Ikzf1      C-A   -2.072634 -1.051466  -4.669578 5.767241e-04
#> C-A.250         Tspan3      C-A   -2.614090 -1.386309  -4.666389 5.798458e-04
#> C-A.251           Tle3      C-A    2.980971  1.575782   4.662656 5.835233e-04
#> C-A.252         Ifitm3      C-A  -20.565057 -4.362123  -4.650544 5.956247e-04
#> C-A.253         Osbpl3      C-A    5.160513  2.367515   4.649642 5.965369e-04
#> C-A.254         Gm1966      C-A   -5.171529 -2.370591  -4.644965 6.012868e-04
#> C-A.255      Serpina3g      C-A   15.981161  3.998300   4.641063 6.052809e-04
#> C-A.256         Ms4a4b      C-A   92.566102  6.532412   4.639620 6.067651e-04
#> C-A.257          Tifab      C-A  -14.088425 -3.816438  -4.638877 6.075304e-04
#> C-A.258         Topbp1      C-A   -2.148635 -1.103421  -4.626676 6.202469e-04
#> C-A.259           Ppt2      C-A   14.228941  3.830756   4.621318 6.259190e-04
#> C-A.260         Zfp422      C-A   -2.392695 -1.258636  -4.610091 6.379805e-04
#> C-A.261           Pten      C-A   -2.114501 -1.080317  -4.608369 6.398517e-04
#> C-A.262            Ltb      C-A  130.361585  7.026375   4.604534 6.440403e-04
#> C-A.263            Ran      C-A   -2.131510 -1.091876  -4.599634 6.494339e-04
#> C-A.264           Jak2      C-A    3.159766  1.659818   4.599110 6.500133e-04
#> C-A.265          Fads2      C-A  -17.275655 -4.110669  -4.599001 6.501345e-04
#> C-A.266         Rundc1      C-A    3.945507  1.980211   4.598055 6.511823e-04
#> C-A.267        Tmem206      C-A   -8.994371 -3.169022  -4.596686 6.527020e-04
#> C-A.268          Hcls1      C-A    2.088695  1.062602   4.591057 6.589885e-04
#> C-A.269           Nme4      C-A   -3.204656 -1.680170  -4.588139 6.622725e-04
#> C-A.270         Rnf166      C-A    2.913262  1.542635   4.584940 6.658921e-04
#> C-A.271          Srsf2      C-A   -2.802122 -1.486520  -4.581404 6.699184e-04
#> C-A.272          Cd24a      C-A   -8.207473 -3.036938  -4.578111 6.736899e-04
#> C-A.273          Abca1      C-A    6.099215  2.608624   4.569201 6.840058e-04
#> C-A.274         Kansl2      C-A   -2.392107 -1.258282  -4.560260 6.945250e-04
#> C-A.275           Ccr2      C-A  160.847624  7.329551   4.559353 6.956008e-04
#> C-A.276          Dirc2      C-A  -14.388973 -3.846892  -4.556759 6.986896e-04
#> C-A.277         Rcbtb2      C-A    2.967017  1.569013   4.542321 7.161413e-04
#> C-A.278          Rgs14      C-A   64.928107  6.020771   4.539337 7.198053e-04
#> C-A.279           Mafk      C-A   -2.866449 -1.519264  -4.532457 7.283261e-04
#> C-A.280           Rhof      C-A   14.097856  3.817404   4.529542 7.319683e-04
#> C-A.281          Cdk19      C-A   -3.333580 -1.737072  -4.524025 7.389144e-04
#> C-A.282         Cldn25      C-A    2.404687  1.265849   4.523104 7.400798e-04
#> C-A.283          U2af1      C-A   -2.486894 -1.314345  -4.520739 7.430835e-04
#> C-A.284        Atp13a2      C-A   -2.872021 -1.522067  -4.520644 7.432045e-04
#> C-A.285        Rps6ka3      C-A   20.248371  4.339734   4.518702 7.456804e-04
#> C-A.286          Snx10      C-A   -2.829555 -1.500575  -4.508636 7.586536e-04
#> C-A.287          Hmgb3      C-A   -3.002381 -1.586107  -4.502396 7.668134e-04
#> C-A.288  4930523C07Rik      C-A    4.502092  2.170595   4.502095 7.672104e-04
#> C-A.289         Il10rb      C-A   -8.836906 -3.143541  -4.499463 7.706815e-04
#> C-A.290          Gnai3      C-A    2.198944  1.136811   4.481729 7.945016e-04
#> C-A.291           Ctsw      C-A   21.305586  4.413160   4.480676 7.959389e-04
#> C-A.292           Bst2      C-A   -4.621401 -2.208330  -4.479713 7.972575e-04
#> C-A.293        Tmem154      C-A   18.627294  4.219346   4.475767 8.026810e-04
#> C-A.294          Cdyl2      C-A   67.132746  6.068945   4.475667 8.028191e-04
#> C-A.295         Sema4a      C-A    9.907828  3.308569   4.447368 8.428667e-04
#> C-A.296          Kdm7a      C-A   13.257511  3.728738   4.445698 8.452938e-04
#> C-A.297          Cstf2      C-A   -2.308821 -1.207156  -4.443742 8.481467e-04
#> C-A.298         Sept11      C-A    4.037908  2.013608   4.443230 8.488958e-04
#> C-A.299           Il7r      C-A    4.983727  2.317225   4.440963 8.522167e-04
#> C-A.300           Bin2      C-A    2.197528  1.135882   4.439683 8.540985e-04
#> C-A.301         Osbpl5      C-A    9.685534  3.275832   4.427854 8.716921e-04
#> C-A.302          Myo1c      C-A   -8.280165 -3.049660  -4.418999 8.851092e-04
#> C-A.303           Tob1      C-A   50.765598  5.665779   4.418326 8.861392e-04
#> C-A.304          Wdr18      C-A   -2.744371 -1.456475  -4.384302 9.397945e-04
#> C-A.305          Anxa6      C-A    3.536980  1.822518   4.377663 9.506534e-04
#> C-A.306        Slc29a1      C-A   -3.092369 -1.628712  -4.374964 9.551049e-04
#> C-A.307        Unc119b      C-A    2.362281  1.240181   4.369313 9.644966e-04
#> C-A.308           Gart      C-A   -2.075635 -1.053553  -4.368995 9.650279e-04
#> C-A.309        Slc11a2      C-A   14.316502  3.839607   4.366475 9.692500e-04
#> C-A.310          Snx14      C-A  142.685579  7.156696   4.362986 9.751249e-04
#> C-A.311         Tspan2      C-A  -65.730903 -6.038500  -4.361559 9.775400e-04
#> C-A.312          Snx30      C-A   -5.676568 -2.505019  -4.361466 9.776976e-04
#> C-A.313           Uap1      C-A   -8.484652 -3.084856  -4.358918 9.820248e-04
#> C-A.314       Colgalt1      C-A   -2.621331 -1.390300  -4.358710 9.823778e-04
#> C-A.315          N4bp3      C-A  -10.439523 -3.383984  -4.356389 9.863377e-04
#> C-A.316         Chchd4      C-A   -3.683870 -1.881222  -4.348727 9.995303e-04
#> C-A.317          Pdcd1      C-A  763.492227  9.576470   4.347209 1.002165e-03
#> C-A.318          Abhd8      C-A    5.769147  2.528358   4.342371 1.010612e-03
#> C-A.319          H2-Aa      C-A   77.015886  6.267084   4.340253 1.014333e-03
#> C-A.320          Nrip1      C-A    2.646015  1.403821   4.339535 1.015598e-03
#> C-A.321          Rrp1b      C-A   -2.570189 -1.361875  -4.337697 1.018844e-03
#> C-A.322           Mcm3      C-A   -2.453113 -1.294614  -4.333175 1.026873e-03
#> C-A.323        Tmem164      C-A    2.067271  1.047728   4.332862 1.027431e-03
#> C-A.324           Cd96      C-A    2.322577  1.215726   4.330869 1.030993e-03
#> C-A.325          Med12      C-A    3.754710  1.908701   4.323772 1.043782e-03
#> C-A.326          Psmb9      C-A   -2.314608 -1.210768  -4.323226 1.044772e-03
#> C-A.327         Tmem43      C-A   -3.911393 -1.967683  -4.321129 1.048588e-03
#> C-A.328           Hat1      C-A   -2.182037 -1.125676  -4.319393 1.051755e-03
#> C-A.329         Tbc1d5      C-A   -4.718855 -2.238437  -4.319025 1.052428e-03
#> C-A.330         Tgfbr2      C-A    6.331739  2.662602   4.309536 1.069936e-03
#> C-A.331         Plxnb2      C-A  -10.621308 -3.408890  -4.305265 1.077915e-03
#> C-A.332           Spi1      C-A  -20.088892 -4.328326  -4.303720 1.080816e-03
#> C-A.333         Eif4e3      C-A    5.821617  2.541420   4.301710 1.084602e-03
#> C-A.334          Slfn2      C-A   73.777186  6.205103   4.299834 1.088149e-03
#> C-A.335          Mknk2      C-A    2.103126  1.072535   4.288220 1.110378e-03
#> C-A.336           Fgl2      C-A    7.767955  2.957535   4.287009 1.112723e-03
#> C-A.337          Capn2      C-A    5.270749  2.398008   4.269021 1.148165e-03
#> C-A.338         Zfp945      C-A   21.531074  4.428348   4.258142 1.170167e-03
#> C-A.339          Enpp4      C-A    4.659565  2.220195   4.256229 1.174081e-03
#> C-A.340         Nhp2l1      C-A   -2.058303 -1.041455  -4.238875 1.210218e-03
#> C-A.341           Sat1      C-A    2.844042  1.507943   4.238343 1.211344e-03
#> C-A.342          Rps17      C-A   -2.130116 -1.090932  -4.236559 1.215126e-03
#> C-A.343          Tra2b      C-A   -2.050048 -1.035658  -4.228265 1.232879e-03
#> C-A.344         Cyfip1      C-A   -2.304323 -1.204343  -4.222423 1.245543e-03
#> C-A.345        Sh3kbp1      C-A    2.464016  1.301012   4.220808 1.249068e-03
#> C-A.346        Il12rb1      C-A   74.898466  6.226864   4.217405 1.256530e-03
#> C-A.347          Cisd3      C-A   16.979260  4.085702   4.214663 1.262574e-03
#> C-A.348          Glrx3      C-A   -2.100272 -1.070576  -4.212880 1.266523e-03
#> C-A.349          Il2rb      C-A 1718.312382 10.746777   4.211999 1.268477e-03
#> C-A.350          Litaf      C-A  -10.749686 -3.426223  -4.201416 1.292211e-03
#> C-A.351        Alox5ap      C-A  -34.386402 -5.103766  -4.197201 1.301791e-03
#> C-A.352         Lpcat3      C-A    2.212992  1.145998   4.189561 1.319348e-03
#> C-A.353          Timp2      C-A    6.625000  2.727920   4.188563 1.321659e-03
#> C-A.354        Fam134b      C-A   -4.338539 -2.117209  -4.187258 1.324687e-03
#> C-A.355       Tmem176a      C-A    2.569416  1.361440   4.186538 1.326360e-03
#> C-A.356          Egln1      C-A   -3.046793 -1.607292  -4.185618 1.328503e-03
#> C-A.357          Mpeg1      C-A  -28.146667 -4.814892  -4.176015 1.351078e-03
#> C-A.358          Mppe1      C-A    4.281402  2.098083   4.169102 1.367575e-03
#> C-A.359           Dtx4      C-A   -4.321368 -2.111488  -4.165133 1.377142e-03
#> C-A.360           Cry1      C-A    2.452663  1.294349   4.157811 1.394973e-03
#> C-A.361  1110034G24Rik      C-A   36.639803  5.195340   4.153497 1.405591e-03
#> C-A.362           Pcna      C-A   -2.141732 -1.098778  -4.148228 1.418672e-03
#> C-A.363         Hivep1      C-A   34.863256  5.123635   4.134892 1.452348e-03
#> C-A.364            Rb1      C-A    4.656263  2.219173   4.132841 1.457603e-03
#> C-A.365          Ccng1      C-A    2.369929  1.244844   4.127424 1.471569e-03
#> C-A.366          Rps25      C-A   -2.054929 -1.039088  -4.124043 1.480356e-03
#> C-A.367         Ruvbl2      C-A   -2.020173 -1.014479  -4.122436 1.484554e-03
#> C-A.368           Inip      C-A   -3.090788 -1.627975  -4.120610 1.489337e-03
#> C-A.369           Zeb2      C-A  -25.328096 -4.662667  -4.120213 1.490377e-03
#> C-A.370           Elf4      C-A    2.542894  1.346471   4.119753 1.491586e-03
#> C-A.371         Tmem71      C-A   15.463137  3.950761   4.118061 1.496040e-03
#> C-A.372         Sqstm1      C-A    2.101706  1.071561   4.112781 1.510024e-03
#> C-A.373          Cdip1      C-A   -3.703599 -1.888928  -4.105308 1.530051e-03
#> C-A.374          Plcb4      C-A    4.135163  2.047944   4.091588 1.567537e-03
#> C-A.375        Pla2g15      C-A   -5.622222 -2.491140  -4.086886 1.580600e-03
#> C-A.376           Rgs3      C-A   49.773705  5.637312   4.085848 1.583499e-03
#> C-A.377          Adpgk      C-A    5.033130  2.331456   4.078882 1.603101e-03
#> C-A.378          Zfp35      C-A   -2.345330 -1.229791  -4.075487 1.612745e-03
#> C-A.379         Impdh1      C-A    2.041595  1.029697   4.075286 1.613317e-03
#> C-A.380         Fbxo33      C-A    7.880359  2.978261   4.071921 1.622937e-03
#> C-A.381           Rfc4      C-A   -2.492754 -1.317741  -4.060320 1.656564e-03
#> C-A.382        Smarcc1      C-A   -2.079646 -1.056338  -4.059435 1.659158e-03
#> C-A.383         Nmral1      C-A   -5.875340 -2.554672  -4.056800 1.666909e-03
#> C-A.384          Cenph      C-A   -6.675015 -2.738771  -4.056485 1.667838e-03
#> C-A.385           Msh6      C-A   -2.601255 -1.379208  -4.055909 1.669541e-03
#> C-A.386          Smad3      C-A    4.396590  2.136385   4.050981 1.684160e-03
#> C-A.387          Mast4      C-A    8.994435  3.169033   4.044557 1.703417e-03
#> C-A.388        Tsc22d1      C-A  -46.631853 -5.543244  -4.043886 1.705442e-03
#> C-A.389         Ogfrl1      C-A   -3.776521 -1.917058  -4.043467 1.706707e-03
#> C-A.390           Rfc2      C-A   -2.335717 -1.223866  -4.041755 1.711888e-03
#> C-A.391         Atp8a1      C-A    4.125118  2.044435   4.040663 1.715201e-03
#> C-A.392          Lmnb1      C-A   -2.199507 -1.137180  -4.035053 1.732325e-03
#> C-A.393         Erlin1      C-A   -4.149142 -2.052813  -4.031275 1.743956e-03
#> C-A.394         Exosc8      C-A   -2.882736 -1.527439  -4.030443 1.746532e-03
#> C-A.395         Nt5dc2      C-A   -2.738753 -1.453519  -4.025207 1.762809e-03
#> C-A.396           H1f0      C-A    2.976397  1.573567   4.021633 1.774014e-03
#> C-A.397         Ndufa4      C-A   -2.184537 -1.127327  -4.019621 1.780353e-03
#> C-A.398          Srpk1      C-A    2.137310  1.095796   4.017363 1.787495e-03
#> C-A.399         Tuba4a      C-A   -2.017452 -1.012534  -4.014840 1.795509e-03
#> C-A.400          Thtpa      C-A   26.397364  4.722322   4.011267 1.806922e-03
#> C-A.401           Rnf6      C-A    2.700214  1.433074   4.004287 1.829439e-03
#> C-A.402         Ruvbl1      C-A   -2.551900 -1.351572  -4.002963 1.833742e-03
#> C-A.403         Cdkn1a      C-A   -6.417589 -2.682031  -3.996654 1.854394e-03
#> C-A.404        Wbscr27      C-A   40.592483  5.343141   3.994355 1.861977e-03
#> C-A.405           Vasp      C-A    3.641042  1.864351   3.993256 1.865617e-03
#> C-A.406          Pole2      C-A   -7.171801 -2.842335  -3.992715 1.867409e-03
#> C-A.407          Cdca5      C-A   -2.741568 -1.455001  -3.988781 1.880501e-03
#> C-A.408          Spsb3      C-A   21.842423  4.449061   3.987176 1.885870e-03
#> C-A.409         Hspbp1      C-A   -2.580426 -1.367609  -3.983100 1.899576e-03
#> C-A.410         Bnip3l      C-A    2.301289  1.202442   3.979407 1.912084e-03
#> C-A.411       Arhgap17      C-A   -4.102460 -2.036489  -3.974902 1.927457e-03
#> C-A.412        Ankrd44      C-A    2.742886  1.455695   3.971565 1.938927e-03
#> C-A.413          Ddit4      C-A   19.149240  4.259215   3.970026 1.944241e-03
#> C-A.414            Cd9      C-A   10.265060  3.359670   3.966088 1.957905e-03
#> C-A.415           Zbp1      C-A   59.196477  5.887439   3.964475 1.963528e-03
#> C-A.416          Yipf3      C-A    2.282986  1.190922   3.960303 1.978156e-03
#> C-A.417          Plcg1      C-A   -2.899855 -1.535981  -3.956219 1.992587e-03
#> C-A.418          Gpsm3      C-A    8.581071  3.101158   3.955766 1.994194e-03
#> C-A.419        Spata13      C-A    3.172266  1.665514   3.954122 2.000036e-03
#> C-A.420          Exoc6      C-A   -3.784783 -1.920211  -3.948285 2.020926e-03
#> C-A.421            Nnt      C-A   -2.018062 -1.012970  -3.945664 2.030379e-03
#> C-A.422          Ccnb1      C-A   -2.248382 -1.168887  -3.945541 2.030827e-03
#> C-A.423         Samhd1      C-A    9.664341  3.272671   3.938204 2.057538e-03
#> C-A.424          Csf1r      C-A  -54.423223 -5.766151  -3.936165 2.065027e-03
#> C-A.425          Golm1      C-A   -3.104947 -1.634569  -3.935897 2.066016e-03
#> C-A.426          Sh2b3      C-A   -3.660402 -1.872002  -3.931467 2.082390e-03
#> C-A.427           Gapt      C-A   -6.139380 -2.618093  -3.929503 2.089692e-03
#> C-A.428           Ly86      C-A  -17.642651 -4.140995  -3.922957 2.114222e-03
#> C-A.429          Tpst2      C-A    2.155977  1.108342   3.922937 2.114296e-03
#> C-A.430           Cyba      C-A    2.155183  1.107810   3.918236 2.132098e-03
#> C-A.431          Sort1      C-A   -3.379932 -1.756994  -3.918095 2.132633e-03
#> C-A.432          Adap1      C-A  -26.050653 -4.703248  -3.916809 2.137530e-03
#> C-A.433          Prkcd      C-A   -2.581168 -1.368024  -3.915962 2.140763e-03
#> C-A.434          Wasf2      C-A   -2.577338 -1.365882  -3.912247 2.154997e-03
#> C-A.435         Tmem64      C-A    6.926344  2.792094   3.908110 2.170965e-03
#> C-A.436          Stat3      C-A    2.392376  1.258444   3.907743 2.172389e-03
#> C-A.437           Tjp3      C-A    8.352535  3.062214   3.902045 2.194594e-03
#> C-A.438           Cast      C-A    3.400921  1.765926   3.900712 2.199824e-03
#> C-A.439          Nfam1      C-A   -2.800269 -1.485566  -3.893068 2.230059e-03
#> C-A.440           Rinl      C-A    3.584561  1.841796   3.887816 2.251080e-03
#> C-A.441            Pbk      C-A   -3.113011 -1.638311  -3.886840 2.255008e-03
#> C-A.442          Ube2h      C-A    3.314046  1.728594   3.886235 2.257449e-03
#> C-A.443          Spidr      C-A   -5.401146 -2.433266  -3.885176 2.261727e-03
#> C-A.444           Gbp2      C-A    6.358521  2.668691   3.883743 2.267525e-03
#> C-A.445         Scarb1      C-A   -4.374477 -2.129111  -3.883617 2.268035e-03
#> C-A.446          Mgea5      C-A    2.351889  1.233820   3.882121 2.274108e-03
#> C-A.447         Plscr3      C-A   -2.086509 -1.061091  -3.877356 2.293560e-03
#> C-A.448          Rab37      C-A    2.570940  1.362296   3.867679 2.333594e-03
#> C-A.449           Crot      C-A    3.063551  1.615205   3.862897 2.353644e-03
#> C-A.450           B9d2      C-A   -5.217426 -2.383338  -3.862812 2.354002e-03
#> C-A.451         Trim44      C-A   -3.627232 -1.858869  -3.859810 2.366678e-03
#> C-A.452        Nsmce4a      C-A   -2.608880 -1.383430  -3.856677 2.379989e-03
#> C-A.453          Ypel3      C-A    2.577000  1.365693   3.854922 2.387476e-03
#> C-A.454          Nop56      C-A   -2.290904 -1.195917  -3.850433 2.406738e-03
#> C-A.455           Phb2      C-A   -2.028928 -1.020718  -3.849617 2.410256e-03
#> C-A.456         Dnajb1      C-A   -2.409966 -1.269013  -3.845131 2.429697e-03
#> C-A.457      Rab11fip1      C-A   25.403049  4.666930   3.836901 2.465785e-03
#> C-A.458           Rrs1      C-A   -2.714568 -1.440723  -3.834454 2.476621e-03
#> C-A.459          Spns2      C-A  -33.183292 -5.052385  -3.833372 2.481427e-03
#> C-A.460          Xlr4b      C-A   23.153572  4.533163   3.831738 2.488701e-03
#> C-A.461           Amd1      C-A   -2.477934 -1.309138  -3.831663 2.489035e-03
#> C-A.462           Tns1      C-A   -6.390884 -2.676015  -3.828694 2.502320e-03
#> C-A.463           Gnal      C-A   17.738575  4.148818   3.827527 2.507558e-03
#> C-A.464         Mtmr10      C-A   20.184449  4.335172   3.824977 2.519050e-03
#> C-A.465          Aurkb      C-A   -2.542698 -1.346360  -3.823480 2.525821e-03
#> C-A.466          Prex1      C-A    2.209945  1.144010   3.822570 2.529942e-03
#> C-A.467            Pnn      C-A   -2.542099 -1.346020  -3.818371 2.549068e-03
#> C-A.468          Ncoa1      C-A    3.325696  1.733656   3.818278 2.549493e-03
#> C-A.469           Pkp3      C-A   15.067407  3.913359   3.812352 2.576736e-03
#> C-A.470         Polr3f      C-A    5.794688  2.534731   3.811072 2.582662e-03
#> C-A.471         H2-Ab1      C-A    2.798695  1.484754   3.800518 2.632049e-03
#> C-A.472  3110082I17Rik      C-A  -21.936850 -4.455285  -3.796893 2.649233e-03
#> C-A.473         Snrpa1      C-A   -2.229524 -1.156736  -3.794948 2.658500e-03
#> C-A.474           Rif1      C-A   -2.307568 -1.206373  -3.791986 2.672680e-03
#> C-A.475         Tbxa2r      C-A   -5.940124 -2.570493  -3.791730 2.673911e-03
#> C-A.476          Ift57      C-A   -6.823186 -2.770446  -3.786069 2.701242e-03
#> C-A.477           Net1      C-A   -8.712168 -3.123032  -3.785027 2.706303e-03
#> C-A.478          Lzts2      C-A   -7.106728 -2.829185  -3.784165 2.710495e-03
#> C-A.479          Herc3      C-A   23.330248  4.544130   3.779974 2.730988e-03
#> C-A.480        Fam129a      C-A   26.177195  4.710239   3.778440 2.738526e-03
#> C-A.481           Ctc1      C-A   -2.614346 -1.386450  -3.777836 2.741499e-03
#> C-A.482         Pi4k2a      C-A    3.423472  1.775460   3.771315 2.773829e-03
#> C-A.483  4833439L19Rik      C-A   -3.588533 -1.843394  -3.764271 2.809185e-03
#> C-A.484          Lims1      C-A   -2.015974 -1.011477  -3.763539 2.812889e-03
#> C-A.485         Syngr2      C-A   -2.165506 -1.114704  -3.763360 2.813792e-03
#> C-A.486          Sf3a1      C-A   -2.253011 -1.171854  -3.762239 2.819474e-03
#> C-A.487        Dennd5a      C-A   -3.192082 -1.674498  -3.759260 2.834621e-03
#> C-A.488          Trap1      C-A   -2.143654 -1.100072  -3.756720 2.847606e-03
#> C-A.489          Ninj1      C-A   -5.234356 -2.388012  -3.753634 2.863461e-03
#> C-A.490           Paox      C-A    2.502596  1.323426   3.745763 2.904320e-03
#> C-A.491        Cbfa2t3      C-A   -5.264569 -2.396316  -3.745596 2.905191e-03
#> C-A.492          Ncoa3      C-A    2.692438  1.428913   3.745213 2.907196e-03
#> C-A.493       AW112010      C-A    2.544434  1.347345   3.741520 2.926588e-03
#> C-A.494           Rara      C-A    5.754319  2.524645   3.740180 2.933659e-03
#> C-A.495          Morc3      C-A   -2.099866 -1.070297  -3.737450 2.948120e-03
#> C-A.496           Ect2      C-A   -3.008266 -1.588932  -3.734541 2.963605e-03
#> C-A.497          Klf10      C-A    2.004813  1.003467   3.731791 2.978326e-03
#> C-A.498          Tgtp1      C-A    4.347462  2.120173   3.727975 2.998873e-03
#> C-A.499        Dpy19l3      C-A    6.155479  2.621871   3.726161 3.008689e-03
#> C-A.500       Zc3hav1l      C-A   -8.098021 -3.017569  -3.719149 3.046953e-03
#> C-A.501          Ddx54      C-A   -2.229374 -1.156639  -3.717542 3.055793e-03
#> C-A.502       Timeless      C-A   -3.383489 -1.758512  -3.715682 3.066062e-03
#> C-A.503         Suclg1      C-A   -2.528355 -1.338199  -3.715065 3.069474e-03
#> C-A.504         Snrpd1      C-A   -2.224251 -1.153320  -3.712962 3.081132e-03
#> C-A.505         Ccdc50      C-A    2.681511  1.423046   3.707786 3.110031e-03
#> C-A.506          Syce2      C-A  -10.617338 -3.408350  -3.703423 3.134602e-03
#> C-A.507          Ddx27      C-A   -2.463754 -1.300858  -3.694504 3.185456e-03
#> C-A.508          Yipf6      C-A   -4.834619 -2.273402  -3.690665 3.207608e-03
#> C-A.509         Pom121      C-A   -5.509900 -2.462026  -3.686159 3.233807e-03
#> C-A.510           Sord      C-A   -4.949726 -2.307349  -3.685073 3.240153e-03
#> C-A.511           Fgd2      C-A  -18.670175 -4.222664  -3.684864 3.241373e-03
#> C-A.512           Mtf2      C-A   -3.988418 -1.995817  -3.682611 3.254591e-03
#> C-A.513           Abi1      C-A    2.794379  1.482528   3.682058 3.257839e-03
#> C-A.514  5430427O19Rik      C-A   -4.489670 -2.166609  -3.680633 3.266238e-03
#> C-A.515         Mif4gd      C-A    4.128497  2.045617   3.671774 3.318926e-03
#> C-A.516          Cytip      C-A    3.265752  1.707415   3.667879 3.342362e-03
#> C-A.517           Ibtk      C-A   -3.253719 -1.702090  -3.665201 3.358581e-03
#> C-A.518        Galnt12      C-A    5.669149  2.503132   3.662330 3.376053e-03
#> C-A.519         Gimap5      C-A    7.351332  2.878006   3.658094 3.402005e-03
#> C-A.520          Wdr34      C-A  -25.345548 -4.663660  -3.655800 3.416140e-03
#> C-A.521         Rnf141      C-A   -5.290160 -2.403311  -3.655564 3.417599e-03
#> C-A.522          Pola2      C-A   -2.197588 -1.135921  -3.654702 3.422932e-03
#> C-A.523    D19Bwg1357e      C-A   -2.167211 -1.115840  -3.649841 3.453151e-03
#> C-A.524           Klf6      C-A    3.770688  1.914828   3.647010 3.470878e-03
#> C-A.525        Fam134c      C-A    2.478431  1.309427   3.646190 3.476032e-03
#> C-A.526           Umps      C-A   -2.529797 -1.339022  -3.644422 3.487165e-03
#> C-A.527            Vcl      C-A  -12.153483 -3.603298  -3.642829 3.497232e-03
#> C-A.528       Ccdc102a      C-A   33.109825  5.049188   3.640484 3.512098e-03
#> C-A.529         Setdb2      C-A   18.543797  4.212865   3.640264 3.513496e-03
#> C-A.530           Tyms      C-A   -2.150558 -1.104711  -3.640083 3.514647e-03
#> C-A.531         Clec2i      C-A    9.604607  3.263727   3.636405 3.538116e-03
#> C-A.532          Erp29      C-A   -2.187075 -1.129003  -3.634905 3.547733e-03
#> C-A.533        Skiv2l2      C-A   -2.048384 -1.034486  -3.631630 3.568829e-03
#> C-A.534         Tsen54      C-A   -4.665364 -2.221990  -3.630439 3.576532e-03
#> C-A.535          Hvcn1      C-A    3.273893  1.711007   3.624512 3.615114e-03
#> C-A.536          Panx1      C-A   -2.607284 -1.382548  -3.624005 3.618429e-03
#> C-A.537       Cyb561a3      C-A   -2.426317 -1.278768  -3.623598 3.621102e-03
#> C-A.538           Rpf2      C-A   -2.967458 -1.569228  -3.622108 3.630885e-03
#> C-A.539        Prelid2      C-A   16.091293  4.008208   3.619312 3.649312e-03
#> C-A.540          Itgam      C-A  -19.642330 -4.295894  -3.614914 3.678501e-03
#> C-A.541         Maged1      C-A    5.948508  2.572528   3.611627 3.700465e-03
#> C-A.542         Zfand6      C-A    2.185191  1.127760   3.609434 3.715199e-03
#> C-A.543           Ncln      C-A   -2.610801 -1.384493  -3.607438 3.728659e-03
#> C-A.544           Sypl      C-A   -3.302497 -1.723557  -3.599441 3.783092e-03
#> C-A.545          Galk2      C-A    2.543822  1.346998   3.596329 3.804491e-03
#> C-A.546          Smap1      C-A    5.680657  2.506058   3.594297 3.818534e-03
#> C-A.547        Dennd1b      C-A    4.716305  2.237657   3.583272 3.895647e-03
#> C-A.548            Sp4      C-A   25.195531  4.655096   3.582758 3.899283e-03
#> C-A.549          Trit1      C-A   -4.419854 -2.143999  -3.582641 3.900108e-03
#> C-A.550        Plekho2      C-A    2.140156  1.097716   3.581189 3.910392e-03
#> C-A.551          Sel1l      C-A    2.459837  1.298562   3.576716 3.942253e-03
#> C-A.552          Nsun2      C-A   -2.101309 -1.071288  -3.575705 3.949495e-03
#> C-A.553          Numa1      C-A   -2.065158 -1.046252  -3.573717 3.963764e-03
#> C-A.554          Trip4      C-A    3.037345  1.602811   3.571929 3.976648e-03
#> C-A.555         Sptlc1      C-A    3.233944  1.693295   3.568264 4.003184e-03
#> C-A.556           Utrn      C-A    7.092137  2.826220   3.566168 4.018443e-03
#> C-A.557          Acss2      C-A    8.928406  3.158403   3.564905 4.027670e-03
#> C-A.558         H2-Eb1      C-A  141.292279  7.142539   3.564300 4.032093e-03
#> C-A.559           Suox      C-A    9.114080  3.188097   3.561749 4.050810e-03
#> C-A.560         Agpat3      C-A    3.297145  1.721217   3.561354 4.053716e-03
#> C-A.561           Cog4      C-A   -4.836034 -2.273824  -3.559194 4.069643e-03
#> C-A.562         Atrnl1      C-A    7.088502  2.825481   3.558900 4.071814e-03
#> C-A.563        Sec14l1      C-A   12.086088  3.595275   3.558567 4.074277e-03
#> C-A.564          Rnf14      C-A    2.869214  1.520655   3.548683 4.148070e-03
#> C-A.565          Etnk1      C-A   -2.003912 -1.002819  -3.545770 4.170080e-03
#> C-A.566           Ctps      C-A   -2.479760 -1.310200  -3.544964 4.176193e-03
#> C-A.567          Dusp2      C-A  -32.759013 -5.033820  -3.542090 4.198050e-03
#> C-A.568            Hlf      C-A    5.605651  2.486882   3.529709 4.293582e-03
#> C-A.569        Creb3l2      C-A   54.482921  5.767732   3.528553 4.302617e-03
#> C-A.570          Gng12      C-A   -2.829050 -1.500318  -3.527558 4.310403e-03
#> C-A.571           Cd48      C-A   -2.783891 -1.477103  -3.527276 4.312618e-03
#> C-A.572         Unc119      C-A   14.603358  3.868228   3.523803 4.339935e-03
#> C-A.573          Tram1      C-A   -2.173038 -1.119713  -3.523801 4.339948e-03
#> C-A.574          Ikbke      C-A   19.604617  4.293122   3.520661 4.364797e-03
#> C-A.575         Klrb1f      C-A    3.775871  1.916809   3.518017 4.385836e-03
#> C-A.576         Arrdc4      C-A   24.935385  4.640123   3.516436 4.398463e-03
#> C-A.577         Vps37b      C-A   25.699708  4.683680   3.509872 4.451294e-03
#> C-A.578          Nabp1      C-A   37.446603  5.226763   3.507828 4.467876e-03
#> C-A.579          Tiam1      C-A   30.295032  4.921009   3.504755 4.492932e-03
#> C-A.580          Naa15      C-A   -4.565645 -2.190819  -3.503918 4.499775e-03
#> C-A.581          Abcb8      C-A   -2.603286 -1.380334  -3.499915 4.532671e-03
#> C-A.582           Pigs      C-A    2.237158  1.161667   3.499655 4.534819e-03
#> C-A.583          Brca1      C-A   -2.746306 -1.457492  -3.496853 4.558008e-03
#> C-A.584           Cdk6      C-A   -2.532136 -1.340355  -3.494886 4.574356e-03
#> C-A.585           Smox      C-A    3.849048  1.944502   3.494881 4.574391e-03
#> C-A.586           Cds2      C-A    2.941936  1.556766   3.494197 4.580090e-03
#> C-A.587         Myo18a      C-A   -2.035523 -1.025400  -3.494111 4.580811e-03
#> C-A.588          Wdr12      C-A   -2.273529 -1.184933  -3.493444 4.586378e-03
#> C-A.589           Arl1      C-A   -2.137112 -1.095663  -3.492325 4.595728e-03
#> C-A.590          Ckap4      C-A  -24.300002 -4.602884  -3.491273 4.604539e-03
#> C-A.591         Asrgl1      C-A   -2.370787 -1.245366  -3.485794 4.650694e-03
#> C-A.592          Ngly1      C-A    2.321868  1.215286   3.484929 4.658031e-03
#> C-A.593         Spata6      C-A    3.445773  1.784827   3.478045 4.716793e-03
#> C-A.594          Itpr2      C-A   13.929863  3.800109   3.473746 4.753875e-03
#> C-A.595         Kif21b      C-A   26.677626  4.737558   3.473083 4.759622e-03
#> C-A.596          Zzef1      C-A   43.750954  5.451243   3.472796 4.762102e-03
#> C-A.597         Rnase6      C-A  -16.194018 -4.017389  -3.471606 4.772443e-03
#> C-A.598          Arl4a      C-A    5.244261  2.390740   3.469064 4.794595e-03
#> C-A.599         Sppl2b      C-A    2.911134  1.541581   3.467474 4.808506e-03
#> C-A.600           Acp1      C-A   -2.588133 -1.371912  -3.466813 4.814298e-03
#> C-A.601          Plin3      C-A    3.467864  1.794048   3.465745 4.823681e-03
#> C-A.602           Ncf2      C-A   -3.239380 -1.695718  -3.465268 4.827872e-03
#> C-A.603          Jade1      C-A   -3.722498 -1.896271  -3.465027 4.829988e-03
#> C-A.604         Ms4a6c      C-A  -10.399235 -3.378406  -3.463154 4.846502e-03
#> C-A.605          Casp9      C-A    3.302006  1.723343   3.460827 4.867105e-03
#> C-A.606         Zfp763      C-A    5.528794  2.466965   3.456808 4.902881e-03
#> C-A.607         Jmjd1c      C-A    2.685988  1.425453   3.449975 4.964340e-03
#> C-A.608         Tgoln1      C-A    2.147354  1.102560   3.444908 5.010414e-03
#> C-A.609         Plxnc1      C-A   -4.565732 -2.190846  -3.443763 5.020883e-03
#> C-A.610        Slc31a1      C-A   -4.189411 -2.066747  -3.443578 5.022584e-03
#> C-A.611          Reps1      C-A   -2.165920 -1.114980  -3.442586 5.031678e-03
#> C-A.612         Smim13      C-A    5.315525  2.410212   3.439132 5.063471e-03
#> C-A.613         Ndufa9      C-A   -2.154970 -1.107668  -3.438487 5.069432e-03
#> C-A.614            Mbp      C-A    2.440092  1.286936   3.437962 5.074289e-03
#> C-A.615        Gpatch3      C-A   30.699102  4.940125   3.437929 5.074593e-03
#> C-A.616         Inpp5f      C-A   -5.884274 -2.556864  -3.434161 5.109589e-03
#> C-A.617        Fam212a      C-A   -9.464488 -3.242524  -3.416243 5.279413e-03
#> C-A.618           Palm      C-A   -3.671076 -1.876203  -3.412022 5.320246e-03
#> C-A.619           Pus7      C-A   -8.611077 -3.106194  -3.410135 5.338604e-03
#> C-A.620           Rrad      C-A   26.332374  4.718766   3.405769 5.381327e-03
#> C-A.621           Chka      C-A   -4.391032 -2.134560  -3.404798 5.390871e-03
#> C-A.622          Sart3      C-A   -3.500618 -1.807609  -3.402598 5.412577e-03
#> C-A.623         Pcmtd1      C-A    3.238362  1.695264   3.402101 5.417489e-03
#> C-A.624          Irgm2      C-A    4.386441  2.133051   3.401854 5.419934e-03
#> C-A.625            Mn1      C-A   -7.153029 -2.838554  -3.401460 5.423835e-03
#> C-A.626          Il4ra      C-A   -2.669714 -1.416685  -3.401333 5.425092e-03
#> C-A.627          Arl11      C-A   -3.537652 -1.822792  -3.400602 5.432336e-03
#> C-A.628          Dusp6      C-A   -3.115722 -1.639566  -3.395993 5.478256e-03
#> C-A.629          Usp24      C-A    3.388168  1.760505   3.395788 5.480309e-03
#> C-A.630        Rapgef3      C-A   30.070205  4.910263   3.392616 5.512148e-03
#> C-A.631          Mtfr2      C-A   -3.126541 -1.644567  -3.391531 5.523085e-03
#> C-A.632           Ier2      C-A    2.440640  1.287260   3.391213 5.526298e-03
#> C-A.633           Pld4      C-A   -5.368732 -2.424581  -3.389557 5.543037e-03
#> C-A.634          Ncor1      C-A    2.383110  1.252845   3.388413 5.554632e-03
#> C-A.635         Vps37a      C-A    4.836660  2.274011   3.384531 5.594169e-03
#> C-A.636           Gga2      C-A   -5.077656 -2.344163  -3.383148 5.608324e-03
#> C-A.637           Cptp      C-A   22.177415  4.471019   3.382803 5.611862e-03
#> C-A.638        St3gal6      C-A   27.661193  4.789791   3.377471 5.666802e-03
#> C-A.639         Ms4a4c      C-A   24.870603  4.636370   3.376743 5.674354e-03
#> C-A.640         Ankib1      C-A   20.312022  4.344262   3.375034 5.692100e-03
#> C-A.641        Laptm4b      C-A    2.392615  1.258588   3.370320 5.741350e-03
#> C-A.642        Akirin1      C-A    2.069110  1.049011   3.369889 5.745879e-03
#> C-A.643        Arfgef2      C-A   -6.714759 -2.747336  -3.369498 5.749989e-03
#> C-A.644         Stoml2      C-A   -2.023644 -1.016956  -3.369121 5.753954e-03
#> C-A.645           Nrgn      C-A    4.170753  2.060308   3.368975 5.755488e-03
#> C-A.646           Cst7      C-A   -2.482775 -1.311953  -3.365138 5.795994e-03
#> C-A.647          Casd1      C-A   -3.032475 -1.600496  -3.363830 5.809870e-03
#> C-A.648           St14      C-A   -5.679473 -2.505757  -3.358862 5.862883e-03
#> C-A.649          Itfg3      C-A    2.848853  1.510381   3.358825 5.863282e-03
#> C-A.650          H2-Q7      C-A    2.057951  1.041208   3.357694 5.875423e-03
#> C-A.651           Ezh1      C-A    2.631024  1.395625   3.357173 5.881027e-03
#> C-A.652          Cyth4      C-A   -2.180788 -1.124850  -3.357033 5.882529e-03
#> C-A.653            Pxk      C-A    3.134989  1.648460   3.355271 5.901508e-03
#> C-A.654        Trmt112      C-A   -2.235758 -1.160764  -3.354173 5.913377e-03
#> C-A.655           Gbp3      C-A   -8.400010 -3.070391  -3.350836 5.949567e-03
#> C-A.656       Slc25a24      C-A    5.647198  2.497535   3.349023 5.969330e-03
#> C-A.657          Eps15      C-A    2.512430  1.329083   3.345338 6.009707e-03
#> C-A.658           Apoe      C-A   32.378194  5.016951   3.341013 6.057440e-03
#> C-A.659  RP23-330G24.3      C-A   26.163408  4.709479   3.338401 6.086450e-03
#> C-A.660          Trpv2      C-A    3.435270  1.780423   3.336743 6.104944e-03
#> C-A.661          Lin54      C-A   -6.777782 -2.760813  -3.332386 6.153808e-03
#> C-A.662         Ifnar2      C-A   -3.830594 -1.937568  -3.332202 6.155875e-03
#> C-A.663         Ankmy2      C-A    2.954366  1.562849   3.326642 6.218836e-03
#> C-A.664          Bub1b      C-A   -2.583549 -1.369354  -3.323575 6.253843e-03
#> C-A.665         Hmbox1      C-A   14.444652  3.852464   3.322759 6.263182e-03
#> C-A.666         Tmbim1      C-A   -5.937054 -2.569747  -3.320853 6.285074e-03
#> C-A.667           Cnr2      C-A   52.518656  5.714758   3.319830 6.296856e-03
#> C-A.668           Phip      C-A    2.067433  1.047840   3.318424 6.313081e-03
#> C-A.669           Msra      C-A   -2.295174 -1.198603  -3.317444 6.324412e-03
#> C-A.670          Tmed4      C-A   -2.340374 -1.226739  -3.316773 6.332191e-03
#> C-A.671         R3hdm1      C-A   -3.174120 -1.666357  -3.314702 6.356242e-03
#> C-A.672         Mgat4a      C-A    2.475811  1.307901   3.312989 6.376212e-03
#> C-A.673         Gtpbp4      C-A   -2.173163 -1.119796  -3.311611 6.392312e-03
#> C-A.674         Crebrf      C-A    6.316936  2.659225   3.307030 6.446153e-03
#> C-A.675        Rsl24d1      C-A   -2.284176 -1.191674  -3.302064 6.505038e-03
#> C-A.676         Zfp263      C-A   -2.099462 -1.070020  -3.298880 6.543088e-03
#> C-A.677        Sipa1l1      C-A    2.394460  1.259700   3.297237 6.562809e-03
#> C-A.678           Mycl      C-A  -19.226431 -4.265019  -3.294112 6.600473e-03
#> C-A.679          Rad18      C-A   -3.985993 -1.994939  -3.293741 6.604957e-03
#> C-A.680          Frmd6      C-A   -4.813669 -2.267137  -3.291623 6.630634e-03
#> C-A.681      Gabarapl1      C-A   -2.762409 -1.465927  -3.290044 6.649845e-03
#> C-A.682          Creb1      C-A    2.203058  1.139508   3.288661 6.666711e-03
#> C-A.683          Zfp64      C-A    4.630724  2.211238   3.281511 6.754621e-03
#> C-A.684         Fam60a      C-A   -2.038271 -1.027346  -3.276143 6.821380e-03
#> C-A.685         C2cd2l      C-A   30.726426  4.941408   3.276114 6.821741e-03
#> C-A.686           Rxra      C-A   20.534229  4.359959   3.275414 6.830506e-03
#> C-A.687         Fam69b      C-A   -3.733595 -1.900566  -3.275179 6.833445e-03
#> C-A.688          Sepp1      C-A    2.192786  1.132765   3.274872 6.837283e-03
#> C-A.689     St6galnac4      C-A   -2.504553 -1.324553  -3.271231 6.883063e-03
#> C-A.690         Ankra2      C-A   12.234673  3.612904   3.268710 6.914935e-03
#> C-A.691          Plod1      C-A    8.974941  3.165902   3.266299 6.945567e-03
#> C-A.692          Btbd2      C-A    2.101530  1.071440   3.264247 6.971729e-03
#> C-A.693           Klf9      C-A   31.336202  4.969758   3.263719 6.978490e-03
#> C-A.694          Psen1      C-A    2.021641  1.015527   3.258319 7.047906e-03
#> C-A.695         Mrps10      C-A   -3.638486 -1.863338  -3.252364 7.125278e-03
#> C-A.696          Rbpms      C-A    2.602151  1.379705   3.247087 7.194547e-03
#> C-A.697           Gatm      C-A  -20.704852 -4.371897  -3.240411 7.283170e-03
#> C-A.698          Erp44      C-A    2.083712  1.059156   3.238647 7.306776e-03
#> C-A.699          Uqcrq      C-A    2.150184  1.104460   3.237505 7.322089e-03
#> C-A.700  1110008F13Rik      C-A   -2.153202 -1.106484  -3.236461 7.336117e-03
#> C-A.701          Sytl1      C-A    4.191572  2.067491   3.235492 7.349171e-03
#> C-A.702         Hectd1      C-A    2.718564  1.442845   3.234984 7.356026e-03
#> C-A.703        Aldh4a1      C-A   21.136305  4.401651   3.234096 7.368008e-03
#> C-A.704         H2-T22      C-A    2.114431  1.080269   3.233777 7.372332e-03
#> C-A.705       Tnfrsf18      C-A   20.669712  4.369446   3.232348 7.391674e-03
#> C-A.706          Dnm1l      C-A   -2.611516 -1.384887  -3.231257 7.406486e-03
#> C-A.707          Brip1      C-A   -5.056027 -2.338004  -3.228212 7.447971e-03
#> C-A.708         Anp32b      C-A   -2.123414 -1.086386  -3.228086 7.449694e-03
#> C-A.709         Golga3      C-A    3.401177  1.766034   3.225879 7.479911e-03
#> C-A.710           H6pd      C-A    4.542613  2.183523   3.219144 7.572906e-03
#> C-A.711           Nfe2      C-A   -2.229307 -1.156596  -3.213008 7.658652e-03
#> C-A.712          H2-Q6      C-A    2.055970  1.039819   3.212284 7.668829e-03
#> C-A.713           Rfc5      C-A   -2.304300 -1.204329  -3.206808 7.746279e-03
#> C-A.714          Rap2c      C-A    6.084326  2.605097   3.206567 7.749711e-03
#> C-A.715         Rnf168      C-A  -12.016697 -3.586969  -3.206358 7.752676e-03
#> C-A.716          Stim2      C-A    4.506294  2.171942   3.204879 7.773744e-03
#> C-A.717           Pkig      C-A   -3.453780 -1.788176  -3.203678 7.790901e-03
#> C-A.718          Tnip1      C-A    3.232385  1.692599   3.196465 7.894730e-03
#> C-A.719         Polr3h      C-A   -3.551069 -1.828253  -3.193687 7.935087e-03
#> C-A.720           Nasp      C-A   -3.071356 -1.618876  -3.192793 7.948122e-03
#> C-A.721            Ahr      C-A    5.812647  2.539195   3.189132 8.001713e-03
#> C-A.722           Ccz1      C-A   -2.042472 -1.030316  -3.187023 8.032750e-03
#> C-A.723         Ticam1      C-A   19.148438  4.259155   3.183361 8.086942e-03
#> C-A.724         Bcl2l1      C-A   11.128195  3.476148   3.179826 8.139602e-03
#> C-A.725         Tuba1a      C-A    2.484514  1.312964   3.178934 8.152947e-03
#> C-A.726           Nsl1      C-A   -4.306704 -2.106584  -3.170732 8.276652e-03
#> C-A.727          Nucb2      C-A   -5.927484 -2.567420  -3.170727 8.276739e-03
#> C-A.728           Cd69      C-A   -5.750617 -2.523717  -3.169450 8.296161e-03
#> C-A.729          Itpr3      C-A    8.312006  3.055197   3.166379 8.343089e-03
#> C-A.730           Orc2      C-A   -5.172324 -2.370813  -3.166149 8.346602e-03
#> C-A.731           Pdpr      C-A   37.561753  5.231192   3.159897 8.442998e-03
#> C-A.732         Il1rap      C-A  -29.396603 -4.877578  -3.159416 8.450463e-03
#> C-A.733         Gemin6      C-A   -3.146937 -1.653948  -3.158925 8.458078e-03
#> C-A.734        Tmem201      C-A   -3.268196 -1.708495  -3.158457 8.465349e-03
#> C-A.735          Uhrf1      C-A   -2.558566 -1.355335  -3.155516 8.511212e-03
#> C-A.736         Dock11      C-A   -2.736210 -1.452179  -3.152838 8.553174e-03
#> C-A.737          Larp7      C-A   -2.400079 -1.263082  -3.152301 8.561620e-03
#> C-A.738           Tia1      C-A   -2.419011 -1.274417  -3.149632 8.603697e-03
#> C-A.739          Irak3      C-A   23.694185  4.566461   3.149038 8.613088e-03
#> C-A.740        Tcp11l2      C-A    3.665098  1.873852   3.147409 8.638901e-03
#> C-A.741         Man2a1      C-A    5.957851  2.574792   3.147237 8.641627e-03
#> C-A.742  2810474O19Rik      C-A    2.184863  1.127543   3.146404 8.654866e-03
#> C-A.743          Dock5      C-A   10.072905  3.332408   3.145384 8.671094e-03
#> C-A.744           Lax1      C-A    2.748367  1.458575   3.144865 8.679373e-03
#> C-A.745        Supv3l1      C-A   -5.445137 -2.444968  -3.138875 8.775412e-03
#> C-A.746          S1pr3      C-A   -8.636122 -3.110384  -3.137701 8.794358e-03
#> C-A.747         Il18r1      C-A   36.431585  5.187118   3.133354 8.864868e-03
#> C-A.748         Dusp23      C-A   11.650695  3.542344   3.132914 8.872033e-03
#> C-A.749        Morf4l2      C-A   -2.399089 -1.262487  -3.127410 8.962215e-03
#> C-A.750          Ascc3      C-A    2.593845  1.375093   3.126237 8.981564e-03
#> C-A.751           Etv3      C-A    3.902660  1.964458   3.125620 8.991744e-03
#> C-A.752         Nudt19      C-A   -2.577567 -1.366010  -3.125515 8.993480e-03
#> C-A.753         Hivep2      C-A   17.332963  4.115446   3.123993 9.018666e-03
#> C-A.754           Ubl4      C-A   -2.105186 -1.073948  -3.119372 9.095578e-03
#> C-A.755          Snrpa      C-A   -2.052619 -1.037466  -3.116356 9.146134e-03
#> C-A.756         Tm6sf1      C-A   -3.539337 -1.823479  -3.112501 9.211161e-03
#> C-A.757           Ier5      C-A   16.896223  4.078629   3.109816 9.256725e-03
#> C-A.758          Arntl      C-A    4.133193  2.047257   3.108678 9.276099e-03
#> C-A.759          Hmces      C-A   -4.132645 -2.047065  -3.107063 9.303680e-03
#> C-A.760         Dennd3      C-A    3.865391  1.950615   3.106804 9.308105e-03
#> C-A.761         Marcks      C-A   -6.085880 -2.605466  -3.104293 9.351166e-03
#> C-A.762          Ttyh3      C-A   -2.525427 -1.336527  -3.101978 9.391040e-03
#> C-A.763         Trim16      C-A   -4.859571 -2.280829  -3.100916 9.409384e-03
#> C-A.764           Plek      C-A   -7.142224 -2.836373  -3.099795 9.428789e-03
#> C-A.765        Tmem30a      C-A    2.052240  1.037200   3.098272 9.455228e-03
#> C-A.766          Hells      C-A   -2.922273 -1.547091  -3.097914 9.461447e-03
#> C-A.767      Secisbp2l      C-A    3.054872  1.611112   3.096384 9.488095e-03
#> C-A.768         Ptger1      C-A    3.268935  1.708821   3.095024 9.511834e-03
#> C-A.769           Tet1      C-A   -4.534742 -2.181021  -3.091179 9.579304e-03
#> C-A.770         Ociad2      C-A   -4.257833 -2.090119  -3.090323 9.594400e-03
#> C-A.771           Pdcl      C-A    2.102163  1.071875   3.089157 9.614980e-03
#> C-A.772            Tec      C-A   -3.469451 -1.794707  -3.088110 9.633505e-03
#> C-A.773  0610010F05Rik      C-A   -6.063562 -2.600166  -3.087862 9.637891e-03
#> C-A.774           Gnl3      C-A   -2.938358 -1.555010  -3.087834 9.638396e-03
#> C-A.775        Ccdc186      C-A    8.889431  3.152091   3.086219 9.667055e-03
#> C-A.776           Gcsh      C-A   -3.108442 -1.636192  -3.084819 9.691954e-03
#> C-A.777            Pnp      C-A   -2.591248 -1.373647  -3.084774 9.692767e-03
#> C-A.778           Gas5      C-A   -2.058596 -1.041661  -3.083875 9.708798e-03
#> C-A.779        Plekhm1      C-A    2.195801  1.134747   3.082477 9.733780e-03
#> C-A.780  C330027C09Rik      C-A   -5.948723 -2.572580  -3.082107 9.740410e-03
#> C-A.781          Ftsj1      C-A   -5.252008 -2.392869  -3.079584 9.785682e-03
#> C-A.782          Whamm      C-A   13.240688  3.726906   3.077212 9.828445e-03
#> C-A.783     St6galnac6      C-A    3.813550  1.931135   3.076160 9.847488e-03
#> C-A.784           Ly6a      C-A    5.185001  2.374344   3.075258 9.863822e-03
#> C-A.785        Gm11223      C-A   -2.907176 -1.539618  -3.072816 9.908217e-03
#> C-A.786          Krit1      C-A    2.955219  1.563265   3.068419 9.988640e-03
#> C-A.787         Chst15      C-A   -4.351695 -2.121577  -3.067560 1.000442e-02
#> C-A.788           Brd1      C-A   -2.381032 -1.251587  -3.065307 1.004594e-02
#> C-A.789           Add1      C-A    2.277545  1.187480   3.061960 1.010796e-02
#> C-A.790           Ctsf      C-A    6.296178  2.654476   3.061722 1.011240e-02
#> C-A.791          Phtf2      C-A    2.397008  1.261235   3.061236 1.012144e-02
#> C-A.792           Mta1      C-A   -2.105098 -1.073887  -3.060740 1.013067e-02
#> C-A.793        St8sia4      C-A   -2.076698 -1.054291  -3.058818 1.016653e-02
#> C-A.794         Gnpda1      C-A    3.033666  1.601062   3.055617 1.022655e-02
#> C-A.795        Cdk2ap2      C-A    2.325433  1.217499   3.047602 1.037839e-02
#> C-A.796           Gse1      C-A    9.906919  3.308436   3.047581 1.037878e-02
#> C-A.797           Slpi      C-A   -6.068282 -2.601288  -3.045755 1.041369e-02
#> C-A.798          Ncoa2      C-A    2.189088  1.130330   3.045261 1.042315e-02
#> C-A.799           Btg2      C-A    2.934431  1.553081   3.044061 1.044618e-02
#> C-A.800          Runx3      C-A    2.051953  1.036998   3.042113 1.048367e-02
#> C-A.801           Cbx1      C-A   -2.062237 -1.044210  -3.038545 1.055269e-02
#> C-A.802         Fkbp1a      C-A   -2.648973 -1.405433  -3.033391 1.065318e-02
#> C-A.803         Pbxip1      C-A    2.112525  1.078968   3.023485 1.084903e-02
#> C-A.804         Zfp777      C-A    3.633196  1.861239   3.021038 1.089798e-02
#> C-A.805  5031425E22Rik      C-A   12.980161  3.698236   3.019732 1.092418e-02
#> C-A.806          Gna13      C-A    2.364291  1.241407   3.018981 1.093927e-02
#> C-A.807          Pold3      C-A   -2.131250 -1.091700  -3.017738 1.096432e-02
#> C-A.808       Slc39a11      C-A   -3.114252 -1.638886  -3.015457 1.101041e-02
#> C-A.809        Dennd4c      C-A   24.928661  4.639733   3.009733 1.112694e-02
#> C-A.810          Nlrc3      C-A    3.356563  1.746985   3.009369 1.113438e-02
#> C-A.811           Clta      C-A   -2.192414 -1.132520  -3.002354 1.127898e-02
#> C-A.812         Trim37      C-A   -2.559953 -1.356117  -2.996939 1.139188e-02
#> C-A.813        Rad54l2      C-A   14.970748  3.904074   2.996677 1.139737e-02
#> C-A.814         Nhlrc3      C-A   15.260214  3.931703   2.995316 1.142592e-02
#> C-A.815          Degs1      C-A    2.151233  1.105164   2.993287 1.146866e-02
#> C-A.816        Il12rb2      C-A    9.857299  3.301192   2.987041 1.160117e-02
#> C-A.817        Ankrd32      C-A   -4.859820 -2.280903  -2.978540 1.178398e-02
#> C-A.818          Hmgn3      C-A    6.812135  2.768107   2.976756 1.182273e-02
#> C-A.819         Pcyt1a      C-A   -3.271051 -1.709754  -2.976484 1.182865e-02
#> C-A.820          Prkch      C-A    2.323944  1.216575   2.973252 1.189917e-02
#> C-A.821          Strn4      C-A   -3.027915 -1.598325  -2.971131 1.194568e-02
#> C-A.822          Rufy1      C-A   -3.163207 -1.661388  -2.967166 1.203312e-02
#> C-A.823        Cyp4f16      C-A    7.047487  2.817109   2.967006 1.203665e-02
#> C-A.824           Narf      C-A    2.233185  1.159103   2.966192 1.205469e-02
#> C-A.825           Pigu      C-A   -3.651839 -1.868623  -2.965297 1.207455e-02
#> C-A.826         Angel1      C-A   22.370846  4.483548   2.961939 1.214937e-02
#> C-A.827          Fnbp4      C-A   -3.865462 -1.950641  -2.960210 1.218806e-02
#> C-A.828         Pdlim1      C-A    2.679495  1.421961   2.959476 1.220454e-02
#> C-A.829         Ikbkap      C-A   -4.118260 -2.042035  -2.958321 1.223048e-02
#> C-A.830         Zbtb24      C-A    2.680483  1.422493   2.955894 1.228522e-02
#> C-A.831         Pik3r5      C-A    6.900950  2.786795   2.954939 1.230681e-02
#> C-A.832           Pkn1      C-A    2.002295  1.001654   2.954530 1.231606e-02
#> C-A.833          P2rx7      C-A    5.007658  2.324136   2.953357 1.234268e-02
#> C-A.834          Emid1      C-A    2.551509  1.351351   2.951560 1.238352e-02
#> C-A.835          Wdhd1      C-A   -3.714624 -1.893216  -2.948801 1.244652e-02
#> C-A.836            Lbh      C-A    8.599977  3.104333   2.943859 1.256018e-02
#> C-A.837          Dnph1      C-A  -18.017872 -4.171357  -2.941776 1.260841e-02
#> C-A.838          Prrc1      C-A    2.263628  1.178637   2.941597 1.261256e-02
#> C-A.839         Tspan4      C-A   -4.376159 -2.129665  -2.938899 1.267528e-02
#> C-A.840           Fut8      C-A    2.432120  1.282214   2.937309 1.271242e-02
#> C-A.841           Pcm1      C-A   -3.379907 -1.756983  -2.935666 1.275088e-02
#> C-A.842        Tsc22d2      C-A    4.215943  2.075855   2.934508 1.277808e-02
#> C-A.843            F8a      C-A    8.332854  3.058811   2.933661 1.279799e-02
#> C-A.844          Cdan1      C-A    2.309925  1.207846   2.929318 1.290064e-02
#> C-A.845          Abhd5      C-A   12.560298  3.650799   2.928834 1.291212e-02
#> C-A.846          Asf1a      C-A   -2.177637 -1.122763  -2.927934 1.293351e-02
#> C-A.847           Pim2      C-A    5.038281  2.332932   2.927052 1.295451e-02
#> C-A.848         Knstrn      C-A   -2.483789 -1.312542  -2.920240 1.311782e-02
#> C-A.849        St6gal1      C-A    4.066926  2.023939   2.918073 1.317021e-02
#> C-A.850          Abtb1      C-A    2.660969  1.411952   2.917100 1.319380e-02
#> C-A.851         Insig2      C-A    4.415644  2.142624   2.915039 1.324391e-02
#> C-A.852           Neu1      C-A    9.127904  3.190284   2.912394 1.330850e-02
#> C-A.853           Acy1      C-A   -5.176160 -2.371882  -2.911289 1.333555e-02
#> C-A.854       Slc25a16      C-A    4.720818  2.239037   2.907281 1.343422e-02
#> C-A.855        Ccdc117      C-A   -2.630435 -1.395302  -2.906614 1.345070e-02
#> C-A.856         Pkmyt1      C-A   -3.519171 -1.815236  -2.905397 1.348085e-02
#> C-A.857         Polr1b      C-A   -3.329149 -1.735154  -2.904054 1.351417e-02
#> C-A.858          Sh2b1      C-A    2.502257  1.323230   2.902215 1.355997e-02
#> C-A.859           Rcn1      C-A   -2.051052 -1.036364  -2.898981 1.364085e-02
#> C-A.860            Zak      C-A  -16.710207 -4.062658  -2.898426 1.365478e-02
#> C-A.861           Cpa3      C-A  158.686237  7.310033   2.897621 1.367499e-02
#> C-A.862           Heca      C-A    2.611727  1.385004   2.897435 1.367967e-02
#> C-A.863          Anxa2      C-A    2.334759  1.223274   2.895954 1.371698e-02
#> C-A.864            Bsn      C-A  -11.995832 -3.584461  -2.895709 1.372316e-02
#> C-A.865        Plekha5      C-A   13.161556  3.718258   2.892843 1.379569e-02
#> C-A.866           Clpb      C-A   -3.746686 -1.905615  -2.891791 1.382238e-02
#> C-A.867          Cd163      C-A   -6.043756 -2.595445  -2.890685 1.385054e-02
#> C-A.868          G6pc3      C-A   -2.556355 -1.354088  -2.889714 1.387528e-02
#> C-A.869          Desi1      C-A   -2.410246 -1.269180  -2.889420 1.388279e-02
#> C-A.870          Cep63      C-A    3.258438  1.704181   2.885327 1.398767e-02
#> C-A.871        Trim12a      C-A    3.683201  1.880960   2.880922 1.410141e-02
#> C-A.872        Bloc1s2      C-A   -2.470767 -1.304959  -2.879085 1.414912e-02
#> C-A.873  2700049A03Rik      C-A   -4.330312 -2.114471  -2.877675 1.418586e-02
#> C-A.874        Gm16039      C-A    5.429684  2.440868   2.874482 1.426938e-02
#> C-A.875           Rpa2      C-A   -2.001515 -1.001092  -2.874135 1.427848e-02
#> C-A.876          Pskh1      C-A    2.156080  1.108411   2.873463 1.429613e-02
#> C-A.877        Herpud1      C-A   -2.151946 -1.105642  -2.872758 1.431470e-02
#> C-A.878           Dctd      C-A   -3.213372 -1.684088  -2.869981 1.438796e-02
#> C-A.879          Kdm4c      C-A    2.877821  1.524977   2.869957 1.438859e-02
#> C-A.880          Kif23      C-A   -2.279780 -1.188895  -2.867605 1.445095e-02
#> C-A.881           Relt      C-A   -3.540655 -1.824016  -2.866404 1.448290e-02
#> C-A.882           Lfng      C-A    2.593270  1.374772   2.861202 1.462205e-02
#> C-A.883         Rnf125      C-A   24.153258  4.594146   2.860426 1.464292e-02
#> C-A.884       Tnfrsf22      C-A   11.111772  3.474017   2.859875 1.465777e-02
#> C-A.885         Cgrrf1      C-A    9.464544  3.242533   2.858553 1.469344e-02
#> C-A.886           Pigl      C-A   11.254332  3.492409   2.858094 1.470583e-02
#> C-A.887           Ipo5      C-A   -2.244156 -1.166173  -2.851729 1.487892e-02
#> C-A.888          Cenpm      C-A   -3.531535 -1.820296  -2.845527 1.504951e-02
#> C-A.889           Gphn      C-A   -4.199447 -2.070200  -2.845390 1.505330e-02
#> C-A.890        Ppip5k2      C-A   -3.116354 -1.639859  -2.844837 1.506860e-02
#> C-A.891          Aplp2      C-A    2.748233  1.458504   2.840138 1.519930e-02
#> C-A.892          Senp6      C-A    2.108840  1.076450   2.839623 1.521371e-02
#> C-A.893            Ski      C-A    3.919352  1.970615   2.839366 1.522089e-02
#> C-A.894           Ice2      C-A  -13.280220 -3.731207  -2.839308 1.522252e-02
#> C-A.895  I830077J02Rik      C-A  -11.409410 -3.512152  -2.837916 1.526151e-02
#> C-A.896           Rmi1      C-A    4.100924  2.035949   2.836928 1.528925e-02
#> C-A.897           Yaf2      C-A    2.667634  1.415561   2.836504 1.530115e-02
#> C-A.898          Tada1      C-A   -3.697858 -1.886690  -2.833273 1.539229e-02
#> C-A.899           Vwa9      C-A   -2.008854 -1.006373  -2.833000 1.540001e-02
#> C-A.900          Pold2      C-A   -2.066885 -1.047458  -2.832592 1.541158e-02
#> C-A.901           Urb2      C-A   -3.590405 -1.844147  -2.830214 1.547906e-02
#> C-A.902           Junb      C-A    3.383840  1.758661   2.829960 1.548629e-02
#> C-A.903         Mettl4      C-A    3.139063  1.650334   2.829044 1.551237e-02
#> C-A.904           Ing4      C-A   -2.623354 -1.391413  -2.828987 1.551399e-02
#> C-A.905       Tmem106c      C-A   -2.566786 -1.359963  -2.828182 1.553696e-02
#> C-A.906          Ptcd1      C-A   -7.638765 -2.933339  -2.827991 1.554243e-02
#> C-A.907            Lat      C-A    3.792928  1.923312   2.823294 1.567715e-02
#> C-A.908         Zfp597      C-A   19.872155  4.312676   2.823266 1.567796e-02
#> C-A.909          Fads1      C-A   -3.613230 -1.853289  -2.821590 1.572631e-02
#> C-A.910         Acadvl      C-A    2.746340  1.457510   2.819936 1.577420e-02
#> C-A.911         Pdlim2      C-A    5.144908  2.363145   2.819670 1.578188e-02
#> C-A.912      Slfn10-ps      C-A    5.330653  2.414312   2.812396 1.599424e-02
#> C-A.913          Thada      C-A   10.528638  3.396247   2.806007 1.618307e-02
#> C-A.914           Poll      C-A   -6.247667 -2.643317  -2.802009 1.630236e-02
#> C-A.915         Zfp213      C-A    7.026556  2.812818   2.801200 1.632661e-02
#> C-A.916  9130401M01Rik      C-A    2.362281  1.240181   2.800695 1.634178e-02
#> C-A.917         Zfp319      C-A    2.948041  1.559757   2.799432 1.637971e-02
#> C-A.918         Thap11      C-A   -4.627235 -2.210150  -2.799056 1.639103e-02
#> C-A.919      Gabarapl2      C-A    2.346211  1.230333   2.796848 1.645765e-02
#> C-A.920          Cdc45      C-A   -2.093362 -1.065822  -2.795668 1.649337e-02
#> C-A.921         Mfsd7b      C-A    3.630732  1.860260   2.791955 1.660622e-02
#> C-A.922           Bub1      C-A   -3.108613 -1.636271  -2.790133 1.666190e-02
#> C-A.923          Dhx58      C-A   -4.801139 -2.263377  -2.789119 1.669295e-02
#> C-A.924        Hnrnpll      C-A    2.558523  1.355311   2.788099 1.672423e-02
#> C-A.925          Parp8      C-A   -2.031090 -1.022254  -2.786147 1.678430e-02
#> C-A.926         B3gnt5      C-A   -7.280975 -2.864132  -2.782326 1.690251e-02
#> C-A.927          Reep4      C-A   -4.369452 -2.127452  -2.782109 1.690922e-02
#> C-A.928          Ptk2b      C-A    2.029906  1.021413   2.779982 1.697539e-02
#> C-A.929          Ddx17      C-A   -2.409111 -1.268501  -2.778979 1.700669e-02
#> C-A.930          H2-Oa      C-A    4.066086  2.023641   2.773846 1.716774e-02
#> C-A.931           Mafg      C-A    4.392480  2.135036   2.773569 1.717647e-02
#> C-A.932          Kdm4a      C-A    3.280502  1.713917   2.773348 1.718345e-02
#> C-A.933          Wdr11      C-A   -8.297858 -3.052739  -2.772992 1.719467e-02
#> C-A.934          Wdr62      C-A    3.956319  1.984159   2.772547 1.720874e-02
#> C-A.935         Tex264      C-A    2.242518  1.165119   2.771749 1.723397e-02
#> C-A.936        Tmem121      C-A   -4.818487 -2.268580  -2.771202 1.725130e-02
#> C-A.937       BC055324      C-A  -12.751600 -3.672606  -2.769680 1.729955e-02
#> C-A.938  A430005L14Rik      C-A   -2.793770 -1.482213  -2.769660 1.730020e-02
#> C-A.939         Tmem9b      C-A    3.223175  1.688482   2.768793 1.732776e-02
#> C-A.940         Ttll12      C-A    2.126620  1.088562   2.768747 1.732921e-02
#> C-A.941          Rock2      C-A   -2.728951 -1.448347  -2.765877 1.742077e-02
#> C-A.942        Psmc3ip      C-A   -4.493748 -2.167919  -2.765671 1.742736e-02
#> C-A.943           Mcat      C-A   -2.140738 -1.098108  -2.765190 1.744274e-02
#> C-A.944        Aldh1b1      C-A   -3.190016 -1.673564  -2.764427 1.746718e-02
#> C-A.945          Hspa2      C-A   -9.594965 -3.262278  -2.762329 1.753459e-02
#> C-A.946          Kntc1      C-A   -4.394922 -2.135838  -2.760198 1.760331e-02
#> C-A.947         Pxylp1      C-A  -17.717332 -4.147089  -2.758507 1.765805e-02
#> C-A.948         Fam63b      C-A    2.378362  1.249968   2.756599 1.771999e-02
#> C-A.949           Pld3      C-A    5.908415  2.562771   2.755327 1.776142e-02
#> C-A.950          Padi4      C-A  -12.654503 -3.661579  -2.754577 1.778587e-02
#> C-A.951         Cdkn2c      C-A    2.698579  1.432200   2.754425 1.779084e-02
#> C-A.952         Fbxo22      C-A   -2.829715 -1.500657  -2.752685 1.784774e-02
#> C-A.953           Sil1      C-A   -2.106816 -1.075064  -2.751864 1.787465e-02
#> C-A.954         Ptdss2      C-A   -3.199262 -1.677739  -2.751350 1.789151e-02
#> C-A.955           Wdr5      C-A   -2.095715 -1.067443  -2.750930 1.790531e-02
#> C-A.956  1700025G04Rik      C-A   -3.388444 -1.760623  -2.750497 1.791955e-02
#> C-A.957          Zfp12      C-A   -4.617662 -2.207162  -2.750004 1.793578e-02
#> C-A.958        Tmem39b      C-A   -3.964840 -1.987263  -2.746669 1.804588e-02
#> C-A.959          Ssbp2      C-A   -2.724344 -1.445909  -2.746396 1.805492e-02
#> C-A.960           Rere      C-A   -5.342373 -2.417481  -2.744252 1.812612e-02
#> C-A.961          Prdx4      C-A   -2.205161 -1.140884  -2.743509 1.815084e-02
#> C-A.962      Uhrf1bp1l      C-A    2.605672  1.381655   2.742361 1.818911e-02
#> C-A.963       Rab3gap1      C-A   -2.084867 -1.059955  -2.741931 1.820345e-02
#> C-A.964           Cd72      C-A   -3.716512 -1.893949  -2.737723 1.834455e-02
#> C-A.965        Irf2bpl      C-A    9.831197  3.297367   2.736900 1.837228e-02
#> C-A.966          Ndrg1      C-A   -8.885583 -3.151466  -2.735015 1.843592e-02
#> C-A.967        Zscan29      C-A    4.782613  2.257799   2.734910 1.843946e-02
#> C-A.968          Spc24      C-A   -2.196936 -1.135493  -2.734651 1.844825e-02
#> C-A.969         Pi4k2b      C-A   -5.454829 -2.447534  -2.734381 1.845738e-02
#> C-A.970        Fam117a      C-A   -2.031189 -1.022324  -2.734185 1.846400e-02
#> C-A.971          Gria3      C-A  -16.416604 -4.037084  -2.733275 1.849488e-02
#> C-A.972        Tnfaip3      C-A   21.407654  4.420055   2.732488 1.852158e-02
#> C-A.973          S1pr4      C-A    7.369424  2.881552   2.731915 1.854105e-02
#> C-A.974          Unc50      C-A   -2.197159 -1.135639  -2.731884 1.854212e-02
#> C-A.975          Zfp68      C-A   -3.637952 -1.863126  -2.730128 1.860193e-02
#> C-A.976          Ssbp1      C-A   -2.080229 -1.056743  -2.729243 1.863218e-02
#> C-A.977          Alas1      C-A   -2.761672 -1.465542  -2.727415 1.869475e-02
#> C-A.978         Bcap29      C-A   -2.491992 -1.317300  -2.726636 1.872149e-02
#> C-A.979           Med1      C-A   -3.671327 -1.876301  -2.725926 1.874588e-02
#> C-A.980         Katnb1      C-A   -3.016257 -1.592760  -2.724556 1.879305e-02
#> C-A.981         Sorbs3      C-A   -4.712493 -2.236491  -2.723206 1.883962e-02
#> C-A.982         Nudcd2      C-A   -2.857076 -1.514540  -2.722248 1.887275e-02
#> C-A.983           Hmbs      C-A   -2.562349 -1.357467  -2.714793 1.913254e-02
#> C-A.984          Vps16      C-A   -2.281957 -1.190272  -2.713767 1.916857e-02
#> C-A.985          Limk1      C-A   -3.431967 -1.779036  -2.712514 1.921268e-02
#> C-A.986       Slc25a13      C-A  -11.918468 -3.575127  -2.710619 1.927956e-02
#> C-A.987          Anxa5      C-A    2.121058  1.084784   2.709616 1.931503e-02
#> C-A.988           Nrp1      C-A    5.980828  2.580345   2.708480 1.935532e-02
#> C-A.989       Trp53i11      C-A   -2.374294 -1.247498  -2.707669 1.938413e-02
#> C-A.990           Pck2      C-A   -2.134717 -1.094045  -2.707484 1.939069e-02
#> C-A.991         Dcaf15      C-A    5.069201  2.341758   2.701794 1.959399e-02
#> C-A.992         Ppp3cc      C-A   -3.147794 -1.654341  -2.695582 1.981839e-02
#> C-A.993          Rab31      C-A   -7.161884 -2.840339  -2.695162 1.983366e-02
#> C-A.994         Lysmd4      C-A    8.165974  3.029625   2.693855 1.988121e-02
#> C-A.995          Syvn1      C-A    2.388784  1.256276   2.692251 1.993974e-02
#> C-A.996        Tmem109      C-A   -2.012285 -1.008835  -2.690705 1.999633e-02
#> C-A.997  2010111I01Rik      C-A   -6.676775 -2.739151  -2.687556 2.011205e-02
#> C-A.998          Nfkb2      C-A    8.226515  3.040281   2.685601 2.018424e-02
#> C-A.999          C2cd3      C-A   -2.869897 -1.520999  -2.684228 2.023508e-02
#> C-A.1000       Ubash3b      C-A    7.144710  2.836875   2.682276 2.030755e-02
#> C-A.1001         Prdm4      C-A   -3.055500 -1.611409  -2.682232 2.030918e-02
#> C-A.1002        Fchsd1      C-A    7.585282  2.923203   2.677941 2.046948e-02
#> C-A.1003        Card11      C-A    6.186227  2.629060   2.677384 2.049040e-02
#> C-A.1004          Lrmp      C-A   -2.152028 -1.105697  -2.676517 2.052297e-02
#> C-A.1005          Sla2      C-A    5.437720  2.443002   2.676488 2.052403e-02
#> C-A.1006         Neil1      C-A    5.869500  2.553238   2.676265 2.053244e-02
#> C-A.1007         Ift22      C-A   -4.814177 -2.267289  -2.675072 2.057734e-02
#> C-A.1008          Nfu1      C-A   -2.088995 -1.062809  -2.674323 2.060559e-02
#> C-A.1009          Lyz2      C-A   -2.246165 -1.167464  -2.672856 2.066104e-02
#> C-A.1010        Zfp113      C-A   -3.398698 -1.764982  -2.672558 2.067233e-02
#> C-A.1011         Gins1      C-A   -3.768906 -1.914146  -2.670099 2.076563e-02
#> C-A.1012         Nxpe3      C-A   -3.723855 -1.896797  -2.670032 2.076816e-02
#> C-A.1013         Mrpl1      C-A   -5.236470 -2.388594  -2.669885 2.077378e-02
#> C-A.1014          Dok1      C-A    4.381710  2.131494   2.668349 2.083230e-02
#> C-A.1015      Suv420h1      C-A    3.570397  1.836085   2.667167 2.087743e-02
#> C-A.1016         Rchy1      C-A    2.505973  1.325371   2.665712 2.093314e-02
#> C-A.1017           Hk3      C-A    3.242963  1.697313   2.665302 2.094885e-02
#> C-A.1018       Angptl2      C-A    7.292760  2.866465   2.665088 2.095709e-02
#> C-A.1019         Abcd3      C-A    2.398695  1.262250   2.664874 2.096529e-02
#> C-A.1020         Mcm10      C-A   -2.722458 -1.444910  -2.663954 2.100062e-02
#> C-A.1021         Rc3h2      C-A    5.735491  2.519917   2.662560 2.105433e-02
#> C-A.1022           Ank      C-A    8.902538  3.154217   2.662052 2.107389e-02
#> C-A.1023        Golgb1      C-A    4.438493  2.150070   2.656328 2.129592e-02
#> C-A.1024         Grina      C-A    3.748811  1.906433   2.656182 2.130163e-02
#> C-A.1025       Ppapdc2      C-A    5.541181  2.470193   2.655065 2.134521e-02
#> C-A.1026        Zfp330      C-A   -2.514787 -1.330436  -2.654030 2.138571e-02
#> C-A.1027         Capn1      C-A    2.263367  1.178471   2.653223 2.141733e-02
#> C-A.1028        Dusp12      C-A    9.617120  3.265605   2.648426 2.160622e-02
#> C-A.1029         Gtse1      C-A   -2.465885 -1.302106  -2.648385 2.160782e-02
#> C-A.1030        Calcrl      C-A   -2.101023 -1.071092  -2.647346 2.164894e-02
#> C-A.1031         Fbxo5      C-A   -2.912315 -1.542166  -2.644971 2.174327e-02
#> C-A.1032         Supt6      C-A   -3.427723 -1.777251  -2.643737 2.179242e-02
#> C-A.1033        Zmynd8      C-A   -2.171213 -1.118501  -2.642868 2.182708e-02
#> C-A.1034         Appl1      C-A    2.388807  1.256290   2.641625 2.187678e-02
#> C-A.1035         Pcgf6      C-A   -6.699824 -2.744123  -2.639199 2.197411e-02
#> C-A.1036        Klhl36      C-A    3.549227  1.827505   2.636553 2.208076e-02
#> C-A.1037         Trub2      C-A  -11.384450 -3.508993  -2.636225 2.209403e-02
#> C-A.1038         Siva1      C-A   -2.063935 -1.045398  -2.635758 2.211291e-02
#> C-A.1039         Enpp5      C-A    6.169725  2.625206   2.631345 2.229214e-02
#> C-A.1040         Mier1      C-A    2.041594  1.029696   2.628502 2.240839e-02
#> C-A.1041          Ccnf      C-A   -2.142208 -1.099099  -2.626580 2.248728e-02
#> C-A.1042         Haus1      C-A   -2.284991 -1.192188  -2.623325 2.262153e-02
#> C-A.1043         Irak2      C-A   44.836068  5.486588   2.621742 2.268711e-02
#> C-A.1044        Nufip1      C-A   -3.593924 -1.845560  -2.620200 2.275121e-02
#> C-A.1045         Usp34      C-A    6.924674  2.791746   2.619330 2.278741e-02
#> C-A.1046 1810022K09Rik      C-A   -2.381129 -1.251646  -2.618074 2.283981e-02
#> C-A.1047       Ppp1r11      C-A    2.416921  1.273170   2.617221 2.287543e-02
#> C-A.1048 1700021K19Rik      C-A   -3.426809 -1.776866  -2.615023 2.296756e-02
#> C-A.1049         Kif15      C-A   -2.242075 -1.164835  -2.613258 2.304179e-02
#> C-A.1050          Ppan      C-A   -2.052038 -1.037057  -2.611775 2.310433e-02
#> C-A.1051         Scmh1      C-A   -4.175073 -2.061801  -2.611647 2.310974e-02
#> C-A.1052         Capn5      C-A   -7.732396 -2.950915  -2.611081 2.313365e-02
#> C-A.1053        Heatr3      C-A   -2.145681 -1.101435  -2.610628 2.315281e-02
#> C-A.1054         Hyal2      C-A   -6.084804 -2.605211  -2.609958 2.318122e-02
#> C-A.1055         Exoc8      C-A    2.619729  1.389418   2.608165 2.325731e-02
#> C-A.1056         Rab4a      C-A    6.375080  2.672443   2.607351 2.329193e-02
#> C-A.1057         Bnip3      C-A   -4.527826 -2.178819  -2.605415 2.337451e-02
#> C-A.1058       Zcchc10      C-A    5.220243  2.384117   2.604532 2.341224e-02
#> C-A.1059          Plk4      C-A   -2.459339 -1.298271  -2.603681 2.344868e-02
#> C-A.1060         C1qbp      C-A   -2.083274 -1.058852  -2.601296 2.355110e-02
#> C-A.1061        Kif20b      C-A   -2.526094 -1.336908  -2.597341 2.372191e-02
#> C-A.1062         Clcn7      C-A    3.554944  1.829827   2.596110 2.377531e-02
#> C-A.1063       Wbscr16      C-A   -4.311249 -2.108106  -2.595140 2.381750e-02
#> C-A.1064          Eny2      C-A   -2.147607 -1.102730  -2.590249 2.403124e-02
#> C-A.1065           Gem      C-A    7.349566  2.877659   2.590096 2.403795e-02
#> C-A.1066         Deaf1      C-A    3.003298  1.586548   2.589831 2.404958e-02
#> C-A.1067         Rbms1      C-A    2.247135  1.168087   2.587835 2.413742e-02
#> C-A.1068         Dhx32      C-A    2.319793  1.213996   2.587488 2.415273e-02
#> C-A.1069          Bptf      C-A    2.263342  1.178455   2.585408 2.424462e-02
#> C-A.1070         Top2a      C-A   -2.602965 -1.380156  -2.582571 2.437055e-02
#> C-A.1071       Slc43a2      C-A    2.054238  1.038603   2.581653 2.441147e-02
#> C-A.1072         Mfsd4      C-A    2.335836  1.223939   2.579076 2.452656e-02
#> C-A.1073         Ppm1a      C-A    2.954166  1.562751   2.578609 2.454752e-02
#> C-A.1074         Ahnak      C-A   11.515678  3.525527   2.578431 2.455546e-02
#> C-A.1075        Arfip1      C-A    5.183852  2.374025   2.576884 2.462493e-02
#> C-A.1076       Hnrnph1      C-A   -2.024526 -1.017584  -2.574156 2.474787e-02
#> C-A.1077       Map3k14      C-A    4.154215  2.054576   2.571992 2.484581e-02
#> C-A.1078        Cdc14a      C-A    3.237208  1.694750   2.570401 2.491805e-02
#> C-A.1079         Acad8      C-A   -3.109398 -1.636635  -2.570260 2.492448e-02
#> C-A.1080       Rps6kb2      C-A    4.012935  2.004658   2.569603 2.495436e-02
#> C-A.1081        Rab3ip      C-A  -18.867899 -4.237862  -2.569024 2.498077e-02
#> C-A.1082      Phospho2      C-A    2.349252  1.232202   2.568893 2.498675e-02
#> C-A.1083         Wdr44      C-A   13.728563  3.779109   2.568392 2.500958e-02
#> C-A.1084           Dut      C-A   -2.033527 -1.023984  -2.568281 2.501463e-02
#> C-A.1085 1110059G10Rik      C-A    3.596160  1.846457   2.567103 2.506850e-02
#> C-A.1086          Pkib      C-A   -4.550984 -2.186178  -2.566306 2.510498e-02
#> C-A.1087        Kif18a      C-A   -6.981280 -2.803492  -2.565348 2.514891e-02
#> C-A.1088       Thumpd3      C-A   -3.342439 -1.740901  -2.563495 2.523405e-02
#> C-A.1089        Mms22l      C-A   -3.212989 -1.683916  -2.563395 2.523867e-02
#> C-A.1090       Galnt11      C-A    2.705560  1.435927   2.563330 2.524168e-02
#> C-A.1091         Myo5a      C-A  -16.044485 -4.004006  -2.562581 2.527616e-02
#> C-A.1092        Tsen34      C-A   -2.137536 -1.095948  -2.561971 2.530434e-02
#> C-A.1093       Arl6ip6      C-A   -4.169700 -2.059944  -2.557095 2.553042e-02
#> C-A.1094          Rgs2      C-A   -4.909687 -2.295631  -2.555279 2.561513e-02
#> C-A.1095         Dcaf5      C-A    2.239844  1.163398   2.555101 2.562343e-02
#> C-A.1096          Nop2      C-A   -2.070566 -1.050025  -2.554255 2.566301e-02
#> C-A.1097     Hist1h2ae      C-A   -3.805904 -1.928239  -2.552323 2.575358e-02
#> C-A.1098        Insig1      C-A   -6.853497 -2.776840  -2.551264 2.580337e-02
#> C-A.1099        Sft2d1      C-A   -3.372546 -1.753838  -2.551096 2.581130e-02
#> C-A.1100         Prtn3      C-A   -3.783108 -1.919572  -2.549184 2.590142e-02
#> C-A.1101         Dyrk3      C-A   -4.960993 -2.310629  -2.548498 2.593384e-02
#> C-A.1102         Smyd5      C-A   -3.347814 -1.743220  -2.543624 2.616535e-02
#> C-A.1103         Zadh2      C-A   12.613558  3.656903   2.543282 2.618166e-02
#> C-A.1104         Tnip2      C-A   11.961885  3.580373   2.540408 2.631916e-02
#> C-A.1105       St3gal2      C-A   -2.907747 -1.539902  -2.539519 2.636184e-02
#> C-A.1106          Mtf1      C-A   -4.233347 -2.081799  -2.539216 2.637641e-02
#> C-A.1107         Mier3      C-A   -2.000520 -1.000375  -2.538536 2.640914e-02
#> C-A.1108          Prcp      C-A   -4.023265 -2.008367  -2.536420 2.651117e-02
#> C-A.1109        Wdr45b      C-A   -2.586265 -1.370870  -2.534347 2.661149e-02
#> C-A.1110          Ppcs      C-A    2.963285  1.567197   2.532875 2.668299e-02
#> C-A.1111        Bckdha      C-A   -2.027156 -1.019457  -2.531765 2.673701e-02
#> C-A.1112          Nqo2      C-A    7.270930  2.862140   2.530888 2.677976e-02
#> C-A.1113        Fam63a      C-A    6.320074  2.659942   2.527118 2.696434e-02
#> C-A.1114         Orai2      C-A    2.018704  1.013429   2.526551 2.699218e-02
#> C-A.1115         Itih5      C-A   -3.112343 -1.638001  -2.520577 2.728749e-02
#> C-A.1116          Ccr7      C-A    2.623172  1.391313   2.520307 2.730090e-02
#> C-A.1117        Kdelc1      C-A   -5.538961 -2.469615  -2.520038 2.731427e-02
#> C-A.1118           Mvd      C-A    9.454615  3.241019   2.519382 2.734691e-02
#> C-A.1119        Ift122      C-A   -3.209967 -1.682559  -2.516510 2.749032e-02
#> C-A.1120        Cops7a      C-A   -2.116310 -1.081551  -2.516250 2.750331e-02
#> C-A.1121         Ttll5      C-A    7.996085  2.999294   2.510098 2.781307e-02
#> C-A.1122          Ulk3      C-A    7.483999  2.903809   2.509831 2.782660e-02
#> C-A.1123         Arl16      C-A   13.403588  3.744547   2.508001 2.791941e-02
#> C-A.1124        Zc3h18      C-A   -2.945886 -1.558702  -2.501272 2.826338e-02
#> C-A.1125        Gpr160      C-A   -3.047617 -1.607682  -2.501134 2.827048e-02
#> C-A.1126       Kansl1l      C-A   -6.020730 -2.589938  -2.499456 2.835693e-02
#> C-A.1127         Acap2      C-A   -2.208427 -1.143019  -2.498618 2.840020e-02
#> C-A.1128        Dyrk1b      C-A    6.589616  2.720194   2.497822 2.844134e-02
#> C-A.1129         Dhodh      C-A   -2.706577 -1.436470  -2.495913 2.854030e-02
#> C-A.1130         Galk1      C-A   -2.637276 -1.399049  -2.495670 2.855288e-02
#> C-A.1131       Gpbp1l1      C-A    2.377491  1.249440   2.495477 2.856290e-02
#> C-A.1132          Dexi      C-A    3.526554  1.818259   2.494159 2.863148e-02
#> C-A.1133         Mmgt1      C-A   -3.366591 -1.751288  -2.492748 2.870502e-02
#> C-A.1134       Pitpnm1      C-A   12.322960  3.623277   2.491533 2.876851e-02
#> C-A.1135         Uqcc1      C-A   -2.126287 -1.088336  -2.491033 2.879469e-02
#> C-A.1136         Spag5      C-A   -8.744934 -3.128448  -2.490161 2.884038e-02
#> C-A.1137         Il1r1      C-A   -2.651034 -1.406555  -2.488667 2.891885e-02
#> C-A.1138          Srrd      C-A    4.341164  2.118082   2.487640 2.897290e-02
#> C-A.1139         Bpnt1      C-A    3.349187  1.743811   2.487089 2.900195e-02
#> C-A.1140          Rab9      C-A   -4.479041 -2.163190  -2.485894 2.906503e-02
#> C-A.1141       Tmem115      C-A    2.175665  1.121457   2.484518 2.913784e-02
#> C-A.1142        Trim59      C-A   -2.081669 -1.057741  -2.481401 2.930338e-02
#> C-A.1143         Dctn4      C-A   -2.398383 -1.262062  -2.481170 2.931568e-02
#> C-A.1144        Zfp952      C-A   16.197367  4.017687   2.479740 2.939196e-02
#> C-A.1145          Glrx      C-A   -2.070609 -1.050055  -2.479714 2.939335e-02
#> C-A.1146        Mfsd12      C-A   -2.869199 -1.520648  -2.478292 2.946942e-02
#> C-A.1147       Bcl2l11      C-A    5.194229  2.376910   2.475991 2.959289e-02
#> C-A.1148         Taok2      C-A    2.711579  1.439133   2.474845 2.965459e-02
#> C-A.1149        Lrrc8b      C-A    3.913514  1.968465   2.474749 2.965973e-02
#> C-A.1150        Dicer1      C-A    2.421095  1.275659   2.474713 2.966169e-02
#> C-A.1151       Pglyrp2      C-A    2.739909  1.454128   2.474354 2.968103e-02
#> C-A.1152 1500011B03Rik      C-A    2.886903  1.529523   2.474207 2.968899e-02
#> C-A.1153       Tubgcp3      C-A   -3.187267 -1.672320  -2.473791 2.971141e-02
#> C-A.1154        Rfxank      C-A    3.382914  1.758267   2.471925 2.981232e-02
#> C-A.1155        Gatsl2      C-A    3.691554  1.884228   2.471362 2.984280e-02
#> C-A.1156         Socs5      C-A    4.847411  2.277214   2.469708 2.993263e-02
#> C-A.1157          Ipmk      C-A    4.057284  2.020514   2.466916 3.008477e-02
#> C-A.1158       Tubgcp4      C-A    8.098592  3.017671   2.465198 3.017882e-02
#> C-A.1159          Ska3      C-A   -3.448571 -1.785999  -2.463633 3.026470e-02
#> C-A.1160         Myo9b      C-A    2.287642  1.193861   2.462344 3.033559e-02
#> C-A.1161       Plekha3      C-A    3.744472  1.904762   2.461067 3.040605e-02
#> C-A.1162         Meis1      C-A   -5.389197 -2.430070  -2.458200 3.056468e-02
#> C-A.1163         Ntmt1      C-A   -3.865615 -1.950698  -2.455442 3.071812e-02
#> C-A.1164      Ankrd13c      C-A    3.813316  1.931046   2.455426 3.071898e-02
#> C-A.1165        Il20ra      C-A    5.831980  2.543986   2.454010 3.079803e-02
#> C-A.1166       Mettl14      C-A   -3.768939 -1.914159  -2.450740 3.098132e-02
#> C-A.1167         Abca3      C-A   -2.083514 -1.059019  -2.449025 3.107788e-02
#> C-A.1168         Gnptg      C-A    3.184882  1.671240   2.447696 3.115294e-02
#> C-A.1169         Stk24      C-A    5.786495  2.532690   2.447218 3.117997e-02
#> C-A.1170       Tmem63a      C-A    2.391608  1.257981   2.445895 3.125486e-02
#> C-A.1171         Kcnk5      C-A   -3.479885 -1.799040  -2.444917 3.131038e-02
#> C-A.1172         Usp20      C-A    7.985567  2.997395   2.444739 3.132046e-02
#> C-A.1173          Nom1      C-A   -3.079413 -1.622655  -2.444213 3.135038e-02
#> C-A.1174          Idh1      C-A   -2.128139 -1.089592  -2.442762 3.143301e-02
#> C-A.1175         Phf20      C-A    2.294357  1.198090   2.442694 3.143686e-02
#> C-A.1176         Paqr7      C-A    3.437282  1.781268   2.441984 3.147740e-02
#> C-A.1177          Prkx      C-A    6.379556  2.673456   2.441749 3.149082e-02
#> C-A.1178         Kif2a      C-A    7.675138  2.940193   2.441642 3.149692e-02
#> C-A.1179         Nup43      C-A   -3.300634 -1.722743  -2.439257 3.163343e-02
#> C-A.1180 1600012H06Rik      C-A    5.677804  2.505333   2.437909 3.171084e-02
#> C-A.1181         Pex13      C-A   -2.878860 -1.525498  -2.437884 3.171229e-02
#> C-A.1182           Ttk      C-A   -5.307324 -2.407985  -2.437550 3.173149e-02
#> C-A.1183         Scrn3      C-A   29.017030  4.858828   2.437185 3.175252e-02
#> C-A.1184         Mapk6      C-A   -2.152075 -1.105728  -2.436005 3.182054e-02
#> C-A.1185          Idnk      C-A    3.783415  1.919689   2.434962 3.188077e-02
#> C-A.1186         Ap1s2      C-A    3.308077  1.725993   2.430705 3.212770e-02
#> C-A.1187       H2-DMb2      C-A    4.503405  2.171016   2.430543 3.213713e-02
#> C-A.1188          Mlkl      C-A   -6.732279 -2.751095  -2.430525 3.213817e-02
#> C-A.1189         Plcl2      C-A   -4.057494 -2.020589  -2.428554 3.225322e-02
#> C-A.1190         Troap      C-A   -5.353582 -2.420505  -2.426473 3.237505e-02
#> C-A.1191         Chpf2      C-A    2.399606  1.262797   2.424369 3.249872e-02
#> C-A.1192        Atp2b1      C-A    2.600598  1.378844   2.417641 3.289719e-02
#> C-A.1193        Abhd11      C-A   -4.215286 -2.075631  -2.415492 3.302550e-02
#> C-A.1194 2810006K23Rik      C-A   -4.005300 -2.001910  -2.413545 3.314211e-02
#> C-A.1195         Ly6c2      C-A    2.798110  1.484453   2.412975 3.317632e-02
#> C-A.1196        Dusp11      C-A    2.815465  1.493373   2.412975 3.317634e-02
#> C-A.1197       Depdc1b      C-A    2.677643  1.420964   2.412843 3.318430e-02
#> C-A.1198          Hagh      C-A    3.767218  1.913500   2.411277 3.327848e-02
#> C-A.1199        Lrrc42      C-A    3.529890  1.819623   2.411120 3.328793e-02
#> C-A.1200          Rit1      C-A    3.508359  1.810796   2.408836 3.342589e-02
#> C-A.1201         Ttc19      C-A    9.879998  3.304511   2.407395 3.351314e-02
#> C-A.1202         Ddias      C-A    2.793704  1.482179   2.406627 3.355979e-02
#> C-A.1203          Mfn1      C-A    6.803711  2.766322   2.405988 3.359862e-02
#> C-A.1204         Sirt2      C-A    2.242263  1.164956   2.405243 3.364394e-02
#> C-A.1205       Elmsan1      C-A    6.009065  2.587141   2.402351 3.382049e-02
#> C-A.1206       Dnajc19      C-A    7.516504  2.910062   2.399758 3.397953e-02
#> C-A.1207         Fem1c      C-A    2.076902  1.054433   2.398779 3.403975e-02
#> C-A.1208           Pcx      C-A   14.316900  3.839647   2.397848 3.409714e-02
#> C-A.1209         Yipf2      C-A    9.324623  3.221045   2.396573 3.417583e-02
#> C-A.1210         Zmym3      C-A   -2.329516 -1.220030  -2.395972 3.421300e-02
#> C-A.1211         Birc3      C-A    2.873541  1.522830   2.395936 3.421521e-02
#> C-A.1212       Rapgef6      C-A    2.017461  1.012541   2.389077 3.464226e-02
#> C-A.1213          Hcst      C-A    4.163782  2.057894   2.388873 3.465500e-02
#> C-A.1214       Sertad1      C-A   14.094210  3.817031   2.386893 3.477927e-02
#> C-A.1215        Zfp608      C-A   -3.742349 -1.903944  -2.384856 3.490755e-02
#> C-A.1216        Rnf157      C-A   -2.208057 -1.142777  -2.382122 3.508049e-02
#> C-A.1217        Zfp629      C-A   -3.298712 -1.721903  -2.375958 3.547333e-02
#> C-A.1218        Nos1ap      C-A   -7.310660 -2.870002  -2.374409 3.557271e-02
#> C-A.1219        Man1a2      C-A    2.208031  1.142761   2.372863 3.567215e-02
#> C-A.1220        Unc45a      C-A    2.234078  1.159680   2.369077 3.591688e-02
#> C-A.1221         Ddx56      C-A   -2.052979 -1.037719  -2.367255 3.603520e-02
#> C-A.1222        Kif20a      C-A   -2.697628 -1.431692  -2.366767 3.606699e-02
#> C-A.1223         Oasl2      C-A   -3.724880 -1.897194  -2.366374 3.609257e-02
#> C-A.1224           Lpp      C-A    4.228406  2.080114   2.366036 3.611460e-02
#> C-A.1225         Cbll1      C-A    3.847484  1.943915   2.365513 3.614875e-02
#> C-A.1226          Lnx2      C-A    4.397328  2.136627   2.365099 3.617572e-02
#> C-A.1227      Mir142hg      C-A    9.090434  3.184349   2.364977 3.618371e-02
#> C-A.1228          Ralb      C-A    4.556671  2.187980   2.363734 3.626496e-02
#> C-A.1229       Ccdc28b      C-A   -4.713905 -2.236923  -2.362419 3.635113e-02
#> C-A.1230      Ankrd13b      C-A    9.316433  3.219778   2.361856 3.638810e-02
#> C-A.1231 4933426M11Rik      C-A    3.133477  1.647765   2.360954 3.644735e-02
#> C-A.1232         Fnip1      C-A   12.022481  3.587663   2.359198 3.656305e-02
#> C-A.1233      Gpatch11      C-A    8.875936  3.149899   2.359192 3.656344e-02
#> C-A.1234        Smurf1      C-A    9.939784  3.313215   2.359121 3.656815e-02
#> C-A.1235        Ifnar1      C-A   -2.388565 -1.256144  -2.357594 3.666899e-02
#> C-A.1236          Rcl1      C-A   -2.711605 -1.439147  -2.356345 3.675172e-02
#> C-A.1237        Zfp599      C-A   13.964659  3.803708   2.347457 3.734554e-02
#> C-A.1238        Fignl1      C-A   -2.574725 -1.364418  -2.346416 3.741568e-02
#> C-A.1239          Lars      C-A   -2.171343 -1.118588  -2.344359 3.755469e-02
#> C-A.1240       Zcchc18      C-A    5.161950  2.367916   2.342897 3.765379e-02
#> C-A.1241        Lztfl1      C-A    4.402256  2.138243   2.342456 3.768368e-02
#> C-A.1242         Mctp2      C-A   -3.772053 -1.915350  -2.340433 3.782132e-02
#> C-A.1243        Ift172      C-A   -4.501313 -2.170346  -2.337857 3.799732e-02
#> C-A.1244        Rab33b      C-A    4.581603  2.195853   2.336867 3.806510e-02
#> C-A.1245         Vamp4      C-A    2.751959  1.460459   2.335436 3.816336e-02
#> C-A.1246         Ifrd2      C-A   -2.472771 -1.306128  -2.332995 3.833154e-02
#> C-A.1247          Haao      C-A   -2.722495 -1.444930  -2.330869 3.847859e-02
#> C-A.1248          Cinp      C-A   -2.814400 -1.492827  -2.327536 3.871014e-02
#> C-A.1249       Slc43a3      C-A   -6.118681 -2.613221  -2.326015 3.881628e-02
#> C-A.1250          Gga1      C-A    2.414015  1.271434   2.324742 3.890533e-02
#> C-A.1251         Pde12      C-A   -2.154342 -1.107248  -2.321581 3.912731e-02
#> C-A.1252         Bod1l      C-A    2.491025  1.316740   2.319436 3.927860e-02
#> C-A.1253         Ttpal      C-A    2.380295  1.251141   2.318435 3.934939e-02
#> C-A.1254         Shprh      C-A    2.094879  1.066867   2.316445 3.949052e-02
#> C-A.1255        Smurf2      C-A    5.060942  2.339406   2.313589 3.969389e-02
#> C-A.1256         Pcf11      C-A    2.963315  1.567212   2.308161 4.008323e-02
#> C-A.1257         Tulp4      C-A   14.330649  3.841032   2.306749 4.018512e-02
#> C-A.1258       Slc35f2      C-A    4.413015  2.141764   2.306374 4.021214e-02
#> C-A.1259         Prim2      C-A   -2.681654 -1.423123  -2.306259 4.022051e-02
#> C-A.1260      Ifi27l2a      C-A  -13.866017 -3.793481  -2.305097 4.030457e-02
#> C-A.1261      Tor1aip1      C-A    2.077390  1.054772   2.304318 4.036105e-02
#> C-A.1262        Gtf3c4      C-A   -2.301724 -1.202715  -2.303950 4.038775e-02
#> C-A.1263        Rad54l      C-A   -3.064688 -1.615740  -2.303257 4.043806e-02
#> C-A.1264       Osbpl1a      C-A   -7.197703 -2.847537  -2.303257 4.043807e-02
#> C-A.1265           Emd      C-A   -2.984248 -1.577368  -2.302863 4.046673e-02
#> C-A.1266        Plscr1      C-A    7.583102  2.922788   2.302652 4.048208e-02
#> C-A.1267        Samd9l      C-A   15.905845  3.991485   2.298792 4.076380e-02
#> C-A.1268      Tmem126b      C-A    2.621165  1.390208   2.296882 4.090389e-02
#> C-A.1269         Ap5z1      C-A    6.162110  2.623424   2.293853 4.112697e-02
#> C-A.1270          Gsap      C-A   10.634415  3.410669   2.292315 4.124072e-02
#> C-A.1271        Slc9a1      C-A    5.089222  2.347445   2.291031 4.133592e-02
#> C-A.1272         Mocs1      C-A    3.433641  1.779739   2.290904 4.134537e-02
#> C-A.1273         Cdk12      C-A    9.552785  3.255921   2.289871 4.142207e-02
#> C-A.1274         Bace1      C-A   -2.663784 -1.413477  -2.288487 4.152517e-02
#> C-A.1275        Lpcat2      C-A    3.345859  1.742376   2.286956 4.163937e-02
#> C-A.1276         Aggf1      C-A   -2.891796 -1.531966  -2.286523 4.167179e-02
#> C-A.1277         Milr1      C-A   -4.608940 -2.204435  -2.284155 4.184925e-02
#> C-A.1278         Extl3      C-A   -2.964812 -1.567941  -2.281922 4.201727e-02
#> C-A.1279         Tmppe      C-A    2.725290  1.446410   2.280667 4.211191e-02
#> C-A.1280        Rsph3b      C-A    5.384087  2.428702   2.280395 4.213251e-02
#> C-A.1281          Nck2      C-A    3.742092  1.903845   2.279736 4.218232e-02
#> C-A.1282       Gpatch8      C-A    2.441826  1.287960   2.279047 4.223446e-02
#> C-A.1283         Nat10      C-A   -2.924485 -1.548183  -2.279008 4.223742e-02
#> C-A.1284        Gatad1      C-A    3.123518  1.643172   2.278811 4.225237e-02
#> C-A.1285         Arrb1      C-A    7.460254  2.899225   2.278589 4.226923e-02
#> C-A.1286       Tbc1d19      C-A    3.659142  1.871505   2.278369 4.228585e-02
#> C-A.1287          Pdxk      C-A   -3.776977 -1.917232  -2.276854 4.240096e-02
#> C-A.1288        Ccdc55      C-A   -2.797466 -1.484120  -2.275750 4.248495e-02
#> C-A.1289          Scly      C-A    3.892027  1.960522   2.275702 4.248863e-02
#> C-A.1290         Mov10      C-A    2.055545  1.039521   2.275455 4.250740e-02
#> C-A.1291        Snapc3      C-A    5.599229  2.485228   2.274939 4.254677e-02
#> C-A.1292 8430419L09Rik      C-A   -3.569412 -1.835687  -2.272695 4.271834e-02
#> C-A.1293         Maml3      C-A   -2.961669 -1.566410  -2.272018 4.277020e-02
#> C-A.1294          Xrn1      C-A   -3.140844 -1.651152  -2.270557 4.288233e-02
#> C-A.1295           Evl      C-A    2.543956  1.347074   2.270046 4.292163e-02
#> C-A.1296        Rassf1      C-A    4.112925  2.040165   2.269764 4.294329e-02
#> C-A.1297          Ass1      C-A   -3.543084 -1.825006  -2.267821 4.309314e-02
#> C-A.1298         Lrch1      C-A    2.364076  1.241277   2.267342 4.313015e-02
#> C-A.1299          Mmab      C-A    2.456565  1.296642   2.267024 4.315471e-02
#> C-A.1300         Gstt2      C-A    3.003960  1.586866   2.260470 4.366439e-02
#> C-A.1301         Snx18      C-A   -2.264376 -1.179114  -2.258550 4.381478e-02
#> C-A.1302        Tmem80      C-A    9.570935  3.258660   2.253587 4.420584e-02
#> C-A.1303         Nsun4      C-A   -4.577396 -2.194527  -2.253057 4.424774e-02
#> C-A.1304        Baiap3      C-A    6.290036  2.653068   2.251865 4.434222e-02
#> C-A.1305         Trak1      C-A   -2.919659 -1.545800  -2.251759 4.435065e-02
#> C-A.1306         Ccsap      C-A    7.679608  2.941033   2.250928 4.441669e-02
#> C-A.1307        Ccp110      C-A   -2.161386 -1.111957  -2.250547 4.444695e-02
#> C-A.1308          Bdh1      C-A    4.770705  2.254202   2.249620 4.452075e-02
#> C-A.1309        Fam64a      C-A   -2.777072 -1.473565  -2.249614 4.452120e-02
#> C-A.1310         Ddx49      C-A   -2.864842 -1.518455  -2.243865 4.498144e-02
#> C-A.1311         Pnpt1      C-A   -2.917825 -1.544893  -2.242947 4.505538e-02
#> C-A.1312       Gm10384      C-A  -11.213551 -3.487171  -2.242331 4.510500e-02
#> C-A.1313         Gpr18      C-A    4.029988  2.010775   2.240585 4.524605e-02
#> C-A.1314          Ganc      C-A    3.455788  1.789015   2.237127 4.552661e-02
#> C-A.1315         Nrde2      C-A    2.293374  1.197472   2.236067 4.561298e-02
#> C-A.1316        Armc10      C-A   -2.717958 -1.442523  -2.230867 4.603877e-02
#> C-A.1317        March7      C-A    2.100494  1.070728   2.230007 4.610953e-02
#> C-A.1318      Selenbp1      C-A   -3.790729 -1.922475  -2.229868 4.612099e-02
#> C-A.1319         Asah1      C-A   -2.254980 -1.173115  -2.228419 4.624056e-02
#> C-A.1320        Snapc2      C-A   -3.277428 -1.712564  -2.227985 4.627638e-02
#> C-A.1321        Man2c1      C-A   -2.514439 -1.330236  -2.223820 4.662186e-02
#> C-A.1322         Atg9a      C-A    3.591156  1.844448   2.220405 4.690689e-02
#> C-A.1323        Ddx26b      C-A   -5.163265 -2.368284  -2.219639 4.697110e-02
#> C-A.1324         Apex1      C-A   -2.231122 -1.157769  -2.219406 4.699060e-02
#> C-A.1325        Sptan1      C-A    2.183932  1.126928   2.219089 4.701719e-02
#> C-A.1326         Cxcr4      C-A   -2.836314 -1.504017  -2.216536 4.723187e-02
#> C-A.1327    D5Ertd579e      C-A    5.626124  2.492141   2.215877 4.728748e-02
#> C-A.1328        Pfkfb3      C-A    7.033982  2.814342   2.213791 4.746374e-02
#> C-A.1329       Ccdc88b      C-A    4.612361  2.205505   2.211070 4.769467e-02
#> C-A.1330        Abhd10      C-A    3.499043  1.806960   2.210787 4.771876e-02
#> C-A.1331         Grwd1      C-A   -2.633190 -1.396812  -2.205297 4.818812e-02
#> C-A.1332         Vezf1      C-A    2.725251  1.446389   2.204849 4.822660e-02
#> C-A.1333          Pdp2      C-A    7.668569  2.938957   2.204232 4.827964e-02
#> C-A.1334         Pank3      C-A    4.013392  2.004822   2.203317 4.835841e-02
#> C-A.1335        N6amt1      C-A    6.699900  2.744140   2.202866 4.839732e-02
#> C-A.1336         Prpf3      C-A   -2.565533 -1.359259  -2.201648 4.850245e-02
#> C-A.1337      Mphosph8      C-A    5.190347  2.375831   2.201638 4.850338e-02
#> C-A.1338           Eed      C-A   -2.103009 -1.072455  -2.201372 4.852631e-02
#> C-A.1339       Arl14ep      C-A    2.632957  1.396684   2.200996 4.855889e-02
#> C-A.1340       St3gal3      C-A    4.356711  2.123239   2.199229 4.871198e-02
#> C-A.1341          Hpse      C-A   -3.623093 -1.857222  -2.197915 4.882612e-02
#> C-A.1342         Eid2b      C-A   10.501766  3.392560   2.197455 4.886611e-02
#> C-A.1343         Clcn3      C-A    3.016921  1.593077   2.194440 4.912925e-02
#> C-A.1344         Qrsl1      C-A    6.985332  2.804329   2.194074 4.916128e-02
#> C-A.1345         Ccbl2      C-A   -4.128281 -2.045541  -2.192932 4.926131e-02
#> C-A.1346         Zfp26      C-A   -5.549324 -2.472312  -2.192501 4.929916e-02
#> C-A.1347         Zmym5      C-A    3.237258  1.694773   2.192344 4.931291e-02
#> C-A.1348         Ep400      C-A   -2.049031 -1.034942  -2.191449 4.939153e-02
#> C-A.1349        Ldoc1l      C-A   -2.650024 -1.406006  -2.191111 4.942125e-02
#> C-A.1350          Map7      C-A   -2.392709 -1.258645  -2.186737 4.980751e-02
#> C-A.1351          Dna2      C-A   -5.320232 -2.411489  -2.185585 4.990971e-02
#> C-A.1352          Cul2      C-A   -2.136858 -1.095491  -2.185491 4.991801e-02
#> C-A.1353         Thoc6      C-A   -2.383210 -1.252906  -2.184712 4.998729e-02
#> B-C.1             Lmo4      B-C  -14.629841 -3.870842 -14.826704 5.946170e-09
#> B-C.2             Cd82      B-C  -11.981857 -3.582780 -10.397522 2.885847e-07
#> B-C.3            H2afy      B-C    3.178350  1.668278  10.346844 3.040030e-07
#> B-C.4           Ifngr1      B-C   -4.318087 -2.110392 -10.266901 3.301570e-07
#> B-C.5         Slc25a12      B-C    8.237681  3.042238  10.089732 3.971659e-07
#> B-C.6            Il2rg      B-C   -5.684932 -2.507143  -9.879782 4.960919e-07
#> B-C.7             Map4      B-C   -7.030684 -2.813665  -9.752136 5.689795e-07
#> B-C.8            Itgb3      B-C  -10.977432 -3.456469  -9.660549 6.283435e-07
#> B-C.9           Ahcyl2      B-C  -12.827720 -3.681193  -9.606961 6.661416e-07
#> B-C.10            Sell      B-C    5.009452  2.324653   9.335240 8.993944e-07
#> B-C.11            Irf8      B-C   12.437638  3.636641   9.290554 9.455332e-07
#> B-C.12        Tmem176b      B-C   -3.517673 -1.814621  -9.275915 9.611968e-07
#> B-C.13            Faah      B-C   -5.131507 -2.359383  -9.261572 9.768142e-07
#> B-C.14            Rora      B-C   -9.330909 -3.222018  -9.233692 1.007956e-06
#> B-C.15           Plac8      B-C   13.360324  3.739883   9.112290 1.156559e-06
#> B-C.16           Nfkb1      B-C    4.477449  2.162677   9.106480 1.164238e-06
#> B-C.17             Emb      B-C   -2.812852 -1.492034  -8.875629 1.518387e-06
#> B-C.18          Shisa5      B-C   -4.285364 -2.099418  -8.552644 2.221504e-06
#> B-C.19            Ctr9      B-C    2.607732  1.382796   8.459175 2.485108e-06
#> B-C.20          Dnajc3      B-C   -3.344088 -1.741613  -8.323623 2.928746e-06
#> B-C.21             Id2      B-C   -2.348624 -1.231816  -8.272192 3.118683e-06
#> B-C.22         Hnrnpa1      B-C    2.727588  1.447626   8.140359 3.668543e-06
#> B-C.23           Rftn1      B-C   -2.605938 -1.381803  -8.105340 3.831476e-06
#> B-C.24           Napsa      B-C    4.497329  2.169068   8.011459 4.307754e-06
#> B-C.25         Tspan13      B-C   -4.068025 -2.024328  -7.856970 5.235081e-06
#> B-C.26          Samsn1      B-C   -6.609430 -2.724526  -7.827347 5.436177e-06
#> B-C.27           Prdx6      B-C   -2.617017 -1.387923  -7.764850 5.888094e-06
#> B-C.28       Serpinb1a      B-C   -4.460892 -2.157332  -7.750205 5.999705e-06
#> B-C.29            Lat2      B-C    5.802524  2.536681   7.717495 6.257252e-06
#> B-C.30           Itga4      B-C    2.962468  1.566799   7.689764 6.484870e-06
#> B-C.31             Pkm      B-C   -2.041998 -1.029982  -7.636515 6.947170e-06
#> B-C.32            Cd53      B-C   -3.060669 -1.613847  -7.337004 1.029796e-05
#> B-C.33           Ptpn6      B-C    3.889165  1.959460   7.267508 1.130026e-05
#> B-C.34            Ccr2      B-C   -4.848963 -2.277676  -7.211128 1.218996e-05
#> B-C.35          Galnt1      B-C   -2.411842 -1.270135  -7.207163 1.225529e-05
#> B-C.36          Il17rb      B-C  -18.909239 -4.241019  -7.152603 1.319318e-05
#> B-C.37          Stk17b      B-C   -5.101359 -2.350882  -7.088493 1.439421e-05
#> B-C.38           Gpr97      B-C    4.668440  2.222941   7.003044 1.617965e-05
#> B-C.39         Plekhb2      B-C   -3.337568 -1.738797  -6.983764 1.661431e-05
#> B-C.40           Itgal      B-C   -3.124347 -1.643555  -6.982573 1.664156e-05
#> B-C.41        Marcksl1      B-C    3.893893  1.961213   6.941915 1.760126e-05
#> B-C.42           Prkcq      B-C   -6.314862 -2.658751  -6.882289 1.911698e-05
#> B-C.43           Gpr56      B-C   26.555813  4.730956   6.844162 2.015873e-05
#> B-C.44        BC035044      B-C   10.398088  3.378246   6.839827 2.028095e-05
#> B-C.45             Cfp      B-C    8.354992  3.062638   6.647305 2.659174e-05
#> B-C.46           Ero1l      B-C   -3.161760 -1.660728  -6.549089 3.059033e-05
#> B-C.47          Rasal3      B-C   -2.205045 -1.140808  -6.540021 3.099058e-05
#> B-C.48            Ncf1      B-C   16.370588  4.033034   6.438657 3.586368e-05
#> B-C.49        Tmem176a      B-C   -5.610672 -2.488174  -6.409640 3.740423e-05
#> B-C.50           Nrros      B-C    2.892289  1.532212   6.357860 4.033100e-05
#> B-C.51           Rogdi      B-C   10.292016  3.363454   6.269420 4.590782e-05
#> B-C.52            Gpx1      B-C    2.472962  1.306240   6.243749 4.767596e-05
#> B-C.53           Gpr65      B-C   -4.953055 -2.308319  -6.236816 4.816580e-05
#> B-C.54            Idh2      B-C    2.821478  1.496451   6.226909 4.887504e-05
#> B-C.55           Stim1      B-C   -3.237752 -1.694992  -6.159777 5.398348e-05
#> B-C.56          Pabpc1      B-C    2.226964  1.155078   6.080460 6.075985e-05
#> B-C.57            Cerk      B-C   -4.582720 -2.196204  -6.069442 6.177045e-05
#> B-C.58            Fli1      B-C   -2.926361 -1.549108  -6.069343 6.177957e-05
#> B-C.59           Sept1      B-C   -2.804881 -1.487940  -6.034426 6.510220e-05
#> B-C.60           Esyt2      B-C   -3.444259 -1.784194  -6.003483 6.820533e-05
#> B-C.61            Snx2      B-C   -2.758691 -1.463984  -5.983246 7.031943e-05
#> B-C.62            Grb2      B-C    2.457600  1.297250   5.902915 7.942228e-05
#> B-C.63            Gclc      B-C   -7.022589 -2.812003  -5.864047 8.426793e-05
#> B-C.64           Esyt1      B-C   -2.037656 -1.026910  -5.854766 8.547073e-05
#> B-C.65         Smpdl3a      B-C   -4.338669 -2.117253  -5.851226 8.593429e-05
#> B-C.66             Ltb      B-C  -48.754237 -5.607456  -5.823340 8.968037e-05
#> B-C.67           Skp1a      B-C    2.093371  1.065828   5.820649 9.005088e-05
#> B-C.68            Eya2      B-C   -8.151008 -3.026978  -5.759238 9.896333e-05
#> B-C.69           Cmtm7      B-C    2.131788  1.092064   5.720584 1.050486e-04
#> B-C.70            Gapt      B-C   15.891042  3.990142   5.695621 1.091880e-04
#> B-C.71           Gata3      B-C   -4.222949 -2.078251  -5.692598 1.097010e-04
#> B-C.72            Aff4      B-C   -3.518382 -1.814912  -5.684975 1.110057e-04
#> B-C.73            Cd34      B-C   18.196817  4.185614   5.673666 1.129715e-04
#> B-C.74            Rbl2      B-C   -3.714155 -1.893034  -5.651821 1.168740e-04
#> B-C.75            Btg1      B-C   -3.293699 -1.719709  -5.630672 1.207883e-04
#> B-C.76          Ndfip1      B-C   -2.413225 -1.270962  -5.625802 1.217092e-04
#> B-C.77          Slc9a9      B-C   -3.422859 -1.775202  -5.616797 1.234315e-04
#> B-C.78          Cdkn1a      B-C   15.368687  3.941922   5.585053 1.297115e-04
#> B-C.79            Lpxn      B-C   -4.300123 -2.104378  -5.560709 1.347563e-04
#> B-C.80       Epb4.1l4b      B-C    4.235775  2.082626   5.557860 1.353600e-04
#> B-C.81         Fam102a      B-C   -4.579991 -2.195345  -5.531001 1.411948e-04
#> B-C.82           Lpar6      B-C   -5.028845 -2.330227  -5.518762 1.439413e-04
#> B-C.83            Lcp2      B-C   -3.703981 -1.889077  -5.497671 1.488076e-04
#> B-C.84           Cmtm6      B-C   -2.364864 -1.241757  -5.486588 1.514341e-04
#> B-C.85          Pdgfrb      B-C    7.120167  2.831911   5.469792 1.555084e-04
#> B-C.86       Serpina3g      B-C   -7.158317 -2.839620  -5.454070 1.594272e-04
#> B-C.87           Cpne2      B-C    4.156962  2.055530   5.391237 1.761606e-04
#> B-C.88           Itm2c      B-C   -3.081061 -1.623427  -5.284110 2.091080e-04
#> B-C.89           Clic4      B-C    3.201533  1.678763   5.270822 2.136267e-04
#> B-C.90           Tgtp2      B-C   -4.626190 -2.209825  -5.257123 2.183933e-04
#> B-C.91          Ppfia1      B-C    4.227823  2.079915   5.238176 2.251717e-04
#> B-C.92             Bid      B-C    5.780014  2.531073   5.224588 2.301688e-04
#> B-C.93        AU040320      B-C   -3.114015 -1.638776  -5.217991 2.326373e-04
#> B-C.94            Sox4      B-C    2.514201  1.330100   5.213952 2.341621e-04
#> B-C.95        Rap1gds1      B-C   -3.606179 -1.850471  -5.197216 2.405940e-04
#> B-C.96         Tmem154      B-C  -20.122145 -4.330712  -5.180480 2.472123e-04
#> B-C.97           Padi2      B-C    2.671399  1.417595   5.149438 2.600005e-04
#> B-C.98           Dapp1      B-C    3.093359  1.629174   5.132116 2.674379e-04
#> B-C.99           Skap2      B-C    4.418413  2.143528   5.103678 2.801375e-04
#> B-C.100           Mycn      B-C    3.700215  1.887609   5.095032 2.841232e-04
#> B-C.101        Ppp1r18      B-C   -2.998680 -1.584328  -5.065992 2.979530e-04
#> B-C.102        Themis2      B-C   50.112151  5.647089   5.061918 2.999493e-04
#> B-C.103  4930523C07Rik      B-C  -10.349199 -3.371447  -5.040761 3.105440e-04
#> B-C.104        Unc93b1      B-C    4.762714  2.251784   5.037085 3.124245e-04
#> B-C.105          Nup85      B-C    2.538809  1.344152   5.017119 3.228501e-04
#> B-C.106           Cd84      B-C   -4.128560 -2.045639  -5.009720 3.268061e-04
#> B-C.107           Ctso      B-C   -5.247724 -2.391692  -5.008136 3.276597e-04
#> B-C.108        Tnfaip1      B-C  -57.095669 -5.835309  -4.981017 3.426410e-04
#> B-C.109            Wls      B-C   -3.368510 -1.752110  -4.977945 3.443828e-04
#> B-C.110          Ncoa4      B-C    2.485732  1.313671   4.969266 3.493543e-04
#> B-C.111         Il10ra      B-C    8.843266  3.144579   4.965192 3.517144e-04
#> B-C.112          Mgat1      B-C    2.612090  1.385205   4.954447 3.580186e-04
#> B-C.113         Rnf145      B-C   -2.431803 -1.282027  -4.952654 3.590824e-04
#> B-C.114          Plcg2      B-C    2.341241  1.227274   4.946388 3.628252e-04
#> B-C.115          Satb1      B-C    4.634266  2.212341   4.931977 3.715890e-04
#> B-C.116           Elp2      B-C    2.772499  1.471187   4.931841 3.716728e-04
#> B-C.117          Prr13      B-C   -2.084972 -1.060028  -4.901969 3.905600e-04
#> B-C.118          Pdcd1      B-C   -4.843750 -2.276124  -4.896712 3.939872e-04
#> B-C.119         Maged1      B-C  -51.457870 -5.685320  -4.877908 4.065051e-04
#> B-C.120         Il27ra      B-C   -5.817526 -2.540406  -4.877600 4.067139e-04
#> B-C.121       Cyb561a3      B-C    3.367219  1.751557   4.836573 4.355133e-04
#> B-C.122         Gpr171      B-C    2.496066  1.319656   4.835973 4.359502e-04
#> B-C.123           Mcm7      B-C    2.042871  1.030598   4.816901 4.500735e-04
#> B-C.124        Spata13      B-C  -11.156302 -3.479787  -4.816833 4.501251e-04
#> B-C.125         Ormdl3      B-C   -6.153965 -2.621516  -4.811926 4.538366e-04
#> B-C.126          Nfam1      B-C    3.703234  1.888786   4.801294 4.619884e-04
#> B-C.127           Cry1      B-C   -3.083545 -1.624590  -4.789443 4.712560e-04
#> B-C.128           Ybx3      B-C    2.433469  1.283015   4.774898 4.828978e-04
#> B-C.129         Ptpn22      B-C   -3.200396 -1.678250  -4.768963 4.877345e-04
#> B-C.130        Tmem260      B-C   -3.897196 -1.962436  -4.746222 5.067414e-04
#> B-C.131        Arhgef3      B-C  -76.592500 -6.259131  -4.738478 5.133899e-04
#> B-C.132          Dusp6      B-C  -61.925088 -5.952452  -4.727320 5.231308e-04
#> B-C.133         Slain2      B-C   -2.248111 -1.168713  -4.724403 5.257092e-04
#> B-C.134        Cbfa2t3      B-C    9.053897  3.178539   4.723160 5.268120e-04
#> B-C.135         Tgfbr2      B-C  -11.047534 -3.465652  -4.715803 5.333877e-04
#> B-C.136          Erp29      B-C    2.804993  1.487997   4.706506 5.418220e-04
#> B-C.137          Furin      B-C   -2.573687 -1.363836  -4.652350 5.938045e-04
#> B-C.138         Pom121      B-C    9.665170  3.272795   4.644944 6.013089e-04
#> B-C.139           Cish      B-C  -60.037828 -5.907800  -4.640123 6.062465e-04
#> B-C.140           Gm2a      B-C   -2.773383 -1.471647  -4.639399 6.069924e-04
#> B-C.141          Cd24a      B-C    9.584826  3.260752   4.629662 6.171090e-04
#> B-C.142        Slco3a1      B-C   -3.699245 -1.887231  -4.617819 6.296520e-04
#> B-C.143           Prr5      B-C   14.415617  3.849561   4.613713 6.340635e-04
#> B-C.144          Zfp35      B-C    2.689792  1.427495   4.593605 6.561354e-04
#> B-C.145          Cnpy3      B-C    2.747108  1.457914   4.593549 6.561974e-04
#> B-C.146         Ltb4r1      B-C  -15.826055 -3.984230  -4.590653 6.594418e-04
#> B-C.147          Runx3      B-C    2.302172  1.202996   4.586643 6.639630e-04
#> B-C.148          Farsa      B-C    2.324653  1.217015   4.583362 6.676857e-04
#> B-C.149           Cd81      B-C   -3.128525 -1.645483  -4.579277 6.723512e-04
#> B-C.150           Tns1      B-C   10.438448  3.383835   4.569024 6.842124e-04
#> B-C.151         Gtpbp4      B-C    2.982686  1.576612   4.559343 6.956130e-04
#> B-C.152            Syk      B-C   14.821392  3.889609   4.547608 7.096981e-04
#> B-C.153         Fbxo33      B-C   -4.191356 -2.067417  -4.541795 7.167863e-04
#> B-C.154          Abce1      B-C    2.424254  1.277541   4.522688 7.406080e-04
#> B-C.155            Dtl      B-C    2.748396  1.458590   4.518473 7.459734e-04
#> B-C.156          Ppm1l      B-C   -8.011987 -3.002160  -4.515681 7.495492e-04
#> B-C.157            Vim      B-C   -2.778650 -1.474384  -4.513843 7.519139e-04
#> B-C.158          Irgm1      B-C   -2.724923 -1.446215  -4.481581 7.947029e-04
#> B-C.159         Camk2d      B-C   -4.365016 -2.125987  -4.474974 8.037748e-04
#> B-C.160        Skiv2l2      B-C    2.447795  1.291483   4.464226 8.187654e-04
#> B-C.161           Ets1      B-C  -32.658401 -5.029382  -4.457525 8.282577e-04
#> B-C.162        Tcrg-C1      B-C   -3.910195 -1.967241  -4.446529 8.440858e-04
#> B-C.163      Trp53inp1      B-C   -4.370352 -2.127750  -4.445390 8.457424e-04
#> B-C.164         Rcbtb2      B-C   -2.823393 -1.497430  -4.442417 8.500845e-04
#> B-C.165          Nfil3      B-C    9.961713  3.316394   4.425935 8.745820e-04
#> B-C.166          Ttyh3      B-C    3.855806  1.947032   4.404241 9.079525e-04
#> B-C.167         Tmem71      B-C  -33.028048 -5.045620  -4.393605 9.247959e-04
#> B-C.168          Snx10      B-C    2.874068  1.523094   4.392922 9.258883e-04
#> B-C.169    D19Bwg1357e      B-C    2.596878  1.376778   4.391597 9.280123e-04
#> B-C.170     St6galnac4      B-C    3.561709  1.832570   4.389891 9.307536e-04
#> B-C.171         Pdlim2      B-C  -27.222747 -4.766741  -4.373612 9.573436e-04
#> B-C.172           Gbp2      B-C  -55.672151 -5.798884  -4.372494 9.591988e-04
#> B-C.173           Mfn2      B-C   -2.228198 -1.155878  -4.349001 9.990552e-04
#> B-C.174         Fam65a      B-C    2.776960  1.473507   4.345633 1.004908e-03
#> B-C.175          Lonp2      B-C   -2.566771 -1.359954  -4.333916 1.025553e-03
#> B-C.176           Cttn      B-C   16.561114  4.049728   4.312861 1.063767e-03
#> B-C.177          Anxa2      B-C    2.506778  1.325834   4.298406 1.090857e-03
#> B-C.178          Clip1      B-C   -3.066192 -1.616448  -4.298374 1.090919e-03
#> B-C.179          Mapk9      B-C   -2.440934 -1.287433  -4.282337 1.121817e-03
#> B-C.180           Sat1      B-C   -2.840049 -1.505916  -4.282224 1.122039e-03
#> B-C.181          Med12      B-C   -3.363576 -1.749996  -4.280234 1.125937e-03
#> B-C.182          Hmga1      B-C    3.024451  1.596673   4.275480 1.135306e-03
#> B-C.183          Abca1      B-C   -4.319563 -2.110885  -4.266440 1.153346e-03
#> B-C.184           Pgls      B-C    2.175636  1.121437   4.256264 1.174010e-03
#> B-C.185         Rundc1      B-C   -3.749049 -1.906525  -4.255545 1.175484e-03
#> B-C.186          Limk1      B-C    7.327300  2.873282   4.237686 1.212736e-03
#> B-C.187         Tyrobp      B-C    2.856939  1.514470   4.237458 1.213219e-03
#> B-C.188            Hlf      B-C  -15.478547 -3.952198  -4.232626 1.223512e-03
#> B-C.189          Papd4      B-C    3.021727  1.595373   4.229097 1.231084e-03
#> B-C.190         Il18r1      B-C  -15.125039 -3.918867  -4.228670 1.232005e-03
#> B-C.191           Zhx1      B-C   -2.049527 -1.035291  -4.228453 1.232472e-03
#> B-C.192        Fam129a      B-C   -7.266627 -2.861286  -4.219956 1.250932e-03
#> B-C.193        Slc35c2      B-C    2.281107  1.189734   4.211815 1.268887e-03
#> B-C.194        Plekha2      B-C   -2.335856 -1.223951  -4.203548 1.287394e-03
#> B-C.195           Cd63      B-C   13.190151  3.721389   4.190869 1.316324e-03
#> B-C.196       AI467606      B-C   -2.121372 -1.084997  -4.187491 1.324146e-03
#> B-C.197            Grn      B-C    2.977277  1.573993   4.177555 1.347430e-03
#> B-C.198           Hcst      B-C  -27.045869 -4.757336  -4.175367 1.352614e-03
#> B-C.199            Myc      B-C    2.114005  1.079979   4.173425 1.357235e-03
#> B-C.200          Abhd8      B-C   -6.017644 -2.589199  -4.171509 1.361808e-03
#> B-C.201         Kansl2      B-C    2.290676  1.195773   4.167417 1.371629e-03
#> B-C.202          Hmgn3      B-C  -36.429726 -5.187044  -4.144621 1.427701e-03
#> B-C.203          Mef2c      B-C    5.837546  2.545362   4.140670 1.437657e-03
#> B-C.204          Plin3      B-C  -14.605283 -3.868418  -4.137609 1.445421e-03
#> B-C.205          Tgtp1      B-C   -8.469459 -3.082270  -4.135341 1.451201e-03
#> B-C.206           Aff3      B-C    2.689209  1.427182   4.126326 1.474418e-03
#> B-C.207           Heca      B-C   -5.982239 -2.580686  -4.114492 1.505479e-03
#> B-C.208          Litaf      B-C   11.790693  3.559577   4.111862 1.512474e-03
#> B-C.209          Plcb4      B-C   -5.579540 -2.480146  -4.110041 1.517336e-03
#> B-C.210         Eif4e3      B-C   -3.880008 -1.956060  -4.102576 1.537440e-03
#> B-C.211            Vcl      B-C   19.722185  4.301748   4.099602 1.545526e-03
#> B-C.212            Npl      B-C  -16.012508 -4.001127  -4.091942 1.566558e-03
#> B-C.213        Fam134c      B-C   -2.854033 -1.513002  -4.089986 1.571973e-03
#> B-C.214           Spi1      B-C   20.008951  4.322574   4.080643 1.598122e-03
#> B-C.215           Ppt2      B-C   -8.317098 -3.056080  -4.063200 1.648150e-03
#> B-C.216          Fnbp1      B-C   -2.084354 -1.059600  -4.057168 1.665826e-03
#> B-C.217             Hp      B-C   11.013201  3.461162   4.055441 1.670923e-03
#> B-C.218        Slc22a3      B-C   14.946572  3.901743   4.052168 1.680626e-03
#> B-C.219         Mgat4a      B-C   -3.932284 -1.975368  -4.050205 1.686474e-03
#> B-C.220        Il12rb1      B-C  -14.611407 -3.869023  -4.048213 1.692427e-03
#> B-C.221  2810474O19Rik      B-C   -2.931394 -1.551587  -4.026943 1.757396e-03
#> B-C.222         Ogfrl1      B-C    4.052057  2.018655   4.026248 1.759562e-03
#> B-C.223          Nrip1      B-C   -2.407049 -1.267265  -4.016445 1.790405e-03
#> B-C.224         Tm7sf3      B-C    2.305864  1.205307   4.013167 1.800845e-03
#> B-C.225          Ikzf2      B-C   -5.365824 -2.423800  -4.008151 1.816938e-03
#> B-C.226           Flt3      B-C   45.013858  5.492297   4.001609 1.838155e-03
#> B-C.227          Spns3      B-C    5.397408  2.432267   3.998899 1.847019e-03
#> B-C.228            Dek      B-C    2.027526  1.019720   3.998276 1.849061e-03
#> B-C.229          Ninj1      B-C    6.414450  2.681326   3.998257 1.849126e-03
#> B-C.230            Fyn      B-C   -2.869707 -1.520903  -3.990587 1.874481e-03
#> B-C.231         Lpcat4      B-C   -2.990366 -1.580322  -3.983975 1.896628e-03
#> B-C.232        Smarcc1      B-C    2.098665  1.069472   3.982725 1.900845e-03
#> B-C.233        Ehbp1l1      B-C   -2.253212 -1.171983  -3.981071 1.906438e-03
#> B-C.234           Fgd2      B-C   27.807132  4.797383   3.980611 1.907998e-03
#> B-C.235           Orc2      B-C    8.910611  3.155524   3.979907 1.910387e-03
#> B-C.236         Spata6      B-C  -26.744561 -4.741174  -3.978464 1.915293e-03
#> B-C.237         Zfp148      B-C   -2.166301 -1.115234  -3.973404 1.932597e-03
#> B-C.238          Fnbp4      B-C    6.542782  2.709904   3.956075 1.993096e-03
#> B-C.239          Ddx27      B-C    2.725046  1.446280   3.953050 2.003857e-03
#> B-C.240          Anxa6      B-C   -2.978167 -1.574424  -3.950420 2.013260e-03
#> B-C.241           Ctps      B-C    2.847284  1.509586   3.946850 2.026096e-03
#> B-C.242          Btaf1      B-C    2.925604  1.548735   3.944924 2.033058e-03
#> B-C.243          Rab44      B-C   71.236964  6.154554   3.942491 2.041885e-03
#> B-C.244             Qk      B-C    2.265504  1.179832   3.942127 2.043211e-03
#> B-C.245          Socs1      B-C  -12.468949 -3.640268  -3.937211 2.061183e-03
#> B-C.246          Ap1g1      B-C   -3.046371 -1.607092  -3.895572 2.220109e-03
#> B-C.247        Tmem63a      B-C   -7.458115 -2.898811  -3.894829 2.223055e-03
#> B-C.248          Rap2c      B-C   -9.924063 -3.310931  -3.885991 2.258433e-03
#> B-C.249        Laptm4b      B-C   -2.882517 -1.527329  -3.885106 2.262010e-03
#> B-C.250          Fmnl3      B-C    2.739051  1.453676   3.883452 2.268704e-03
#> B-C.251          Runx1      B-C    2.138876  1.096853   3.879558 2.284548e-03
#> B-C.252         Plxnb2      B-C    9.592742  3.261943   3.878500 2.288872e-03
#> B-C.253        Tubgcp3      B-C    6.356395  2.668209   3.876006 2.299102e-03
#> B-C.254        Galnt12      B-C  -15.151196 -3.921360  -3.873843 2.308013e-03
#> B-C.255          Spns2      B-C   41.441983  5.373021   3.871596 2.317303e-03
#> B-C.256          Kdm7a      B-C   -2.855521 -1.513754  -3.863498 2.351115e-03
#> B-C.257          Tifab      B-C   10.321432  3.367571   3.855775 2.383831e-03
#> B-C.258         H2-DMa      B-C    2.729016  1.448381   3.841727 2.444558e-03
#> B-C.259         Map2k1      B-C   -2.104587 -1.073537  -3.839435 2.454612e-03
#> B-C.260       Leprotl1      B-C   -2.296603 -1.199502  -3.829942 2.496726e-03
#> B-C.261         Tardbp      B-C    2.138212  1.096405   3.827778 2.506431e-03
#> B-C.262          Pcgf6      B-C   18.758038  4.229437   3.824390 2.521702e-03
#> B-C.263         Trip12      B-C   -2.192576 -1.132627  -3.824307 2.522078e-03
#> B-C.264          Adap1      B-C   28.850512  4.850525   3.819600 2.543454e-03
#> B-C.265          Anxa1      B-C   -3.698329 -1.886874  -3.818159 2.550037e-03
#> B-C.266           Palm      B-C    4.594710  2.199974   3.817887 2.551282e-03
#> B-C.267          Ckap5      B-C   -2.897123 -1.534621  -3.814009 2.569089e-03
#> B-C.268        Alox5ap      B-C   29.085165  4.862212   3.811742 2.579557e-03
#> B-C.269        Dpy19l3      B-C  -19.331558 -4.272886  -3.805927 2.606618e-03
#> B-C.270         Ms4a6c      B-C   15.205954  3.926564   3.798245 2.642808e-03
#> B-C.271        Zdhhc18      B-C   -2.656855 -1.409720  -3.795066 2.657939e-03
#> B-C.272          Itgam      B-C   27.968831  4.805748   3.787698 2.693345e-03
#> B-C.273          Timp2      B-C   -8.188134 -3.033535  -3.784958 2.706636e-03
#> B-C.274           Dok2      B-C   -2.821816 -1.496624  -3.779677 2.732442e-03
#> B-C.275          Capn2      B-C   -3.462001 -1.791606  -3.772864 2.766112e-03
#> B-C.276         Tspan3      B-C    2.257946  1.175011   3.767579 2.792523e-03
#> B-C.277           Ncf2      B-C    3.818783  1.933113   3.763251 2.814346e-03
#> B-C.278          Wdfy4      B-C    4.131136  2.046539   3.761422 2.823617e-03
#> B-C.279        Golph3l      B-C   -3.091723 -1.628411  -3.755408 2.854334e-03
#> B-C.280          Prtn3      B-C    7.845251  2.971820   3.753157 2.865923e-03
#> B-C.281          Itpr1      B-C    4.456550  2.155927   3.751583 2.874050e-03
#> B-C.282           Tifa      B-C    3.306820  1.725444   3.750033 2.882078e-03
#> B-C.283            Lyn      B-C    2.608699  1.383330   3.747956 2.892877e-03
#> B-C.284          U2af1      B-C    2.203335  1.139689   3.740687 2.930981e-03
#> B-C.285          Ncoa1      B-C   -3.326542 -1.734023  -3.739806 2.935636e-03
#> B-C.286        S100a10      B-C   -2.216925 -1.148560  -3.728806 2.994384e-03
#> B-C.287            Bsn      B-C   30.170710  4.915077   3.728330 2.996952e-03
#> B-C.288         Mtmr10      B-C  -20.504330 -4.357857  -3.727438 3.001776e-03
#> B-C.289          Tpd52      B-C    2.106390  1.074772   3.724514 3.017634e-03
#> B-C.290           Tbl3      B-C    2.329270  1.219878   3.723953 3.020687e-03
#> B-C.291         Suclg1      B-C    2.626347  1.393057   3.723824 3.021389e-03
#> B-C.292          Iigp1      B-C  -60.114535 -5.909642  -3.711373 3.089974e-03
#> B-C.293         Map3k5      B-C   -5.829545 -2.543383  -3.707849 3.109676e-03
#> B-C.294          Gcnt1      B-C   -2.557657 -1.354823  -3.707176 3.113453e-03
#> B-C.295       Tbc1d10c      B-C   -2.063055 -1.044782  -3.705157 3.124814e-03
#> B-C.296            Evl      B-C  -14.445642 -3.852562  -3.704470 3.128686e-03
#> B-C.297          Mppe1      B-C   -3.359639 -1.748306  -3.699622 3.156172e-03
#> B-C.298           Nrp1      B-C   -6.218330 -2.636527  -3.697694 3.167170e-03
#> B-C.299          Mast4      B-C   -6.492566 -2.698789  -3.695035 3.182404e-03
#> B-C.300         Cnot6l      B-C   -5.335300 -2.415569  -3.692392 3.197623e-03
#> B-C.301         Il20ra      B-C  -43.935260 -5.457307  -3.684089 3.245913e-03
#> B-C.302         Ddx26b      B-C   16.685587  4.060531   3.673273 3.309949e-03
#> B-C.303         Cited2      B-C   -2.788159 -1.479313  -3.654487 3.424263e-03
#> B-C.304           Rpf2      B-C    3.140454  1.650973   3.637627 3.530304e-03
#> B-C.305        Snrnp70      B-C    2.508890  1.327049   3.633992 3.553605e-03
#> B-C.306         Tm6sf1      B-C    4.635960  2.212868   3.626338 3.603183e-03
#> B-C.307            Hk1      B-C   -2.157847 -1.109593  -3.622357 3.629248e-03
#> B-C.308          Uhrf2      B-C   -2.841593 -1.506700  -3.621750 3.633240e-03
#> B-C.309         Glipr2      B-C   -8.896753 -3.153279  -3.619072 3.650899e-03
#> B-C.310           Ccr7      B-C   -8.216745 -3.038567  -3.617893 3.658707e-03
#> B-C.311           Ralb      B-C  -16.247781 -4.022171  -3.616614 3.667191e-03
#> B-C.312        Sh3kbp1      B-C   -2.126857 -1.088723  -3.615360 3.675529e-03
#> B-C.313           Tob1      B-C  -25.260287 -4.658799  -3.612526 3.694447e-03
#> B-C.314        Hnrnpdl      B-C    3.684087  1.881307   3.609768 3.712948e-03
#> B-C.315         Ankmy2      B-C   -3.145133 -1.653121  -3.599773 3.780815e-03
#> B-C.316            Rb1      B-C   -2.809729 -1.490431  -3.597632 3.795518e-03
#> B-C.317           Hes6      B-C    2.034511  1.024682   3.596836 3.801001e-03
#> B-C.318          Ddx54      B-C    2.226436  1.154736   3.596548 3.802981e-03
#> B-C.319           Inip      B-C    2.810884  1.491024   3.587319 3.867158e-03
#> B-C.320  I830077J02Rik      B-C   25.296190  4.660848   3.574676 3.956877e-03
#> B-C.321            Fes      B-C    2.162962  1.113008   3.560410 4.060667e-03
#> B-C.322        Dennd1b      B-C   -3.158059 -1.659038  -3.550928 4.131189e-03
#> B-C.323            Cpq      B-C  -16.692939 -4.061166  -3.536749 4.238996e-03
#> B-C.324           Nfe2      B-C    2.484683  1.313062   3.536706 4.239326e-03
#> B-C.325         Sh3bp5      B-C  -19.250368 -4.266814  -3.535579 4.248014e-03
#> B-C.326          Gstm1      B-C   -4.377231 -2.130019  -3.533117 4.267068e-03
#> B-C.327        St6gal1      B-C  -28.905675 -4.853281  -3.525497 4.326589e-03
#> B-C.328         Ruvbl1      B-C    2.369427  1.244538   3.524876 4.331473e-03
#> B-C.329           Pkn1      B-C   -2.327274 -1.218641  -3.522810 4.347776e-03
#> B-C.330          Nlrc3      B-C   -2.872027 -1.522069  -3.522119 4.353240e-03
#> B-C.331            Srm      B-C    2.038866  1.027767   3.521801 4.355760e-03
#> B-C.332         Frmd4b      B-C  -33.508526 -5.066456  -3.521716 4.356429e-03
#> B-C.333          Nsun2      B-C    2.120885  1.084667   3.517868 4.387020e-03
#> B-C.334         Nfkbia      B-C   -2.049409 -1.035208  -3.513582 4.421352e-03
#> B-C.335           Ppan      B-C    2.673696  1.418835   3.513041 4.425708e-03
#> B-C.336         Samhd1      B-C   -3.265603 -1.707349  -3.500644 4.526668e-03
#> B-C.337           Cd74      B-C   -2.236787 -1.161428  -3.496868 4.557884e-03
#> B-C.338         Bcl11a      B-C    8.142013  3.025386   3.493350 4.587161e-03
#> B-C.339       Cdc42ep3      B-C   -2.428119 -1.279839  -3.489975 4.615432e-03
#> B-C.340           Tti1      B-C    2.665053  1.414164   3.478631 4.711758e-03
#> B-C.341           Ldb1      B-C   -2.319508 -1.213819  -3.474893 4.743948e-03
#> B-C.342           Pmm2      B-C   -4.149928 -2.053086  -3.461294 4.862963e-03
#> B-C.343           Yaf2      B-C   -4.177355 -2.062590  -3.460973 4.865810e-03
#> B-C.344          Mdfic      B-C   -7.556741 -2.917764  -3.460298 4.871798e-03
#> B-C.345          Tmed3      B-C    2.642818  1.402077   3.460131 4.873283e-03
#> B-C.346           Cog4      B-C    5.107459  2.352606   3.452847 4.938409e-03
#> B-C.347       Colgalt1      B-C    2.227296  1.155293   3.452354 4.942857e-03
#> B-C.348           Cpa3      B-C    3.663112  1.873070   3.451489 4.950653e-03
#> B-C.349           Pus7      B-C   10.090281  3.334894   3.448762 4.975329e-03
#> B-C.350        Phf20l1      B-C   -2.797416 -1.484095  -3.448321 4.979329e-03
#> B-C.351         Atrnl1      B-C  -10.651302 -3.412958  -3.447495 4.986835e-03
#> B-C.352          Ints1      B-C    2.147117  1.102401   3.446882 4.992414e-03
#> B-C.353         Map4k2      B-C   -4.873171 -2.284861  -3.443989 5.018814e-03
#> B-C.354         Dopey1      B-C  -52.159229 -5.704851  -3.439726 5.057990e-03
#> B-C.355           Rbpj      B-C    2.134684  1.094023   3.437209 5.081264e-03
#> B-C.356           Ctsh      B-C    6.091160  2.606717   3.434850 5.103175e-03
#> B-C.357         Rnf166      B-C   -2.189337 -1.130494  -3.434287 5.108413e-03
#> B-C.358          Birc3      B-C  -33.392708 -5.061461  -3.430883 5.140237e-03
#> B-C.359          Ckap4      B-C   27.249380  4.768152   3.429616 5.152128e-03
#> B-C.360          Pycr1      B-C    5.488306  2.456361   3.419138 5.251584e-03
#> B-C.361          Nudt4      B-C   -2.687580 -1.426308  -3.414878 5.292576e-03
#> B-C.362         Kif13a      B-C  -26.936414 -4.751486  -3.406350 5.375620e-03
#> B-C.363            Pnn      B-C    2.380320  1.251156   3.394755 5.490662e-03
#> B-C.364          H2-Oa      B-C   -6.905313 -2.787707  -3.393398 5.504284e-03
#> B-C.365         Mettl4      B-C   -8.330647 -3.058429  -3.385891 5.580289e-03
#> B-C.366           Klf3      B-C   -9.653893 -3.271111  -3.384436 5.595143e-03
#> B-C.367           Ly86      B-C   13.745055  3.780841   3.378720 5.653882e-03
#> B-C.368           Nme4      B-C    2.489467  1.315837   3.376757 5.674202e-03
#> B-C.369        St3gal6      B-C  -77.208213 -6.270682  -3.376325 5.678686e-03
#> B-C.370           Rara      B-C   -3.101549 -1.632989  -3.376067 5.681363e-03
#> B-C.371            Myb      B-C    2.826428  1.498980   3.374354 5.699175e-03
#> B-C.372           Il7r      B-C  -34.847668 -5.122990  -3.371520 5.728771e-03
#> B-C.373           Cds2      B-C   -2.922578 -1.547242  -3.364891 5.798620e-03
#> B-C.374           Snrk      B-C   -2.312934 -1.209724  -3.356440 5.888905e-03
#> B-C.375           Rgs2      B-C    8.949999  3.161887   3.348811 5.971643e-03
#> B-C.376          Galk2      B-C   -2.394455 -1.259697  -3.347880 5.981823e-03
#> B-C.377           Aspm      B-C  -11.338212 -3.503121  -3.347215 5.989106e-03
#> B-C.378        Arl14ep      B-C   -6.741327 -2.753033  -3.343803 6.026596e-03
#> B-C.379        Atp13a2      B-C    2.275319  1.186069   3.339331 6.076107e-03
#> B-C.380          Panx1      B-C    2.521527  1.334298   3.339030 6.079456e-03
#> B-C.381           Rfc2      B-C    2.074603  1.052835   3.333395 6.142460e-03
#> B-C.382  D030056L22Rik      B-C   -5.146560 -2.363608  -3.331696 6.161587e-03
#> B-C.383        Fam117a      B-C    2.413448  1.271096   3.326617 6.219119e-03
#> B-C.384          Hmgb3      B-C    2.374902  1.247868   3.325894 6.227356e-03
#> B-C.385         Tex264      B-C   -3.020336 -1.594709  -3.322869 6.261929e-03
#> B-C.386         Sigirr      B-C   -2.750564 -1.459727  -3.322533 6.265776e-03
#> B-C.387         Lmbrd1      B-C  -24.444103 -4.611415  -3.316474 6.335656e-03
#> B-C.388          Sbno1      B-C    2.265913  1.180092   3.314563 6.357865e-03
#> B-C.389        Wbscr16      B-C    6.978128  2.802840   3.307218 6.443940e-03
#> B-C.390          Wdr18      B-C    2.235072  1.160322   3.305948 6.458943e-03
#> B-C.391            Btk      B-C    6.849694  2.776040   3.301249 6.514764e-03
#> B-C.392           Sdc1      B-C -138.474097 -7.113472  -3.297510 6.559522e-03
#> B-C.393         Zfp868      B-C    4.705636  2.234390   3.291117 6.636786e-03
#> B-C.394        Plekhm1      B-C   -2.045439 -1.032410  -3.289108 6.661250e-03
#> B-C.395          Fbxo4      B-C  -10.733426 -3.424039  -3.288130 6.673204e-03
#> B-C.396        Slc11a2      B-C   -4.804393 -2.264354  -3.286275 6.695913e-03
#> B-C.397          Enpp4      B-C   -2.795720 -1.483220  -3.274145 6.846402e-03
#> B-C.398         Rnase6      B-C   16.110306  4.009912   3.273643 6.852700e-03
#> B-C.399        Ccdc186      B-C  -11.934882 -3.577112  -3.268017 6.923732e-03
#> B-C.400        Poglut1      B-C   -3.265432 -1.707274  -3.266541 6.942478e-03
#> B-C.401         B3gnt1      B-C  -18.293959 -4.193295  -3.262000 7.000515e-03
#> B-C.402           Accs      B-C  -18.407709 -4.202238  -3.260225 7.023325e-03
#> B-C.403       Slc25a13      B-C   22.935730  4.519525   3.260034 7.025789e-03
#> B-C.404         Gm7694      B-C    5.315758  2.410276   3.256586 7.070338e-03
#> B-C.405         Arrdc4      B-C   -6.061298 -2.599627  -3.252426 7.124469e-03
#> B-C.406         Rabac1      B-C   -2.281243 -1.189820  -3.249359 7.164639e-03
#> B-C.407           Arl1      B-C    2.079694  1.056371   3.243712 7.239218e-03
#> B-C.408           Hars      B-C    2.006632  1.004776   3.243622 7.240414e-03
#> B-C.409          Arl4a      B-C   -3.814842 -1.931623  -3.240110 7.287190e-03
#> B-C.410          Trip4      B-C   -6.095414 -2.607724  -3.239402 7.296663e-03
#> B-C.411  2700049A03Rik      B-C    5.587573  2.482222   3.237554 7.321431e-03
#> B-C.412          Irgm2      B-C   -2.784737 -1.477541  -3.237287 7.325022e-03
#> B-C.413          Pde4b      B-C   -3.335737 -1.738006  -3.236046 7.341713e-03
#> B-C.414          Fbxl5      B-C   -2.774111 -1.472026  -3.230947 7.410694e-03
#> B-C.415           Ier2      B-C   -2.277734 -1.187600  -3.230748 7.413395e-03
#> B-C.416        Slc29a3      B-C    3.660327  1.871972   3.227059 7.463742e-03
#> B-C.417         Ldoc1l      B-C    4.335055  2.116050   3.227030 7.464128e-03
#> B-C.418            Mbp      B-C   -2.280221 -1.189174  -3.226301 7.474122e-03
#> B-C.419          Mrps7      B-C    2.169233  1.117185   3.224427 7.499867e-03
#> B-C.420          Dgcr2      B-C    2.090355  1.063748   3.223558 7.511830e-03
#> B-C.421         Pycard      B-C    2.124131  1.086873   3.219413 7.569169e-03
#> B-C.422           Ncln      B-C    2.450353  1.292990   3.217286 7.598768e-03
#> B-C.423        Zc3hav1      B-C   -2.044810 -1.031967  -3.214917 7.631873e-03
#> B-C.424          Cops2      B-C    2.670917  1.417335   3.214101 7.643299e-03
#> B-C.425         Pi4k2a      B-C   -2.788944 -1.479719  -3.211550 7.679166e-03
#> B-C.426        Dennd4c      B-C   -2.895766 -1.533945  -3.208276 7.725447e-03
#> B-C.427        Nsmce4a      B-C    2.310646  1.208296   3.207824 7.731844e-03
#> B-C.428          Wdr35      B-C  -31.040315 -4.956071  -3.202708 7.804790e-03
#> B-C.429        Rtn4ip1      B-C    6.087647  2.605885   3.202191 7.812199e-03
#> B-C.430         Tmem64      B-C   -4.039462 -2.014163  -3.199626 7.849057e-03
#> B-C.431          Hspa2      B-C   15.760724  3.978262   3.199275 7.854120e-03
#> B-C.432           Gga2      B-C    5.093882  2.348765   3.198637 7.863313e-03
#> B-C.433       Arhgap17      B-C    3.354436  1.746070   3.195594 7.907370e-03
#> B-C.434            Fn1      B-C  -30.767293 -4.943326  -3.195111 7.914376e-03
#> B-C.435          Ccnl2      B-C   -2.471354 -1.305302  -3.192454 7.953074e-03
#> B-C.436      Uhrf1bp1l      B-C   -3.401814 -1.766304  -3.191213 7.971213e-03
#> B-C.437          Gna12      B-C    2.020724  1.014872   3.191070 7.973295e-03
#> B-C.438          Inpp1      B-C  -12.643757 -3.660353  -3.179970 8.137451e-03
#> B-C.439          Naa15      B-C    4.307975  2.107010   3.173310 8.237573e-03
#> B-C.440         Agpat3      B-C   -2.806988 -1.489023  -3.171157 8.270197e-03
#> B-C.441          Cdip1      B-C    2.931383  1.551581   3.158965 8.457466e-03
#> B-C.442        Smarca2      B-C   -3.062142 -1.614541  -3.155440 8.512392e-03
#> B-C.443           Pld4      B-C    5.282471  2.401213   3.154783 8.522680e-03
#> B-C.444         Apol7e      B-C   -9.272105 -3.212897  -3.154750 8.523195e-03
#> B-C.445           Icmt      B-C   -2.936740 -1.554216  -3.151102 8.580489e-03
#> B-C.446          H2-Aa      B-C   -2.247877 -1.168563  -3.148010 8.629360e-03
#> B-C.447          Anks1      B-C    7.428560  2.893083   3.142189 8.722147e-03
#> B-C.448          Tor4a      B-C   -2.026237 -1.018803  -3.139194 8.770265e-03
#> B-C.449         Tspan2      B-C   23.313947  4.543121   3.137285 8.801075e-03
#> B-C.450          Ngly1      B-C   -2.145886 -1.101574  -3.135682 8.827038e-03
#> B-C.451          Arrb1      B-C  -22.906316 -4.517674  -3.132783 8.874170e-03
#> B-C.452           Lax1      B-C   -3.024990 -1.596930  -3.124328 9.013110e-03
#> B-C.453         Rnf122      B-C    6.082409  2.604643   3.122743 9.039407e-03
#> B-C.454        Cdk2ap2      B-C   -2.256208 -1.173900  -3.120885 9.070323e-03
#> B-C.455          Ssbp3      B-C   -2.064577 -1.045846  -3.120414 9.078175e-03
#> B-C.456       Trappc13      B-C    2.288753  1.194562   3.120320 9.079746e-03
#> B-C.457         Heatr1      B-C    2.743225  1.455873   3.117932 9.119669e-03
#> B-C.458  0610010F05Rik      B-C    6.804293  2.766445   3.115300 9.163897e-03
#> B-C.459          Nr3c1      B-C   -2.540979 -1.345384  -3.115151 9.166402e-03
#> B-C.460          Sgpl1      B-C   -3.157546 -1.658804  -3.114160 9.183110e-03
#> B-C.461         Scarb1      B-C    3.543987  1.825373   3.112617 9.209190e-03
#> B-C.462           Tle6      B-C   30.998338  4.954119   3.111335 9.230909e-03
#> B-C.463          Hipk2      B-C    3.351371  1.744751   3.111329 9.231014e-03
#> B-C.464         Parp10      B-C   -4.091218 -2.032530  -3.109745 9.257923e-03
#> B-C.465           Gatb      B-C    7.580396  2.922273   3.109610 9.260216e-03
#> B-C.466          Mapk6      B-C    2.720018  1.443616   3.106349 9.315891e-03
#> B-C.467        Tmem206      B-C    4.889249  2.289613   3.105005 9.338934e-03
#> B-C.468           Mfn1      B-C  -32.882947 -5.039268  -3.104768 9.342995e-03
#> B-C.469         Plxnc1      B-C    4.275110  2.095962   3.101782 9.394424e-03
#> B-C.470         Calcrl      B-C    2.439136  1.286370   3.100547 9.415764e-03
#> B-C.471          Casc3      B-C    2.904242  1.538162   3.097894 9.461795e-03
#> B-C.472        Pglyrp2      B-C   -3.990028 -1.996399  -3.095069 9.511049e-03
#> B-C.473          Spidr      B-C    4.215477  2.075696   3.094990 9.512442e-03
#> B-C.474          Ptcd1      B-C   10.314836  3.366649   3.092496 9.556152e-03
#> B-C.475            Cd9      B-C   -3.868526 -1.951784  -3.092178 9.561728e-03
#> B-C.476        Osbpl1a      B-C   16.536356  4.047569   3.091986 9.565110e-03
#> B-C.477           Pdcl      B-C   -2.128230 -1.089654  -3.091799 9.568395e-03
#> B-C.478            Eng      B-C   -2.751316 -1.460122  -3.091447 9.574589e-03
#> B-C.479         Tuba1a      B-C   -2.334946 -1.223389  -3.090686 9.587998e-03
#> B-C.480          Ppm1m      B-C    3.858623  1.948086   3.075136 9.866027e-03
#> B-C.481           Ezh1      B-C   -2.209540 -1.143746  -3.072174 9.919903e-03
#> B-C.482          Nup43      B-C    4.784912  2.258492   3.070265 9.954789e-03
#> B-C.483          Unc50      B-C    2.497842  1.320682   3.070004 9.959563e-03
#> B-C.484         Klhl42      B-C   -3.566587 -1.834544  -3.068272 9.991341e-03
#> B-C.485          Prkcd      B-C    2.180234  1.124483   3.062020 1.010684e-02
#> B-C.486         Rfxank      B-C   -5.441477 -2.443998  -3.061593 1.011478e-02
#> B-C.487        Rhobtb3      B-C  -43.505984 -5.443142  -3.060019 1.014410e-02
#> B-C.488         Gimap6      B-C   -2.353051 -1.234532  -3.055746 1.022412e-02
#> B-C.489         Dcaf13      B-C    2.157846  1.109592   3.055670 1.022554e-02
#> B-C.490           Bin3      B-C    2.332883  1.222114   3.054258 1.025214e-02
#> B-C.491         Tuba1c      B-C   -2.196873 -1.135452  -3.050799 1.031755e-02
#> B-C.492         Zfp629      B-C    4.801395  2.263454   3.050443 1.032431e-02
#> B-C.493         R3hdm1      B-C    3.081107  1.623449   3.050084 1.033111e-02
#> B-C.494      Rab11fip1      B-C   -8.363276 -3.064068  -3.048499 1.036127e-02
#> B-C.495          Sesn1      B-C   -2.219477 -1.150220  -3.048047 1.036989e-02
#> B-C.496         Phf21a      B-C  -17.011981 -4.088479  -3.046807 1.039356e-02
#> B-C.497          Pus7l      B-C    6.830962  2.772089   3.041731 1.049103e-02
#> B-C.498           Pcm1      B-C    3.763279  1.911990   3.041081 1.050359e-02
#> B-C.499           Ier5      B-C   -9.700434 -3.278049  -3.039546 1.053327e-02
#> B-C.500         Pea15a      B-C  -10.343830 -3.370699  -3.032791 1.066494e-02
#> B-C.501          Spsb3      B-C  -10.773822 -3.429458  -3.032487 1.067090e-02
#> B-C.502         Rab27a      B-C    2.106872  1.075103   3.025761 1.080373e-02
#> B-C.503           Gps2      B-C    3.195652  1.676110   3.025266 1.081355e-02
#> B-C.504           Mios      B-C    3.626823  1.858706   3.025151 1.081586e-02
#> B-C.505           Zbp1      B-C  -21.772631 -4.444444  -3.023517 1.084841e-02
#> B-C.506          Srsf2      B-C    2.058615  1.041674   3.021663 1.088546e-02
#> B-C.507         Zfp808      B-C   21.751622  4.443051   3.018793 1.094306e-02
#> B-C.508           Ly6a      B-C  -21.609783 -4.433613  -3.016089 1.099762e-02
#> B-C.509          Etfdh      B-C   -2.183560 -1.126682  -3.013004 1.106019e-02
#> B-C.510         Zcchc2      B-C    5.946836  2.572122   3.010555 1.111013e-02
#> B-C.511          Actr6      B-C   -4.451971 -2.154444  -3.008503 1.115215e-02
#> B-C.512           Cast      B-C   -2.331861 -1.221482  -3.004410 1.123641e-02
#> B-C.513           Sord      B-C    4.028207  2.010138   3.002261 1.128091e-02
#> B-C.514        Mettl13      B-C    3.607953  1.851180   3.000096 1.132591e-02
#> B-C.515           Ctsz      B-C    2.130540  1.091219   2.999574 1.133679e-02
#> B-C.516        Gramd1b      B-C   -3.742390 -1.903960  -2.996036 1.141081e-02
#> B-C.517          Parp8      B-C    2.187233  1.129107   2.993446 1.146529e-02
#> B-C.518        Irf2bpl      B-C  -14.739220 -3.881588  -2.990434 1.152899e-02
#> B-C.519         Pik3ca      B-C   -2.857501 -1.514754  -2.985021 1.164435e-02
#> B-C.520         Atp1a3      B-C   11.391535  3.509890   2.972546 1.191462e-02
#> B-C.521          Senp2      B-C    2.995883  1.582981   2.968753 1.199805e-02
#> B-C.522         Rassf4      B-C    5.984755  2.581292   2.968668 1.199991e-02
#> B-C.523        Heatr5b      B-C    6.497311  2.699843   2.968564 1.200221e-02
#> B-C.524          Trpv2      B-C   -2.862004 -1.517026  -2.967285 1.203049e-02
#> B-C.525          S1pr3      B-C    8.555851  3.096911   2.965694 1.206574e-02
#> B-C.526           Crem      B-C  -19.748420 -4.303665  -2.962299 1.214132e-02
#> B-C.527          Enpp5      B-C  -20.418116 -4.351778  -2.961983 1.214837e-02
#> B-C.528          G6pc3      B-C    2.716828  1.441923   2.960882 1.217300e-02
#> B-C.529          Dhx32      B-C   -2.553311 -1.352370  -2.959713 1.219922e-02
#> B-C.530       Mapkapk3      B-C    3.350041  1.744179   2.957684 1.224483e-02
#> B-C.531          Sesn2      B-C    5.915900  2.564598   2.957439 1.225034e-02
#> B-C.532          Fads2      B-C    7.172988  2.842574   2.957125 1.225743e-02
#> B-C.533          Cytip      B-C   -2.483686 -1.312483  -2.954332 1.232056e-02
#> B-C.534         Pcmtd1      B-C   -2.910264 -1.541150  -2.953696 1.233497e-02
#> B-C.535         Gatsl2      B-C   -4.825424 -2.270656  -2.953677 1.233541e-02
#> B-C.536           Tia1      B-C    2.373547  1.247044   2.952749 1.235649e-02
#> B-C.537          Mocs1      B-C   -9.513843 -3.250028  -2.950804 1.240076e-02
#> B-C.538           Acy1      B-C    5.832538  2.544124   2.949666 1.242675e-02
#> B-C.539       Zc3hav1l      B-C    5.805765  2.537486   2.940491 1.263823e-02
#> B-C.540          Hspe1      B-C    2.055287  1.039340   2.939997 1.264972e-02
#> B-C.541  3110082I17Rik      B-C   12.693801  3.666052   2.929925 1.288624e-02
#> B-C.542          Cyth4      B-C    2.021193  1.015207   2.929289 1.290131e-02
#> B-C.543         Atp2b1      B-C   -3.459739 -1.790663  -2.927780 1.293717e-02
#> B-C.544         Pi4k2b      B-C    6.678741  2.739576   2.927272 1.294926e-02
#> B-C.545           Plek      B-C    7.198470  2.847690   2.925414 1.299358e-02
#> B-C.546        Akirin2      B-C   -2.098302 -1.069222  -2.925163 1.299960e-02
#> B-C.547          Lzts2      B-C    5.084571  2.346126   2.924214 1.302229e-02
#> B-C.548          Lsm11      B-C  -16.359660 -4.032071  -2.923672 1.303530e-02
#> B-C.549           Gfm1      B-C    2.195077  1.134271   2.923616 1.303664e-02
#> B-C.550         Trafd1      B-C    2.082919  1.058607   2.921085 1.309745e-02
#> B-C.551           Mmab      B-C   -3.705514 -1.889674  -2.913229 1.328806e-02
#> B-C.552           Snx9      B-C   12.651585  3.661246   2.912286 1.331113e-02
#> B-C.553        Map3k14      B-C   -6.796055 -2.764697  -2.910687 1.335033e-02
#> B-C.554           Smc5      B-C   -2.105645 -1.074263  -2.904541 1.350207e-02
#> B-C.555          Nxpe3      B-C    4.476550  2.162387   2.902219 1.355985e-02
#> B-C.556       Arhgap25      B-C   -3.929542 -1.974361  -2.901420 1.357980e-02
#> B-C.557            Tec      B-C    3.438442  1.781755   2.900490 1.360304e-02
#> B-C.558           Lgmn      B-C   -6.057130 -2.598634  -2.898527 1.365223e-02
#> B-C.559          Rad18      B-C    3.625057  1.858004   2.897758 1.367156e-02
#> B-C.560        Sec14l1      B-C   -2.705014 -1.435636  -2.895184 1.373642e-02
#> B-C.561          Capn5      B-C   11.804139  3.561221   2.894990 1.374132e-02
#> B-C.562           Ctsf      B-C   -6.167440 -2.624672  -2.892291 1.380968e-02
#> B-C.563           Vwa8      B-C    5.188759  2.375389   2.889489 1.388104e-02
#> B-C.564           Amz2      B-C    5.215872  2.382908   2.889301 1.388583e-02
#> B-C.565        Aldh3a2      B-C   -3.015625 -1.592457  -2.889199 1.388843e-02
#> B-C.566          Oasl2      B-C  -16.985707 -4.086249  -2.889080 1.389146e-02
#> B-C.567        Slc33a1      B-C  -11.922053 -3.575561  -2.888593 1.390392e-02
#> B-C.568           Dok1      B-C   -4.964856 -2.311752  -2.888201 1.391394e-02
#> B-C.569         Dock11      B-C    2.616011  1.387369   2.883759 1.402805e-02
#> B-C.570        Adprhl2      B-C   -2.400376 -1.263260  -2.880816 1.410417e-02
#> B-C.571         Neurl3      B-C  -27.449559 -4.778711  -2.880304 1.411746e-02
#> B-C.572           Klf6      B-C   -2.152613 -1.106089  -2.874389 1.427182e-02
#> B-C.573         Tespa1      B-C   -2.380230 -1.251101  -2.867373 1.445710e-02
#> B-C.574          Padi4      B-C   16.454135  4.040378   2.867211 1.446141e-02
#> B-C.575          Mlycd      B-C    4.056187  2.020124   2.861861 1.460435e-02
#> B-C.576        Arl6ip6      B-C    5.299301  2.405802   2.855950 1.476392e-02
#> B-C.577        Angptl4      B-C    9.072015  3.181423   2.855796 1.476810e-02
#> B-C.578          Cisd3      B-C   -4.345116 -2.119395  -2.855700 1.477071e-02
#> B-C.579            Gsn      B-C    2.250617  1.170321   2.853461 1.483163e-02
#> B-C.580           Msra      B-C    2.110376  1.077500   2.849672 1.493528e-02
#> B-C.581          Sirt3      B-C    2.437363  1.285321   2.848002 1.498119e-02
#> B-C.582           Crot      B-C   -2.161768 -1.112211  -2.844489 1.507824e-02
#> B-C.583            Nmi      B-C   -3.729699 -1.899059  -2.840240 1.519645e-02
#> B-C.584           Ibtk      B-C    2.641217  1.401203   2.835832 1.532006e-02
#> B-C.585         Zbtb11      B-C    5.797753  2.535494   2.832083 1.542600e-02
#> B-C.586           Fdps      B-C    2.224931  1.153761   2.831000 1.545672e-02
#> B-C.587         Incenp      B-C   -2.145830 -1.101536  -2.829283 1.550558e-02
#> B-C.588         Lgals8      B-C   -3.270456 -1.709492  -2.825789 1.560545e-02
#> B-C.589           Pld3      B-C   -5.116417 -2.355134  -2.823510 1.567092e-02
#> B-C.590           Utrn      B-C  -12.579210 -3.652969  -2.820430 1.575988e-02
#> B-C.591          Zfp68      B-C    4.048546  2.017404   2.818495 1.581600e-02
#> B-C.592        Pcyox1l      B-C    6.243959  2.642461   2.818205 1.582443e-02
#> B-C.593           Ing4      B-C    2.739710  1.454023   2.817094 1.585678e-02
#> B-C.594          Rcor2      B-C    5.144226  2.362954   2.813560 1.596007e-02
#> B-C.595        Gm10384      B-C   22.986933  4.522742   2.807504 1.613863e-02
#> B-C.596          Mki67      B-C   -2.655325 -1.408888  -2.806067 1.618129e-02
#> B-C.597           Acp1      B-C    2.253455  1.172138   2.805916 1.618579e-02
#> B-C.598         Nos1ap      B-C   11.660099  3.543508   2.805155 1.620844e-02
#> B-C.599         Pgm2l1      B-C   -9.012980 -3.172004  -2.805126 1.620929e-02
#> B-C.600  9130401M01Rik      B-C   -2.448860 -1.292110  -2.802744 1.628038e-02
#> B-C.601         Bcl2l1      B-C   -5.731367 -2.518879  -2.802030 1.630175e-02
#> B-C.602           Zzz3      B-C    2.169165  1.117140   2.800982 1.633315e-02
#> B-C.603        Rapgef6      B-C   -2.367114 -1.243129  -2.798690 1.640207e-02
#> B-C.604           Pisd      B-C    2.886099  1.529121   2.798630 1.640388e-02
#> B-C.605        Psmc3ip      B-C    4.865725  2.282655   2.796953 1.645446e-02
#> B-C.606          Acsl4      B-C   -2.020917 -1.015010  -2.793111 1.657101e-02
#> B-C.607        Slc44a1      B-C   -3.814843 -1.931624  -2.790113 1.666250e-02
#> B-C.608        Arhgap4      B-C   -2.585385 -1.370379  -2.788870 1.670058e-02
#> B-C.609         Lztfl1      B-C   -7.936827 -2.988562  -2.788414 1.671456e-02
#> B-C.610      Slfn10-ps      B-C   -8.379265 -3.066824  -2.787497 1.674276e-02
#> B-C.611          Abcd3      B-C   -2.875752 -1.523939  -2.783471 1.686699e-02
#> B-C.612         Gpr160      B-C    3.625761  1.858284   2.780501 1.695922e-02
#> B-C.613            Lxn      B-C   -5.506642 -2.461173  -2.779022 1.700535e-02
#> B-C.614          Cdca2      B-C   -2.256611 -1.174158  -2.778408 1.702455e-02
#> B-C.615         Hmgxb4      B-C    4.151171  2.053518   2.776807 1.707467e-02
#> B-C.616         Zfp516      B-C   -7.217975 -2.851594  -2.775449 1.711728e-02
#> B-C.617          Pde6d      B-C  -12.129732 -3.600476  -2.771805 1.723219e-02
#> B-C.618        Slc38a7      B-C  -22.811743 -4.511705  -2.770976 1.725844e-02
#> B-C.619           Hmbs      B-C    2.705588  1.435942   2.769674 1.729974e-02
#> B-C.620          Cers5      B-C   -2.360530 -1.239111  -2.769096 1.731813e-02
#> B-C.621          Rnf44      B-C    2.258276  1.175222   2.767199 1.737853e-02
#> B-C.622  1110059G10Rik      B-C   -4.489990 -2.166712  -2.765030 1.744788e-02
#> B-C.623          Ehmt2      B-C   -3.036440 -1.602381  -2.753712 1.781415e-02
#> B-C.624          Adpgk      B-C   -2.373015 -1.246721  -2.753127 1.783329e-02
#> B-C.625         Trim44      B-C    2.664708  1.413978   2.752012 1.786979e-02
#> B-C.626           Bcor      B-C    3.207447  1.681426   2.751970 1.787119e-02
#> B-C.627          Egln3      B-C    3.884552  1.957748   2.750200 1.792934e-02
#> B-C.628         Fchsd1      B-C   -8.839925 -3.144034  -2.749614 1.794861e-02
#> B-C.629        Slc10a7      B-C  -15.172089 -3.923348  -2.747161 1.802958e-02
#> B-C.630         Qtrtd1      B-C  -14.998328 -3.906730  -2.745172 1.809553e-02
#> B-C.631          Lmtk2      B-C   -3.575700 -1.838226  -2.744668 1.811228e-02
#> B-C.632         Gtf3c2      B-C    3.884832  1.957852   2.743372 1.815538e-02
#> B-C.633          Casd1      B-C    2.611960  1.385133   2.742845 1.817295e-02
#> B-C.634          Tdrd7      B-C  -21.366892 -4.417305  -2.738875 1.830581e-02
#> B-C.635          Phf20      B-C   -2.509290 -1.327279  -2.735737 1.841151e-02
#> B-C.636         Med12l      B-C   13.521722  3.757207   2.735030 1.843543e-02
#> B-C.637          Ncor1      B-C   -2.009902 -1.007125  -2.730715 1.858191e-02
#> B-C.638         Tbc1d5      B-C    2.862747  1.517400   2.723012 1.884632e-02
#> B-C.639         Erlin1      B-C    2.800856  1.485868   2.721761 1.888964e-02
#> B-C.640          Siah2      B-C   -7.023149 -2.812118  -2.721392 1.890241e-02
#> B-C.641           Rit1      B-C   -4.497012 -2.168967  -2.719963 1.895200e-02
#> B-C.642         Atp8b2      B-C   -4.575461 -2.193917  -2.717939 1.902250e-02
#> B-C.643          Exoc8      B-C   -3.538235 -1.823030  -2.717809 1.902702e-02
#> B-C.644       Tnfrsf18      B-C   -8.787357 -3.135429  -2.708794 1.934417e-02
#> B-C.645         Polr1b      B-C    3.264147  1.706706   2.708547 1.935292e-02
#> B-C.646          Mipep      B-C    3.694847  1.885515   2.706463 1.942701e-02
#> B-C.647           Tlk1      B-C   -2.026018 -1.018647  -2.705834 1.944942e-02
#> B-C.648          Uqcc1      B-C    2.329850  1.220237   2.703628 1.952823e-02
#> B-C.649         Mrps27      B-C   -2.679579 -1.422006  -2.698870 1.969931e-02
#> B-C.650           Lfng      B-C   -2.425565 -1.278321  -2.698790 1.970220e-02
#> B-C.651        Tmem192      B-C   -2.390291 -1.257186  -2.696910 1.977022e-02
#> B-C.652          Zbed4      B-C    5.773331  2.529404   2.696491 1.978539e-02
#> B-C.653         Insig2      B-C   -3.520739 -1.815878  -2.695348 1.982689e-02
#> B-C.654           Dna2      B-C    8.863618  3.147896   2.694289 1.986540e-02
#> B-C.655         Armcx5      B-C  -10.488543 -3.390742  -2.691047 1.998379e-02
#> B-C.656          Cela1      B-C   -7.214579 -2.850915  -2.689491 2.004086e-02
#> B-C.657          Terf1      B-C   -3.686310 -1.882177  -2.687101 2.012882e-02
#> B-C.658          Dzip3      B-C    9.456177  3.241257   2.684798 2.021392e-02
#> B-C.659          P4ha1      B-C   -3.132897 -1.647497  -2.683926 2.024628e-02
#> B-C.660         Pik3cg      B-C   -2.099414 -1.069987  -2.682826 2.028709e-02
#> B-C.661  2610008E11Rik      B-C   11.111765  3.474016   2.680309 2.038086e-02
#> B-C.662            Aen      B-C    3.925332  1.972815   2.676306 2.053087e-02
#> B-C.663         Rassf2      B-C   -2.215320 -1.147515  -2.676241 2.053332e-02
#> B-C.664          Sart3      B-C    2.857977  1.514994   2.674378 2.060353e-02
#> B-C.665           Dntt      B-C    8.343972  3.060734   2.672167 2.068712e-02
#> B-C.666      Eif4enif1      B-C    2.024983  1.017910   2.665976 2.092302e-02
#> B-C.667          Stag1      B-C    2.334709  1.223243   2.665860 2.092746e-02
#> B-C.668         Gnptab      B-C    2.232590  1.158718   2.665192 2.095308e-02
#> B-C.669           Pigu      B-C    3.452972  1.787839   2.664885 2.096488e-02
#> B-C.670        Slc35f2      B-C   -8.720143 -3.124352  -2.664105 2.099484e-02
#> B-C.671         Ifitm2      B-C    4.713413  2.236772   2.660908 2.111810e-02
#> B-C.672           Nck2      B-C   -5.062391 -2.339819  -2.660737 2.112471e-02
#> B-C.673          Grina      B-C   -4.388943 -2.133873  -2.656001 2.130866e-02
#> B-C.674         Atp8a1      B-C   -5.771193 -2.528869  -2.654779 2.135641e-02
#> B-C.675          Itpr3      B-C   -3.795936 -1.924456  -2.653485 2.140705e-02
#> B-C.676          Golm1      B-C    2.254520  1.172820   2.653331 2.141308e-02
#> B-C.677         Map4k5      B-C   -6.553636 -2.712296  -2.652662 2.143931e-02
#> B-C.678           Gnl3      B-C    2.648472  1.405160   2.651870 2.147043e-02
#> B-C.679          Mrpl1      B-C    5.612293  2.488590   2.651710 2.147673e-02
#> B-C.680         Zfp213      B-C   -6.178500 -2.627257  -2.644231 2.177271e-02
#> B-C.681         Zfp398      B-C   -7.496938 -2.906301  -2.644051 2.177987e-02
#> B-C.682          Abcb9      B-C   -3.454098 -1.788309  -2.642743 2.183209e-02
#> B-C.683         Chchd4      B-C    2.341670  1.227538   2.642438 2.184429e-02
#> B-C.684          Itpr2      B-C   -3.554712 -1.829733  -2.639706 2.195377e-02
#> B-C.685          Prdx4      B-C   -6.681762 -2.740229  -2.638801 2.199014e-02
#> B-C.686          Plod3      B-C   19.528562  4.287514   2.637125 2.205768e-02
#> B-C.687  4833439L19Rik      B-C    2.608272  1.383094   2.633259 2.221423e-02
#> B-C.688          Rab4a      B-C   -8.276824 -3.049077  -2.632650 2.223900e-02
#> B-C.689          Ccnt1      B-C   -2.025298 -1.018134  -2.632261 2.225482e-02
#> B-C.690         Zfp319      B-C   -6.319052 -2.659708  -2.631475 2.228685e-02
#> B-C.691        Prkar2b      B-C    5.401429  2.433341   2.628804 2.239601e-02
#> B-C.692         Rad54l      B-C    3.792153  1.923017   2.628794 2.239642e-02
#> B-C.693          Arntl      B-C   -4.689624 -2.229472  -2.625774 2.252048e-02
#> B-C.694          Pear1      B-C   -2.488108 -1.315049  -2.625294 2.254025e-02
#> B-C.695          Uckl1      B-C    4.360060  2.124348   2.625144 2.254642e-02
#> B-C.696         Slain1      B-C   -2.206845 -1.141985  -2.624661 2.256634e-02
#> B-C.697        Rasgrp2      B-C    2.477474  1.308870   2.624618 2.256813e-02
#> B-C.698          Sfxn1      B-C    2.120478  1.084389   2.623064 2.263233e-02
#> B-C.699          Tasp1      B-C   -5.155765 -2.366187  -2.621056 2.271561e-02
#> B-C.700         Mrpl57      B-C    2.374644  1.247711   2.619952 2.276149e-02
#> B-C.701           Urb2      B-C    3.511270  1.811993   2.617818 2.285048e-02
#> B-C.702        Pitpnm1      B-C   -5.597450 -2.484770  -2.610883 2.314203e-02
#> B-C.703         Pik3r1      B-C   -2.235226 -1.160421  -2.598383 2.367678e-02
#> B-C.704          Ints4      B-C    2.054329  1.038667   2.594049 2.386501e-02
#> B-C.705         Ociad2      B-C    3.646774  1.866621   2.592815 2.391886e-02
#> B-C.706       Slc4a1ap      B-C   -3.236656 -1.694504  -2.592755 2.392147e-02
#> B-C.707        Fam120b      B-C   -2.688714 -1.426917  -2.586625 2.419084e-02
#> B-C.708          Clcn3      B-C   -4.256261 -2.089587  -2.586259 2.420698e-02
#> B-C.709          Cenpl      B-C   -2.360621 -1.239167  -2.580300 2.447182e-02
#> B-C.710          Gng12      B-C    2.237222  1.161709   2.579930 2.448835e-02
#> B-C.711         Eif2b5      B-C    2.055170  1.039258   2.579874 2.449086e-02
#> B-C.712          Trub2      B-C   12.389175  3.631008   2.578658 2.454532e-02
#> B-C.713  1110008L16Rik      B-C    6.254396  2.644871   2.576425 2.464555e-02
#> B-C.714           Dlg1      B-C   -2.968602 -1.569784  -2.574900 2.471426e-02
#> B-C.715            Mmd      B-C   -3.582738 -1.841063  -2.571970 2.484681e-02
#> B-C.716            Mvd      B-C   -9.359795 -3.226477  -2.570871 2.489668e-02
#> B-C.717          Gins3      B-C   -6.288281 -2.652666  -2.569578 2.495553e-02
#> B-C.718           Ccnh      B-C    2.010602  1.007628   2.569033 2.498037e-02
#> B-C.719          Lrch1      B-C   -2.869321 -1.520709  -2.568181 2.501922e-02
#> B-C.720          Thoc6      B-C    2.868391  1.520242   2.568057 2.502488e-02
#> B-C.721           Tmc8      B-C  -13.937430 -3.800893  -2.567673 2.504244e-02
#> B-C.722         Ms4a4b      B-C   -2.933035 -1.552394  -2.567026 2.507203e-02
#> B-C.723        Tmem39b      B-C    3.869106  1.952000   2.565777 2.512924e-02
#> B-C.724       Vkorc1l1      B-C    4.331370  2.114824   2.563661 2.522644e-02
#> B-C.725         Mllt10      B-C    4.182977  2.064530   2.562307 2.528882e-02
#> B-C.726           Abi2      B-C   -3.090888 -1.628021  -2.555391 2.560987e-02
#> B-C.727          Nlrx1      B-C   -4.796688 -2.262039  -2.553829 2.568293e-02
#> B-C.728         Zbtb22      B-C   -3.799634 -1.925860  -2.553515 2.569765e-02
#> B-C.729     Csgalnact2      B-C    5.239477  2.389423   2.553163 2.571418e-02
#> B-C.730          Cenph      B-C    3.701355  1.888054   2.553095 2.571734e-02
#> B-C.731          Dtwd1      B-C   -4.750370 -2.248040  -2.552537 2.574354e-02
#> B-C.732          Nfkb2      B-C   -7.075421 -2.822816  -2.551699 2.578291e-02
#> B-C.733           Ppcs      B-C   -2.543629 -1.346888  -2.551438 2.579521e-02
#> B-C.734          Fbxw4      B-C    4.543095  2.183675   2.551232 2.580486e-02
#> B-C.735        Tbc1d19      B-C   -5.067924 -2.341395  -2.550328 2.584743e-02
#> B-C.736           Numb      B-C   -9.067969 -3.180779  -2.549516 2.588575e-02
#> B-C.737         Elmod3      B-C   -4.276996 -2.096598  -2.548918 2.591398e-02
#> B-C.738         Cdkn2c      B-C   -6.210087 -2.634614  -2.548441 2.593656e-02
#> B-C.739          Wdr34      B-C   11.386887  3.509302   2.548048 2.595512e-02
#> B-C.740         Camk1d      B-C   -2.090871 -1.064104  -2.545244 2.608818e-02
#> B-C.741         Dgcr14      B-C   -2.508384 -1.326758  -2.542913 2.619928e-02
#> B-C.742        Hnrnpll      B-C   -2.202152 -1.138914  -2.538776 2.639756e-02
#> B-C.743         Lrrc8b      B-C   -4.331837 -2.114979  -2.538458 2.641288e-02
#> B-C.744            Ahr      B-C   -2.748316 -1.458548  -2.537187 2.647412e-02
#> B-C.745        Zcchc18      B-C   -4.337832 -2.116974  -2.534274 2.661503e-02
#> B-C.746          Qtrt1      B-C    4.552944  2.186800   2.531667 2.674179e-02
#> B-C.747         Zfp182      B-C   16.365464  4.032583   2.531545 2.674772e-02
#> B-C.748           Pkp3      B-C   -3.640303 -1.864058  -2.528104 2.691594e-02
#> B-C.749          Alcam      B-C   25.643506  4.680522   2.523814 2.712707e-02
#> B-C.750          Krit1      B-C   -2.205841 -1.141329  -2.521797 2.722690e-02
#> B-C.751          Xlr4b      B-C   -5.924194 -2.566619  -2.518248 2.740343e-02
#> B-C.752          Gna13      B-C   -2.045242 -1.032271  -2.515505 2.754068e-02
#> B-C.753           Eri2      B-C   -8.496445 -3.086859  -2.514961 2.756797e-02
#> B-C.754          P2rx4      B-C    2.520190  1.333533   2.514061 2.761318e-02
#> B-C.755          Anxa4      B-C   -3.185898 -1.671700  -2.513188 2.765708e-02
#> B-C.756           Zeb2      B-C    8.038041  3.006844   2.512996 2.766676e-02
#> B-C.757          Hmces      B-C    3.412842  1.770974   2.510711 2.778205e-02
#> B-C.758          Syce2      B-C    5.569943  2.477663   2.510532 2.779113e-02
#> B-C.759           Als2      B-C   13.407083  3.744924   2.510364 2.779959e-02
#> B-C.760         Rnf157      B-C    2.362678  1.240423   2.509997 2.781817e-02
#> B-C.761           Gatm      B-C   12.439546  3.636862   2.504884 2.807827e-02
#> B-C.762           Mycl      B-C   10.940074  3.451551   2.504115 2.811755e-02
#> B-C.763           Gnal      B-C   -4.269140 -2.093945  -2.502702 2.818996e-02
#> B-C.764           Sik1      B-C    4.413106  2.141794   2.498421 2.841038e-02
#> B-C.765           Rcl1      B-C    3.035806  1.602080   2.498192 2.842219e-02
#> B-C.766           Smox      B-C   -2.339888 -1.226439  -2.495741 2.854921e-02
#> B-C.767         Mboat1      B-C   -6.818880 -2.769535  -2.494945 2.859056e-02
#> B-C.768        Rad54l2      B-C   -4.783496 -2.258065  -2.493005 2.869161e-02
#> B-C.769        Rps6ka3      B-C   -3.821803 -1.934253  -2.492763 2.870423e-02
#> B-C.770          Ift57      B-C    3.913200  1.968349   2.489245 2.888846e-02
#> B-C.771          Ints3      B-C   -2.665066 -1.414172  -2.488757 2.891411e-02
#> B-C.772          Ddx17      B-C    2.282521  1.190628   2.483215 2.920690e-02
#> B-C.773        Zdhhc15      B-C    8.420827  3.073962   2.482403 2.925003e-02
#> B-C.774         Zfp275      B-C   -3.425063 -1.776131  -2.482229 2.925933e-02
#> B-C.775            Pxk      B-C   -5.170703 -2.370360  -2.480776 2.933669e-02
#> B-C.776          Zmym3      B-C    2.478226  1.309308   2.478385 2.946445e-02
#> B-C.777          Adck3      B-C   -2.425139 -1.278067  -2.477699 2.950119e-02
#> B-C.778        Depdc1b      B-C   -2.383964 -1.253363  -2.475736 2.960661e-02
#> B-C.779           Clk4      B-C   -5.154297 -2.365776  -2.474458 2.967542e-02
#> B-C.780           Mcat      B-C    2.030296  1.021690   2.469592 2.993894e-02
#> B-C.781         Dolpp1      B-C   -4.456512 -2.155915  -2.468725 2.998608e-02
#> B-C.782         Sapcd2      B-C   -9.391967 -3.231427  -2.467401 3.005829e-02
#> B-C.783          Nr2c1      B-C    6.344441  2.665493   2.462825 3.030913e-02
#> B-C.784          Btbd6      B-C   -5.354573 -2.420772  -2.462218 3.034253e-02
#> B-C.785          Scmh1      B-C    4.123488  2.043865   2.459838 3.047396e-02
#> B-C.786          Myo1c      B-C    3.577493  1.838949   2.457963 3.057787e-02
#> B-C.787       Timeless      B-C    2.375529  1.248249   2.457589 3.059864e-02
#> B-C.788           Ltbr      B-C    4.271951  2.094895   2.457545 3.060107e-02
#> B-C.789          Phka2      B-C    2.704337  1.435275   2.453104 3.084869e-02
#> B-C.790            Gem      B-C   -8.753210 -3.129812  -2.452170 3.090104e-02
#> B-C.791         Fam13b      B-C   -2.308223 -1.206782  -2.452155 3.090190e-02
#> B-C.792          Mmgt1      B-C    3.499240  1.807042   2.451656 3.092990e-02
#> B-C.793         Sft2d1      B-C    3.378670  1.756456   2.450524 3.099346e-02
#> B-C.794  6030458C11Rik      B-C    5.228134  2.386296   2.447520 3.116286e-02
#> B-C.795        Arfgap3      B-C   -5.578456 -2.479866  -2.446849 3.120081e-02
#> B-C.796          Mex3c      B-C    3.375110  1.754935   2.444820 3.131586e-02
#> B-C.797          Arl11      B-C   -6.456885 -2.690838  -2.444106 3.135644e-02
#> B-C.798            Lnp      B-C   -3.589693 -1.843861  -2.441753 3.149060e-02
#> B-C.799         Abhd15      B-C   -3.517664 -1.814618  -2.439795 3.160257e-02
#> B-C.800          P2rx7      B-C    2.115648  1.081099   2.437164 3.175371e-02
#> B-C.801        Ndufaf3      B-C   -2.881445 -1.526792  -2.434890 3.188488e-02
#> B-C.802           Tep1      B-C    9.360297  3.226554   2.433844 3.194542e-02
#> B-C.803          Ypel1      B-C   11.713876  3.550147   2.432932 3.199825e-02
#> B-C.804          Hyal2      B-C    6.206706  2.633828   2.431670 3.207156e-02
#> B-C.805          Itih5      B-C    3.157821  1.658929   2.431550 3.207857e-02
#> B-C.806          Csf1r      B-C   13.702184  3.776334   2.431406 3.208691e-02
#> B-C.807          Prdm4      B-C    2.919961  1.545949   2.430877 3.211771e-02
#> B-C.808         Zfp142      B-C   -3.320660 -1.731470  -2.424098 3.251466e-02
#> B-C.809         Scamp1      B-C    4.524531  2.177768   2.412011 3.323430e-02
#> B-C.810          Stk10      B-C   -2.212238 -1.145507  -2.409493 3.338613e-02
#> B-C.811         Il1rap      B-C   15.827473  3.984359   2.408862 3.342430e-02
#> B-C.812         Ssx2ip      B-C   -3.461175 -1.791262  -2.407568 3.350267e-02
#> B-C.813         Lysmd4      B-C   -7.047036 -2.817017  -2.406338 3.357734e-02
#> B-C.814          Paqr8      B-C    5.475654  2.453031   2.405013 3.365798e-02
#> B-C.815       Slc25a16      B-C   -3.375504 -1.755103  -2.404836 3.366875e-02
#> B-C.816          Ocel1      B-C   -6.274058 -2.649399  -2.403676 3.373948e-02
#> B-C.817          Mpeg1      B-C    7.663093  2.937927   2.403253 3.376534e-02
#> B-C.818          Itfg3      B-C   -2.054905 -1.039072  -2.402954 3.378359e-02
#> B-C.819           Igtp      B-C   -2.256521 -1.174100  -2.401655 3.386306e-02
#> B-C.820        Zfp322a      B-C   -2.246013 -1.167366  -2.401436 3.387651e-02
#> B-C.821         Clec2i      B-C   -2.821445 -1.496434  -2.395458 3.424481e-02
#> B-C.822          Ap1g2      B-C   -2.088635 -1.062560  -2.393470 3.436814e-02
#> B-C.823           Bms1      B-C    3.108698  1.636310   2.392305 3.444067e-02
#> B-C.824          Thada      B-C   -3.489633 -1.803075  -2.387962 3.471212e-02
#> B-C.825          Ncoa7      B-C   -5.874060 -2.554358  -2.386842 3.478247e-02
#> B-C.826         Tbc1d7      B-C    6.886079  2.783683   2.385867 3.484384e-02
#> B-C.827         Nusap1      B-C   -2.092114 -1.064961  -2.384656 3.492017e-02
#> B-C.828        Dennd2d      B-C   -4.261183 -2.091254  -2.384073 3.495699e-02
#> B-C.829         Oxnad1      B-C   15.534884  3.957440   2.383486 3.499409e-02
#> B-C.830           Sypl      B-C    2.338354  1.225493   2.382746 3.504096e-02
#> B-C.831          Dnm1l      B-C    2.109357  1.076803   2.381665 3.510948e-02
#> B-C.832           Emc8      B-C    2.017864  1.012829   2.380583 3.517818e-02
#> B-C.833         Fbxo42      B-C    3.060682  1.613853   2.379992 3.521576e-02
#> B-C.834           Crat      B-C    3.516666  1.814208   2.379879 3.522291e-02
#> B-C.835          Hace1      B-C   11.333524  3.502525   2.378164 3.533227e-02
#> B-C.836            F8a      B-C   -3.696993 -1.886352  -2.375665 3.549209e-02
#> B-C.837        Bcl2l11      B-C  -14.035368 -3.810995  -2.374790 3.554826e-02
#> B-C.838          Rab31      B-C    6.383018  2.674239   2.370256 3.584050e-02
#> B-C.839         Nudcd2      B-C    2.621808  1.390562   2.368641 3.594517e-02
#> B-C.840          Crtc3      B-C   -3.910248 -1.967260  -2.366434 3.608865e-02
#> B-C.841          Ttc19      B-C    4.674680  2.224868   2.366101 3.611036e-02
#> B-C.842          Eltd1      B-C   -2.277947 -1.187734  -2.364346 3.622491e-02
#> B-C.843          Nabp1      B-C   -3.925760 -1.972972  -2.362261 3.636152e-02
#> B-C.844           Ganc      B-C   -3.903702 -1.964843  -2.360382 3.648503e-02
#> B-C.845           Aff1      B-C    3.579273  1.839667   2.359927 3.651498e-02
#> B-C.846          Nsun4      B-C    5.479542  2.454055   2.358922 3.658124e-02
#> B-C.847          Mylip      B-C  -12.450846 -3.638172  -2.358298 3.662249e-02
#> B-C.848          Plod1      B-C   -3.307432 -1.725711  -2.356229 3.675945e-02
#> B-C.849          Sh2b3      B-C    2.305816  1.205277   2.355512 3.680698e-02
#> B-C.850          Msto1      B-C    2.275110  1.185937   2.354662 3.686346e-02
#> B-C.851         Zfp839      B-C   10.856826  3.440531   2.354109 3.690024e-02
#> B-C.852          Neil1      B-C   -5.329383 -2.413969  -2.353520 3.693949e-02
#> B-C.853          Stap1      B-C   -3.611091 -1.852435  -2.351440 3.707833e-02
#> B-C.854         Polr3h      B-C    2.712545  1.439647   2.348881 3.724975e-02
#> B-C.855         Smurf2      B-C   -8.210055 -3.037392  -2.347521 3.734126e-02
#> B-C.856            Htt      B-C    3.694788  1.885492   2.343267 3.762866e-02
#> B-C.857           Taf1      B-C   -2.743481 -1.456007  -2.341082 3.777710e-02
#> B-C.858           Ern1      B-C   -3.697775 -1.886657  -2.340834 3.779403e-02
#> B-C.859            Zak      B-C   11.289886  3.496959   2.340749 3.779981e-02
#> B-C.860         Man2c1      B-C    2.706664  1.436516   2.340185 3.783824e-02
#> B-C.861         Ifnar2      B-C    2.755631  1.462383   2.338754 3.793593e-02
#> B-C.862          Nim1k      B-C  -14.920826 -3.899255  -2.337993 3.798795e-02
#> B-C.863           Med1      B-C    3.243157  1.697399   2.336455 3.809340e-02
#> B-C.864          Rpl22      B-C    2.053503  1.038087   2.335978 3.812613e-02
#> B-C.865         Zfp599      B-C   -9.809824 -3.294227  -2.334441 3.823185e-02
#> B-C.866         Cdc25b      B-C   -3.367857 -1.751831  -2.333502 3.829656e-02
#> B-C.867          Rc3h2      B-C   -3.618776 -1.855502  -2.333266 3.831284e-02
#> B-C.868          Naa16      B-C    3.015665  1.592476   2.333064 3.832679e-02
#> B-C.869          Stim2      B-C   -2.196216 -1.135020  -2.332606 3.835837e-02
#> B-C.870       B3galnt2      B-C    3.214670  1.684671   2.332462 3.836835e-02
#> B-C.871            Tk2      B-C   -3.057732 -1.612462  -2.329356 3.858353e-02
#> B-C.872          Rint1      B-C   -4.539662 -2.182585  -2.329046 3.860505e-02
#> B-C.873          Herc3      B-C   -2.440975 -1.287458  -2.326634 3.877308e-02
#> B-C.874          Grwd1      B-C    2.921747  1.546831   2.326215 3.880230e-02
#> B-C.875           Klf9      B-C   -8.522317 -3.091246  -2.324862 3.889693e-02
#> B-C.876         Armc10      B-C    2.991848  1.581037   2.324585 3.891632e-02
#> B-C.877        Ubash3a      B-C    4.549161  2.185601   2.324553 3.891859e-02
#> B-C.878          Ascc3      B-C   -4.169970 -2.060037  -2.324536 3.891977e-02
#> B-C.879          Nat10      B-C    3.150767  1.655703   2.319984 3.923993e-02
#> B-C.880          Txlna      B-C    2.382806  1.252661   2.317139 3.944128e-02
#> B-C.881         Mfsd7b      B-C   -3.361878 -1.749267  -2.316095 3.951544e-02
#> B-C.882           Tjp3      B-C   -2.363306 -1.240806  -2.315491 3.955835e-02
#> B-C.883           Suox      B-C   -2.350043 -1.232687  -2.315382 3.956609e-02
#> B-C.884          Ptplb      B-C    2.905006  1.538541   2.313617 3.969195e-02
#> B-C.885          Alas1      B-C    2.519616  1.333204   2.313474 3.970215e-02
#> B-C.886         Tspan4      B-C    3.456466  1.789298   2.312291 3.978668e-02
#> B-C.887          Basp1      B-C   -3.431748 -1.778943  -2.310974 3.988102e-02
#> B-C.888          Rptor      B-C   -2.031312 -1.022412  -2.309659 3.997539e-02
#> B-C.889          Pole2      B-C    3.509435  1.811239   2.308443 4.006292e-02
#> B-C.890         Snapc2      B-C    3.669490  1.875580   2.306623 4.019416e-02
#> B-C.891          Bace1      B-C    2.808193  1.489642   2.306344 4.021438e-02
#> B-C.892          Smco4      B-C   -9.337324 -3.223009  -2.300364 4.064883e-02
#> B-C.893       Slc39a11      B-C    2.519127  1.332924   2.300358 4.064931e-02
#> B-C.894           Rdm1      B-C   -6.178195 -2.627185  -2.298125 4.081269e-02
#> B-C.895          Ift20      B-C   -2.447463 -1.291287  -2.297984 4.082297e-02
#> B-C.896           Bst2      B-C    2.363826  1.241124   2.297576 4.085290e-02
#> B-C.897           Dgkd      B-C    2.152626  1.106098   2.296536 4.092934e-02
#> B-C.898         Osbpl2      B-C    2.870372  1.521238   2.295462 4.100838e-02
#> B-C.899         Vps37b      B-C   -3.611950 -1.852778  -2.294356 4.108984e-02
#> B-C.900           Tmx4      B-C    2.035926  1.025685   2.289752 4.143094e-02
#> B-C.901         Thap11      B-C    3.787532  1.921258   2.289186 4.147303e-02
#> B-C.902         Abhd11      B-C    4.195024  2.068679   2.287796 4.157670e-02
#> B-C.903         H2-Eb1      B-C   -4.046151 -2.016550  -2.286308 4.168787e-02
#> B-C.904          Senp7      B-C   -4.178470 -2.062975  -2.282604 4.196583e-02
#> B-C.905          Lcorl      B-C   -8.594556 -3.103423  -2.281720 4.203244e-02
#> B-C.906           Selo      B-C   -3.559381 -1.831626  -2.278951 4.224175e-02
#> B-C.907          Adnp2      B-C   -3.415313 -1.772018  -2.278680 4.226233e-02
#> B-C.908          Prkd2      B-C   -5.739426 -2.520906  -2.275992 4.246650e-02
#> B-C.909         Uap1l1      B-C    2.887552  1.529847   2.272031 4.276922e-02
#> B-C.910           Twf1      B-C    3.653775  1.869388   2.269239 4.298376e-02
#> B-C.911           Pak1      B-C   11.103646  3.472962   2.268372 4.305060e-02
#> B-C.912          Frat2      B-C    2.626557  1.393173   2.267567 4.311276e-02
#> B-C.913          Xylt2      B-C   -3.038781 -1.603493  -2.267344 4.313000e-02
#> B-C.914         Gtf2e1      B-C    2.276526  1.186834   2.266356 4.320640e-02
#> B-C.915         Tsen54      B-C    2.833165  1.502415   2.264850 4.332313e-02
#> B-C.916         Zfp788      B-C   -7.827098 -2.968478  -2.262596 4.349846e-02
#> B-C.917        Supv3l1      B-C    3.745872  1.905302   2.260362 4.367285e-02
#> B-C.918           Phyh      B-C   -2.197944 -1.136155  -2.259092 4.377227e-02
#> B-C.919         Exoc6b      B-C   -2.106245 -1.074673  -2.259059 4.377486e-02
#> B-C.920         Man1a2      B-C   -2.061748 -1.043868  -2.253136 4.424150e-02
#> B-C.921         Hs6st1      B-C   -4.176103 -2.062157  -2.251013 4.440992e-02
#> B-C.922        Akr1c13      B-C   -9.109877 -3.187432  -2.250620 4.444115e-02
#> B-C.923  2310022A10Rik      B-C   -2.209999 -1.144045  -2.249460 4.453347e-02
#> B-C.924  2510003E04Rik      B-C   -2.348851 -1.231955  -2.248515 4.460888e-02
#> B-C.925           Shq1      B-C    2.540750  1.345254   2.247818 4.466456e-02
#> B-C.926          Usp34      B-C   -3.671332 -1.876304  -2.247411 4.469710e-02
#> B-C.927          Naa60      B-C   -2.130946 -1.091494  -2.245955 4.481365e-02
#> B-C.928          Gria3      B-C   11.417696  3.513200   2.245577 4.484392e-02
#> B-C.929           Mkl1      B-C    2.376064  1.248574   2.241840 4.514467e-02
#> B-C.930          Acap2      B-C   -3.154804 -1.657551  -2.240665 4.523957e-02
#> B-C.931           Msh3      B-C    4.092422  2.032955   2.240356 4.526458e-02
#> B-C.932           Sgcb      B-C    4.622437  2.208654   2.239679 4.531946e-02
#> B-C.933          Klhl2      B-C   -4.270335 -2.094349  -2.238890 4.538336e-02
#> B-C.934          Ltbp3      B-C   -5.167079 -2.369349  -2.238650 4.540285e-02
#> B-C.935          Akip1      B-C   -2.955510 -1.563407  -2.235304 4.567519e-02
#> B-C.936       Selenbp1      B-C    4.155274  2.054944   2.233425 4.582888e-02
#> B-C.937         Nudt18      B-C    4.646329  2.216091   2.230637 4.605769e-02
#> B-C.938           Selt      B-C   -2.994248 -1.582194  -2.229520 4.614972e-02
#> B-C.939           Lyz2      B-C   -3.977001 -1.991681  -2.228023 4.627325e-02
#> B-C.940          Tmub1      B-C   -3.611179 -1.852470  -2.227181 4.634290e-02
#> B-C.941           H6pd      B-C   -2.013905 -1.009996  -2.219585 4.697562e-02
#> B-C.942          Synj2      B-C   -5.585019 -2.481562  -2.210050 4.778151e-02
#> B-C.943          Kpna3      B-C    3.420525  1.774218   2.209985 4.778701e-02
#> B-C.944            Atm      B-C   -4.666163 -2.222237  -2.209408 4.783620e-02
#> B-C.945           Idua      B-C    2.780311  1.475246   2.208797 4.788835e-02
#> B-C.946         Zfp219      B-C    2.592818  1.374521   2.206091 4.811995e-02
#> B-C.947          Kntc1      B-C    3.580769  1.840270   2.204356 4.826897e-02
#> B-C.948           Gbp3      B-C    4.529871  2.179470   2.203021 4.838395e-02
#> B-C.949        Zfyve26      B-C    5.020509  2.327834   2.202405 4.843708e-02
#> B-C.950         Angpt1      B-C    2.831445  1.501539   2.200953 4.856258e-02
#> B-C.951          Urgcp      B-C   -2.112064 -1.078654  -2.200850 4.857153e-02
#> B-C.952           Prkx      B-C   -5.299816 -2.405942  -2.200537 4.859859e-02
#> B-C.953        Msantd4      B-C    6.672624  2.738254   2.197316 4.887821e-02
#> B-C.954          Stat4      B-C   -6.647083 -2.732721  -2.197263 4.888282e-02
#> B-C.955         Gm4759      B-C   -2.516545 -1.331445  -2.196153 4.897956e-02
#> B-C.956        Aldh7a1      B-C   -3.759412 -1.910507  -2.194654 4.911055e-02
#> B-C.957         Tecpr1      B-C  -17.692562 -4.145071  -2.194584 4.911661e-02
#> B-C.958           Lap3      B-C    2.577123  1.365761   2.194482 4.912558e-02
#> B-C.959         Zfp952      B-C  -10.258947 -3.358811  -2.192873 4.926651e-02
#> B-C.960           Tsr3      B-C   -3.225202 -1.689389  -2.190942 4.943612e-02
#> B-C.961           Nbas      B-C    4.434447  2.148754   2.190797 4.944894e-02
#> B-C.962        Slc50a1      B-C   -2.368075 -1.243715  -2.189213 4.958850e-02
#>               adjpval
#> B-A.1    2.509467e-05
#> B-A.2    1.096354e-04
#> B-A.3    2.894973e-04
#> B-A.4    3.040401e-04
#> B-A.5    3.232894e-04
#> B-A.6    3.232894e-04
#> B-A.7    6.272189e-04
#> B-A.8    9.557830e-04
#> B-A.9    9.557830e-04
#> B-A.10   1.184587e-03
#> B-A.11   1.960952e-03
#> B-A.12   2.359380e-03
#> B-A.13   3.427435e-03
#> B-A.14   3.427435e-03
#> B-A.15   3.427435e-03
#> B-A.16   3.427435e-03
#> B-A.17   3.589338e-03
#> B-A.18   4.101032e-03
#> B-A.19   4.101032e-03
#> B-A.20   4.101032e-03
#> B-A.21   4.511768e-03
#> B-A.22   5.953528e-03
#> B-A.23   7.056101e-03
#> B-A.24   7.529738e-03
#> B-A.25   7.858119e-03
#> B-A.26   8.045012e-03
#> B-A.27   8.355734e-03
#> B-A.28   9.465061e-03
#> B-A.29   9.465061e-03
#> B-A.30   9.465061e-03
#> B-A.31   1.055430e-02
#> B-A.32   1.055430e-02
#> B-A.33   1.256456e-02
#> B-A.34   1.367065e-02
#> B-A.35   1.367065e-02
#> B-A.36   1.647570e-02
#> B-A.37   1.647570e-02
#> B-A.38   1.650174e-02
#> B-A.39   1.650174e-02
#> B-A.40   1.650174e-02
#> B-A.41   1.650174e-02
#> B-A.42   1.925930e-02
#> B-A.43   2.125234e-02
#> B-A.44   2.527489e-02
#> B-A.45   2.533508e-02
#> B-A.46   2.650962e-02
#> B-A.47   2.674545e-02
#> B-A.48   2.857544e-02
#> B-A.49   3.073929e-02
#> B-A.50   3.073929e-02
#> B-A.51   3.073929e-02
#> B-A.52   3.089323e-02
#> B-A.53   3.115306e-02
#> B-A.54   3.209999e-02
#> B-A.55   3.209999e-02
#> B-A.56   3.712646e-02
#> B-A.57   3.738807e-02
#> B-A.58   3.785099e-02
#> B-A.59   4.194622e-02
#> B-A.60   4.430649e-02
#> B-A.61   4.691216e-02
#> B-A.62   4.767408e-02
#> B-A.63   4.767408e-02
#> B-A.64   4.767408e-02
#> B-A.65   4.767408e-02
#> B-A.66   4.767408e-02
#> B-A.67   4.826498e-02
#> B-A.68   6.293500e-02
#> B-A.69   6.300934e-02
#> B-A.70   6.306932e-02
#> B-A.71   6.306932e-02
#> B-A.72   6.306932e-02
#> B-A.73   6.306932e-02
#> B-A.74   6.306932e-02
#> B-A.75   6.306932e-02
#> B-A.76   6.313834e-02
#> B-A.77   6.313834e-02
#> B-A.78   6.505258e-02
#> B-A.79   6.505258e-02
#> B-A.80   6.505258e-02
#> B-A.81   6.568366e-02
#> B-A.82   6.839293e-02
#> B-A.83   6.839293e-02
#> B-A.84   6.945283e-02
#> B-A.85   6.945283e-02
#> B-A.86   7.031347e-02
#> B-A.87   7.104384e-02
#> B-A.88   7.104384e-02
#> B-A.89   7.104384e-02
#> B-A.90   7.104384e-02
#> B-A.91   7.223211e-02
#> B-A.92   7.666712e-02
#> B-A.93   8.145382e-02
#> B-A.94   8.206530e-02
#> B-A.95   8.227777e-02
#> B-A.96   8.227777e-02
#> B-A.97   8.494191e-02
#> B-A.98   8.648331e-02
#> B-A.99   8.692134e-02
#> B-A.100  8.828377e-02
#> B-A.101  8.828377e-02
#> B-A.102  8.832066e-02
#> B-A.103  9.155275e-02
#> B-A.104  9.340631e-02
#> B-A.105  9.340631e-02
#> B-A.106  9.375351e-02
#> B-A.107  9.375351e-02
#> B-A.108  9.375351e-02
#> B-A.109  9.825303e-02
#> B-A.110  9.825303e-02
#> B-A.111  9.825303e-02
#> B-A.112  9.842597e-02
#> B-A.113  9.842597e-02
#> B-A.114  9.842597e-02
#> B-A.115  1.017161e-01
#> B-A.116  1.021385e-01
#> B-A.117  1.021385e-01
#> B-A.118  1.025346e-01
#> B-A.119  1.027956e-01
#> B-A.120  1.060138e-01
#> B-A.121  1.064092e-01
#> B-A.122  1.064092e-01
#> B-A.123  1.071722e-01
#> B-A.124  1.083124e-01
#> B-A.125  1.083124e-01
#> B-A.126  1.083124e-01
#> B-A.127  1.083124e-01
#> B-A.128  1.083124e-01
#> B-A.129  1.099492e-01
#> B-A.130  1.103313e-01
#> B-A.131  1.108567e-01
#> B-A.132  1.108567e-01
#> B-A.133  1.108567e-01
#> B-A.134  1.134822e-01
#> B-A.135  1.143531e-01
#> B-A.136  1.187621e-01
#> B-A.137  1.211446e-01
#> B-A.138  1.211446e-01
#> B-A.139  1.211446e-01
#> B-A.140  1.211446e-01
#> B-A.141  1.211446e-01
#> B-A.142  1.211446e-01
#> B-A.143  1.247963e-01
#> B-A.144  1.247963e-01
#> B-A.145  1.247963e-01
#> B-A.146  1.265964e-01
#> B-A.147  1.273871e-01
#> B-A.148  1.324227e-01
#> B-A.149  1.336615e-01
#> B-A.150  1.338010e-01
#> B-A.151  1.338010e-01
#> B-A.152  1.360341e-01
#> B-A.153  1.362746e-01
#> B-A.154  1.370948e-01
#> B-A.155  1.387349e-01
#> B-A.156  1.387349e-01
#> B-A.157  1.387349e-01
#> B-A.158  1.387349e-01
#> B-A.159  1.394606e-01
#> B-A.160  1.407266e-01
#> B-A.161  1.407266e-01
#> B-A.162  1.407266e-01
#> B-A.163  1.419877e-01
#> B-A.164  1.419877e-01
#> B-A.165  1.419877e-01
#> B-A.166  1.419877e-01
#> B-A.167  1.419877e-01
#> B-A.168  1.440416e-01
#> B-A.169  1.440416e-01
#> B-A.170  1.450489e-01
#> B-A.171  1.468849e-01
#> B-A.172  1.473105e-01
#> B-A.173  1.473297e-01
#> B-A.174  1.481587e-01
#> B-A.175  1.484332e-01
#> B-A.176  1.484332e-01
#> B-A.177  1.508783e-01
#> B-A.178  1.508783e-01
#> B-A.179  1.508783e-01
#> B-A.180  1.508783e-01
#> B-A.181  1.508783e-01
#> B-A.182  1.508783e-01
#> B-A.183  1.517359e-01
#> B-A.184  1.518985e-01
#> B-A.185  1.557526e-01
#> B-A.186  1.557526e-01
#> B-A.187  1.557526e-01
#> B-A.188  1.585569e-01
#> B-A.189  1.585569e-01
#> B-A.190  1.598197e-01
#> B-A.191  1.609137e-01
#> B-A.192  1.609137e-01
#> B-A.193  1.609137e-01
#> B-A.194  1.625534e-01
#> B-A.195  1.716982e-01
#> B-A.196  1.716982e-01
#> B-A.197  1.716982e-01
#> B-A.198  1.716982e-01
#> B-A.199  1.716982e-01
#> B-A.200  1.717815e-01
#> B-A.201  1.727633e-01
#> B-A.202  1.727633e-01
#> B-A.203  1.727633e-01
#> B-A.204  1.727633e-01
#> B-A.205  1.727633e-01
#> B-A.206  1.727633e-01
#> B-A.207  1.727633e-01
#> B-A.208  1.727633e-01
#> B-A.209  1.727633e-01
#> B-A.210  1.729102e-01
#> B-A.211  1.750593e-01
#> B-A.212  1.750593e-01
#> B-A.213  1.750593e-01
#> B-A.214  1.750593e-01
#> B-A.215  1.843651e-01
#> B-A.216  1.843651e-01
#> B-A.217  1.847131e-01
#> B-A.218  1.847131e-01
#> B-A.219  1.847131e-01
#> B-A.220  1.866423e-01
#> B-A.221  1.866423e-01
#> B-A.222  1.866423e-01
#> B-A.223  1.866423e-01
#> B-A.224  1.881756e-01
#> B-A.225  1.881756e-01
#> B-A.226  1.886735e-01
#> B-A.227  1.888425e-01
#> B-A.228  1.904785e-01
#> B-A.229  1.904785e-01
#> B-A.230  1.910933e-01
#> B-A.231  1.915260e-01
#> B-A.232  1.951582e-01
#> B-A.233  1.967938e-01
#> B-A.234  1.978077e-01
#> B-A.235  1.978077e-01
#> B-A.236  1.982172e-01
#> B-A.237  1.998484e-01
#> B-A.238  2.011802e-01
#> B-A.239  2.033918e-01
#> B-A.240  2.033918e-01
#> B-A.241  2.044134e-01
#> B-A.242  2.080894e-01
#> B-A.243  2.080894e-01
#> B-A.244  2.080894e-01
#> B-A.245  2.094292e-01
#> B-A.246  2.094292e-01
#> B-A.247  2.094292e-01
#> B-A.248  2.103721e-01
#> B-A.249  2.103721e-01
#> B-A.250  2.116339e-01
#> B-A.251  2.130456e-01
#> B-A.252  2.133670e-01
#> B-A.253  2.133670e-01
#> B-A.254  2.133670e-01
#> B-A.255  2.136021e-01
#> B-A.256  2.151750e-01
#> B-A.257  2.151750e-01
#> B-A.258  2.151750e-01
#> B-A.259  2.166044e-01
#> B-A.260  2.166044e-01
#> B-A.261  2.190061e-01
#> B-A.262  2.238778e-01
#> B-A.263  2.259522e-01
#> B-A.264  2.261134e-01
#> B-A.265  2.281033e-01
#> B-A.266  2.294201e-01
#> B-A.267  2.297850e-01
#> B-A.268  2.297850e-01
#> B-A.269  2.297850e-01
#> B-A.270  2.314053e-01
#> B-A.271  2.314053e-01
#> B-A.272  2.314053e-01
#> B-A.273  2.314053e-01
#> B-A.274  2.314053e-01
#> B-A.275  2.314053e-01
#> B-A.276  2.314053e-01
#> B-A.277  2.314053e-01
#> B-A.278  2.314053e-01
#> B-A.279  2.314053e-01
#> B-A.280  2.314053e-01
#> B-A.281  2.333309e-01
#> B-A.282  2.346004e-01
#> B-A.283  2.364200e-01
#> B-A.284  2.364200e-01
#> B-A.285  2.364200e-01
#> B-A.286  2.364200e-01
#> B-A.287  2.364200e-01
#> B-A.288  2.390670e-01
#> B-A.289  2.390670e-01
#> B-A.290  2.390670e-01
#> B-A.291  2.390670e-01
#> B-A.292  2.410836e-01
#> B-A.293  2.410836e-01
#> B-A.294  2.432183e-01
#> B-A.295  2.433953e-01
#> B-A.296  2.437908e-01
#> B-A.297  2.439351e-01
#> B-A.298  2.486941e-01
#> B-A.299  2.486941e-01
#> B-A.300  2.489212e-01
#> B-A.301  2.496917e-01
#> B-A.302  2.496917e-01
#> B-A.303  2.506273e-01
#> B-A.304  2.506273e-01
#> B-A.305  2.517608e-01
#> B-A.306  2.524787e-01
#> B-A.307  2.528427e-01
#> B-A.308  2.581862e-01
#> B-A.309  2.581862e-01
#> B-A.310  2.581862e-01
#> B-A.311  2.581862e-01
#> B-A.312  2.581862e-01
#> B-A.313  2.581862e-01
#> B-A.314  2.581862e-01
#> B-A.315  2.581862e-01
#> B-A.316  2.581862e-01
#> B-A.317  2.587564e-01
#> B-A.318  2.604686e-01
#> B-A.319  2.621808e-01
#> B-A.320  2.621808e-01
#> B-A.321  2.652134e-01
#> B-A.322  2.656583e-01
#> B-A.323  2.666350e-01
#> B-A.324  2.668849e-01
#> B-A.325  2.668849e-01
#> B-A.326  2.668849e-01
#> B-A.327  2.681472e-01
#> B-A.328  2.682243e-01
#> B-A.329  2.705014e-01
#> B-A.330  2.725465e-01
#> B-A.331  2.725465e-01
#> B-A.332  2.725465e-01
#> B-A.333  2.725465e-01
#> B-A.334  2.756918e-01
#> B-A.335  2.756918e-01
#> B-A.336  2.756918e-01
#> B-A.337  2.756918e-01
#> B-A.338  2.758340e-01
#> B-A.339  2.762163e-01
#> B-A.340  2.762163e-01
#> B-A.341  2.808536e-01
#> B-A.342  2.835070e-01
#> B-A.343  2.865282e-01
#> B-A.344  2.866064e-01
#> B-A.345  2.888206e-01
#> B-A.346  2.909781e-01
#> B-A.347  2.914846e-01
#> B-A.348  2.927597e-01
#> B-A.349  2.950326e-01
#> B-A.350  2.953156e-01
#> B-A.351  2.957245e-01
#> B-A.352  2.958138e-01
#> B-A.353  2.958138e-01
#> B-A.354  2.958138e-01
#> B-A.355  2.958138e-01
#> B-A.356  2.975615e-01
#> B-A.357  2.981010e-01
#> B-A.358  2.981010e-01
#> B-A.359  2.999486e-01
#> B-A.360  3.007124e-01
#> B-A.361  3.034817e-01
#> B-A.362  3.034817e-01
#> B-A.363  3.034817e-01
#> B-A.364  3.034817e-01
#> B-A.365  3.034817e-01
#> B-A.366  3.038593e-01
#> B-A.367  3.050556e-01
#> B-A.368  3.070344e-01
#> B-A.369  3.076495e-01
#> B-A.370  3.090215e-01
#> B-A.371  3.090215e-01
#> B-A.372  3.100313e-01
#> B-A.373  3.103354e-01
#> B-A.374  3.103354e-01
#> B-A.375  3.114728e-01
#> B-A.376  3.117285e-01
#> B-A.377  3.117285e-01
#> B-A.378  3.117285e-01
#> B-A.379  3.117285e-01
#> B-A.380  3.117285e-01
#> B-A.381  3.135317e-01
#> B-A.382  3.135317e-01
#> B-A.383  3.173913e-01
#> B-A.384  3.189445e-01
#> B-A.385  3.190365e-01
#> B-A.386  3.201145e-01
#> B-A.387  3.201145e-01
#> B-A.388  3.201145e-01
#> B-A.389  3.227399e-01
#> B-A.390  3.235756e-01
#> B-A.391  3.252229e-01
#> B-A.392  3.268006e-01
#> B-A.393  3.270473e-01
#> B-A.394  3.270473e-01
#> B-A.395  3.270473e-01
#> B-A.396  3.270473e-01
#> B-A.397  3.272167e-01
#> B-A.398  3.279091e-01
#> B-A.399  3.292835e-01
#> B-A.400  3.292835e-01
#> B-A.401  3.299339e-01
#> B-A.402  3.299339e-01
#> B-A.403  3.299339e-01
#> B-A.404  3.299339e-01
#> B-A.405  3.321444e-01
#> B-A.406  3.322928e-01
#> B-A.407  3.322928e-01
#> B-A.408  3.322928e-01
#> B-A.409  3.343926e-01
#> B-A.410  3.343926e-01
#> B-A.411  3.343926e-01
#> B-A.412  3.343926e-01
#> B-A.413  3.343926e-01
#> B-A.414  3.343926e-01
#> B-A.415  3.343926e-01
#> B-A.416  3.343926e-01
#> B-A.417  3.343926e-01
#> B-A.418  3.349667e-01
#> B-A.419  3.354262e-01
#> B-A.420  3.354262e-01
#> B-A.421  3.358365e-01
#> B-A.422  3.371430e-01
#> B-A.423  3.371430e-01
#> B-A.424  3.371430e-01
#> B-A.425  3.371430e-01
#> B-A.426  3.371430e-01
#> B-A.427  3.374487e-01
#> B-A.428  3.440167e-01
#> B-A.429  3.440167e-01
#> B-A.430  3.440167e-01
#> B-A.431  3.440167e-01
#> B-A.432  3.440167e-01
#> B-A.433  3.440167e-01
#> B-A.434  3.440167e-01
#> B-A.435  3.445724e-01
#> B-A.436  3.453656e-01
#> B-A.437  3.462322e-01
#> B-A.438  3.462322e-01
#> B-A.439  3.462322e-01
#> B-A.440  3.464674e-01
#> B-A.441  3.464674e-01
#> B-A.442  3.471781e-01
#> B-A.443  3.478400e-01
#> B-A.444  3.484105e-01
#> B-A.445  3.494257e-01
#> B-A.446  3.494257e-01
#> B-A.447  3.507546e-01
#> B-A.448  3.520026e-01
#> B-A.449  3.520026e-01
#> B-A.450  3.536486e-01
#> B-A.451  3.536486e-01
#> B-A.452  3.564904e-01
#> B-A.453  3.594076e-01
#> B-A.454  3.608412e-01
#> B-A.455  3.612695e-01
#> B-A.456  3.664231e-01
#> B-A.457  3.664231e-01
#> B-A.458  3.665333e-01
#> B-A.459  3.665333e-01
#> B-A.460  3.666846e-01
#> B-A.461  3.666846e-01
#> B-A.462  3.710535e-01
#> B-A.463  3.718161e-01
#> B-A.464  3.718161e-01
#> B-A.465  3.728706e-01
#> B-A.466  3.733668e-01
#> B-A.467  3.733668e-01
#> B-A.468  3.735636e-01
#> B-A.469  3.735636e-01
#> B-A.470  3.744408e-01
#> B-A.471  3.744408e-01
#> B-A.472  3.744408e-01
#> B-A.473  3.744408e-01
#> B-A.474  3.745958e-01
#> B-A.475  3.768817e-01
#> B-A.476  3.770076e-01
#> B-A.477  3.772288e-01
#> B-A.478  3.772288e-01
#> B-A.479  3.772288e-01
#> B-A.480  3.773389e-01
#> B-A.481  3.773389e-01
#> B-A.482  3.786923e-01
#> B-A.483  3.790512e-01
#> B-A.484  3.790512e-01
#> B-A.485  3.809561e-01
#> B-A.486  3.853751e-01
#> B-A.487  3.853751e-01
#> B-A.488  3.864493e-01
#> B-A.489  3.864493e-01
#> B-A.490  3.867227e-01
#> B-A.491  3.897293e-01
#> B-A.492  3.931201e-01
#> B-A.493  3.931201e-01
#> B-A.494  3.931201e-01
#> B-A.495  3.931201e-01
#> B-A.496  3.931201e-01
#> B-A.497  3.933274e-01
#> B-A.498  3.935770e-01
#> B-A.499  3.935770e-01
#> B-A.500  3.935770e-01
#> B-A.501  3.954103e-01
#> B-A.502  3.961877e-01
#> B-A.503  3.968263e-01
#> B-A.504  3.987463e-01
#> B-A.505  4.019541e-01
#> B-A.506  4.022873e-01
#> B-A.507  4.022873e-01
#> B-A.508  4.022873e-01
#> B-A.509  4.045898e-01
#> B-A.510  4.045898e-01
#> B-A.511  4.050968e-01
#> B-A.512  4.050968e-01
#> B-A.513  4.050968e-01
#> B-A.514  4.050968e-01
#> B-A.515  4.050968e-01
#> B-A.516  4.050968e-01
#> B-A.517  4.050968e-01
#> B-A.518  4.050968e-01
#> B-A.519  4.050968e-01
#> B-A.520  4.050968e-01
#> B-A.521  4.050968e-01
#> B-A.522  4.050968e-01
#> B-A.523  4.071532e-01
#> B-A.524  4.074347e-01
#> B-A.525  4.095430e-01
#> B-A.526  4.118374e-01
#> B-A.527  4.138602e-01
#> B-A.528  4.140937e-01
#> B-A.529  4.140937e-01
#> B-A.530  4.146766e-01
#> B-A.531  4.146766e-01
#> B-A.532  4.146766e-01
#> B-A.533  4.151399e-01
#> B-A.534  4.151399e-01
#> B-A.535  4.151399e-01
#> B-A.536  4.157985e-01
#> B-A.537  4.173150e-01
#> B-A.538  4.173150e-01
#> B-A.539  4.173150e-01
#> B-A.540  4.179057e-01
#> B-A.541  4.183056e-01
#> B-A.542  4.188114e-01
#> B-A.543  4.188114e-01
#> B-A.544  4.188114e-01
#> B-A.545  4.193647e-01
#> B-A.546  4.193898e-01
#> B-A.547  4.224819e-01
#> B-A.548  4.224819e-01
#> B-A.549  4.224819e-01
#> B-A.550  4.259498e-01
#> B-A.551  4.259498e-01
#> B-A.552  4.259498e-01
#> B-A.553  4.259498e-01
#> B-A.554  4.259498e-01
#> B-A.555  4.259498e-01
#> B-A.556  4.259498e-01
#> B-A.557  4.259498e-01
#> B-A.558  4.259498e-01
#> B-A.559  4.279879e-01
#> B-A.560  4.304997e-01
#> B-A.561  4.312304e-01
#> B-A.562  4.312304e-01
#> B-A.563  4.321032e-01
#> B-A.564  4.345011e-01
#> B-A.565  4.345011e-01
#> B-A.566  4.345011e-01
#> B-A.567  4.345011e-01
#> B-A.568  4.345011e-01
#> B-A.569  4.345011e-01
#> B-A.570  4.345011e-01
#> B-A.571  4.345011e-01
#> B-A.572  4.345011e-01
#> B-A.573  4.345011e-01
#> B-A.574  4.345011e-01
#> B-A.575  4.345011e-01
#> B-A.576  4.345011e-01
#> B-A.577  4.346855e-01
#> B-A.578  4.346855e-01
#> B-A.579  4.346855e-01
#> B-A.580  4.368225e-01
#> B-A.581  4.384295e-01
#> B-A.582  4.396172e-01
#> B-A.583  4.405233e-01
#> B-A.584  4.405496e-01
#> B-A.585  4.407212e-01
#> B-A.586  4.407212e-01
#> B-A.587  4.407653e-01
#> B-A.588  4.441033e-01
#> B-A.589  4.441033e-01
#> B-A.590  4.454603e-01
#> B-A.591  4.462953e-01
#> B-A.592  4.464718e-01
#> B-A.593  4.469287e-01
#> B-A.594  4.473969e-01
#> B-A.595  4.473969e-01
#> B-A.596  4.476050e-01
#> B-A.597  4.476050e-01
#> B-A.598  4.476050e-01
#> B-A.599  4.476050e-01
#> B-A.600  4.476050e-01
#> B-A.601  4.476050e-01
#> B-A.602  4.476050e-01
#> B-A.603  4.476050e-01
#> B-A.604  4.476050e-01
#> B-A.605  4.476050e-01
#> B-A.606  4.488276e-01
#> B-A.607  4.488276e-01
#> B-A.608  4.488276e-01
#> B-A.609  4.504149e-01
#> B-A.610  4.509440e-01
#> B-A.611  4.509440e-01
#> B-A.612  4.509440e-01
#> C-A.1    1.718259e-05
#> C-A.2    3.279939e-05
#> C-A.3    1.815498e-04
#> C-A.4    1.815498e-04
#> C-A.5    1.815498e-04
#> C-A.6    1.815498e-04
#> C-A.7    1.815498e-04
#> C-A.8    1.815498e-04
#> C-A.9    1.815498e-04
#> C-A.10   1.815498e-04
#> C-A.11   1.815498e-04
#> C-A.12   1.831896e-04
#> C-A.13   1.831896e-04
#> C-A.14   1.831896e-04
#> C-A.15   1.887658e-04
#> C-A.16   1.887658e-04
#> C-A.17   2.033990e-04
#> C-A.18   2.055845e-04
#> C-A.19   2.200977e-04
#> C-A.20   2.200977e-04
#> C-A.21   2.200977e-04
#> C-A.22   2.200977e-04
#> C-A.23   2.407833e-04
#> C-A.24   2.422509e-04
#> C-A.25   2.639199e-04
#> C-A.26   3.278056e-04
#> C-A.27   3.278056e-04
#> C-A.28   3.278056e-04
#> C-A.29   3.300273e-04
#> C-A.30   3.300273e-04
#> C-A.31   3.765305e-04
#> C-A.32   3.939147e-04
#> C-A.33   3.939147e-04
#> C-A.34   4.906170e-04
#> C-A.35   5.060033e-04
#> C-A.36   5.060033e-04
#> C-A.37   5.060033e-04
#> C-A.38   5.369857e-04
#> C-A.39   5.427159e-04
#> C-A.40   5.762145e-04
#> C-A.41   6.640546e-04
#> C-A.42   6.816507e-04
#> C-A.43   7.045020e-04
#> C-A.44   7.846914e-04
#> C-A.45   7.870199e-04
#> C-A.46   7.870199e-04
#> C-A.47   7.870199e-04
#> C-A.48   8.131678e-04
#> C-A.49   8.344280e-04
#> C-A.50   8.625042e-04
#> C-A.51   9.149130e-04
#> C-A.52   1.015646e-03
#> C-A.53   1.032445e-03
#> C-A.54   1.037394e-03
#> C-A.55   1.037394e-03
#> C-A.56   1.121097e-03
#> C-A.57   1.121097e-03
#> C-A.58   1.121097e-03
#> C-A.59   1.143406e-03
#> C-A.60   1.196626e-03
#> C-A.61   1.205943e-03
#> C-A.62   1.254471e-03
#> C-A.63   1.406694e-03
#> C-A.64   1.477762e-03
#> C-A.65   1.530270e-03
#> C-A.66   1.575513e-03
#> C-A.67   1.575513e-03
#> C-A.68   1.686300e-03
#> C-A.69   1.717674e-03
#> C-A.70   1.762318e-03
#> C-A.71   1.789702e-03
#> C-A.72   1.883478e-03
#> C-A.73   1.885235e-03
#> C-A.74   1.885235e-03
#> C-A.75   1.988306e-03
#> C-A.76   2.038561e-03
#> C-A.77   2.085430e-03
#> C-A.78   2.126150e-03
#> C-A.79   2.212272e-03
#> C-A.80   2.217770e-03
#> C-A.81   2.217770e-03
#> C-A.82   2.234850e-03
#> C-A.83   2.234850e-03
#> C-A.84   2.243080e-03
#> C-A.85   2.335970e-03
#> C-A.86   2.572550e-03
#> C-A.87   2.572550e-03
#> C-A.88   2.692151e-03
#> C-A.89   2.845881e-03
#> C-A.90   2.845881e-03
#> C-A.91   2.912242e-03
#> C-A.92   2.941129e-03
#> C-A.93   2.974422e-03
#> C-A.94   3.027335e-03
#> C-A.95   3.176517e-03
#> C-A.96   3.274867e-03
#> C-A.97   3.274867e-03
#> C-A.98   3.290704e-03
#> C-A.99   3.364131e-03
#> C-A.100  3.624524e-03
#> C-A.101  3.628027e-03
#> C-A.102  3.628027e-03
#> C-A.103  3.712518e-03
#> C-A.104  3.712518e-03
#> C-A.105  3.766136e-03
#> C-A.106  3.832854e-03
#> C-A.107  3.832854e-03
#> C-A.108  3.849367e-03
#> C-A.109  4.025755e-03
#> C-A.110  4.036032e-03
#> C-A.111  4.036032e-03
#> C-A.112  4.081274e-03
#> C-A.113  4.203205e-03
#> C-A.114  4.203205e-03
#> C-A.115  4.630894e-03
#> C-A.116  4.680054e-03
#> C-A.117  4.742680e-03
#> C-A.118  4.742680e-03
#> C-A.119  4.816168e-03
#> C-A.120  4.816168e-03
#> C-A.121  4.816168e-03
#> C-A.122  4.845678e-03
#> C-A.123  4.887516e-03
#> C-A.124  4.939061e-03
#> C-A.125  4.939061e-03
#> C-A.126  5.312497e-03
#> C-A.127  5.366938e-03
#> C-A.128  5.388614e-03
#> C-A.129  5.468321e-03
#> C-A.130  5.648372e-03
#> C-A.131  5.706438e-03
#> C-A.132  5.826870e-03
#> C-A.133  5.931743e-03
#> C-A.134  6.030486e-03
#> C-A.135  6.030486e-03
#> C-A.136  6.030486e-03
#> C-A.137  6.030486e-03
#> C-A.138  6.030486e-03
#> C-A.139  6.030486e-03
#> C-A.140  6.087636e-03
#> C-A.141  6.136018e-03
#> C-A.142  6.404411e-03
#> C-A.143  6.435965e-03
#> C-A.144  6.542071e-03
#> C-A.145  6.553380e-03
#> C-A.146  6.579603e-03
#> C-A.147  6.656114e-03
#> C-A.148  6.667405e-03
#> C-A.149  6.671419e-03
#> C-A.150  6.714501e-03
#> C-A.151  6.926400e-03
#> C-A.152  6.975433e-03
#> C-A.153  7.045186e-03
#> C-A.154  7.581792e-03
#> C-A.155  7.581792e-03
#> C-A.156  7.588601e-03
#> C-A.157  7.610737e-03
#> C-A.158  7.799036e-03
#> C-A.159  7.843401e-03
#> C-A.160  7.919330e-03
#> C-A.161  7.919330e-03
#> C-A.162  7.972817e-03
#> C-A.163  7.972817e-03
#> C-A.164  8.136623e-03
#> C-A.165  8.371165e-03
#> C-A.166  8.371165e-03
#> C-A.167  8.371165e-03
#> C-A.168  8.371165e-03
#> C-A.169  8.611691e-03
#> C-A.170  8.684923e-03
#> C-A.171  9.026274e-03
#> C-A.172  9.026274e-03
#> C-A.173  9.292014e-03
#> C-A.174  9.292014e-03
#> C-A.175  9.292014e-03
#> C-A.176  9.292014e-03
#> C-A.177  9.406443e-03
#> C-A.178  9.731288e-03
#> C-A.179  9.762948e-03
#> C-A.180  9.815448e-03
#> C-A.181  9.815448e-03
#> C-A.182  9.901433e-03
#> C-A.183  9.923733e-03
#> C-A.184  1.018336e-02
#> C-A.185  1.018336e-02
#> C-A.186  1.021033e-02
#> C-A.187  1.021033e-02
#> C-A.188  1.021033e-02
#> C-A.189  1.021033e-02
#> C-A.190  1.021033e-02
#> C-A.191  1.021033e-02
#> C-A.192  1.021033e-02
#> C-A.193  1.044737e-02
#> C-A.194  1.044737e-02
#> C-A.195  1.044737e-02
#> C-A.196  1.054484e-02
#> C-A.197  1.054484e-02
#> C-A.198  1.054484e-02
#> C-A.199  1.054484e-02
#> C-A.200  1.065144e-02
#> C-A.201  1.075073e-02
#> C-A.202  1.088902e-02
#> C-A.203  1.094139e-02
#> C-A.204  1.115781e-02
#> C-A.205  1.123375e-02
#> C-A.206  1.130056e-02
#> C-A.207  1.141754e-02
#> C-A.208  1.150689e-02
#> C-A.209  1.150689e-02
#> C-A.210  1.181101e-02
#> C-A.211  1.198447e-02
#> C-A.212  1.204785e-02
#> C-A.213  1.204785e-02
#> C-A.214  1.204785e-02
#> C-A.215  1.215223e-02
#> C-A.216  1.237698e-02
#> C-A.217  1.246689e-02
#> C-A.218  1.246689e-02
#> C-A.219  1.250966e-02
#> C-A.220  1.267522e-02
#> C-A.221  1.281068e-02
#> C-A.222  1.309799e-02
#> C-A.223  1.319734e-02
#> C-A.224  1.319734e-02
#> C-A.225  1.319734e-02
#> C-A.226  1.319734e-02
#> C-A.227  1.325318e-02
#> C-A.228  1.330943e-02
#> C-A.229  1.352698e-02
#> C-A.230  1.365600e-02
#> C-A.231  1.365600e-02
#> C-A.232  1.371302e-02
#> C-A.233  1.377035e-02
#> C-A.234  1.395588e-02
#> C-A.235  1.396980e-02
#> C-A.236  1.456312e-02
#> C-A.237  1.456312e-02
#> C-A.238  1.456312e-02
#> C-A.239  1.558280e-02
#> C-A.240  1.568646e-02
#> C-A.241  1.568646e-02
#> C-A.242  1.568646e-02
#> C-A.243  1.596037e-02
#> C-A.244  1.599803e-02
#> C-A.245  1.622834e-02
#> C-A.246  1.643026e-02
#> C-A.247  1.643515e-02
#> C-A.248  1.665789e-02
#> C-A.249  1.665789e-02
#> C-A.250  1.668737e-02
#> C-A.251  1.673258e-02
#> C-A.252  1.698313e-02
#> C-A.253  1.698313e-02
#> C-A.254  1.699160e-02
#> C-A.255  1.699160e-02
#> C-A.256  1.699160e-02
#> C-A.257  1.699160e-02
#> C-A.258  1.728639e-02
#> C-A.259  1.738348e-02
#> C-A.260  1.763405e-02
#> C-A.261  1.763405e-02
#> C-A.262  1.763405e-02
#> C-A.263  1.763405e-02
#> C-A.264  1.763405e-02
#> C-A.265  1.763405e-02
#> C-A.266  1.763405e-02
#> C-A.267  1.763405e-02
#> C-A.268  1.774355e-02
#> C-A.269  1.777173e-02
#> C-A.270  1.780869e-02
#> C-A.271  1.785625e-02
#> C-A.272  1.789381e-02
#> C-A.273  1.805003e-02
#> C-A.274  1.817486e-02
#> C-A.275  1.817486e-02
#> C-A.276  1.819571e-02
#> C-A.277  1.852870e-02
#> C-A.278  1.856303e-02
#> C-A.279  1.872199e-02
#> C-A.280  1.875492e-02
#> C-A.281  1.880023e-02
#> C-A.282  1.880023e-02
#> C-A.283  1.880023e-02
#> C-A.284  1.880023e-02
#> C-A.285  1.880298e-02
#> C-A.286  1.900942e-02
#> C-A.287  1.910330e-02
#> C-A.288  1.910330e-02
#> C-A.289  1.912976e-02
#> C-A.290  1.942520e-02
#> C-A.291  1.942520e-02
#> C-A.292  1.942520e-02
#> C-A.293  1.944144e-02
#> C-A.294  1.944144e-02
#> C-A.295  2.007131e-02
#> C-A.296  2.007131e-02
#> C-A.297  2.007131e-02
#> C-A.298  2.007131e-02
#> C-A.299  2.007131e-02
#> C-A.300  2.007131e-02
#> C-A.301  2.042434e-02
#> C-A.302  2.058071e-02
#> C-A.303  2.058071e-02
#> C-A.304  2.176323e-02
#> C-A.305  2.195070e-02
#> C-A.306  2.198956e-02
#> C-A.307  2.202648e-02
#> C-A.308  2.202648e-02
#> C-A.309  2.204245e-02
#> C-A.310  2.204245e-02
#> C-A.311  2.204245e-02
#> C-A.312  2.204245e-02
#> C-A.313  2.204245e-02
#> C-A.314  2.204245e-02
#> C-A.315  2.206896e-02
#> C-A.316  2.223519e-02
#> C-A.317  2.223519e-02
#> C-A.318  2.234597e-02
#> C-A.319  2.234597e-02
#> C-A.320  2.234597e-02
#> C-A.321  2.235546e-02
#> C-A.322  2.237479e-02
#> C-A.323  2.237479e-02
#> C-A.324  2.237479e-02
#> C-A.325  2.247159e-02
#> C-A.326  2.247159e-02
#> C-A.327  2.247159e-02
#> C-A.328  2.247159e-02
#> C-A.329  2.247159e-02
#> C-A.330  2.278418e-02
#> C-A.331  2.289272e-02
#> C-A.332  2.289312e-02
#> C-A.333  2.291222e-02
#> C-A.334  2.292618e-02
#> C-A.335  2.325885e-02
#> C-A.336  2.325885e-02
#> C-A.337  2.387401e-02
#> C-A.338  2.414192e-02
#> C-A.339  2.415991e-02
#> C-A.340  2.473446e-02
#> C-A.341  2.473446e-02
#> C-A.342  2.474806e-02
#> C-A.343  2.498152e-02
#> C-A.344  2.515132e-02
#> C-A.345  2.515132e-02
#> C-A.346  2.520358e-02
#> C-A.347  2.525191e-02
#> C-A.348  2.525191e-02
#> C-A.349  2.525191e-02
#> C-A.350  2.559608e-02
#> C-A.351  2.572171e-02
#> C-A.352  2.586347e-02
#> C-A.353  2.586347e-02
#> C-A.354  2.586347e-02
#> C-A.355  2.586347e-02
#> C-A.356  2.586347e-02
#> C-A.357  2.623867e-02
#> C-A.358  2.649426e-02
#> C-A.359  2.661470e-02
#> C-A.360  2.682874e-02
#> C-A.361  2.696765e-02
#> C-A.362  2.715305e-02
#> C-A.363  2.766426e-02
#> C-A.364  2.769793e-02
#> C-A.365  2.782796e-02
#> C-A.366  2.782796e-02
#> C-A.367  2.782796e-02
#> C-A.368  2.782796e-02
#> C-A.369  2.782796e-02
#> C-A.370  2.782796e-02
#> C-A.371  2.782915e-02
#> C-A.372  2.802364e-02
#> C-A.373  2.832913e-02
#> C-A.374  2.895569e-02
#> C-A.375  2.911513e-02
#> C-A.376  2.911513e-02
#> C-A.377  2.940747e-02
#> C-A.378  2.945880e-02
#> C-A.379  2.945880e-02
#> C-A.380  2.956648e-02
#> C-A.381  2.993490e-02
#> C-A.382  2.993490e-02
#> C-A.383  2.993490e-02
#> C-A.384  2.993490e-02
#> C-A.385  2.993490e-02
#> C-A.386  3.012902e-02
#> C-A.387  3.027520e-02
#> C-A.388  3.027520e-02
#> C-A.389  3.027520e-02
#> C-A.390  3.027520e-02
#> C-A.391  3.027520e-02
#> C-A.392  3.044216e-02
#> C-A.393  3.055661e-02
#> C-A.394  3.055661e-02
#> C-A.395  3.070612e-02
#> C-A.396  3.083369e-02
#> C-A.397  3.087630e-02
#> C-A.398  3.093263e-02
#> C-A.399  3.100375e-02
#> C-A.400  3.113315e-02
#> C-A.401  3.125625e-02
#> C-A.402  3.125625e-02
#> C-A.403  3.154058e-02
#> C-A.404  3.155921e-02
#> C-A.405  3.155921e-02
#> C-A.406  3.155921e-02
#> C-A.407  3.167094e-02
#> C-A.408  3.167094e-02
#> C-A.409  3.183193e-02
#> C-A.410  3.197407e-02
#> C-A.411  3.216343e-02
#> C-A.412  3.222992e-02
#> C-A.413  3.224030e-02
#> C-A.414  3.239924e-02
#> C-A.415  3.242475e-02
#> C-A.416  3.259853e-02
#> C-A.417  3.272703e-02
#> C-A.418  3.272703e-02
#> C-A.419  3.275524e-02
#> C-A.420  3.298624e-02
#> C-A.421  3.298745e-02
#> C-A.422  3.298745e-02
#> C-A.423  3.317208e-02
#> C-A.424  3.317208e-02
#> C-A.425  3.317208e-02
#> C-A.426  3.334763e-02
#> C-A.427  3.339723e-02
#> C-A.428  3.365502e-02
#> C-A.429  3.365502e-02
#> C-A.430  3.380532e-02
#> C-A.431  3.380532e-02
#> C-A.432  3.380532e-02
#> C-A.433  3.380532e-02
#> C-A.434  3.396258e-02
#> C-A.435  3.403409e-02
#> C-A.436  3.403409e-02
#> C-A.437  3.431430e-02
#> C-A.438  3.432849e-02
#> C-A.439  3.466410e-02
#> C-A.440  3.484526e-02
#> C-A.441  3.484526e-02
#> C-A.442  3.484526e-02
#> C-A.443  3.484526e-02
#> C-A.444  3.484526e-02
#> C-A.445  3.484526e-02
#> C-A.446  3.487111e-02
#> C-A.447  3.510163e-02
#> C-A.448  3.557723e-02
#> C-A.449  3.575112e-02
#> C-A.450  3.575112e-02
#> C-A.451  3.587505e-02
#> C-A.452  3.600811e-02
#> C-A.453  3.605271e-02
#> C-A.454  3.625884e-02
#> C-A.455  3.625884e-02
#> C-A.456  3.640702e-02
#> C-A.457  3.674621e-02
#> C-A.458  3.681640e-02
#> C-A.459  3.681640e-02
#> C-A.460  3.681640e-02
#> C-A.461  3.681640e-02
#> C-A.462  3.694411e-02
#> C-A.463  3.695275e-02
#> C-A.464  3.705336e-02
#> C-A.465  3.707626e-02
#> C-A.466  3.707626e-02
#> C-A.467  3.722540e-02
#> C-A.468  3.722540e-02
#> C-A.469  3.755415e-02
#> C-A.470  3.757158e-02
#> C-A.471  3.822004e-02
#> C-A.472  3.832389e-02
#> C-A.473  3.832389e-02
#> C-A.474  3.840665e-02
#> C-A.475  3.840665e-02
#> C-A.476  3.872206e-02
#> C-A.477  3.872206e-02
#> C-A.478  3.872206e-02
#> C-A.479  3.894477e-02
#> C-A.480  3.895479e-02
#> C-A.481  3.895479e-02
#> C-A.482  3.927366e-02
#> C-A.483  3.948754e-02
#> C-A.484  3.948754e-02
#> C-A.485  3.948754e-02
#> C-A.486  3.949749e-02
#> C-A.487  3.957011e-02
#> C-A.488  3.968164e-02
#> C-A.489  3.983269e-02
#> C-A.490  4.022971e-02
#> C-A.491  4.022971e-02
#> C-A.492  4.022971e-02
#> C-A.493  4.031498e-02
#> C-A.494  4.031498e-02
#> C-A.495  4.044373e-02
#> C-A.496  4.058606e-02
#> C-A.497  4.071746e-02
#> C-A.498  4.092791e-02
#> C-A.499  4.099145e-02
#> C-A.500  4.138195e-02
#> C-A.501  4.142008e-02
#> C-A.502  4.144874e-02
#> C-A.503  4.144874e-02
#> C-A.504  4.148040e-02
#> C-A.505  4.179861e-02
#> C-A.506  4.205768e-02
#> C-A.507  4.266793e-02
#> C-A.508  4.289230e-02
#> C-A.509  4.312601e-02
#> C-A.510  4.312601e-02
#> C-A.511  4.312601e-02
#> C-A.512  4.320036e-02
#> C-A.513  4.320036e-02
#> C-A.514  4.323954e-02
#> C-A.515  4.379108e-02
#> C-A.516  4.395428e-02
#> C-A.517  4.409457e-02
#> C-A.518  4.425080e-02
#> C-A.519  4.442541e-02
#> C-A.520  4.442541e-02
#> C-A.521  4.442541e-02
#> C-A.522  4.442541e-02
#> C-A.523  4.474450e-02
#> C-A.524  4.479636e-02
#> C-A.525  4.479636e-02
#> C-A.526  4.481034e-02
#> C-A.527  4.481034e-02
#> C-A.528  4.481034e-02
#> C-A.529  4.481034e-02
#> C-A.530  4.481034e-02
#> C-A.531  4.503727e-02
#> C-A.532  4.508743e-02
#> C-A.533  4.528308e-02
#> C-A.534  4.530844e-02
#> C-A.535  4.558227e-02
#> C-A.536  4.558227e-02
#> C-A.537  4.558227e-02
#> C-A.538  4.562029e-02
#> C-A.539  4.562489e-02
#> C-A.540  4.586866e-02
#> C-A.541  4.607021e-02
#> C-A.542  4.618127e-02
#> C-A.543  4.624323e-02
#> C-A.544  4.673266e-02
#> C-A.545  4.685128e-02
#> C-A.546  4.695141e-02
#> C-A.547  4.758611e-02
#> C-A.548  4.758611e-02
#> C-A.549  4.758611e-02
#> C-A.550  4.763841e-02
#> C-A.551  4.795300e-02
#> C-A.552  4.796765e-02
#> C-A.553  4.806745e-02
#> C-A.554  4.815017e-02
#> C-A.555  4.839770e-02
#> C-A.556  4.850835e-02
#> C-A.557  4.852563e-02
#> C-A.558  4.852563e-02
#> C-A.559  4.863847e-02
#> C-A.560  4.863847e-02
#> C-A.561  4.866464e-02
#> C-A.562  4.866464e-02
#> C-A.563  4.866464e-02
#> C-A.564  4.939748e-02
#> C-A.565  4.958371e-02
#> C-A.566  4.958371e-02
#> C-A.567  4.976883e-02
#> C-A.568  5.052378e-02
#> C-A.569  5.052378e-02
#> C-A.570  5.052378e-02
#> C-A.571  5.052378e-02
#> C-A.572  5.069442e-02
#> C-A.573  5.069442e-02
#> C-A.574  5.090981e-02
#> C-A.575  5.108020e-02
#> C-A.576  5.115226e-02
#> C-A.577  5.161552e-02
#> C-A.578  5.173227e-02
#> C-A.579  5.194666e-02
#> C-A.580  5.195016e-02
#> C-A.581  5.219140e-02
#> C-A.582  5.219140e-02
#> C-A.583  5.219140e-02
#> C-A.584  5.219140e-02
#> C-A.585  5.219140e-02
#> C-A.586  5.219140e-02
#> C-A.587  5.219140e-02
#> C-A.588  5.219140e-02
#> C-A.589  5.222299e-02
#> C-A.590  5.224836e-02
#> C-A.591  5.262174e-02
#> C-A.592  5.262979e-02
#> C-A.593  5.314254e-02
#> C-A.594  5.331823e-02
#> C-A.595  5.331823e-02
#> C-A.596  5.331823e-02
#> C-A.597  5.331823e-02
#> C-A.598  5.331823e-02
#> C-A.599  5.331823e-02
#> C-A.600  5.331823e-02
#> C-A.601  5.331823e-02
#> C-A.602  5.331823e-02
#> C-A.603  5.331823e-02
#> C-A.604  5.331823e-02
#> C-A.605  5.347084e-02
#> C-A.606  5.378948e-02
#> C-A.607  5.419968e-02
#> C-A.608  5.435202e-02
#> C-A.609  5.435202e-02
#> C-A.610  5.435202e-02
#> C-A.611  5.437635e-02
#> C-A.612  5.446959e-02
#> C-A.613  5.446959e-02
#> C-A.614  5.446959e-02
#> C-A.615  5.446959e-02
#> C-A.616  5.477121e-02
#> C-A.617  5.643927e-02
#> C-A.618  5.679934e-02
#> C-A.619  5.691883e-02
#> C-A.620  5.722072e-02
#> C-A.621  5.722644e-02
#> C-A.622  5.722644e-02
#> C-A.623  5.722644e-02
#> C-A.624  5.722644e-02
#> C-A.625  5.722644e-02
#> C-A.626  5.722644e-02
#> C-A.627  5.722685e-02
#> C-A.628  5.750343e-02
#> C-A.629  5.750343e-02
#> C-A.630  5.775708e-02
#> C-A.631  5.775708e-02
#> C-A.632  5.775708e-02
#> C-A.633  5.785590e-02
#> C-A.634  5.790084e-02
#> C-A.635  5.823655e-02
#> C-A.636  5.826800e-02
#> C-A.637  5.826800e-02
#> C-A.638  5.876163e-02
#> C-A.639  5.876323e-02
#> C-A.640  5.887024e-02
#> C-A.641  5.914080e-02
#> C-A.642  5.914080e-02
#> C-A.643  5.914080e-02
#> C-A.644  5.914080e-02
#> C-A.645  5.914080e-02
#> C-A.646  5.948008e-02
#> C-A.647  5.954554e-02
#> C-A.648  5.988431e-02
#> C-A.649  5.988431e-02
#> C-A.650  5.988431e-02
#> C-A.651  5.988431e-02
#> C-A.652  5.988431e-02
#> C-A.653  5.991065e-02
#> C-A.654  5.991065e-02
#> C-A.655  6.012393e-02
#> C-A.656  6.017054e-02
#> C-A.657  6.041593e-02
#> C-A.658  6.075031e-02
#> C-A.659  6.096428e-02
#> C-A.660  6.099568e-02
#> C-A.661  6.135020e-02
#> C-A.662  6.135020e-02
#> C-A.663  6.190002e-02
#> C-A.664  6.217055e-02
#> C-A.665  6.218557e-02
#> C-A.666  6.228634e-02
#> C-A.667  6.228634e-02
#> C-A.668  6.236916e-02
#> C-A.669  6.237816e-02
#> C-A.670  6.237816e-02
#> C-A.671  6.245319e-02
#> C-A.672  6.245319e-02
#> C-A.673  6.252973e-02
#> C-A.674  6.290146e-02
#> C-A.675  6.339818e-02
#> C-A.676  6.369087e-02
#> C-A.677  6.374175e-02
#> C-A.678  6.390126e-02
#> C-A.679  6.390126e-02
#> C-A.680  6.394771e-02
#> C-A.681  6.394771e-02
#> C-A.682  6.397724e-02
#> C-A.683  6.471888e-02
#> C-A.684  6.504016e-02
#> C-A.685  6.504016e-02
#> C-A.686  6.504016e-02
#> C-A.687  6.504016e-02
#> C-A.688  6.504016e-02
#> C-A.689  6.539733e-02
#> C-A.690  6.562166e-02
#> C-A.691  6.575523e-02
#> C-A.692  6.590980e-02
#> C-A.693  6.590980e-02
#> C-A.694  6.648637e-02
#> C-A.695  6.713652e-02
#> C-A.696  6.762874e-02
#> C-A.697  6.821960e-02
#> C-A.698  6.836010e-02
#> C-A.699  6.842277e-02
#> C-A.700  6.847330e-02
#> C-A.701  6.848940e-02
#> C-A.702  6.848940e-02
#> C-A.703  6.848940e-02
#> C-A.704  6.848940e-02
#> C-A.705  6.853006e-02
#> C-A.706  6.856610e-02
#> C-A.707  6.880572e-02
#> C-A.708  6.880572e-02
#> C-A.709  6.900457e-02
#> C-A.710  6.978143e-02
#> C-A.711  7.040818e-02
#> C-A.712  7.042024e-02
#> C-A.713  7.086249e-02
#> C-A.714  7.086249e-02
#> C-A.715  7.086249e-02
#> C-A.716  7.097339e-02
#> C-A.717  7.104837e-02
#> C-A.718  7.158429e-02
#> C-A.719  7.186818e-02
#> C-A.720  7.190425e-02
#> C-A.721  7.222422e-02
#> C-A.722  7.234029e-02
#> C-A.723  7.274584e-02
#> C-A.724  7.305407e-02
#> C-A.725  7.309126e-02
#> C-A.726  7.378467e-02
#> C-A.727  7.378467e-02
#> C-A.728  7.384272e-02
#> C-A.729  7.407493e-02
#> C-A.730  7.407493e-02
#> C-A.731  7.471140e-02
#> C-A.732  7.471140e-02
#> C-A.733  7.471140e-02
#> C-A.734  7.471140e-02
#> C-A.735  7.494962e-02
#> C-A.736  7.514359e-02
#> C-A.737  7.514359e-02
#> C-A.738  7.542525e-02
#> C-A.739  7.542525e-02
#> C-A.740  7.542906e-02
#> C-A.741  7.542906e-02
#> C-A.742  7.546169e-02
#> C-A.743  7.550960e-02
#> C-A.744  7.550960e-02
#> C-A.745  7.617825e-02
#> C-A.746  7.618985e-02
#> C-A.747  7.643228e-02
#> C-A.748  7.643228e-02
#> C-A.749  7.712554e-02
#> C-A.750  7.714386e-02
#> C-A.751  7.714386e-02
#> C-A.752  7.714386e-02
#> C-A.753  7.727644e-02
#> C-A.754  7.785148e-02
#> C-A.755  7.819994e-02
#> C-A.756  7.858674e-02
#> C-A.757  7.880618e-02
#> C-A.758  7.882119e-02
#> C-A.759  7.882119e-02
#> C-A.760  7.882119e-02
#> C-A.761  7.900490e-02
#> C-A.762  7.917239e-02
#> C-A.763  7.917239e-02
#> C-A.764  7.924958e-02
#> C-A.765  7.924958e-02
#> C-A.766  7.924958e-02
#> C-A.767  7.924958e-02
#> C-A.768  7.936187e-02
#> C-A.769  7.954939e-02
#> C-A.770  7.954939e-02
#> C-A.771  7.958189e-02
#> C-A.772  7.958189e-02
#> C-A.773  7.958189e-02
#> C-A.774  7.958189e-02
#> C-A.775  7.965598e-02
#> C-A.776  7.969943e-02
#> C-A.777  7.969943e-02
#> C-A.778  7.974869e-02
#> C-A.779  7.981388e-02
#> C-A.780  7.981388e-02
#> C-A.781  7.999253e-02
#> C-A.782  8.015127e-02
#> C-A.783  8.022420e-02
#> C-A.784  8.027493e-02
#> C-A.785  8.047133e-02
#> C-A.786  8.095894e-02
#> C-A.787  8.100418e-02
#> C-A.788  8.125756e-02
#> C-A.789  8.152774e-02
#> C-A.790  8.152774e-02
#> C-A.791  8.152774e-02
#> C-A.792  8.152774e-02
#> C-A.793  8.165091e-02
#> C-A.794  8.196719e-02
#> C-A.795  8.259252e-02
#> C-A.796  8.259252e-02
#> C-A.797  8.270835e-02
#> C-A.798  8.270835e-02
#> C-A.799  8.280837e-02
#> C-A.800  8.299514e-02
#> C-A.801  8.323733e-02
#> C-A.802  8.361486e-02
#> C-A.803  8.490035e-02
#> C-A.804  8.511567e-02
#> C-A.805  8.515649e-02
#> C-A.806  8.515649e-02
#> C-A.807  8.515649e-02
#> C-A.808  8.515649e-02
#> C-A.809  8.578116e-02
#> C-A.810  8.578116e-02
#> C-A.811  8.643209e-02
#> C-A.812  8.713119e-02
#> C-A.813  8.713119e-02
#> C-A.814  8.726550e-02
#> C-A.815  8.750774e-02
#> C-A.816  8.809571e-02
#> C-A.817  8.931314e-02
#> C-A.818  8.948089e-02
#> C-A.819  8.948089e-02
#> C-A.820  8.975791e-02
#> C-A.821  8.987868e-02
#> C-A.822  9.030897e-02
#> C-A.823  9.030897e-02
#> C-A.824  9.030897e-02
#> C-A.825  9.030897e-02
#> C-A.826  9.069777e-02
#> C-A.827  9.090120e-02
#> C-A.828  9.093867e-02
#> C-A.829  9.104655e-02
#> C-A.830  9.136842e-02
#> C-A.831  9.142663e-02
#> C-A.832  9.142663e-02
#> C-A.833  9.153863e-02
#> C-A.834  9.175591e-02
#> C-A.835  9.213674e-02
#> C-A.836  9.280510e-02
#> C-A.837  9.301909e-02
#> C-A.838  9.301909e-02
#> C-A.839  9.339496e-02
#> C-A.840  9.356566e-02
#> C-A.841  9.369123e-02
#> C-A.842  9.380432e-02
#> C-A.843  9.386372e-02
#> C-A.844  9.452624e-02
#> C-A.845  9.452624e-02
#> C-A.846  9.459567e-02
#> C-A.847  9.466209e-02
#> C-A.848  9.567940e-02
#> C-A.849  9.592007e-02
#> C-A.850  9.592007e-02
#> C-A.851  9.615758e-02
#> C-A.852  9.642784e-02
#> C-A.853  9.647020e-02
#> C-A.854  9.668683e-02
#> C-A.855  9.668683e-02
#> C-A.856  9.672846e-02
#> C-A.857  9.679268e-02
#> C-A.858  9.703320e-02
#> C-A.859  9.710681e-02
#> C-A.860  9.710681e-02
#> C-A.861  9.710681e-02
#> C-A.862  9.710681e-02
#> C-A.863  9.715714e-02
#> C-A.864  9.715714e-02
#> C-A.865  9.749038e-02
#> C-A.866  9.759218e-02
#> C-A.867  9.767140e-02
#> C-A.868  9.767140e-02
#> C-A.869  9.767140e-02
#> C-A.870  9.822732e-02
#> C-A.871  9.866780e-02
#> C-A.872  9.884476e-02
#> C-A.873  9.892741e-02
#> C-A.874  9.926208e-02
#> C-A.875  9.926208e-02
#> C-A.876  9.926208e-02
#> C-A.877  9.930274e-02
#> C-A.878  9.964127e-02
#> C-A.879  9.964127e-02
#> C-A.880  9.998596e-02
#> C-A.881  1.000328e-01
#> C-A.882  1.006438e-01
#> C-A.883  1.007002e-01
#> C-A.884  1.007151e-01
#> C-A.885  1.007506e-01
#> C-A.886  1.007506e-01
#> C-A.887  1.015320e-01
#> C-A.888  1.024579e-01
#> C-A.889  1.024579e-01
#> C-A.890  1.024742e-01
#> C-A.891  1.029919e-01
#> C-A.892  1.029919e-01
#> C-A.893  1.029919e-01
#> C-A.894  1.029919e-01
#> C-A.895  1.031678e-01
#> C-A.896  1.032600e-01
#> C-A.897  1.032600e-01
#> C-A.898  1.037408e-01
#> C-A.899  1.037408e-01
#> C-A.900  1.037408e-01
#> C-A.901  1.038297e-01
#> C-A.902  1.038297e-01
#> C-A.903  1.038297e-01
#> C-A.904  1.038297e-01
#> C-A.905  1.038297e-01
#> C-A.906  1.038297e-01
#> C-A.907  1.043840e-01
#> C-A.908  1.043840e-01
#> C-A.909  1.046181e-01
#> C-A.910  1.047247e-01
#> C-A.911  1.047247e-01
#> C-A.912  1.057804e-01
#> C-A.913  1.068513e-01
#> C-A.914  1.073712e-01
#> C-A.915  1.074418e-01
#> C-A.916  1.074526e-01
#> C-A.917  1.075983e-01
#> C-A.918  1.075983e-01
#> C-A.919  1.077684e-01
#> C-A.920  1.078245e-01
#> C-A.921  1.082949e-01
#> C-A.922  1.085689e-01
#> C-A.923  1.086820e-01
#> C-A.924  1.087965e-01
#> C-A.925  1.090088e-01
#> C-A.926  1.093509e-01
#> C-A.927  1.093509e-01
#> C-A.928  1.096224e-01
#> C-A.929  1.097353e-01
#> C-A.930  1.103209e-01
#> C-A.931  1.103209e-01
#> C-A.932  1.103209e-01
#> C-A.933  1.103209e-01
#> C-A.934  1.103221e-01
#> C-A.935  1.103947e-01
#> C-A.936  1.104089e-01
#> C-A.937  1.104089e-01
#> C-A.938  1.104089e-01
#> C-A.939  1.104089e-01
#> C-A.940  1.104089e-01
#> C-A.941  1.107277e-01
#> C-A.942  1.107277e-01
#> C-A.943  1.107277e-01
#> C-A.944  1.107277e-01
#> C-A.945  1.109779e-01
#> C-A.946  1.113241e-01
#> C-A.947  1.114928e-01
#> C-A.948  1.117244e-01
#> C-A.949  1.117980e-01
#> C-A.950  1.117980e-01
#> C-A.951  1.117980e-01
#> C-A.952  1.119784e-01
#> C-A.953  1.119816e-01
#> C-A.954  1.119816e-01
#> C-A.955  1.119816e-01
#> C-A.956  1.119816e-01
#> C-A.957  1.119816e-01
#> C-A.958  1.123905e-01
#> C-A.959  1.123905e-01
#> C-A.960  1.125690e-01
#> C-A.961  1.126344e-01
#> C-A.962  1.127839e-01
#> C-A.963  1.127847e-01
#> C-A.964  1.132170e-01
#> C-A.965  1.132170e-01
#> C-A.966  1.132170e-01
#> C-A.967  1.132170e-01
#> C-A.968  1.132170e-01
#> C-A.969  1.132170e-01
#> C-A.970  1.132170e-01
#> C-A.971  1.132651e-01
#> C-A.972  1.132924e-01
#> C-A.973  1.132924e-01
#> C-A.974  1.132924e-01
#> C-A.975  1.134832e-01
#> C-A.976  1.135805e-01
#> C-A.977  1.138746e-01
#> C-A.978  1.139500e-01
#> C-A.979  1.140111e-01
#> C-A.980  1.142106e-01
#> C-A.981  1.144061e-01
#> C-A.982  1.145197e-01
#> C-A.983  1.155664e-01
#> C-A.984  1.156960e-01
#> C-A.985  1.158742e-01
#> C-A.986  1.160282e-01
#> C-A.987  1.161387e-01
#> C-A.988  1.162930e-01
#> C-A.989  1.163295e-01
#> C-A.990  1.163295e-01
#> C-A.991  1.172834e-01
#> C-A.992  1.182701e-01
#> C-A.993  1.182724e-01
#> C-A.994  1.183782e-01
#> C-A.995  1.186377e-01
#> C-A.996  1.187964e-01
#> C-A.997  1.190387e-01
#> C-A.998  1.192726e-01
#> C-A.999  1.194110e-01
#> C-A.1000 1.196705e-01
#> C-A.1001 1.196705e-01
#> C-A.1002 1.202582e-01
#> C-A.1003 1.202722e-01
#> C-A.1004 1.202722e-01
#> C-A.1005 1.202722e-01
#> C-A.1006 1.202722e-01
#> C-A.1007 1.204464e-01
#> C-A.1008 1.205229e-01
#> C-A.1009 1.207355e-01
#> C-A.1010 1.207355e-01
#> C-A.1011 1.210610e-01
#> C-A.1012 1.210610e-01
#> C-A.1013 1.210610e-01
#> C-A.1014 1.211908e-01
#> C-A.1015 1.212267e-01
#> C-A.1016 1.212267e-01
#> C-A.1017 1.212267e-01
#> C-A.1018 1.212267e-01
#> C-A.1019 1.212267e-01
#> C-A.1020 1.212267e-01
#> C-A.1021 1.212970e-01
#> C-A.1022 1.212970e-01
#> C-A.1023 1.218871e-01
#> C-A.1024 1.218871e-01
#> C-A.1025 1.218871e-01
#> C-A.1026 1.220307e-01
#> C-A.1027 1.221234e-01
#> C-A.1028 1.229448e-01
#> C-A.1029 1.229448e-01
#> C-A.1030 1.230025e-01
#> C-A.1031 1.234502e-01
#> C-A.1032 1.235526e-01
#> C-A.1033 1.235727e-01
#> C-A.1034 1.236778e-01
#> C-A.1035 1.241397e-01
#> C-A.1036 1.245514e-01
#> C-A.1037 1.245514e-01
#> C-A.1038 1.245694e-01
#> C-A.1039 1.253124e-01
#> C-A.1040 1.256990e-01
#> C-A.1041 1.260525e-01
#> C-A.1042 1.266264e-01
#> C-A.1043 1.268147e-01
#> C-A.1044 1.269051e-01
#> C-A.1045 1.270178e-01
#> C-A.1046 1.271315e-01
#> C-A.1047 1.271515e-01
#> C-A.1048 1.275743e-01
#> C-A.1049 1.278079e-01
#> C-A.1050 1.280060e-01
#> C-A.1051 1.280060e-01
#> C-A.1052 1.280448e-01
#> C-A.1053 1.280448e-01
#> C-A.1054 1.280448e-01
#> C-A.1055 1.282048e-01
#> C-A.1056 1.282635e-01
#> C-A.1057 1.284932e-01
#> C-A.1058 1.286054e-01
#> C-A.1059 1.287166e-01
#> C-A.1060 1.291003e-01
#> C-A.1061 1.297680e-01
#> C-A.1062 1.299706e-01
#> C-A.1063 1.300273e-01
#> C-A.1064 1.310191e-01
#> C-A.1065 1.310191e-01
#> C-A.1066 1.310191e-01
#> C-A.1067 1.313108e-01
#> C-A.1068 1.313108e-01
#> C-A.1069 1.316302e-01
#> C-A.1070 1.322236e-01
#> C-A.1071 1.322656e-01
#> C-A.1072 1.325027e-01
#> C-A.1073 1.325027e-01
#> C-A.1074 1.325027e-01
#> C-A.1075 1.327618e-01
#> C-A.1076 1.332694e-01
#> C-A.1077 1.337061e-01
#> C-A.1078 1.338574e-01
#> C-A.1079 1.338574e-01
#> C-A.1080 1.338890e-01
#> C-A.1081 1.338890e-01
#> C-A.1082 1.338890e-01
#> C-A.1083 1.338890e-01
#> C-A.1084 1.338890e-01
#> C-A.1085 1.340869e-01
#> C-A.1086 1.341917e-01
#> C-A.1087 1.343361e-01
#> C-A.1088 1.344699e-01
#> C-A.1089 1.344699e-01
#> C-A.1090 1.344699e-01
#> C-A.1091 1.345634e-01
#> C-A.1092 1.345999e-01
#> C-A.1093 1.353726e-01
#> C-A.1094 1.356846e-01
#> C-A.1095 1.356846e-01
#> C-A.1096 1.358036e-01
#> C-A.1097 1.361016e-01
#> C-A.1098 1.361349e-01
#> C-A.1099 1.361349e-01
#> C-A.1100 1.364290e-01
#> C-A.1101 1.365093e-01
#> C-A.1102 1.373586e-01
#> C-A.1103 1.373586e-01
#> C-A.1104 1.378978e-01
#> C-A.1105 1.380157e-01
#> C-A.1106 1.380157e-01
#> C-A.1107 1.380960e-01
#> C-A.1108 1.385383e-01
#> C-A.1109 1.388798e-01
#> C-A.1110 1.391615e-01
#> C-A.1111 1.392604e-01
#> C-A.1112 1.393916e-01
#> C-A.1113 1.400770e-01
#> C-A.1114 1.401300e-01
#> C-A.1115 1.410645e-01
#> C-A.1116 1.410645e-01
#> C-A.1117 1.410645e-01
#> C-A.1118 1.411414e-01
#> C-A.1119 1.417571e-01
#> C-A.1120 1.417571e-01
#> C-A.1121 1.426899e-01
#> C-A.1122 1.426899e-01
#> C-A.1123 1.428810e-01
#> C-A.1124 1.440362e-01
#> C-A.1125 1.440362e-01
#> C-A.1126 1.442915e-01
#> C-A.1127 1.444192e-01
#> C-A.1128 1.444435e-01
#> C-A.1129 1.446111e-01
#> C-A.1130 1.446111e-01
#> C-A.1131 1.446111e-01
#> C-A.1132 1.448534e-01
#> C-A.1133 1.449485e-01
#> C-A.1134 1.450846e-01
#> C-A.1135 1.451245e-01
#> C-A.1136 1.452626e-01
#> C-A.1137 1.455655e-01
#> C-A.1138 1.457453e-01
#> C-A.1139 1.457991e-01
#> C-A.1140 1.459314e-01
#> C-A.1141 1.462046e-01
#> C-A.1142 1.469113e-01
#> C-A.1143 1.469113e-01
#> C-A.1144 1.471149e-01
#> C-A.1145 1.471149e-01
#> C-A.1146 1.474028e-01
#> C-A.1147 1.474795e-01
#> C-A.1148 1.474795e-01
#> C-A.1149 1.474795e-01
#> C-A.1150 1.474795e-01
#> C-A.1151 1.474795e-01
#> C-A.1152 1.474795e-01
#> C-A.1153 1.474986e-01
#> C-A.1154 1.479071e-01
#> C-A.1155 1.479659e-01
#> C-A.1156 1.481647e-01
#> C-A.1157 1.483323e-01
#> C-A.1158 1.487037e-01
#> C-A.1159 1.488499e-01
#> C-A.1160 1.491062e-01
#> C-A.1161 1.493601e-01
#> C-A.1162 1.497688e-01
#> C-A.1163 1.503394e-01
#> C-A.1164 1.503394e-01
#> C-A.1165 1.506214e-01
#> C-A.1166 1.511576e-01
#> C-A.1167 1.515357e-01
#> C-A.1168 1.517154e-01
#> C-A.1169 1.517540e-01
#> C-A.1170 1.520253e-01
#> C-A.1171 1.520651e-01
#> C-A.1172 1.520651e-01
#> C-A.1173 1.521173e-01
#> C-A.1174 1.521198e-01
#> C-A.1175 1.521198e-01
#> C-A.1176 1.521198e-01
#> C-A.1177 1.521198e-01
#> C-A.1178 1.521198e-01
#> C-A.1179 1.524662e-01
#> C-A.1180 1.525685e-01
#> C-A.1181 1.525685e-01
#> C-A.1182 1.525685e-01
#> C-A.1183 1.525773e-01
#> C-A.1184 1.526272e-01
#> C-A.1185 1.527638e-01
#> C-A.1186 1.535942e-01
#> C-A.1187 1.535942e-01
#> C-A.1188 1.535942e-01
#> C-A.1189 1.539587e-01
#> C-A.1190 1.543548e-01
#> C-A.1191 1.548514e-01
#> C-A.1192 1.561879e-01
#> C-A.1193 1.565164e-01
#> C-A.1194 1.567891e-01
#> C-A.1195 1.567891e-01
#> C-A.1196 1.567891e-01
#> C-A.1197 1.567891e-01
#> C-A.1198 1.570107e-01
#> C-A.1199 1.570107e-01
#> C-A.1200 1.574744e-01
#> C-A.1201 1.577919e-01
#> C-A.1202 1.578655e-01
#> C-A.1203 1.578655e-01
#> C-A.1204 1.579102e-01
#> C-A.1205 1.585694e-01
#> C-A.1206 1.590450e-01
#> C-A.1207 1.592330e-01
#> C-A.1208 1.594076e-01
#> C-A.1209 1.596777e-01
#> C-A.1210 1.596777e-01
#> C-A.1211 1.596777e-01
#> C-A.1212 1.609735e-01
#> C-A.1213 1.609735e-01
#> C-A.1214 1.614563e-01
#> C-A.1215 1.618626e-01
#> C-A.1216 1.623802e-01
#> C-A.1217 1.637215e-01
#> C-A.1218 1.640848e-01
#> C-A.1219 1.643526e-01
#> C-A.1220 1.650902e-01
#> C-A.1221 1.653667e-01
#> C-A.1222 1.653667e-01
#> C-A.1223 1.653667e-01
#> C-A.1224 1.653667e-01
#> C-A.1225 1.653667e-01
#> C-A.1226 1.653667e-01
#> C-A.1227 1.653667e-01
#> C-A.1228 1.656427e-01
#> C-A.1229 1.658166e-01
#> C-A.1230 1.658166e-01
#> C-A.1231 1.658166e-01
#> C-A.1232 1.658638e-01
#> C-A.1233 1.658638e-01
#> C-A.1234 1.658638e-01
#> C-A.1235 1.659611e-01
#> C-A.1236 1.661462e-01
#> C-A.1237 1.681608e-01
#> C-A.1238 1.683811e-01
#> C-A.1239 1.689111e-01
#> C-A.1240 1.692609e-01
#> C-A.1241 1.692995e-01
#> C-A.1242 1.698218e-01
#> C-A.1243 1.703418e-01
#> C-A.1244 1.704346e-01
#> C-A.1245 1.707783e-01
#> C-A.1246 1.713379e-01
#> C-A.1247 1.718984e-01
#> C-A.1248 1.724479e-01
#> C-A.1249 1.727270e-01
#> C-A.1250 1.729295e-01
#> C-A.1251 1.735278e-01
#> C-A.1252 1.739074e-01
#> C-A.1253 1.741238e-01
#> C-A.1254 1.745538e-01
#> C-A.1255 1.750631e-01
#> C-A.1256 1.764862e-01
#> C-A.1257 1.766754e-01
#> C-A.1258 1.766754e-01
#> C-A.1259 1.766754e-01
#> C-A.1260 1.766754e-01
#> C-A.1261 1.766754e-01
#> C-A.1262 1.766754e-01
#> C-A.1263 1.766754e-01
#> C-A.1264 1.766754e-01
#> C-A.1265 1.766754e-01
#> C-A.1266 1.766754e-01
#> C-A.1267 1.775147e-01
#> C-A.1268 1.779297e-01
#> C-A.1269 1.787043e-01
#> C-A.1270 1.789050e-01
#> C-A.1271 1.791385e-01
#> C-A.1272 1.791385e-01
#> C-A.1273 1.793000e-01
#> C-A.1274 1.796484e-01
#> C-A.1275 1.799464e-01
#> C-A.1276 1.799886e-01
#> C-A.1277 1.805587e-01
#> C-A.1278 1.808906e-01
#> C-A.1279 1.811903e-01
#> C-A.1280 1.811903e-01
#> C-A.1281 1.812609e-01
#> C-A.1282 1.812609e-01
#> C-A.1283 1.812609e-01
#> C-A.1284 1.812609e-01
#> C-A.1285 1.812609e-01
#> C-A.1286 1.812609e-01
#> C-A.1287 1.815584e-01
#> C-A.1288 1.815954e-01
#> C-A.1289 1.815954e-01
#> C-A.1290 1.815954e-01
#> C-A.1291 1.815954e-01
#> C-A.1292 1.817846e-01
#> C-A.1293 1.818649e-01
#> C-A.1294 1.821207e-01
#> C-A.1295 1.821207e-01
#> C-A.1296 1.821207e-01
#> C-A.1297 1.824567e-01
#> C-A.1298 1.824767e-01
#> C-A.1299 1.824767e-01
#> C-A.1300 1.837977e-01
#> C-A.1301 1.843331e-01
#> C-A.1302 1.852942e-01
#> C-A.1303 1.852942e-01
#> C-A.1304 1.852942e-01
#> C-A.1305 1.852942e-01
#> C-A.1306 1.852942e-01
#> C-A.1307 1.852942e-01
#> C-A.1308 1.853417e-01
#> C-A.1309 1.853417e-01
#> C-A.1310 1.869218e-01
#> C-A.1311 1.870752e-01
#> C-A.1312 1.871834e-01
#> C-A.1313 1.876707e-01
#> C-A.1314 1.886374e-01
#> C-A.1315 1.888967e-01
#> C-A.1316 1.900655e-01
#> C-A.1317 1.902072e-01
#> C-A.1318 1.902072e-01
#> C-A.1319 1.906013e-01
#> C-A.1320 1.906501e-01
#> C-A.1321 1.917181e-01
#> C-A.1322 1.926481e-01
#> C-A.1323 1.927026e-01
#> C-A.1324 1.927026e-01
#> C-A.1325 1.927026e-01
#> C-A.1326 1.929979e-01
#> C-A.1327 1.929979e-01
#> C-A.1328 1.932976e-01
#> C-A.1329 1.938773e-01
#> C-A.1330 1.938773e-01
#> C-A.1331 1.951852e-01
#> C-A.1332 1.952415e-01
#> C-A.1333 1.953567e-01
#> C-A.1334 1.955758e-01
#> C-A.1335 1.956335e-01
#> C-A.1336 1.957914e-01
#> C-A.1337 1.957914e-01
#> C-A.1338 1.957914e-01
#> C-A.1339 1.957914e-01
#> C-A.1340 1.962676e-01
#> C-A.1341 1.963906e-01
#> C-A.1342 1.964289e-01
#> C-A.1343 1.971873e-01
#> C-A.1344 1.972162e-01
#> C-A.1345 1.972268e-01
#> C-A.1346 1.972268e-01
#> C-A.1347 1.972268e-01
#> C-A.1348 1.974418e-01
#> C-A.1349 1.974613e-01
#> C-A.1350 1.986050e-01
#> C-A.1351 1.987463e-01
#> C-A.1352 1.987463e-01
#> C-A.1353 1.989224e-01
#> B-C.1    4.723043e-05
#> B-C.2    5.718712e-04
#> B-C.3    5.718712e-04
#> B-C.4    5.718712e-04
#> B-C.5    5.718712e-04
#> B-C.6    5.718712e-04
#> B-C.7    5.718712e-04
#> B-C.8    5.718712e-04
#> B-C.9    5.718712e-04
#> B-C.10   5.718712e-04
#> B-C.11   5.718712e-04
#> B-C.12   5.718712e-04
#> B-C.13   5.718712e-04
#> B-C.14   5.718712e-04
#> B-C.15   5.779713e-04
#> B-C.16   5.779713e-04
#> B-C.17   7.094441e-04
#> B-C.18   9.803005e-04
#> B-C.19   1.038906e-03
#> B-C.20   1.163151e-03
#> B-C.21   1.179605e-03
#> B-C.22   1.323192e-03
#> B-C.23   1.323192e-03
#> B-C.24   1.425687e-03
#> B-C.25   1.660752e-03
#> B-C.26   1.660752e-03
#> B-C.27   1.701988e-03
#> B-C.28   1.701988e-03
#> B-C.29   1.713840e-03
#> B-C.30   1.716977e-03
#> B-C.31   1.780044e-03
#> B-C.32   2.556146e-03
#> B-C.33   2.719939e-03
#> B-C.34   2.781251e-03
#> B-C.35   2.781251e-03
#> B-C.36   2.910929e-03
#> B-C.37   3.090088e-03
#> B-C.38   3.304598e-03
#> B-C.39   3.304598e-03
#> B-C.40   3.304598e-03
#> B-C.41   3.409923e-03
#> B-C.42   3.615386e-03
#> B-C.43   3.661173e-03
#> B-C.44   3.661173e-03
#> B-C.45   4.693737e-03
#> B-C.46   5.237407e-03
#> B-C.47   5.237407e-03
#> B-C.48   5.934691e-03
#> B-C.49   6.063302e-03
#> B-C.50   6.406983e-03
#> B-C.51   7.058445e-03
#> B-C.52   7.058445e-03
#> B-C.53   7.058445e-03
#> B-C.54   7.058445e-03
#> B-C.55   7.656978e-03
#> B-C.56   8.317205e-03
#> B-C.57   8.317205e-03
#> B-C.58   8.317205e-03
#> B-C.59   8.618446e-03
#> B-C.60   8.881229e-03
#> B-C.61   9.008826e-03
#> B-C.62   1.001351e-02
#> B-C.63   1.034206e-02
#> B-C.64   1.034206e-02
#> B-C.65   1.034206e-02
#> B-C.66   1.051874e-02
#> B-C.67   1.051874e-02
#> B-C.68   1.139226e-02
#> B-C.69   1.175213e-02
#> B-C.70   1.191511e-02
#> B-C.71   1.191511e-02
#> B-C.72   1.191511e-02
#> B-C.73   1.196444e-02
#> B-C.74   1.205624e-02
#> B-C.75   1.210391e-02
#> B-C.76   1.210391e-02
#> B-C.77   1.210391e-02
#> B-C.78   1.256461e-02
#> B-C.79   1.279957e-02
#> B-C.80   1.279957e-02
#> B-C.81   1.319424e-02
#> B-C.82   1.329449e-02
#> B-C.83   1.347684e-02
#> B-C.84   1.351507e-02
#> B-C.85   1.372448e-02
#> B-C.86   1.391571e-02
#> B-C.87   1.504563e-02
#> B-C.88   1.766962e-02
#> B-C.89   1.786144e-02
#> B-C.90   1.806977e-02
#> B-C.91   1.843855e-02
#> B-C.92   1.859949e-02
#> B-C.93   1.859949e-02
#> B-C.94   1.859949e-02
#> B-C.95   1.892117e-02
#> B-C.96   1.925105e-02
#> B-C.97   2.005033e-02
#> B-C.98   2.023104e-02
#> B-C.99   2.099181e-02
#> B-C.100  2.109150e-02
#> B-C.101  2.185777e-02
#> B-C.102  2.185777e-02
#> B-C.103  2.235665e-02
#> B-C.104  2.235665e-02
#> B-C.105  2.282983e-02
#> B-C.106  2.282983e-02
#> B-C.107  2.282983e-02
#> B-C.108  2.337976e-02
#> B-C.109  2.337976e-02
#> B-C.110  2.347620e-02
#> B-C.111  2.347620e-02
#> B-C.112  2.357183e-02
#> B-C.113  2.357183e-02
#> B-C.114  2.362230e-02
#> B-C.115  2.380804e-02
#> B-C.116  2.380804e-02
#> B-C.117  2.481775e-02
#> B-C.118  2.483683e-02
#> B-C.119  2.523850e-02
#> B-C.120  2.523850e-02
#> B-C.121  2.663656e-02
#> B-C.122  2.663656e-02
#> B-C.123  2.708594e-02
#> B-C.124  2.708594e-02
#> B-C.125  2.710394e-02
#> B-C.126  2.738488e-02
#> B-C.127  2.752343e-02
#> B-C.128  2.799750e-02
#> B-C.129  2.807301e-02
#> B-C.130  2.895717e-02
#> B-C.131  2.912754e-02
#> B-C.132  2.926201e-02
#> B-C.133  2.926201e-02
#> B-C.134  2.926201e-02
#> B-C.135  2.942152e-02
#> B-C.136  2.968064e-02
#> B-C.137  3.214227e-02
#> B-C.138  3.214227e-02
#> B-C.139  3.214227e-02
#> B-C.140  3.214227e-02
#> B-C.141  3.246157e-02
#> B-C.142  3.268840e-02
#> B-C.143  3.270368e-02
#> B-C.144  3.294054e-02
#> B-C.145  3.294054e-02
#> B-C.146  3.294054e-02
#> B-C.147  3.294054e-02
#> B-C.148  3.294054e-02
#> B-C.149  3.296596e-02
#> B-C.150  3.334171e-02
#> B-C.151  3.369058e-02
#> B-C.152  3.395863e-02
#> B-C.153  3.409242e-02
#> B-C.154  3.468185e-02
#> B-C.155  3.468185e-02
#> B-C.156  3.468185e-02
#> B-C.157  3.468185e-02
#> B-C.158  3.627773e-02
#> B-C.159  3.648219e-02
#> B-C.160  3.695144e-02
#> B-C.161  3.716865e-02
#> B-C.162  3.730509e-02
#> B-C.163  3.730509e-02
#> B-C.164  3.730509e-02
#> B-C.165  3.796068e-02
#> B-C.166  3.919493e-02
#> B-C.167  3.932434e-02
#> B-C.168  3.932434e-02
#> B-C.169  3.932434e-02
#> B-C.170  3.932434e-02
#> B-C.171  3.968186e-02
#> B-C.172  3.968186e-02
#> B-C.173  4.090461e-02
#> B-C.174  4.093325e-02
#> B-C.175  4.156105e-02
#> B-C.176  4.289088e-02
#> B-C.177  4.354355e-02
#> B-C.178  4.354355e-02
#> B-C.179  4.419065e-02
#> B-C.180  4.419065e-02
#> B-C.181  4.419065e-02
#> B-C.182  4.420458e-02
#> B-C.183  4.468793e-02
#> B-C.184  4.510564e-02
#> B-C.185  4.510564e-02
#> B-C.186  4.553268e-02
#> B-C.187  4.553268e-02
#> B-C.188  4.553268e-02
#> B-C.189  4.553268e-02
#> B-C.190  4.553268e-02
#> B-C.191  4.553268e-02
#> B-C.192  4.600070e-02
#> B-C.193  4.644595e-02
#> B-C.194  4.690720e-02
#> B-C.195  4.752528e-02
#> B-C.196  4.759136e-02
#> B-C.197  4.786212e-02
#> B-C.198  4.786212e-02
#> B-C.199  4.786212e-02
#> B-C.200  4.786212e-02
#> B-C.201  4.799494e-02
#> B-C.202  4.930533e-02
#> B-C.203  4.943424e-02
#> B-C.204  4.947165e-02
#> B-C.205  4.947165e-02
#> B-C.206  5.004831e-02
#> B-C.207  5.085316e-02
#> B-C.208  5.085316e-02
#> B-C.209  5.085316e-02
#> B-C.210  5.131045e-02
#> B-C.211  5.136449e-02
#> B-C.212  5.180989e-02
#> B-C.213  5.180989e-02
#> B-C.214  5.223822e-02
#> B-C.215  5.365269e-02
#> B-C.216  5.373337e-02
#> B-C.217  5.373337e-02
#> B-C.218  5.377181e-02
#> B-C.219  5.377181e-02
#> B-C.220  5.377181e-02
#> B-C.221  5.502442e-02
#> B-C.222  5.502442e-02
#> B-C.223  5.576936e-02
#> B-C.224  5.587544e-02
#> B-C.225  5.615542e-02
#> B-C.226  5.627436e-02
#> B-C.227  5.627436e-02
#> B-C.228  5.627436e-02
#> B-C.229  5.627436e-02
#> B-C.230  5.676556e-02
#> B-C.231  5.676556e-02
#> B-C.232  5.676556e-02
#> B-C.233  5.676556e-02
#> B-C.234  5.676556e-02
#> B-C.235  5.676556e-02
#> B-C.236  5.676556e-02
#> B-C.237  5.706548e-02
#> B-C.238  5.863393e-02
#> B-C.239  5.873297e-02
#> B-C.240  5.879164e-02
#> B-C.241  5.880154e-02
#> B-C.242  5.880154e-02
#> B-C.243  5.880154e-02
#> B-C.244  5.880154e-02
#> B-C.245  5.910461e-02
#> B-C.246  6.283888e-02
#> B-C.247  6.283888e-02
#> B-C.248  6.345183e-02
#> B-C.249  6.345183e-02
#> B-C.250  6.345183e-02
#> B-C.251  6.356823e-02
#> B-C.252  6.356823e-02
#> B-C.253  6.362984e-02
#> B-C.254  6.365467e-02
#> B-C.255  6.368974e-02
#> B-C.256  6.439623e-02
#> B-C.257  6.462380e-02
#> B-C.258  6.604465e-02
#> B-C.259  6.609147e-02
#> B-C.260  6.677271e-02
#> B-C.261  6.677622e-02
#> B-C.262  6.677622e-02
#> B-C.263  6.677622e-02
#> B-C.264  6.688064e-02
#> B-C.265  6.688064e-02
#> B-C.266  6.688064e-02
#> B-C.267  6.712591e-02
#> B-C.268  6.717843e-02
#> B-C.269  6.766134e-02
#> B-C.270  6.837727e-02
#> B-C.271  6.854549e-02
#> B-C.272  6.922861e-02
#> B-C.273  6.922861e-02
#> B-C.274  6.956342e-02
#> B-C.275  7.019562e-02
#> B-C.276  7.064015e-02
#> B-C.277  7.096619e-02
#> B-C.278  7.097466e-02
#> B-C.279  7.152043e-02
#> B-C.280  7.153858e-02
#> B-C.281  7.153858e-02
#> B-C.282  7.153858e-02
#> B-C.283  7.158294e-02
#> B-C.284  7.219119e-02
#> B-C.285  7.219119e-02
#> B-C.286  7.272392e-02
#> B-C.287  7.272392e-02
#> B-C.288  7.272392e-02
#> B-C.289  7.272392e-02
#> B-C.290  7.272392e-02
#> B-C.291  7.272392e-02
#> B-C.292  7.370470e-02
#> B-C.293  7.374229e-02
#> B-C.294  7.374229e-02
#> B-C.295  7.374229e-02
#> B-C.296  7.374229e-02
#> B-C.297  7.417004e-02
#> B-C.298  7.420895e-02
#> B-C.299  7.434657e-02
#> B-C.300  7.448305e-02
#> B-C.301  7.538681e-02
#> B-C.302  7.615399e-02
#> B-C.303  7.793388e-02
#> B-C.304  7.943685e-02
#> B-C.305  7.973527e-02
#> B-C.306  8.039349e-02
#> B-C.307  8.059906e-02
#> B-C.308  8.059906e-02
#> B-C.309  8.059906e-02
#> B-C.310  8.059906e-02
#> B-C.311  8.059906e-02
#> B-C.312  8.059906e-02
#> B-C.313  8.059906e-02
#> B-C.314  8.059906e-02
#> B-C.315  8.111532e-02
#> B-C.316  8.111532e-02
#> B-C.317  8.111532e-02
#> B-C.318  8.111532e-02
#> B-C.319  8.213058e-02
#> B-C.320  8.358902e-02
#> B-C.321  8.555405e-02
#> B-C.322  8.680962e-02
#> B-C.323  8.809915e-02
#> B-C.324  8.809915e-02
#> B-C.325  8.809915e-02
#> B-C.326  8.826386e-02
#> B-C.327  8.872594e-02
#> B-C.328  8.872594e-02
#> B-C.329  8.872594e-02
#> B-C.330  8.872594e-02
#> B-C.331  8.872594e-02
#> B-C.332  8.872594e-02
#> B-C.333  8.889515e-02
#> B-C.334  8.922183e-02
#> B-C.335  8.922183e-02
#> B-C.336  9.079626e-02
#> B-C.337  9.119212e-02
#> B-C.338  9.154729e-02
#> B-C.339  9.165093e-02
#> B-C.340  9.309824e-02
#> B-C.341  9.350168e-02
#> B-C.342  9.464140e-02
#> B-C.343  9.464140e-02
#> B-C.344  9.464140e-02
#> B-C.345  9.464140e-02
#> B-C.346  9.464140e-02
#> B-C.347  9.464140e-02
#> B-C.348  9.464140e-02
#> B-C.349  9.464140e-02
#> B-C.350  9.464140e-02
#> B-C.351  9.464140e-02
#> B-C.352  9.464140e-02
#> B-C.353  9.491534e-02
#> B-C.354  9.520288e-02
#> B-C.355  9.539243e-02
#> B-C.356  9.539243e-02
#> B-C.357  9.539243e-02
#> B-C.358  9.539243e-02
#> B-C.359  9.539243e-02
#> B-C.360  9.633563e-02
#> B-C.361  9.663052e-02
#> B-C.362  9.748528e-02
#> B-C.363  9.889417e-02
#> B-C.364  9.891522e-02
#> B-C.365  9.960502e-02
#> B-C.366  9.964623e-02
#> B-C.367  9.999874e-02
#> B-C.368  9.999874e-02
#> B-C.369  9.999874e-02
#> B-C.370  9.999874e-02
#> B-C.371  9.999874e-02
#> B-C.372  1.002283e-01
#> B-C.373  1.010512e-01
#> B-C.374  1.023536e-01
#> B-C.375  1.034162e-01
#> B-C.376  1.034162e-01
#> B-C.377  1.034162e-01
#> B-C.378  1.038065e-01
#> B-C.379  1.038065e-01
#> B-C.380  1.038065e-01
#> B-C.381  1.044744e-01
#> B-C.382  1.045758e-01
#> B-C.383  1.052423e-01
#> B-C.384  1.052423e-01
#> B-C.385  1.054429e-01
#> B-C.386  1.054429e-01
#> B-C.387  1.061690e-01
#> B-C.388  1.063169e-01
#> B-C.389  1.075299e-01
#> B-C.390  1.075543e-01
#> B-C.391  1.082568e-01
#> B-C.392  1.087730e-01
#> B-C.393  1.093693e-01
#> B-C.394  1.095150e-01
#> B-C.395  1.095150e-01
#> B-C.396  1.096611e-01
#> B-C.397  1.115009e-01
#> B-C.398  1.115009e-01
#> B-C.399  1.122351e-01
#> B-C.400  1.123098e-01
#> B-C.401  1.129673e-01
#> B-C.402  1.129673e-01
#> B-C.403  1.129673e-01
#> B-C.404  1.134539e-01
#> B-C.405  1.140921e-01
#> B-C.406  1.145045e-01
#> B-C.407  1.150212e-01
#> B-C.408  1.150212e-01
#> B-C.409  1.154418e-01
#> B-C.410  1.154418e-01
#> B-C.411  1.154418e-01
#> B-C.412  1.154418e-01
#> B-C.413  1.154757e-01
#> B-C.414  1.161432e-01
#> B-C.415  1.161432e-01
#> B-C.416  1.161780e-01
#> B-C.417  1.161780e-01
#> B-C.418  1.161780e-01
#> B-C.419  1.163089e-01
#> B-C.420  1.163089e-01
#> B-C.421  1.168054e-01
#> B-C.422  1.169710e-01
#> B-C.423  1.172022e-01
#> B-C.424  1.172022e-01
#> B-C.425  1.175253e-01
#> B-C.426  1.178772e-01
#> B-C.427  1.178772e-01
#> B-C.428  1.184204e-01
#> B-C.429  1.184204e-01
#> B-C.430  1.185167e-01
#> B-C.431  1.185167e-01
#> B-C.432  1.185167e-01
#> B-C.433  1.186111e-01
#> B-C.434  1.186111e-01
#> B-C.435  1.188215e-01
#> B-C.436  1.188215e-01
#> B-C.437  1.188215e-01
#> B-C.438  1.205891e-01
#> B-C.439  1.213934e-01
#> B-C.440  1.214381e-01
#> B-C.441  1.234883e-01
#> B-C.442  1.235397e-01
#> B-C.443  1.235397e-01
#> B-C.444  1.235397e-01
#> B-C.445  1.238583e-01
#> B-C.446  1.241721e-01
#> B-C.447  1.252803e-01
#> B-C.448  1.257441e-01
#> B-C.449  1.257540e-01
#> B-C.450  1.258764e-01
#> B-C.451  1.260958e-01
#> B-C.452  1.269726e-01
#> B-C.453  1.269726e-01
#> B-C.454  1.269726e-01
#> B-C.455  1.269726e-01
#> B-C.456  1.269726e-01
#> B-C.457  1.269944e-01
#> B-C.458  1.269944e-01
#> B-C.459  1.269944e-01
#> B-C.460  1.269944e-01
#> B-C.461  1.269944e-01
#> B-C.462  1.269944e-01
#> B-C.463  1.269944e-01
#> B-C.464  1.269944e-01
#> B-C.465  1.269944e-01
#> B-C.466  1.269944e-01
#> B-C.467  1.269944e-01
#> B-C.468  1.269944e-01
#> B-C.469  1.271208e-01
#> B-C.470  1.271410e-01
#> B-C.471  1.271410e-01
#> B-C.472  1.271410e-01
#> B-C.473  1.271410e-01
#> B-C.474  1.271410e-01
#> B-C.475  1.271410e-01
#> B-C.476  1.271410e-01
#> B-C.477  1.271410e-01
#> B-C.478  1.271410e-01
#> B-C.479  1.271410e-01
#> B-C.480  1.294743e-01
#> B-C.481  1.294743e-01
#> B-C.482  1.294743e-01
#> B-C.483  1.294743e-01
#> B-C.484  1.296079e-01
#> B-C.485  1.304249e-01
#> B-C.486  1.304249e-01
#> B-C.487  1.305909e-01
#> B-C.488  1.310863e-01
#> B-C.489  1.310863e-01
#> B-C.490  1.310863e-01
#> B-C.491  1.310863e-01
#> B-C.492  1.310863e-01
#> B-C.493  1.310863e-01
#> B-C.494  1.311593e-01
#> B-C.495  1.311593e-01
#> B-C.496  1.312496e-01
#> B-C.497  1.320095e-01
#> B-C.498  1.320095e-01
#> B-C.499  1.320868e-01
#> B-C.500  1.330596e-01
#> B-C.501  1.330596e-01
#> B-C.502  1.338167e-01
#> B-C.503  1.338167e-01
#> B-C.504  1.338167e-01
#> B-C.505  1.340107e-01
#> B-C.506  1.342597e-01
#> B-C.507  1.347608e-01
#> B-C.508  1.352231e-01
#> B-C.509  1.354572e-01
#> B-C.510  1.354572e-01
#> B-C.511  1.356531e-01
#> B-C.512  1.362608e-01
#> B-C.513  1.363840e-01
#> B-C.514  1.366436e-01
#> B-C.515  1.366436e-01
#> B-C.516  1.369125e-01
#> B-C.517  1.371518e-01
#> B-C.518  1.376750e-01
#> B-C.519  1.384597e-01
#> B-C.520  1.399968e-01
#> B-C.521  1.401964e-01
#> B-C.522  1.401964e-01
#> B-C.523  1.401964e-01
#> B-C.524  1.403204e-01
#> B-C.525  1.404224e-01
#> B-C.526  1.408679e-01
#> B-C.527  1.408679e-01
#> B-C.528  1.409477e-01
#> B-C.529  1.410457e-01
#> B-C.530  1.411025e-01
#> B-C.531  1.411025e-01
#> B-C.532  1.411025e-01
#> B-C.533  1.411803e-01
#> B-C.534  1.411803e-01
#> B-C.535  1.411803e-01
#> B-C.536  1.411803e-01
#> B-C.537  1.411803e-01
#> B-C.538  1.412099e-01
#> B-C.539  1.428056e-01
#> B-C.540  1.428056e-01
#> B-C.541  1.440195e-01
#> B-C.542  1.440195e-01
#> B-C.543  1.440195e-01
#> B-C.544  1.440195e-01
#> B-C.545  1.440195e-01
#> B-C.546  1.440195e-01
#> B-C.547  1.440195e-01
#> B-C.548  1.440195e-01
#> B-C.549  1.440195e-01
#> B-C.550  1.444904e-01
#> B-C.551  1.460507e-01
#> B-C.552  1.460507e-01
#> B-C.553  1.460507e-01
#> B-C.554  1.467126e-01
#> B-C.555  1.471392e-01
#> B-C.556  1.471546e-01
#> B-C.557  1.472057e-01
#> B-C.558  1.475370e-01
#> B-C.559  1.475451e-01
#> B-C.560  1.476868e-01
#> B-C.561  1.476868e-01
#> B-C.562  1.477519e-01
#> B-C.563  1.477519e-01
#> B-C.564  1.477519e-01
#> B-C.565  1.477519e-01
#> B-C.566  1.477519e-01
#> B-C.567  1.477519e-01
#> B-C.568  1.477519e-01
#> B-C.569  1.487648e-01
#> B-C.570  1.488781e-01
#> B-C.571  1.488781e-01
#> B-C.572  1.497504e-01
#> B-C.573  1.509421e-01
#> B-C.574  1.509421e-01
#> B-C.575  1.520345e-01
#> B-C.576  1.527653e-01
#> B-C.577  1.527653e-01
#> B-C.578  1.527653e-01
#> B-C.579  1.529969e-01
#> B-C.580  1.538663e-01
#> B-C.581  1.540819e-01
#> B-C.582  1.545374e-01
#> B-C.583  1.553480e-01
#> B-C.584  1.560435e-01
#> B-C.585  1.560435e-01
#> B-C.586  1.560435e-01
#> B-C.587  1.562324e-01
#> B-C.588  1.569039e-01
#> B-C.589  1.571671e-01
#> B-C.590  1.576583e-01
#> B-C.591  1.579064e-01
#> B-C.592  1.579064e-01
#> B-C.593  1.580306e-01
#> B-C.594  1.582773e-01
#> B-C.595  1.593445e-01
#> B-C.596  1.593445e-01
#> B-C.597  1.593445e-01
#> B-C.598  1.593445e-01
#> B-C.599  1.593445e-01
#> B-C.600  1.595747e-01
#> B-C.601  1.595747e-01
#> B-C.602  1.595747e-01
#> B-C.603  1.598724e-01
#> B-C.604  1.598724e-01
#> B-C.605  1.599728e-01
#> B-C.606  1.605165e-01
#> B-C.607  1.608144e-01
#> B-C.608  1.609257e-01
#> B-C.609  1.609257e-01
#> B-C.610  1.610021e-01
#> B-C.611  1.618050e-01
#> B-C.612  1.624935e-01
#> B-C.613  1.627269e-01
#> B-C.614  1.627269e-01
#> B-C.615  1.630097e-01
#> B-C.616  1.632204e-01
#> B-C.617  1.641190e-01
#> B-C.618  1.641722e-01
#> B-C.619  1.643464e-01
#> B-C.620  1.643464e-01
#> B-C.621  1.645264e-01
#> B-C.622  1.647295e-01
#> B-C.623  1.668047e-01
#> B-C.624  1.668047e-01
#> B-C.625  1.668047e-01
#> B-C.626  1.668047e-01
#> B-C.627  1.670580e-01
#> B-C.628  1.670580e-01
#> B-C.629  1.674550e-01
#> B-C.630  1.676758e-01
#> B-C.631  1.676758e-01
#> B-C.632  1.678462e-01
#> B-C.633  1.678462e-01
#> B-C.634  1.686810e-01
#> B-C.635  1.694585e-01
#> B-C.636  1.694822e-01
#> B-C.637  1.702733e-01
#> B-C.638  1.721810e-01
#> B-C.639  1.721810e-01
#> B-C.640  1.721810e-01
#> B-C.641  1.724350e-01
#> B-C.642  1.725247e-01
#> B-C.643  1.725247e-01
#> B-C.644  1.744839e-01
#> B-C.645  1.744839e-01
#> B-C.646  1.749532e-01
#> B-C.647  1.749567e-01
#> B-C.648  1.752686e-01
#> B-C.649  1.760344e-01
#> B-C.650  1.760344e-01
#> B-C.651  1.763550e-01
#> B-C.652  1.763550e-01
#> B-C.653  1.763550e-01
#> B-C.654  1.763814e-01
#> B-C.655  1.767608e-01
#> B-C.656  1.770685e-01
#> B-C.657  1.776480e-01
#> B-C.658  1.780035e-01
#> B-C.659  1.780910e-01
#> B-C.660  1.782526e-01
#> B-C.661  1.788787e-01
#> B-C.662  1.798194e-01
#> B-C.663  1.798194e-01
#> B-C.664  1.802355e-01
#> B-C.665  1.803708e-01
#> B-C.666  1.817948e-01
#> B-C.667  1.817948e-01
#> B-C.668  1.817948e-01
#> B-C.669  1.817948e-01
#> B-C.670  1.818561e-01
#> B-C.671  1.823843e-01
#> B-C.672  1.823843e-01
#> B-C.673  1.835734e-01
#> B-C.674  1.836272e-01
#> B-C.675  1.836272e-01
#> B-C.676  1.836272e-01
#> B-C.677  1.836272e-01
#> B-C.678  1.836272e-01
#> B-C.679  1.836272e-01
#> B-C.680  1.851874e-01
#> B-C.681  1.851874e-01
#> B-C.682  1.851874e-01
#> B-C.683  1.851874e-01
#> B-C.684  1.851874e-01
#> B-C.685  1.851874e-01
#> B-C.686  1.851874e-01
#> B-C.687  1.856828e-01
#> B-C.688  1.856828e-01
#> B-C.689  1.856828e-01
#> B-C.690  1.857550e-01
#> B-C.691  1.862772e-01
#> B-C.692  1.862772e-01
#> B-C.693  1.865334e-01
#> B-C.694  1.865334e-01
#> B-C.695  1.865334e-01
#> B-C.696  1.865334e-01
#> B-C.697  1.865334e-01
#> B-C.698  1.866756e-01
#> B-C.699  1.871682e-01
#> B-C.700  1.873518e-01
#> B-C.701  1.876213e-01
#> B-C.702  1.892348e-01
#> B-C.703  1.924920e-01
#> B-C.704  1.934307e-01
#> B-C.705  1.934911e-01
#> B-C.706  1.934911e-01
#> B-C.707  1.946114e-01
#> B-C.708  1.946114e-01
#> B-C.709  1.959022e-01
#> B-C.710  1.959022e-01
#> B-C.711  1.959022e-01
#> B-C.712  1.961403e-01
#> B-C.713  1.967434e-01
#> B-C.714  1.970282e-01
#> B-C.715  1.977537e-01
#> B-C.716  1.978271e-01
#> B-C.717  1.979225e-01
#> B-C.718  1.979225e-01
#> B-C.719  1.979225e-01
#> B-C.720  1.979225e-01
#> B-C.721  1.979225e-01
#> B-C.722  1.979594e-01
#> B-C.723  1.982140e-01
#> B-C.724  1.983897e-01
#> B-C.725  1.985335e-01
#> B-C.726  1.998210e-01
#> B-C.727  1.998210e-01
#> B-C.728  1.998210e-01
#> B-C.729  1.998210e-01
#> B-C.730  1.998210e-01
#> B-C.731  1.998210e-01
#> B-C.732  1.998210e-01
#> B-C.733  1.998210e-01
#> B-C.734  1.998210e-01
#> B-C.735  1.998210e-01
#> B-C.736  1.998210e-01
#> B-C.737  1.998210e-01
#> B-C.738  1.998210e-01
#> B-C.739  1.998210e-01
#> B-C.740  2.004047e-01
#> B-C.741  2.010637e-01
#> B-C.742  2.019225e-01
#> B-C.743  2.019225e-01
#> B-C.744  2.021961e-01
#> B-C.745  2.030770e-01
#> B-C.746  2.033346e-01
#> B-C.747  2.033346e-01
#> B-C.748  2.038435e-01
#> B-C.749  2.050146e-01
#> B-C.750  2.054091e-01
#> B-C.751  2.061191e-01
#> B-C.752  2.061191e-01
#> B-C.753  2.061191e-01
#> B-C.754  2.061191e-01
#> B-C.755  2.061191e-01
#> B-C.756  2.061191e-01
#> B-C.757  2.061191e-01
#> B-C.758  2.061191e-01
#> B-C.759  2.061191e-01
#> B-C.760  2.061191e-01
#> B-C.761  2.070226e-01
#> B-C.762  2.070226e-01
#> B-C.763  2.071350e-01
#> B-C.764  2.080714e-01
#> B-C.765  2.080714e-01
#> B-C.766  2.086167e-01
#> B-C.767  2.087269e-01
#> B-C.768  2.089805e-01
#> B-C.769  2.089805e-01
#> B-C.770  2.097395e-01
#> B-C.771  2.097395e-01
#> B-C.772  2.111832e-01
#> B-C.773  2.111832e-01
#> B-C.774  2.111832e-01
#> B-C.775  2.114531e-01
#> B-C.776  2.119893e-01
#> B-C.777  2.120615e-01
#> B-C.778  2.126268e-01
#> B-C.779  2.127364e-01
#> B-C.780  2.136612e-01
#> B-C.781  2.138056e-01
#> B-C.782  2.139364e-01
#> B-C.783  2.147595e-01
#> B-C.784  2.148046e-01
#> B-C.785  2.154224e-01
#> B-C.786  2.154760e-01
#> B-C.787  2.154760e-01
#> B-C.788  2.154760e-01
#> B-C.789  2.164548e-01
#> B-C.790  2.164548e-01
#> B-C.791  2.164548e-01
#> B-C.792  2.164548e-01
#> B-C.793  2.167087e-01
#> B-C.794  2.173930e-01
#> B-C.795  2.173930e-01
#> B-C.796  2.177135e-01
#> B-C.797  2.177135e-01
#> B-C.798  2.182220e-01
#> B-C.799  2.184676e-01
#> B-C.800  2.189217e-01
#> B-C.801  2.191675e-01
#> B-C.802  2.191675e-01
#> B-C.803  2.191675e-01
#> B-C.804  2.191675e-01
#> B-C.805  2.191675e-01
#> B-C.806  2.191675e-01
#> B-C.807  2.191675e-01
#> B-C.808  2.207384e-01
#> B-C.809  2.239016e-01
#> B-C.810  2.243775e-01
#> B-C.811  2.243775e-01
#> B-C.812  2.243775e-01
#> B-C.813  2.246881e-01
#> B-C.814  2.247319e-01
#> B-C.815  2.247319e-01
#> B-C.816  2.247429e-01
#> B-C.817  2.247429e-01
#> B-C.818  2.247429e-01
#> B-C.819  2.247962e-01
#> B-C.820  2.247962e-01
#> B-C.821  2.262949e-01
#> B-C.822  2.269211e-01
#> B-C.823  2.270226e-01
#> B-C.824  2.278664e-01
#> B-C.825  2.281397e-01
#> B-C.826  2.283536e-01
#> B-C.827  2.286652e-01
#> B-C.828  2.287178e-01
#> B-C.829  2.287721e-01
#> B-C.830  2.288901e-01
#> B-C.831  2.291492e-01
#> B-C.832  2.293242e-01
#> B-C.833  2.293242e-01
#> B-C.834  2.293242e-01
#> B-C.835  2.295238e-01
#> B-C.836  2.303098e-01
#> B-C.837  2.303098e-01
#> B-C.838  2.319292e-01
#> B-C.839  2.319354e-01
#> B-C.840  2.323628e-01
#> B-C.841  2.323628e-01
#> B-C.842  2.327949e-01
#> B-C.843  2.331070e-01
#> B-C.844  2.332738e-01
#> B-C.845  2.332738e-01
#> B-C.846  2.332738e-01
#> B-C.847  2.332738e-01
#> B-C.848  2.335842e-01
#> B-C.849  2.336993e-01
#> B-C.850  2.338710e-01
#> B-C.851  2.339175e-01
#> B-C.852  2.339796e-01
#> B-C.853  2.346719e-01
#> B-C.854  2.353817e-01
#> B-C.855  2.356901e-01
#> B-C.856  2.367854e-01
#> B-C.857  2.367854e-01
#> B-C.858  2.367854e-01
#> B-C.859  2.367854e-01
#> B-C.860  2.368393e-01
#> B-C.861  2.372638e-01
#> B-C.862  2.373519e-01
#> B-C.863  2.373519e-01
#> B-C.864  2.373519e-01
#> B-C.865  2.373519e-01
#> B-C.866  2.373519e-01
#> B-C.867  2.373519e-01
#> B-C.868  2.373519e-01
#> B-C.869  2.373519e-01
#> B-C.870  2.373519e-01
#> B-C.871  2.379829e-01
#> B-C.872  2.379829e-01
#> B-C.873  2.379829e-01
#> B-C.874  2.379829e-01
#> B-C.875  2.379829e-01
#> B-C.876  2.379829e-01
#> B-C.877  2.379829e-01
#> B-C.878  2.379829e-01
#> B-C.879  2.394890e-01
#> B-C.880  2.401666e-01
#> B-C.881  2.401666e-01
#> B-C.882  2.401666e-01
#> B-C.883  2.401666e-01
#> B-C.884  2.403614e-01
#> B-C.885  2.403614e-01
#> B-C.886  2.405065e-01
#> B-C.887  2.408935e-01
#> B-C.888  2.412801e-01
#> B-C.889  2.416155e-01
#> B-C.890  2.418038e-01
#> B-C.891  2.418038e-01
#> B-C.892  2.434973e-01
#> B-C.893  2.434973e-01
#> B-C.894  2.437976e-01
#> B-C.895  2.437976e-01
#> B-C.896  2.437976e-01
#> B-C.897  2.438873e-01
#> B-C.898  2.439922e-01
#> B-C.899  2.442939e-01
#> B-C.900  2.455866e-01
#> B-C.901  2.456330e-01
#> B-C.902  2.458253e-01
#> B-C.903  2.458253e-01
#> B-C.904  2.463670e-01
#> B-C.905  2.465758e-01
#> B-C.906  2.475588e-01
#> B-C.907  2.475588e-01
#> B-C.908  2.480696e-01
#> B-C.909  2.489780e-01
#> B-C.910  2.498772e-01
#> B-C.911  2.498772e-01
#> B-C.912  2.498772e-01
#> B-C.913  2.498772e-01
#> B-C.914  2.501373e-01
#> B-C.915  2.505499e-01
#> B-C.916  2.505499e-01
#> B-C.917  2.513721e-01
#> B-C.918  2.515946e-01
#> B-C.919  2.515946e-01
#> B-C.920  2.534179e-01
#> B-C.921  2.537715e-01
#> B-C.922  2.537715e-01
#> B-C.923  2.539335e-01
#> B-C.924  2.540105e-01
#> B-C.925  2.541337e-01
#> B-C.926  2.541368e-01
#> B-C.927  2.546071e-01
#> B-C.928  2.546071e-01
#> B-C.929  2.554808e-01
#> B-C.930  2.554808e-01
#> B-C.931  2.554808e-01
#> B-C.932  2.554808e-01
#> B-C.933  2.555881e-01
#> B-C.934  2.555881e-01
#> B-C.935  2.569323e-01
#> B-C.936  2.569323e-01
#> B-C.937  2.579945e-01
#> B-C.938  2.581358e-01
#> B-C.939  2.581358e-01
#> B-C.940  2.581358e-01
#> B-C.941  2.605783e-01
#> B-C.942  2.634089e-01
#> B-C.943  2.634089e-01
#> B-C.944  2.634368e-01
#> B-C.945  2.634368e-01
#> B-C.946  2.640366e-01
#> B-C.947  2.645966e-01
#> B-C.948  2.650439e-01
#> B-C.949  2.651227e-01
#> B-C.950  2.651227e-01
#> B-C.951  2.651227e-01
#> B-C.952  2.651227e-01
#> B-C.953  2.657606e-01
#> B-C.954  2.657606e-01
#> B-C.955  2.658069e-01
#> B-C.956  2.658069e-01
#> B-C.957  2.658069e-01
#> B-C.958  2.658069e-01
#> B-C.959  2.663879e-01
#> B-C.960  2.668295e-01
#> B-C.961  2.668295e-01
#> B-C.962  2.674008e-01
```
