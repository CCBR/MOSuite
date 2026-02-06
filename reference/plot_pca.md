# Perform and plot a Principal Components Analysis

Perform and plot a Principal Components Analysis

## Usage

``` r
plot_pca(moo_counts, principal_components = c(1, 2), ...)
```

## Arguments

- moo_counts:

  counts dataframe or `multiOmicDataSet` containing `count_type` &
  `sub_count_type` in the counts slot

- principal_components:

  vector with numbered principal components to plot. Use 2 for a 2D pca
  with ggplot, or 3 for a 3D pca with plotly. (Default: `c(1,2)`)

- ...:

  additional arguments forwarded to method (see Details below)

## Value

PCA plot (2D or 3D depending on the number of `principal_components`)

## Details

See the low-level function docs for additional arguments depending on
whether you're plotting 2 or 3 PCs:

- [plot_pca_2d](https://ccbr.github.io/MOSuite/reference/plot_pca_2d.md) -
  used when there are **2** principal components

- [plot_pca_3d](https://ccbr.github.io/MOSuite/reference/plot_pca_3d.md) -
  used when there are **3** principal components

## Methods

|                                                                          |                    |
|--------------------------------------------------------------------------|--------------------|
| link to docs                                                             | class              |
| [plot_pca_moo](https://ccbr.github.io/MOSuite/reference/plot_pca_moo.md) | `multiOmicDataSet` |
| [plot_pca_dat](https://ccbr.github.io/MOSuite/reference/plot_pca_dat.md) | `data.frame`       |

## See also

Other plotters:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`print_or_save_plot()`](https://ccbr.github.io/MOSuite/reference/print_or_save_plot.md)

Other PCA functions:
[`calc_pca()`](https://ccbr.github.io/MOSuite/reference/calc_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_2d.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_3d.md)

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)

## Examples

``` r
# multiOmicDataSet
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = nidap_raw_counts,
    "clean" = nidap_clean_raw_counts
  )
)
plot_pca(moo, count_type = "clean", principal_components = c(1, 2))
#> Saving 6.67 x 6.67 in image


# 3D
plot_pca(moo, count_type = "clean", principal_components = c(1, 2, 3))

{"x":{"visdat":{"1c6d1be3dacf":["function () ","plotlyVisDat"]},"cur_data":"1c6d1be3dacf","attrs":{"1c6d1be3dacf":{"x":{},"y":{},"z":{},"mode":"markers","marker":{"size":8},"hoverinfo":"text","text":{},"color":{},"size":24,"colors":["#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter3d"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"PC1"},"yaxis":{"title":"PC2"},"zaxis":{"title":"PC3"}},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-40.624166881644534,-56.213316043361083,-69.107071102044515],"y":[25.229726861915882,6.1338577161232672,-21.895234510691218],"z":[-9.1146448284550949,-29.509826980510688,-36.358202516253556],"mode":"markers","marker":{"color":"rgba(89,84,214,1)","size":8,"sizemode":"area","line":{"color":"rgba(89,84,214,1)"}},"hoverinfo":["text","text","text"],"text":["A1","A2","A3"],"type":"scatter3d","name":"A","textfont":{"color":"rgba(89,84,214,1)","size":55},"error_y":{"color":"rgba(89,84,214,1)","width":55},"error_x":{"color":"rgba(89,84,214,1)","width":55},"line":{"color":"rgba(89,84,214,1)","width":55},"frame":null},{"x":[-36.166025121574357,-25.865255255389851,-9.6232450176944013],"y":[7.8050429797890049,-11.213808049471595,9.3272469604228885],"z":[18.215438247891754,30.86264186787793,58.980343929227793],"mode":"markers","marker":{"color":"rgba(166,163,121,1)","size":8,"sizemode":"area","line":{"color":"rgba(166,163,121,1)"}},"hoverinfo":["text","text","text"],"text":["B1","B2","B3"],"type":"scatter3d","name":"B","textfont":{"color":"rgba(166,163,121,1)","size":55},"error_y":{"color":"rgba(166,163,121,1)","width":55},"error_x":{"color":"rgba(166,163,121,1)","width":55},"line":{"color":"rgba(166,163,121,1)","width":55},"frame":null},{"x":[74.334557668027145,85.044222680901512,78.220299072780236],"y":[-86.72868022299312,117.99234054350583,-46.650492278600858],"z":[-73.882686673689975,-31.798500917774614,72.605437871686476],"mode":"markers","marker":{"color":"rgba(135,133,0,1)","size":8,"sizemode":"area","line":{"color":"rgba(135,133,0,1)"}},"hoverinfo":["text","text","text"],"text":["C1","C2","C3"],"type":"scatter3d","name":"C","textfont":{"color":"rgba(135,133,0,1)","size":55},"error_y":{"color":"rgba(135,133,0,1)","width":55},"error_x":{"color":"rgba(135,133,0,1)","width":55},"line":{"color":"rgba(135,133,0,1)","width":55},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
# dataframe
plot_pca(nidap_clean_raw_counts,
  sample_metadata = nidap_sample_metadata,
  principal_components = c(1, 2)
)
#> Saving 6.67 x 6.67 in image

```
