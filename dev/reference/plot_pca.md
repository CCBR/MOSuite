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

PCA is sensitive to feature scale. For raw, cleaned, filtered, or
CPM-like count data, setting `log_transform = TRUE` compresses
high-count features before PCA so ordination is less dominated by
abundance differences. Do not enable this for already
normalized/log-scale data, such as voom/logCPM `norm` counts or
batch-corrected normalized counts, because that would apply a second log
transform. The default PCA transform uses the original MOSuite PCA
scaling, `log(x + 0.5)`.

See the low-level function docs for additional arguments depending on
whether you're plotting 2 or 3 PCs:

- [`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md) -
  used when there are **2** principal components

- [`plot_pca_3d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_3d.md) -
  used when there are **3** principal components

## See also

- [`plot_pca.multiOmicDataSet()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.multiOmicDataSet.md)

- [`plot_pca.data.frame()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.data.frame.md)

Other plotters:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`print_or_save_plot()`](https://ccbr.github.io/MOSuite/dev/reference/print_or_save_plot.md)

Other PCA functions:
[`calc_pca()`](https://ccbr.github.io/MOSuite/dev/reference/calc_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_3d.md)

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md)

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


# 3D
plot_pca(moo, count_type = "clean", principal_components = c(1, 2, 3))

{"x":{"visdat":{"1ac72457b5e6":["function () ","plotlyVisDat"]},"cur_data":"1ac72457b5e6","attrs":{"1ac72457b5e6":{"x":{},"y":{},"z":{},"mode":"markers","marker":{"size":8},"hoverinfo":"text","text":{},"color":{},"size":24,"colors":["#5954d6","#e1562c","#b80058"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter3d"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"PC1"},"yaxis":{"title":"PC2"},"zaxis":{"title":"PC3"}},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-40.624166881644548,-56.213316043361097,-69.107071102044515],"y":[25.229726861915786,6.1338577161232601,-21.895234510691349],"z":[-9.1146448284552086,-29.509826980510748,-36.358202516253542],"mode":"markers","marker":{"color":"rgba(89,84,214,1)","size":8,"sizemode":"area","line":{"color":"rgba(89,84,214,1)"}},"hoverinfo":["text","text","text"],"text":["Label: A1<br>Group: A","Label: A2<br>Group: A","Label: A3<br>Group: A"],"type":"scatter3d","name":"A","textfont":{"color":"rgba(89,84,214,1)","size":55},"error_y":{"color":"rgba(89,84,214,1)","width":55},"error_x":{"color":"rgba(89,84,214,1)","width":55},"line":{"color":"rgba(89,84,214,1)","width":55},"frame":null},{"x":[-36.16602512157435,-25.865255255389862,-9.6232450176944173],"y":[7.8050429797890217,-11.213808049471634,9.3272469604228387],"z":[18.215438247891832,30.862641867877986,58.980343929227899],"mode":"markers","marker":{"color":"rgba(225,86,44,1)","size":8,"sizemode":"area","line":{"color":"rgba(225,86,44,1)"}},"hoverinfo":["text","text","text"],"text":["Label: B1<br>Group: B","Label: B2<br>Group: B","Label: B3<br>Group: B"],"type":"scatter3d","name":"B","textfont":{"color":"rgba(225,86,44,1)","size":55},"error_y":{"color":"rgba(225,86,44,1)","width":55},"error_x":{"color":"rgba(225,86,44,1)","width":55},"line":{"color":"rgba(225,86,44,1)","width":55},"frame":null},{"x":[74.334557668027173,85.044222680901271,78.220299072780293],"y":[-86.72868022299302,117.99234054350596,-46.650492278600879],"z":[-73.882686673690046,-31.798500917774529,72.605437871686334],"mode":"markers","marker":{"color":"rgba(184,0,88,1)","size":8,"sizemode":"area","line":{"color":"rgba(184,0,88,1)"}},"hoverinfo":["text","text","text"],"text":["Label: C1<br>Group: C","Label: C2<br>Group: C","Label: C3<br>Group: C"],"type":"scatter3d","name":"C","textfont":{"color":"rgba(184,0,88,1)","size":55},"error_y":{"color":"rgba(184,0,88,1)","width":55},"error_x":{"color":"rgba(184,0,88,1)","width":55},"line":{"color":"rgba(184,0,88,1)","width":55},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
# dataframe
plot_pca(nidap_clean_raw_counts,
  sample_metadata = nidap_sample_metadata,
  principal_components = c(1, 2)
)

```
