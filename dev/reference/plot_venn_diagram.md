# Plot a venn diagram, UpSet plot, or table of intersections

generates Venn diagram of intersections across a series of sets (e.g.,
intersections of significant genes across tested contrasts). This Venn
diagram is available for up to five sets; Intersection plot is available
for any number of sets. Specific sets can be selected for the
visualizations and the returned dataset may include all (default) or
specified intersections.

## Usage

``` r
plot_venn_diagram(
  diff_summary_dat,
  feature_id_colname = NULL,
  contrasts_colname = "Contrast",
  select_contrasts = c(),
  plot_type = "Venn diagram",
  intersection_ids = c(),
  venn_force_unique = TRUE,
  venn_numbers_format = "raw",
  venn_significant_digits = 2,
  venn_fill_colors = c("darkgoldenrod2", "darkolivegreen2", "mediumpurple3",
    "darkorange2", "lightgreen"),
  venn_fill_transparency = 0.2,
  venn_border_colors = "fill colors",
  venn_font_size_for_category_names = 3,
  venn_category_names_distance = c(),
  venn_category_names_position = c(),
  venn_font_size_for_counts = 6,
  venn_outer_margin = 0,
  intersections_order = "degree",
  display_empty_intersections = FALSE,
  intersection_bar_color = "steelblue4",
  intersection_point_size = 2.2,
  intersection_line_width = 0.7,
  table_font_size = 0.7,
  table_content = "all intersections",
  graphics_device = grDevices::png,
  dpi = 300,
  image_width = 4000,
  image_height = 3000,
  plot_filename = "venn_diagram.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff"
)
```

## Arguments

- diff_summary_dat:

  Summarized differential expression analysis

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- contrasts_colname:

  Name of the column in `diff_summary_dat` that contains the contrast
  names (default: "Contrast")

- select_contrasts:

  A vector of contrast names to select for the plot. If empty, all
  contrasts are used.

- plot_type:

  Type of plot to generate: "Venn diagram" or "Intersection plot".
  Default: "Venn diagram"

- intersection_ids:

  A vector of intersection IDs to select for the plot. If empty, all
  intersections are used.

- venn_force_unique:

  If TRUE, forces unique elements in the Venn diagram. Default: TRUE

- venn_numbers_format:

  Format for the numbers in the Venn diagram. Options: "raw", "percent",
  "raw-percent", "percent-raw". Default: "raw"

- venn_significant_digits:

  Number of significant digits for the Venn diagram numbers. Default: 2

- venn_fill_colors:

  A vector of colors to fill the Venn diagram categories. Default:
  c("darkgoldenrod2", "darkolivegreen2", "mediumpurple3", "darkorange2",
  "lightgreen")

- venn_fill_transparency:

  Transparency level for the Venn diagram fill colors. Default: 0.2

- venn_border_colors:

  Colors for the borders of the Venn diagram categories. Default: "fill
  colors" (uses the same colors as `venn_fill_colors`)

- venn_font_size_for_category_names:

  Font size for the category names in the Venn diagram. Default: 3

- venn_category_names_distance:

  Distance of the category names from the Venn diagram circles. Default:
  c()

- venn_category_names_position:

  Position of the category names in the Venn diagram. Default: c()

- venn_font_size_for_counts:

  Font size for the counts in the Venn diagram. Default: 6

- venn_outer_margin:

  Outer margin for the Venn diagram. Default: 0

- intersections_order:

  Order of the intersections in the plot. Default: "by size"

- display_empty_intersections:

  If TRUE, displays empty intersections in the plot. Default: FALSE

- intersection_bar_color:

  Color for the intersection bars in the plot. Default: "lightgray"

- intersection_point_size:

  Size of the points in the intersection plot. Default: 2

- intersection_line_width:

  Width of the lines in the intersection plot. Default: 0.5

- table_font_size:

  Font size for the table in the plot. Default: 3

- table_content:

  Content of the table in the plot. Default: NULL

- graphics_device:

  passed to `ggsave(device)`. Default:
  [`grDevices::png`](https://rdrr.io/r/grDevices/png.html)

- dpi:

  dots-per-inch of the output image (see `ggsave()`) - only used if
  save_plots is TRUE

- image_width:

  output image width in pixels - only used if save_plots is TRUE

- image_height:

  output image height in pixels - only used if save_plots is TRUE

- plot_filename:

  plot output filename - only used if save_plots is TRUE

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- plots_subdir:

  subdirectory in where plots will be saved if `save_plots` is `TRUE`

## Examples

``` r
plot_venn_diagram(nidap_volcano_summary_dat, print_plots = TRUE)
#> All intersections: 1:7,c(1, 2, 3, 4, 5, 6, 7),c(163, 237, 518, 780, 225, 379, 766),c("Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes")
#> Intersections returned: 1:7,c(1, 2, 3, 4, 5, 6, 7),c(163, 237, 518, 780, 225, 379, 766)

#>                Gene      Intersection Id Size
#> 1              Dntt (B-A ∩ B-C ∩ C-A)  1  163
#> 2              Flt3 (B-A ∩ B-C ∩ C-A)  1  163
#> 3              Perp (B-A ∩ B-C ∩ C-A)  1  163
#> 4              Chad (B-A ∩ B-C ∩ C-A)  1  163
#> 5              Tox2 (B-A ∩ B-C ∩ C-A)  1  163
#> 6               Id2 (B-A ∩ B-C ∩ C-A)  1  163
#> 7              Rxrg (B-A ∩ B-C ∩ C-A)  1  163
#> 8              Emp1 (B-A ∩ B-C ∩ C-A)  1  163
#> 9              Rag2 (B-A ∩ B-C ∩ C-A)  1  163
#> 10           Nlrp1a (B-A ∩ B-C ∩ C-A)  1  163
#> 11          Trbj1-1 (B-A ∩ B-C ∩ C-A)  1  163
#> 12            Adam8 (B-A ∩ B-C ∩ C-A)  1  163
#> 13            Eltd1 (B-A ∩ B-C ∩ C-A)  1  163
#> 14          Tcrg-C4 (B-A ∩ B-C ∩ C-A)  1  163
#> 15            Trbc1 (B-A ∩ B-C ∩ C-A)  1  163
#> 16              Txk (B-A ∩ B-C ∩ C-A)  1  163
#> 17           S100a4 (B-A ∩ B-C ∩ C-A)  1  163
#> 18            Spon2 (B-A ∩ B-C ∩ C-A)  1  163
#> 19            Cass4 (B-A ∩ B-C ∩ C-A)  1  163
#> 20             Cmah (B-A ∩ B-C ∩ C-A)  1  163
#> 21            Dusp6 (B-A ∩ B-C ∩ C-A)  1  163
#> 22         Serpine2 (B-A ∩ B-C ∩ C-A)  1  163
#> 23            Runx3 (B-A ∩ B-C ∩ C-A)  1  163
#> 24            Lpar1 (B-A ∩ B-C ∩ C-A)  1  163
#> 25             Nav1 (B-A ∩ B-C ∩ C-A)  1  163
#> 26             Lig4 (B-A ∩ B-C ∩ C-A)  1  163
#> 27              Myc (B-A ∩ B-C ∩ C-A)  1  163
#> 28            Mef2c (B-A ∩ B-C ∩ C-A)  1  163
#> 29            Prr13 (B-A ∩ B-C ∩ C-A)  1  163
#> 30            Anxa2 (B-A ∩ B-C ∩ C-A)  1  163
#> 31         Trbv13-2 (B-A ∩ B-C ∩ C-A)  1  163
#> 32            Trat1 (B-A ∩ B-C ∩ C-A)  1  163
#> 33           Spock2 (B-A ∩ B-C ∩ C-A)  1  163
#> 34            Aqp11 (B-A ∩ B-C ∩ C-A)  1  163
#> 35             Zpbp (B-A ∩ B-C ∩ C-A)  1  163
#> 36             Rarg (B-A ∩ B-C ∩ C-A)  1  163
#> 37             Cd34 (B-A ∩ B-C ∩ C-A)  1  163
#> 38           Apcdd1 (B-A ∩ B-C ∩ C-A)  1  163
#> 39             Ahsg (B-A ∩ B-C ∩ C-A)  1  163
#> 40          Afap1l1 (B-A ∩ B-C ∩ C-A)  1  163
#> 41            Esyt1 (B-A ∩ B-C ∩ C-A)  1  163
#> 42            Mgat1 (B-A ∩ B-C ∩ C-A)  1  163
#> 43              Ckb (B-A ∩ B-C ∩ C-A)  1  163
#> 44              Itk (B-A ∩ B-C ∩ C-A)  1  163
#> 45           Ifitm2 (B-A ∩ B-C ∩ C-A)  1  163
#> 46          Tcrg-C1 (B-A ∩ B-C ∩ C-A)  1  163
#> 47            Oasl2 (B-A ∩ B-C ∩ C-A)  1  163
#> 48          Hspbap1 (B-A ∩ B-C ∩ C-A)  1  163
#> 49             Icos (B-A ∩ B-C ∩ C-A)  1  163
#> 50            Kcng1 (B-A ∩ B-C ∩ C-A)  1  163
#> 51             Rora (B-A ∩ B-C ∩ C-A)  1  163
#> 52            Gata3 (B-A ∩ B-C ∩ C-A)  1  163
#> 53          Tmem119 (B-A ∩ B-C ∩ C-A)  1  163
#> 54            Il12a (B-A ∩ B-C ∩ C-A)  1  163
#> 55            Arl11 (B-A ∩ B-C ∩ C-A)  1  163
#> 56             Thy1 (B-A ∩ B-C ∩ C-A)  1  163
#> 57            P2rx7 (B-A ∩ B-C ∩ C-A)  1  163
#> 58          Tcrg-V4 (B-A ∩ B-C ∩ C-A)  1  163
#> 59           Cox6a2 (B-A ∩ B-C ∩ C-A)  1  163
#> 60             Lmo4 (B-A ∩ B-C ∩ C-A)  1  163
#> 61            Fnbp1 (B-A ∩ B-C ∩ C-A)  1  163
#> 62             Cd74 (B-A ∩ B-C ∩ C-A)  1  163
#> 63            Ero1l (B-A ∩ B-C ∩ C-A)  1  163
#> 64          Tcrg-V1 (B-A ∩ B-C ∩ C-A)  1  163
#> 65            Ikzf2 (B-A ∩ B-C ∩ C-A)  1  163
#> 66            Prdx4 (B-A ∩ B-C ∩ C-A)  1  163
#> 67             Shpk (B-A ∩ B-C ∩ C-A)  1  163
#> 68             Pask (B-A ∩ B-C ∩ C-A)  1  163
#> 69            Sh2d5 (B-A ∩ B-C ∩ C-A)  1  163
#> 70          Emilin2 (B-A ∩ B-C ∩ C-A)  1  163
#> 71            Mpeg1 (B-A ∩ B-C ∩ C-A)  1  163
#> 72            Itgb3 (B-A ∩ B-C ∩ C-A)  1  163
#> 73            Trbc2 (B-A ∩ B-C ∩ C-A)  1  163
#> 74            H2-Aa (B-A ∩ B-C ∩ C-A)  1  163
#> 75           Gpr183 (B-A ∩ B-C ∩ C-A)  1  163
#> 76             Cpa3 (B-A ∩ B-C ∩ C-A)  1  163
#> 77             Snx9 (B-A ∩ B-C ∩ C-A)  1  163
#> 78             Gapt (B-A ∩ B-C ∩ C-A)  1  163
#> 79             Cttn (B-A ∩ B-C ∩ C-A)  1  163
#> 80             Lef1 (B-A ∩ B-C ∩ C-A)  1  163
#> 81           Tom1l1 (B-A ∩ B-C ∩ C-A)  1  163
#> 82             Zeb2 (B-A ∩ B-C ∩ C-A)  1  163
#> 83          Tcrg-V5 (B-A ∩ B-C ∩ C-A)  1  163
#> 84            Pdcd1 (B-A ∩ B-C ∩ C-A)  1  163
#> 85           Ccp110 (B-A ∩ B-C ∩ C-A)  1  163
#> 86           Tspan2 (B-A ∩ B-C ∩ C-A)  1  163
#> 87             Ccr2 (B-A ∩ B-C ∩ C-A)  1  163
#> 88           Ifitm1 (B-A ∩ B-C ∩ C-A)  1  163
#> 89            Abcb9 (B-A ∩ B-C ∩ C-A)  1  163
#> 90          Pik3ap1 (B-A ∩ B-C ∩ C-A)  1  163
#> 91            Cela1 (B-A ∩ B-C ∩ C-A)  1  163
#> 92           Gnptab (B-A ∩ B-C ∩ C-A)  1  163
#> 93          Fam124b (B-A ∩ B-C ∩ C-A)  1  163
#> 94           Cdkn1a (B-A ∩ B-C ∩ C-A)  1  163
#> 95           Rassf4 (B-A ∩ B-C ∩ C-A)  1  163
#> 96             Eya2 (B-A ∩ B-C ∩ C-A)  1  163
#> 97          Slc15a1 (B-A ∩ B-C ∩ C-A)  1  163
#> 98              Btk (B-A ∩ B-C ∩ C-A)  1  163
#> 99            Ssbp3 (B-A ∩ B-C ∩ C-A)  1  163
#> 100          Gprin3 (B-A ∩ B-C ∩ C-A)  1  163
#> 101           Stag1 (B-A ∩ B-C ∩ C-A)  1  163
#> 102          Podnl1 (B-A ∩ B-C ∩ C-A)  1  163
#> 103           Trgj4 (B-A ∩ B-C ∩ C-A)  1  163
#> 104          Wfdc17 (B-A ∩ B-C ∩ C-A)  1  163
#> 105           Bcas1 (B-A ∩ B-C ∩ C-A)  1  163
#> 106         Gramd1b (B-A ∩ B-C ∩ C-A)  1  163
#> 107         Gm26517 (B-A ∩ B-C ∩ C-A)  1  163
#> 108         Slc29a3 (B-A ∩ B-C ∩ C-A)  1  163
#> 109           Trgj1 (B-A ∩ B-C ∩ C-A)  1  163
#> 110            Fsd1 (B-A ∩ B-C ∩ C-A)  1  163
#> 111            Ccr9 (B-A ∩ B-C ∩ C-A)  1  163
#> 112         Siglecg (B-A ∩ B-C ∩ C-A)  1  163
#> 113         Slc16a3 (B-A ∩ B-C ∩ C-A)  1  163
#> 114            Faah (B-A ∩ B-C ∩ C-A)  1  163
#> 115            Aff4 (B-A ∩ B-C ∩ C-A)  1  163
#> 116           Csf1r (B-A ∩ B-C ∩ C-A)  1  163
#> 117            Dll1 (B-A ∩ B-C ∩ C-A)  1  163
#> 118           Cd226 (B-A ∩ B-C ∩ C-A)  1  163
#> 119           Gcnt1 (B-A ∩ B-C ∩ C-A)  1  163
#> 120          Rbpms2 (B-A ∩ B-C ∩ C-A)  1  163
#> 121         B4galt4 (B-A ∩ B-C ∩ C-A)  1  163
#> 122            Asb2 (B-A ∩ B-C ∩ C-A)  1  163
#> 123          Lpcat4 (B-A ∩ B-C ∩ C-A)  1  163
#> 124           Furin (B-A ∩ B-C ∩ C-A)  1  163
#> 125            Esr1 (B-A ∩ B-C ∩ C-A)  1  163
#> 126            Gbp2 (B-A ∩ B-C ∩ C-A)  1  163
#> 127           Ighj1 (B-A ∩ B-C ∩ C-A)  1  163
#> 128           Stim1 (B-A ∩ B-C ∩ C-A)  1  163
#> 129          Maged1 (B-A ∩ B-C ∩ C-A)  1  163
#> 130          Spata6 (B-A ∩ B-C ∩ C-A)  1  163
#> 131            Lpxn (B-A ∩ B-C ∩ C-A)  1  163
#> 132           Itm2c (B-A ∩ B-C ∩ C-A)  1  163
#> 133           Hoxa9 (B-A ∩ B-C ∩ C-A)  1  163
#> 134         Ehbp1l1 (B-A ∩ B-C ∩ C-A)  1  163
#> 135         Trbj2-1 (B-A ∩ B-C ∩ C-A)  1  163
#> 136        Tmem176a (B-A ∩ B-C ∩ C-A)  1  163
#> 137           Myo1c (B-A ∩ B-C ∩ C-A)  1  163
#> 138           Fads2 (B-A ∩ B-C ∩ C-A)  1  163
#> 139          Pdlim2 (B-A ∩ B-C ∩ C-A)  1  163
#> 140         Zfp36l1 (B-A ∩ B-C ∩ C-A)  1  163
#> 141         Slc39a4 (B-A ∩ B-C ∩ C-A)  1  163
#> 142           Anks1 (B-A ∩ B-C ∩ C-A)  1  163
#> 143          Cd300a (B-A ∩ B-C ∩ C-A)  1  163
#> 144            Ice2 (B-A ∩ B-C ∩ C-A)  1  163
#> 145          Rnase4 (B-A ∩ B-C ∩ C-A)  1  163
#> 146           Kdm7a (B-A ∩ B-C ∩ C-A)  1  163
#> 147           Limk1 (B-A ∩ B-C ∩ C-A)  1  163
#> 148          Stxbp6 (B-A ∩ B-C ∩ C-A)  1  163
#> 149            Rorc (B-A ∩ B-C ∩ C-A)  1  163
#> 150          Lrrc52 (B-A ∩ B-C ∩ C-A)  1  163
#> 151         Spata13 (B-A ∩ B-C ∩ C-A)  1  163
#> 152           Stk39 (B-A ∩ B-C ∩ C-A)  1  163
#> 153           Myl10 (B-A ∩ B-C ∩ C-A)  1  163
#> 154            Drc7 (B-A ∩ B-C ∩ C-A)  1  163
#> 155          Mageh1 (B-A ∩ B-C ∩ C-A)  1  163
#> 156         Map3k12 (B-A ∩ B-C ∩ C-A)  1  163
#> 157            Ccr8 (B-A ∩ B-C ∩ C-A)  1  163
#> 158           Prtn3 (B-A ∩ B-C ∩ C-A)  1  163
#> 159           Ikzf3 (B-A ∩ B-C ∩ C-A)  1  163
#> 160       Serpinb6b (B-A ∩ B-C ∩ C-A)  1  163
#> 161         Dennd4c (B-A ∩ B-C ∩ C-A)  1  163
#> 162           Bach1 (B-A ∩ B-C ∩ C-A)  1  163
#> 163          Ahcyl2 (B-A ∩ B-C ∩ C-A)  1  163
#> 164         Tspan13       (B-A ∩ B-C)  2  237
#> 165           Nfkb1       (B-A ∩ B-C)  2  237
#> 166           Itga4       (B-A ∩ B-C)  2  237
#> 167           Skp1a       (B-A ∩ B-C)  2  237
#> 168          Samsn1       (B-A ∩ B-C)  2  237
#> 169       Serpinb1a       (B-A ∩ B-C)  2  237
#> 170             Emb       (B-A ∩ B-C)  2  237
#> 171          Cyp2j6       (B-A ∩ B-C)  2  237
#> 172            Grb2       (B-A ∩ B-C)  2  237
#> 173            Nrp2       (B-A ∩ B-C)  2  237
#> 174          Lrrc49       (B-A ∩ B-C)  2  237
#> 175          Itga2b       (B-A ∩ B-C)  2  237
#> 176             Cfp       (B-A ∩ B-C)  2  237
#> 177        Tmem176b       (B-A ∩ B-C)  2  237
#> 178           Fgf11       (B-A ∩ B-C)  2  237
#> 179           Mcf2l       (B-A ∩ B-C)  2  237
#> 180           Syne2       (B-A ∩ B-C)  2  237
#> 181           Cpne2       (B-A ∩ B-C)  2  237
#> 182           Mdga1       (B-A ∩ B-C)  2  237
#> 183           Trpv4       (B-A ∩ B-C)  2  237
#> 184           Padi2       (B-A ∩ B-C)  2  237
#> 185            Mycn       (B-A ∩ B-C)  2  237
#> 186          Zfp105       (B-A ∩ B-C)  2  237
#> 187           Rogdi       (B-A ∩ B-C)  2  237
#> 188            Rnd3       (B-A ∩ B-C)  2  237
#> 189          Clec4d       (B-A ∩ B-C)  2  237
#> 190            Cd81       (B-A ∩ B-C)  2  237
#> 191            Rgs1       (B-A ∩ B-C)  2  237
#> 192           Usp54       (B-A ∩ B-C)  2  237
#> 193           Cpeb4       (B-A ∩ B-C)  2  237
#> 194          Gpr141       (B-A ∩ B-C)  2  237
#> 195            Vill       (B-A ∩ B-C)  2  237
#> 196         Napepld       (B-A ∩ B-C)  2  237
#> 197          Cdkn1c       (B-A ∩ B-C)  2  237
#> 198          Fam26f       (B-A ∩ B-C)  2  237
#> 199         Gm15518       (B-A ∩ B-C)  2  237
#> 200           Adrb2       (B-A ∩ B-C)  2  237
#> 201           Itsn1       (B-A ∩ B-C)  2  237
#> 202         Clec4a2       (B-A ∩ B-C)  2  237
#> 203          Cyp2j9       (B-A ∩ B-C)  2  237
#> 204          Tyrobp       (B-A ∩ B-C)  2  237
#> 205          Tas1r1       (B-A ∩ B-C)  2  237
#> 206            Aif1       (B-A ∩ B-C)  2  237
#> 207           Iigp1       (B-A ∩ B-C)  2  237
#> 208          Gm4955       (B-A ∩ B-C)  2  237
#> 209         Col23a1       (B-A ∩ B-C)  2  237
#> 210          Ptger2       (B-A ∩ B-C)  2  237
#> 211        Fam171a2       (B-A ∩ B-C)  2  237
#> 212           Mpzl1       (B-A ∩ B-C)  2  237
#> 213           Tmtc3       (B-A ∩ B-C)  2  237
#> 214         Siglecf       (B-A ∩ B-C)  2  237
#> 215           Sytl2       (B-A ∩ B-C)  2  237
#> 216             Vit       (B-A ∩ B-C)  2  237
#> 217           Bfsp2       (B-A ∩ B-C)  2  237
#> 218         Arhgef3       (B-A ∩ B-C)  2  237
#> 219            Bbs1       (B-A ∩ B-C)  2  237
#> 220            Rbpj       (B-A ∩ B-C)  2  237
#> 221           Zmat1       (B-A ∩ B-C)  2  237
#> 222            Cd53       (B-A ∩ B-C)  2  237
#> 223          Il17re       (B-A ∩ B-C)  2  237
#> 224            Ehd3       (B-A ∩ B-C)  2  237
#> 225            Accs       (B-A ∩ B-C)  2  237
#> 226         Tnfaip1       (B-A ∩ B-C)  2  237
#> 227              Hp       (B-A ∩ B-C)  2  237
#> 228         Gm20257       (B-A ∩ B-C)  2  237
#> 229           Tgtp2       (B-A ∩ B-C)  2  237
#> 230            Cish       (B-A ∩ B-C)  2  237
#> 231             Ak8       (B-A ∩ B-C)  2  237
#> 232          Cnrip1       (B-A ∩ B-C)  2  237
#> 233             Npl       (B-A ∩ B-C)  2  237
#> 234          Man1c1       (B-A ∩ B-C)  2  237
#> 235           Nfil3       (B-A ∩ B-C)  2  237
#> 236          Rwdd2b       (B-A ∩ B-C)  2  237
#> 237             Gal       (B-A ∩ B-C)  2  237
#> 238           Rab44       (B-A ∩ B-C)  2  237
#> 239            Ctso       (B-A ∩ B-C)  2  237
#> 240           Klhl4       (B-A ∩ B-C)  2  237
#> 241         Trim34b       (B-A ∩ B-C)  2  237
#> 242         Pik3c2a       (B-A ∩ B-C)  2  237
#> 243             Gsn       (B-A ∩ B-C)  2  237
#> 244          Chst10       (B-A ∩ B-C)  2  237
#> 245            Yes1       (B-A ∩ B-C)  2  237
#> 246            Gm2a       (B-A ∩ B-C)  2  237
#> 247            Klf3       (B-A ∩ B-C)  2  237
#> 248          Sh3bp5       (B-A ∩ B-C)  2  237
#> 249            Sdc1       (B-A ∩ B-C)  2  237
#> 250           Ckap5       (B-A ∩ B-C)  2  237
#> 251        Arhgap12       (B-A ∩ B-C)  2  237
#> 252             Cpq       (B-A ∩ B-C)  2  237
#> 253          Ankrd6       (B-A ∩ B-C)  2  237
#> 254          Prdm10       (B-A ∩ B-C)  2  237
#> 255         Camsap2       (B-A ∩ B-C)  2  237
#> 256         Ccdc176       (B-A ∩ B-C)  2  237
#> 257            Polq       (B-A ∩ B-C)  2  237
#> 258          Camk1d       (B-A ∩ B-C)  2  237
#> 259          Trim58       (B-A ∩ B-C)  2  237
#> 260         Rhobtb3       (B-A ∩ B-C)  2  237
#> 261          Map4k5       (B-A ∩ B-C)  2  237
#> 262          Gm5141       (B-A ∩ B-C)  2  237
#> 263           Vipr2       (B-A ∩ B-C)  2  237
#> 264          Ttc30b       (B-A ∩ B-C)  2  237
#> 265           Dapp1       (B-A ∩ B-C)  2  237
#> 266          Ppfia1       (B-A ∩ B-C)  2  237
#> 267           Kalrn       (B-A ∩ B-C)  2  237
#> 268          Cyp2r1       (B-A ∩ B-C)  2  237
#> 269          Camk2d       (B-A ∩ B-C)  2  237
#> 270            Eps8       (B-A ∩ B-C)  2  237
#> 271         Gm16194       (B-A ∩ B-C)  2  237
#> 272           Ttc19       (B-A ∩ B-C)  2  237
#> 273         Ppp1r18       (B-A ∩ B-C)  2  237
#> 274            Abi2       (B-A ∩ B-C)  2  237
#> 275         Rps6ka2       (B-A ∩ B-C)  2  237
#> 276           Fbxo4       (B-A ∩ B-C)  2  237
#> 277            Cd63       (B-A ∩ B-C)  2  237
#> 278         Gm11423       (B-A ∩ B-C)  2  237
#> 279           Clip1       (B-A ∩ B-C)  2  237
#> 280          Slc9a2       (B-A ∩ B-C)  2  237
#> 281          Lgalsl       (B-A ∩ B-C)  2  237
#> 282            Lgmn       (B-A ∩ B-C)  2  237
#> 283          Nynrin       (B-A ∩ B-C)  2  237
#> 284           Gpsm2       (B-A ∩ B-C)  2  237
#> 285         Aldh7a1       (B-A ∩ B-C)  2  237
#> 286           Rcan3       (B-A ∩ B-C)  2  237
#> 287         Bcl2l14       (B-A ∩ B-C)  2  237
#> 288           Mllt4       (B-A ∩ B-C)  2  237
#> 289          Nusap1       (B-A ∩ B-C)  2  237
#> 290            Aspm       (B-A ∩ B-C)  2  237
#> 291          Plxnd1       (B-A ∩ B-C)  2  237
#> 292          Dopey1       (B-A ∩ B-C)  2  237
#> 293           Sesn2       (B-A ∩ B-C)  2  237
#> 294          Kif13a       (B-A ∩ B-C)  2  237
#> 295          Rassf2       (B-A ∩ B-C)  2  237
#> 296         Slc22a3       (B-A ∩ B-C)  2  237
#> 297          B3gnt1       (B-A ∩ B-C)  2  237
#> 298      D3Ertd254e       (B-A ∩ B-C)  2  237
#> 299          Zfp868       (B-A ∩ B-C)  2  237
#> 300        Arhgef11       (B-A ∩ B-C)  2  237
#> 301          Frmd4b       (B-A ∩ B-C)  2  237
#> 302            Smtn       (B-A ∩ B-C)  2  237
#> 303           Gna12       (B-A ∩ B-C)  2  237
#> 304            Il15       (B-A ∩ B-C)  2  237
#> 305          Lmbrd1       (B-A ∩ B-C)  2  237
#> 306         Gm16299       (B-A ∩ B-C)  2  237
#> 307          Sapcd2       (B-A ∩ B-C)  2  237
#> 308          Vps13d       (B-A ∩ B-C)  2  237
#> 309          Atp8b2       (B-A ∩ B-C)  2  237
#> 310         Gm26656       (B-A ∩ B-C)  2  237
#> 311            Gpx7       (B-A ∩ B-C)  2  237
#> 312           Tube1       (B-A ∩ B-C)  2  237
#> 313           Il1r2       (B-A ∩ B-C)  2  237
#> 314           Snx33       (B-A ∩ B-C)  2  237
#> 315          Map3k5       (B-A ∩ B-C)  2  237
#> 316            Gch1       (B-A ∩ B-C)  2  237
#> 317           Sgpl1       (B-A ∩ B-C)  2  237
#> 318         Golph3l       (B-A ∩ B-C)  2  237
#> 319           Hipk2       (B-A ∩ B-C)  2  237
#> 320          Mboat1       (B-A ∩ B-C)  2  237
#> 321          Pik3cg       (B-A ∩ B-C)  2  237
#> 322            Bcor       (B-A ∩ B-C)  2  237
#> 323          Armcx5       (B-A ∩ B-C)  2  237
#> 324          Dusp22       (B-A ∩ B-C)  2  237
#> 325           Cdk20       (B-A ∩ B-C)  2  237
#> 326            Sdc4       (B-A ∩ B-C)  2  237
#> 327   D130017N08Rik       (B-A ∩ B-C)  2  237
#> 328           Cdca2       (B-A ∩ B-C)  2  237
#> 329         Kbtbd11       (B-A ∩ B-C)  2  237
#> 330          Sigirr       (B-A ∩ B-C)  2  237
#> 331         Cyp4f18       (B-A ∩ B-C)  2  237
#> 332          Incenp       (B-A ∩ B-C)  2  237
#> 333   A630033H20Rik       (B-A ∩ B-C)  2  237
#> 334          Zfp931       (B-A ∩ B-C)  2  237
#> 335           Terf1       (B-A ∩ B-C)  2  237
#> 336           Ces2g       (B-A ∩ B-C)  2  237
#> 337         Plekhg3       (B-A ∩ B-C)  2  237
#> 338           Mki67       (B-A ∩ B-C)  2  237
#> 339         Smarca2       (B-A ∩ B-C)  2  237
#> 340          Nhlrc1       (B-A ∩ B-C)  2  237
#> 341           Cers5       (B-A ∩ B-C)  2  237
#> 342           Smug1       (B-A ∩ B-C)  2  237
#> 343        Ankrd13d       (B-A ∩ B-C)  2  237
#> 344          Pmepa1       (B-A ∩ B-C)  2  237
#> 345         Sec61a2       (B-A ∩ B-C)  2  237
#> 346         Fam229b       (B-A ∩ B-C)  2  237
#> 347           Gins3       (B-A ∩ B-C)  2  237
#> 348             Lss       (B-A ∩ B-C)  2  237
#> 349          Srgap2       (B-A ∩ B-C)  2  237
#> 350           Tdrd7       (B-A ∩ B-C)  2  237
#> 351            Eri2       (B-A ∩ B-C)  2  237
#> 352           Wdr35       (B-A ∩ B-C)  2  237
#> 353          Map2k6       (B-A ∩ B-C)  2  237
#> 354   E230001N04Rik       (B-A ∩ B-C)  2  237
#> 355            Ispd       (B-A ∩ B-C)  2  237
#> 356           Cpne5       (B-A ∩ B-C)  2  237
#> 357           Krt83       (B-A ∩ B-C)  2  237
#> 358          Zfp398       (B-A ∩ B-C)  2  237
#> 359          Prss16       (B-A ∩ B-C)  2  237
#> 360           Matn2       (B-A ∩ B-C)  2  237
#> 361           Inpp1       (B-A ∩ B-C)  2  237
#> 362          Pex11a       (B-A ∩ B-C)  2  237
#> 363           Klhl2       (B-A ∩ B-C)  2  237
#> 364           Siah2       (B-A ∩ B-C)  2  237
#> 365          Pik3ca       (B-A ∩ B-C)  2  237
#> 366         Tmem243       (B-A ∩ B-C)  2  237
#> 367          Trim62       (B-A ∩ B-C)  2  237
#> 368           Pycr1       (B-A ∩ B-C)  2  237
#> 369            Rgcc       (B-A ∩ B-C)  2  237
#> 370           Itpr1       (B-A ∩ B-C)  2  237
#> 371            Epg5       (B-A ∩ B-C)  2  237
#> 372           Anxa4       (B-A ∩ B-C)  2  237
#> 373           P2rx4       (B-A ∩ B-C)  2  237
#> 374          Zfp449       (B-A ∩ B-C)  2  237
#> 375         Gm20667       (B-A ∩ B-C)  2  237
#> 376            Zhx3       (B-A ∩ B-C)  2  237
#> 377           Prkd2       (B-A ∩ B-C)  2  237
#> 378           Stox1       (B-A ∩ B-C)  2  237
#> 379          Map3k8       (B-A ∩ B-C)  2  237
#> 380          Zfp568       (B-A ∩ B-C)  2  237
#> 381          Zfp808       (B-A ∩ B-C)  2  237
#> 382            Als2       (B-A ∩ B-C)  2  237
#> 383           Naprt       (B-A ∩ B-C)  2  237
#> 384          Ankle1       (B-A ∩ B-C)  2  237
#> 385             Evl       (B-A ∩ B-C)  2  237
#> 386          Parp10       (B-A ∩ B-C)  2  237
#> 387             Tst       (B-A ∩ B-C)  2  237
#> 388            Mbd5       (B-A ∩ B-C)  2  237
#> 389         Tmem260       (B-A ∩ B-C)  2  237
#> 390         Slc33a1       (B-A ∩ B-C)  2  237
#> 391           Gstm1       (B-A ∩ B-C)  2  237
#> 392           Nim1k       (B-A ∩ B-C)  2  237
#> 393            Hcst       (B-A ∩ B-C)  2  237
#> 394          Med12l       (B-A ∩ B-C)  2  237
#> 395         Ppip5k1       (B-A ∩ B-C)  2  237
#> 396   9930104L06Rik       (B-A ∩ B-C)  2  237
#> 397            Layn       (B-A ∩ B-C)  2  237
#> 398          Glipr2       (B-A ∩ B-C)  2  237
#> 399             Gk5       (B-A ∩ B-C)  2  237
#> 400          Zfp788       (B-A ∩ B-C)  2  237
#> 401          Zbtb16       (B-A ∩ C-A)  3  518
#> 402          Tmsb4x       (B-A ∩ C-A)  3  518
#> 403           Tapt1       (B-A ∩ C-A)  3  518
#> 404           Itgb7       (B-A ∩ C-A)  3  518
#> 405            Bcl2       (B-A ∩ C-A)  3  518
#> 406            Cnn3       (B-A ∩ C-A)  3  518
#> 407          Il17ra       (B-A ∩ C-A)  3  518
#> 408            Ccl5       (B-A ∩ C-A)  3  518
#> 409            Ighm       (B-A ∩ C-A)  3  518
#> 410           Csrp1       (B-A ∩ C-A)  3  518
#> 411           Xrcc6       (B-A ∩ C-A)  3  518
#> 412   A630014C17Rik       (B-A ∩ C-A)  3  518
#> 413            Cnn2       (B-A ∩ C-A)  3  518
#> 414           Paqr5       (B-A ∩ C-A)  3  518
#> 415            Add3       (B-A ∩ C-A)  3  518
#> 416           Pgam2       (B-A ∩ C-A)  3  518
#> 417           Tuba8       (B-A ∩ C-A)  3  518
#> 418          Notch1       (B-A ∩ C-A)  3  518
#> 419          Fcer1g       (B-A ∩ C-A)  3  518
#> 420             Tox       (B-A ∩ C-A)  3  518
#> 421   9230105E05Rik       (B-A ∩ C-A)  3  518
#> 422          Lgals9       (B-A ∩ C-A)  3  518
#> 423            Tlr7       (B-A ∩ C-A)  3  518
#> 424         Siglech       (B-A ∩ C-A)  3  518
#> 425          Ptp4a3       (B-A ∩ C-A)  3  518
#> 426          Col9a3       (B-A ∩ C-A)  3  518
#> 427          Ttc39c       (B-A ∩ C-A)  3  518
#> 428            Dbn1       (B-A ∩ C-A)  3  518
#> 429            Bin2       (B-A ∩ C-A)  3  518
#> 430           Ephb2       (B-A ∩ C-A)  3  518
#> 431         Tnfrsf9       (B-A ∩ C-A)  3  518
#> 432          Papss2       (B-A ∩ C-A)  3  518
#> 433        BC021614       (B-A ∩ C-A)  3  518
#> 434            Ddx4       (B-A ∩ C-A)  3  518
#> 435           Ramp1       (B-A ∩ C-A)  3  518
#> 436          Tmem51       (B-A ∩ C-A)  3  518
#> 437           Eomes       (B-A ∩ C-A)  3  518
#> 438           Batf3       (B-A ∩ C-A)  3  518
#> 439          Gm5111       (B-A ∩ C-A)  3  518
#> 440            Polm       (B-A ∩ C-A)  3  518
#> 441          Megf10       (B-A ∩ C-A)  3  518
#> 442            Jak1       (B-A ∩ C-A)  3  518
#> 443          Gpr114       (B-A ∩ C-A)  3  518
#> 444          Fam71b       (B-A ∩ C-A)  3  518
#> 445           Sytl4       (B-A ∩ C-A)  3  518
#> 446           Fgf15       (B-A ∩ C-A)  3  518
#> 447           Gcnt2       (B-A ∩ C-A)  3  518
#> 448   8430408G22Rik       (B-A ∩ C-A)  3  518
#> 449             Fv1       (B-A ∩ C-A)  3  518
#> 450         Fam189b       (B-A ∩ C-A)  3  518
#> 451         Gorasp1       (B-A ∩ C-A)  3  518
#> 452            Fgf3       (B-A ∩ C-A)  3  518
#> 453        Ighv1-77       (B-A ∩ C-A)  3  518
#> 454            Cd52       (B-A ∩ C-A)  3  518
#> 455         Plekha7       (B-A ∩ C-A)  3  518
#> 456          Atp2a3       (B-A ∩ C-A)  3  518
#> 457        Serpini1       (B-A ∩ C-A)  3  518
#> 458         Gm16565       (B-A ∩ C-A)  3  518
#> 459           Egfl7       (B-A ∩ C-A)  3  518
#> 460           Rapsn       (B-A ∩ C-A)  3  518
#> 461            Tcf4       (B-A ∩ C-A)  3  518
#> 462         Tacstd2       (B-A ∩ C-A)  3  518
#> 463            Rgl1       (B-A ∩ C-A)  3  518
#> 464           Htra3       (B-A ∩ C-A)  3  518
#> 465          Smim14       (B-A ∩ C-A)  3  518
#> 466           Hnf4a       (B-A ∩ C-A)  3  518
#> 467           Gfod1       (B-A ∩ C-A)  3  518
#> 468            Cd33       (B-A ∩ C-A)  3  518
#> 469          Fam73b       (B-A ∩ C-A)  3  518
#> 470            Hes1       (B-A ∩ C-A)  3  518
#> 471         Clec12a       (B-A ∩ C-A)  3  518
#> 472           Ap3s1       (B-A ∩ C-A)  3  518
#> 473       D8Ertd82e       (B-A ∩ C-A)  3  518
#> 474        Epb4.1l3       (B-A ∩ C-A)  3  518
#> 475           Nol4l       (B-A ∩ C-A)  3  518
#> 476         Gm16710       (B-A ∩ C-A)  3  518
#> 477           Tnni1       (B-A ∩ C-A)  3  518
#> 478           Emid1       (B-A ∩ C-A)  3  518
#> 479           Gp1bb       (B-A ∩ C-A)  3  518
#> 480            Cd93       (B-A ∩ C-A)  3  518
#> 481            Tcf7       (B-A ∩ C-A)  3  518
#> 482            Fhl2       (B-A ∩ C-A)  3  518
#> 483            Myof       (B-A ∩ C-A)  3  518
#> 484            Myh9       (B-A ∩ C-A)  3  518
#> 485          Nccrp1       (B-A ∩ C-A)  3  518
#> 486            Nek6       (B-A ∩ C-A)  3  518
#> 487           Flt3l       (B-A ∩ C-A)  3  518
#> 488            Blnk       (B-A ∩ C-A)  3  518
#> 489        Rap1gap2       (B-A ∩ C-A)  3  518
#> 490            Mc5r       (B-A ∩ C-A)  3  518
#> 491            Cdc6       (B-A ∩ C-A)  3  518
#> 492            Cd96       (B-A ∩ C-A)  3  518
#> 493            Bin1       (B-A ∩ C-A)  3  518
#> 494            Ctsc       (B-A ∩ C-A)  3  518
#> 495          Gm5960       (B-A ∩ C-A)  3  518
#> 496        Tmem120b       (B-A ∩ C-A)  3  518
#> 497         Trbj1-2       (B-A ∩ C-A)  3  518
#> 498            Chn2       (B-A ∩ C-A)  3  518
#> 499           Cd160       (B-A ∩ C-A)  3  518
#> 500         Fam134b       (B-A ∩ C-A)  3  518
#> 501            Ttpa       (B-A ∩ C-A)  3  518
#> 502            Ugcg       (B-A ∩ C-A)  3  518
#> 503          Malat1       (B-A ∩ C-A)  3  518
#> 504         Tcrg-C2       (B-A ∩ C-A)  3  518
#> 505             Cd7       (B-A ∩ C-A)  3  518
#> 506        Trbv12-2       (B-A ∩ C-A)  3  518
#> 507            Lsm3       (B-A ∩ C-A)  3  518
#> 508            Fgl2       (B-A ∩ C-A)  3  518
#> 509            Fgd4       (B-A ∩ C-A)  3  518
#> 510   1700040L02Rik       (B-A ∩ C-A)  3  518
#> 511          Spint1       (B-A ∩ C-A)  3  518
#> 512            Rfc4       (B-A ∩ C-A)  3  518
#> 513           Azin2       (B-A ∩ C-A)  3  518
#> 514         Col19a1       (B-A ∩ C-A)  3  518
#> 515           Cxcr5       (B-A ∩ C-A)  3  518
#> 516           Runx2       (B-A ∩ C-A)  3  518
#> 517          Gm4208       (B-A ∩ C-A)  3  518
#> 518           Rnf43       (B-A ∩ C-A)  3  518
#> 519           Trbd1       (B-A ∩ C-A)  3  518
#> 520           Cdk19       (B-A ∩ C-A)  3  518
#> 521           Gpr68       (B-A ∩ C-A)  3  518
#> 522            Lmo2       (B-A ∩ C-A)  3  518
#> 523        Ighv1-82       (B-A ∩ C-A)  3  518
#> 524           Epha2       (B-A ∩ C-A)  3  518
#> 525        Cacna2d2       (B-A ∩ C-A)  3  518
#> 526         Tmem109       (B-A ∩ C-A)  3  518
#> 527           Itgb2       (B-A ∩ C-A)  3  518
#> 528           Prrc1       (B-A ∩ C-A)  3  518
#> 529            Rgs8       (B-A ∩ C-A)  3  518
#> 530           Smad3       (B-A ∩ C-A)  3  518
#> 531           Mecom       (B-A ∩ C-A)  3  518
#> 532            Tle1       (B-A ∩ C-A)  3  518
#> 533          Mif4gd       (B-A ∩ C-A)  3  518
#> 534   RP24-350F17.2       (B-A ∩ C-A)  3  518
#> 535           Il2rb       (B-A ∩ C-A)  3  518
#> 536            Nrgn       (B-A ∩ C-A)  3  518
#> 537          Acot11       (B-A ∩ C-A)  3  518
#> 538   2010300C02Rik       (B-A ∩ C-A)  3  518
#> 539            Tyw5       (B-A ∩ C-A)  3  518
#> 540            Uaca       (B-A ∩ C-A)  3  518
#> 541         Trbj1-5       (B-A ∩ C-A)  3  518
#> 542           Clip4       (B-A ∩ C-A)  3  518
#> 543         Slc35d3       (B-A ∩ C-A)  3  518
#> 544           Acap2       (B-A ∩ C-A)  3  518
#> 545           Mier3       (B-A ∩ C-A)  3  518
#> 546          Thnsl2       (B-A ∩ C-A)  3  518
#> 547          Serac1       (B-A ∩ C-A)  3  518
#> 548            Nkg7       (B-A ∩ C-A)  3  518
#> 549          Fbxo17       (B-A ∩ C-A)  3  518
#> 550            Rec8       (B-A ∩ C-A)  3  518
#> 551          Cryba4       (B-A ∩ C-A)  3  518
#> 552           Lims1       (B-A ∩ C-A)  3  518
#> 553          Mfhas1       (B-A ∩ C-A)  3  518
#> 554          Pcp4l1       (B-A ∩ C-A)  3  518
#> 555            Cd38       (B-A ∩ C-A)  3  518
#> 556            Tle3       (B-A ∩ C-A)  3  518
#> 557            Lrp1       (B-A ∩ C-A)  3  518
#> 558            Chka       (B-A ∩ C-A)  3  518
#> 559            St14       (B-A ∩ C-A)  3  518
#> 560           Robo4       (B-A ∩ C-A)  3  518
#> 561           Hvcn1       (B-A ∩ C-A)  3  518
#> 562         Map3k19       (B-A ∩ C-A)  3  518
#> 563          Ifitm3       (B-A ∩ C-A)  3  518
#> 564           Vegfa       (B-A ∩ C-A)  3  518
#> 565          Ndufa4       (B-A ∩ C-A)  3  518
#> 566          Srgap3       (B-A ∩ C-A)  3  518
#> 567          Tspyl4       (B-A ∩ C-A)  3  518
#> 568           Lipt2       (B-A ∩ C-A)  3  518
#> 569            Tymp       (B-A ∩ C-A)  3  518
#> 570           Hif1a       (B-A ∩ C-A)  3  518
#> 571            Gse1       (B-A ∩ C-A)  3  518
#> 572           Pvrl4       (B-A ∩ C-A)  3  518
#> 573           Spag6       (B-A ∩ C-A)  3  518
#> 574          Man2a1       (B-A ∩ C-A)  3  518
#> 575            Dtx4       (B-A ∩ C-A)  3  518
#> 576           Egln1       (B-A ∩ C-A)  3  518
#> 577          Zfp526       (B-A ∩ C-A)  3  518
#> 578           Dirc2       (B-A ∩ C-A)  3  518
#> 579          Lonrf1       (B-A ∩ C-A)  3  518
#> 580          Rhbdl3       (B-A ∩ C-A)  3  518
#> 581             Erg       (B-A ∩ C-A)  3  518
#> 582          Il10rb       (B-A ∩ C-A)  3  518
#> 583         Tsc22d1       (B-A ∩ C-A)  3  518
#> 584           Dcaf7       (B-A ∩ C-A)  3  518
#> 585           Lmnb1       (B-A ∩ C-A)  3  518
#> 586   2900076A07Rik       (B-A ∩ C-A)  3  518
#> 587            Eya1       (B-A ∩ C-A)  3  518
#> 588            Gab1       (B-A ∩ C-A)  3  518
#> 589        Trp53i11       (B-A ∩ C-A)  3  518
#> 590           Phtf2       (B-A ∩ C-A)  3  518
#> 591           Susd1       (B-A ∩ C-A)  3  518
#> 592            Rag1       (B-A ∩ C-A)  3  518
#> 593           Rnf14       (B-A ∩ C-A)  3  518
#> 594          P2ry13       (B-A ∩ C-A)  3  518
#> 595         Rasgrp3       (B-A ∩ C-A)  3  518
#> 596           Rbpms       (B-A ∩ C-A)  3  518
#> 597          Ankib1       (B-A ∩ C-A)  3  518
#> 598            Ubl4       (B-A ∩ C-A)  3  518
#> 599            Tab2       (B-A ∩ C-A)  3  518
#> 600          Plxdc1       (B-A ∩ C-A)  3  518
#> 601         Plekhf1       (B-A ∩ C-A)  3  518
#> 602         Cd300lb       (B-A ∩ C-A)  3  518
#> 603          Entpd7       (B-A ∩ C-A)  3  518
#> 604           Tpst2       (B-A ∩ C-A)  3  518
#> 605             Xdh       (B-A ∩ C-A)  3  518
#> 606           Plcg1       (B-A ∩ C-A)  3  518
#> 607          Kctd11       (B-A ∩ C-A)  3  518
#> 608           Itgax       (B-A ∩ C-A)  3  518
#> 609           Arap3       (B-A ∩ C-A)  3  518
#> 610            Cdk1       (B-A ∩ C-A)  3  518
#> 611         Slc37a2       (B-A ∩ C-A)  3  518
#> 612            Exo1       (B-A ∩ C-A)  3  518
#> 613         Gm26597       (B-A ∩ C-A)  3  518
#> 614         Trbj1-3       (B-A ∩ C-A)  3  518
#> 615            Pygm       (B-A ∩ C-A)  3  518
#> 616         Hdgfrp3       (B-A ∩ C-A)  3  518
#> 617         Trbj1-4       (B-A ∩ C-A)  3  518
#> 618          Ptger3       (B-A ∩ C-A)  3  518
#> 619           Gfra1       (B-A ∩ C-A)  3  518
#> 620          Tcf7l2       (B-A ∩ C-A)  3  518
#> 621           Ndrg2       (B-A ∩ C-A)  3  518
#> 622            Rxra       (B-A ∩ C-A)  3  518
#> 623          Nsmce1       (B-A ∩ C-A)  3  518
#> 624            Sdc2       (B-A ∩ C-A)  3  518
#> 625           Stat3       (B-A ∩ C-A)  3  518
#> 626          Ifitm6       (B-A ∩ C-A)  3  518
#> 627            Lyz2       (B-A ∩ C-A)  3  518
#> 628           Cdyl2       (B-A ∩ C-A)  3  518
#> 629            Fan1       (B-A ∩ C-A)  3  518
#> 630          Cldn25       (B-A ∩ C-A)  3  518
#> 631           Glis2       (B-A ∩ C-A)  3  518
#> 632           Snx14       (B-A ∩ C-A)  3  518
#> 633         Il12rb2       (B-A ∩ C-A)  3  518
#> 634           Camkv       (B-A ∩ C-A)  3  518
#> 635          Mcoln3       (B-A ∩ C-A)  3  518
#> 636            Pim2       (B-A ∩ C-A)  3  518
#> 637           Sepp1       (B-A ∩ C-A)  3  518
#> 638          Arpp21       (B-A ∩ C-A)  3  518
#> 639            Trac       (B-A ∩ C-A)  3  518
#> 640             Nlk       (B-A ∩ C-A)  3  518
#> 641          Jmjd1c       (B-A ∩ C-A)  3  518
#> 642           Rab37       (B-A ∩ C-A)  3  518
#> 643          Chd3os       (B-A ∩ C-A)  3  518
#> 644          Gas2l3       (B-A ∩ C-A)  3  518
#> 645        Rab3gap1       (B-A ∩ C-A)  3  518
#> 646           Thtpa       (B-A ∩ C-A)  3  518
#> 647   9130221H12Rik       (B-A ∩ C-A)  3  518
#> 648           Srxn1       (B-A ∩ C-A)  3  518
#> 649         Herpud1       (B-A ∩ C-A)  3  518
#> 650        Ccdc102a       (B-A ∩ C-A)  3  518
#> 651            Cd3g       (B-A ∩ C-A)  3  518
#> 652         Fam212a       (B-A ∩ C-A)  3  518
#> 653          Sema4a       (B-A ∩ C-A)  3  518
#> 654        AA467197       (B-A ∩ C-A)  3  518
#> 655         Slc18a1       (B-A ∩ C-A)  3  518
#> 656          Gpr155       (B-A ∩ C-A)  3  518
#> 657   4930579K19Rik       (B-A ∩ C-A)  3  518
#> 658          Phlpp1       (B-A ∩ C-A)  3  518
#> 659        Bcas1os2       (B-A ∩ C-A)  3  518
#> 660           Ddit4       (B-A ∩ C-A)  3  518
#> 661            Styx       (B-A ∩ C-A)  3  518
#> 662          Trim59       (B-A ∩ C-A)  3  518
#> 663            Mafk       (B-A ∩ C-A)  3  518
#> 664            Ctsw       (B-A ∩ C-A)  3  518
#> 665            Ctc1       (B-A ∩ C-A)  3  518
#> 666          Tmem8b       (B-A ∩ C-A)  3  518
#> 667           Sez6l       (B-A ∩ C-A)  3  518
#> 668         Tmem87a       (B-A ∩ C-A)  3  518
#> 669            Fggy       (B-A ∩ C-A)  3  518
#> 670   9430091E24Rik       (B-A ∩ C-A)  3  518
#> 671          Hivep1       (B-A ∩ C-A)  3  518
#> 672           Ttc13       (B-A ∩ C-A)  3  518
#> 673           S1pr2       (B-A ∩ C-A)  3  518
#> 674             Pbk       (B-A ∩ C-A)  3  518
#> 675         Trbj1-6       (B-A ∩ C-A)  3  518
#> 676          Fam69b       (B-A ∩ C-A)  3  518
#> 677          Tmem43       (B-A ∩ C-A)  3  518
#> 678         Gm11944       (B-A ∩ C-A)  3  518
#> 679   1110008P14Rik       (B-A ∩ C-A)  3  518
#> 680           Snx30       (B-A ∩ C-A)  3  518
#> 681            Cds1       (B-A ∩ C-A)  3  518
#> 682            Uap1       (B-A ∩ C-A)  3  518
#> 683   9330136K24Rik       (B-A ∩ C-A)  3  518
#> 684           Efcc1       (B-A ∩ C-A)  3  518
#> 685           Tram1       (B-A ∩ C-A)  3  518
#> 686         Ankrd44       (B-A ∩ C-A)  3  518
#> 687          Adssl1       (B-A ∩ C-A)  3  518
#> 688   2610015P09Rik       (B-A ∩ C-A)  3  518
#> 689            Rinl       (B-A ∩ C-A)  3  518
#> 690           Fmnl2       (B-A ∩ C-A)  3  518
#> 691          Entpd1       (B-A ∩ C-A)  3  518
#> 692           Iffo2       (B-A ∩ C-A)  3  518
#> 693            Gphn       (B-A ∩ C-A)  3  518
#> 694           H2-Ob       (B-A ∩ C-A)  3  518
#> 695          Ms4a4c       (B-A ∩ C-A)  3  518
#> 696         Ppapdc2       (B-A ∩ C-A)  3  518
#> 697            Tha1       (B-A ∩ C-A)  3  518
#> 698          B3gnt5       (B-A ∩ C-A)  3  518
#> 699            Sdc3       (B-A ∩ C-A)  3  518
#> 700           Ifi44       (B-A ∩ C-A)  3  518
#> 701            Btg2       (B-A ∩ C-A)  3  518
#> 702            Naaa       (B-A ∩ C-A)  3  518
#> 703          Lgals3       (B-A ∩ C-A)  3  518
#> 704   C330011M18Rik       (B-A ∩ C-A)  3  518
#> 705         Aldh4a1       (B-A ∩ C-A)  3  518
#> 706          Nudt19       (B-A ∩ C-A)  3  518
#> 707           Cd302       (B-A ∩ C-A)  3  518
#> 708           Asap1       (B-A ∩ C-A)  3  518
#> 709        Arhgef10       (B-A ∩ C-A)  3  518
#> 710           Mllt3       (B-A ∩ C-A)  3  518
#> 711           Kif23       (B-A ∩ C-A)  3  518
#> 712          Gm1966       (B-A ∩ C-A)  3  518
#> 713           Slfn2       (B-A ∩ C-A)  3  518
#> 714          Hspbp1       (B-A ∩ C-A)  3  518
#> 715          Klrb1f       (B-A ∩ C-A)  3  518
#> 716           Fads3       (B-A ∩ C-A)  3  518
#> 717            Amd1       (B-A ∩ C-A)  3  518
#> 718         Trbj1-7       (B-A ∩ C-A)  3  518
#> 719         Slc41a3       (B-A ∩ C-A)  3  518
#> 720          Osbpl3       (B-A ∩ C-A)  3  518
#> 721            Bbc3       (B-A ∩ C-A)  3  518
#> 722        Ighv1-23       (B-A ∩ C-A)  3  518
#> 723            Apoe       (B-A ∩ C-A)  3  518
#> 724          Sept11       (B-A ∩ C-A)  3  518
#> 725         Gm10612       (B-A ∩ C-A)  3  518
#> 726           Nucb2       (B-A ∩ C-A)  3  518
#> 727            Elf4       (B-A ∩ C-A)  3  518
#> 728          Shisa8       (B-A ∩ C-A)  3  518
#> 729           Ggta1       (B-A ∩ C-A)  3  518
#> 730           Frmd6       (B-A ∩ C-A)  3  518
#> 731           Trim8       (B-A ∩ C-A)  3  518
#> 732          Havcr2       (B-A ∩ C-A)  3  518
#> 733            Klf7       (B-A ∩ C-A)  3  518
#> 734         Slc29a1       (B-A ∩ C-A)  3  518
#> 735            Tpk1       (B-A ∩ C-A)  3  518
#> 736            Rgs3       (B-A ∩ C-A)  3  518
#> 737           Vldlr       (B-A ∩ C-A)  3  518
#> 738           Il6st       (B-A ∩ C-A)  3  518
#> 739           Prss2       (B-A ∩ C-A)  3  518
#> 740          Gnpda1       (B-A ∩ C-A)  3  518
#> 741           Gpsm3       (B-A ∩ C-A)  3  518
#> 742          Cacnb2       (B-A ∩ C-A)  3  518
#> 743            Vasp       (B-A ∩ C-A)  3  518
#> 744            Fut8       (B-A ∩ C-A)  3  518
#> 745            Glrx       (B-A ∩ C-A)  3  518
#> 746           Mrvi1       (B-A ∩ C-A)  3  518
#> 747            Rrs1       (B-A ∩ C-A)  3  518
#> 748           Ttbk2       (B-A ∩ C-A)  3  518
#> 749           Kdm4a       (B-A ∩ C-A)  3  518
#> 750         Gm15590       (B-A ∩ C-A)  3  518
#> 751           Acss2       (B-A ∩ C-A)  3  518
#> 752         Creb3l2       (B-A ∩ C-A)  3  518
#> 753            Cd48       (B-A ∩ C-A)  3  518
#> 754            Pdpr       (B-A ∩ C-A)  3  518
#> 755            Rhof       (B-A ∩ C-A)  3  518
#> 756          Osbpl5       (B-A ∩ C-A)  3  518
#> 757           Llgl2       (B-A ∩ C-A)  3  518
#> 758          Cuedc1       (B-A ∩ C-A)  3  518
#> 759            Il18       (B-A ∩ C-A)  3  518
#> 760         Aldh1b1       (B-A ∩ C-A)  3  518
#> 761           Sort1       (B-A ∩ C-A)  3  518
#> 762          Lingo4       (B-A ∩ C-A)  3  518
#> 763           Klrk1       (B-A ∩ C-A)  3  518
#> 764         Clec10a       (B-A ∩ C-A)  3  518
#> 765          Ms4a4b       (B-A ∩ C-A)  3  518
#> 766         Fam179a       (B-A ∩ C-A)  3  518
#> 767         Fam129b       (B-A ∩ C-A)  3  518
#> 768            Zfp2       (B-A ∩ C-A)  3  518
#> 769         Slc18a2       (B-A ∩ C-A)  3  518
#> 770            Net1       (B-A ∩ C-A)  3  518
#> 771         Ighd4-1       (B-A ∩ C-A)  3  518
#> 772           Kif15       (B-A ∩ C-A)  3  518
#> 773         Dennd5a       (B-A ∩ C-A)  3  518
#> 774          Shisa2       (B-A ∩ C-A)  3  518
#> 775          Tbxa2r       (B-A ∩ C-A)  3  518
#> 776          Golga3       (B-A ∩ C-A)  3  518
#> 777             C8g       (B-A ∩ C-A)  3  518
#> 778   4632428N05Rik       (B-A ∩ C-A)  3  518
#> 779           Prune       (B-A ∩ C-A)  3  518
#> 780            Btla       (B-A ∩ C-A)  3  518
#> 781            Sgce       (B-A ∩ C-A)  3  518
#> 782       Hist3h2ba       (B-A ∩ C-A)  3  518
#> 783            Ect2       (B-A ∩ C-A)  3  518
#> 784          Mrps10       (B-A ∩ C-A)  3  518
#> 785          Celsr1       (B-A ∩ C-A)  3  518
#> 786          Ift172       (B-A ∩ C-A)  3  518
#> 787           Whamm       (B-A ∩ C-A)  3  518
#> 788          Chst15       (B-A ∩ C-A)  3  518
#> 789          Gm8822       (B-A ∩ C-A)  3  518
#> 790         C1qtnf6       (B-A ∩ C-A)  3  518
#> 791         Pla2g15       (B-A ∩ C-A)  3  518
#> 792          Pkmyt1       (B-A ∩ C-A)  3  518
#> 793           Fcgr3       (B-A ∩ C-A)  3  518
#> 794          Fam64a       (B-A ∩ C-A)  3  518
#> 795          Knstrn       (B-A ∩ C-A)  3  518
#> 796          Amica1       (B-A ∩ C-A)  3  518
#> 797           Gtse1       (B-A ∩ C-A)  3  518
#> 798            B9d2       (B-A ∩ C-A)  3  518
#> 799       Gabarapl1       (B-A ∩ C-A)  3  518
#> 800           Ly6c2       (B-A ∩ C-A)  3  518
#> 801          Rab27b       (B-A ∩ C-A)  3  518
#> 802           Klri2       (B-A ∩ C-A)  3  518
#> 803           Prkch       (B-A ∩ C-A)  3  518
#> 804         Mterf1b       (B-A ∩ C-A)  3  518
#> 805           Ccnd1       (B-A ∩ C-A)  3  518
#> 806         Ugt1a7c       (B-A ∩ C-A)  3  518
#> 807           Capn1       (B-A ∩ C-A)  3  518
#> 808          Zfp945       (B-A ∩ C-A)  3  518
#> 809           Chil3       (B-A ∩ C-A)  3  518
#> 810            Pfkp       (B-A ∩ C-A)  3  518
#> 811         Rundc3a       (B-A ∩ C-A)  3  518
#> 812           Adam2       (B-A ∩ C-A)  3  518
#> 813             Hk3       (B-A ∩ C-A)  3  518
#> 814            Pck2       (B-A ∩ C-A)  3  518
#> 815            Jak2       (B-A ∩ C-A)  3  518
#> 816           Uhrf1       (B-A ∩ C-A)  3  518
#> 817           Abhd5       (B-A ∩ C-A)  3  518
#> 818            Lrmp       (B-A ∩ C-A)  3  518
#> 819           Dock5       (B-A ∩ C-A)  3  518
#> 820         Tcp11l2       (B-A ∩ C-A)  3  518
#> 821   5430435G22Rik       (B-A ∩ C-A)  3  518
#> 822         Pla2g4a       (B-A ∩ C-A)  3  518
#> 823          Phlda3       (B-A ∩ C-A)  3  518
#> 824           Zufsp       (B-A ∩ C-A)  3  518
#> 825            Cd69       (B-A ∩ C-A)  3  518
#> 826   0610009L18Rik       (B-A ∩ C-A)  3  518
#> 827          Exoc3l       (B-A ∩ C-A)  3  518
#> 828         Ankrd32       (B-A ∩ C-A)  3  518
#> 829        Slc25a24       (B-A ∩ C-A)  3  518
#> 830            Mafg       (B-A ∩ C-A)  3  518
#> 831          H2-Eb1       (B-A ∩ C-A)  3  518
#> 832   RP23-330G24.3       (B-A ∩ C-A)  3  518
#> 833           Kank2       (B-A ∩ C-A)  3  518
#> 834           Rgs14       (B-A ∩ C-A)  3  518
#> 835            Ccnf       (B-A ∩ C-A)  3  518
#> 836             Mn1       (B-A ∩ C-A)  3  518
#> 837          Ccdc30       (B-A ∩ C-A)  3  518
#> 838            Art4       (B-A ∩ C-A)  3  518
#> 839         Ankrd37       (B-A ∩ C-A)  3  518
#> 840        Tnfrsf22       (B-A ∩ C-A)  3  518
#> 841         Zdhhc14       (B-A ∩ C-A)  3  518
#> 842           Tsacc       (B-A ∩ C-A)  3  518
#> 843            Cptp       (B-A ∩ C-A)  3  518
#> 844          Unc119       (B-A ∩ C-A)  3  518
#> 845          Hmbox1       (B-A ∩ C-A)  3  518
#> 846   A730017L22Rik       (B-A ∩ C-A)  3  518
#> 847           Kif3c       (B-A ∩ C-A)  3  518
#> 848   1110034G24Rik       (B-A ∩ C-A)  3  518
#> 849         Wbscr27       (B-A ∩ C-A)  3  518
#> 850           Fnip1       (B-A ∩ C-A)  3  518
#> 851          Zfp628       (B-A ∩ C-A)  3  518
#> 852           Arl4c       (B-A ∩ C-A)  3  518
#> 853          Vps37a       (B-A ∩ C-A)  3  518
#> 854           Pskh1       (B-A ∩ C-A)  3  518
#> 855           Ifi30       (B-A ∩ C-A)  3  518
#> 856            Nsmf       (B-A ∩ C-A)  3  518
#> 857         Tmem38a       (B-A ∩ C-A)  3  518
#> 858           Bub1b       (B-A ∩ C-A)  3  518
#> 859          Fam43a       (B-A ∩ C-A)  3  518
#> 860            Emc9       (B-A ∩ C-A)  3  518
#> 861          Dusp11       (B-A ∩ C-A)  3  518
#> 862          Gimap5       (B-A ∩ C-A)  3  518
#> 863           Rufy1       (B-A ∩ C-A)  3  518
#> 864           Simc1       (B-A ∩ C-A)  3  518
#> 865            Nab2       (B-A ∩ C-A)  3  518
#> 866           Cmss1       (B-A ∩ C-A)  3  518
#> 867          Gm5526       (B-A ∩ C-A)  3  518
#> 868         Zscan29       (B-A ∩ C-A)  3  518
#> 869            Fbp1       (B-A ∩ C-A)  3  518
#> 870          Rasal1       (B-A ∩ C-A)  3  518
#> 871          Hspb11       (B-A ∩ C-A)  3  518
#> 872   5031414D18Rik       (B-A ∩ C-A)  3  518
#> 873          Gm4890       (B-A ∩ C-A)  3  518
#> 874           Atad5       (B-A ∩ C-A)  3  518
#> 875            Cd28       (B-A ∩ C-A)  3  518
#> 876          Fbrsl1       (B-A ∩ C-A)  3  518
#> 877           Rdh13       (B-A ∩ C-A)  3  518
#> 878         Elmsan1       (B-A ∩ C-A)  3  518
#> 879           Brip1       (B-A ∩ C-A)  3  518
#> 880          Sptlc1       (B-A ∩ C-A)  3  518
#> 881           Cmtr1       (B-A ∩ C-A)  3  518
#> 882           Brca2       (B-A ∩ C-A)  3  518
#> 883           Zfp64       (B-A ∩ C-A)  3  518
#> 884            Aatk       (B-A ∩ C-A)  3  518
#> 885          Rnf141       (B-A ∩ C-A)  3  518
#> 886            Etv3       (B-A ∩ C-A)  3  518
#> 887           Cox17       (B-A ∩ C-A)  3  518
#> 888          Zmynd8       (B-A ∩ C-A)  3  518
#> 889          Gm9825       (B-A ∩ C-A)  3  518
#> 890           Irak3       (B-A ∩ C-A)  3  518
#> 891           Zfp40       (B-A ∩ C-A)  3  518
#> 892          Tmbim1       (B-A ∩ C-A)  3  518
#> 893          Trim16       (B-A ∩ C-A)  3  518
#> 894        Slc16a13       (B-A ∩ C-A)  3  518
#> 895            Ccl9       (B-A ∩ C-A)  3  518
#> 896           Dstyk       (B-A ∩ C-A)  3  518
#> 897            Ell2       (B-A ∩ C-A)  3  518
#> 898          Zfp608       (B-A ∩ C-A)  3  518
#> 899            Ass1       (B-A ∩ C-A)  3  518
#> 900            Rrad       (B-A ∩ C-A)  3  518
#> 901           Wdr47       (B-A ∩ C-A)  3  518
#> 902          Sorbs3       (B-A ∩ C-A)  3  518
#> 903           Mfsd4       (B-A ∩ C-A)  3  518
#> 904          Nmral1       (B-A ∩ C-A)  3  518
#> 905          Gm5637       (B-A ∩ C-A)  3  518
#> 906          Zfp763       (B-A ∩ C-A)  3  518
#> 907           Trit1       (B-A ∩ C-A)  3  518
#> 908          Kif21b       (B-A ∩ C-A)  3  518
#> 909            Cnr2       (B-A ∩ C-A)  3  518
#> 910           N4bp3       (B-A ∩ C-A)  3  518
#> 911           Lynx1       (B-A ∩ C-A)  3  518
#> 912             Bcr       (B-A ∩ C-A)  3  518
#> 913            Mdc1       (B-A ∩ C-A)  3  518
#> 914      St6galnac6       (B-A ∩ C-A)  3  518
#> 915           Herc3       (B-A ∩ C-A)  3  518
#> 916          Polr3f       (B-A ∩ C-A)  3  518
#> 917          Acadvl       (B-A ∩ C-A)  3  518
#> 918          Zfp947       (B-A ∩ C-A)  3  518
#> 919           Il2ra       (B-C ∩ C-A)  4  780
#> 920            Cd82       (B-C ∩ C-A)  4  780
#> 921           H2afy       (B-C ∩ C-A)  4  780
#> 922          Ifngr1       (B-C ∩ C-A)  4  780
#> 923        Slc25a12       (B-C ∩ C-A)  4  780
#> 924          Il1rl1       (B-C ∩ C-A)  4  780
#> 925           Il2rg       (B-C ∩ C-A)  4  780
#> 926            Map4       (B-C ∩ C-A)  4  780
#> 927           Tktl1       (B-C ∩ C-A)  4  780
#> 928           Spsb1       (B-C ∩ C-A)  4  780
#> 929          Fam20a       (B-C ∩ C-A)  4  780
#> 930           Camk4       (B-C ∩ C-A)  4  780
#> 931         Ccdc184       (B-C ∩ C-A)  4  780
#> 932            Sell       (B-C ∩ C-A)  4  780
#> 933            Tnik       (B-C ∩ C-A)  4  780
#> 934           Pcsk1       (B-C ∩ C-A)  4  780
#> 935            Irf8       (B-C ∩ C-A)  4  780
#> 936           Plac8       (B-C ∩ C-A)  4  780
#> 937            Ctr9       (B-C ∩ C-A)  4  780
#> 938          Shisa5       (B-C ∩ C-A)  4  780
#> 939         Hnrnpa1       (B-C ∩ C-A)  4  780
#> 940          Coro2b       (B-C ∩ C-A)  4  780
#> 941          Dnajc3       (B-C ∩ C-A)  4  780
#> 942           Rftn1       (B-C ∩ C-A)  4  780
#> 943           Napsa       (B-C ∩ C-A)  4  780
#> 944          Cercam       (B-C ∩ C-A)  4  780
#> 945            Bmp7       (B-C ∩ C-A)  4  780
#> 946           Gm973       (B-C ∩ C-A)  4  780
#> 947            Lat2       (B-C ∩ C-A)  4  780
#> 948          Il17rb       (B-C ∩ C-A)  4  780
#> 949            Bex6       (B-C ∩ C-A)  4  780
#> 950            Cap2       (B-C ∩ C-A)  4  780
#> 951           Ptgir       (B-C ∩ C-A)  4  780
#> 952           Ptpn6       (B-C ∩ C-A)  4  780
#> 953            Nt5e       (B-C ∩ C-A)  4  780
#> 954           Pydc3       (B-C ∩ C-A)  4  780
#> 955           Nmur1       (B-C ∩ C-A)  4  780
#> 956           Gpr97       (B-C ∩ C-A)  4  780
#> 957           Prkcq       (B-C ∩ C-A)  4  780
#> 958        Marcksl1       (B-C ∩ C-A)  4  780
#> 959          Galnt1       (B-C ∩ C-A)  4  780
#> 960          Stk17b       (B-C ∩ C-A)  4  780
#> 961         Cysltr1       (B-C ∩ C-A)  4  780
#> 962           Gpr56       (B-C ∩ C-A)  4  780
#> 963         Plekhb2       (B-C ∩ C-A)  4  780
#> 964        BC035044       (B-C ∩ C-A)  4  780
#> 965          Adam12       (B-C ∩ C-A)  4  780
#> 966           Itgal       (B-C ∩ C-A)  4  780
#> 967            Xcr1       (B-C ∩ C-A)  4  780
#> 968          Zfp184       (B-C ∩ C-A)  4  780
#> 969           Anks6       (B-C ∩ C-A)  4  780
#> 970          Rasal3       (B-C ∩ C-A)  4  780
#> 971         B3galt2       (B-C ∩ C-A)  4  780
#> 972           Xlr4a       (B-C ∩ C-A)  4  780
#> 973           Nrros       (B-C ∩ C-A)  4  780
#> 974            Ncf1       (B-C ∩ C-A)  4  780
#> 975            Idh2       (B-C ∩ C-A)  4  780
#> 976           Gpr65       (B-C ∩ C-A)  4  780
#> 977         Themis2       (B-C ∩ C-A)  4  780
#> 978          Ppfia4       (B-C ∩ C-A)  4  780
#> 979          Gpr146       (B-C ∩ C-A)  4  780
#> 980            Mtbp       (B-C ∩ C-A)  4  780
#> 981          Fbxo27       (B-C ∩ C-A)  4  780
#> 982        Serpinf1       (B-C ∩ C-A)  4  780
#> 983            Gpx1       (B-C ∩ C-A)  4  780
#> 984            Nod1       (B-C ∩ C-A)  4  780
#> 985             Ltb       (B-C ∩ C-A)  4  780
#> 986      AC161246.1       (B-C ∩ C-A)  4  780
#> 987            Fli1       (B-C ∩ C-A)  4  780
#> 988           Sept1       (B-C ∩ C-A)  4  780
#> 989            Snx2       (B-C ∩ C-A)  4  780
#> 990         Tnfsf11       (B-C ∩ C-A)  4  780
#> 991            Cerk       (B-C ∩ C-A)  4  780
#> 992         Hnrnpab       (B-C ∩ C-A)  4  780
#> 993          Prss57       (B-C ∩ C-A)  4  780
#> 994           Esyt2       (B-C ∩ C-A)  4  780
#> 995          Ccdc92       (B-C ∩ C-A)  4  780
#> 996         Smpdl3a       (B-C ∩ C-A)  4  780
#> 997           Itgae       (B-C ∩ C-A)  4  780
#> 998           Rab19       (B-C ∩ C-A)  4  780
#> 999          Cxcl10       (B-C ∩ C-A)  4  780
#> 1000          Lppr3       (B-C ∩ C-A)  4  780
#> 1001         Slc4a8       (B-C ∩ C-A)  4  780
#> 1002        Gm27252       (B-C ∩ C-A)  4  780
#> 1003          Csf3r       (B-C ∩ C-A)  4  780
#> 1004           Upp1       (B-C ∩ C-A)  4  780
#> 1005        Slco4a1       (B-C ∩ C-A)  4  780
#> 1006      Epb4.1l4b       (B-C ∩ C-A)  4  780
#> 1007           Gclc       (B-C ∩ C-A)  4  780
#> 1008         Ndfip1       (B-C ∩ C-A)  4  780
#> 1009  4930506M07Rik       (B-C ∩ C-A)  4  780
#> 1010        Fam102a       (B-C ∩ C-A)  4  780
#> 1011          Smim5       (B-C ∩ C-A)  4  780
#> 1012           Btg1       (B-C ∩ C-A)  4  780
#> 1013          Xlr4c       (B-C ∩ C-A)  4  780
#> 1014        Gm16271       (B-C ∩ C-A)  4  780
#> 1015  9030619P08Rik       (B-C ∩ C-A)  4  780
#> 1016           Rbl2       (B-C ∩ C-A)  4  780
#> 1017          Aldoc       (B-C ∩ C-A)  4  780
#> 1018         Slc9a9       (B-C ∩ C-A)  4  780
#> 1019          Lpar6       (B-C ∩ C-A)  4  780
#> 1020      Serpina3g       (B-C ∩ C-A)  4  780
#> 1021           Lcp2       (B-C ∩ C-A)  4  780
#> 1022         Il27ra       (B-C ∩ C-A)  4  780
#> 1023          Skap1       (B-C ∩ C-A)  4  780
#> 1024          Cmtm6       (B-C ∩ C-A)  4  780
#> 1025          Grhl2       (B-C ∩ C-A)  4  780
#> 1026         Pdgfrb       (B-C ∩ C-A)  4  780
#> 1027  1700113H08Rik       (B-C ∩ C-A)  4  780
#> 1028            Bid       (B-C ∩ C-A)  4  780
#> 1029          Clic4       (B-C ∩ C-A)  4  780
#> 1030           Sox4       (B-C ∩ C-A)  4  780
#> 1031           Prr5       (B-C ∩ C-A)  4  780
#> 1032        Aldh3a1       (B-C ∩ C-A)  4  780
#> 1033       AU040320       (B-C ∩ C-A)  4  780
#> 1034          Skap2       (B-C ∩ C-A)  4  780
#> 1035            Fah       (B-C ∩ C-A)  4  780
#> 1036          S1pr1       (B-C ∩ C-A)  4  780
#> 1037       Rap1gds1       (B-C ∩ C-A)  4  780
#> 1038          Trdv4       (B-C ∩ C-A)  4  780
#> 1039          Nup85       (B-C ∩ C-A)  4  780
#> 1040          Sirpa       (B-C ∩ C-A)  4  780
#> 1041            Maf       (B-C ∩ C-A)  4  780
#> 1042          Myo1e       (B-C ∩ C-A)  4  780
#> 1043         Steap3       (B-C ∩ C-A)  4  780
#> 1044          Zbtb4       (B-C ∩ C-A)  4  780
#> 1045        Tmem154       (B-C ∩ C-A)  4  780
#> 1046         Sh2d2a       (B-C ∩ C-A)  4  780
#> 1047  4930430F08Rik       (B-C ∩ C-A)  4  780
#> 1048           Cd84       (B-C ∩ C-A)  4  780
#> 1049  C030017K20Rik       (B-C ∩ C-A)  4  780
#> 1050         Ormdl3       (B-C ∩ C-A)  4  780
#> 1051        Unc93b1       (B-C ∩ C-A)  4  780
#> 1052          Adap1       (B-C ∩ C-A)  4  780
#> 1053         Themis       (B-C ∩ C-A)  4  780
#> 1054          Muc13       (B-C ∩ C-A)  4  780
#> 1055          Plcg2       (B-C ∩ C-A)  4  780
#> 1056            Vcl       (B-C ∩ C-A)  4  780
#> 1057           Elp2       (B-C ∩ C-A)  4  780
#> 1058       Cyb561a3       (B-C ∩ C-A)  4  780
#> 1059         Galnt3       (B-C ∩ C-A)  4  780
#> 1060            Ak4       (B-C ∩ C-A)  4  780
#> 1061           Mcm7       (B-C ∩ C-A)  4  780
#> 1062         Tcf7l1       (B-C ∩ C-A)  4  780
#> 1063  4930523C07Rik       (B-C ∩ C-A)  4  780
#> 1064         Ddx26b       (B-C ∩ C-A)  4  780
#> 1065            Wls       (B-C ∩ C-A)  4  780
#> 1066            Bik       (B-C ∩ C-A)  4  780
#> 1067         Gpr171       (B-C ∩ C-A)  4  780
#> 1068          Nfam1       (B-C ∩ C-A)  4  780
#> 1069          Atxn1       (B-C ∩ C-A)  4  780
#> 1070        Ctnnal1       (B-C ∩ C-A)  4  780
#> 1071          Satb1       (B-C ∩ C-A)  4  780
#> 1072          Rab38       (B-C ∩ C-A)  4  780
#> 1073           Ybx3       (B-C ∩ C-A)  4  780
#> 1074          Prrt1       (B-C ∩ C-A)  4  780
#> 1075           Cry1       (B-C ∩ C-A)  4  780
#> 1076          Cd24a       (B-C ∩ C-A)  4  780
#> 1077         Rnf145       (B-C ∩ C-A)  4  780
#> 1078          Hmgn3       (B-C ∩ C-A)  4  780
#> 1079           Arg1       (B-C ∩ C-A)  4  780
#> 1080  2010016I18Rik       (B-C ∩ C-A)  4  780
#> 1081          Erp29       (B-C ∩ C-A)  4  780
#> 1082          Sept8       (B-C ∩ C-A)  4  780
#> 1083         Il10ra       (B-C ∩ C-A)  4  780
#> 1084        Cbfa2t3       (B-C ∩ C-A)  4  780
#> 1085           Hps4       (B-C ∩ C-A)  4  780
#> 1086          Cxcr3       (B-C ∩ C-A)  4  780
#> 1087            Bsn       (B-C ∩ C-A)  4  780
#> 1088         Pdlim5       (B-C ∩ C-A)  4  780
#> 1089          Zfp35       (B-C ∩ C-A)  4  780
#> 1090          Cnpy3       (B-C ∩ C-A)  4  780
#> 1091          Bspry       (B-C ∩ C-A)  4  780
#> 1092         Slc7a8       (B-C ∩ C-A)  4  780
#> 1093         Tgfbr2       (B-C ∩ C-A)  4  780
#> 1094          Prmt1       (B-C ∩ C-A)  4  780
#> 1095         Ptpn22       (B-C ∩ C-A)  4  780
#> 1096            Syk       (B-C ∩ C-A)  4  780
#> 1097           Cd97       (B-C ∩ C-A)  4  780
#> 1098         Gtpbp4       (B-C ∩ C-A)  4  780
#> 1099         Slain2       (B-C ∩ C-A)  4  780
#> 1100  1700012B07Rik       (B-C ∩ C-A)  4  780
#> 1101           Csf2       (B-C ∩ C-A)  4  780
#> 1102        Col18a1       (B-C ∩ C-A)  4  780
#> 1103       Ppp1r16b       (B-C ∩ C-A)  4  780
#> 1104            Dtl       (B-C ∩ C-A)  4  780
#> 1105        Gm15417       (B-C ∩ C-A)  4  780
#> 1106           Igf2       (B-C ∩ C-A)  4  780
#> 1107          Spns2       (B-C ∩ C-A)  4  780
#> 1108        Dpy19l3       (B-C ∩ C-A)  4  780
#> 1109         Sec31b       (B-C ∩ C-A)  4  780
#> 1110        Skiv2l2       (B-C ∩ C-A)  4  780
#> 1111          Abhd8       (B-C ∩ C-A)  4  780
#> 1112         P2ry10       (B-C ∩ C-A)  4  780
#> 1113          Snx10       (B-C ∩ C-A)  4  780
#> 1114        Gm11613       (B-C ∩ C-A)  4  780
#> 1115          Kank3       (B-C ∩ C-A)  4  780
#> 1116          Sytl3       (B-C ∩ C-A)  4  780
#> 1117           Grap       (B-C ∩ C-A)  4  780
#> 1118          Ttyh3       (B-C ∩ C-A)  4  780
#> 1119          Pcgf6       (B-C ∩ C-A)  4  780
#> 1120           Ets1       (B-C ∩ C-A)  4  780
#> 1121          Litaf       (B-C ∩ C-A)  4  780
#> 1122          Abca1       (B-C ∩ C-A)  4  780
#> 1123          Calca       (B-C ∩ C-A)  4  780
#> 1124          Itgam       (B-C ∩ C-A)  4  780
#> 1125     St6galnac4       (B-C ∩ C-A)  4  780
#> 1126           Cd68       (B-C ∩ C-A)  4  780
#> 1127         Fam84b       (B-C ∩ C-A)  4  780
#> 1128         Tmem71       (B-C ∩ C-A)  4  780
#> 1129    D19Bwg1357e       (B-C ∩ C-A)  4  780
#> 1130          Gfi1b       (B-C ∩ C-A)  4  780
#> 1131         Fam65a       (B-C ∩ C-A)  4  780
#> 1132          Pdia5       (B-C ∩ C-A)  4  780
#> 1133           Pgls       (B-C ∩ C-A)  4  780
#> 1134          Ppm1l       (B-C ∩ C-A)  4  780
#> 1135           Fgd2       (B-C ∩ C-A)  4  780
#> 1136          Cd8b1       (B-C ∩ C-A)  4  780
#> 1137          Impg1       (B-C ∩ C-A)  4  780
#> 1138         Rcbtb2       (B-C ∩ C-A)  4  780
#> 1139          Hmga1       (B-C ∩ C-A)  4  780
#> 1140          Plin3       (B-C ∩ C-A)  4  780
#> 1141      Trp53inp1       (B-C ∩ C-A)  4  780
#> 1142         Rundc1       (B-C ∩ C-A)  4  780
#> 1143         Hs3st1       (B-C ∩ C-A)  4  780
#> 1144        Gm26740       (B-C ∩ C-A)  4  780
#> 1145          Msrb3       (B-C ∩ C-A)  4  780
#> 1146         Kansl2       (B-C ∩ C-A)  4  780
#> 1147          Med12       (B-C ∩ C-A)  4  780
#> 1148          Socs1       (B-C ∩ C-A)  4  780
#> 1149          Tgtp1       (B-C ∩ C-A)  4  780
#> 1150          Lonp2       (B-C ∩ C-A)  4  780
#> 1151           Sat1       (B-C ∩ C-A)  4  780
#> 1152         Eif4e3       (B-C ∩ C-A)  4  780
#> 1153          Prkcb       (B-C ∩ C-A)  4  780
#> 1154           Orc2       (B-C ∩ C-A)  4  780
#> 1155         Fbxo33       (B-C ∩ C-A)  4  780
#> 1156        Cyp4f17       (B-C ∩ C-A)  4  780
#> 1157          Rap2a       (B-C ∩ C-A)  4  780
#> 1158         Cyb5rl       (B-C ∩ C-A)  4  780
#> 1159         Ltb4r1       (B-C ∩ C-A)  4  780
#> 1160         Fam83f       (B-C ∩ C-A)  4  780
#> 1161         Ms4a6c       (B-C ∩ C-A)  4  780
#> 1162            Zp1       (B-C ∩ C-A)  4  780
#> 1163          Traf1       (B-C ∩ C-A)  4  780
#> 1164           Palm       (B-C ∩ C-A)  4  780
#> 1165           Tns1       (B-C ∩ C-A)  4  780
#> 1166        Plekha2       (B-C ∩ C-A)  4  780
#> 1167           Grk5       (B-C ∩ C-A)  4  780
#> 1168         Klrb1b       (B-C ∩ C-A)  4  780
#> 1169  5830444B04Rik       (B-C ∩ C-A)  4  780
#> 1170        Tmem140       (B-C ∩ C-A)  4  780
#> 1171         Pom121       (B-C ∩ C-A)  4  780
#> 1172            Dek       (B-C ∩ C-A)  4  780
#> 1173          Spns3       (B-C ∩ C-A)  4  780
#> 1174         Gm7160       (B-C ∩ C-A)  4  780
#> 1175          Ddx27       (B-C ∩ C-A)  4  780
#> 1176           Ctps       (B-C ∩ C-A)  4  780
#> 1177          Nupl2       (B-C ∩ C-A)  4  780
#> 1178         Il18r1       (B-C ∩ C-A)  4  780
#> 1179       AA415398       (B-C ∩ C-A)  4  780
#> 1180         Lepre1       (B-C ∩ C-A)  4  780
#> 1181           Sgk1       (B-C ∩ C-A)  4  780
#> 1182           Spi1       (B-C ∩ C-A)  4  780
#> 1183        Smarcc1       (B-C ∩ C-A)  4  780
#> 1184         Mgat4a       (B-C ∩ C-A)  4  780
#> 1185        Cd200r1       (B-C ∩ C-A)  4  780
#> 1186          Myo7a       (B-C ∩ C-A)  4  780
#> 1187        Fam129a       (B-C ∩ C-A)  4  780
#> 1188        Gm13572       (B-C ∩ C-A)  4  780
#> 1189        Fam134c       (B-C ∩ C-A)  4  780
#> 1190            Hlf       (B-C ∩ C-A)  4  780
#> 1191           Irf6       (B-C ∩ C-A)  4  780
#> 1192     St6galnac3       (B-C ∩ C-A)  4  780
#> 1193       Traf3ip2       (B-C ∩ C-A)  4  780
#> 1194          Bach2       (B-C ∩ C-A)  4  780
#> 1195          Cnnm2       (B-C ∩ C-A)  4  780
#> 1196          Plcb4       (B-C ∩ C-A)  4  780
#> 1197  2810474O19Rik       (B-C ∩ C-A)  4  780
#> 1198       Slc25a13       (B-C ∩ C-A)  4  780
#> 1199         Zfp575       (B-C ∩ C-A)  4  780
#> 1200        Tubgcp3       (B-C ∩ C-A)  4  780
#> 1201          Nrip1       (B-C ∩ C-A)  4  780
#> 1202        Galnt12       (B-C ∩ C-A)  4  780
#> 1203          Ninj1       (B-C ∩ C-A)  4  780
#> 1204         Ogfrl1       (B-C ∩ C-A)  4  780
#> 1205  I830077J02Rik       (B-C ∩ C-A)  4  780
#> 1206         Homer1       (B-C ∩ C-A)  4  780
#> 1207           Ncf2       (B-C ∩ C-A)  4  780
#> 1208           Tifa       (B-C ∩ C-A)  4  780
#> 1209        Laptm4b       (B-C ∩ C-A)  4  780
#> 1210           Nek3       (B-C ∩ C-A)  4  780
#> 1211          U2af1       (B-C ∩ C-A)  4  780
#> 1212         Tspan3       (B-C ∩ C-A)  4  780
#> 1213          Tifab       (B-C ∩ C-A)  4  780
#> 1214       BC026585       (B-C ∩ C-A)  4  780
#> 1215         H2-DMa       (B-C ∩ C-A)  4  780
#> 1216          Anxa6       (B-C ∩ C-A)  4  780
#> 1217        Slc27a6       (B-C ∩ C-A)  4  780
#> 1218          Tpd52       (B-C ∩ C-A)  4  780
#> 1219          Rap2c       (B-C ∩ C-A)  4  780
#> 1220           Tob1       (B-C ∩ C-A)  4  780
#> 1221           Heca       (B-C ∩ C-A)  4  780
#> 1222        Plekho1       (B-C ∩ C-A)  4  780
#> 1223          Anxa3       (B-C ∩ C-A)  4  780
#> 1224        Alox5ap       (B-C ∩ C-A)  4  780
#> 1225          Cxcr6       (B-C ∩ C-A)  4  780
#> 1226          S1pr3       (B-C ∩ C-A)  4  780
#> 1227           Rpf2       (B-C ∩ C-A)  4  780
#> 1228        Tmem63a       (B-C ∩ C-A)  4  780
#> 1229          Ptcd1       (B-C ∩ C-A)  4  780
#> 1230           Gatb       (B-C ∩ C-A)  4  780
#> 1231           Pus7       (B-C ∩ C-A)  4  780
#> 1232         Trip12       (B-C ∩ C-A)  4  780
#> 1233  8030462N17Rik       (B-C ∩ C-A)  4  780
#> 1234          Capn2       (B-C ∩ C-A)  4  780
#> 1235           Nrp1       (B-C ∩ C-A)  4  780
#> 1236           Hhex       (B-C ∩ C-A)  4  780
#> 1237          Capn5       (B-C ∩ C-A)  4  780
#> 1238          Afap1       (B-C ∩ C-A)  4  780
#> 1239           Hes6       (B-C ∩ C-A)  4  780
#> 1240         Suclg1       (B-C ∩ C-A)  4  780
#> 1241          Ncoa1       (B-C ∩ C-A)  4  780
#> 1242         Zcchc3       (B-C ∩ C-A)  4  780
#> 1243      Tnfrsf10b       (B-C ∩ C-A)  4  780
#> 1244         Tm6sf1       (B-C ∩ C-A)  4  780
#> 1245        Hnrnpdl       (B-C ∩ C-A)  4  780
#> 1246         Plxnb2       (B-C ∩ C-A)  4  780
#> 1247      Serpina3f       (B-C ∩ C-A)  4  780
#> 1248          Timp2       (B-C ∩ C-A)  4  780
#> 1249          Mppe1       (B-C ∩ C-A)  4  780
#> 1250           Dok2       (B-C ∩ C-A)  4  780
#> 1251           Ppan       (B-C ∩ C-A)  4  780
#> 1252          Il6ra       (B-C ∩ C-A)  4  780
#> 1253        Col16a1       (B-C ∩ C-A)  4  780
#> 1254          Plod3       (B-C ∩ C-A)  4  780
#> 1255  4933424M12Rik       (B-C ∩ C-A)  4  780
#> 1256           Nfe2       (B-C ∩ C-A)  4  780
#> 1257        St6gal1       (B-C ∩ C-A)  4  780
#> 1258          Ddx54       (B-C ∩ C-A)  4  780
#> 1259        S100a10       (B-C ∩ C-A)  4  780
#> 1260         Inpp4b       (B-C ∩ C-A)  4  780
#> 1261        Il12rb1       (B-C ∩ C-A)  4  780
#> 1262           Inip       (B-C ∩ C-A)  4  780
#> 1263          Ints1       (B-C ∩ C-A)  4  780
#> 1264            Rb1       (B-C ∩ C-A)  4  780
#> 1265           Coq4       (B-C ∩ C-A)  4  780
#> 1266            Srm       (B-C ∩ C-A)  4  780
#> 1267         Gpr132       (B-C ∩ C-A)  4  780
#> 1268           Ppt2       (B-C ∩ C-A)  4  780
#> 1269        St3gal6       (B-C ∩ C-A)  4  780
#> 1270         Ruvbl1       (B-C ∩ C-A)  4  780
#> 1271  3110082I17Rik       (B-C ∩ C-A)  4  780
#> 1272           Nme4       (B-C ∩ C-A)  4  780
#> 1273        Sh3kbp1       (B-C ∩ C-A)  4  780
#> 1274           Mlf1       (B-C ∩ C-A)  4  780
#> 1275           Aqp3       (B-C ∩ C-A)  4  780
#> 1276        Zdhhc15       (B-C ∩ C-A)  4  780
#> 1277        Ccdc186       (B-C ∩ C-A)  4  780
#> 1278           Dna2       (B-C ∩ C-A)  4  780
#> 1279           Tns4       (B-C ∩ C-A)  4  780
#> 1280           Ier5       (B-C ∩ C-A)  4  780
#> 1281          Ube2o       (B-C ∩ C-A)  4  780
#> 1282          Nsun2       (B-C ∩ C-A)  4  780
#> 1283          Eva1b       (B-C ∩ C-A)  4  780
#> 1284           Yaf2       (B-C ∩ C-A)  4  780
#> 1285       Zc3hav1l       (B-C ∩ C-A)  4  780
#> 1286            Pnn       (B-C ∩ C-A)  4  780
#> 1287          Fnbp4       (B-C ∩ C-A)  4  780
#> 1288         Ankmy2       (B-C ∩ C-A)  4  780
#> 1289          Tmed3       (B-C ∩ C-A)  4  780
#> 1290        Fam117a       (B-C ∩ C-A)  4  780
#> 1291         Pi4k2b       (B-C ∩ C-A)  4  780
#> 1292       Colgalt1       (B-C ∩ C-A)  4  780
#> 1293         Samhd1       (B-C ∩ C-A)  4  780
#> 1294        Dennd1b       (B-C ∩ C-A)  4  780
#> 1295          Birc3       (B-C ∩ C-A)  4  780
#> 1296           Mfn1       (B-C ∩ C-A)  4  780
#> 1297          Hspa2       (B-C ∩ C-A)  4  780
#> 1298          Trub2       (B-C ∩ C-A)  4  780
#> 1299          Nlrc3       (B-C ∩ C-A)  4  780
#> 1300           Gja1       (B-C ∩ C-A)  4  780
#> 1301         Bcl11a       (B-C ∩ C-A)  4  780
#> 1302          Hmgb3       (B-C ∩ C-A)  4  780
#> 1303          Panx1       (B-C ∩ C-A)  4  780
#> 1304           Ccr7       (B-C ∩ C-A)  4  780
#> 1305          Apbb2       (B-C ∩ C-A)  4  780
#> 1306          Mast4       (B-C ∩ C-A)  4  780
#> 1307         Rnase6       (B-C ∩ C-A)  4  780
#> 1308         H2afy2       (B-C ∩ C-A)  4  780
#> 1309          Wdr18       (B-C ∩ C-A)  4  780
#> 1310          Pde4a       (B-C ∩ C-A)  4  780
#> 1311         Rnf166       (B-C ∩ C-A)  4  780
#> 1312           Ncr1       (B-C ∩ C-A)  4  780
#> 1313         Mtmr10       (B-C ∩ C-A)  4  780
#> 1314         Tmem64       (B-C ∩ C-A)  4  780
#> 1315           Ncln       (B-C ∩ C-A)  4  780
#> 1316          Ypel1       (B-C ∩ C-A)  4  780
#> 1317          Wdr34       (B-C ∩ C-A)  4  780
#> 1318          Stab2       (B-C ∩ C-A)  4  780
#> 1319           Rfc2       (B-C ∩ C-A)  4  780
#> 1320         Zfp953       (B-C ∩ C-A)  4  780
#> 1321            Myb       (B-C ∩ C-A)  4  780
#> 1322        Atp13a2       (B-C ∩ C-A)  4  780
#> 1323        Nkiras1       (B-C ∩ C-A)  4  780
#> 1324         Elovl7       (B-C ∩ C-A)  4  780
#> 1325        Tmem206       (B-C ∩ C-A)  4  780
#> 1326          Tcea2       (B-C ∩ C-A)  4  780
#> 1327           Arl1       (B-C ∩ C-A)  4  780
#> 1328          Ckap4       (B-C ∩ C-A)  4  780
#> 1329        Arl6ip6       (B-C ∩ C-A)  4  780
#> 1330           Ly86       (B-C ∩ C-A)  4  780
#> 1331        Gm10384       (B-C ∩ C-A)  4  780
#> 1332           Cds2       (B-C ∩ C-A)  4  780
#> 1333          Galk2       (B-C ∩ C-A)  4  780
#> 1334          Arl4a       (B-C ∩ C-A)  4  780
#> 1335  2700049A03Rik       (B-C ∩ C-A)  4  780
#> 1336  B230217O12Rik       (B-C ∩ C-A)  4  780
#> 1337        Slc11a2       (B-C ∩ C-A)  4  780
#> 1338          Hace1       (B-C ∩ C-A)  4  780
#> 1339          Cdip1       (B-C ∩ C-A)  4  780
#> 1340        Slc6a13       (B-C ∩ C-A)  4  780
#> 1341         Mblac2       (B-C ∩ C-A)  4  780
#> 1342           Cog4       (B-C ∩ C-A)  4  780
#> 1343         Tex264       (B-C ∩ C-A)  4  780
#> 1344         Zbtb49       (B-C ∩ C-A)  4  780
#> 1345          Enpp4       (B-C ∩ C-A)  4  780
#> 1346          Zfp74       (B-C ∩ C-A)  4  780
#> 1347        Wbscr16       (B-C ∩ C-A)  4  780
#> 1348         Klhl30       (B-C ∩ C-A)  4  780
#> 1349           Rgs2       (B-C ∩ C-A)  4  780
#> 1350           Acp5       (B-C ∩ C-A)  4  780
#> 1351          Alcam       (B-C ∩ C-A)  4  780
#> 1352          Hyal2       (B-C ∩ C-A)  4  780
#> 1353           Mycl       (B-C ∩ C-A)  4  780
#> 1354           Rara       (B-C ∩ C-A)  4  780
#> 1355          Unc50       (B-C ∩ C-A)  4  780
#> 1356        Nsmce4a       (B-C ∩ C-A)  4  780
#> 1357         Arrdc4       (B-C ∩ C-A)  4  780
#> 1358           Ier2       (B-C ∩ C-A)  4  780
#> 1359         Atrnl1       (B-C ∩ C-A)  4  780
#> 1360           Dkk3       (B-C ∩ C-A)  4  780
#> 1361         Kcnk12       (B-C ∩ C-A)  4  780
#> 1362           Gga2       (B-C ∩ C-A)  4  780
#> 1363          Irgm2       (B-C ∩ C-A)  4  780
#> 1364          Wdfy4       (B-C ∩ C-A)  4  780
#> 1365         Trerf1       (B-C ∩ C-A)  4  780
#> 1366          Mapk6       (B-C ∩ C-A)  4  780
#> 1367         Prss34       (B-C ∩ C-A)  4  780
#> 1368           Il7r       (B-C ∩ C-A)  4  780
#> 1369         Plxnc1       (B-C ∩ C-A)  4  780
#> 1370         Gimap8       (B-C ∩ C-A)  4  780
#> 1371           Acy1       (B-C ∩ C-A)  4  780
#> 1372         Mettl4       (B-C ∩ C-A)  4  780
#> 1373  0610010F05Rik       (B-C ∩ C-A)  4  780
#> 1374          H2-Oa       (B-C ∩ C-A)  4  780
#> 1375            Mbp       (B-C ∩ C-A)  4  780
#> 1376          Trip4       (B-C ∩ C-A)  4  780
#> 1377         Snrpa1       (B-C ∩ C-A)  4  780
#> 1378        Pcyox1l       (B-C ∩ C-A)  4  780
#> 1379           Ly6a       (B-C ∩ C-A)  4  780
#> 1380          Nsun4       (B-C ∩ C-A)  4  780
#> 1381         Zbtb42       (B-C ∩ C-A)  4  780
#> 1382         Zfp202       (B-C ∩ C-A)  4  780
#> 1383         R3hdm1       (B-C ∩ C-A)  4  780
#> 1384      Uhrf1bp1l       (B-C ∩ C-A)  4  780
#> 1385         Agpat3       (B-C ∩ C-A)  4  780
#> 1386        Hsd17b7       (B-C ∩ C-A)  4  780
#> 1387         Pi4k2a       (B-C ∩ C-A)  4  780
#> 1388         Calcrl       (B-C ∩ C-A)  4  780
#> 1389          Lzts2       (B-C ∩ C-A)  4  780
#> 1390          Prkcd       (B-C ∩ C-A)  4  780
#> 1391        Pglyrp2       (B-C ∩ C-A)  4  780
#> 1392          Rrp1b       (B-C ∩ C-A)  4  780
#> 1393         Zbtb11       (B-C ∩ C-A)  4  780
#> 1394           Dgke       (B-C ∩ C-A)  4  780
#> 1395         Hmgxb4       (B-C ∩ C-A)  4  780
#> 1396           Mmp8       (B-C ∩ C-A)  4  780
#> 1397          Mcph1       (B-C ∩ C-A)  4  780
#> 1398           Lax1       (B-C ∩ C-A)  4  780
#> 1399          G6pc3       (B-C ∩ C-A)  4  780
#> 1400        Msantd4       (B-C ∩ C-A)  4  780
#> 1401        Gm22107       (B-C ∩ C-A)  4  780
#> 1402          Ngly1       (B-C ∩ C-A)  4  780
#> 1403           Myo6       (B-C ∩ C-A)  4  780
#> 1404        Cdk2ap2       (B-C ∩ C-A)  4  780
#> 1405           Pcm1       (B-C ∩ C-A)  4  780
#> 1406          Lrrk1       (B-C ∩ C-A)  4  780
#> 1407          Nr1d1       (B-C ∩ C-A)  4  780
#> 1408        Osbpl1a       (B-C ∩ C-A)  4  780
#> 1409          Syne1       (B-C ∩ C-A)  4  780
#> 1410          Parp8       (B-C ∩ C-A)  4  780
#> 1411          Srsf2       (B-C ∩ C-A)  4  780
#> 1412         Zfp629       (B-C ∩ C-A)  4  780
#> 1413          Ddx55       (B-C ∩ C-A)  4  780
#> 1414          Ercc6       (B-C ∩ C-A)  4  780
#> 1415           Tia1       (B-C ∩ C-A)  4  780
#> 1416          Qtrt1       (B-C ∩ C-A)  4  780
#> 1417           Pdcl       (B-C ∩ C-A)  4  780
#> 1418          Lphn1       (B-C ∩ C-A)  4  780
#> 1419          Krt10       (B-C ∩ C-A)  4  780
#> 1420        Psmc3ip       (B-C ∩ C-A)  4  780
#> 1421          Cbwd1       (B-C ∩ C-A)  4  780
#> 1422         Ldoc1l       (B-C ∩ C-A)  4  780
#> 1423         Tuba1a       (B-C ∩ C-A)  4  780
#> 1424           Msra       (B-C ∩ C-A)  4  780
#> 1425           Phf1       (B-C ∩ C-A)  4  780
#> 1426         Ehhadh       (B-C ∩ C-A)  4  780
#> 1427          Dusp1       (B-C ∩ C-A)  4  780
#> 1428          Dcp1b       (B-C ∩ C-A)  4  780
#> 1429  2610020H08Rik       (B-C ∩ C-A)  4  780
#> 1430           Pld4       (B-C ∩ C-A)  4  780
#> 1431           Ezh1       (B-C ∩ C-A)  4  780
#> 1432          Lace1       (B-C ∩ C-A)  4  780
#> 1433           Klf6       (B-C ∩ C-A)  4  780
#> 1434          Anxa9       (B-C ∩ C-A)  4  780
#> 1435  1110008L16Rik       (B-C ∩ C-A)  4  780
#> 1436      Rab11fip1       (B-C ∩ C-A)  4  780
#> 1437          Nxpe3       (B-C ∩ C-A)  4  780
#> 1438          Sirt3       (B-C ∩ C-A)  4  780
#> 1439         Pcmtd1       (B-C ∩ C-A)  4  780
#> 1440           Utrn       (B-C ∩ C-A)  4  780
#> 1441           Klf5       (B-C ∩ C-A)  4  780
#> 1442           Acp1       (B-C ∩ C-A)  4  780
#> 1443          Smpd5       (B-C ∩ C-A)  4  780
#> 1444            Zak       (B-C ∩ C-A)  4  780
#> 1445          Syce2       (B-C ∩ C-A)  4  780
#> 1446           Espn       (B-C ∩ C-A)  4  780
#> 1447          Spidr       (B-C ∩ C-A)  4  780
#> 1448         Il1rap       (B-C ∩ C-A)  4  780
#> 1449          Cyth4       (B-C ∩ C-A)  4  780
#> 1450           Pld3       (B-C ∩ C-A)  4  780
#> 1451         Zfp563       (B-C ∩ C-A)  4  780
#> 1452          Cytip       (B-C ∩ C-A)  4  780
#> 1453           Cast       (B-C ∩ C-A)  4  780
#> 1454      Slfn10-ps       (B-C ∩ C-A)  4  780
#> 1455           Fdps       (B-C ∩ C-A)  4  780
#> 1456           Ing4       (B-C ∩ C-A)  4  780
#> 1457          Adpgk       (B-C ∩ C-A)  4  780
#> 1458        Sh2d1b1       (B-C ∩ C-A)  4  780
#> 1459        Sec14l1       (B-C ∩ C-A)  4  780
#> 1460          Rad18       (B-C ∩ C-A)  4  780
#> 1461           Ibtk       (B-C ∩ C-A)  4  780
#> 1462         Zfp395       (B-C ∩ C-A)  4  780
#> 1463        Gm16740       (B-C ∩ C-A)  4  780
#> 1464           Def8       (B-C ∩ C-A)  4  780
#> 1465          Trpv2       (B-C ∩ C-A)  4  780
#> 1466         Dnaja4       (B-C ∩ C-A)  4  780
#> 1467          Dhx32       (B-C ∩ C-A)  4  780
#> 1468         Dock11       (B-C ∩ C-A)  4  780
#> 1469          Mrpl1       (B-C ∩ C-A)  4  780
#> 1470           Hmbs       (B-C ∩ C-A)  4  780
#> 1471        Slc35f2       (B-C ∩ C-A)  4  780
#> 1472          Trib2       (B-C ∩ C-A)  4  780
#> 1473          Nup43       (B-C ∩ C-A)  4  780
#> 1474         Atp2b1       (B-C ∩ C-A)  4  780
#> 1475           Plek       (B-C ∩ C-A)  4  780
#> 1476         Gpr160       (B-C ∩ C-A)  4  780
#> 1477        Slc35d2       (B-C ∩ C-A)  4  780
#> 1478       Arhgap17       (B-C ∩ C-A)  4  780
#> 1479         S100a6       (B-C ∩ C-A)  4  780
#> 1480         Gatsl2       (B-C ∩ C-A)  4  780
#> 1481         Erlin1       (B-C ∩ C-A)  4  780
#> 1482           Tyw3       (B-C ∩ C-A)  4  780
#> 1483          Padi4       (B-C ∩ C-A)  4  780
#> 1484        Gprasp2       (B-C ∩ C-A)  4  780
#> 1485          Uqcc1       (B-C ∩ C-A)  4  780
#> 1486        Rasgrp1       (B-C ∩ C-A)  4  780
#> 1487         Kdelc1       (B-C ∩ C-A)  4  780
#> 1488         Zfp934       (B-C ∩ C-A)  4  780
#> 1489         Hrsp12       (B-C ∩ C-A)  4  780
#> 1490        Tcp11l1       (B-C ∩ C-A)  4  780
#> 1491          Naa15       (B-C ∩ C-A)  4  780
#> 1492          Usp46       (B-C ∩ C-A)  4  780
#> 1493            Srl       (B-C ∩ C-A)  4  780
#> 1494         Rnf168       (B-C ∩ C-A)  4  780
#> 1495           Ctsf       (B-C ∩ C-A)  4  780
#> 1496  2010111I01Rik       (B-C ∩ C-A)  4  780
#> 1497           Aak1       (B-C ∩ C-A)  4  780
#> 1498         Scarb1       (B-C ∩ C-A)  4  780
#> 1499         Chchd4       (B-C ∩ C-A)  4  780
#> 1500          Brwd3       (B-C ∩ C-A)  4  780
#> 1501           Ecm1       (B-C ∩ C-A)  4  780
#> 1502         Trim44       (B-C ∩ C-A)  4  780
#> 1503         Fchsd1       (B-C ∩ C-A)  4  780
#> 1504         Bcl2l1       (B-C ∩ C-A)  4  780
#> 1505           Sord       (B-C ∩ C-A)  4  780
#> 1506            Cd9       (B-C ∩ C-A)  4  780
#> 1507  9130401M01Rik       (B-C ∩ C-A)  4  780
#> 1508         Lztfl1       (B-C ∩ C-A)  4  780
#> 1509           Crot       (B-C ∩ C-A)  4  780
#> 1510          Zfp68       (B-C ∩ C-A)  4  780
#> 1511           Msh6       (B-C ∩ C-A)  4  780
#> 1512           Gatm       (B-C ∩ C-A)  4  780
#> 1513         Asrgl1       (B-C ∩ C-A)  4  780
#> 1514         Cdkn2c       (B-C ∩ C-A)  4  780
#> 1515          Nicn1       (B-C ∩ C-A)  4  780
#> 1516          Casd1       (B-C ∩ C-A)  4  780
#> 1517       Tnfrsf18       (B-C ∩ C-A)  4  780
#> 1518        Rasgrp2       (B-C ∩ C-A)  4  780
#> 1519         Polr1b       (B-C ∩ C-A)  4  780
#> 1520          Rab3a       (B-C ∩ C-A)  4  780
#> 1521          Xlr4b       (B-C ∩ C-A)  4  780
#> 1522         Apol7e       (B-C ∩ C-A)  4  780
#> 1523          Abcd3       (B-C ∩ C-A)  4  780
#> 1524            Tec       (B-C ∩ C-A)  4  780
#> 1525           Nav2       (B-C ∩ C-A)  4  780
#> 1526        Zfp280c       (B-C ∩ C-A)  4  780
#> 1527        Rps6ka3       (B-C ∩ C-A)  4  780
#> 1528        Tnfsf12       (B-C ∩ C-A)  4  780
#> 1529        Tmsb15a       (B-C ∩ C-A)  4  780
#> 1530          Itpr2       (B-C ∩ C-A)  4  780
#> 1531          Dusp2       (B-C ∩ C-A)  4  780
#> 1532         Bcl11b       (B-C ∩ C-A)  4  780
#> 1533          Mdfic       (B-C ∩ C-A)  4  780
#> 1534       Selenbp1       (B-C ∩ C-A)  4  780
#> 1535           Vwa9       (B-C ∩ C-A)  4  780
#> 1536        Gm11427       (B-C ∩ C-A)  4  780
#> 1537         Mpped2       (B-C ∩ C-A)  4  780
#> 1538           Gnl3       (B-C ∩ C-A)  4  780
#> 1539         Atp8a1       (B-C ∩ C-A)  4  780
#> 1540         Plxdc2       (B-C ∩ C-A)  4  780
#> 1541        Gm15738       (B-C ∩ C-A)  4  780
#> 1542          Golm1       (B-C ∩ C-A)  4  780
#> 1543        Supv3l1       (B-C ∩ C-A)  4  780
#> 1544         Gimap3       (B-C ∩ C-A)  4  780
#> 1545        Tbc1d23       (B-C ∩ C-A)  4  780
#> 1546          Scmh1       (B-C ∩ C-A)  4  780
#> 1547          Phf20       (B-C ∩ C-A)  4  780
#> 1548          Sfxn1       (B-C ∩ C-A)  4  780
#> 1549          Sart3       (B-C ∩ C-A)  4  780
#> 1550         Nos1ap       (B-C ∩ C-A)  4  780
#> 1551 9930111J21Rik1       (B-C ∩ C-A)  4  780
#> 1552  1700112E06Rik       (B-C ∩ C-A)  4  780
#> 1553         Ociad2       (B-C ∩ C-A)  4  780
#> 1554         Rad54l       (B-C ∩ C-A)  4  780
#> 1555         Ccdc77       (B-C ∩ C-A)  4  780
#> 1556       Tctex1d2       (B-C ∩ C-A)  4  780
#> 1557          Nfkb2       (B-C ∩ C-A)  4  780
#> 1558           Lfng       (B-C ∩ C-A)  4  780
#> 1559          Gng12       (B-C ∩ C-A)  4  780
#> 1560          Exoc8       (B-C ∩ C-A)  4  780
#> 1561            Gem       (B-C ∩ C-A)  4  780
#> 1562           Zbp1       (B-C ∩ C-A)  4  780
#> 1563           Clta       (B-C ∩ C-A)  4  780
#> 1564          Spsb3       (B-C ∩ C-A)  4  780
#> 1565          Ddx17       (B-C ∩ C-A)  4  780
#> 1566           Snx8       (B-C ∩ C-A)  4  780
#> 1567          Zmym3       (B-C ∩ C-A)  4  780
#> 1568           Mcat       (B-C ∩ C-A)  4  780
#> 1569          Pole2       (B-C ∩ C-A)  4  780
#> 1570        Tmem39b       (B-C ∩ C-A)  4  780
#> 1571          Atp4a       (B-C ∩ C-A)  4  780
#> 1572          Fnip2       (B-C ∩ C-A)  4  780
#> 1573          Ift57       (B-C ∩ C-A)  4  780
#> 1574            Pxk       (B-C ∩ C-A)  4  780
#> 1575          Grina       (B-C ∩ C-A)  4  780
#> 1576         Rnf157       (B-C ∩ C-A)  4  780
#> 1577        Gm28529       (B-C ∩ C-A)  4  780
#> 1578           Mlkl       (B-C ∩ C-A)  4  780
#> 1579  4833439L19Rik       (B-C ∩ C-A)  4  780
#> 1580         Bckdhb       (B-C ∩ C-A)  4  780
#> 1581           E2f7       (B-C ∩ C-A)  4  780
#> 1582        Ubash3a       (B-C ∩ C-A)  4  780
#> 1583       Slc26a11       (B-C ∩ C-A)  4  780
#> 1584        Fastkd2       (B-C ∩ C-A)  4  780
#> 1585         D2hgdh       (B-C ∩ C-A)  4  780
#> 1586         Insig2       (B-C ∩ C-A)  4  780
#> 1587        Bcl2l11       (B-C ∩ C-A)  4  780
#> 1588         Dusp10       (B-C ∩ C-A)  4  780
#> 1589           Asb1       (B-C ∩ C-A)  4  780
#> 1590          Lacc1       (B-C ∩ C-A)  4  780
#> 1591         Tnrc6a       (B-C ∩ C-A)  4  780
#> 1592        Cysltr2       (B-C ∩ C-A)  4  780
#> 1593          Sh2b3       (B-C ∩ C-A)  4  780
#> 1594          Arntl       (B-C ∩ C-A)  4  780
#> 1595          Gria3       (B-C ∩ C-A)  4  780
#> 1596          Hmces       (B-C ∩ C-A)  4  780
#> 1597        Zfp119b       (B-C ∩ C-A)  4  780
#> 1598          Krit1       (B-C ∩ C-A)  4  780
#> 1599           Zfp3       (B-C ∩ C-A)  4  780
#> 1600          Myo5a       (B-C ∩ C-A)  4  780
#> 1601         Nudt22       (B-C ∩ C-A)  4  780
#> 1602         Wdsub1       (B-C ∩ C-A)  4  780
#> 1603           Urb2       (B-C ∩ C-A)  4  780
#> 1604          Plod2       (B-C ∩ C-A)  4  780
#> 1605           Pigu       (B-C ∩ C-A)  4  780
#> 1606            Ahr       (B-C ∩ C-A)  4  780
#> 1607          Mmgt1       (B-C ∩ C-A)  4  780
#> 1608         Sft2d1       (B-C ∩ C-A)  4  780
#> 1609       Gimap1os       (B-C ∩ C-A)  4  780
#> 1610            Cd5       (B-C ∩ C-A)  4  780
#> 1611         Zswim3       (B-C ∩ C-A)  4  780
#> 1612          Ric8b       (B-C ∩ C-A)  4  780
#> 1613           Rhob       (B-C ∩ C-A)  4  780
#> 1614           Smox       (B-C ∩ C-A)  4  780
#> 1615          Rpl22       (B-C ∩ C-A)  4  780
#> 1616          Itih5       (B-C ∩ C-A)  4  780
#> 1617          Thoc6       (B-C ∩ C-A)  4  780
#> 1618        Hnrnpll       (B-C ∩ C-A)  4  780
#> 1619       Timeless       (B-C ∩ C-A)  4  780
#> 1620        Slc52a2       (B-C ∩ C-A)  4  780
#> 1621          Troap       (B-C ∩ C-A)  4  780
#> 1622         Tbc1d5       (B-C ∩ C-A)  4  780
#> 1623          Cenph       (B-C ∩ C-A)  4  780
#> 1624        Zcchc18       (B-C ∩ C-A)  4  780
#> 1625           Sypl       (B-C ∩ C-A)  4  780
#> 1626          Dkkl1       (B-C ∩ C-A)  4  780
#> 1627         Man2c1       (B-C ∩ C-A)  4  780
#> 1628          Gna13       (B-C ∩ C-A)  4  780
#> 1629           Gab3       (B-C ∩ C-A)  4  780
#> 1630  1700030K09Rik       (B-C ∩ C-A)  4  780
#> 1631         Ccl27a       (B-C ∩ C-A)  4  780
#> 1632           Ppcs       (B-C ∩ C-A)  4  780
#> 1633       BC028528       (B-C ∩ C-A)  4  780
#> 1634         Ccdc66       (B-C ∩ C-A)  4  780
#> 1635         Nudcd2       (B-C ∩ C-A)  4  780
#> 1636          Bace1       (B-C ∩ C-A)  4  780
#> 1637          Prdm4       (B-C ∩ C-A)  4  780
#> 1638          Sh2b2       (B-C ∩ C-A)  4  780
#> 1639           Bst2       (B-C ∩ C-A)  4  780
#> 1640          Pank1       (B-C ∩ C-A)  4  780
#> 1641  1810055G02Rik       (B-C ∩ C-A)  4  780
#> 1642        Irf2bpl       (B-C ∩ C-A)  4  780
#> 1643           Pkp3       (B-C ∩ C-A)  4  780
#> 1644          Larp7       (B-C ∩ C-A)  4  780
#> 1645         Abhd11       (B-C ∩ C-A)  4  780
#> 1646          Dnm1l       (B-C ∩ C-A)  4  780
#> 1647          Rab4a       (B-C ∩ C-A)  4  780
#> 1648  5430405H02Rik       (B-C ∩ C-A)  4  780
#> 1649           Rcl1       (B-C ∩ C-A)  4  780
#> 1650           Atn1       (B-C ∩ C-A)  4  780
#> 1651          Rtkn2       (B-C ∩ C-A)  4  780
#> 1652           Rhoc       (B-C ∩ C-A)  4  780
#> 1653         Klhl23       (B-C ∩ C-A)  4  780
#> 1654          Rab31       (B-C ∩ C-A)  4  780
#> 1655          Stim2       (B-C ∩ C-A)  4  780
#> 1656         Igfbp7       (B-C ∩ C-A)  4  780
#> 1657           Galm       (B-C ∩ C-A)  4  780
#> 1658          Itfg3       (B-C ∩ C-A)  4  780
#> 1659          Ascc3       (B-C ∩ C-A)  4  780
#> 1660           Klf2       (B-C ∩ C-A)  4  780
#> 1661          Dhx58       (B-C ∩ C-A)  4  780
#> 1662           Prkx       (B-C ∩ C-A)  4  780
#> 1663         Tspan4       (B-C ∩ C-A)  4  780
#> 1664         Rab3ip       (B-C ∩ C-A)  4  780
#> 1665         Zfp213       (B-C ∩ C-A)  4  780
#> 1666           Med1       (B-C ∩ C-A)  4  780
#> 1667  4933404O12Rik       (B-C ∩ C-A)  4  780
#> 1668          Thada       (B-C ∩ C-A)  4  780
#> 1669       Slc25a16       (B-C ∩ C-A)  4  780
#> 1670          Lair1       (B-C ∩ C-A)  4  780
#> 1671  6820431F20Rik       (B-C ∩ C-A)  4  780
#> 1672         Tsen54       (B-C ∩ C-A)  4  780
#> 1673          Mtfr2       (B-C ∩ C-A)  4  780
#> 1674        Gadd45a       (B-C ∩ C-A)  4  780
#> 1675        Slc37a3       (B-C ∩ C-A)  4  780
#> 1676         Nt5dc3       (B-C ∩ C-A)  4  780
#> 1677         Snapc2       (B-C ∩ C-A)  4  780
#> 1678         Kif18a       (B-C ∩ C-A)  4  780
#> 1679          Itpr3       (B-C ∩ C-A)  4  780
#> 1680           Lmna       (B-C ∩ C-A)  4  780
#> 1681         Zfp319       (B-C ∩ C-A)  4  780
#> 1682           Tlr2       (B-C ∩ C-A)  4  780
#> 1683          Alas1       (B-C ∩ C-A)  4  780
#> 1684           Hpse       (B-C ∩ C-A)  4  780
#> 1685         Lrrc8b       (B-C ∩ C-A)  4  780
#> 1686           Pex3       (B-C ∩ C-A)  4  780
#> 1687           Gnal       (B-C ∩ C-A)  4  780
#> 1688           Cd3d       (B-C ∩ C-A)  4  780
#> 1689        Thumpd3       (B-C ∩ C-A)  4  780
#> 1690          Neil1       (B-C ∩ C-A)  4  780
#> 1691         Pxylp1       (B-C ∩ C-A)  4  780
#> 1692         Clec2i       (B-C ∩ C-A)  4  780
#> 1693          Kntc1       (B-C ∩ C-A)  4  780
#> 1694           Chd7       (B-C ∩ C-A)  4  780
#> 1695        Gm28053       (B-C ∩ C-A)  4  780
#> 1696         Mfsd7b       (B-C ∩ C-A)  4  780
#> 1697          Supt6       (B-C ∩ C-A)  4  780
#> 1698          Zfp26       (B-C ∩ C-A)  4  780
#> 1699         Tagln2             (B-A)  5  225
#> 1700       Ankrd13a             (B-A)  5  225
#> 1701          Coro7             (B-A)  5  225
#> 1702           Dbnl             (B-A)  5  225
#> 1703            Xpc             (B-A)  5  225
#> 1704           Il16             (B-A)  5  225
#> 1705            Gyg             (B-A)  5  225
#> 1706           Fut7             (B-A)  5  225
#> 1707           Bbs2             (B-A)  5  225
#> 1708           Jrkl             (B-A)  5  225
#> 1709        Pip4k2a             (B-A)  5  225
#> 1710           Polh             (B-A)  5  225
#> 1711         Spice1             (B-A)  5  225
#> 1712        Zc3hav1             (B-A)  5  225
#> 1713        Ppp2r5a             (B-A)  5  225
#> 1714          Fgf13             (B-A)  5  225
#> 1715         Shcbp1             (B-A)  5  225
#> 1716           Vav3             (B-A)  5  225
#> 1717          Aldh2             (B-A)  5  225
#> 1718  1700026L06Rik             (B-A)  5  225
#> 1719           Tprg             (B-A)  5  225
#> 1720         Cyb5r3             (B-A)  5  225
#> 1721           Fhit             (B-A)  5  225
#> 1722           Rin3             (B-A)  5  225
#> 1723            Ivd             (B-A)  5  225
#> 1724   RP24-242N1.2             (B-A)  5  225
#> 1725        Tsc22d3             (B-A)  5  225
#> 1726           Wee1             (B-A)  5  225
#> 1727        Laptm4a             (B-A)  5  225
#> 1728          Abca7             (B-A)  5  225
#> 1729           Jag2             (B-A)  5  225
#> 1730          Cdc20             (B-A)  5  225
#> 1731          Ccnb2             (B-A)  5  225
#> 1732          Kcnk5             (B-A)  5  225
#> 1733  1110059E24Rik             (B-A)  5  225
#> 1734         Ss18l1             (B-A)  5  225
#> 1735          Cep89             (B-A)  5  225
#> 1736          Wdr59             (B-A)  5  225
#> 1737  D630008O14Rik             (B-A)  5  225
#> 1738         Fbxo10             (B-A)  5  225
#> 1739         Ctdspl             (B-A)  5  225
#> 1740          Abhd4             (B-A)  5  225
#> 1741         Cobll1             (B-A)  5  225
#> 1742  1600020E01Rik             (B-A)  5  225
#> 1743           Rfc3             (B-A)  5  225
#> 1744          Itpkc             (B-A)  5  225
#> 1745          Disc1             (B-A)  5  225
#> 1746  4833418N02Rik             (B-A)  5  225
#> 1747           Capg             (B-A)  5  225
#> 1748          Ypel5             (B-A)  5  225
#> 1749         Tspan6             (B-A)  5  225
#> 1750          Paqr8             (B-A)  5  225
#> 1751         Fam57a             (B-A)  5  225
#> 1752           Sdk1             (B-A)  5  225
#> 1753          Ftsj2             (B-A)  5  225
#> 1754          Klhl5             (B-A)  5  225
#> 1755       Rasgef1b             (B-A)  5  225
#> 1756       Slc22a17             (B-A)  5  225
#> 1757           Prc1             (B-A)  5  225
#> 1758        Plekhm2             (B-A)  5  225
#> 1759          Stap1             (B-A)  5  225
#> 1760           Rfx1             (B-A)  5  225
#> 1761          Uvrag             (B-A)  5  225
#> 1762          Abcd2             (B-A)  5  225
#> 1763           Gbp5             (B-A)  5  225
#> 1764  1110038B12Rik             (B-A)  5  225
#> 1765         Hoxa10             (B-A)  5  225
#> 1766         Fam72a             (B-A)  5  225
#> 1767        Tmem194             (B-A)  5  225
#> 1768           Rffl             (B-A)  5  225
#> 1769          Otud3             (B-A)  5  225
#> 1770           Dgkz             (B-A)  5  225
#> 1771        Abhd17b             (B-A)  5  225
#> 1772         Kif20b             (B-A)  5  225
#> 1773         Scamp3             (B-A)  5  225
#> 1774           Mzt1             (B-A)  5  225
#> 1775           Klf1             (B-A)  5  225
#> 1776  A630001G21Rik             (B-A)  5  225
#> 1777          Arrb2             (B-A)  5  225
#> 1778           Brd8             (B-A)  5  225
#> 1779          Rnf13             (B-A)  5  225
#> 1780          Prr11             (B-A)  5  225
#> 1781         Cep85l             (B-A)  5  225
#> 1782            Kmo             (B-A)  5  225
#> 1783         Zfp664             (B-A)  5  225
#> 1784        Nckipsd             (B-A)  5  225
#> 1785          Birc5             (B-A)  5  225
#> 1786          Xrcc2             (B-A)  5  225
#> 1787          Dgat1             (B-A)  5  225
#> 1788  2810414N06Rik             (B-A)  5  225
#> 1789         Erlin2             (B-A)  5  225
#> 1790            Fuk             (B-A)  5  225
#> 1791       Ighv1-75             (B-A)  5  225
#> 1792          Rbm44             (B-A)  5  225
#> 1793            Tk1             (B-A)  5  225
#> 1794          Zfpm1             (B-A)  5  225
#> 1795         Zfp236             (B-A)  5  225
#> 1796         Cep152             (B-A)  5  225
#> 1797  9530077C05Rik             (B-A)  5  225
#> 1798        Tmem55a             (B-A)  5  225
#> 1799          Tmtc1             (B-A)  5  225
#> 1800            Pml             (B-A)  5  225
#> 1801  D230025D16Rik             (B-A)  5  225
#> 1802           B9d1             (B-A)  5  225
#> 1803          Aim1l             (B-A)  5  225
#> 1804          P2rx1             (B-A)  5  225
#> 1805          Chtf8             (B-A)  5  225
#> 1806  2310075C17Rik             (B-A)  5  225
#> 1807        Gm16973             (B-A)  5  225
#> 1808          Ptgs1             (B-A)  5  225
#> 1809        Cyp11a1             (B-A)  5  225
#> 1810        Eif2ak2             (B-A)  5  225
#> 1811         Sh3bp4             (B-A)  5  225
#> 1812      Bambi-ps1             (B-A)  5  225
#> 1813         Cd2bp2             (B-A)  5  225
#> 1814          Rai14             (B-A)  5  225
#> 1815           Gpd2             (B-A)  5  225
#> 1816         Glt8d1             (B-A)  5  225
#> 1817          Atg9a             (B-A)  5  225
#> 1818  RP23-242C19.8             (B-A)  5  225
#> 1819          Akip1             (B-A)  5  225
#> 1820          Slfn8             (B-A)  5  225
#> 1821           Tpx2             (B-A)  5  225
#> 1822       Mettl7a1             (B-A)  5  225
#> 1823           Acp2             (B-A)  5  225
#> 1824  4933407K13Rik             (B-A)  5  225
#> 1825         Rhbdf2             (B-A)  5  225
#> 1826         Fancd2             (B-A)  5  225
#> 1827       Sh3bgrl2             (B-A)  5  225
#> 1828          Socs2             (B-A)  5  225
#> 1829         Zfp692             (B-A)  5  225
#> 1830        Zdhhc21             (B-A)  5  225
#> 1831         Sh3tc1             (B-A)  5  225
#> 1832        Gm14230             (B-A)  5  225
#> 1833       Fam160b2             (B-A)  5  225
#> 1834         Gm8995             (B-A)  5  225
#> 1835         Ckap2l             (B-A)  5  225
#> 1836       AB124611             (B-A)  5  225
#> 1837         Ddx19b             (B-A)  5  225
#> 1838          Wdr90             (B-A)  5  225
#> 1839         Pik3cb             (B-A)  5  225
#> 1840         Strip2             (B-A)  5  225
#> 1841           Bmp1             (B-A)  5  225
#> 1842       Marveld2             (B-A)  5  225
#> 1843         Man2a2             (B-A)  5  225
#> 1844           Nptn             (B-A)  5  225
#> 1845         Mroh2a             (B-A)  5  225
#> 1846          Pdcd7             (B-A)  5  225
#> 1847         Ahctf1             (B-A)  5  225
#> 1848         Atp1a1             (B-A)  5  225
#> 1849           Ficd             (B-A)  5  225
#> 1850         Zfp110             (B-A)  5  225
#> 1851       Ighv1-83             (B-A)  5  225
#> 1852           Car3             (B-A)  5  225
#> 1853         Qtrtd1             (B-A)  5  225
#> 1854         Fam45a             (B-A)  5  225
#> 1855           Fgd3             (B-A)  5  225
#> 1856         Zdhhc7             (B-A)  5  225
#> 1857         Clstn3             (B-A)  5  225
#> 1858          Pcgf5             (B-A)  5  225
#> 1859        Tmem63b             (B-A)  5  225
#> 1860         Zfp109             (B-A)  5  225
#> 1861           E4f1             (B-A)  5  225
#> 1862         Zfp951             (B-A)  5  225
#> 1863          Rsph9             (B-A)  5  225
#> 1864         Dcaf10             (B-A)  5  225
#> 1865            Smo             (B-A)  5  225
#> 1866       Hist1h4i             (B-A)  5  225
#> 1867         Zfp429             (B-A)  5  225
#> 1868         Tmem80             (B-A)  5  225
#> 1869         Trim65             (B-A)  5  225
#> 1870        Fam210a             (B-A)  5  225
#> 1871          Kif2c             (B-A)  5  225
#> 1872           Cln8             (B-A)  5  225
#> 1873         C77080             (B-A)  5  225
#> 1874         Supt16             (B-A)  5  225
#> 1875        Tmem199             (B-A)  5  225
#> 1876          Tmem9             (B-A)  5  225
#> 1877  A230072C01Rik             (B-A)  5  225
#> 1878        St3gal1             (B-A)  5  225
#> 1879       Arhgap26             (B-A)  5  225
#> 1880      Rab11fip3             (B-A)  5  225
#> 1881          Stau2             (B-A)  5  225
#> 1882         Hexim2             (B-A)  5  225
#> 1883          Ppil4             (B-A)  5  225
#> 1884           Irf5             (B-A)  5  225
#> 1885           Rtp4             (B-A)  5  225
#> 1886           Sos2             (B-A)  5  225
#> 1887          Vars2             (B-A)  5  225
#> 1888          Ndor1             (B-A)  5  225
#> 1889         Shkbp1             (B-A)  5  225
#> 1890           Mavs             (B-A)  5  225
#> 1891        Gm14018             (B-A)  5  225
#> 1892        Ube2cbp             (B-A)  5  225
#> 1893          Casc4             (B-A)  5  225
#> 1894         Klhl12             (B-A)  5  225
#> 1895         Mfsd11             (B-A)  5  225
#> 1896         Dnase1             (B-A)  5  225
#> 1897           Pcnt             (B-A)  5  225
#> 1898           Sc5d             (B-A)  5  225
#> 1899        Slc27a1             (B-A)  5  225
#> 1900           Pogk             (B-A)  5  225
#> 1901         Zdhhc8             (B-A)  5  225
#> 1902            Nin             (B-A)  5  225
#> 1903          Neil3             (B-A)  5  225
#> 1904         Ptpn18             (B-A)  5  225
#> 1905          Dscr3             (B-A)  5  225
#> 1906          Snx27             (B-A)  5  225
#> 1907            App             (B-A)  5  225
#> 1908         Nufip2             (B-A)  5  225
#> 1909         P2ry12             (B-A)  5  225
#> 1910          Dhx29             (B-A)  5  225
#> 1911          Kif11             (B-A)  5  225
#> 1912            Shf             (B-A)  5  225
#> 1913        Tmprss3             (B-A)  5  225
#> 1914          Prkce             (B-A)  5  225
#> 1915        Zkscan3             (B-A)  5  225
#> 1916        Jakmip1             (B-A)  5  225
#> 1917        Tubgcp4             (B-A)  5  225
#> 1918           Pcnx             (B-A)  5  225
#> 1919        Dennd2c             (B-A)  5  225
#> 1920     D5Ertd579e             (B-A)  5  225
#> 1921        Rtn4rl1             (B-A)  5  225
#> 1922          Otub2             (B-A)  5  225
#> 1923        Trmt61b             (B-A)  5  225
#> 1924            Pkm             (B-C)  6  379
#> 1925          Prdx6             (B-C)  6  379
#> 1926         Pabpc1             (B-C)  6  379
#> 1927          Cmtm7             (B-C)  6  379
#> 1928          Ncoa4             (B-C)  6  379
#> 1929          Farsa             (B-C)  6  379
#> 1930        Slco3a1             (B-C)  6  379
#> 1931          Abce1             (B-C)  6  379
#> 1932            Vim             (B-C)  6  379
#> 1933          Irgm1             (B-C)  6  379
#> 1934          Papd4             (B-C)  6  379
#> 1935           Mfn2             (B-C)  6  379
#> 1936        Slc35c2             (B-C)  6  379
#> 1937          Mapk9             (B-C)  6  379
#> 1938           Zhx1             (B-C)  6  379
#> 1939            Grn             (B-C)  6  379
#> 1940           Aff3             (B-C)  6  379
#> 1941         Tm7sf3             (B-C)  6  379
#> 1942       AI467606             (B-C)  6  379
#> 1943           Pid1             (B-C)  6  379
#> 1944          Btaf1             (B-C)  6  379
#> 1945          Rmdn2             (B-C)  6  379
#> 1946          Runx1             (B-C)  6  379
#> 1947          Bbs12             (B-C)  6  379
#> 1948             Qk             (B-C)  6  379
#> 1949         Zfp148             (B-C)  6  379
#> 1950          Fmnl3             (B-C)  6  379
#> 1951            Fyn             (B-C)  6  379
#> 1952         Tardbp             (B-C)  6  379
#> 1953      Tnfaip8l1             (B-C)  6  379
#> 1954           Cbx7             (B-C)  6  379
#> 1955           Tbl3             (B-C)  6  379
#> 1956         Map2k1             (B-C)  6  379
#> 1957          Ap1g1             (B-C)  6  379
#> 1958            Lyn             (B-C)  6  379
#> 1959       Leprotl1             (B-C)  6  379
#> 1960          Anxa1             (B-C)  6  379
#> 1961        Tmem177             (B-C)  6  379
#> 1962        Zdhhc18             (B-C)  6  379
#> 1963        Snrnp70             (B-C)  6  379
#> 1964         Zfp128             (B-C)  6  379
#> 1965          Lyrm1             (B-C)  6  379
#> 1966         Rasl12             (B-C)  6  379
#> 1967       Tbc1d10c             (B-C)  6  379
#> 1968            Fes             (B-C)  6  379
#> 1969         Cited2             (B-C)  6  379
#> 1970         Cnot6l             (B-C)  6  379
#> 1971           Tle6             (B-C)  6  379
#> 1972          Uhrf2             (B-C)  6  379
#> 1973          Hemk1             (B-C)  6  379
#> 1974            Fn1             (B-C)  6  379
#> 1975           Tti1             (B-C)  6  379
#> 1976         Trim46             (B-C)  6  379
#> 1977  2610008E11Rik             (B-C)  6  379
#> 1978         Nipal3             (B-C)  6  379
#> 1979  1300002E11Rik             (B-C)  6  379
#> 1980            Hk1             (B-C)  6  379
#> 1981        Arl14ep             (B-C)  6  379
#> 1982         Nfkbia             (B-C)  6  379
#> 1983           Pkn1             (B-C)  6  379
#> 1984       Cdc42ep3             (B-C)  6  379
#> 1985           Ctsh             (B-C)  6  379
#> 1986         Il20ra             (B-C)  6  379
#> 1987         Zbtb37             (B-C)  6  379
#> 1988           Ldb1             (B-C)  6  379
#> 1989           Ralb             (B-C)  6  379
#> 1990           Tep1             (B-C)  6  379
#> 1991        Phf20l1             (B-C)  6  379
#> 1992         Pycard             (B-C)  6  379
#> 1993          Sbno1             (B-C)  6  379
#> 1994           Hars             (B-C)  6  379
#> 1995          Arrb1             (B-C)  6  379
#> 1996          Nudt4             (B-C)  6  379
#> 1997         Gtpbp3             (B-C)  6  379
#> 1998           Snrk             (B-C)  6  379
#> 1999          Mrps7             (B-C)  6  379
#> 2000           Gusb             (B-C)  6  379
#> 2001         Phf21a             (B-C)  6  379
#> 2002           Pmm2             (B-C)  6  379
#> 2003         Zfp182             (B-C)  6  379
#> 2004          Dgcr2             (B-C)  6  379
#> 2005          Pus7l             (B-C)  6  379
#> 2006         Gm7694             (B-C)  6  379
#> 2007         Rnf217             (B-C)  6  379
#> 2008  B130055M24Rik             (B-C)  6  379
#> 2009          Cops2             (B-C)  6  379
#> 2010           Csf1             (B-C)  6  379
#> 2011         Rabac1             (B-C)  6  379
#> 2012       Trappc13             (B-C)  6  379
#> 2013          Enpp5             (B-C)  6  379
#> 2014           Gps2             (B-C)  6  379
#> 2015        Poglut1             (B-C)  6  379
#> 2016          Mocs1             (B-C)  6  379
#> 2017          Pde4b             (B-C)  6  379
#> 2018           Bin3             (B-C)  6  379
#> 2019          Fbxl5             (B-C)  6  379
#> 2020        Tnfsf14             (B-C)  6  379
#> 2021          Ccnl2             (B-C)  6  379
#> 2022         Heatr1             (B-C)  6  379
#> 2023           Crem             (B-C)  6  379
#> 2024         Map4k2             (B-C)  6  379
#> 2025         Fam49b             (B-C)  6  379
#> 2026      Hist1h2bj             (B-C)  6  379
#> 2027         Dcaf13             (B-C)  6  379
#> 2028            Mvd             (B-C)  6  379
#> 2029        Rtn4ip1             (B-C)  6  379
#> 2030         Rnf122             (B-C)  6  379
#> 2031          Hspe1             (B-C)  6  379
#> 2032         Map3k6             (B-C)  6  379
#> 2033          Efna4             (B-C)  6  379
#> 2034           Icmt             (B-C)  6  379
#> 2035          Tor4a             (B-C)  6  379
#> 2036         Trafd1             (B-C)  6  379
#> 2037         Inafm1             (B-C)  6  379
#> 2038         Ccrn4l             (B-C)  6  379
#> 2039            Eng             (B-C)  6  379
#> 2040        Zscan20             (B-C)  6  379
#> 2041         Rab27a             (B-C)  6  379
#> 2042        Angptl4             (B-C)  6  379
#> 2043         Kbtbd8             (B-C)  6  379
#> 2044           Ier3             (B-C)  6  379
#> 2045          Nr3c1             (B-C)  6  379
#> 2046           Gfm1             (B-C)  6  379
#> 2047           Ctsz             (B-C)  6  379
#> 2048           Aqp1             (B-C)  6  379
#> 2049          Casc3             (B-C)  6  379
#> 2050        Heatr5b             (B-C)  6  379
#> 2051         Klhl42             (B-C)  6  379
#> 2052          Senp2             (B-C)  6  379
#> 2053           Vwa8             (B-C)  6  379
#> 2054          Lsm11             (B-C)  6  379
#> 2055         Gimap6             (B-C)  6  379
#> 2056         Tuba1c             (B-C)  6  379
#> 2057          Sesn1             (B-C)  6  379
#> 2058           Mios             (B-C)  6  379
#> 2059          Abcg3             (B-C)  6  379
#> 2060        Fam173b             (B-C)  6  379
#> 2061          Etfdh             (B-C)  6  379
#> 2062          Myo10             (B-C)  6  379
#> 2063          Ppm1m             (B-C)  6  379
#> 2064       Arhgap25             (B-C)  6  379
#> 2065         Phldb1             (B-C)  6  379
#> 2066       Epb4.1l1             (B-C)  6  379
#> 2067       BC033916             (B-C)  6  379
#> 2068         Atp1a3             (B-C)  6  379
#> 2069           Pisd             (B-C)  6  379
#> 2070          Tmcc2             (B-C)  6  379
#> 2071          Dzip3             (B-C)  6  379
#> 2072           Amz2             (B-C)  6  379
#> 2073           Ssr3             (B-C)  6  379
#> 2074          Egln3             (B-C)  6  379
#> 2075          Ribc1             (B-C)  6  379
#> 2076          Rnf44             (B-C)  6  379
#> 2077           Numb             (B-C)  6  379
#> 2078        Akirin2             (B-C)  6  379
#> 2079          Cebpa             (B-C)  6  379
#> 2080         Zfp617             (B-C)  6  379
#> 2081          Rcor2             (B-C)  6  379
#> 2082  3830403N18Rik             (B-C)  6  379
#> 2083          Nr2c1             (B-C)  6  379
#> 2084        Aldh3a2             (B-C)  6  379
#> 2085           Zzz3             (B-C)  6  379
#> 2086           Tab3             (B-C)  6  379
#> 2087           Smc5             (B-C)  6  379
#> 2088           Pak1             (B-C)  6  379
#> 2089         Srd5a1             (B-C)  6  379
#> 2090          Fancf             (B-C)  6  379
#> 2091        Adprhl2             (B-C)  6  379
#> 2092           Otos             (B-C)  6  379
#> 2093           Gale             (B-C)  6  379
#> 2094           Lrba             (B-C)  6  379
#> 2095         Tbc1d7             (B-C)  6  379
#> 2096         Lgals8             (B-C)  6  379
#> 2097         Zcchc2             (B-C)  6  379
#> 2098          Ift88             (B-C)  6  379
#> 2099         Pgm2l1             (B-C)  6  379
#> 2100         Tespa1             (B-C)  6  379
#> 2101          Fbxw4             (B-C)  6  379
#> 2102          Large             (B-C)  6  379
#> 2103             Xk             (B-C)  6  379
#> 2104            Nmi             (B-C)  6  379
#> 2105        Mettl13             (B-C)  6  379
#> 2106         Mrpl57             (B-C)  6  379
#> 2107          Mipep             (B-C)  6  379
#> 2108  D130020L05Rik             (B-C)  6  379
#> 2109          Mlycd             (B-C)  6  379
#> 2110         Zbtb34             (B-C)  6  379
#> 2111          Foxo1             (B-C)  6  379
#> 2112        Dync1h1             (B-C)  6  379
#> 2113          Ehmt2             (B-C)  6  379
#> 2114         Gtf3c2             (B-C)  6  379
#> 2115          Actr6             (B-C)  6  379
#> 2116           Ldlr             (B-C)  6  379
#> 2117        Rapgef6             (B-C)  6  379
#> 2118            Aen             (B-C)  6  379
#> 2119          Slfn5             (B-C)  6  379
#> 2120          Zbed4             (B-C)  6  379
#> 2121      Eif4enif1             (B-C)  6  379
#> 2122        Arhgap4             (B-C)  6  379
#> 2123        Slc44a1             (B-C)  6  379
#> 2124          Daam1             (B-C)  6  379
#> 2125          Lcorl             (B-C)  6  379
#> 2126        Tbc1d31             (B-C)  6  379
#> 2127           Per3             (B-C)  6  379
#> 2128         Rfxank             (B-C)  6  379
#> 2129        Zfyve26             (B-C)  6  379
#> 2130          Ier5l             (B-C)  6  379
#> 2131          Copz2             (B-C)  6  379
#> 2132        Map3k14             (B-C)  6  379
#> 2133          Kcnh2             (B-C)  6  379
#> 2134        Proser3             (B-C)  6  379
#> 2135         Scamp1             (B-C)  6  379
#> 2136        Zfp955b             (B-C)  6  379
#> 2137         Trim45             (B-C)  6  379
#> 2138           Ccnh             (B-C)  6  379
#> 2139           Coa4             (B-C)  6  379
#> 2140          P4ha1             (B-C)  6  379
#> 2141          Ints4             (B-C)  6  379
#> 2142           Nbas             (B-C)  6  379
#> 2143         Gm4070             (B-C)  6  379
#> 2144  D030056L22Rik             (B-C)  6  379
#> 2145           Sik1             (B-C)  6  379
#> 2146          Nlrx1             (B-C)  6  379
#> 2147           Dok1             (B-C)  6  379
#> 2148          Cd274             (B-C)  6  379
#> 2149         Eif2b5             (B-C)  6  379
#> 2150          Nsun6             (B-C)  6  379
#> 2151        Prkar2b             (B-C)  6  379
#> 2152           Rras             (B-C)  6  379
#> 2153        Clec16a             (B-C)  6  379
#> 2154          Taf4a             (B-C)  6  379
#> 2155         Slain1             (B-C)  6  379
#> 2156        Dync2h1             (B-C)  6  379
#> 2157         Mllt10             (B-C)  6  379
#> 2158          Ift80             (B-C)  6  379
#> 2159         Tm7sf2             (B-C)  6  379
#> 2160         Neurl3             (B-C)  6  379
#> 2161           Aff1             (B-C)  6  379
#> 2162        Fam120b             (B-C)  6  379
#> 2163          Pear1             (B-C)  6  379
#> 2164         Smim24             (B-C)  6  379
#> 2165          Uckl1             (B-C)  6  379
#> 2166       Mapkapk3             (B-C)  6  379
#> 2167       Vkorc1l1             (B-C)  6  379
#> 2168           Ltbr             (B-C)  6  379
#> 2169         Cluap1             (B-C)  6  379
#> 2170         Zswim7             (B-C)  6  379
#> 2171         Mrps27             (B-C)  6  379
#> 2172  6030458C11Rik             (B-C)  6  379
#> 2173        Slc38a7             (B-C)  6  379
#> 2174            Lxn             (B-C)  6  379
#> 2175          Clcn3             (B-C)  6  379
#> 2176         Pik3r1             (B-C)  6  379
#> 2177          Cd247             (B-C)  6  379
#> 2178         Smurf2             (B-C)  6  379
#> 2179         Zfp839             (B-C)  6  379
#> 2180        Slc10a7             (B-C)  6  379
#> 2181          Phf13             (B-C)  6  379
#> 2182        Slc31a2             (B-C)  6  379
#> 2183  1600002H07Rik             (B-C)  6  379
#> 2184         Nudt18             (B-C)  6  379
#> 2185          Msto1             (B-C)  6  379
#> 2186         Plxna3             (B-C)  6  379
#> 2187       Slc25a27             (B-C)  6  379
#> 2188          Dnpep             (B-C)  6  379
#> 2189         Dgcr14             (B-C)  6  379
#> 2190     Csgalnact2             (B-C)  6  379
#> 2191          Pomt2             (B-C)  6  379
#> 2192           Emc8             (B-C)  6  379
#> 2193           Mmab             (B-C)  6  379
#> 2194           Rem2             (B-C)  6  379
#> 2195          Ocel1             (B-C)  6  379
#> 2196         Pea15a             (B-C)  6  379
#> 2197          Pde6d             (B-C)  6  379
#> 2198        Cwf19l1             (B-C)  6  379
#> 2199          Tmco6             (B-C)  6  379
#> 2200         Fbxo42             (B-C)  6  379
#> 2201       Zkscan14             (B-C)  6  379
#> 2202          Xylt1             (B-C)  6  379
#> 2203          Lrch1             (B-C)  6  379
#> 2204         Mboat2             (B-C)  6  379
#> 2205        Tmem192             (B-C)  6  379
#> 2206          Phka2             (B-C)  6  379
#> 2207          Adck3             (B-C)  6  379
#> 2208          Cenpl             (B-C)  6  379
#> 2209          Prdm5             (B-C)  6  379
#> 2210          Tnip3             (B-C)  6  379
#> 2211           Dlg1             (B-C)  6  379
#> 2212         Dolpp1             (B-C)  6  379
#> 2213          Ulbp1             (B-C)  6  379
#> 2214          Wdpcp             (B-C)  6  379
#> 2215         Polr3a             (B-C)  6  379
#> 2216          Ints3             (B-C)  6  379
#> 2217           Aco1             (B-C)  6  379
#> 2218          Ptplb             (B-C)  6  379
#> 2219          Naa16             (B-C)  6  379
#> 2220          Rbm15             (B-C)  6  379
#> 2221        Fastkd1             (B-C)  6  379
#> 2222         Fam13b             (B-C)  6  379
#> 2223          Fbln1             (B-C)  6  379
#> 2224          Nudt8             (B-C)  6  379
#> 2225         Zbtb18             (B-C)  6  379
#> 2226         Gtf2e1             (B-C)  6  379
#> 2227          Clybl             (B-C)  6  379
#> 2228           Crat             (B-C)  6  379
#> 2229           Msh3             (B-C)  6  379
#> 2230       B3galnt2             (B-C)  6  379
#> 2231         Uap1l1             (B-C)  6  379
#> 2232       Pisd-ps1             (B-C)  6  379
#> 2233           Bms1             (B-C)  6  379
#> 2234            Met             (B-C)  6  379
#> 2235          Mfsd9             (B-C)  6  379
#> 2236          Aamdc             (B-C)  6  379
#> 2237        Depdc1b             (B-C)  6  379
#> 2238           Wwox             (B-C)  6  379
#> 2239           Sgcb             (B-C)  6  379
#> 2240          Tomm5             (B-C)  6  379
#> 2241        Akr1c13             (B-C)  6  379
#> 2242           Rhag             (B-C)  6  379
#> 2243        Prkar2a             (B-C)  6  379
#> 2244         Nfkbid             (B-C)  6  379
#> 2245           Shq1             (B-C)  6  379
#> 2246           Idua             (B-C)  6  379
#> 2247           Igtp             (B-C)  6  379
#> 2248            Mmd             (B-C)  6  379
#> 2249           Tmx4             (B-C)  6  379
#> 2250           Islr             (B-C)  6  379
#> 2251          Cdk10             (B-C)  6  379
#> 2252       Slc4a1ap             (B-C)  6  379
#> 2253        Ccdc122             (B-C)  6  379
#> 2254          Frat2             (B-C)  6  379
#> 2255         Metrnl             (B-C)  6  379
#> 2256          Mex3c             (B-C)  6  379
#> 2257         Angpt1             (B-C)  6  379
#> 2258           Nck2             (B-C)  6  379
#> 2259            Htt             (B-C)  6  379
#> 2260          Nat10             (B-C)  6  379
#> 2261         Zfp219             (B-C)  6  379
#> 2262        Zkscan5             (B-C)  6  379
#> 2263        Zfp322a             (B-C)  6  379
#> 2264          Stk10             (B-C)  6  379
#> 2265           Lap3             (B-C)  6  379
#> 2266         Tbc1d2             (B-C)  6  379
#> 2267           Dgkd             (B-C)  6  379
#> 2268          Txlna             (B-C)  6  379
#> 2269  D430042O09Rik             (B-C)  6  379
#> 2270          Rras2             (B-C)  6  379
#> 2271          Cspp1             (B-C)  6  379
#> 2272       AW554918             (B-C)  6  379
#> 2273          Ap1g2             (B-C)  6  379
#> 2274        Pitpnm1             (B-C)  6  379
#> 2275           Taf1             (B-C)  6  379
#> 2276          Dapk1             (B-C)  6  379
#> 2277          Smim3             (B-C)  6  379
#> 2278           Atl3             (B-C)  6  379
#> 2279       Mapk8ip3             (B-C)  6  379
#> 2280        N4bp2l1             (B-C)  6  379
#> 2281       BC002059             (B-C)  6  379
#> 2282         Oxnad1             (B-C)  6  379
#> 2283         Cdc25b             (B-C)  6  379
#> 2284          Pot1b             (B-C)  6  379
#> 2285           Cdk8             (B-C)  6  379
#> 2286          Ap4e1             (B-C)  6  379
#> 2287           Ccl3             (B-C)  6  379
#> 2288        Gnpnat1             (B-C)  6  379
#> 2289          Zbtb3             (B-C)  6  379
#> 2290          Zfp59             (B-C)  6  379
#> 2291           Twf1             (B-C)  6  379
#> 2292          Kctd9             (B-C)  6  379
#> 2293           Ece1             (B-C)  6  379
#> 2294         Zfp738             (B-C)  6  379
#> 2295           Gnb4             (B-C)  6  379
#> 2296       AA987161             (B-C)  6  379
#> 2297         Zbtb22             (B-C)  6  379
#> 2298         Heatr2             (B-C)  6  379
#> 2299          Senp7             (B-C)  6  379
#> 2300            Tk2             (B-C)  6  379
#> 2301          Jade3             (B-C)  6  379
#> 2302       Cdc42ep4             (B-C)  6  379
#> 2303          Stmn1             (C-A)  7  766
#> 2304           Mcm6             (C-A)  7  766
#> 2305         Zfp706             (C-A)  7  766
#> 2306          Actr3             (C-A)  7  766
#> 2307            Msn             (C-A)  7  766
#> 2308          Srsf7             (C-A)  7  766
#> 2309          Paics             (C-A)  7  766
#> 2310          Srsf3             (C-A)  7  766
#> 2311           Cd47             (C-A)  7  766
#> 2312          Hspa9             (C-A)  7  766
#> 2313           Ly6e             (C-A)  7  766
#> 2314          Nup62             (C-A)  7  766
#> 2315          Aldoa             (C-A)  7  766
#> 2316          Hmgn2             (C-A)  7  766
#> 2317           Flna             (C-A)  7  766
#> 2318          Dnmt1             (C-A)  7  766
#> 2319          Cdca7             (C-A)  7  766
#> 2320           Hn1l             (C-A)  7  766
#> 2321         Anapc5             (C-A)  7  766
#> 2322           Akna             (C-A)  7  766
#> 2323          Phgdh             (C-A)  7  766
#> 2324          Idh3g             (C-A)  7  766
#> 2325          Naa50             (C-A)  7  766
#> 2326          Myo1g             (C-A)  7  766
#> 2327            Sla             (C-A)  7  766
#> 2328          Cops4             (C-A)  7  766
#> 2329           Sdhd             (C-A)  7  766
#> 2330         Nfatc3             (C-A)  7  766
#> 2331       Tmem229b             (C-A)  7  766
#> 2332           Npm1             (C-A)  7  766
#> 2333           Cdk4             (C-A)  7  766
#> 2334           Rrm1             (C-A)  7  766
#> 2335        Fam111a             (C-A)  7  766
#> 2336          Hmha1             (C-A)  7  766
#> 2337           Mcm2             (C-A)  7  766
#> 2338            Fh1             (C-A)  7  766
#> 2339          Rpl32             (C-A)  7  766
#> 2340          Ikzf1             (C-A)  7  766
#> 2341          Ap3d1             (C-A)  7  766
#> 2342        Erbb2ip             (C-A)  7  766
#> 2343           Pten             (C-A)  7  766
#> 2344         Topbp1             (C-A)  7  766
#> 2345            Ran             (C-A)  7  766
#> 2346         Zfp422             (C-A)  7  766
#> 2347          Gpr25             (C-A)  7  766
#> 2348            Erh             (C-A)  7  766
#> 2349         Rps27a             (C-A)  7  766
#> 2350  2310034G01Rik             (C-A)  7  766
#> 2351         Clint1             (C-A)  7  766
#> 2352          Hcls1             (C-A)  7  766
#> 2353          Cstf2             (C-A)  7  766
#> 2354           Mcm3             (C-A)  7  766
#> 2355           Gart             (C-A)  7  766
#> 2356           Hat1             (C-A)  7  766
#> 2357          Psmb9             (C-A)  7  766
#> 2358          Gnai3             (C-A)  7  766
#> 2359         Nhp2l1             (C-A)  7  766
#> 2360          Rps17             (C-A)  7  766
#> 2361        Unc119b             (C-A)  7  766
#> 2362          Glrx3             (C-A)  7  766
#> 2363          Tra2b             (C-A)  7  766
#> 2364         Cyfip1             (C-A)  7  766
#> 2365        Tmem164             (C-A)  7  766
#> 2366         Ruvbl2             (C-A)  7  766
#> 2367          Rps25             (C-A)  7  766
#> 2368         Cited4             (C-A)  7  766
#> 2369           Pcna             (C-A)  7  766
#> 2370         Il31ra             (C-A)  7  766
#> 2371          Mknk2             (C-A)  7  766
#> 2372         Exosc8             (C-A)  7  766
#> 2373          Il21r             (C-A)  7  766
#> 2374         Lpcat3             (C-A)  7  766
#> 2375          Ccng1             (C-A)  7  766
#> 2376         Tuba4a             (C-A)  7  766
#> 2377          Cdca5             (C-A)  7  766
#> 2378          Trgv2             (C-A)  7  766
#> 2379         Tmem65             (C-A)  7  766
#> 2380          Cisd3             (C-A)  7  766
#> 2381         Nt5dc2             (C-A)  7  766
#> 2382            Nnt             (C-A)  7  766
#> 2383         Impdh1             (C-A)  7  766
#> 2384         Gm5601             (C-A)  7  766
#> 2385          Ccnb1             (C-A)  7  766
#> 2386         Sqstm1             (C-A)  7  766
#> 2387         Plscr3             (C-A)  7  766
#> 2388          Srpk1             (C-A)  7  766
#> 2389           H1f0             (C-A)  7  766
#> 2390          Yipf3             (C-A)  7  766
#> 2391          Wasf2             (C-A)  7  766
#> 2392          Sstr2             (C-A)  7  766
#> 2393           Rnf6             (C-A)  7  766
#> 2394           Phb2             (C-A)  7  766
#> 2395          Kif17             (C-A)  7  766
#> 2396          Eno1b             (C-A)  7  766
#> 2397         Bnip3l             (C-A)  7  766
#> 2398          Nop56             (C-A)  7  766
#> 2399         Spink4             (C-A)  7  766
#> 2400           Cyba             (C-A)  7  766
#> 2401         Dnajb1             (C-A)  7  766
#> 2402            Mdk             (C-A)  7  766
#> 2403         Snrpd1             (C-A)  7  766
#> 2404          Aurkb             (C-A)  7  766
#> 2405          Morc3             (C-A)  7  766
#> 2406          Mturn             (C-A)  7  766
#> 2407           Rif1             (C-A)  7  766
#> 2408        Gm13910             (C-A)  7  766
#> 2409          mt-Tp             (C-A)  7  766
#> 2410            Sp4             (C-A)  7  766
#> 2411         Syngr2             (C-A)  7  766
#> 2412          Trap1             (C-A)  7  766
#> 2413        Arfgef2             (C-A)  7  766
#> 2414          Sf3a1             (C-A)  7  766
#> 2415          Exoc6             (C-A)  7  766
#> 2416          Ube2h             (C-A)  7  766
#> 2417          Ypel3             (C-A)  7  766
#> 2418        Tbc1d12             (C-A)  7  766
#> 2419          Mgea5             (C-A)  7  766
#> 2420          Ncoa3             (C-A)  7  766
#> 2421           Tjp3             (C-A)  7  766
#> 2422          Pola2             (C-A)  7  766
#> 2423           Tyms             (C-A)  7  766
#> 2424           Paox             (C-A)  7  766
#> 2425         H2-Ab1             (C-A)  7  766
#> 2426         Gm5560             (C-A)  7  766
#> 2427         Phf11c             (C-A)  7  766
#> 2428         Vps37b             (C-A)  7  766
#> 2429          Dfna5             (C-A)  7  766
#> 2430       AW112010             (C-A)  7  766
#> 2431         Gm4735             (C-A)  7  766
#> 2432        Gm12696             (C-A)  7  766
#> 2433          Abcb8             (C-A)  7  766
#> 2434         Ccdc50             (C-A)  7  766
#> 2435  6330419J24Rik             (C-A)  7  766
#> 2436       Scp2-ps2             (C-A)  7  766
#> 2437           Umps             (C-A)  7  766
#> 2438           Evpl             (C-A)  7  766
#> 2439         Gm7353             (C-A)  7  766
#> 2440          Prex1             (C-A)  7  766
#> 2441           Abi1             (C-A)  7  766
#> 2442        Gm14176             (C-A)  7  766
#> 2443          Etnk1             (C-A)  7  766
#> 2444        Gm26804             (C-A)  7  766
#> 2445        Gm10480             (C-A)  7  766
#> 2446          Numa1             (C-A)  7  766
#> 2447         Myo18a             (C-A)  7  766
#> 2448         Gm5499             (C-A)  7  766
#> 2449           Pnck             (C-A)  7  766
#> 2450         Gm6055             (C-A)  7  766
#> 2451          Yipf6             (C-A)  7  766
#> 2452        Gm10241             (C-A)  7  766
#> 2453          Reps1             (C-A)  7  766
#> 2454         Zfand6             (C-A)  7  766
#> 2455          Wdr12             (C-A)  7  766
#> 2456           Cdk6             (C-A)  7  766
#> 2457        Plekho2             (C-A)  7  766
#> 2458          Kif24             (C-A)  7  766
#> 2459         Inpp5f             (C-A)  7  766
#> 2460        Trmt112             (C-A)  7  766
#> 2461           Cst7             (C-A)  7  766
#> 2462          Osgep             (C-A)  7  766
#> 2463         Gm1840             (C-A)  7  766
#> 2464          Jade1             (C-A)  7  766
#> 2465         Pdlim7             (C-A)  7  766
#> 2466         Ndufa9             (C-A)  7  766
#> 2467          Sel1l             (C-A)  7  766
#> 2468           Suox             (C-A)  7  766
#> 2469         Stoml2             (C-A)  7  766
#> 2470          Tmed4             (C-A)  7  766
#> 2471         Zfp263             (C-A)  7  766
#> 2472         Vcpkmt             (C-A)  7  766
#> 2473         Setdb2             (C-A)  7  766
#> 2474         Gm8355             (C-A)  7  766
#> 2475           Cd86             (C-A)  7  766
#> 2476         Gm8399             (C-A)  7  766
#> 2477        Gm15542             (C-A)  7  766
#> 2478            Ftx             (C-A)  7  766
#> 2479          Il4ra             (C-A)  7  766
#> 2480        Gm18588             (C-A)  7  766
#> 2481         Sppl2b             (C-A)  7  766
#> 2482         Tbc1d8             (C-A)  7  766
#> 2483           Pigs             (C-A)  7  766
#> 2484         Chrnb1             (C-A)  7  766
#> 2485          Alox5             (C-A)  7  766
#> 2486         Fam60a             (C-A)  7  766
#> 2487          Casp9             (C-A)  7  766
#> 2488        Gm12183             (C-A)  7  766
#> 2489        Gm21399             (C-A)  7  766
#> 2490        Gm10093             (C-A)  7  766
#> 2491  1700056E22Rik             (C-A)  7  766
#> 2492        Krtcap3             (C-A)  7  766
#> 2493  5430427O19Rik             (C-A)  7  766
#> 2494     D7Ertd715e             (C-A)  7  766
#> 2495        Gm10524             (C-A)  7  766
#> 2496         Nanos3             (C-A)  7  766
#> 2497           Nsl1             (C-A)  7  766
#> 2498           Poll             (C-A)  7  766
#> 2499         Tgoln1             (C-A)  7  766
#> 2500        Rsl24d1             (C-A)  7  766
#> 2501          Usp24             (C-A)  7  766
#> 2502         Smim13             (C-A)  7  766
#> 2503  1110008F13Rik             (C-A)  7  766
#> 2504        Akirin1             (C-A)  7  766
#> 2505           Pkig             (C-A)  7  766
#> 2506           Thbd             (C-A)  7  766
#> 2507       BC055324             (C-A)  7  766
#> 2508           Ccz1             (C-A)  7  766
#> 2509      Vdac3-ps1             (C-A)  7  766
#> 2510         Gm5913             (C-A)  7  766
#> 2511        Prelid2             (C-A)  7  766
#> 2512        Sipa1l1             (C-A)  7  766
#> 2513          Eps15             (C-A)  7  766
#> 2514         Anp32b             (C-A)  7  766
#> 2515           Phip             (C-A)  7  766
#> 2516          H2-Q7             (C-A)  7  766
#> 2517         Tmem98             (C-A)  7  766
#> 2518          Btbd2             (C-A)  7  766
#> 2519          Sytl1             (C-A)  7  766
#> 2520          Ncor1             (C-A)  7  766
#> 2521          Abhd6             (C-A)  7  766
#> 2522          Ikbke             (C-A)  7  766
#> 2523        Slc31a1             (C-A)  7  766
#> 2524          Rgs10             (C-A)  7  766
#> 2525          Nabp1             (C-A)  7  766
#> 2526  C330027C09Rik             (C-A)  7  766
#> 2527           Mtf2             (C-A)  7  766
#> 2528          Ndrg1             (C-A)  7  766
#> 2529          Wdr11             (C-A)  7  766
#> 2530           Unkl             (C-A)  7  766
#> 2531          Smap1             (C-A)  7  766
#> 2532         Crebrf             (C-A)  7  766
#> 2533           Nasp             (C-A)  7  766
#> 2534           Gbp3             (C-A)  7  766
#> 2535          Creb1             (C-A)  7  766
#> 2536          Zzef1             (C-A)  7  766
#> 2537           Cbx1             (C-A)  7  766
#> 2538          Uqcrq             (C-A)  7  766
#> 2539           Rfc5             (C-A)  7  766
#> 2540         Gm4149             (C-A)  7  766
#> 2541         Gemin6             (C-A)  7  766
#> 2542        Gm26781             (C-A)  7  766
#> 2543           Mta1             (C-A)  7  766
#> 2544       Ifi27l2a             (C-A)  7  766
#> 2545        Tmem201             (C-A)  7  766
#> 2546          Tcam1             (C-A)  7  766
#> 2547         H2-T22             (C-A)  7  766
#> 2548        Gm12184             (C-A)  7  766
#> 2549          Tnip1             (C-A)  7  766
#> 2550          Snrpa             (C-A)  7  766
#> 2551        Gm24366             (C-A)  7  766
#> 2552        Gm11868             (C-A)  7  766
#> 2553       BC068157             (C-A)  7  766
#> 2554         Cuedc2             (C-A)  7  766
#> 2555          Dnph1             (C-A)  7  766
#> 2556  5830415F09Rik             (C-A)  7  766
#> 2557           Gas5             (C-A)  7  766
#> 2558        Gm15772             (C-A)  7  766
#> 2559          Erp44             (C-A)  7  766
#> 2560           Gcsh             (C-A)  7  766
#> 2561        Tnfaip3             (C-A)  7  766
#> 2562         Gm9755             (C-A)  7  766
#> 2563        Gm12657             (C-A)  7  766
#> 2564         Hectd1             (C-A)  7  766
#> 2565         Gm6382             (C-A)  7  766
#> 2566         Gm9855             (C-A)  7  766
#> 2567        Gm11223             (C-A)  7  766
#> 2568           Brd1             (C-A)  7  766
#> 2569        Gm15210             (C-A)  7  766
#> 2570          Taf9b             (C-A)  7  766
#> 2571           H6pd             (C-A)  7  766
#> 2572         Gm8325             (C-A)  7  766
#> 2573           Zxdb             (C-A)  7  766
#> 2574         Csf2rb             (C-A)  7  766
#> 2575        Gm13436             (C-A)  7  766
#> 2576         Zfp773             (C-A)  7  766
#> 2577          H2-Q6             (C-A)  7  766
#> 2578           Adnp             (C-A)  7  766
#> 2579          Hells             (C-A)  7  766
#> 2580        Morf4l2             (C-A)  7  766
#> 2581        St8sia4             (C-A)  7  766
#> 2582            Pnp             (C-A)  7  766
#> 2583        Gm12174             (C-A)  7  766
#> 2584      Pgam1-ps2             (C-A)  7  766
#> 2585        Gm15501             (C-A)  7  766
#> 2586          Asf1a             (C-A)  7  766
#> 2587            Ttl             (C-A)  7  766
#> 2588         Zfp777             (C-A)  7  766
#> 2589           Tet1             (C-A)  7  766
#> 2590         Gm9844             (C-A)  7  766
#> 2591         Angel1             (C-A)  7  766
#> 2592        Gm20605             (C-A)  7  766
#> 2593         Pcyt1a             (C-A)  7  766
#> 2594         Fbxl12             (C-A)  7  766
#> 2595          Pold3             (C-A)  7  766
#> 2596          Hoxa7             (C-A)  7  766
#> 2597       BC002163             (C-A)  7  766
#> 2598         Gm9625             (C-A)  7  766
#> 2599        Plekhm1             (C-A)  7  766
#> 2600         Ptger1             (C-A)  7  766
#> 2601         Dennd3             (C-A)  7  766
#> 2602        Gpatch3             (C-A)  7  766
#> 2603        Tmem30a             (C-A)  7  766
#> 2604        Gm13835             (C-A)  7  766
#> 2605      Secisbp2l             (C-A)  7  766
#> 2606           Ighd             (C-A)  7  766
#> 2607          Ypel2             (C-A)  7  766
#> 2608          Ncoa2             (C-A)  7  766
#> 2609         Fkbp1a             (C-A)  7  766
#> 2610         Gm6180             (C-A)  7  766
#> 2611          Btbd3             (C-A)  7  766
#> 2612         Trim37             (C-A)  7  766
#> 2613         Marcks             (C-A)  7  766
#> 2614           Rcn1             (C-A)  7  766
#> 2615           Dctd             (C-A)  7  766
#> 2616        Bloc1s2             (C-A)  7  766
#> 2617           Lifr             (C-A)  7  766
#> 2618         Gm7488             (C-A)  7  766
#> 2619           Add1             (C-A)  7  766
#> 2620  4930402H24Rik             (C-A)  7  766
#> 2621          Ftsj1             (C-A)  7  766
#> 2622           Myl4             (C-A)  7  766
#> 2623        Gm10053             (C-A)  7  766
#> 2624          Acot2             (C-A)  7  766
#> 2625         Gm3604             (C-A)  7  766
#> 2626           Slpi             (C-A)  7  766
#> 2627         Gm7816             (C-A)  7  766
#> 2628         Pbxip1             (C-A)  7  766
#> 2629            Hyi             (C-A)  7  766
#> 2630         Ticam1             (C-A)  7  766
#> 2631  9130020K20Rik             (C-A)  7  766
#> 2632          Desi1             (C-A)  7  766
#> 2633         Nanos1             (C-A)  7  766
#> 2634          Degs1             (C-A)  7  766
#> 2635        Gm12346             (C-A)  7  766
#> 2636          Dyrk3             (C-A)  7  766
#> 2637        Gm11560             (C-A)  7  766
#> 2638           Sil1             (C-A)  7  766
#> 2639        Gm14150             (C-A)  7  766
#> 2640         Gm9531             (C-A)  7  766
#> 2641       Tmem106c             (C-A)  7  766
#> 2642           Relt             (C-A)  7  766
#> 2643          Cks1b             (C-A)  7  766
#> 2644          Scn1b             (C-A)  7  766
#> 2645           Clpb             (C-A)  7  766
#> 2646          Reep4             (C-A)  7  766
#> 2647  A430005L14Rik             (C-A)  7  766
#> 2648           Narf             (C-A)  7  766
#> 2649        Slc2a10             (C-A)  7  766
#> 2650      Hmgb1-ps5             (C-A)  7  766
#> 2651          Gna11             (C-A)  7  766
#> 2652         Gm3531             (C-A)  7  766
#> 2653           Rpa2             (C-A)  7  766
#> 2654          Plod1             (C-A)  7  766
#> 2655          Tada1             (C-A)  7  766
#> 2656          Snhg4             (C-A)  7  766
#> 2657  6330549D23Rik             (C-A)  7  766
#> 2658          Kdm6b             (C-A)  7  766
#> 2659          Ift22             (C-A)  7  766
#> 2660          Cdan1             (C-A)  7  766
#> 2661         Pdlim1             (C-A)  7  766
#> 2662         Rps3a3             (C-A)  7  766
#> 2663          Kantr             (C-A)  7  766
#> 2664         Zbtb24             (C-A)  7  766
#> 2665        Gm12263             (C-A)  7  766
#> 2666           Vrk1             (C-A)  7  766
#> 2667          Pold2             (C-A)  7  766
#> 2668      Rps18-ps3             (C-A)  7  766
#> 2669           Neu1             (C-A)  7  766
#> 2670           Cd72             (C-A)  7  766
#> 2671         Zfp146             (C-A)  7  766
#> 2672          Rock2             (C-A)  7  766
#> 2673            Ski             (C-A)  7  766
#> 2674          Ssbp1             (C-A)  7  766
#> 2675  D930048G16Rik             (C-A)  7  766
#> 2676         Insig1             (C-A)  7  766
#> 2677          Cep63             (C-A)  7  766
#> 2678        Ppip5k2             (C-A)  7  766
#> 2679       Slc39a11             (C-A)  7  766
#> 2680        Rapgef3             (C-A)  7  766
#> 2681          Cdc45             (C-A)  7  766
#> 2682           Ipo5             (C-A)  7  766
#> 2683          Brca1             (C-A)  7  766
#> 2684           Nme2             (C-A)  7  766
#> 2685          Abtb1             (C-A)  7  766
#> 2686       Tmsb15b1             (C-A)  7  766
#> 2687           Klf9             (C-A)  7  766
#> 2688          Lin54             (C-A)  7  766
#> 2689          Spc24             (C-A)  7  766
#> 2690        Gm12669             (C-A)  7  766
#> 2691          Hdac9             (C-A)  7  766
#> 2692        Gm10131             (C-A)  7  766
#> 2693      Kansl2-ps             (C-A)  7  766
#> 2694          Kdm4c             (C-A)  7  766
#> 2695          Sh2b1             (C-A)  7  766
#> 2696          Mcm10             (C-A)  7  766
#> 2697           Nfu1             (C-A)  7  766
#> 2698          Yif1b             (C-A)  7  766
#> 2699       Tpt1-ps3             (C-A)  7  766
#> 2700          Zfp12             (C-A)  7  766
#> 2701         Polr3h             (C-A)  7  766
#> 2702         Necap1             (C-A)  7  766
#> 2703        Gm10698             (C-A)  7  766
#> 2704        Ccdc117             (C-A)  7  766
#> 2705     Eif5al3-ps             (C-A)  7  766
#> 2706          Aplp2             (C-A)  7  766
#> 2707       Arhgap32             (C-A)  7  766
#> 2708         Ankra2             (C-A)  7  766
#> 2709         Gm9385             (C-A)  7  766
#> 2710          Cd163             (C-A)  7  766
#> 2711         Mrpl17             (C-A)  7  766
#> 2712          Stau1             (C-A)  7  766
#> 2713  1810022K09Rik             (C-A)  7  766
#> 2714         Zfp113             (C-A)  7  766
#> 2715  1700025G04Rik             (C-A)  7  766
#> 2716           Wdr5             (C-A)  7  766
#> 2717      Smt3h2-ps             (C-A)  7  766
#> 2718          Vps16             (C-A)  7  766
#> 2719          Senp6             (C-A)  7  766
#> 2720  2410004P03Rik             (C-A)  7  766
#> 2721          C2cd3             (C-A)  7  766
#> 2722          Gins1             (C-A)  7  766
#> 2723        Gm10073             (C-A)  7  766
#> 2724         Zfp330             (C-A)  7  766
#> 2725        Plekha5             (C-A)  7  766
#> 2726        Gm10146             (C-A)  7  766
#> 2727           Eny2             (C-A)  7  766
#> 2728         Zfp346             (C-A)  7  766
#> 2729        Gm12355             (C-A)  7  766
#> 2730          F2rl3             (C-A)  7  766
#> 2731      Gabarapl2             (C-A)  7  766
#> 2732         Bcap29             (C-A)  7  766
#> 2733         Ifnar2             (C-A)  7  766
#> 2734          Tiam1             (C-A)  7  766
#> 2735          Haus1             (C-A)  7  766
#> 2736         Parpbp             (C-A)  7  766
#> 2737         Ttll12             (C-A)  7  766
#> 2738          Strn4             (C-A)  7  766
#> 2739          Siva1             (C-A)  7  766
#> 2740           Rere             (C-A)  7  766
#> 2741         Tmem9b             (C-A)  7  766
#> 2742          Cenpm             (C-A)  7  766
#> 2743          Peg12             (C-A)  7  766
#> 2744          Itpkb             (C-A)  7  766
#> 2745        Kansl1l             (C-A)  7  766
#> 2746       Hist1h1c             (C-A)  7  766
#> 2747        Gm14586             (C-A)  7  766
#> 2748       AI504432             (C-A)  7  766
#> 2749  D330023K18Rik             (C-A)  7  766
#> 2750         Tsen34             (C-A)  7  766
#> 2751        Angptl2             (C-A)  7  766
#> 2752        Gm10269             (C-A)  7  766
#> 2753        Gm11478             (C-A)  7  766
#> 2754        Tmem121             (C-A)  7  766
#> 2755         Thap11             (C-A)  7  766
#> 2756        Slc43a1             (C-A)  7  766
#> 2757         Cpped1             (C-A)  7  766
#> 2758           Gcat             (C-A)  7  766
#> 2759        Fam188b             (C-A)  7  766
#> 2760        St3gal2             (C-A)  7  766
#> 2761         Fam63b             (C-A)  7  766
#> 2762         Rnf125             (C-A)  7  766
#> 2763         Ift122             (C-A)  7  766
#> 2764         Ikbkap             (C-A)  7  766
#> 2765        Mypopos             (C-A)  7  766
#> 2766        Gm26623             (C-A)  7  766
#> 2767         Camk2g             (C-A)  7  766
#> 2768         Ptdss2             (C-A)  7  766
#> 2769          Wdr62             (C-A)  7  766
#> 2770         Mms22l             (C-A)  7  766
#> 2771        Gm13493             (C-A)  7  766
#> 2772         Ppp6r2             (C-A)  7  766
#> 2773            Lat             (C-A)  7  766
#> 2774          Ntmt1             (C-A)  7  766
#> 2775          Prmt2             (C-A)  7  766
#> 2776        Tmem107             (C-A)  7  766
#> 2777         Gm7332             (C-A)  7  766
#> 2778          Chpt1             (C-A)  7  766
#> 2779      Rpl31-ps8             (C-A)  7  766
#> 2780        Trim12a             (C-A)  7  766
#> 2781          Cwc27             (C-A)  7  766
#> 2782         Heatr3             (C-A)  7  766
#> 2783         Ppp3cc             (C-A)  7  766
#> 2784        Gm13212             (C-A)  7  766
#> 2785          Ssbp2             (C-A)  7  766
#> 2786          Acad8             (C-A)  7  766
#> 2787        Gm17383             (C-A)  7  766
#> 2788        Gm24187             (C-A)  7  766
#> 2789           Prcp             (C-A)  7  766
#> 2790         Rps3a2             (C-A)  7  766
#> 2791            Lbh             (C-A)  7  766
#> 2792          Fads1             (C-A)  7  766
#> 2793          Anxa5             (C-A)  7  766
#> 2794         Cx3cr1             (C-A)  7  766
#> 2795        Gramd1c             (C-A)  7  766
#> 2796      Hmgb1-ps1             (C-A)  7  766
#> 2797        Cyp4f16             (C-A)  7  766
#> 2798           Ptms             (C-A)  7  766
#> 2799         Mfsd12             (C-A)  7  766
#> 2800        Gm21975             (C-A)  7  766
#> 2801       Suv420h1             (C-A)  7  766
#> 2802           Plk4             (C-A)  7  766
#> 2803  A830080D01Rik             (C-A)  7  766
#> 2804          Rchy1             (C-A)  7  766
#> 2805           Pkib             (C-A)  7  766
#> 2806          Syvn1             (C-A)  7  766
#> 2807            Dut             (C-A)  7  766
#> 2808          C1qbp             (C-A)  7  766
#> 2809          Rmdn1             (C-A)  7  766
#> 2810         Pik3r5             (C-A)  7  766
#> 2811            F8a             (C-A)  7  766
#> 2812           Sla2             (C-A)  7  766
#> 2813          Enkd1             (C-A)  7  766
#> 2814          Dctn4             (C-A)  7  766
#> 2815        Zcchc10             (C-A)  7  766
#> 2816         Inpp5a             (C-A)  7  766
#> 2817           Bub1             (C-A)  7  766
#> 2818        Hnrnph1             (C-A)  7  766
#> 2819         Katnb1             (C-A)  7  766
#> 2820          Wdhd1             (C-A)  7  766
#> 2821           Nqo2             (C-A)  7  766
#> 2822        Ccdc28b             (C-A)  7  766
#> 2823           Nop2             (C-A)  7  766
#> 2824          Fcgrt             (C-A)  7  766
#> 2825          Alg14             (C-A)  7  766
#> 2826        Ubash3b             (C-A)  7  766
#> 2827        Rad54l2             (C-A)  7  766
#> 2828          Meis1             (C-A)  7  766
#> 2829         Cgrrf1             (C-A)  7  766
#> 2830         Dusp23             (C-A)  7  766
#> 2831         Cops7a             (C-A)  7  766
#> 2832          Top2a             (C-A)  7  766
#> 2833  2810006K23Rik             (C-A)  7  766
#> 2834         Bckdha             (C-A)  7  766
#> 2835  2810408I11Rik             (C-A)  7  766
#> 2836          Mrpl9             (C-A)  7  766
#> 2837        Gm10443             (C-A)  7  766
#> 2838  D130058E05Rik             (C-A)  7  766
#> 2839          Hexdc             (C-A)  7  766
#> 2840  1700021K19Rik             (C-A)  7  766
#> 2841        Gm16039             (C-A)  7  766
#> 2842        Tsc22d2             (C-A)  7  766
#> 2843          Deaf1             (C-A)  7  766
#> 2844          Tulp4             (C-A)  7  766
#> 2845          Aph1c             (C-A)  7  766
#> 2846           Mtf1             (C-A)  7  766
#> 2847        Map3k10             (C-A)  7  766
#> 2848         Gm5835             (C-A)  7  766
#> 2849          Appl1             (C-A)  7  766
#> 2850          Slfn3             (C-A)  7  766
#> 2851           Skp2             (C-A)  7  766
#> 2852         C2cd2l             (C-A)  7  766
#> 2853    Tmem181b-ps             (C-A)  7  766
#> 2854           Bptf             (C-A)  7  766
#> 2855        Gm14303             (C-A)  7  766
#> 2856          Rbms1             (C-A)  7  766
#> 2857         Nufip1             (C-A)  7  766
#> 2858        Gm24270             (C-A)  7  766
#> 2859           Rab9             (C-A)  7  766
#> 2860       Gpatch11             (C-A)  7  766
#> 2861          Smyd5             (C-A)  7  766
#> 2862        Ppp1r11             (C-A)  7  766
#> 2863          Maml3             (C-A)  7  766
#> 2864          Wdr31             (C-A)  7  766
#> 2865          Galk1             (C-A)  7  766
#> 2866         Armcx4             (C-A)  7  766
#> 2867        Galnt11             (C-A)  7  766
#> 2868      Hmga1-rs1             (C-A)  7  766
#> 2869         Fbxo22             (C-A)  7  766
#> 2870      Fam120aos             (C-A)  7  766
#> 2871          Wdr44             (C-A)  7  766
#> 2872        Slc43a2             (C-A)  7  766
#> 2873          Zadh2             (C-A)  7  766
#> 2874         Gm6136             (C-A)  7  766
#> 2875          Usp34             (C-A)  7  766
#> 2876         Inpp4a             (C-A)  7  766
#> 2877          Ppm1a             (C-A)  7  766
#> 2878       Tpm3-rs7             (C-A)  7  766
#> 2879          Cxx1a             (C-A)  7  766
#> 2880         Pou2f2             (C-A)  7  766
#> 2881         Engase             (C-A)  7  766
#> 2882         Camk2a             (C-A)  7  766
#> 2883           Emc2             (C-A)  7  766
#> 2884          Dcaf5             (C-A)  7  766
#> 2885     Rpl17-ps10             (C-A)  7  766
#> 2886       Phospho2             (C-A)  7  766
#> 2887      Hist1h2ae             (C-A)  7  766
#> 2888          Pex13             (C-A)  7  766
#> 2889           Ska3             (C-A)  7  766
#> 2890         Gm8430             (C-A)  7  766
#> 2891          Fbxo5             (C-A)  7  766
#> 2892         Zfand3             (C-A)  7  766
#> 2893            Mr1             (C-A)  7  766
#> 2894           Sv2a             (C-A)  7  766
#> 2895           Xrn1             (C-A)  7  766
#> 2896           Junb             (C-A)  7  766
#> 2897       Pcna-ps2             (C-A)  7  766
#> 2898         Zc3h18             (C-A)  7  766
#> 2899         Card11             (C-A)  7  766
#> 2900          Dhodh             (C-A)  7  766
#> 2901           Epc2             (C-A)  7  766
#> 2902          Spag5             (C-A)  7  766
#> 2903          Ifrd2             (C-A)  7  766
#> 2904          Rc3h2             (C-A)  7  766
#> 2905            Ttk             (C-A)  7  766
#> 2906      Nutf2-ps1             (C-A)  7  766
#> 2907         Gm7931             (C-A)  7  766
#> 2908          Clcn7             (C-A)  7  766
#> 2909        Mettl14             (C-A)  7  766
#> 2910      Rps16-ps2             (C-A)  7  766
#> 2911        Gm14494             (C-A)  7  766
#> 2912        Pomgnt2             (C-A)  7  766
#> 2913         Utp14a             (C-A)  7  766
#> 2914     Rpl21-ps15             (C-A)  7  766
#> 2915           Skil             (C-A)  7  766
#> 2916         Exosc5             (C-A)  7  766
#> 2917        mt-Atp6             (C-A)  7  766
#> 2918           Ulk3             (C-A)  7  766
#> 2919           Cinp             (C-A)  7  766
#> 2920          Plcl2             (C-A)  7  766
#> 2921        Tmem115             (C-A)  7  766
#> 2922          Pde12             (C-A)  7  766
#> 2923          Pnpt1             (C-A)  7  766
#> 2924            Apc             (C-A)  7  766
#> 2925          Myo9b             (C-A)  7  766
#> 2926          Apex1             (C-A)  7  766
#> 2927         Fam63a             (C-A)  7  766
#> 2928          Bnip3             (C-A)  7  766
#> 2929         Gm4202             (C-A)  7  766
#> 2930        Gpbp1l1             (C-A)  7  766
#> 2931         Gtf3c4             (C-A)  7  766
#> 2932          Rad52             (C-A)  7  766
#> 2933        Slc16a7             (C-A)  7  766
#> 2934            Lbp             (C-A)  7  766
#> 2935        Gm21769             (C-A)  7  766
#> 2936           Nom1             (C-A)  7  766
#> 2937        Gm13341             (C-A)  7  766
#> 2938          Ddx49             (C-A)  7  766
#> 2939           Rmi1             (C-A)  7  766
#> 2940         Dicer1             (C-A)  7  766
#> 2941          Pcf11             (C-A)  7  766
#> 2942          Snx18             (C-A)  7  766
#> 2943          Aasdh             (C-A)  7  766
#> 2944           Sun1             (C-A)  7  766
#> 2945         Slc2a6             (C-A)  7  766
#> 2946         Hivep2             (C-A)  7  766
#> 2947          Kif2a             (C-A)  7  766
#> 2948         Gm2423             (C-A)  7  766
#> 2949          Socs5             (C-A)  7  766
#> 2950          Ddx56             (C-A)  7  766
#> 2951      Rps24-ps3             (C-A)  7  766
#> 2952         Cep104             (C-A)  7  766
#> 2953          Ccbl2             (C-A)  7  766
#> 2954          Jmjd8             (C-A)  7  766
#> 2955            Pcx             (C-A)  7  766
#> 2956        Tmsb15l             (C-A)  7  766
#> 2957         Ctnna1             (C-A)  7  766
#> 2958            Ank             (C-A)  7  766
#> 2959        Gm21887             (C-A)  7  766
#> 2960      Hmgn2-ps1             (C-A)  7  766
#> 2961           Zhx2             (C-A)  7  766
#> 2962          Rcan1             (C-A)  7  766
#> 2963          Sfxn2             (C-A)  7  766
#> 2964         Gm6863             (C-A)  7  766
#> 2965            Emd             (C-A)  7  766
#> 2966          Taok2             (C-A)  7  766
#> 2967          Meaf6             (C-A)  7  766
#> 2968        Slc27a2             (C-A)  7  766
#> 2969          Stk24             (C-A)  7  766
#> 2970           Ipmk             (C-A)  7  766
#> 2971          Nthl1             (C-A)  7  766
#> 2972         Gm3699             (C-A)  7  766
#> 2973          Tlr12             (C-A)  7  766
#> 2974          Sirt2             (C-A)  7  766
#> 2975          S1pr4             (C-A)  7  766
#> 2976          Shprh             (C-A)  7  766
#> 2977        Foxred2             (C-A)  7  766
#> 2978        Atp6v0b             (C-A)  7  766
#> 2979         Gm7846             (C-A)  7  766
#> 2980  1600012H06Rik             (C-A)  7  766
#> 2981         Pnpla6             (C-A)  7  766
#> 2982       BC068281             (C-A)  7  766
#> 2983          Gnptg             (C-A)  7  766
#> 2984          Prdm2             (C-A)  7  766
#> 2985           Idh1             (C-A)  7  766
#> 2986           Exd2             (C-A)  7  766
#> 2987         Wdr45b             (C-A)  7  766
#> 2988         Zfp599             (C-A)  7  766
#> 2989         Man1a2             (C-A)  7  766
#> 2990          Hmgn5             (C-A)  7  766
#> 2991          Fem1c             (C-A)  7  766
#> 2992           Dexi             (C-A)  7  766
#> 2993          Prpf3             (C-A)  7  766
#> 2994           Pdxk             (C-A)  7  766
#> 2995           Uxs1             (C-A)  7  766
#> 2996          Ttpal             (C-A)  7  766
#> 2997       Gtf2ird1             (C-A)  7  766
#> 2998         Gm3608             (C-A)  7  766
#> 2999         Dusp12             (C-A)  7  766
#> 3000          Scrn3             (C-A)  7  766
#> 3001          Il1r1             (C-A)  7  766
#> 3002         Cfap20             (C-A)  7  766
#> 3003         Iqgap3             (C-A)  7  766
#> 3004        Camsap3             (C-A)  7  766
#> 3005           Eldr             (C-A)  7  766
#> 3006         Lrsam1             (C-A)  7  766
#> 3007            Mcu             (C-A)  7  766
#> 3008           Pan3             (C-A)  7  766
#> 3009         Golgb1             (C-A)  7  766
#> 3010           Map7             (C-A)  7  766
#> 3011         Cdc14a             (C-A)  7  766
#> 3012          Chpf2             (C-A)  7  766
#> 3013         Unc45a             (C-A)  7  766
#> 3014           Cbx2             (C-A)  7  766
#> 3015  B430212C06Rik             (C-A)  7  766
#> 3016         Gm9843             (C-A)  7  766
#> 3017          Maml1             (C-A)  7  766
#> 3018          Trak1             (C-A)  7  766
#> 3019          Usp20             (C-A)  7  766
#> 3020          Ahnak             (C-A)  7  766
#> 3021        Slc43a3             (C-A)  7  766
#> 3022           Lnx2             (C-A)  7  766
#> 3023          Nudt6             (C-A)  7  766
#> 3024         Arfip1             (C-A)  7  766
#> 3025  8430419L09Rik             (C-A)  7  766
#> 3026        Plekha3             (C-A)  7  766
#> 3027          Anks3             (C-A)  7  766
#> 3028        Slc30a5             (C-A)  7  766
#> 3029          Tmppe             (C-A)  7  766
#> 3030          Taf4b             (C-A)  7  766
#> 3031          Zc3h3             (C-A)  7  766
#> 3032          Mctp2             (C-A)  7  766
#> 3033  5031425E22Rik             (C-A)  7  766
#> 3034          Extl3             (C-A)  7  766
#> 3035        Gm14597             (C-A)  7  766
#> 3036           Rpia             (C-A)  7  766
#> 3037        Ralgps2             (C-A)  7  766
#> 3038           Lars             (C-A)  7  766
#> 3039         Plscr1             (C-A)  7  766
#> 3040            Axl             (C-A)  7  766
#> 3041           Ly6d             (C-A)  7  766
#> 3042         Hexim1             (C-A)  7  766
#> 3043         Spire1             (C-A)  7  766
#> 3044         Ogfod3             (C-A)  7  766
#> 3045          Med19             (C-A)  7  766
#> 3046          Fkbp7             (C-A)  7  766
#> 3047          Asah1             (C-A)  7  766
#> 3048        Ccdc88b             (C-A)  7  766
#> 3049         Alg10b             (C-A)  7  766
#> 3050          Atg4c             (C-A)  7  766
#> 3051      Rps10-ps2             (C-A)  7  766
#> 3052         Zfp583             (C-A)  7  766
#> 3053         Crybg3             (C-A)  7  766
#> 3054          Efnb1             (C-A)  7  766
#> 3055            Mt1             (C-A)  7  766
#> 3056       Tor1aip1             (C-A)  7  766
#> 3057         Ptpn23             (C-A)  7  766
#> 3058         Ccdc55             (C-A)  7  766
#> 3059  4933426M11Rik             (C-A)  7  766
#> 3060  D630044L22Rik             (C-A)  7  766
#> 3061          Cxcr4             (C-A)  7  766
#> 3062          Neat1             (C-A)  7  766
#> 3063            Ddt             (C-A)  7  766
#> 3064        Sertad1             (C-A)  7  766
#> 3065          Msrb1             (C-A)  7  766
#> 3066         Gm8524             (C-A)  7  766
#> 3067        Nckap5l             (C-A)  7  766
#> 3068        Rps6kb2             (C-A)  7  766
```
