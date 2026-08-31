#' Plot a venn diagram, UpSet plot, or table of intersections
#'
#' Generates Venn diagram of intersections across a series of sets (e.g., intersections of significant genes across
#' tested contrasts). This Venn diagram is available for up to five sets; Intersection plot is available for any number
#' of sets. Specific sets can be selected for the visualizations and the returned dataset may include all (default) or
#' specified intersections.
#' An S7 generic for `multiOmicDataSet` inputs.
#'
#' @param moo_diff_summary_dat multiOmicDataSet containing per-contrast differential expression results.
#' @param feature_id_colname Feature identifier column in each contrast result. Defaults to the first column.
#' @param signif_colname Significance column in each contrast result. Default: `"adjpval"`.
#' @param signif_threshold Features pass the significance threshold when their value is less than this cutoff
#'   (exclusive). Default: `0.05`.
#' @param change_colname Log fold-change column in each contrast result. Default: `"logFC"`.
#' @param change_threshold Features pass the fold-change threshold when their absolute value is greater than this
#'   cutoff (exclusive). Default: `1.0`.
#' @param select_contrasts Contrast names to include. If empty, all available contrasts are used.
#' @inheritParams plot_venn_sets
#'
#' @export
plot_venn_diagram <- S7::new_generic(
  "plot_venn_diagram",
  "moo_diff_summary_dat",
  function(
    moo_diff_summary_dat,
    feature_id_colname = NULL,
    signif_colname = "adjpval",
    signif_threshold = 0.05,
    change_colname = "logFC",
    change_threshold = 1.0,
    select_contrasts = c(),
    plot_type = "Venn diagram",
    intersection_ids = c(),
    venn_force_unique = TRUE,
    venn_numbers_format = "raw",
    venn_significant_digits = 2,
    venn_fill_colors = c(
      "darkgoldenrod2",
      "darkolivegreen2",
      "mediumpurple3",
      "darkorange2",
      "lightgreen"
    ),
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
  ) {
    return(S7::S7_dispatch())
  }
)

#' @rdname plot_venn_diagram
S7::method(plot_venn_diagram, multiOmicDataSet) <- function(
  moo_diff_summary_dat,
  feature_id_colname = NULL,
  signif_colname = "adjpval",
  signif_threshold = 0.05,
  change_colname = "logFC",
  change_threshold = 1.0,
  select_contrasts = c(),
  plot_type = "Venn diagram",
  intersection_ids = c(),
  venn_force_unique = TRUE,
  venn_numbers_format = "raw",
  venn_significant_digits = 2,
  venn_fill_colors = c(
    "darkgoldenrod2",
    "darkolivegreen2",
    "mediumpurple3",
    "darkorange2",
    "lightgreen"
  ),
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
) {
  diff_results <- moo_diff_summary_dat@analyses$diff
  if (!is.list(diff_results) || length(diff_results) == 0) {
    stop(
      "No per-contrast differential expression results found in moo@analyses$diff"
    )
  }
  if (length(select_contrasts) > 0) {
    missing_contrasts <- setdiff(select_contrasts, names(diff_results))
    if (length(missing_contrasts) > 0) {
      stop(glue::glue(
        "Selected contrasts not found: {paste(missing_contrasts, collapse = ', ')}"
      ))
    }
    diff_results <- diff_results[select_contrasts]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- names(diff_results[[1]])[1]
  }

  venn_sets <- lapply(diff_results, function(diff_result) {
    required_columns <- c(feature_id_colname, signif_colname, change_colname)
    missing_columns <- setdiff(required_columns, names(diff_result))
    if (length(missing_columns) > 0) {
      stop(glue::glue(
        "Required DEG columns not found: {paste(missing_columns, collapse = ', ')}"
      ))
    }
    is_significant <- !is.na(diff_result[[signif_colname]]) &
      diff_result[[signif_colname]] < signif_threshold
    has_change <- !is.na(diff_result[[change_colname]]) &
      abs(diff_result[[change_colname]]) > change_threshold
    return(unique(diff_result[[feature_id_colname]][
      is_significant & has_change
    ]))
  })

  return(
    plot_venn_sets(
      venn_sets = venn_sets,
      plot_type = plot_type,
      intersection_ids = intersection_ids,
      venn_force_unique = venn_force_unique,
      venn_numbers_format = venn_numbers_format,
      venn_significant_digits = venn_significant_digits,
      venn_fill_colors = venn_fill_colors,
      venn_fill_transparency = venn_fill_transparency,
      venn_border_colors = venn_border_colors,
      venn_font_size_for_category_names = venn_font_size_for_category_names,
      venn_category_names_distance = venn_category_names_distance,
      venn_category_names_position = venn_category_names_position,
      venn_font_size_for_counts = venn_font_size_for_counts,
      venn_outer_margin = venn_outer_margin,
      intersections_order = intersections_order,
      display_empty_intersections = display_empty_intersections,
      intersection_bar_color = intersection_bar_color,
      intersection_point_size = intersection_point_size,
      intersection_line_width = intersection_line_width,
      table_font_size = table_font_size,
      table_content = table_content,
      graphics_device = graphics_device,
      dpi = dpi,
      image_width = image_width,
      image_height = image_height,
      plot_filename = plot_filename,
      print_plots = print_plots,
      save_plots = save_plots,
      plots_subdir = plots_subdir
    )
  )
}

#' Build Venn or UpSet intersections from feature sets
#'
#' @param venn_sets Named list of feature identifier vectors, one per contrast.
#' @param plot_type Type of plot to generate: "Venn diagram" or "Intersection plot". Default: "Venn diagram"
#' @param intersection_ids A vector of intersection IDs to select for the plot. If empty, all intersections are used.
#' @param venn_force_unique If TRUE, forces unique elements in the Venn diagram. Default: TRUE
#' @param venn_numbers_format Format for the numbers in the Venn diagram. Options: "raw", "percent", "raw-percent",
#'   "percent-raw". Default: "raw"
#' @param venn_significant_digits Number of significant digits for the Venn diagram numbers. Default: 2
#' @param venn_fill_colors A vector of colors to fill the Venn diagram categories. Default: c("darkgoldenrod2",
#'   "darkolivegreen2", "mediumpurple3", "darkorange2", "lightgreen")
#' @param venn_fill_transparency Transparency level for the Venn diagram fill colors. Default: 0.2
#' @param venn_border_colors Colors for the borders of the Venn diagram categories. Default: "fill colors" (uses the
#'   same colors as `venn_fill_colors`)
#' @param venn_font_size_for_category_names Font size for the category names in the Venn diagram. Default: 3
#' @param venn_category_names_distance Distance of the category names from the Venn diagram circles. Default: c()
#' @param venn_category_names_position Position of the category names in the Venn diagram. Default: c()
#' @param venn_font_size_for_counts Font size for the counts in the Venn diagram. Default: 6
#' @param venn_outer_margin Outer margin for the Venn diagram. Default: 0
#' @param intersections_order Order of the intersections in the plot. Default: "by size"
#' @param display_empty_intersections If TRUE, displays empty intersections in the plot. Default: FALSE
#' @param intersection_bar_color Color for the intersection bars in the plot. Default: "lightgray"
#' @param intersection_point_size Size of the points in the intersection plot. Default: 2
#' @param intersection_line_width Width of the lines in the intersection plot. Default: 0.5
#' @param table_font_size Font size for the table in the plot. Default: 3
#' @param table_content Content of the table in the plot. Default: NULL
#'
#' @keywords internal
#'
plot_venn_sets <- function(
  venn_sets,
  plot_type = "Venn diagram",
  intersection_ids = c(),
  venn_force_unique = TRUE,
  venn_numbers_format = "raw",
  venn_significant_digits = 2,
  venn_fill_colors = c(
    "darkgoldenrod2",
    "darkolivegreen2",
    "mediumpurple3",
    "darkorange2",
    "lightgreen"
  ),
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
) {
  Freq <- Gene <- Id <- Size <- Var1 <- NULL
  abort_packages_not_installed(c(
    "VennDiagram",
    "gridExtra",
    "patchwork",
    "UpSetR"
  ))

  if (
    !is.list(venn_sets) || length(venn_sets) == 0 || is.null(names(venn_sets))
  ) {
    stop("venn_sets must be a non-empty named list")
  }
  vlist <- venn_sets
  num_categories <- length(vlist)
  if (num_categories == 0) {
    stop("Zero categories found")
  }

  # generate upset object

  # modify UpSetR function (keep gene names as rownames of intersection matrix)
  fromList <- function(input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
      return(as.vector(match(elements, x)))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = FALSE))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
  }

  if (num_categories > 1) {
    sets <- fromList(vlist)

    Intersection <- sets

    # generate intersection frequency table and gene list (all intersections for the output dataset/table not the plot)
    intersection_matrix <- Intersection
    Intersection <- sapply(colnames(intersection_matrix), function(x) {
      return(ifelse(intersection_matrix[, x] == 1, x, "{}"))
    })
    rownames(Intersection) <- rownames(sets)
    Intersection <- apply(Intersection, 1, function(x) {
      return(sprintf("(%s)", paste(x, collapse = " ")))
    })
    tab <- table(Intersection)
    tab <- tab[order(tab)]
    nn <- stringr::str_count(names(tab), pattern = "\\{\\}")
    tab <- tab[order(nn, decreasing = FALSE)]
    names(tab) <- gsub("\\{\\} | \\{\\}|\\{\\} |\\{\\} \\{\\}", "", names(tab))
    names(tab) <- sub("\\( ", "(", names(tab))
    names(tab) <- gsub(" ", " \u2229 ", names(tab))
    tab <- tab[names(tab) != "()"] |>
      data.frame() |>
      dplyr::rename("Intersection" = Var1, "Size" = Freq) |>
      tibble::rownames_to_column("Id") |>
      dplyr::mutate(Id = as.numeric(Id)) |>
      dplyr::select(Intersection, Id, Size)
    Intersection <- gsub(
      "\\{\\} | \\{\\}|\\{\\} |\\{\\} \\{\\}",
      "",
      Intersection
    )
    Intersection <- sub("\\( ", "(", Intersection)
    Intersection <- gsub(" ", " \u2229 ", Intersection)
    Intersection <- data.frame(Intersection) |>
      tibble::rownames_to_column("Gene") |>
      dplyr::inner_join(tab, by = c(Intersection = "Intersection")) |>
      dplyr::select(Gene, Intersection, Id, Size) |>
      dplyr::arrange(Id)
  } else if (num_categories == 1) {
    Intersection <- data.frame(
      Gene = vlist[[1]],
      Intersection = sprintf("(%s)", names(vlist)),
      Id = 1,
      Size = length(vlist[[1]])
    )
    tab <- table(Intersection$Intersection)
    tab <- data.frame(Id = 1, tab) |>
      dplyr::rename(Intersection = Var1, Size = Freq) |>
      dplyr::select(Intersection, Id, Size)
  }

  # returned intersections

  if (length(intersection_ids) > 0) {
    intersection_ids <- sort(as.numeric(intersection_ids))
    tabsel <- tab[tab$Id %in% intersection_ids, ]
    Intersectionsel <- Intersection[Intersection$Id %in% intersection_ids, ]
  } else {
    tabsel <- tab
    Intersectionsel <- Intersection
  }
  tab$"Return" <- ifelse(
    tab$Intersection %in% tabsel$Intersection,
    "Yes",
    "\u2014"
  )

  if (intersections_order == "freq") {
    tab <- tab |> dplyr::arrange(-Size)
    tabsel <- tabsel |> dplyr::arrange(-Size)
  }

  message(glue::glue("All intersections: {paste(tab, collapse=',')}"))
  message(glue::glue("\nIntersections returned: {paste(tabsel, collapse=',')}"))

  ### PH: End Create venn Table from DEG table

  ### PH: START This section sets immage output parameters and includes error check

  # Logic Error Check ====
  if (num_categories == 1) {
    if (plot_type == "Intersection plot") {
      plot_type <- "Venn diagram"
      cat(
        "\nIntersection plot not available for a single contrast, the Venn diagram generated instead"
      )
    }
  } else if (num_categories > 5) {
    plot_type <- "Intersection plot"
    cat(
      "\nVenn diagram available for up to 5 contrasts, the Intersection plot generated instead"
    )
  }

  ### PH: End This section sets immage output parameters and includes error check

  # DO PLOT ====

  ### PH: START Create Intersection plot

  # Intersection Plot

  if (plot_type == "Intersection plot") {
    # do plot
    empty <- display_empty_intersections
    if (empty) {
      keepEmpty <- "on"
    } else {
      keepEmpty <- NULL
    }

    barcol <- intersection_bar_color

    pSet <- UpSetR::upset(
      sets,
      nsets = num_categories,
      sets = names(vlist),
      order.by = intersections_order,
      nintersects = NA,
      text.scale = 2,
      empty.intersections = keepEmpty,
      matrix.color = barcol,
      main.bar.color = barcol,
      sets.bar.color = barcol,
      point.size = intersection_point_size,
      line.size = intersection_line_width
    )

    output_plot <- pSet
    ### PH: End Create Intersection plot

    ### PH: START Create Venn Diagram
  } else if (plot_type == "Venn diagram") {
    # Venn diagram

    ## If venn fill color param empty upon template upgrade,
    ## then fill it with the default colors.
    if (length(venn_fill_colors) == 0) {
      venn_fill_colors <- c(
        "darkgoldenrod2",
        "darkolivegreen2",
        "mediumpurple3",
        "darkorange2",
        "lightgreen"
      )
    }

    color_border <- venn_border_colors
    if (color_border != "black") {
      color_border <- venn_fill_colors[1:num_categories]
    }

    print_mode <- venn_numbers_format
    if (print_mode == "raw-percent") {
      print_mode <- c("raw", "percent")
    } else if (print_mode == "percent-raw") {
      print_mode <- c("percent", "raw")
    }

    distance <- venn_category_names_distance
    position <- venn_category_names_position

    if (is.null(distance) && is.null(position)) {
      vobj <- VennDiagram::venn.diagram(
        vlist,
        file = NULL,
        force_unique = venn_force_unique,
        print.mode = print_mode,
        sigdigs = venn_significant_digits,
        margin = venn_outer_margin,
        main = "",
        cat.cex = venn_font_size_for_category_names,
        cex = venn_font_size_for_counts,
        main.cex = 3,
        fill = venn_fill_colors[1:num_categories],
        alpha = venn_fill_transparency,
        col = color_border
      )
    } else if (!is.null(distance) && is.null(position)) {
      distance <- as.numeric(distance)

      vobj <- VennDiagram::venn.diagram(
        vlist,
        file = NULL,
        force_unique = venn_force_unique,
        print.mode = print_mode,
        sigdigs = venn_significant_digits,
        margin = venn_outer_margin,
        main = "",
        cat.cex = venn_font_size_for_category_names,
        cex = venn_font_size_for_counts,
        main.cex = 3,
        fill = venn_fill_colors[1:num_categories],
        alpha = venn_fill_transparency,
        col = color_border,
        cat.dist = distance
      )
    } else if (is.null(distance) && !is.null(position)) {
      position <- as.numeric(position)

      vobj <- VennDiagram::venn.diagram(
        vlist,
        file = NULL,
        force_unique = venn_force_unique,
        print.mode = print_mode,
        sigdigs = venn_significant_digits,
        margin = venn_outer_margin,
        main = "",
        cat.cex = venn_font_size_for_category_names,
        cex = venn_font_size_for_counts,
        main.cex = 3,
        fill = venn_fill_colors[1:num_categories],
        alpha = venn_fill_transparency,
        col = color_border,
        cat.pos = position
      )
    } else {
      distance <- as.numeric(distance)
      position <- as.numeric(position)

      vobj <- VennDiagram::venn.diagram(
        vlist,
        file = NULL,
        force_unique = venn_force_unique,
        print.mode = print_mode,
        sigdigs = venn_significant_digits,
        margin = venn_outer_margin,
        main = "",
        cat.cex = venn_font_size_for_category_names,
        cex = venn_font_size_for_counts,
        main.cex = 3,
        fill = venn_fill_colors[1:num_categories],
        alpha = venn_fill_transparency,
        col = color_border,
        cat.dist = distance,
        cat.pos = position
      )
    }

    pVenn <- patchwork::wrap_elements(grid::gTree(children = vobj))
    output_plot <- pVenn

    ### PH: END Create Venn Diagram

    ### PH: START Create figure from Intersect table
    ## this might be a NIDAP specific function. as long as we have access to the table in the MOobject
  } else {
    font_size_table <- table_font_size
    table_content <- table_content
    if (table_content == "all intersections") {
      pTab <- patchwork::wrap_elements(gridExtra::tableGrob(
        tab,
        rows = NULL,
        theme = gridExtra::ttheme_default(
          core = list(fg_params = list(cex = font_size_table)),
          colhead = list(fg_params = list(cex = font_size_table)),
          rowhead = list(fg_params = list(cex = font_size_table))
        )
      ))
    } else {
      pTab <- patchwork::wrap_elements(gridExtra::tableGrob(
        tabsel,
        rows = NULL,
        theme = gridExtra::ttheme_default(
          core = list(fg_params = list(cex = font_size_table)),
          colhead = list(fg_params = list(cex = font_size_table)),
          rowhead = list(fg_params = list(cex = font_size_table))
        )
      ))
    }
    output_plot <- pTab
  }

  print_or_save_plot(
    output_plot,
    filename = file.path(plots_subdir, plot_filename),
    device = graphics_device,
    dpi = dpi,
    width = image_width,
    height = image_height,
    units = "px",
    print_plots = print_plots,
    save_plots = save_plots
  )
  ### PH: END Create figure from Intersect table

  # SAVE DATASET ====
  return(Intersectionsel)
}
