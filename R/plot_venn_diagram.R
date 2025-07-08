# Venn Diagram [CCBR] (1f11bf8e-ddf5-43ba-aca0-10ba466df2b4): v569
#' generates Venn diagram of intersections across a series of sets (e.g., intersections of significant genes across tested contrasts). This Venn diagram is available for up to five sets; Intersection plot is available for any number of sets. Specific sets can be selected for the visualizations and the returned dataset may include all (default) or specified intersections.
plot_venn_diagram <- function(VolcanoSummary) {

  ## --------- ##
  ## Libraries ##
  ## --------- ##

  library(VennDiagram)
  library(gridExtra)
  library(patchwork)
  library(UpSetR)
  library(dplyr)
  library(tibble)

  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  ### PH:
  # Input - DEG table from Volcano Summary, I think we need to make this function more generic The input should be the Limma DEG table and maybe be used with the DEG Gene List Template
  # Output - Venn Diagram Figure + Venn table
  # Purpose - compare DEGS from differetn Comparisons


  # Basic parameters
  input_dataset <- VolcanoSummary
  elements_column <- 'Gene'
  categories_column <- 'Contrast'

  # Advanced parameters
  selected_categories = c()
  select_plot_type <- 'Venn diagram'
  intersection_ids = c()

  # Venn Diagram parameters
  venn_force_unique  <- TRUE
  venn_numbers_format <- 'raw'
  venn_significant_digits <- 2
  venn_fill_colors = c("darkgoldenrod2","darkolivegreen2","mediumpurple3","darkorange2","lightgreen")
  venn_fill_transparency <- 0.2
  venn_border_colors <- 'fill colors'
  venn_font_size_for_category_names <- 3
  venn_category_names_distance <- c()
  venn_category_names_position <- c()
  venn_font_size_for_counts <- 6
  venn_outer_margin <- 0

  # Intersection Plot parameters
  intersections_order <- 'degree'
  display_empty_intersections <- FALSE
  intersection_bar_color <- "steelblue4"
  intersection_point_size <- 2.2
  intersection_line_width <-  0.7

  # Table parameters
  table_font_size <- 0.7
  table_content <- "all intersections"

  # Image parameters
  image_output_format <- 'png'
  image_resolution <- 300
  image_width <- 4000
  image_height <- 3000

  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##

  ## --------- ##
  ## Functions ##
  ## --------- ##




  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##


  ### PH: Create venn Table from DEG table

  # SET INPUT ====

  # select required columns
  set_elements <- input_dataset[, elements_column]
  set_names <-    input_dataset[, categories_column]

  # prepare format - R list
  vlist = split(set_elements, set_names)
  if(!is.null(selected_categories)){
    vlist = vlist[selected_categories]
  }
  num_categories = length(vlist)

  # generate upset object

  # modify UpSetR function (keep gene names as rownames of intersection matrix)
  fromList <- function (input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
  }


  if(num_categories > 1) {

    sets = fromList(vlist)

    if(!is.null(selected_categories)){
      Intersection = sets[,match(selected_categories, colnames(sets))]
    } else {
      Intersection = sets
    }

    # generate intersection frequency table and gene list (all intersections for the output dataset/table not the plot)
    Intersection = sapply(colnames(Intersection), function(x){ifelse(Intersection[,x]==1, x, "{}")})
    rownames(Intersection) = rownames(sets)
    Intersection = apply(Intersection, 1, function(x) sprintf("(%s)", paste(x, collapse=' ')))
    tab = table(Intersection)
    tab = tab[order(tab)]
    nn = stringr::str_count(names(tab), pattern = "\\{\\}")
    tab = tab[order(nn, decreasing=FALSE)]
    names(tab) = gsub("\\{\\} | \\{\\}|\\{\\} |\\{\\} \\{\\}","", names(tab))
    names(tab) = sub("\\( ","(", names(tab))
    names(tab) = gsub(" "," ∩ ", names(tab))
    tab = tab[names(tab) != "()"] %>% data.frame() %>% dplyr::rename("Intersection"=Var1, "Size"=Freq) %>% tibble::rownames_to_column('Id') %>% dplyr::mutate(Id=as.numeric(Id))  %>% dplyr::select(Intersection, Id, Size)
    Intersection = gsub("\\{\\} | \\{\\}|\\{\\} |\\{\\} \\{\\}","", Intersection)
    Intersection = sub("\\( ","(", Intersection)
    Intersection = gsub(" "," ∩ ", Intersection)
    Intersection = data.frame(Intersection) %>% tibble::rownames_to_column("Gene") %>% dplyr::inner_join(tab, by=c(Intersection="Intersection")) %>% dplyr::select(Gene, Intersection, Id, Size) %>% dplyr::arrange(Id)

  } else if (num_categories == 1) {
    Intersection = data.frame(Gene=vlist[[1]], Intersection = sprintf("(%s)", names(vlist)), Id = 1, Size = length(vlist[[1]]))
    tab = table(Intersection$Intersection)
    tab = data.frame(Id=1, tab) %>% dplyr::rename(Intersection=Var1, Size=Freq) %>% dplyr::select(Intersection, Id, Size)
  }

  # returned intersections

  if (!is.null(intersection_ids) ) {
    intersection_ids = sort(as.numeric(intersection_ids))
    tabsel = tab[tab$Id %in% intersection_ids,]
    Intersectionsel = Intersection[Intersection$Id %in% intersection_ids,]
  } else {
    tabsel = tab
    Intersectionsel = Intersection
  }
  tab$"Return" = ifelse(tab$Intersection %in% tabsel$Intersection, "Yes", "—")

  if(intersections_order == 'freq'){
    tab = tab %>% dplyr::arrange(-Size)
    tabsel = tabsel %>% dplyr::arrange(-Size)
  }


  # screen log
  cat('All intersections\n')
  print(tab)
  cat('\nIntersections returned\n')
  print(tabsel)

  ### PH: End Create venn Table from DEG table


  ### PH: START This section sets immage output parmeters and includes error check

  # SET IMAGE ====

  if(image_output_format == 'png') {
    png(filename=graphicsFile, width=image_width, height=image_height, units="px", pointsize=4, bg="white", res=image_resolution, type="cairo")
  } else {
    svglite::svglite(file=graphicsFile, width=round(image_width/image_resolution,digits=2), height=round(image_height/image_resolution,digits=2), pointsize=1, bg="white")
  }

  # Logic Error Check ====
  if(num_categories == 1) {
    if(select_plot_type == "Intersection plot") {
      select_plot_type = 'Venn diagram'
      cat("\nIntersection plot not available for a single contrast, the Venn diagram genereated instead")
    }
  } else if(num_categories > 5) {
    select_plot_type = 'Intersection plot'
    cat("\nVenn diagram available for up to 5 contrasts, the Intersection plot genereated instead")
  }

  ### PH: End This section sets immage output parmeters and includes error check

  # DO PLOT ====

  ### PH: START Create Intersection plot

  # Intersection Plot

  if(select_plot_type == 'Intersection plot') {

    # do plot
    empty = display_empty_intersections
    if(empty) {keepEmpty='on'} else {keepEmpty=NULL}

    barcol = intersection_bar_color

    pSet = upset(sets,
                 nsets = num_categories,
                 sets = selected_categories,
                 order.by = intersections_order,
                 nintersects = NA,
                 text.scale = 2,
                 empty.intersections = keepEmpty,
                 matrix.color = barcol, main.bar.color = barcol, sets.bar.color = barcol,
                 point.size =intersection_point_size, line.size = intersection_line_width)

    print(pSet)
    ### PH: End Create Intersection plot

    ### PH: START Create Venn Diagram

  } else if (select_plot_type == 'Venn diagram') {
    # Venn diagram

    ## If venn fill color param empty upon template upgrade,
    ## then fill it with the default colors.
    if (length(venn_fill_colors) == 0) {
      venn_fill_colors <- c("darkgoldenrod2","darkolivegreen2","mediumpurple3","darkorange2","lightgreen")
    }

    color_border = venn_border_colors
    if(color_border != 'black') { color_border = venn_fill_colors[1:num_categories] }

    print_mode = venn_numbers_format
    if(print_mode == 'raw-percent') {
      print_mode = c('raw','percent')
    } else if(print_mode == 'percent-raw'){
      print_mode = c('percent','raw')
    }

    distance = venn_category_names_distance
    position = venn_category_names_position

    if( is.null(distance) & is.null(position) ) {

      vobj = venn.diagram( vlist, file=NULL, force_unique = venn_force_unique, print.mode = print_mode, sigdigs = venn_significant_digits, margin=venn_outer_margin,  main = '', cat.cex=venn_font_size_for_category_names, cex=venn_font_size_for_counts,  main.cex=3, fill=venn_fill_colors[1:num_categories], alpha=venn_fill_transparency, col=color_border )

    } else if ( !is.null(distance) & is.null(position) ) {

      distance = as.numeric(distance)

      vobj = venn.diagram( vlist, file=NULL, force_unique = venn_force_unique, print.mode = print_mode, sigdigs = venn_significant_digits, margin=venn_outer_margin,  main = '', cat.cex=venn_font_size_for_category_names, cex=venn_font_size_for_counts,  main.cex=3, fill=venn_fill_colors[1:num_categories], alpha=venn_fill_transparency, col=color_border, cat.dist = distance)

    } else if ( is.null(distance) & !is.null(position) ) {

      position = as.numeric(position)

      vobj = venn.diagram( vlist, file=NULL, force_unique = venn_force_unique, print.mode = print_mode, sigdigs = venn_significant_digits, margin=venn_outer_margin,  main = '', cat.cex=venn_font_size_for_category_names, cex=venn_font_size_for_counts,  main.cex=3, fill=venn_fill_colors[1:num_categories], alpha=venn_fill_transparency, col=color_border, cat.pos = position)

    } else {

      distance = as.numeric(distance)
      position = as.numeric(position)

      vobj = venn.diagram( vlist, file=NULL, force_unique = venn_force_unique, print.mode = print_mode, sigdigs = venn_significant_digits, margin=venn_outer_margin,  main = '', cat.cex=venn_font_size_for_category_names, cex=venn_font_size_for_counts,  main.cex=3, fill=venn_fill_colors[1:num_categories], alpha=venn_fill_transparency, col=color_border, cat.dist = distance, cat.pos = position)

    }

    pVenn = wrap_elements( gTree(children=vobj) )
    print(pVenn)

    ### PH: END Create Venn Diagram

    ### PH: START Create figure from Intersect table
    ## this might be a NIDAP specific function. as long as we have access to the table in the MOobject
  } else {

    font_size_table = table_font_size
    table_content = table_content
    if(table_content == "all intersections"){
      pTab = wrap_elements( tableGrob(tab, rows=NULL, theme = ttheme_default(core=list(fg_params=list(cex=font_size_table)), colhead = list(fg_params=list(cex = font_size_table)), rowhead=list(fg_params=list(cex= font_size_table)))) )
    } else {
      pTab = wrap_elements( tableGrob(tabsel, rows=NULL, theme = ttheme_default(core=list(fg_params=list(cex=font_size_table)), colhead = list(fg_params=list(cex = font_size_table)), rowhead=list(fg_params=list(cex= font_size_table)))) )
    }
    print(pTab)
  }

  ### PH: END Create figure from Intersect table

  # SAVE DATASET ====
  return(Intersectionsel)

}

## ---------------------------- ##
## Global Imports and Functions ##
## ---------------------------- ##

## Functions defined here will be available to call in the code for any table.

## --------------- ##
## End of Template ##
## --------------- ##
