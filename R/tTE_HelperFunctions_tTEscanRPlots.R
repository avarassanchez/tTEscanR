SavePlot <- function(plot, save_format, out_name, out_directory){

  ###
  # CALL: Each plot generating function
  # DESCRIPTION: This function is called every time the user specifies that the plot needs to be exported into a .png or .pdf format.
  ###

  if (save_format == "png" || save_format == "pdf"){
    output_file <- GetOutputName(action = "plot", out_name = out_name, out_directory = out_directory, save_format = save_format)
    ggplot2::ggsave(output_file, plot = plot, width = 8, height = 6)
  } else {
    message("The input `save_format` is not recognized.\n", "Supported formats: png or pdf.\n", "The plot was generated but not stored.")
  }
}

GetAnnotData <- function(permut_data, sig_data){

  ###
  # CALL: tTE_PermutationPlot()
  # DESCRIPTION: This function merges the permutation and the significance table to generate the permutation plot. Therefore, the naming of the columns is crucial.
  ###

  # Check that the expected columns are present in the data
  CheckValueInData(param = "permut_data", observed = colnames(permut_data), expected = c("codon", "freq", "group", "aa"))
  CheckValueInData(param = "sig_data", observed = colnames(sig_data), expected = c("codon", "group", "sig_adj", "p_val_adj"))

  # Generate the annot_data
  annot_data <- permut_data %>% dplyr::left_join(sig_data, by = c("group", "codon"))
  annot_data <- annot_data %>% dplyr::mutate(label = ifelse(.data$sig_adj, paste("p =", signif(.data$p_val_adj, 2), "*"), paste("p =", signif(.data$p_val_adj, 2))))

  # Generate the anot_data_labels
  annot_data_labels <- annot_data %>%
    dplyr::group_by(.data$codon) %>%
    dplyr::summarise(combined_label = paste(.data$label, collapse = " | "), .groups = "drop")
  annot_data_labels <- annot_data_labels %>% dplyr::distinct(.data$codon, .data$combined_label)

  return(list(annot_data, annot_data_labels))
}

CheckValueInData <- function(param, observed, expected){

  ###
  # CALL: Multiple
  # DESCRIPTION: This function takes the input data and variables to check that all of them are actually present in the data.
  # It does not have a return value except if it reports and error.
  ###

  present <- all(observed %in% expected)
  if (isFALSE(present)) stop(paste("Please, specify a proper", param, "parameter.\n"), paste("Supported formats:", paste(expected, collapse = ", ")))
}

CheckDataInLongFormat <- function(data) {

  ###
  # CALL: Multiple
  # DESCRIPTION: This function checks if the input data is in long format. If the data that is not in long format will not be computed.
  ###

  num_cols <- sapply(data, is.numeric) # One numeric column
  cat_cols <- sapply(data, is.character) | sapply(data, is.factor) # One categorical column

  is_long <- any(cat_cols) & any(num_cols) # Will return TRUE if the data is in long format
  if(isFALSE(is_long)) stop("The input data needs to be in long format, containing at least one numeric and one categorical colummn.")
}

GetOutputName <- function(action, out_name, out_directory, save_format){

  ###
  # CALL: SavePlot().
  # DESCRIPTION: This function is used to save as .pdf and .png plots and matrices as .rds files.
  # Requires an out_name and out_directory but if there are not given default values are used.
  ###

  # Use user-defined values.
  if (!is.null(out_name) && !is.null(out_directory)) output_file <- paste(out_directory, out_name, sep = "/")

  if (is.null(out_directory)){ # Use current working directory
    message("- Parameter `out_directory` has not been specified.\n", "- The plot will be stored in the current directory.")
    out_directory <- getwd()
  }

  if (is.null(out_name)){ # Use default output name
    message("- Parameter `out_name` has not been specified.\n", "- A standard name will be given.")
    if (action == "plot") out_name <- "distribution_plot" else out_name <- "tRNA_expression_matrix"
  }

  output_file <- paste(getwd(), out_name, sep = "/") # Combine output directory and output name

  # Check if the output name contains already the format extension
  if (!(save_format %in% output_file)) output_file <- paste(output_file, save_format, sep = ".")
  message(paste("- The generated", action, "will be save in:", output_file))
  return (output_file) # Returns the output name.
}

GenerateDistPlot <- function(level, target, data, x_axis_col, y_axis_col, color_palette, add_titles, show_legend, ncols, bar_position, facet_col){

  if (level == "jitter"){
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], color = .data[[target]]))
  } else if (level == "dot"){
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], group = .data[[target]], color = .data[[target]]))
  } else {
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], fill = .data[[target]]))
  }

  if(!is.null(facet_col)) plot <- plot + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), ncol = ncols)

  if (level == "jitter") plot <- plot + ggplot2::geom_jitter(size = 0.5) + ggplot2::theme_bw()
  if (level == "barplot") plot <- plot + ggplot2::geom_bar(stat = "identity", position = bar_position) +  ggplot2::theme_bw() # ENABLE TO CUSTOMIZE stat
  if (level == "boxplot") plot <- plot + ggplot2::geom_boxplot() + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (level == "dot") plot <- plot + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text = ggplot2::element_blank())

  # Customize the colors
  if (is.null(color_palette) && length(unique(data[[target]])) < 36) color_palette <- gradual_groups_35
  if (!(is.null(color_palette))){
    if (level == "jitter" || (level == "dot")) plot <- plot + ggplot2::scale_color_manual(values = color_palette) else  plot <- plot + ggplot2::scale_fill_manual(values = color_palette)
  }

  if(isTRUE(add_titles)) plot <- plot + ggplot2::labs(x = x_axis_col, y = paste(x_axis_col, "usage counts")) # + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 15))

  plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90), legend.position = show_legend)
  return(plot)
}

DrawBarCountsPlot <- function(data, var_numerical, var_categorical, var_color, color_palette, show_legend) {

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = stats::reorder(.data[[var_categorical]], .data[[var_numerical]]), y = .data[[var_numerical]], fill = .data[[var_color]])) +
    ggplot2::geom_col(color = "black")

  if (var_categorical != var_color) plot <- plot + ggplot2::geom_text(ggplot2::aes(label = .data[[var_color]]), hjust = -0.3, size = 3)
  if (!(is.null(color_palette)) && (length(color_palette) == length(unique(data[[var_color]])))) plot <- plot + ggplot2::scale_fill_manual(values = color_palette)
  if (length(unique(data[[var_color]])) < 36) plot <- plot + ggplot2::scale_fill_manual(values = gradual_groups_35)

  plot <- plot + ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::theme(legend.position = show_legend) + ggplot2::labs(x = "", y = "")

  return(plot)
}

DrawDonutPlot <- function(data, var_numerical, var_categorical, color_palette, show_legend){

  # Compute percentages
  total <- sum(data[[var_numerical]], na.rm = TRUE)
  if (total == 0) stop("Error: Sum of var_numerical is zero, cannot compute fractions.")
  fraction <- data[[var_numerical]] / total

  # Compute the cumulative percentages and the position of the label
  ymax <- cumsum(fraction)
  ymin <- c(0, utils::head(ymax, -1))
  labelPosition <- (ymax + ymin) / 2

  if (show_legend == "none") label <- paste(data[[var_categorical]], ":", round(data[[var_numerical]], 2)) else label <- paste0(round(data[[var_numerical]], 2))

  # Adjusting the sizes automatically based on the number of categories
  n_cat <- length(unique(data[[var_categorical]]))
  circle_min <- 2 + log10(n_cat) # Radius scales with log of categories
  circle_max <- circle_min + 1
  label_size <- max(1.5, 6 - 0.1 * n_cat) # Label size decreases as n_cat grows
  label_position <- 4 + 0.2 * n_cat # Label position shifts outward with n_cat

  plot <- ggplot2::ggplot(data, ggplot2::aes(ymax = ymax, ymin = ymin, xmax = circle_max, xmin = circle_min, fill = .data[[var_categorical]])) +
    ggplot2::geom_rect() +
    ggplot2::geom_label(x = label_position, ggplot2::aes(y = labelPosition, label = label), size = label_size, nudge_x = 2, nudge_y = 0) +
    ggplot2::coord_polar(theta = "y") + ggplot2::xlim(c(-1, label_position)) + ggplot2::theme_void()

  if (!(is.null(color_palette)) && (length(color_palette) == length(unique(data[[var_categorical]])))) plot <- plot + ggplot2::scale_fill_manual(values = color_palette)
  if (length(unique(data[[var_categorical]])) < 36) plot <- plot + ggplot2::scale_fill_manual(values = gradual_groups_35)

  plot <- plot + ggplot2::theme(legend.position = show_legend)
  return(plot)
}

DrawRadarPlot <- function(data, show_legend, var_color, var_categorical, var_numerical){

  var_color_sym <- rlang::sym(var_color)

  radar_data <- data %>%
    dplyr::select(
      .data[[var_color]],
      .data[[var_categorical]],
      .data[[var_numerical]]
    ) %>%
    # Reshape the data to a wide format
    tidyr::pivot_wider(
      names_from = .data[[var_categorical]],
      values_from = .data[[var_numerical]],
      values_fill = 0
    ) %>%
    # Add a normalization step
    # dplyr::mutate(total_codons = rowSums(dplyr::select(., -!!var_color_sym))) %>%
    # dplyr::mutate(dplyr::across(-c(total_codons, (.data[[var_color]])), ~ .x / total_codons)) %>%
    # dplyr::select(-total_codons)
    dplyr::mutate(total_codons = rowSums(dplyr::select(., -!!var_color_sym))) %>%
    dplyr::mutate(dplyr::across(-c("total_codons", var_color), ~ .x / .data$total_codons)) %>%
    dplyr::select(-"total_codons")


  plot <- ggradar::ggradar(radar_data,
                           background.circle.colour = "white",
                           grid.min = 0, grid.mid = 0.5, grid.max = 1,
                           gridline.min.linetype = 1, gridline.mid.linetype = 1, gridline.max.linetype = 1,
                           group.line.width = 0.8, group.point.size = 2,
                           legend.position = show_legend)
  return(plot)
}

GenerateProportionPlot <- function(level, data, var_numerical, var_categorical, var_color, color_palette, show_legend){

  if (level == "donut") plot <- DrawDonutPlot(data = data, var_numerical = var_numerical, var_categorical = var_categorical, color_palette = color_palette, show_legend = show_legend)
  if (level == "radar") plot <- DrawRadarPlot(data = data, show_legend = show_legend, var_color = var_color, var_categorical = var_categorical, var_numerical = var_numerical)
  if (level == "bar") plot <- DrawBarCountsPlot(data = data, var_numerical = var_numerical, var_categorical = var_categorical, var_color = var_color, color_palette = color_palette, show_legend = show_legend)
  return(plot)
}

SignificanceSymbol <- function(pvalue) {
  dplyr::case_when(is.na(pvalue) ~ NA_character_,
                   pvalue <= 0.001 ~ "***",
                   pvalue <= 0.01  ~ "**",
                   pvalue <= 0.05  ~ "*",
                   TRUE ~ "ns")
}

Compute_tTE_Significance <- function(merged, x_col = "class", target_col, group_col = NULL) {

  if(length(unique(merged[[target_col]])) != 2) return(NULL)  # skip if not 2 levels

  target_levels <- sort(unique(merged[[target_col]]))
  print(target_levels)
  if (!is.null(group_col)) {
    sig_table <- merged %>%
      dplyr::group_by(.data[[x_col]], .data[[group_col]]) %>%
      dplyr::summarise(
        p_value = tryCatch(stats::wilcox.test(tTE ~ .data[[target_col]])$p.value, error = function(e) NA_real_), .groups = "drop")
  } else {
    sig_table <- merged %>%
      dplyr::group_by(.data[[x_col]]) %>%
      dplyr::summarise(
        p_value = tryCatch(stats::wilcox.test(tTE ~ .data[[target_col]])$p.value, error = function(e) NA_real_),
        .groups = "drop"
      )
  }

  if (nrow(sig_table) > 0) sig_table <- sig_table %>%
    dplyr::mutate(p_signif = SignificanceSymbol(.data$p_value), group1 = target_levels[1], group2 = target_levels[2])

  return(sig_table)
}



