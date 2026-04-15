#' @importFrom rlang :=

GetSafeColorScale <- function(n_colors, color_palette = NULL, aes_type = "fill") {

  use_custom <- !is.null(color_palette) && (length(color_palette) >= n_colors)
  use_internal <- n_colors <= length(gradual_groups_35)

  if (use_custom) {
    if (aes_type == "fill") return(ggplot2::scale_fill_manual(values = color_palette))
    if (aes_type == "color") return(ggplot2::scale_color_manual(values = color_palette))
  }

  if (use_internal) {
    if (!is.null(color_palette)) warning("The provided 'color_palette' does not have enough colors for ", n_colors, " categories. Defaulting to internal palette.")
    if (aes_type == "fill") return(ggplot2::scale_fill_manual(values = gradual_groups_35))
    if (aes_type == "color") return(ggplot2::scale_color_manual(values = gradual_groups_35))
  }

  if (!is.null(color_palette)) warning("The provided 'color_palette' does not have enough colors for ", n_colors, " categories. Defaulting to viridis palette.")

  if (aes_type == "fill") return(ggplot2::scale_fill_viridis_d(option = "viridis"))
  if (aes_type == "color") return(ggplot2::scale_color_viridis_d(option = "viridis"))
}

SavePlot <- function(plot, save_format, out_name, out_directory, width = 8, height = 6){

  ###
  # CALL: Each plot generating function
  # DESCRIPTION: This function is called every time the user specifies that the plot needs to be exported into a .png or .pdf format.
  ###

  save_format <- tolower(save_format)

  if (save_format %in% c("png", "pdf")){
    output_file <- GetOutputName(action = "plot", out_name = out_name, out_directory = out_directory, save_format = save_format)
    ggplot2::ggsave(filename = output_file, plot = plot, width = width, height = height)
  } else {
    warning("The input 'save_format' ('", save_format, "') is not recognized.\n",
            "Supported formats: png or pdf.\n", "The plot was generated but not stored.")
  }
}

GetOutputName <- function(action, out_name, out_directory, save_format){

  ###
  # CALL: SavePlot().
  # DESCRIPTION: This function is used to save as .pdf and .png plots and matrices as .rds files.
  # Requires an out_name and out_directory but if there are not given default values are used.
  ###

  if (is.null(out_directory)){
    message("- Parameter 'out_directory' has not been specified.\n  The file will be stored in the current working directory.")
    out_directory <- getwd()
  }

  if (is.null(out_name)){
    message("- Parameter 'out_name' has not been specified.\n  A standard name will be used.")
    out_name <- ifelse(action == "plot", "distribution_plot", "tRNA_expression_matrix")
  }

  # Check if the output name already contains the format extension
  if (!endsWith(tolower(out_name), paste0(".", tolower(save_format)))) out_name <- paste0(out_name, ".", save_format)

  output_file <- file.path(out_directory, out_name)
  message(paste("- The generated", action, "will be save in:", output_file))

  return (output_file) # Returns the output name.
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
  labels_df <- sig_data %>% dplyr::mutate(label = ifelse(.data$sig_adj, paste("p =", signif(.data$p_val_adj, 2), "*"), paste("p =", signif(.data$p_val_adj, 2))))

  # annot_data <- permut_data %>% dplyr::left_join(sig_data, by = c("group", "codon"))
  # annot_data <- annot_data %>% dplyr::mutate(label = ifelse(.data$sig_adj, paste("p =", signif(.data$p_val_adj, 2), "*"), paste("p =", signif(.data$p_val_adj, 2))))

  # Generate the anot_data_labels
  annot_data_labels <- labels_df %>%
    dplyr::group_by(.data$codon) %>%
    dplyr::summarise(combined_label = paste(.data$label, collapse = " | "), .groups = "drop")
  # annot_data_labels <- annot_data_labels %>% dplyr::distinct(.data$codon, .data$combined_label)

  annot_data <- permut_data %>% dplyr::left_join(labels_df, by = c("group", "codon"))

  return(list(annot_data, annot_data_labels))
}

CheckValueInData <- function(param, observed, expected){

  ###
  # CALL: Multiple
  # DESCRIPTION: This function takes the input data and variables to check that all of them are actually present in the data.
  # It does not have a return value except if it reports and error.
  ###

  if (is.null(observed)) return()

  present <- all(observed %in% expected)
  if (!present) {
    missing_vals <- observed[!observed %in% expected]

    stop("Invalid '", param, "' provided: '", paste(missing_vals, collapse = "', '"), "'.\n",
         "Available columns/options are: '", paste(expected, collapse = "', '"), "'.")
  }
}

CheckDataInLongFormat <- function(data) {

  ###
  # CALL: Multiple
  # DESCRIPTION: This function checks if the input data is in long format. If the data that is not in long format will not be computed.
  ###

  if (!is.data.frame(data)) stop("The input 'data' must be a data.frame or tibble.")

  num_cols <- vapply(data, is.numeric, FUN.VALUE = logical(1)) # One numeric column
  cat_cols <- vapply(data, function(x) is.character(x) || is.factor(x), FUN.VALUE = logical(1)) # One categorical column

  is_long <- any(cat_cols) & any(num_cols) # Will return TRUE if the data is in long format
  if(!is_long) stop("The input data needs to be in long format, containing at least one numeric and one categorical colummn.")
}

GenerateDistPlot <- function(level, target, data, x_axis_col, y_axis_col, color_palette, add_titles, show_legend, ncols, bar_position, facet_col, add_stats){

  if (level %in% c("jitter", "dot")) {
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], color = .data[[target]]))
  } else {
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], fill = .data[[target]]))
  }

  if (level == "jitter") plot <- plot + ggplot2::geom_jitter(size = 0.5)
  if (level == "barplot") plot <- plot + ggplot2::geom_col(stat = "identity", position = bar_position)  # Previously used geom_bar instead of geom_col
  if (level == "boxplot") plot <- plot + ggplot2::geom_boxplot()
  if (level == "dot") plot <- plot + ggplot2::aes(group = .data[[target]]) + ggplot2::geom_point()

  plot <- plot + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = show_legend)

  if (!is.null(facet_col)) plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_col]]), ncol = ncols) # stats::as.formula(paste("~", facet_col))
  if (level == "dot") plot <- plot + ggplot2::theme(strip.text = ggplot2::element_blank())

  n_colors <- length(unique(data[[target]]))
  current_aes <- if (level %in% c("jitter", "dot")) "color" else "fill"
  plot <- plot + GetSafeColorScale(n_colors = n_colors, color_palette = color_palette, aes_type = current_aes) # Customize the colors

  if(isTRUE(add_titles)) plot <- plot + ggplot2::labs(x = x_axis_col, y = paste(x_axis_col, "usage counts")) # + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 15))

  sig_table <- NULL
  if (isTRUE(add_stats)) { # Statistical Logic
    sig_table <- Compute_Boxplot_Significance(data, x_col = x_axis_col, y_col = y_axis_col, target_col = target, group_col = facet_col)

    if (!is.null(sig_table) && nrow(sig_table) > 0) {
      y_max <- max(data[[y_axis_col]], na.rm = TRUE) * 1.05 # Position starts slightly above the highest data point
      plot <- plot + ggpubr::stat_pvalue_manual(sig_table, x = x_axis_col, label = "p_signif", y.position = y_max, hide.ns = TRUE)
    }
  }

  if (isTRUE(add_stats)) return(list(plot = plot, stats = sig_table)) else return(list(plot = plot)) # Return results
}

DrawBarCountsPlot <- function(data, var_numerical, var_categorical, var_color, color_palette, show_legend, order, x_limits, facet_col) {

  # Set a specific order of the features (if available)
  data[[var_categorical]] <- if (!is.null(order)) factor(data[[var_categorical]], levels = order) else stats::reorder(data[[var_categorical]], data[[var_numerical]])

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[var_numerical]], y = .data[[var_categorical]], fill = .data[[var_color]])) + ggplot2::geom_col(color = "black") # + aes_order

  if (!is.null(x_limits)) plot <- plot + ggplot2::scale_x_continuous(limits = x_limits)
  if (!is.null(facet_col)) plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_col]]))

  n_colors <- length(unique(data[[var_color]]))
  plot <- plot + GetSafeColorScale(n_colors = n_colors, color_palette = color_palette, aes_type = "fill")
  plot <- plot + ggplot2::theme_bw() + ggplot2::theme(legend.position = show_legend) + ggplot2::labs(x = "", y = "")
  return(plot)
}

DrawDonutPlot <- function(data, var_numerical, var_categorical, color_palette, show_legend){

  # Compute percentages
  total <- sum(data[[var_numerical]], na.rm = TRUE)
  if (total == 0) stop("Error: Sum of '", var_numerical, "' is zero, cannot compute fractions.")

  data$fraction <- data[[var_numerical]] / total

  # Compute the cumulative percentages and the position of the label
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, utils::head(data$ymax, -1))
  data$labelPosition <- (data$ymax + data$ymin) / 2

  if (show_legend == "none") data$label <- paste(data[[var_categorical]], ":", round(data[[var_numerical]], 2)) else data$label <- as.character(round(data[[var_numerical]], 2))

  # Adjusting the sizes automatically based on the number of categories
  n_cat <- length(unique(data[[var_categorical]]))
  circle_min <- 1 + log10(n_cat) * 0.5 # Radius scales with log of categories
  donut_width <- ifelse(n_cat <= 3, 2.5, 1.5)
  circle_max <- circle_min + donut_width # before + 1
  label_size <- min(5, max(1.5, 6 - 0.1 * n_cat)) # Label size decreases as n_cat grows
  label_position <- 4 + 0.2 * n_cat # Label position shifts outward with n_cat

  plot <- ggplot2::ggplot(data, ggplot2::aes(ymax = .data$ymax, ymin = .data$ymin, xmax = circle_max, xmin = circle_min, fill = .data[[var_categorical]])) + ggplot2::geom_rect() +
    ggplot2::geom_label(x = label_position, ggplot2::aes(y = .data$labelPosition, label = .data$label), size = label_size, nudge_x = 2, nudge_y = 0) +
    ggplot2::coord_polar(theta = "y") + ggplot2::xlim(c(-1, label_position)) + ggplot2::theme_void()

  plot <- plot + GetSafeColorScale(n_colors = n_cat, color_palette = color_palette, aes_type = "fill")
  plot <- plot + ggplot2::theme(legend.position = show_legend)
  return(plot)
}

ComputeRings <- function(data, var_color, var_categorical, var_numerical, normalize, n_rings, zoom) {

  # Standardize dummy groups if needed
  if (is.null(var_color) || var_color == var_categorical) {
    data$.__dummy_group__ <- "all"
    var_color <- ".__dummy_group__"
  }

  if (!is.numeric(data[[var_numerical]])) data[[var_numerical]] <- as.numeric(as.character(data[[var_numerical]])) # Ensure numerical column is numeric

  # Pivot data to get wide format for ggradar
  tmp <- data %>%
    dplyr::select(dplyr::all_of(c(var_color, var_categorical, var_numerical))) %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of(var_color), names_from = dplyr::all_of(var_categorical),
                       values_from = dplyr::all_of(var_numerical), values_fn = sum, values_fill = 0)

  num_cols <- setdiff(colnames(tmp), var_color)

  if (isTRUE(normalize)) { # Apply normalization if requested
    row_totals <- rowSums(tmp[, num_cols, drop = FALSE], na.rm = TRUE) # Ensure row_totals is a vector for division
    row_totals[row_totals == 0] <- 1 # Prevent division by zero
    tmp[num_cols] <- tmp[num_cols] / row_totals
  }

  all_numeric_values <- unlist(tmp[num_cols]) # Extract all numeric values to calculate distribution

  if (isTRUE(zoom)) {
    max_val <- stats::quantile(all_numeric_values, 0.95, na.rm = TRUE) # Calculate the 95th percentile
    max_val <- as.numeric(max_val) + 0.00001 # Apply a tiny buffer to avoid "Value > grid.max" floating point errors
  } else {
    max_val <- if (isTRUE(normalize)) 1.00001 else max(all_numeric_values, na.rm = TRUE) + 0.00001
  }

  ring_values <- seq(0, max_val, length.out = n_rings) # Generate equally spaced ring values
  labels_rings <- if (isTRUE(normalize)) paste0(round(ring_values * 100, 1), "%") else as.character(round(ring_values, 2)) # Create labels based on normalization

  return(list(max_val = max_val, labels_rings = labels_rings))
}

DrawRadarPlot <- function(data, var_color, var_categorical, var_numerical, normalize, zoom, add_titles, title, show_legend, global_max_val = NULL, labels_rings = NULL, color_palette) {

  if (is.null(var_color) || var_color == var_categorical) {
    data$.__dummy_group__ <- "all"
    var_color <- ".__dummy_group__."
  }

  if (!is.numeric(data[[var_numerical]])) data[[var_numerical]] <- as.numeric(as.character(data[[var_numerical]]))

  p <- data %>%
    dplyr::select(dplyr::all_of(c(var_color, var_categorical, var_numerical))) %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of(var_color), names_from = dplyr::all_of(var_categorical), values_from = dplyr::all_of(var_numerical), values_fn = sum, values_fill = 0)

  if (isTRUE(normalize)) {
    num_cols <- setdiff(colnames(p), var_color)
    p[num_cols] <- p[num_cols] / rowSums(p[num_cols], na.rm = TRUE)
  }

  num_cols <- setdiff(colnames(p), var_color)
  p[num_cols] <- lapply(p[num_cols], function(x) ifelse(x > global_max_val, global_max_val, x))

  n_groups <- nrow(p)
  plot_colors <- color_palette
  if (!is.null(color_palette)) {
    if (length(color_palette) < n_groups) {
      warning("Not enough colors in 'color_palette' for radar plot. Defaulting to ggradar's internal palette.")
      plot_colors <- NULL
    } else {
      plot_colors <- color_palette[1:n_groups]
    }
  }

  radar_plot <- ggradar::ggradar(p, background.circle.colour = "white", group.colours = plot_colors, values.radar = labels_rings, grid.min = 0, grid.mid = global_max_val/2, grid.max = global_max_val,
                                 gridline.min.linetype = 1, gridline.mid.linetype = 1, gridline.max.linetype = 1, group.line.width = 0.8, group.point.size = 2,
                                 legend.position = show_legend)

  if (isTRUE(add_titles)) radar_plot <- radar_plot + ggplot2::ggtitle(title)

  return(radar_plot)
}

GenerateProportionPlot <- function(level, data, var_numerical, var_categorical, var_color, color_palette, show_legend, order, x_limits, facet_col, normalize, zoom, add_titles, n_rings){

  if (level == "donut") {
    plot <- DrawDonutPlot(data = data, var_numerical = var_numerical, var_categorical = var_categorical, color_palette = color_palette, show_legend = show_legend)
  } else if (level == "bar") {
    plot <- DrawBarCountsPlot(data = data, var_numerical = var_numerical, var_categorical = var_categorical, var_color = var_color, color_palette = color_palette,
                              show_legend = show_legend, order = order, x_limits = x_limits, facet_col = facet_col)
  } else if (level == "radar") {
    conditions <- sort(unique(data[[var_categorical]])) # Replace long category names with numeric indices
    condition_index <- stats::setNames(seq_along(conditions), conditions)
    data[[var_categorical]] <- condition_index[data[[var_categorical]]]

    legend_table <- data.frame(Index = seq_along(conditions), Condition = conditions) # Save legend table

    rings_info <- ComputeRings(data = data, var_color = var_color, var_categorical = var_categorical, var_numerical = var_numerical, normalize, n_rings = n_rings, zoom = zoom) # Compute global rings once

    if(is.null(facet_col)){
      title_plot <- paste("Radar plot by", var_categorical)
      plot <- DrawRadarPlot(data = data, show_legend = show_legend, var_color = var_color, var_categorical = var_categorical, var_numerical = var_numerical,
                            normalize = normalize, zoom = zoom, add_titles = add_titles, title = title_plot, global_max_val = rings_info$max_val,
                            labels_rings = rings_info$labels_rings, color_palette = color_palette)
    } else {
      data_split <- split(data, data[[facet_col]])
      plots <- list()
      for (group_name in names(data_split)) {
        data_subset <- data_split[[group_name]]
        plots[[group_name]] <- DrawRadarPlot(data = data_subset, var_color = var_color, var_categorical = var_categorical, var_numerical = var_numerical,
                                             normalize = normalize, zoom = zoom, add_titles = add_titles, title = group_name, show_legend = show_legend, global_max_val = rings_info$max_val,
                                             labels_rings = rings_info$labels_rings, color_palette = color_palette)
      }
      plot <- patchwork::wrap_plots(plots)
    }
    return(list(plot = plot, legend_table = legend_table))
  }
  return(list(plot = plot))
}

SignificanceSymbol <- function(pvalue) {

  if (!is.numeric(pvalue)) return(rep(NA_character_, length(pvalue)))

  symbols <- as.character(stats::symnum(pvalue, corr = FALSE, na = FALSE, cutpoints = c(-Inf, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))
  symbols[is.na(pvalue)] <- NA_character_

  return(symbols)
}

Compute_tTE_Significance <- function(merged, x_col = "class", target_col, group_col = NULL) {

  target_levels <- sort(unique(merged[[target_col]]))
  if(length(target_levels) < 2) return(NULL) # We need at least 2 levels to compare anything
  comparisons <- utils::combn(target_levels, 2) # Generate all possible pairs: 2 columns, N rows

  group_vars <- if (!is.null(group_col)) group_col else NULL # Define the grouping variables (facets only)

  run_pairwise_stats <- function(df) {
    purrr::map_df(seq_len(ncol(comparisons)), function(i) { # Loop through each pair in the 'comparisons' matrix
      g1_name <- comparisons[1, i]
      g2_name <- comparisons[2, i]

      p_val <- tryCatch({
        val1 <- df$tTE[df[[target_col]] == g1_name]
        val2 <- df$tTE[df[[target_col]] == g2_name]
        stats::wilcox.test(val1, val2)$p.value
      }, error = function(e) NA_real_)

      data.frame(group1 = g1_name, group2 = g2_name, p_value = p_val, comparison = paste0(g1_name, "_vs_", g2_name))
    })
  }

  if (is.null(group_vars)) {
    sig_table <- run_pairwise_stats(merged)
  } else {
    sig_table <- merged %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::group_modify(~ run_pairwise_stats(.x)) %>%
      dplyr::ungroup()
  }

  if (nrow(sig_table) > 0) { # Final formatting
    sig_table <- sig_table %>%
      dplyr::mutate(p_signif = SignificanceSymbol(.data$p_value), !!x_col := .data$group2) %>%
      dplyr::filter(!is.na(.data$p_value))
  }

  return(sig_table)
}

Compute_tTE_Significance_simple <- function(merged, x_col = "class", target_col, group_col = NULL) {

  target_levels <- sort(unique(merged[[target_col]]))
  if (length(target_levels) != 2) return(NULL)

  group_vars <- if (!is.null(group_col)) group_col else NULL

  # If there are no group_vars, we run a single test on the whole dataset
  if (is.null(group_vars)) {
    sig_table <- data.frame(temp_id = 1) %>%
      dplyr::mutate(
        p_value = tryCatch({
          g1 <- merged$tTE[merged[[target_col]] == target_levels[1]]
          g2 <- merged$tTE[merged[[target_col]] == target_levels[2]]
          stats::wilcox.test(g1, g2)$p.value
        }, error = function(e) NA_real_)
      )
  } else {
    # If we have facets, group by those facets only
    sig_table <- merged %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        p_value = tryCatch({
          g1 <- .data$tTE[.data[[target_col]] == target_levels[1]]
          g2 <- .data$tTE[.data[[target_col]] == target_levels[2]]
          stats::wilcox.test(g1, g2)$p.value
        }, error = function(e) NA_real_),
        .groups = "drop")
  }

  # Add columns back for the plot to recognize where to put the stars
  if (nrow(sig_table) > 0) {
    sig_table <- sig_table %>%
      dplyr::mutate(p_signif = SignificanceSymbol(.data$p_value), group1 = target_levels[1], group2 = target_levels[2],
                    comparison = paste0(target_levels[1], "_vs_", target_levels[2]))
  }

  return(sig_table)
}

Compute_Boxplot_Significance <- function(merged, x_col, y_col, target_col, group_col = NULL) {

  if(length(unique(merged[[target_col]])) != 2) return(NULL)
  target_levels <- sort(unique(merged[[target_col]]))

  form <- stats::as.formula(paste(y_col, "~", target_col)) # Create formula dynamically for the test

  sig_table <- merged %>% # Perform the Wilcoxon test
    dplyr::group_by(dplyr::across(dplyr::all_of(c(x_col, group_col)))) %>%
    dplyr::summarise(
      p_value = tryCatch({
        stats::wilcox.test(form, data = dplyr::pick(dplyr::everything()))$p.value
      }, error = function(e) NA_real_), .groups = "drop")

  sig_table <- sig_table %>% dplyr::filter(!is.na(.data$p_value)) # Filter out NA values (failed tests)

  if (nrow(sig_table) > 0){
    sig_table <- sig_table %>%
      dplyr::mutate(p_value = stats::p.adjust(.data$p_value, method = "BH")) %>%
      dplyr::mutate(p_signif = SignificanceSymbol(.data$p_value), group1 = as.character(target_levels[1]), group2 = as.character(target_levels[2]))
  } else {
    return(NULL) # All tests failed
  }

  return(sig_table)
}
