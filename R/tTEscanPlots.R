# All the functions in this file are called by the user

#' Distribution Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply
#' @description
#' This function generates a visualization of codon-anticodon usage or amino acid demand-supply distributions across conditions.
#' The input data (\code{data}) is expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' Based on the selected \code{plot} parameter the user can select the most convenient layout to display the data.
#' It is crucial to ensure consistency between the name of the columns in \code{data} and the parameters \code{x_axis_col}, \code{y_axis_col} and \code{condition_col}.
#'
#' @param data A long format table. This format can be obtained using \code{\link{TransformFormat}}.
#' @param plot Either \code{"jitter"} (default), \code{"barplot"}, \code{"boxplot"} or \code{"dot"} to indicate the type of plot to generate.
#' @param bar_position Either \code{"dodge"} (default), \code{"stack"} or \code{"fill"} to indicate the display of the plot if \code{plot} is \code{"barplot"}.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param condition_col Name of the categorical variable to group the data points.
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param targeted_arg Optional; a vector defining key feature clusters to highlight or label.
#' @param ncols Numeric; number of columns for arranging panels. Defaults to 1.
#' @param facet_col Optional; name of the categorical variable to divide the plot into different panels. Required if \code{ncols} bigger than 1.
#' @param add_stats Logical; if \code{TRUE}, performs a statistical analysis based on the available parameters. Defaults to \code{FALSE}.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested distribution plot. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_DistributionPlot <- function(data, plot = "jitter", bar_position = "dodge", x_axis_col, y_axis_col, condition_col, ncols = 1, facet_col = NULL, add_stats = FALSE,
                                 color_palette = NULL, targeted_arg = NULL, save_format = NULL, out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE) {

  ###
  # DESCRIPTION: This function takes a long format data expressing the usage of different features across a set of conditions.
  ###

  # Checking the input parameters
  CheckDataInLongFormat(data) # Data in long format.
  CheckValueInData(param = "plot", observed = plot, expected = c("jitter", "boxplot", "barplot", "dot")) # Valid plot layout
  CheckValueInData(param = "_col", observed = c(x_axis_col, y_axis_col, condition_col), expected = colnames(data)) # All columns present in data

  if (!is.null(facet_col)){
    CheckValueInData(param = "_col", observed = c(facet_col), expected = colnames(data))
    if (is.numeric(data[[facet_col]])) stop(paste("Error in 'data': Column", facet_col, "is not categorical."))
  }

  if (!show_legend %in% c("none", "top", "bottom", "right", "left")) stop("Please specify a suitable 'show_legend' parameter.\nSupported formats: none, top, bottom, right or left.")
  if (!(is.numeric(data[[y_axis_col]]))) stop(paste("Error in 'data': Column", y_axis_col, "is not numeric."))
  if (is.numeric(data[[x_axis_col]])) stop(paste("Error in 'data': Column", x_axis_col, "is not categorical."))
  if (is.numeric(data[[condition_col]])) stop(paste("Error in 'data': Column", condition_col, "is not categorical."))
  if (!bar_position %in% c("stack", "dodge", "fill")) stop("Please specify a suitable 'bar_position' parameter.\nSupported formats: stack, dodge or fill.")

  target_col_name <- condition_col

  # TARGETED APPROACH - compare specific conditions against the rest
  if (!(is.null(targeted_arg))){
    data$target <- "other" # Initialize all the values of the target column as other

    for (arg in targeted_arg){ # Iterate over the values in targeted_arg
      index <- which(data[[condition_col]] == arg) # Check which entry match the value in targeted_arg
      if (length(index) == 0) stop(paste("The argument", arg, "has not been found in the categories defined in the data."))
      data$target_group[index] <- arg # Update the value of the target column
    }
    target_col_name <- "target_group"
    CheckDataFrame(data = data) # Ensure that the update data is suitable for the downstream analysis
  }

  # Adjust the colors of the plot
  if (is.list(color_palette)) {
    warning("The input 'color_palette' needs to be a vector.\n",
            "The default 'color_palette' will be used instead.")
    color_palette <- NULL
  }

  num_target_categories <- length(unique(data[[target_col_name]]))

  if (!(is.null(color_palette)) && (length(color_palette) < num_target_categories)){
    stop("The number of colors in 'color_palette' does not correspond to the number of categories needed.\n",
         "Note: if the 'targeted_arg' is given, an extra color is required.\n",
         paste("Number of colors given:", length(color_palette), "\n"), paste("Number of colors needed:", num_target_categories))
  }

  # Customization
  plot <- GenerateDistPlot(level = plot, target = target_col_name, data = data, x_axis_col = x_axis_col, y_axis_col = y_axis_col, color_palette = color_palette,
                           add_titles = add_titles, show_legend = show_legend, ncols = ncols, facet_col = facet_col, bar_position = bar_position, add_stats = add_stats)
  if (!(is.null(save_format))) SavePlot(plot = plot$plot, save_format = save_format, out_name = out_name, out_directory = out_directory)

  return(plot) # Returns the generated plot and exports a file if enabled.
}

#' Distribution Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply Compared to a Target
#' @description
#' This function generates a distribution plot to compare distribution usages between a targeted set of conditions against the rest.
#' A common application of this plot is to compare mean usages across conditions.
#' Both input data sources (\code{target_data} and \code{overall_data}) are expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' It is crucial to ensure consistency between the name of the columns in \code{target_data} and \code{overall_data}, together with the parameters \code{x_axis_col} and \code{y_axis_col}.
#'
#' @param target_data A long format table of a condition of interest. This format can be obtained using \code{\link{TransformFormat}}.
#' @param overall_data A long format table of the whole set of conditions (with or without the data in \code{target_data}. This format can be obtained using \code{\link{TransformFormat}}.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param show_difference Logical; if \code{TRUE}, displays the difference between \code{overall_data} and \code{target_data}.
#' @param color_palette Optional; a vector of color codes (min. 2 codes and max. 3 codes) to customize plot appearance. Colors for (i) data, (ii) target value, and (iii) difference bar.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested distribution plot. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_PlotTargetComparison <- function(target_data, overall_data, x_axis_col, y_axis_col, show_difference = TRUE, color_palette = NULL,
                                     save_format = NULL, out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE){

  ###
  # DESCRIPTION: This function takes two long format data and expressing the usage of different features across a set of conditions and compares them.
  ###

  # Check that the names given correspond to column names.
  CheckValueInData(param = "x_axis_col", observed = x_axis_col, expected = colnames(overall_data)) # features
  CheckValueInData(param = "y_axis_col", observed = y_axis_col, expected = colnames(overall_data)) # values
  CheckValueInData(param = "x_axis_col", observed = x_axis_col, expected = colnames(target_data))
  CheckValueInData(param = "y_axis_col", observed = y_axis_col, expected = colnames(target_data))

  if (!(is.numeric(overall_data[[y_axis_col]]))) stop(paste("Error in 'overall_data': Column", y_axis_col, "is not numeric."))
  if (!(is.numeric(target_data[[y_axis_col]]))) stop(paste("Error in 'target_data': Column", y_axis_col, "is not numeric."))
  if (is.numeric(overall_data[[x_axis_col]])) stop(paste("Error in 'overall_data': Column", x_axis_col, "is not categorical."))
  if (is.numeric(target_data[[x_axis_col]])) stop(paste("Error in 'target_data': Column", x_axis_col, "is not categorical."))

  target_subset <- target_data[, c(x_axis_col, y_axis_col), drop = FALSE]
  colnames(target_subset)[2] <- "target_value"

  plot_data <- merge(overall_data, target_subset, by = x_axis_col, all.x = TRUE)

  # Adjust the colors of the plot
  if (is.list(color_palette)) {
    warning("The input 'color_palette' needs to be a vector.\n",
            "The default 'color_palette' will be used instead.")
    color_palette <- NULL
  }

  if (!is.null(color_palette) && ((isFALSE(show_difference) && length(color_palette) < 2) || (isTRUE(show_difference) && length(color_palette) < 3))) {
    colors_needed <- ifelse(isTRUE(show_difference), 3, 2)

    stop("The number of colors in 'color_palette' does not correspond to the number of graphical elements.\n",
         "Note: If 'show_difference' is TRUE, an extra color is required for the segments.\n",
         "Number of colors given: ", length(color_palette), "\n",
         "Number of colors needed: ", colors_needed)
  } else {
    color_palette <- c("#bababa", "#b2182b", "#fddbc7") # Use default color_palette
  }

  # Generate plot
  plot <- ggplot2::ggplot(data = plot_data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]])) +
    ggplot2::geom_col(fill = color_palette[1]) +
    ggplot2::theme_bw()

  if (isTRUE(show_difference)){
    plot <- plot + ggplot2::geom_segment(mapping = ggplot2::aes(xend = .data[[x_axis_col]], yend = .data$target_value), color = color_palette[3]) # Include difference bars between the targeted point and the rest of the data
  }

  plot <- plot + ggplot2::geom_point(mapping = ggplot2::aes(y = .data$target_value), color = color_palette[2], size = 2) +
    ggplot2::theme(ggplot2::labs(x = x_axis_col, y = y_axis_col), axis.text.x = ggplot2::element_text(angle = 90), legend.position = show_legend) # Arrange elements in the plot
  if (!(is.null(save_format))) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory) # Save the ggplot

  return(plot) # Returns the generated plot and exports a file if enabled.
}

#' Proportion Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply
#' @description
#' This function generates a proportion plot to compare codon-anticodon usage or amino acid demand-supply frequencies across conditions.
#' The input data (\code{data}) is expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' Based on the selected \code{plot} parameter the user can select the most convenient layout to display the data.
#' It is crucial to ensure consistency between the name of the columns in \code{data} and the parameters \code{var_numerical}, \code{var_categorical} and \code{var_color}.
#'
#' @param data A long format table. This format can be obtained using \code{\link{TransformFormat}}.
#' @param plot Either \code{"bar"} (default), \code{"donut"} or \code{"radar"} to indicate the type of plot to generate.
#' @param var_numerical Name of the numerical variable to reflect in the plot.
#' @param var_categorical Name of the categorical variable to reflect in the plot.
#' @param var_color Optional; name of the categorical variable to group the data points when coloring. Required if \code{plot} is \code{"bar"}.
#' @param facet_col Optional; name of the categorical variable to divide the plot into different panels.
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param num_limits Optional; a vector with the upper and lower ranges of the values in \code{var_numerical}.
#' @param num_rings Optional; a number specifying the amount of rings to display if \code{plot} is \code{"radar"}. Defaults to 5.
#' @param order Optional; a vector of the levels to organize the data, based on the \code{var_categorical}.
#' @param normalize Logical; if \code{TRUE}, normalizes the data in order to display the relative contribution of the \code{var_categorical}.  Defaults to \code{TRUE}.
#' @param zoom Logical; if \code{TRUE}, centers the plot display to the values in \code{data}.Defaults to \code{TRUE}.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested proportion plot. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_ProportionPlot <- function(data, plot = "bar", var_numerical, var_categorical, var_color = NULL, facet_col = NULL, color_palette = NULL, num_limits = NULL, num_rings = 5,
                               save_format = NULL,  out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE, order = NULL, normalize = TRUE, zoom = FALSE){

  ###
  # DESCRIPTION: This function takes long format table to compare proportions of each condition with respect to the total values of the data.
  ###

  # Checking the type of plot and the input variables
  CheckValueInData(param = "plot", observed = plot, expected = c("donut", "radar", "bar"))
  if (is.null(var_color)) var_color <- var_categorical

  CheckValueInData(param = "var_", observed = c(var_categorical, var_numerical, var_color), expected = colnames(data))
  if (!is.null(facet_col)) CheckValueInData(param = "_col", observed = c(facet_col), expected = colnames(data))

  if (!(is.numeric(data[[var_numerical]]))) stop(paste("Error in 'data': Column", var_numerical, "is not numeric."))
  if (is.numeric(data[[var_categorical]])) stop(paste("Error in 'data': Column", var_categorical, "is not categorical."))
  if (is.numeric(data[[var_color]])) stop(paste("Error in 'data': Column", var_color, "is not categorical."))
  if (!is.null(facet_col) && is.numeric(data[[facet_col]])) stop(paste("Error in 'data': Column", facet_col, "is not categorical."))

  # Generate plot
  plot <- GenerateProportionPlot(data = data, level = plot, var_numerical = var_numerical, var_categorical = var_categorical, var_color = var_color,
                                 color_palette = color_palette, show_legend = show_legend, order = order, x_limits = num_limits, facet_col = facet_col,
                                 normalize = normalize, zoom = zoom, add_titles = add_titles, n_rings = num_rings)

  # Save the ggplot
  if(!is.null(save_format)) SavePlot(plot = plot$plot, save_format = save_format, out_name = out_name, out_directory = out_directory)
  return(plot) # Returns the generated plot and exports a file if enabled.
}

#' Violin Plot Displaying tTE Scores
#' @description
#' This function generates a violin plot to visualize the distribution of tTE scores across different condition. The tTE scores should be obtained using \code{\link{Compute_tTE}}.
#' The plot allows the user to easily compare the tTE distribution between different conditions, with the option to highlight specific feature clusters based on the provided \code{targets}.
#'
#' @param data A tTE results table obtained from \code{\link{Compute_tTE}}.
#' @param metadata A table with additional information regarding the conditions in \code{data}.
#' @param index_col Name of the categorical variable that links the conditions in \code{data} with the \code{metadata}.
#' @param class_col Name of the categorical variable to reflect in the plot.
#' @param add_stats Logical; if \code{TRUE}, performs a statistical analysis based on the available parameters. Defaults to \code{TRUE}.
#' @param target_col Optional; name of the categorical variable to perform the statistical comparison (the most specific level). Used if \code{add_stats} is \code{TRUE}.
#' @param facet_col Optional; name of the categorical variable separate the plot into different panels.
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend in the plot.
#' @param show_outliers Logical; if \code{TRUE}, labels the outlier data points based on \code{index_col}.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the tTE scores. If \code{save_format} is provided, the plot will also be saved to the specified location. If \code{add_stats} reports a table with the statitical measures summarized.
#' @export

tTE_ScoresPlot <- function(data, metadata, class_col, index_col, target_col = NULL, facet_col = NULL, color_palette = NULL, save_format = NULL,
                           out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE, add_stats = TRUE, show_outliers = FALSE) {

  required_cols <- c("condition", "tTE", "p_value", "neg_log10_tTE_p_value") # Extracted when computing the tTE through tTEscanR
  CheckValueInData("data columns", colnames(data), required_cols) # Modified CheckValueInData handles this!
  CheckValueInData("index_col", index_col, colnames(metadata))
  CheckValueInData("class_col", class_col, colnames(metadata))

  if (isFALSE(add_stats) && !is.null(target_col)){
    warning("The 'target_col' will be ignored as 'add_stats' is not enabled.\n Check the documentation of the function for more details.")
    add_stats <- NULL
  }
  if (!is.null(target_col)) CheckValueInData("target_col", target_col, colnames(metadata))
  if (!is.null(facet_col)) CheckValueInData("facet_col", facet_col, colnames(metadata))

  merged <- dplyr::left_join(data, metadata, by = stats::setNames(index_col, "condition"))
  if (nrow(merged) == 0) stop("Error: Merged data has 0 rows. Check if 'index_col' matches 'data$condition'.")

  merged$class <- merged[[class_col]]
  na_classes <- sum(is.na(merged$class))
  if (na_classes > 0) warning(paste(na_classes, "rows in data could not be matched to a group in metadata!"))

  group_counts <- table(merged$class)
  if (any(group_counts < 2)) warning("One or more groups have fewer than 2 samples. Statistics (p-values) will return NA.")

  merged$label <- merged[["condition"]]
  merged$target <- if (!is.null(target_col)) merged[[target_col]] else merged$class

  if (isTRUE(show_outliers)){ # Label outliers
    group_vars <- c("class", "target")
    if (!is.null(facet_col)) group_vars <- c(group_vars, facet_col)

    merged <- merged %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::mutate(q1 = stats::quantile(.data$tTE, 0.25, na.rm = TRUE),
                    q3 = stats::quantile(.data$tTE, 0.75, na.rm = TRUE),
                    iqr = .data$q3 - .data$q1,
                    is_extreme = .data$tTE < (.data$q1 - 1.5 * .data$iqr) | .data$tTE > (.data$q3 + 1.5 * .data$iqr)) %>%
      dplyr::ungroup()
  }

  fill_label <- if (!is.null(target_col)) target_col else "Group"

  # plot <- ggplot2::ggplot(merged, ggplot2::aes(x = .data$class, y = .data$tTE, fill = .data$target)) +
  #   ggplot2::geom_violin(trim = FALSE, scale = "width", width = 0.6, position = ggplot2::position_dodge(width = 0.9), alpha = 0.7) +
  #   ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.08, dodge.width = 0.7), size = 0.8, alpha = 0.5) +
  #   ggplot2::theme_bw() +
  #   ggplot2::scale_x_discrete(expand = c(0.02, 0)) +
  #   ggplot2::labs(x = class_col, y = "tTE", fill = fill_label)

  plot <- ggplot2::ggplot(merged,
    ggplot2::aes(x = .data$class, y = .data$tTE, fill = .data$target)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", width = 0.7, position = ggplot2::position_dodge(width = 0.7) , alpha = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(dodge.width = 0.7, jitter.width = 0.03), size = 0.8, alpha = 0.5) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.4)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = class_col, y = "tTE", fill = fill_label)

  if (isTRUE(show_outliers)){ # Label outliers
    plot <- plot + ggrepel::geom_text_repel(data = dplyr::filter(merged, .data$is_extreme), ggplot2::aes(x = .data$class, y = .data$tTE, label = .data$label),
                                            position = ggplot2::position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15),
                                            size = 3, max.overlaps = 15, show.legend = FALSE)
  }

  if (length(unique(merged$class)) > 5) plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  if (!is.null(facet_col)) plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_col]]), scales = "free_x")

  n_colors <- length(unique(merged$target))
  plot <- plot + GetSafeColorScale(n_colors = n_colors, color_palette = color_palette, aes_type = "fill")

  if (add_titles) {
    title_text <- paste("tTE scores by", class_col)
    if (!is.null(target_col)) title_text <- paste(title_text, "highlighting", target_col)
    if (!is.null(facet_col)) title_text <- paste(title_text, "faceted by", facet_col)
    plot <- plot + ggplot2::labs(title = title_text)
  }

  plot <- plot + ggplot2::theme(legend.position = show_legend)

  sig_table <- NULL
  if (isTRUE(add_stats)){
    if (!is.null(target_col)) {
      stat_target <- target_col
      if (stat_target != class_col) {
        actual_group_col <- unique(c("class", facet_col))
        actual_x_col <- "dummy_x"
      } else {
        actual_group_col <- facet_col
        actual_x_col <- "class"
      }
    } else {
      stat_target <- class_col
      actual_group_col <- facet_col
      actual_x_col <- "class"
    }

    variances <- tapply(merged$tTE, merged[[stat_target]], stats::var, na.rm = TRUE)
    if (any(variances == 0, na.rm = TRUE) || any(is.na(variances))) warning("One or more groups have zero variance (all values are identical). Stats will return NA.")

    sig_table <- Compute_tTE_Significance(merged, x_col = actual_x_col, target_col = stat_target, group_col = actual_group_col)

    if (!is.null(sig_table) && nrow(sig_table) > 0) {
      y_max <- max(merged$tTE, na.rm = TRUE) * 1.15  # stars above violins

      if (stat_target != class_col) {
        plot <- plot + ggpubr::stat_pvalue_manual(
          sig_table, x = "class", y.position = y_max, label = "p_signif", size = 7, hide.ns = TRUE, step.increase = 0.12)
      } else {
        plot <- plot + ggpubr::stat_pvalue_manual(
          sig_table, y.position = y_max, label = "p_signif", size = 7, hide.ns = TRUE, step.increase = 0.12)
      }
    } else {
      message("The statistics could not be computed, check the selected comparisons.")
    }
  }

  if (!is.null(save_format)) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory)
  return(list(plot = plot, stats = sig_table))
}

#' Correlation Plot: Exonic Codon Background - Mean Codon Usage
#' @description
#' This function generates a visualization to correlate different parmaeters of the \code{data}.
#' The input data (\code{data}) is expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' A common usage of this kind of plot is to represent the codon frequencies correlations of the exonic background and the mean codon usage across conditions.
#' Based on the selected \code{plot} parameter the user can select the most convenient layout to display the data.
#' It is crucial to ensure consistency between the name of the columns in \code{data} and the parameters \code{x_axis_col}, \code{y_axis_col} and \code{condition_col}.
#'
#' @param data A long format table. This format can be obtained using \code{\link{TransformFormat}}.
#' @param plot Either \code{"MeanCodonUsage"} (default) or \code{"PoolDiversity"} to indicate the \code{data} source.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param condition_col Name of the categorical variable to group the data points.
#' @param label_col Name of the categorical variable to label the data points.
#' @param extra_val Optional; variable with additional information to include in the plot (e.g. correlation value).
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param targeted_arg Optional; a vector defining key feature clusters to highlight or label.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the correlation. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_CorrelationPlot <- function(data, plot = "MeanCodonUsage", x_axis_col, y_axis_col, condition_col, label_col = NULL, extra_val = NULL, color_palette = NULL,
                                targeted_arg = NULL, save_format = NULL, out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE) {

  ###
  # DESCRIPTION: This function takes a long format data with two parameters to correlate.
  # It is possible to perform a more targeted approach to compare specific conditions against the rest
  ###

  plot_type <- plot

  # Checking the input parameters
  if (!(plot_type %in% c("MeanCodonUsage", "PoolDiversity"))) stop("Please use MeanCodonUsage or PoolDiversity as 'plot'.") # Suitable plot
  CheckDataInLongFormat(data) # Data in long format.

  expected_cols <- c(x_axis_col, y_axis_col, condition_col, label_col)
  CheckValueInData("data columns", expected_cols, colnames(data))

  # if (is.null(label_col)){
  #   CheckValueInData(param = "_col", expected = colnames(data), observed = c(x_axis_col, y_axis_col, condition_col)) # Columns present in data
  # } else {
  #   CheckValueInData(param = "_col", expected = colnames(data), observed = c(x_axis_col, y_axis_col, condition_col, label_col)) # Columns present in data
  #   if (is.numeric(data[[label_col]])) stop(paste("Error in `data`: Column", label_col, "is not categorical."))
  # }

  if (!(is.numeric(data[[y_axis_col]]))) stop(paste("Error in 'data': Column", y_axis_col, "is not numeric."))
  if (!(is.numeric(data[[x_axis_col]]))) stop(paste("Error in 'data': Column", x_axis_col, "is not numeric"))
  if (!is.character(data[[condition_col]]) && !is.factor(data[[condition_col]])) stop(paste("Error in 'data': Column", condition_col, "is not categorical."))
  if (!is.null(label_col) && is.numeric(data[[label_col]])) stop(paste("Error in 'data': Column", label_col, "is not categorical."))

  plot_color_col <- condition_col

  # TARGETED APPROACH
  if (!is.null(targeted_arg)){
    data$target <- "other" # Initialize all the values of the target column as other

    missing_targets <- targeted_arg[!targeted_arg %in% data[[condition_col]]]
    if (length(missing_targets) > 0) stop("The following arguments in 'targeted_arg' were not found in the data: ", paste(missing_targets, collapse = ", "))

    data$.__target__ <- ifelse(data[[condition_col]] %in% targeted_arg, as.character(data[[condition_col]]), "other")
    plot_color_col <- ".__target__" # Update the pointer to use our new column

    # for (i in 1:length(targeted_arg)){ # Iterate over the values in targeted_arg
    #   index <- which(data[[condition_col]] == targeted_arg[i]) # Check which entry match the value in targeted_arg
    #   if (length(index) == 0) stop(paste("The argument", targeted_arg[i], "has not been found in the categories defined in the data."))
    #  data$target[index] <- targeted_arg[i] # Update the value of the target column
    # }
  }

  # Generate the plot
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], color = .data[[plot_color_col]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_bw()

  # Adjust the colors of the plot
  n_colors <- length(unique(data[[plot_color_col]]))
  plot <- plot + GetSafeColorScale(n_colors = n_colors, color_palette = color_palette, aes_type = "color")

  # if (is.null(label_col)){
  #   plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], color = .data[[target]]))
  # } else {
  #   plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], label = .data[[label_col]], color = .data[[target]]))
  # }

  if (plot_type == "MeanCodonUsage") {
    plot <- plot + ggplot2::geom_abline(slope = 1, intercept = 0)
    if (isTRUE(add_titles)) {
      plot <- plot + ggplot2::labs(title = "Codon frequencies of exonic background \n vs. mean codon usage across conditions",
                                   x = "Mean codon usage across conditions", y = "Codon usage of exonic background")
    }
  } else if (plot_type == "PoolDiversity" && isTRUE(add_titles)) {
    plot <- plot + ggplot2::labs(title = "Conditions' correlation to mean codon usage vs. \n codon pool diversity",
                                 x = "Codon pool diversity", y = "Correlation")
  }

  if (!is.null(label_col)) plot <- plot + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_col]]), max.overlaps = 10)

  if (!(is.null(extra_val))){
    plot <- plot + ggplot2::annotate("text", x = Inf, y = -Inf, label = sprintf("rho == %.3f", extra_val),
                                     parse = TRUE, hjust = 1.2, vjust = -0.2)
  }

  plot <- plot + ggplot2::theme(legend.position = show_legend)

  # Save the ggplot
  if(!is.null(save_format)) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory)
  return(plot) # Returns the generated plot and exports a file if enabled
}

#' Permutation Plot
#' @description
#' This function generates a plot to compare the baseline codon exonic background against the current codon usage.
#' For a better interpretation, the codons are colored by amino acid.
#'
#' @param permut_data A table with the codons and their frequencies after computing all the permutations. Output from \code{\link{GetPermutationDist}}.
#' @param sig_data A table with the codon exonic background and their significance level before (p-value) and after the correction (p-adjusted value). Output from \code{\link{ObtainSignificance}}
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend in the plot.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return Permutation plot.
#' @export

tTE_PermutationPlot <- function(permut_data, sig_data, color_palette = NULL, save_format = NULL, out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE){

  ###
  # DESCRIPTION: This function takes the output tables of GetPermutationDist() and ObtainSignificance().
  # Generates a histogram per each codon and displays the background (expected) and the observed signal colored by groups and amino acids.
  ###

  # Processing the annot_data parameter
  annot_data_results <- GetAnnotData(permut_data = permut_data, sig_data = sig_data)
  annot_data <- annot_data_results[[1]]
  annot_data_labels <- annot_data_results[[2]]

  # Defining the size of the bars in the histogram based on the data
  bin_width <- (max(annot_data$freq) - min(annot_data$freq)) / 25
  if (is.null(bin_width)) stop("The bind width could not be calculated.\n", "Revise the `freq` column in the data.")

  # Defining the plot background
  plot <- ggplot2::ggplot(permut_data, ggplot2::aes(x = .data$freq, fill = .data$aa)) +
    ggplot2::geom_histogram(binwidth = bin_width, color = "black") +
    ggplot2::geom_vline(data = annot_data, ggplot2::aes(xintercept = .data$freq, color = .data$group), linetype = "dashed") +
    ggplot2::geom_text(data = annot_data_labels, ggplot2::aes(x = -Inf, y = Inf, label = .data$combined_label), hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(permut_data$codon), scales = "fixed")

  n_aa <- length(unique(permut_data$aa))
  if (is.null(color_palette) && exists("aa_colors")) color_palette <- aa_colors
  plot <- plot + GetSafeColorScale(n_colors = n_aa, color_palette = color_palette, aes_type = "fill")

  plot <- plot + ggplot2::theme_bw() + ggplot2::theme(legend.position = show_legend) # Including the legend
  if (isTRUE(add_titles)) plot <- plot + ggplot2::labs(title = "Histogram of Codon Frequencies in Genes of Interest", x = "Frequency", y = "Count") # Adding standard titles to the plot
  if (!is.null(save_format)) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory) # Save the ggplot

  return(plot)
}
