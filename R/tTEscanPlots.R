# All the functions in this file are called by the user

#' Distribution Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply
#' @description
#' This function generates a visualization of codon-anticodon usage or amino acid demand-supply distributions across conditions.
#' The input data (\code{data}) is expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' Based on the selected \code{plot} parameter the user can select the most convenient layout to display the data.
#' It is crucial to ensure consistency between the name of the columns in \code{data} and the parameters \code{x_axis_col}, \code{y_axis_col} and \code{condition_col}.
#'
#' @param data A long format table. This format can be obtained using \code{\link{DataToLongFormat}}.
#' @param plot Either \code{"jitter"} (default), \code{"barplot"}, \code{"boxplot"} or \code{"line"} to indicate the type of plot to generate.
#' @param bar_position Either \code{"dodge"} (default), \code{"stack"} or \code{"fill"} to indicate the display of the plot if \code{plot} is \code{line}.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param condition_col Name of the categorical variable to group the data points.
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param targeted_arg Optional; a vector defining key feature clusters to highlight or label.
#' @param ncols Numeric; number of columns for arranging panels. Defaults to 1.
#' @param facet_col Optional; name of the categorical variable to divide the plot into different panels. Required if \code{ncols} bigger than 1.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested distribution plot. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_DistributionPlot <- function(data, plot = "jitter", bar_position = "dodge", x_axis_col, y_axis_col, condition_col, ncols = 1, facet_col = NULL,
                                 color_palette = NULL, targeted_arg = NULL, save_format = NULL, out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE) {

  ###
  # DESCRIPTION: This function takes a long format data expressing the usage of different features across a set of conditions.
  ###

  # Checking the input parameters
  CheckDataInLongFormat(data) # Data in long format.
  CheckValueInData(param = "plot", observed = plot, expected = c("jitter", "boxplot", "barplot", "line")) # Valid plot layout
  CheckValueInData(param = "_col", observed = c(x_axis_col, y_axis_col, condition_col), expected = colnames(data)) # All columns present in data

  if (!is.null(facet_col)){
    CheckValueInData(param = "_col", observed = c(facet_col), expected = colnames(data))
    if (is.numeric(data[[facet_col]])) stop(paste("Error in `data`: Column", facet_col, "is not categorical."))
  }

  if (!show_legend %in% c("none", "top", "bottom", "right", "left")) stop("Please specify a suitable `show_legend` parameter.\nSupported formats: none, top, bottom, right or left.")
  if (!(is.numeric(data[[y_axis_col]]))) stop(paste("Error in `data`: Column", y_axis_col, "is not numeric."))
  if (is.numeric(data[[x_axis_col]])) stop(paste("Error in `data`: Column", x_axis_col, "is not categorical."))
  if (is.numeric(data[[condition_col]])) stop(paste("Error in `data`: Column", condition_col, "is not categorical."))
  if (!bar_position %in% c("stack", "dodge", "fill")) stop("Please specify a suitable `bar_position` parameter.\nSupported formats: stack, dodge or fill.")

  target <- condition_col

  # TARGETED APPROACH - compare specific conditions against the rest
  if (!(is.null(targeted_arg))){
    data$target <- "other" # Initialize all the values of the target column as other

    for (i in 1:length(targeted_arg)){ # Iterate over the values in targeted_arg
      # Check which entry match the value in targeted_arg
      index <- which(data[[condition_col]] == targeted_arg[i])
      if (length(index) == 0) stop(paste("The argument", targeted_arg[i], "has not been found in the categories defined in the data."))
      data$target[index] <- targeted_arg[i] # Update the value of the target column
    }
    CheckDataFrame(data = data) # Ensure that the update data is suitable for the downstream analysis
  }

  # Adjust the colors of the plot
  if (!(is.null(color_palette)) && (length(color_palette) != (length(unique(data[[target]]))))){
    stop("The number of colors in `color_palette` does not correspond to the number of categories needed.\n",
         "Note: if the `targeted_arg` is given, an extra color is required.\n",
         paste("Number of colors given:", length(color_palette), "\n"), paste("Number of colors needed:", length(unique(data[[target]]))))
  }

  # Customization
  plot <- GenerateDistPlot(level = plot, target = target, data = data, x_axis_col = x_axis_col, y_axis_col = y_axis_col, color_palette = color_palette,
                           add_titles = add_titles, show_legend = show_legend, ncols = ncols, facet_col = facet_col, bar_position = bar_position)
  if (!(is.null(save_format))) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory)

  return(plot) # Returns the generated plot and exports a file if enabled.
}

#' Distribution Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply Compared to the Mean
#' @description
#' This function generates a distribution plot to compare distribution usages between a targeted set of conditions against the rest.
#' A common application of this plot is to compare mean usages across conditions.
#' Both input data sources (\code{target_data} and \code{overall_data}) are expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' It is crucial to ensure consistency between the name of the columns in \code{target_data} and \code{overall_data}, together with the parameters \code{x_axis_col} and \code{y_axis_col}.
#'
#' @param target_data A long format table of a condition of interest. This format can be obtained using \code{\link{DataToLongFormat}}.
#' @param overall_data A long format table of the whole set of conditions (with or without the data in \code{target_data}. This format can be obtained using \code{\link{DataToLongFormat}}.
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

tTE_CompareTargetToMean <- function(target_data, overall_data, x_axis_col, y_axis_col, show_difference = TRUE, color_palette = NULL,
                                    save_format = NULL, out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE){

  ###
  # DESCRIPTION: This function takes two long format data and expressing the usage of different features across a set of conditions and compares them.
  ###

  # Check that the names given correspond to column names.
  CheckValueInData(param = "x_axis_col", observed = x_axis_col, expected = colnames(overall_data)) # features
  CheckValueInData(param = "y_axis_col", observed = y_axis_col, expected = colnames(overall_data)) # values
  CheckValueInData(param = "x_axis_col", observed = x_axis_col, expected = colnames(target_data))
  CheckValueInData(param = "y_axis_col", observed = y_axis_col, expected = colnames(target_data))
  if (!(is.numeric(overall_data[[y_axis_col]]))) stop(paste("Error in `overall_data`: Column", y_axis_col, "is not numeric."))
  if (!(is.numeric(target_data[[y_axis_col]]))) stop(paste("Error in `target_data`: Column", y_axis_col, "is not numeric."))
  if (is.numeric(overall_data[[x_axis_col]])) stop(paste("Error in `overall_data`: Column", x_axis_col, "is not categorical."))
  if (is.numeric(target_data[[x_axis_col]])) stop(paste("Error in `target_data`: Column", x_axis_col, "is not categorical."))

  # Extract the data from the target condition - Way to combine both data inputs into a single to be called when generating the plot
  overall_data$target <- target_data[[y_axis_col]]

  # Adjust the colors of the plot
  if (!is.null(color_palette)){
    if (isTRUE(is.list(color_palette))) stop("`color_palette` needs to be a vector.")
    if (length(color_palette) < 2) stop("Wrong number of colors in `color_palette`.\n", "Required colors: one for the data points, one for the mean bars and (if applicable) an extra color to show the difference.\n")

    if (isTRUE(show_difference) && (length(color_palette) < 3)){
      message("- The difference can not be plotted as `color_palette` contains just two colors.")
      show_difference <- FALSE
    }

  } else {
    color_palette <- c("#CCCCCC", "#B10DA1", "blue") # Use default color_palette
  }

  # Generate plot
  plot <- ggplot2::ggplot(data = overall_data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]])) +
    ggplot2::geom_bar(stat = "identity", fill = color_palette[1]) +
    ggplot2::geom_point(mapping = ggplot2::aes(y = .data$target), color = color_palette[2]) +
    ggplot2::theme_bw()

  # Customization:
  if (isTRUE(show_difference)) plot <- plot + ggplot2::geom_segment(mapping = ggplot2::aes(xend = .data[[x_axis_col]], yend = .data$target), color = color_palette[3]) # Include difference bars between the targeted point and the rest of the data
  plot <- plot + ggplot2::theme(ggplot2::labs(x = x_axis_col, y = y_axis_col), axis.text.x = ggplot2::element_text(angle = 90), legend.position = show_legend) # Arrange elements in the plot
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
#' @param data A long format table. This format can be obtained using \code{\link{DataToLongFormat}}.
#' @param plot Either \code{"bar"} (default), \code{"donut"} or \code{"radar"} to indicate the type of plot to generate.
#' @param var_numerical Name of the numerical variable to reflect in the plot.
#' @param var_categorical Name of the categorical variable to reflect in the plot.
#' @param var_color Optional; name of the categorical variable to group the data points when coloring. Required if \code{plot} is \code{"bar"}.
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested proportion plot. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_ProportionPlot <- function(data, plot = "bar", var_numerical, var_categorical, var_color = NULL, color_palette = NULL,
                               save_format = NULL,  out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE){

  ###
  # DESCRIPTION: This function takes long format table to compare proportions of each condition with respect to the total values of the data.
  ###

  # Checking the type of plot and the input variables
  CheckValueInData(param = "plot", observed = plot, expected = c("donut", "radar", "bar"))
  if (is.null(var_color)) var_color <- var_categorical
  CheckValueInData(param = "var_", observed = c(var_categorical, var_numerical, var_color), expected = colnames(data))
  if (!(is.numeric(data[[var_numerical]]))) stop(paste("Error in `data`: Column", var_numerical, "is not numeric."))
  if (is.numeric(data[[var_categorical]])) stop(paste("Error in `data`: Column", var_categorical, "is not categorical."))
  if (is.numeric(data[[var_color]])) stop(paste("Error in `data`: Column", var_color, "is not categorical."))

  # Generate plot
  plot <- GenerateProportionPlot(data = data, level = plot, var_numerical = var_numerical, var_categorical = var_categorical, var_color = var_color, color_palette = color_palette, show_legend = show_legend)

  # Save the ggplot
  if(!is.null(save_format)) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory)
  return(plot) # Returns the generated plot and exports a file if enabled.
}

#' Violin Plot Displaying tTE Scores
#' @description
#' This function generates a violin plot to visualize the distribution of tTE scores across different condition. The tTE scores should be obtained using \code{\link{Compute_tTE}}.
#' The plot allows the user to easily compare the tTE distribution between different conditions, with the option to highlight specific feature clusters based on the provided \code{targets}.
#'
#' @param data A tTE results table obtained from \code{\link{Compute_tTE}}.
#' @param conditions A character vector specifying condition labels for samples.
#' @param name_sep A string delimiter used in \code{conditions} to separate the sample group labels.
#' @param targets A table defining key feature clusters to highlight. The first column should contain the conditions to select, and the second column the labels to use in the comparisons.
#' @param color_palette Optional; a vector of color codes to customize plot appearance.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format} specified).
#' @param out_directory Optional; path to the directory where the plot will be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"}, \code{"bottom"}, \code{"right"} or \code{"left"} to specify the position of the legend in the plot.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot. Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the tTE scores. If \code{save_format} is provided, the plot will also be saved to the specified location.
#' @export

tTE_ScoresPlot <- function(data, conditions, name_sep, targets, color_palette = NULL, save_format = NULL,
                           out_name = NULL, out_directory = NULL, show_legend = "none", add_titles = TRUE){

  ###
  # DESCRIPTION: This function takes the output table of Compute_tTE(), which contains the tTE scores nd the significant scores.
  # It performs a targeted approach displaying the tTE scores of a selected conditions against all the others.
  ###

  # Check that all the required columns are in the input data
  CheckValueInData(param = "data", observed = colnames(data), expected = c("condition", "tTE", "p_value", "neg_log10_tTE_p_value"))
  if (is.numeric(data$condition)) stop(paste("Error in `data`: Column condition is not categorical."))
  if (!(is.numeric(data$tTE))) stop(paste("Error in `data`: Column tTE is not categorical."))
  if (!(is.numeric(data$p_value))) stop(paste("Error in `data`: Column p_value is not categorical."))
  if (!(is.numeric(data$neg_log10_tTE_p_value))) stop(paste("Error in `data`: Column neg_log10_tTE_p_value is not categorical."))

  data <- tidyr::separate(data, .data$condition, into = conditions, sep = name_sep)
  data$class <- "other" # Initialize the class column - all values set to "other"

  # Iterates over the number of targets (groups of interest to compare)
  for (i in 1:nrow(targets)) data$class[grep(targets[i,1], data[[conditions[2]]])] <- targets[i,2]

  # Evaluates the color_palette if specified
  if(!(is.null(color_palette)) && (length(color_palette) != length(unique(data$class)))){
    stop("Wrong number of colors in `color_palette`.\n",
         "Required colors: one per condition in `targets` plus an extra for the non-targeted groups.\n",
         paste("Number of colors needed:", length(unique(data$class))))
  }

  # Defines the plot background
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data$class, y = .data$tTE, fill = .data$class)) +
    ggplot2::theme_bw() + ggplot2::geom_violin() + ggplot2::geom_jitter(shape = 16, position = ggplot2::position_jitter(0.2), size = 3)

  # Customization
  if (!(is.null(color_palette))) plot <- plot + ggplot2::scale_fill_manual(values = color_palette) # Incorporates the color_palette (if any)
  plot <- plot + ggplot2::theme(legend.position = show_legend) # Incorporates the legend
  if (isTRUE(add_titles)) plot + ggplot2::labs(title = "Theoretical translation efficiencies (tTE)")
  plot + ggplot2::labs(x = "Condition classification", y = "tTE") # Adding standard titles to the plot
  if (!is.null(save_format)) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory) # Save the ggplot

  return(plot) # Returns the generated violin plot and exports a file if enabled.
}

#' Correlation Plot: Exonic Codon Background - Mean Codon Usage
#' @description
#' This function generates a visualization to correlate different parmaeters of the \code{data}.
#' The input data (\code{data}) is expected to be in long format and contain a minimum information regarding the features, the conditions they belong to and their usage counts.
#' A common usage of this kind of plot is to represent the codon frequencies correlations of the exonic background and the mean codon usage across conditions.
#' Based on the selected \code{plot} parameter the user can select the most convenient layout to display the data.
#' It is crucial to ensure consistency between the name of the columns in \code{data} and the parameters \code{x_axis_col}, \code{y_axis_col} and \code{condition_col}.
#'
#' @param data A long format table. This format can be obtained using \code{\link{DataToLongFormat}}.
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
  if (!(plot_type %in% c("MeanCodonUsage", "PoolDiversity"))) stop("Please use MeanCodonUsage or PoolDiversity as `plot`.") # Suitable plot
  CheckDataInLongFormat(data) # Data in long format.
  if (is.null(label_col)){
    CheckValueInData(param = "_col", expected = colnames(data), observed = c(x_axis_col, y_axis_col, condition_col)) # Columns present in data
  } else {
    CheckValueInData(param = "_col", expected = colnames(data), observed = c(x_axis_col, y_axis_col, condition_col, label_col)) # Columns present in data
    if (is.numeric(data[[label_col]])) stop(paste("Error in `data`: Column", label_col, "is not categorical."))
  }

  if (!(is.numeric(data[[y_axis_col]]))) stop(paste("Error in `data`: Column", y_axis_col, "is not numeric."))
  if (!(is.numeric(data[[x_axis_col]]))) stop(paste("Error in `data`: Column", x_axis_col, "is not numeric"))
  if (is.numeric(data[[condition_col]])) stop(paste("Error in `data`: Column", condition_col, "is not categorical."))

  target <- condition_col

  # TARGETED APPROACH
  if (!is.null(targeted_arg)){
    data$target <- "other" # Initialize all the values of the target column as other

    for (i in 1:length(targeted_arg)){ # Iterate over the values in targeted_arg
      # Check which entry match the value in targeted_arg
      index <- which(data[[condition_col]] == targeted_arg[i])
      if (length(index) == 0) stop(paste("The argument", targeted_arg[i], "has not been found in the categories defined in the data."))

      data$target[index] <- targeted_arg[i] # Update the value of the target column
    }
    CheckDataFrame(data = data) # Ensure that the update data is suitable for the downstream analysis
  }

  # Adjust the colors of the plot
  if(!(is.null(color_palette)) && (length(color_palette) != (length(unique(data[[target]]))))){
    stop("The number of colors in `color_palette` does not correspond to the number of categories needed.\n",
         "Note: if the `targeted_arg` is given, an extra color is required.\n",
         paste("Number of colors given:", length(color_palette), "\n"), paste("Number of colors needed:", length(unique(data[[target]]))))
  }
  if (is.null(color_palette) && (length(unique(data[[target]])) < 36)) color_palette <- gradual_groups_35

  # Generate the plot
  if (is.null(label_col)){
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], color = .data[[target]]))
  } else {
    plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = .data[[x_axis_col]], y = .data[[y_axis_col]], label = .data[[label_col]], color = .data[[target]]))
  }

  # Customization:
  if (!(is.null(color_palette))) plot <- plot + ggplot2::scale_color_manual(values = color_palette) # Considering the user's color_palette (if any)
  plot <- plot + ggplot2::theme_bw() + ggplot2::geom_point(size = 3)
  if (plot_type == "MeanCodonUsage") plot <- plot + ggplot2::geom_abline(slope = 1, intercept = 0)
  if (!is.null(label_col)) plot <- plot + ggrepel::geom_text_repel(max.overlaps = 10)

  if (isTRUE(add_titles)){
    if(plot_type == "MeanCodonUsage") plot <- plot + ggplot2::labs(title = "Codon frequencies of exonic background \n vs. mean codon usage across conditions", x = "Mean codon usage across conditions", y = "Codon usage of exonic background")
    if(plot_type == "PoolDiversity") plot <- plot + ggplot2::labs(title = "Conditions' correlation to mean codon usage vs. \n codon pool diversity", x = "Codon pool diversity", y = "Correlation")
  }

  if (!(is.null(extra_val))){
    if (plot_type == "MeanCodonUsage") coord_values <- list(x = 0.047, y = 0.001)
    if (plot_type == "PoolDiversity") coord_values <- list(x = 0.95, y = 0.85)
    if (!is.null(label_col)) plot <- plot + ggplot2::annotate("text", x = coord_values$x, y = coord_values$y, label = sprintf("rho == %f", extra_val), parse = TRUE)
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

  # Check the feasibility of the color_palette
  if (!(is.null(color_palette)) && (length(color_palette) != length(unique(annot_data$aa)))) stop("The number of colors in `color_palette` does not correspond to the number of amino acids.\n",
                                                                                                 paste("Number of colors needed:", length(unique(annot_data$aa))))
  # Defining the size of the bars in the histogram based on the data
  bin_width <- (max(annot_data$freq) - min(annot_data$freq)) / 25
  if (is.null(bin_width)) stop("The bind width could not be calculated.\n", "Revise the `freq` column in the data.")

  # Defining the plot background
  plot <- ggplot2::ggplot(permut_data, ggplot2::aes(x = .data$freq, fill = .data$aa)) +
    ggplot2::geom_histogram(binwidth = bin_width, color = "black") +
    ggplot2::geom_vline(data = annot_data, ggplot2::aes(xintercept = .data$freq, color = .data$group), linetype = "dashed") +
    ggplot2::geom_text(data = annot_data_labels, ggplot2::aes(x = -Inf, y = Inf, label = .data$combined_label), hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ codon, scales = "fixed")

  # Customization:
  if (!(is.null(color_palette))) plot <- plot + ggplot2::scale_fill_manual(values = aa_colors) # Considering the user's color_palette (if any)
  plot <- plot + ggplot2::theme_bw() + ggplot2::theme(legend.position = show_legend) # Including the legend
  if (isTRUE(add_titles)) plot <- plot + ggplot2::labs(title = "Histogram of Codon Frequencies in Genes of Interest", x = "Frequency", y = "Count") # Adding standard titles to the plot
  if (!is.null(save_format)) SavePlot(plot = plot, save_format = save_format, out_name = out_name, out_directory = out_directory) # Save the ggplot

  return(plot)
}
