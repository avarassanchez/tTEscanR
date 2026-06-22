#' Distribution Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply
#' @description
#' This function generates a visualization of codon-anticodon usage or amino
#' acid demand-supply distributions across conditions. The input data
#' (\code{data}) is expected to be in long format and contain a minimum
#' information regarding the features, the conditions they belong to and their
#' usage counts. Based on the selected \code{plot} parameter the user can
#' select the most convenient layout to display the data. It is crucial to
#' ensure consistency between the name of the columns in \code{data} and the
#' parameters \code{x_axis_col}, \code{y_axis_col} and \code{condition_col}.
#'
#' @param data A long format table. This format can be obtained using
#'     \code{\link{transformFormat}}.
#' @param plot Either \code{"jitter"} (default), \code{"barplot"},
#'     \code{"boxplot"} or \code{"dot"} to indicate the type of plot to
#'     generate.
#' @param bar_position Either \code{"dodge"} (default), \code{"stack"} or
#'     \code{"fill"} to indicate the display of the plot if \code{plot} is
#'     \code{"barplot"}.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param condition_col Name of the categorical variable to group the data
#'     points.
#' @param color_palette Optional; a vector of color codes to customize plot
#'     appearance.
#' @param targeted_arg Optional; a vector defining key feature clusters to
#'     highlight or label.
#' @param ncols Numeric; number of columns for arranging panels. Defaults to 1.
#' @param facet_col Optional; name of the categorical variable to divide the
#'     plot into different panels. Required if \code{ncols} bigger than 1.
#' @param add_stats Logical; if \code{TRUE}, performs a statistical analysis
#'     based on the available parameters. Defaults to \code{FALSE}.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify
#'     the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format}
#'     specified).
#' @param out_directory Optional; path to the directory where the plot will
#'     be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot.
#'     Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested distribution plot.
#'     If \code{save_format} is provided, the plot will also be saved to the
#'     specified location.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#'
#' # Define the object and compute the codon usage
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data,
#'     assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
#'
#' # Transform the data
#' long_cu <- transformFormat(
#'     data = getAssay(tTEscanR_obj, "CodonUsage"), normalize = TRUE,
#'     rownames_to_column = "codon", names_to = "condition", values_to = "usage"
#' )
#' long_cu <- long_cu |>
#'     tidyr::separate(condition, into = c("tissue", "cell_type"), sep = "-")
#'
#' # Generate the plot
#' codon_usage_plot <- plotDistribution(
#'     data = long_cu, plot = "jitter", x_axis_col = "codon",
#'     y_axis_col = "usage", condition_col = "tissue", show_legend = "right",
#'     add_titles = FALSE
#' )
plotDistribution <- function(data, plot = "jitter", bar_position = "dodge",
    x_axis_col, y_axis_col, condition_col, color_palette = NULL, ncols = 1,
    facet_col = NULL, add_stats = FALSE, targeted_arg = NULL,
    save_format = NULL, out_name = NULL, out_directory = NULL,
    show_legend = "none", add_titles = TRUE, verbose = TRUE) {
    generalChecksDistPlot( # Checking the input parameters
        data = data, x_axis_col = x_axis_col, y_axis_col = y_axis_col,
        condition_col = condition_col, facet_col = facet_col,
        bar_position = bar_position, plot = plot, show_legend = show_legend
    )
    if (!(is.null(targeted_arg))) { # TARGETED - specific condition vs the rest
        data$target_group <- "other" # Initialize values of the target column
        for (arg in targeted_arg) { # Iterate and check match with targeted_arg
            index <- which(data[[condition_col]] == arg)
            if (length(index) == 0) {
                stop("The argument ", arg, " has not been found in the data.")
            }
            data$target_group[index] <- arg # Update value of the target column
        }
        condition_col <- "target_group"
        checkDataFrame(data = data) # Ensure the data is suitable
    }
    if (is.list(color_palette)) color_palette <- unlist(color_palette)
    num_target <- length(unique(data[[condition_col]]))
    if (!(is.null(color_palette)) && (length(color_palette) < num_target)) {
        stop(
            "The number of colors in 'color_palette' does not correspond to ",
            "the number of categories needed.\n", "Note: if 'targeted_arg' is ",
            "given, an extra color is required.\n", "Number of colors given: ",
            length(color_palette), "\n", "Number of colors needed: ", num_target
        )
    }
    plot <- generateDistPlot( # Customization
        level = plot, add_titles = add_titles, facet = facet_col, data = data,
        add_stats = add_stats, x_axis = x_axis_col, y_axis = y_axis_col,
        color = color_palette, target = condition_col, ncols = ncols,
        show_legend = show_legend, bar = bar_position
    )
    if (!(is.null(save_format))) {
        savePlot(
            plot = plot$plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory, verbose = verbose
        )
    }
    return(plot$plot) # Returns the generated plot & exports a file if enabled
}

generalChecksDistPlot <- function(data, plot, x_axis_col, y_axis_col,
    condition_col, facet_col, bar_position, show_legend) {
    checkDataInLongFormat(data) # Data in long format.

    # Valid plot layout
    checkValueInData(
        param = "plot", observed = plot,
        expected = c("jitter", "boxplot", "barplot", "dot")
    )

    checkValueInData(
        param = "_col", observed = c(x_axis_col, y_axis_col, condition_col),
        expected = colnames(data)
    ) # All columns present in data

    if (!is.null(facet_col)) {
        checkValueInData(
            param = "_col", observed = c(facet_col), expected = colnames(data)
        )
        if (is.numeric(data[[facet_col]])) {
            stop("The 'data' column ", facet_col, " is not categorical.")
        }
    }

    if (!show_legend %in% c("none", "top", "bottom", "right", "left")) {
        stop(
            "Specify a suitable 'show_legend' parameter.\n",
            "Supported formats: none, top, bottom, right or left."
        )
    }
    if (!(is.numeric(data[[y_axis_col]]))) {
        stop("The 'data' column ", y_axis_col, " is not numeric.")
    }
    if (is.numeric(data[[x_axis_col]])) {
        stop("The 'data' column ", x_axis_col, " is not categorical.")
    }
    if (is.numeric(data[[condition_col]])) {
        stop("The 'data' column ", condition_col, " is not categorical.")
    }
    if (!bar_position %in% c("stack", "dodge", "fill")) {
        stop(
            "Specify a suitable 'bar_position' parameter.\n",
            "Supported formats: stack, dodge or fill."
        )
    }

    return(invisible(NULL))
}

#' Distribution Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply
#' Compared to a Target
#' @description
#' This function generates a distribution plot to compare distribution usages
#' between a targeted set of conditions against the rest. A common application
#' of this plot is to compare mean usages across conditions. Both input data
#' sources (\code{target_data} and \code{overall_data}) are expected to be in
#' long format and contain a minimum information regarding the features, the
#' conditions they belong to and their usage counts. It is crucial to ensure
#' consistency between the name of the columns in \code{target_data} and
#' \code{overall_data}, together with the parameters \code{x_axis_col} and
#' \code{y_axis_col}.
#'
#' @param target_data A long format table of a condition of interest. This
#'     format can be obtained using \code{\link{transformFormat}}.
#' @param overall_data A long format table of the whole set of conditions
#'     (with or without the data in \code{target_data}. This format can be
#'     obtained using \code{\link{transformFormat}}.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param show_difference Logical; if \code{TRUE}, displays the difference
#'     between \code{overall_data} and \code{target_data}.
#' @param color_palette Optional; a vector of color codes (min. 2 codes and max.
#'     3 codes) to customize plot appearance. Colors for (i) data, (ii) target
#'     value, and (iii) difference bar.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify
#'     the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format}
#'     specified).
#' @param out_directory Optional; path to the directory where the plot will
#'     be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot.
#'     Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested distribution plot.
#'     If \code{save_format} is provided, the plot will also be saved to the
#'     specified location.
#' @export
plotTargetComparison <- function(target_data, overall_data, x_axis_col,
    y_axis_col, color_palette = NULL, show_difference = TRUE,
    save_format = NULL, out_name = NULL, add_titles = TRUE,
    out_directory = NULL, show_legend = "none", verbose = TRUE) {
    generalChecksPlotTargetComp( # Check names correspond to column names
        x_axis = x_axis_col, y_axis = y_axis_col,
        target = target_data, data = overall_data
    )
    target_subset <- target_data[, c(x_axis_col, y_axis_col), drop = FALSE]
    colnames(target_subset)[2] <- "target_value"
    plot_data <- merge(
        overall_data, target_subset,
        by = x_axis_col, all.x = TRUE
    )

    color_palette <- selectColorPalette(
        palette = color_palette, diff = show_difference
    )
    plot <- ggplot2::ggplot(data = plot_data, mapping = ggplot2::aes(
        x = .data[[x_axis_col]], y = .data[[y_axis_col]]
    )) +
        ggplot2::geom_col(fill = color_palette[1]) +
        ggplot2::theme_bw()

    if (isTRUE(show_difference)) { # Add diff. bars: targeted point-data
        plot <- plot + ggplot2::geom_segment(
            mapping = ggplot2::aes(
                xend = .data[[x_axis_col]], yend = .data$target_value
            ),
            color = color_palette[3]
        )
    }
    plot <- plot + ggplot2::geom_point(mapping = ggplot2::aes(
        y = .data$target_value
    ), color = color_palette[2], size = 2) + ggplot2::labs(
        x = x_axis_col, y = y_axis_col) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
            legend.position = show_legend
        ) # Arrange elements in the plot
    if (!(is.null(save_format))) {
        savePlot(
            plot = plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory, verbose = verbose
        )
    } # Save the ggplot
    return(plot) # Returns the generated plot and exports a file if enabled.
}

selectColorPalette <- function(palette, diff) {
    # The color palette is a list - vector like required
    if (is.list(palette)) palette <- unlist(palette)

    # Use default color_palette
    if (is.null(palette)) palette <- c("#bababa", "#b2182b", "#fddbc7")

    # Validate lenght of the oclor palette
    colors_needed <- ifelse(isTRUE(diff), 3, 2)

    if (length(palette) < colors_needed) {
        stop(
            "The number of colors in 'color_palette' does not correspond to ",
            "the number of graphical elements.\n", "Note: If 'show_difference'",
            " is TRUE, an extra color is required for the segments.\n",
            "Number of colors given: ", length(palette), "\n",
            "Number of colors needed: ", colors_needed
        )
    }

    return(palette)
}

generalChecksPlotTargetComp <- function(x_axis, y_axis, target, data) {
    checkValueInData(
        param = "x_axis_col", observed = x_axis, expected = colnames(data)
    ) # features
    checkValueInData(
        param = "y_axis_col", observed = y_axis, expected = colnames(data)
    ) # values
    checkValueInData(
        param = "x_axis_col", observed = x_axis, expected = colnames(target)
    )
    checkValueInData(
        param = "y_axis_col", observed = y_axis, expected = colnames(target)
    )

    if (!(is.numeric(data[[y_axis]]))) {
        stop("The 'overall_data' column ", y_axis, " is not numeric.")
    }

    if (!(is.numeric(target[[y_axis]]))) {
        stop("The 'target_data' column ", y_axis, " is not numeric.")
    }

    if (is.numeric(data[[x_axis]])) {
        stop("The 'overall_data' column ", x_axis, " is not categorical.")
    }

    if (is.numeric(target[[x_axis]])) {
        stop("The 'target_data' column ", x_axis, " is not categorical.")
    }

    return(invisible(NULL))
}

#' Proportion Plot of Codon/Anticodon Usage or Amino Acid Demand/Supply
#' @description
#' This function generates a proportion plot to compare codon-anticodon usage
#' or amino acid demand-supply frequencies across conditions. The input data
#' (\code{data}) is expected to be in long format and contain a minimum
#' information regarding the features, the conditions they belong to and their
#' usage counts. Based on the selected \code{plot} parameter the user can
#' select the most convenient layout to display the data. It is crucial to
#' ensure consistency between the name of the columns in \code{data} and the
#' parameters \code{var_numerical}, \code{var_categorical} and \code{var_color}.
#'
#' @param data A long format table. This format can be obtained using
#'     \code{\link{transformFormat}}.
#' @param plot Either \code{"bar"} (default), \code{"donut"} or \code{"radar"}
#'     to indicate the type of plot to generate.
#' @param var_numerical Name of the numerical variable to reflect in the plot.
#' @param var_categorical Name of the categorical variable to reflect in the
#'     plot.
#' @param var_color Optional; name of the categorical variable to group the
#'     data points when coloring. Required if \code{plot} is \code{"bar"}.
#' @param facet_col Optional; name of the categorical variable to divide the
#'     plot into different panels.
#' @param color_palette Optional; a vector of color codes to customize plot
#'     appearance.
#' @param num_limits Optional; a vector with the upper and lower ranges of the
#'     values in \code{var_numerical}.
#' @param num_rings Optional; a number specifying the amount of rings to
#'     display if \code{plot} is \code{"radar"}. Defaults to 5.
#' @param order Optional; a vector of the levels to organize the data, based
#'     on the \code{var_categorical}.
#' @param normalize Logical; if \code{TRUE}, normalizes the data in order to
#'     display the relative contribution of the \code{var_categorical}.
#'     Defaults to \code{TRUE}.
#' @param zoom Logical; if \code{TRUE}, centers the plot display to the values
#'     in \code{data}.Defaults to \code{TRUE}.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify
#'     the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format}
#'     specified).
#' @param out_directory Optional; path to the directory where the plot will be
#'     saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot.
#'     Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the requested proportion plot.
#'     If \code{save_format} is provided, the plot will also be saved to the
#'     specified location.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#'
#' # Define the object and compute the codon usage
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data, assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = TRUE, reduce = 1000
#' )
#'
#' # Compute and extract the mean codon usage
#' additional_metrics <- getMetadata(
#'     tTEscanR_obj, "CodonUsage_AdditionalMetrics"
#' )
#' mean_codon_usage <- additional_metrics$MeanCodonUsage
#' mean_codon_usage$codon <- mean_codon_usage$feature
#'
#' # Translate the codons to amino acids
#' mean_codon_usage <- featuresToAA(
#'     data = mean_codon_usage, position = "feature",
#'     notation_from = "codon", notation_to = "aa", verbose = FALSE
#' )
#'
#' # Generate the plot
#' plotProportion(
#'     data = mean_codon_usage, plot = "bar",
#'     var_numerical = "mean_usage_across_conditions",
#'     var_categorical = "codon", var_color = "feature", show_legend = "none"
#' )
plotProportion <- function(data, plot = "bar", var_numerical, var_categorical,
    var_color = NULL, facet_col = NULL, color_palette = NULL, num_limits = NULL,
    num_rings = 5, save_format = NULL, out_name = NULL,
    out_directory = NULL, zoom = FALSE, show_legend = "none",
    add_titles = TRUE, order = NULL, normalize = TRUE, verbose = TRUE) {
    # Checking the type of plot and the input variables
    if (is.null(var_color)) var_color <- var_categorical
    generalChecksProportionPlot(
        plot = plot, data = data, var_color = var_color, facet_col = facet_col,
        var_categorical = var_categorical, var_numerical = var_numerical
    )

    # Generate plot
    plot <- generateProportionPlot(
        data = data, level = plot, var_numerical = var_numerical,
        var_categorical = var_categorical, var_color = var_color,
        color_palette = color_palette, show_legend = show_legend, order = order,
        x_limits = num_limits, facet_col = facet_col, normalize = normalize,
        zoom = zoom, add_titles = add_titles, n_rings = num_rings
    )

    # Save the ggplot
    if (!is.null(save_format)) {
        savePlot(
            plot = plot$plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory, verbose = verbose
        )
    }

    return(plot) # Returns the generated plot and exports a file if enabled.
}

generalChecksProportionPlot <- function(plot, data, var_color, var_categorical,
    var_numerical, facet_col) {
    checkValueInData(
        param = "plot", observed = plot, expected = c("donut", "radar", "bar")
    )

    checkValueInData(
        param = "var_",
        observed = c(var_categorical, var_numerical, var_color),
        expected = colnames(data)
    )

    if (!is.null(facet_col)) {
        checkValueInData(
            param = "_col", observed = c(facet_col), expected = colnames(data)
        )
    }

    if (!(is.numeric(data[[var_numerical]]))) {
        stop("The 'data' column ", var_numerical, " is not numeric.")
    }

    if (is.numeric(data[[var_categorical]])) {
        stop("The 'data' column ", var_categorical, " is not categorical.")
    }

    if (is.numeric(data[[var_color]])) {
        stop("The 'data' column ", var_color, " is not categorical.")
    }

    if (!is.null(facet_col) && is.numeric(data[[facet_col]])) {
        stop("The 'data' column ", facet_col, " is not categorical.")
    }

    return(invisible(NULL))
}

#' Violin Plot Displaying tTE Scores
#' @description
#' This function generates a violin plot to visualize the distribution of tTE
#' scores across different condition. The tTE scores should be obtained using
#' \code{\link{computeTheoreticalTE}}. The plot allows the user to easily
#' compare the tTE distribution between different conditions, with the option
#' to highlight specific feature clusters based on the provided \code{targets}.
#'
#' @param data A tTE results table obtained from
#'      \code{\link{computeTheoreticalTE}}.
#' @param metadata A table with additional information regarding the conditions
#'     in \code{data}.
#' @param score_col Name of the numerical variable that contains the tTE scores
#'      in the \code{data}. Defaults to \code{"tTE"}.
#' @param cond_col Name of the categorical variable that contains the conditions
#'      in the \code{data}. Defaults to \code{"condition"}.
#' @param pval_col Name of the numerical variable that contains the significance
#'      scores in \code{data}. Defaults to \code{"p_value"}.
#' @param index_col Name of the categorical variable that links the conditions
#'     in \code{data} with the \code{metadata}.
#' @param class_col Name of the categorical variable to reflect in the plot.
#' @param add_stats Logical; if \code{TRUE}, performs a statistical analysis
#'     based on the available parameters. Defaults to \code{TRUE}.
#' @param target_col Optional; name of the categorical variable to perform the
#'     statistical comparison (the most specific level). Used if
#'     \code{add_stats} is \code{TRUE}.
#' @param facet_col Optional; name of the categorical variable separate the
#'     plot into different panels.
#' @param color_palette Optional; a vector of color codes to customize plot
#'     appearance.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify
#'     the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format}
#'     specified).
#' @param out_directory Optional; path to the directory where the plot will be
#'     saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend in the plot.
#' @param show_outliers Logical; if \code{TRUE}, labels the outlier data points
#'     based on \code{index_col}.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot.
#'     Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the tTE scores. If
#'     \code{save_format} is provided, the plot will also be saved to the
#'     specified location. If \code{add_stats} reports a table with the
#'     statitical measures summarized.
#' @export
#'
#' @examples
#' data(
#'     default_tTEscanR_mRNA_data, default_tTEscanR_tRNA_data,
#'     default_tTEscanR_metadata
#' )
#'
#' # Define the tTEscanR object
#' tTEscanR_obj <- createObject(
#'     counts = list(
#'         mRNA = default_tTEscanR_mRNA_data, tRNA = default_tTEscanR_tRNA_data
#'     ),
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#'
#' # Compute the codon and anticodon usage
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 10000
#' )
#' tTEscanR_obj <- computeAnticodonUsage(object = tTEscanR_obj)
#'
#' # Compute the theoretical translation efficiency (tTE scores)
#' tTEscanR_obj <- computeTheoreticalTE(
#'     object = tTEscanR_obj, level = "codon", compute_significance = TRUE
#' )
#' tTEresults_codon <- getMetadata(tTEscanR_obj, "tTEresults_codon")
#' conditions_metadata <- getMetadata(tTEscanR_obj, "ConditionsLabels")
#'
#' # Visualize the tTE scores
#' plotTEscore(
#'     data = tTEresults_codon, metadata = conditions_metadata,
#'     index_col = "conditions", class_col = "tissue", add_stats = TRUE
#' )
plotTEscore <- function(data, metadata, class_col, index_col, target_col = NULL,
    score_col = "tTE", cond_col = "condition", pval_col = "p_value",
    facet_col = NULL, color_palette = NULL, save_format = NULL,
    out_name = NULL, add_stats = TRUE, out_directory = NULL, verbose = TRUE,
    show_legend = "none", add_titles = TRUE, show_outliers = FALSE) {
    generalChecksScoresPlot(
        data = data, meta = metadata, stats = add_stats, facet = facet_col,
        class = class_col, target = target_col, index = index_col,
        score_col = score_col, cond_col = cond_col, pval_col = pval_col
    )
    merged <- defineMergedData(
        data = data, meta = metadata, index = index_col, facet = facet_col,
        class = class_col, target = target_col, outliers = show_outliers,
        cond = cond_col, score = score_col
    )
    fill_label <- if (!is.null(target_col)) target_col else "Group"
    plot <- ggplot2::ggplot(merged, ggplot2::aes(
        x = .data$class, y = .data[[score_col]], fill = .data$target
    )) + ggplot2::geom_violin(
            trim = FALSE, scale = "width", width = 0.7,
            position = ggplot2::position_dodge(width = 0.7), alpha = 0.7
        ) + ggplot2::geom_point(position = ggplot2::position_jitterdodge(
            dodge.width = 0.7, jitter.width = 0.03
        ), size = 0.8, alpha = 0.5) +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.4)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = class_col, y = "tTE", fill = fill_label)
    plot <- additionalCustomization(
        plot = plot, merged = merged, facet = facet_col, colors = color_palette,
        add_titles = add_titles, class = class_col, target = target_col,
        legend = show_legend, outliers = show_outliers, score = score_col
    )
    sig_table <- NULL
    if (isTRUE(add_stats)) {
        results_stat <- addStats(
            add_stats = add_stats, target = target_col, score = score_col,
            plot = plot, class = class_col, facet = facet_col, merged = merged
        )
        plot <- results_stat$plot
        sig_table <- results_stat$sig_table
    }
    if (!is.null(save_format)) {
        savePlot(
            plot = plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory, verbose = verbose
        )
    }
    return(list(plot = plot, stats = sig_table))
}

additionalCustomization <- function(plot, merged, outliers, facet, colors,
    add_titles, class, target, legend, score) {
    if (isTRUE(outliers)) { # Label outliers
        plot <- plot + ggrepel::geom_text_repel(
            data = dplyr::filter(merged, .data$is_extreme),
            ggplot2::aes(x=.data$class, y=.data[[score]], label=.data$label),
            position = ggplot2::position_jitterdodge(
                dodge.width = 0.8, jitter.width = 0.15
            ), size = 3, max.overlaps = 15, show.legend = FALSE
        )
    }
    if (length(unique(merged$class)) > 5) {
        plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90, hjust = 1, vjust = 0.5
        ))
    }
    if (!is.null(facet)) {
        plot <- plot + ggplot2::facet_wrap(
            ggplot2::vars(.data[[facet]]),
            scales = "free_x"
        )
    }
    if (!is.null(colors) && !is.null(names(colors))) {
        if (!any(names(colors) %in% unique(merged$target))) {
            warning(
                "Palette names do not match data categories. ",
                "Applying sequentially."
            )
            colors <- unname(colors) # Clean strip: forces sequential matching
        }
    }
    unique_categories <- unique(merged$target)
    checked_palette <- checkPaletteNames(
        color_palette = colors, actual_categories = unique_categories
    )
    plot <- plot + getSafeColorScale(
        n_colors = length(unique_categories), aes_type = "fill",
        color_palette = checked_palette
    )
    if (add_titles) {
        title <- paste("tTE scores by", class)
        if (!is.null(target)) {
            title <- paste(title, "highlighting", target)
        }
        if (!is.null(facet)) title <- paste(title, "faceted by", facet)
        plot <- plot + ggplot2::labs(title = title)
    }
    plot <- plot + ggplot2::theme(legend.position = legend)
    return(plot)
}

generalChecksScoresPlot <- function(data, meta, stats, target, facet,
    index, class, score_col, cond_col, pval_col) {
    # Extracted when computing the tTE through tTEscanR
    # required_cols <- c("condition", "tTE", "p_value", "neg_log10_tTE_p_value")
    required_cols <- c(score_col, cond_col, pval_col)

    checkValueInData("data columns", required_cols, colnames(data))
    checkValueInData("index_col", index, colnames(meta))
    checkValueInData("class_col", class, colnames(meta))

    if (isFALSE(stats) && !is.null(target)) {
        warning(
            "The 'target_col' will be ignored as 'add_stats' is not ",
            "enabled. \n Check the documentation of the function for more",
            " details."
        )
    }

    if (!is.null(target)) {
        checkValueInData("target_col", target, colnames(meta))
    }

    if (!is.null(facet)) {
        checkValueInData("facet_col", facet, colnames(meta))
    }

    return(invisible(NULL))
}

defineMergedData <- function(data, meta, index, class, target,
    outliers, facet, cond, score) {
    merged <- dplyr::left_join(
        data, meta,
        by = stats::setNames(index, cond)
    )
    if (nrow(merged) == 0) {
        stop(
            "The merged data ('data' and 'metadata') has 0 rows. ",
            "Check if 'index_col' matches 'cond_col'."
        )
    }
    merged$class <- merged[[class]]
    na_classes <- sum(is.na(merged$class))
    if (na_classes > 0) {
        warning(
            na_classes, " rows in data could not be matched to a group ",
            "in metadata!"
        )
    }
    group_counts <- table(merged$class)
    if (any(group_counts < 2)) {
        warning(
            "One or more groups have fewer than 2 samples. Statistics ",
            "(p-values) will return NA."
        )
    }
    merged$label <- merged[[cond]]
    if (!is.null(target)) {
        merged$target <- merged[[target]]
    } else {
        merged$target <- merged$class
    }
    if (isTRUE(outliers)) { # Label outliers
        group_vars <- c("class", "target")
        if (!is.null(facet)) group_vars <- c(group_vars, facet)
        merged <- merged %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
            dplyr::mutate(
                q1 = stats::quantile(.data[[score]], 0.25, na.rm = TRUE),
                q3 = stats::quantile(.data[[score]], 0.75, na.rm = TRUE),
                iqr = .data$q3 - .data$q1,
                is_extreme = .data[[score]] < (.data$q1 - 1.5 * .data$iqr) |
                    .data[[score]] > (.data$q3 + 1.5 * .data$iqr)
            ) %>%
            dplyr::ungroup()
    }
    return(merged)
}

addStats <- function(plot, add_stats, target, class, score,
    facet, merged) {
    if (!is.null(target)) {
        stat_target <- target
        if (stat_target != class) {
            group_col <- unique(c("class", facet))
            x_col <- "dummy_x"
        } else {
            group_col <- facet
            x_col <- "class"
        }
    } else {
        stat_target <- class
        group_col <- facet
        x_col <- "class"
    }
    var <- tapply(
        merged[[score]], merged[[stat_target]], stats::var, na.rm = TRUE
    )
    if (any(var == 0, na.rm = TRUE) || any(is.na(var))) {
        warning(
            "One or more groups have zero variance (all values are ",
            "identical). Stats will return NA."
        )
    }
    sig_table <- computeTEsignificance(
        merged,
        x_col = x_col, target = stat_target, group = group_col, score = score
    )
    if (!is.null(sig_table) && nrow(sig_table) > 0) {
        y_max <- max(merged[[score]], na.rm = TRUE) * 1.15 # stars above violins

        if (stat_target != class) {
            plot <- plot + ggpubr::stat_pvalue_manual(
                sig_table,
                x = "class", y.position = y_max, label = "p_signif",
                size = 7, hide.ns = TRUE, step.increase = 0.12
            )
        } else {
            plot <- plot + ggpubr::stat_pvalue_manual(
                sig_table,
                y.position = y_max, label = "p_signif", size = 7,
                hide.ns = TRUE, step.increase = 0.12
            )
        }
    } else {
        warning("The statistics could not be computed, check the comparisons.")
    }
    return(list(plot = plot, sig_table = sig_table))
}

#' Correlation Plot: Exonic Codon Background - Mean Codon Usage
#' @description
#' This function generates a visualization to correlate different parameters of
#' the \code{data}. The input data (\code{data}) is expected to be in long
#' format and contain a minimum information regarding the features, the
#' conditions they belong to and their usage counts. A common usage of this
#' kind of plot is to represent the codon frequencies correlations of the
#' exonic background and the mean codon usage across conditions. Based on the
#' selected \code{plot} parameter the user can select the most convenient layout
#' to display the data. It is crucial to ensure consistency between the name of
#' the columns in \code{data} and the parameters \code{x_axis_col},
#' \code{y_axis_col} and \code{condition_col}.
#'
#' @param data A long format table. This format can be obtained using
#'     \code{\link{transformFormat}}.
#' @param plot Either \code{"MeanCodonUsage"} (default) or
#'     \code{"PoolDiversity"} to indicate the \code{data} source.
#' @param x_axis_col Name of the categorical variable to reflect in the plot.
#' @param y_axis_col Name of the numerical variable to reflect in the plot.
#' @param condition_col Name of the categorical variable to group the data
#'     points.
#' @param label_col Name of the categorical variable to label the data points.
#' @param extra_val Optional; variable with additional information to include
#'     in the plot (e.g. correlation value).
#' @param color_palette Optional; a vector of color codes to customize plot
#'     appearance.
#' @param targeted_arg Optional; a vector defining key feature clusters to
#'     highlight or label.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify
#'     the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format}
#'     specified).
#' @param out_directory Optional; path to the directory where the plot will
#'     be saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot.
#'     Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the correlation. If
#'     \code{save_format} is provided, the plot will also be saved to the
#'     specified location.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#'
#' # Define the object and compute the codon usage
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data, assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = TRUE, reduce = 1000
#' )
#'
#' # Compute and extract the mean codon usage
#' additional_metrics <- getMetadata(
#'     tTEscanR_obj,
#'     "CodonUsage_AdditionalMetrics"
#' )
#' mean_codon_usage <- additional_metrics$MeanCodonUsage
#' exonic_background <- additional_metrics$CodonExonicBackground
#' exonic_background <- as.data.frame(exonic_background)
#' correlation_mean_background <- cbind(mean_codon_usage, exonic_background)
#'
#' plotCorrelation(
#'     data = correlation_mean_background, plot = "MeanCodonUsage",
#'     x_axis_col = "mean_usage_across_conditions",
#'     y_axis_col = "exonic_background", condition_col = "feature",
#'     extra_val = additional_metrics$MeanCodonCorr, add_titles = TRUE,
#'     show_legend = "none"
#' )
plotCorrelation <- function(data, plot = "MeanCodonUsage", x_axis_col,
    y_axis_col, condition_col, extra_val = NULL, label_col = NULL,
    color_palette = NULL, out_name = NULL, show_legend = "none",
    targeted_arg = NULL, save_format = NULL, out_directory = NULL,
    add_titles = TRUE, verbose = TRUE) {
    generalChecksCorrelationPlot( # Checking the input parameters
        plot = plot, x_axis_col = x_axis_col, y_axis_col = y_axis_col,
        label_col = label_col, data = data, condition_col = condition_col
    )
    if (!is.null(targeted_arg)) { # TARGETED APPROACH
        targeted_res <- targetedCorrelation(
            data = data, target = targeted_arg, cond = condition_col
        )
        data <- targeted_res$data
        condition_col <- targeted_res$condition_col
    }
    p <- ggplot2::ggplot(data, ggplot2::aes( # Generate the plot
        x = .data[[x_axis_col]], y = .data[[y_axis_col]],
        color = .data[[condition_col]]
    )) + ggplot2::geom_point(size = 3) + ggplot2::theme_bw()
    unique_categories <- unique(data[[condition_col]])
    checked_palette <- checkPaletteNames(
        color_palette = color_palette, actual_categories = unique_categories
    )
    p <- p + getSafeColorScale( # Adjust the colors of the plot
        n_colors = length(unique_categories), color_palette = checked_palette,
        aes_type = "color"
    )
    p <- addTitleCorrelation(p = p, plot = plot, add_titles = add_titles)
    if (!is.null(label_col)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(
            label = .data[[label_col]]
        ), max.overlaps = 10)
    }
    if (!(is.null(extra_val))) {
        p <- p + ggplot2::annotate(
            "text",
            x = Inf, y = -Inf, hjust = 1.2, vjust = -0.2,
            label = sprintf("rho == %.3f", extra_val), parse = TRUE
        )
    }
    p <- p + ggplot2::theme(legend.position = show_legend)
    if (!is.null(save_format)) {
        savePlot( # Save the ggplot
            plot = p, save_format = save_format, out_name = out_name,
            out_directory = out_directory, verbose = verbose
        )
    }
    return(p) # Returns the generated plot and exports a file if enabled
}

targetedCorrelation <- function(data, target, cond) {
    data$target <- "other" # Initialize all values of the target column
    missing <- target[!target %in% data[[cond]]]

    if (length(missing) > 0) {
        stop(
            "Arguments in 'targeted_arg' not in 'data': ",
            paste(missing, collapse = ", ")
        )
    }

    data$.__target__ <- ifelse(
        data[[cond]] %in% target, as.character(data[[cond]]), "other"
    )

    cond <- ".__target__" # Update pointer to use new column

    return(list(data = data, condition_col = cond))
}

addTitleCorrelation <- function(p, plot, add_titles) {
    if (plot == "MeanCodonUsage") {
        p <- p + ggplot2::geom_abline(slope = 1, intercept = 0)
        if (isTRUE(add_titles)) {
            p <- p + ggplot2::labs(
                title = "Codon frequencies of exonic background \n vs.
                mean codon usage across conditions",
                x = "Mean codon usage across conditions",
                y = "Codon usage of exonic background"
            )
        }
    } else if (plot == "PoolDiversity" && isTRUE(add_titles)) {
        p <- p + ggplot2::labs(
            title = "Correlation mean codon usage vs. \n codon pool diversity",
            x = "Codon pool diversity", y = "Correlation"
        )
    }

    return(p)
}

generalChecksCorrelationPlot <- function(plot, data, x_axis_col, y_axis_col,
    condition_col, label_col) {
    # Suitable plot
    if (!(plot %in% c("MeanCodonUsage", "PoolDiversity"))) {
        stop("Please use MeanCodonUsage or PoolDiversity as 'plot'.")
    }

    checkDataInLongFormat(data) # Data in long format.

    expected_cols <- c(x_axis_col, y_axis_col, condition_col, label_col)
    checkValueInData("data columns", expected_cols, colnames(data))

    if (!(is.numeric(data[[y_axis_col]]))) {
        stop("The 'data' column ", y_axis_col, " is not numeric.")
    }

    if (!(is.numeric(data[[x_axis_col]]))) {
        stop("The 'data' column ", x_axis_col, " is not numeric")
    }

    if (!is.character(data[[condition_col]]) &&
        !is.factor(data[[condition_col]])) {
        stop("The 'data' column ", condition_col, " is not categorical.")
    }

    if (!is.null(label_col) && is.numeric(data[[label_col]])) {
        stop("The 'data' column ", label_col, " is not categorical.")
    }

    return(invisible(NULL))
}

#' Permutation Plot
#' @description
#' This function generates a plot to compare the baseline codon exonic
#' background against the current codon usage. For a better interpretation, the
#' codons are colored by amino acid.
#'
#' @param permut_data A table with the codons and their frequencies after
#'     computing all the permutations. Output from
#'     \code{\link{getPermutationDist}}.
#' @param sig_data A table with the codon exonic background and their
#'     significance level before (p-value) and after the correction (p-adjusted
#'     value). Output from \code{\link{obtainSignificance}}
#' @param color_palette Optional; a vector of color codes to customize plot
#'     appearance.
#' @param save_format Optional; either \code{"png"} or \code{"pdf"} to specify
#'     the format to save the plot.
#' @param out_name Optional; name for the saved plot (if \code{save_format}
#'     specified).
#' @param out_directory Optional; path to the directory where the plot will be
#'     saved (if \code{save_format} specified).
#' @param show_legend Either \code{"none"} (default), \code{"top"},
#'     \code{"bottom"}, \code{"right"} or \code{"left"} to specify the
#'     position of the legend in the plot.
#' @param add_titles Logical; if \code{TRUE}, includes titles in the plot.
#'     Defaults to \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, displays information messages.
#'     Defaults to \code{TRUE}.
#'
#' @return Permutation plot.
#' @export
#'
#' @examples
#' data(default_tTEscanR_mRNA_data, default_tTEscanR_metadata)
#' selected_genes <- default_tTEscanR_mRNA_data[1:20, ]
#' permutation_test <- getPermutationDist(
#'     n_permut = 100, target_data = selected_genes, species = "hg38"
#' ) # Generate table with codon and freq
#' tTEscanR_obj <- createObject(
#'     counts = default_tTEscanR_mRNA_data, assay = "mRNA",
#'     meta.data = list(default_tTEscanR_metadata, "tissue"),
#'     meta.data.ids = list("ConditionsLabels", "CorrectionFactor")
#' )
#' tTEscanR_obj <- computeCodonUsage(
#'     object = tTEscanR_obj, species = "hg38",
#'     additional_metrics = FALSE, reduce = 1000
#' )
#'
#' codon_usage <- getAssay(tTEscanR_obj, "CodonUsage")
#' codon_background <- rowSums(codon_usage) / sum(rowSums(codon_usage))
#' codons_to_AA <- featuresToAA(
#'     data = names(codon_background),
#'     notation_from = "codon", notation_to = "aa"
#' )
#' codon_background <- data.frame(
#'     group = codons_to_AA, codon = names(codon_background),
#'     freq = as.numeric(codon_background), row.names = NULL
#' )
#' significance <- obtainSignificance(
#'     dist = permutation_test, value = codon_background
#' )
#'
#' plotPermutation(permut_data = permutation_test, sig_data = significance)
plotPermutation <- function(permut_data, sig_data, color_palette = NULL,
    save_format = NULL, out_name = NULL, out_directory = NULL,
    show_legend = "none", add_titles = TRUE, verbose = TRUE) {
    annot_results <- getAnnotData(permut_data, sig_data)
    observed_data <- annot_results$sig_data
    label_data <- annot_results$annot_data_labels
    val_range <- range(permut_data$freq, na.rm = TRUE)
    bin_width <- (val_range[2] - val_range[1]) / 25 # Size of the bins
    if (is.infinite(bin_width) || bin_width == 0) bin_width <- 0.001
    plot <- ggplot2::ggplot() + # Background of the histogram (permut_data)
        ggplot2::geom_histogram(
            data = permut_data, ggplot2::aes(x = .data$freq),
            binwidth = bin_width, fill = "lightgrey", color = "black"
        ) +
        ggplot2::geom_vline( # The observed values (observed_data / sig_data)
            data = observed_data,
            ggplot2::aes(xintercept = .data$freq, color = .data$group),
            linetype = "dashed", linewidth = 0.8
        ) +
        ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(x = -Inf, y = Inf, label = .data$combined_label),
            hjust = -0.1, vjust = 1.2, size = 3
        ) + # Facet by codon
        ggplot2::facet_wrap(ggplot2::vars(.data$codon), scales = "free_y") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = show_legend)
    unique_groups <- unique(observed_data$group)
    checked_palette <- checkPaletteNames(
        color_palette = color_palette, actual_categories = unique_groups
    )
    plot <- plot + getSafeColorScale(
        n_colors = length(unique_groups), color_palette = checked_palette,
        aes_type = "color"
    )
    if (isTRUE(add_titles)) { # Adding standard titles to the plot
        plot <- plot + ggplot2::labs(
            title = "Codon Frequency Permutation Test",
            subtitle = "Grey bars = Null Distribution; Dashed lines = Observed",
            x = "Frequency", y = "Count"
        )
    }
    if (!is.null(save_format)) { # Save the ggplot
        savePlot(
            plot = plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory, verbose = verbose
        )
    }
    return(plot)
}
