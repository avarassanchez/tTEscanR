#' @importFrom rlang :=

checkPaletteNames <- function(color_palette, actual_categories) {
    if (is.list(color_palette)) color_palette <- unlist(color_palette)

    if (!is.null(color_palette) && !is.null(names(color_palette))) {
        if (!any(names(color_palette) %in% unique(actual_categories))) {
            warning(
                "Provided 'color_palette' names do not match data categories. ",
                "Applying colors sequentially instead."
            )
            color_palette <- unname(color_palette)
        }
    }
    return(color_palette)
}

getSafeColorScale <- function(n_colors, color_palette = NULL,
    aes_type = "fill") {
    if (is.list(color_palette)) color_palette <- unlist(color_palette)

    use_custom <- !is.null(color_palette) && (length(color_palette) >= n_colors)
    use_internal <- n_colors <= length(gradual_groups_35)

    if (use_custom) {
        if (aes_type == "fill") {
            return(ggplot2::scale_fill_manual(values = color_palette))
        }
        if (aes_type == "color") {
            return(ggplot2::scale_color_manual(values = color_palette))
        }
    }

    if (use_internal) {
        if (!is.null(color_palette)) {
            warning(
                "The provided 'color_palette' does not have enough colors for ",
                n_colors, " categories. Defaulting to internal palette."
            )
        }
        selected_colors <- unname(gradual_groups_35[seq_len(n_colors)])
        if (aes_type == "fill") {
            return(ggplot2::scale_fill_manual(values = selected_colors))
        }
        if (aes_type == "color") {
            return(ggplot2::scale_color_manual(values = selected_colors))
        }
    }

    if (!is.null(color_palette)) {
        warning(
            "The provided 'color_palette' does not have enough colors for ",
            n_colors, " categories. Defaulting to viridis palette."
        )
    }
    if (aes_type == "fill") {
        return(ggplot2::scale_fill_viridis_d(option = "viridis"))
    }
    if (aes_type == "color") {
        return(ggplot2::scale_color_viridis_d(option = "viridis"))
    }
}

savePlot <- function(plot, save_format, out_name, out_directory, width = 8,
    height = 6, verbose) {
    save_format <- tolower(save_format)

    if (save_format %in% c("png", "pdf")) {
        output_file <- getOutputName(
            action = "plot", out_name = out_name,
            out_directory = out_directory,
            save_format = save_format,
            verbose = verbose
        )
        ggplot2::ggsave(
            filename = output_file, plot = plot,
            width = width, height = height
        )
    } else {
        warning(
            "The input 'save_format' ('", save_format,
            "') is not recognized.\n", "Supported formats: png or pdf.\n",
            "The plot was generated but not stored."
        )
    }
}

getOutputName <- function(action, out_name, out_directory, save_format,
    verbose) {
    if (is.null(out_directory)) {
        if (verbose) {
            message(
                "- Parameter 'out_directory' has not been specified.\n",
                "The file will be stored in the current working directory."
            )
        }
        out_directory <- getwd()
    }

    if (is.null(out_name)) {
        if (verbose) {
            message(
                "- Parameter 'out_name' has not been specified.\n",
                "A standard name will be used."
            )
        }

        if (action == "plot") {
            out_name <- "distribution_plot"
        } else {
            out_name <- "tRNA_expression_matrix"
        }
    }

    ## Check if the output name already contains the format extension
    if (!endsWith(tolower(out_name), paste0(".", tolower(save_format)))) {
        out_name <- paste0(out_name, ".", save_format)
    }

    output_file <- file.path(out_directory, out_name)
    if (verbose) {
        message("- The generated ", action, " will be save in: ", output_file)
    }
    return(output_file) # Returns the output name.
}

getAnnotData <- function(permut_data, sig_data) {
    ## Check that the expected columns are present in the data
    checkValueInData(
        param = "permut_data", observed = colnames(permut_data),
        expected = c("codon", "freq")
    )

    checkValueInData(
        param = "sig_data", observed = colnames(sig_data),
        expected = c(
            "group", "codon", "freq", "p_val", "tail", "p_val_adj", "sig_adj"
        )
    )

    labels_df <- sig_data %>%
        dplyr::mutate(
            label = ifelse(
                .data$sig_adj, paste0("p=", signif(.data$p_val_adj, 2), "*"),
                paste0("p=", signif(.data$p_val_adj, 2))
            )
        )

    annot_data_labels <- labels_df %>%
        dplyr::group_by(.data$codon) %>%
        dplyr::summarise(
            combined_label = paste(.data$label, collapse = " | "),
            .groups = "drop"
        )

    return(list(sig_data = sig_data, annot_data_labels = annot_data_labels))
}

checkValueInData <- function(param, observed, expected) {
    if (is.null(observed)) {
        return()
    }
    present <- all(observed %in% expected)

    if (!present) {
        missing_vals <- observed[!observed %in% expected]

        stop(
            "Invalid '", param, "' provided: '",
            paste(missing_vals, collapse = "', '"), "'.\n",
            "Available columns/options are: '",
            paste(expected, collapse = "', '"), "'."
        )
    }
}

checkDataInLongFormat <- function(data) {
    if (!is.data.frame(data)) {
        stop("The input 'data' must be a data.frame or tibble.")
    }

    ## One numeric column
    num_cols <- vapply(data, is.numeric, FUN.VALUE = logical(1))

    ## One categorical column
    cat_cols <- vapply(
        data, function(x) is.character(x) || is.factor(x),
        FUN.VALUE = logical(1)
    )

    ## Will return TRUE if the data is in long format
    is_long <- any(cat_cols) & any(num_cols)

    if (!is_long) {
        stop(
            "The input data needs to be in long format, ",
            "containing at least one numeric and one categorical colummn."
        )
    }
}

generateDistPlot <- function(level, target, data, x_axis, y_axis, color,
    add_titles, show_legend, ncols, bar, facet, add_stats) {
    if (level %in% c("jitter", "dot")) {
        p <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(
            x = .data[[x_axis]], y = .data[[y_axis]], color = .data[[target]]
        ))
    } else {
        p <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(
            x = .data[[x_axis]], y = .data[[y_axis]], fill = .data[[target]]
        ))
    }
    if (level == "jitter") p <- p + ggplot2::geom_jitter(size = 0.5)
    if (level == "barplot") {
        p <- p + ggplot2::geom_col(stat = "identity", position = bar)
    }
    if (level == "boxplot") p <- p + ggplot2::geom_boxplot()
    if (level == "dot") {
        p <- p + ggplot2::aes(group = .data[[target]]) + ggplot2::geom_point()
    }
    p <- p + ggplot2::theme_bw() + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = show_legend
    )
    if (!is.null(facet)) {
        p <-p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet]]), ncol = ncols)
    }
    if (level == "dot") {
        p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())
    }
    current_aes <- if (level %in% c("jitter", "dot")) "color" else "fill"
    unique_cat <- unique(data[[target]])
    palette <- checkPaletteNames(color, unique_cat)
    p <- p + getSafeColorScale(
        n_colors = length(unique_cat),
        color_palette = palette, aes_type = current_aes
    )
    if (isTRUE(add_titles)) {
        p <- p + ggplot2::labs(x = x_axis, y = paste(x_axis, "usage counts"))
    }
    if (isTRUE(add_stats)) { # Statistical Logic
        stats_results <- addSignificanceDist(
            data = data, plot = p, x_axis_col = x_axis, y_axis_col = y_axis,
            target = target, facet_col = facet
        )
        return(list(plot = stats_results$plot, stats = stats_results$sig_table))
    } else {
        return(list(plot = p))
    }
}

addSignificanceDist <- function(data, plot, x_axis_col, y_axis_col, target,
    facet_col) {
    sig_table <- computeBoxplotSignificance(
        data,
        x_col = x_axis_col, y_col = y_axis_col, target_col = target,
        group_col = facet_col
    )

    if (!is.null(sig_table) && nrow(sig_table) > 0) {
        ## Extract the position of the significance info
        y_max <- max(data[[y_axis_col]], na.rm = TRUE) * 1.05
        plot <- plot + ggpubr::stat_pvalue_manual(
            sig_table,
            x = x_axis_col, label = "p_signif",
            y.position = y_max, hide.ns = TRUE
        )
    }

    return(list(sig_table = sig_table, plot = plot))
}

drawBarCountsPlot <- function(data, var_numerical, var_categorical, var_color,
    color_palette, show_legend, order, x_limits, facet_col) {
    ## Set a specific order of the features (if available)
    if (!is.null(order)) {
        data[[var_categorical]] <- factor(
            data[[var_categorical]],
            levels = order
        )
    } else {
        data[[var_categorical]] <- stats::reorder(
            data[[var_categorical]], data[[var_numerical]]
        )
    }

    plot <- ggplot2::ggplot(data, ggplot2::aes(
        x = .data[[var_numerical]], y = .data[[var_categorical]],
        fill = .data[[var_color]]
    )) +
        ggplot2::geom_col(color = "black")

    if (!is.null(x_limits)) {
        plot <- plot + ggplot2::scale_x_continuous(limits = x_limits)
    }

    if (!is.null(facet_col)) {
        plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_col]]))
    }

    n_colors <- length(unique(data[[var_color]]))
    plot <- plot + getSafeColorScale(
        n_colors = n_colors, color_palette = color_palette, aes_type = "fill"
    )

    plot <- plot + ggplot2::theme_bw() + ggplot2::theme(
        legend.position = show_legend
    ) + ggplot2::labs(x = "", y = "")

    return(plot)
}

drawDonutPlot <- function(data, var_numerical, var_categorical, color_palette,
    show_legend) {
    total <- sum(data[[var_numerical]], na.rm = TRUE) # Compute percentages
    if (total == 0) {
        stop(
            "Sum of '", var_numerical, "' is zero, cannot compute fractions."
        )
    }
    ## Compute the cumulative percentages and the position of the label
    data$fraction <- data[[var_numerical]] / total
    data$ymax <- cumsum(data$fraction)
    data$ymin <- c(0, utils::head(data$ymax, -1))
    data$labelPosition <- (data$ymax + data$ymin) / 2

    if (show_legend == "none") {
        data$label <- paste(
            data[[var_categorical]], ":", round(data[[var_numerical]], 2)
        )
    } else {
        data$label <- as.character(round(data[[var_numerical]], 2))
    }
    ## Adjusting the sizes automatically based on the number of categories
    n_cat <- length(unique(data[[var_categorical]]))
    circle_min <- 1 + log10(n_cat) * 0.5 # Radius scales with log of categories
    donut_width <- ifelse(n_cat <= 3, 2.5, 1.5)
    circle_max <- circle_min + donut_width # before + 1
    label_size <- min(5, max(1.5, 6 - 0.1 * n_cat)) # Size decreases n_cat grows
    label_position <- 4 + 0.2 * n_cat # Label position shifts outward with n_cat

    plot <- ggplot2::ggplot(data, ggplot2::aes(
        ymax = .data$ymax, ymin = .data$ymin, xmax = circle_max,
        xmin = circle_min, fill = .data[[var_categorical]]
    )) +
        ggplot2::geom_rect() +
        ggplot2::geom_label(
            x = label_position, ggplot2::aes(
                y = .data$labelPosition, label = .data$label
            ), size = label_size, nudge_x = 2, nudge_y = 0
        ) +
        ggplot2::coord_polar(theta = "y") +
        ggplot2::xlim(c(-1, label_position)) +
        ggplot2::theme_void()
    plot <- plot + getSafeColorScale(
        n_colors = n_cat, color_palette = color_palette, aes_type = "fill"
    )
    plot <- plot + ggplot2::theme(legend.position = show_legend)
    return(plot)
}

computeRings <- function(data, var_color, var_cat, var_num,
    norm, n_rings, zoom) {
    ## Standardize dummy groups if needed
    if (is.null(var_color) || var_color == var_cat) {
        data$.__dummy_group__ <- "all"
        var_color <- ".__dummy_group__"
    }
    if (!is.numeric(data[[var_num]])) { # Ensure column is numeric
        data[[var_num]] <- as.numeric(as.character(data[[var_num]]))
    }

    tmp <- data %>% # Pivot data to get wide format for ggradar
        dplyr::select(dplyr::all_of(c(
            var_color, var_cat, var_num
        ))) %>% tidyr::pivot_wider(
            id_cols = dplyr::all_of(var_color), names_from = dplyr::all_of(
                var_cat
            ), values_from = dplyr::all_of(var_num),
            values_fn = sum, values_fill = 0
        )
    num_cols <- setdiff(colnames(tmp), var_color)

    if (isTRUE(norm)) { # Apply normalization if requested
        ## Ensure is a vector for division
        row_totals <- rowSums(tmp[, num_cols, drop = FALSE], na.rm = TRUE)
        row_totals[row_totals == 0] <- 1 # Prevent division by zero
        tmp[num_cols] <- tmp[num_cols] / row_totals
    }
    all_numeric_values <- unlist(tmp[num_cols]) # Get numbers to calculate dist.
    if (isTRUE(zoom)) {
        ## Calculate the 95th percentile
        max_val <- stats::quantile(all_numeric_values, 0.95, na.rm = TRUE)
        max_val <- as.numeric(max_val) + 0.00001
    } else {
        if (isTRUE(norm)) {
            max_val <- 1.00001
        } else {
            max_val <- max(all_numeric_values, na.rm = TRUE) + 0.00001
        }
    }
    ring_values <- seq(0, max_val, length.out = n_rings) # Equally spaced ring
    if (isTRUE(norm)) { # Create labels based on normalization
        labels_rings <- paste0(round(ring_values * 100, 1), "%")
    } else {
        labels_rings <- as.character(round(ring_values, 2))
    }
    return(list(max_val = max_val, labels_rings = labels_rings))
}

drawRadarPlot <- function(data, var_color, var_categorical, var_numerical,
    normalize, zoom, title, add_titles, show_legend, global_max_val = NULL,
    labels_rings = NULL, color_palette) {
    if (is.null(var_color) || var_color == var_categorical) {
        data$.__dummy_group__. <- "all"
        var_color <- ".__dummy_group__."
    }
    if (!is.numeric(data[[var_numerical]])) {
        data[[var_numerical]] <- as.numeric(as.character(data[[var_numerical]]))
    }
    p <- data %>%
        dplyr::select(dplyr::all_of(
            c(var_color, var_categorical, var_numerical)
        )) %>%
        tidyr::pivot_wider(
            id_cols = dplyr::all_of(var_color), names_from =
                dplyr::all_of(var_categorical), values_from = dplyr::all_of(
                var_numerical
            ), values_fn = sum, values_fill = 0
        )
    if (isTRUE(normalize)) {
        num_cols <- setdiff(colnames(p), var_color)
        p[num_cols] <- p[num_cols] / rowSums(p[num_cols], na.rm = TRUE)
    }
    num_cols <- setdiff(colnames(p), var_color)
    p[num_cols] <- lapply(p[num_cols], function(x) {
        ifelse(x > global_max_val, global_max_val, x)
    })
    if (!is.null(color_palette)) {
        n_groups <- nrow(p)
        if (length(color_palette) < n_groups) {
            warning("Uncomplete 'color_palette', using internal palette.")
            color_palette <- NULL
        } else {
            color_palette <- color_palette[seq_len(n_groups)]
        }
    }
    radar_plot <- ggradar::ggradar(
        p,
        background.circle.colour = "white", group.colours = color_palette,
        grid.min = 0, values.radar = labels_rings, grid.max = global_max_val,
        grid.mid = global_max_val / 2, gridline.min.linetype = 1,
        gridline.mid.linetype = 1, gridline.max.linetype = 1,
        group.line.width = 0.8, group.point.size = 2,
        legend.position = show_legend
    )
    if (isTRUE(add_titles)) radar_plot <- radar_plot + ggplot2::ggtitle(title)
    return(radar_plot)
}

generateProportionPlot <- function(level, data, var_numerical, var_categorical,
    var_color, color_palette, zoom, show_legend, order, x_limits, facet_col,
    normalize, add_titles, n_rings) {
    base_args <- list(
        data = data, var_numerical = var_numerical, show_legend = show_legend,
        var_categorical = var_categorical, color_palette = color_palette
    )
    res <- list(plot = NULL, legend = NULL)
    if (level == "donut") res$plot <- do.call(drawDonutPlot, base_args)
    if (level == "bar") {
        bar_args <- c(base_args, list(
            var_color = var_color, order = order, x_limits = x_limits,
            facet_col = facet_col
        ))
        res$plot <- do.call(drawBarCountsPlot, bar_args)
    }
    if (level == "radar") {
        cond <- sort(unique(data[[var_categorical]]))
        condition_index <- stats::setNames(seq_along(cond), cond)
        data[[var_categorical]] <- condition_index[data[[var_categorical]]]
        res$legend <- data.frame(Index = seq_along(cond), Condition = cond)
        rings_info <- computeRings(
            data = data, var_color = var_color, n_rings = n_rings, zoom = zoom,
            norm = normalize, var_cat = var_categorical, var_num = var_numerical
        )
        radar_args <- c(base_args[names(base_args) != "data"], list(
            var_color = var_color, add_titles = add_titles,
            labels_rings = rings_info$labels_rings, normalize = normalize,
            zoom = zoom, global_max_val = rings_info$max_val
        ))
        if (is.null(facet_col)) {
            res$plot <- do.call(drawRadarPlot, c(list(
                data = data, title = paste("Radar plot by", var_categorical)
            ), radar_args))
        } else {
            data_split <- split(data, data[[facet_col]])
            plots_list <- lapply(names(data_split), function(name) {
                do.call(drawRadarPlot, c(
                    list(data = data_split[[name]], title = name), radar_args
                ))
            })
            res$plot <- patchwork::wrap_plots(plots_list)
        }
    }
    if (is.null(res$legend)) {
        return(list(plot = res$plot))
    } else {
        return(res)
    }
}

significanceSymbol <- function(pvalue) {
    if (!is.numeric(pvalue)) {
        return(rep(NA_character_, length(pvalue)))
    }

    symbols <- as.character(stats::symnum(
        pvalue,
        corr = FALSE, na = FALSE,
        cutpoints = c(-Inf, 0.001, 0.01, 0.05, Inf),
        symbols = c("***", "**", "*", "ns")
    ))

    symbols[is.na(pvalue)] <- NA_character_

    return(symbols)
}

computeTEsignificance <- function(merged, x_col = "class", target, score,
    group = NULL) {
    target_levels <- sort(unique(merged[[target]]))
    if (length(target_levels) < 2) { # Need at least 2 levels to compare
        return(NULL)
    }

    ## Generate all possible pairs: 2 columns, N rows
    comparisons <- utils::combn(target_levels, 2)
    ## Define the grouping variables (facets only)
    group_vars <- if (!is.null(group)) group else NULL

    run_pairwise_stats <- function(df) {
        ## Loop through each pair in the 'comparisons' matrix
        purrr::map_df(seq_len(ncol(comparisons)), function(i) {
            g1_name <- comparisons[1, i]
            g2_name <- comparisons[2, i]
            p_val <- tryCatch(
                {
                    val1 <- df[[score]][df[[target]] == g1_name]
                    val2 <- df[[score]][df[[target]] == g2_name]
                    stats::wilcox.test(val1, val2)$p.value
                },
                error = function(e) NA_real_
            )
            data.frame(
                group1 = g1_name, group2 = g2_name, p_value = p_val,
                comparison = paste0(g1_name, "_vs_", g2_name)
            )
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
            dplyr::mutate(
                p_signif = significanceSymbol(.data$p_value),
                !!x_col := .data$group2
            ) %>%
            dplyr::filter(!is.na(.data$p_value))
    }
    return(sig_table)
}

computeBoxplotSignificance <- function(merged, x_col, y_col, target_col,
    group_col = NULL) {
    if (length(unique(merged[[target_col]])) != 2) {
        return(NULL)
    }
    target_levels <- sort(unique(merged[[target_col]]))

    ## Create formula dynamically
    form <- stats::as.formula(paste(y_col, "~", target_col))

    sig_table <- merged %>% # Perform the Wilcoxon test
        dplyr::group_by(dplyr::across(dplyr::all_of(c(x_col, group_col)))) %>%
        dplyr::summarise(
            p_value = tryCatch(
                {
                    stats::wilcox.test(
                        form,
                        data = dplyr::pick(dplyr::everything())
                    )$p.value
                },
                error = function(e) NA_real_
            ), .groups = "drop"
        )

    ## Filter out NA values (failed tests)
    sig_table <- sig_table %>% dplyr::filter(!is.na(.data$p_value))

    if (nrow(sig_table) > 0) {
        sig_table <- sig_table %>%
            dplyr::mutate(
                p_value = stats::p.adjust(.data$p_value, method = "BH")
            ) %>%
            dplyr::mutate(
                p_signif = significanceSymbol(.data$p_value),
                group1 = as.character(target_levels[1]),
                group2 = as.character(target_levels[2])
            )
    } else {
        return(NULL) # All tests failed
    }
    return(sig_table)
}
