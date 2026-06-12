RunIterations <- function(num_iter, data, anticodon, supply, cuts, slope, rho) {
    cutoff_res <- list() # Store optimal results of each iteration

    pb <- utils::txtProgressBar(min = 0, max = num_iter, style = 3) # Timer

    for (i in seq_len(num_iter)) { # Obtain the optimal cutoff per iteration
        out <- Iterate_tRNACutoff(
            data = data, anticodon = anticodon, supply = supply, cutoffs = cuts,
            slope_threshold = slope, rho_threshold = rho
        )

        cutoff_res[[i]] <- out
        utils::setTxtProgressBar(pb, i) # Increase the progress bar
    }

    close(pb)

    return(cutoff_res)
}

TransformCounts <- function(data) {
    data_long <- TransformFormat(
        data = data, normalize = FALSE, rownames_to_column = "features",
        names_to = "conditions", values_to = "counts"
    )

    data_long <- dplyr::filter(data_long, .data$counts > 0)

    # Transforms every count into an independent row
    data_long <- tidyr::uncount(data_long, weights = .data$counts)

    return(data_long)
}

ReferenceObject <- function(data, compute_aa) {
    # Create object
    ref <- CreateObject(counts = data, assay = "tRNA", verbose = FALSE)

    # Compute anticodon usage
    ref <- ComputeAnticodonUsage(ref, verbose = FALSE)
    anticodon <- GetReference(getAssay(ref, "AnticodonUsage"))

    # Compute amino acid supply if applicable
    if (isTRUE(compute_aa)) {
        ref <- ComputeAAUsage(ref, level = "supply", verbose = FALSE)
        supply <- GetReference(getAssay(ref, "AASupply"))
    } else {
        supply <- NULL
    }

    return(list(anticodon = anticodon, supply = supply))
}

GetReference <- function(assay_slot) {
    value <- rowSums(assay_slot)
    return(value / sum(value))
}

CutoffMatrix <- function(data, cutoff) {
    data <- dplyr::slice_sample(data, n = cutoff)
    data <- data %>%
        dplyr::count(.data$features, .data$conditions) %>%
        tidyr::pivot_wider(
            names_from = .data$conditions,
            values_from = .data$n, values_fill = 0
        )

    matrix <- data %>%
        tibble::column_to_rownames("features") %>%
        as.matrix()

    return(matrix)
}

FilteringCutoffs <- function(data, cutoff, aa_supply) {
    data_above <- CutoffMatrix(data = data, cutoff = cutoff)

    object_above <- CreateObject(
        counts = data_above, assay = "tRNA", verbose = FALSE
    )
    object_above <- ComputeAnticodonUsage(
        object = object_above, verbose = FALSE
    )

    total_anticodon <- rowSums(getAssay(object_above, "AnticodonUsage"))
    total_anticodon <- total_anticodon / sum(total_anticodon)

    if (!is.null(aa_supply)) {
        object_above <- ComputeAAUsage(
            object = object_above, level = "supply", verbose = FALSE
        )
        total_supply <- rowSums(getAssay(object_above, "AASupply"))
        total_supply <- total_supply / sum(total_supply)
    } else {
        total_supply <- NULL
    }

    return(list(anticodon = total_anticodon, supply = total_supply))
}

ComputeCorrelations <- function(data, ref_anticodon, ref_supply, cutoffs) {
    cor_results <- data.frame(
        cutoff = cutoffs, anticodon_spearman = NA_real_,
        supply_spearman = NA_real_,
        total_anticodon = NA_real_, total_supply = NA_real_
    )

    for (i in seq_along(cutoffs)) {
        cut <- cutoffs[i]
        obj <- FilteringCutoffs(
            data = data, cutoff = cut, aa_supply = ref_supply
        )

        if (length(obj$anticodon) == 0) next
        if (!is.null(ref_supply) && length(obj$supply) == 0) next

        anticodon_vec <- obj$anticodon[names(ref_anticodon)]
        anticodon_vec[is.na(anticodon_vec)] <- 0

        cor_results$anticodon_spearman[i] <- stats::cor(
            anticodon_vec, ref_anticodon,
            method = "spearman"
        )
        cor_results$total_anticodon[i] <- sum(obj$anticodon)

        if (!is.null(ref_supply)) {
            supply_vec <- obj$supply[names(ref_supply)]
            supply_vec[is.na(supply_vec)] <- 0

            cor_results$supply_spearman[i] <- stats::cor(
                supply_vec, ref_supply,
                method = "spearman"
            )
            cor_results$total_supply[i] <- sum(obj$supply)
        }
    }

    return(cor_results)
}

SelectionCutoff <- function(
    cor_long, slope_threshold = 0.001, rho_threshold = 0.9
) {
    # Retain those cutoff that give a high correlation score
    cor_filtered <- cor_long %>% dplyr::filter(
        abs(.data$spearman_corr) >= rho_threshold
    )

    # Compute the slope on the remaining points
    cor_stability <- cor_filtered %>%
        dplyr::group_by(.data$type) %>%
        dplyr::arrange(.data$cutoff, .by_group = TRUE) %>%
        dplyr::mutate(
            delta_corr = abs(
                .data$spearman_corr - dplyr::lag(.data$spearman_corr)
            ),
            slope_ratio = .data$delta_corr / abs(
                .data$cutoff - dplyr::lag(.data$cutoff)
            )
        )

    stable_points <- cor_stability %>%
        dplyr::filter(
            is.finite(.data$slope_ratio), .data$slope_ratio < slope_threshold
        )

    if (nrow(stable_points) == 0) {
        warning(
            "No stable points below the slope threshold = ", slope_threshold
        )
        if ("supply_spearman" %in% cor_long$type) {
            return(tibble::tibble(
                type = c("anticodon_spearman", "supply_spearman"),
                optimal_cutoff = NA_real_
            ))
        } else {
            return(tibble::tibble(
                type = c("anticodon_spearman"), optimal_cutoff = NA_real_
            ))
        }
    }
    optimal_cutoff <- stable_points %>%
        dplyr::group_by(.data$type) %>%
        dplyr::summarise(optimal_cutoff = min(.data$cutoff), .groups = "drop")

    return(optimal_cutoff)
}

Iterate_tRNACutoff <- function(
    data, anticodon, supply, cutoffs, slope_threshold, rho_threshold
) {
    cor_results <- ComputeCorrelations(
        data = data, ref_anticodon = anticodon,
        ref_supply = supply, cutoffs = cutoffs
    )

    if (is.null(supply)) {
        cor_long <- cor_results %>% tidyr::pivot_longer(
            cols = c("anticodon_spearman"), names_to = "type",
            values_to = "spearman_corr"
        )
    } else {
        cor_long <- cor_results %>% tidyr::pivot_longer(
            cols = c("anticodon_spearman", "supply_spearman"),
            names_to = "type", values_to = "spearman_corr"
        )
    }
    optimal_cutoff <- SelectionCutoff(
        cor_long = cor_long, slope_threshold = slope_threshold,
        rho_threshold = rho_threshold
    )

    return(optimal_cutoff)
}

CorrelationCutoffPlot <- function(
    data, add_titles = TRUE, save_format = NULL, out_name = NULL,
    out_directory = NULL, show_legend = "bottom"
) {
    cor_long <- data %>%
        tidyr::pivot_longer(
            cols = c("anticodon_spearman", "supply_spearman"),
            names_to = "type", values_to = "spearman_corr"
        )

    # Defining the plot background
    plot <- ggplot2::ggplot(cor_long, ggplot2::aes(
        x = .data$cutoff, y = .data$spearman_corr, color = .data$type
    )) +
        ggplot2::geom_point(size = 2, alpha = 0.8) +
        ggplot2::scale_color_manual(
            values = c(
                "anticodon_spearman" = "#034e7b", "supply_spearman" = "#a6bddb"
            ),
            labels = c("Anticodon", "Supply")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = show_legend)

    if (isTRUE(add_titles)) {
        plot <- plot + ggplot2::labs(
            x = "Sample size (number of cuts)", y = "Spearman Correlation",
            color = "Feature Type",
            title = "Correlation vs Reference Across Cutoff Thresholds"
        )
    }

    if (!is.null(save_format)) { # Save the ggplot
        SavePlot(
            plot = plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory
        )
    }

    return(plot)
}

SelectionCutoffPlot <- function(
    data, add_titles = TRUE, save_format = NULL, out_name = NULL,
    out_directory = NULL, show_legend = "bottom"
) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(
        x = .data$optimal_cutoff, fill = .data$type
    )) +
        ggplot2::geom_bar(position = "identity") +
        ggplot2::facet_wrap(~type, scales = "free") +
        ggplot2::theme_bw() + # coord_cartesian(ylim = c(0, 150)) +
        ggplot2::scale_fill_manual(values = c("#034e7b", "#a6bddb")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90, vjust = 0.5, hjust = 1
        ))

    if (isTRUE(add_titles)) { # Add titles
        plot <- plot + ggplot2::labs(
            title = "Distribution of Cutoff Values", x = "Cutoff", y = "Density"
        )
    }

    if (!is.null(save_format)) { # Save the ggplot
        SavePlot(
            plot = plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory
        )
    }

    return(plot)
}

CumulativeSelectionCutoffPlot <- function(
    data, cutoffs, add_titles = TRUE, save_format = NULL, out_name = NULL,
    out_directory = NULL, show_legend = "bottom"
) {
    conditions_per_threshold <- c()

    for (i in cutoffs) {
        conditions <- length(which(data$total_counts > i))
        conditions_per_threshold <- c(conditions_per_threshold, conditions)
    }

    conditions_per_threshold_table <- data.frame(
        cuts = cutoffs, conditions = conditions_per_threshold
    )

    plot <- ggplot2::ggplot(conditions_per_threshold_table, ggplot2::aes(
        x = .data$cuts, y = .data$conditions
    )) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(
            xintercept = 5000, linetype = "dashed",
            color = "#636363", linewidth = 0.8
        ) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::theme_bw()

    if (isTRUE(add_titles)) { # Add titles
        plot <- plot + ggplot2::labs(
            title = "Distribution of Cutoff Values", x = "Cutoff", y = "Density"
        )
    }

    if (!is.null(save_format)) { # Save the ggplot
        SavePlot(
            plot = plot, save_format = save_format, out_name = out_name,
            out_directory = out_directory
        )
    }

    return(plot)
}
