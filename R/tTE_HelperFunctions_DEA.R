computeSizeCorrection <- function(data, metadata, batch = NULL, reduce = 100,
    verbose = FALSE) {
    checkDataFrame(data = data) # Evaluate the input parameter: 'data'
    checkDataFrame(data = metadata) # Evaluate parameters: 'metadata' & 'batch'

    if (!(batch %in% colnames(metadata))) {
        stop("The correction factor ('batch') was not found in the 'metadata'.")
    }

    if (verbose) {
        message(
            "- The 'data', 'metadata' and 'batch' parameters have been ",
            "properly loaded.\n", "- Size-correcting the counts using DESeq2.",
            "\n- Filtering (if necessary) the 'data' and 'metadata' for common",
            " features."
        )
    }

    ## Check consistency in conditions included in the data and metadata
    filtered <- filterByMetadata(
        data = data, metadata = metadata, verbose = verbose
    )
    DESeq2_run <- computeDESeq2(
        data = filtered$data, metadata = filtered$metadata,
        batch = batch, reduce = reduce, verbose = verbose
    )

    ## Extracting the normalized data from the DESeq2 object
    if (verbose) message("- Extracting the size-corrected counts.")
    size_corrected_output_matrix <- DESeq2::counts(
        DESeq2_run,
        normalized = TRUE
    )
    checkDataFrame(data = size_corrected_output_matrix)

    return(size_corrected_output_matrix) # Return the size-corrected matrix
}

computeAllPairwiseComp <- function(dds, factor_name, padj_threshold = 0.05,
    verbose) {
    if (!factor_name %in% colnames(SummarizedExperiment::colData(dds))) {
        stop("Factor '", factor_name, "' not found in colData.")
    }

    factor_data <- SummarizedExperiment::colData(dds)[[factor_name]]
    if (!is.factor(factor_data)) {
        stop(
            "The column '", factor_name,
            "' must be a factor in the DESeqDataSet."
        )
    }

    levels_factor <- levels(factor_data)
    if (length(levels_factor) < 2) {
        stop("Factor must have at least 2 levels to compute contrasts.")
    }

    ## Generate all pairwise combinations
    combn_matrix <- utils::combn(levels_factor, 2)
    results_list <- list()

    for (i in seq_len(ncol(combn_matrix))) {
        level1 <- combn_matrix[1, i]
        level2 <- combn_matrix[2, i]
        contrast_name <- paste0(level1, "_vs_", level2)

        if (verbose) message("Extracting contrasts: ", contrast_name)
        res <- DESeq2::results(
            dds,
            contrast = c(factor_name, level1, level2),
            alpha = padj_threshold
        )
        results_list[[contrast_name]] <- res
    }
    return(results_list)
}

getDesignFormula <- function(batch, condition, verbose) {
    ## Dynamic definition of the design formula
    if (!is.null(batch) && !is.null(condition)) {
        if (batch == condition) {
            design_formula <- stats::as.formula(paste("~", condition))
            if (verbose) {
                message(
                    "- Batch and Condition are identical. ",
                    "Building single-factor model: ~ ", condition
                )
            }
        } else {
            design_formula <- stats::as.formula(
                paste("~", batch, "+", condition) # condition
            )
            if (verbose) {
                message(
                    "- Building model with batch correction: ~ ",
                    batch, " + ", condition
                )
            }
        }
    } else if (!is.null(condition)) {
        design_formula <- stats::as.formula(paste("~", condition))
        if (verbose) message("- Building model: ~ ", condition)
    } else if (!is.null(batch)) {
        design_formula <- stats::as.formula(paste("~", batch))
        if (verbose) {
            message("- Building exploratory model based on batch: ~ ", batch)
        }
    } else {
        design_formula <- stats::as.formula("~ 1")
        if (verbose) {
            message(
                "- No factors provided. Building intercept-only ",
                "exploratory model: ~ 1"
            )
        }
    }
    return(design_formula)
}

computeDESeq2 <- function(data, metadata, condition = NULL, batch = NULL,
    reference = NULL, reduce = 100, verbose = TRUE) {
    data <- as.matrix(data)
    data <- checkIntegerLength(data = data, reduce = reduce, verbose = verbose)
    storage.mode(data) <- "integer" # Force data as integer
    valid <- (Matrix::colSums(data)) > 0 # Remove empty samples
    if (sum(!valid) > 0) { # Prevent DESeq2 crash
        if (verbose) {
            message("- Removing ", sum(!valid), " samples with zero counts.")
        }
        data <- data[, valid, drop = FALSE]
        metadata <- metadata[valid, , drop = FALSE]
    }
    col <- as.data.frame(metadata) # Ensures a standard data.frame
    if (ncol(data) != nrow(col)) {
        stop(
            "Dimension mismatch: Count matrix has ", ncol(data),
            " samples, but metadata has ", nrow(col), " rows."
        )
    }
    if (!is.null(condition)) {
        col <- setConditionParam(
            cond = condition, col = col, ref = reference, verbose = verbose
        )
    }
    if (!is.null(batch)) {
        if (!batch %in% colnames(col)) {
            stop("The batch factor '", batch, "' was not found in 'metadata'.")
        }
        col[[batch]] <- factor(col[[batch]])
    }
    design_formula <- getDesignFormula(
        batch = batch, condition = condition, verbose = verbose
    )
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = data, colData = col, design = design_formula
    )
    if (verbose) message("- Running DESeq2 pipeline...")
    dds <- DESeq2::DESeq(dds)
    return(dds)
}

setConditionParam <- function(cond, col, ref, verbose) {
    if (!cond %in% colnames(col)) {
        stop(
            "The condition factor '", cond, "' was not found in 'metadata'."
        )
    }
    col[[cond]] <- factor(col[[cond]])

    ## Apply reference level to primary condition
    if (!is.null(ref)) {
        if (!(ref %in% levels(col[[cond]]))) {
            stop(
                "The reference level '", ref,
                "' was not found in the 'condition' column."
            )
        }
        col[[cond]] <- stats::relevel(col[[cond]], ref = ref)
        if (verbose) {
            message("- Set reference level to: ", ref, " in ", cond, ".")
        }
    }

    return(col)
}

targetedApproach <- function(DESeq2, verbose = TRUE, fc, padj, sig, condition) {
    if (verbose) {
        message("- Executing a targeted analysis based on '", condition, "'.")
    }
    col_data <- SummarizedExperiment::colData(DESeq2)
    levels_factor <- levels(col_data[[condition]])
    combs <- utils::combn(levels_factor, 2, simplify = FALSE)
    results.list <- list()

    if (verbose) {
        message(
            "    Levels found in 'target' factor: ",
            paste(levels_factor, collapse = ", "),
            "\n    Number of pairwise combinations: ", length(combs)
        )
    }
    for (cmb in combs) {
        level1 <- cmb[1]
        level2 <- cmb[2]
        contrast_name <- paste(level1, "vs", level2, sep = "_")

        if (verbose) message("Processing contrast: ", contrast_name)
        res <- DESeq2::results(
            DESeq2, contrast = c(condition, level1, level2), alpha = padj
        )
        res_df <- as.data.frame(res)

        target_results <- tibble::tibble(
            feature = rownames(res_df), log2_FC = res_df$log2FoldChange,
            p_value = res_df$pvalue,
            neg_log10_adjusted_p_value = -log10(res_df$padj)
        ) %>%
            dplyr::filter(
                !is.na(.data$log2_FC) &
                    !is.na(.data$neg_log10_adjusted_p_value)
            )

        results.list[[paste0("results_", contrast_name)]] <- target_results
        if (verbose) message("- Generating volcano plot.")
        volcano_plot <- generateVolcanoPlot(
            data = target_results, fc_threshold = fc,
            padj_threshold = padj, label_significant = sig
        ) +
            ggplot2::ggtitle(contrast_name)
        results.list[[paste0("volcano_", contrast_name)]] <- volcano_plot
    }
    return(results.list)
}

makeScatterPlot <- function(df, x_col, y_col, color_factor, shape_factor = NULL,
    label_factor = NULL, median_size = 6, highlight_median = FALSE,
    color_palette = NULL, show_legend = "none", point_alpha = 0.5, title = NULL,
    x_axis = NULL, y_axis = NULL) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data[[x_col]], y = .data[[y_col]]
    ))
    if (!is.null(shape_factor)) {
        p <- p + ggplot2::geom_point(ggplot2::aes(
            color = .data[[color_factor]], shape = .data[[shape_factor]]
        ), alpha = if (highlight_median) point_alpha else 1, size = 3)
    } else {
        p <- p + ggplot2::geom_point(
            ggplot2::aes(color = .data[[color_factor]]),
            alpha = if (highlight_median) point_alpha else 1, size = 3
        )
    }
    if (highlight_median) { # Median points and labels
        p <- highlightMedian(
            df = df, color = color_factor, median = median_size,
            x_col = x_col, y_col = y_col, p = p
        )
    }
    if (!is.null(label_factor)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(
            label = .data[[label_factor]]
        ), size = 2.8, show.legend = FALSE)
    }
    if (!is.null(color_palette)) {
        p <- p + ggplot2::scale_color_manual(values = color_palette)
    } else p <- p + ggplot2::scale_color_viridis_d(option = "viridis")
    p <- p + ggplot2::guides(
        color = ggplot2::guide_legend(), shape = ggplot2::guide_legend()
    ) + ggplot2::theme_bw()
    p <- p + ggplot2::labs(
        title = title, color = color_factor, shape = shape_factor,
        x = if (!is.null(x_axis)) {
            sprintf("%s (%.1f%%)", x_col, x_axis)
        } else {
            x_col
        },
        y = if (!is.null(y_axis)) {
            sprintf("%s (%.1f%%)", y_col, y_axis)
        } else {
            y_col
        }
    )
    return(p + ggplot2::theme(legend.position = show_legend))
}

highlightMedian <- function(df, color, x_col, y_col, p, median) {
    medians <- df %>%
        dplyr::group_by(.data[[color]]) %>%
        dplyr::summarise(
            dplyr::across(
                .cols = dplyr::all_of(c(x_col, y_col)),
                .fns = ~ stats::median(.x, na.rm = TRUE)
            ),
            .groups = "drop"
        )

    p <- p + ggplot2::geom_point(
        data = medians, ggplot2::aes(
            x = .data[[x_col]], y = .data[[y_col]],
            color = .data[[color]]
        ),
        size = median, inherit.aes = FALSE
    ) +
        ggrepel::geom_text_repel(
            data = medians, ggplot2::aes(
                x = .data[[x_col]], y = .data[[y_col]],
                label = .data[[color]]
            ),
            fontface = "bold", size = 3, show.legend = FALSE,
            inherit.aes = FALSE
        )

    return(p)
}

runDimReduct <- function(vst, metadata, condition = NULL, show_legend, numPC,
    color = condition, shape = NULL, label = NULL, dim_reduct, median,
    palette = NULL, scale = TRUE, verbose) {
    plots <- list()
    M <- t(vst) # Transpose vst to samples x features
    var <- NULL
    if (dim_reduct == "PCA") {
        low_var <- which(matrixStats::colVars(M, useNames = TRUE) < 1e-5)
        nlow <- length(low_var)
        if (nlow > 0) { # Remove zero variance columns
            M <- M[, -low_var, drop = FALSE]
            if (verbose) message("- Removed ", nlow, " low-variance features.")
        }
        pca <- stats::prcomp(M, scale. = scale, center = TRUE)
        coords <- pca$x[, seq_len(min(numPC, ncol(pca$x))), drop = FALSE]
        colnames(coords) <- paste0("PC", seq_len(ncol(coords)))
        var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
        plots[["ElbowPlot"]] <- produceElbowPlot(data = pca, variance = var)
    } else if (dim_reduct == "UMAP") {
        coords <- tryCatch(
            {
                uwot::umap(M)
            },
            error = function(e) {
                warning("UMAP failed: ", e$message, "\nUsing PCA instead.")
                pca <- stats::prcomp(M, scale. = scale, center = TRUE)$x
                pca[, seq_len(min(2, ncol(pca))), drop = FALSE]
            }
        )
        colnames(coords) <- c("Dim1", "Dim2")
    } else if (dim_reduct == "tSNE") {
        perp <- min(30, max(1, floor((nrow(M) - 1) / 3))) # Dynamic perplexity
        if (verbose) message("- Running t-SNE with perplexity = ", perp)
        coords <- Rtsne::Rtsne(M, perplexity = perp, check_duplicates = FALSE)$Y
        colnames(coords) <- c("Dim1", "Dim2")
    }
    rownames(coords) <- rownames(M) # Restore sample names before merging
    df <- cbind(as.data.frame(coords), metadata)
    common_args <- list(
        df = df, color_factor = color, shape_factor = shape,
        label_factor = label, color_palette = palette,
        highlight_median = median, show_legend = show_legend
    )
    plots <- generateReductPlot(
        plots = plots, dim_reduct = dim_reduct, coords = coords, numPC = numPC,
        common_args = common_args, var = var
    )
    return(plots)
}

generateReductPlot <- function(plots, dim_reduct, coords, numPC,
    common_args, var) {
    if (dim_reduct == "PCA") {
        pc_names <- colnames(coords)[seq_len(numPC)]
        pc_pairs <- utils::combn(pc_names, 2, simplify = FALSE)

        for (pair in pc_pairs) {
            x_pc <- pair[1]
            y_pc <- pair[2]
            pca_args <- c(common_args, list(
                x_col = x_pc, y_col = y_pc,
                title = sprintf("%s vs %s", x_pc, y_pc),
                x_axis = var[which(pc_names == x_pc)],
                y_axis = var[which(pc_names == y_pc)]
            ))
            plots[[paste0(x_pc, "_vs_", y_pc)]] <- do.call(
                makeScatterPlot, pca_args
            )
        }
    } else {
        other_args <- c(common_args, list(
            x_col = colnames(coords)[1], y_col = colnames(coords)[2],
            title = dim_reduct
        ))
        plots[[dim_reduct]] <- do.call(makeScatterPlot, other_args)
    }

    return(plots)
}

generateVolcanoPlot <- function(data, fc_threshold, padj_threshold,
    label_significant) {
    ## Extract the statistical significant data
    sig_data <- data %>% dplyr::filter(
        abs(.data$log2_FC) > fc_threshold,
        .data$neg_log10_adjusted_p_value > -log10(padj_threshold)
    )

    volcano_plot <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_point(
            data = data, mapping = ggplot2::aes(
                x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value
            ), color = "black",
            alpha = 0.5
        ) +
        ggplot2::geom_point(
            data = sig_data, mapping = ggplot2::aes(
                x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value
            ),
            color = "darkred", size = 3
        )

    if (isTRUE(label_significant)) {
        volcano_plot <- volcano_plot +
            ggrepel::geom_text_repel(data = sig_data, ggplot2::aes(
                x = .data$log2_FC, y = .data$neg_log10_adjusted_p_value,
                label = .data$feature
            ), size = 2, max.overlaps = 20)
    }

    volcano_plot <- volcano_plot + ggplot2::geom_hline(
        yintercept = -log10(padj_threshold), linetype = "dotted"
    ) +
        ggplot2::geom_vline(
            xintercept = c(-fc_threshold, fc_threshold),
            linetype = "dotted"
        ) +
        ggplot2::labs(
            x = expression(log[2] ~ fold ~ change),
            y = expression(-log[10] ~ adjusted ~ p ~ value)
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(
            hjust = 0.5, size = 20
        ))

    return(volcano_plot)
}

heatmapFontSize <- function(data, base_size = 10, shrink_factor = 1.5,
    min_size = 2) {
    fontsize_row <- max(min_size, base_size - shrink_factor * log10(nrow(data)))
    fontsize_col <- max(min_size, base_size - shrink_factor * log10(ncol(data)))

    return(list(fontsize_row = fontsize_row, fontsize_col = fontsize_col))
}

produceHeatmapDiffExp <- function(data) {
    Euc_Dists <- stats::dist(t(data)) # Compute Euclidean distances
    dist_matrix <- as.matrix(Euc_Dists)

    ## Sort nodes in hierarchical clustering
    dend <- dendsort::dendsort(stats::as.dendrogram(stats::hclust(Euc_Dists)))

    ## Control the size of the labels in the heatmap
    font_size <- heatmapFontSize(data = dist_matrix)

    ## Produce heatmap that has been hierarchically clustered
    heatmap_grop <- grid::grid.grabExpr({
        Dist_heatmap <- ComplexHeatmap::Heatmap(
            matrix = dist_matrix, name = "Euclidean\nDistance",
            col = viridis::inferno(100), cluster_columns = dend,
            cluster_rows = dend, show_row_dend = TRUE, show_column_dend = TRUE,
            border = FALSE,
            row_names_gp = grid::gpar(fontsize = font_size$fontsize_row),
            column_names_gp = grid::gpar(fontsize = font_size$fontsize_col)
        )
        ComplexHeatmap::draw(Dist_heatmap)
    })

    return(heatmap_grop)
}

produceElbowPlot <- function(data, variance) {
    ## To keep the output cleaner the number of PCs will be limited to 50
    numPC <- min(length(variance), 50)
    ## Generates table with each PC and their variance
    elbow_tibble <- tibble::tibble(
        PC = seq_len(numPC),
        variance_explained = variance[seq_len(numPC)]
    )

    ## Generate plot - reorder to add PCs in the x-axis & variance in the y-axis
    elbow_plot <- ggplot2::ggplot(data = elbow_tibble, ggplot2::aes(
        x = .data$PC, y = .data$variance_explained
    )) +
        ggplot2::theme_bw() +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
        ggplot2::geom_line(
            mapping = ggplot2::aes(group = 1), color = "black",
            linetype = "dashed", linewidth = 0.8
        ) +
        ggplot2::geom_point(color = "black", size = 2) +
        ggplot2::labs(
            title = "Elbow Plot: Variance Explained per PC",
            x = "Principal Component", y = "Variance Explained (%)"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5,
                face = "bold", size = 14
            ),
            axis.text.x = ggplot2::element_text(angle = 0)
        )

    return(elbow_plot) # Outputs the elbow plot
}
