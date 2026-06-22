mock_plot_data <- data.frame(
    SampleID = paste0("S", 1:12),
    Expression = c(1.2, 1.5, 3.4, 3.8, 1.1, 0.9, 4.5, 4.2, 2.1, 2.5, 5.0, 4.8),
    Codon = rep(c("AUU", "AUC", "AUA"), each = 4),
    Condition = rep(c("Control", "Treatment", "Mutant"), times = 4),
    Tissue = rep(c("Liver", "Kidney"), each = 6),
    stringsAsFactors = FALSE
)

mock_overall <- data.frame(
    Codon = c("AUU", "AUC", "AUA", "AUG"),
    Frequency = c(0.4, 0.3, 0.2, 0.1),
    stringsAsFactors = FALSE
)

mock_target <- data.frame(
    Codon = c("AUU", "AUC", "AUA", "AUG"),
    Frequency = c(0.1, 0.5, 0.1, 0.3), # Different values to show comparison
    stringsAsFactors = FALSE
)

mock_prop_data <- data.frame(
    AminoAcid = rep(c("Ala", "Arg", "Asn"), each = 4),
    Codon = c("GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AAU", "AAC", "AAA", "AAG"),
    Frequency = c(0.4, 0.3, 0.2, 0.1, 0.5, 0.2, 0.2, 0.1, 0.6, 0.4, 0.0, 0.0),
    Group = rep(c("TypeA", "TypeB"), times = 6),
    stringsAsFactors = FALSE
)

mock_te_scores <- data.frame(
    SampleID = paste0("Sample_", 1:12),
    tTE = c(0.45, 0.48, 0.52, 0.50, 0.71, 0.75, 0.72, 0.78, 0.31, 0.29, 0.35, 0.33),
    stringsAsFactors = FALSE
)

mock_metadata_scores <- data.frame(
    SampleID = paste0("Sample_", 1:12),
    CellType = rep(c("T-Cell", "B-Cell"), each = 6),
    Treatment = rep(c("Control", "Drug_A", "Drug_B"), times = 4),
    TissueOrigin = rep(c("Spleen", "Blood"), each = 6),
    stringsAsFactors = FALSE
)

mock_corr_data <- data.frame(
    Codon = c("UUU", "UUC", "UUA", "UUG", "CUU", "CUC"),
    Usage_Control = c(0.12, 0.45, 0.22, 0.15, 0.33, 0.55),
    Usage_Treatment = c(0.14, 0.42, 0.20, 0.18, 0.30, 0.59),
    Group = c("Phe", "Phe", "Leu", "Leu", "Leu", "Leu"),
    stringsAsFactors = FALSE
)

set.seed(42)
mock_permut_data <- data.frame(
    codon = rep(c("AUG", "UAA"), each = 100),
    freq = c(rnorm(100, mean = 0.25, sd = 0.03), rnorm(100, mean = 0.10, sd = 0.01)),
    stringsAsFactors = FALSE
)

mock_sig_data <- data.frame(
    codon = c("AUG", "UAA"),
    freq = c(0.34, 0.11),
    group = c("Treatment_A", "Treatment_A"),
    sig_adj = c(0.002, 0.045),
    p_val_adj = c(0.002, 0.045), # Added to prevent the mutate evaluation crash
    stringsAsFactors = FALSE
)

mock_plot <- ggplot2::ggplot() + ggplot2::geom_blank()

test_that("Testing parameters from plotDistribution", {
    expect_error(
        p <- plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Codon",
            condition_col = "Condition",
            targeted_arg = "Treatment",
            verbose = FALSE
        ),
        "The 'data' column Codon is not numeric."
    )

    expect_error(
        p <- plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Expression",
            y_axis_col = "Expression",
            condition_col = "Condition",
            targeted_arg = "Treatment",
            verbose = FALSE
        ),
        "The 'data' column Expression is not categorical."
    )

    expect_error(
        p <- plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Codon",
            targeted_arg = "Treatment", # Treatment belongs to the Condition column
            verbose = FALSE
        ),
        "The argument Treatment has not been found in the data."
    )

    expect_error(
        p <- plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Expression",
            verbose = FALSE
        ),
        "The 'data' column Expression is not categorical."
    )
})

test_that("plotDistribution returns a valid ggplot object", {
    p <- plotDistribution(
        data = mock_plot_data,
        plot = "jitter",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_label(p$mapping$x), "Codon")
    expect_equal(rlang::as_label(p$mapping$y), "Expression")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "barplot",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_label(p$mapping$x), "Codon")
    expect_equal(rlang::as_label(p$mapping$y), "Expression")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "boxplot",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_label(p$mapping$x), "Codon")
    expect_equal(rlang::as_label(p$mapping$y), "Expression")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "dot",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_label(p$mapping$x), "Codon")
    expect_equal(rlang::as_label(p$mapping$y), "Expression")

    expect_error(
        p <- plotDistribution(
            data = mock_plot_data,
            plot = "bar",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            verbose = FALSE
        ),
        "Invalid 'plot' provided: 'bar'"
    )
})

test_that("plotDistribution handles targeted_arg logic correctly", {

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "jitter",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        targeted_arg = "Treatment",
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "jitter",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        targeted_arg = c("Treatment"),
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "jitter",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        targeted_arg = list("Treatment"),
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")

    # Case B: Invalid targeted argument (should trigger the explicit stop() error)
    expect_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            targeted_arg = c("NonExistentCondition"),
            verbose = FALSE
        ),
        "The argument NonExistentCondition has not been found in the data"
    )
})

test_that("plotDistribution handels different color palettes", {
    insufficient_colors <- c("lightblue", "orange") # Only 2 colors provided
    expect_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            color_palette = insufficient_colors,
            verbose = FALSE
        ),
        "The number of colors in 'color_palette' does not correspond to the number of categories needed"
    )

    extra_colors <- c("lightblue", "orange", "lightpink", "lightgreen") # Additional colors provided
    expect_no_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            color_palette = extra_colors,
            verbose = FALSE
        )
    )

    expect_no_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "jitter",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            color_palette = list(extra_colors),
            verbose = FALSE
        )
    )

    named_colors <- c(Control = "lightblue", Treatment = "orange", Mutant = "lightpink")
    expect_no_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "barplot",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            color_palette = named_colors,
            verbose = FALSE
        )
    )

    expect_warning(
        plotDistribution(
            data = mock_plot_data,
            plot = "barplot",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Codon",
            color_palette = named_colors,
            verbose = FALSE
        ),
        "Provided 'color_palette' names do not match data categories"
    )
})

test_that("plotDistribution applies faceting parameters", {
    p <- plotDistribution(
        data = mock_plot_data,
        plot = "barplot",
        bar_position = "dodge",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        facet_col = "Tissue",
        verbose = FALSE
    )
    expect_s3_class(p$facet, "FacetWrap")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "barplot",
        bar_position = "dodge",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        facet_col = "Tissue",
        ncols = 2,
        verbose = FALSE
    )
    expect_s3_class(p$facet, "FacetWrap")

    expect_error(
        p <- plotDistribution(
            data = mock_plot_data,
            plot = "barplot",
            bar_position = "dodge",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            facet_col = "Expression",
            verbose = FALSE
        ),
        "The 'data' column Expression is not categorical."
    )
})

test_that("plotDistribution testing the display of the legend", {
    p <- plotDistribution(
        data = mock_plot_data,
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        show_legend = "none",
        verbose = FALSE
    )
    expect_equal(p$theme$legend.position, "none")

    p <- plotDistribution(
        data = mock_plot_data,
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        show_legend = "bottom",
        verbose = FALSE
    )
    expect_equal(p$theme$legend.position, "bottom")

    # Invalid string position
    expect_error(
        plotDistribution(
            data = mock_plot_data,
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            show_legend = "corner",
            verbose = FALSE
        ),
        "Specify a suitable 'show_legend' parameter"
    )

    # Incompatible logical data type
    expect_error(
        plotDistribution(
            data = mock_plot_data,
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            show_legend = FALSE,
            verbose = FALSE
        ),
        "Specify a suitable 'show_legend' parameter"
    )
})

test_that("plotDistribution checking the different bar parameters", {
    p <- plotDistribution(
        data = mock_plot_data,
        plot = "barplot",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        verbose = FALSE
    )
    pos <- p$layers[[1]]$position
    expect_s3_class(pos, "PositionDodge")

    p <- plotDistribution(
        data = mock_plot_data,
        plot = "barplot",
        bar_position = "stack",
        x_axis_col = "Codon",
        y_axis_col = "Expression",
        condition_col = "Condition",
        verbose = FALSE
    )
    pos <- p$layers[[1]]$position
    expect_s3_class(pos, "PositionStack")

    expect_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "barplot",
            bar_position = "doge",
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            verbose = FALSE
        ),
        "Specify a suitable 'bar_position' parameter."
    )

    expect_error(
        plotDistribution(
            data = mock_plot_data,
            plot = "barplot",
            bar_position = 123,
            x_axis_col = "Codon",
            y_axis_col = "Expression",
            condition_col = "Condition",
            verbose = FALSE
        ),
        "Specify a suitable 'bar_position' parameter."
    )
})

test_that("plotTargetComparison general checks", {
    expect_no_error(
        plotTargetComparison(
            target_data = mock_target,
            overall_data = mock_target,
            x_axis_col = "Codon",
            y_axis_col = "Frequency",
            verbose = FALSE
        )
    )

    expect_no_error(
        plotTargetComparison(
            target_data = mock_target,
            overall_data = mock_overall,
            x_axis_col = "Codon",
            y_axis_col = "Frequency",
            verbose = FALSE
        )
    )

    expect_error(
        plotTargetComparison(
            target_data = mock_plot_data,
            overall_data = mock_overall,
            x_axis_col = "Codon",
            y_axis_col = "Frequency",
            verbose = FALSE
        ),
        "Invalid 'y_axis_col' provided: 'Frequency'."
    )
})

test_that("plotTargetComparison throws an error for missing columns", {
    expect_error(
        plotTargetComparison(
            target_data = mock_target,
            overall_data = mock_overall,
            x_axis_col = "NonExistentColumn",
            y_axis_col = "Frequency",
            verbose = FALSE
        ),
        "Invalid 'x_axis_col' provided: 'NonExistentColumn'."
    )

    expect_error(
        plotTargetComparison(
            target_data = mock_target,
            overall_data = mock_overall,
            x_axis_col = "Codon",
            y_axis_col = "NonExistentColumn",
            verbose = FALSE
        ),
        "Invalid 'y_axis_col' provided: 'NonExistentColumn'."
    )
})

test_that("plotTargetComparison builds all 3 layers when show_difference = TRUE", {
    p <- plotTargetComparison(
        target_data = mock_target,
        overall_data = mock_overall,
        x_axis_col = "Codon",
        y_axis_col = "Frequency",
        color_palette = c("grey", "lightpink", "lightblue"),
        show_difference = TRUE,
        verbose = FALSE
    )

    expect_s3_class(p, "ggplot")

    expect_equal(rlang::as_label(p$mapping$x), "Codon")
    expect_equal(rlang::as_label(p$mapping$y), "Frequency")

    expect_length(p$layers, 3)
    expect_s3_class(p$layers[[1]]$geom, "GeomCol")
    expect_s3_class(p$layers[[2]]$geom, "GeomSegment")
    expect_s3_class(p$layers[[3]]$geom, "GeomPoint")
})

test_that("plotTargetComparison skips geom_segment when show_difference = FALSE", {
    p <- plotTargetComparison(
        target_data = mock_target,
        overall_data = mock_overall,
        x_axis_col = "Codon",
        y_axis_col = "Frequency",
        color_palette = c("grey", "lightpink", "lightblue"),
        show_difference = FALSE,
        verbose = FALSE
    )

    expect_length(p$layers, 2)
    expect_s3_class(p$layers[[1]]$geom, "GeomCol")
    expect_s3_class(p$layers[[2]]$geom, "GeomPoint")
})

test_that("Testing parameters for plotTargetComparison", {
    expect_error(
        plotTargetComparison(
            target_data = mock_target,
            overall_data = mock_overall,
            x_axis_col = "Frequency",
            y_axis_col = "Frequency",
            verbose = FALSE
        ),
        "The 'overall_data' column Frequency is not categorical."
    )

    expect_error(
        plotTargetComparison(
            target_data = mock_target,
            overall_data = mock_overall,
            x_axis_col = "Codon",
            y_axis_col = "Codon",
            verbose = FALSE
        ),
        "The 'overall_data' column Codon is not numeric"
    )
})

test_that("selectColorPalette returns default colors when palette is NULL", {
    default_palette <- selectColorPalette(palette = NULL, diff = TRUE)
    expect_equal(default_palette, c("#bababa", "#b2182b", "#fddbc7"))
})

test_that("selectColorPalette enforces correct color counts based on diff", {
    two_colors <- c("#111111", "#222222")

    expect_silent( # 2 colors given and 2 required
        res <- selectColorPalette(palette = two_colors, diff = FALSE)
    )
    expect_equal(length(res), 2)

    expect_error( # 2 colors given and 3 required
        selectColorPalette(palette = two_colors, diff = TRUE),
        "The number of colors in 'color_palette' does not correspond to the number of graphical elements"
    )
})

test_that("plotProportion executes successfully with standard layouts", {
    # BAR
    res <- plotProportion(
        data = mock_prop_data,
        plot = "bar",
        var_numerical = "Frequency",
        var_categorical = "Codon",
        verbose = FALSE
    )
    expect_type(res, "list")
    expect_s3_class(res$plot, "ggplot")
    expect_equal(rlang::as_label(res$plot$mapping$x), "Frequency")
    expect_equal(rlang::as_label(res$plot$mapping$y), "Codon")

    # DONUT
    res <- plotProportion(
        data = mock_prop_data,
        plot = "donut",
        var_numerical = "Frequency",
        var_categorical = "Codon",
        verbose = FALSE
    )
    expect_false(is.null(res))

    # RADAR
    res <- plotProportion(
        data = mock_prop_data,
        plot = "radar",
        var_numerical = "Frequency",
        var_categorical = "Codon",
        verbose = FALSE
    )
    expect_false(is.null(res))
})

test_that("Targeted checks for the radar layout of plotProportion", {
    # CASE 1
    mock_labels <- c("0%", "50%", "100%")
    p <- drawRadarPlot(
        data = mock_prop_data,
        var_color = NULL, # Tests the fallback dummy assignment
        var_categorical = "Codon",
        var_numerical = "Frequency",
        normalize = TRUE,
        zoom = FALSE,
        title = "Test Radar",
        add_titles = TRUE,
        show_legend = "none",
        global_max_val = 1.0,
        labels_rings = mock_labels,
        color_palette = c("#FF0000")
    )
    expect_s3_class(p, "ggplot")
    expect_true(length(p$layers) > 0)

    # CASE 2
    res_single <- generateProportionPlot(
        level = "radar", data = mock_prop_data,
        var_numerical = "Frequency", var_categorical = "Codon",
        var_color = "Group", color_palette = c("#111111", "#222222"),
        zoom = FALSE, show_legend = "none", order = NULL, x_limits = NULL,
        facet_col = NULL, normalize = TRUE, add_titles = TRUE, n_rings = 3
    )
    expect_type(res_single, "list")
    expect_s3_class(res_single$plot, "ggplot")
    expect_s3_class(res_single$legend, "data.frame")

    # CASE 3
    res_faceted <- generateProportionPlot(
        level = "radar", data = mock_prop_data,
        var_numerical = "Frequency", var_categorical = "Codon",
        var_color = "Group", color_palette = c("#111111", "#222222"),
        zoom = FALSE, show_legend = "none", order = NULL, x_limits = NULL,
        facet_col = "Group", # Faceting enabled!
        normalize = TRUE, add_titles = TRUE, n_rings = 3
    )
    expect_s3_class(res_faceted$plot, "patchwork")
})

test_that("plotTEscore runs successfully and returns a list with a plot and statistics table", {
    expect_error(
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            index_col = "SampleID",
            class_col = "CellType",
            target_col = "Treatment",
            add_stats = TRUE,
            verbose = FALSE
        ), # Defaults parameters for score_col and pval_col
        "Invalid 'data columns' provided: 'condition', 'p_value'."
    )

    expect_error(
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            cond_col = "condition",
            score_col = "tTE",
            pval_col = "tTE", # Use the same for the example
            index_col = "SampleID",
            class_col = "CellType",
            target_col = "Treatment",
            add_stats = TRUE,
            verbose = FALSE
        ), # Incorrect cond_col
        "Invalid 'data columns' provided: 'condition'."
    )

    expect_error(
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            cond_col = "SampleID",
            score_col = "tTE",
            pval_col = "tTE", # Use the same for the example
            index_col = "conditions",
            class_col = "CellType",
            target_col = "Treatment",
            add_stats = TRUE,
            verbose = FALSE
        ), # Incorrect index_col
        "Invalid 'index_col' provided: 'conditions'."
    )

    res <- plotTEscore(
        data = mock_te_scores,
        metadata = mock_metadata_scores,
        cond_col = "SampleID",
        score_col = "tTE",
        pval_col = "tTE", # Use the same for the example
        index_col = "SampleID",
        class_col = "CellType",
        target_col = "Treatment",
        add_stats = TRUE,
        verbose = FALSE
    )

    expect_type(res, "list")
    expect_named(res, c("plot", "stats"))

    expect_s3_class(res$plot, "ggplot")
    expect_equal(rlang::as_label(res$plot$mapping$x), "class")
    expect_equal(rlang::as_label(res$plot$mapping$y), "tTE")
    expect_equal(rlang::as_label(res$plot$mapping$fill), "target")

    expect_s3_class(res$plot$layers[[1]]$geom, "GeomViolin")
    expect_s3_class(res$plot$layers[[2]]$geom, "GeomPoint")

    colnames(mock_te_scores) <- c("SampleID", "tTEscore")
    expect_no_error({
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            cond_col = "SampleID",
            score_col = "tTEscore",
            pval_col = "tTEscore", # Use the same for the example
            index_col = "SampleID",
            class_col = "CellType",
            target_col = "Treatment",
            add_stats = TRUE,
            verbose = FALSE
        )
    })
})

test_that("plotTEscore skips stat overlays and returns NULL for stats when add_stats is FALSE", {

    expect_warning(
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            cond_col = "SampleID",
            score_col = "tTE",
            pval_col = "tTE", # Use the same for the example
            index_col = "SampleID",
            class_col = "CellType",
            target_col = "Treatment",
            add_stats = FALSE,
            verbose = FALSE
        ),
        "The 'target_col' will be ignored as 'add_stats' is not enabled."
    )

    res <- plotTEscore(
        data = mock_te_scores,
        metadata = mock_metadata_scores,
        cond_col = "SampleID",
        score_col = "tTE",
        pval_col = "tTE", # Use the same for the example
        index_col = "SampleID",
        class_col = "CellType",
        target_col = NULL, # No warning as targer_col is set to NULL
        verbose = FALSE
    )
    expect_s3_class(res$plot, "ggplot")
    expect_equal(res$plot$labels$fill, "Group")

    res <- plotTEscore(
        data = mock_te_scores,
        metadata = mock_metadata_scores,
        cond_col = "SampleID",
        score_col = "tTE",
        pval_col = "tTE", # Use the same for the example
        index_col = "SampleID",
        class_col = "CellType",
        add_stats = FALSE,
        verbose = FALSE
    )
    expect_type(res, "list")
    expect_s3_class(res$plot, "ggplot")
    expect_null(res$stats)

})

test_that("plotTEscore successfully applies different modes of color_palette", {
    my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A")
    list_colors <- list(ctrl = "#000000", drugA = "#FFFFFF", drugB = "#FF0000")
    list_colors_correct <- list(Control = "#000000", Drug_A = "#FFFFFF", Drug_B = "#FF0000")

    res <- plotTEscore(
        data = mock_te_scores,
        metadata = mock_metadata_scores,
        cond_col = "SampleID",
        score_col = "tTE",
        pval_col = "tTE", # Use the same for the example
        index_col = "SampleID",
        class_col = "CellType",
        color_palette = my_colors,
        add_stats = FALSE,
        verbose = FALSE
    )
    ggplot_built <- ggplot2::ggplot_build(res$plot)
    plot_colors <- unique(ggplot_built$data[[1]]$fill)
    expect_true(all(plot_colors %in% my_colors))

    res <- plotTEscore(
        data = mock_te_scores,
        metadata = mock_metadata_scores,
        cond_col = "SampleID",
        score_col = "tTE",
        pval_col = "tTE", # Use the same for the example
        index_col = "SampleID",
        class_col = "Treatment",
        color_palette = my_colors,
        add_stats = FALSE,
        verbose = FALSE
    )
    ggplot_built <- ggplot2::ggplot_build(res$plot)
    plot_colors <- unique(ggplot_built$data[[1]]$fill)
    expect_true(all(plot_colors %in% my_colors))

    expect_no_error(
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            cond_col = "SampleID",
            score_col = "tTE",
            pval_col = "tTE", # Use the same for the example
            index_col = "SampleID",
            class_col = "Treatment",
            color_palette = list_colors_correct,
            add_stats = FALSE,
            verbose = FALSE
        )
    )
    expect_s3_class(res$plot, "ggplot")

    expect_warning(
        res <- plotTEscore(
            data = mock_te_scores,
            metadata = mock_metadata_scores,
            cond_col = "SampleID",
            score_col = "tTE",
            pval_col = "tTE", # Use the same for the example
            index_col = "SampleID",
            class_col = "Treatment",
            color_palette = list_colors,
            add_stats = FALSE,
            verbose = FALSE
        ),
        "Palette names do not match data categories"
    )
    expect_s3_class(res$plot, "ggplot")
})

test_that("generalChecksProportionPlot passes valid inputs silently", {
    expect_error(
        generalChecksProportionPlot(
            plot = "bar",
            data = mock_prop_data,
            var_categorical = "Codon",
            var_numerical = "Frequency",
            var_color = "Group",
            facet_col = "Tissue"
        ),
        "Invalid '_col' provided: 'Tissue'."
    )

    expect_no_error(
        generalChecksProportionPlot(
            plot = "bar",
            data = mock_prop_data,
            var_categorical = "Codon",
            var_numerical = "Frequency",
            var_color = "Group",
            facet_col = "AminoAcid"
        )
    )

    expect_no_error(
        generalChecksProportionPlot(
            plot = "donut",
            data = mock_prop_data,
            var_categorical = "Codon",
            var_numerical = "Frequency",
            var_color = "Codon", # Using the categorical feature as its own color map
            facet_col = NULL
        )
    )

    expect_no_error(
        generalChecksProportionPlot(
            plot = "radar",
            data = mock_prop_data,
            var_categorical = "Codon",
            var_numerical = "Frequency",
            var_color = "Group",
            facet_col = NULL
        )
    )
})

test_that("plotCorrelation runs successfully with basic settings", {
    p <- plotCorrelation(
        data = mock_corr_data,
        plot = "MeanCodonUsage",
        x_axis_col = "Usage_Control",
        y_axis_col = "Usage_Treatment",
        condition_col = "Group",
        add_titles = FALSE,
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")
    expect_equal(rlang::as_label(p$mapping$x), "Usage_Control")
    expect_equal(rlang::as_label(p$mapping$y), "Usage_Treatment")
    expect_equal(rlang::as_label(p$mapping$colour), "Group")
    expect_s3_class(p$layers[[1]]$geom, "GeomPoint")

    p <- plotCorrelation(
        data = mock_corr_data,
        x_axis_col = "Usage_Control",
        y_axis_col = "Usage_Treatment",
        condition_col = "Group",
        extra_val = 0.854, # Triggers annotation block
        add_titles = FALSE,
        verbose = FALSE
    )
    ann_layer <- p$layers[[length(p$layers)]]
    expect_s3_class(ann_layer$geom, "GeomText")
    expect_equal(ann_layer$aes_params$label, "rho == 0.854")
})

test_that("plotCorrelation safely uses the additional parameters", {
    p <- plotCorrelation(
        data = mock_corr_data,
        x_axis_col = "Usage_Control",
        y_axis_col = "Usage_Treatment",
        condition_col = "Group",
        label_col = "Codon", # Triggers ggrepel injection
        add_titles = FALSE,
        verbose = FALSE
    )

    geoms <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_true("GeomTextRepel" %in% geoms)

    expect_silent({
        p <- plotCorrelation(
            data = mock_corr_data,
            x_axis_col = "Usage_Control",
            y_axis_col = "Usage_Treatment",
            condition_col = "Group",
            targeted_arg = "Phe", # Triggers targeted subset routing
            add_titles = FALSE,
            verbose = FALSE
        )
    })
    expect_s3_class(p, "ggplot")
})

test_that("plotPermutation builds a multi-layered empirical null visualization", {
    p <- plotPermutation(
        permut_data = mock_permut_data,
        sig_data = mock_sig_data,
        color_palette = c(Treatment_A = "#FF0000"),
        add_titles = TRUE,
        verbose = FALSE
    )
    expect_s3_class(p, "ggplot")
    geoms <- sapply(p$layers, function(l) class(l$geom)[1])
    expect_true("GeomBar" %in% geoms)
    expect_true("GeomVline" %in% geoms)
    expect_true("GeomText" %in% geoms)
})

test_that("plotPermutation successfully handles flat, zero-variance permutation vectors", {
    flat_permut_data <- data.frame(
        codon = rep("AUG", 10),
        freq = rep(0.25, 10),
        stringsAsFactors = FALSE
    )

    expect_silent({
        p <- plotPermutation(
            permut_data = flat_permut_data,
            sig_data = mock_sig_data[1, ],
            add_titles = FALSE,
            verbose = FALSE
        )
    })

    expect_s3_class(p, "ggplot")
})

test_that("plotPermutation isolates codons inside individual free-scaled facets", {
    p <- plotPermutation(
        permut_data = mock_permut_data,
        sig_data = mock_sig_data,
        add_titles = FALSE,
        verbose = FALSE
    )

    expect_s3_class(p$facet, "FacetWrap")
    expect_true(p$facet$params$free$y)  # Confirms independent y-axes per codon
    expect_false(p$facet$params$free$x) # Confirms shared x-axis baseline
})

test_that("savePlot successfully writes files inside isolated environments", {
    virtual_dir <- withr::local_tempdir()
    expect_silent({
        savePlot(
            plot = mock_plot,
            save_format = "Pdf", # Mix case to check tolower()
            out_name = "test_run",
            out_directory = virtual_dir,
            verbose = FALSE
        )
    })
    expect_true(file.exists(file.path(virtual_dir, "test_run.pdf")))

    savePlot(
        plot = mock_plot,
        save_format = "png",
        out_name = "test_run",
        out_directory = virtual_dir,
        verbose = FALSE
    )
    expect_true(file.exists(file.path(virtual_dir, "test_run.png")))

    virtual_dir <- withr::local_tempdir()
    expect_warning(
        savePlot(
            plot = mock_plot,
            save_format = "tiff", # Unsupported format
            out_name = "test_corrupt",
            out_directory = virtual_dir,
            verbose = FALSE
        ),
        "is not recognized"
    )
    expect_false(file.exists(file.path(virtual_dir, "test_corrupt.tiff")))
})

test_that("getOutputName computes correct fallback paths and formats strings", {
    res_path1 <- getOutputName(
        action = "plot",
        out_name = "my_visualization",
        out_directory = "/mock/dir",
        save_format = "png",
        verbose = FALSE
    )
    expect_equal(basename(res_path1), "my_visualization.png")

    res_path2 <- getOutputName(
        action = "plot",
        out_name = "my_visualization.png",
        out_directory = "/mock/dir",
        save_format = "png",
        verbose = FALSE
    )
    expect_equal(basename(res_path2), "my_visualization.png")

    expect_message(
        res_plot <- getOutputName(
            action = "plot",
            out_name = NULL, # Triggers fallback to "distribution_plot"
            out_directory = "/mock/dir",
            save_format = "pdf",
            verbose = TRUE
        ),
        "A standard name will be used"
    )
    expect_equal(basename(res_plot), "distribution_plot.pdf")

    res_matrix <- getOutputName(
        action = "matrix",
        out_name = NULL, # Triggers fallback to "tRNA_expression_matrix"
        out_directory = "/mock/dir",
        save_format = "pdf",
        verbose = FALSE
    )
    expect_equal(basename(res_matrix), "tRNA_expression_matrix.pdf")
})

test_that("addSignificanceDist accurately maps and appends significance bars to distribution plots", {
    base_canvas <- ggplot2::ggplot(
        data = mock_plot_data,
        mapping = ggplot2::aes(x = Condition, y = Expression)
    ) + ggplot2::geom_boxplot()

    res <- addSignificanceDist(
        data = mock_plot_data,
        plot = base_canvas,
        x_axis_col = "Condition",
        y_axis_col = "Expression",
        target = "Tissue",
        facet_col = NULL
    )

    expect_type(res, "list")
    expect_named(res, c("sig_table", "plot"))
    expect_s3_class(res$plot, "ggplot")

    if (!is.null(res$sig_table)) {
        expect_true("group1" %in% colnames(res$sig_table))
        expect_true("group2" %in% colnames(res$sig_table))
        expect_equal(res$sig_table$group1[1], "Kidney")
        expect_equal(res$sig_table$group2[1], "Liver")
    }
})

test_that("addSignificanceDist handles target groups with more than 2 levels by returning the base plot", {

    base_canvas <- ggplot2::ggplot(
        data = mock_plot_data,
        mapping = ggplot2::aes(x = Condition, y = Expression)
    ) + ggplot2::geom_boxplot()

    res <- addSignificanceDist(
        data = mock_plot_data,
        plot = base_canvas,
        x_axis_col = "Condition",
        y_axis_col = "Expression",
        target = "Codon", # 3 levels -> fails
        facet_col = NULL
    )

    expect_null(res$sig_table)
    expect_equal(length(res$plot$layers), length(base_canvas$layers))
})
