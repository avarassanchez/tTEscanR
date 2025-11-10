TransformCounts <- function(data) {
  data_long <- DataToLongFormat(data = data, normalize = FALSE, rownames_to_column = "features",
                                names_to = "conditions",values_to = "counts")

  data_long <- dplyr::filter(data_long, .data$counts > 0)
  data_long <- tidyr::uncount(data_long, weights = .data$counts)

  return(data_long)
}

CutoffMatrix <- function(data, cutoff) {
  data <- dplyr::slice_sample(data, n = cutoff)

  data <- data %>% dplyr::count(.data$features, .data$conditions) %>%
    tidyr::pivot_wider(names_from = .data$conditions, values_from = .data$n, values_fill = 0)

  matrix <- data %>% tibble::column_to_rownames("features") %>% as.matrix()

  return(matrix)
}

FilteringCutoffs <- function(data, cutoff) {
  data_above <- CutoffMatrix(data = data, cutoff = cutoff)

  object_above <- suppressMessages(Create_tTEscanR_Object(counts = data_above, assay = "tRNA", verbose = FALSE))
  object_above <- suppressMessages(ComputeAnticodonUsage(object = object_above, verbose = FALSE))
  object_above <- suppressMessages(ComputeAAUsage(object = object_above, level = "supply", verbose = FALSE))

  total_anticodon <- rowSums(object_above@assays$AnticodonUsage)
  total_supply <- rowSums(object_above@assays$AASupply)

  # Normalize
  total_anticodon <- total_anticodon / sum(total_anticodon)
  total_supply <- total_supply / sum(total_supply)

  return(list(anticodon = total_anticodon, supply = total_supply))
}

ComputeCorrelations <- function(data, ref_anticodon, ref_supply, cutoffs) {

  cor_results <- data.frame(cutoff = cutoffs, anticodon_spearman = NA_real_,
                            supply_spearman = NA_real_, total_anticodon = NA_real_,
                            total_supply = NA_real_)

  for (i in seq_along(cutoffs)) {
    cut <- cutoffs[i]
    obj <- FilteringCutoffs(data = data, cutoff = cut)

    if (length(obj$anticodon) == 0 || length(obj$supply) == 0)
      next

    anticodon_vec <- obj$anticodon[names(ref_anticodon)]
    anticodon_vec[is.na(anticodon_vec)] <- 0

    supply_vec <- obj$supply[names(ref_supply)]
    supply_vec[is.na(supply_vec)] <- 0

    cor_results$anticodon_spearman[i] <- suppressWarnings(stats::cor(anticodon_vec, ref_anticodon, method = "spearman"))
    cor_results$supply_spearman[i] <- suppressWarnings(stats::cor(supply_vec, ref_supply, method = "spearman"))

    cor_results$total_anticodon[i] <- sum(obj$anticodon)
    cor_results$total_supply[i] <- sum(obj$supply)
  }
  return(cor_results)
}

Selection_Cutoff <- function(data, slope_threshold = 0.001) {

  cor_long <- data %>%
    tidyr::pivot_longer(cols = c("anticodon_spearman", "supply_spearman"), names_to = "type", values_to = "spearman_corr")

  cor_stability <- cor_long %>%
    dplyr::group_by(.data$type) %>%
    dplyr::arrange(.data$cutoff, .by_group = TRUE) %>%
    dplyr::mutate(delta_corr = abs(.data$spearman_corr - dplyr::lag(.data$spearman_corr)),
                  slope_ratio = .data$delta_corr / abs(.data$cutoff - dplyr::lag(.data$cutoff)))

  stable_points <- cor_stability %>%
    dplyr::filter(is.finite(.data$slope_ratio), .data$slope_ratio < slope_threshold)

  if (nrow(stable_points) == 0) {
    message("No stable points found below slope threshold = ", slope_threshold)
    return(NULL)
  }

  optimal_cutoff <- stable_points %>%
    dplyr::group_by(.data$type) %>%
    dplyr::summarise(optimal_cutoff = min(.data$cutoff), .groups = "drop")

  return(optimal_cutoff)
}
