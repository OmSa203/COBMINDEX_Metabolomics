#======================== ENVIRONMENT & LIBRARIES ========================
#----------------Cleaning the environment
rm(list = ls()) 
#----------------Clear all plots
graphics.off()
#----------------Clear console
cat("\014") 

suppressPackageStartupMessages({
  library(tidyverse) # For data manipulation (dplyr) and plotting (ggplot2)
  library(vegan)     # For diversity analysis and ordination methods
  library(scales)
  library(ggrepel)   # For non-overlapping labels in ggplot
  library(ggpubr)    # For adding p-values to plots
  library(rstatix)   # For statistical test annotations
})


#======================== HELPER FUNCTION ========================

# Helper function for creating a directory and returning a full file path.
check_dir_and_save <- function(directory, filename) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    message("Created directory: ", directory)
  }
  file.path(directory, filename)
}

#======================== KRUSKAL-WALLIS FEATURE OUTPUTS (NEW) ========================

#' Create a boxplot for a single feature from a Kruskal-Wallis test.
plot_significant_kw_feature <- function(feature_data, kw_results_row, title, grouping_column) {
  # Generate pairwise comparisons using Dunn's test for p-value annotation
  stat_test <- feature_data %>%
    dunn_test(as.formula(paste("intensity ~", grouping_column)), p.adjust.method = "bh") %>%
    add_xy_position(x = grouping_column)
  
  # Overall KW p-value for the subtitle
  kw_p_val <- scales::pvalue(kw_results_row$kw_p_adj_bh)
  
  # Create the plot
  p <- ggplot(feature_data, aes_string(x = grouping_column, y = "intensity")) +
    geom_boxplot(aes_string(fill = grouping_column), outlier.shape = NA, alpha = 0.6, show.legend = FALSE) +
    geom_jitter(aes_string(color = grouping_column), width = 0.15, height = 0, size = 2, show.legend = FALSE) +
    stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    scale_y_continuous(labels = scales::scientific_format(digits=2)) +
    labs(
      title = title,
      subtitle = paste("Kruskal-Wallis adjusted p-value:", kw_p_val),
      y = "Normalized Intensity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1) # Rotate labels for readability
    )
  return(p)
}


#' Calculate summary statistics for a single feature and save to CSV.
calculate_and_save_kw_stats <- function(feature_data, output_path, grouping_column) {
  stats_df <- feature_data %>%
    group_by(!!sym(grouping_column)) %>%
    summarise(
      n_samples = n(),
      n_zeros = sum(intensity == 0, na.rm = TRUE),
      zero_percent = mean(intensity == 0, na.rm = TRUE) * 100,
      mean_intensity = mean(intensity, na.rm = TRUE),
      median_intensity = median(intensity, na.rm = TRUE),
      sd_intensity = sd(intensity, na.rm = TRUE),
      cv_percent = (sd(intensity, na.rm = TRUE) / mean(intensity, na.rm = TRUE)) * 100,
      .groups = 'drop'
    )
  write.csv(stats_df, file = output_path, row.names = FALSE)
}


#' Main function to generate outputs for all features from a KW test.
generate_kw_feature_outputs <- function(kw_results, feature_table, metadata, config, mode, result_dir) {
  
  features_to_process <- kw_results %>%
    distinct(m.z, retention_time, .keep_all = TRUE) %>%
    arrange(kw_p_adj_bh)
  
  if (nrow(features_to_process) == 0) {
    cat("No features found in Kruskal-Wallis results to generate individual outputs.\n")
    return()
  }
  
  cat(paste("Found", nrow(features_to_process), "features from Kruskal-Wallis test. Generating outputs for all...\n"))
  
  output_dir <- file.path(result_dir, config$folder_name, "All_Feature_Outputs")
  grouping_col <- config$grouping_column %||% "Group"
  
  clean_sample_name <- function(x, ion_mode) gsub(paste0("_", ion_mode, "\\.mzML$"), "", x)
  
  for (i in 1:nrow(features_to_process)) {
    feature_info <- features_to_process[i, ]
    mz <- feature_info$m.z
    rt <- feature_info$retention_time
    
    # Prepare data for this one feature
    feature_data_long <- feature_table %>%
      filter(row.m.z == mz, row.retention.time == rt) %>%
      select(-row.ID) %>%
      pivot_longer(
        cols = matches("S\\d+_.*\\.mzML"),
        names_to = "filename",
        values_to = "intensity"
      ) %>%
      mutate(filename = clean_sample_name(filename, ion_mode = mode)) %>%
      inner_join(metadata, by = "filename")
    
    # Apply filters based on the analysis config
    if (grouping_col == "Group") {
      feature_data_long <- feature_data_long %>%
        filter(Timepoint == config$timepoint & !!sym(grouping_col) %in% config$groups)
    } else {
      feature_data_long <- feature_data_long %>%
        filter(!!sym(grouping_col) %in% config$groups)
    }
    
    if (nrow(feature_data_long) == 0) next
    
    plot_title <- sprintf("Feature at m/z: %.4f, RT: %.2f", mz, rt)
    feature_plot <- plot_significant_kw_feature(feature_data_long, feature_info, plot_title, grouping_col)
    
    p_val_str <- format(feature_info$kw_p_adj_bh, scientific = TRUE, digits = 2)
    base_filename <- sprintf("feature_%.4f_%.2f_pval_%s", mz, rt, p_val_str)
    
    plot_path <- check_dir_and_save(output_dir, paste0(base_filename, ".png"))
    ggsave(plot_path, feature_plot, width = 8, height = 7, dpi = 300)
    
    stats_path <- check_dir_and_save(output_dir, paste0(base_filename, "_statistics.csv"))
    calculate_and_save_kw_stats(feature_data_long, stats_path, grouping_col)
  }
  cat("Finished generating individual Kruskal-Wallis feature outputs.\n")
}


#======================== WILCOXON FEATURE OUTPUTS ========================
# This section is for the Wilcoxon test outputs and remains as previously defined.
#' Create a boxplot for a single significant feature from a Wilcoxon test.
plot_significant_wilcoxon_feature <- function(feature_data, result_row, title) {
  p_val_adj <- scales::pvalue(result_row$p_adjusted_BH)
  feature_data <- feature_data %>%
    mutate(Group_Timepoint = paste(Group, "T", Timepoint, sep=""))
  p <- ggplot(feature_data, aes(x = Group_Timepoint, y = intensity)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.6) +
    geom_jitter(aes(color = Group), width = 0.15, height = 0, size = 2.5) +
    scale_y_continuous(labels = scales::scientific_format(digits=2)) +
    labs(
      title = title,
      subtitle = paste("Wilcoxon Test, Adjusted p-value:", p_val_adj),
      x = "Comparison Groups",
      y = "Normalized Intensity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  return(p)
}

#' Calculate summary statistics for a single feature from a Wilcoxon test.
calculate_and_save_wilcoxon_stats <- function(feature_data, output_path) {
  stats_df <- feature_data %>%
    group_by(Group, Timepoint) %>%
    summarise(
      n_samples = n(),
      n_zeros = sum(intensity == 0, na.rm = TRUE),
      zero_percent = mean(intensity == 0, na.rm = TRUE) * 100,
      mean_intensity = mean(intensity, na.rm = TRUE),
      median_intensity = median(intensity, na.rm = TRUE),
      sd_intensity = sd(intensity, na.rm = TRUE),
      cv_percent = (sd(intensity, na.rm = TRUE) / mean(intensity, na.rm = TRUE)) * 100,
      .groups = 'drop'
    )
  write.csv(stats_df, file = output_path, row.names = FALSE)
}

#' Main function to generate outputs for ALL features from a Wilcoxon test.
generate_wilcoxon_feature_outputs <- function(wilcoxon_results, feature_table, metadata, config, mode, result_dir) {
  features_to_process <- wilcoxon_results %>%
    arrange(p_value)
  if (nrow(features_to_process) == 0) {
    cat("No features found in Wilcoxon results to generate individual outputs.\n")
    return()
  }
  cat(paste("Found", nrow(features_to_process), "features from Wilcoxon test. Generating outputs for all...\n"))
  output_dir <- file.path(result_dir, config$folder_name, "All_Feature_Outputs")
  clean_sample_name <- function(x, ion_mode) gsub(paste0("_", ion_mode, "\\.mzML$"), "", x)
  for (i in 1:nrow(features_to_process)) {
    feature_info <- features_to_process[i, ]
    mz <- feature_info$row.m.z
    rt <- feature_info$row.retention.time
    feature_data_long <- feature_table %>%
      filter(row.m.z == mz, row.retention.time == rt) %>%
      select(-row.ID) %>%
      pivot_longer(
        cols = matches("S\\d+_.*\\.mzML"),
        names_to = "filename",
        values_to = "intensity"
      ) %>%
      mutate(filename = clean_sample_name(filename, ion_mode = mode)) %>%
      inner_join(
        metadata %>%
          mutate(filename = clean_sample_name(filename, ion_mode = mode)),
        by = "filename"
      ) %>%
      filter(
        (Group == config$group1 & Timepoint == config$timepoint1) |
          (Group == config$group2 & Timepoint == config$timepoint2)
      )
    if (nrow(feature_data_long) == 0) next
    plot_title <- sprintf("Feature at m/z: %.4f, RT: %.2f", mz, rt)
    feature_plot <- plot_significant_wilcoxon_feature(feature_data_long, feature_info, plot_title)
    p_val_str <- format(feature_info$p_adjusted_BH, scientific = TRUE, digits = 2)
    base_filename <- sprintf("feature_%.4f_%.2f_pval_%s", mz, rt, p_val_str)
    plot_path <- check_dir_and_save(output_dir, paste0(base_filename, ".png"))
    ggsave(plot_path, feature_plot, width = 7, height = 6, dpi = 300)
    stats_path <- check_dir_and_save(output_dir, paste0(base_filename, "_statistics.csv"))
    calculate_and_save_wilcoxon_stats(feature_data_long, stats_path)
  }
  cat("Finished generating individual Wilcoxon feature outputs.\n")
}


#======================== KRUSKAL-WALLIS CORE & WRAPPERS (MODIFIED) ========================

perform_kruskal_posthoc_analysis <- function(feature_table, metadata, config, mode) {
  
  grouping_col <- config$grouping_column %||% "Group"
  
  clean_sample_name <- function(x, ion_mode) gsub(paste0("_", ion_mode, "\\.mzML$"), "", x)
  
  # Clean metadata filenames and filter for valid samples BEFORE the join
  metadata_cleaned <- metadata %>%
    mutate(filename = clean_sample_name(filename, ion_mode = mode)) %>%
    filter(Sample_Type == "Sample")
  
  prepared_data <- feature_table %>%
    select(-row.ID) %>%
    pivot_longer(
      cols = matches("S\\d+_.*\\.mzML"),
      names_to = "filename",
      values_to = "intensity"
    ) %>%
    mutate(
      filename = clean_sample_name(filename, ion_mode = mode),
      m.z = row.m.z,
      retention_time = row.retention.time
    ) %>%
    inner_join(metadata_cleaned, by = "filename") # Join with the cleaned and filtered metadata
  
  # Apply filters based on the analysis config
  if (grouping_col == "Group") {
    prepared_data <- prepared_data %>%
      filter(Timepoint == config$timepoint & !!sym(grouping_col) %in% config$groups)
  } else {
    prepared_data <- prepared_data %>%
      filter(!!sym(grouping_col) %in% config$groups)
  }
  
  if (nrow(prepared_data) == 0) {
    warning("No data found for the specified filters. Returning an empty data frame.")
    return(tibble())
  }
  
  analysis_results <- prepared_data %>%
    group_by(m.z, retention_time) %>%
    nest() %>%
    mutate(
      kw_test_result = map(data, ~ tryCatch(kruskal.test(as.formula(paste("intensity ~", grouping_col)), data = .x), error = function(e) NULL)),
      statistic = map_dbl(kw_test_result, ~ if(is.null(.x)) NA_real_ else .x$statistic),
      p.value   = map_dbl(kw_test_result, ~ if(is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    select(-kw_test_result) %>%
    ungroup() %>%
    filter(!is.na(p.value)) %>%
    mutate(kw_p_adj_bh = p.adjust(p.value, method = "BH"))
  
  final_results <- analysis_results %>%
    mutate(
      posthoc_test = map(data, ~ {
        pw_test <- pairwise.wilcox.test(.x$intensity, .x[[grouping_col]], p.adjust.method = "BH")
        if (!is.null(pw_test$p.value)) {
          pw_test$p.value %>%
            as.data.frame() %>%
            rownames_to_column(var = "group1") %>%
            pivot_longer(cols = -group1, names_to = "group2", values_to = "posthoc_p_adj") %>%
            filter(!is.na(posthoc_p_adj))
        } else { NULL }
      })
    ) %>%
    select(-data) %>%
    unnest(posthoc_test, keep_empty = TRUE) %>%
    select(
      m.z, retention_time,
      kw_statistic = statistic, kw_p_value = p.value, kw_p_adj_bh,
      group1, group2, posthoc_p_adj
    ) %>%
    arrange(kw_p_adj_bh, m.z, retention_time)
  
  return(final_results)
}

run_and_save_kw_analysis <- function(feature_table, metadata, result_dir, analysis_config, mode) {
  cat("--------------------------------------------------\n")
  cat("Running Kruskal-Wallis Analysis:", analysis_config$title, "\n")
  
  results_df <- perform_kruskal_posthoc_analysis(
    feature_table = feature_table,
    metadata = metadata,
    config = analysis_config,
    mode = mode
  )
  
  if (nrow(results_df) > 0) {
    analysis_specific_dir <- file.path(result_dir, analysis_config$folder_name)
    results_path <- check_dir_and_save(analysis_specific_dir, paste0("KW_results_", mode, analysis_config$file_suffix, ".csv"))
    write.csv(results_df, file = results_path, row.names = FALSE)
    cat("Successfully saved KW results for:", analysis_config$title, "\n")
    
    # Generate and save plots/stats for all features
    generate_kw_feature_outputs(
      kw_results = results_df,
      feature_table = feature_table,
      metadata = metadata,
      config = analysis_config,
      mode = mode,
      result_dir = result_dir
    )
  } else {
    cat("Skipping file save for", analysis_config$title, "as no results were generated.\n")
  }
}

run_kw_analyses_for_mode <- function(mode, data_dir, result_dir, analysis_list) {
  cat("\n==================================================\n")
  cat("STARTING KRUSKAL-WALLIS ANALYSIS FOR MODE:", toupper(mode), "\n")
  cat("==================================================\n")
  
  mode_string <- ifelse(mode %in% c("Positive", "POS", "pos"), "Positive", "Negative")
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", mode_string, ".csv"))
  metadata_file <- file.path(data_dir, paste0(mode_string, "_metadata_w_subjects.csv"))
  
  if (!file.exists(feature_table_file) || !file.exists(metadata_file)) {
    warning("Data files not found for KW mode '", mode, "'. Skipping this mode.", immediate. = TRUE)
    return()
  }
  
  feature_table <- read.csv(feature_table_file)
  metadata <- read.csv(metadata_file)
  
  # Create the combined Group_Timepoint column for the new analysis
  metadata <- metadata %>%
    mutate(Group_Timepoint = case_when(
      Group == "Healthy" ~ "Healthy",
      TRUE ~ paste0(Group, Timepoint)
    ))
  
  for (analysis_config in analysis_list) {
    run_and_save_kw_analysis(
      feature_table = feature_table,
      metadata = metadata,
      result_dir = result_dir,
      analysis_config = analysis_config,
      mode = mode
    )
  }
}


#======================== WILCOXON TEST CORE & WRAPPERS ========================
# This section remains as previously defined.
perform_wilcoxon_analysis <- function(feature_table, metadata, config) {
  group1_samples <- metadata %>%
    filter(Group == config$group1, Timepoint == config$timepoint1) %>%
    arrange(Subject)
  group2_samples <- metadata %>%
    filter(Group == config$group2, Timepoint == config$timepoint2) %>%
    arrange(Subject)
  if (config$paired) {
    paired_subjects <- intersect(group1_samples$Subject, group2_samples$Subject)
    if (length(paired_subjects) == 0) {
      warning("No subjects found with both conditions for paired analysis: ", config$title)
      return(NULL)
    }
    sample_pairs <- data.frame(Subject = paired_subjects)
    sample_pairs$Sample1 <- sapply(paired_subjects, function(s) group1_samples$filename[group1_samples$Subject == s][1])
    sample_pairs$Sample2 <- sapply(paired_subjects, function(s) group2_samples$filename[group2_samples$Subject == s][1])
    sample1_cols <- sample_pairs$Sample1
    sample2_cols <- sample_pairs$Sample2
  } else {
    sample1_cols <- group1_samples$filename
    sample2_cols <- group2_samples$filename
  }
  if (!all(c(sample1_cols, sample2_cols) %in% colnames(feature_table))) {
    warning("One or more sample columns not found in feature table for analysis: ", config$title)
    return(NULL)
  }
  group1_data <- feature_table[, sample1_cols, drop = FALSE]
  group2_data <- feature_table[, sample2_cols, drop = FALSE]
  results_list <- lapply(1:nrow(feature_table), function(i) {
    g1_vals <- as.numeric(group1_data[i, ])
    g2_vals <- as.numeric(group2_data[i, ])
    if(all(is.na(g1_vals)) || all(is.na(g2_vals))) return(NULL)
    test_res <- tryCatch(
      wilcox.test(g1_vals, g2_vals, paired = config$paired, exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(test_res)) return(NULL)
    data.frame(
      p_value = test_res$p.value,
      statistic = test_res$statistic
    )
  })
  results_df <- bind_rows(results_list)
  if(nrow(results_df) == 0) return(NULL)
  results_df <- cbind(feature_table[, c("row.ID", "row.m.z", "row.retention.time")], results_df)
  results_df$p_adjusted_BH <- p.adjust(results_df$p_value, method = "BH")
  results_df$mean_group1 <- rowMeans(group1_data, na.rm = TRUE)
  results_df$mean_group2 <- rowMeans(group2_data, na.rm = TRUE)
  results_df$fold_change <- results_df$mean_group2 / results_df$mean_group1
  results_df$log2_fold_change <- log2(results_df$fold_change)
  results_df$significance <- ifelse(results_df$p_adjusted_BH < 0.05, "Significant (FDR < 0.05)", "Not Significant")
  results_df <- results_df %>%
    arrange(p_value) %>%
    mutate(label = if_else(row_number() <= 10, as.character(row.ID), NA_character_))
  comparison_name <- paste0(config$group1, config$timepoint1, "_vs_", config$group2, config$timepoint2)
  volcano_plot <- ggplot(results_df, aes(x = log2_fold_change, y = -log10(p_value))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = c("Significant (FDR < 0.05)" = "red", "Not Significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste("Volcano Plot:", comparison_name, ifelse(config$paired, "(Paired)", "(Unpaired)")),
      x = paste("Log2 Fold Change (", config$group2, config$timepoint2, "/", config$group1, config$timepoint1, ")", sep=""),
      y = "-Log10(P-value)"
    ) +
    theme_minimal()
  p_value_hist <- ggplot(results_df, aes(x = p_value)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    labs(
      title = paste("P-value Distribution:", comparison_name),
      x = "P-value", y = "Frequency"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5))
  return(list(
    results_df = arrange(results_df, p_value),
    volcano_plot = volcano_plot,
    p_value_hist = p_value_hist
  ))
}

run_and_save_wilcoxon_analysis <- function(feature_table, metadata, result_dir, analysis_config, mode) {
  cat("--------------------------------------------------\n")
  cat("Running Wilcoxon Analysis:", analysis_config$title, "\n")
  results <- perform_wilcoxon_analysis(
    feature_table = feature_table,
    metadata = metadata,
    config = analysis_config
  )
  if (!is.null(results)) {
    analysis_specific_dir <- file.path(result_dir, analysis_config$folder_name)
    results_path <- check_dir_and_save(analysis_specific_dir, paste0("Wilcoxon_results_", mode, analysis_config$file_suffix, ".csv"))
    write.csv(results$results_df, file = results_path, row.names = FALSE)
    volcano_path <- check_dir_and_save(analysis_specific_dir, paste0("Volcano_plot_", mode, analysis_config$file_suffix, ".png"))
    ggsave(volcano_path, results$volcano_plot, width = 10, height = 8, dpi = 300)
    pval_hist_path <- check_dir_and_save(analysis_specific_dir, paste0("Pval_hist_", mode, analysis_config$file_suffix, ".png"))
    ggsave(pval_hist_path, results$p_value_hist, width = 8, height = 6, dpi = 300)
    cat("Successfully saved Wilcoxon results for:", analysis_config$title, "\n")
    generate_wilcoxon_feature_outputs(
      wilcoxon_results = results$results_df,
      feature_table = feature_table,
      metadata = metadata,
      config = analysis_config,
      mode = mode,
      result_dir = result_dir
    )
  } else {
    cat("Skipping Wilcoxon analysis for", analysis_config$title, "due to data issues.\n")
  }
}

run_wilcoxon_analyses_for_mode <- function(mode, data_dir, result_dir, analysis_list) {
  cat("\n==================================================\n")
  cat("STARTING WILCOXON ANALYSIS FOR MODE:", toupper(mode), "\n")
  cat("==================================================\n")
  mode_string <- ifelse(mode %in% c("Positive", "POS", "pos"), "Positive", "Negative")
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", mode_string, ".csv"))
  metadata_file <- file.path(data_dir, paste0(mode_string, "_metadata_w_subjects.csv"))
  if (!file.exists(feature_table_file) || !file.exists(metadata_file)) {
    warning("Data files not found for Wilcoxon mode '", mode, "'. Skipping this mode.", immediate. = TRUE)
    return()
  }
  feature_table <- read.csv(feature_table_file)
  metadata <- read.csv(metadata_file)
  for (analysis_config in analysis_list) {
    run_and_save_wilcoxon_analysis(
      feature_table = feature_table,
      metadata = metadata,
      result_dir = result_dir,
      analysis_config = analysis_config,
      mode = mode
    )
  }
}


#======================== INTERACTION (CHANGE-SCORE) CORE & WRAPPERS ========================

perform_interaction_analysis <- function(feature_table, metadata, config) {
  get_change_scores <- function(group_name) {
    t1_samples <- metadata %>% filter(Group == group_name, Timepoint == config$timepoint1)
    t2_samples <- metadata %>% filter(Group == group_name, Timepoint == config$timepoint2)
    paired_subjects <- intersect(t1_samples$Subject, t2_samples$Subject)
    if (length(paired_subjects) == 0) {
      warning("No paired subjects found for group: ", group_name)
      return(NULL)
    }
    t1_paired_files <- sapply(paired_subjects, function(s) t1_samples$filename[t1_samples$Subject == s][1])
    t2_paired_files <- sapply(paired_subjects, function(s) t2_samples$filename[t2_samples$Subject == s][1])
    t1_data <- feature_table[, t1_paired_files, drop = FALSE]
    t2_data <- feature_table[, t2_paired_files, drop = FALSE]
    change_scores <- t2_data - t1_data
    return(change_scores)
  }
  group1_changes <- get_change_scores(config$group1)
  group2_changes <- get_change_scores(config$group2)
  if (is.null(group1_changes) || is.null(group2_changes)) {
    warning("Could not calculate change scores for one or both groups for analysis: ", config$title)
    return(NULL)
  }
  results_list <- lapply(1:nrow(feature_table), function(i) {
    g1_diffs <- as.numeric(group1_changes[i, ])
    g2_diffs <- as.numeric(group2_changes[i, ])
    if(all(is.na(g1_diffs)) || all(is.na(g2_diffs))) return(NULL)
    test_res <- tryCatch(
      wilcox.test(g1_diffs, g2_diffs, paired = FALSE, exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(test_res)) return(NULL)
    data.frame(p_value = test_res$p.value, statistic = test_res$statistic)
  })
  results_df <- bind_rows(results_list)
  if(nrow(results_df) == 0) return(NULL)
  results_df <- cbind(feature_table[, c("row.ID", "row.m.z", "row.retention.time")], results_df)
  results_df$p_adjusted_BH <- p.adjust(results_df$p_value, method = "BH")
  results_df$mean_change_group1 <- rowMeans(group1_changes, na.rm = TRUE)
  results_df$mean_change_group2 <- rowMeans(group2_changes, na.rm = TRUE)
  epsilon <- 1e-9
  ratio <- (results_df$mean_change_group2 + epsilon) / (results_df$mean_change_group1 + epsilon)
  results_df$log2_fold_change_of_change <- log2(abs(ratio)) * sign(ratio)
  results_df$significance <- ifelse(results_df$p_adjusted_BH < 0.05, "Significant (FDR < 0.05)", "Not Significant")
  results_df <- results_df %>%
    arrange(p_value) %>%
    mutate(label = if_else(row_number() <= 10, as.character(row.ID), NA_character_))
  volcano_plot <- ggplot(results_df, aes(x = log2_fold_change_of_change, y = -log10(p_value))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = c("Significant (FDR < 0.05)" = "purple", "Not Significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_text_repel(aes(label = label), max.overlaps = 15) +
    labs(
      title = paste("Interaction Volcano Plot:", config$title),
      x = "Log2(FC)",
      y = "-Log10(P-value)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom")
  return(list(results_df = arrange(results_df, p_value), volcano_plot = volcano_plot))
}

run_and_save_interaction_analysis <- function(feature_table, metadata, result_dir, analysis_config, mode) {
  cat("--------------------------------------------------\n")
  cat("Running Interaction Analysis:", analysis_config$title, "\n")
  results <- perform_interaction_analysis(
    feature_table = feature_table, metadata = metadata, config = analysis_config
  )
  if (!is.null(results)) {
    analysis_specific_dir <- file.path(result_dir, analysis_config$folder_name)
    results_path <- check_dir_and_save(analysis_specific_dir, paste0("Interaction_results_", mode, analysis_config$file_suffix, ".csv"))
    write.csv(results$results_df, file = results_path, row.names = FALSE)
    volcano_path <- check_dir_and_save(analysis_specific_dir, paste0("Interaction_volcano_", mode, analysis_config$file_suffix, ".png"))
    ggsave(volcano_path, results$volcano_plot, width = 10, height = 8, dpi = 300)
    cat("Successfully saved Interaction results for:", analysis_config$title, "\n")
  } else {
    cat("Skipping Interaction analysis for", analysis_config$title, "due to data issues.\n")
  }
}

run_interaction_analyses_for_mode <- function(mode, data_dir, result_dir, analysis_list) {
  cat("\n==================================================\n")
  cat("STARTING INTERACTION ANALYSIS FOR MODE:", toupper(mode), "\n")
  cat("==================================================\n")
  mode_string <- ifelse(mode %in% c("Positive", "POS", "pos"), "Positive", "Negative")
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", mode_string, ".csv"))
  metadata_file <- file.path(data_dir, paste0(mode_string, "_metadata_w_subjects.csv"))
  if (!file.exists(feature_table_file) || !file.exists(metadata_file)) {
    warning("Data files not found for Interaction mode '", mode, "'. Skipping this mode.", immediate. = TRUE)
    return()
  }
  feature_table <- read.csv(feature_table_file)
  metadata <- read.csv(metadata_file)
  for (analysis_config in analysis_list) {
    run_and_save_interaction_analysis(
      feature_table = feature_table, metadata = metadata, result_dir = result_dir,
      analysis_config = analysis_config, mode = mode
    )
  }
}


#======================== MAIN SCRIPT ========================

# --- 1. Configuration ---
data_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data"
result_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Univariate"

# --- 2. Define Kruskal-Wallis Analyses ---
kw_analysis_list <- list(
  
  list(title = "All Groups Comparison",
       grouping_column = "Group_Timepoint",
       groups = c("Healthy", "TR1", "TR2", "PR1", "PR2"),
       file_suffix = "_All_Groups", folder_name = "All_Groups_Comparison"),
  
  
  list(title = "Healthy vs TR vs PR at Timepoint 1",
       groups = c("Healthy", "TR", "PR"), timepoint = 1,
       file_suffix = "_H_TR_PR_T1", folder_name = "H_vs_TR_vs_PR_T1"),
  
  list(title = "Healthy vs TR vs PR at Timepoint 2",
       groups = c("Healthy", "TR", "PR"), timepoint = 2,
       file_suffix = "_H_TR_PR_T2", folder_name = "H_vs_TR_vs_PR_T2")
)

# --- 3. Define Wilcoxon Analyses ---
wilcoxon_analysis_list <- list(
  list(title = "Paired TR Timepoint 1 vs 2",
       group1 = "TR", timepoint1 = 1, group2 = "TR", timepoint2 = 2, paired = TRUE,
       file_suffix = "_TR1_vs_TR2_paired", folder_name = "TR1_vs_TR2_paired"),
  
  list(title = "Paired PR Timepoint 1 vs 2",
       group1 = "PR", timepoint1 = 1, group2 = "PR", timepoint2 = 2, paired = TRUE,
       file_suffix = "_PR1_vs_PR2_paired", folder_name = "PR1_vs_PR2_paired"),
  
  list(title = "Unpaired Healthy vs TR at Timepoint 1",
       group1 = "Healthy", timepoint1 = 1, group2 = "TR", timepoint2 = 1, paired = FALSE,
       file_suffix = "_H1_vs_TR1_unpaired", folder_name = "H1_vs_TR1_unpaired")
)

# --- 4. Define Interaction (Change-Score) Analyses ---
interaction_analysis_list <- list(
  list(title = "Interaction of Change (T2-T1) between TR and PR groups",
       group1 = "TR", group2 = "PR", timepoint1 = 1, timepoint2 = 2,
       file_suffix = "_Interaction_TR_vs_PR", folder_name = "Interaction_TR_vs_PR")
)


# --- 5. Execute All Analyses ---

# START logging to a file
log_file_path <- file.path(result_dir, "analysis_console_log.md")
sink(log_file_path, type = "output", split = TRUE) # split=TRUE shows output in console AND file

cat("# Univariate Analysis Log\n\n")
cat("This file contains the console output from the analysis script run on:", format(Sys.time()), "\n\n")

# Define mode-specific result directories
result_dir_neg <- file.path(result_dir, "Negative")
result_dir_pos <- file.path(result_dir, "Positive")

# Run for Negative Mode
cat("\n\n\n<<<<< RUNNING ALL ANALYSES FOR NEGATIVE MODE >>>>>\n")
run_kw_analyses_for_mode("neg", data_dir, result_dir_neg, kw_analysis_list)
run_wilcoxon_analyses_for_mode("neg", data_dir, result_dir_neg, wilcoxon_analysis_list)
run_interaction_analyses_for_mode("neg", data_dir, result_dir_neg, interaction_analysis_list)

# Run for Positive Mode
cat("\n\n\n<<<<< RUNNING ALL ANALYSES FOR POSITIVE MODE >>>>>\n")
run_kw_analyses_for_mode("pos", data_dir, result_dir_pos, kw_analysis_list)
run_wilcoxon_analyses_for_mode("pos", data_dir, result_dir_pos, wilcoxon_analysis_list)
run_interaction_analyses_for_mode("pos", data_dir, result_dir_pos, interaction_analysis_list)


cat("\nAll analyses completed successfully!\n")
sessionInfo()

# STOP logging to a file and return output to the console
sink()

