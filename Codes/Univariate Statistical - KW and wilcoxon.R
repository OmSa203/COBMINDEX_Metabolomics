#======================== ENVIRONMENT & LIBRARIES ========================
#----------------Cleaning the environment
rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse) # For data manipulation (dplyr) and plotting (ggplot2)
  library(vegan)     # For diversity analysis and ordination methods
  library(scales)
  library(ggrepel)   # For non-overlapping labels in ggplot
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


#======================== KRUSKAL-WALLIS CORE & WRAPPERS ========================

# The core Kruskal-Wallis function.
perform_kruskal_posthoc_analysis <- function(feature_table, metadata, groups_to_compare,
                                             timepoint_to_analyze, mode,
                                             return_all_features = TRUE, p_adj_threshold = 0.05) {
  # --- Input Validation ---
  required_cols <- c("filename", "Sample_Type", "Timepoint", "Group")
  missing_cols <- setdiff(required_cols, names(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("The following required columns are missing from the metadata dataframe:", paste(missing_cols, collapse = ", ")))
  }
  if (!mode %in% c("pos", "neg")) {
    stop("The 'mode' argument must be either 'pos' or 'neg'.")
  }
  
  # --- Helper function to clean sample names dynamically based on mode ---
  clean_sample_name <- function(x, ion_mode) {
    pattern <- paste0("_", ion_mode, "\\.mzML$")
    gsub(pattern, "", x)
  }
  
  # --- Data Preparation and Analysis Pipeline ---
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
    inner_join(
      metadata %>%
        mutate(filename = clean_sample_name(filename, ion_mode = mode)) %>%
        filter(Sample_Type == "Sample"),
      by = "filename"
    ) %>%
    filter(Timepoint == timepoint_to_analyze & Group %in% groups_to_compare)
  
  if (nrow(prepared_data) == 0) {
    warning(paste("No data found for the specified filters (Timepoint:", timepoint_to_analyze,
                  ", Groups:", paste(groups_to_compare, collapse=", "), "). Returning an empty data frame."))
    return(tibble())
  }
  
  analysis_results <- prepared_data %>%
    group_by(m.z, retention_time) %>%
    nest() %>%
    mutate(
      kw_test_result = map(data, ~ tryCatch(kruskal.test(intensity ~ Group, data = .x), error = function(e) NULL)),
      statistic = map_dbl(kw_test_result, ~ if(is.null(.x)) NA_real_ else .x$statistic),
      p.value   = map_dbl(kw_test_result, ~ if(is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    select(-kw_test_result) %>%
    ungroup() %>%
    filter(!is.na(p.value)) %>%
    mutate(kw_p_adj_bh = p.adjust(p.value, method = "BH"))
  
  if (!return_all_features) {
    analysis_results <- analysis_results %>%
      filter(kw_p_adj_bh < p_adj_threshold)
    
    if (nrow(analysis_results) == 0) {
      warning("No features were significant after FDR correction. Returning an empty data frame.")
      return(tibble())
    }
  }
  
  final_results <- analysis_results %>%
    mutate(
      posthoc_test = map(data, ~ {
        pw_test <- pairwise.wilcox.test(.x$intensity, .x$Group, p.adjust.method = "BH")
        if (!is.null(pw_test$p.value)) {
          pw_test$p.value %>%
            as.data.frame() %>%
            rownames_to_column(var = "group1") %>%
            pivot_longer(cols = -group1, names_to = "group2", values_to = "posthoc_p_adj") %>%
            filter(!is.na(posthoc_p_adj))
        } else {
          NULL
        }
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

# Wrapper function to run a KW analysis and save the results.
run_and_save_kw_analysis <- function(feature_table, metadata, result_dir, analysis_config, mode) {
  cat("--------------------------------------------------\n")
  cat("Running Kruskal-Wallis Analysis:", analysis_config$title, "\n")
  
  results_df <- perform_kruskal_posthoc_analysis(
    feature_table = feature_table,
    metadata = metadata,
    groups_to_compare = analysis_config$groups,
    timepoint_to_analyze = analysis_config$timepoint,
    mode = mode,
    return_all_features = analysis_config$return_all,
    p_adj_threshold = analysis_config$p_adj_thresh
  )
  
  if (nrow(results_df) > 0) {
    analysis_specific_dir <- file.path(result_dir, analysis_config$folder_name)
    results_path <- check_dir_and_save(analysis_specific_dir, paste0("KW_results_", mode, analysis_config$file_suffix, ".csv"))
    write.csv(results_df, file = results_path, row.names = FALSE)
    cat("Successfully saved KW results for:", analysis_config$title, "\n")
  } else {
    cat("Skipping file save for", analysis_config$title, "as no results were generated.\n")
  }
}

# General function to run all KW analyses for a specific mode (pos/neg).
run_kw_analyses_for_mode <- function(mode, data_dir, result_dir, analysis_list) {
  cat("\n==================================================\n")
  cat("STARTING KRUSKAL-WALLIS ANALYSIS FOR MODE:", toupper(mode), "\n")
  cat("==================================================\n")
  
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", ifelse(mode =="Positive"| mode == "POS" | mode == "pos", "Positive", "Negative"), ".csv"))
  metadata_file <- file.path(data_dir, paste0(toupper(mode), "metadata_w_subjects.csv"))
  
  if (!file.exists(feature_table_file) || !file.exists(metadata_file)) {
    warning("Data files not found for KW mode '", mode, "'. Skipping this mode.", immediate. = TRUE)
    return()
  }
  
  feature_table <- read.csv(feature_table_file)
  metadata <- read.csv(metadata_file)
  
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

# The core Wilcoxon test function, adapted from your script.
perform_wilcoxon_analysis <- function(feature_table, metadata, config) {
  
  # --- 1. Filter metadata for specified groups and timepoints ---
  group1_samples <- metadata %>%
    filter(Group == config$group1, Timepoint == config$timepoint1) %>%
    arrange(Subject)
  
  group2_samples <- metadata %>%
    filter(Group == config$group2, Timepoint == config$timepoint2) %>%
    arrange(Subject)
  
  # --- 2. Determine sample pairing strategy ---
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
  
  # --- 3. Extract and verify data ---
  if (!all(c(sample1_cols, sample2_cols) %in% colnames(feature_table))) {
    warning("One or more sample columns not found in feature table for analysis: ", config$title)
    return(NULL)
  }
  group1_data <- feature_table[, sample1_cols, drop = FALSE]
  group2_data <- feature_table[, sample2_cols, drop = FALSE]
  
  # --- 4. Perform Wilcoxon test for each metabolite ---
  results_list <- lapply(1:nrow(feature_table), function(i) {
    g1_vals <- as.numeric(group1_data[i, ])
    g2_vals <- as.numeric(group2_data[i, ])
    
    # Skip if there's no data
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
  
  # --- 5. Calculate adjusted p-values, fold change, and effect size ---
  results_df$p_adjusted_BH <- p.adjust(results_df$p_value, method = "BH")
  results_df$mean_group1 <- rowMeans(group1_data, na.rm = TRUE)
  results_df$mean_group2 <- rowMeans(group2_data, na.rm = TRUE)
  results_df$fold_change <- results_df$mean_group2 / results_df$mean_group1
  results_df$log2_fold_change <- log2(results_df$fold_change)
  
  # --- 6. Generate Plots ---
  results_df$significance <- ifelse(results_df$p_adjusted_BH < 0.05, "Significant (FDR < 0.05)", "Not Significant")

  # Add a column for labels, identifying the top 10 most significant features by p-value
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

# Wrapper to run and save a single Wilcoxon analysis
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
    
    # Save results table
    results_path <- check_dir_and_save(analysis_specific_dir, paste0("Wilcoxon_results_", mode, analysis_config$file_suffix, ".csv"))
    write.csv(results$results_df, file = results_path, row.names = FALSE)
    
    # Save volcano plot
    volcano_path <- check_dir_and_save(analysis_specific_dir, paste0("Volcano_plot_", mode, analysis_config$file_suffix, ".png"))
    ggsave(volcano_path, results$volcano_plot, width = 10, height = 8, dpi = 300)
    
    # Save p-value histogram
    pval_hist_path <- check_dir_and_save(analysis_specific_dir, paste0("Pval_hist_", mode, analysis_config$file_suffix, ".png"))
    ggsave(pval_hist_path, results$p_value_hist, width = 8, height = 6, dpi = 300)
    
    cat("Successfully saved Wilcoxon results for:", analysis_config$title, "\n")
  } else {
    cat("Skipping Wilcoxon analysis for", analysis_config$title, "due to data issues.\n")
  }
}

# General function to run all Wilcoxon analyses for a specific mode
run_wilcoxon_analyses_for_mode <- function(mode, data_dir, result_dir, analysis_list) {
  cat("\n==================================================\n")
  cat("STARTING WILCOXON ANALYSIS FOR MODE:", toupper(mode), "\n")
  cat("==================================================\n")
  
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", ifelse(mode =="Positive"| mode == "POS" | mode == "pos", "Positive", "Negative"), ".csv"))
  metadata_file <- file.path(data_dir, paste0(toupper(mode), "metadata_w_subjects.csv"))
  
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

# Function to test for interaction effects (e.g., is the change over time different between groups?)
perform_interaction_analysis <- function(feature_table, metadata, config) {
  
  # --- 1. Helper function to get change scores for one group ---
  get_change_scores <- function(group_name) {
    # Filter metadata for the two timepoints within the specified group
    t1_samples <- metadata %>% filter(Group == group_name, Timepoint == config$timepoint1)
    t2_samples <- metadata %>% filter(Group == group_name, Timepoint == config$timepoint2)
    
    # Find subjects present at both timepoints
    paired_subjects <- intersect(t1_samples$Subject, t2_samples$Subject)
    if (length(paired_subjects) == 0) {
      warning("No paired subjects found for group: ", group_name)
      return(NULL)
    }
    
    # Get filenames for paired samples
    t1_paired_files <- sapply(paired_subjects, function(s) t1_samples$filename[t1_samples$Subject == s][1])
    t2_paired_files <- sapply(paired_subjects, function(s) t2_samples$filename[t2_samples$Subject == s][1])
    
    # Extract data and calculate difference (T2 - T1)
    t1_data <- feature_table[, t1_paired_files, drop = FALSE]
    t2_data <- feature_table[, t2_paired_files, drop = FALSE]
    change_scores <- t2_data - t1_data
    return(change_scores)
  }
  
  # --- 2. Get change scores for both groups ---
  group1_changes <- get_change_scores(config$group1)
  group2_changes <- get_change_scores(config$group2)
  
  if (is.null(group1_changes) || is.null(group2_changes)) {
    warning("Could not calculate change scores for one or both groups for analysis: ", config$title)
    return(NULL)
  }
  
  # --- 3. Perform Mann-Whitney U test on the change scores ---
  results_list <- lapply(1:nrow(feature_table), function(i) {
    g1_diffs <- as.numeric(group1_changes[i, ])
    g2_diffs <- as.numeric(group2_changes[i, ])
    
    if(all(is.na(g1_diffs)) || all(is.na(g2_diffs))) return(NULL)
    
    # Use unpaired Wilcoxon test to compare the two sets of differences
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
  
  # --- 4. Calculate adjusted p-values and mean changes ---
  results_df$p_adjusted_BH <- p.adjust(results_df$p_value, method = "BH")
  results_df$mean_change_group1 <- rowMeans(group1_changes, na.rm = TRUE)
  results_df$mean_change_group2 <- rowMeans(group2_changes, na.rm = TRUE)
  
  # This represents the difference of the differences
  # A small epsilon is added to avoid division by zero
  epsilon <- 1e-9 
  ratio <- (results_df$mean_change_group2 + epsilon) / (results_df$mean_change_group1 + epsilon)
  results_df$log2_fold_change_of_change <- log2(abs(ratio)) * sign(ratio)
  
  # --- 5. Generate Plots ---
  results_df$significance <- ifelse(results_df$p_adjusted_BH < 0.05, "Significant (FDR < 0.05)", "Not Significant")
  
  # Add a column for labels, identifying the top 10 most significant features by p-value
  results_df <- results_df %>%
    arrange(p_value) %>%
    mutate(label = if_else(row_number() <= 10, as.character(row.ID), NA_character_))
  
  volcano_plot <- ggplot(results_df, aes(x = log2_fold_change_of_change, y = -log10(p_value))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = c("Significant (FDR < 0.05)" = "purple", "Not Significant" = "gray")) +
    # Add reference lines for clarity
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    # Add labels to the top features
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

# Wrapper to run and save a single Interaction analysis
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

# General function to run all Interaction analyses for a specific mode
run_interaction_analyses_for_mode <- function(mode, data_dir, result_dir, analysis_list) {
  MODE <- ifelse(mode =="Positive"| mode == "POS" | mode == "pos", "Positive", "Negative")
  cat("\n==================================================\n")
  cat("STARTING INTERACTION ANALYSIS FOR MODE:", toupper(MODE), "\n")
  cat("==================================================\n")
  
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", ifelse(mode =="Positive"| mode == "POS" | mode == "pos", "Positive", "Negative"), ".csv"))
  metadata_file <- file.path(data_dir, paste0(toupper(mode), "metadata_w_subjects.csv"))
  
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
  list(title = "Healthy vs TR vs PR at Timepoint 1",
       groups = c("Healthy", "TR", "PR"), timepoint = 1, return_all = TRUE, p_adj_thresh = 0.05,
       file_suffix = "_H_TR_PR_T1", folder_name = "H_vs_TR_vs_PR_T1"),
  
  list(title = "Healthy vs TR vs PR at Timepoint 2",
       groups = c("Healthy", "TR", "PR"), timepoint = 2, return_all = TRUE, p_adj_thresh = 0.05,
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
# Run for Negative Mode
run_kw_analyses_for_mode("neg", data_dir, result_dir, kw_analysis_list)
run_wilcoxon_analyses_for_mode("neg", data_dir, result_dir, wilcoxon_analysis_list)
run_interaction_analyses_for_mode("neg", data_dir, result_dir, interaction_analysis_list)

# Run for Positive Mode
run_kw_analyses_for_mode("pos", data_dir, result_dir, kw_analysis_list)
run_wilcoxon_analyses_for_mode("pos", data_dir, result_dir, wilcoxon_analysis_list)
run_interaction_analyses_for_mode("pos", data_dir, result_dir, interaction_analysis_list)

cat("\nAll analyses completed successfully!\n")
sessionInfo()
