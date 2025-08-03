#----------------Cleaning the environment
rm(list=ls())

#----------------Loading Libraries
suppressPackageStartupMessages({
  library(tidyverse) 
  library(ggplot2)
  library(tidyr)
})


#======================== FUNCTIONS ========================

# Function for creating a directory and returning a full file path.
check_dir_and_save <- function(directory, filename) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    message("Created directory: ", directory)
  }
  file.path(directory, filename)
}

# The core PCA function. It is now more focused on just the analysis.
perform_pca <- function(feature_table, metadata, plot_title, groups_to_include = NULL){
  
  #=========== 1. Align Metadata and Feature Table ====================
  if (is.null(groups_to_include)) {
    metadata_filtered <- metadata %>%
      filter(!Group %in% c("IGNORE", "INGNORE", "Blank"))
  } else {
    metadata_filtered <- metadata %>%
      filter(Group %in% groups_to_include)
  }
  
  samples_to_keep <- metadata_filtered$filename
  
  feature_table_filtered <- feature_table %>%
    select(1:3, any_of(samples_to_keep))
  
  # --- Sanity Check ---
  cat("--------------------------------------------------\n")
  cat("Running Analysis for Groups:", paste(if(is.null(groups_to_include)) "All" else groups_to_include, collapse=", "), "\n")
  cat("Number of samples in filtered metadata:", nrow(metadata_filtered), "\n")
  cat("Number of samples in filtered feature table:", ncol(feature_table_filtered) - 3, "\n")
  
  if (nrow(metadata_filtered) != (ncol(feature_table_filtered) - 3)) {
    stop("ERROR: Sample mismatch after filtering. Check sample names.")
  }
  
  #=========== 2. Prepare Data for PCA ====================
  Groups <- metadata_filtered$Group
  feature_ids <- paste(
    feature_table_filtered$row.ID,
    round(feature_table_filtered$row.m.z, digits = 4),
    round(feature_table_filtered$row.retention.time, digits = 2),
    sep = "_"
  )
  peak_data <- as.matrix(feature_table_filtered[, -c(1:3)])
  peak_data_t <- t(peak_data)
  scaled_data <- scale(log1p(peak_data_t), center = TRUE, scale = TRUE)
  colnames(scaled_data) <- feature_ids
  
  #=========== 3. Perform PCA ====================
  pca_result <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  pca_df <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2]
  )
  pca_df$sample_name <- rownames(peak_data_t)
  pca_df$Groups <- Groups
  
  #=========== 4. Process and Save Results ====================
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pc1_var <- round(var_explained[1] * 100, 1)
  pc2_var <- round(var_explained[2] * 100, 1)
  
  # Create a summary of component importance (proportion of variance, etc.)
  pca_importance_df <- as.data.frame(summary(pca_result)$importance) %>%
    t() %>%
    as_tibble(rownames = "principal_component")
  
  # Create a summary of feature loadings
  pca_loadings_summary <- pca_result$rotation %>%
    as_tibble(rownames = "feature") %>%
    pivot_longer(
      cols = starts_with("PC"),
      names_to = "principal_component",
      values_to = "loading"
    )
  
  # Create a table of eigenvalues for each principal component
  eigenvalues_df <- tibble(
    principal_component = colnames(pca_result$rotation),
    eigenvalue = pca_result$sdev^2
  )
  
  # Combine loadings, eigenvalues, and importance into a single file for saving.
  pca_for_save <- pca_loadings_summary %>%
    left_join(eigenvalues_df, by = "principal_component") %>%
    left_join(pca_importance_df, by = "principal_component")
  
  
  #=========== 5. Create Plot ====================
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Groups)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, type = "t") +
    theme_bw() +
    labs(
      title = plot_title,
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)"),
      color = "Groups"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  print(pca_plot)
  
  return(list(
    plot = pca_plot, 
    pca_scores = pca_df, 
    pca_loadings = pca_for_save
  ))
}


# General wrapper function to run an analysis and save the results
run_and_save_analysis <- function(analysis_function, feature_table, metadata, result_dir, analysis_config, mode, analysis_type) {
  
  # Run the specified analysis function (e.g., perform_pca)
  results <- analysis_function(
    feature_table = feature_table,
    metadata = metadata,
    plot_title = analysis_config$title,
    groups_to_include = analysis_config$groups
  )
  
  # --- Saving Results ---
  file_suffix <- analysis_config$file_suffix
  folder_name <- analysis_config$folder_name # Get the folder name for the subdirectory
  
  # Create the specific subdirectory path for this analysis
  analysis_specific_dir <- file.path(result_dir, folder_name)
  
  # Save the plot with a dynamic filename based on mode and analysis type
  plot_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_", mode, file_suffix, ".png"))
  ggsave(plot_path, results$plot, width = 10, height = 8, dpi = 300)
  
  # Save the scores
  scores_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_scores_", mode, file_suffix, ".csv"))
  write.csv(results$pca_scores, file = scores_path, row.names = FALSE)
  
  # Save the loadings
  loadings_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_loadings_", mode, file_suffix, ".csv"))
  write.csv(results$pca_loadings, file = loadings_path, row.names = FALSE)
  
  cat("Successfully saved results for analysis:", analysis_config$title, "\n")
}


#It now takes the list of analyses as an argument.
run_analysis_for_mode <- function(mode, data_dir, result_dir, analysis_list, analysis_function, analysis_type) {
  
  cat("\n==================================================\n")
  cat("STARTING", toupper(analysis_type), "ANALYSIS FOR MODE:", mode, "\n")
  cat("==================================================\n")
  
  # --- 1. Construct dynamic file paths ---
  # NOTE: Assumes metadata files are named like "Negative_metaData.csv"
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", mode, ".csv"))
  metadata_file <- file.path(data_dir, paste0(mode, "_metaData.csv"))
  
  # --- 2. Load Data ---
  if (!file.exists(feature_table_file)) {
    warning("Feature table file not found for mode '", mode, "'. Skipping. Path: ", feature_table_file)
    return()
  }
  if (!file.exists(metadata_file)) {
    warning("Metadata file not found for mode '", mode, "'. Skipping. Path: ", metadata_file)
    return()
  }
  
  feature_table <- read.csv(feature_table_file)
  metadata <- read.csv(metadata_file)
  
  # --- 3. Execute All Analyses in a Loop ---
  # The loop now iterates over the list of analyses passed into the function.
  for (analysis_config in analysis_list) {
    run_and_save_analysis(
      analysis_function = analysis_function,
      feature_table = feature_table,
      metadata = metadata,
      result_dir = result_dir,
      analysis_config = analysis_config,
      mode = mode,
      analysis_type = analysis_type
    )
  }
}


#======================== MAIN SCRIPT ========================

# --- 1. Configuration ---
data_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data"
result_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate"

# --- 2. Define the list of PCA analyses you want to run ---
# A 'folder_name' key has been added to create a sub-directory for each analysis.
pca_analysis_list <- list(
  list(title = "PCA of Serum Metabolites (All Groups)", groups = NULL, file_suffix = "_All", folder_name = "All_Groups"),
  list(title = "PCA of Serum Metabolites (H, TR1, TR2)", groups = c("Healthy", "TR1", "TR2"), file_suffix = "_H_TR1_TR2", folder_name = "H_TR1_TR2"),
  list(title = "PCA of Serum Metabolites (TR1, TR2)", groups = c("TR1", "TR2"), file_suffix = "_TR1_TR2", folder_name = "TR1_TR2"),
  list(title = "PCA of Serum Metabolites (H, PR1, PR2)", groups = c("Healthy", "PR1", "PR2"), file_suffix = "_H_PR1_PR2", folder_name = "H_PR1_PR2"),
  list(title = "PCA of Serum Metabolites (PR1, PR2)", groups = c("PR1", "PR2"), file_suffix = "_PR1_PR2", folder_name = "PR1_PR2"),
  list(title = "PCA of Serum Metabolites (TR2, PR2)", groups = c("TR2", "PR2"), file_suffix = "_TR2_PR2", folder_name = "TR2_PR2")
)


# --- 3. Execute Analyses for Each Mode ---
{
  # Negative Mode
  run_analysis_for_mode(
    mode = "Negative", 
    data_dir = data_dir, 
    result_dir = result_dir,
    analysis_list = pca_analysis_list,
    analysis_function = perform_pca,
    analysis_type = "PCA"
  )
  
  # Positive Mode
  run_analysis_for_mode(
    mode = "Positive", 
    data_dir = data_dir, 
    result_dir = result_dir,
    analysis_list = pca_analysis_list,
    analysis_function = perform_pca,
    analysis_type = "PCA"
  )
  
  cat("\nAll analyses completed successfully!\n")
  sessionInfo()
}