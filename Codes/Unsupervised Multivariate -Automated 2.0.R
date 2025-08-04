#----------------Cleaning the environment
rm(list=ls())

#----------------Loading Libraries
suppressPackageStartupMessages({
  library(tidyverse)  
  library(ggplot2)
  library(tidyr)
  library(pheatmap)
  library(vegan)
})


#======================== HELPER FUNCTIONS ========================

# Function for creating a directory and returning a full file path.
check_dir_and_save <- function(directory, filename) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    message("Created directory: ", directory)
  }
  file.path(directory, filename)
}

#======================== ANALYSIS FUNCTIONS ========================

# The core PCA function.
perform_pca <- function(feature_table, metadata, plot_title, groups_to_include = NULL){
  
  #=========== 1. Pre-process Metadata & Align with Feature Table ====================
  if ("Timepoint" %in% colnames(metadata) && "Group" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Group = if_else(
        Group %in% c("Healthy", "Blank", "IGNORE", "INGNORE"), 
        as.character(Group), 
        paste0(Group, Timepoint)
      ))
  }
  
  if (is.null(groups_to_include)) {
    metadata_to_use <- metadata %>%
      filter(!Group %in% c("IGNORE", "INGNORE", "Blank"))
  } else {
    metadata_to_use <- metadata %>%
      filter(Group %in% groups_to_include)
  }
  
  # ROBUST ALIGNMENT: Find common samples and filter both tables
  samples_in_metadata <- metadata_to_use$filename
  samples_in_feature_table <- colnames(feature_table)
  common_samples <- intersect(samples_in_metadata, samples_in_feature_table)
  
  if (length(common_samples) < 3) {
    stop("Fewer than 3 common samples found for this group.")
  }
  
  metadata_final <- metadata_to_use %>% filter(filename %in% common_samples)
  feature_table_final <- feature_table %>% select(1:3, all_of(common_samples))
  
  #=========== 2. Prepare Data for PCA ====================
  peak_data <- as.matrix(feature_table_final[, -c(1:3)])
  peak_data_t <- t(peak_data)
  scaled_data <- scale(log1p(peak_data_t), center = TRUE, scale = TRUE)
  
  #=========== 3. Perform PCA ====================
  pca_result <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  pca_df <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2],
    sample_name = rownames(scaled_data)
  )
  
  ordered_metadata <- metadata_final[match(pca_df$sample_name, metadata_final$filename), ]
  pca_df$Groups <- ordered_metadata$Group
  
  #=========== 4. Process and Save Results ====================
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pc1_var <- round(var_explained[1] * 100, 1)
  pc2_var <- round(var_explained[2] * 100, 1)
  
  pca_importance_df <- as.data.frame(summary(pca_result)$importance) %>%
    t() %>%
    as_tibble(rownames = "principal_component")
  
  pca_loadings_summary <- pca_result$rotation %>%
    as_tibble(rownames = "feature") %>%
    pivot_longer(
      cols = starts_with("PC"),
      names_to = "principal_component",
      values_to = "loading"
    )
  
  eigenvalues_df <- tibble(
    principal_component = colnames(pca_result$rotation),
    eigenvalue = pca_result$sdev^2
  )
  
  pca_for_save <- pca_loadings_summary %>%
    left_join(eigenvalues_df, by = "principal_component") %>%
    left_join(pca_importance_df, by = "principal_component")
  
  
  #=========== 5. Create Plot ====================
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Groups)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, type = "norm") +
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

# Function for Hierarchical Clustering Analysis (HCA)
perform_hca <- function(feature_table, metadata, plot_title, groups_to_include = NULL) {
  
  #=========== 1. Pre-process Metadata & Align with Feature Table ====================
  if ("Timepoint" %in% colnames(metadata) && "Group" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Group = if_else(
        Group %in% c("Healthy", "Blank", "IGNORE", "INGNORE"), 
        as.character(Group), 
        paste0(Group, Timepoint)
      ))
  }
  
  if (is.null(groups_to_include)) {
    metadata_to_use <- metadata %>%
      filter(!Group %in% c("IGNORE", "INGNORE", "Blank"))
  } else {
    metadata_to_use <- metadata %>%
      filter(Group %in% groups_to_include)
  }
  
  # ROBUST ALIGNMENT: Find common samples and filter both tables
  samples_in_metadata <- metadata_to_use$filename
  samples_in_feature_table <- colnames(feature_table)
  common_samples <- intersect(samples_in_metadata, samples_in_feature_table)
  
  if (length(common_samples) < 3) {
    stop("Fewer than 3 common samples found for this group.")
  }
  
  metadata_final <- metadata_to_use %>% filter(filename %in% common_samples)
  feature_table_final <- feature_table %>% select(1:3, all_of(common_samples))
  
  #=========== 2. Prepare Data for HCA ====================
  peak_data <- as.matrix(feature_table_final[, -c(1:3)])
  peak_data_t <- t(peak_data)
  scaled_data <- scale(log1p(peak_data_t), center = TRUE, scale = TRUE)
  
  top_var_features <- head(order(apply(scaled_data, 2, var), decreasing = TRUE), 50)
  scaled_data_subset <- scaled_data[, top_var_features]
  
  #=========== 3. Create Annotation for Plot ====================
  final_sample_names <- rownames(scaled_data_subset)
  ordered_metadata <- metadata_final[match(final_sample_names, metadata_final$filename), ]
  
  annotation_row <- data.frame(
    Group = ordered_metadata$Group,
    row.names = final_sample_names
  )
  
  #=========== 4. Perform HCA and Create Plot ====================
  hca_plot <- pheatmap(
    scaled_data_subset,
    annotation_row = annotation_row,
    main = plot_title,
    show_colnames = FALSE,
    fontsize_row = 8
  )
  
  #=========== 5. Package Results ====================
  sample_order_df <- tibble(
    sample_name = rownames(scaled_data_subset)[hca_plot$tree_row$order],
    cluster_order = 1:nrow(scaled_data_subset)
  )
  
  feature_order_df <- tibble(
    feature_name = colnames(scaled_data_subset)[hca_plot$tree_col$order],
    cluster_order = 1:ncol(scaled_data_subset)
  )
  
  # CLARITY FIX: Return list with descriptive names for HCA results.
  return(list(
    plot = hca_plot,
    sample_clustering = sample_order_df,
    feature_clustering = feature_order_df 
  ))
}

# UPDATED Function for Principal Coordinates Analysis (PCoA / MDS)
perform_pcoa <- function(feature_table, metadata, plot_title, groups_to_include = NULL, distance_method = "euclidean") {
  
  #=========== 1. Pre-process Metadata & Align with Feature Table ====================
  if ("Timepoint" %in% colnames(metadata) && "Group" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Group = if_else(
        Group %in% c("Healthy", "Blank", "IGNORE", "INGNORE"), 
        as.character(Group), 
        paste0(Group, Timepoint)
      ))
  }
  
  if (is.null(groups_to_include)) {
    metadata_to_use <- metadata %>%
      filter(!Group %in% c("IGNORE", "INGNORE", "Blank"))
  } else {
    metadata_to_use <- metadata %>%
      filter(Group %in% groups_to_include)
  }
  
  # ROBUST ALIGNMENT: Find common samples and filter both tables
  samples_in_metadata <- metadata_to_use$filename
  samples_in_feature_table <- colnames(feature_table)
  common_samples <- intersect(samples_in_metadata, samples_in_feature_table)
  
  if (length(common_samples) < 3) {
    stop("Fewer than 3 common samples found for this group.")
  }
  
  metadata_final <- metadata_to_use %>% filter(filename %in% common_samples)
  feature_table_final <- feature_table %>% select(1:3, all_of(common_samples))
  
  #=========== 2. Prepare Data for PCoA ====================
  peak_data <- as.matrix(feature_table_final[, -c(1:3)])
  peak_data_t <- t(peak_data)
  data_for_dist <- log1p(peak_data_t)
  
  #=========== 3. Perform PCoA & PERMANOVA ====================
  dist_matrix <- dist(data_for_dist, method = distance_method)
  
  # PERMANOVA
  ordered_metadata_for_adonis <- metadata_final[match(rownames(as.matrix(dist_matrix)), metadata_final$filename), ]
  adon_res <- vegan::adonis2(dist_matrix ~ Group, data = ordered_metadata_for_adonis)
  r_sq <- adon_res$R2[1]
  p_val <- adon_res$`Pr(>F)`[1]
  
  # PCoA
  pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE, add = TRUE)
  pcoa_df <- as.data.frame(pcoa_result$points)
  colnames(pcoa_df) <- c("PCo1", "PCo2")
  pcoa_df$sample_name <- rownames(pcoa_result$points)
  
  ordered_metadata <- metadata_final[match(pcoa_df$sample_name, metadata_final$filename), ]
  pcoa_df$Groups <- ordered_metadata$Group
  
  #=========== 4. Process and Save Results ====================
  variance <- round(pcoa_result$eig*100/sum(pcoa_result$eig),1)
  pco1_var <- variance[1]
  pco2_var <- variance[2]
  
  permanova_df <- tibble(
    R2 = r_sq,
    p_value = p_val
  )
  
  #=========== 5. Create Plot ====================
  plot_subtitle <- paste0("Metric: ", distance_method, " (p=", round(p_val, 4), ", R2=", round(r_sq, 4), ")")
  
  pcoa_plot <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Groups)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, type = "norm") +
    theme_bw() +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = paste0("PCo1 (", pco1_var, "%)"),
      y = paste0("PCo2 (", pco2_var, "%)"),
      color = "Groups"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  print(pcoa_plot)
  
  # CLARITY FIX: Return list with descriptive names for PCoA results.
  return(list(
    plot = pcoa_plot,
    pcoa_scores = pcoa_df,
    permanova_results = permanova_df
  ))
}


#======================== WRAPPER FUNCTIONS ========================

# General wrapper function to run an analysis and save the results
run_and_save_analysis <- function(analysis_function, feature_table, metadata, result_dir, analysis_config, mode, analysis_type) {
  
  tryCatch({
    args <- list(
      feature_table = feature_table,
      metadata = metadata,
      plot_title = analysis_config$title,
      groups_to_include = analysis_config$groups
    )
    
    if (!is.null(analysis_config$distance_method)) {
      args$distance_method <- analysis_config$distance_method
    }
    
    results <- do.call(analysis_function, args)
    
    file_suffix <- analysis_config$file_suffix
    folder_name <- analysis_config$folder_name
    analysis_specific_dir <- file.path(result_dir, folder_name)
    
    # --- 1. SAVE PLOT (BUG FIX) ---
    # The 'ggsave' function can correctly handle both ggplot objects (from PCA/PCoA) 
    # and pheatmap objects (from HCA). This removes the need for the incorrect
    # 'pheatmap(results$plot$gtable, ...)' call that was causing the error.
    plot_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_", mode, file_suffix, ".png"))
    ggsave(plot_path, results$plot, width = 10, height = 8, dpi = 300)
    
    # --- 2. SAVE DATA TABLES (CLARITY FIX) ---
    # This block now saves the data tables with descriptive filenames based on the 
    # analysis type, using the clearer returned list names.
    if (analysis_type == "PCA") {
      scores_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_scores_", mode, file_suffix, ".csv"))
      write.csv(results$pca_scores, file = scores_path, row.names = FALSE)
      
      loadings_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_loadings_", mode, file_suffix, ".csv"))
      write.csv(results$pca_loadings, file = loadings_path, row.names = FALSE)
      
    } else if (analysis_type == "HCA") {
      sample_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_sample_clustering_", mode, file_suffix, ".csv"))
      write.csv(results$sample_clustering, file = sample_path, row.names = FALSE)
      
      feature_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_feature_clustering_", mode, file_suffix, ".csv"))
      write.csv(results$feature_clustering, file = feature_path, row.names = FALSE)
      
    } else if (analysis_type == "PCoA") {
      scores_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_scores_", mode, file_suffix, ".csv"))
      write.csv(results$pcoa_scores, file = scores_path, row.names = FALSE)
      
      permanova_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_permanova_", mode, file_suffix, ".csv"))
      write.csv(results$permanova_results, file = permanova_path, row.names = FALSE)
    }
    
    cat("Successfully saved results for analysis:", analysis_config$title, "\n")
    
  }, error = function(e) {
    cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("ERROR: The following analysis failed and was skipped:\n")
    cat("  - Analysis Type:", analysis_type, "\n")
    cat("  - Mode:", mode, "\n")
    cat("  - Title:", analysis_config$title, "\n")
    cat("  - R Error Message:", conditionMessage(e), "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
  })
}

# Main loop function that handles data loading
run_analysis_for_mode <- function(mode, data_dir, result_dir, analysis_list, analysis_function, analysis_type) {
  
  cat("\n==================================================\n")
  cat("STARTING", toupper(analysis_type), "ANALYSIS FOR MODE:", mode, "\n")
  cat("==================================================\n")
  
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", mode, ".csv"))
  metadata_file <- file.path(data_dir, paste0(mode, "_metaData.csv"))
  
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
# !!! IMPORTANT: UPDATE THESE PATHS TO MATCH YOUR COMPUTER !!!
data_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data"
base_result_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate"

# Create separate main result folders for each analysis type
pca_result_dir <- file.path(base_result_dir, "PCA")
hca_result_dir <- file.path(base_result_dir, "HCA")
pcoa_result_dir <- file.path(base_result_dir, "PCoA")

# --- 2. Define the list of analyses to run ---
# This list is for PCA and HCA
analysis_definitions <- list(
  list(title = "Analysis of Serum Metabolites (All Groups)", groups = NULL, file_suffix = "_All", folder_name = "All_Groups"),
  list(title = "Analysis of Serum Metabolites (H, TR1, TR2)", groups = c("Healthy", "TR1", "TR2"), file_suffix = "_H_TR1_TR2", folder_name = "H_TR1_TR2"),
  list(title = "Analysis of Serum Metabolites (TR1, TR2)", groups = c("TR1", "TR2"), file_suffix = "_TR1_TR2", folder_name = "TR1_TR2"),
  list(title = "Analysis of Serum Metabolites (H, PR1, PR2)", groups = c("Healthy", "PR1", "PR2"), file_suffix = "_H_PR1_PR2", folder_name = "H_PR1_PR2"),
  list(title = "Analysis of Serum Metabolites (PR1, PR2)", groups = c("PR1", "PR2"), file_suffix = "_PR1_PR2", folder_name = "PR1_PR2"),
  list(title = "Analysis of Serum Metabolites (TR2, PR2)", groups = c("TR2", "PR2"), file_suffix = "_TR2_PR2", folder_name = "TR2_PR2")
)

# NEW: Specific list for PCoA to define distance metrics
pcoa_analysis_definitions <- list(
  list(title = "PCoA of All Groups", groups = NULL, file_suffix = "_All_euc", folder_name = "All_Groups_Euclidean", distance_method = "euclidean"),
  list(title = "PCoA of H, TR1, TR2", groups = c("Healthy", "TR1", "TR2"), file_suffix = "_H_TR1_TR2_euc", folder_name = "H_TR1_TR2_Euclidean", distance_method = "euclidean"),
  list(title = "PCoA of All Groups", groups = NULL, file_suffix = "_All_man", folder_name = "All_Groups_Manhattan", distance_method = "manhattan")
)


# --- 3. Execute Analyses for Each Mode and Type ---
for (current_mode in c("Negative", "Positive")) {
  
  # Run PCA
  run_analysis_for_mode(
    mode = current_mode,
    data_dir = data_dir,
    result_dir = pca_result_dir,
    analysis_list = analysis_definitions,
    analysis_function = perform_pca,
    analysis_type = "PCA"
  )
  
  # Run HCA
  run_analysis_for_mode(
    mode = current_mode,
    data_dir = data_dir,
    result_dir = hca_result_dir,
    analysis_list = analysis_definitions,
    analysis_function = perform_hca,
    analysis_type = "HCA"
  )
  
  # Run PCoA
  run_analysis_for_mode(
    mode = current_mode,
    data_dir = data_dir,
    result_dir = pcoa_result_dir,
    analysis_list = pcoa_analysis_definitions,
    analysis_function = perform_pcoa, 
    analysis_type = "PCoA"
  )
}

cat("\nAll analyses completed successfully!\n")
sessionInfo()
