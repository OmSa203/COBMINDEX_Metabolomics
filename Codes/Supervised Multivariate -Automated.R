#----------------Cleaning the environment
rm(list=ls())

#----------------Loading Libraries
# Ensure you have these packages installed: install.packages(c("ropls", "randomForest", "pals", "ggrepel"))
suppressPackageStartupMessages({
  library(tidyverse)  
  library(ggplot2)
  library(tidyr)
  library(ropls)
  library(randomForest)
  library(pals) # For color palettes
  library(ggrepel) # For non-overlapping plot labels
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

#======================== SUPERVISED ANALYSIS FUNCTIONS ========================

#' Perform OPLS-DA on feature data
#'
#' @param full_feature_table A data frame with features in rows and samples in columns. The first column must be a unique feature identifier.
#' @param metadata A data frame with sample metadata.
#' @param plot_title The main title for the plot.
#' @param groups_to_include A character vector of EXACTLY TWO groups to compare.
#' @return A list containing the S-plot, Scores plot, model summary, and VIP scores table, or NULL if model is not significant.
perform_oplsda <- function(full_feature_table, metadata, plot_title, groups_to_include) {
  
  #=========== 1. Validate Input & Pre-process Data ====================
  if (length(groups_to_include) != 2) {
    stop("OPLS-DA requires exactly two groups for comparison.")
  }
  
  # **IMPROVED**: Create a composite feature identifier in the format ID_mz_retention-time
  feature_identifiers <- paste(
    full_feature_table[[1]], 
    round(full_feature_table[[2]], 4), # Assumes column 2 is m/z
    round(full_feature_table[[3]], 2), # Assumes column 3 is retention time
    sep = "_"
  )
  feature_table <- full_feature_table[, -c(1:3)]
  rownames(feature_table) <- feature_identifiers
  
  if ("Timepoint" %in% colnames(metadata) && "Group" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Group = if_else(
        Group %in% c("Healthy", "Blank", "IGNORE", "INGNORE"), 
        as.character(Group), 
        paste0(Group, Timepoint)
      ))
  }
  
  metadata_to_use <- metadata %>%
    filter(Group %in% groups_to_include)
  
  # ROBUST ALIGNMENT
  common_samples <- intersect(metadata_to_use$filename, colnames(feature_table))
  if (length(common_samples) < 5) {
    stop("Fewer than 5 common samples found for the specified groups.")
  }
  
  metadata_final <- metadata_to_use %>% filter(filename %in% common_samples)
  feature_table_final <- feature_table %>% select(all_of(common_samples))
  
  # Order data consistently
  metadata_final <- metadata_final[match(colnames(feature_table_final), metadata_final$filename), ]
  
  #=========== 2. Prepare Data for OPLS-DA ====================
  peak_data <- as.matrix(feature_table_final)
  peak_data_t <- t(peak_data)
  scaled_data <- scale(log1p(peak_data_t), center = TRUE, scale = TRUE)
  
  # Response vector (the groups)
  y_vector <- as.factor(metadata_final$Group)
  
  #=========== 3. Perform OPLS-DA ====================
  set.seed(123)
  opls_model <- opls(x = scaled_data, y = y_vector, predI = 1, orthoI = NA)
  
  #=========== 4. Check for Model Validity ====================
  if (nrow(opls_model@summaryDF) == 0 || opls_model@summaryDF[, "pre"] < 1) {
    cat("INFO: OPLS-DA model for '", plot_title, "' was not significant or could not be built. Skipping.\n", sep="")
    return(NULL)
  }
  
  #=========== 5. Extract Results for Plotting & Saving ====================
  model_summary <- as_tibble(opls_model@summaryDF, rownames = "Metric")
  
  vip_scores <- as_tibble(opls_model@vipVn, rownames = "feature") %>%
    rename(VIP = value)
  
  #=========== 6. Create S-Plot with Labels ====================
  loadings_p1 <- opls_model@loadingMN %>% as_tibble(rownames = "feature") %>% select(feature, p1)
  corr_pcorr1 <- opls_model@corrMN %>% as_tibble(rownames = "feature") %>% select(feature, pcorr1)
  
  s_plot_data <- left_join(loadings_p1, corr_pcorr1, by = "feature") %>%
    left_join(vip_scores, by = "feature") %>%
    arrange(desc(VIP)) %>%
    mutate(feature_label = if_else(row_number() <= 10, feature, "")) # Label top 10 VIP features
  
  s_plot <- ggplot(s_plot_data, aes(x = p1, y = pcorr1)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(aes(label = feature_label), max.overlaps = 15, size = 3) +
    labs(
      title = paste("S-Plot:", plot_title),
      subtitle = paste0("R2X=", round(model_summary$`R2X(cum)`, 2), 
                        ", R2Y=", round(model_summary$`R2Y(cum)`, 2), 
                        ", Q2=", round(model_summary$`Q2(cum)`, 2)),
      x = "Modeled Covariance p(1)",
      y = "Modeled Correlation p(corr)[1]"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  print(s_plot)
  
  #=========== 7. Create Scores Plot with Labels ====================
  scores_data <- opls_model@scoreMN %>%
    as_tibble(rownames = "sample_id") %>%
    left_join(metadata_final %>% select(filename, Group), by = c("sample_id" = "filename"))
  
  scores_plot <- ggplot(scores_data, aes(x = t1, y = to1, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(aes(label = sample_id), max.overlaps = 15, size = 2.5) +
    stat_ellipse(type = "t", level = 0.95, linetype = 2, geom = "path") +
    labs(
      title = paste("Scores Plot:", plot_title),
      subtitle = paste0("R2X=", round(model_summary$`R2X(cum)`, 2), 
                        ", R2Y=", round(model_summary$`R2Y(cum)`, 2), 
                        ", Q2=", round(model_summary$`Q2(cum)`, 2)),
      x = "t1 (Predictive Component)",
      y = "to1 (First Orthogonal Component)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  print(scores_plot)
  
  #=========== 8. Return All Results ====================
  return(list(
    s_plot = s_plot,
    scores_plot = scores_plot,
    model_summary = model_summary,
    vip_scores = vip_scores
  ))
}

#' Perform Random Forest classification
#'
#' @param full_feature_table A data frame with features in rows and samples in columns. The first column must be a unique feature identifier.
#' @param metadata A data frame with sample metadata.
#' @param plot_title The main title for the plot.
#' @param groups_to_include A character vector of groups to include in the model.
#' @return A list containing the feature importance plot, confusion matrix, and full importance table.
perform_random_forest <- function(full_feature_table, metadata, plot_title, groups_to_include) {
  
  #=========== 1. Pre-process Data ====================
  # **IMPROVED**: Create a composite feature identifier in the format ID_mz_retention-time
  feature_identifiers <- paste(
    full_feature_table[[1]], 
    round(full_feature_table[[2]], 4), # Assumes column 2 is m/z
    round(full_feature_table[[3]], 2), # Assumes column 3 is retention time
    sep = "_"
  )
  feature_table <- full_feature_table[, -c(1:3)]
  rownames(feature_table) <- feature_identifiers
  
  if ("Timepoint" %in% colnames(metadata) && "Group" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Group = if_else(
        Group %in% c("Healthy", "Blank", "IGNORE", "INGNORE"), 
        as.character(Group), 
        paste0(Group, Timepoint)
      ))
  }
  
  metadata_to_use <- metadata %>%
    filter(Group %in% groups_to_include)
  
  # ROBUST ALIGNMENT
  common_samples <- intersect(metadata_to_use$filename, colnames(feature_table))
  if (length(common_samples) < 5) {
    stop("Fewer than 5 common samples found for the specified groups.")
  }
  
  metadata_final <- metadata_to_use %>% filter(filename %in% common_samples)
  feature_table_final <- feature_table %>% select(all_of(common_samples))
  
  # Order data consistently
  metadata_final <- metadata_final[match(colnames(feature_table_final), metadata_final$filename), ]
  
  #=========== 2. Prepare Data for Random Forest ====================
  # Transposing the data now carries over the meaningful feature names as column names
  peak_data_t <- t(as.matrix(feature_table_final))
  
  # Response vector must be a factor for classification
  y_vector <- as.factor(metadata_final$Group)
  
  #=========== 3. Perform Random Forest ====================
  set.seed(123) # For reproducible results
  rf_model <- randomForest(x = peak_data_t, y = y_vector, importance = TRUE)
  
  #=========== 4. Extract Results for Plotting & Saving ====================
  # The 'feature' column now contains your meaningful identifiers
  importance_df <- as.data.frame(importance(rf_model)) %>%
    as_tibble(rownames = "feature") %>%
    arrange(desc(MeanDecreaseAccuracy))
  
  confusion_matrix <- rf_model$confusion %>%
    as.data.frame()
  
  #=========== 5. Create Feature Importance Plot ====================
  top_20_features <- head(importance_df, 20)
  
  importance_plot <- ggplot(top_20_features, aes(x = reorder(feature, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +
    labs(
      title = plot_title,
      subtitle = "Top 20 Most Important Features",
      x = "Feature",
      y = "Mean Decrease in Accuracy"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  print(importance_plot)
  
  return(list(
    plot = importance_plot,
    confusion_matrix = confusion_matrix,
    importance_scores = importance_df
  ))
}


#======================== WRAPPER FUNCTION ========================

# General wrapper to run a supervised analysis and save results
run_and_save_supervised_analysis <- function(analysis_function, feature_table, metadata, result_dir, analysis_config, mode, analysis_type) {
  
  tryCatch({
    # **CHANGED**: Pass the full feature table to the analysis functions
    results <- analysis_function(
      full_feature_table = feature_table,
      metadata = metadata,
      plot_title = analysis_config$title,
      groups_to_include = analysis_config$groups
    )
    
    if (is.null(results)) {
      return()
    }
    
    file_suffix <- analysis_config$file_suffix
    folder_name <- analysis_config$folder_name
    analysis_specific_dir <- file.path(result_dir, folder_name)
    
    # --- 1. SAVE PLOTS ---
    if (analysis_type == "OPLSDA") {
      s_plot_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_S-plot_", mode, file_suffix, ".png"))
      ggsave(s_plot_path, results$s_plot, width = 10, height = 8, dpi = 300)
      
      scores_plot_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_Scores-plot_", mode, file_suffix, ".png"))
      ggsave(scores_plot_path, results$scores_plot, width = 10, height = 8, dpi = 300)
    } else { # For Random Forest
      plot_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_plot_", mode, file_suffix, ".png"))
      ggsave(plot_path, results$plot, width = 10, height = 8, dpi = 300)
    }
    
    # --- 2. SAVE DATA TABLES ---
    if (analysis_type == "OPLSDA") {
      summary_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_summary_", mode, file_suffix, ".csv"))
      write.csv(results$model_summary, file = summary_path, row.names = FALSE)
      
      vip_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_VIPscores_", mode, file_suffix, ".csv"))
      write.csv(results$vip_scores, file = vip_path, row.names = FALSE)
      
    } else if (analysis_type == "RF") {
      confusion_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_confusion_matrix_", mode, file_suffix, ".csv"))
      write.csv(results$confusion_matrix, file = confusion_path, row.names = TRUE)
      
      importance_path <- check_dir_and_save(analysis_specific_dir, paste0(analysis_type, "_importance_", mode, file_suffix, ".csv"))
      write.csv(results$importance_scores, file = importance_path, row.names = FALSE)
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

#======================== MAIN SCRIPT ========================

# --- 1. Configuration ---
# !!! IMPORTANT: UPDATE THESE PATHS TO MATCH YOUR COMPUTER !!!
data_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data"
base_result_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Supervised_Multivariate"

# Create separate main result folders for each analysis type
oplsda_result_dir <- file.path(base_result_dir, "OPLSDA")
rf_result_dir <- file.path(base_result_dir, "RandomForest")

# --- 2. Define the list of analyses to run ---
# This list defines pairwise comparisons for OPLS-DA and can also be used for Random Forest.
supervised_analysis_definitions <- list(
  list(title = "OPLS-DA of Healthy vs TR1", groups = c("Healthy", "TR1"), file_suffix = "_H_vs_TR1", folder_name = "H_vs_TR1"),
  list(title = "OPLS-DA of Healthy vs PR1", groups = c("Healthy", "PR1"), file_suffix = "_H_vs_PR1", folder_name = "H_vs_PR1"),
  list(title = "OPLS-DA of TR1 vs TR2", groups = c("TR1", "TR2"), file_suffix = "_TR1_vs_TR2", folder_name = "TR1_vs_TR2"),
  list(title = "OPLS-DA of PR1 vs PR2", groups = c("PR1", "PR2"), file_suffix = "_PR1_vs_PR2", folder_name = "PR1_vs_PR2"),
  list(title = "OPLS-DA of TR2 vs PR2", groups = c("TR2", "PR2"), file_suffix = "_TR2_vs_PR2", folder_name = "TR2_vs_PR2"),
  # This one is for Random Forest, including more groups
  list(title = "Random Forest of H, TR1, TR2", groups = c("Healthy", "TR1", "TR2"), file_suffix = "_H_TR1_TR2", folder_name = "H_TR1_TR2_Classifier")
)

# --- 3. Execute Analyses for Each Mode and Type ---
for (current_mode in c("Negative", "Positive")) {
  
  # --- Load Data for the current mode ---
  feature_table_file <- file.path(data_dir, paste0("Imputed_TIC_Normalised_table_", current_mode, ".csv"))
  metadata_file <- file.path(data_dir, paste0(current_mode, "_metaData.csv"))
  
  if (!file.exists(feature_table_file) || !file.exists(metadata_file)) {
    warning("Data files not found for mode '", current_mode, "'. Skipping.")
    next
  }
  feature_table <- read.csv(feature_table_file)
  metadata <- read.csv(metadata_file)
  
  # --- Run OPLS-DA Analyses ---
  cat("\n==================================================\n")
  cat("STARTING OPLS-DA ANALYSIS FOR MODE:", current_mode, "\n")
  cat("==================================================\n")
  for (config in supervised_analysis_definitions) {
    # OPLS-DA is strictly for 2 groups
    if(length(config$groups) == 2) {
      run_and_save_supervised_analysis(
        analysis_function = perform_oplsda,
        feature_table = feature_table,
        metadata = metadata,
        result_dir = oplsda_result_dir,
        analysis_config = config,
        mode = current_mode,
        analysis_type = "OPLSDA"
      )
    }
  }
  
  # --- Run Random Forest Analyses ---
  cat("\n==================================================\n")
  cat("STARTING RANDOM FOREST ANALYSIS FOR MODE:", current_mode, "\n")
  cat("==================================================\n")
  for (config in supervised_analysis_definitions) {
    # Update title for RF context
    rf_config <- config
    rf_config$title <- gsub("OPLS-DA", "Random Forest", rf_config$title)
    
    run_and_save_supervised_analysis(
      analysis_function = perform_random_forest,
      feature_table = feature_table,
      metadata = metadata,
      result_dir = rf_result_dir,
      analysis_config = rf_config,
      mode = current_mode,
      analysis_type = "RF"
    )
  }
}

cat("\nAll supervised analyses completed successfully!\n")
sessionInfo()

