#----------------Cleaning the environment
rm(list=ls())

suppressPackageStartupMessages({
  library(tidyverse) # For data manipulation (dplyr) and plotting (ggplot2)
  library(vegan)       # For diversity analysis and ordination methods
  library(scales)
  })


#============Function============

#Function for saving. if the directory doesn't exist it will create new.
check_dir_and_save <- function(directory, filename) {
  # Construct full directory path
  full_dir <- directory
  
  # Create directory if it doesn't exist
  if (!dir.exists(full_dir)) {
    dir.create(full_dir, recursive = TRUE)
    message("Created directory: ", full_dir)
  }
  
  # Construct the full file path
  final_file_path <- file.path(full_dir, filename)
  
  
  return(final_file_path)
}
Preform_PCA <- function(feature_table, metadata, plot_title, groups_to_include = NULL){
  
  #=========== 1. Align Metadata and Feature Table ====================
  if (is.null(groups_to_include)) {
    # If no groups are specified, remove IGNORE/Blank by default.
    metadata_filtered <- metadata %>%
      filter(!Group %in% c("IGNORE", "INGNORE", "Blank"))
  } else {
    # Otherwise, filter for the user-specified groups.
    metadata_filtered <- metadata %>%
      filter(Group %in% groups_to_include)
  }
  
  # Extract the final list of sample names to keep.
  samples_to_keep <- metadata_filtered$filename
  
  # Filter the feature_table to keep ONLY the columns that might not be in the feature table for some reason.
  feature_table_filtered <- feature_table %>%
    select(1:3, any_of(samples_to_keep))
  
  # --- Sanity Check ---
  cat("Number of samples in filtered metadata:", nrow(metadata_filtered), "\n")
  cat("Number of samples in filtered feature table:", ncol(feature_table_filtered) - 3, "\n")
  
  if (nrow(metadata_filtered) != (ncol(feature_table_filtered) - 3)) {
    stop("ERROR: Sample mismatch after filtering. Check sample names.")
  }
  
  #=========== 2. Prepare Data for PCA ====================
  
  # Extract the final, ordered group information for the plot.
  Groups <- metadata_filtered$Group
  
  # Create unique feature IDs from the FILTERED feature table.
  feature_ids <- paste(
    feature_table_filtered$row.ID,
    round(feature_table_filtered$row.m.z, digits = 4),
    round(feature_table_filtered$row.retention.time, digits = 2),
    sep = "_"
  )
  
  # Create the numeric matrix for PCA, starting from the filtered data.
  peak_data <- as.matrix(feature_table_filtered[, -c(1:3)])
  
  # Transpose so samples are rows and features are columns.
  peak_data_t <- t(peak_data)
  
  # Log transform and scale the data.
  scaled_data <- scale(log1p(peak_data_t), center = TRUE, scale = TRUE)
  colnames(scaled_data) <- feature_ids
  
  
  #=========== 3. Perform PCA ====================
  
  pca_result <- prcomp(scaled_data, center = FALSE, scale. = FALSE)
  
  # Create data frame with PCA scores.
  pca_df <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2]
  )
  
  # Add sample names and groups.
  pca_df$sample_name <- rownames(peak_data_t)
  pca_df$Groups <- Groups
  
  #=========== 4. Process and Save Results ====================
  
  # Calculate variance explained for plot labels.
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pc1_var <- round(var_explained[1] * 100, 1)
  pc2_var <- round(var_explained[2] * 100, 1)
  
  # Create a summary of component importance.
  pca_importance_df <- as.data.frame(summary(pca_result)$importance) %>%
    t() %>%
    as_tibble(rownames = "principal_component")
  
  # Create a summary of feature loadings.
  pca_loadings_summary <- pca_result$rotation %>%
    as_tibble(rownames = "feature") %>%
    pivot_longer(
      cols = starts_with("PC"),
      names_to = "principal_component",
      values_to = "loading"
    )
  
  # Combine loadings and importance into a single file for saving.
  pca_for_save <- left_join(pca_loadings_summary, pca_importance_df, by = "principal_component")
  
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
  
  # Print the plot to the plot pane.
  print(pca_plot)
  
  # Return a named list containing all the important results.
  return(list(
    plot = pca_plot, 
    pca_scores = pca_df, 
    pca_loadings = pca_for_save
  ))
}

#===============Negative==========
#shows working directory
getwd()
#set the R working directory
setwd("C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data")
getwd()
#show files in directory
dir()
#Reading R Data Files+ naming the files

#Quant_Table output
feature_table <- read.csv("Imputed_TIC_Normalised_table_Negative.csv")
head(feature_table,2)

#metadata
metadata <- read.csv("NEGmetaData.csv")
head(metadata, 4)



#=====Overall -Negative =======
PCA_resluts <- Preform_PCA(feature_table,metadata, plot_title = "PCA of Serum Metabolites - Negative")
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Negative.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Negative_H_TR1_TR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Negative",c("Healthy", "TR1","TR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_H_TR1_TR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Negative_H_TR1_TR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_H_TR1_TR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Negative_TR1_TR2=======================================================

PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Negative",c("TR1","TR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_TR1_TR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Negative_TR1_TR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_TR1_TR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Negative_H_PR1_PR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Negative",c("Healthy", "PR1","PR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_H_PR1_PR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Negative_H_PR1_PR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_H_PR1_PR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Negative_PR1_PR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Negatvie",c("PR1","PR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_PR1_PR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Negative_PR1_PR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_PR1_PR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Negative_TR2_PR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Negative",c("TR2","PR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_TR2_PR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Negative_TR2_PR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative_TR2_PR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)



#===============Positive ==========

#shows working directory
getwd()
#set the R working directory
setwd("C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data")
getwd()
#show files in directory
dir()
#Reading R Data Files+ naming the files

#Quant_Table output
feature_table <- read.csv("Imputed_TIC_Normalised_table_Positive.csv")
head(feature_table,2)

#metadata
metadata <- read.csv("POSmetaData.csv")
head(metadata, 4)

#====Overall -Positive=========
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Positive")
PCA_resluts$pca_scores
PCA_resluts$pca_loadings
# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Positive.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Positive_H_TR1_TR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Positive",c("Healthy", "TR1","TR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_H_TR1_TR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Positive_H_TR1_TR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_H_TR1_TR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Positive_TR1_TR2=======================================================

PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Positive",c("TR1","TR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_TR1_TR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Positive_TR1_TR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_TR1_TR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Positive_H_PR1_PR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Positive",c("Healthy", "PR1","PR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_H_PR1_PR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Positive_H_PR1_PR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_H_PR1_PR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Positive_PR1_PR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Positive",c("PR1","PR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_PR1_PR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Positive_PR1_PR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_PR1_PR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)

#=====Positive_TR2_PR2=======================================================
PCA_resluts <- Preform_PCA(feature_table,metadata,plot_title = "PCA of Serum Metabolites - Positive",c("TR2","PR2"))
PCA_resluts$pca_scores
PCA_resluts$pca_loadings

# Saving Results

#Save the plot
plot_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_TR2_PR2.png")
ggsave(plot_to_save, PCA_resluts$plot, width = 10, height = 8, dpi = 300)


#The data for the PCA plot 
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_PC1andPC2_Positive_TR2_PR2.csv")
write.csv(PCA_resluts$pca_scores, file = final_file_path, row.names = FALSE)

#The data for deeper analyzation of the results
final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Positive_TR2_PR2.csv")
write.csv(PCA_resluts$pca_loadings, file = final_file_path, row.names = FALSE)
