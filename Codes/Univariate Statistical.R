#----------------Cleaning the environment
rm(list=ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(cowplot)
  library(ggrepel)
  library(ComplexHeatmap)
})


#===========================***Functions***========================================

# Function to perform ANOVA and extract results
perform_anova <- function(feature_data, groups) {
  tryCatch({
    # Perform ANOVA
    anova_result <- aov(feature_data ~ groups)
    anova_summary <- summary(anova_result)[[1]]
    
    # Extract statistics
    f_value <- anova_summary$"F value"[1]
    p_value <- anova_summary$"Pr(>F)"[1]
    
    return(c(f_value, p_value))
  }, error = function(e) {
    return(c(NA, NA))
  })
}

# Function to perform Friedman and extract results
perform_friedman <- function(feature_data, groups) {
  tryCatch({
    # Perform Friedman
    friedman_result <- friedman.test(feature_data ~ groups)
    friedman_summary <- summary(friedman_result)[[1]]
    
    # Extract statistics
    f_value <- friedman_summary$"F value"[1]
    p_value <- friedman_summary$"Pr(>F)"[1]
    
    return(c(f_value, p_value))
  }, error = function(e) {
    return(c(NA, NA))
  })
}

# Function to create violin plots - p-adjusted
create_violin_plots <- function(feature_id, data, anova_results) {
  # Get feature statistics
  feature_stats <- anova_results[anova_results$Feature_ID == feature_id, ]
  
  # Create violin plot
  p <- ggplot(data[data$Feature_ID == feature_id, ], 
              aes(x = Group, y = Intensity, fill = Group, color = Group)) +
    geom_violin(trim = FALSE, alpha = 1) +
    geom_jitter(width = 0.2, size = 2, alpha = 1) +
    labs(
      title = paste("Peak Area of", feature_id),
      subtitle = paste("FDR adjusted p-value:", 
                       format(feature_stats$FDR_adjusted_P, scientific = TRUE, digits = 3)),
      x = "Group",
      y = "Peak Area"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    scale_fill_manual(values = c("Healthy" = "#6d8abf", "TR1" = "#ccbb9c", "TR2" = "#ebe8e1")) +
    scale_color_manual(values = c("Healthy" = "black", "TR1" = "black", "TR2" = "black"))
  
  # Add statistical comparisons
  p <- p + stat_compare_means(
    comparisons = list(
      c("Healthy", "TR1"),
      c("Healthy", "TR2"),
      c("TR1", "TR2")
    ),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  )
  
  return(p)
}

# Function to create violin plots - P-value
create_violin_plots_p_value <- function(feature_id, data, anova_results) {
  # Get feature statistics
  feature_stats <- anova_results[anova_results$Feature_ID == feature_id, ]
  
  # Create violin plot
  p <- ggplot(data[data$Feature_ID == feature_id, ], 
              aes(x = Group, y = `Peak Area`, fill = Group, color = Group)) +
    geom_violin(trim = FALSE, alpha = 1) +
    geom_jitter(width = 0.2, size = 2, alpha = 1) +
    labs(
      title = paste("Peak Area of", feature_id),
      subtitle = paste("p-value:", 
                       format(feature_stats$P_value, scientific = TRUE, digits = 3)),
      x = "Group",
      y = "Peak Area"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    scale_fill_manual(values = c("Healthy" = "#6d8abf", "TR1" = "#ccbb9c", "TR2" = "#ebe8e1")) +
    scale_color_manual(values = c("Healthy" = "black", "TR1" = "black", "TR2" = "black"))
  
  # Add statistical comparisons
  p <- p + stat_compare_means(
    comparisons = list(
      c("Healthy", "TR1"),
      c("Healthy", "TR2"),
      c("TR1", "TR2")
    ),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE
  )
  
  return(p)
}


#===========================For Negative========================================

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


# Create unique feature IDs combining m/z and retention time
feature_ids <- paste(
  feature_table$row.ID,
  round(feature_table$row.m.z, digits = 4),
  round(feature_table$row.retention.time, digits = 2),
  sep = "_"
)
head(feature_ids,4)

# Select only the intensity columns (those ending with mzML)
intensity_cols <- grep("mzML$", colnames(feature_table))
intensity_data <- feature_table[, intensity_cols]
rownames(intensity_data) <- feature_ids

# Transpose the data so samples are rows
intensity_data_t <- as.data.frame(t(intensity_data))

# Clean up sample names to match metadata
rownames(intensity_data_t) <- gsub("-", ".", rownames(intensity_data_t))

# Filter metadata to keep only the three groups of interest
metadata_filtered <- metadata %>%
  filter(Group %in% c("Healthy", "TR1", "TR2"))

# Make sure sample order matches between metadata and intensity data
intensity_data_t <- intensity_data_t[metadata_filtered$filename, ]


# Initialize results data frame
anova_results <- data.frame(
  Feature_ID = colnames(intensity_data_t),
  F_value = NA,
  P_value = NA,
  FDR_adjusted_P = NA,
  Significant = NA
)
# Perform ANOVA for each feature
for(i in 1:ncol(intensity_data_t)) {
  result <- perform_anova(intensity_data_t[,i], metadata_filtered$Group)
  anova_results$F_value[i] <- result[1]
  anova_results$P_value[i] <- result[2]
}


data <- as.data.frame(t(intensity_data))

data$filename <- rownames(data)

# 3. (Optional but good practice) Remove the now-redundant row names
rownames(data) <- NULL

# 4. Reorder the columns to put SampleID first
data <- data[, c("filename", setdiff(names(data), "filename"))]

Data <- merge(metadata_filtered,data,by.x="filename",by.y="filename")

library(dplyr)

# Check the structure of your data
# This will count how many times each subject appears in each group.
# A perfect dataset for Friedman's test will have a "count" of 1 for every row.
problem_check <- Data %>%
  group_by(sample.ID, Group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(sample.ID)

 # View the results.
# Look for counts that are not equal to 1, or for subjects who don't appear
# for all three groups.
print(problem_check, n = 50) # Print more rows to see the full picture
