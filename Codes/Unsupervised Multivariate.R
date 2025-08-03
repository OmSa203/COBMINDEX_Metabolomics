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


metadata <- subset(metadata,Group %in% c("Healthy","CD","UC"))
# Remove the first three columns (ID, m/z, RT) and convert to matrix
peak_data <- as.matrix(feature_table[, -c(1:3)])
head(peak_data)
# Transpose so samples are rows and features are columns
peak_data_t <- t(peak_data)
head(peak_data_t)
# Create sample names
sample_names <- rownames(peak_data_t)
metadata$sample_clean <- gsub(".mzML$", "", metadata$filename)


Groups <- filter(metadata, Group != "IGNORE" & Group != "Blank")
Groups <- Groups$Group

# Create unique feature IDs combining m/z and retention time
feature_ids <- paste(
  feature_table$row.ID,
  round(feature_table$row.m.z, digits = 4),
  round(feature_table$row.retention.time, digits = 2),
  sep = "_"
)

# Log transform and scale the data
scaled_data <- scale(log1p(peak_data_t))
colnames(scaled_data) <- feature_ids
# Perform PCA
pca_result <- prcomp(scaled_data, scale = TRUE)

# Combine loadings and eigenvalues into a single tibble
pca_summary <- pca_result$rotation %>%
  as_tibble(rownames = "feature") %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "principal_component",
    values_to = "loading"
  ) %>%
  mutate(
    eigenvalue = rep(pca_result$sdev^2, each = nrow(pca_result$rotation))
  )

# --- 3. View the Result ---
print(pca_summary)

# Create data frame with PCA results
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2]
  #Sample_clean = NA
)

# Add sample names as a column
pca_df$sample_clean <- sample_names
pca_df$Groups <- Groups
head(pca_df)

# Print diagnostic information
cat("Dimensions of PCA dataframe:", dim(pca_df), "\n")
cat("Number of NA values in Groups:", sum(is.na(pca_df$Groups)), "\n")

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)
total_var <- round(sum(var_explained[1:2]) * 100, 1)

# Define custom colors for 'Groups'
#custom_colors <- c("Healthy" = "#6d8abf", "TR1" = "#ccbb9c", "TR2" = "#ebe8e1", "PR1"= black, )

# Create PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Groups)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(
    title = " PCA of Serum Metabolites - Negative",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Groups",
    shape = "Groups"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) 
  #+ scale_color_manual(values = custom_colors)

pca_plot

# Add confidence ellipses (95%)
pca_plot <- pca_plot + 
  stat_ellipse(level = 0.95, type = "t")

pca_plot

# Save the plot

setwd()
file_to_save <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative.png")
ggsave(file_to_save, pca_plot, width = 10, height = 8, dpi = 300)

final_file_path <- check_dir_and_save(
  "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Unsupervised Multivariate/",
  "serum_pca_Negative.csv")
write.csv(pca_df, file = final_file_path, row.names = FALSE)
# Print variance information
cat("\nExplained variance:\n")
cat(paste("PC1:", pc1_var, "%\n"))
cat(paste("PC2:", pc2_var, "%\n"))
cat(paste("\nTotal variance explained by PC1 and PC2:", total_var, "%\n"))

# Print first few rows of final data to check
print(head(pca_df))



