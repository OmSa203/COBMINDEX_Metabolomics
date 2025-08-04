# Load required libraries

if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")


library(tidyverse)
library(tidyr)
library(reshape2)
library(stats)
library(ggplot2)
library(rstatix)
library(corrplot)

# Read the data
setwd("/Users/batelhagbi/Michal/Serum_20.2.25/")
dir()
metabolites <- read.csv("Features_with_high_zeros_Positive.csv")
metadata <- read.csv("POSmetaData.csv")

# Function to clean sample names
clean_sample_name <- function(x) {
  gsub("_pos\\.mzML$", "", x)
}

# Prepare metadata first
metadata_clean <- metadata %>%
  mutate(filename = clean_sample_name(filename)) %>%
  filter(Sample_Type == "Sample")  # Remove blanks

# Reshape metabolites data to long format
metabolites_long <- metabolites %>%
  select(-row.ID) %>%
  pivot_longer(
    cols = matches("S\\d+_pos\\.mzML"),
    names_to = "filename",
    values_to = "intensity"
  ) %>%
  mutate(
    filename = clean_sample_name(filename),
    m.z = row.m.z,
    retention_time = row.retention.time
  )

# Merge with metadata
merged_data <- metabolites_long %>%
  inner_join(metadata_clean, by = "filename")  # Changed to inner_join to ensure only matched samples

# Function to perform statistical analysis for each feature
analyze_feature <- function(data) {
  # Ensure we have the Group column
  if(!"Group" %in% names(data)) {
    stop("Group column not found in data")
  }
  
  # Kruskal-Wallis test between all groups
  tryCatch({
    kw_test <- kruskal.test(intensity ~ Group, data = data)
  }, error = function(e) {
    kw_test <- list(statistic = NA, p.value = NA, method = "Failed to compute")
  })
  
  # Modified pairwise Wilcoxon test with appropriate handling of ties
  if(length(unique(data$Group)) > 1) {
    pairwise_test <- tryCatch({
      pairwise.wilcox.test(
        data$intensity,
        data$Group,
        p.adjust.method = "BH",
        exact = FALSE  # Use approximation for ties
      )
    }, error = function(e) {
      NULL
    })
  } else {
    pairwise_test <- NULL
  }
  
  # Calculate comprehensive group statistics
  group_stats <- data %>%
    group_by(Group) %>%
    summarise(
      n_samples = n(),
      n_zeros = sum(intensity == 0),
      zero_percent = round(100 * sum(intensity == 0) / n(), 1),
      mean_intensity = mean(intensity),
      median_intensity = median(intensity),
      sd_intensity = sd(intensity),
      cv_percent = round(100 * sd_intensity / mean_intensity, 1),
      .groups = 'drop'
    ) %>%
    arrange(desc(median_intensity))
  
  # Calculate fold changes between groups
  if(length(unique(data$Group)) > 1) {
    group_medians <- group_stats$median_intensity
    group_names <- group_stats$Group
    
    # Create fold change matrix
    fc_matrix <- matrix(NA, 
                        nrow = length(group_names), 
                        ncol = length(group_names),
                        dimnames = list(group_names, group_names))
    
    for(i in 1:length(group_names)) {
      for(j in 1:length(group_names)) {
        if(i != j && group_medians[j] != 0) {
          fc_matrix[i,j] <- log2(group_medians[i] / group_medians[j])
        }
      }
    }
  } else {
    fc_matrix <- NULL
  }
  
  list(
    kw_test = kw_test,
    pairwise_test = pairwise_test,
    group_stats = group_stats,
    fold_changes = fc_matrix
  )
}

# Analyze each feature
features_analysis <- merged_data %>%
  group_by(m.z, retention_time) %>%
  nest() %>%
  mutate(
    analysis = map(data, safely(analyze_feature))
  ) %>%
  mutate(
    analysis = map(analysis, ~.$result)  # Extract successful results
  )

# Function to create visualization for a feature
plot_feature <- function(feature_data, mz, rt) {
  # Add small constant to handle zeros in log scale
  min_non_zero <- min(feature_data$intensity[feature_data$intensity > 0])
  plotting_data <- feature_data %>%
    mutate(plot_intensity = intensity + min_non_zero/10)
  
  ggplot(plotting_data, aes(x = Group, y = plot_intensity)) +
    geom_boxplot(aes(fill = Group)) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_y_log10() +
    labs(
      title = sprintf("Feature m/z %.4f, RT %.2f", mz, rt),
      y = "Intensity (log scale)",
      x = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 11)
    )
}

# Function to print analysis results for a feature
print_feature_analysis <- function(mz, rt, analysis_results) {
  cat(sprintf("\nFeature Analysis: m/z %.4f, RT %.2f\n", mz, rt))
  cat("\nKruskal-Wallis test:\n")
  print(analysis_results$kw_test)
  
  cat("\nGroup Statistics:\n")
  print(analysis_results$group_stats)
  
  if(!is.null(analysis_results$pairwise_test)) {
    cat("\nPairwise Wilcoxon tests (BH-adjusted p-values):\n")
    print(analysis_results$pairwise_test)
  }
  
  if(!is.null(analysis_results$fold_changes)) {
    cat("\nLog2 Fold Changes between groups:\n")
    print(round(analysis_results$fold_changes, 2))
  }
}

# Function to find features with significant differences
find_significant_features <- function(features_analysis, p_threshold = 0.05) {
  features_analysis %>%
    mutate(
      kw_pvalue = map_dbl(analysis, ~.$kw_test$p.value)
    ) %>%
    filter(!is.na(kw_pvalue) & kw_pvalue < p_threshold) %>%
    arrange(kw_pvalue)
}

# Get significant features
significant_features <- find_significant_features(features_analysis)

# Print summary of significant features
cat(sprintf("\nFound %d significant features (p < 0.05)\n", nrow(significant_features)))

# Export results with more detailed information
if(nrow(significant_features) > 0) {
  # Create detailed results dataframe
  detailed_results <- significant_features %>%
    mutate(
      group_stats = map_chr(analysis, ~paste(capture.output(print(.$group_stats)), collapse="\n")),
      kw_pvalue = map_dbl(analysis, ~.$kw_test$p.value)
    ) %>%
    select(m.z, retention_time, kw_pvalue, group_stats)
  
  # Create results directory if it doesn't exist
  dir.create("metabolomics_results", showWarnings = FALSE)
  
  # Save all significant features details
  write.csv(
    detailed_results,
    file.path("metabolomics_results", "significant_features_detailed.csv")
  )
  
  # Save summary statistics for all features
  all_stats <- significant_features %>%
    mutate(
      group_stats = map(analysis, ~.$group_stats),
      fold_changes = map(analysis, ~.$fold_changes)
    ) %>%
    select(m.z, retention_time, kw_pvalue, group_stats, fold_changes)
  
  # Save complete results workbook
  write.csv(
    all_stats %>% select(m.z, retention_time, kw_pvalue),
    file.path("metabolomics_results", "all_significant_features_summary.csv")
  )
  
  # Create separate directory for plots
  dir.create(file.path("metabolomics_results", "plots"), showWarnings = FALSE)
  
  # Analyze top 5 features
  top_features <- significant_features %>% 
    slice_head(n = 5)
  
  for(i in 1:nrow(top_features)) {
    feature <- top_features[i,]
    cat(sprintf("\n\nAnalyzing feature %d of 5:", i))
    print_feature_analysis(
      feature$m.z,
      feature$retention_time,
      feature$analysis[[1]]
    )
    
    # Create and save plot with more details
    p <- plot_feature(feature$data[[1]], feature$m.z, feature$retention_time)
    
    # Add p-value annotation
    p <- p + labs(
      subtitle = sprintf("Kruskal-Wallis p-value: %.2e", 
                         feature$analysis[[1]]$kw_test$p.value)
    )
    
    # Save plot with descriptive filename
    plot_filename <- sprintf(
      "feature_%d_mz%.4f_rt%.2f_pval%.2e.png",
      i,
      feature$m.z,
      feature$retention_time,
      feature$analysis[[1]]$kw_test$p.value
    )
    
    ggsave(
      file.path("metabolomics_results", "plots", plot_filename),
      p,
      width = 8,
      height = 6,
      dpi = 300
    )
    
    # Save feature-specific statistics
    feature_stats <- feature$analysis[[1]]$group_stats
    stats_filename <- sprintf(
      "feature_%d_mz%.4f_rt%.2f_statistics.csv",
      i,
      feature$m.z,
      feature$retention_time
    )
    write.csv(
      feature_stats,
      file.path("metabolomics_results", "plots", stats_filename)
    )
  }
}
#================================
"key aspects of how the script analyzes significance in your metabolomics data:

1. Statistical Tests Used:
- Primary test: Kruskal-Wallis test (`kruskal.test`)
   - Tests if there are any significant differences between your 5 groups (Healthy, TR1, TR2, PR1, PR2)
   - Non-parametric test chosen because metabolomics data often isn't normally distributed
   - P-value < 0.05 is considered significant

2. Secondary Analysis:
- Pairwise Wilcoxon tests (`pairwise.wilcox.test`)
   - Compares each pair of groups (e.g., TR1 vs TR2, TR1 vs PR1, etc.)
   - Uses Benjamini-Hochberg correction (BH) to adjust for multiple comparisons
   - Helps identify which specific groups differ from each other

3. What Makes a Feature Significant:
```r
find_significant_features <- function(features_analysis, p_threshold = 0.05) {
  features_analysis %>%
    mutate(
      kw_pvalue = map_dbl(analysis, ~.$kw_test$p.value)
    ) %>%
    filter(!is.na(kw_pvalue) & kw_pvalue < p_threshold) %>%
    arrange(kw_pvalue)
}
```
- A feature is considered significant if:
   - Kruskal-Wallis p-value < 0.05
   - The test successfully ran (no NA values)
   - The data contains valid measurements

4. Additional Metrics Calculated:
- For each group:
   - Number of samples (n_samples)
   - Number and percentage of zeros (n_zeros, zero_percent)
   - Mean and median intensity
   - Standard deviation (sd_intensity)
   - Coefficient of variation (cv_percent)
- Between groups:
   - Log2 fold changes between all group pairs
   - Shows magnitude of differences between groups

For example, in the output you'll see:
```
Feature Analysis: m/z 185.0233, RT 0.5246604

Kruskal-Wallis test:
p-value = 0.003  # Significant if < 0.05

Group Statistics:
Group    n_samples    n_zeros    mean_intensity    median_intensity
TR1      10          2          1000.5            950.3
TR2      10          3          800.2             750.1
...
```"