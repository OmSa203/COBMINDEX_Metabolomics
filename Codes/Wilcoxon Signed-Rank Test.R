# Paired Wilcoxon Signed-Rank Test for TR1 vs TR2 Metabolomics Data
# Author: Analysis script for Crohn's disease metabolomics study
# Date: August 2025

# --- 1. Configuration ---
data_dir <- "/Users/batelhagbi/Documents/GitHub/COBMINDEX_Metabolomics/Data"
result_dir <- "/Users/batelhagbi/Documents/GitHub/COBMINDEX_Metabolomics/results/Univariate"

# --- 2. Analysis Configuration ---
# Define which comparison to perform
# Options for groups: "TR", "PR", "Healthy"
# Options for timepoints: 1, 2
# For paired analysis: set group1 = group2, timepoint1 ≠ timepoint2
# For unpaired analysis between groups: set group1 ≠ group2, timepoint1 = timepoint2

group1 <- "TR"        # First group (e.g., "TR", "PR", "Healthy")
timepoint1 <- 1       # First timepoint (1 = before treatment, 2 = after treatment)
group2 <- "TR"        # Second group 
timepoint2 <- 2       # Second timepoint

# Analysis type
paired_analysis <- TRUE  # Set to TRUE for paired analysis (same subjects), FALSE for unpaired

print(paste("Analysis Configuration:"))
print(paste("Comparison:", paste0(group1, timepoint1), "vs", paste0(group2, timepoint2)))
print(paste("Analysis type:", ifelse(paired_analysis, "Paired", "Unpaired")))

# Create results directory if it doesn't exist
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

# Set working directory to data directory
setwd(data_dir)

# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# Read the data files
print("Loading data files...")
metadata <- read_csv("POSmetadata_w_subjects.csv")
metabolomics_data <- read_csv("Imputed_TIC_Normalised_table_Positive.csv")

# Clean column names in metadata (remove spaces)
colnames(metadata) <- trimws(colnames(metadata))

# Display data structure
print("Metadata structure:")
print(head(metadata))
print(paste("Number of metabolites:", nrow(metabolomics_data)))
print(paste("Number of samples in metabolomics data:", ncol(metabolomics_data) - 3))

# Filter metadata for specified groups and timepoints
group1_samples <- metadata %>%
  filter(Group == group1, Timepoint == timepoint1) %>%
  arrange(Subject)

group2_samples <- metadata %>%
  filter(Group == group2, Timepoint == timepoint2) %>%
  arrange(Subject)

print(paste("Group 1 samples (", group1, timepoint1, "):", sep=""))
print(group1_samples)
print(paste("Group 2 samples (", group2, timepoint2, "):", sep=""))
print(group2_samples)

# Determine sample pairing strategy
if (paired_analysis) {
  # For paired analysis, match samples by subject
  if (group1 == group2) {
    print("Performing paired analysis within the same group across timepoints")
  } else {
    print("Performing paired analysis between different groups at different timepoints")
  }
  
  # Find subjects that have both conditions
  subjects_group1 <- group1_samples$Subject
  subjects_group2 <- group2_samples$Subject
  paired_subjects <- intersect(subjects_group1, subjects_group2)
  
  if (length(paired_subjects) == 0) {
    stop("No subjects found with both conditions for paired analysis!")
  }
  
  # Create paired sample mapping
  sample_pairs <- data.frame(
    Subject = paired_subjects,
    stringsAsFactors = FALSE
  )
  
  # Get sample names for each subject
  sample_pairs$Sample1 <- sapply(paired_subjects, function(subj) {
    group1_samples$filename[group1_samples$Subject == subj]
  })
  
  sample_pairs$Sample2 <- sapply(paired_subjects, function(subj) {
    group2_samples$filename[group2_samples$Subject == subj]
  })
  
  # Keep .mzML extension for column matching
  sample_pairs$Sample1_col <- sample_pairs$Sample1
  sample_pairs$Sample2_col <- sample_pairs$Sample2
  
  print(paste("Found", length(paired_subjects), "paired subjects:"))
  print(sample_pairs)
  
  # Extract sample column names
  sample1_columns <- sample_pairs$Sample1_col
  sample2_columns <- sample_pairs$Sample2_col
  
} else {
  # For unpaired analysis
  print("Performing unpaired analysis between groups")
  
  sample1_columns <- group1_samples$filename
  sample2_columns <- group2_samples$filename
  
  print(paste("Group 1 samples:", length(sample1_columns)))
  print(paste("Group 2 samples:", length(sample2_columns)))
}

# Verify all sample columns exist in the metabolomics data
missing_sample1 <- sample1_columns[!sample1_columns %in% colnames(metabolomics_data)]
missing_sample2 <- sample2_columns[!sample2_columns %in% colnames(metabolomics_data)]

if(length(missing_sample1) > 0) {
  stop(paste("Missing Group 1 samples:", paste(missing_sample1, collapse = ", ")))
}
if(length(missing_sample2) > 0) {
  stop(paste("Missing Group 2 samples:", paste(missing_sample2, collapse = ", ")))
}

print("All sample columns found in metabolomics data. Proceeding with analysis...")

# Extract data for the two groups
group1_data <- metabolomics_data[, sample1_columns]
group2_data <- metabolomics_data[, sample2_columns]

# Ensure column order matches for paired analysis
if (paired_analysis) {
  print("Verifying paired data dimensions:")
  print(paste("Group 1 data dimensions:", nrow(group1_data), "x", ncol(group1_data)))
  print(paste("Group 2 data dimensions:", nrow(group2_data), "x", ncol(group2_data)))
  
  if (ncol(group1_data) != ncol(group2_data)) {
    stop("Unequal number of samples between groups for paired analysis!")
  }
}

# Function to perform Wilcoxon test for each metabolite
perform_wilcoxon_test <- function(group1_values, group2_values, paired = TRUE) {
  # Perform Wilcoxon test (paired or unpaired based on configuration)
  test_result <- wilcox.test(group1_values, group2_values, paired = paired, exact = FALSE)
  
  return(list(
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    n_samples = ifelse(paired, length(group1_values), length(group1_values) + length(group2_values))
  ))
}

print(paste("Performing", ifelse(paired_analysis, "paired", "unpaired"), "Wilcoxon tests..."))
print("This may take a few moments...")

# Perform the test for each metabolite (row)
results <- data.frame(
  row_ID = metabolomics_data$row.ID,
  mz = metabolomics_data$row.m.z,
  retention_time = metabolomics_data$row.retention.time,
  p_value = numeric(nrow(metabolomics_data)),
  statistic = numeric(nrow(metabolomics_data)),
  n_samples = rep(ifelse(paired_analysis, ncol(group1_data), ncol(group1_data) + ncol(group2_data)), nrow(metabolomics_data)),
  stringsAsFactors = FALSE
)

# Loop through each metabolite
for(i in 1:nrow(metabolomics_data)) {
  if(i %% 100 == 0) {
    print(paste("Processing metabolite", i, "of", nrow(metabolomics_data)))
  }
  
  group1_values <- as.numeric(group1_data[i, ])
  group2_values <- as.numeric(group2_data[i, ])
  
  test_result <- perform_wilcoxon_test(group1_values, group2_values, paired = paired_analysis)
  
  results$p_value[i] <- test_result$p_value
  results$statistic[i] <- test_result$statistic
}


# Calculate adjusted p-values using Benjamini-Hochberg (FDR) correction
results$p_adjusted_BH <- p.adjust(results$p_value, method = "BH")

# Calculate adjusted p-values using Bonferroni correction
results$p_adjusted_Bonferroni <- p.adjust(results$p_value, method = "bonferroni")

# Calculate effect size for Wilcoxon test
calculate_effect_size <- function(group1_values, group2_values, paired = TRUE) {
  if (paired) {
    # For paired analysis: rank-biserial correlation
    differences <- group1_values - group2_values
    n <- length(differences)
    
    test_result <- wilcox.test(group1_values, group2_values, paired = TRUE, exact = FALSE)
    W <- test_result$statistic
    r <- (2 * W) / (n * (n + 1)) - 0.5
    
    return(r)
  } else {
    # For unpaired analysis: rank-biserial correlation (Mann-Whitney U)
    n1 <- length(group1_values)
    n2 <- length(group2_values)
    
    test_result <- wilcox.test(group1_values, group2_values, paired = FALSE, exact = FALSE)
    U <- test_result$statistic
    r <- (2 * U) / (n1 * n2) - 1
    
    return(r)
  }
}

print("Calculating effect sizes...")
results$effect_size <- numeric(nrow(results))

for(i in 1:nrow(metabolomics_data)) {
  group1_values <- as.numeric(group1_data[i, ])
  group2_values <- as.numeric(group2_data[i, ])
  results$effect_size[i] <- calculate_effect_size(group1_values, group2_values, paired = paired_analysis)
}

# Add mean values for both groups for each metabolite
results$mean_group1 <- rowMeans(group1_data, na.rm = TRUE)
results$mean_group2 <- rowMeans(group2_data, na.rm = TRUE)
results$fold_change <- results$mean_group2 / results$mean_group1
results$log2_fold_change <- log2(results$fold_change)

# Sort results by p-value
results <- results[order(results$p_value), ]

# Summary statistics
print("Analysis completed!")
print("=== SUMMARY STATISTICS ===")
print(paste("Total metabolites analyzed:", nrow(results)))
print(paste("Metabolites with valid p-values:", sum(!is.na(results$p_value))))
print(paste("Significant metabolites (p < 0.05):", sum(results$p_value < 0.05, na.rm = TRUE)))
print(paste("Significant metabolites (FDR < 0.05):", sum(results$p_adjusted_BH < 0.05, na.rm = TRUE)))
print(paste("Significant metabolites (Bonferroni < 0.05):", sum(results$p_adjusted_Bonferroni < 0.05, na.rm = TRUE)))

# Display top 20 most significant results
print("=== TOP 20 MOST SIGNIFICANT METABOLITES ===")
top_results <- head(results, 20)
print(top_results[, c("row_ID", "mz", "retention_time", "p_value", "p_adjusted_BH", 
                      "effect_size", "mean_group1", "mean_group2", "log2_fold_change")])

# Create a more detailed view of significant metabolites
if(sum(results$p_adjusted_BH < 0.05, na.rm = TRUE) > 0) {
  print("=== SIGNIFICANT METABOLITES (FDR < 0.05) - DETAILED VIEW ===")
  significant_detailed <- results[results$p_adjusted_BH < 0.05 & !is.na(results$p_adjusted_BH), ]
  significant_detailed <- significant_detailed[order(significant_detailed$p_adjusted_BH), ]
  
  print(significant_detailed[, c("row_ID", "mz", "retention_time", "p_value", "p_adjusted_BH", 
                                 "effect_size", "fold_change", "log2_fold_change")])
}

# Save results to CSV
comparison_name <- paste0(group1, timepoint1, "_vs_", group2, timepoint2)
results_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_wilcoxon_results.csv")
results_file <- file.path(result_dir, results_filename)
write_csv(results, results_file)
print(paste("Results saved to:", results_file))

# Create volcano plot
print("Creating volcano plot...")
results$significance <- ifelse(results$p_adjusted_BH < 0.05, "Significant (FDR < 0.05)", "Not Significant")

volcano_plot <- ggplot(results, aes(x = log2_fold_change, y = -log10(p_value))) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(values = c("Significant (FDR < 0.05)" = "red", "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = paste("Volcano Plot:", comparison_name, ifelse(paired_analysis, "(Paired)", "(Unpaired)"), "Wilcoxon Test"),
    x = paste("Log2 Fold Change (", group2, timepoint2, "/", group1, timepoint1, ")", sep=""),
    y = "-Log10(P-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

volcano_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_volcano_plot.png")
volcano_plot_file <- file.path(result_dir, volcano_filename)
ggsave(volcano_plot_file, volcano_plot, width = 10, height = 8, dpi = 300)
print(paste("Volcano plot saved to:", volcano_plot_file))

# Create histogram of p-values
print("Creating p-value distribution plot...")
p_value_hist <- ggplot(results, aes(x = p_value)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(
    title = paste("P-value Distribution:", comparison_name, ifelse(paired_analysis, "(Paired)", "(Unpaired)"), "Wilcoxon Test"),
    x = "P-value",
    y = "Frequency"
  ) +
  theme_minimal()

pvalue_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_pvalue_distribution.png")
pvalue_plot_file <- file.path(result_dir, pvalue_filename)
ggsave(pvalue_plot_file, p_value_hist, width = 8, height = 6, dpi = 300)
print(paste("P-value distribution plot saved to:", pvalue_plot_file))

# Generate detailed summary report
print("=== DETAILED RESULTS SUMMARY ===")
cat("\nAnalysis:", ifelse(paired_analysis, "Paired", "Unpaired"), "Wilcoxon Test\n")
cat("Comparison:", comparison_name, "\n")
if (paired_analysis) {
  cat("Number of paired subjects:", ncol(group1_data), "\n")
} else {
  cat("Group 1 sample size:", ncol(group1_data), "\n")
  cat("Group 2 sample size:", ncol(group2_data), "\n")
}
cat("Number of metabolites tested:", nrow(results), "\n")
cat("Significance threshold: p < 0.05\n")
cat("Multiple testing correction: Benjamini-Hochberg (FDR)\n\n")

cat("Results Summary:\n")
cat("- Metabolites with p < 0.05:", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("- Metabolites with FDR < 0.05:", sum(results$p_adjusted_BH < 0.05, na.rm = TRUE), "\n")
cat("- Metabolites with Bonferroni < 0.05:", sum(results$p_adjusted_Bonferroni < 0.05, na.rm = TRUE), "\n")

if(sum(results$p_adjusted_BH < 0.05, na.rm = TRUE) > 0) {
  cat("\nTop 10 Significant metabolites (FDR < 0.05) with full identification:\n")
  significant_metabolites <- results[results$p_adjusted_BH < 0.05 & !is.na(results$p_adjusted_BH), ]
  significant_metabolites <- significant_metabolites[order(significant_metabolites$p_adjusted_BH), ]
  
  for(i in 1:min(10, nrow(significant_metabolites))) {
    cat(sprintf("- Row ID %s | m/z %.4f | RT %.2f min | p = %.2e | FDR = %.2e | FC = %.2f\n",
                significant_metabolites$row_ID[i],
                significant_metabolites$mz[i],
                significant_metabolites$retention_time[i],
                significant_metabolites$p_value[i],
                significant_metabolites$p_adjusted_BH[i],
                significant_metabolites$fold_change[i]))
  }
  
  # Also save a focused significant results file
  significant_results_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_significant_metabolites.csv")
  significant_results_file <- file.path(result_dir, significant_results_filename)
  write_csv(significant_metabolites, significant_results_file)
  cat(paste("\nSignificant metabolites saved separately to:", significant_results_file, "\n"))
}

print("Analysis complete! Check the generated files for detailed results.")
print("Files generated:")
print(paste("1.", results_file, "- Complete results table with metabolite IDs"))
print(paste("2.", volcano_plot_file, "- Volcano plot visualization"))
print(paste("3.", pvalue_plot_file, "- P-value distribution plot"))
if(exists("significant_results_file")) {
  print(paste("4.", significant_results_file, "- Significant metabolites only"))
}