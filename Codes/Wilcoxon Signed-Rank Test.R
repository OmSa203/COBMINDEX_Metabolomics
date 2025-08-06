# Paired Wilcoxon Signed-Rank Test for TR1 vs TR2 Metabolomics Data
# Author: Analysis script for Crohn's disease metabolomics study
# Date: August 2025

# --- 1. Configuration ---
data_dir <- "/Users/batelhagbi/Documents/GitHub/COBMINDEX_Metabolomics/Data"
base_result_dir <- "/Users/batelhagbi/Documents/GitHub/COBMINDEX_Metabolomics/results/Univariate"

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
metadata <- read_csv("Negative_metadata_w_subjects.csv")
metabolomics_data <- read_csv("Imputed_TIC_Normalised_table_Negative.csv")

# Clean column names in metadata (remove spaces)
colnames(metadata) <- trimws(colnames(metadata))

# Automatically detect polarity from the loaded file names
# Extract polarity from the actual file paths used above
metadata_filename <- "Negative_metadata_w_subjects.csv"  # Update this when you change the file
data_filename <- "Imputed_TIC_Normalised_table_Negative.csv"  # Update this when you change the file

polarity <- "Unknown"  # Default

# Check for polarity indicators in the actual filenames being used
all_filenames <- c(metadata_filename, data_filename)

if (any(grepl("POS|Positive", all_filenames, ignore.case = TRUE))) {
  polarity <- "Positive"
} else if (any(grepl("NEG|Negative", all_filenames, ignore.case = TRUE))) {
  polarity <- "Negative"
}

print(paste("Detected polarity:", polarity))

# Create organized directory structure
comparison_name <- paste0(group1, timepoint1, "_vs_", group2, timepoint2)
analysis_type <- ifelse(paired_analysis, "paired", "unpaired")
comparison_folder <- paste0(comparison_name, "_", analysis_type)

# Create the full directory path: base/TR1_vs_TR2_paired/POS/
result_dir <- file.path(base_result_dir, comparison_folder, polarity)

# Create results directory structure if it doesn't exist
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
  print(paste("Created directory structure:", result_dir))
} else {
  print(paste("Using existing directory:", result_dir))
}

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
print(paste("Metabolites with large effect size (|effect| > 0.3):", sum(abs(results$effect_size) > 0.3, na.rm = TRUE)))

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

# Find metabolites with meaningful changes (good effect size + statistical support)
meaningful_changes <- results[abs(results$effect_size) > 0.3 & results$p_value < 0.05, ]

if(nrow(meaningful_changes) > 0) {
  cat("\n=== MEANINGFUL CHANGES (|Effect Size| > 0.3 & p < 0.05) ===\n")
  cat("Found", nrow(meaningful_changes), "metabolites with meaningful changes:\n")
  
  meaningful_changes <- meaningful_changes[order(meaningful_changes$p_value), ]
  
  for(i in 1:min(20, nrow(meaningful_changes))) {
    cat(sprintf("%d. Row ID %s | m/z %.4f | RT %.2f min | p = %.3f | Effect = %.3f | FC = %.2f\n",
                i,
                meaningful_changes$row_ID[i],
                meaningful_changes$mz[i],
                meaningful_changes$retention_time[i],
                meaningful_changes$p_value[i],
                meaningful_changes$effect_size[i],
                meaningful_changes$fold_change[i]))
  }
  
  # Save meaningful changes to a separate file
  meaningful_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_meaningful_changes.csv")
  meaningful_file <- file.path(result_dir, meaningful_filename)
  write_csv(meaningful_changes, meaningful_file)
  cat(paste("\nMeaningful changes saved to:", meaningful_file, "\n"))
  
  # Create subfolder for individual metabolite plots within the organized structure
  individual_plots_dir <- file.path(result_dir, "individual_metabolite_plots")
  if (!dir.exists(individual_plots_dir)) {
    dir.create(individual_plots_dir, recursive = TRUE)
  }
  print(paste("Created directory for individual plots:", individual_plots_dir))
  
  # Create individual plots for meaningful changes - ALL metabolites
  print("Creating individual plots for ALL meaningful changes...")
  
  # First create individual plots for each metabolite
  print("Creating separate plot files for each meaningful metabolite...")
  
  for (i in 1:nrow(meaningful_changes)) {
    metabolite_idx <- which(metabolomics_data$row.ID == meaningful_changes$row_ID[i])
    
    group1_vals <- as.numeric(group1_data[metabolite_idx, ])
    group2_vals <- as.numeric(group2_data[metabolite_idx, ])
    
    # Create data for this specific metabolite
    if (paired_analysis) {
      individual_data <- data.frame(
        Subject = rep(1:length(group1_vals), 2),
        Group = rep(c(paste0(group1, timepoint1), paste0(group2, timepoint2)), each = length(group1_vals)),
        Value = c(group1_vals, group2_vals),
        stringsAsFactors = FALSE
      )
    } else {
      individual_data <- data.frame(
        Subject = c(1:length(group1_vals), 1:length(group2_vals)),
        Group = c(rep(paste0(group1, timepoint1), length(group1_vals)), 
                  rep(paste0(group2, timepoint2), length(group2_vals))),
        Value = c(group1_vals, group2_vals),
        stringsAsFactors = FALSE
      )
    }
    
    # Create individual plot
    metabolite_title <- paste("m/z", round(meaningful_changes$mz[i], 4), "| RT", round(meaningful_changes$retention_time[i], 2), "min")
    metabolite_subtitle <- paste("Row ID:", meaningful_changes$row_ID[i], "| p =", round(meaningful_changes$p_value[i], 4), 
                                 "| Effect Size =", round(meaningful_changes$effect_size[i], 3), 
                                 "| Fold Change =", round(meaningful_changes$fold_change[i], 2))
    
    individual_plot <- ggplot(individual_data, aes(x = Group, y = Value)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, fill = "lightblue") +
      geom_point(aes(color = Group), size = 3, alpha = 0.8, position = position_jitter(width = 0.15)) +
      {if(paired_analysis) geom_line(aes(group = Subject), alpha = 0.6, color = "gray50", size = 0.8)} +
      labs(
        title = metabolite_title,
        subtitle = metabolite_subtitle,
        x = paste("Treatment Group (", comparison_name, ")", sep=""),
        y = "Normalized Intensity",
        color = "Group"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray30"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
      ) +
      scale_color_manual(values = c("steelblue", "darkred"))
    
    # Add significance annotation
    y_max <- max(individual_data$Value) * 1.1
    if (meaningful_changes$p_value[i] < 0.001) {
      sig_text <- "***"
    } else if (meaningful_changes$p_value[i] < 0.01) {
      sig_text <- "**"
    } else if (meaningful_changes$p_value[i] < 0.05) {
      sig_text <- "*"
    }
    
    individual_plot <- individual_plot +
      annotate("text", x = 1.5, y = y_max, label = sig_text, size = 6, color = "red")
    
    # Save individual plot
    safe_mz <- round(meaningful_changes$mz[i], 4)
    safe_rt <- round(meaningful_changes$retention_time[i], 2)
    individual_filename <- paste0("Row_", meaningful_changes$row_ID[i], "_mz_", safe_mz, "_RT_", safe_rt, ".png")
    individual_file_path <- file.path(individual_plots_dir, individual_filename)
    
    ggsave(individual_file_path, individual_plot, width = 8, height = 6, dpi = 300)
    
    if (i %% 10 == 0) {
      print(paste("Saved individual plot", i, "of", nrow(meaningful_changes)))
    }
  }
  
  cat(paste("Individual plots for all", nrow(meaningful_changes), "meaningful metabolites saved to:", individual_plots_dir, "\n"))
  
  # Now create the combined plot with all metabolites
  print("Creating combined plot with all meaningful metabolites...")
  
  # Prepare data for plotting ALL meaningful metabolites in combined plot
  plot_data <- data.frame()
  
  for (i in 1:nrow(meaningful_changes)) {  # Plot ALL meaningful metabolites
    metabolite_idx <- which(metabolomics_data$row.ID == meaningful_changes$row_ID[i])
    
    group1_vals <- as.numeric(group1_data[metabolite_idx, ])
    group2_vals <- as.numeric(group2_data[metabolite_idx, ])
    
    # Create paired data for plotting
    if (paired_analysis) {
      temp_data <- data.frame(
        Subject = rep(1:length(group1_vals), 2),
        Group = rep(c(paste0(group1, timepoint1), paste0(group2, timepoint2)), each = length(group1_vals)),
        Value = c(group1_vals, group2_vals),
        Metabolite = paste("m/z", round(meaningful_changes$mz[i], 4), 
                           "(p =", round(meaningful_changes$p_value[i], 3), 
                           ", Effect =", round(meaningful_changes$effect_size[i], 3), ")"),
        MetaboliteID = meaningful_changes$row_ID[i],
        stringsAsFactors = FALSE
      )
    } else {
      temp_data <- data.frame(
        Subject = c(1:length(group1_vals), 1:length(group2_vals)),
        Group = c(rep(paste0(group1, timepoint1), length(group1_vals)), 
                  rep(paste0(group2, timepoint2), length(group2_vals))),
        Value = c(group1_vals, group2_vals),
        Metabolite = paste("m/z", round(meaningful_changes$mz[i], 4), 
                           "(p =", round(meaningful_changes$p_value[i], 3), 
                           ", Effect =", round(meaningful_changes$effect_size[i], 3), ")"),
        MetaboliteID = meaningful_changes$row_ID[i],
        stringsAsFactors = FALSE
      )
    }
    
    plot_data <- rbind(plot_data, temp_data)
  }
  
  # Calculate optimal plot dimensions based on number of metabolites
  n_metabolites <- nrow(meaningful_changes)
  n_cols <- min(6, ceiling(sqrt(n_metabolites)))  # Max 6 columns
  n_rows <- ceiling(n_metabolites / n_cols)
  
  # Calculate plot dimensions
  plot_width <- max(15, n_cols * 2.5)
  plot_height <- max(10, n_rows * 2)
  
  print(paste("Creating plot with", n_metabolites, "metabolites in", n_rows, "rows and", n_cols, "columns"))
  
  # Create the plot
  meaningful_plot <- ggplot(plot_data, aes(x = Group, y = Value)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_point(aes(color = Group), size = 1.5, alpha = 0.8, position = position_jitter(width = 0.2)) +
    {if(paired_analysis) geom_line(aes(group = Subject), alpha = 0.5, color = "gray")} +
    facet_wrap(~ Metabolite, scales = "free_y", ncol = n_cols) +
    labs(
      title = paste("ALL Meaningful Metabolite Changes:", comparison_name),
      subtitle = paste("All", n_metabolites, "metabolites with Effect Size > 0.3 & p < 0.05 (", ifelse(paired_analysis, "Paired", "Unpaired"), "Analysis )"),
      x = "Group",
      y = "Normalized Intensity",
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 7),
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("steelblue", "darkred"))
  
  # Save the meaningful changes plot
  meaningful_plot_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_ALL_meaningful_changes_plot.png")
  meaningful_plot_file <- file.path(result_dir, meaningful_plot_filename)
  ggsave(meaningful_plot_file, meaningful_plot, width = plot_width, height = plot_height, dpi = 300)
  cat(paste("ALL meaningful changes plot (", n_metabolites, "metabolites) saved to:", meaningful_plot_file, "\n"))
  
  # Also create a focused plot with just the top 12 most significant
  if (n_metabolites > 12) {
    top12_data <- plot_data[plot_data$MetaboliteID %in% meaningful_changes$row_ID[1:12], ]
    
    top12_plot <- ggplot(top12_data, aes(x = Group, y = Value)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_point(aes(color = Group), size = 2, alpha = 0.8) +
      {if(paired_analysis) geom_line(aes(group = Subject), alpha = 0.5, color = "gray")} +
      facet_wrap(~ Metabolite, scales = "free_y", ncol = 3) +
      labs(
        title = paste("Top 12 Most Significant Metabolites:", comparison_name),
        subtitle = paste("Lowest p-values among meaningful changes (Effect Size > 0.3 & p < 0.05)"),
        x = "Group",
        y = "Normalized Intensity",
        color = "Group"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        legend.position = "bottom"
      ) +
      scale_color_manual(values = c("steelblue", "darkred"))
    
    top12_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_top12_meaningful_changes.png")
    top12_file <- file.path(result_dir, top12_filename)
    ggsave(top12_file, top12_plot, width = 15, height = 12, dpi = 300)
    cat(paste("Top 12 meaningful changes plot saved to:", top12_file, "\n"))
  }
  
} else {
  cat("\nNo metabolites found with meaningful changes (|Effect Size| > 0.3 & p < 0.05)\n")
  cat("Consider relaxing criteria or trying different group comparisons.\n")
}

# Create effect size vs p-value plot for all results
print("Creating overall effect size visualization...")

results$meaningful <- ifelse(abs(results$effect_size) > 0.3 & results$p_value < 0.05, 
                             "Meaningful Change", "Other")

effect_plot <- ggplot(results, aes(x = effect_size, y = -log10(p_value))) +
  geom_point(aes(color = meaningful), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Meaningful Change" = "red", "Other" = "gray70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "darkgreen", alpha = 0.7) +
  labs(
    title = paste("Effect Size vs Statistical Significance:", comparison_name),
    subtitle = paste("Meaningful changes:", sum(abs(results$effect_size) > 0.3 & results$p_value < 0.05, na.rm = TRUE), "metabolites with |Effect Size| > 0.3 & p < 0.05"),
    x = "Effect Size (Rank-biserial correlation)",
    y = "-Log10(P-value)",
    color = "Classification"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  annotate("text", x = 0.4, y = max(-log10(results$p_value)) * 0.9, 
           label = "Large Effect\n& Significant", size = 3, color = "darkgreen") +
  annotate("text", x = -0.4, y = max(-log10(results$p_value)) * 0.9, 
           label = "Large Effect\n& Significant", size = 3, color = "darkgreen")

effect_plot_filename <- paste0(comparison_name, "_", ifelse(paired_analysis, "paired", "unpaired"), "_effect_size_plot.png")
effect_plot_file <- file.path(result_dir, effect_plot_filename)
ggsave(effect_plot_file, effect_plot, width = 10, height = 8, dpi = 300)
print(paste("Effect size plot saved to:", effect_plot_file))

# Handle FDR significant metabolites if any exist
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

# Final summary
print("Analysis complete! Check the generated files for detailed results.")
print(paste("All results saved in organized directory structure:"))
print(paste("Main directory:", result_dir))
print("")
print("Files generated:")
print(paste("1.", results_file, "- Complete results table with metabolite IDs"))
print(paste("2.", volcano_plot_file, "- Volcano plot visualization"))
print(paste("3.", pvalue_plot_file, "- P-value distribution plot"))
if(exists("meaningful_file")) {
  print(paste("4.", meaningful_file, "- Meaningful changes (Effect > 0.3 & p < 0.05)"))
}
if(exists("individual_plots_dir")) {
  print(paste("5.", individual_plots_dir, "/ - Individual plots for each meaningful metabolite"))
}
if(exists("meaningful_plot_file")) {
  print(paste("6.", meaningful_plot_file, "- ALL meaningful metabolites combined plot"))
}
if(exists("top12_file")) {
  print(paste("7.", top12_file, "- Top 12 meaningful metabolites plot"))
}
if(exists("effect_plot_file")) {
  print(paste("8.", effect_plot_file, "- Effect size vs p-value visualization"))
}
if(exists("significant_results_file")) {
  print(paste("9.", significant_results_file, "- FDR significant metabolites only"))
}

print("")
print(paste("Directory structure created:"))
print(paste("└── ", base_result_dir))
print(paste("    └── ", comparison_folder))
print(paste("        └── ", polarity))
print(paste("            ├── Results files (.csv, .png)"))
print(paste("            └── individual_metabolite_plots/"))

