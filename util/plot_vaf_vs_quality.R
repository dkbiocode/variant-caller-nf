#!/usr/bin/env Rscript
# Plot VAF vs Somatic Score for different sample groups

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the variant data
variants <- read.delim("results/varscan/all_variants_summary.tsv", stringsAsFactors = FALSE)

# Add sample groupings
variants <- variants %>%
  mutate(
    # Create tissue type (biopsy vs plasma)
    tissue_type = ifelse(grepl("biopsy", sample_type), "Biopsy", "Plasma"),
    # Create timepoint grouping
    timepoint = case_when(
      grepl("biopsy", sample_type) ~ "Biopsy",
      grepl("baseline", sample_type) ~ "Baseline",
      grepl("tp2", sample_type) ~ "TP2",
      grepl("tp3", sample_type) ~ "TP3",
      grepl("tp4", sample_type) ~ "TP4",
      grepl("tp5", sample_type) ~ "TP5",
      grepl("tp6", sample_type) ~ "TP6",
      TRUE ~ "Other"
    )
  )

# Filter for somatic variants only
somatic <- variants %>% filter(somatic_status == "Somatic")

# Convert somatic_score to numeric (handle any non-numeric values)
somatic$somatic_score <- as.numeric(somatic$somatic_score)

# Remove NA values
somatic <- somatic %>% filter(!is.na(vaf) & !is.na(somatic_score))

cat(sprintf("Total somatic variants: %d\n", nrow(somatic)))
cat(sprintf("Biopsy variants: %d\n", sum(somatic$tissue_type == "Biopsy")))
cat(sprintf("Plasma variants: %d\n", sum(somatic$tissue_type == "Plasma")))

# ============================================================================
# Plot 1: VAF vs Quality - All samples
# ============================================================================
p1 <- ggplot(somatic, aes(x = somatic_score, y = vaf)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(
    title = "VAF vs Somatic Score - All Samples",
    x = "Somatic Score (SSC)",
    y = "Variant Allele Frequency (%)"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("vaf_vs_quality_all.png", p1, width = 8, height = 6, dpi = 300)
cat("Saved: vaf_vs_quality_all.png\n")

# ============================================================================
# Plot 2: VAF vs Quality - Biopsy vs Plasma
# ============================================================================
p2 <- ggplot(somatic, aes(x = somatic_score, y = vaf, color = tissue_type)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("Biopsy" = "#E41A1C", "Plasma" = "#377EB8")) +
  labs(
    title = "VAF vs Somatic Score - Biopsy vs Plasma",
    x = "Somatic Score (SSC)",
    y = "Variant Allele Frequency (%)",
    color = "Tissue Type"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom"
  )

ggsave("vaf_vs_quality_tissue.png", p2, width = 8, height = 6, dpi = 300)
cat("Saved: vaf_vs_quality_tissue.png\n")

# ============================================================================
# Plot 3: VAF vs Quality - Faceted by timepoint
# ============================================================================
p3 <- ggplot(somatic, aes(x = somatic_score, y = vaf)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#377EB8") +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 0.8) +
  facet_wrap(~timepoint, scales = "free_y") +
  labs(
    title = "VAF vs Somatic Score by Timepoint",
    x = "Somatic Score (SSC)",
    y = "Variant Allele Frequency (%)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    strip.background = element_rect(fill = "lightgray")
  )

ggsave("vaf_vs_quality_timepoints.png", p3, width = 12, height = 8, dpi = 300)
cat("Saved: vaf_vs_quality_timepoints.png\n")

# ============================================================================
# Plot 4: Quality threshold analysis
# ============================================================================
# Add quality categories
somatic <- somatic %>%
  mutate(
    quality_tier = case_when(
      somatic_score < 20 ~ "Low (<20)",
      somatic_score < 40 ~ "Medium (20-40)",
      somatic_score < 100 ~ "High (40-100)",
      TRUE ~ "Very High (≥100)"
    ),
    quality_tier = factor(quality_tier,
                         levels = c("Low (<20)", "Medium (20-40)",
                                   "High (40-100)", "Very High (≥100)"))
  )

p4 <- ggplot(somatic, aes(x = quality_tier, y = vaf, fill = tissue_type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("Biopsy" = "#E41A1C", "Plasma" = "#377EB8")) +
  labs(
    title = "VAF Distribution by Quality Tier",
    x = "Somatic Score Tier",
    y = "Variant Allele Frequency (%)",
    fill = "Tissue Type"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("vaf_by_quality_tier.png", p4, width = 8, height = 6, dpi = 300)
cat("Saved: vaf_by_quality_tier.png\n")

# ============================================================================
# Summary statistics
# ============================================================================
cat("\n=== Summary Statistics ===\n")

summary_stats <- somatic %>%
  group_by(tissue_type, quality_tier) %>%
  summarise(
    n_variants = n(),
    mean_vaf = mean(vaf, na.rm = TRUE),
    median_vaf = median(vaf, na.rm = TRUE),
    mean_ssc = mean(somatic_score, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

write.table(summary_stats, "vaf_quality_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved: vaf_quality_summary.tsv\n")
