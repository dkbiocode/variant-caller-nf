#!/usr/bin/env Rscript
# Compare VarScan vs LoFreq somatic variant calls
# VarScan: results/varscan/all_variants_summary.tsv
# LoFreq:  results/lofreq/lofreq/all_lofreq_variants.tsv

library(ggplot2)
library(dplyr)
library(tidyr)

RESISTANCE_GENES <- c('KRAS', 'MAP2K1', 'TP53', 'APC', 'PTEN', 'BRAF', 'NRAS', 'PIK3CA')

# =============================================================================
# Load data
# =============================================================================
varscan <- read.delim("results/varscan/all_variants_summary.tsv",
                      stringsAsFactors = FALSE) %>%
  filter(somatic_status == "Somatic") %>%
  mutate(caller = "VarScan",
         vaf = as.numeric(vaf))

lofreq <- read.delim("results/lofreq/lofreq/all_lofreq_variants.tsv",
                     stringsAsFactors = FALSE) %>%
  mutate(caller = "LoFreq",
         vaf = as.numeric(vaf))

cat(sprintf("VarScan somatic variants: %d\n", nrow(varscan)))
cat(sprintf("LoFreq variants:          %d\n", nrow(lofreq)))

# Combine for joint plots (align common columns)
common_cols <- c("patient", "sample_type", "chr", "pos", "ref", "alt",
                 "vaf", "gene", "effect", "impact", "protein_change", "caller")
combined <- bind_rows(
  varscan %>% select(any_of(common_cols)),
  lofreq  %>% select(any_of(common_cols))
) %>%
  mutate(
    tissue_type = ifelse(grepl("biopsy", sample_type), "Biopsy", "Plasma"),
    timepoint = case_when(
      grepl("biopsy",    sample_type) ~ "Biopsy",
      grepl("baseline",  sample_type) ~ "Baseline",
      grepl("tp2",       sample_type) ~ "TP2",
      grepl("tp3",       sample_type) ~ "TP3",
      grepl("tp4",       sample_type) ~ "TP4",
      grepl("tp5",       sample_type) ~ "TP5",
      grepl("tp6",       sample_type) ~ "TP6",
      TRUE ~ "Other"
    )
  )

# =============================================================================
# Plot 1: Variant counts per caller, split by tissue type
# =============================================================================
count_summary <- combined %>%
  group_by(caller, tissue_type) %>%
  summarise(n_variants = n(), .groups = "drop")

p1 <- ggplot(count_summary, aes(x = caller, y = n_variants, fill = tissue_type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Biopsy" = "#E41A1C", "Plasma" = "#377EB8")) +
  labs(
    title = "Somatic Variant Counts by Caller",
    x = "Caller",
    y = "Number of Variants",
    fill = "Tissue Type"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("caller_variant_counts.png", p1, width = 7, height = 5, dpi = 300)
cat("Saved: caller_variant_counts.png\n")

# =============================================================================
# Plot 2: VAF distributions per caller (density), biopsy vs plasma
# =============================================================================
p2 <- combined %>%
  filter(!is.na(vaf)) %>%
  ggplot(aes(x = vaf, fill = caller, color = caller)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~tissue_type, scales = "free_y") +
  scale_fill_manual(values  = c("VarScan" = "#FF7F00", "LoFreq" = "#4DAF4A")) +
  scale_color_manual(values = c("VarScan" = "#FF7F00", "LoFreq" = "#4DAF4A")) +
  labs(
    title = "VAF Distribution by Caller and Tissue Type",
    x = "Variant Allele Frequency (%)",
    y = "Density",
    fill = "Caller", color = "Caller"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("caller_vaf_distribution.png", p2, width = 10, height = 5, dpi = 300)
cat("Saved: caller_vaf_distribution.png\n")

# =============================================================================
# Plot 3: Low-VAF variant detection (VAF < 5%)
# LoFreq is designed for low-VAF; check if it finds more low-VAF variants
# =============================================================================
low_vaf_summary <- combined %>%
  filter(!is.na(vaf)) %>%
  mutate(vaf_tier = case_when(
    vaf < 1   ~ "<1%",
    vaf < 2   ~ "1-2%",
    vaf < 5   ~ "2-5%",
    vaf < 10  ~ "5-10%",
    vaf < 20  ~ "10-20%",
    TRUE      ~ ">=20%"
  ),
  vaf_tier = factor(vaf_tier,
                    levels = c("<1%","1-2%","2-5%","5-10%","10-20%",">=20%"))) %>%
  group_by(caller, vaf_tier, tissue_type) %>%
  summarise(n_variants = n(), .groups = "drop")

p3 <- ggplot(low_vaf_summary,
             aes(x = vaf_tier, y = n_variants, fill = caller)) +
  geom_col(position = "dodge") +
  facet_wrap(~tissue_type) +
  scale_fill_manual(values = c("VarScan" = "#FF7F00", "LoFreq" = "#4DAF4A")) +
  labs(
    title = "Variant Counts by VAF Tier",
    x = "VAF Tier",
    y = "Number of Variants",
    fill = "Caller"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("caller_vaf_tiers.png", p3, width = 10, height = 5, dpi = 300)
cat("Saved: caller_vaf_tiers.png\n")

# =============================================================================
# Plot 4: Per-sample variant counts — caller comparison
# =============================================================================
per_sample <- combined %>%
  group_by(patient, sample_type, caller, tissue_type) %>%
  summarise(n_variants = n(), .groups = "drop")

p4 <- ggplot(per_sample, aes(x = caller, y = n_variants, fill = caller)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
  facet_wrap(~tissue_type) +
  scale_fill_manual(values = c("VarScan" = "#FF7F00", "LoFreq" = "#4DAF4A")) +
  labs(
    title = "Per-Sample Variant Count Distribution by Caller",
    x = "Caller",
    y = "Variants per Sample",
    fill = "Caller"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none")

ggsave("caller_per_sample_counts.png", p4, width = 8, height = 5, dpi = 300)
cat("Saved: caller_per_sample_counts.png\n")

# =============================================================================
# Plot 5: Resistance gene variants detected by each caller
# =============================================================================
resist_summary <- combined %>%
  filter(!is.na(gene), gene %in% RESISTANCE_GENES) %>%
  group_by(gene, caller, tissue_type) %>%
  summarise(n_variants = n(), .groups = "drop")

p5 <- ggplot(resist_summary,
             aes(x = gene, y = n_variants, fill = caller)) +
  geom_col(position = "dodge") +
  facet_wrap(~tissue_type) +
  scale_fill_manual(values = c("VarScan" = "#FF7F00", "LoFreq" = "#4DAF4A")) +
  labs(
    title = "Resistance Gene Variant Counts by Caller",
    x = "Gene",
    y = "Number of Variants",
    fill = "Caller"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("caller_resistance_genes.png", p5, width = 10, height = 5, dpi = 300)
cat("Saved: caller_resistance_genes.png\n")

# =============================================================================
# Plot 6: Concordance — variant overlap at chr:pos level per sample
# Find variants called by both callers vs caller-unique
# =============================================================================
# Create a key for matching: patient + sample_type + chr + pos + ref + alt
varscan_keys <- varscan %>%
  mutate(key = paste(patient, sample_type, chr, pos, ref, alt, sep = ":")) %>%
  pull(key)

lofreq_keys <- lofreq %>%
  mutate(key = paste(patient, sample_type, chr, pos, ref, alt, sep = ":")) %>%
  pull(key)

n_shared      <- length(intersect(varscan_keys, lofreq_keys))
n_varscan_only <- length(setdiff(varscan_keys, lofreq_keys))
n_lofreq_only  <- length(setdiff(lofreq_keys, varscan_keys))

cat(sprintf("\n=== Variant Concordance ===\n"))
cat(sprintf("Shared (both callers):  %d\n", n_shared))
cat(sprintf("VarScan only:           %d\n", n_varscan_only))
cat(sprintf("LoFreq only:            %d\n", n_lofreq_only))
cat(sprintf("VarScan Jaccard:        %.3f\n",
            n_shared / (n_shared + n_varscan_only + n_lofreq_only)))

concordance_df <- data.frame(
  category = c("Both callers", "VarScan only", "LoFreq only"),
  count    = c(n_shared, n_varscan_only, n_lofreq_only)
)

p6 <- ggplot(concordance_df, aes(x = category, y = count, fill = category)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Both callers" = "#984EA3",
    "VarScan only" = "#FF7F00",
    "LoFreq only"  = "#4DAF4A"
  )) +
  labs(
    title = "Variant Caller Concordance (chr:pos:ref:alt)",
    x = NULL,
    y = "Number of Variants"
  ) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none")

ggsave("caller_concordance.png", p6, width = 6, height = 5, dpi = 300)
cat("Saved: caller_concordance.png\n")

# =============================================================================
# Summary table
# =============================================================================
cat("\n=== Summary Statistics by Caller ===\n")
summary_stats <- combined %>%
  filter(!is.na(vaf)) %>%
  group_by(caller, tissue_type) %>%
  summarise(
    n_variants   = n(),
    median_vaf   = median(vaf),
    mean_vaf     = mean(vaf),
    pct_below_5  = mean(vaf < 5) * 100,
    pct_below_10 = mean(vaf < 10) * 100,
    .groups = "drop"
  )
print(summary_stats)
write.table(summary_stats, "caller_comparison_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved: caller_comparison_summary.tsv\n")
