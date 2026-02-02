#!/usr/bin/env python3
"""
Extract VarScan somatic scores and create MultiQC custom content for visualization.

VarScan uses SSC (Somatic Score) field instead of standard QUAL field.
This script extracts SSC values from the all_variants_summary.tsv file
and creates a MultiQC-compatible custom content file.
"""

import sys
import pandas as pd
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: extract_somatic_scores.py <all_variants_summary.tsv> [<output_multiqc.yaml>]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "somatic_scores_mqc.yaml"

    # Read the variant summary
    try:
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    # Filter for somatic variants only
    somatic_df = df[df['somatic_status'] == 'Somatic'].copy()

    if len(somatic_df) == 0:
        print("Warning: No somatic variants found in input file", file=sys.stderr)
        # Create empty output file for MultiQC
        with open(output_file, 'w') as f:
            f.write("# MultiQC custom content: VarScan Somatic Scores\n")
            f.write("# No somatic variants found\n")
        return

    # Convert somatic_score to numeric (VarScan SSC ranges 0-255)
    somatic_df['somatic_score'] = pd.to_numeric(somatic_df['somatic_score'], errors='coerce')

    # Remove NA values
    somatic_df = somatic_df.dropna(subset=['somatic_score'])

    if len(somatic_df) == 0:
        print("Warning: No valid somatic scores found", file=sys.stderr)
        with open(output_file, 'w') as f:
            f.write("# MultiQC custom content: VarScan Somatic Scores\n")
            f.write("# No valid somatic scores found\n")
        return

    # Create binned data for histogram (MultiQC format)
    # Bin SSC scores into 20 bins from 0-255
    bins = range(0, 260, 13)  # Creates bins: 0-13, 13-26, ..., 247-260
    somatic_df['ssc_bin'] = pd.cut(somatic_df['somatic_score'], bins=bins, right=False)

    # Count variants per bin
    bin_counts = somatic_df['ssc_bin'].value_counts().sort_index()

    # Also get per-sample statistics
    sample_stats = somatic_df.groupby(['patient', 'sample_type']).agg({
        'somatic_score': ['count', 'mean', 'median', 'min', 'max']
    }).reset_index()
    sample_stats.columns = ['patient', 'sample_type', 'variant_count', 'mean_ssc', 'median_ssc', 'min_ssc', 'max_ssc']
    sample_stats['sample'] = sample_stats['patient'] + '_' + sample_stats['sample_type']

    # Write MultiQC custom content YAML
    with open(output_file, 'w') as f:
        f.write("# MultiQC custom content file for VarScan Somatic Scores\n")
        f.write("# This file provides custom visualizations for SSC (Somatic Score) distribution\n\n")

        # Section 1: Overall somatic score distribution (histogram)
        f.write("id: 'varscan_somatic_scores'\n")
        f.write("section_name: 'VarScan Somatic Scores'\n")
        f.write("description: 'Distribution of VarScan SSC (Somatic Score) values for all somatic variants. VarScan uses SSC (0-255) instead of QUAL field.'\n")
        f.write("plot_type: 'bargraph'\n")
        f.write("pconfig:\n")
        f.write("    id: 'varscan_somatic_score_dist'\n")
        f.write("    title: 'VarScan Somatic Score Distribution'\n")
        f.write("    xlab: 'Somatic Score (SSC)'\n")
        f.write("    ylab: 'Number of Variants'\n")
        f.write("data:\n")
        f.write("    'All Somatic Variants':\n")
        for interval, count in bin_counts.items():
            # Format bin label as "min-max"
            bin_label = f"{int(interval.left)}-{int(interval.right)}"
            f.write(f"        '{bin_label}': {count}\n")

        f.write("\n---\n\n")

        # Section 2: Per-sample statistics (table)
        f.write("id: 'varscan_somatic_scores_table'\n")
        f.write("section_name: 'VarScan Somatic Scores by Sample'\n")
        f.write("description: 'Summary statistics of VarScan SSC values per sample'\n")
        f.write("plot_type: 'table'\n")
        f.write("pconfig:\n")
        f.write("    id: 'varscan_ssc_table'\n")
        f.write("    title: 'VarScan SSC Statistics'\n")
        f.write("    col1_header: 'Sample'\n")
        f.write("headers:\n")
        f.write("    variant_count:\n")
        f.write("        title: 'Variants'\n")
        f.write("        description: 'Number of somatic variants'\n")
        f.write("        format: '{:,.0f}'\n")
        f.write("    mean_ssc:\n")
        f.write("        title: 'Mean SSC'\n")
        f.write("        description: 'Mean somatic score'\n")
        f.write("        format: '{:,.1f}'\n")
        f.write("    median_ssc:\n")
        f.write("        title: 'Median SSC'\n")
        f.write("        description: 'Median somatic score'\n")
        f.write("        format: '{:,.1f}'\n")
        f.write("    min_ssc:\n")
        f.write("        title: 'Min SSC'\n")
        f.write("        description: 'Minimum somatic score'\n")
        f.write("        format: '{:,.0f}'\n")
        f.write("    max_ssc:\n")
        f.write("        title: 'Max SSC'\n")
        f.write("        description: 'Maximum somatic score'\n")
        f.write("        format: '{:,.0f}'\n")
        f.write("data:\n")
        for _, row in sample_stats.iterrows():
            f.write(f"    '{row['sample']}':\n")
            f.write(f"        variant_count: {row['variant_count']}\n")
            f.write(f"        mean_ssc: {row['mean_ssc']:.1f}\n")
            f.write(f"        median_ssc: {row['median_ssc']:.1f}\n")
            f.write(f"        min_ssc: {row['min_ssc']:.0f}\n")
            f.write(f"        max_ssc: {row['max_ssc']:.0f}\n")

    # Print summary
    print(f"\n=== VarScan Somatic Score Summary ===")
    print(f"Total somatic variants: {len(somatic_df)}")
    print(f"SSC range: {somatic_df['somatic_score'].min():.0f} - {somatic_df['somatic_score'].max():.0f}")
    print(f"Mean SSC: {somatic_df['somatic_score'].mean():.1f}")
    print(f"Median SSC: {somatic_df['somatic_score'].median():.1f}")
    print(f"\nCreated MultiQC custom content: {output_file}")

if __name__ == '__main__':
    main()
