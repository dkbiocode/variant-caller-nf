#!/usr/bin/env python3
"""
Combine multiple patient sample lists into a single CSV for batch processing.
"""

import pandas as pd
import sys
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: combine_patient_lists.py <num_patients> [output.csv]")
        print("  num_patients: Number of top patients to include (max 22)")
        print("  output.csv: Output filename (default: multi_patient_samples.csv)")
        sys.exit(1)

    num_patients = int(sys.argv[1])
    output_file = sys.argv[2] if len(sys.argv) > 2 else "sample_lists/multi_patient_samples.csv"

    # Get available patient CSV files
    patient_files = sorted(Path('sample_lists').glob('C_fw*_samples.csv'))

    if len(patient_files) == 0:
        print("Error: No patient CSV files found. Run select_patients.py first.")
        sys.exit(1)

    # Limit to requested number
    patient_files = patient_files[:num_patients]

    print(f"Combining {len(patient_files)} patients:")
    print("-" * 50)

    all_data = []
    total_comparisons = 0

    for csv_file in patient_files:
        df = pd.read_csv(csv_file)
        patient = df['patient'].iloc[0]
        n_comparisons = len(df)
        total_comparisons += n_comparisons

        print(f"  {patient}: {n_comparisons} comparisons")
        all_data.append(df)

    # Combine all dataframes
    combined = pd.concat(all_data, ignore_index=True)

    # Save
    combined.to_csv(output_file, index=False)

    print("-" * 50)
    print(f"Combined CSV saved: {output_file}")
    print(f"Total patients: {len(patient_files)}")
    print(f"Total comparisons: {total_comparisons}")
    print()
    print("To run pipeline with this CSV:")
    print(f"  nextflow run main.nf --samples_csv {output_file}")

if __name__ == '__main__':
    main()
