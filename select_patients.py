#!/usr/bin/env python3
"""
Select good candidate patients for variant analysis.

We want patients with:
1. PBMC/normal sample (for tumor-normal comparison)
2. Biopsy sample (tissue tumor)
3. Multiple plasma timepoints (for tracking resistance)
"""

import pandas as pd
import sys

def main():
    # Read parsed metadata
    df = pd.read_csv('SraRunTable_parsed.csv')

    print("=" * 70)
    print("Identifying Patients with Complete Sample Sets")
    print("=" * 70)
    print()

    # Group by patient
    patient_groups = df.groupby('patient_id')

    # Find patients with complete sets
    complete_patients = []

    for patient_id, group in patient_groups:
        if patient_id is None or pd.isna(patient_id):
            continue

        categories = set(group['sample_category'].dropna())
        n_samples = len(group)
        n_plasma = len(group[group['sample_category'] == 'plasma'])
        n_pbmc = len(group[group['sample_category'] == 'PBMC'])
        n_biopsy = len(group[group['sample_category'] == 'biopsy'])

        # Check if patient has PBMC, biopsy, and at least 2 plasma timepoints
        has_normal = n_pbmc > 0
        has_biopsy = n_biopsy > 0
        has_plasma = n_plasma >= 2

        score = 0
        if has_normal: score += 10
        if has_biopsy: score += 10
        score += n_plasma * 5  # More plasma = better

        complete_patients.append({
            'patient_id': patient_id,
            'total_samples': n_samples,
            'pbmc': n_pbmc,
            'biopsy': n_biopsy,
            'plasma': n_plasma,
            'categories': ', '.join(sorted(categories)),
            'has_complete_set': has_normal and has_biopsy and has_plasma,
            'score': score
        })

    # Convert to DataFrame and sort by score
    summary_df = pd.DataFrame(complete_patients)
    summary_df = summary_df.sort_values('score', ascending=False)

    print("Top 20 patients ranked by completeness:")
    print("-" * 70)
    print(summary_df.head(20).to_string(index=False))
    print()

    # Show patients with complete sets
    complete = summary_df[summary_df['has_complete_set'] == True]

    print()
    print("=" * 70)
    print(f"Patients with COMPLETE sets (PBMC + Biopsy + 2+ Plasma): {len(complete)}")
    print("=" * 70)
    print()

    if len(complete) > 0:
        print(complete.to_string(index=False))
        print()

        # Create sample lists for top 5 complete patients
        print()
        print("=" * 70)
        print("Creating sample CSV files for top patients...")
        print("=" * 70)
        print()

        for idx, row in complete.head(5).iterrows():
            patient = row['patient_id']
            patient_data = df[df['patient_id'] == patient].copy()

            # Create CSV with tumor-normal pairs
            output_rows = []

            # Get normal SRR
            normal_samples = patient_data[patient_data['sample_category'] == 'PBMC']
            if len(normal_samples) == 0:
                print(f"Warning: {patient} has no PBMC sample, skipping")
                continue
            normal_srr = normal_samples.iloc[0]['Run']

            # Get biopsy
            biopsy_samples = patient_data[patient_data['sample_category'] == 'biopsy']
            if len(biopsy_samples) > 0:
                for _, sample in biopsy_samples.iterrows():
                    output_rows.append({
                        'patient': patient,
                        'tumor_srr': sample['Run'],
                        'normal_srr': normal_srr,
                        'sample_type': 'biopsy'
                    })

            # Get plasma samples
            plasma_samples = patient_data[patient_data['sample_category'] == 'plasma'].sort_values('timepoint')
            for idx2, sample in plasma_samples.iterrows():
                tp = sample['timepoint']
                if tp == 1:
                    sample_type = 'plasma_baseline'
                else:
                    sample_type = f'plasma_tp{tp}'

                output_rows.append({
                    'patient': patient,
                    'tumor_srr': sample['Run'],
                    'normal_srr': normal_srr,
                    'sample_type': sample_type
                })

            # Save to CSV
            output_df = pd.DataFrame(output_rows)
            filename = f"sample_lists/{patient}_samples.csv"
            output_df.to_csv(filename, index=False)
            print(f"Created: {filename} ({len(output_rows)} comparisons)")

    else:
        print("No patients found with complete sets.")
        print("Showing patients with at least PBMC + plasma:")
        partial = summary_df[(summary_df['pbmc'] > 0) & (summary_df['plasma'] >= 2)]
        print(partial.head(10).to_string(index=False))

if __name__ == '__main__':
    main()
