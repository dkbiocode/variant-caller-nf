#!/usr/bin/env python3
"""
Parse SRA metadata to extract patient IDs and timepoints from sample names.

The sample naming appears to follow these patterns:
- CTC### or CTC###-# (plasma samples with patient ID and timepoint)
- PBMC_CTC### or PBMC_C_fw### (normal blood controls)
- FFPE_C_fw### (tissue biopsy samples)
- C_fw###-# (plasma samples with patient ID and timepoint)

This script extracts:
1. Base patient ID (e.g., CTC362, C_fw012, etc.)
2. Timepoint (e.g., baseline=1, second=-2, etc.)
"""

import pandas as pd
import re
import sys

def parse_sample_name(sample_name):
    """
    Extract patient ID and timepoint from sample name.

    Returns: (patient_id, timepoint, sample_category)
    """

    # Pattern 1: CTC###-# (patient with timepoint)
    match = re.match(r'(CTC\d+)-(\d+)', sample_name)
    if match:
        return match.group(1), int(match.group(2)), 'plasma'

    # Pattern 2: CTC### (patient baseline)
    match = re.match(r'(CTC\d+)$', sample_name)
    if match:
        return match.group(1), 1, 'plasma'

    # Pattern 3: PBMC_CTC### or PBMC_C_fw### (normal blood)
    match = re.match(r'PBMC_(CTC\d+|C_fw\d+)', sample_name)
    if match:
        return match.group(1), 0, 'PBMC'

    # Pattern 4: FFPE_C_fw### (tissue biopsy)
    match = re.match(r'FFPE_(C_fw\d+)', sample_name)
    if match:
        return match.group(1), 0, 'biopsy'

    # Pattern 5: C_fw###-# (patient with timepoint)
    match = re.match(r'(C_fw\d+)-(\d+)', sample_name)
    if match:
        return match.group(1), int(match.group(2)), 'plasma'

    # Pattern 6: C_fw### (patient baseline)
    match = re.match(r'(C_fw\d+)$', sample_name)
    if match:
        return match.group(1), 1, 'plasma'

    # Pattern 7: PBMC_CTC### for PBMC samples
    match = re.match(r'PBMC_CTC(\d+)', sample_name)
    if match:
        return f'CTC{match.group(1)}', 0, 'PBMC'

    # Unknown pattern
    return None, None, 'unknown'

def main():
    # Read the SRA metadata
    input_file = 'SraRunTable_all.csv'
    output_file = 'SraRunTable_parsed.csv'

    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file)

    print(f"Found {len(df)} samples")

    # Parse sample names
    print("Parsing sample names...")
    parsed = df['Sample Name'].apply(parse_sample_name)

    # Add new columns
    df['patient_id'] = [x[0] for x in parsed]
    df['timepoint'] = [x[1] for x in parsed]
    df['sample_category'] = [x[2] for x in parsed]

    # Reorder columns to put new ones near the front
    cols = list(df.columns)
    # Move Run, Sample Name, tissue, patient_id, timepoint, sample_category to front
    priority_cols = ['Run', 'Sample Name', 'tissue', 'patient_id', 'timepoint', 'sample_category']
    other_cols = [c for c in cols if c not in priority_cols]
    df = df[priority_cols + other_cols]

    # Save parsed version
    print(f"Writing to {output_file}...")
    df.to_csv(output_file, index=False)

    # Print summary statistics
    print("\n=== Summary ===")
    print(f"Total samples: {len(df)}")
    print(f"\nSample categories:")
    print(df['sample_category'].value_counts())

    print(f"\nSamples per patient (top 20):")
    patient_counts = df['patient_id'].value_counts().head(20)
    print(patient_counts)

    print(f"\nTimepoint distribution:")
    print(df['timepoint'].value_counts().sort_index())

    # Find patients with multiple sample types
    print("\n=== Patients with complete sets (PBMC + biopsy + plasma) ===")
    patient_groups = df.groupby('patient_id')

    complete_patients = []
    for patient_id, group in patient_groups:
        if patient_id is None:
            continue
        categories = set(group['sample_category'])
        if 'PBMC' in categories and 'biopsy' in categories and 'plasma' in categories:
            complete_patients.append({
                'patient_id': patient_id,
                'n_samples': len(group),
                'categories': ', '.join(sorted(categories)),
                'timepoints': ', '.join(map(str, sorted(group['timepoint'].unique())))
            })

    if complete_patients:
        complete_df = pd.DataFrame(complete_patients)
        print(f"\nFound {len(complete_df)} patients with PBMC + biopsy + plasma:")
        print(complete_df.to_string(index=False))
    else:
        print("No patients found with all three sample types (PBMC, biopsy, plasma)")

    print(f"\n✓ Parsed data saved to {output_file}")

if __name__ == '__main__':
    main()
