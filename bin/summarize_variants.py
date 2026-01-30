#!/usr/bin/env python3
"""
Extract variant information from VarScan VCF files and create a summary table
similar to Table 4 from the publication.

This script focuses on key resistance genes: KRAS, MAP2K1, TP53, APC, PTEN
"""

import sys
import pandas as pd
from pathlib import Path

# Key genes from Table 4 (resistance mutations)
RESISTANCE_GENES = ['KRAS', 'MAP2K1', 'TP53', 'APC', 'PTEN', 'BRAF', 'NRAS', 'PIK3CA']

def parse_vcf_line(line, patient, sample_type):
    """Parse a VCF data line and extract key information."""
    fields = line.strip().split('\t')

    if len(fields) < 10:
        return None

    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    info = fields[7]
    format_fields = fields[8].split(':')
    tumor_data = fields[10].split(':')

    # Parse INFO field
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value

    # Extract variant allele frequency from tumor sample
    vaf = None
    if 'FREQ' in format_fields:
        freq_idx = format_fields.index('FREQ')
        if freq_idx < len(tumor_data):
            freq_str = tumor_data[freq_idx].rstrip('%')
            try:
                vaf = float(freq_str)
            except ValueError:
                vaf = None

    # Parse somatic p-value
    ss_score = info_dict.get('SS', '')
    somatic_status = {
        '1': 'Germline',
        '2': 'Somatic',
        '3': 'LOH',
        '5': 'Unknown'
    }.get(ss_score, 'Unknown')

    # Get somatic p-value
    ssc = info_dict.get('SSC', 'NA')

    variant = {
        'patient': patient,
        'sample_type': sample_type,
        'chr': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'vaf': vaf,
        'somatic_status': somatic_status,
        'somatic_score': ssc,
        'info': info_dict
    }

    return variant

def parse_vcf_file(vcf_path, patient, sample_type):
    """Parse a VCF file and extract variants."""
    variants = []

    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                variant = parse_vcf_line(line, patient, sample_type)
                if variant:
                    variants.append(variant)
    except FileNotFoundError:
        print(f"Warning: VCF file not found: {vcf_path}", file=sys.stderr)

    return variants

def main():
    if len(sys.argv) < 4:
        print("Usage: summarize_variants.py <patient> <sample_type> <snp.vcf> <indel.vcf> [<output.tsv>]")
        sys.exit(1)

    patient = sys.argv[1]
    sample_type = sys.argv[2]
    snp_vcf = sys.argv[3]
    indel_vcf = sys.argv[4]
    output_file = sys.argv[5] if len(sys.argv) > 5 else f"{patient}_{sample_type}_variants.tsv"

    # Parse both SNP and INDEL VCF files
    all_variants = []
    all_variants.extend(parse_vcf_file(snp_vcf, patient, sample_type))
    all_variants.extend(parse_vcf_file(indel_vcf, patient, sample_type))

    if not all_variants:
        # Create empty output file
        pd.DataFrame(columns=[
            'patient', 'sample_type', 'chr', 'pos', 'ref', 'alt',
            'vaf', 'somatic_status', 'somatic_score'
        ]).to_csv(output_file, sep='\t', index=False)
        print(f"No variants found. Created empty file: {output_file}")
        return

    # Convert to DataFrame
    df = pd.DataFrame(all_variants)

    # Sort by chromosome and position
    df['chr_num'] = df['chr'].str.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')
    df['chr_num'] = pd.to_numeric(df['chr_num'], errors='coerce')
    df = df.sort_values(['chr_num', 'pos'])
    df = df.drop('chr_num', axis=1)

    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)

    # Print summary
    print(f"\n=== Variant Summary for {patient} - {sample_type} ===")
    print(f"Total variants: {len(df)}")
    print(f"Somatic variants: {len(df[df['somatic_status'] == 'Somatic'])}")

    if 'vaf' in df.columns:
        somatic_df = df[df['somatic_status'] == 'Somatic']
        if len(somatic_df) > 0:
            print(f"\nTop somatic variants by VAF:")
            top_variants = somatic_df.nlargest(10, 'vaf')[['chr', 'pos', 'ref', 'alt', 'vaf']]
            print(top_variants.to_string(index=False))

    print(f"\nSaved to: {output_file}")

if __name__ == '__main__':
    main()
