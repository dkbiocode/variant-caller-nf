#!/usr/bin/env python3
"""
Extract variant information from LoFreq VCF files and create a summary table.
LoFreq uses QUAL field for quality scores (not SSC like VarScan).
"""

import sys
import pandas as pd
from pathlib import Path

# Key genes from Table 4 (resistance mutations)
RESISTANCE_GENES = ['KRAS', 'MAP2K1', 'TP53', 'APC', 'PTEN', 'BRAF', 'NRAS', 'PIK3CA']

def parse_snpeff_annotation(ann_string):
    """Parse SnpEff ANN field - same as VarScan version."""
    if not ann_string:
        return {'gene': None, 'effect': None, 'impact': None, 'protein_change': None}

    annotations = ann_string.split(',')
    if not annotations:
        return {'gene': None, 'effect': None, 'impact': None, 'protein_change': None}

    fields = annotations[0].split('|')
    if len(fields) < 11:
        return {'gene': None, 'effect': None, 'impact': None, 'protein_change': None}

    return {
        'gene': fields[3] if fields[3] else None,
        'effect': fields[1] if fields[1] else None,
        'impact': fields[2] if fields[2] else None,
        'protein_change': fields[10] if fields[10] else None
    }

def parse_vcf_line(line, patient, sample_type):
    """Parse a LoFreq VCF data line and extract key information."""
    fields = line.strip().split('\t')

    # LoFreq VCFs have 8 columns (no FORMAT/sample columns)
    if len(fields) < 8:
        return None

    chrom = fields[0]
    pos = fields[1]
    ref = fields[3]
    alt = fields[4]
    qual = fields[5]
    info = fields[7]

    # Parse INFO field
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value

    # Extract VAF from LoFreq
    # LoFreq stores allele frequency in AF field
    vaf = None
    af = info_dict.get('AF', None)
    if af:
        try:
            vaf = float(af) * 100  # Convert to percentage
        except ValueError:
            vaf = None

    # Get quality score
    quality_score = None
    try:
        quality_score = float(qual) if qual != '.' else None
    except ValueError:
        quality_score = None

    # Get depth (DP field)
    depth = info_dict.get('DP', 'NA')

    # Parse SnpEff annotation if present
    ann_info = parse_snpeff_annotation(info_dict.get('ANN', ''))

    variant = {
        'patient': patient,
        'sample_type': sample_type,
        'chr': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'vaf': vaf,
        'quality': quality_score,
        'depth': depth,
        'gene': ann_info['gene'],
        'effect': ann_info['effect'],
        'impact': ann_info['impact'],
        'protein_change': ann_info['protein_change'],
        'caller': 'LoFreq'
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
        print("Usage: summarize_lofreq_variants.py <patient> <sample_type> <snv.vcf> <indel.vcf> [<output.tsv>]")
        sys.exit(1)

    patient = sys.argv[1]
    sample_type = sys.argv[2]
    snv_vcf = sys.argv[3]
    indel_vcf = sys.argv[4]
    output_file = sys.argv[5] if len(sys.argv) > 5 else f"{patient}_{sample_type}_lofreq_variants.tsv"

    # Parse both SNV and INDEL VCF files
    all_variants = []
    all_variants.extend(parse_vcf_file(snv_vcf, patient, sample_type))
    all_variants.extend(parse_vcf_file(indel_vcf, patient, sample_type))

    if not all_variants:
        # Create empty output file
        pd.DataFrame(columns=[
            'patient', 'sample_type', 'chr', 'pos', 'ref', 'alt',
            'vaf', 'quality', 'depth', 'gene', 'effect', 'impact', 'protein_change', 'caller'
        ]).to_csv(output_file, sep='\t', index=False)
        print(f"No variants found. Created empty file: {output_file}")
        return

    # Convert to DataFrame
    df = pd.DataFrame(all_variants)

    # Reorder columns for better readability
    column_order = [
        'patient', 'sample_type', 'chr', 'pos', 'ref', 'alt', 'vaf', 'quality', 'depth',
        'gene', 'effect', 'impact', 'protein_change', 'caller'
    ]
    existing_cols = [col for col in column_order if col in df.columns]
    df = df[existing_cols]

    # Sort by chromosome and position
    df['chr_num'] = df['chr'].str.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')
    df['chr_num'] = pd.to_numeric(df['chr_num'], errors='coerce')
    df = df.sort_values(['chr_num', 'pos'])
    df = df.drop('chr_num', axis=1)

    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)

    # Print summary
    print(f"\n=== LoFreq Variant Summary for {patient} - {sample_type} ===")
    print(f"Total variants: {len(df)}")

    # Show resistance gene variants if any
    if 'gene' in df.columns:
        resistance_vars = df[df['gene'].isin(RESISTANCE_GENES)]
        if len(resistance_vars) > 0:
            print(f"\nResistance gene variants ({len(resistance_vars)} found):")
            print(resistance_vars[['gene', 'chr', 'pos', 'effect', 'protein_change', 'vaf']].to_string(index=False))

    if 'vaf' in df.columns:
        valid_vaf = df[df['vaf'].notna()]
        if len(valid_vaf) > 0:
            print(f"\nTop variants by VAF:")
            display_cols = ['chr', 'pos', 'gene', 'effect', 'protein_change', 'vaf'] if 'gene' in df.columns else ['chr', 'pos', 'ref', 'alt', 'vaf']
            top_variants = valid_vaf.nlargest(10, 'vaf')[display_cols]
            print(top_variants.to_string(index=False))

    print(f"\nSaved to: {output_file}")

if __name__ == '__main__':
    main()
