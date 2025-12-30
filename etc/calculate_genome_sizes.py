#!/usr/bin/env python3
"""
Calculate assembled base counts for reference genomes.
"""

import pandas as pd
import os
from pathlib import Path

# Define paths to .fai files
base_dir = Path('..')
fai_files = {
    'GRCh37': base_dir / 'resources/references/GRCh37.fa.gz.fai',
    'GRCh38': base_dir / 'resources/references/GRCh38.fa.gz.fai',
    'CHM13v2.0': base_dir / 'resources/references/CHM13v2.0.fa.gz.fai',
    'hg002v1.1': base_dir / 'data/hg002v1.1.mat_Y_EBV_MT.fasta.gz.fai'
}

# Verify all files exist
print("Checking for .fai files...")
for ref, path in fai_files.items():
    if not path.exists():
        print(f"Warning: {ref} file not found at {path}")
    else:
        print(f"Found {ref}: {path}")

def read_fai_file(fai_path):
    """Read a FASTA index (.fai) file and return a DataFrame with chromosome lengths."""
    df = pd.read_csv(
        fai_path,
        sep='\t',
        header=None,
        names=['chrom', 'length', 'offset', 'linebases', 'linewidth']
    )
    return df[['chrom', 'length']]

def normalize_chrom_name(chrom):
    """Normalize chromosome names to standard format (chr1, chr2, ..., chrX, chrY)."""
    chrom = str(chrom)

    # Handle hg002 format (chr1_MATERNAL, chr1_PATERNAL)
    if '_MATERNAL' in chrom or '_PATERNAL' in chrom:
        base_chrom = chrom.split('_')[0]
        return base_chrom

    # Handle GRCh37 format (may not have 'chr' prefix)
    if not chrom.startswith('chr'):
        # Convert numeric chromosomes
        if chrom.isdigit() or chrom in ['X', 'Y']:
            return f"chr{chrom}"

    return chrom

def get_chrom_type(chrom):
    """Classify chromosome type."""
    chrom_normalized = normalize_chrom_name(chrom)

    # Extract number/letter after 'chr'
    if chrom_normalized.startswith('chr'):
        chrom_id = chrom_normalized[3:]

        # Autosomes (1-22)
        if chrom_id.isdigit() and 1 <= int(chrom_id) <= 22:
            return 'autosome'
        # Sex chromosomes
        elif chrom_id in ['X', 'Y']:
            return 'sex'

    return 'other'

# Test the functions
print("\nTesting normalization:")
test_chroms = ['chr1', '1', 'chr1_MATERNAL', 'chr1_PATERNAL', 'chrX', 'X', 'chrY', 'chrM', 'chr22']
for tc in test_chroms:
    print(f"{tc:20s} -> {normalize_chrom_name(tc):10s} -> {get_chrom_type(tc)}")

# Read all .fai files and process them
print("\n" + "="*60)
print("Processing reference genomes...")
print("="*60)
all_data = []

for ref_name, fai_path in fai_files.items():
    print(f"\nProcessing {ref_name}...")

    # Read the .fai file
    df = read_fai_file(fai_path)

    # Add normalized chromosome name and type
    df['chrom_normalized'] = df['chrom'].apply(normalize_chrom_name)
    df['chrom_type'] = df['chrom'].apply(get_chrom_type)
    df['reference'] = ref_name

    # For hg002, sum maternal and paternal chromosomes
    if ref_name == 'hg002v1.1':
        # Group by normalized chromosome name and sum lengths
        df_grouped = df.groupby(['chrom_normalized', 'chrom_type', 'reference'], as_index=False)['length'].sum()
        df_grouped = df_grouped.rename(columns={'chrom_normalized': 'chrom'})
        all_data.append(df_grouped)

        print(f"  Total chromosomes (after merging maternal/paternal): {len(df_grouped)}")
        print(f"  Sample chromosomes: {', '.join(df_grouped['chrom'].head(5).tolist())}")
    else:
        # For other references, use normalized names
        df_clean = df[['chrom_normalized', 'length', 'chrom_type', 'reference']].copy()
        df_clean = df_clean.rename(columns={'chrom_normalized': 'chrom'})
        all_data.append(df_clean)

        print(f"  Total chromosomes: {len(df_clean)}")
        print(f"  Sample chromosomes: {', '.join(df_clean['chrom'].head(5).tolist())}")

# Combine all data
combined_df = pd.concat(all_data, ignore_index=True)
print(f"\nTotal rows in combined dataset: {len(combined_df)}")
print(f"References: {combined_df['reference'].unique().tolist()}")

# Display summary by reference
print("\n" + "="*60)
print("Chromosome counts by reference and type:")
print("="*60)
summary = combined_df.groupby(['reference', 'chrom_type']).size().unstack(fill_value=0)
print(summary)

# Create pivot table for detailed chromosome data
display_df = combined_df.pivot(index='chrom', columns='reference', values='length')
display_df = display_df[['GRCh37', 'GRCh38', 'CHM13v2.0', 'hg002v1.1']]

print("\n" + "="*60)
print("First 25 chromosomes:")
print("="*60)
print(display_df.head(25).to_string())

# Calculate totals for each reference
print("\n" + "="*60)
print("Calculating totals...")
print("="*60)
totals = []

for ref_name in ['GRCh37', 'GRCh38', 'CHM13v2.0', 'hg002v1.1']:
    ref_data = combined_df[combined_df['reference'] == ref_name]

    # Autosome total (chr1-22)
    autosome_total = ref_data[ref_data['chrom_type'] == 'autosome']['length'].sum()

    # Autosome + sex chromosomes (chr1-22, X, Y)
    autosome_sex_total = ref_data[ref_data['chrom_type'].isin(['autosome', 'sex'])]['length'].sum()

    totals.append({
        'reference': ref_name,
        'autosome_bases': autosome_total,
        'autosome_sex_bases': autosome_sex_total
    })

    print(f"\n{ref_name}:")
    print(f"  Autosomes (chr1-22): {autosome_total:,} bp")
    print(f"  Autosomes + Sex (chr1-22,X,Y): {autosome_sex_total:,} bp")

totals_df = pd.DataFrame(totals)
print("\n" + "="*60)
print("Summary table:")
print("="*60)
print(totals_df.to_string(index=False))

# Save the detailed chromosome data to TSV
output_file = 'reference_genome_sizes.tsv'
display_df.to_csv(output_file, sep='\t')
print(f"\nSaved detailed chromosome data to: {output_file}")

# Also save the totals
totals_file = 'reference_genome_totals.tsv'
totals_df.to_csv(totals_file, sep='\t', index=False)
print(f"Saved totals to: {totals_file}")

# Display the totals in a format ready for the config file
print("\n" + "="*60)
print("YAML format for config file:")
print("="*60)
print("\ngenome_sizes:")
for _, row in totals_df.iterrows():
    print(f"  {row['reference']}:")
    print(f"    autosome_bases: {int(row['autosome_bases'])}")
    print(f"    autosome_sex_bases: {int(row['autosome_sex_bases'])}")
