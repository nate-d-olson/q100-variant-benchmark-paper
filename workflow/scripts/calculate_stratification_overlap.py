#!/usr/bin/env python3
"""
Calculate percent overlap of benchmark regions with stratification intervals.

For each interval in the stratification BED files, calculate what percentage
of that interval is covered by the benchmark regions. This allows visualization
of how much of difficult genomic contexts are included in the benchmark.

Outputs a table with:
- stratification: Name of the stratification (TR, HP, SD, etc.)
- chrom: Chromosome
- start: Interval start position
- end: Interval end position
- length: Length of the interval
- overlap_bp: Bases overlapping with benchmark regions
- percent_overlap: Percentage of interval covered by benchmark
"""

import sys
import subprocess
import pandas as pd
from pathlib import Path


def run_bedtools_intersect(benchmark_bed, strat_bed, output_file):
    """
    Run bedtools intersect to find overlaps between benchmark and stratification.
    
    Args:
        benchmark_bed: Path to benchmark regions BED file
        strat_bed: Path to stratification BED file
        output_file: Path for temporary output
        
    Returns:
        Path to output file
    """
    cmd = [
        "bedtools", "intersect",
        "-a", str(strat_bed),
        "-b", str(benchmark_bed),
        "-wa",  # Write original A entry
        "-wb",  # Write original B entry
    ]
    
    with open(output_file, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
    
    return output_file


def calculate_overlap_per_interval(benchmark_bed, strat_bed, strat_name):
    """
    Calculate overlap percentage for each interval in stratification.
    
    Args:
        benchmark_bed: Path to benchmark regions BED file
        strat_bed: Path to stratification BED file (can be gzipped)
        strat_name: Name of the stratification
        
    Returns:
        DataFrame with overlap statistics per interval
    """
    # Read stratification intervals
    if str(strat_bed).endswith('.gz'):
        strat_df = pd.read_csv(
            strat_bed, 
            sep='\t', 
            header=None, 
            names=['chrom', 'start', 'end'],
            compression='gzip'
        )
    else:
        strat_df = pd.read_csv(
            strat_bed, 
            sep='\t', 
            header=None, 
            names=['chrom', 'start', 'end']
        )
    
    # Calculate original interval lengths
    strat_df['length'] = strat_df['end'] - strat_df['start']
    strat_df['stratification'] = strat_name
    
    # Create temporary file for intersection results
    tmp_intersect = Path("tmp_intersect.bed")
    
    try:
        # Run bedtools intersect to get overlaps
        run_bedtools_intersect(benchmark_bed, strat_bed, tmp_intersect)
        
        # Read intersection results if file is not empty
        if tmp_intersect.stat().st_size > 0:
            # Columns: strat_chrom, strat_start, strat_end, bench_chrom, bench_start, bench_end
            intersect_df = pd.read_csv(
                tmp_intersect,
                sep='\t',
                header=None,
                names=['strat_chrom', 'strat_start', 'strat_end', 
                       'bench_chrom', 'bench_start', 'bench_end']
            )
            
            # Calculate overlap for each intersection
            # Overlap is the intersection of the two intervals
            intersect_df['overlap_start'] = intersect_df[['strat_start', 'bench_start']].max(axis=1)
            intersect_df['overlap_end'] = intersect_df[['strat_end', 'bench_end']].min(axis=1)
            intersect_df['overlap_bp'] = intersect_df['overlap_end'] - intersect_df['overlap_start']
            
            # Sum overlaps per stratification interval (handle multi-overlaps)
            overlap_summary = intersect_df.groupby(
                ['strat_chrom', 'strat_start', 'strat_end']
            )['overlap_bp'].sum().reset_index()
            
            # Merge with original intervals
            result_df = strat_df.merge(
                overlap_summary,
                left_on=['chrom', 'start', 'end'],
                right_on=['strat_chrom', 'strat_start', 'strat_end'],
                how='left'
            )
            
            # Fill NaN with 0 (no overlap)
            result_df['overlap_bp'] = result_df['overlap_bp'].fillna(0)
        else:
            # No overlaps found
            result_df = strat_df.copy()
            result_df['overlap_bp'] = 0
    
    finally:
        # Clean up temporary file
        if tmp_intersect.exists():
            tmp_intersect.unlink()
    
    # Calculate percent overlap
    result_df['percent_overlap'] = (result_df['overlap_bp'] / result_df['length'] * 100).round(2)
    
    # Select and order columns
    result_df = result_df[[
        'stratification', 'chrom', 'start', 'end', 
        'length', 'overlap_bp', 'percent_overlap'
    ]]
    
    return result_df


def main():
    """Main function to process all stratifications for a benchmark."""
    # Get arguments from Snakemake
    benchmark_bed = snakemake.input.benchmark_bed
    strat_beds = snakemake.input.strat_beds
    strat_names = snakemake.params.strat_names
    output_file = snakemake.output.table
    
    # Process each stratification
    all_results = []
    
    for strat_bed, strat_name in zip(strat_beds, strat_names):
        print(f"Processing {strat_name}...", file=sys.stderr)
        result_df = calculate_overlap_per_interval(benchmark_bed, strat_bed, strat_name)
        all_results.append(result_df)
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Sort by stratification and position
    combined_df = combined_df.sort_values(['stratification', 'chrom', 'start'])
    
    # Write output
    combined_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Wrote {len(combined_df)} intervals to {output_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
