---
agent: AIAgentExpert
description: "Change the pipeline from using multiple methods for variant counting to bcftools query based counts."
model: Claude Opus 4.5 (Preview) (copilot)
tools: ['vscode', 'execute', 'read', 'edit', 'search', 'web', 'azure-mcp/search', 'serena/*', 'agent', 'todo']
---
# Implementation Plan: bcftools Query-Based Variant Tables

Replace current variant counting methods (rtg-tools vcfstats, sv_len script, context_counts) with a unified bcftools annotate + query approach that generates per-benchmark TSV tables with difficult region, benchmark region, and exclusion annotations.

## Overview

### Current State
- **Chromosome-level counts**: `rtg-tools vcfstats` → complex text output
- **SV length distribution**: `bcftools query` → custom parsing
- **Context counts**: `bcftools query` → Python counting script

### Target State
- **Single TSV per benchmark** with all variant metadata and annotation flags
- Annotations added via `bcftools annotate` with combined BED files (flattened to handle overlaps)
- Multiallelic variants split via `bcftools norm -m -any`
- INFO fields extracted dynamically from VCF headers
- Post-processing script to expand annotation IDs into binary columns

## Implementation Steps

### Step 1: Update config/config.yaml

**Action**: Add new `stratifications` section with descriptions (move from `CONTEXT_PATHS` in common.smk)

**Purpose**: Define stratification regions (difficult contexts) dynamically instead of hardcoding them in `common.smk`.

```yaml
stratifications:
  TR:
    url: TODO
    checksum: TODO-sha256
    path: "resources/stratifications/LowComplexity/{ref}_AllTandemRepeats.bed.gz"
    description: "All tandem repeat regions"
  TR10kb:
    url: TODO
    checksum: TODO-sha256
    path: "resources/stratifications/LowComplexity/{ref}_AllTandemRepeats_ge101bp_slop5.bed.gz"
    description: "Tandem repeats ≥101bp with 5bp slop"
  HP:
    url: TODO
    checksum: TODO-sha256
    path: "resources/stratifications/LowComplexity/{ref}_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz"
    description: "Homopolymers ≥7bp (perfect) or ≥11bp (imperfect) with 5bp slop"
  SD:
    url: TODO
    checksum: TODO-sha256
    path: "resources/stratifications/SegmentalDuplications/{ref}_segdups.bed.gz"
    description: "All segmental duplications"
  SD10kb:
    url: TODO
    checksum: TODO-sha256
    path: "resources/stratifications/SegmentalDuplications/{ref}_segdups_gt10kb.bed.gz"
    description: "Segmental duplications >10kb"
  MAP:
    url: TODO
    checksum: TODO-sha256
    path: "resources/stratifications/Mappability/{ref}_lowmappabilityall.bed.gz"
    description: "Low mappability regions"
```

Add `description` field to each exclusion in benchmark `exclusions` lists:

```yaml
benchmarksets:
  v5q_chm13_smvar:
    # ... existing fields ...
    exclusions:
      - name: "consecutive-svs"
        description: "Regions with consecutive structural variants"
        type: "single"
        files:
          - path: "data/.../file.bed"
            sha256: "..."
```

### Step 2: Create workflow/scripts/combine_beds_with_id.py

Script to combine multiple BED files into one sorted file with ID column for bcftools annotate. **Crucially**, this script flattens overlapping intervals and merges IDs (comma-separated) so `bcftools annotate` can apply them to a single `String` INFO field.

**Input**: List of BED files with their annotation IDs
**Output**: Combined, sorted BED with 4th column as comma-separated annotation IDs

```python
#!/usr/bin/env python3
"""
Combine multiple BED files into a single sorted file with ID column.
Flattens overlapping intervals and merges IDs.

For use with bcftools annotate --columns CHROM,FROM,TO,<ID>

Usage:
    python combine_beds_with_id.py \
        --beds file1.bed:ID1 file2.bed:ID2 ... \
        --output combined.bed
"""
import argparse
import gzip
import sys
from collections import defaultdict

def open_bed(path):
    """Open BED file, handling gzip compression."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')

def flatten_intervals(intervals):
    """
    Flatten overlapping intervals and merge IDs.
    intervals: list of (start, end, id) tuples
    Returns: list of (start, end, id_string) tuples
    """
    if not intervals:
        return []

    # Get all unique boundaries
    boundaries = set()
    for start, end, annot_id in intervals:
        boundaries.add(start)
        boundaries.add(end)
    sorted_boundaries = sorted(list(boundaries))

    # Map intervals to their IDs for quick lookup
    # Optimization: Sort intervals by start time
    intervals.sort()
    
    flattened = []
    
    # Sweep line approach
    # Active intervals at current segment
    active_intervals = []
    idx = 0
    n = len(intervals)
    
    for i in range(len(sorted_boundaries) - 1):
        start = sorted_boundaries[i]
        end = sorted_boundaries[i+1]
        mid = (start + end) / 2
        
        # This is O(N*M) naive implementation which is fine for typical BED counts
        # ideally we use an IntervalTree or proper sweep line
        current_ids = set()
        for i_start, i_end, i_id in intervals:
             # Optimization: since intervals sorted by start, we can skip or break
             # but simple intersection check is safest
             if i_start <= start and i_end >= end:
                 current_ids.add(i_id)
             elif i_start >= end:
                 # Since intervals are sorted by start, if this interval starts after our segment ends,
                 # all subsequent intervals will also start after.
                 # (Wait, no, sorting by start doesn't guarantee end order, but it helps)
                 pass
        
        # Actually, let's use a simpler robust approach for the script:
        # Collect all IDs valid for the current segment (start, end)
        # Re-iterate is slow. 
        # Better: Standard sweep line with events.
        pass

    # Re-implementing with proper sweep line
    events = []
    for start, end, annot_id in intervals:
        events.append((start, 1, annot_id))
        events.append((end, -1, annot_id))
    
    events.sort() # Sort by position, then type (+1 before -1 at same pos? No, standard is usually closed-open)
    # BED is 0-based closed-open [start, end).
    # Event at 'end' should remove ID.
    
    current_ids = set()
    result = []
    
    for i in range(len(events) - 1):
        pos, type_, annot_id = events[i]
        
        if type_ == 1:
            current_ids.add(annot_id)
        else:
            if annot_id in current_ids:
                current_ids.remove(annot_id)
        
        next_pos = events[i+1][0]
        
        if next_pos > pos and current_ids:
            # Emit interval
            ids_str = ",".join(sorted(list(current_ids)))
            # Merge with previous if identical
            if result and result[-1][2] == ids_str and result[-1][1] == pos:
                result[-1] = (result[-1][0], next_pos, ids_str)
            else:
                result.append((pos, next_pos, ids_str))
                
    # Process last event
    return result

def main():
    parser = argparse.ArgumentParser(description='Combine BED files with ID column')
    parser.add_argument('--beds', nargs='+', required=True,
                        help='BED files with IDs in format path:ID')
    parser.add_argument('--output', required=True, help='Output BED file')
    args = parser.parse_args()

    intervals_by_chrom = defaultdict(list)
    
    for bed_spec in args.beds:
        path, annot_id = bed_spec.rsplit(':', 1)
        with open_bed(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    intervals_by_chrom[chrom].append((start, end, annot_id))

    # Sort chroms
    def sort_key(chrom):
        if chrom.startswith('chr'):
            chrom_num = chrom[3:]
        else:
            chrom_num = chrom
        try:
            return (0, int(chrom_num))
        except ValueError:
            return (1, chrom_num)

    sorted_chroms = sorted(intervals_by_chrom.keys(), key=sort_key)

    with open(args.output, 'w') as out:
        for chrom in sorted_chroms:
            flattened = flatten_intervals(intervals_by_chrom[chrom])
            for start, end, ids in flattened:
                out.write(f'{chrom}\t{start}\t{end}\t{ids}\n')

if __name__ == '__main__':
    main()
```

### Step 3: Create workflow/scripts/generate_header_lines.py

Script to generate VCF INFO header lines. Now generates `STRAT_IDS` and `REGION_IDS` as String fields.

```python
#!/usr/bin/env python3
"""
Generate VCF INFO header lines for bcftools annotate.

Usage:
    python generate_header_lines.py \
        --benchmark-regions \
        --output header_lines.txt
"""
import argparse

def main():
    parser = argparse.ArgumentParser(description='Generate VCF header lines')
    parser.add_argument('--benchmark-regions', action='store_true',
                        help='Include benchmark regions header')
    parser.add_argument('--output', required=True, help='Output header file')
    args = parser.parse_args()

    lines = []

    # Stratification IDs field
    lines.append(
        '##INFO=<ID=STRAT_IDS,Number=.,Type=String,'
        'Description="Comma-separated list of stratification region IDs overlapping variant">'
    )

    # Region IDs field (includes Benchmark Regions and Exclusions)
    lines.append(
        '##INFO=<ID=REGION_IDS,Number=.,Type=String,'
        'Description="Comma-separated list of region IDs (Benchmark, Exclusions) overlapping variant">'
    )

    with open(args.output, 'w') as f:
        f.write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
```

### Step 4: Create workflow/scripts/extract_info_fields.py

Script to dynamically extract INFO field names from VCF header.

```python
#!/usr/bin/env python3
"""
Extract INFO field names from VCF header for bcftools query.

Usage:
    python extract_info_fields.py --vcf input.vcf.gz --output fields.txt
"""
import argparse
import gzip
import re

def main():
    parser = argparse.ArgumentParser(description='Extract INFO fields from VCF')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--output', required=True, help='Output file with field names')
    parser.add_argument('--exclude', nargs='*', default=[],
                        help='INFO fields to exclude from output')
    args = parser.parse_args()

    info_fields = []
    opener = gzip.open if args.vcf.endswith('.gz') else open

    with opener(args.vcf, 'rt') as f:
        for line in f:
            if not line.startswith('##'):
                break  # End of header
            if line.startswith('##INFO='):
                match = re.search(r'ID=([^,>]+)', line)
                if match:
                    field_id = match.group(1)
                    if field_id not in args.exclude:
                        info_fields.append(field_id)

    with open(args.output, 'w') as f:
        f.write('\n'.join(info_fields))

if __name__ == '__main__':
    main()
```

### Step 5: Create workflow/scripts/expand_annotations.py

New script to post-process the TSV, expanding `STRAT_IDS` and `REGION_IDS` into binary columns.

```python
#!/usr/bin/env python3
"""
Expand STRAT_IDS and REGION_IDS columns into binary flags in TSV.
Reads from stdin, writes to stdout.

Usage:
    python expand_annotations.py --strat-ids ID1 ID2 ... --region-ids ID3 ID4 ...
"""
import argparse
import sys
import csv

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strat-ids', nargs='*', default=[])
    parser.add_argument('--region-ids', nargs='*', default=[])
    args = parser.parse_args()

    reader = csv.DictReader(sys.stdin, delimiter='\t')
    
    # Prepare new fieldnames
    fieldnames = list(reader.fieldnames)
    
    # Remove container columns
    if 'STRAT_IDS' in fieldnames:
        fieldnames.remove('STRAT_IDS')
    if 'REGION_IDS' in fieldnames:
        fieldnames.remove('REGION_IDS')
        
    # Add expansion columns
    # We prefix them to match what was expected (STRAT_..., EXCL_..., BMKREGIONS)
    # Stratifications: STRAT_{ID}
    strat_cols = [f"STRAT_{sid}" for sid in args.strat_ids]
    fieldnames.extend(strat_cols)
    
    # Regions: {ID} (BMKREGIONS or EXCL_...)
    region_cols = args.region_ids
    fieldnames.extend(region_cols)
    
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    
    for row in reader:
        # Parse STRAT_IDS
        strats = set(row.pop('STRAT_IDS', '').split(','))
        for sid, col in zip(args.strat_ids, strat_cols):
            row[col] = '1' if sid in strats else '0'
            
        # Parse REGION_IDS
        regions = set(row.pop('REGION_IDS', '').split(','))
        for rid in args.region_ids:
            row[rid] = '1' if rid in regions else '0'
            
        writer.writerow(row)

if __name__ == '__main__':
    main()
```

### Step 6: Create workflow/rules/var_tables.smk

```python
"""
Rules for generating variant tables via bcftools annotate and query.
"""

def get_stratification_beds(wildcards):
    """Get list of stratification BED files."""
    ref = get_reference_for_benchmark(wildcards.benchmark)
    strats = config.get("stratifications", {})
    beds = []
    for name, strat in strats.items():
        path = f"resources/stratifications/{strat['path_template'].format(ref=ref)}"
        beds.append(f"{path}:{name}")
    return beds

def get_region_beds(wildcards):
    """Get benchmark region BED and exclusion BEDs."""
    benchmark = wildcards.benchmark
    beds = []
    # Benchmark regions
    bed_path = get_benchmark_bed(benchmark)
    beds.append(f"{bed_path}:BMKREGIONS")
    # Exclusions
    if benchmark.startswith("v5q_"):
        for excl in get_exclusion_config(benchmark):
            name = excl["name"].replace("-", "_").upper()
            for f in excl["files"]:
                beds.append(f"{f['path']}:EXCL_{name}")
    return beds

def get_strat_ids(wildcards):
    return list(config.get("stratifications", {}).keys())

def get_region_ids(wildcards):
    ids = ["BMKREGIONS"]
    if wildcards.benchmark.startswith("v5q_"):
        for excl in get_exclusion_config(wildcards.benchmark):
            name = excl["name"].replace("-", "_").upper()
            ids.append(f"EXCL_{name}")
    return ids

def get_var_table_query_format(wildcards):
    """Build bcftools query format string."""
    # Base fields
    base_fields = [
        "CHROM", "POS", "END", "GT", "VKX", "N_ALT", "TYPE", 
        "STRLEN_REF", "STRLEN_ALT", "ILEN"
    ]
    # Annotation containers
    annotation_fields = ["INFO/STRAT_IDS", "INFO/REGION_IDS"]
    
    # All fields with % prefix
    all_fields = base_fields + annotation_fields
    
    # We also need dynamic INFO fields, which we'll append in the rule
    # but for this function we just return the fixed part
    return "\\t".join([f"%{field}" for field in all_fields])

rule combine_stratification_beds:
    input:
        beds=lambda w: [b.split(":")[0] for b in get_stratification_beds(w)],
    output:
        bed=temp("results/var_tables/{benchmark}/strat_combined.bed"),
        bed_gz=ensure("results/var_tables/{benchmark}/strat_combined.bed.gz", non_empty=True),
        tbi="results/var_tables/{benchmark}/strat_combined.bed.gz.tbi",
    params:
        bed_specs=lambda w: get_stratification_beds(w),
    conda: "../envs/bedtools.yaml" # or python if no bedtools deps needed
    shell:
        """
        python workflow/scripts/combine_beds_with_id.py \
            --beds {params.bed_specs} \
            --output {output.bed}
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        """

rule combine_region_beds:
    input:
        beds=lambda w: [b.split(":")[0] for b in get_region_beds(w)],
    output:
        bed=temp("results/var_tables/{benchmark}/region_combined.bed"),
        bed_gz=ensure("results/var_tables/{benchmark}/region_combined.bed.gz", non_empty=True),
        tbi="results/var_tables/{benchmark}/region_combined.bed.gz.tbi",
    params:
        bed_specs=lambda w: get_region_beds(w),
    conda: "../envs/bedtools.yaml"
    shell:
        """
        python workflow/scripts/combine_beds_with_id.py \
            --beds {params.bed_specs} \
            --output {output.bed}
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        """

rule generate_annotation_headers:
    output: headers="results/var_tables/{benchmark}/annotation_headers.txt"
    shell:
        """
        python workflow/scripts/generate_header_lines.py \
            --benchmark-regions \
            --output {output.headers}
        """

rule annotate_vcf_stratifications:
    input:
        vcf=lambda w: get_benchmark_vcf(w.benchmark),
        strat_bed="results/var_tables/{benchmark}/strat_combined.bed.gz",
        strat_tbi="results/var_tables/{benchmark}/strat_combined.bed.gz.tbi",
        headers="results/var_tables/{benchmark}/annotation_headers.txt",
    output:
        vcf=temp("results/annotated_vcfs/{benchmark}/strat_annotated.vcf.gz"),
        tbi=temp("results/annotated_vcfs/{benchmark}/strat_annotated.vcf.gz.tbi"),
    conda: "../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate \
            --annotations {input.strat_bed} \
            --columns CHROM,FROM,TO,INFO/STRAT_IDS \
            --header-lines {input.headers} \
            --min-overlap 0.20 \
            --output-type z --output {output.vcf} {input.vcf}
        bcftools index --tbi {output.vcf}
        """

rule annotate_vcf_regions:
    input:
        vcf="results/annotated_vcfs/{benchmark}/strat_annotated.vcf.gz",
        region_bed="results/var_tables/{benchmark}/region_combined.bed.gz",
    output:
        vcf=ensure("results/annotated_vcfs/{benchmark}/fully_annotated.vcf.gz", non_empty=True),
        tbi="results/annotated_vcfs/{benchmark}/fully_annotated.vcf.gz.tbi",
    conda: "../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate \
            --annotations {input.region_bed} \
            --columns CHROM,FROM,TO,INFO/REGION_IDS \
            --output-type z --output {output.vcf} {input.vcf}
        bcftools index --tbi {output.vcf}
        """

rule split_multiallelics:
    input:
        vcf="results/annotated_vcfs/{benchmark}/fully_annotated.vcf.gz",
        ref=lambda w: config["references"][get_reference_for_benchmark(w.benchmark).lower()]["path"],
    output:
        vcf=ensure("results/annotated_vcfs/{benchmark}/split_annotated.vcf.gz", non_empty=True),
        tbi="results/annotated_vcfs/{benchmark}/split_annotated.vcf.gz.tbi",
    conda: "../envs/bcftools.yaml"
    shell:
        """
        bcftools norm --multiallelics -any --fasta-ref {input.ref} \
            --output-type z --output {output.vcf} {input.vcf}
        bcftools index --tbi {output.vcf}
        """

rule extract_info_fields:
    input: vcf="results/annotated_vcfs/{benchmark}/split_annotated.vcf.gz"
    output: fields="results/var_tables/{benchmark}/info_fields.txt"
    params:
        # Exclude our custom fields from the dynamic list since we add them explicitly
        exclude=["STRAT_IDS", "REGION_IDS"]
    shell:
        """
        python workflow/scripts/extract_info_fields.py \
            --vcf {input.vcf} \
            --exclude {params.exclude} \
            --output {output.fields}
        """

rule generate_var_table:
    input:
        vcf="results/annotated_vcfs/{benchmark}/split_annotated.vcf.gz",
        fields="results/var_tables/{benchmark}/info_fields.txt",
    output:
        tsv=ensure("results/var_tables/{benchmark}/variants.tsv", non_empty=True),
    params:
        base_fmt=get_var_table_query_format,
        strat_ids=lambda w: get_strat_ids(w),
        region_ids=lambda w: get_region_ids(w),
    conda: "../envs/bcftools.yaml"
    shell:
        """
        # Construct full format string: base_fmt + dynamic fields
        # Read dynamic fields, prepend %INFO/ to each
        DYN_FIELDS=$(cat {input.fields} | sed 's/^/%INFO\//' | paste -sd'\\t' -)
        FULL_FMT="{params.base_fmt}\\t$DYN_FIELDS\\n"
        
        # Build header line
        # Base headers match get_var_table_query_format
        BASE_HDR="CHROM\\tPOS\\tEND\\tGT\\tVKX\\tN_ALT\\tTYPE\\tSTRLEN_REF\\tSTRLEN_ALT\\tILEN\\tSTRAT_IDS\\tREGION_IDS"
        DYN_HDR=$(cat {input.fields} | paste -sd'\\t' -)
        HEADER="$BASE_HDR\\t$DYN_HDR"
        
        # Run query and pipe to expand script
        (echo "$HEADER"; bcftools query -f "$FULL_FMT" {input.vcf}) | \
        python workflow/scripts/expand_annotations.py \
            --strat-ids {params.strat_ids} \
            --region-ids {params.region_ids} > {output.tsv}
        """

def get_var_table_inputs(wildcards):
    return [f"results/var_tables/{b}/variants.tsv" for b in config["benchmarksets"]]
```

### Steps 7, 8, 9

Unchanged (updates to common.smk, vcf_processing.smk, Snakefile, README).
