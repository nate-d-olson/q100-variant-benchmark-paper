#!/usr/bin/env python3
"""
Generate YAML config entries for draft benchmark sets with exclusions.
Calculates SHA256 checksums for all BED files.
"""

import os
import hashlib
import json
from pathlib import Path


def calculate_sha256(filepath):
    """Calculate SHA256 checksum of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def get_exclusion_name(filename):
    """
    Derive exclusion name from filename.
    Handles both single files and start/end pairs.
    """
    name = filename.replace("_sorted.bed", "")
    # For start/end pairs, remove the _start_sorted or _end_sorted suffix
    if name.endswith("_start"):
        name = name[:-6]  # Remove "_start"
    elif name.endswith("_end"):
        name = name[:-4]  # Remove "_end"
    return name


def process_draft_set(set_dir, set_name):
    """Process a single draft benchmark set directory."""
    # Find dip bed (benchmark.bed)
    dip_bed = None
    for f in os.listdir(set_dir):
        if f.endswith(".benchmark.bed"):
            dip_bed = os.path.join(set_dir, f)
            break

    if not dip_bed:
        print(f"Warning: No benchmark bed found in {set_name}")
        return None

    # Get exclusion BED files from intersections directory
    intersections_dir = os.path.join(set_dir, "exclusions", "intersections")
    if not os.path.exists(intersections_dir):
        print(f"Warning: No intersections directory in {set_name}")
        return None

    bed_files = sorted([f for f in os.listdir(intersections_dir) if f.endswith(".bed")])

    # Group files: identify pairs vs singles
    exclusions = []
    processed = set()

    for f in bed_files:
        if f in processed:
            continue

        filepath = os.path.join(intersections_dir, f)

        # Check if this is a start file of a pair
        if "_start_sorted.bed" in f:
            base_name = f.replace("_start_sorted.bed", "")
            end_file = f.replace("_start_sorted.bed", "_end_sorted.bed")

            if end_file in bed_files:
                # It's a pair
                end_filepath = os.path.join(intersections_dir, end_file)
                exclusion_name = get_exclusion_name(f)

                exclusions.append({
                    "name": exclusion_name,
                    "type": "pair",
                    "files": [
                        {
                            "path": filepath,
                            "sha256": calculate_sha256(filepath)
                        },
                        {
                            "path": end_filepath,
                            "sha256": calculate_sha256(end_filepath)
                        }
                    ]
                })
                processed.add(f)
                processed.add(end_file)
                continue

        # Skip end files (handled with their start pair)
        if "_end_sorted.bed" in f:
            continue

        # Single file
        exclusion_name = get_exclusion_name(f)
        exclusions.append({
            "name": exclusion_name,
            "type": "single",
            "files": [
                {
                    "path": filepath,
                    "sha256": calculate_sha256(filepath)
                }
            ]
        })
        processed.add(f)

    return {
        "dip_bed": {
            "path": dip_bed,
            "sha256": calculate_sha256(dip_bed)
        },
        "exclusions": exclusions
    }


def main():
    root_dir = "data/20250117_v0.020_HG002Q100v1.1/draft_benchmarksets"
    
    if not os.path.exists(root_dir):
        print(f"Error: Directory not found: {root_dir}")
        return

    sets = sorted([d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d))])

    # Map set directory names to config benchmark IDs
    set_to_id = {
        "CHM13_HG002-T2TQ100v1.1-dipz2k_smvar-excluded": "v5q_chm13_smvar",
        "CHM13_HG002-T2TQ100v1.1-dipz2k_stvar-excluded": "v5q_chm13_stvar",
        "GRCh37_HG002-T2TQ100v1.1-dipz2k_smvar-excluded": "v5q_grch37_smvar",
        "GRCh37_HG002-T2TQ100v1.1-dipz2k_stvar-excluded": "v5q_grch37_stvar",
        "GRCh38_HG002-T2TQ100v1.1-dipz2k_smvar-excluded": "v5q_grch38_smvar",
        "GRCh38_HG002-T2TQ100v1.1-dipz2k_stvar-excluded": "v5q_grch38_stvar",
    }

    config_data = {}

    for s in sets:
        set_dir = os.path.join(root_dir, s)
        benchmark_id = set_to_id.get(s, s)
        
        print(f"Processing {s} -> {benchmark_id}...", file=__import__('sys').stderr)
        result = process_draft_set(set_dir, s)
        
        if result:
            config_data[benchmark_id] = result

    # Output JSON for easier parsing
    import json
    print(json.dumps(config_data, indent=2))


if __name__ == "__main__":
    main()
