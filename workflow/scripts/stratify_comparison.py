import pandas as pd
import subprocess
import os
import tempfile
import logging
import gzip
from typing import List, Dict, Any

# Setup logging
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def run_bedtools_intersect_vcf(vcf_file: str, bed_file: str) -> int:
    """
    Count variants in VCF overlapping with BED regions.

    Uses bedtools intersect with -u flag to identify unique overlapping
    variant records (excluding duplicates).

    Args:
        vcf_file: Path to input VCF file (can be gzipped)
        bed_file: Path to BED file defining regions of interest

    Returns:
        Number of variants overlapping with BED regions

    Raises:
        FileNotFoundError: If input files don't exist
        RuntimeError: If bedtools fails

    Example:
        >>> count = run_bedtools_intersect_vcf("variants.vcf.gz", "regions.bed")
        >>> print(f"Found {count} variants in regions")
    """
    # Validate inputs
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"BED file not found: {bed_file}")

    try:
        # Use subprocess to pipe
        # bedtools intersect -u -a {vcf} -b {bed}
        p1 = subprocess.Popen(
            ["bedtools", "intersect", "-u", "-a", vcf_file, "-b", bed_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        # Grep to exclude header (though bedtools output shouldn't have header if only overlapping lines?
        # VCF entries in bedtools output still have headers? No, intersect usually strips header if not -header.
        # But let's be safe: count lines.
        # Wait, bedtools intersect on VCF returns VCF lines.
        output, err = p1.communicate()
        if p1.returncode != 0:
            error_msg = (
                f"bedtools intersect failed:\n"
                f"  VCF: {vcf_file}\n"
                f"  BED: {bed_file}\n"
                f"  Error: {err.decode()}"
            )
            logging.error(error_msg)
            raise RuntimeError(error_msg)

        # Count non-header lines
        count = 0
        for line in output.splitlines():
            if not line.startswith(b'#'):
                count += 1
        return count
    except FileNotFoundError:
        raise
    except RuntimeError:
        raise
    except Exception as e:
        raise RuntimeError(
            f"Unexpected error in intersect_vcf:\n"
            f"  VCF: {vcf_file}\n"
            f"  BED: {bed_file}\n"
            f"  Error: {str(e)}"
        ) from e

def run_bedtools_intersect_bed_bases(bed_a, bed_b):
    """
    Calculate number of bases in intersection of bed_a and bed_b.
    bedtools intersect -a bed_a -b bed_b | awk '{sum+=$3-$2} END {print sum}'
    """
    try:
        p1 = subprocess.Popen(
            ["bedtools", "intersect", "-a", bed_a, "-b", bed_b],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        output, err = p1.communicate()
        if p1.returncode != 0:
            logging.error(f"bedtools intersect beds failed: {err.decode()}")
            return 0
        
        total_bases = 0
        for line in output.splitlines():
            parts = line.decode().split('\t')
            if len(parts) >= 3:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    total_bases += (end - start)
                except ValueError:
                    continue
        return total_bases
    except Exception as e:
        logging.error(f"Error in intersect_bed_bases: {e}")
        return 0

def run_bedtools_subtract_bed_bases(bed_a, bed_b):
    """
    Calculate bases in bed_a NOT in bed_b.
    bedtools subtract -a bed_a -b bed_b
    """
    try:
        p1 = subprocess.Popen(
            ["bedtools", "subtract", "-a", bed_a, "-b", bed_b],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        output, err = p1.communicate()
        
        total_bases = 0
        for line in output.splitlines():
            parts = line.decode().split('\t')
            if len(parts) >= 3:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    total_bases += (end - start)
                except ValueError:
                    continue
        return total_bases
    except Exception as e:
        logging.error(f"Error in subtract_bed_bases: {e}")
        return 0

def main():
    logging.info("Starting stratification analysis")
    
    # Inputs
    tp_vcf = snakemake.input.tp
    fp_vcf = snakemake.input.fp
    fn_vcf = snakemake.input.fn
    
    new_bed = snakemake.input.new_bed
    old_bed = snakemake.input.old_bed
    
    strat_beds = snakemake.input.strat_beds
    # strat_beds is a list. The rule should probably pass names/IDs too?
    # Or we parse filename.
    # We can get the ID from the filename or wildcards if passed, but snakemake input is just files.
    # We'll use os.path.basename and strip extension.
    
    # Global Counts (No Stratification)
    logging.info("Calculating global stats")

    # Variants
    # Just count lines in input VCFs (excluding header)
    def count_vcf(f):
        open_func = gzip.open if f.endswith('.gz') else open
        c = 0
        try:
            with open_func(f, 'rb') as fh:
                for line in fh:
                    if not line.startswith(b'#'):
                        c += 1
        except Exception as e:
            logging.error(f"Error counting lines in {f}: {e}")
            return 0
        return c

    tp_count = count_vcf(tp_vcf)
    fp_count = count_vcf(fp_vcf)
    fn_count = count_vcf(fn_vcf)
    
    var_results = [{
        "Stratification": "Global",
        "Shared_Variants": tp_count,
        "Unique_New_Variants": fp_count,
        "Unique_Old_Variants": fn_count
    }]
    
    # Regions
    # Intersection of New & Old
    shared_bases = run_bedtools_intersect_bed_bases(new_bed, old_bed)
    # Unique New (New - Old)
    unique_new_bases = run_bedtools_subtract_bed_bases(new_bed, old_bed)
    # Unique Old (Old - New)
    unique_old_bases = run_bedtools_subtract_bed_bases(old_bed, new_bed)
    
    region_results = [{
        "Stratification": "Global",
        "Shared_Bases": shared_bases,
        "Unique_New_Bases": unique_new_bases,
        "Unique_Old_Bases": unique_old_bases
    }]
    
    # Stratified Counts
    # To be efficient, we should create temporary BED files for:
    # 1. Global Shared Regions (BED)
    # 2. Global Unique New Regions (BED)
    # 3. Global Unique Old Regions (BED)
    # Then intersect these with each strat bed.
    
    with tempfile.TemporaryDirectory() as tmpdir:
        shared_bed_tmp = os.path.join(tmpdir, "shared.bed")
        unique_new_bed_tmp = os.path.join(tmpdir, "unique_new.bed")
        unique_old_bed_tmp = os.path.join(tmpdir, "unique_old.bed")
        
        # Generate these temp beds
        with open(shared_bed_tmp, 'w') as f:
            subprocess.run(["bedtools", "intersect", "-a", new_bed, "-b", old_bed], stdout=f)
        
        with open(unique_new_bed_tmp, 'w') as f:
            subprocess.run(["bedtools", "subtract", "-a", new_bed, "-b", old_bed], stdout=f)
            
        with open(unique_old_bed_tmp, 'w') as f:
            subprocess.run(["bedtools", "subtract", "-a", old_bed, "-b", new_bed], stdout=f)
            
        logging.info("Generated temp global region beds")
        
        for strat_bed in strat_beds:
            strat_name = os.path.basename(strat_bed).replace(".bed.gz", "").replace(".bed", "")
            # Remove ref prefix if present (e.g. GRCh38_)
            if "_" in strat_name:
                parts = strat_name.split("_", 1)
                # If first part is a known ref, strip it
                if parts[0] in ["GRCh37", "GRCh38", "CHM13v2.0"]:
                    strat_name = parts[1]
            
            logging.info(f"Processing {strat_name}")
            
            # Variant Intersections
            s_tp = run_bedtools_intersect_vcf(tp_vcf, strat_bed)
            s_fp = run_bedtools_intersect_vcf(fp_vcf, strat_bed)
            s_fn = run_bedtools_intersect_vcf(fn_vcf, strat_bed)
            
            var_results.append({
                "Stratification": strat_name,
                "Shared_Variants": s_tp,
                "Unique_New_Variants": s_fp,
                "Unique_Old_Variants": s_fn
            })
            
            # Region Intersections
            # Intersect Strat Bed with the Temp Global Beds
            s_shared = run_bedtools_intersect_bed_bases(shared_bed_tmp, strat_bed)
            s_unique_new = run_bedtools_intersect_bed_bases(unique_new_bed_tmp, strat_bed)
            s_unique_old = run_bedtools_intersect_bed_bases(unique_old_bed_tmp, strat_bed)
            
            region_results.append({
                "Stratification": strat_name,
                "Shared_Bases": s_shared,
                "Unique_New_Bases": s_unique_new,
                "Unique_Old_Bases": s_unique_old
            })

    # Save outputs
    pd.DataFrame(var_results).to_csv(snakemake.output.variants_csv, index=False)
    pd.DataFrame(region_results).to_csv(snakemake.output.regions_csv, index=False)
    
    logging.info("Analysis complete")

if __name__ == "__main__":
    main()
