#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from Claude (Anthropic)
# for code generation and debugging.
"""
Annotate old benchmark variants/regions with v5.0q exclusion status.

For variants in the old benchmark regions that are NOT in v5.0q benchmark
regions, determines whether they are:
- not_in_dipbed: region not in v5.0q dip.bed (assembly didn't cover it)
- excluded: in dip.bed but removed by one or more exclusions
- in_v5_dipbed_not_excluded: in dip.bed but not in any exclusion
"""

import csv
import gzip
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, List


def _run_cmd(cmd: str, log) -> subprocess.CompletedProcess:
    """Run a shell command and log errors."""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )
    if result.returncode != 0:
        log.write(f"ERROR running: {cmd}\n{result.stderr}\n")
    return result


def _compute_bed_size(bed_path: str, log) -> int:
    """Compute total bp in a BED file."""
    result = _run_cmd(
        f"awk '{{sum+=$3-$2}} END {{print sum+0}}' {bed_path}", log
    )
    return int(result.stdout.strip() or 0)


def _compute_variant_size(ref: str, alt: str, svlen: str | None) -> int:
    """Compute variant size."""
    if svlen and svlen != ".":
        try:
            return abs(int(svlen))
        except ValueError:
            pass
    return max(len(ref), len(alt)) - 1


def analyze_regions(
    old_bed: str,
    new_benchmark_bed: str,
    new_dip_bed: str,
    exclusion_beds: List[str],
    exclusion_names: List[str],
    output_csv: str,
    log,
) -> None:
    """Analyze old-only regions: are they in dip.bed? Which exclusions?"""
    with tempfile.TemporaryDirectory() as tmpdir:
        old_only = os.path.join(tmpdir, "old_only.bed")
        _run_cmd(
            f"bedtools subtract -a {old_bed} -b {new_benchmark_bed} "
            f"| bedtools sort -i - | bedtools merge -i - > {old_only}",
            log,
        )

        old_only_size = _compute_bed_size(old_only, log)
        log.write(f"Old-only region size: {old_only_size:,} bp\n")

        # Not in dip.bed
        not_in_dip = os.path.join(tmpdir, "not_in_dip.bed")
        _run_cmd(
            f"bedtools subtract -a {old_only} -b {new_dip_bed} > {not_in_dip}",
            log,
        )
        not_in_dip_bp = _compute_bed_size(not_in_dip, log)

        # In dip.bed
        in_dip = os.path.join(tmpdir, "in_dip.bed")
        _run_cmd(
            f"bedtools intersect -a {old_only} -b {new_dip_bed} "
            f"| bedtools sort -i - | bedtools merge -i - > {in_dip}",
            log,
        )
        in_dip_bp = _compute_bed_size(in_dip, log)

        # In dip.bed and in any exclusion
        excluded_bp = 0
        if exclusion_beds:
            all_excl = os.path.join(tmpdir, "all_excl.bed")
            beds_str = " ".join(exclusion_beds)
            _run_cmd(
                f"cat {beds_str} | bedtools sort -i - | bedtools merge -i - > {all_excl}",
                log,
            )
            in_excl = os.path.join(tmpdir, "in_excl.bed")
            _run_cmd(
                f"bedtools intersect -a {in_dip} -b {all_excl} "
                f"| bedtools sort -i - | bedtools merge -i - > {in_excl}",
                log,
            )
            excluded_bp = _compute_bed_size(in_excl, log)

        not_excluded_bp = in_dip_bp - excluded_bp

        with open(output_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["category", "bases_bp", "pct_of_old_only"])
            total = old_only_size if old_only_size > 0 else 1
            writer.writerow(["not_in_dipbed", not_in_dip_bp, f"{not_in_dip_bp / total * 100:.2f}"])
            writer.writerow(["excluded", excluded_bp, f"{excluded_bp / total * 100:.2f}"])
            writer.writerow(["in_v5_dipbed_not_excluded", not_excluded_bp, f"{not_excluded_bp / total * 100:.2f}"])

        log.write(f"Region analysis: not_in_dipbed={not_in_dip_bp}, excluded={excluded_bp}, not_excluded={not_excluded_bp}\n")


def analyze_variants(
    old_vcf: str,
    old_bed: str,
    new_benchmark_bed: str,
    new_dip_bed: str,
    exclusion_beds: List[str],
    exclusion_names: List[str],
    comp_type: str,
    output_tsv: str,
    output_summary: str,
    log,
) -> None:
    """Analyze old-only variants: are they in dip.bed? Which exclusions?"""
    is_smvar = comp_type == "smvar"

    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract variants from old VCF within old benchmark regions,
        # then filter to those NOT in v5.0q benchmark regions
        old_only_vcf = os.path.join(tmpdir, "old_only.vcf")
        _run_cmd(
            f"bcftools view -T {old_bed} {old_vcf} "
            f"| bcftools view -T ^{new_benchmark_bed} "
            f"| bcftools view --exclude 'GT=\"0/0\" || GT=\"./.\"' "
            f"> {old_only_vcf}",
            log,
        )

        # Parse variants and create position BED
        variants = []
        positions_bed = os.path.join(tmpdir, "positions.bed")

        open_func = gzip.open if old_only_vcf.endswith(".gz") else open
        with open(old_only_vcf, "r") as vcf_f, open(positions_bed, "w") as bed_f:
            for line in vcf_f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4].split(",")[0]  # First alt allele

                # Parse INFO for SVLEN, SVTYPE, END
                info = parts[7]
                svlen = None
                svtype = None
                end = pos + len(ref) - 1

                for field in info.split(";"):
                    if field.startswith("SVLEN="):
                        svlen = field.split("=", 1)[1]
                    elif field.startswith("SVTYPE="):
                        svtype = field.split("=", 1)[1]
                    elif field.startswith("END="):
                        try:
                            end = int(field.split("=", 1)[1])
                        except ValueError:
                            pass

                var_size = _compute_variant_size(ref, alt, svlen)

                # Size filter
                if is_smvar and var_size >= 50:
                    continue
                if not is_smvar and var_size < 50:
                    continue

                # Determine variant type
                if var_size == 0:
                    var_type = "SNP"
                elif is_smvar:
                    var_type = "INDEL"
                else:
                    if svtype and "INS" in svtype.upper():
                        var_type = "SV_INS"
                    elif svtype and "DEL" in svtype.upper():
                        var_type = "SV_DEL"
                    elif len(alt) > len(ref):
                        var_type = "SV_INS"
                    else:
                        var_type = "SV_DEL"

                # BED interval for overlap checking
                bed_start = pos - 1  # VCF is 1-based, BED is 0-based
                bed_end = max(pos, end)

                variants.append({
                    "chrom": chrom,
                    "pos": pos,
                    "end": end,
                    "ref": ref,
                    "alt": alt,
                    "var_type": var_type,
                    "bed_start": bed_start,
                    "bed_end": bed_end,
                })
                bed_f.write(f"{chrom}\t{bed_start}\t{bed_end}\n")

        log.write(f"Extracted {len(variants)} size-filtered old-only variants\n")

        if not variants:
            # Write empty outputs
            with open(output_tsv, "w") as f:
                f.write("chrom\tpos\tend\tref\talt\tvar_type\tin_dip_bed\texclusion_ids\tstatus\n")
            with open(output_summary, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["status", "variant_count", "pct_of_total"])
            log.write("No variants to analyze\n")
            return

        # Check dip.bed overlap
        dip_overlap_out = os.path.join(tmpdir, "dip_overlap.bed")
        _run_cmd(
            f"bedtools intersect -a {positions_bed} -b {new_dip_bed} -c > {dip_overlap_out}",
            log,
        )
        with open(dip_overlap_out, "r") as f:
            for i, line in enumerate(f):
                parts = line.strip().split("\t")
                if i < len(variants):
                    variants[i]["in_dip_bed"] = int(parts[-1]) > 0

        # Check each exclusion overlap
        for excl_bed, excl_name in zip(exclusion_beds, exclusion_names):
            excl_out = os.path.join(tmpdir, f"excl_{excl_name}.bed")
            _run_cmd(
                f"bedtools intersect -a {positions_bed} -b {excl_bed} -c > {excl_out}",
                log,
            )
            with open(excl_out, "r") as f:
                for i, line in enumerate(f):
                    parts = line.strip().split("\t")
                    if i < len(variants):
                        if int(parts[-1]) > 0:
                            existing = variants[i].get("exclusion_ids", [])
                            existing.append(excl_name)
                            variants[i]["exclusion_ids"] = existing

        # Classify and write output
        status_counts: Counter = Counter()

        with open(output_tsv, "w") as f:
            f.write("chrom\tpos\tend\tref\talt\tvar_type\tin_dip_bed\texclusion_ids\tstatus\n")
            for v in variants:
                in_dip = v.get("in_dip_bed", False)
                excl_ids = v.get("exclusion_ids", [])

                if not in_dip:
                    status = "not_in_dipbed"
                    excl_str = ""
                elif excl_ids:
                    status = "excluded"
                    excl_str = ",".join(excl_ids)
                else:
                    status = "in_v5_dipbed_not_excluded"
                    excl_str = ""

                status_counts[status] += 1
                f.write(
                    f"{v['chrom']}\t{v['pos']}\t{v['end']}\t{v['ref']}\t{v['alt']}\t"
                    f"{v['var_type']}\t{in_dip}\t{excl_str}\t{status}\n"
                )

        # Write summary
        total = len(variants)
        with open(output_summary, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["status", "variant_count", "pct_of_total"])
            for status in ["not_in_dipbed", "excluded", "in_v5_dipbed_not_excluded"]:
                count = status_counts.get(status, 0)
                pct = count / total * 100 if total > 0 else 0
                writer.writerow([status, count, f"{pct:.2f}"])

        log.write(f"Status counts: {dict(status_counts)}\n")


def main():
    comp_id = snakemake.wildcards.comp_id
    comp_type = snakemake.params.comp_type
    excl_names = list(snakemake.params.excl_names)

    log_path = str(snakemake.log[0])

    with open(log_path, "w") as log:
        log.write(f"Annotating old benchmark status for {comp_id}\n")
        log.write(f"Comparison type: {comp_type}\n\n")

        # Region analysis
        log.write("=== Region Analysis ===\n")
        analyze_regions(
            old_bed=str(snakemake.input.old_bed),
            new_benchmark_bed=str(snakemake.input.new_benchmark_bed),
            new_dip_bed=str(snakemake.input.new_dip_bed),
            exclusion_beds=[str(p) for p in snakemake.input.exclusion_beds],
            exclusion_names=excl_names,
            output_csv=str(snakemake.output.regions_csv),
            log=log,
        )

        # Variant analysis
        log.write("\n=== Variant Analysis ===\n")
        analyze_variants(
            old_vcf=str(snakemake.input.old_vcf),
            old_bed=str(snakemake.input.old_bed),
            new_benchmark_bed=str(snakemake.input.new_benchmark_bed),
            new_dip_bed=str(snakemake.input.new_dip_bed),
            exclusion_beds=[str(p) for p in snakemake.input.exclusion_beds],
            exclusion_names=excl_names,
            comp_type=comp_type,
            output_tsv=str(snakemake.output.variants_tsv),
            output_summary=str(snakemake.output.summary_csv),
            log=log,
        )

        log.write("\nAnalysis complete\n")


if __name__ == "__main__":
    main()
