#!/usr/bin/env python3
"""Create a compact GRCh38 test/debug dataset from full v5.0q files."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, TextIO


@dataclass(frozen=True)
class Region:
    """Genomic region definition."""

    chrom: str
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    label: str

    @property
    def length_bp(self) -> int:
        return self.end - self.start

    @property
    def id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


REGIONS: tuple[Region, ...] = (
    Region("chr1", 125_000_000, 125_200_000, "overlap_rich_1"),
    Region("chr4", 131_600_000, 131_800_000, "overlap_rich_2"),
    Region("chr7", 62_400_000, 62_600_000, "overlap_rich_3"),
    Region("chr9", 134_832_420, 134_842_674, "no_overlap_control"),
    Region("chr17", 26_800_000, 27_000_000, "overlap_rich_4"),
    Region("chrX", 200_000, 400_000, "overlap_rich_X"),
    Region("chrY", 10_200_000, 10_400_000, "overlap_rich_Y"),
)

NO_OVERLAP_LABEL = "no_overlap_control"

BENCHMARKSETS_DIR = Path("resources/benchmarksets")
EXCLUSIONS_DIR = Path("resources/exclusions/v5.0q_GRCh38_stvar")
STRATIFICATIONS_DIR = Path("resources/stratifications")

INPUT_VCFS = (
    BENCHMARKSETS_DIR / "v5.0q_GRCh38_smvar_benchmark.vcf.gz",
    BENCHMARKSETS_DIR / "v5.0q_GRCh38_stvar_benchmark.vcf.gz",
)
INPUT_BEDS = (
    BENCHMARKSETS_DIR / "v5.0q_GRCh38_smvar_benchmark.bed",
    BENCHMARKSETS_DIR / "v5.0q_GRCh38_stvar_benchmark.bed",
    BENCHMARKSETS_DIR / "v5.0q_GRCh38_smvar_dip.bed",
    BENCHMARKSETS_DIR / "v5.0q_GRCh38_stvar_dip.bed",
)

CONTEXT_FILES = tuple(
    STRATIFICATIONS_DIR / f"GRCh38_{name}.bed.gz"
    for name in ("HP", "MAP", "SD", "SD10kb", "TR", "TR10kb")
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compact GRCh38 debug subset dataset."
    )
    parser.add_argument(
        "--output-root",
        default="tests/fixtures/grch38_debug_subset",
        help="Directory where subset files are written.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove output-root before writing new files.",
    )
    return parser.parse_args()


def open_text(path: Path, mode: str) -> TextIO:
    """Open text file with transparent gzip handling."""
    if "b" in mode:
        raise ValueError("Binary mode is not supported in open_text().")
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def region_by_chrom(regions: Iterable[Region]) -> Dict[str, Region]:
    out: Dict[str, Region] = {}
    for region in regions:
        if region.chrom in out:
            raise ValueError(f"Duplicate region for chromosome: {region.chrom}")
        out[region.chrom] = region
    return out


def overlaps(region: Region, start0: int, end0: int) -> bool:
    return end0 > region.start and start0 < region.end


def clip_to_region(region: Region, start0: int, end0: int) -> tuple[int, int] | None:
    clipped_start = max(region.start, start0)
    clipped_end = min(region.end, end0)
    if clipped_end <= clipped_start:
        return None
    return clipped_start, clipped_end


def write_regions_bed(path: Path, regions: Iterable[Region]) -> None:
    ensure_parent(path)
    with open_text(path, "wt") as handle:
        for region in regions:
            handle.write(
                f"{region.chrom}\t{region.start}\t{region.end}\t{region.label}\n"
            )


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if info == ".":
        return out
    for token in info.split(";"):
        if not token:
            continue
        if "=" not in token:
            out[token] = "True"
            continue
        key, value = token.split("=", 1)
        out[key] = value
    return out


def parse_vcf_interval(fields: list[str]) -> tuple[int, int, dict[str, str]]:
    """Return 0-based [start, end) interval and INFO dict for a VCF record."""
    pos1 = int(fields[1])
    ref = fields[3]
    info = parse_info(fields[7])
    start0 = pos1 - 1

    end1 = None
    if "END" in info:
        raw = info["END"].split(",")[0]
        try:
            end1 = int(float(raw))
        except ValueError:
            end1 = None
    if end1 is None:
        end1 = pos1 + max(1, len(ref)) - 1

    end0 = max(start0 + 1, end1)
    return start0, end0, info


def parse_svlen(info: dict[str, str]) -> int | None:
    if "SVLEN" not in info:
        return None
    raw = info["SVLEN"].split(",")[0]
    try:
        return abs(int(float(raw)))
    except ValueError:
        return None


def first_alt(alt_field: str) -> str:
    return alt_field.split(",")[0]


def allele_length(allele: str) -> int:
    if allele in {"", "."}:
        return 0
    if allele.startswith("<") and allele.endswith(">"):
        return 0
    return len(allele)


def classify_variant(ref: str, alt: str, info: dict[str, str]) -> str:
    if "SVTYPE" in info and info["SVTYPE"] not in {"", "."}:
        return info["SVTYPE"].split(",")[0]

    ref_len = allele_length(ref)
    alt_len = allele_length(first_alt(alt))
    if ref_len == 1 and alt_len == 1:
        return "SNP"
    if ref_len == alt_len and ref_len > 1:
        return "MNP"
    if alt_len > ref_len:
        return "INS"
    if ref_len > alt_len:
        return "DEL"
    return "COMPLEX"


def variant_size_bp(ref: str, alt: str, start0: int, end0: int, info: dict[str, str]) -> int:
    svlen = parse_svlen(info)
    if svlen is not None:
        return max(1, svlen)

    ref_len = allele_length(ref)
    alt_len = allele_length(first_alt(alt))
    delta = abs(alt_len - ref_len)
    span = max(1, end0 - start0)
    return max(1, delta, span)


def size_bin(size_bp: int) -> str:
    if size_bp == 1:
        return "1"
    if size_bp <= 5:
        return "2-5"
    if size_bp <= 20:
        return "6-20"
    if size_bp <= 50:
        return "21-50"
    if size_bp <= 200:
        return "51-200"
    if size_bp <= 1000:
        return "201-1000"
    if size_bp <= 5000:
        return "1001-5000"
    return ">5000"


def filter_vcf(
    input_path: Path,
    output_path: Path,
    regions_map: Dict[str, Region],
) -> dict[str, object]:
    """Filter VCF records to region subset."""
    ensure_parent(output_path)
    kept = 0
    type_counts: Counter[str] = Counter()
    size_bins: Counter[str] = Counter()
    per_region_counts: Counter[str] = Counter()
    large_variants = 0

    with open_text(input_path, "rt") as src, open_text(output_path, "wt") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            region = regions_map.get(chrom)
            if region is None:
                continue

            start0, end0, info = parse_vcf_interval(fields)
            if not overlaps(region, start0, end0):
                continue

            dst.write(line)
            kept += 1

            ref = fields[3]
            alt = fields[4]
            vtype = classify_variant(ref, alt, info)
            type_counts[vtype] += 1

            var_size = variant_size_bp(ref, alt, start0, end0, info)
            size_bucket = size_bin(var_size)
            size_bins[size_bucket] += 1
            if var_size > 5000:
                large_variants += 1

            per_region_counts[region.label] += 1

    return {
        "input": str(input_path),
        "output": str(output_path),
        "variants_kept": kept,
        "variant_type_counts": dict(sorted(type_counts.items())),
        "size_bin_counts": dict(
            (bin_name, size_bins.get(bin_name, 0))
            for bin_name in (
                "1",
                "2-5",
                "6-20",
                "21-50",
                "51-200",
                "201-1000",
                "1001-5000",
                ">5000",
            )
        ),
        "variants_gt_5kb": large_variants,
        "per_region_variant_counts": dict(sorted(per_region_counts.items())),
    }


def clip_bed(
    input_path: Path,
    output_path: Path,
    regions_map: Dict[str, Region],
) -> dict[str, object]:
    """Clip BED intervals to region subset."""
    ensure_parent(output_path)
    interval_count = 0
    bp_count = 0
    per_region_bp: Counter[str] = Counter()

    with open_text(input_path, "rt") as src, open_text(output_path, "wt") as dst:
        for line in src:
            if not line.strip():
                continue
            if line.startswith(("#", "track", "browser")):
                dst.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            chrom = parts[0]
            region = regions_map.get(chrom)
            if region is None:
                continue

            try:
                start0 = int(parts[1])
                end0 = int(parts[2])
            except ValueError:
                continue

            clipped = clip_to_region(region, start0, end0)
            if clipped is None:
                continue
            clipped_start, clipped_end = clipped

            clipped_parts = [chrom, str(clipped_start), str(clipped_end), *parts[3:]]
            dst.write("\t".join(clipped_parts) + "\n")

            interval_count += 1
            clipped_bp = clipped_end - clipped_start
            bp_count += clipped_bp
            per_region_bp[region.label] += clipped_bp

    return {
        "input": str(input_path),
        "output": str(output_path),
        "intervals_kept": interval_count,
        "bp_kept": bp_count,
        "per_region_bp": dict(sorted(per_region_bp.items())),
    }


def canonical_exclusion_name(path: Path) -> str:
    return re.sub(r"_[01]$", "", path.stem)


def write_readme(output_root: Path) -> None:
    readme_path = output_root / "README.md"
    ensure_parent(readme_path)
    readme = """# GRCh38 Debug Subset

Compact subset extracted from GRCh38 `v5.0q` benchmark resources for testing/debugging.

## Design Constraints

- One region per chromosome
- Exactly 7 chromosomes: 5 autosomes + `chrX` + `chrY`
- Includes overlap-rich regions spanning multiple exclusion/context tracks
- Includes one zero-overlap control region (`no_overlap_control`)
- Captures variant size bins from 1bp through 1-5kb with zero >5kb variants

## Regeneration

```bash
python scripts/create_grch38_debug_subset.py --force
```
"""
    readme_path.write_text(readme, encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    regions_map = region_by_chrom(REGIONS)

    if args.force and output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    write_regions_bed(output_root / "regions.bed", REGIONS)

    benchmark_dir = output_root / "benchmarksets"
    contexts_dir = output_root / "stratifications"
    exclusions_dir = output_root / "exclusions" / "v5.0q_GRCh38_stvar"

    bed_stats: list[dict[str, object]] = []
    vcf_stats: list[dict[str, object]] = []
    context_stats: list[dict[str, object]] = []
    exclusion_stats: list[dict[str, object]] = []

    for bed_file in INPUT_BEDS:
        output_path = benchmark_dir / bed_file.name
        bed_stats.append(clip_bed(bed_file, output_path, regions_map))

    for vcf_file in INPUT_VCFS:
        output_path = benchmark_dir / vcf_file.name
        vcf_stats.append(filter_vcf(vcf_file, output_path, regions_map))

    for context_file in CONTEXT_FILES:
        output_name = context_file.name
        output_path = contexts_dir / output_name
        context_stats.append(clip_bed(context_file, output_path, regions_map))

    exclusion_files = sorted(EXCLUSIONS_DIR.glob("*.bed"))
    for exclusion_file in exclusion_files:
        output_path = exclusions_dir / exclusion_file.name
        exclusion_stats.append(clip_bed(exclusion_file, output_path, regions_map))

    context_types_overlapping = sorted(
        stat["output"]
        .split("/")[-1]
        .replace(".bed.gz", "")
        .replace(".bed", "")
        .replace("GRCh38_", "")
        for stat in context_stats
        if stat["bp_kept"] > 0
    )
    exclusion_types_overlapping = sorted(
        {
            canonical_exclusion_name(Path(stat["output"]))
            for stat in exclusion_stats
            if stat["bp_kept"] > 0
        }
    )

    no_overlap_context_bp = sum(
        int(stat["per_region_bp"].get(NO_OVERLAP_LABEL, 0)) for stat in context_stats
    )
    no_overlap_exclusion_bp = sum(
        int(stat["per_region_bp"].get(NO_OVERLAP_LABEL, 0)) for stat in exclusion_stats
    )

    stvar_stats = next(
        (
            stat
            for stat in vcf_stats
            if Path(stat["input"]).name == "v5.0q_GRCh38_stvar_benchmark.vcf.gz"
        ),
        None,
    )
    if stvar_stats is None:
        available_inputs = ", ".join(
            Path(s.get("input", "<missing>")).name for s in vcf_stats
        )
        raise RuntimeError(
            "Expected stvar VCF stats for 'v5.0q_GRCh38_stvar_benchmark.vcf.gz' "
            f"not found in vcf_stats. Available inputs: {available_inputs}"
        )

    summary = {
        "dataset_name": "grch38_debug_subset_v5.0q",
        "source": {
            "benchmarksets_dir": str(BENCHMARKSETS_DIR),
            "exclusions_dir": str(EXCLUSIONS_DIR),
            "stratifications_dir": str(STRATIFICATIONS_DIR),
        },
        "regions": [
            {
                "chrom": region.chrom,
                "start": region.start,
                "end": region.end,
                "label": region.label,
                "length_bp": region.length_bp,
            }
            for region in REGIONS
        ],
        "constraint_checks": {
            "one_region_per_chromosome": len({r.chrom for r in REGIONS}) == len(REGIONS),
            "autosome_count": sum(1 for r in REGIONS if r.chrom[3:].isdigit()),
            "includes_chrX": any(r.chrom == "chrX" for r in REGIONS),
            "includes_chrY": any(r.chrom == "chrY" for r in REGIONS),
            "overlapping_context_types": context_types_overlapping,
            "n_overlapping_context_types": len(context_types_overlapping),
            "overlapping_exclusion_types": exclusion_types_overlapping,
            "n_overlapping_exclusion_types": len(exclusion_types_overlapping),
            "no_overlap_region_label": NO_OVERLAP_LABEL,
            "no_overlap_region_context_bp": no_overlap_context_bp,
            "no_overlap_region_exclusion_bp": no_overlap_exclusion_bp,
            "stvar_variant_types": stvar_stats["variant_type_counts"],
            "stvar_size_bins": stvar_stats["size_bin_counts"],
            "stvar_variants_gt_5kb": stvar_stats["variants_gt_5kb"],
        },
        "files": {
            "regions_bed": str(output_root / "regions.bed"),
            "benchmark_beds": [stat["output"] for stat in bed_stats],
            "benchmark_vcfs": [stat["output"] for stat in vcf_stats],
            "context_beds": [stat["output"] for stat in context_stats],
            "exclusion_beds": [stat["output"] for stat in exclusion_stats],
        },
        "stats": {
            "benchmark_beds": bed_stats,
            "benchmark_vcfs": vcf_stats,
            "context_beds": context_stats,
            "exclusion_beds": exclusion_stats,
        },
    }

    summary_path = output_root / "subset_summary.json"
    ensure_parent(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    write_readme(output_root)
    print(f"Wrote GRCh38 debug subset to: {output_root}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
