#!/usr/bin/env bash
# Run SV use-case benchmarking with GIAB v5.0q recommended truvari parameters.
# Requires: truvari conda env (workflow/envs/truvari.yaml)
# AI Disclosure: Developed with assistance from Claude Opus 4.6 (Anthropic).
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${PROJ_DIR}/results/use_case/stvar"
CALLSET_DIR="${PROJ_DIR}/results/use_case/callsets"
BENCH_VCF="${PROJ_DIR}/resources/benchmarksets/v5.0q_GRCh38_stvar_benchmark.vcf.gz"
BENCH_BED="${PROJ_DIR}/resources/benchmarksets/v5.0q_GRCh38_stvar_benchmark.bed"
REF="${PROJ_DIR}/resources/references/GRCh38.fa.gz"
STRAT_DIR="${PROJ_DIR}/resources/stratifications"
THREADS=4

CALLSET_NAMES="ont-sniffles ont-verkko"
SNIFFLES_SRC="/Users/nolson/Google Drive/Shared drives/BBD_Human_Genomics/HG002-Q100v1.1-Variant-Benchmark-Evaluations/evaluations/stvar/callsets/ont-sniffles/v5hac40x_snifflesPhased_v2.5.3tr.pass.vcf.gz"
VERKKO_SRC="/Users/nolson/Google Drive/Shared drives/BBD_Human_Genomics/HG002-Q100v1.1-Variant-Benchmark-Evaluations/evaluations/stvar/callsets/ont-verkko/twoFCverkko.dipcall.dip.svs.splitmulti.vcf.gz"
SNIFFLES_VCF="${CALLSET_DIR}/ont-sniffles.vcf.gz"
VERKKO_VCF="${CALLSET_DIR}/ont-verkko.vcf.gz"

for callset in ${CALLSET_NAMES}; do
    case "${callset}" in
        ont-sniffles) src_vcf="${SNIFFLES_SRC}"; vcf="${SNIFFLES_VCF}" ;;
        ont-verkko)   src_vcf="${VERKKO_SRC}"; vcf="${VERKKO_VCF}" ;;
    esac
    run_dir="${OUTDIR}/${callset}_v5.0q"
    filtered_vcf="${OUTDIR}/${callset}_filtered.vcf.gz"

    echo "=== Processing ${callset} ==="

    # Step 0: Copy callset from Google Drive if not present
    if [ ! -f "${vcf}" ]; then
        echo "Copying ${callset} from Google Drive..."
        cp "${src_vcf}" "${vcf}"
        cp "${src_vcf}.tbi" "${vcf}.tbi"
    else
        echo "Callset ${vcf} already exists, skipping copy"
    fi

    # Step 1: Filter ALT="*" variants (per GIAB README)
    echo "Filtering ALT=* variants..."
    bcftools view --exclude 'ALT="*"' --output-type z --output "${filtered_vcf}" "${vcf}"
    bcftools index --tbi "${filtered_vcf}"

    # Step 2: Remove previous output (truvari bench refuses to overwrite)
    rm -rf "${run_dir}"

    # Step 3: Run truvari bench with GIAB-recommended parameters
    echo "Running truvari bench..."
    truvari bench \
        -b "${BENCH_VCF}" \
        -c "${filtered_vcf}" \
        -o "${run_dir}" \
        -f "${REF}" \
        --includebed "${BENCH_BED}" \
        --pick ac \
        --passonly \
        -r 2000 \
        -C 5000

    # Step 4: Run truvari refine (mirrors benchmark_comparisons.smk)
    echo "Running truvari refine..."
    truvari refine \
        --reference "${REF}" \
        --use-original-vcfs \
        --threads "${THREADS}" \
        "${run_dir}"

    # Step 5: Run truvari stratify for each genomic context
    # truvari stratify takes positional args: BED BENCH_DIR [-o OUTPUT]
    echo "Running truvari stratify..."
    for context in HP TR SD MAP; do
        echo "  Stratifying ${context}..."
        truvari stratify \
            "${STRAT_DIR}/GRCh38_${context}.bed.gz" \
            "${run_dir}" \
            -o "${run_dir}/stratify_${context}.txt"
    done

    # Cleanup filtered VCF
    rm -f "${filtered_vcf}" "${filtered_vcf}.tbi"

    echo "=== Done ${callset} ==="
done

echo "All SV use-case benchmarks complete."
echo "Run scripts/extract_sv_metrics.py to generate CSV metrics."
