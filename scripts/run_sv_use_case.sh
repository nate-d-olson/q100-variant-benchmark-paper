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
THREADS=16

CALLSET_NAMES="baylor_ont dragen"
BAYLOR_ONT_VCF="${CALLSET_DIR}/hg002_ONT_BCM_GRCh38.vcf.gz"
DRAGEN_VCF="${CALLSET_DIR}/DRAGEN_4.4.x_NIST_HG002.sv.vcf.gz"

for callset in ${CALLSET_NAMES}; do
    case "${callset}" in
        baylor_ont) vcf="${BAYLOR_ONT_VCF}" ;;
        dragen)     vcf="${DRAGEN_VCF}" ;;
    esac
    run_dir="${OUTDIR}/${callset}_v5.0q"
    filtered_vcf="${OUTDIR}/${callset}_filtered.vcf.gz"

    echo "=== Processing ${callset} ==="

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
        --debug \
        "${run_dir}"

    # Cleanup filtered VCF
    rm -f "${filtered_vcf}" "${filtered_vcf}.tbi"

    echo "=== Done ${callset} ==="
done

echo "All SV use-case benchmarks complete."
echo "Run scripts/extract_sv_metrics.py to generate CSV metrics."
