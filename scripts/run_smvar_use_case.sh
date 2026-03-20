#!/usr/bin/env bash
# Run small variant use-case benchmarking with rtg vcfeval.
# Produces stratified ROC files using --roc-regions and --roc-subset.
# Requires: rtg-tools (available in hlienv conda env)
# AI Disclosure: Developed with assistance from Claude Opus 4.6 (Anthropic).
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${PROJ_DIR}/results/use_case/smvar"
CALLSET="${PROJ_DIR}/results/use_case/callsets/RN_HG002_Illumina_PacBio_Oxford.vcf.gz"
SDF="${PROJ_DIR}/results/use_case/GRCh38.sdf"
STRAT_DIR="${PROJ_DIR}/resources/stratifications"
THREADS=8

BENCH_VERSIONS="v4.2.1 v5.0q"

for bench in ${BENCH_VERSIONS}; do
    bench_vcf="${PROJ_DIR}/resources/benchmarksets/${bench}_GRCh38_smvar_benchmark.vcf.gz"
    bench_bed="${PROJ_DIR}/resources/benchmarksets/${bench}_GRCh38_smvar_benchmark.bed"
    run_dir="${OUTDIR}/roche_${bench}"

    echo "=== Processing ${bench} ==="

    # Remove previous output (vcfeval refuses to overwrite)
    if [ -d "${run_dir}" ]; then
        echo "Removing previous output: ${run_dir}"
        rm -rf "${run_dir}"
    fi

    rtg vcfeval \
        --calls "${CALLSET}" \
        --baseline "${bench_vcf}" \
        --bed-regions "${bench_bed}" \
        --template "${SDF}" \
        --output "${run_dir}" \
        --vcf-score-field=QUAL \
        --threads "${THREADS}" \
        --roc-regions "HP=${STRAT_DIR}/GRCh38_HP.bed.gz" \
        --roc-regions "TR=${STRAT_DIR}/GRCh38_TR.bed.gz" \
        --roc-regions "SD=${STRAT_DIR}/GRCh38_SD.bed.gz" \
        --roc-regions "MAP=${STRAT_DIR}/GRCh38_MAP.bed.gz" \
        --roc-subset snp,non-snp

    echo "=== Done ${bench} ==="
done

echo "All smvar vcfeval runs complete."
echo "Run scripts/extract_smvar_metrics.py to generate CSV metrics."
