#!/usr/bin/env bash
# Build per-reference annotation BEDs (benchmark regions + "large excluded
# regions") for the SVbyEye same-scale figure.
#
# Reconstruction note: reimplements a script lost from a deleted scratch
# worktree (scratch/ideogram-explore/scripts/prep_svbyeye_beds.sh, 2026-06-17
# session); see scripts/make_svbyeye_samescale.R header for full context,
# including why these 5 categories (and not the other 7 exclusion
# categories) make up the red track.
#
# AI Disclosure: Developed with assistance from Claude (Anthropic).
#
# Usage: bash scripts/prep_svbyeye_beds.sh <REF>
#   REF: GRCh38 | GRCh37 | CHM13v2.0
#
# Outputs (plain bedtools merge; the >=10kb size filter on the exclusion
# track is applied at plot time in R, not here):
#   svbyeye/<REF>/v5_benchmark_all.bed  - union of v5.0q smvar + stvar regions
#   svbyeye/<REF>/excl_large_all.bed    - union of {segdups, satellites,
#                                          tandem-repeats, flanks, gaps}

set -euo pipefail

REF="${1:?Usage: prep_svbyeye_beds.sh <REF>}"
DATA_ROOT="${Q100_DATA_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BMK_DIR="$DATA_ROOT/resources/benchmarksets"
EXCL_DIR="$DATA_ROOT/resources/exclusions/v5.0q_${REF}_stvar"

OUT_DIR="results/svbyeye/${REF}"
mkdir -p "$OUT_DIR"

# GRCh37 benchmark/exclusion BEDs use bare contig names (1, 2, ..., X); add
# the chr prefix so downstream R code can assume canonical chr-prefixed names
# regardless of reference (the PAF target column keeps the bare name --
# handled separately in make_svbyeye_samescale.R's ref_paf_name()).
add_chr_prefix() {
  if [[ "$REF" == "GRCh37" ]]; then
    sed -E 's/^([0-9XYM]+)\t/chr\1\t/'
  else
    cat
  fi
}

echo "[prep_svbyeye_beds] building v5_benchmark_all.bed ($REF)..."
cat "$BMK_DIR/v5.0q_${REF}_smvar_benchmark.bed" "$BMK_DIR/v5.0q_${REF}_stvar_benchmark.bed" \
  | add_chr_prefix \
  | sort -k1,1 -k2,2n \
  | bedtools merge -d 10000 -i - \
  >"$OUT_DIR/v5_benchmark_all.bed"
echo "[prep_svbyeye_beds]   -> $(wc -l <"$OUT_DIR/v5_benchmark_all.bed" | tr -d ' ') intervals"

echo "[prep_svbyeye_beds] building excl_large_all.bed ($REF)..."
excl_files=()
for cat in segdups satellites tandem-repeats flanks gaps; do
  for f in "$EXCL_DIR"/${cat}_*.bed; do
    [[ -e "$f" ]] && excl_files+=("$f")
  done
done
if [[ ${#excl_files[@]} -eq 0 ]]; then
  echo "[prep_svbyeye_beds] ERROR: no exclusion BEDs found under $EXCL_DIR for segdups/satellites/tandem-repeats/flanks/gaps" >&2
  exit 1
fi
cat "${excl_files[@]}" \
  | add_chr_prefix \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  >"$OUT_DIR/excl_large_all.bed"
echo "[prep_svbyeye_beds]   -> $(wc -l <"$OUT_DIR/excl_large_all.bed" | tr -d ' ') intervals (>=10kb filter applied at plot time)"
