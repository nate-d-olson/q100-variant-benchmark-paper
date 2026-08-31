#!/usr/bin/env bash
# Generate per-chromosome minimap2 asm5 PAFs (REF vs HG002 maternal/paternal)
# for the SVbyEye same-scale figure.
#
# Reconstruction note: reimplements a script lost from a deleted scratch
# worktree (scratch/ideogram-explore/scripts/prep_svbyeye_pafs.sh, 2026-06-17
# session); see scripts/make_svbyeye_samescale.R header for full context.
#
# AI Disclosure: Developed with assistance from Claude (Anthropic).
#
# Usage: bash scripts/prep_svbyeye_pafs.sh <REF> <chrom1,chrom2,...>
#   REF: GRCh38 | GRCh37 | CHM13v2.0
#   chroms: comma-separated, e.g. chr6,chr8,chr15,chr20,chrX (default: all
#           autosomes + chrX + chrY)
#
# Requires minimap2 + samtools on PATH (e.g. `mamba activate hlienv`).
# chrX has no paternal-origin homolog in a male sample (HG002 is XY); the
# PAT alignment for chrX is skipped rather than fabricated.

set -euo pipefail

REF="${1:?Usage: prep_svbyeye_pafs.sh <REF> [chroms]}"
CHROMS="${2:-$(seq 1 22 | sed 's/^/chr/' | paste -sd, -),chrX,chrY}"

DATA_ROOT="${Q100_DATA_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
REF_FASTA="$DATA_ROOT/resources/references/${REF}.fa.gz"
MAT_FASTA="$DATA_ROOT/resources/references/HG002_mat.fa.gz"
PAT_FASTA="$DATA_ROOT/resources/references/HG002_pat.fa.gz"

OUT_ROOT="results/svbyeye/${REF}"
mkdir -p "$OUT_ROOT"

IFS=',' read -ra CHROM_LIST <<<"$CHROMS"

for chrom in "${CHROM_LIST[@]}"; do
  outdir="$OUT_ROOT/$chrom"
  mkdir -p "$outdir"

  # GRCh37 FASTA uses bare contig names (1, 2, ..., X) despite chr-prefixed
  # benchmark BEDs; strip the prefix only for the reference extraction.
  if [[ "$REF" == "GRCh37" ]]; then
    ref_contig="${chrom#chr}"
  else
    ref_contig="$chrom"
  fi

  ref_fa="$outdir/ref.fa"
  if [[ ! -s "$ref_fa" ]]; then
    samtools faidx "$REF_FASTA" "$ref_contig" | sed "s/^>.*/>${chrom}/" >"$ref_fa"
  fi

  for hap in mat pat; do
    hap_long="MATERNAL"
    [[ "$hap" == "pat" ]] && hap_long="PATERNAL"
    hap_contig="${chrom}_${hap_long}"
    hap_fasta="$MAT_FASTA"
    [[ "$hap" == "pat" ]] && hap_fasta="$PAT_FASTA"

    paf_out="$outdir/ref_${hap}.named.paf"
    if [[ "$chrom" == "chrX" && "$hap" == "pat" ]]; then
      echo "[prep_svbyeye_pafs] chrX has no paternal homolog in a male sample -- skipping ref_pat for chrX"
      continue
    fi
    if [[ "$chrom" == "chrY" && "$hap" == "mat" ]]; then
      echo "[prep_svbyeye_pafs] chrY has no maternal homolog in a male sample -- skipping ref_mat for chrY"
      continue
    fi
    if [[ -s "$paf_out" ]]; then
      echo "[prep_svbyeye_pafs] $paf_out already exists, skipping"
      continue
    fi

    hap_fa="$outdir/${hap}.fa"
    samtools faidx "$hap_fasta" "$hap_contig" >"$hap_fa"

    echo "[prep_svbyeye_pafs] aligning $chrom ref vs $hap ($REF)..."
    minimap2 -cx asm5 --eqx -t "${THREADS:-8}" "$ref_fa" "$hap_fa" 2>"$outdir/${hap}.minimap2.log" \
      | awk -v h="HG002_$(echo "$hap" | tr '[:lower:]' '[:upper:]')" 'BEGIN{OFS="\t"} {$1=h; print}' \
        >"$paf_out"
    rm -f "$hap_fa"
    echo "[prep_svbyeye_pafs]   -> $paf_out ($(wc -l <"$paf_out" | tr -d ' ') alignment records)"
  done
done
