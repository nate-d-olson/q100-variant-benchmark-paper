#!/usr/bin/env bash
# Calculate total bases in each dip.bed file and overall total.
# AI Disclosure: Script developed with assistance from Claude (Anthropic).
set -euo pipefail

DIP_DIR="resources/benchmarksets"
GRAND_TOTAL=0

printf "%-45s %15s\n" "File" "Total Bases"
printf "%-45s %15s\n" "----" "-----------"

for bed_file in "${DIP_DIR}"/*_dip.bed; do
    if [[ ! -f "${bed_file}" ]]; then
        echo "No dip.bed files found in ${DIP_DIR}" >&2
        exit 1
    fi
    total=$(awk '{sum += $3 - $2} END {print sum + 0}' "${bed_file}")
    GRAND_TOTAL=$((GRAND_TOTAL + total))
    printf "%-45s %'15d\n" "$(basename "${bed_file}")" "${total}"
done

printf "%-45s %15s\n" "" "==============="
printf "%-45s %'15d\n" "TOTAL" "${GRAND_TOTAL}"
