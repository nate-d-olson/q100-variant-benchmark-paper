#!/bin/bash

# Script to check S3 object status and restore from Glacier if needed
# Based on download failures in logs/downloads/exclusions/

BUCKET="giab-data"
PROFILE="commercial"
RESTORE_TIER="Standard"  # Options: Standard (3-5 hrs), Bulk (5-12 hrs)

# List of S3 keys that failed to download (403 Forbidden errors)
KEYS=(
    "defrabb_runs/20250117_v0.020_HG002Q100v1.1/resources/exclusions/CHM13v2.0/HG002Q100-errors_sorted.bed"
    "defrabb_runs/20250117_v0.020_HG002Q100v1.1/resources/exclusions/GRCh37/HG002Q100-errors_sorted.bed"
    "defrabb_runs/20250117_v0.020_HG002Q100v1.1/resources/exclusions/GRCh38/HG002Q100-errors_sorted.bed"
    "defrabb_runs/20250117_v0.020_HG002Q100v1.1/results/draft_benchmarksets/CHM13_HG002-T2TQ100v1.1-dipz2k_smvar-excluded/exclusions/intersections/CHM13v2.0_HG2-T2TQ100-V1.1_smvar_dipcall-z2k_svs-and-simple-repeats_sorted.bed"
    "defrabb_runs/20250117_v0.020_HG002Q100v1.1/results/draft_benchmarksets/GRCh37_HG002-T2TQ100v1.1-dipz2k_smvar-excluded/exclusions/intersections/GRCh37_HG2-T2TQ100-V1.1_smvar_dipcall-z2k_svs-and-simple-repeats_sorted.bed"
    "defrabb_runs/20250117_v0.020_HG002Q100v1.1/results/draft_benchmarksets/GRCh38_HG002-T2TQ100v1.1-dipz2k_smvar-excluded/exclusions/intersections/GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z2k_svs-and-simple-repeats_sorted.bed"
)

echo "========================================="
echo "S3 Glacier Restore Script"
echo "Bucket: ${BUCKET}"
echo "Profile: ${PROFILE}"
echo "Restore Tier: ${RESTORE_TIER} (${RESTORE_DAYS} days)"
echo "========================================="
echo ""

# Function to check object status and restore if needed
check_and_restore() {
    local key="$1"

    echo "Checking: ${key}"

    # Get object metadata
    local metadata
    if ! metadata=$(aws s3api head-object \
        --bucket "${BUCKET}" \
        --key "${key}" \
        --profile "${PROFILE}" 2>&1); then
        echo "  ✗ ERROR: Cannot access object"
        echo "    ${metadata}"
        echo ""
        return 1
    fi

    # Extract storage class and archive status
    local storage_class
    local archive_status
    local restore_status

    storage_class=$(echo "${metadata}" | grep -o '"StorageClass": "[^"]*"' | cut -d'"' -f4)
    archive_status=$(echo "${metadata}" | grep -o '"ArchiveStatus": "[^"]*"' | cut -d'"' -f4)
    restore_status=$(echo "${metadata}" | grep -o '"Restore": "[^"]*"' | cut -d'"' -f4)

    echo "  Storage Class: ${storage_class:-STANDARD}"
    echo "  Archive Status: ${archive_status:-N/A}"

    # Check if restore is already in progress or completed
    if [[ -n "${restore_status}" ]]; then
        if echo "${restore_status}" | grep -q "ongoing-request=\"true\""; then
            echo "  ✓ Restore already in progress"
        elif echo "${restore_status}" | grep -q "ongoing-request=\"false\""; then
            echo "  ✓ Restore completed - object is accessible"
            # Extract expiry date if available
            expiry=$(echo "${restore_status}" | grep -o 'expiry-date="[^"]*"' | cut -d'"' -f2)
            if [[ -n "${expiry}" ]]; then
                echo "    Available until: ${expiry}"
            fi
        fi
        echo ""
        return 0
    fi

    # Check if object needs restoration
    if [[ "${archive_status}" == "ARCHIVE_ACCESS" ]] || [[ "${archive_status}" == "DEEP_ARCHIVE_ACCESS" ]]; then
        echo "  ⚠ Object is archived (${archive_status})"
        echo "  → Initiating restore request..."

        # Choose appropriate tier based on archive status
        local tier="${RESTORE_TIER}"
        if [[ "${archive_status}" == "ARCHIVE_ACCESS" ]]; then
            # ARCHIVE_ACCESS supports Expedited retrieval (1-5 minutes)
            tier="Expedited"
        fi
        # DEEP_ARCHIVE_ACCESS only supports Standard (12 hrs) or Bulk (48 hrs)

        # Initiate restore
        if aws s3api restore-object \
            --bucket "${BUCKET}" \
            --key "${key}" \
            --restore-request "{\"GlacierJobParameters\":{\"Tier\":\"${tier}\"}}" \
            --profile "${PROFILE}" 2>&1; then
            echo "  ✓ Restore request submitted successfully (${tier} tier)"

            if [[ "${archive_status}" == "DEEP_ARCHIVE_ACCESS" ]]; then
                echo "    Estimated time: 12 hours (Standard tier)"
            else
                echo "    Estimated time: 1-5 minutes (Expedited tier)"
            fi
        else
            echo "  ✗ Failed to submit restore request"
        fi
    else
        echo "  ✓ Object is accessible (not archived)"
    fi

    echo ""
}

# Process each key
for key in "${KEYS[@]}"; do
    check_and_restore "${key}"
done

echo "========================================="
echo "Summary:"
echo "  Total objects checked: ${#KEYS[@]}"
echo ""
echo "To check restore status later, run:"
echo "  aws s3api head-object --bucket ${BUCKET} --key <key> --profile ${PROFILE}"
echo ""
echo "After restoration completes, retry your downloads."
echo "========================================="
