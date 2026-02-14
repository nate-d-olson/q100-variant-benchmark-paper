"""
Rules for downloading and caching remote input files.

This module implements robust download rules following Snakemake 8+ best practices:
- Uses wget with retries for network stability
- Validates checksums (SHA256 or MD5) from config
- Outputs to resources/ directory for consistent local paths
- Proper logging for debugging failed downloads
- Supports local file paths for test fixtures (file:// URLs, absolute/relative paths)

Download types handled:
- Benchmark VCF and BED files
- Reference genome files
- Stratification BED files
- Exclusion BED files

For local fixture testing, use file:// URLs or absolute/relative paths in config.
Example: "file://$PWD/tests/fixtures/grch38_debug_subset/benchmarksets/..."
"""

import os


# ============================================================================
# Benchmark File Downloads (VCF, BED, dip_bed)
# ============================================================================


rule download_benchmark_vcf:
    """
    Download benchmark VCF file with checksum validation.

    Downloads VCF files from URLs specified in config["benchmarksets"][benchmark]["vcf"]
    and validates against SHA256 checksum using Snakemake's ensure() function.
    
    Supports both remote URLs and local file paths (including file:// URLs).
    """
    output:
        vcf=ensure(
            "resources/benchmarksets/{benchmark}_benchmark.vcf.gz",
            non_empty=True,
            sha256=lambda w: config["benchmarksets"][w.benchmark]["vcf"]["sha256"],
        ),
    params:
        url=lambda w: config["benchmarksets"][w.benchmark]["vcf"]["url"],
        is_local=lambda w: is_local_file(
            config["benchmarksets"][w.benchmark]["vcf"]["url"]
        ),
        local_path=lambda w: (
            resolve_file_url(config["benchmarksets"][w.benchmark]["vcf"]["url"])
            if is_local_file(config["benchmarksets"][w.benchmark]["vcf"]["url"])
            else ""
        ),
    log:
        "logs/download_benchmark_vcf/{benchmark}.log",
    retries: 3
    resources:
        mem_mb=512,
        runtime=30,  # 30 minutes
    conda:
        "../envs/downloads.yaml"
    shell:
        """
        exec 2> {log}
        set -euo pipefail

        echo "[$(date)] Downloading benchmark VCF for {wildcards.benchmark}" | tee -a {log}
        echo "URL: {params.url}" | tee -a {log}

        if [ "{params.is_local}" = "True" ]; then
            echo "Local file detected: {params.local_path}" | tee -a {log}
            cp "{params.local_path}" {output.vcf} 2>&1 | tee -a {log}
        else
            echo "Remote URL detected, using wget" | tee -a {log}
            wget --no-verbose -O {output.vcf} "{params.url}" 2>&1 | tee -a {log}
        fi

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """


rule download_benchmark_bed:
    """
    Download benchmark BED file with checksum validation.

    Downloads BED files from URLs specified in config["benchmarksets"][benchmark]["bed"]
    and validates against SHA256 checksum using Snakemake's ensure() function.
    
    Supports both remote URLs and local file paths (including file:// URLs).
    """
    output:
        bed=ensure(
            "resources/benchmarksets/{benchmark}_benchmark.bed",
            non_empty=True,
            sha256=lambda w: config["benchmarksets"][w.benchmark]["bed"]["sha256"],
        ),
    params:
        url=lambda w: config["benchmarksets"][w.benchmark]["bed"]["url"],
        is_local=lambda w: is_local_file(
            config["benchmarksets"][w.benchmark]["bed"]["url"]
        ),
        local_path=lambda w: (
            resolve_file_url(config["benchmarksets"][w.benchmark]["bed"]["url"])
            if is_local_file(config["benchmarksets"][w.benchmark]["bed"]["url"])
            else ""
        ),
    log:
        "logs/download_benchmark_bed/{benchmark}.log",
    retries: 3
    resources:
        mem_mb=512,
        runtime=30,  # 30 minutes
    conda:
        "../envs/downloads.yaml"
    shell:
        """
        exec 2> {log}
        set -euo pipefail

        echo "[$(date)] Downloading benchmark BED for {wildcards.benchmark}" | tee -a {log}
        echo "URL: {params.url}" | tee -a {log}

        if [ "{params.is_local}" = "True" ]; then
            echo "Local file detected: {params.local_path}" | tee -a {log}
            cp "{params.local_path}" {output.bed} 2>&1 | tee -a {log}
        else
            echo "Remote URL detected, using wget" | tee -a {log}
            wget --no-verbose -O {output.bed} "{params.url}" 2>&1 | tee -a {log}
        fi

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """


rule download_benchmark_dip_bed:
    """
    Download dipcall BED file with checksum validation.

    Downloads dip.bed files from URLs specified in
    config["benchmarksets"][benchmark]["dip_bed"] and validates against SHA256 checksum
    using Snakemake's ensure() function.

    Only applies to benchmarks that have dip_bed configured (v5q benchmarks).
    
    Supports both remote URLs and local file paths (including file:// URLs).
    """
    output:
        bed=ensure(
            "resources/benchmarksets/{benchmark}_dip.bed",
            non_empty=True,
        ),
    params:
        url=lambda w: config["benchmarksets"][w.benchmark]
        .get("dip_bed", {})
        .get("url", ""),
        sha256=lambda w: config["benchmarksets"][w.benchmark]
        .get("dip_bed", {})
        .get("sha256", ""),
        is_local=lambda w: is_local_file(
            config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("url", "")
        ),
        local_path=lambda w: (
            resolve_file_url(
                config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("url", "")
            )
            if is_local_file(
                config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("url", "")
            )
            else ""
        ),
    log:
        "logs/download_benchmark_dip_bed/{benchmark}.log",
    retries: 3
    resources:
        mem_mb=512,
        runtime=30,  # 30 minutes
    conda:
        "../envs/downloads.yaml"
    shell:
        """
        exec 2> {log}
        set -euo pipefail

        echo "[$(date)] Downloading dip.bed for {wildcards.benchmark}" | tee -a {log}
        echo "URL: {params.url}" | tee -a {log}
        echo "Expected SHA256: {params.sha256}" | tee -a {log}

        if [ "{params.is_local}" = "True" ]; then
            echo "Local file detected: {params.local_path}" | tee -a {log}
            # Copy to temporary location first
            cp "{params.local_path}" {output.bed}.tmp 2>&1 | tee -a {log}
        else
            echo "Remote URL detected, using wget" | tee -a {log}
            # Download to temporary location
            wget --no-verbose -O {output.bed}.tmp "{params.url}" 2>&1 | tee -a {log}
        fi

        # Validate SHA256 checksum
        echo "[$(date)] Validating checksum..." | tee -a {log}
        echo "{params.sha256}  {output.bed}.tmp" | sha256sum -c - 2>&1 | tee -a {log}

        # Move to final location after validation
        mv {output.bed}.tmp {output.bed}

        echo "[$(date)] Download and validation completed successfully" | tee -a {log}
        """


# ============================================================================
# Reference Genome Downloads and Preparation
# ============================================================================


rule prepare_reference:
    """
    Download and prepare reference genome for bcftools.

    Downloads reference from URL, validates checksum, converts to bgzip format,
    and creates faidx and gzi indices for random access.

    Note: Checksum validation is done in shell (not via ensure()) because
    the checksum in config is for the original downloaded gzip file, not the
    final bgzip-converted output. The file is validated before conversion.
    
    Supports both remote URLs and local file paths (including file:// URLs).
    """
    output:
        ref=ensure("resources/references/{ref_name}.fa.gz", non_empty=True),
        fai=ensure("resources/references/{ref_name}.fa.gz.fai", non_empty=True),
        gzi=ensure("resources/references/{ref_name}.fa.gz.gzi", non_empty=True),
    params:
        url=lambda w: config["references"][w.ref_name]["url"],
        checksum=lambda w: get_reference_checksum(w.ref_name),
        checksum_type="sha256",  ## hardcoded using sha256 in config
        is_local=lambda w: is_local_file(config["references"][w.ref_name]["url"]),
        local_path=lambda w: (
            resolve_file_url(config["references"][w.ref_name]["url"])
            if is_local_file(config["references"][w.ref_name]["url"])
            else ""
        ),
    log:
        "logs/prepare_reference/{ref_name}.log",
    retries: 3
    resources:
        mem_mb=2048,
        runtime=30,  # 30 minutes
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        exec 2> {log}
        set -euo pipefail

        echo "[$(date)] Downloading and preparing reference {wildcards.ref_name}" | tee -a {log}
        echo "URL: {params.url}" | tee -a {log}

        # Create temp file for download
        TMPFILE=$(mktemp)
        trap "rm -f $TMPFILE" EXIT

        if [ "{params.is_local}" = "True" ]; then
            echo "Local file detected: {params.local_path}" | tee -a {log}
            cp "{params.local_path}" "$TMPFILE" 2>&1 | tee -a {log}
        else
            echo "Remote URL detected, using wget" | tee -a {log}
            wget --no-verbose -O "$TMPFILE" "{params.url}" 2>&1 | tee -a {log}
        fi

        # Validate checksum
        echo "[$(date)] Validating {params.checksum_type} checksum..." | tee -a {log}
        if [ "{params.checksum_type}" = "md5" ]; then
            echo "{params.checksum}  $TMPFILE" | md5sum -c - 2>&1 | tee -a {log}
        else
            echo "{params.checksum}  $TMPFILE" | sha256sum -c - 2>&1 | tee -a {log}
        fi

        # Convert gzip to bgzip
        echo "[$(date)] Converting to bgzip format..." | tee -a {log}
        gunzip -c "$TMPFILE" | bgzip -c > {output.ref}

        # Index the reference
        echo "[$(date)] Creating faidx index..." | tee -a {log}
        samtools faidx {output.ref}

        # Verify all outputs were created
        for f in {output.ref} {output.fai} {output.gzi}; do
            if [ ! -s "$f" ]; then
                echo "ERROR: Output file $f is missing or empty" | tee -a {log}
                exit 1
            fi
        done

        echo "[$(date)] Reference preparation completed successfully" | tee -a {log}
        """


# ============================================================================
# Stratification File Downloads
# ============================================================================


rule download_stratification:
    """
    Download stratification BED file with URL templating.

    Stratification URLs in config use {ref} template which gets replaced
    based on the reference genome (GRCh37, GRCh38, CHM13).

    Note: GIAB stratifications URL format uses {ref}@all pattern which needs
    special handling.
    
    Supports both remote URLs and local file paths (including file:// URLs).
    """
    output:
        bed=ensure(
            "resources/stratifications/{ref}_{strat_name}.bed.gz",
            sha256=lambda w: get_stratification_sha256(w),
            non_empty=True,
        ),
    params:
        url=lambda w: get_stratification_url(w),
        is_local=lambda w: is_local_file(get_stratification_url(w)),
        local_path=lambda w: (
            resolve_file_url(get_stratification_url(w))
            if is_local_file(get_stratification_url(w))
            else ""
        ),
    log:
        "logs/download_stratification/{ref}_{strat_name}.log",
    wildcard_constraints:
        ref="GRCh37|GRCh38|CHM13v2.0",
        strat_name="TR|TR10kb|HP|SD|SD10kb|MAP",
    retries: 3
    resources:
        mem_mb=256,
        runtime=30,  # 30 minutes
    conda:
        "../envs/downloads.yaml"
    shell:
        """
        exec 2> {log}
        set -euo pipefail

        echo "[$(date)] Downloading stratification {wildcards.strat_name} for {wildcards.ref}" | tee -a {log}
        echo "URL: {params.url}" | tee -a {log}

        if [ "{params.is_local}" = "True" ]; then
            echo "Local file detected: {params.local_path}" | tee -a {log}
            cp "{params.local_path}" {output.bed} 2>&1 | tee -a {log}
        else
            echo "Remote URL detected, using wget" | tee -a {log}
            wget --no-verbose -O {output.bed} "{params.url}" 2>&1 | tee -a {log}
        fi

        # Verify non-empty
        if [ ! -s {output.bed} ]; then
            echo "ERROR: Downloaded file is empty" | tee -a {log}
            exit 1
        fi

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """


# ============================================================================
# Exclusion File Downloads
# ============================================================================


rule download_exclusion:
    """
    Download exclusion BED file with checksum validation.

    Downloads exclusion BED files specified in
    config["benchmarksets"][benchmark]["exclusions"][exclusion_name]["files"]
    and validates against SHA256 checksum using Snakemake's ensure() function.
    
    Supports both remote URLs and local file paths (including file:// URLs).
    """
    output:
        bed=ensure(
            "resources/exclusions/{benchmark}/{exclusion_name}_{file_idx}.bed",
            non_empty=True,
            sha256=lambda w: get_exclusion_file_checksum(
                w.benchmark, w.exclusion_name, int(w.file_idx)
            ),
        ),
    params:
        url=lambda w: get_exclusion_file_url(
            w.benchmark, w.exclusion_name, int(w.file_idx)
        ),
        is_local=lambda w: is_local_file(
            get_exclusion_file_url(w.benchmark, w.exclusion_name, int(w.file_idx))
        ),
        local_path=lambda w: (
            resolve_file_url(
                get_exclusion_file_url(w.benchmark, w.exclusion_name, int(w.file_idx))
            )
            if is_local_file(
                get_exclusion_file_url(w.benchmark, w.exclusion_name, int(w.file_idx))
            )
            else ""
        ),
    log:
        "logs/download_exclusion/{benchmark}_{exclusion_name}_{file_idx}.log",
    retries: 3
    resources:
        mem_mb=256,
        runtime=30,  # 30 minutes
    conda:
        "../envs/downloads.yaml"
    shell:
        """
        exec 2> {log}
        set -euo pipefail

        echo "[$(date)] Downloading exclusion {wildcards.exclusion_name} file {wildcards.file_idx} for {wildcards.benchmark}" | tee -a {log}
        echo "URL: {params.url}" | tee -a {log}

        if [ "{params.is_local}" = "True" ]; then
            echo "Local file detected: {params.local_path}" | tee -a {log}
            cp "{params.local_path}" {output.bed} 2>&1 | tee -a {log}
        else
            echo "Remote URL detected, using wget" | tee -a {log}
            wget --no-verbose -O {output.bed} "{params.url}" 2>&1 | tee -a {log}
        fi

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """
