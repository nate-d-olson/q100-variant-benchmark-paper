"""
Rules for downloading and caching remote input files.

This module implements robust download rules following Snakemake 8+ best practices:
- Uses wget with retries for network stability
- Validates checksums (SHA256 or MD5) from config
- Outputs to resources/ directory for consistent local paths
- Proper logging for debugging failed downloads

Download types handled:
- Benchmark VCF and BED files
- Reference genome files
- Stratification BED files
- Exclusion BED files
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
    """
    output:
        vcf=ensure(
            "resources/benchmarksets/{benchmark}_benchmark.vcf.gz",
            non_empty=True,
            sha256=lambda w: config["benchmarksets"][w.benchmark]["vcf"]["sha256"],
        ),
    params:
        url=lambda w: config["benchmarksets"][w.benchmark]["vcf"]["url"],
    log:
        "logs/downloads/{benchmark}_vcf.log",
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

        # Download with wget (checksum validated by Snakemake ensure())
        wget --no-verbose -O {output.vcf} "{params.url}" 2>&1 | tee -a {log}

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """


rule download_benchmark_bed:
    """
    Download benchmark BED file with checksum validation.

    Downloads BED files from URLs specified in config["benchmarksets"][benchmark]["bed"]
    and validates against SHA256 checksum using Snakemake's ensure() function.
    """
    output:
        bed=ensure(
            "resources/benchmarksets/{benchmark}_benchmark.bed",
            non_empty=True,
            sha256=lambda w: config["benchmarksets"][w.benchmark]["bed"]["sha256"],
        ),
    params:
        url=lambda w: config["benchmarksets"][w.benchmark]["bed"]["url"],
    log:
        "logs/downloads/{benchmark}_bed.log",
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

        # Download with wget (checksum validated by Snakemake ensure())
        wget --no-verbose -O {output.bed} "{params.url}" 2>&1 | tee -a {log}

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """


rule download_benchmark_dip_bed:
    """
    Download dipcall BED file with checksum validation.

    Downloads dip.bed files from URLs specified in
    config["benchmarksets"][benchmark]["dip_bed"] and validates against SHA256 checksum
    using Snakemake's ensure() function.

    Only applies to benchmarks that have dip_bed configured (v5q benchmarks).
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
    log:
        "logs/downloads/{benchmark}_dip_bed.log",
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

        # Download to temporary location
        wget --no-verbose -O {output.bed}.tmp "{params.url}" 2>&1 | tee -a {log}

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
    """
    output:
        ref=ensure("resources/references/{ref_name}.fa.gz", non_empty=True),
        fai=ensure("resources/references/{ref_name}.fa.gz.fai", non_empty=True),
        gzi=ensure("resources/references/{ref_name}.fa.gz.gzi", non_empty=True),
    params:
        url=lambda w: config["references"][w.ref_name]["url"],
        checksum=lambda w: get_reference_checksum(w.ref_name),
        checksum_type="sha256",  ## hardcoded using sha256 in config
    log:
        "logs/references/{ref_name}_prepare.log",
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

        # Download with wget
        wget --no-verbose -O "$TMPFILE" "{params.url}" 2>&1 | tee -a {log}

        # Validate checksum if configured (skipped when no checksum is provided)
        if [ -n "{params.checksum}" ]; then
            echo "[$(date)] Validating {params.checksum_type} checksum..." | tee -a {log}
            if [ "{params.checksum_type}" = "md5" ]; then
                echo "{params.checksum}  $TMPFILE" | md5sum -c - 2>&1 | tee -a {log}
            else
                echo "{params.checksum}  $TMPFILE" | sha256sum -c - 2>&1 | tee -a {log}
            fi
        else
            echo "[$(date)] No checksum configured; skipping validation" | tee -a {log}
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
    """
    output:
        bed=ensure(
            "resources/stratifications/{ref}_{strat_name}.bed.gz",
            sha256=lambda w: get_stratification_sha256(w),
            non_empty=True,
        ),
    params:
        url=lambda w: get_stratification_url(w),
    log:
        "logs/downloads/stratifications/{ref}_{strat_name}.log",
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

        # Download with wget
        wget --no-verbose -O {output.bed} "{params.url}" 2>&1 | tee -a {log}

        # Verify non-empty (stratifications don't have checksums in config)
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
    log:
        "logs/downloads/exclusions/{benchmark}_{exclusion_name}_{file_idx}.log",
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

        # Download with wget (checksum validated by Snakemake ensure())
        wget --no-verbose -O {output.bed} "{params.url}" 2>&1 | tee -a {log}

        echo "[$(date)] Download completed successfully" | tee -a {log}
        """
